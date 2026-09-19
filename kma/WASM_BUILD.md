# Building KMA for WebAssembly (OpenPathogen)

This documents how the `kma.js` / `kma.wasm` artifacts in this directory are produced
from upstream [KMA](https://bitbucket.org/genomicepidemiology/kma) **v1.5.1**. KMA is a
native, POSIX, multi-threaded C99 program. Running it unmodified in a browser is
impossible because it relies on `fork()`, `pipe()`, System V shared memory
(`shmget`/`shmat`), and `mmap`-style file access. The WASM port replaces each of those
with a browser-safe equivalent, guarded behind `#ifdef __EMSCRIPTEN__` so the same tree
still builds natively.

The build is **wasm64 (Memory64) + pthreads**:

- **Memory64** lifts the wasm32 4 GB linear-memory ceiling. That ceiling — not the
  database index — is what limited sample size: reads, the gzip k-mer-anker
  intermediate, and outputs all live in MEMFS (= linear memory), so large metagenomes
  used to run out of address space well before anything else.
- **pthreads** re-enable KMA's native within-stage parallelism (`-t N`), which scales
  near-linearly with core count.

Browsers without Memory64 (Safari, and every iOS browser) get a **wasm32 + pthreads**
build of the same sources, `kma32.js` / `kma32.wasm`; see §3f.

## 1. Toolchain

| Tool | Version | Notes |
|------|---------|-------|
| Emscripten (`emcc`/`emar`) | latest (3.1.5x+ for solid Memory64 + pthreads) | Install via [emsdk](https://emscripten.org/docs/getting_started/downloads.html) |
| zlib | provided by `-sUSE_ZLIB=1` | Emscripten ports build it |

```bash
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk && ./emsdk install latest && ./emsdk activate latest
source ./emsdk_env.sh        # puts emcc/emar on PATH
emcc --version               # sanity check
```

## 2. Build

From this `kma/` directory, with `emsdk_env.sh` sourced:

```bash
emmake make -f Makefile.wasm clean      # REQUIRED when flags change — see note
emmake make -f Makefile.wasm
```

This compiles every `.c` to a `.o`, archives them into `libkma.a`, then links
`wasm_wrapper.c` + `main.c` against that archive, emitting `kma.js` + `kma.wasm`.

> **Always `clean` first if you change build flags.** Memory64 and `-pthread` change
> code generation (pointer width, atomics), so every object must be rebuilt with the
> same flags. Make only rebuilds `.o`s when the `.c` is newer, so stale objects would
> otherwise be linked silently.

Modern Emscripten **inlines the pthread worker** into `kma.js`; there is no separate
`kma.worker.js` artifact to ship.

### Key compile/link flags (`Makefile.wasm`)

Flags are split into `COMMON` (must match on every compile **and** the link, because
they affect codegen/target/threading) and link-only flags.

| Flag | Why |
|------|-----|
| `-O3 -std=c99 -include stdint.h` | Same optimisation as native; force fixed-width int types |
| `-pthread` | POSIX threads via Web Workers + `SharedArrayBuffer`. **Compile + link.** |
| `-sMEMORY64=1` | wasm64 linear memory (>4 GB). Changes pointer width → **compile + link.** Implies `-sWASM_BIGINT`. |
| `-sUSE_ZLIB=1` | KMA reads/writes gzip |
| `-sALLOW_MEMORY_GROWTH=1` | Linear memory grows on demand |
| `-sINITIAL_MEMORY=67108864` | Start at 64 MB |
| `-sMAXIMUM_MEMORY=17179869184` | Cap at 16 GB (tunable; only meaningful with Memory64) |
| `-sPTHREAD_POOL_SIZE='Math.max((navigator.hardwareConcurrency\|\|4),4)'` | **Pre-warm** the worker pool at module init: KMA spawns its threads in the middle of a blocking `kma_main()`, when no new Workers can be created, so they must already exist. Sized to the device. |
| `-sFORCE_FILESYSTEM=1` | Keep MEMFS so JS can stage input/DB files |
| `-sEXIT_RUNTIME=0` | KMA's `main` returns per run; keep the module alive to read outputs |
| `-lworkerfs.js` `--no-entry` | WorkerFS support; module is a library of exports, not a `main()` |
| `EXPORTED_FUNCTIONS` / `EXPORTED_RUNTIME_METHODS` | `_kma_run`, … `_malloc`/`_free`; `ccall`, `cwrap`, `UTF8ToString`, `stringToUTF8`, `FS`, … |

We deliberately do **not** use `-sPROXY_TO_PTHREAD` — there is no `main()`; the JS side
loads the module inside its own Web Worker and calls `kma_run` via `ccall` (see §4).

## 3. The source changes (vs. upstream)

All port edits are wrapped in `#ifdef __EMSCRIPTEN__`, so apart from the upstream fixes
in `patches/` (§3d and §3e, applied to every build) the native build is byte-for-byte
upstream.

### 3a. Two new files
- **`wasm_shm_stubs.{c,h}`** — stubs for System V shared memory. `shmget` always returns
  `-1`, so KMA loads the database from disk into normal memory every run — correct for a
  browser. The header also defines the `_SYS_IPC_H` / `_SYS_SHM_H` guards so the real
  system headers are never pulled in.
- **`wasm_wrapper.c`** — the JS entry points. Each `kma_*_run(const char *args)` splits a
  space-separated string into a C `argv[]` (`argv[0]` = `"kma"`) and calls the matching
  upstream `*_main()`.

### 3b. Cross-stage pipeline (`kmapipe.c` + `wasm_ringpipe.c`)
Native KMA pipelines its stages (read-conversion `-s1`/`-s2` → mapping → scoring) by
running a child stage in a **thread** (`kmaPipeThread`, the native default) connected to
the parent by a `pipe()`, which **streams** with backpressure at ~constant memory (~160 MB
even for GB inputs).

The first WASM port used **`kmaPipeWasm()`**, which ran the child stage *sequentially*
in-process and handed its entire output to the parent through a `tmpfile()` in MEMFS. That
materialises the WHOLE stage intermediate in linear memory, so large samples (more than a
few hundred MB) died with `RangeError: Array buffer allocation failed` (MEMFS
`expandFileStorage`). Neither Memory64 nor input staging helps — the intermediate itself is
the problem.

Emscripten's own `pipe()` (PIPEFS) can't replace it: its reads are hard-coded non-blocking
(`EAGAIN` on empty) and its buffers are per-worker JS objects proxied to the runtime thread,
so a blocking read would deadlock the proxied write.

The fix (**`wasm_ringpipe.c`**, bound to `kmaPipe` for `__EMSCRIPTEN__`) keeps native's
threaded model but replaces the transport with a **bounded ring buffer in WASM linear
memory**. Under `-pthread` the linear memory is a `SharedArrayBuffer`, so the ring is
directly R/W from the producer pthread and the consumer runtime thread with no FS proxy, and
WebAssembly atomic wait/notify give real blocking + backpressure. The two ends are wrapped as
`FILE*` via `fopencookie()`, so KMA's `fread`/`fwrite` pipe sites are unchanged. Peak memory
drops from *whole intermediate* to a fixed 16 MiB.

**Key subtlety:** the producer stage reads its *input* through FS calls that Emscripten
proxies to the runtime (main) thread. If the consumer hard-blocked there (futex) it would
stall those and deadlock, so when the consumer is the main runtime thread and the ring is
empty it calls `emscripten_main_thread_process_queued_calls()` (and only naps ≤0.5 ms)
instead of blocking — servicing the producer's proxied reads so it can make progress. A
pthread consumer (nested pipeline) blocks normally, since the main thread stays free.

`PTHREAD_POOL_SIZE` is raised to `2*max(hwConcurrency,4)+2` so the concurrent producer +
consumer stages (each may use `-t N` within-stage threads) both have workers. (The
within-stage worker threads are upstream KMA code, untouched.) `kmaPipeWasm` is kept in
`kmapipe.c` as an unused fallback.

### 3c. Reading native 64-bit index files (`hashmapkma.c`)
`hashMapKMA_readHeader()` reads the header optimistically, sanity-checks it, and if that
fails re-reads the trailing fields as 8-byte values — auto-detecting the on-disk word
size. Under wasm64 `size_t` is 8 bytes and now matches the native index format, so this
heuristic is effectively a no-op, but it is kept (harmless, and keeps wasm32 buildable).

> The bundled `kma_index_*.{comp,length,seq}.b` / `.name` files come from a native
> 64-bit `kma index`. The WASM build **maps reads against pre-built indices**; it does
> not build them.

### 3d. Two upstream memory bugs, fixed in both builds (`patches/`)

Found by the 50-run native-vs-browser benchmark (September 2026), where two MinION
Nanopore runs crashed KMA in the browser at the VFDB step (`RuntimeError: memory access
out of bounds` in `malloc`, or a hang at `-t 1`) while native KMA completed. Native KMA
built with AddressSanitizer reports both; they are plain C bugs, not port issues.

- **`patches/0001-seqmenttree-use-after-free.patch`** — `rcpSeqmentTree()` copies a
  node into the grown array but writes the copied children's new addresses into the
  *old* node (`branch->branch[i] = …` instead of `dest->branch[i] = …`), and
  `resizeSeqmentTree()` then frees the old array. The next insertion walks freed memory.
  It triggers when one read collects more anker segments than the tree's first 64
  slots (long reads against a redundant database such as VFDB). Native macOS malloc
  leaves the freed block intact, so native KMA usually survives; Emscripten's malloc
  reuses it and corrupts its own heap.
- **`patches/0002-strjoin-nul-terminate.patch`** — `strjoin()` (the command line
  written to `.mapstat` as `## command`) never NUL-terminates the string, and sizes it
  with a different quoting rule than it writes with. The fix terminates it and counts
  quotes for `-i…` and `-o…` alike.

With both applied, native KMA is clean under ASan and its output is MD5-identical to
the unpatched build on the same input; the browser build completes the two failing
runs (at `-t 10` and `-t 1`) with VFDB `.res` MD5-identical to native. Apply to a fresh
upstream tree with `patch -p1 < patches/000N-*.patch` from `kma/`.

### 3e. printf formats for 64-bit counters (`patches/0003-wasm32-printf-formats.patch`)

Several `fprintf` calls pass `uint64_t` values with `%lu`/`%ld`, which is only right
where `long` is 64-bit. On wasm32 the varargs are misread and every later field is
garbage (`.mapstat` columns shift; `.res` is unaffected). The patch uses `PRIu64`
instead, which leaves native and wasm64 output unchanged; `-Wformat` is now clean on
native, wasm64 and wasm32.

### 3f. The wasm32 fallback build (`kma32.js` / `kma32.wasm`)

Safari as of 27.0, and therefore **every iOS browser** (all are WebKit), has no
Memory64: `kma.wasm` fails to compile there (`Memory64 is not enabled`). For those
browsers the same sources are also built as **wasm32 + pthreads** with
[`Makefile.wasm32`](Makefile.wasm32) (objects in `build32/`, so it never mixes with the
wasm64 objects):

```bash
emmake make -f Makefile.wasm32 clean
emmake make -f Makefile.wasm32
```

Differences from the wasm64 build:

- **Memory:** `MAXIMUM_MEMORY` is 4 GB (the wasm32 ceiling), but a shared memory
  reserves its maximum up front and phones refuse large reservations, so the runner
  creates the memory itself (`Module.wasmMemory`) with the largest maximum the device
  accepts: 4 GB, then 2 GB, 1 GB, 512 MB. The log says which one it got.
- **Threads:** `-t` is capped at 4, and the pool size is read from `Module.kmaPoolSize`
  (`-sPTHREAD_POOL_SIZE='Module.kmaPoolSize||10'`), which the runner sets to
  `2*max(t,4)+2`, the same rule as the wasm64 build.
- **Index files:** the native 64-bit indexes load unchanged; the header fields are read
  as `uint64_t` in both builds (§3c).

`js/engine.js` picks the build: wasm64 when the browser validates a Memory64 module and
accepts a 64-bit memory descriptor, wasm32 otherwise. `?engine=wasm32` forces the
fallback for testing on desktop.

Checked on the example data (`example2_*.fastq.gz`, fastp + ResFinder + CARD + VFDB)
in headless Chrome with both builds and in Playwright WebKit (Safari 26.5 engine, which
selects wasm32 by itself): `.res`, `.aln`, `.fsa`, `.mat.gz`, `.vcf.gz` and fastp's
`report.json` are byte-identical across all three; `.mapstat` differs only in the
`## command` line (`-t 4` vs `-t 10`) and `.frag.gz` only in line order, the same
thread-order effect seen between two native runs.

Known wasm32 limits, none of which the bundled databases reach: index allocations and
seeks are computed in 64 bits and truncated to 32-bit `size_t`/`long` (fine below 2 GB
per index file), and the per-template score sum `read_score` is a 32-bit `long` in
`runkma.c`/`spltdb.c` (it would need billions of aligned bases on one template).

## 4. How the browser drives it

KMA runs **off the main thread** because `pthread_join` blocks, and you cannot block the
main browser thread. The runtime ([`js/wasm-runtime.js`](../js/wasm-runtime.js) +
[`js/kma-runner.worker.js`](../js/kma-runner.worker.js)) works per run:

1. `runAnalysis()` loads the database's four index files (IndexedDB cache, else fetch)
   and reads the input FASTQs into `ArrayBuffer`s.
2. It spawns a **fresh** `kma-runner.worker.js` (a classic Web Worker) and transfers the
   bytes + a prebuilt MEMFS layout + the argument string (`-ipe`/`-i`, `-t_db`, `-o`,
   `-ef`, `-t <hardwareConcurrency>`, thresholds, `-bcNano` for Nanopore).
3. The worker `importScripts('../kma/kma.js')` with `Module.mainScriptUrlOrBlob` set to
   kma.js's URL — **required** so KMA's pthread pool spawns from kma.js itself, not from
   the runner script. It writes the files into MEMFS and calls
   `Module.ccall('kma_run', 'number', ['string'], [args])`.
4. Outputs (`.res`, `.aln`, `.fsa`, `.frag.gz`, …) are read back and posted to the main
   thread; the worker is then **terminated**.

Spawning a fresh worker per run (and terminating it) replaces the old whole-linear-memory
snapshot/restore hack: it gives clean KMA globals every run *and* frees the run's
(possibly multi-GB) memory back to the OS — both of which the snapshot approach could not
do once Memory64 made the heap large.

### Cross-origin isolation
`SharedArrayBuffer` (and therefore the whole pthread/shared-memory module) only loads when
the page is **cross-origin isolated** (`COOP: same-origin` + `COEP: require-corp`). The
site is static, so [`coi-serviceworker.min.js`](../coi-serviceworker.min.js) (vendored at
the repo root, loaded first in `index.html`) injects those headers via a service worker,
reloading once on first visit. This works on `localhost` and any HTTPS static host
(Render, GitHub Pages). All app assets are same-origin (no Google Fonts, no off-site
database), so the default `require-corp` mode is sufficient.

## 5. Verification

- Build compiles + links clean; the log shows `wasm64-emscripten` **`-mt`** system libs.
- Native vs. WASM agree exactly. Example data (`example2_*.fastq.gz`, ResFinder,
  `-ID 80 -mrc 0.8`) yields identical `.res` (blaOXA-181 98.25/98.62, blaTEM-158
  93.96/94.08, blaCTX-M-184 96.12/96.12) from native single-thread `kma` and from the
  10-thread wasm64 build in headless Chrome (`crossOriginIsolated === true`, ~0.4 s).
- A handy Node smoke test of `ccall('kma_run', …)` is **not** representative — Node's
  worker_threads pthread bootstrap can hang at pool init even though the browser path
  works. Validate in a real browser.

## 6. Tuning / future work

- **Thread count** = `navigator.hardwareConcurrency`, bounded by the baked pool size.
  Raise both (the `PTHREAD_POOL_SIZE` expression + the cap) together if desired.
- **`MAXIMUM_MEMORY`** (16 GB) only caps growth; bump it for very large samples.
- **WorkerFS** could back huge input FASTQs by `Blob` instead of copying them into MEMFS,
  cutting peak memory further. Inputs are currently written into MEMFS.
- **`-msimd128`** (WASM SIMD) may speed KMA's inner alignment loops if they vectorize.
