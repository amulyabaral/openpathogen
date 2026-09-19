# Building fastp for WebAssembly (OpenPathogen)

This documents how `fastp.js` / `fastp.wasm` in this directory are produced
from upstream [fastp](https://github.com/OpenGene/fastp) **v0.22.0** — the
same pattern as `kma/WASM_BUILD.md`. The patch (`wasm-port.patch`) is guarded
by `#ifdef __EMSCRIPTEN__` so the tree still builds natively.

## Why v0.22.0

v0.23+ hardwires ISA-L (`isa-l/igzip_lib.h`) for gzip I/O; ISA-L is
assembly-centric and does not port. v0.22.0 reads gzip through zlib, which
Emscripten provides via `-sUSE_ZLIB=1`.

## The build is wasm64 (Memory64) + pthreads

- **Memory64** lifts the wasm32 4 GB linear-memory ceiling, so full-size
  samples (inputs staged in MEMFS, gzipped outputs, duplication buffers) fit.
- **pthreads** run fastp's *native* producer/consumer/writer `std::thread`
  pipeline unmodified, so **the whole sample is processed** — there is no
  `--reads_to_process` cap any more. The wrapper still passes `--thread 1`
  (one *consumer*): with more consumers fastp's filtered output ordering is
  nondeterministic, and in the browser the proxied MEMFS writes serialise on
  the runtime thread anyway, so extra consumers buy no wall-clock win.

The earlier wasm32 build ran the pipeline sequentially on one thread and
capped processing at 400k reads, because a worker blocked inside `callMain`
cannot dispatch Emscripten's proxied FS operations and the pipeline deadlocks.
The fix is not "call main differently" but "don't block the runtime thread at
all": `fastp_run()` (in the patch's `src/wasm_wrapper.cpp`) spawns **one**
pthread that runs `main(argc, argv)`; the whole pipeline then executes on
pthreads while the worker's runtime thread sits in the JS event loop
servicing the proxied MEMFS reads/writes. fastp's own bounded pack queues
live in shared linear memory between real pthreads, so no ring-buffer
transport (as KMA's `wasm_ringpipe.c` needs) is required.

`fastp_state()` / `fastp_code()` poll the run (idle / running / exit code).
The wrapper also overrides `exit()`: fastp calls `exit(-1)` on error paths,
and the real Emscripten `exit()` would tear the runtime down and *drop the
proxied stderr prints* (the error message vanishes). The override records the
code, flushes stdio, and ends the calling thread — the module stays alive
(`EXIT_RUNTIME=0`) so outputs remain readable.

## Patch contents (`wasm-port.patch`)

1. **`src/wasm_wrapper.cpp` (new)** — `fastp_run` / `fastp_state` /
   `fastp_code` and the `exit()` override described above. The JS side calls
   `fastp_run` via `ccall` after staging inputs into MEMFS.
2. **`src/evaluator.cpp`** — adapter-detection `count*size` compared in
   `long long` under `#ifdef __EMSCRIPTEN__`. A no-op on this wasm64 (LP64)
   build and on native LP64, but protects the latent ILP32/LLP64 overflow
   (also present in 64-bit Windows builds) if the tree is ever built wasm32
   again.

The duplication buffer stays at the native `1L << 29` — under wasm64 the
`mBufLenInBits` shift no longer overflows, and native-size buffers keep the
reported duplication rate byte-identical to native fastp on full samples.

## Build

From a fresh upstream clone (or the pristine tree in
`benchmark/fastp_val/fastp-src/`), with the repo's `emsdk/` activated:

```bash
git clone https://github.com/OpenGene/fastp && cd fastp
git checkout v0.22.0
git apply /path/to/openpathogen/fastp/wasm-port.patch
source /path/to/openpathogen/emsdk/emsdk_env.sh
em++ -O2 -std=c++11 -pthread -sMEMORY64=1 -sUSE_ZLIB=1 \
     -sINVOKE_RUN=0 -sEXIT_RUNTIME=0 -sFORCE_FILESYSTEM=1 \
     -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=67108864 \
     -sMAXIMUM_MEMORY=17179869184 \
     -sSTACK_SIZE=1048576 -sDEFAULT_PTHREAD_STACK_SIZE=1048576 \
     -sPTHREAD_POOL_SIZE='Math.max((navigator.hardwareConcurrency||4),4)+12' \
     -sEXPORTED_FUNCTIONS='["_fastp_run","_fastp_state","_fastp_code","_malloc","_free"]' \
     -sEXPORTED_RUNTIME_METHODS='["ccall","cwrap","UTF8ToString","stringToUTF8","FS"]' \
     -o fastp.js src/*.cpp
```

Notes:

- The pool is pre-warmed and oversized (+12): fastp spawns 1 producer,
  `--thread` consumers and up to 7 writers mid-run, when no new workers can
  be created. Raise the expression if `--thread` goes above 4.
- `-pthread` needs `SharedArrayBuffer`, i.e. a cross-origin-isolated page —
  the site already forces this via `coi-serviceworker.min.js` for KMA. fastp
  therefore no longer works without cross-origin isolation (the old wasm32
  build's one advantage).
- 1 MiB stacks: fastp's option parsing exceeds the wasm default.
- No `-sPROXY_TO_PTHREAD` — the module runs inside its own worker and is
  driven via `ccall`, exactly like the KMA build.

## The wasm32 fallback build (`fastp32.js` / `fastp32.wasm`)

Safari as of 27.0, and so every iOS browser (all are WebKit), has no Memory64, so
`fastp.wasm` cannot load there. Those browsers get a wasm32 + pthreads build of the
same tree. On top of `wasm-port.patch`, apply `wasm32.patch`:

1. **`long` → `long long`** in every type (not in strings or comments). fastp keeps
   read, base and Q20/Q30 totals in `long`, which is 32-bit on wasm32 and would
   overflow past 2^31 bases (a few hundred MB of gzipped reads). The rewrite is
   mechanical: `wasm32-widen-long.py` reproduces it from the wasm64 tree. It makes
   the evaluator fix above redundant for this build.
2. **Duplication buffer 2 × 64 MiB** instead of 2 × 512 MiB, under
   `#if defined(__EMSCRIPTEN__) && !defined(__wasm64__)`: a 1 GB bitmap does not fit
   in a phone's browser tab, and `(1<<29)<<3` overflows a 32-bit `size_t`. Only the
   reported duplication rate changes, and only on large samples (it is an estimate
   in either build).

```bash
git apply /path/to/openpathogen/fastp/wasm-port.patch
git apply /path/to/openpathogen/fastp/wasm32.patch
em++ -O2 -std=c++11 -pthread -sUSE_ZLIB=1 \
     -sINVOKE_RUN=0 -sEXIT_RUNTIME=0 -sFORCE_FILESYSTEM=1 \
     -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=67108864 \
     -sMAXIMUM_MEMORY=4294967296 \
     -sSTACK_SIZE=1048576 -sDEFAULT_PTHREAD_STACK_SIZE=1048576 \
     -sPTHREAD_POOL_SIZE=12 \
     -sEXPORTED_FUNCTIONS='["_fastp_run","_fastp_state","_fastp_code","_malloc","_free"]' \
     -sEXPORTED_RUNTIME_METHODS='["ccall","cwrap","UTF8ToString","stringToUTF8","FS"]' \
     -o fastp32.js src/*.cpp
```

- No `-sMEMORY64`; the ceiling is 4 GB. A shared memory reserves its maximum up
  front and phones refuse large reservations, so the runner creates the memory
  itself (`Module.wasmMemory`) with the largest maximum the device accepts: 4 GB,
  then 2 GB, 1 GB, 512 MB.
- The pool is a fixed 12, since `--thread` is always 1 (1 wrapper + 1 producer +
  1 consumer + up to 7 writers = 10).

`js/engine.js` picks the build per browser; `?engine=wasm32` forces it on desktop.
On the example data, fastp32's `report.json` is byte-identical to fastp's, and the
filtered reads give the same KMA results.

## Gotchas encoded in `js/fastp-runner.worker.js`

- `Module.mainScriptUrlOrBlob` must point at `fastp.js` itself, or the
  pthread pool spawns the wrong file (same as the KMA runner).
- fastp only `gzopen`s paths **ending in `.gz`** — staged input names must
  match whether the incoming bytes are actually gzip. Outputs are named
  `filtered*.fastq.gz` (fastp gzips by extension), keeping full-sample output
  memory in MEMFS manageable; KMA gz-opens them unchanged.
- The worker is fresh-per-run and terminated afterwards, so the run's grown
  heap is returned to the OS between samples.
- `test.html` / `test-worker.js` in this directory drive the module directly
  (`ccall fastp_run` + state polling); they must be served cross-origin
  isolated (e.g. `node benchmark/harness/server.mjs`).
