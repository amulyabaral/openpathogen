# Building fastp for WebAssembly (OpenPathogen)

This documents how `fastp.js` / `fastp.wasm` in this directory are produced
from upstream [fastp](https://github.com/OpenGene/fastp) **v0.22.0** — the
same pattern as `kma/WASM_BUILD.md`. Two patches are required, both guarded
by `#ifdef __EMSCRIPTEN__` so the tree still builds natively
(`fastp/wasm-port.patch`).

## Why v0.22.0

v0.23+ hardwires ISA-L (`isa-l/igzip_lib.h`) for gzip I/O; ISA-L is
assembly-centric and does not port. v0.22.0 reads gzip through zlib, which
Emscripten provides via `-sUSE_ZLIB=1`.

## Patch 1 — sequential execution (peprocessor.cpp, seprocessor.cpp)

fastp runs a producer/consumer/writer `std::thread` pipeline. Two problems
in a browser worker: (a) without `-pthread`, `std::thread` construction
aborts in `-fno-exceptions` mode; (b) with `-pthread`, the main thread
blocked inside `callMain` cannot dispatch proxied operations and the
pipeline deadlocks or exits abnormally.

The patch runs the stages back-to-back on one thread:
`producerTask(); consumerTask(each); writeTask(each writer);`.
This is safe **only because the wrapper caps processing with
`--reads_to_process 400000`**, below `PACK_IN_MEM_LIMIT` (500 packs × 1000
reads), so the producer never blocks on the bounded pack queue or writer
buffers. QC statistics on the first 400k reads are unchanged to any
practical precision (fastp itself samples the first reads for adapter
detection). The wrapper must keep that cap — remove it and the sequential
build deadlocks.

## Patch 2 — duplication buffer overflow (duplicate.cpp)

`Duplicate` sets `mBufLenInBytes = 1L<<29` and computes
`mBufLenInBits = mBufLenInBytes << 3` = 2³². On wasm32, `long`/`size_t` are
32-bit, so the shift overflows to **0**, and `statRead`'s
`positions[i] % mBufLenInBits` traps with "remainder by zero". The 2 GiB
buffer also cannot live in a browser heap. The patch uses 64 MiB buffers
(with the bits computation cast to `uint64`), which keeps duplication
estimation accurate at our read caps.

## Build

```bash
git clone https://github.com/OpenGene/fastp && cd fastp
git checkout v0.22.0
git apply /path/to/openpathogen/fastp/wasm-port.patch
em++ -O2 -std=c++11 -sUSE_ZLIB=1 -sEXIT_RUNTIME=1 -sINVOKE_RUN=0 \
     -sFORCE_FILESYSTEM=1 -sALLOW_MEMORY_GROWTH=1 \
     -sEXPORTED_RUNTIME_METHODS='["FS","callMain"]' \
     -o fastp.js src/*.cpp
```

No pthreads, no SharedArrayBuffer needed — fastp runs even without
cross-origin isolation (KMA still requires it).

## Gotchas encoded in `js/fastp-runner.worker.js`

- fastp only `gzopen`s paths **ending in `.gz`** — staged input names must
  match whether the incoming bytes are actually gzip.
- Outputs are plain FASTQ in MEMFS; they are re-packaged as `File` objects
  so the KMA runner can stream them into its own worker unchanged.
- `test.html` / `test-worker.js` in this directory are a flag-bisect harness
  used to find the wasm32 overflow; keep them for future debugging.
