# Human read removal: Deacon in WebAssembly

`deacon_wasm.js` / `deacon_wasm_bg.wasm` are the `deacon-wasm` crate of
[Deacon](https://github.com/bede/deacon) **0.18.0** (tag `0.18.0`, commit
`ffcb766`), unmodified, built by `build.sh`. Deacon removes host reads by
counting the minimizers a read (or read pair) shares with a host index; with
the default thresholds a read is human when at least 2 of its distinct
minimizers, and at least 1% of them, are in the index.

Built with rustc 1.97.1, `wasm32-unknown-unknown` with SIMD128 and
reference types (Deacon's own `.cargo/config-wasm.toml`), and
wasm-bindgen-cli 0.2.128 (`--target web`). The module is 32-bit and
single-threaded, so the same file runs in every browser including Safari on
iPhone (SIMD128 needs Safari 16.4 or later).

## The index

`databases/deacon/panhuman-1.k31w61.b16.pidx` (250,609,707 bytes,
SHA-256 pinned in `js/integrity.js`) is Deacon's panhuman-1 index —
HPRC year 1 ∪ CHM13v2.0 ∪ GRCh38.p14, minus minimizers shared with FDA-ARGOS
bacteria and RefSeq viruses; CC0 — built at k=31, w=61 (the compact
variant Deacon's browser demo recommends; 111,300,492 minimizers) and frozen
with `deacon index freeze -b 16` into a binary fuse filter with 16-bit
fingerprints. The exact index is 890 MB; the filter is 239 MiB and loads in
under a second. `build.sh --index` downloads the exact index (SHA-256
checked), builds native Deacon from the same tag and freezes it; freezing is
deterministic, so the output matches the pin.

A fuse filter has no false negatives, and a false-positive rate of 2^-16
per minimizer lookup. With the 2-hit threshold a lookup error can only
change the call for a read that already has exactly one genuine human hit.

## Validation (2026-09-25)

Browser engine (under Node, streaming 1 MiB chunks as the worker does)
against native Deacon 0.18.0 (`deacon filter -d --ordered`), paired reads:

| Run | Reads | Removed as human | Browser + frozen index vs native + exact index |
| --- | ---: | ---: | --- |
| example2 (bundled test pair) | 100,000 | 14 | identical output |
| SRR10341524, S. aureus (the app's example) | 515,164 | 26 | identical output |
| SRR849803, K. pneumoniae | 14,195,034 | 250 | identical output (MD5) |
| SRR13449732, human low-pass WGS, 2×150 | 512,196 | 504,504 (98.5%) | identical output |
| ERR13720561, cervicovaginal metagenome | 1,017,074 | 150,012 (14.7%) | 3 more pairs removed (native: 150,006 reads) |

The 3 pairs were each one hit short of the threshold in the exact index
(`--debug`: 1 hit) and reached 2 hits through a fingerprint collision. Native
Deacon on the frozen index produces byte-identical output to the browser
engine, so the difference is the index approximation, not the port.

The browser engine filters about 15 MB/s of gzipped input on one Apple M4
core (SRR849803, 1.4 GB: 86 s; native Deacon on 8 threads: 9 s), including
decompressing the input and gzipping the kept reads.
