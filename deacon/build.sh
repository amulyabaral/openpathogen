#!/usr/bin/env bash
# build.sh — rebuild the human read filter from Deacon's 0.18.0 release:
#   deacon_wasm.js + deacon_wasm_bg.wasm   the browser engine (always)
#   ../databases/deacon/panhuman-1.k31w61.b16.pidx   the human index (--index)
# See WASM_BUILD.md. Needs git, curl, a Rust toolchain with the
# wasm32-unknown-unknown target, and wasm-bindgen-cli 0.2.128 (the version
# Deacon's Cargo.lock pins; `cargo install wasm-bindgen-cli --version 0.2.128`).
set -euo pipefail
cd "$(dirname "$0")"

TAG=0.18.0
COMMIT=ffcb7663dd3edf0bbc90fd6efe368c34f6f78301
SRC=build-src

if [ ! -d "$SRC" ]; then
  git clone --quiet --branch "$TAG" --depth 1 https://github.com/bede/deacon.git "$SRC"
fi
[ "$(git -C "$SRC" rev-parse HEAD)" = "$COMMIT" ] || { echo "$SRC is not deacon $TAG ($COMMIT)" >&2; exit 1; }

# 1. The engine: the deacon-wasm crate with the flags of Deacon's own
#    .cargo/config-wasm.toml (SIMD128 for the minimizer kernel).
(cd "$SRC" && RUSTFLAGS="-C target-feature=+reference-types,+simd128" \
  cargo build --locked --release --target wasm32-unknown-unknown -p deacon-wasm)
wasm-bindgen --target web --out-dir . "$SRC/target/wasm32-unknown-unknown/release/deacon_wasm.wasm"
rm -f deacon_wasm.d.ts deacon_wasm_bg.wasm.d.ts
echo "engine: $(wc -c < deacon_wasm_bg.wasm | tr -d ' ') bytes (deacon_wasm_bg.wasm)"

[ "${1:-}" = "--index" ] || exit 0

# 2. The index: panhuman-1 at k=31, w=61 (the compact build Deacon's
#    browser demo recommends), frozen into a binary fuse filter with 16-bit
#    fingerprints. Freezing is deterministic, so the SHA-256 pinned in
#    js/integrity.js is reproducible.
IDX_URL=https://objectstorage.uk-london-1.oraclecloud.com/n/lrbvkel2wjot/b/human-genome-bucket/o/deacon/3/panhuman-1.k31w61.idx
IDX_SHA256=a3abdaf4138936edb8bf6838663158e786678dcbaf0afb8acd0cb8b0641a680a
IDX="$SRC/panhuman-1.k31w61.idx"
OUT=../databases/deacon/panhuman-1.k31w61.b16.pidx

if ! echo "$IDX_SHA256  $IDX" | shasum -a 256 -c --status 2>/dev/null; then
  curl -L --fail --retry 5 -C - -o "$IDX" "$IDX_URL"
  echo "$IDX_SHA256  $IDX" | shasum -a 256 -c
fi
(cd "$SRC" && cargo build --locked --release)
mkdir -p "$(dirname "$OUT")"
"$SRC/target/release/deacon" index freeze -b 16 "$IDX" -o "$OUT"
echo "index: $OUT"
shasum -a 256 "$OUT"
