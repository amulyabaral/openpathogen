#!/usr/bin/env bash
# Assemble dist/ — everything Cloudflare Pages serves. Files over ~24 MiB
# (Pages' hard cap is 25 MiB) ship as `<name>.partNN` chunks instead;
# js/assets.js fetches the parts and reassembles them in the browser.
set -euo pipefail
cd "$(dirname "$0")/.."

rm -rf dist
mkdir -p dist/css dist/js/vendor dist/kma dist/fastp dist/hvprof dist/databases dist/resfinder_db

cp index.html about.html coi-serviceworker.min.js logo_dark.png logo_clear.png dist/
cp css/style.css dist/css/
cp -R js/. dist/js/
cp kma/kma.js kma/kma.wasm dist/kma/
cp fastp/fastp.js fastp/fastp.wasm dist/fastp/
cp hvprof/hvprof.wasm hvprof/hv_patho_v1.hvd dist/hvprof/
cp resfinder_db/phenotypes.txt dist/resfinder_db/

LIMIT=$((24 * 1024 * 1024))
CHUNK=$((16 * 1024 * 1024))

ship() { # repo-relative file path; chunked into dist/ if over the limit
  local f=$1 size
  size=$(wc -c < "$f" | tr -d ' ')
  mkdir -p "dist/$(dirname "$f")"
  if [ "$size" -gt "$LIMIT" ]; then
    split -b $CHUNK -d -a 2 "$f" "dist/$f.part"
    echo "  chunked $f → $(( (size + CHUNK - 1) / CHUNK )) parts"
  else
    cp "$f" "dist/$f"
  fi
}

for f in databases/*; do ship "$f"; done
ship SRR21386014_sub_1.fastq.gz
ship SRR21386014_sub_2.fastq.gz

echo "dist assembled: $(du -sh dist | cut -f1), $(find dist -type f | wc -l | tr -d ' ') files"
