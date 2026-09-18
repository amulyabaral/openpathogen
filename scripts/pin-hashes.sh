#!/usr/bin/env bash
# pin-hashes.sh — regenerate the pinned SHA-256 table in js/integrity.js from
# the local database files.
#
# Run this after rebuilding or re-publishing ANY database file (KMA indexes
# under databases/, CABBAGE snapshots under databases/cabbage/) and deploy
# the updated js/integrity.js together with the new files. Clients refuse any
# downloadable file whose hash is not pinned or does not match — shipping new
# bytes without new pins therefore breaks downloads, by design.
set -euo pipefail
cd "$(dirname "$0")/.."

INDEX_FILES=(
  databases/kma_index_resfinder_2_6_0.comp.b
  databases/kma_index_resfinder_2_6_0.length.b
  databases/kma_index_resfinder_2_6_0.name
  databases/kma_index_resfinder_2_6_0.seq.b
  databases/kma_index_card_4_0_1_homolog.comp.b
  databases/kma_index_card_4_0_1_homolog.length.b
  databases/kma_index_card_4_0_1_homolog.name
  databases/kma_index_card_4_0_1_homolog.seq.b
  databases/VFDB_setA_nt.fas.gz.comp.b
  databases/VFDB_setA_nt.fas.gz.length.b
  databases/VFDB_setA_nt.fas.gz.name
  databases/VFDB_setA_nt.fas.gz.seq.b
)
CABBAGE_FILES=(
  databases/cabbage/manifest.json
  databases/cabbage/associations.json.gz
  databases/cabbage/combined.cbg.gz
  databases/cabbage/genotypes.cbg.gz
  databases/cabbage/phenotypes.cbg.gz
)

row() { printf "  '%s': '%s',\n" "$1" "$(shasum -a 256 "$1" | cut -d' ' -f1)"; }

table_file=$(mktemp)
{
  echo "  // KMA indexes"
  for f in "${INDEX_FILES[@]}"; do row "$f"; done
  echo "  // CABBAGE snapshots + manifest"
  for f in "${CABBAGE_FILES[@]}"; do row "$f"; done
} > "$table_file"

# Replace the table between `const PINNED_SHA256 = {` and `};` (newline-free
# -v values are not portable across awks, so the table goes via a temp file).
awk -v tblfile="$table_file" '
  /^const PINNED_SHA256 = \{$/ { print; while ((getline line < tblfile) > 0) print line; close(tblfile); skip=1; next }
  skip && /^\};$/ { skip=0; print; next }
  skip { next }
  { print }
' js/integrity.js > js/integrity.js.tmp && mv js/integrity.js.tmp js/integrity.js
rm -f "$table_file"

echo "Pinned $((${#INDEX_FILES[@]} + ${#CABBAGE_FILES[@]})) files in js/integrity.js."
echo "Deploy js/integrity.js together with the new files."
