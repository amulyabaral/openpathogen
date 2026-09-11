/* integrity.js — pinned SHA-256 checksums for every database/index file the
 * app can download (KMA indexes, CABBAGE snapshots, the CABBAGE manifest).
 *
 * Downloads — whether same-origin or from Zenodo — are verified against these
 * pins before they are used or cached, so a compromised or wrongly-updated
 * remote record cannot silently change which genes get called. Bytes read
 * back from the cache are verified on every load too: a mismatch evicts the
 * entry and forces a fresh download.
 *
 * Updating the pins: after rebuilding or re-deploying any database, run
 * scripts/pin-hashes.sh — it regenerates the table below from the local
 * files. Bump the app.js/css version query in index.html when shipping.
 */

const PINNED_SHA256 = {
  // KMA indexes
  'databases/kma_index_resfinder_2_6_0.comp.b': '3689e42fb2154d5b69bbe9b52b433c725fb9bcc9470c40c1b568234050e853e6',
  'databases/kma_index_resfinder_2_6_0.length.b': 'e04a9d692c8d60b4e3c32da60d376ea4ec007c2a886d7d3078804fe6d5d5cda5',
  'databases/kma_index_resfinder_2_6_0.name': '55ab63f78601dc151eccaf86dd734b70b47989e3574bf4b0df7593670f125a0f',
  'databases/kma_index_resfinder_2_6_0.seq.b': 'b1d9ad0fb3c7c4da1d2fb2025e3ca0d1a363e313038ca38bbd5f7107c4fd9974',
  'databases/kma_index_card_4_0_1_homolog.comp.b': '3c1312a1ac35c7dbb486d19721926f3ba1e67b37853ed3462218f6fb2018a298',
  'databases/kma_index_card_4_0_1_homolog.length.b': '53c6a38600c76d0d2d268ece5ac5673ea47de1f3ae2dabae82baf748f6ccd944',
  'databases/kma_index_card_4_0_1_homolog.name': 'b4685ca3c34d661abe81f3a20c152dea1ce5e3a748f8c5e72b39b059aadd9ee8',
  'databases/kma_index_card_4_0_1_homolog.seq.b': 'fa703455dde98e7cef597ec6f2e21b7d15fd93ac94b59949b68314bf2f73941c',
  'databases/VFDB_setA_nt.fas.gz.comp.b': '96ffc8bbba7e580b476d4f32e874d6049f99314b8965b8969eba4050d7ca884b',
  'databases/VFDB_setA_nt.fas.gz.length.b': 'bf2a39f5cf2c5bb5cc5130279c72ca3c4553e99c2ac4a46e6048aa61d68682c5',
  'databases/VFDB_setA_nt.fas.gz.name': '48d09fa30ca97adff28ede5cb5c0e4e31a7b0c217de5761141df62f1cdbf8317',
  'databases/VFDB_setA_nt.fas.gz.seq.b': 'd5b1567d6e59fbfdec0ae5613cb558a6206df6dbc13ae9eeba54eb8a42087520',
  // CABBAGE snapshots + manifest
  'databases/cabbage/manifest.json': '2a8f5577c4eccaa2cbeb5b86338c9684150071b25d58da5e1c63d6d3880fa0d3',
  'databases/cabbage/associations.json.gz': 'a638e157174dc23e52cd733e3edda43d1176fbdb0aad518446f58faf380aa84b',
  'databases/cabbage/combined.cbg.gz': 'a29d44bb3dd1fa1aef81e0aa9cb1dd96fd49da755da5b5f6ec159616be403892',
  'databases/cabbage/genotypes.cbg.gz': 'b47c242a9e3daa2141a27cc47f34281b1c5116ddb7b4a558fd006f087f598bb9',
  'databases/cabbage/phenotypes.cbg.gz': 'e28b866d5c7dd584a5f17f8cd848ed0353728fabd95518acf38ab8eaac1c380a',
};

export function isPinned(path) {
  return Object.prototype.hasOwnProperty.call(PINNED_SHA256, path);
}

export async function sha256Hex(bytes) {
  if (!crypto?.subtle) {
    throw new Error('integrity checks need a secure context (https or localhost)');
  }
  const buf = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes);
  const digest = await crypto.subtle.digest('SHA-256', buf);
  return [...new Uint8Array(digest)].map(b => b.toString(16).padStart(2, '0')).join('');
}

// Throws a descriptive Error when the bytes do not match the pin. Unpinned
// paths pass through (fetch-run FASTQs from ENA are intentionally unpinned —
// they are reads, not call-defining references). Every load is hashed (a few
// hundred ms for the largest index, once per run — the cost of not trusting
// anything, including this session's own cache).
export async function verifyPinned(path, bytes) {
  const expected = PINNED_SHA256[path];
  if (!expected) return;
  const got = await sha256Hex(bytes);
  if (got !== expected) {
    throw new Error(
      `${path} failed its integrity check (expected SHA-256 ${expected.slice(0, 12)}…, `
      + `got ${got.slice(0, 12)}…). The download is corrupted or was modified in transit — `
      + 'it has not been used or cached.');
  }
}
