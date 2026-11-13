/* integrity.js — pinned SHA-256 checksums for every database/index file the
 * app can download (KMA indexes, the CABBAGE association table and manifest,
 * the Deacon human index).
 *
 * Downloads — whether same-origin or from Zenodo — are verified against these
 * pins before they are used or cached, so a compromised or wrongly-updated
 * remote record cannot silently change which genes get called. Bytes read
 * back from the cache are verified on every load too: a mismatch evicts the
 * entry and forces a fresh download.
 *
 * Updating the pins: after rebuilding or re-deploying any database, run
 * scripts/pin-hashes.sh — it regenerates the table below from the local
 * files.
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
  // CABBAGE association table + manifest
  'databases/cabbage/manifest.json': '061bc62b86580834baab69ddc4cde23349aec3282504e54ce6c9202d28a9e5b1',
  'databases/cabbage/associations.json.gz': 'ecb9e230f588ea414bf574b8fe421eb0b986f993646a54b084ad70236d677166',
  // Deacon human index (deacon/build.sh --index)
  'databases/deacon/panhuman-1.k31w61.b16.pidx': '6ffda39e0bf5c6d6241b66bb32fb5173b793cbc98bca974d8078ccc475c98acd',
  // PointFinder references and catalogues (scripts/build-pointfinder.py)
  'databases/pointfinder/manifest.json': 'eb04d0e3df5a7945048b77b1a1d9602b75ecca268abd337134dc01de4818d24c',
  'databases/pointfinder/campylobacter/campylobacter.comp.b.gz': 'cbc4131e24f0091b632007a54dfd4e030e64ea2c925460522c3d4165ba00cf1e',
  'databases/pointfinder/campylobacter/campylobacter.length.b': 'f24c40b0026f23f9708bc8c689689a0eacfe57677d4db40a48fa297b692369f9',
  'databases/pointfinder/campylobacter/campylobacter.name': 'c0819b1995deba391c4b14778a30922490e015023e6c98355bc70a563082a820',
  'databases/pointfinder/campylobacter/campylobacter.seq.b': 'a463a5f6ad28f86b8322018f36f74a768f3dea32e5dd92f95e610c7970cebb12',
  'databases/pointfinder/campylobacter/catalogue.json.gz': 'fbbcfbcc50b5c2c22db3be5ebaf86bb0e0b6b4a33907c526ddc7fbaee32b7cf9',
  'databases/pointfinder/enterococcus_faecalis/catalogue.json.gz': '61fa7f97ee56ef7291e9db5e601c1086a2c04d3f5149960aaaf9f2227576c6d6',
  'databases/pointfinder/enterococcus_faecalis/enterococcus_faecalis.comp.b.gz': 'f2aaa8b8ddbfb8891f34e1d9b2a7d326f60e5f26bee2f7556529e67bd52fac2b',
  'databases/pointfinder/enterococcus_faecalis/enterococcus_faecalis.length.b': 'fdcff5521e0a3bd1bce11a8f485f3a183c4bb3ade0c2ef6f5a44b536072e1737',
  'databases/pointfinder/enterococcus_faecalis/enterococcus_faecalis.name': '305b3a732bd26cbf3f86b8f842ffe98a3f4b16d3683a901c546c630817bcbd1e',
  'databases/pointfinder/enterococcus_faecalis/enterococcus_faecalis.seq.b': '80bcfbda7376ffcd325b95d9bbc4f1a9608d259ea939d12c239e7f13c42f793f',
  'databases/pointfinder/enterococcus_faecium/catalogue.json.gz': 'b78b9cde918a7f3e3330feaca00e8e8b32c5ad07b74e8d77a60b2543effd5efc',
  'databases/pointfinder/enterococcus_faecium/enterococcus_faecium.comp.b.gz': 'd1753837389fb219817742c771816ac3b146719985134ee69fdf7237dd913324',
  'databases/pointfinder/enterococcus_faecium/enterococcus_faecium.length.b': 'fe5cf6313f0110acd021bbae00d1851cf3717417d6e6b5d368ad3965d31d21f3',
  'databases/pointfinder/enterococcus_faecium/enterococcus_faecium.name': 'e2f8a9d300d4722892d9bda0c3fe51f216fb7eb73626e9d9143c1966c45ea113',
  'databases/pointfinder/enterococcus_faecium/enterococcus_faecium.seq.b': '537460e8ec30fe6c822cc8701afe0f1e2bff5055732a56e0fcbe43e66cf7c2e9',
  'databases/pointfinder/escherichia_coli/catalogue.json.gz': '91d41f571e3ffeb20eed4737ec7fbbe63cfce5d6360d987de2a136876db12deb',
  'databases/pointfinder/escherichia_coli/escherichia_coli.comp.b.gz': '0bc9b5d36d7b575cb03c182c3be2457d79dbab72056a2650e4bebf3eb4460392',
  'databases/pointfinder/escherichia_coli/escherichia_coli.length.b': '86a50a54cc63cb50431c79e6a282bd2bb22bbf3851bc19859ad4e55e57f3a66f',
  'databases/pointfinder/escherichia_coli/escherichia_coli.name': '4443fba86096f24dab069597dc575401081d5de0cef15a61bca5bfc2fc66b3ca',
  'databases/pointfinder/escherichia_coli/escherichia_coli.seq.b': 'c2c984159dc817e1d6de83daeb7c2491dae2814ac0d0d3e629386371160cea68',
  'databases/pointfinder/helicobacter_pylori/catalogue.json.gz': '951ee207200275fe9b1ddcb9dd9297fb32dc45324641630392a1e2370cc2e7c3',
  'databases/pointfinder/helicobacter_pylori/helicobacter_pylori.comp.b.gz': '7fa0b325d308eed3d4d7b3f27baeb581ec8c849b565acda76162a98b9f58280a',
  'databases/pointfinder/helicobacter_pylori/helicobacter_pylori.length.b': 'b2faf0c2378ee26098cde54485e6afcf7b70eac092b64183e0ced6d92373228f',
  'databases/pointfinder/helicobacter_pylori/helicobacter_pylori.name': '665cb2a438fdb99c3903fb2ab0ce492b9811894e2766aa1bbb6c8fb6ffa3c7ca',
  'databases/pointfinder/helicobacter_pylori/helicobacter_pylori.seq.b': 'f3d09cd4a8c8fccf4cbbf81c22433bba848a01d081d714df7c6cafafeb7b4f40',
  'databases/pointfinder/klebsiella/catalogue.json.gz': '78d597a17cb4048b1f8ddce77dd3db63f0c03cd604289d4f5063137456ee6156',
  'databases/pointfinder/klebsiella/klebsiella.comp.b.gz': 'ec846f03abc153628cb9346f4767e82bf3067b4378bd5a260f9a07042d195454',
  'databases/pointfinder/klebsiella/klebsiella.length.b': '5ecb6c3ed3b816a93195a3c79d5e169e1ef57e1c9e722693fd36a7d7c55f8f28',
  'databases/pointfinder/klebsiella/klebsiella.name': 'a44058b4bc613b2a9dfd3f2fcc8f167c47ae2a69d91d8d80be54baa5872c4171',
  'databases/pointfinder/klebsiella/klebsiella.seq.b': 'd9f1d501d76cb32c04c8ea444ba4d25b3a5149cd40928c4c63a271f1eb9beecf',
  'databases/pointfinder/mycobacterium_tuberculosis/catalogue.json.gz': 'a98b85d27907d27f6c56e7106be7047ffb05bd3cfb174c1ffdcc3daf9cbbab72',
  'databases/pointfinder/mycobacterium_tuberculosis/mycobacterium_tuberculosis.comp.b.gz': '61e1ed0ad37a5d2daa8e3544adafa7b790dec85e569f9bb909f756d60dd1fd1c',
  'databases/pointfinder/mycobacterium_tuberculosis/mycobacterium_tuberculosis.length.b': '0f7e960f4a67d8335bb760bb38bcd3f21571b5e36256d6aa1269162b76d3a76f',
  'databases/pointfinder/mycobacterium_tuberculosis/mycobacterium_tuberculosis.name': '61f1ac0c4682e202a3c8bfba7c1bde487a91aad80480cc0f6236949d422d30bb',
  'databases/pointfinder/mycobacterium_tuberculosis/mycobacterium_tuberculosis.seq.b': '3aba6e4f8777341840d6de8c347f303a0ed37c57037339e68888fc120158cf91',
  'databases/pointfinder/neisseria_gonorrhoeae/catalogue.json.gz': 'fc8c52b6240869c81a92df6fbfcbd0594015e9865f3f56bda4c88ea6afc30725',
  'databases/pointfinder/neisseria_gonorrhoeae/neisseria_gonorrhoeae.comp.b.gz': 'dbf30645c41ba2506575dcc71365669be2dd121dcd38f89cfa999fdf2af7547f',
  'databases/pointfinder/neisseria_gonorrhoeae/neisseria_gonorrhoeae.length.b': '94477c81ab2891e15c17eb968cf66cedd6e33743b25479704dfe87da4d7b59c7',
  'databases/pointfinder/neisseria_gonorrhoeae/neisseria_gonorrhoeae.name': 'cef0aec83294f95d0fac59c3b5782ec0bef901ad44764a671418d534ade8d432',
  'databases/pointfinder/neisseria_gonorrhoeae/neisseria_gonorrhoeae.seq.b': 'a29ab65eaaddf7c9bcd1470fd8a68b4284d2c95777f619025074ab38d88498f4',
  'databases/pointfinder/salmonella/catalogue.json.gz': 'f887248f099f629c55b16f612edfb893714315d5c0474fd657855d2b6532965a',
  'databases/pointfinder/salmonella/salmonella.comp.b.gz': '9cf0b2ce63d36da39a58d3d26602d4e45606cd41fb405c05504a2fcb20c17ec6',
  'databases/pointfinder/salmonella/salmonella.length.b': '0a03ab784eec9653692d2dc769b3dab3acfe89e06f5cba41b77397ed7b47879e',
  'databases/pointfinder/salmonella/salmonella.name': 'e2c70e2a41acb6fa54e9ed8dcb358634bed937ac83ce5becb87cd9f46074ab3b',
  'databases/pointfinder/salmonella/salmonella.seq.b': 'f320cbc6331ecdd38550604baeef453e1de4e884cc8666f9104aac5ae721ffb3',
  'databases/pointfinder/staphylococcus_aureus/catalogue.json.gz': '0f25d32763306ba1f76499ccb9e8310efbf9697bf2f583e46e79cbfb68e68315',
  'databases/pointfinder/staphylococcus_aureus/staphylococcus_aureus.comp.b.gz': 'c98b7c686f06f6d79826d7b7457f0cc3373143231cd0d58afcb7dd14dad7ec9f',
  'databases/pointfinder/staphylococcus_aureus/staphylococcus_aureus.length.b': '6afc86b19687ac4491744ac7cb061568b9f7a61f88a1cb2f9195d5a11fdd0e09',
  'databases/pointfinder/staphylococcus_aureus/staphylococcus_aureus.name': '35ee452b183ab5af14198af1d32e45be85a50a6af2f25eefce5e3d545d76e1d2',
  'databases/pointfinder/staphylococcus_aureus/staphylococcus_aureus.seq.b': 'dd60c0d4828ab1bedc11dacaf59707d6faec5b7830be5f5769773190d16efff5',
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
