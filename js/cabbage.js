/* cabbage.js — CABBAGE, the comprehensive AMR genotype–phenotype database
 * (Dickens et al., Nucleic Acids Research 2026, doi:10.1093/nar/gkag780),
 * served entirely client-side from a compact snapshot of the EMBL-EBI AMR
 * portal (https://www.ebi.ac.uk/amr) at a pinned release.
 *
 * The snapshot (CBG1, see scripts/build-cabbage.py) stores each portal view
 * as column dictionaries plus varint row indices. This module decodes it
 * and powers the "Phenotype prediction (CABBAGE)" card in the results tab:
 *   - the precomputed gene→phenotype association table (join of ~165k
 *     isolates that have both a genome and an antibiogram), and
 *   - the full phenotypes view, for country / year / isolation-source
 *     filtered species antibiograms.
 *
 * Column data is streamed straight off the varint bytes (one varint per
 * row, so a column scan is a single sequential walk), which keeps the
 * ~1.7M-row view at a few MB instead of hundreds.
 *
 * Files are fetched same-origin first (databases/cabbage/), then from
 * Zenodo (CABBAGE_BASE), and cached in IndexedDB like the KMA indexes.
 */

import { cacheDBFile, getCachedDBFile, deleteDBFile } from './db.js';
import { fetchAssetWithProgress } from './assets.js';
import { verifyPinned } from './integrity.js';

// Snapshot files: same-origin databases/cabbage/ first (local dev / any
// static host that ships them), then Zenodo (set by
// scripts/deploy-cabbage-zenodo.sh; the /api/records/…/files/<name>/content
// form is the CORS-enabled one).
const CABBAGE_BASE = '';

const viewsByKey = new Map(); // decoded snapshots, kept for the session

// ── Snapshot acquisition (network or IndexedDB) ──

let manifestPromise = null;
export function getManifest() {
  if (!manifestPromise) {
    manifestPromise = (async () => {
      const tryUrls = CABBAGE_BASE
        ? ['databases/cabbage/manifest.json', CABBAGE_BASE + 'manifest.json/content']
        : ['databases/cabbage/manifest.json'];
      for (const url of tryUrls) {
        let bytes = null;
        try {
          const resp = await fetch(url);
          if (resp.ok) {
            const ct = resp.headers.get('content-type') || '';
            if (!ct.includes('text/html')) bytes = new Uint8Array(await resp.arrayBuffer());
          }
        } catch (_) { /* try next */ }
        if (bytes) {
          // Checksum-verify before trusting the file list (throws with a clear
          // message if the manifest does not match the pinned release).
          await verifyPinned('databases/cabbage/manifest.json', bytes);
          return JSON.parse(new TextDecoder().decode(bytes));
        }
      }
      throw new Error('CABBAGE manifest not found — run scripts/fetch-cabbage.py and scripts/build-cabbage.py');
    })();
  }
  return manifestPromise;
}

export async function loadView(key, onProgress) {
  if (viewsByKey.has(key)) return viewsByKey.get(key);
  const manifest = await getManifest();
  const info = manifest.views[key];
  if (!info) throw new Error(`Unknown CABBAGE view: ${key}`);
  const gz = await fetchSnapshotBytes(info.file, info.bytes, onProgress);
  const bytes = window.fflate.gunzipSync(gz);
  const view = decodeSnapshot(bytes);
  if (view.rowCount !== info.rows) {
    throw new Error(`Snapshot corrupt: ${key} has ${view.rowCount} rows, manifest says ${info.rows}`);
  }
  viewsByKey.set(key, view);
  return view;
}

export function isViewLoaded(key) { return viewsByKey.has(key); }

async function fetchSnapshotBytes(file, expectedBytes, onProgress) {
  const cacheKey = 'cabbage:' + file;
  const pinPath = `databases/cabbage/${file}`;
  let data = await getCachedDBFile(cacheKey);
  if (data && data.byteLength) {
    try {
      await verifyPinned(pinPath, data);
      return data;
    } catch (_) {
      // Cached copy does not match the pinned release — evict, re-download.
      try { await deleteDBFile(cacheKey); } catch (_) {}
    }
  }
  const tryPaths = CABBAGE_BASE
    ? [`databases/cabbage/${file}`, `${CABBAGE_BASE}${file}/content`]
    : [`databases/cabbage/${file}`];
  let lastErr = null;
  for (const p of tryPaths) {
    try {
      data = await fetchAssetWithProgress(p, expectedBytes, onProgress);
      if (data.length) {
        // Verify against the pinned hash before caching or using the bytes.
        await verifyPinned(pinPath, data);
        try { await cacheDBFile(cacheKey, data); } catch (_) {}
        return data;
      }
    } catch (err) { lastErr = err; }
  }
  throw lastErr || new Error('snapshot download failed');
}

// ── CBG1 decoder ──

const TD = new TextDecoder('utf-8');

function decodeSnapshot(u8) {
  let p = 0;
  const magic = String.fromCharCode(u8[p++], u8[p++], u8[p++], u8[p++]);
  if (magic !== 'CBG1') throw new Error('Not a CBG1 snapshot');
  const version = u8[p++];
  if (version !== 1) throw new Error(`Unsupported CBG version ${version}`);

  const str = () => {
    const len = u8[p] | (u8[p + 1] << 8);
    p += 2;
    const s = TD.decode(u8.subarray(p, p + len));
    p += len;
    return s;
  };

  const release = str();
  const viewId = u8[p++];
  const name = str();
  const rowCount = (u8[p] | (u8[p + 1] << 8) | (u8[p + 2] << 16) | (u8[p + 3] << 24)) >>> 0;
  p += 4;
  const nCols = u8[p] | (u8[p + 1] << 8);
  p += 2;

  const columns = [];
  for (let c = 0; c < nCols; c++) {
    const id = str();
    const kind = u8[p++] === 0 ? 's' : 'n';
    if (kind === 's') {
      const nDict = (u8[p] | (u8[p + 1] << 8) | (u8[p + 2] << 16) | (u8[p + 3] << 24)) >>> 0;
      p += 4;
      const dict = new Array(nDict + 1); // index 0 = null
      for (let i = 1; i <= nDict; i++) dict[i] = str();
      columns.push({ id, kind, dict, varints: null, keys: null, vals: null });
    } else {
      const vals = new Float64Array(rowCount);
      for (let i = 0; i < rowCount; i++) {
        // little-endian f64, read byte-wise (DataView-free, avoids alignment)
        const b = [];
        for (let k = 0; k < 8; k++) b[k] = u8[p++];
        F64BYTES.set(b);
        vals[i] = F64VIEW[0];
      }
      columns.push({ id, kind, dict: null, varints: null, keys: null, vals });
    }
  }

  // Remaining bytes are the concatenated varint blocks, in column order.
  for (const col of columns) {
    if (col.kind !== 's') continue;
    let end = p, n = 0;
    while (n < rowCount && end < u8.length) {
      let b;
      do { b = u8[end++]; } while (b & 0x80);
      n++;
    }
    col.varints = u8.subarray(p, end);
    p = end;
  }

  return { id: viewId, name, release, rowCount, columns, colById: new Map(columns.map(c => [c.id, c])) };
}

const F64BYTES = new Uint8Array(8);
const F64VIEW = new Float64Array(F64BYTES.buffer);

// ── Column access ──

// Cached Uint32Array of dictionary indices (random access).
function keys(col) {
  if (col.keys) return col.keys;
  const v = col.varints;
  if (!v) return new Uint32Array(0);
  let n = 0;
  for (let i = 0; i < v.length; i++) if (!(v[i] & 0x80)) n++;
  const arr = new Uint32Array(n);
  let p = 0, i = 0;
  while (p < v.length) {
    let shift = 0, val = 0, b;
    do {
      b = v[p++];
      val |= (b & 0x7F) << shift;
      shift += 7;
    } while (b & 0x80);
    arr[i++] = val >>> 0;
  }
  col.keys = arr;
  return arr;
}

// ── Gene→phenotype associations (results-card engine) ──
//
// The association table is precomputed at snapshot-build time by joining the
// genotypes and phenotypes views on BioSample (~165k isolates have both a
// genome and an antibiogram — far richer than the portal's own merged view).
// Each row: for (species, gene), how many isolates carrying that gene were
// recorded R / I / S per antibiotic, phenotypes taken from the updated
// (2025 CLSI/EUCAST) breakpoint reinterpretation first.

let associationsCache = null;

export async function loadAssociations(onProgress) {
  if (associationsCache) return associationsCache;
  let expectedBytes = 0;
  try { expectedBytes = (await getManifest()).associations?.bytes || 0; } catch (_) {}
  const gz = await fetchSnapshotBytes('associations.json.gz', expectedBytes, onProgress);
  const parsed = JSON.parse(window.fflate.strFromU8(window.fflate.gunzipSync(gz)));
  // Lookup: "species\x00gene" -> rows; plus per-species row lists.
  const bySpeciesGene = new Map();
  const speciesSet = new Set();
  for (const row of parsed.rows) {
    speciesSet.add(row.species);
    const key = row.species + '\x00' + row.gene;
    let list = bySpeciesGene.get(key);
    if (!list) bySpeciesGene.set(key, list = []);
    list.push(row);
  }
  associationsCache = {
    release: null, // filled by callers from the manifest
    species: parsed.species.filter(s => speciesSet.has(s)),
    bySpeciesGene,
    background: parsed.background,
  };
  return associationsCache;
}

// Canonical gene symbol — mirrors canon_gene() in scripts/build-cabbage.py.
// Strips quotes and parentheses and trailing _N row suffixes ("vanA_2" ->
// "vanA", "erm(C)" -> "ermC", "ermC'" -> "ermC") so AMRFinderPlus and
// ResFinder/CARD spelling variants unify.
export function geneCanon(g) {
  let s = String(g).replace(/['\u2019()]/g, '');
  for (;;) {
    const m = s.match(/^(.*)_(\d+)$/);
    if (!m) break;
    s = m[1];
  }
  return s;
}

// Extract a gene symbol from a KMA template name, per source database:
// ResFinder "mecA_1_U72714" -> "mecA";
// CARD "gb|HQ845196.1|+|0-861|ARO:3001109|SHV-52 [Klebsiella pneumoniae]" -> "SHV-52".
// Returned as detected (canonicalisation happens at match time).
export function templateGene(template, db) {
  if (!template) return null;
  if (db === 'card_homolog') {
    const bars = template.split('|');
    if (bars.length >= 6) {
      const gene = bars[5].split(' [')[0].trim();
      if (gene) return gene;
    }
    return template;
  }
  return template.split('_')[0];
}

// Match a detected gene against the association vocabulary. Exact canonical
// match first; then a single trailing digit is dropped as a fallback
// ("ermA1" vs ResFinder's "ermA") — the matched symbol is always shown, so
// a loose match is visible to the user.
export function matchAssociationGene(assoc, species, gene) {
  const canon = geneCanon(gene);
  const key = species + '\x00' + canon;
  if (assoc.bySpeciesGene.has(key)) return canon;
  const stripped = canon.replace(/\d+$/, '');
  if (stripped && stripped !== canon && assoc.bySpeciesGene.has(species + '\x00' + stripped)) {
    return stripped;
  }
  return null;
}

// Predicted phenotypic resistance for one species from the detected genes.
// Returns one entry per antibiotic with >=1 informative gene (n >= minN),
// plus the species background rate where available.
export function predictPhenotypes(assoc, species, detectedGenes, minN = 20) {
  const byAntibiotic = new Map();
  for (const g of detectedGenes) {
    const matched = matchAssociationGene(assoc, species, g.name);
    if (!matched) continue;
    for (const row of assoc.bySpeciesGene.get(species + '\x00' + matched) || []) {
      const n = row.r + row.i + row.s;
      if (n < minN) continue;
      let e = byAntibiotic.get(row.antibiotic);
      if (!e) byAntibiotic.set(row.antibiotic, e = { antibiotic: row.antibiotic, genes: [] });
      e.genes.push({ gene: matched, detectedAs: g.name, db: g.db, r: row.r, i: row.i, s: row.s, n });
    }
  }
  const bg = new Map(assoc.background
    .filter(b => b.species === species)
    .map(b => [b.antibiotic, b]));
  const out = [...byAntibiotic.values()].map(e => {
    e.genes.sort((a, b) => b.n - a.n);
    const strongest = e.genes.reduce((acc, g) => Math.max(acc, (g.r + g.i) / g.n), 0);
    const evidence = e.genes.reduce((acc, g) => Math.max(acc, g.n), 0);
    const b = bg.get(e.antibiotic);
    return {
      ...e,
      rate: strongest,
      evidence,
      verdict: strongest >= 0.85 ? 'resistant' : strongest >= 0.4 ? 'uncertain' : 'susceptible',
      background: b ? { r: b.r, i: b.i, s: b.s, n: b.r + b.i + b.s } : null,
    };
  });
  // Most evidence first, then rate: a 100%-resistant n=21 fluke must not
  // outrank a 97%-resistant n=1,791 signal.
  out.sort((a, b) => b.evidence - a.evidence || b.rate - a.rate);
  return out;
}

// Which detected genes found no CABBAGE data for this species (so the card
// can list them with their ResFinder drug class instead).
export function unmatchedGenes(assoc, species, detectedGenes) {
  return detectedGenes.filter(g => !matchAssociationGene(assoc, species, g.name));
}

// Instant background antibiogram for a species (all 1.7M AST records,
// aggregated at build time — no phenotype-view download needed).
export function backgroundAntibiogram(assoc, species) {
  return assoc.background
    .filter(b => b.species === species)
    .map(b => ({ antibiotic: b.antibiotic, r: b.r, i: b.i, s: b.s, n: b.r + b.i + b.s }))
    .sort((a, b) => b.n - a.n);
}

// Fuzzy-resolve an organism string ("S. aureus USA300_TCH1516" from ENA)
// against the species vocabulary; exact match first, then abbreviated
// genus + species epithet.
export function resolveSpecies(speciesList, organism) {
  if (!organism) return null;
  const q = String(organism).toLowerCase().trim();
  for (const s of speciesList) if (s.toLowerCase() === q) return s;
  const [g, e] = q.split(/\s+/);
  if (!g || !e) return null;
  const abbrev = g.endsWith('.') ? g.slice(0, -1) : null;
  const cands = [];
  for (const s of speciesList) {
    const toks = s.toLowerCase().split(/\s+/);
    const genusOk = abbrev ? toks[0].startsWith(abbrev) : toks[0] === g;
    if (genusOk && toks[1] === e) cands.push(s);
  }
  return cands.length === 1 ? cands[0] : null;
}

// Runtime-filtered antibiogram from the full phenotypes view (country /
// year / isolation-source filters). Requires the 5.7 MB view to be loaded.
export function antibiogram(view, species, { country = null, sourceCategory = null, yearFrom = null, yearTo = null } = {}) {
  const spCol = view.colById.get('phenotype-species');
  const abCol = view.colById.get('phenotype-antibiotic_name');
  const phCol = view.colById.get('phenotype-resistance_phenotype');
  const clsiCol = view.colById.get('phenotype-Updated_phenotype_CLSI');
  const eucastCol = view.colById.get('phenotype-Updated_phenotype_EUCAST');
  const ctryCol = country && view.colById.get('phenotype-country');
  const srcCol = sourceCategory && view.colById.get('phenotype-isolation_source_category');
  const yrCol = (yearFrom != null || yearTo != null) && view.colById.get('phenotype-collection_year');
  if (!spCol || !abCol || !phCol) return [];

  const spIdx = keys(spCol), abIdx = keys(abCol), phIdx = keys(phCol);
  const clsiIdx = clsiCol ? keys(clsiCol) : null;
  const eucastIdx = eucastCol ? keys(eucastCol) : null;
  const ctryIdx = ctryCol ? keys(ctryCol) : null;
  const srcIdx = srcCol ? keys(srcCol) : null;
  const yrV = yrCol ? yrCol.vals : null;

  const speciesIdx = new Set();
  for (let i = 1; i < spCol.dict.length; i++) {
    if (spCol.dict[i] === species) speciesIdx.add(i);
  }
  const ctrySel = ctryIdx ? new Set(ctryCol.dict.map((v, i) => (v === country ? i : -1)).filter(i => i > 0)) : null;
  const srcSel = srcIdx ? new Set(srcCol.dict.map((v, i) => (v === sourceCategory ? i : -1)).filter(i => i > 0)) : null;

  // Phenotype bucket by dictionary index, preferring updated breakpoints.
  const bucketOf = (phDict, di) => {
    if (!di) return null;
    const ph = phDict[di];
    if (!ph) return null;
    if (ph === 'resistant' || ph === 'non-susceptible') return 'r';
    if (ph.includes('intermediate')) return 'i';
    if (ph.startsWith('susceptible')) return 's';
    return null;
  };
  const clsiBucket = clsiCol ? new Array(clsiCol.dict.length) : null;
  if (clsiBucket) for (let i = 1; i < clsiCol.dict.length; i++) clsiBucket[i] = bucketOf(clsiCol.dict, i);
  const eucastBucket = eucastCol ? new Array(eucastCol.dict.length) : null;
  if (eucastBucket) for (let i = 1; i < eucastCol.dict.length; i++) eucastBucket[i] = bucketOf(eucastCol.dict, i);
  const phBucket = new Array(phCol.dict.length);
  for (let i = 1; i < phCol.dict.length; i++) phBucket[i] = bucketOf(phCol.dict, i);

  const counts = new Map(); // antibiotic -> [r, i, s]
  const n = view.rowCount;
  for (let r = 0; r < n; r++) {
    if (!speciesIdx.has(spIdx[r])) continue;
    if (ctrySel && !ctrySel.has(ctryIdx[r])) continue;
    if (srcSel && !srcSel.has(srcIdx[r])) continue;
    if (yrV) {
      const y = yrV[r];
      if (Number.isNaN(y)) continue;
      if (yearFrom != null && y < yearFrom) continue;
      if (yearTo != null && y > yearTo) continue;
    }
    const b = (clsiBucket && clsiBucket[clsiIdx[r]]) || (eucastBucket && eucastBucket[eucastIdx[r]]) || phBucket[phIdx[r]];
    if (!b) continue;
    const ab = abCol.dict[abIdx[r]];
    let e = counts.get(ab);
    if (!e) counts.set(ab, e = [0, 0, 0]);
    e['ris'.indexOf(b)]++;
  }
  return [...counts.entries()]
    .map(([antibiotic, e]) => ({ antibiotic, r: e[0], i: e[1], s: e[2], n: e[0] + e[1] + e[2] }))
    .sort((a, b) => b.n - a.n);
}

export function fmtInt(n) {
  if (n == null) return '—';
  if (n >= 1e6) return (n / 1e6).toFixed(2) + 'M';
  if (n >= 1e3) return (n / 1e3).toFixed(1) + 'k';
  return String(n);
}
