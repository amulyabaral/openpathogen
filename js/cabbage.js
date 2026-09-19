/* cabbage.js — tested-resistance rates from CABBAGE, the AMR genotype–
 * phenotype database (Dickens et al., Nucleic Acids Research 2026,
 * doi:10.1093/nar/gkag780), served client-side from a snapshot of the
 * EMBL-EBI AMR portal at a pinned release.
 *
 * The report uses one thing from it: for a detected gene, how often isolates
 * of the same species carrying that gene tested resistant to the antibiotics
 * the gene acts on. The association table behind this was precomputed by
 * scripts/build-cabbage.py (the genotypes and phenotypes views joined on the
 * isolate); each row is (species, gene, antibiotic, R, I, S). Genes are
 * AMRFinderPlus element symbols, each with the antibiotics the portal links
 * it to.
 *
 * Files are fetched same-origin (databases/cabbage/), then from Zenodo when
 * CABBAGE_BASE is set, verified against pinned hashes and cached in
 * IndexedDB like the KMA indexes.
 */

import { cacheDBFile, getCachedDBFile, deleteDBFile } from './db.js';
import { fetchAssetWithProgress } from './assets.js';
import { verifyPinned } from './integrity.js';
import { antibioticClasses, linkedAntibiotics } from './antibiotic-classes.js';

// Remote copy of databases/cabbage/ (the /api/records/…/files/<name>/content
// form is the CORS-enabled one). Empty: same-origin only.
const CABBAGE_BASE = '';

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
    // A failed attempt (offline) must not stick for the whole session.
    manifestPromise.catch(() => { manifestPromise = null; });
  }
  return manifestPromise;
}

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

// ── Association table ──

let associationsPromise = null;

// Loads once per session (about 190 KB). Concurrent callers share the load.
export function loadAssociations() {
  if (!associationsPromise) {
    associationsPromise = (async () => {
      let expectedBytes = 0;
      try { expectedBytes = (await getManifest()).associations?.bytes || 0; } catch (_) {}
      const gz = await fetchSnapshotBytes('associations.json.gz', expectedBytes);
      const parsed = JSON.parse(window.fflate.strFromU8(window.fflate.gunzipSync(gz)));
      // Lookup: "species\x00gene" -> rows.
      const bySpeciesGene = new Map();
      const speciesSet = new Set();
      for (const row of parsed.rows) {
        speciesSet.add(row.species);
        const key = row.species + '\x00' + row.gene;
        let list = bySpeciesGene.get(key);
        if (!list) bySpeciesGene.set(key, list = []);
        list.push(row);
      }
      const genes = new Map();
      for (const g of parsed.genes || []) {
        genes.set(g.gene, { symbol: g.symbol || g.gene, cls: g.class || '', links: new Set(g.links || []) });
      }
      return {
        minN: parsed.min_n || 20,
        species: parsed.species.filter(s => speciesSet.has(s)),
        bySpeciesGene,
        genes,
        keyIndex: new Map(),
      };
    })();
    associationsPromise.catch(() => { associationsPromise = null; });
  }
  return associationsPromise;
}

// ── Matching detected genes to CABBAGE symbols ──
//
// ResFinder, CARD and AMRFinderPlus spell the same gene differently. The
// lookup tries, in order: a known alias, the canonical symbol, the symbol
// without its trailing allele letter (aph(3')-IIIa vs aph(3')-III), the
// symbol without a trailing number (ermA1 vs ermA), and a "bla" prefix for
// beta-lactamases named without it (CARD SHV-52 vs blaSHV-52). Fusion
// symbols in CABBAGE (aac(6')-Ie/aph(2'')-Ia) are also indexed by each
// half. The report names the matched CABBAGE symbol, so a loose match is
// visible to the reader.

// Canonical gene symbol — mirrors canon_gene() in scripts/build-cabbage.py.
// Strips quotes and parentheses and trailing _N row suffixes ("vanA_2" ->
// "vanA", "erm(C)" -> "ermC", "ermC'" -> "ermC").
function geneCanon(g) {
  let s = String(g).replace(/['’()]/g, '');
  for (;;) {
    const m = s.match(/^(.*)_(\d+)$/);
    if (!m) break;
    s = m[1];
  }
  return s;
}

const ALIASES = {
  'aac6-aph2': 'aac6-Ie/aph2-Ia',         // ResFinder aac(6')-aph(2'')
  'aac6-ie-aph2-ia': 'aac6-Ie/aph2-Ia',   // CARD AAC(6')-Ie-APH(2'')-Ia
};

// Canonical symbol without a trailing allele letter, lower-cased: after a
// digit ("blaTEM-1B" -> "blatem-1") or after an upper-case roman numeral
// ("aph3-IIIa" -> "aph3-iii"). The allele letter after a numeral is lower
// case, so "aph3-VI" keeps its I.
function looseKey(c) {
  return c.replace(/(\d)[A-Za-z]$|([IVX])[a-z]$/, '$1$2').toLowerCase();
}

// Candidate lookup keys for a detected gene, most specific first.
function geneKeys(name) {
  const c = geneCanon(name);
  const lc = c.toLowerCase();
  const out = [];
  const add = (k) => { if (k && !out.includes(k)) out.push(k); };
  if (ALIASES[lc]) add(ALIASES[lc].toLowerCase());
  add(lc);
  add(looseKey(c));
  const d = c.replace(/\d+$/, '');
  if (d && d !== c) { add(d.toLowerCase()); add(looseKey(d)); }
  if (!/^bla/i.test(c)) { add('bla' + lc); add('bla' + looseKey(c)); }
  return out;
}

function speciesIndex(assoc, species) {
  let idx = assoc.keyIndex.get(species);
  if (idx) return idx;
  idx = new Map();
  const prefix = species + '\x00';
  const reg = (k, g) => {
    if (!k) return;
    let l = idx.get(k);
    if (!l) idx.set(k, l = []);
    if (!l.includes(g)) l.push(g);
  };
  for (const key of assoc.bySpeciesGene.keys()) {
    if (!key.startsWith(prefix)) continue;
    const g = key.slice(prefix.length);
    reg(g.toLowerCase(), g);
    reg(looseKey(g), g);
    if (g.includes('/')) {
      for (const part of g.split('/')) { reg(part.toLowerCase(), g); reg(looseKey(part), g); }
    }
  }
  assoc.keyIndex.set(species, idx);
  return idx;
}

function evidence(assoc, species, gene) {
  return (assoc.bySpeciesGene.get(species + '\x00' + gene) || [])
    .reduce((n, r) => Math.max(n, r.r + r.i + r.s), 0);
}

// The CABBAGE gene (canonical symbol) a detected gene corresponds to for
// this species, or null. Ambiguous loose matches resolve to the gene with
// the most tested isolates.
function matchAssociationGene(assoc, species, gene) {
  const idx = speciesIndex(assoc, species);
  for (const k of geneKeys(gene)) {
    const cands = idx.get(k);
    if (!cands || !cands.length) continue;
    if (cands.length === 1) return cands[0];
    return cands.slice().sort((a, b) => evidence(assoc, species, b) - evidence(assoc, species, a))[0];
  }
  return null;
}

function geneInfo(assoc, gene) {
  return assoc.genes.get(gene) || { symbol: gene, cls: '', links: new Set() };
}

// Antibiotic names as CABBAGE writes them: lower case, combinations joined
// with "-" (ResFinder writes "Amoxicillin+Clavulanic acid").
function cabbageName(ab) {
  const k = String(ab).trim().toLowerCase().replace(/\s*\+\s*/g, '-');
  return k === 'rifampicin' ? 'rifampin' : k;
}

// A portal link that names a drug (mecA → methicillin, ceftaroline) rather
// than a class (blaZ → "beta-lactam antibiotic", blaCTX-M-15 → "cephalosporin").
const isDrugLink = (l) => !!antibioticClasses(l) || String(l).toLowerCase() === 'kanamycin a';

// Tested-resistance rates for one detected gene in one species, for the
// antibiotics the gene is known to act on: those ResFinder lists for it
// (resfinderAntibiotics) and those the portal links to it by drug name.
// Class-level links are not used: "beta-lactam antibiotic" would put
// methicillin under blaZ and ceftazidime-avibactam under CTX-M-15, rates that
// come from other genes the same isolates carry. names: the gene's names,
// tried in order. Drug-linked antibiotics come first (mecA: methicillin,
// ceftaroline), then by number of isolates tested.
export function linkedEvidence(assoc, species, names, resfinderAntibiotics) {
  const minN = assoc.minN || 20;
  const listed = new Set((resfinderAntibiotics || []).map(cabbageName));
  for (const name of names) {
    const matched = matchAssociationGene(assoc, species, name);
    if (!matched) continue;
    const info = geneInfo(assoc, matched);
    const drugLinks = [...info.links].filter(isDrugLink);
    const rows = [];
    for (const r of assoc.bySpeciesGene.get(species + '\x00' + matched) || []) {
      const n = r.r + r.i + r.s;
      if (n < minN) continue;
      const ab = r.antibiotic.toLowerCase();
      const linked = drugLinks.some(l => linkedAntibiotics(l).has(ab));
      if (linked || listed.has(ab)) rows.push({ antibiotic: r.antibiotic, r: r.r, i: r.i, s: r.s, n, linked });
    }
    rows.sort((a, b) => b.linked - a.linked || b.n - a.n);
    return { symbol: info.symbol, detectedAs: name, rows };
  }
  return null;
}

// Fuzzy-resolve an organism string ("Staphylococcus aureus subsp. aureus
// USA300_TCH1516" from ENA, or "K. pneumoniae" typed by hand) against the
// species vocabulary: exact match first, then genus (or its initial) plus
// species epithet.
export function resolveSpecies(speciesList, organism) {
  if (!organism) return null;
  const q = String(organism).toLowerCase().trim().replace(/\s+/g, ' ');
  for (const s of speciesList) if (s.toLowerCase() === q) return s;
  const [g, e] = q.split(' ');
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
