/* pointfinder.js — chromosomal point mutations with PointFinder (CGE), read
 * from the reads the way the rest of the app reads genes.
 *
 * PointFinder holds, for each species, wild-type reference sequences of the
 * genes whose mutations confer resistance (gyrA, parC, 23S rRNA, promoters
 * and more) and a curated table of those mutations. The reads are mapped to
 * the species' references with KMA, with the same command as for acquired
 * genes, and every curated position is read from KMA's per-base counts
 * (.mat.gz): the base most reads carry there, how many do, and whether the
 * position has reads enough to tell. scripts/build-pointfinder.py prepares
 * the index and a catalogue of the rows, checked against the references.
 *
 * Positions (catalogue rows keep phenotypes.txt's coordinates): codon n of a
 * coding template is template bases P+3n-2 .. P+3n, nucleotide p of an RNA
 * gene (or of the gene part of a promoter template) is base P+p, promoter
 * position -k is base P+1-k, where P is the promoter length (else 0).
 *
 * Where this differs from ResFinder 4.7.2, on purpose:
 *   - a position the reads do not reach is "not determined": KMA's .aln
 *     writes it as "-", the same as a deletion, and ResFinder reads it as one
 *   - the copies of a multi-copy rRNA locus (E. coli 16S rrsB/C/H, N.
 *     gonorrhoeae 23S a1–a4) are each read from the pooled reads of all the
 *     copies, since KMA divides the reads between them; a mutation that only
 *     some of the reads carry is reported with its share, because a few
 *     mutated copies can be enough
 *   - insertions and deletions are compared as sequence around the site, so
 *     an event that reads place differently in a repeat still matches
 *   - required mutations are checked once every position has been read,
 *     against any change at the required position that the row allows
 *     (ResFinder only counts partners that are rows themselves, which leaves
 *     some requirements unsatisfiable), on whichever template of the gene
 *     was read (C. coli gyrA rows require T86I of the C. jejuni template)
 *   - a gene with a template per species (Campylobacter gyrA) is read on the
 *     template the reads match; the other's rows are not counted
 *   - a row PointFinder marks as under development is reported as found only
 *     when AMRFinderPlus curates the same mutation (build-pointfinder.py);
 *     otherwise it is listed apart, "use with caution" (every Klebsiella
 *     row is marked so: its porin and acrR rows read differences between
 *     lineages, its gyrA and parC rows are well established)
 */

import { cacheDBFile, getCachedDBFile, deleteDBFile } from './db.js';
import { fetchAssetWithProgress } from './assets.js';
import { verifyPinned } from './integrity.js';
import { parseResFile } from './resfinder-db.js';
import { parseTaxon } from './genes.js';
import { matrixOf, alignTemplates } from './evidence.js';

const BASE = 'databases/pointfinder/';
const BASES = 'ACGT';
const STOPS = new Set(['TAA', 'TAG', 'TGA']);

// The calls a row's evaluation rests on: a position counts as read from
// this depth, a base as called when more than half the reads carry it, and
// a minority as present from this share, with at least 3 reads. An rRNA
// gene has several copies (one of E. coli's seven 23S copies is 14% of the
// reads); in a single-copy gene a share under 20% was, in the truth set, a
// context error of the sequencer (pmrB V161G at 10 and 14% in two
// unrelated E. coli).
export const POINT_CALL = {
  minDepth: { illumina: 5, nanopore: 10 },
  minorityShare: { illumina: 0.2, nanopore: 0.25 },
  minorityShareRna: { illumina: 0.1, nanopore: 0.25 },
  minorityReads: 3,
  // ResFinder's gene filter, on the .res row of a template (percent)
  minTemplateIdentity: 80,
  minTemplateCoverage: 1,
};

// ── Species ──

// PointFinder's species by database folder, with the sample species that use
// each: ResFinder's rule (the genus_species folder, else the genus folder),
// plus Shigella, which is E. coli genetically (AMRFinderPlus files both under
// Escherichia). The Klebsiella folder is the K. pneumoniae complex only
// (strict): its genes are K. pneumoniae's, 79 to 89% identical to those of
// K. oxytoca and K. aerogenes, too far for the reads to be read on them.
const SPECIES = [
  { folder: 'campylobacter', genus: 'campylobacter', species: ['jejuni', 'coli', 'lari', 'upsaliensis', 'fetus'] },
  { folder: 'enterococcus_faecalis', genus: 'enterococcus', epithet: 'faecalis' },
  { folder: 'enterococcus_faecium', genus: 'enterococcus', epithet: 'faecium' },
  { folder: 'escherichia_coli', genus: 'escherichia', epithet: 'coli' },
  { folder: 'escherichia_coli', genus: 'shigella', species: ['sonnei', 'flexneri', 'dysenteriae', 'boydii'] },
  { folder: 'helicobacter_pylori', genus: 'helicobacter', epithet: 'pylori' },
  { folder: 'klebsiella', genus: 'klebsiella', species: ['pneumoniae', 'quasipneumoniae', 'variicola', 'quasivariicola', 'africana'], strict: true },
  { folder: 'mycobacterium_tuberculosis', genus: 'mycobacterium', epithet: 'tuberculosis' },
  { folder: 'neisseria_gonorrhoeae', genus: 'neisseria', epithet: 'gonorrhoeae' },
  { folder: 'salmonella', genus: 'salmonella', species: ['enterica', 'bongori'] },
  { folder: 'staphylococcus_aureus', genus: 'staphylococcus', epithet: 'aureus' },
];

// The PointFinder folder for a species as typed or from run metadata
// ("Escherichia coli", "E. coli", "Klebsiella pneumoniae subsp. ozaenae"),
// or null. A whole genus folder takes any species of the genus when the
// genus is spelled out (unless it is strict), and its known species when it
// is abbreviated (an initial alone cannot tell Salmonella from Shigella).
export function pointfinderSpecies(text) {
  const t = parseTaxon(text);
  if (!t) return null;
  const hits = SPECIES.filter(s => (t.genus ? t.genus === s.genus : t.initial === s.genus[0])
    && (s.epithet ? s.epithet === t.epithet : (t.genus && !s.strict) || s.species.includes(t.epithet)));
  if (!hits.length) return null;
  const specific = hits.filter(s => s.epithet);
  const pick = specific.length ? specific : hits;
  return new Set(pick.map(s => s.folder)).size === 1 ? pick[0].folder : null;
}

// The species PointFinder covers, for the input page's note.
export const POINTFINDER_COVERS = 'Campylobacter, E. coli and Shigella, Enterococcus faecalis and faecium, H. pylori, the K. pneumoniae complex, M. tuberculosis, N. gonorrhoeae, Salmonella and S. aureus';

// ── Files ──

let manifestPromise = null;

// databases/pointfinder/manifest.json: version, commit and the species.
export function loadPointFinderManifest() {
  if (!manifestPromise) {
    manifestPromise = (async () => {
      const path = BASE + 'manifest.json';
      const bytes = await fetchPinned(path);
      return JSON.parse(new TextDecoder().decode(bytes));
    })();
    manifestPromise.catch(() => { manifestPromise = null; });
  }
  return manifestPromise;
}

// Fetch a pinned file: the IndexedDB cache when its bytes still match the
// pin, else the network (same origin only), verified before it is cached.
async function fetchPinned(path) {
  let data = await getCachedDBFile(path);
  if (data && data.length) {
    try { await verifyPinned(path, data); return data; } catch (_) { try { await deleteDBFile(path); } catch (__) {} }
  }
  data = await fetchAssetWithProgress(path, 0);
  await verifyPinned(path, data);
  if (data.length) { try { await cacheDBFile(path, data); } catch (_) {} }
  return data;
}

// The KMA index of a species, as wasm-runtime.js runAnalysis takes it. Its
// hash table ships gzipped (scripts/build-pointfinder.py).
export function pointfinderDatabase(folder, manifest) {
  const sp = manifest?.species?.find(s => s.folder === folder);
  return {
    prefix: `${BASE}${folder}/${folder}`,
    shortName: 'pointfinder',
    displayName: `PointFinder ${manifest?.version || ''} (${sp?.label || folder})`.replace('  ', ' '),
    sameOriginOnly: true,
    gzipped: ['.comp.b'],
    species: folder,
  };
}

const catalogues = new Map();

// A species' catalogue: templates (with sequences) and the checked rows.
export function loadCatalogue(folder) {
  if (!catalogues.has(folder)) {
    const p = (async () => {
      const gz = await fetchPinned(`${BASE}${folder}/catalogue.json.gz`);
      return JSON.parse(globalThis.fflate.strFromU8(globalThis.fflate.gunzipSync(gz)));
    })();
    p.catch(() => catalogues.delete(folder));
    catalogues.set(folder, p);
  }
  return catalogues.get(folder);
}

// ── Reading positions ──

const CODE = (() => {
  const aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
  const map = new Map();
  let i = 0;
  for (const x of 'TCAG') for (const y of 'TCAG') for (const z of 'TCAG') map.set(x + y + z, aa[i++]);
  return map;
})();
const translate = (codon) => CODE.get(codon) || 'X';

// 0-based template index of a row's nucleotide position.
const tindex = (t, pos) => (pos > 0 ? t.promoter + pos - 1 : t.promoter + pos);

// Counts at one position of a view { counts, n } (A C G T N -) and the call:
// the base (or "-" for a deletion) more than half the reads carry, when the
// depth is enough. reason says why there is no call.
function callAt(view, i, minDepth) {
  if (i < 0 || i >= view.n) return { base: null, depth: 0, reason: 'outside' };
  const c = view.counts;
  const k = i * 6;
  const depth = c[k] + c[k + 1] + c[k + 2] + c[k + 3] + c[k + 4] + c[k + 5];
  if (depth < minDepth) return { base: null, depth, reason: depth ? 'low depth' : 'no reads' };
  let best = 0;
  for (const j of [1, 2, 3, 5]) if (c[k + j] > c[k + best]) best = j;
  const share = c[k + best] / depth;
  if (share <= 0.5) return { base: null, depth, share, reason: 'mixed' };
  return { base: best === 5 ? '-' : BASES[best], depth, share };
}

const countOf = (view, i, base) => {
  const j = base === '-' ? 5 : BASES.indexOf(base);
  return j < 0 || i < 0 || i >= view.n ? 0 : view.counts[i * 6 + j];
};
const depthAt = (view, i) => {
  if (i < 0 || i >= view.n) return 0;
  let d = 0;
  for (let j = 0; j < 6; j++) d += view.counts[i * 6 + j];
  return d;
};

// Bases inserted after position i that most reads passing it carry, from
// the matrix's insertion rows (a row's A C G T N count the reads carrying a
// base there, its "-" the reads without).
function insertedAfter(view, i) {
  const rows = view.ins?.get(i);
  if (!rows) return '';
  let s = '';
  for (const r of rows) {
    const withBase = r[0] + r[1] + r[2] + r[3] + r[4];
    if (withBase <= r[5]) continue;
    let best = 0;
    for (let j = 1; j < 4; j++) if (r[j] > r[best]) best = j;
    s += BASES[best];
  }
  return s;
}

// The sample's sequence over template indices [from, to), as the reads have
// it (deletions out, insertions in), or null when some position is not read.
function sampleSeq(view, from, to, minDepth) {
  let s = '';
  for (let i = Math.max(0, from); i < Math.min(view.n, to); i++) {
    const c = callAt(view, i, minDepth);
    if (!c.base) return null;
    if (c.base !== '-') s += c.base;
    s += insertedAfter(view, i);
  }
  return s;
}

// ── One species run ──

// evaluatePointFinder(result, catalogue, { nanopore }) reads every row of the
// catalogue on the KMA result of the species' index. It returns
//   { species, label, genes: [gene], found: [mutation], minority: [mutation],
//     needs: [mutation], caution: [mutation], other: [change],
//     otherChanges: [change at a listed codon, not in the table],
//     undetermined: [{ row, reason }],
//     calls: [{ row, status, read, depth, share, reason }], screened, read,
//     speciesMatch: { matched, total } }
// A mutation is { gene, template, locus, name, protein, symbols, codon, kind,
// region, pos, ref, alt, depth, share, antibiotics, classes, pmids, notes,
// mechanism, status: 'found' | 'minority' | 'needs' | 'caution', requires,
// rows }. calls
// has every row, as read: 'mutation', 'minority', 'other' (another change
// there), 'wildtype' or 'not read'.
export function evaluatePointFinder(result, cat, { nanopore = false } = {}) {
  const mode = nanopore ? 'nanopore' : 'illumina';
  const minDepth = POINT_CALL.minDepth[mode];
  const res = new Map();
  if (result?.resTable) for (const r of parseResFile(result.resTable).rows) res.set(String(r.Template).trim(), r);
  const tmpl = new Map(cat.templates.map(t => [t.name, t]));
  const geneClasses = classesOfGenes(cat.rows);

  // The reads of each template: its own matrix, or for a copy of a
  // multi-copy rRNA locus the pooled matrices of all copies, mapped onto it.
  const own = new Map();
  for (const t of cat.templates) {
    let m = null;
    try { m = res.has(t.name) ? matrixOf(result, t.name) : null; } catch (_) { m = null; }
    own.set(t.name, m);
  }
  const loci = new Map();
  for (const t of cat.templates) {
    if (!loci.has(t.locus)) loci.set(t.locus, []);
    loci.get(t.locus).push(t);
  }
  const views = new Map();
  for (const t of cat.templates) {
    const copies = t.rna ? loci.get(t.locus) : [t];
    const n = t.length;
    const counts = new Uint32Array(n * 6);
    let any = false;
    for (const c of copies) {
      const m = own.get(c.name);
      if (!m) continue;
      any = true;
      const map = c === t ? null : alignTemplates(t.seq, c.seq)?.map;
      if (c !== t && !map) continue;
      for (let i = 0; i < n; i++) {
        const j = map ? map[i] : i;
        if (j < 0 || j >= m.n) continue;
        for (let k = 0; k < 6; k++) counts[i * 6 + k] += m.counts[j * 6 + k];
      }
    }
    const m = own.get(t.name);
    views.set(t.name, any ? { n, counts, ins: m?.ins || new Map(), pooled: copies.length > 1 } : null);
  }

  // ResFinder's gene filter, per template; a copy of an rRNA locus passes
  // with its locus.
  const pass = (t) => {
    const r = res.get(t.name);
    return !!r && Number(r.Template_Identity) >= POINT_CALL.minTemplateIdentity && Number(r.Template_Coverage) >= POINT_CALL.minTemplateCoverage;
  };
  const passes = new Map(cat.templates.map(t => [t.name, t.rna ? loci.get(t.locus).some(pass) : pass(t)]));

  // A coding gene with a template per species (Campylobacter gyrA) is read
  // on the templates that pass; when none does, on the one the reads match
  // best, which then counts as not found.
  const byGene = new Map();
  for (const t of cat.templates) if (!t.rna) { if (!byGene.has(t.gene)) byGene.set(t.gene, []); byGene.get(t.gene).push(t); }
  const ti = (t) => Number(res.get(t.name)?.Template_Identity) || 0;
  const used = (t) => {
    const all = t.rna ? [t] : byGene.get(t.gene);
    if (all.length < 2 || passes.get(t.name)) return true;
    if (all.some(o => passes.get(o.name))) return false;
    return all.slice().sort((a, b) => ti(b) - ti(a))[0] === t;
  };

  const genes = cat.templates.filter(t => t.screened && used(t)).map(t => {
    const r = res.get(t.name);
    return {
      template: t.name, gene: t.gene, locus: t.locus, rna: t.rna, promoter: t.promoter, length: t.length,
      identity: r ? Number(r.Template_Identity) : null, coverage: r ? Number(r.Template_Coverage) : null,
      depth: r ? Number(r.Depth) : null, read: passes.get(t.name), rows: 0, determined: 0, found: 0,
    };
  });
  const geneOf = new Map(genes.map(g => [g.template, g]));

  // ── Evaluate rows ──
  const hits = []; // { row, t, alt, depth, share, status: 'found' | 'minority' }
  const undetermined = [];
  const other = [];
  const calls = [];
  const aaAt = new Map(); // `${gene}:${codon}` -> amino acid read (for requirements)
  for (const row of cat.rows) {
    const t = tmpl.get(row.template);
    const g = geneOf.get(row.template);
    if (!t || !g) continue; // an alternative template the reads did not choose
    g.rows++;
    const view = views.get(t.name);
    if (!passes.get(t.name) || !view) {
      undetermined.push({ row, reason: 'gene not found' });
      calls.push({ row, status: 'not read', reason: 'gene not found' });
      continue;
    }
    const out = evaluateRow(row, t, view, minDepth, mode);
    calls.push({
      row, status: CALL_STATUS[out.status], read: out.read || '', depth: out.depth ?? null, share: out.share ?? null,
      reason: out.reason || (out.status === 'other' ? out.change : ''),
    });
    if (out.status === 'undetermined') { undetermined.push({ row, reason: out.reason }); continue; }
    g.determined++;
    if (out.aa) aaAt.set(`${t.gene}:${row.pos}`, out.aa);
    if (out.status === 'hit' || out.status === 'minority') {
      hits.push({ row, t, alt: out.alt, depth: out.depth, share: out.share, status: out.status === 'hit' ? 'found' : 'minority', codon: out.codon, haplotype: out.haplotype, at: out.at });
    } else if (out.status === 'other') {
      other.push({ row, t, gene: t.gene, template: t.name, kind: row.kind, change: out.change, depth: out.depth, share: out.share });
    }
  }

  // Required mutations: a row with requirements counts only when one of its
  // groups is present in full (members: a found change at that position
  // that the requirement allows).
  const present = (m) => {
    const a = aaAt.get(`${tmpl.get(m.template)?.gene || m.template}:${m.pos}`);
    return !!a && a !== m.ref && m.alts.includes(a);
  };
  const satisfied = (row) => !row.requires?.length || row.requires.some(group => group.every(present));

  // One mutation per (template, position, kind, change), with the
  // antibiotics of every row whose conditions hold. Copies of an rRNA locus
  // that report the same site (on the pooled reads) are one mutation.
  const byKey = new Map();
  for (const h of hits) {
    // Indel rows that give the same sequence (an insertion written at four
    // places of a repeat) are one mutation, named by the first row.
    const key = h.haplotype ? `${h.t.name}|${h.row.kind}|hap:${h.haplotype}` : `${h.t.name}|${h.row.kind}|${h.row.pos}|${h.alt}`;
    if (!byKey.has(key)) byKey.set(key, { h, rows: [] });
    byKey.get(key).rows.push(h.row);
  }
  const mutations = [];
  for (const { h, rows } of byKey.values()) {
    const ok = rows.filter(satisfied);
    const use = ok.length ? ok : rows;
    // Under development in PointFinder and not curated by AMRFinderPlus.
    const caution = use.every(r => r.caution && !r.amrfinder);
    const antibiotics = uniq(use.flatMap(r => r.antibiotics));
    const altIdx = (r) => Math.max(0, r.alts.indexOf(h.alt));
    const symbols = h.row.symbols ? [h.row.symbols[altIdx(h.row)]].filter(Boolean)
      : indelSymbols(h.row, h.alt, h.t, h.at);
    const mut = {
      gene: h.t.gene, template: h.t.name, locus: h.t.locus, accession: h.t.accession,
      kind: h.row.kind, region: h.row.region, pos: h.row.pos, ref: h.row.ref, alt: h.alt,
      name: mutationName(h.row, h.alt, h.t),
      protein: proteinName(h.row, h.alt, h.t, h.at),
      symbols,
      codon: h.codon || null,
      depth: h.depth, share: h.share,
      status: h.status === 'minority' ? 'minority' : !ok.length ? 'needs' : caution ? 'caution' : 'found',
      antibiotics, classes: uniq(use.flatMap(r => rowClasses(r, geneClasses))),
      pmids: uniq(use.flatMap(r => r.pmids[altIdx(r)] || r.pmids.flat())).sort(),
      notes: uniq(use.map(r => r.notes).filter(Boolean)),
      mechanism: use.map(r => r.mechanism).find(Boolean) || '',
      requires: ok.length ? null : rows.map(r => requirementText(r.requires, tmpl)).filter(Boolean)[0] || null,
      fixed: use.map(r => r.fixed).find(Boolean) || null,
      rows: rows.map(r => r.line),
    };
    mutations.push(mut);
  }
  // rRNA copies of one locus report one site each on the pooled reads.
  const seen = new Map();
  const deduped = [];
  for (const m of mutations.sort((a, b) => b.depth - a.depth)) {
    const t = tmpl.get(m.template);
    if (!t.rna || loci.get(t.locus).length < 2) { deduped.push(m); continue; }
    const canon = canonicalSite(m, t, loci.get(t.locus));
    const key = `${t.locus}|${m.kind}|${canon}|${m.alt}`;
    if (seen.has(key)) continue;
    seen.set(key, true);
    deduped.push(m);
  }
  for (const m of deduped) if (m.status === 'found') geneOf.get(m.template).found++;
  // Mutations PointFinder lists as a set that acts together (E. faecium's
  // nineteen pbp5 mutations): each knows the set's size, so the report can
  // say how many of it the reads carry.
  const sets = new Map();
  for (const r of cat.rows) {
    if (!TOGETHER_RE.test(r.notes)) continue;
    const k = `${r.template}|${r.notes}`;
    if (!sets.has(k)) sets.set(k, new Set());
    sets.get(k).add(r.pos);
  }
  for (const m of deduped) {
    const note = m.notes.find(n => TOGETHER_RE.test(n));
    if (note && sets.has(`${m.template}|${note}`)) m.together = { note, size: sets.get(`${m.template}|${note}`).size };
  }
  for (const m of deduped.filter(x => x.status === 'needs' || x.status === 'caution')) {
    const why = m.status === 'needs' ? `counted only with ${m.requires}` : 'under development in PointFinder, not in AMRFinderPlus';
    for (const c of calls) if (m.rows.includes(c.row.line)) c.reason = why;
  }

  // Changes at a listed codon that PointFinder's table does not hold (gyrA
  // S83I in K. pneumoniae, where it lists S83L and S83F): one per site, not
  // where a mutation of the table was found, nor at sites whose rows are all
  // to be used with caution, with the antibiotics of the site's rows as
  // context and an AMRFinderPlus symbol for CABBAGE.
  const atSite = new Set(deduped.map(m => `${m.template}:${m.kind}:${m.pos}`));
  const otherChanges = [];
  const seenOther = new Set();
  for (const o of other) {
    if (o.kind !== 'aa_sub') continue;
    const key = `${o.template}:aa_sub:${o.row.pos}`;
    const aa = o.change.slice(o.row.ref.length + String(o.row.pos).length);
    if (atSite.has(key) || seenOther.has(key + aa)) continue;
    seenOther.add(key + aa);
    const siteRows = cat.rows.filter(r => r.template === o.template && r.kind === 'aa_sub' && r.pos === o.row.pos);
    if (siteRows.every(r => r.caution && !r.amrfinder)) continue; // a site the table itself is unsure of
    otherChanges.push({
      gene: o.gene, template: o.template, locus: o.t.locus, region: o.row.region, kind: 'aa_sub', pos: o.row.pos, ref: o.row.ref, alt: aa,
      name: `p.${o.row.ref}${o.row.pos}${aa}`, protein: null, depth: o.depth, share: o.share,
      symbols: o.t.symbolGene && aa !== '*' ? [`${o.t.symbolGene}_${o.row.ref}${o.row.pos}${aa}`] : [],
      listed: uniq(siteRows.flatMap(r => r.alts)), antibiotics: uniq(siteRows.flatMap(r => r.antibiotics)),
      classes: uniq(siteRows.flatMap(r => rowClasses(r, geneClasses))), status: 'other',
    });
  }

  // How many of the species' coding genes the reads match: a check that the
  // species is right (rRNA is too conserved to tell).
  const coding = genes.filter(g => !g.rna);
  const speciesMatch = { matched: coding.filter(g => g.read).length, total: coding.length };

  const order = (a, b) => a.gene.localeCompare(b.gene) || a.pos - b.pos;
  return {
    species: cat.species, label: cat.label,
    genes,
    found: deduped.filter(m => m.status === 'found').sort(order),
    needs: deduped.filter(m => m.status === 'needs').sort(order),
    minority: deduped.filter(m => m.status === 'minority').sort(order),
    caution: deduped.filter(m => m.status === 'caution').sort(order),
    other: other.sort((a, b) => a.t.gene.localeCompare(b.t.gene) || a.row.pos - b.row.pos),
    otherChanges: otherChanges.sort(order),
    undetermined,
    calls,
    screened: calls.length,
    read: calls.length - undetermined.length,
    speciesMatch,
    dropped: cat.dropped?.length || 0,
    nanopore,
  };
}

const uniq = (a) => [...new Set(a)];

const CALL_STATUS = { hit: 'mutation', minority: 'minority', other: 'other', wildtype: 'wildtype', undetermined: 'not read' };
const TOGETHER_RE = /present simultaneously|must be present together/i;

// PointFinder's drug classes as the report's class tokens (genes.js
// CLASS_GROUPS); a class it leaves empty or "unspecified" is read from the
// antibiotics (Cephalosporins, Ethambutol), and a row with neither takes
// the classes of its gene's other rows (S. aureus grlA D79V).
const CLASS_TOKEN = {
  rifampicin: 'rifamycin', carbapenem: 'beta-lactam', cephalosporins: 'beta-lactam', fluoroquinolone: 'quinolone',
};
const classToken = (s) => CLASS_TOKEN[s] || s;
const ownClasses = (row) => {
  const col = String(row.class || '').split(',').map(x => x.trim().toLowerCase()).filter(x => x && x !== 'unspecified');
  return uniq((col.length ? col : row.antibiotics.map(a => a.trim().toLowerCase()).filter(Boolean)).map(classToken));
};
function classesOfGenes(rows) {
  const out = new Map();
  for (const r of rows) {
    const c = ownClasses(r);
    if (!c.length) continue;
    if (!out.has(r.template)) out.set(r.template, new Map());
    const m = out.get(r.template);
    for (const x of c) m.set(x, (m.get(x) || 0) + 1);
  }
  return new Map([...out].map(([t, m]) => [t, [...m].sort((a, b) => b[1] - a[1]).slice(0, 1).map(([x]) => x)]));
}
const rowClasses = (row, byTemplate) => {
  const own = ownClasses(row);
  return own.length ? own : byTemplate.get(row.template) || [];
};

// AMRFinderPlus symbols of an in-frame indel in a coding part, for the
// CABBAGE look-up, at every placement that gives the same protein (in a
// repeat the tools place it differently): an insertion of X after residue
// c is <gene>_<aa c><c><aa c><X> (penA_D346DD), a deletion of the residues
// R from c is <gene>_<R><c>del (rpsE_V27del).
function indelSymbols(row, alt, t, at) {
  if (!t.symbolGene || row.region !== 'coding' || at == null || alt.length % 3 || !['ins', 'del'].includes(row.kind)) return [];
  const rel = at - t.promoter;
  if (rel < 0 || rel % 3) return [];
  const prot = translateSeq(t.seq.slice(t.promoter));
  const n = alt.length / 3;
  const out = [];
  if (row.kind === 'ins') {
    const c = rel / 3;
    const mutant = prot.slice(0, c) + translateSeq(alt) + prot.slice(c);
    for (let k = Math.max(1, c - n - 3); k <= c + n + 3; k++) {
      const x = mutant.slice(k, k + n);
      if (prot.slice(0, k) + x + prot.slice(k) === mutant) out.push(`${t.symbolGene}_${prot[k - 1]}${k}${prot[k - 1]}${x}`);
    }
  } else {
    const c = rel / 3 + 1;
    const mutant = prot.slice(0, c - 1) + prot.slice(c - 1 + n);
    for (let k = Math.max(1, c - n - 3); k <= c + n + 3; k++) {
      if (prot.slice(0, k - 1) + prot.slice(k - 1 + n) === mutant) out.push(`${t.symbolGene}_${prot.slice(k - 1, k - 1 + n)}${k}del`);
    }
  }
  return uniq(out);
}

const translateSeq = (s) => {
  let p = '';
  for (let k = 0; k + 3 <= s.length; k += 3) p += translate(s.slice(k, k + 3));
  return p;
};

// A mutation as the report and tables name it: the gene (an rRNA copy by
// its locus, a promoter change as "ampC promoter", a coding change on a
// promoter template by the gene) and the change, as a protein change when
// there is one.
// The changes PointFinder lists at the codon of another change: "S83L and
// S83F".
export function listedAt(m) {
  const all = (m.listed || []).map(a => `${m.ref}${m.pos}${a}`);
  return all.length > 1 ? `${all.slice(0, -1).join(', ')} and ${all[all.length - 1]}` : all[0] || '';
}

export function mutationLabel(m) {
  const gene = m.region === 'rna' && m.locus ? m.locus
    : m.region === 'coding' ? m.gene.replace(/-promoter$/, '') // katG R463L on the katG-promoter template
      : m.gene.replace(/-promoter$/, ' promoter');
  return `${gene} ${m.protein || m.name}`;
}

// A row's required mutations in words: "gyrA S83L, W, A or V", groups
// joined by "or", members of a group by "and".
function requirementText(requires, tmpl) {
  if (!requires?.length) return '';
  const alts = (a) => (a.length > 1 ? `${a.slice(0, -1).join(', ')} or ${a[a.length - 1]}` : a[0]);
  const one = (m) => `${tmpl.get(m.template)?.gene || m.template} ${m.ref}${m.pos}${alts(m.alts)}`;
  return requires.map(group => group.map(one).join(' and ')).join(', or ');
}

// A site of an rRNA copy in the coordinates of the locus's first copy, so
// that the same site reported on several copies is recognised.
function canonicalSite(m, t, copies) {
  const first = copies[0];
  const i = tindex(t, m.pos);
  if (first === t) return i;
  const map = alignTemplates(first.seq, t.seq)?.map;
  const j = map ? map.indexOf(i) : -1;
  return j < 0 ? `${t.name}:${m.pos}` : j;
}

// Reads one row. Returns { status: 'hit' | 'minority' | 'wildtype' |
// 'other' | 'undetermined', alt, depth, share, aa, codon, reason }.
function evaluateRow(row, t, view, minDepth, mode) {
  switch (row.kind) {
    case 'aa_sub': return evalCodon(row, t, view, minDepth, mode);
    case 'nuc_sub': return evalBase(row, t, view, minDepth, mode);
    case 'stop_any': return evalStop(row, t, view, minDepth);
    case 'del': case 'ins': return evalIndel(row, t, view, minDepth);
    case 'aa_ins': case 'aa_del': return evalAaIndel(row, t, view, minDepth);
    default: return { status: 'undetermined', reason: 'unknown row kind' };
  }
}

function evalCodon(row, t, view, minDepth, mode) {
  const i = t.promoter + 3 * (row.pos - 1);
  const calls = [0, 1, 2].map(k => callAt(view, i + k, minDepth));
  const depth = Math.min(...calls.map(c => c.depth || 0));
  if (calls.some(c => !c.base)) return { status: 'undetermined', reason: calls.find(c => !c.base).reason, depth };
  if (calls.some(c => c.base === '-') || [0, 1].some(k => insertedAfter(view, i + k))) {
    return { status: 'undetermined', reason: 'indel in the codon', depth };
  }
  const codon = calls.map(c => c.base).join('');
  const aa = translate(codon);
  const refCodon = t.seq.slice(i, i + 3);
  const share = Math.min(...calls.map((c, k) => (c.base === refCodon[k] ? 1 : c.share)));
  const read = `${aa} (${codon.toLowerCase()})`;
  if (aa !== row.ref && row.alts.includes(aa)) return { status: 'hit', alt: aa, depth, share, aa, read, codon: `${refCodon.toLowerCase()}>${codon.toLowerCase()}` };
  // A minority of the reads with an allowed change, read base by base (a
  // single-base change; reads are not phased within a codon in the matrix),
  // where most reads carry the reference codon: a base the reads share with
  // another codon (ATC, Ile, against TCC) is not a minority of TTC (Phe).
  const minority = codon === refCodon ? minorityCodon(row, t, view, i, refCodon, mode) : null;
  if (minority) return { status: 'minority', ...minority, aa, read };
  if (aa !== row.ref) return { status: 'other', change: `${row.ref}${row.pos}${aa}`, depth, share, aa, read };
  return { status: 'wildtype', depth, share, aa, read };
}

function minorityCodon(row, t, view, i, refCodon, mode) {
  for (let k = 0; k < 3; k++) {
    for (const b of BASES) {
      if (b === refCodon[k]) continue;
      const codon = refCodon.slice(0, k) + b + refCodon.slice(k + 1);
      const aa = translate(codon);
      if (!row.alts.includes(aa)) continue;
      const n = countOf(view, i + k, b), d = depthAt(view, i + k);
      if (n >= POINT_CALL.minorityReads && d && n / d >= POINT_CALL.minorityShare[mode] && n / d < 0.5) {
        return { alt: aa, depth: d, share: n / d, codon: `${refCodon.toLowerCase()}>${codon.toLowerCase()}` };
      }
    }
  }
  return null;
}

function evalBase(row, t, view, minDepth, mode) {
  const i = tindex(t, row.pos);
  const c = callAt(view, i, minDepth);
  const d = depthAt(view, i);
  const read = c.base === '-' ? 'deleted' : c.base || '';
  if (c.base && row.alts.includes(c.base)) return { status: 'hit', alt: c.base, depth: c.depth, share: c.share, read };
  // Any allowed base on a minority of the reads (mixed rRNA copies).
  const minShare = (t.rna ? POINT_CALL.minorityShareRna : POINT_CALL.minorityShare)[mode];
  if (d >= minDepth) {
    for (const b of row.alts) {
      const n = countOf(view, i, b);
      if (n >= POINT_CALL.minorityReads && n / d >= minShare && n / d < 0.5) return { status: 'minority', alt: b, depth: d, share: n / d, read };
    }
  }
  if (!c.base) return { status: 'undetermined', reason: c.reason, depth: d };
  if (c.base !== row.ref) return { status: 'other', change: `${row.ref}${row.pos}${c.base === '-' ? 'del' : c.base}`, depth: c.depth, share: c.share, read };
  return { status: 'wildtype', depth: c.depth, share: c.share, read };
}

// "Any premature stop in this gene" (M. tuberculosis katG, ethA, pncA, gidB):
// translate the coding part as the reads have it.
function evalStop(row, t, view, minDepth) {
  const start = t.promoter;
  const seq = sampleSeq(view, start, t.length, minDepth);
  if (seq == null) return { status: 'undetermined', reason: 'part of the gene not read' };
  const codons = Math.floor((t.length - start) / 3) - 1; // without the template's stop
  for (let c = 0; c + 3 <= seq.length && c / 3 < codons; c += 3) {
    if (STOPS.has(seq.slice(c, c + 3))) {
      const n = c / 3 + 1;
      const refAa = translate(t.seq.slice(start + c, start + c + 3));
      return { status: 'hit', alt: `${refAa}${n}*`, depth: depthAt(view, start + c), share: callAt(view, start + c, minDepth).share || 1, read: `stop at codon ${n}` };
    }
  }
  return { status: 'wildtype', read: 'no stop' };
}

// An indel row matches when the sample's sequence around the site equals the
// reference with the indel applied, at any placement that gives the same
// sequence (repeats). WINDOW bases are compared on each side.
const WINDOW = 25;

function evalIndel(row, t, view, minDepth) {
  const i = tindex(t, row.pos);
  const from = i - WINDOW, to = i + WINDOW + Math.max(...row.alts.map(a => a.length), 1) + 1;
  const off = Math.max(0, from);
  const ref = t.seq.slice(off, Math.min(t.length, to));
  const mutants = [];
  for (const alt of row.alts) {
    for (const shift of [0, -1, 1]) {
      if (row.kind === 'del') {
        const s = i + shift - off;
        if (ref.slice(s, s + alt.length) !== alt) continue;
        mutants.push({ alt, at: s + off, seq: ref.slice(0, s) + ref.slice(s + alt.length) });
      } else {
        const s = i + 1 + shift - off; // after the stated base
        mutants.push({ alt, at: s + off, seq: ref.slice(0, s) + alt + ref.slice(s) });
      }
    }
  }
  // The whole template with the change applied identifies the event, so
  // rows that write one insertion at several places of a repeat match the
  // same one.
  const whole = (m) => t.seq.slice(0, off) + m.seq + t.seq.slice(off + ref.length);
  const sample = sampleSeq(view, from, to, minDepth);
  if (sample != null) {
    const m = mutants.find(x => x.seq === sample);
    if (m) return { status: 'hit', alt: m.alt, depth: depthAt(view, i), share: indelSupport(row, view, i, m.alt), haplotype: whole(m), at: m.at, read: `${row.kind === 'del' ? 'deleted' : 'inserted'} ${m.alt}` };
  }
  // An insertion that the reads place at different spots of a repeat has
  // no majority at any one spot, so the consensus misses it. Each spot's
  // inserted bases, applied to the reference, are checked against the
  // row's sequence; the reads of the spots that give it are added up.
  if (row.kind === 'ins') {
    const found = splitInsertion(view, t, ref, off, i, mutants, minDepth);
    if (found) return { ...found, haplotype: t.seq.slice(0, off) + found.haplotype + t.seq.slice(off + ref.length), read: `inserted ${found.alt}, placed differently by the reads` };
  }
  if (sample == null) return { status: 'undetermined', reason: 'site not read', depth: depthAt(view, i) };
  return sample === ref ? { status: 'wildtype', depth: depthAt(view, i), read: 'as the reference' }
    : { status: 'other', change: 'another change within 25 bases', depth: depthAt(view, i), read: 'differs near the site' };
}

function splitInsertion(view, t, ref, off, i, mutants, minDepth) {
  const span = Math.max(...mutants.map(m => m.alt.length), 1) + 3;
  const depth = depthAt(view, i);
  if (depth < minDepth) return null;
  let support = 0, best = null;
  for (let k = i - span; k <= i + span; k++) {
    const rows = view.ins?.get(k);
    if (!rows?.length) continue;
    // the inserted string at this spot: each column's most common base,
    // with the reads that carry the whole string (the smallest count)
    let s = '', n = Infinity;
    for (const r of rows) {
      let b = 0;
      for (let j = 1; j < 4; j++) if (r[j] > r[b]) b = j;
      if (!r[b]) { s = ''; break; }
      s += BASES[b];
      n = Math.min(n, r[b]);
    }
    if (!s) continue;
    const at = k + 1 - off;
    const seq = ref.slice(0, at) + s + ref.slice(at);
    const m = mutants.find(x => x.seq === seq);
    if (!m) continue;
    support += n;
    best = best || { m, seq };
  }
  if (!best || support / depth < 0.5) return null;
  return { status: 'hit', alt: best.m.alt, depth, share: Math.min(1, support / depth), haplotype: best.seq, at: best.m.at };
}

// The share of reads with the indel: for a deletion, the gap counts over the
// deleted bases; for an insertion, the reads with inserted bases at the site
// and its neighbours (reads place an insertion in a repeat differently).
function indelSupport(row, view, i, alt) {
  const d = depthAt(view, i) || 1;
  if (row.kind === 'del') {
    let gaps = 0;
    for (let k = -1; k <= alt.length; k++) gaps += countOf(view, i + k, '-');
    return Math.min(1, gaps / alt.length / d);
  }
  let inserted = 0;
  for (let k = i - alt.length - 2; k <= i + alt.length + 2; k++) {
    for (const r of view.ins?.get(k) || []) inserted += r[0] + r[1] + r[2] + r[3] + r[4];
  }
  return Math.min(1, inserted / alt.length / d);
}

// Amino-acid insertions (ompK36 131 ins SG, porA 139 ins D): translate the
// sample around the codon and look for the inserted residues after codon
// pos-1 or pos.
function evalAaIndel(row, t, view, minDepth) {
  const c0 = t.promoter + 3 * (row.pos - 1);
  const from = c0 - 3 * 8, to = c0 + 3 * 8 + 3 * Math.max(...row.alts.map(a => a.length), 1);
  const sample = sampleSeq(view, from, to, minDepth);
  if (sample == null) return { status: 'undetermined', reason: 'site not read', depth: depthAt(view, c0) };
  const ref = t.seq.slice(Math.max(0, from), Math.min(t.length, to));
  if ((sample.length - ref.length) % 3 !== 0) return { status: 'other', change: 'a frameshift near the site', depth: depthAt(view, c0), read: 'frameshift' };
  const frame = (Math.max(0, from) - t.promoter) % 3;
  const prot = (s) => { let p = ''; for (let k = (3 - frame) % 3; k + 3 <= s.length; k += 3) p += translate(s.slice(k, k + 3)); return p; };
  const pr = prot(ref), ps = prot(sample);
  const at = Math.floor((c0 - Math.max(0, from) - ((3 - frame) % 3)) / 3);
  for (const alt of row.alts) {
    if (row.kind === 'aa_ins') {
      for (const k of [at, at + 1]) {
        if (pr.slice(0, k) + alt + pr.slice(k) === ps) return { status: 'hit', alt, depth: depthAt(view, c0), share: 1, read: `inserted ${alt}` };
      }
    } else if (alt === 'DEL') {
      if (pr.slice(0, at) + pr.slice(at + 1) === ps) return { status: 'hit', alt: 'del', depth: depthAt(view, c0), share: 1, read: 'codon deleted' };
    }
  }
  return ps === pr ? { status: 'wildtype', depth: depthAt(view, c0), read: 'as the reference' }
    : { status: 'other', change: 'another change within 8 codons', depth: depthAt(view, c0), read: 'differs near the site' };
}

// ResFinder's notation: p.S83L, r.2058A>G, n.-42C>T, p.V27del, n.-57delA,
// n.-15_-14insGT, p.131insSG; a stop anywhere in the gene, p.Q30*.
export function mutationName(row, alt, t) {
  const nuc = row.region === 'rna' ? 'r' : 'n';
  switch (row.kind) {
    case 'aa_sub': return `p.${row.ref}${row.pos}${alt}`;
    case 'nuc_sub': return `${nuc}.${row.pos}${row.ref}>${alt}`;
    case 'stop_any': return `p.${alt}`;
    case 'del': return `${nuc}.${row.pos}${alt.length > 1 ? `_${row.pos + alt.length - 1}` : ''}del${alt}`;
    case 'ins': {
      const next = row.pos === -1 ? 1 : row.pos + 1;
      return `${nuc}.${row.pos}_${next}ins${alt}`;
    }
    case 'aa_ins': return `p.${row.pos}ins${alt}`;
    case 'aa_del': return `p.${translate(t.seq.slice(t.promoter + 3 * (row.pos - 1), t.promoter + 3 * row.pos))}${row.pos}del`;
    default: return `${row.pos}${alt}`;
  }
}

// The protein-level name of an in-frame indel in a coding part, when it sits
// between whole codons (rpsE n.78_80delGTT is p.V27del; penA's GAC after
// codon 345 is p.R345_G346insD), else null. at: template index where the
// deleted bases start or the inserted ones go.
function proteinName(row, alt, t, at) {
  if (at == null || row.region !== 'coding' || alt.length % 3) return null;
  const rel = at - t.promoter;
  if (rel < 0 || rel % 3) return null;
  const aaAt = (c) => translate(t.seq.slice(t.promoter + 3 * (c - 1), t.promoter + 3 * c));
  if (row.kind === 'del') {
    const c1 = rel / 3 + 1, c2 = c1 + alt.length / 3 - 1;
    return c1 === c2 ? `p.${aaAt(c1)}${c1}del` : `p.${aaAt(c1)}${c1}_${aaAt(c2)}${c2}del`;
  }
  if (row.kind === 'ins') {
    const c = rel / 3; // inserted after codon c
    let ins = '';
    for (let k = 0; k + 3 <= alt.length; k += 3) ins += translate(alt.slice(k, k + 3));
    return c >= 1 ? `p.${aaAt(c)}${c}_${aaAt(c + 1)}${c + 1}ins${ins}` : null;
  }
  return null;
}

// ── Output ──

// A row's site, without the change: p.S83, r.2058A, n.-42C, n.-13 (indel).
export function siteName(row) {
  const nuc = row.region === 'rna' ? 'r' : 'n';
  switch (row.kind) {
    case 'aa_sub': return `p.${row.ref}${row.pos}`;
    case 'nuc_sub': return `${nuc}.${row.pos}${row.ref}`;
    case 'stop_any': return 'any stop codon';
    case 'aa_ins': case 'aa_del': return `p.${row.pos}`;
    default: return `${nuc}.${row.pos}`;
  }
}

// pointfinder_results.tsv in the ZIP: the mutations found, in some of the
// reads, without the mutation PointFinder requires with them, or to use with
// caution, and the other changes at listed codons, in the columns of
// ResFinder's PointFinder_results.txt and more.
export function pointfinderTsv(ev) {
  const lines = [['gene', 'template', 'mutation', 'nucleotide_change', 'amino_acid_change', 'status', 'depth', 'read_share',
    'resistance', 'classes', 'requires', 'pmids', 'notes'].join('\t')];
  for (const m of [...ev.found, ...ev.needs, ...ev.minority, ...(ev.caution || []), ...(ev.otherChanges || [])]) {
    const aa = m.kind === 'aa_sub' ? `${m.ref}>${m.alt}` : m.protein || '';
    const notes = m.status === 'other' ? [`not in the table, which lists ${listedAt(m)} at this codon`] : m.notes;
    lines.push([m.gene, m.template, m.name, m.codon || '', aa, m.status, m.depth, share(m.share),
      m.antibiotics.join(', '), m.classes.join(', '), m.requires || '', (m.pmids || []).join(', '), notes.join(' ')].map(tsvCell).join('\t'));
  }
  return lines.join('\n') + '\n';
}

// pointfinder_positions.tsv in the ZIP: every row of PointFinder's table for
// the species as the reads have it, with the reason when it could not be read
// (a position the reads do not reach is not wild type).
export function pointfinderPositionsTsv(ev) {
  const lines = [['line', 'gene', 'template', 'site', 'wild_type', 'resistant', 'reads_show', 'depth', 'read_share', 'result', 'reason', 'resistance'].join('\t')];
  for (const c of [...ev.calls].sort((a, b) => a.row.line - b.row.line)) {
    const r = c.row;
    lines.push([r.line, r.gene, r.template, siteName(r), r.ref, r.alts.join(','), c.read || '', c.depth ?? '', share(c.share),
      c.status, c.reason || '', r.antibiotics.join(', ')].map(tsvCell).join('\t'));
  }
  return lines.join('\n') + '\n';
}

const share = (v) => (v != null ? Number(v).toFixed(3) : '');
const tsvCell = (v) => String(v ?? '').replace(/[\t\r\n]+/g, ' ');
