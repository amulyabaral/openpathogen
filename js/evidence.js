/* evidence.js — what the reads say about a hit, beyond KMA's .res numbers.
 *
 * Every .res number describes KMA's consensus rather than the reads under
 * it: Template_Coverage counts inserted bases as coverage, the identity of a
 * gene cut short includes the ragged ends of the reads that cross the cut,
 * and when the sample's allele is not in the database KMA can share its
 * reads between two sister templates, so one gene shows up as two partly
 * covered hits. KMA's per-base counts (.mat.gz: A, C, G, T, N and gaps under
 * every template position) and its alignment (.aln) answer what the .res
 * cannot:
 *
 *   coverage     which positions of the gene have reads at all. A deletion
 *                in the sample still has reads across it, so only a stretch
 *                that no read reaches counts as uncovered.
 *   gaps         whether thin coverage can explain such a stretch. Under
 *                Poisson read starts (Lander–Waterman) a stretch of l bases
 *                gets no read with probability exp(-(d/r)(l + r - 1)) at
 *                depth d and read length r. Against the gene's own depth, a
 *                long stretch at high depth cannot be missed by chance: that
 *                part of the gene is missing from the sample, or too
 *                different from the template to align.
 *   sisters      whether two templates of one gene are one allele whose
 *                reads KMA split, or two alleles. At the positions where the
 *                templates differ, the pooled reads of one allele show one
 *                base and those of two alleles show both. Split reads are
 *                pooled into one gene, named by its closest template.
 *   translation  the consensus translated in frame, for a stop codon or a
 *                frameshift in a gene whose coverage looks complete.
 *   mixed bases  positions where the reads of one template disagree beyond
 *                sequencing error: a second, related allele or a mixed sample.
 *
 * Everything is computed from the KMA output already in memory, lazily and
 * once per result. Without a matrix (an older result, a failed run) the
 * functions return null and the callers keep KMA's numbers.
 */

import { parseAlnBlock, alnBlocksOf, alignmentCheck } from './aln-check.js';

export const EVIDENCE = {
  // Depth is divided by this before the gap test, an allowance for coverage
  // that is less even than Poisson (GC bias, library effects).
  dispersion: 2,
  // A stretch with no reads is "missing" when its chance under thin
  // coverage, times the number of places it could start, is below this, and
  // when it is at least minMissing bases or minMissingShare of the gene long
  // (a base or two at a gene's end is usually where a template's start or end
  // is annotated differently, not a missing part).
  missingP: 1e-3,
  minMissing: 10,
  minMissingShare: 0.01,
  // KMA -mrs under -cge: a read needs half its length on the template, so
  // reads hanging over a gene's ends reach less far into it.
  mrs: 0.5,
  // Read length when .mapstat gives none (Illumina).
  defaultReadLength: 150,
  // A site counts for the sister and mixed-base tests from this depth.
  minSiteDepth: 8,
  // Sequencing error per read and base when the template's own reads do not
  // give a usable estimate, and the floor of that estimate.
  errorFloor: { illumina: 0.002, nanopore: 0.03 },
  // A site is mixed when the minor base is at least this share of the reads
  // and its binomial tail against the error rate is below mixedP.
  mixedShare: { illumina: 0.05, nanopore: 0.2 },
  mixedP: 1e-4,
  // Templates align as sisters from this identity.
  sisterIdentity: 0.8,
  // Two sets of reads agree on a base when this share of each carries it
  // (the sister test's second-allele check).
  agree: 0.8,
  // Translation calls need this depth and this share of reads at the bases
  // that make the stop codon or the frameshift.
  stopDepth: 5,
  stopShare: 0.8,
  // A protein cut to less than this share of its length is likely inactive
  // (Kleborate's threshold for truncated genes).
  truncatedProtein: 0.9,
};

const BASES = 'ACGT';
const COL = { A: 0, C: 1, G: 2, T: 3, N: 4, '-': 5 };
const STOPS = new Set(['TAA', 'TAG', 'TGA']);
const STARTS = new Set(['ATG', 'GTG', 'TTG']);

// ── KMA output text ──

function fflate() {
  const f = globalThis.fflate;
  if (!f) throw new Error('fflate is not loaded');
  return f;
}

function fileText(result, ext) {
  const info = result?.files?.[ext];
  if (!info) return null;
  if (typeof info.data === 'string') return info.data;
  const u8 = info.data instanceof Uint8Array ? info.data : new Uint8Array(info.data);
  return fflate().strFromU8(ext.endsWith('.gz') ? fflate().gunzipSync(u8) : u8);
}

// Offsets of each "# name" block of a .mat or .aln text, so a template's
// block is parsed only when it is needed.
function indexBlocks(text) {
  const out = new Map();
  let at = text.startsWith('#') ? 0 : text.indexOf('\n#');
  if (at > 0) at += 1;
  while (at >= 0 && at < text.length) {
    const eol = text.indexOf('\n', at);
    const name = text.slice(at + 1, eol < 0 ? text.length : eol).trim();
    const next = eol < 0 ? -1 : text.indexOf('\n#', eol);
    out.set(name, [eol < 0 ? text.length : eol + 1, next < 0 ? text.length : next + 1]);
    at = next < 0 ? -1 : next + 1;
  }
  return out;
}

// One template's block of the matrix: n template positions, counts[i*6+k]
// for k = A C G T N -, and the insertion rows after each position
// (ins.get(i), i = -1 before the first). Rows whose reference is "-" are
// insertions: their A C G T count reads carrying an inserted base there, and
// their "-" reads that pass without one.
export function parseMatrixBlock(body) {
  const ref = [];
  const counts = [];
  const ins = new Map();
  for (const line of body.split('\n')) {
    if (!line || line[0] === '#') continue;
    const p = line.split('\t');
    if (p.length < 7) continue;
    const c = [+p[1] || 0, +p[2] || 0, +p[3] || 0, +p[4] || 0, +p[5] || 0, +p[6] || 0];
    if (p[0] === '-') {
      const at = ref.length - 1;
      if (!ins.has(at)) ins.set(at, []);
      ins.get(at).push(c);
    } else {
      ref.push(p[0].toUpperCase());
      counts.push(...c);
    }
  }
  return { ref: ref.join(''), n: ref.length, counts: Uint32Array.from(counts), ins };
}

function parseMapstat(text) {
  const out = new Map();
  if (!text) return out;
  let cols = null;
  for (const line of text.split('\n')) {
    if (!line) continue;
    if (line.startsWith('# refSequence')) { cols = line.slice(2).split('\t'); continue; }
    if (line[0] === '#' || !cols) continue;
    const p = line.split('\t');
    const row = {};
    cols.forEach((c, i) => { row[c] = p[i]; });
    out.set(String(p[0]).trim(), row);
  }
  return out;
}

// Per result: the matrix index, parsed blocks, .mapstat rows and the read
// type, built on first use.
function sourceOf(result) {
  if (result._evidence) return result._evidence;
  const src = { mat: null, matIndex: null, blocks: new Map(), mapstat: new Map(), aln: null, cache: new Map() };
  try {
    const text = fileText(result, '.mat.gz');
    if (text) { src.mat = text; src.matIndex = indexBlocks(text); }
  } catch (_) { /* no matrix: evidence is unavailable */ }
  try { src.mapstat = parseMapstat(fileText(result, '.mapstat')); } catch (_) {}
  try { src.aln = alnBlocksOf(result); } catch (_) { src.aln = new Map(); }
  src.nanopore = /(^|\s)-ont(\s|$)|-bcNano/.test(result.command || '');
  result._evidence = src;
  return src;
}

export function matrixOf(result, template) {
  const src = sourceOf(result);
  const key = String(template).trim();
  if (src.blocks.has(key)) return src.blocks.get(key);
  let m = null;
  const span = src.matIndex?.get(key);
  if (span) {
    try { m = parseMatrixBlock(src.mat.slice(span[0], span[1])); } catch (_) { m = null; }
    if (m && !m.n) m = null;
  }
  src.blocks.set(key, m);
  return m;
}

export function hasEvidence(result) {
  return !!sourceOf(result).matIndex;
}

// ── Per-position summaries ──

// Reads over each template position (bases, Ns and deletions), the
// positions some read reaches, and the stretches no read reaches (1-based,
// inclusive; edge marks a stretch at the start or end of the template).
export function depthOf(m) {
  const depth = new Uint32Array(m.n);
  for (let i = 0; i < m.n; i++) {
    let d = 0;
    for (let k = 0; k < 6; k++) d += m.counts[i * 6 + k];
    depth[i] = d;
  }
  return depth;
}

export function gapsOf(depth) {
  const gaps = [];
  let start = -1;
  for (let i = 0; i <= depth.length; i++) {
    const empty = i < depth.length && depth[i] === 0;
    if (empty && start < 0) start = i;
    if (!empty && start >= 0) {
      gaps.push({ start: start + 1, end: i, length: i - start, edge: start === 0 ? 'start' : i === depth.length ? 'end' : null });
      start = -1;
    }
  }
  return gaps;
}

function medianOf(values) {
  if (!values.length) return 0;
  const v = Float64Array.from(values).sort();
  const h = v.length >> 1;
  return v.length % 2 ? v[h] : (v[h - 1] + v[h]) / 2;
}

function median(depth) {
  const v = [];
  for (let i = 0; i < depth.length; i++) if (depth[i]) v.push(depth[i]);
  return medianOf(v);
}

// log10 of the chance that a stretch like this has no read at the given
// depth, times the n places it could start (so a gene with many positions
// does not produce "significant" gaps by sheer length).
export function gapLog10P(gap, n, depth, readLength, o = EVIDENCE) {
  const d = depth / o.dispersion;
  const r = Math.max(1, readLength);
  if (!(d > 0)) return 0;
  const starts = gap.edge ? gap.length + (1 - o.mrs) * r : gap.length + r - 1;
  return Math.min(0, (-(d / r) * starts) / Math.LN10 + Math.log10(Math.max(1, n)));
}

// Mean aligned read length on the template from .mapstat, or a default.
function readLengthOf(src, templates, n) {
  let bp = 0, reads = 0;
  for (const t of templates) {
    const row = src.mapstat.get(String(t).trim());
    if (!row) continue;
    bp += Number(row.bpTotal) || 0;
    reads += Number(row.readCount) || 0;
  }
  const r = reads ? bp / reads : (src.nanopore ? n : EVIDENCE.defaultReadLength);
  return Math.max(30, Math.min(n, r));
}

// The consensus at each template position from the .aln: base[i] is the
// consensus character under template position i ("-" when it has none) and
// insert[i] the bases KMA placed after it (insert[-1]: before the first).
export function consensusColumns(alnBlock) {
  const { t, q } = parseAlnBlock(alnBlock);
  if (!t.length || t.length !== q.length) return null;
  const base = [];
  const insert = new Map();
  for (let c = 0; c < t.length; c++) {
    if (t[c] === '-') {
      if (q[c] !== '-') {
        const at = base.length - 1;
        insert.set(at, (insert.get(at) || '') + q[c]);
      }
    } else {
      base.push(q[c]);
    }
  }
  return { base, insert, n: base.length };
}

// ── Sequencing error and mixed bases ──

const binomTailCache = new Map();

// P(X >= k) for X ~ Binomial(n, p), in log10 (exact sum, k small).
function log10BinomTail(k, n, p) {
  if (k <= 0) return 0;
  if (p <= 0) return -Infinity;
  const key = `${k},${n},${p.toFixed(6)}`;
  if (binomTailCache.has(key)) return binomTailCache.get(key);
  let logTerm = 0;
  // log C(n,k) p^k (1-p)^(n-k), then the ratio of successive terms.
  for (let i = 0; i < k; i++) logTerm += Math.log(n - i) - Math.log(i + 1);
  logTerm += k * Math.log(p) + (n - k) * Math.log1p(-p);
  let sum = 0;
  let term = 0; // relative to logTerm
  for (let x = k; x <= n; x++) {
    sum += Math.exp(term);
    term += Math.log(n - x) - Math.log(x + 1) + Math.log(p) - Math.log1p(-p);
    if (term < -40) break;
  }
  const out = (logTerm + Math.log(sum)) / Math.LN10;
  if (binomTailCache.size > 5000) binomTailCache.clear();
  binomTailCache.set(key, out);
  return out;
}

// Per-base disagreement among the reads of a template: the share of reads
// whose base is not the top base, over positions with enough depth. It
// includes real minor alleles, which only makes the mixed-base test stricter.
function errorRate(counts, n, skip, nanopore) {
  let minor = 0, total = 0;
  for (let i = 0; i < n; i++) {
    if (skip?.has(i)) continue;
    let top = 0, sum = 0;
    for (let k = 0; k < 4; k++) { const c = counts[i * 6 + k]; sum += c; if (c > top) top = c; }
    if (sum < EVIDENCE.minSiteDepth) continue;
    minor += sum - top;
    total += sum;
  }
  const floor = EVIDENCE.errorFloor[nanopore ? 'nanopore' : 'illumina'];
  return Math.max(floor, total ? minor / total : floor);
}

// A site's two most common bases and whether the second is beyond error.
function siteMix(c, err, nanopore) {
  const order = [0, 1, 2, 3].sort((a, b) => c[b] - c[a]);
  const top = c[order[0]], second = c[order[1]];
  const depth = c[0] + c[1] + c[2] + c[3];
  const share = depth ? second / depth : 0;
  const mixed = depth >= EVIDENCE.minSiteDepth && second >= 3
    && share >= EVIDENCE.mixedShare[nanopore ? 'nanopore' : 'illumina']
    && log10BinomTail(second, depth, err / 3) < Math.log10(EVIDENCE.mixedP);
  return { top: BASES[order[0]], second: BASES[order[1]], topN: top, secondN: second, depth, share, mixed };
}

// Positions of one template whose reads carry a second base beyond error,
// over the whole gene (so the tail probability is divided by its length),
// leaving out the positions in skip (next to stretches with no reads or in
// poorly aligned windows, where the reads crossing a junction disagree).
export function mixedSites(m, nanopore, skip = null) {
  const err = errorRate(m.counts, m.n, skip, nanopore);
  const out = [];
  const limit = Math.log10(EVIDENCE.mixedP / Math.max(1, m.n));
  const minShare = EVIDENCE.mixedShare[nanopore ? 'nanopore' : 'illumina'];
  for (let i = 0; i < m.n; i++) {
    if (skip?.has(i)) continue;
    const c = [m.counts[i * 6], m.counts[i * 6 + 1], m.counts[i * 6 + 2], m.counts[i * 6 + 3]];
    const depth = c[0] + c[1] + c[2] + c[3];
    if (depth < EVIDENCE.minSiteDepth) continue;
    const s = siteMix(c, err, nanopore);
    if (s.secondN >= 3 && s.share >= minShare && log10BinomTail(s.secondN, depth, err / 3) < limit) {
      out.push({ pos: i + 1, ...s, mixed: true });
    }
  }
  return { sites: out, error: err };
}

// ── Aligning two templates ──

// Global alignment of two template sequences within a band (sister alleles
// differ by substitutions and a few indels). Returns map[i] = position in b
// aligned to position i of a, or -1, and the identity over aligned columns.
export function alignTemplates(a, b) {
  if (a.length === b.length) {
    let same = 0;
    for (let i = 0; i < a.length; i++) if (a[i] === b[i]) same++;
    if (same / a.length >= EVIDENCE.sisterIdentity) {
      return { map: Int32Array.from({ length: a.length }, (_, i) => i), identity: same / a.length, aligned: a.length };
    }
  }
  const n = a.length, m = b.length;
  const band = Math.abs(n - m) + 64;
  if (n * (2 * band + 1) > 4e7) return null;
  const W = 2 * band + 1;
  const NEG = -1e9;
  const score = new Float64Array((n + 1) * W).fill(NEG);
  const move = new Uint8Array((n + 1) * W); // 1 diagonal, 2 up (gap in b), 3 left (gap in a)
  const diag = (i) => Math.round(i * m / Math.max(1, n));
  const at = (i, j) => { const k = j - diag(i) + band; return k >= 0 && k < W ? i * W + k : -1; };
  const GAP = -2;
  for (let i = 0; i <= n; i++) {
    const c = diag(i);
    for (let j = Math.max(0, c - band); j <= Math.min(m, c + band); j++) {
      const k = at(i, j);
      if (i === 0 && j === 0) { score[k] = 0; continue; }
      let best = NEG, mv = 0;
      if (i > 0 && j > 0) {
        const p = at(i - 1, j - 1);
        if (p >= 0 && score[p] > NEG) { const s = score[p] + (a[i - 1] === b[j - 1] ? 1 : -1); if (s > best) { best = s; mv = 1; } }
      }
      if (i > 0) {
        const p = at(i - 1, j);
        if (p >= 0 && score[p] > NEG) { const s = score[p] + GAP; if (s > best) { best = s; mv = 2; } }
      }
      if (j > 0) {
        const p = at(i, j - 1);
        if (p >= 0 && score[p] > NEG) { const s = score[p] + GAP; if (s > best) { best = s; mv = 3; } }
      }
      score[k] = best;
      move[k] = mv;
    }
  }
  if (at(n, m) < 0 || score[at(n, m)] <= NEG) return null;
  const map = new Int32Array(n).fill(-1);
  let i = n, j = m, same = 0, cols = 0;
  while (i > 0 || j > 0) {
    const mv = move[at(i, j)];
    if (mv === 1) { map[i - 1] = j - 1; if (a[i - 1] === b[j - 1]) same++; i--; j--; }
    else if (mv === 2) i--;
    else j--;
    cols++;
  }
  return { map, identity: cols ? same / cols : 0, aligned: map.reduce((x, v) => x + (v >= 0 ? 1 : 0), 0) };
}

// ── One gene: a template, or sister templates whose reads are pooled ──

// Pools the matrices of the members (the anchor first) in the anchor's
// coordinates. maps[k][i] = member k's position under anchor position i.
function pooledCounts(members, maps) {
  const n = members[0].m.n;
  const counts = new Uint32Array(n * 6);
  for (let k = 0; k < members.length; k++) {
    const { m } = members[k];
    const map = maps[k];
    for (let i = 0; i < n; i++) {
      const j = map ? map[i] : i;
      if (j < 0) continue;
      for (let c = 0; c < 6; c++) counts[i * 6 + c] += m.counts[j * 6 + c];
    }
  }
  return counts;
}

// The pooled consensus in anchor coordinates: at each position, the member
// with the most reads there gives its KMA consensus base (and the bases KMA
// inserted after it), so the bases are KMA's own calls wherever one exists.
function pooledConsensus(members, maps, alnOf) {
  const n = members[0].m.n;
  const cols = members.map(x => consensusColumns(alnOf(x.template)));
  const base = new Array(n).fill('-');
  const insert = new Map();
  for (let i = 0; i < n; i++) {
    let best = -1, bestDepth = -1;
    for (let k = 0; k < members.length; k++) {
      const j = maps[k] ? maps[k][i] : i;
      if (j < 0 || !cols[k]) continue;
      let d = 0;
      for (let c = 0; c < 6; c++) d += members[k].m.counts[j * 6 + c];
      if (d > bestDepth) { bestDepth = d; best = k; }
    }
    if (best < 0 || bestDepth <= 0) continue;
    const j = maps[best] ? maps[best][i] : i;
    base[i] = cols[best].base[j] ?? '-';
    const ins = cols[best].insert.get(j);
    if (ins) insert.set(i, ins);
  }
  const lead = cols[0]?.insert.get(-1);
  if (lead) insert.set(-1, lead);
  return { base, insert, n };
}

// The reads that cross the edge of a missing stretch (the junction of a
// truncation) are aligned base by base into it, so next to a stretch with no
// reads the consensus is junk: mismatches and inserted bases. From the edge,
// walk into the covered part until the consensus matches the template at 11
// of 12 positions; the positions walked over belong to the missing stretch.
// Returns how far the stretch extends to the left and to the right.
const JUNK_WINDOW = 12, JUNK_MATCH = 11, EDGE_MARGIN = 5;

function junkAt(cons, ref, i) {
  const q = String(cons.base[i] || '-').toUpperCase();
  return !(q === ref[i] && !cons.insert.has(i));
}

export function edgeJunk(gap, cons, ref, maxWalk) {
  const n = ref.length;
  const walk = (from, step) => {
    // from: first covered position (0-based) next to the gap
    for (let k = 0; k <= maxWalk; k++) {
      const s = from + step * k;
      let good = 0, seen = 0;
      for (let w = 0; w < JUNK_WINDOW; w++) {
        const i = s + step * w;
        if (i < 0 || i >= n) break;
        seen++;
        if (!junkAt(cons, ref, i)) good++;
      }
      if (seen < JUNK_WINDOW) return k;
      if (good >= JUNK_MATCH) {
        // Trim back to the first matching base of the window.
        let t = k;
        while (t < k + JUNK_WINDOW && junkAt(cons, ref, from + step * t)) t++;
        return t;
      }
    }
    return maxWalk;
  };
  const left = gap.start > 1 ? walk(gap.start - 2, -1) : 0; // the covered base before the gap
  const right = gap.end < n ? walk(gap.end, 1) : 0;          // the covered base after it
  return { left, right };
}

// Differences between a consensus (anchor coordinates) and a template
// mapped onto the anchor, over positions with a called base, leaving out
// the junk at the edges of missing stretches (skip).
function differencesTo(cons, depth, seq, map, skip = null) {
  let diffs = 0, compared = 0;
  for (let i = 0; i < cons.n; i++) {
    if (!depth[i] || skip?.has(i)) continue;
    const j = map ? map[i] : i;
    if (j < 0) { continue; }
    const q = String(cons.base[i] || '-').toUpperCase();
    if (!BASES.includes(q)) continue;
    compared++;
    if (q !== seq[j]) diffs++;
  }
  return { diffs, compared };
}

// ── Translation ──

const CODE = (() => {
  const aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
  const map = new Map();
  let i = 0;
  for (const x of 'TCAG') for (const y of 'TCAG') for (const z of 'TCAG') map.set(x + y + z, aa[i++]);
  return map;
})();

// A template whose frame-1 translation is a clean coding sequence: whole
// codons, a stop at the end and none before it (ResFinder's van operons and
// a few others are not, and are skipped).
export function codingTemplate(seq) {
  if (seq.length < 9 || seq.length % 3) return false;
  for (let i = 0; i < seq.length - 3; i += 3) if (STOPS.has(seq.slice(i, i + 3))) return false;
  return STOPS.has(seq.slice(-3));
}

// The homopolymer run of the template around position i, [from, to]: an
// indel in a run can be aligned to any of its positions, and the reads that
// carry it are spread over them.
function runAround(seq, i) {
  let from = i, to = i;
  while (from > 0 && seq[from - 1] === seq[i]) from--;
  while (to + 1 < seq.length && seq[to + 1] === seq[i]) to++;
  return [from, to];
}

// Share of the reads at position i carrying a deletion (summed over the
// homopolymer run), and of the reads passing position i that carry an
// insertion after it (insertion rows of the run, ins: the matrix's).
function deletionSupport(seq, counts, depth, i) {
  const [from, to] = runAround(seq, i);
  let del = 0;
  for (let k = from; k <= to; k++) del += counts[k * 6 + 5];
  return depth[i] ? Math.min(1, del / depth[i]) : 0;
}

function insertionSupport(seq, ins, i) {
  if (!ins) return 0;
  const [from, to] = i >= 0 ? runAround(seq, i) : [-1, -1];
  let best = 0;
  for (let k = from - 1; k <= to; k++) {
    for (const c of ins.get(k) || []) {
      const with_ = c[0] + c[1] + c[2] + c[3] + c[4];
      const total = with_ + c[5];
      if (total) best = Math.max(best, with_ / total);
    }
  }
  return best;
}

// Translates the consensus of a gene in the template's frame: uncovered
// positions become N (the frame is kept), deletions in the sample drop the
// base, inserted bases are added. Returns the first stop codon before the
// template's own and what caused it, or that the protein is complete.
// support (per stop) is the smallest share of reads behind the bases or the
// indel that make it; confident needs EVIDENCE.stopShare of the reads at
// EVIDENCE.stopDepth, and is never set for Nanopore reads, whose consensus
// can carry homopolymer errors.
export function translationCheck(seq, cons, counts, depth, { nanopore = false, ins = null } = {}) {
  if (!codingTemplate(seq)) return { status: 'na', reason: 'not a coding template' };
  const codons = seq.length / 3 - 1; // protein length without the stop
  const bases = []; // { b, tpos (0-based template position or -1 for inserted), share }
  const events = []; // indels: { tpos, kind, length, share }
  const share = (i, k) => (depth[i] ? counts[i * 6 + k] / depth[i] : 0);
  const lead = cons.insert.get(-1);
  if (lead) for (const ch of lead) bases.push({ b: ch.toUpperCase(), tpos: -1, share: insertionSupport(seq, ins, -1) });
  for (let i = 0; i < cons.n; i++) {
    const Q = String(cons.base[i] || '-').toUpperCase();
    if (BASES.includes(Q)) {
      bases.push({ b: Q, tpos: i, share: depth[i] >= EVIDENCE.stopDepth ? share(i, COL[Q]) : 0 });
    } else if (depth[i] && deletionSupport(seq, counts, depth, i) >= 0.5) {
      events.push({ tpos: i, kind: 'del', length: 1, share: depth[i] >= EVIDENCE.stopDepth ? deletionSupport(seq, counts, depth, i) : 0 });
    } else {
      bases.push({ b: 'N', tpos: i, share: 0 });
    }
    const inserted = cons.insert.get(i);
    if (inserted) {
      const s = depth[i] >= EVIDENCE.stopDepth ? insertionSupport(seq, ins, i) : 0;
      events.push({ tpos: i, kind: 'ins', length: inserted.length, share: s });
      for (const ch of inserted) bases.push({ b: ch.toUpperCase(), tpos: -1, share: s });
    }
  }
  // Merge consecutive deletions into one event.
  const merged = [];
  for (const e of events) {
    const last = merged[merged.length - 1];
    if (last && last.kind === 'del' && e.kind === 'del' && e.tpos === last.tpos + last.length) {
      last.length++;
      last.share = Math.min(last.share, e.share);
    } else merged.push({ ...e });
  }
  const triplet = (k) => bases[k].b + bases[k + 1].b + bases[k + 2].b;
  // A stop or frameshift near the start can move the start codon instead of
  // ending the protein, as in alleles annotated from an upstream start
  // (aadA2 is aadA2b with 39 more bases): a start codon after it from which
  // the sample's sequence reads, without a stop, to the template's own stop
  // codon, giving at least EVIDENCE.truncatedProtein of the protein.
  const laterStart = (from) => {
    for (let p = Math.max(0, from); p + 2 < bases.length; p++) {
      if (!STARTS.has(triplet(p))) continue;
      let k = p;
      while (k + 2 < bases.length && !STOPS.has(triplet(k))) k += 3;
      if (k + 2 >= bases.length) return null;
      const length = (k - p) / 3;
      if (length < EVIDENCE.truncatedProtein * codons) return null; // later starts only give shorter proteins
      if (bases[k].tpos === seq.length - 3) return { codon: Math.floor(Math.max(0, bases[p].tpos) / 3) + 1, templatePosition: bases[p].tpos + 1, length };
    }
    return null;
  };
  for (let c = 0; c + 2 < bases.length; c += 3) {
    const codon = triplet(c);
    if (!STOPS.has(codon)) continue;
    const index = c / 3; // 0-based codon in the sample's protein
    if (index >= codons) return { status: 'complete', codons };
    const firstT = bases[c].tpos >= 0 ? bases[c].tpos : bases[c + 1].tpos >= 0 ? bases[c + 1].tpos : bases[c + 2].tpos;
    const before = merged.filter(e => e.tpos < (firstT >= 0 ? firstT : Infinity));
    const shift = before.reduce((s, e) => s + (e.kind === 'ins' ? e.length : -e.length), 0);
    const cause = shift % 3 !== 0 ? before.find(e => e.length % 3 !== 0) : null;
    const support = cause ? cause.share : Math.min(bases[c].share, bases[c + 1].share, bases[c + 2].share);
    const from = cause ? bases.findIndex(x => x.tpos > cause.tpos) : c + 3;
    const start = from >= 0 ? laterStart(from) : null;
    if (start) return { status: 'complete', codons, start, startCause: cause ? 'frameshift' : 'stop', startSupport: support };
    return {
      status: 'premature_stop',
      codons,
      stopCodon: index + 1,
      ratio: index / codons,
      templatePosition: firstT >= 0 ? firstT + 1 : null,
      cause: cause ? 'frameshift' : 'substitution',
      indel: cause ? { position: cause.tpos + 1, kind: cause.kind, length: cause.length } : null,
      support,
      confident: !nanopore && support >= EVIDENCE.stopShare,
    };
  }
  // No stop at all in the consensus: complete unless the reads end early.
  return { status: bases.some(x => x.b === 'N') ? 'incomplete' : 'complete', codons };
}

// ── Evidence for a group of hits ──

// Describes the reads of one gene: members are { template, row } with the
// anchor first; pooled when there is more than one. Returns null when the
// matrix is missing.
export function geneEvidence(result, members, maps = null) {
  const src = sourceOf(result);
  const ms = members.map(x => ({ ...x, m: matrixOf(result, x.template) }));
  if (!ms.length || ms.some(x => !x.m)) return null;
  const anchor = ms[0];
  const n = anchor.m.n;
  const pooled = ms.length > 1;
  const counts = pooled ? pooledCounts(ms, maps) : anchor.m.counts;
  const view = { n, counts, ins: anchor.m.ins, ref: anchor.m.ref };
  const depth = depthOf(view);
  const gaps = gapsOf(depth);
  let covered = 0;
  for (let i = 0; i < n; i++) if (depth[i]) covered++;
  const med = median(depth);
  const readLength = readLengthOf(src, ms.map(x => x.template), n);
  const alnOf = (t) => src.aln.get(String(t).trim()) || src.aln.get(t) || '';
  const cons = pooled ? pooledConsensus(ms, maps, alnOf) : consensusColumns(alnOf(anchor.template));

  // Stretches with no reads, and whether thin coverage explains them. A
  // missing stretch extends over the junk the reads crossing its edges leave
  // (extent: [start, end] with those positions, 1-based).
  const junk = new Set();
  const scored = gaps.map(g => {
    const log10P = gapLog10P(g, n, med, readLength);
    const missing = log10P < Math.log10(EVIDENCE.missingP)
      && g.length >= Math.min(EVIDENCE.minMissing, Math.max(1, Math.ceil(EVIDENCE.minMissingShare * n)));
    let extent = [g.start, g.end];
    if (missing && cons) {
      const e = edgeJunk(g, cons, anchor.m.ref, Math.round(readLength));
      extent = [g.start - e.left, g.end + e.right];
      // The edge is known to a few bases (a junk base matches the template
      // one time in four), so allele comparisons leave EDGE_MARGIN more out.
      for (let i = Math.max(0, extent[0] - 1 - EDGE_MARGIN); i < g.start - 1; i++) junk.add(i);
      for (let i = g.end; i < Math.min(n, extent[1] + EDGE_MARGIN); i++) junk.add(i);
    }
    return { ...g, log10P, missing, extent };
  });

  // Covered stretches that align poorly (the long-read trace of a missing
  // part, or a divergent region), from the .aln with the uncovered
  // positions left to the gap test.
  let poor = null;
  if (!pooled) {
    const spanned = new Uint8Array(n);
    for (let i = 0; i < n; i++) spanned[i] = depth[i] ? 1 : 0;
    poor = alignmentCheck(alnOf(anchor.template), { spanned });
  }

  // The consensus against the anchor's template (pooled consensus differs
  // from each member at the positions it took from the others), leaving out
  // the junk at missing stretches and, on long reads, the stretches that
  // align poorly (flanking sequence aligned across a missing part).
  if (poor?.flagged) for (const s of poor.segments) for (let i = s.start - 1; i < s.end; i++) junk.add(i);
  const vs = cons ? differencesTo(cons, depth, anchor.m.ref, null, junk) : { diffs: 0, compared: 0 };
  // The translation of a gene with a missing or poorly aligned part says
  // nothing new (its protein is already known to be incomplete), and the
  // junk at a truncation's edge would read as a frameshift.
  const missing = scored.filter(g => g.missing);
  const partial = missing.length || poor?.flagged || covered < 0.9 * n;
  const translation = !cons ? { status: 'na', reason: 'no alignment' }
    : partial ? { status: 'na', reason: 'partial' }
      : translationCheck(anchor.m.ref, cons, counts, depth, { nanopore: src.nanopore, ins: anchor.m.ins });

  // Mixed bases, away from gap edges, poorly aligned windows and the gene's
  // own ends (reads hanging over an end align their overhang base by base).
  const skip = new Set();
  for (let i = 0; i < Math.min(10, n); i++) { skip.add(i); skip.add(n - 1 - i); }
  const margin = Math.round(readLength / 2);
  for (const g of scored) for (let i = g.start - 1 - margin; i < g.end + margin; i++) if (i >= 0 && i < n) skip.add(i);
  if (poor) {
    for (const s of poor.segments) for (let i = s.start - 1 - 10; i < s.end + 10; i++) if (i >= 0 && i < n) skip.add(i);
  }
  const mixed = mixedSites(view, src.nanopore, skip);

  return {
    n,
    covered,
    breadth: n ? covered / n : 0,
    meanDepth: depth.reduce((s, d) => s + d, 0) / Math.max(1, n),
    medianDepth: med,
    readLength,
    gaps: scored,
    junk,
    missing,
    thin: scored.filter(g => !g.missing),
    poor,
    differences: vs.diffs,
    compared: vs.compared,
    identity: vs.compared ? 1 - vs.diffs / vs.compared : null,
    translation,
    mixed,
    consensus: cons,
    depth,
    nanopore: src.nanopore,
  };
}

// ── Sister templates ──

// Whether the reads of `other` are reads of the same allele as the cluster
// (anchor coordinates) or of a second allele. At each position where the two
// templates differ and the pooled reads are deep enough, the site is mixed
// when both bases are there beyond sequencing error. One allele: no mixed
// site (or one, among five or more). Two alleles: at least two mixed sites,
// half of those that can be read or more; or, wherever they are, at least
// two positions where the reads of `other` and the cluster's own each agree
// on a base (EVIDENCE.agree of them) and the bases differ: when KMA spreads
// two alleles over many templates, they need not differ where the
// templates do (a chromosomal and a plasmid SHV).
export function sisterTest(result, clusterMembers, clusterMaps, other) {
  const src = sourceOf(result);
  const anchorM = matrixOf(result, clusterMembers[0].template);
  const otherM = matrixOf(result, other.template);
  if (!anchorM || !otherM) return { relation: 'unknown' };
  const aln = alignTemplates(anchorM.ref, otherM.ref);
  if (!aln || aln.identity < EVIDENCE.sisterIdentity || aln.aligned < 0.5 * Math.min(anchorM.n, otherM.n)) {
    return { relation: 'unrelated', identity: aln?.identity ?? 0 };
  }
  const ms = clusterMembers.map(x => ({ ...x, m: matrixOf(result, x.template) }));
  const counts = ms.length > 1 ? pooledCounts(ms, clusterMaps) : anchorM.counts;
  const disc = [];
  for (let i = 0; i < anchorM.n; i++) {
    const j = aln.map[i];
    if (j < 0) continue;
    const a = anchorM.ref[i], b = otherM.ref[j];
    if (a !== b && BASES.includes(a) && BASES.includes(b)) disc.push([i, j]);
  }
  const skip = new Set(disc.map(([i]) => i));
  const err = errorRate(counts, anchorM.n, skip, src.nanopore);
  let evaluable = 0, mixedN = 0, otherBase = 0;
  const sites = [];
  for (const [i, j] of disc) {
    const c = [0, 1, 2, 3].map(k => counts[i * 6 + k] + otherM.counts[j * 6 + k]);
    const depth = c[0] + c[1] + c[2] + c[3];
    if (depth < EVIDENCE.minSiteDepth) continue;
    evaluable++;
    const s = siteMix(c, err, src.nanopore);
    const ai = COL[anchorM.ref[i]], bi = COL[otherM.ref[j]];
    const mixed = s.mixed && ((s.top === anchorM.ref[i] && s.second === otherM.ref[j]) || (s.top === otherM.ref[j] && s.second === anchorM.ref[i]));
    if (mixed) mixedN++;
    if (c[bi] > c[ai]) otherBase++;
    sites.push({ pos: i + 1, a: anchorM.ref[i], b: otherM.ref[j], na: c[ai], nb: c[bi], depth, mixed });
  }
  // Positions where the two sets of reads each agree on a different base.
  let disagree = 0;
  for (let i = 0; i < anchorM.n; i++) {
    const j = aln.map[i];
    if (j < 0) continue;
    const a = majority(counts, i), b = majority(otherM.counts, j);
    if (a && b && a.base !== b.base) disagree++;
  }
  let relation;
  if (disagree >= 2 || (mixedN >= 2 && mixedN >= 0.5 * evaluable)) relation = 'two';
  else if (mixedN === 0 || (mixedN === 1 && evaluable >= 5)) relation = 'one';
  else relation = 'unclear';
  const shares = sites.filter(s => s.mixed).map(s => Math.min(s.na, s.nb) / Math.max(1, s.na + s.nb));
  return {
    relation, identity: aln.identity, map: aln.map, sites, evaluable, mixed: mixedN, disagree,
    otherBase, discriminating: disc.length, minorShare: shares.length ? medianOf(shares) : null,
  };
}

// The base most reads at position i carry, when the reads are deep enough
// and EVIDENCE.agree of them carry it, or null.
function majority(counts, i) {
  let depth = 0, best = 0;
  for (let k = 0; k < 4; k++) {
    depth += counts[i * 6 + k];
    if (counts[i * 6 + k] > counts[i * 6 + best]) best = k;
  }
  if (depth < EVIDENCE.minSiteDepth || counts[i * 6 + best] < EVIDENCE.agree * depth) return null;
  return { base: BASES[best], depth };
}

// One gene from sister templates whose reads are one allele (sisterTest
// "one"): the pooled evidence in the coordinates of the template closest to
// the pooled consensus, which then becomes the anchor (first member).
// Returns { members, maps, evidence, closest, anchorCovered } or null.
export function poolSisters(result, members) {
  const align = (anchor) => members.map(x => {
    if (x === anchor) return null;
    const a = matrixOf(result, anchor.template), b = matrixOf(result, x.template);
    return a && b ? alignTemplates(a.ref, b.ref)?.map || null : null;
  });
  const build = (anchor) => {
    const ordered = [anchor, ...members.filter(x => x !== anchor)];
    const maps = align(anchor);
    const orderedMaps = [null, ...members.filter(x => x !== anchor).map(x => maps[members.indexOf(x)])];
    if (orderedMaps.slice(1).some(m => !m)) return null;
    return { members: ordered, maps: orderedMaps, evidence: geneEvidence(result, ordered, orderedMaps) };
  };
  let pool = build(members[0]);
  if (!pool?.evidence) return null;
  // The closest template: the one the pooled consensus aligns to best,
  // matches minus mismatches over the positions it shares with the anchor
  // (a shorter template that is a part of another matches without a
  // difference, but leaves out sequence the reads have). Ties keep KMA's
  // best-scoring template.
  // A template the reads match exactly over its whole length names the
  // allele, as in assembly-based tools (AMRFinderPlus's EXACT), when it is
  // nearly as long as the gene, or is a whole coding sequence itself (aadA1
  // inside ant(3'')-Ia, which ResFinder annotates from an upstream start):
  // the longest such one.
  const ev = pool.evidence;
  const scores = pool.members.map((x, k) => {
    const m = matrixOf(result, x.template);
    const d = differencesTo(ev.consensus, ev.depth, m.ref, pool.maps[k], ev.junk);
    return { template: x.template, length: m.n, coding: codingTemplate(m.ref), diffs: d.diffs, compared: d.compared, score: d.compared - 2 * d.diffs };
  });
  const exact = scores.filter(s => s.diffs === 0 && s.compared >= s.length && (s.length >= 0.9 * ev.n || s.coding))
    .sort((a, b) => b.length - a.length);
  let best = exact[0] || scores.reduce((a, b) => (b.score > a.score ? b : a), scores[0]);
  if (best.template !== pool.members[0].template) {
    const re = build(pool.members.find(x => x.template === best.template));
    if (re?.evidence) pool = re;
  }
  // A closest template whose frame the sample breaks (a frameshift or stop
  // in ant(3'')-Ia's upstream start) gives way to a sister nearly as close
  // whose frame gives the whole protein (aadA1, two differences): the
  // database annotates the same gene from two starts.
  const broken = (e) => e?.translation?.status === 'premature_stop' && e.translation.confident;
  if (broken(pool.evidence)) {
    const alts = scores.filter(s => s.template !== best.template && s.coding && s.compared && s.diffs / s.compared <= 0.01)
      .sort((a, b) => b.score - a.score);
    for (const a of alts) {
      const re = build(pool.members.find(x => x.template === a.template));
      if (re?.evidence && re.evidence.translation?.status === 'complete') { pool = re; best = a; break; }
    }
  }
  const own = matrixOf(result, pool.members[0].template);
  let anchorCovered = 0;
  for (const d of depthOf(own)) if (d) anchorCovered++;
  const ties = scores.filter(s => s.score === best.score).map(s => s.template);
  return {
    ...pool,
    anchorCovered,
    closest: { template: pool.members[0].template, differences: pool.evidence.differences, compared: pool.evidence.compared, ties },
    scores,
  };
}

export { fileText as evidenceFileText };
