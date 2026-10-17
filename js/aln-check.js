/* aln-check.js — a check of each called gene's alignment, from KMA's .aln
 * output (the template against the read consensus).
 *
 * KMA's identity and coverage figures can hide a gene that is only partly
 * there. With long reads, the sequence flanking a truncated gene is
 * gap-aligned across the missing part at 55–75% identity, and KMA's base
 * caller breaks ties toward the template base, so Template_Coverage reaches
 * 100% and Template_Identity stays above 90%. With short reads the missing
 * part is simply uncovered. Either way the alignment shows a stretch of
 * mismatches, gaps or uncovered positions, which is what this scans for:
 * windows of the template whose identity falls below a floor, counted as a
 * fraction of the gene. The part outside those windows is the "good part":
 * when it alone would pass the sliders, a hit that fails them as a whole is
 * a partial gene (partialByAlignment).
 */

import { normaliseThresholds, PARTIAL_MIN_COVERAGE } from './thresholds.js';

export const ALN_CHECK = {
  window: 60,        // template positions per window
  step: 20,
  minIdentity: 0.85, // a window below this is poor
  flagFraction: 0.1, // flag a gene when this much of it lies in poor windows
};

// Alignment blocks of a .aln file by template name: "# name" headers, then
// template / match / query lines in groups.
export function indexAln(text) {
  const map = new Map();
  let name = null, buf = [];
  for (const line of String(text || '').split('\n')) {
    if (line[0] === '#') {
      if (name) map.set(name, buf.join('\n').replace(/\n+$/, ''));
      name = line.replace(/^#\s*/, '').trim();
      buf = [];
    } else if (name !== null) {
      buf.push(line);
    }
  }
  if (name) map.set(name, buf.join('\n').replace(/\n+$/, ''));
  return map;
}

// Per result (runAnalysis output), the .aln index, parsed once.
export function alnBlocksOf(result) {
  if (result._alnIndex) return result._alnIndex;
  let map = new Map();
  try {
    const info = result.files?.['.aln'];
    if (info) {
      const text = typeof info.data === 'string' ? info.data : new TextDecoder().decode(info.data);
      map = indexAln(text);
    }
  } catch (_) { /* no alignment: the check is skipped */ }
  result._alnIndex = map;
  return map;
}

// The template and consensus rows of one block as two strings of the same
// length, gaps included.
export function parseAlnBlock(block) {
  let t = '', q = '';
  for (const line of String(block || '').split('\n')) {
    if (line.startsWith('template:')) t += line.slice(line.indexOf('\t') + 1).trim();
    else if (line.startsWith('query:')) q += line.slice(line.indexOf('\t') + 1).trim();
  }
  return { t, q };
}

// Returns null when there is no usable alignment, else
//   { length, poorBases, poorFraction, uncovered, uncoveredFraction,
//     uncoveredInPoor, segments: [{start, end}], flagged,
//     goodFraction, goodIdentity }
// with positions 1-based in template coordinates; goodFraction is the share
// of the template outside poor windows and goodIdentity the identity there.
export function alignmentCheck(block, opts = {}) {
  const o = { ...ALN_CHECK, ...opts };
  const { t, q } = parseAlnBlock(block);
  if (!t.length || t.length !== q.length) return null;

  // Per template position: identical to the consensus, covered at all, and
  // the number of consensus insertions just before it (gaps in the template).
  const ok = [], covered = [], ins = [];
  let pending = 0;
  for (let i = 0; i < t.length; i++) {
    const tc = t[i];
    if (tc === '-') { pending++; continue; }
    const qc = q[i];
    covered.push(qc !== '-' ? 1 : 0);
    ok.push(qc !== '-' && qc.toUpperCase() === tc.toUpperCase() ? 1 : 0);
    ins.push(pending);
    pending = 0;
  }
  const n = ok.length;
  if (!n) return null;

  const w = Math.min(o.window, n);
  const poor = new Uint8Array(n);
  const starts = [];
  for (let s = 0; s + w <= n; s += o.step) starts.push(s);
  if (starts[starts.length - 1] !== n - w) starts.push(n - w);
  for (const s of starts) {
    let m = 0, extra = 0;
    for (let i = s; i < s + w; i++) { m += ok[i]; extra += ins[i]; }
    if (m / (w + extra) < o.minIdentity) poor.fill(1, s, s + w);
  }

  let poorBases = 0, uncovered = 0, uncoveredInPoor = 0;
  let goodBases = 0, goodMatches = 0, goodIns = 0;
  const segments = [];
  for (let i = 0; i < n; i++) {
    if (!covered[i]) uncovered++;
    if (!poor[i]) { goodBases++; goodMatches += ok[i]; goodIns += ins[i]; continue; }
    poorBases++;
    if (!covered[i]) uncoveredInPoor++;
    const last = segments[segments.length - 1];
    if (last && last.end === i) last.end = i + 1;
    else segments.push({ start: i + 1, end: i + 1 });
  }
  return {
    length: n,
    poorBases,
    poorFraction: poorBases / n,
    uncovered,
    uncoveredFraction: uncovered / n,
    uncoveredInPoor,
    segments,
    flagged: poorBases / n >= o.flagFraction,
    goodFraction: goodBases / n,
    goodIdentity: goodBases ? goodMatches / (goodBases + goodIns) : 0,
  };
}

// A hit that fails the sliders as a whole but whose well-aligned part alone
// would pass them: at least the partial floor of the template outside the
// poor windows, at the identity the slider asks. A gene that is uniformly
// divergent has no poor windows and is not partial.
export function partialByAlignment(check, th) {
  if (!check || !check.flagged) return false;
  const t = normaliseThresholds(th);
  const round4 = (x) => Math.round(x * 10000) / 10000;
  return round4(check.goodFraction) >= round4(Math.min(t.coverage, PARTIAL_MIN_COVERAGE))
    && round4(check.goodIdentity) >= round4(t.identity);
}

// One sentence for a flagged check, or null.
//   "29% of the gene (positions 589–828) aligns at under 85% identity."
//   "26% of the gene (positions 614–828) is not covered by reads."
export function alignmentWarning(check, opts = {}) {
  if (!check || !check.flagged) return null;
  const o = { ...ALN_CHECK, ...opts };
  const share = Math.round(100 * check.poorFraction);
  const where = check.segments.slice(0, 3).map(s => `${s.start}–${s.end}`).join(', ')
    + (check.segments.length > 3 ? ', …' : '');
  const mostlyUncovered = check.uncoveredInPoor >= 0.5 * check.poorBases;
  return mostlyUncovered
    ? `${share}% of the gene (positions ${where}) is not covered by reads.`
    : `${share}% of the gene (positions ${where}) aligns at under ${Math.round(100 * o.minIdentity)}% identity.`;
}
