/* thresholds.js — which KMA hits count, from the identity and coverage
 * sliders, applied in the browser after KMA has run.
 *
 * KMA's only gene-level filter is -ID, on Template_Identity: identical
 * bases divided by the whole template length. That bounds coverage too: a
 * gene 90% identical over its whole length is at least 90% covered, so a
 * coverage threshold below the identity threshold can never bind. KMA has
 * no gene coverage filter at all; -mrc, which earlier builds passed the
 * coverage slider to, is "Minimum query coverage", a per-read check.
 *
 * So KMA runs with a low -ID floor and the two sliders apply here, to two
 * independent numbers of each .res row:
 *   identity  Query_Identity: identical bases over the aligned part
 *   coverage  Template_Coverage: the aligned part over the template length
 * Template_Identity equals identity × coverage, which is how the -ID floor
 * follows from the sliders.
 *
 * This differs from ResFinder's read mode on purpose. ResFinder keeps hits
 * with Template_Identity ≥ 90% and Template_Coverage ≥ 60%, so it never
 * reports a gene that is less than 90% covered. Here a hit that passes the
 * identity slider but not the coverage slider is reported as partial, down
 * to PARTIAL_MIN_COVERAGE of the gene, so a truncated or thinly covered
 * gene is visible instead of silently absent.
 */

export const DEFAULT_THRESHOLDS = { identity: 0.9, coverage: 0.6 };

// Hits below the coverage slider that still cover at least this fraction
// of the gene are listed as partial.
export const PARTIAL_MIN_COVERAGE = 0.3;

// Percent with two decimals, free of binary rounding (100 * 0.9 is
// 90.00000000000001 in JavaScript).
const pct = (frac) => Math.round(frac * 10000) / 100;

// Accepts {identity, coverage} (fractions of 1) and the former
// {id_threshold, mrc} shape.
export function normaliseThresholds(th) {
  const pick = (a, b, fallback) => {
    const v = th?.[a] ?? th?.[b];
    return Number.isFinite(v) ? Math.min(1, Math.max(0, v)) : fallback;
  };
  return {
    identity: pick('identity', 'id_threshold', DEFAULT_THRESHOLDS.identity),
    coverage: pick('coverage', 'mrc', DEFAULT_THRESHOLDS.coverage),
  };
}

// The -ID (percent) KMA runs with: just under the lowest Template_Identity a
// listed hit can have. -ID 0 means something else to KMA (it prints every
// template with zero identity), so the floor is at least 1.
export function kmaIdentityFloor(th) {
  const t = normaliseThresholds(th);
  const lowestCoverage = Math.min(t.coverage, PARTIAL_MIN_COVERAGE);
  return Math.max(1, Math.floor(pct(lowestCoverage * t.identity)) - 1);
}

export const num = (v) => { const x = parseFloat(v); return Number.isFinite(x) ? x : null; };

// 'called' passes both sliders; 'partial' passes identity but covers less of
// the gene than the coverage slider asks, and at least PARTIAL_MIN_COVERAGE;
// 'below' is everything else, which the gene tables hide.
export function hitStatus(identity, coverage, th) {
  const t = normaliseThresholds(th);
  if (identity == null || coverage == null) return 'below';
  if (identity < pct(t.identity)) return 'below';
  if (coverage >= pct(t.coverage)) return 'called';
  return coverage >= pct(Math.min(t.coverage, PARTIAL_MIN_COVERAGE)) ? 'partial' : 'below';
}

export function classifyRow(row, th) {
  return hitStatus(num(row.Query_Identity), num(row.Template_Coverage), th);
}

export function splitRows(rows, th) {
  const out = { called: [], partial: [], below: [] };
  for (const row of rows) out[classifyRow(row, th)].push(row);
  return out;
}

// "identity ≥ 90% over the aligned part, coverage ≥ 60% of the gene"
export function describeThresholds(th) {
  const t = normaliseThresholds(th);
  return `identity ≥ ${Math.round(pct(t.identity))}% over the aligned part, coverage ≥ ${Math.round(pct(t.coverage))}% of the gene`;
}
