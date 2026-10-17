/* hits.js — the hits of one KMA result sorted into called, partial and
 * hidden, so the report, the gene tables and the hAMRonization export agree.
 *
 *   called   passes both sliders (thresholds.js)
 *   partial  passes the identity slider but not the coverage slider, or
 *            fails them as a whole while its well-aligned part alone would
 *            pass (aln-check.js): a truncated gene on long reads, where the
 *            flanking sequence aligned across the missing part drags the
 *            whole-gene identity under the bar
 *   hidden   below the sliders, or a partial hit of a gene that is already
 *            called in the same database ("shadowed"): with a redundant
 *            database, reads of a called allele also land on sister alleles
 *            at a fraction of their coverage and depth, and those are not
 *            second copies of the gene
 * Every row stays in the .res file in the ZIP.
 */

import { parseResFile } from './resfinder-db.js';
import { geneName, sameVariantFamily } from './genes.js';
import { classifyRow, normaliseThresholds } from './thresholds.js';
import { alnBlocksOf, alignmentCheck, alignmentWarning, partialByAlignment } from './aln-check.js';

const blaFamily = (name) => {
  const m = String(name).toLowerCase().match(/^bla([a-z]+(?:-[a-z]+)?)-?\d/);
  return m ? m[1] : null;
};

// The same gene, or alleles of one family (sul1 and sul1; blaCTX-M-15 and
// blaCTX-M-33; aph(3')-III and aph(3')-IIIa). sul1 and sul2 stay apart.
function sameGene(a, b) {
  if (a.toLowerCase() === b.toLowerCase() || sameVariantFamily(a, b)) return true;
  const fa = blaFamily(a);
  return !!fa && fa === blaFamily(b);
}

// { called: [hit], partial: [hit], shadowed, below } for one runAnalysis
// result, cached on it per thresholds. A hit is
//   { row, name, status, by: 'coverage' | 'alignment', check, warning }.
export function classifyResult(result, thresholds) {
  const th = normaliseThresholds(thresholds);
  const key = JSON.stringify(th);
  if (result._hits && result._hits.key === key) return result._hits;
  const out = { key, called: [], partial: [], shadowed: 0, below: 0 };
  if (result.exitCode === 0 && result.resTable) {
    const blocks = alnBlocksOf(result);
    const partials = [];
    for (const row of parseResFile(result.resTable).rows) {
      let status = classifyRow(row, th);
      const check = alignmentCheck(blocks.get(row.Template));
      let by = 'coverage';
      if (status === 'below' && partialByAlignment(check, th)) { status = 'partial'; by = 'alignment'; }
      if (status === 'below') { out.below++; continue; }
      const hit = { row, name: geneName(row.Template, result.database), status, by, check, warning: alignmentWarning(check) };
      if (status === 'called') out.called.push(hit); else partials.push(hit);
    }
    // One partial hit per gene: the one with the most of the gene well
    // aligned, unless the gene is already called.
    const extent = (h) => h.check ? h.check.goodFraction : parseFloat(h.row.Template_Coverage) / 100 || 0;
    partials.sort((a, b) => extent(b) - extent(a));
    for (const p of partials) {
      if (out.called.some(c => sameGene(c.name, p.name)) || out.partial.some(k => sameGene(k.name, p.name))) out.shadowed++;
      else out.partial.push(p);
    }
  }
  result._hits = out;
  return out;
}
