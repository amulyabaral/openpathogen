/* hits.js — the genes of one KMA result: called, partial and hidden, so the
 * report, the gene tables and the hAMRonization export agree.
 *
 *   called   passes both sliders (thresholds.js)
 *   partial  passes the identity slider but covers less of the gene than
 *            the coverage slider asks; or fails as a whole while its
 *            well-aligned part alone would pass (aln-check.js: a truncated
 *            gene on long reads, whose flanks are aligned across the missing
 *            part); or is complete but its protein is cut short by a stop
 *            codon or a frameshift that the reads clearly carry
 *   hidden   below the sliders, or the reads of a gene that KMA placed on a
 *            sister allele
 *
 * KMA reports templates; a gene can be several. When the sample's allele is
 * not in the database, KMA can divide its reads between two or more sister
 * templates, and with a redundant database a few reads of every gene land on
 * its sister alleles. The templates of one gene are grouped by name, and
 * evidence.js decides from the reads whether two of them are one allele
 * (their reads are pooled into one gene, named by the closest template) or
 * two alleles (both stay). The pooled reads, rather than each template's,
 * then decide whether a stretch of the gene is missing (evidence.js).
 *
 * Each gene's numbers are KMA's .res numbers for its template, unless the
 * pooled reads cover more of the gene than that template's own, when they
 * are the pooled coverage, the identity of the pooled consensus to the
 * closest template and the summed depth. Results without a depth matrix
 * fall back to KMA's numbers alone. Every .res row stays in the ZIP.
 */

import { parseResFile } from './resfinder-db.js';
import { alleleLabel, geneName, sameVariantFamily } from './genes.js';
import { classifyRow, hitStatus, normaliseThresholds, num } from './thresholds.js';
import { alnBlocksOf, alignmentCheck, alignmentWarning, partialByAlignment } from './aln-check.js';
import { EVIDENCE, hasEvidence, geneEvidence, sisterTest, poolSisters, matrixOf, alignTemplates } from './evidence.js';

// A beta-lactamase family: ResFinder's "blaCTX-M-15", CARD's "CTX-M-15".
const blaFamily = (name) => {
  const m = String(name).toLowerCase().match(/^bla([a-z]+(?:-[a-z]+)?)-?\d/);
  if (m) return m[1];
  const card = String(name).match(/^([A-Z]{2,5}(?:-[A-Z]{1,2})?)-\d+[A-Za-z]?$/);
  return card ? card[1].toLowerCase() : null;
};

// Names ResFinder gives one gene twice: ant(3'')-Ia (X02340) is aadA1 with
// 180 bases of upstream sequence.
const ALIASES = { "ant(3'')-ia": 'aadA1' };
const alias = (name) => ALIASES[String(name).toLowerCase()] || name;

// The same gene, or alleles of one family (sul1 and sul1; blaCTX-M-15 and
// blaCTX-M-33; aph(3')-III and aph(3')-IIIa). sul1 and sul2 stay apart.
export function sameGene(a, b) {
  a = alias(a);
  b = alias(b);
  if (a.toLowerCase() === b.toLowerCase() || sameVariantFamily(a, b)) return true;
  const fa = blaFamily(a);
  return !!fa && fa === blaFamily(b);
}

const pct = (frac) => Math.round(frac * 10000) / 100;

// { called: [hit], partial: [hit], merged, below } for one runAnalysis
// result, cached on it per thresholds. A hit is one gene:
//   { row, template, name, status, by, numbers: {identity, coverage, depth,
//     source}, members: [{template, row}], evidence, closest, related,
//     warnings: [text], notes: [text], warning, check }
// row is its template's .res row (the closest allele when reads were
// pooled), by says what made a hit partial: 'coverage', 'alignment',
// 'stop' or 'frameshift'.
export function classifyResult(result, thresholds) {
  const th = normaliseThresholds(thresholds);
  const key = JSON.stringify(th);
  if (result._hits && result._hits.key === key) return result._hits;
  const out = { key, called: [], partial: [], merged: 0, below: 0 };
  if (result.exitCode === 0 && result.resTable) {
    const rows = parseResFile(result.resTable).rows;
    let done = false;
    if (hasEvidence(result)) {
      try {
        classifyWithEvidence(result, rows, th, out);
        done = true;
      } catch (err) {
        // The results must show even if the reads cannot be read this way.
        console.warn('[hits] reading the evidence failed; using KMA\'s numbers:', err);
        Object.assign(out, { called: [], partial: [], merged: 0, below: 0 });
      }
    }
    if (!done) classifyByNumbers(result, rows, th, out);
  }
  out.shadowed = out.merged; // the name earlier callers used
  result._hits = out;
  return out;
}

// ── With the depth matrix ──

function classifyWithEvidence(result, rows, th, out) {
  const blocks = alnBlocksOf(result);
  const idMin = pct(th.identity);
  // Candidates: every template at the identity threshold, whatever its
  // coverage (a sister holding a small part of the reads still completes
  // its gene), and those whose well-aligned part passes on its own.
  const cands = [];
  for (const row of rows) {
    const qi = num(row.Query_Identity);
    const check = alignmentCheck(blocks.get(row.Template));
    const rescue = !(qi != null && qi >= idMin) && partialByAlignment(check, th);
    if ((qi != null && qi >= idMin) || rescue) {
      cands.push({ template: row.Template, row, name: geneName(row.Template, result.database), score: num(row.Score) || 0, check, rescue });
    } else {
      out.below++;
    }
  }
  cands.sort((a, b) => b.score - a.score || a.template.localeCompare(b.template));

  // Templates of one gene, by name; within a gene, sisters whose reads are
  // one allele form one cluster, best-scoring template first.
  const families = [];
  for (const c of cands) {
    let fam = families.find(f => sameGene(f[0].name, c.name));
    if (!fam) families.push(fam = []);
    fam.push(c);
  }
  const clusters = [];
  for (const fam of families) {
    const local = [];
    for (const c of fam) {
      const links = [];
      let home = null;
      for (const cl of local) {
        const t = sisterTest(result, cl.members, cl.maps, c);
        if (t.relation === 'one') { home = cl; cl.members.push(c); cl.maps.push(t.map); cl.tests.push({ template: c.template, test: t }); break; }
        if (t.relation === 'two' || t.relation === 'unclear') links.push({ cl, test: t });
      }
      if (!home) {
        home = { members: [c], maps: [null], tests: [], related: [] };
        local.push(home);
      }
      for (const { cl, test } of links) {
        cl.related.push({ cluster: home, relation: test.relation, test, other: c.template });
        home.related.push({ cluster: cl, relation: test.relation, test, other: cl.members[0].template });
      }
    }
    clusters.push(...local);
  }
  mergeAcrossNames(result, clusters);

  const genes = [];
  for (const cl of clusters) {
    const gene = buildGene(result, cl, th);
    if (!gene) continue;
    if (gene.status === 'below') { out.below += gene.members.length; continue; }
    out.merged += gene.members.length - 1;
    genes.push(gene);
  }
  // Second alleles name each other by the anchor each ended up with.
  for (const g of genes) {
    g.related = g.cluster.related
      .map(r => ({ relation: r.relation, test: r.test, other: genes.find(x => x.cluster === r.cluster) }))
      .filter(r => r.other);
  }
  for (const g of genes) delete g.cluster;
  flagLowDepth(genes);
  for (const g of genes) {
    Object.assign(g, describeGene(g));
    (g.status === 'called' ? out.called : out.partial).push(g);
  }
}

// Templates the database names as different genes can be one locus whose
// reads KMA spread over them (Serratia's chromosomal AmpC as blaSRT-2,
// blaSRT-1 and blaSST-1). A cluster joins a higher-scoring one of another
// name when their templates share most short words (KMER-mers, sampled every
// STEP bases: about 88% identity or more) and the sister test finds one
// allele at 90% identity or more.
const KMER = 12, STEP = 6, SHARED = 0.2, MERGE_IDENTITY = 0.9;
function mergeAcrossNames(result, clusters) {
  if (clusters.length < 2) return;
  const info = clusters.map(cl => {
    const m = matrixOf(result, cl.members[0].template);
    if (!m) return null;
    const all = new Set(), sample = [];
    for (let i = 0; i + KMER <= m.ref.length; i++) {
      const w = m.ref.slice(i, i + KMER);
      all.add(w);
      if (i % STEP === 0) sample.push(w);
    }
    return { m, all, sample, score: cl.members[0].score };
  });
  const order = clusters.map((cl, k) => k).sort((a, b) => (info[a]?.score ?? 0) - (info[b]?.score ?? 0));
  const gone = new Set();
  for (const k of order) {
    const b = clusters[k], ib = info[k];
    if (!ib || !ib.sample.length) continue;
    let target = null;
    for (let a = 0; a < clusters.length && !target; a++) {
      const ia = info[a];
      if (a === k || gone.has(a) || !ia || ia.score < ib.score) continue;
      if (sameGene(clusters[a].members[0].name, b.members[0].name)) continue; // tested already
      const len = ib.m.n / ia.m.n;
      if (len < 0.8 || len > 1.25) continue;
      let shared = 0;
      for (const w of ib.sample) if (ia.all.has(w)) shared++;
      if (shared < SHARED * ib.sample.length) continue;
      const t = sisterTest(result, clusters[a].members, clusters[a].maps, b.members[0]);
      if (t.relation === 'one' && t.identity >= MERGE_IDENTITY) target = { a, t };
    }
    if (!target) continue;
    const A = clusters[target.a];
    for (const x of b.members) {
      const map = x === b.members[0] ? target.t.map : alignTemplates(info[target.a].m.ref, matrixOf(result, x.template)?.ref || '')?.map;
      if (!map) continue;
      A.members.push(x);
      A.maps.push(map);
    }
    A.tests.push({ template: b.members[0].template, test: target.t }, ...b.tests);
    // Links to the merged cluster now point at the one it joined.
    for (const cl of clusters) for (const r of cl.related) if (r.cluster === b) r.cluster = A;
    for (const r of b.related) if (r.cluster !== A && !A.related.some(x => x.cluster === r.cluster)) A.related.push(r);
    A.related = A.related.filter(r => r.cluster !== A);
    gone.add(k);
  }
  for (const k of [...gone].sort((x, y) => y - x)) clusters.splice(k, 1);
}

function buildGene(result, cl, th) {
  let members = cl.members;
  let ev = null, closest = null, numbers = null;
  const anchorRow = () => members[0].row;
  const kmaNumbers = () => ({
    identity: num(anchorRow().Query_Identity), coverage: num(anchorRow().Template_Coverage),
    depth: num(anchorRow().Depth), source: 'kma',
  });
  if (members.length > 1) {
    const pool = poolSisters(result, members);
    if (pool) {
      members = pool.members;
      ev = pool.evidence;
      closest = pool.closest;
      numbers = kmaNumbers();
      if (ev.covered > pool.anchorCovered + 0.01 * ev.n) {
        // Two decimals, as KMA writes its own.
        const r2 = (x) => Math.round(x * 100) / 100;
        numbers = {
          identity: ev.identity != null ? r2(100 * ev.identity) : numbers.identity,
          coverage: r2(100 * ev.breadth),
          depth: r2(members.reduce((s, x) => s + (num(x.row.Depth) || 0), 0)),
          source: 'pooled',
        };
      }
    }
  }
  if (!ev) {
    ev = geneEvidence(result, [members[0]]);
    numbers = kmaNumbers();
    if (ev) closest = { template: members[0].template, differences: ev.differences, compared: ev.compared, ties: [members[0].template] };
  }
  const anchor = members[0];
  let status = hitStatus(numbers.identity, numbers.coverage, th);
  let by = 'coverage';
  if (status === 'below' && anchor.rescue) { status = 'partial'; by = 'alignment'; }
  const tr = ev?.translation;
  if (status === 'called' && tr?.status === 'premature_stop' && tr.confident && tr.ratio < EVIDENCE.truncatedProtein) {
    status = 'partial';
    by = tr.cause === 'frameshift' ? 'frameshift' : 'stop';
  }
  return {
    row: anchor.row, template: anchor.template, name: anchor.name, db: result.database,
    status, by, numbers, evidence: ev, closest,
    members: members.map(x => ({ template: x.template, row: x.row })),
    pooled: numbers.source === 'pooled',
    tests: cl.tests,
    check: ev?.poor || anchor.check,
    cluster: cl,
    lowDepth: null,
  };
}

// Genes at under a fifth of the median depth of the others called in the
// same run: a gene in a minority of the cells, or from contamination or
// index hopping. Needs three called genes for a median.
const VERY_LOW_DEPTH = 0.05;

function flagLowDepth(genes) {
  const depths = genes.filter(g => g.status === 'called' && g.numbers.depth > 0).map(g => g.numbers.depth).sort((a, b) => a - b);
  if (depths.length < 3) return;
  const med = depths[depths.length >> 1];
  for (const g of genes) {
    if (g.numbers.depth != null && g.numbers.depth < 0.2 * med) g.lowDepth = { depth: g.numbers.depth, median: med };
  }
}

// ── Without a depth matrix: KMA's numbers alone ──

function classifyByNumbers(result, rows, th, out) {
  const blocks = alnBlocksOf(result);
  const partials = [];
  const asGene = (row, status, by, check) => {
    const g = {
      row, template: row.Template, name: geneName(row.Template, result.database), db: result.database, status, by,
      numbers: { identity: num(row.Query_Identity), coverage: num(row.Template_Coverage), depth: num(row.Depth), source: 'kma' },
      members: [{ template: row.Template, row }], evidence: null, closest: null, related: [], pooled: false, tests: [],
      check, lowDepth: null,
    };
    return Object.assign(g, describeGene(g));
  };
  for (const row of rows) {
    let status = classifyRow(row, th);
    const check = alignmentCheck(blocks.get(row.Template));
    let by = 'coverage';
    if (status === 'below' && partialByAlignment(check, th)) { status = 'partial'; by = 'alignment'; }
    if (status === 'below') { out.below++; continue; }
    const g = asGene(row, status, by, check);
    if (status === 'called') out.called.push(g); else partials.push(g);
  }
  // One partial hit per gene, the one with the most of the gene well
  // aligned, unless the gene is already called.
  const extent = (h) => h.check ? h.check.goodFraction : (h.numbers.coverage || 0) / 100;
  partials.sort((a, b) => extent(b) - extent(a));
  for (const p of partials) {
    if (out.called.some(c => sameGene(c.name, p.name)) || out.partial.some(k => sameGene(k.name, p.name))) out.merged++;
    else out.partial.push(p);
  }
}

// ── What the evidence says, in words ──

const fmt = (n) => Number(n).toLocaleString('en-US');
const share = (x) => `${Math.round(100 * x)}%`;
const depthText = (d) => `${d < 10 ? d.toFixed(1) : Math.round(d)}×`;
const range = (a, b) => (a === b ? `position ${a}` : `positions ${a}–${b}`);

// Missing stretches, with the junk the reads crossing their edges leave
// (evidence.js edgeJunk): [{ start, end, gap }].
function missingStretches(g) {
  const ev = g.evidence;
  if (!ev) return [];
  return ev.missing.map(gap => ({ start: gap.extent[0], end: gap.extent[1], gap }));
}

// The compact phrase for a partial gene's line: "the first 238 bases are
// missing", "thinly covered at 2×, so it may be complete".
export function partialPhrase(g) {
  const ev = g.evidence;
  if (g.by === 'alignment') {
    const segs = g.check?.segments || [];
    if (!segs.length) return null;
    const where = segs.slice(0, 2).map(x => range(x.start, x.end)).join(' and ') + (segs.length > 2 ? ', …' : '');
    return `${where} align${segs.length === 1 && segs[0].start === segs[0].end ? 's' : ''} at under 85% identity`;
  }
  if (g.by === 'stop' || g.by === 'frameshift') {
    const t = ev.translation;
    return `${g.by === 'frameshift' ? 'a frameshift' : 'a stop codon'} cuts the protein to ${share(t.ratio)} of its length`;
  }
  if (!ev) return null;
  const miss = missingStretches(g);
  if (miss.length) {
    const n = ev.n;
    const m = miss[0];
    const len = m.end - m.start + 1;
    const where = m.start === 1 ? `the first ${fmt(len)} bases are missing`
      : m.end === n ? `the last ${fmt(len)} bases are missing`
        : `${range(m.start, m.end)} are missing`;
    return miss.length > 1 ? `${where}, and ${miss.length - 1} more stretch${miss.length > 2 ? 'es' : ''}` : where;
  }
  if (ev.thin.length && ev.covered < ev.n) {
    return `thinly covered at ${depthText(ev.medianDepth)}, so the gene may be complete`;
  }
  return null;
}

// Mixed sites packed into one stretch: at least JUNCTION.sites of them
// within JUNCTION.span bases, holding most of the gene's mixed sites, or
// null. Two alleles differ along the whole gene instead.
const JUNCTION = { sites: 5, span: 60, share: 0.8 };
function junctionOf(sites) {
  if (sites.length < JUNCTION.sites) return null;
  let best = null;
  for (let i = 0, j = 0; j < sites.length; j++) {
    while (sites[j].pos - sites[i].pos > JUNCTION.span) i++;
    if (!best || j - i + 1 > best.count) best = { count: j - i + 1, from: i, to: j };
  }
  if (best.count < JUNCTION.sites || best.count < JUNCTION.share * sites.length) return null;
  const run = sites.slice(best.from, best.to + 1);
  return { sites: run.length, start: run[0].pos, end: run[run.length - 1].pos,
    lo: Math.min(...run.map(s => s.share)), hi: Math.max(...run.map(s => s.share)) };
}

// Warnings (the "check alignment" mark) and notes for one gene.
export function describeGene(g) {
  const warnings = [];
  const notes = [];
  const ev = g.evidence;
  const label = (t) => alleleLabel(t, g.db);
  if (!ev) {
    const w = g.status === 'called' ? alignmentWarning(g.check) : null;
    if (w) warnings.push(w);
    return { warnings, notes, warning: warnings[0] || null };
  }
  const n = ev.n;
  for (const m of missingStretches(g)) {
    const len = m.end - m.start + 1;
    const g0 = m.gap;
    const edgeNote = m.end > g0.end || m.start < g0.start
      ? ` (no read reaches ${range(g0.start, g0.end)}, and the reads at the edge only partly match)` : '';
    warnings.push(`${cap(range(m.start, m.end))} (${share(len / n)} of the gene) ${len === 1 ? 'is' : 'are'} not in the reads${edgeNote}, while the rest of the gene is covered at about ${depthText(ev.medianDepth)}. A stretch this long is not missed by chance at that depth, so this part of the gene is missing from the sample, or too different to align (its reads may have gone to a related sequence in the database).`);
  }
  if (!ev.missing.length && ev.thin.length && ev.covered < n) {
    const unc = n - ev.covered;
    notes.push(`${share(unc / n)} of the gene has no reads, in ${ev.thin.length} stretch${ev.thin.length > 1 ? 'es' : ''}. At a depth of ${depthText(ev.medianDepth)} gaps like these are expected, so the gene may be complete.`);
  }
  // Covered stretches that align poorly, apart from the edges of missing ones.
  if (ev.poor?.flagged) {
    const miss = missingStretches(g);
    const segs = ev.poor.segments.filter(s => !miss.some(m => s.start >= m.start && s.end <= m.end));
    const bases = segs.reduce((x, s) => x + s.end - s.start + 1, 0);
    if (segs.length && bases / n >= 0.1) {
      const where = segs.slice(0, 3).map(s => `${s.start}–${s.end}`).join(', ') + (segs.length > 3 ? ', …' : '');
      warnings.push(`${share(bases / n)} of the gene (positions ${where}) aligns at under 85% identity: the gene may be truncated or interrupted there, or the reads come from a related sequence.`);
    }
  }
  const t = ev.translation;
  if (t?.status === 'complete' && t.start) {
    notes.push(`${t.startCause === 'frameshift' ? 'A frameshift' : 'A stop codon'} near the start of the gene is followed by a start codon at position ${fmt(t.start.templatePosition)}, from which the protein is complete (${share(t.start.length / t.codons)} of its length): the allele likely starts there, as some alleles in the database do.`);
  }
  if (t?.status === 'premature_stop' && t.ratio < 1) {
    const how = t.cause === 'frameshift' && t.indel
      ? `${t.indel.kind === 'del' ? (t.indel.length === 1 ? 'A deleted base' : `${t.indel.length} deleted bases`) : (t.indel.length === 1 ? 'An inserted base' : `${t.indel.length} inserted bases`)} at position ${t.indel.position} shifts the reading frame, and the protein stops at codon ${fmt(t.stopCodon)} of ${fmt(t.codons)}`
      : `A stop codon at codon ${fmt(t.stopCodon)} of ${fmt(t.codons)}${t.templatePosition ? ` (position ${t.templatePosition})` : ''} cuts the protein short`;
    const len = `, ${share(t.ratio)} of its length`;
    const sup = ` ${share(t.support)} of the reads there carry it.`;
    if (t.confident && t.ratio < EVIDENCE.truncatedProtein) {
      warnings.push(`${how}${len}, so the gene is probably not functional.${sup}`);
    } else if (ev.nanopore) {
      warnings.push(`${how}${len}.${sup} On Nanopore reads this can also be a sequencing error, typically in a run of one base.`);
    } else if (!t.confident) {
      warnings.push(`${how}${len}, but only ${share(t.support)} of the reads there carry it: the gene may be interrupted here, or the reads are mixed.`);
    } else {
      notes.push(`${how}${len}, near its end.${sup}`);
    }
  }
  // Mixed bases within one template: spread along the gene, two alleles on
  // one template; packed into one short stretch, the reads of two different
  // sequences meeting there (a gene in two partial copies whose reads
  // together cover it, or interrupted at that point).
  const mixed = ev.mixed?.sites || [];
  const junction = junctionOf(mixed);
  if (junction) {
    warnings.push(`The reads disagree at ${junction.sites} positions between ${junction.start} and ${junction.end}, where ${share(junction.lo)}–${share(junction.hi)} of them carry other bases, while the rest of the gene reads as one sequence. Reads of two different sequences meet there: the gene may be in two partial copies that together cover it, or be interrupted at that point.`);
  } else if (!g.pooled && (mixed.length >= 2 || (mixed.length === 1 && mixed[0].share >= 0.2))) {
    const lo = Math.min(...mixed.map(s => s.share)), hi = Math.max(...mixed.map(s => s.share));
    notes.push(`The reads disagree at ${mixed.length} position${mixed.length > 1 ? 's' : ''} (${mixed.slice(0, 4).map(s => s.pos).join(', ')}${mixed.length > 4 ? ', …' : ''}), where ${share(lo)}${hi - lo > 0.02 ? `–${share(hi)}` : ''} of them carry a second base. A second, closely related allele or a mixed sample is likely.`);
  }
  // Reads pooled from sister templates.
  if (g.pooled) {
    const others = g.members.slice(1).map(x => label(x.template));
    const test = (g.tests || []).find(x => x.test?.discriminating);
    const disc = test?.test;
    notes.push(`KMA divided the reads of this gene between ${others.length + 1} database alleles (${[label(g.template), ...others].join(', ')}). ${disc?.evaluable ? `At the ${disc.evaluable} position${disc.evaluable > 1 ? 's' : ''} where they differ and the reads can tell, the reads show one base each time, so they are one gene` : 'The reads do not tell them apart, so they are one gene'}, shown here with the reads of all of them.`);
  }
  // The closest allele, when the sample's is not an exact match.
  const c = g.closest;
  if (c && c.compared) {
    const ties = (c.ties || []).filter(x => x !== c.template).map(label);
    if (c.differences > 0) {
      notes.push(`Closest allele in the database: ${label(c.template)}, ${fmt(c.differences)} difference${c.differences > 1 ? 's' : ''} in ${fmt(c.compared)} compared bases${ties.length ? ` (as close: ${ties.slice(0, 3).join(', ')}${ties.length > 3 ? ', …' : ''})` : ''}. The allele in this sample is not in the database, and what it does may differ.${ev.nanopore ? ' On Nanopore reads some of the differences can be sequencing errors.' : ''}`);
    } else if (ties.length) {
      notes.push(`The reads match ${label(c.template)} and ${ties.length === 1 ? ties[0] : `${ties.length} other alleles`} equally well over the part they cover, so the allele cannot be told.`);
    }
  }
  for (const r of g.related || []) {
    const other = label(r.other.template);
    const t2 = r.test;
    if (r.relation === 'two') {
      notes.push(`A second allele of this gene is present too: ${other}. At the ${t2.mixed} position${t2.mixed > 1 ? 's' : ''} where the two differ, the reads show both bases${t2.minorShare != null ? ` (the less common one in about ${share(t2.minorShare)} of them)` : ''}.`);
    } else {
      notes.push(`At ${t2.mixed} of the ${t2.evaluable} positions where this allele and ${other} differ, some reads carry ${other}'s base${t2.minorShare != null ? ` (about ${share(t2.minorShare)} of them)` : ''}, so a minority of a second allele may be present.`);
    }
  }
  // A gene at a twentieth of the others' depth or less is most often
  // contamination (in the truth set, the calls with no locus in the genome
  // sat at 0.027 and 0.04 of the median; true chromosomal genes under a
  // high-copy plasmid's at 0.13 and 0.19), so it is marked to check.
  if (g.lowDepth && g.status === 'called' && g.lowDepth.depth < VERY_LOW_DEPTH * g.lowDepth.median) {
    warnings.push(`Its depth (${depthText(g.lowDepth.depth)}) is under a twentieth of that of the other genes found (median ${depthText(g.lowDepth.median)}): it most likely comes from contamination, or from a small minority of the cells.`);
  } else if (g.lowDepth) {
    notes.push(`Its depth (${depthText(g.lowDepth.depth)}) is far below that of the other genes found (median ${depthText(g.lowDepth.median)}), so it may be in a minority of the cells, or come from contamination.`);
  }
  return { warnings, notes, warning: warnings[0] || null };
}

const cap = (s) => s.charAt(0).toUpperCase() + s.slice(1);
