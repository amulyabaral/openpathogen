/* report.js — the results report, read top to bottom:
 *
 *   Resistance genes found, grouped by drug class, each with the antibiotics
 *   it affects (ResFinder's phenotype table) and, when the species is known,
 *   how often isolates of that species carrying the gene tested resistant
 *   (CABBAGE), next to the rate among all sequenced isolates of the species.
 *   Then the classes with no acquired gene, the genes ResFinder lists as
 *   natural to the species, partial genes, CARD matches that have no drug
 *   class, and virulence factors grouped by VFDB category.
 *
 * Which hits count is decided from the sliders (hits.js): a hit is called
 * when it passes both, partial when it passes identity but not coverage or
 * when only its well-aligned part would pass. Each called gene's alignment
 * is checked for a stretch that is missing or aligns poorly (aln-check.js),
 * the trace a truncated gene leaves even when its identity and coverage
 * look complete.
 *
 * Everything is derived from the KMA results already in memory. Changing the
 * species or the breakpoint standard re-renders the body; nothing about the
 * sample leaves the browser.
 */

import { parseResFile } from './resfinder-db.js';
import { describeHit, classGroup, CLASS_GROUPS, antibioticPhrase, sameVariantFamily, intrinsicNote } from './genes.js';
import { loadAssociations, getManifest, linkedEvidence, resolveSpecies, STANDARDS, DEFAULT_STANDARD } from './cabbage.js';
import { antibioticClasses } from './antibiotic-classes.js';
import { describeThresholds, num } from './thresholds.js';
import { classifyResult } from './hits.js';

const AMR_DBS = new Set(['resfinder', 'card_homolog']);
export const DB_NAMES = { resfinder: 'ResFinder 2.6.0', card_homolog: 'CARD 4.0.1', vfdb_core: 'VFDB set A' };
const READ_TYPES = { paired: 'Illumina paired-end', single: 'single-end', nanopore: 'Oxford Nanopore' };

// VFDB categories, most clinically telling first; others follow alphabetically.
const VF_ORDER = ['Exotoxin', 'Exoenzyme', 'Adherence', 'Invasion', 'Immune modulation', 'Effector delivery system',
  'Biofilm', 'Nutritional/Metabolic factor', 'Motility', 'Stress survival', 'Regulation',
  'Antimicrobial activity/Competitive advantage', 'Post-translational modification', 'Others'];

const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;')
  .replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const fmtN = (n) => Number(n).toLocaleString('en-US');
const pct = (num, den) => (100 * num / den).toFixed(1) + '%';
const cap1 = (s) => s.charAt(0).toUpperCase() + s.slice(1);

// "Staphylococcus aureus" → "S. aureus"
const shortSpecies = (s) => s.replace(/^([A-Z])[a-z]+ /, '$1. ');

// ── Model (species-independent) ──

// The same beta-lactamase family or van cluster, whose alleles the two
// databases name differently (blaTEM-1B / TEM-219, VanHAX / vanA).
function familyKey(gene) {
  const g = gene.toLowerCase();
  let m = g.match(/^bla([a-z]+(?:-[a-z]+)?)-?\d/);
  if (m) return 'bla:' + m[1];
  m = g.match(/^van(?:h)?([a-z])/);
  return m ? 'van:' + m[1] : null;
}

function buildModel(results, thresholds) {
  const m = {
    amrRan: [], vfRan: false, failed: [], amrFailed: false, vfFailed: false,
    genes: new Map(), otherCard: [], vf: [], partial: [],
  };
  const cardHits = [];
  results.forEach((r, section) => {
    const label = DB_NAMES[r.database] || r.dbLabel || r.database;
    if (r.exitCode !== 0) {
      m.failed.push(label);
      if (AMR_DBS.has(r.database)) m.amrFailed = true; else m.vfFailed = true;
      return;
    }
    if (AMR_DBS.has(r.database)) m.amrRan.push(label); else m.vfRan = true;
    if (!r.resTable) return;
    const sorted = classifyResult(r, thresholds);
    for (const h of [...sorted.called, ...sorted.partial]) {
      const row = h.row;
      const hit = describeHit(row.Template, r.database, row);
      Object.assign(hit, {
        section, dbLabel: label, status: h.status, by: h.by, check: h.check, warning: h.warning,
        identity: num(row.Query_Identity), templateIdentity: num(row.Template_Identity),
        coverage: num(row.Template_Coverage), depth: num(row.Depth),
      });
      if (h.status === 'partial') { m.partial.push(hit); continue; }
      if (r.database === 'vfdb_core') m.vf.push(hit);
      else if (!hit.classes.length) m.otherCard.push(hit);
      else if (r.database === 'card_homolog') cardHits.push(hit);
      else addGene(m, hit);
    }
  });
  // CARD calls with the same name as a ResFinder call merge first, so a
  // variant-name merge can never take a CARD call that has its own match.
  const rest = cardHits.filter(hit => {
    const e = m.genes.get(hit.gene.toLowerCase());
    if (e) e.hits.push(hit);
    return !e;
  });
  for (const hit of rest) addGene(m, hit);
  for (const e of m.genes.values()) {
    // Annotation from a hit with a ResFinder entry when there is one.
    const a = e.hits.find(h => h.entry) || e.hits[0];
    Object.assign(e, {
      classes: a.classes, byClass: a.byClass, antibiotics: a.antibiotics,
      spectrumKnown: a.spectrumKnown, notes: a.notes, tag: a.tag,
    });
    e.best = e.hits.slice().sort((x, y) => (y.identity ?? 0) - (x.identity ?? 0) || (y.coverage ?? 0) - (x.coverage ?? 0))[0];
    e.warning = e.best.warning || null;
  }
  return m;
}

// One entry per gene: ResFinder and CARD calls of the same gene merge. A
// CARD allele that only shares a family with a ResFinder call (TEM-219 next
// to blaTEM-1B, aph(3')-IIIa next to aph(3')-III) joins that entry if it has
// no CARD call yet: the two databases named the same locus differently.
function addGene(m, hit) {
  const key = hit.gene.toLowerCase();
  let e = m.genes.get(key);
  if (!e && hit.db === 'card_homolog') {
    const fam = familyKey(hit.gene);
    for (const f of m.genes.values()) {
      if (f.hits.some(h => h.db === 'card_homolog')) continue;
      if ((fam && familyKey(f.gene) === fam) || sameVariantFamily(f.gene, hit.gene)) { e = f; break; }
    }
  }
  if (!e) {
    e = { gene: hit.gene, hits: [] };
    m.genes.set(key, e);
  }
  e.hits.push(hit);
}

// Class groups with their genes; a gene with several classes appears in
// each, with that class's antibiotics. Identity and sources are shown once,
// where the gene first appears; main marks the gene's own (first) class.
// Entries in exclude (the species' intrinsic genes) are left out.
function groupGenes(m, exclude = new Set()) {
  const groups = new Map();
  for (const e of m.genes.values()) {
    if (exclude.has(e)) continue;
    const seen = new Map();
    for (const token of e.classes) {
      const g = classGroup(token);
      if (!groups.has(g.key)) groups.set(g.key, { g, items: [] });
      if (seen.has(g.key)) { seen.get(g.key).abs.push(...(e.byClass.get(token) || [])); continue; }
      const item = { e, abs: (e.byClass.get(token) || []).slice(), main: token === e.classes[0] };
      seen.set(g.key, item);
      groups.get(g.key).items.push(item);
    }
  }
  const order = (key) => { const i = CLASS_GROUPS.findIndex(g => g.key === key); return i < 0 ? CLASS_GROUPS.length : i; };
  const list = [...groups.values()].sort((a, b) => order(a.g.key) - order(b.g.key) || a.g.label.localeCompare(b.g.label));
  const shown = new Set();
  for (const grp of list) {
    grp.items.sort((a, b) => a.e.gene.localeCompare(b.e.gene, 'en', { sensitivity: 'base' }));
    for (const item of grp.items) {
      item.first = !shown.has(item.e);
      shown.add(item.e);
    }
  }
  return list;
}

// ── Rendering ──

const geneBtn = (hit, text) =>
  `<button type="button" class="gene-btn" data-sec="${hit.section}" data-template="${esc(hit.template)}" title="Open the gene viewer">${esc(text ?? hit.gene)}</button>`;

function hitMeta(e) {
  const b = e.best;
  const bits = [];
  if (b.identity != null) bits.push(`${b.identity.toFixed(1)}% identity`);
  if (b.coverage != null) bits.push(`${b.coverage.toFixed(1)}% coverage`);
  if (b.depth != null) bits.push(`${b.depth.toFixed(b.depth < 10 ? 1 : 0)}× depth`);
  const dbs = [];
  for (const h of e.hits) {
    const label = h.dbLabel.split(' ')[0];
    const t = h.name !== e.gene ? `${label} (${h.name})` : label;
    if (!dbs.includes(t)) dbs.push(t);
  }
  bits.push(dbs.join(', '));
  return bits.join(' · ');
}

// "In CABBAGE, S. aureus with mecA tested resistant to methicillin in 97.7%
// of 1,792 isolates, …, and to ceftaroline in 1.0% of 208. Across all
// sequenced S. aureus with a test result, 55.2% of 12,345 were resistant to
// methicillin, …" The four antibiotics with the most isolates tested; the
// CSV has all of them.
function evidenceText(species, ev, rows) {
  const top = rows.slice(0, 4);
  const parts = top.map((r, i) => {
    const iTxt = r.i / r.n >= 0.05 ? ` (${pct(r.i, r.n)} intermediate)` : '';
    return i === 0
      ? `to ${r.antibiotic} in ${pct(r.r, r.n)} of ${fmtN(r.n)} isolates${iTxt}`
      : `to ${r.antibiotic} in ${pct(r.r, r.n)} of ${fmtN(r.n)}${iTxt}`;
  });
  const join = (items) => items.length > 1 ? items.slice(0, -1).join(', ') + ', and ' + items[items.length - 1] : items[0];
  const sp = esc(shortSpecies(species));
  let text = `In CABBAGE, <i>${sp}</i> with ${esc(ev.symbol)} tested resistant ${esc(join(parts))}.`;
  const bg = top.filter(r => r.background && r.background.n >= 20)
    .map((r, i) => `${pct(r.background.r, r.background.n)} of ${fmtN(r.background.n)}${i === 0 ? ' were resistant' : ''} to ${r.antibiotic}`);
  if (bg.length) text += ` Across all sequenced <i>${sp}</i> with a test result, ${esc(join(bg))}.`;
  return text;
}

// Evidence rows that belong to a class group; rows whose antibiotic has no
// known class go with the gene's main class.
function rowsForGroup(rows, grp, item) {
  return rows.filter(r => {
    const cls = antibioticClasses(r.antibiotic);
    if (!cls) return item.main;
    return grp.g.tokens.some(t => cls.has(t));
  });
}

const ALN_ADVICE = 'The gene may be truncated or interrupted in this sample, or the reads there may come from a related sequence. Open the gene to see the alignment.';

// "73.7% of the gene at 97.6% identity · 25.7× depth · ResFinder", or for a
// hit that is partial by its alignment, the well-aligned part and what is
// wrong with the rest.
const fmtPartial = (h) => {
  const depth = h.depth != null ? `${h.depth.toFixed(h.depth < 10 ? 1 : 0)}× depth` : null;
  const db = h.dbLabel.split(' ')[0];
  if (h.by === 'alignment' && h.check) {
    return [`${Math.round(100 * h.check.goodFraction)}% of the gene at ${(100 * h.check.goodIdentity).toFixed(1)}% identity (${h.warning})`, depth, db]
      .filter(Boolean).join(' · ');
  }
  return [
    h.coverage != null ? `${h.coverage.toFixed(1)}% of the gene` + (h.identity != null ? ` at ${h.identity.toFixed(1)}% identity` : '') : null,
    depth, db,
  ].filter(Boolean).join(' · ');
};

export function mountReport(el, ctx) {
  // ctx: { results, readType, inputNames, sampleName, thresholds,
  //        qc: {tone, text} | null, host: {stats, error} | null (human read
  //        removal), species, onSpeciesChange, onCsv, openHit, log }
  const model = buildModel(ctx.results, ctx.thresholds);
  const st = { species: ctx.species || '', standard: DEFAULT_STANDARD, assoc: null, assocErr: null, release: null, loading: true };

  renderShell();
  renderBody();
  loadAssociations()
    .then(async (a) => {
      st.assoc = a;
      if (!a.standards.includes(st.standard)) st.standard = a.defaultStandard;
      try { st.release = (await getManifest()).release; } catch (_) {}
    })
    .catch((err) => { st.assocErr = err; ctx.log?.('CABBAGE data unavailable: ' + err.message, 'warn'); })
    .finally(() => { st.loading = false; renderBody(); });

  function resolved() {
    return st.assoc && st.species ? resolveSpecies(st.assoc.species, st.species) : null;
  }

  function renderShell() {
    const sample = ctx.sampleName || 'sample';
    const dbs = ctx.results.map(r => DB_NAMES[r.database] || r.dbLabel).join(', ');
    const th = ctx.thresholds;
    const when = new Date().toLocaleString('en-GB', { dateStyle: 'medium', timeStyle: 'short' });
    const standardOptions = Object.entries(STANDARDS)
      .map(([k, s]) => `<option value="${k}" title="${esc(s.detail)}"${k === st.standard ? ' selected' : ''}>${esc(s.label)}</option>`).join('');
    const facts = [
      ['Species', `<input type="text" class="rep-species-input" id="rep-species" list="species-list" autocomplete="off" spellcheck="false" placeholder="Enter species (optional)" value="${esc(st.species)}" aria-label="Species"><span class="print-only rep-species-print"></span><span class="rep-species-status"></span>`],
      ['Reads', esc([READ_TYPES[ctx.readType], (ctx.inputNames || []).join(', ')].filter(Boolean).join(' · '))],
      ['Databases', esc(dbs) + (th ? ` · ${esc(describeThresholds(th))}` : '')],
      ['Breakpoints', `<select id="rep-standard" class="rep-standard" aria-label="Breakpoint standard for the CABBAGE rates">${standardOptions}</select><span class="print-only rep-standard-print"></span>`, 'rep-fact-standard'],
    ];
    if (ctx.qc) facts.push(['Read quality', `<span class="qc-${esc(ctx.qc.tone)}">${esc(ctx.qc.text)}</span>`]);
    if (ctx.host) facts.push(['Human reads', hostFact(ctx.host)]);
    facts.push(['Analyzed', esc(when)]);
    el.innerHTML = `
      <article class="report" aria-labelledby="rep-title">
        <header class="report-head">
          <h2 class="report-title" id="rep-title">Results for ${esc(sample)}</h2>
          ${model.failed.length ? `<p class="rep-alert">${esc(model.failed.join(' and '))} failed. Open the Logs panel below for details.</p>` : ''}
          <dl class="report-facts">${facts.map(([k, v, id]) => `<div${id ? ` id="${id}" hidden` : ''}><dt>${k}</dt><dd>${v}</dd></div>`).join('')}</dl>
        </header>
        <div class="report-body"></div>
      </article>`;
    const input = el.querySelector('#rep-species');
    input.addEventListener('change', () => {
      st.species = input.value.trim();
      ctx.onSpeciesChange?.(st.species);
      renderBody();
    });
    el.querySelector('#rep-standard').addEventListener('change', (e) => {
      st.standard = e.target.value;
      renderBody();
    });
    el.addEventListener('click', (e) => {
      const b = e.target.closest('.gene-btn');
      if (b) ctx.openHit?.(Number(b.dataset.sec), b.dataset.template);
    });
  }

  // Called by the app when the species is changed on the input page.
  function setSpecies(text) {
    st.species = text || '';
    const input = el.querySelector('#rep-species');
    if (input) input.value = st.species;
    renderBody();
  }

  function renderBody() {
    const body = el.querySelector('.report-body');
    const species = resolved();
    const printEl = el.querySelector('.rep-species-print');
    if (printEl) printEl.textContent = st.species || 'not given';
    const statusEl = el.querySelector('.rep-species-status');
    if (statusEl) {
      statusEl.textContent = st.species && st.assoc && !species ? 'not in CABBAGE' : '';
    }
    // The breakpoint choice matters only once CABBAGE rates can be shown.
    const stdRow = el.querySelector('#rep-fact-standard');
    if (stdRow) {
      stdRow.hidden = !(st.assoc && st.assoc.standards.length > 1 && species);
      const print = stdRow.querySelector('.rep-standard-print');
      if (print) print.textContent = STANDARDS[st.standard]?.label || st.standard;
    }

    // Genes ResFinder lists as natural to the sample's species are set
    // apart and not counted as acquired resistance. This needs only the
    // species as typed, not CABBAGE.
    const intrinsic = new Map();
    if (st.species) {
      for (const e of model.genes.values()) {
        const note = intrinsicNote(e.notes, st.species);
        if (note) intrinsic.set(e, note);
      }
    }
    const acquired = [...model.genes.values()].filter(e => !intrinsic.has(e));
    const groups = groupGenes(model, new Set(intrinsic.keys()));

    const evByGene = new Map();
    if (species) {
      for (const e of acquired) {
        const ev = linkedEvidence(st.assoc, species, [e.gene, ...e.hits.map(h => h.name)], e.antibiotics, st.standard);
        if (ev && ev.rows.length) evByGene.set(e, ev);
      }
    }

    let html = '';
    if (model.amrRan.length) {
      html += '<section class="rep-sec"><h3 class="rep-sec-title">Resistance genes found</h3>';
      if (!acquired.length) {
        html += `<p class="rep-empty">No acquired resistance genes were found with ${esc(model.amrRan.join(' and '))}.${intrinsic.size ? ` ${intrinsic.size === 1 ? 'One gene natural to the species is' : `${intrinsic.size} genes natural to the species are`} listed below.` : ''}</p>`;
      } else {
        html += '<div class="rep-classes">';
        for (const grp of groups) {
          html += `<div class="rep-class"><h4 class="rep-class-name">${esc(grp.g.label)}</h4><div class="rep-genes">`;
          for (const item of grp.items) html += geneRow(grp, item, species, evByGene.get(item.e));
          html += '</div></div>';
        }
        html += '</div>';
      }
      html += cabbageStatus(species, evByGene, acquired.length > 0);
      html += '</section>';

      const found = new Set(groups.map(g => g.g.key));
      const missing = CLASS_GROUPS.filter(g => !found.has(g.key));
      html += `<section class="rep-sec"><h3 class="rep-sec-title">${acquired.length ? 'No acquired resistance gene found for' : 'Classes checked'}</h3>
        <p class="rep-none">${missing.map(g => esc(g.hint ? `${g.label} (${g.hint})` : g.label)).join(' · ')}</p>
        <p class="rep-caveat">Resistance can still come from chromosomal mutations, which are not checked, or be intrinsic to the species.</p>
      </section>`;

      if (intrinsic.size) html += intrinsicSection(intrinsic, species);
      if (model.partial.length) html += partialSection();
      if (model.otherCard.length) html += otherCardSection();
    } else if (model.amrFailed) {
      html += '<section class="rep-sec"><h3 class="rep-sec-title">Resistance genes</h3><p class="rep-empty">The resistance database run failed. Open the Logs panel below for details.</p></section>';
    } else {
      html += '<section class="rep-sec"><p class="rep-empty">Resistance genes were not screened: select ResFinder or CARD on the input page to include them.</p></section>';
      if (model.partial.length) html += partialSection();
    }

    html += virulenceSection();
    html += footer(evByGene.size > 0);
    body.innerHTML = html;
    emitCsv(species, evByGene);
  }

  function geneRow(grp, item, species, ev) {
    const e = item.e;
    const drugs = antibioticPhrase(grp.g.key, item.abs);
    const drugTxt = drugs
      ? cap1(drugs) + (e.spectrumKnown ? '' : ' (not all antibiotics are listed for this allele)')
      : 'Antibiotics not listed for this allele';
    const closest = grp.g.key === 'beta-lactam' && e.best.identity != null && e.best.identity < 100
      ? `<p class="rep-note">Closest known allele (${e.best.identity.toFixed(1)}% identity). The spectrum of the allele in this sample may differ.</p>` : '';
    const warning = e.warning ? `<p class="rep-note rep-warn">Alignment check: ${esc(e.warning)} ${ALN_ADVICE}</p>` : '';
    const rows = ev ? rowsForGroup(ev.rows, grp, item) : [];
    return `
      <div class="rep-gene">
        <div class="rep-gene-id">${geneBtn(e.best, e.gene)}${e.tag ? `<span class="rep-tag">${esc(e.tag)}</span>` : ''}${e.warning ? '<span class="rep-tag rep-tag-warn">check alignment</span>' : ''}</div>
        <div class="rep-gene-info">
          <p class="rep-drugs">${esc(drugTxt)}</p>
          ${item.first ? `<p class="rep-meta">${esc(hitMeta(e))}</p>` : ''}
          ${item.first ? warning : ''}
          ${rows.length ? `<p class="rep-evidence">${evidenceText(species, ev, rows)}</p>` : ''}
          ${e.notes.map(n => `<p class="rep-note">ResFinder note: ${esc(n)}</p>`).join('')}
          ${item.first ? closest : ''}
        </div>
      </div>`;
  }

  // Genes ResFinder marks "Natural in" or "Chromosomal in" the sample's
  // species: shown with their antibiotics and caveats, outside the classes.
  function intrinsicSection(intrinsic, species) {
    const name = species || st.species;
    let html = `<section class="rep-sec"><h3 class="rep-sec-title">Genes natural to <i>${esc(shortSpecies(name))}</i></h3>
      <p class="rep-caveat">ResFinder lists these as part of the species' own genome, so they are not counted as acquired resistance above. They can still matter: read each note.</p>
      <div class="rep-genes">`;
    const entries = [...intrinsic.keys()].sort((a, b) => a.gene.localeCompare(b.gene, 'en', { sensitivity: 'base' }));
    for (const e of entries) {
      const groupsOf = [...new Set(e.classes.map(t => classGroup(t).key))];
      const drugs = groupsOf.map(key => antibioticPhrase(key, e.antibiotics)).filter(Boolean).join('; ');
      const drugTxt = drugs ? cap1(drugs) : 'Antibiotics not listed for this allele';
      const warning = e.warning ? `<p class="rep-note rep-warn">Alignment check: ${esc(e.warning)} ${ALN_ADVICE}</p>` : '';
      html += `
      <div class="rep-gene">
        <div class="rep-gene-id">${geneBtn(e.best, e.gene)}<span class="rep-tag rep-tag-muted">intrinsic</span></div>
        <div class="rep-gene-info">
          <p class="rep-drugs">${esc(drugTxt)}</p>
          <p class="rep-meta">${esc(hitMeta(e))}</p>
          ${warning}
          ${e.notes.map(n => `<p class="rep-note">ResFinder note: ${esc(n)}</p>`).join('')}
        </div>
      </div>`;
    }
    return html + '</div></section>';
  }

  // Hits that pass the identity threshold but cover less of the gene than
  // the coverage threshold asks.
  function partialSection() {
    const hits = model.partial.slice().sort((a, b) => (b.coverage ?? 0) - (a.coverage ?? 0) || a.name.localeCompare(b.name));
    const shown = hits.slice(0, 15);
    const more = hits.length > shown.length ? `<p class="rep-caveat">${hits.length - shown.length} more in the gene tables.</p>` : '';
    return `<section class="rep-sec"><h3 class="rep-sec-title">Partial genes</h3>
      <p class="rep-caveat">Hits whose well-aligned part reaches the identity threshold but which cover less of the gene than the coverage threshold asks, or whose remainder aligns poorly. A gene can be partial because it is truncated or interrupted in this sample, because coverage is thin, or because only part of it is shared with a related gene. They are not counted above.</p>
      <div class="rep-partials">${shown.map(h => `<p class="rep-partial">${geneBtn(h, h.name)}<span class="rep-tag rep-tag-muted">partial</span><span class="rep-meta">${esc(fmtPartial(h))}</span></p>`).join('')}</div>
      ${more}
    </section>`;
  }

  function cabbageStatus(species, evByGene, anyAcquired) {
    if (!anyAcquired) return '';
    let text = '';
    if (!st.species) text = 'Enter the species above to add tested-resistance rates from CABBAGE.';
    else if (st.loading) text = 'Loading CABBAGE data…';
    else if (st.assocErr) text = `CABBAGE data could not be loaded (${st.assocErr.message}).`;
    else if (!species) text = `CABBAGE has no data for “${st.species}”. It covers ${st.assoc.species.length} species, listed in the species box.`;
    else if (!evByGene.size) text = `CABBAGE has fewer than ${st.assoc.minN} tested ${species} isolates carrying these genes.`;
    return text ? `<p class="rep-cabbage-note">${esc(text)}</p>` : '';
  }

  function otherCardSection() {
    const hits = model.otherCard.slice().sort((a, b) => a.name.localeCompare(b.name));
    const shown = hits.slice(0, 12).map(h => geneBtn(h, h.name)).join(', ');
    const more = hits.length > 12 ? ` and ${hits.length - 12} more in the CARD table` : '';
    return `<section class="rep-sec"><h3 class="rep-sec-title">Other CARD matches</h3>
      <p class="rep-other">${shown}${esc(more)}</p>
      <p class="rep-caveat">ResFinder has no drug class for these genes. Many are efflux pumps, regulators or genes found in most isolates of a species, so their presence alone does not mean resistance. Each gene's CARD entry is linked in the gene viewer.</p>
    </section>`;
  }

  function virulenceSection() {
    if (model.vfFailed) {
      return '<section class="rep-sec"><h3 class="rep-sec-title">Virulence factors</h3><p class="rep-empty">The VFDB run failed. Open the Logs panel below for details.</p></section>';
    }
    if (!model.vfRan) {
      return model.amrRan.length
        ? '<section class="rep-sec"><h3 class="rep-sec-title">Virulence factors</h3><p class="rep-empty">Not screened. Select VFDB on the input page to include virulence factors.</p></section>'
        : '';
    }
    if (!model.vf.length) {
      return '<section class="rep-sec"><h3 class="rep-sec-title">Virulence factors</h3><p class="rep-empty">No virulence factors were found with VFDB.</p></section>';
    }
    const cats = new Map();
    for (const h of model.vf) {
      const cat = h.vf?.category || 'Other';
      const factor = h.vf?.factor || h.name;
      if (!cats.has(cat)) cats.set(cat, new Map());
      const f = cats.get(cat);
      if (!f.has(factor)) f.set(factor, []);
      if (!f.get(factor).some(x => x.name === h.name)) f.get(factor).push(h);
    }
    const rank = (c) => { const i = VF_ORDER.indexOf(c); return i < 0 ? VF_ORDER.length : i; };
    const ordered = [...cats.entries()].sort((a, b) => rank(a[0]) - rank(b[0]) || a[0].localeCompare(b[0]));
    let html = '<section class="rep-sec"><h3 class="rep-sec-title">Virulence factors</h3><div class="rep-vf">';
    for (const [cat, factors] of ordered) {
      const items = [...factors.entries()].map(([factor, hits]) => {
        const shown = hits.slice(0, 6).map(h => geneBtn(h, h.name)).join(', ');
        const more = hits.length > 6 ? `, +${hits.length - 6} more` : '';
        return `<span class="rep-factor">${esc(factor)} (${shown}${more})</span>`;
      });
      html += `<div class="rep-vf-row"><h4 class="rep-vf-cat">${esc(cat)}</h4><p class="rep-vf-list">${items.join('<span class="rep-sep"> · </span>')}</p></div>`;
    }
    return html + '</div></section>';
  }

  function footer(withCabbage) {
    const th = ctx.thresholds ? describeThresholds(ctx.thresholds) : null;
    const lines = [
      `Genes were detected by mapping reads with KMA, with ResFinder's settings.${th ? ` A gene is called when it reaches ${th}; a hit at the identity threshold but below the coverage threshold is listed as partial.` : ''} Drug classes and antibiotics come from the ResFinder phenotype table. This is a prediction from the genome, not a susceptibility test. Confirm results that matter for treatment in the laboratory.`,
    ];
    if (ctx.host?.stats) {
      lines.push('Human reads were removed before gene detection with Deacon 0.18.0 and its panhuman-1 index (k=31, w=61), on this device.');
    }
    if (withCabbage) {
      const s = STANDARDS[st.standard] || { label: st.standard, detail: '' };
      lines.push(`Tested-resistance rates are from CABBAGE (Dickens et al., NAR 2026${st.release ? `, release ${esc(st.release)}` : ''}): isolates of the same species that carry the gene, tested against antibiotics the gene is known to affect, read with ${esc(s.detail || s.label)}, shown when at least ${st.assoc.minN} were tested. The rate across all sequenced isolates of the species is given for comparison. They describe other isolates, not this sample, and include the effect of any other genes those isolates carry. To explore the data, use the <a href="https://www.ebi.ac.uk/amr" target="_blank" rel="noopener">EMBL-EBI AMR portal</a>.`);
    }
    lines.push('Generated by openpathogen, openpathogen.org.');
    return `<footer class="rep-foot">${lines.map(l => `<p>${l}</p>`).join('')}</footer>`;
  }

  // cabbage_predictions.csv in the ZIP: every linked rate behind the report,
  // with the species background and the breakpoint standard.
  function emitCsv(species, evByGene) {
    if (!ctx.onCsv) return;
    if (!species || !evByGene.size) { ctx.onCsv('cabbage_predictions.csv', null); return; }
    const lines = ['species,gene,cabbage_symbol,gene_database,antibiotic,resistant,intermediate,susceptible,n,resistant_rate,non_susceptible_rate,background_n,background_resistant_rate,breakpoints'];
    for (const [e, ev] of evByGene) {
      const dbs = [...new Set(e.hits.map(h => h.dbLabel.split(' ')[0]))].join('+');
      for (const r of ev.rows) {
        const bg = r.background;
        lines.push([species, e.gene, ev.symbol, dbs, r.antibiotic, r.r, r.i, r.s, r.n,
          (r.r / r.n).toFixed(4), ((r.r + r.i) / r.n).toFixed(4),
          bg ? bg.n : '', bg ? (bg.r / bg.n).toFixed(4) : '', st.standard].map(csvCell).join(','));
      }
    }
    ctx.onCsv('cabbage_predictions.csv', lines.join('\n') + '\n');
  }

  return { setSpecies };
}

// "504,504 of 512,196 reads removed (98.5%)", or the failure.
function hostFact({ stats, error }) {
  if (!stats) return `<span class="qc-warn">Removal failed${error ? ` (${esc(error)})` : ''}. Results are from all reads.</span>`;
  const removed = stats.readsIn - stats.readsOut;
  const p = stats.readsIn ? removed / stats.readsIn : 0;
  const share = p > 0 && p < 0.001 ? '<0.1%' : (100 * p).toFixed(1) + '%';
  return esc(`${fmtN(removed)} of ${fmtN(stats.readsIn)} reads removed (${share}) with Deacon`);
}

function csvCell(v) {
  const s = String(v);
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}
