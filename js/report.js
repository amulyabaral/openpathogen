/* report.js — the results report, read top to bottom:
 *
 *   Resistance genes found, grouped by drug class, each with the antibiotics
 *   it affects (ResFinder's phenotype table) and, when the species is known,
 *   how often isolates of that species carrying the gene tested resistant
 *   (CABBAGE). Then the classes with no acquired gene, CARD matches that have
 *   no drug class, and virulence factors grouped by VFDB category.
 *
 * Everything is derived from the KMA results already in memory. Changing the
 * species re-renders the body; nothing about the sample leaves the browser.
 */

import { parseResFile } from './resfinder-db.js';
import { describeHit, classGroup, CLASS_GROUPS, antibioticPhrase, sameVariantFamily } from './genes.js';
import { loadAssociations, getManifest, linkedEvidence, resolveSpecies } from './cabbage.js';
import { antibioticClasses } from './antibiotic-classes.js';

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
const num = (v) => { const x = parseFloat(v); return Number.isFinite(x) ? x : null; };
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

function buildModel(results) {
  const m = { amrRan: [], vfRan: false, failed: [], amrFailed: false, vfFailed: false, genes: new Map(), otherCard: [], vf: [] };
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
    for (const row of parseResFile(r.resTable).rows) {
      const hit = describeHit(row.Template, r.database, row);
      Object.assign(hit, {
        section, dbLabel: label,
        identity: num(row.Template_Identity), coverage: num(row.Template_Coverage), depth: num(row.Depth),
      });
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
function groupGenes(m) {
  const groups = new Map();
  for (const e of m.genes.values()) {
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
// of 1,792 isolates, …, and to ceftaroline in 1.0% of 208." The four
// antibiotics with the most isolates tested; the CSV has all of them.
function evidenceText(species, ev, rows) {
  const parts = rows.slice(0, 4).map((r, i) => {
    const iTxt = r.i / r.n >= 0.05 ? ` (${pct(r.i, r.n)} intermediate)` : '';
    return i === 0
      ? `to ${r.antibiotic} in ${pct(r.r, r.n)} of ${fmtN(r.n)} isolates${iTxt}`
      : `to ${r.antibiotic} in ${pct(r.r, r.n)} of ${fmtN(r.n)}${iTxt}`;
  });
  const joined = parts.length > 1 ? parts.slice(0, -1).join(', ') + ', and ' + parts[parts.length - 1] : parts[0];
  return `In CABBAGE, <i>${esc(shortSpecies(species))}</i> with ${esc(ev.symbol)} tested resistant ${esc(joined)}.`;
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

export function mountReport(el, ctx) {
  // ctx: { results, readType, inputNames, sampleName, thresholds,
  //        qc: {tone, text} | null, species, onSpeciesChange, onCsv, openHit, log }
  const model = buildModel(ctx.results);
  const groups = groupGenes(model);
  const st = { species: ctx.species || '', assoc: null, assocErr: null, release: null, loading: true };
  const hasAmrGenes = model.genes.size > 0;

  renderShell();
  renderBody();
  loadAssociations()
    .then(async (a) => {
      st.assoc = a;
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
    const facts = [
      ['Species', `<input type="text" class="rep-species-input" id="rep-species" list="species-list" autocomplete="off" spellcheck="false" placeholder="Enter species (optional)" value="${esc(st.species)}" aria-label="Species"><span class="print-only rep-species-print"></span><span class="rep-species-status"></span>`],
      ['Reads', esc([READ_TYPES[ctx.readType], (ctx.inputNames || []).join(', ')].filter(Boolean).join(' · '))],
      ['Databases', esc(dbs) + (th ? ` · identity ≥ ${Math.round(th.id_threshold * 100)}%, coverage ≥ ${Math.round(th.mrc * 100)}%` : '')],
    ];
    if (ctx.qc) facts.push(['Read quality', `<span class="qc-${esc(ctx.qc.tone)}">${esc(ctx.qc.text)}</span>`]);
    facts.push(['Analyzed', esc(when)]);
    el.innerHTML = `
      <article class="report" aria-labelledby="rep-title">
        <header class="report-head">
          <h2 class="report-title" id="rep-title">Results for ${esc(sample)}</h2>
          ${model.failed.length ? `<p class="rep-alert">${esc(model.failed.join(' and '))} failed. Open the Logs panel below for details.</p>` : ''}
          <dl class="report-facts">${facts.map(([k, v]) => `<div><dt>${k}</dt><dd>${v}</dd></div>`).join('')}</dl>
        </header>
        <div class="report-body"></div>
      </article>`;
    const input = el.querySelector('#rep-species');
    input.addEventListener('change', () => {
      st.species = input.value.trim();
      ctx.onSpeciesChange?.(st.species);
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

    const evByGene = new Map();
    if (species) {
      for (const e of model.genes.values()) {
        const ev = linkedEvidence(st.assoc, species, [e.gene, ...e.hits.map(h => h.name)], e.antibiotics);
        if (ev && ev.rows.length) evByGene.set(e, ev);
      }
    }

    let html = '';
    if (model.amrRan.length) {
      html += '<section class="rep-sec"><h3 class="rep-sec-title">Resistance genes found</h3>';
      if (!hasAmrGenes) {
        html += `<p class="rep-empty">No acquired resistance genes were found with ${esc(model.amrRan.join(' and '))}.</p>`;
      } else {
        html += '<div class="rep-classes">';
        for (const grp of groups) {
          html += `<div class="rep-class"><h4 class="rep-class-name">${esc(grp.g.label)}</h4><div class="rep-genes">`;
          for (const item of grp.items) html += geneRow(grp, item, species, evByGene.get(item.e));
          html += '</div></div>';
        }
        html += '</div>';
      }
      html += cabbageStatus(species, evByGene);
      html += '</section>';

      const found = new Set(groups.map(g => g.g.key));
      const missing = CLASS_GROUPS.filter(g => !found.has(g.key));
      html += `<section class="rep-sec"><h3 class="rep-sec-title">${hasAmrGenes ? 'No acquired resistance gene found for' : 'Classes checked'}</h3>
        <p class="rep-none">${missing.map(g => esc(g.hint ? `${g.label} (${g.hint})` : g.label)).join(' · ')}</p>
        <p class="rep-caveat">Resistance can still come from chromosomal mutations, which are not checked, or be intrinsic to the species.</p>
      </section>`;

      if (model.otherCard.length) html += otherCardSection();
    } else if (model.amrFailed) {
      html += '<section class="rep-sec"><h3 class="rep-sec-title">Resistance genes</h3><p class="rep-empty">The resistance database run failed. Open the Logs panel below for details.</p></section>';
    } else {
      html += '<section class="rep-sec"><p class="rep-empty">Resistance genes were not screened: select ResFinder or CARD on the input page to include them.</p></section>';
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
    const rows = ev ? rowsForGroup(ev.rows, grp, item) : [];
    return `
      <div class="rep-gene">
        <div class="rep-gene-id">${geneBtn(e.best, e.gene)}${e.tag ? `<span class="rep-tag">${esc(e.tag)}</span>` : ''}</div>
        <div class="rep-gene-info">
          <p class="rep-drugs">${esc(drugTxt)}</p>
          ${item.first ? `<p class="rep-meta">${esc(hitMeta(e))}</p>` : ''}
          ${rows.length ? `<p class="rep-evidence">${evidenceText(species, ev, rows)}</p>` : ''}
          ${e.notes.map(n => `<p class="rep-note">ResFinder note: ${esc(n)}</p>`).join('')}
          ${item.first ? closest : ''}
        </div>
      </div>`;
  }

  function cabbageStatus(species, evByGene) {
    if (!hasAmrGenes) return '';
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
    const lines = [
      'Genes were detected by mapping reads with KMA. Drug classes and antibiotics come from the ResFinder phenotype table. This is a prediction from the genome, not a susceptibility test. Confirm results that matter for treatment in the laboratory.',
    ];
    if (withCabbage) {
      lines.push(`Tested-resistance rates are from CABBAGE (Dickens et al., NAR 2026${st.release ? `, release ${esc(st.release)}` : ''}): isolates of the same species that carry the gene, tested with 2025 CLSI or EUCAST breakpoints against antibiotics the gene is known to affect, shown when at least ${st.assoc.minN} were tested. They describe other isolates, not this sample, and include the effect of any other genes those isolates carry. To explore the data, use the <a href="https://www.ebi.ac.uk/amr" target="_blank" rel="noopener">EMBL-EBI AMR portal</a>.`);
    }
    lines.push('Generated by openpathogen, openpathogen.org.');
    return `<footer class="rep-foot">${lines.map(l => `<p>${l}</p>`).join('')}</footer>`;
  }

  // cabbage_predictions.csv in the ZIP: every linked rate behind the report.
  function emitCsv(species, evByGene) {
    if (!ctx.onCsv) return;
    if (!species || !evByGene.size) { ctx.onCsv('cabbage_predictions.csv', null); return; }
    const lines = ['species,gene,cabbage_symbol,gene_database,antibiotic,resistant,intermediate,susceptible,n,resistant_rate,non_susceptible_rate'];
    for (const [e, ev] of evByGene) {
      const dbs = [...new Set(e.hits.map(h => h.dbLabel.split(' ')[0]))].join('+');
      for (const r of ev.rows) {
        lines.push([species, e.gene, ev.symbol, dbs, r.antibiotic, r.r, r.i, r.s, r.n,
          (r.r / r.n).toFixed(4), ((r.r + r.i) / r.n).toFixed(4)].map(csvCell).join(','));
      }
    }
    ctx.onCsv('cabbage_predictions.csv', lines.join('\n') + '\n');
  }

  return { setSpecies };
}

function csvCell(v) {
  const s = String(v);
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}
