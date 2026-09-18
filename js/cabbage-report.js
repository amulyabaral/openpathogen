/* cabbage-report.js — the "Phenotype associations (CABBAGE)" card that mounts
 * into the analysis results. Two parts:
 *
 *  1. Gene→phenotype associations: for the sample's species, every detected
 *     resistance gene is matched to its CABBAGE (AMRFinderPlus) symbol and
 *     looked up in the precomputed association table: "S. aureus isolates
 *     carrying mecA: 97.7% methicillin resistant (n=1,791)", next to the
 *     rate in all sequenced isolates of the species. Antibiotics the portal
 *     links to the gene are listed first; other antibiotics are shown as
 *     co-occurrence, because the gene does not act on them.
 *
 *  2. Species antibiogram: per-antibiotic resistance rates across all of
 *     CABBAGE's AST records for the selected species, optionally refined by
 *     country / year / isolation source (loads the phenotype snapshot on
 *     demand).
 *
 * The association table loads automatically with the card; nothing about
 * the sample is sent anywhere.
 */

import {
  loadAssociations, loadView, isViewLoaded, getManifest,
  predictPhenotypes, unmatchedGenes, backgroundAntibiogram, antibiogram,
  resolveSpecies, templateGene, geneInfo,
} from './cabbage.js';
import { lookupGenePhenotype } from './resfinder-db.js';
import { antibioticClasses, geneClasses, sameClass, linkedAntibiotics } from './antibiotic-classes.js';

const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;')
  .replace(/>/g, '&gt;').replace(/"/g, '&quot;');

// Exact counts with separators ("1,791"), never rounded to "1.8k".
const fmtN = (n) => Number(n).toLocaleString('en-US');
const pct = (num, den) => (100 * num / den).toFixed(1) + '%';

// Only the AMR databases feed the card: virulence factors have no antibiotic
// phenotype, and VFDB template names are not gene symbols.
const AMR_DBS = { resfinder: 'ResFinder', card_homolog: 'CARD' };

const PRED_PREVIEW = 12;
const TITLE = 'Phenotype associations (CABBAGE)';

const GROUPS = {
  match: {
    title: 'Antibiotics the detected genes act on',
    note: 'CABBAGE links these antibiotics to the gene, so the gene itself can explain the rate.',
  },
  unknown: {
    title: 'Antibiotics with no known link to the detected genes',
    note: 'CABBAGE records no target antibiotic for these genes. The rate may or may not be caused by the gene.',
  },
  cooc: {
    title: 'Other antibiotics: co-occurrence only',
    note: 'The detected genes do not act on these antibiotics. These rates show what else the same lineages carry.',
  },
};

// Extract gene symbols (and their ResFinder drug classes, used as a fallback
// when CABBAGE has no link for a gene) from the run's KMA results.
export function collectDetectedGenes(results, parseRes) {
  const out = [];
  const seen = new Set();
  for (const r of results || []) {
    if (!r.resTable || !AMR_DBS[r.database]) continue;
    const rows = parseRes(r.resTable).rows;
    for (const row of rows) {
      const gene = templateGene(row.Template, r.database);
      if (!gene) continue;
      const key = gene.toLowerCase();
      if (seen.has(key)) continue;
      seen.add(key);
      const pheno = row._phenotype || lookupGenePhenotype(gene) || lookupGenePhenotype('bla' + gene);
      out.push({ name: gene, db: AMR_DBS[r.database], classes: geneClasses(pheno?.class) });
    }
  }
  return out;
}

export function mountCabbageReport(container, { organism, detectedGenes, onCsv, log }) {
  const state = {
    assoc: null,
    release: null,
    species: null,
    phenotypeView: null,       // lazy; enables filtered antibiogram
    filters: { country: '', sourceCategory: '', yearFrom: '', yearTo: '' },
    refineOpen: false,
  };
  const classesByGene = new Map((detectedGenes || []).map(g => [g.name, g.classes || new Set()]));

  const push = (msg, level) => log?.(msg, level);

  renderLoading();
  (async () => {
    try {
      state.assoc = await loadAssociations((got, total) => {
        const p = total ? Math.round(100 * got / total) : 0;
        setStatus(`Loading CABBAGE data… ${p ? p + '%' : ''}`);
      });
      try { state.release = (await getManifest()).release; } catch (_) {}
      state.species = resolveSpecies(state.assoc.species, organism);
      push(`CABBAGE ${state.release || ''}: ${fmtN(state.assoc.bySpeciesGene.size)} gene–antibiotic associations loaded`, 'ok');
      render();
      emitCsv();
    } catch (err) {
      renderError(err);
    }
  })();

  function setStatus(text) {
    const el = container.querySelector('.cabb-status');
    if (el) el.textContent = text;
  }

  function renderLoading() {
    container.innerHTML = `
      <div class="card">
        <div class="card-title comp-title"><span>${TITLE}</span><span class="comp-src">AMR genotype–phenotype database</span></div>
        <div class="card-body"><p class="cabb-status opt-note">Loading CABBAGE data…</p></div>
      </div>`;
  }

  function renderError(err) {
    container.innerHTML = `
      <div class="card">
        <div class="card-title comp-title"><span>${TITLE}</span><span class="comp-src">AMR genotype–phenotype database</span></div>
        <div class="card-body">
          <p class="opt-note">Could not load the CABBAGE data (${esc(err.message)}). The gene results above are not affected.</p>
        </div>
      </div>`;
    push('CABBAGE data unavailable: ' + err.message, 'warn');
  }

  function render() {
    const { assoc } = state;
    const options = assoc.species
      .map(s => `<option value="${esc(s)}"${s === state.species ? ' selected' : ''}>${esc(s)}</option>`)
      .join('');

    container.innerHTML = `
      <div class="card cabb-card" id="cabb-card">
        <div class="card-title comp-title">
          <span>${TITLE}</span>
          <span class="comp-src">${state.release ? 'release ' + esc(state.release) + ' · ' : ''}Dickens et al., <a href="https://doi.org/10.1093/nar/gkag780" target="_blank" rel="noopener">NAR 2026</a></span>
        </div>
        <div class="card-body">
          <div class="cabb-head-row">
            <label class="cabb-species-label" for="cabb-species">Species</label>
            <select id="cabb-species">${state.species ? '' : '<option value="">— select species —</option>'}${options}</select>
            <span class="cabb-isolates" id="cabb-isolates"></span>
          </div>
          <div id="cabb-predictions"></div>
          <div class="cabb-ab-head">
            <span class="cabb-ab-title">Species antibiogram: resistance rates in all CABBAGE records</span>
            <button class="cabb-refine-btn" id="cabb-refine" type="button">Filter by country / year / source</button>
          </div>
          <div id="cabb-abfilters" hidden></div>
          <div id="cabb-antibiogram"></div>
          <p class="cabb-footnote">These are associations from CABBAGE, not a susceptibility test. For each gene, the dark bar is the share of isolates of this species carrying the gene that tested resistant, and the grey bar the share that tested intermediate, using 2025 CLSI or EUCAST breakpoints. The background is all sequenced isolates of the species. Only associations with at least ${state.assoc.minN} isolates are shown. Confirm important results in the laboratory. Snapshot of the <a href="https://www.ebi.ac.uk/amr" target="_blank" rel="noopener">EMBL-EBI AMR portal</a>, computed in your browser.</p>
        </div>
      </div>`;

    container.querySelector('#cabb-species').addEventListener('change', (e) => {
      state.species = e.target.value || null;
      renderPredictions();
      renderAntibiogram();
      emitCsv();
    });
    container.querySelector('#cabb-refine').addEventListener('click', () => {
      state.refineOpen = !state.refineOpen;
      if (state.refineOpen) enableFilters();
      else { container.querySelector('#cabb-abfilters').hidden = true; renderAntibiogram(); }
    });

    renderPredictions();
    renderAntibiogram();
  }

  // ── Section 1: gene → phenotype associations ──

  // true: CABBAGE links the gene to this antibiotic; false: the gene has
  // links but not this one; null: CABBAGE has no links for the gene, so
  // fall back to ResFinder's drug class where known.
  function linkStatus(g, antibiotic) {
    const info = geneInfo(state.assoc, g.gene);
    const ab = antibiotic.toLowerCase();
    if (info.links.size) {
      for (const l of info.links) if (linkedAntibiotics(l).has(ab)) return true;
      return false;
    }
    const abCls = antibioticClasses(antibiotic);
    const gCls = classesByGene.get(g.detectedAs);
    if (!abCls || !gCls || !gCls.size) return null;
    return sameClass(gCls, abCls);
  }

  function classify(p) {
    const status = p.genes.map(g => linkStatus(g, p.antibiotic));
    const linked = p.genes.filter((g, i) => status[i] === true);
    if (linked.length) return { kind: 'match', genes: linked };
    return { kind: status.some(s => s === false) ? 'cooc' : 'unknown', genes: [] };
  }

  // For a linked antibiotic only the linked genes count: the bar, the
  // verdict and the gene list come from them, so a co-occurring gene cannot
  // set the verdict. Other rows use every detected gene.
  function groupedPredictions() {
    const preds = predictPhenotypes(state.assoc, state.species, detectedGenes);
    const order = { match: 0, unknown: 1, cooc: 2 };
    const byRate = (a, b) => b.r / b.n - a.r / a.n;
    return preds
      .map(p => {
        const cls = classify(p);
        cls.genes = cls.genes.slice().sort(byRate);
        const pool = (cls.kind === 'match' ? cls.genes : p.genes).slice().sort(byRate);
        const driver = pool[0];
        const rate = driver.r / driver.n;
        const rateNS = pool.reduce((m, g) => Math.max(m, (g.r + g.i) / g.n), 0);
        const verdict = rate >= 0.85 ? 'resistant' : rateNS >= 0.4 ? 'uncertain' : 'susceptible';
        return { ...p, cls, pool, driver, rate, rateNS, verdict };
      })
      .sort((a, b) => order[a.cls.kind] - order[b.cls.kind]);
  }

  // Stacked bar: resistant (dark) then intermediate (grey), label = % resistant.
  function bars(r, i, n, extraClass, title) {
    const showI = i / n >= 0.005;
    return `
      <div class="cabb-bar-row${extraClass ? ' ' + extraClass : ''}" title="${esc(title)}">
        <span class="cabb-bar-track"><span class="cabb-bar-fill" style="width:${Math.round(100 * r / n)}%"></span><span class="cabb-bar-fill i" style="width:${Math.round(100 * i / n)}%"></span></span>
        <span class="cabb-bar-val">${pct(r, n)}${showI ? `<small class="cabb-i">+${pct(i, n)} I</small>` : ''}</span>
      </div>`;
  }

  function renderPredictions() {
    const area = container.querySelector('#cabb-predictions');
    const isoEl = container.querySelector('#cabb-isolates');
    if (isoEl) isoEl.textContent = '';
    if (!state.species) {
      area.innerHTML = detectedGenes.length
        ? `<p class="opt-note">Select the sample's species to look up the ${detectedGenes.length} detected resistance gene${detectedGenes.length === 1 ? '' : 's'} in CABBAGE.</p>`
        : '';
      return;
    }

    const preds = groupedPredictions();
    const missing = unmatchedGenes(state.assoc, state.species, detectedGenes);

    let html = '';
    if (!detectedGenes.length) {
      html += `<p class="opt-note">No resistance genes were detected, so there is nothing to look up. The antibiogram below shows the background resistance rates of this species.</p>`;
    } else if (!preds.length) {
      html += `<p class="opt-note">None of the detected genes has CABBAGE data for <b>${esc(state.species)}</b> with at least ${state.assoc.minN} isolates.</p>`;
    } else {
      html += '<div class="cabb-preds">';
      let lastKind = null;
      for (const p of preds) {
        if (p.cls.kind !== lastKind) {
          lastKind = p.cls.kind;
          html += `<div class="cabb-group-title" data-group="${p.cls.kind}">${esc(GROUPS[p.cls.kind].title)}</div>
            <p class="cabb-group-note" data-group="${p.cls.kind}">${esc(GROUPS[p.cls.kind].note)}</p>`;
        }
        const geneChips = p.pool.map(g => {
          const detectedNote = g.detectedAs.toLowerCase() !== g.symbol.toLowerCase() ? ` (detected as ${esc(g.detectedAs)})` : '';
          const iTxt = g.i / g.n >= 0.005 ? `, ${pct(g.i, g.n)} intermediate` : '';
          return `<span class="cabb-gene" title="${esc(g.db)} hit${detectedNote}">${esc(g.symbol)} <b>${pct(g.r, g.n)}</b> resistant${iTxt} (n=${fmtN(g.n)})</span>`;
        }).join('');
        const bgTxt = p.background
          ? `all sequenced ${state.species} isolates: ${pct(p.background.r, p.background.n)} resistant (n=${fmtN(p.background.n)})`
          : 'no background data';
        const verdictTxt = p.verdict === 'resistant' ? 'predicted resistant'
          : p.verdict === 'uncertain' ? 'uncertain' : 'likely susceptible';
        const tag = p.cls.kind === 'match'
          ? `<span class="cabb-tag cabb-tag-match" title="CABBAGE links this antibiotic to the gene">linked: ${esc(p.cls.genes.map(g => g.symbol).join(', '))}</span>`
          : p.cls.kind === 'cooc'
            ? '<span class="cabb-tag" title="the detected genes do not act on this antibiotic">co-occurrence</span>'
            : '';
        html += `
          <div class="cabb-pred cabb-verdict-${p.verdict}${p.cls.kind === 'cooc' ? ' cabb-cooc' : ''}" data-group="${p.cls.kind}">
            <span class="cabb-verdict-pill" title="based on the strongest gene association for this antibiotic">${verdictTxt}</span>
            <span class="cabb-ab-cell"><span class="cabb-ab-name" title="${esc(p.antibiotic)}">${esc(p.antibiotic)}</span>${tag}</span>
            <div class="cabb-bars">
              ${bars(p.driver.r, p.driver.i, p.driver.n, '', 'isolates of this species carrying ' + p.driver.symbol)}
              ${p.background ? bars(p.background.r, p.background.i, p.background.n, 'bg', 'all sequenced isolates of this species') : ''}
            </div>
            <span class="cabb-genes">${geneChips}</span>
            <span class="cabb-bg-note">${esc(bgTxt)}</span>
          </div>`;
      }
      html += '</div>';
      if (preds.length > PRED_PREVIEW) {
        html += `<button class="cabb-preds-more" type="button">Show all ${preds.length} antibiotics</button>`;
      }
    }

    if (missing.length) {
      html += `<p class="cabb-missing">No CABBAGE data for ${missing.map(g => `<code title="${esc(g.db)} hit">${esc(g.name)}</code>`).join(', ')}: fewer than ${state.assoc.minN} tested ${esc(state.species)} isolates carry the gene, or the name has no match in CABBAGE.</p>`;
    }
    area.innerHTML = html;

    // Cap the visible rows until the user expands (same pattern as the gene
    // tables); a group heading whose rows are all hidden hides with them.
    const rows = [...area.querySelectorAll('.cabb-pred')];
    const extra = rows.slice(PRED_PREVIEW);
    extra.forEach(el => el.classList.add('row-limited'));
    const hiddenHeads = [];
    for (const head of area.querySelectorAll('.cabb-group-title')) {
      const kind = head.dataset.group;
      if (rows.some(r => r.dataset.group === kind && !r.classList.contains('row-limited'))) continue;
      const note = head.nextElementSibling;
      head.hidden = true; note.hidden = true;
      hiddenHeads.push(head, note);
    }
    area.querySelector('.cabb-preds-more')?.addEventListener('click', (e) => {
      extra.forEach(el => el.classList.remove('row-limited'));
      hiddenHeads.forEach(el => { el.hidden = false; });
      e.target.remove();
    });
  }

  // ── Section 2: species antibiogram (background or filtered) ──

  function currentFilters() {
    const f = state.filters;
    const active = !!(state.phenotypeView && (f.country || f.sourceCategory || f.yearFrom || f.yearTo));
    return {
      active,
      opts: {
        country: f.country || null,
        sourceCategory: f.sourceCategory || null,
        yearFrom: f.yearFrom ? Number(f.yearFrom) : null,
        yearTo: f.yearTo ? Number(f.yearTo) : null,
      },
    };
  }

  function antibiogramRows() {
    const { active, opts } = currentFilters();
    return {
      filtered: active,
      rows: active
        ? antibiogram(state.phenotypeView, state.species, opts)
        : backgroundAntibiogram(state.assoc, state.species),
    };
  }

  function renderAntibiogram() {
    const area = container.querySelector('#cabb-antibiogram');
    if (!state.species) { area.innerHTML = ''; return; }
    const { filtered, rows } = antibiogramRows();

    const isoEl = container.querySelector('#cabb-isolates');
    if (isoEl) {
      const total = rows.reduce((n, r) => Math.max(n, r.n), 0);
      isoEl.textContent = total ? `up to ${fmtN(total)} records per antibiotic` : '';
    }

    if (!rows.length) { area.innerHTML = '<p class="opt-note">No AST records for this species' + (filtered ? ' with these filters' : '') + '.</p>'; return; }

    const max = Math.max(...rows.map(r => (r.r + r.i) / r.n), 0.001);
    let html = '<div class="cabb-abtable">';
    for (const r of rows.slice(0, 25)) {
      const rns = (r.r + r.i) / r.n;
      const w = Math.max(1, Math.round(100 * rns / max));
      html += `
        <div class="cabb-abrow">
          <span class="cabb-ab-name" title="${esc(r.antibiotic)}">${esc(r.antibiotic)}</span>
          <span class="cabb-bar-track"><span class="cabb-bar-fill${filtered ? '' : ' bg'}" style="width:${w}%"></span></span>
          <span class="cabb-bar-val">${(100 * rns).toFixed(1)}%</span>
          <span class="cabb-bar-n">n=${fmtN(r.n)}</span>
        </div>`;
    }
    html += '</div>';
    if (rows.length > 25) html += `<p class="opt-note">${rows.length - 25} more antibiotics are in the results ZIP.</p>`;
    area.innerHTML = html;
  }

  // ── Filters (require the full phenotype view) ──

  async function enableFilters() {
    const box = container.querySelector('#cabb-abfilters');
    box.hidden = false;
    const btn = container.querySelector('#cabb-refine');
    if (state.phenotypeView) { buildFilterControls(box); return; }
    box.innerHTML = '<p class="opt-note">Loading the full CABBAGE phenotype table (5.7 MB, once)…</p>';
    btn.disabled = true;
    try {
      if (!isViewLoaded('phenotypes')) {
        state.phenotypeView = await loadView('phenotypes', (got, total) => {
          const p = total ? Math.round(100 * got / total) : 0;
          box.innerHTML = `<p class="opt-note">Downloading CABBAGE phenotypes… ${p}%</p>`;
        });
        push('CABBAGE phenotype table loaded; antibiogram filters enabled', 'ok');
      } else {
        state.phenotypeView = await loadView('phenotypes');
      }
      buildFilterControls(box);
    } catch (err) {
      box.innerHTML = `<p class="opt-note">Could not load the phenotype table (${esc(err.message)}). The unfiltered antibiogram is shown.</p>`;
      state.refineOpen = false;
    } finally {
      btn.disabled = false;
    }
  }

  function buildFilterControls(box) {
    const view = state.phenotypeView;
    const opts = (colId, label) => {
      const col = view.colById.get(colId);
      if (!col) return '';
      const values = col.dict.slice(1).filter(Boolean).sort((a, b) => a.localeCompare(b));
      return `<select data-filter="${label}" aria-label="${label === 'country' ? 'Country' : 'Isolation source'}"><option value="">all ${label === 'country' ? 'countries' : label === 'sourceCategory' ? 'source categories' : ''}</option>${values.map(v => `<option value="${esc(v)}"${state.filters[label] === v ? ' selected' : ''}>${esc(v)}</option>`).join('')}</select>`;
    };
    let years = [];
    const yr = view.colById.get('phenotype-collection_year');
    if (yr) {
      const set = new Set();
      for (let i = 0; i < yr.vals.length; i++) if (!Number.isNaN(yr.vals[i])) set.add(yr.vals[i]);
      years = [...set].sort((a, b) => a - b);
    }
    const yearOpts = (sel, key) => `<option value="">${sel}</option>` + years.map(y => `<option value="${y}"${String(state.filters[key]) === String(y) ? ' selected' : ''}>${y}</option>`).join('');

    box.innerHTML = `
      <div class="cabb-filters">
        ${opts('phenotype-country', 'country')}
        ${opts('phenotype-isolation_source_category', 'sourceCategory')}
        <select data-filter="yearFrom" aria-label="Year from">${yearOpts('year from', 'yearFrom')}</select>
        <select data-filter="yearTo" aria-label="Year to">${yearOpts('year to', 'yearTo')}</select>
        <button type="button" class="cabb-clear-filters" id="cabb-clear-filters">clear</button>
      </div>`;
    box.querySelectorAll('select[data-filter]').forEach(sel => {
      sel.addEventListener('change', () => {
        state.filters[sel.dataset.filter] = sel.value;
        renderAntibiogram();
        emitCsv();
      });
    });
    box.querySelector('#cabb-clear-filters')?.addEventListener('click', () => {
      state.filters = { country: '', sourceCategory: '', yearFrom: '', yearTo: '' };
      buildFilterControls(box);
      renderAntibiogram();
      emitCsv();
    });
  }

  // ── CSV for the results ZIP ──

  function emitCsv() {
    if (!onCsv || !state.species) return;
    const preds = groupedPredictions();
    const lines = ['species,gene,detected_as,gene_database,antibiotic,link,resistant,intermediate,susceptible,n,resistant_rate,non_susceptible_rate'];
    for (const p of preds) {
      for (const g of p.genes) {
        const st = linkStatus(g, p.antibiotic);
        const link = st === true ? 'linked' : st === false ? 'co-occurrence' : 'unknown';
        lines.push([state.species, g.symbol, g.detectedAs, g.db, JSON.stringify(p.antibiotic), link, g.r, g.i, g.s, g.n, (g.r / g.n).toFixed(4), ((g.r + g.i) / g.n).toFixed(4)].join(','));
      }
    }
    const bgSeq = (state.assoc.backgroundSeq || []).filter(b => b.species === state.species);
    if (bgSeq.length) {
      lines.push('');
      lines.push('# background: all sequenced isolates of the species with an antibiogram');
      lines.push('species,antibiotic,resistant,intermediate,susceptible,n,resistant_rate,non_susceptible_rate');
      for (const b of bgSeq) {
        const n = b.r + b.i + b.s;
        lines.push([state.species, JSON.stringify(b.antibiotic), b.r, b.i, b.s, n, (b.r / n).toFixed(4), ((b.r + b.i) / n).toFixed(4)].join(','));
      }
    }
    const { filtered, rows } = antibiogramRows();
    lines.push('');
    lines.push('# species antibiogram (all CABBAGE AST records' + (filtered ? ', filtered' : '') + ')');
    lines.push('species,antibiotic,resistant,intermediate,susceptible,n,resistant_rate,non_susceptible_rate');
    for (const r of rows) {
      lines.push([state.species, JSON.stringify(r.antibiotic), r.r, r.i, r.s, r.n, (r.r / r.n).toFixed(4), ((r.r + r.i) / r.n).toFixed(4)].join(','));
    }
    onCsv('cabbage_predictions.csv', lines.join('\n') + '\n');
  }
}
