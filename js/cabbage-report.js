/* cabbage-report.js — the "Phenotype prediction (CABBAGE)" card that mounts
 * into the analysis results. Two functions:
 *
 *  1. Gene→phenotype prediction: for the sample's species, every detected
 *     resistance gene is looked up in CABBAGE's precomputed association
 *     table (join of ~165k isolates that have both a genome and an
 *     antibiogram): "S. aureus isolates carrying mecA: 97.7% methicillin
 *     resistant (n=1,791)". This turns genotype calls into an empirically
 *     supported phenotype prediction — the bridge the "genotype ≠
 *     phenotype" caveat always needed.
 *
 *  2. Species antibiogram: per-antibiotic resistance rates across all of
 *     CABBAGE's 1.7M AST records for the selected species, optionally
 *     refined by country / year / isolation source (loads the phenotype
 *     snapshot on demand), as context for how typical the predicted profile
 *     is.
 *
 * The association table is ~0.3 MB and loads automatically with the card;
 * nothing about the sample is sent anywhere.
 */

import {
  loadAssociations, loadView, isViewLoaded, getManifest,
  predictPhenotypes, unmatchedGenes, backgroundAntibiogram, antibiogram,
  resolveSpecies, templateGene, fmtInt,
} from './cabbage.js';

const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;')
  .replace(/>/g, '&gt;').replace(/"/g, '&quot;');

const nextFrame = () => new Promise((r) => requestAnimationFrame(() => setTimeout(r, 0)));

// Extract gene symbols from the run's KMA results (per database).
export function collectDetectedGenes(results, parseRes) {
  const out = [];
  const seen = new Set();
  for (const r of results || []) {
    if (!r.resTable) continue;
    const rows = parseRes(r.resTable).rows;
    for (const row of rows) {
      const gene = templateGene(row.Template, r.database);
      if (!gene) continue;
      const key = gene.toLowerCase();
      if (seen.has(key)) continue;
      seen.add(key);
      out.push({ name: gene, db: r.database === 'card_homolog' ? 'CARD' : 'ResFinder' });
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

  const push = (msg, level) => log?.(msg, level);

  renderLoading();
  (async () => {
    try {
      state.assoc = await loadAssociations((got, total) => {
        const pct = total ? Math.round(100 * got / total) : 0;
        setStatus(`Loading CABBAGE phenotype data… ${pct ? pct + '%' : ''}`);
      });
      try { state.release = (await getManifest()).release; } catch (_) {}
      state.species = resolveSpecies(state.assoc.species, organism);
      push(`CABBAGE ${state.release || ''}: ${fmtInt(state.assoc.bySpeciesGene.size)} gene–phenotype associations loaded`, 'ok');
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
        <div class="card-title comp-title"><span>Phenotype prediction (CABBAGE)</span><span class="comp-src">AMR genotype–phenotype database</span></div>
        <div class="card-body"><p class="cabb-status opt-note">Loading CABBAGE phenotype data…</p></div>
      </div>`;
  }

  function renderError(err) {
    container.innerHTML = `
      <div class="card">
        <div class="card-title comp-title"><span>Phenotype prediction (CABBAGE)</span><span class="comp-src">AMR genotype–phenotype database</span></div>
        <div class="card-body">
          <p class="opt-note">Could not load the CABBAGE data (${esc(err.message)}), so phenotype predictions are unavailable for this run. The gene-detection results above are unaffected.</p>
        </div>
      </div>`;
    push('CABBAGE data unavailable: ' + err.message, 'warn');
  }

  function render() {
    const { assoc } = state;
    const options = assoc.species
      .map(s => `<option value="${esc(s)}"${s === state.species ? ' selected' : ''}>${esc(s)}</option>`)
      .join('');

    let head = `
      <div class="card cabb-card" id="cabb-card">
        <div class="card-title comp-title">
          <span>Phenotype prediction (CABBAGE)</span>
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
            <span class="cabb-ab-title">Species antibiogram — resistance rates from all CABBAGE records</span>
            <button class="cabb-refine-btn" id="cabb-refine" type="button">Refine by country / year / source</button>
          </div>
          <div id="cabb-abfilters" hidden></div>
          <div id="cabb-antibiogram"></div>
          <p class="cabb-footnote">Predictions are associations, not susceptibility testing: for each gene they show how often CABBAGE isolates of this species carrying that gene were recorded non-susceptible (updated 2025 CLSI/EUCAST breakpoints first). Always confirm critical results with a laboratory antibiogram. Snapshot of the <a href="https://www.ebi.ac.uk/amr" target="_blank" rel="noopener">EMBL-EBI AMR portal</a>, computed locally in your browser.</p>
        </div>
      </div>`;

    container.innerHTML = head;
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

  // ── Section 1: gene → phenotype predictions ──

  function renderPredictions() {
    const area = container.querySelector('#cabb-predictions');
    const isoEl = container.querySelector('#cabb-isolates');
    if (isoEl) isoEl.textContent = '';
    if (!state.species) {
      area.innerHTML = detectedGenes.length
        ? `<p class="opt-note">Select the sample's species above to interpret the ${detectedGenes.length} detected gene${detectedGenes.length === 1 ? '' : 's'} against CABBAGE's genotype–phenotype pairs.</p>`
        : '';
      return;
    }

    const preds = predictPhenotypes(state.assoc, state.species, detectedGenes);
    const missing = unmatchedGenes(state.assoc, state.species, detectedGenes);

    let html = '';
    if (!detectedGenes.length) {
      html += `<p class="opt-note">No resistance genes were detected in this run, so no phenotype predictions apply. The antibiogram below shows this species' background resistance rates.</p>`;
    } else if (!preds.length) {
      html += `<p class="opt-note">None of the detected genes have phenotype data for <b>${esc(state.species)}</b> in CABBAGE (each association needs ≥20 isolates).</p>`;
    } else {
      const PRED_PREVIEW = 12;
      html += '<div class="cabb-preds">';
      for (const p of preds) {
        const geneChips = p.genes.map(g => `
          <span class="cabb-gene" title="${esc(g.db)} hit">${esc(g.gene)}
            <b>${(100 * (g.r + g.i) / g.n).toFixed(1)}%</b> non-susceptible (n=${fmtInt(g.n)})
          </span>`).join('');
        const bgTxt = p.background
          ? `background ${(100 * (p.background.r + p.background.i) / p.background.n).toFixed(1)}% (n=${fmtInt(p.background.n)})`
          : 'no background data';
        html += `
          <div class="cabb-pred cabb-verdict-${p.verdict}">
            <span class="cabb-verdict-pill" title="strongest gene association for this antibiotic">${p.verdict === 'resistant' ? 'predicted resistant' : p.verdict === 'uncertain' ? 'uncertain' : 'not resistant'}</span>
            <span class="cabb-ab-name" title="${esc(p.antibiotic)}">${esc(p.antibiotic)}</span>
            <div class="cabb-bars">
              <div class="cabb-bar-row" title="isolates of this species carrying the detected gene(s)">
                <span class="cabb-bar-track"><span class="cabb-bar-fill gene" style="width:${Math.round(100 * p.rate)}%"></span></span>
                <span class="cabb-bar-val">${(100 * p.rate).toFixed(1)}%</span>
              </div>
              ${p.background ? `
              <div class="cabb-bar-row bg" title="all isolates of this species in CABBAGE">
                <span class="cabb-bar-track"><span class="cabb-bar-fill bg" style="width:${Math.round(100 * (p.background.r + p.background.i) / p.background.n)}%"></span></span>
                <span class="cabb-bar-val">${(100 * (p.background.r + p.background.i) / p.background.n).toFixed(1)}%</span>
              </div>` : ''}
            </div>
            <span class="cabb-genes">${geneChips}</span>
            <span class="cabb-bg-note">${esc(bgTxt)}</span>
          </div>`;
      }
      html += '</div>';
      if (preds.length > PRED_PREVIEW) {
        html += `<button class="cabb-preds-more" type="button">show all ${preds.length} antibiotics</button>`;
      }
    }

    if (missing.length) {
      html += `<p class="cabb-missing">No CABBAGE phenotype link for: ${missing.map(g => `<code>${esc(g.name)}</code>`).join(', ')} — detected by ${missing[0].db === 'CARD' ? 'CARD' : 'ResFinder'} but no isolates in the database connect ${esc(state.species)} carriers of these genes to an AST result.</p>`;
    }
    area.innerHTML = html;
    // Cap the visible rows at 12 until the user expands (same pattern as the
    // gene tables).
    const extra = [...area.querySelectorAll('.cabb-pred')].slice(12);
    extra.forEach(el => el.classList.add('row-limited'));
    area.querySelector('.cabb-preds-more')?.addEventListener('click', (e) => {
      extra.forEach(el => el.classList.remove('row-limited'));
      e.target.remove();
    });
  }

  // ── Section 2: species antibiogram (background or filtered) ──

  function renderAntibiogram() {
    const area = container.querySelector('#cabb-antibiogram');
    if (!state.species) { area.innerHTML = ''; return; }
    const f = state.filters;
    const filtered = state.phenotypeView && (f.country || f.sourceCategory || f.yearFrom || f.yearTo);

    let rows;
    if (filtered) {
      rows = antibiogram(state.phenotypeView, state.species, {
        country: f.country || null,
        sourceCategory: f.sourceCategory || null,
        yearFrom: f.yearFrom ? Number(f.yearFrom) : null,
        yearTo: f.yearTo ? Number(f.yearTo) : null,
      });
    } else {
      rows = backgroundAntibiogram(state.assoc, state.species);
    }

    const isoEl = container.querySelector('#cabb-isolates');
    if (isoEl) {
      const total = rows.reduce((n, r) => Math.max(n, r.n), 0);
      isoEl.textContent = total ? `up to ${fmtInt(total)} records per antibiotic` : '';
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
          <span class="cabb-bar-n">n=${fmtInt(r.n)}</span>
        </div>`;
    }
    html += '</div>';
    if (rows.length > 25) html += `<p class="opt-note">+ ${rows.length - 25} more antibiotics — all of them are included in the results ZIP.</p>`;
    area.innerHTML = html;
  }

  // ── Filters (require the full phenotype view) ──

  async function enableFilters() {
    const box = container.querySelector('#cabb-abfilters');
    box.hidden = false;
    const btn = container.querySelector('#cabb-refine');
    if (state.phenotypeView) { buildFilterControls(box); return; }
    box.innerHTML = '<p class="opt-note">Loading the full CABBAGE phenotype snapshot (5.7 MB, one time)…</p>';
    btn.disabled = true;
    try {
      if (!isViewLoaded('phenotypes')) {
        state.phenotypeView = await loadView('phenotypes', (got, total) => {
          const pct = total ? Math.round(100 * got / total) : 0;
          box.innerHTML = `<p class="opt-note">Downloading CABBAGE phenotypes… ${pct}%</p>`;
        });
        push('CABBAGE phenotype view loaded — antibiogram filters enabled', 'ok');
      } else {
        state.phenotypeView = await loadView('phenotypes');
      }
      buildFilterControls(box);
    } catch (err) {
      box.innerHTML = `<p class="opt-note">Could not load the phenotype snapshot (${esc(err.message)}); the unfiltered antibiogram is shown.</p>`;
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
      return `<select data-filter="${label}"><option value="">all ${label === 'country' ? 'countries' : label === 'sourceCategory' ? 'source categories' : ''}</option>${values.map(v => `<option value="${esc(v)}"${state.filters[label] === v ? ' selected' : ''}>${esc(v)}</option>`).join('')}</select>`;
    };
    let years = [];
    const yr = view.colById.get('phenotype-collection_year');
    if (yr) {
      const set = new Set();
      for (let i = 0; i < yr.vals.length; i++) if (!Number.isNaN(yr.vals[i])) set.add(yr.vals[i]);
      years = [...set].sort((a, b) => a - b);
    }
    const yearOpts = (sel) => `<option value="">${sel}</option>` + years.map(y => `<option value="${y}"${String(state.filters.yearFrom) === String(y) || String(state.filters.yearTo) === String(y) ? ' selected' : ''}>${y}</option>`).join('');

    box.innerHTML = `
      <div class="cabb-filters">
        ${opts('phenotype-country', 'country')}
        ${opts('phenotype-isolation_source_category', 'sourceCategory')}
        <select data-filter="yearFrom">${yearOpts('year from')}</select>
        <select data-filter="yearTo">${yearOpts('year to')}</select>
        <button type="button" class="cabb-clear-filters" id="cabb-clear-filters">clear</button>
      </div>`;
    box.querySelectorAll('select[data-filter]').forEach(sel => {
      sel.addEventListener('change', () => {
        state.filters[sel.dataset.filter] = sel.value;
        box.querySelectorAll('select[data-filter]').forEach(s => {
          if (s !== sel) s.value = state.filters[s.dataset.filter] || '';
        });
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
    const preds = predictPhenotypes(state.assoc, state.species, detectedGenes);
    const lines = ['species,gene,gene_database,antibiotic,resistant,intermediate,susceptible,n,non_susceptible_rate'];
    for (const p of preds) {
      for (const g of p.genes) {
        lines.push([state.species, g.gene, g.db, p.antibiotic, g.r, g.i, g.s, g.n, ((g.r + g.i) / g.n).toFixed(4)].join(','));
      }
    }
    lines.push('');
    lines.push('# species antibiogram (all CABBAGE AST records' + (state.filters.country || state.filters.sourceCategory || state.filters.yearFrom || state.filters.yearTo ? ', filtered' : '') + ')');
    lines.push('species,antibiotic,resistant,intermediate,susceptible,n,non_susceptible_rate');
    const f = state.filters;
    const filtered = state.phenotypeView && (f.country || f.sourceCategory || f.yearFrom || f.yearTo);
    const rows = filtered
      ? antibiogram(state.phenotypeView, state.species, {
        country: f.country || null, sourceCategory: f.sourceCategory || null,
        yearFrom: f.yearFrom ? Number(f.yearFrom) : null, yearTo: f.yearTo ? Number(f.yearTo) : null,
      })
      : backgroundAntibiogram(state.assoc, state.species);
    for (const r of rows) {
      lines.push([state.species, JSON.stringify(r.antibiotic), r.r, r.i, r.s, r.n, ((r.r + r.i) / r.n).toFixed(4)].join(','));
    }
    onCsv('cabbage_predictions.csv', lines.join('\n') + '\n');
  }
}
