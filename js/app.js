import { initWasm, runAnalysis, setLogCallback, formatBytes } from './wasm-runtime.js';
import { loadPhenotypesDB, parseResFile } from './resfinder-db.js';
import { cacheDBFile, getCachedDBFile } from './db.js';
import { openGeneViewer } from './gene-view.js';
import { fetchAssetWithProgress } from './assets.js';
import {
  runComprehensive, summariseQc, qcVerdict, COMPREHENSIVE_DBS,
} from './comprehensive.js';

const state = {
  readType: 'paired',
  mode: localStorage.getItem('op-mode') || 'simple',
  files: { r1: null, r2: null },
  wasmReady: false,
  running: false,
};

const term = {
  el: null,
  lines: [],
  push(msg, level) {
    if (!this.el) this.el = document.getElementById('terminal');
    const line = document.createElement('div');
    line.className = 'term-line term-' + (level || 'info');
    line.textContent = msg;
    this.el.appendChild(line);
    this.el.scrollTop = this.el.scrollHeight;
    this.lines.push(msg);
  },
  clear() {
    if (!this.el) this.el = document.getElementById('terminal');
    this.el.innerHTML = '';
    this.lines = [];
  },
  // Plain-text transcript of everything logged this run — bundled with the
  // downloads so the exact command and outcome are recoverable.
  text() {
    return this.lines.join('\n') + '\n';
  },
};

document.addEventListener('DOMContentLoaded', async () => {
  setupReadTypeToggle();
  setupUploadSlots();
  setupSliders();
  setupModeToggle();
  document.getElementById('btn-run').addEventListener('click', doRun);
  document.getElementById('btn-example').addEventListener('click', () => loadExample('sub'));
  document.getElementById('btn-example-full').addEventListener('click', () => loadExample('full'));
  document.getElementById('btn-back').addEventListener('click', () => {
    document.getElementById('workspace').classList.remove('show-results');
  });
  document.getElementById('db-select').addEventListener('change', applyMode);
  setupLogsToggle();
  applyReadType();

  const logToTerm = (msg, level) => term.push(msg, level);
  setLogCallback(logToTerm);

  try {
    await Promise.all([initWasm(), loadPhenotypesDB()]);
    state.wasmReady = true;
    setStatus('Ready', 'ok');
    updateRunButton();
  } catch (err) {
    setStatus('Error', 'error');
    term.push('Init failed: ' + err.message, 'error');
  }
});

function setStatus(text, kind) {
  const el = document.getElementById('wasm-status');
  if (!el) return;
  el.textContent = text;
  el.className = 'badge badge-' + kind;
}

// ── Logs panel ──

function setupLogsToggle() {
  const card = document.getElementById('terminal-card');
  const toggle = document.getElementById('terminal-toggle');
  toggle.addEventListener('click', () => {
    const collapsed = card.classList.toggle('collapsed');
    toggle.setAttribute('aria-expanded', String(!collapsed));
  });
}

// Errors land in the Logs panel; surface it instead of failing silently
// when the panel happens to be collapsed.
function expandLogs() {
  document.getElementById('terminal-card').classList.remove('collapsed');
  document.getElementById('terminal-toggle').setAttribute('aria-expanded', 'true');
}

// ── Mode (Simple / Advanced) ──
//
// Simple (default): one click → fastp QC → ResFinder + CARD + VFDB →
// plain-language synthesis. Advanced: the original single-database UI.

function setupModeToggle() {
  document.querySelectorAll('#mode-toggle .seg-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      state.mode = btn.dataset.mode;
      localStorage.setItem('op-mode', state.mode);
      document.querySelectorAll('#mode-toggle .seg-btn').forEach(b =>
        b.classList.toggle('active', b === btn));
      applyMode();
    });
  });
  document.querySelectorAll('#mode-toggle .seg-btn').forEach(b =>
    b.classList.toggle('active', b.dataset.mode === state.mode));
  applyMode();
}

function applyMode() {
  const simple = state.mode === 'simple';
  const kmaOpts = document.getElementById('kma-options');
  if (kmaOpts) kmaOpts.style.display = simple ? 'none' : '';
  // the hint explains both modes, so it stays visible in either
  if (!state.running) {
    document.getElementById('btn-run').textContent = simple ? 'Analyze sample' : 'Run analysis';
  }
}

// ── Read type ──

function setupReadTypeToggle() {
  document.querySelectorAll('#read-type .seg-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      document.querySelectorAll('#read-type .seg-btn').forEach(b => b.classList.toggle('active', b === btn));
      state.readType = btn.dataset.type;
      applyReadType();
    });
  });
}

function applyReadType() {
  const slotR2 = document.querySelector('[data-slot="r2"]');
  const slotR1Title = document.querySelector('[data-slot="r1"] .slot-title');
  if (state.readType === 'paired') {
    slotR2.style.display = '';
    slotR1Title.textContent = 'Forward reads (R1)';
  } else {
    slotR2.style.display = 'none';
    state.files.r2 = null;
    setSlotName('r2', 'No file');
    slotR1Title.textContent = state.readType === 'nanopore' ? 'Nanopore reads' : 'Reads';
  }
  updateRunButton();
}

// ── Uploads ──

function setupUploadSlots() {
  document.querySelectorAll('[data-pick]').forEach(btn => {
    const slot = btn.dataset.pick;
    const input = document.querySelector(`[data-input="${slot}"]`);
    btn.addEventListener('click', () => input.click());
    input.addEventListener('change', () => {
      const f = input.files[0];
      if (f) {
        state.files[slot] = f;
        setSlotName(slot, `${f.name} (${formatBytes(f.size)})`);
      } else {
        state.files[slot] = null;
        setSlotName(slot, 'No file');
      }
      input.value = '';
      updateRunButton();
    });
  });
}

function setSlotName(slot, text) {
  const el = document.querySelector(`[data-name="${slot}"]`);
  if (el) el.textContent = text;
}

// ── Sliders ──

function setupSliders() {
  const id = document.getElementById('slider-id');
  const cov = document.getElementById('slider-cov');
  id.addEventListener('input', () => { document.getElementById('slider-id-val').textContent = id.value; });
  cov.addEventListener('input', () => { document.getElementById('slider-cov-val').textContent = cov.value; });
}

function updateRunButton() {
  const btn = document.getElementById('btn-run');
  const haveR1 = !!state.files.r1;
  const needR2 = state.readType === 'paired';
  const haveR2 = !!state.files.r2;
  const ready = state.wasmReady && !state.running && haveR1 && (!needR2 || haveR2);
  btn.disabled = !ready;
}

// ── Example data ──
//
// S. aureus JKD6159 — Australian CA-MRSA ST93-IV (Chua et al. 2010/2011);
// Illumina PE resequencing SRR21386014 (Wick et al. 2023).
//   sub  — 300 000 read pairs (~51 MB) shipped with the app for a fast demo.
//   full — the complete public run (~570 MB), fetched from ENA on demand.
const EXAMPLE = {
  accession: 'SRR21386014',
  organism: 'S. aureus JKD6159',
  sub: {
    note: '300 000 read-pair subsample of SRR21386014 (~51 MB), shipped with the app',
    files: [
      { slot: 'r1', name: 'SRR21386014_sub_1.fastq.gz', bytes: 26214400, url: 'SRR21386014_sub_1.fastq.gz' },
      { slot: 'r2', name: 'SRR21386014_sub_2.fastq.gz', bytes: 27262976, url: 'SRR21386014_sub_2.fastq.gz' },
    ],
  },
  full: {
    note: 'full SRR21386014 run from ENA (~3.4 million PE reads, ~570 MB gzipped)',
    files: [
      {
        slot: 'r1', name: 'SRR21386014_1.fastq.gz', bytes: 292999499,
        url: 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR213/014/SRR21386014/SRR21386014_1.fastq.gz',
      },
      {
        slot: 'r2', name: 'SRR21386014_2.fastq.gz', bytes: 303106337,
        url: 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR213/014/SRR21386014/SRR21386014_2.fastq.gz',
      },
    ],
  },
};

async function loadExample(kind) {
  const specSet = EXAMPLE[kind];
  if (!specSet) return;
  const btn = document.getElementById(kind === 'full' ? 'btn-example-full' : 'btn-example');
  const other = document.getElementById(kind === 'full' ? 'btn-example' : 'btn-example-full');
  if (btn.disabled) return;
  btn.disabled = true;
  other.disabled = true;
  const original = btn.textContent;
  const fromEna = kind === 'full';
  try {
    document.querySelectorAll('#read-type .seg-btn').forEach(b => {
      b.classList.toggle('active', b.dataset.type === 'paired');
    });
    state.readType = 'paired';
    applyReadType();
    term.push(`Example isolate: ${EXAMPLE.organism} · ENA/SRA ${EXAMPLE.accession}`, 'info');
    term.push(specSet.note, 'info');
    term.push('Wick et al. Microbiol Resour Announc. 2023;12:e01129-22.  doi:10.1128/mra.01129-22', 'info');

    for (const spec of specSet.files) {
      const cacheKey = 'example:' + spec.name;
      let data = fromEna ? await getCachedDBFile(cacheKey) : null;
      if (data) {
        term.push(`${spec.name} (cached, ${formatBytes(data.byteLength)})`, 'info');
        btn.textContent = `Cached ${spec.name}`;
      } else {
        term.push(fromEna ? `Downloading ${spec.name} from ENA…` : `Loading ${spec.name}…`, 'progress');
        data = await fetchAssetWithProgress(spec.url, spec.bytes, (got, total) => {
          const pct = total ? Math.min(100, Math.round(100 * got / total)) : 0;
          btn.textContent = `${spec.name}  ${formatBytes(got)}${total ? ' / ' + formatBytes(total) : ''}  (${pct}%)`;
        });
        term.push(`${spec.name} (${formatBytes(data.byteLength)})`, 'ok');
        if (fromEna) {
          try { await cacheDBFile(cacheKey, data); } catch (_) {
            term.push('Could not cache the FASTQ in IndexedDB; next load will re-download.', 'warn');
          }
        }
      }
      const file = new File([data], spec.name, { type: 'application/gzip' });
      state.files[spec.slot] = file;
      setSlotName(spec.slot, `${file.name} (${formatBytes(file.size)})`);
    }
    updateRunButton();
    term.push('Ready. Choose ResFinder (mecA, blaZ) or VFDB (PVL, hla, ica…) and run.', 'ok');
  } catch (err) {
    term.push('Failed to load example: ' + err.message, 'error');
    expandLogs();
    if (fromEna) term.push('ENA must be reachable from this browser (CORS). Try the ~50 MB subsample, or upload your own FASTQ.', 'info');
  } finally {
    btn.textContent = original;
    btn.disabled = false;
    other.disabled = false;
  }
}

// ── Run ──

async function doRun() {
  if (state.running || !state.wasmReady || !state.files.r1) return;
  state.running = true;
  updateRunButton();
  setRunning(true);
  term.clear();
  document.getElementById('results').innerHTML = '';

  const files = [state.files.r1];
  if (state.readType === 'paired' && state.files.r2) files.push(state.files.r2);

  try {
    if (state.mode === 'simple') {
      await doRunSimple(files);
    } else {
      const dbKey = document.getElementById('db-select').value;
      if (dbKey === 'hv_profile') {
        // research tool, not exposed in the UI — kept for future work
        await runHvProfile(files);
      } else {
        const config = {
          id_threshold: parseInt(document.getElementById('slider-id').value) / 100,
          mrc: parseInt(document.getElementById('slider-cov').value) / 100,
          nanopore: state.readType === 'nanopore',
        };
        const result = await runAnalysis(files, dbKey, config);
        result.files = result.files || {};
        result.files['.log'] = { data: term.text(), binary: false };
        renderResults(result);
      }
    }
    document.getElementById('workspace').classList.add('show-results');
  } catch (err) {
    term.push('Error: ' + err.message, 'error');
  } finally {
    state.running = false;
    setRunning(false);
    updateRunButton();
  }
}

// ── Simple mode: the comprehensive pipeline ──

async function doRunSimple(files) {
  const statusText = document.querySelector('#run-status .run-status-text');
  const onStep = (i, n, text) => {
    if (statusText) statusText.textContent = `Step ${i} of ${n}: ${text}`;
    term.push(`[${i}/${n}] ${text}`, 'progress');
  };
  const { setFastpLog } = await import('./comprehensive.js');
  setFastpLog((t) => term.push(t, 'info'));
  const result = await runComprehensive(files, state.readType, { onStep });
  if (statusText) statusText.textContent = 'Analyzing… this may take a moment.';
  renderComprehensive(result);
}

function renderComprehensive({ qc, results }) {
  const area = document.getElementById('results');
  const metrics = summariseQc(qc);
  const verdict = qcVerdict(metrics);

  const anyOk = results.some(r => r.exitCode === 0);
  let html = `
    <div class="result-header ${anyOk ? 'ok' : 'error'}">
      <strong>${anyOk ? 'Analysis complete' : 'Analysis failed'}</strong>
      <span class="result-meta">Simple mode · QC + ResFinder + CARD + VFDB</span>
    </div>`;

  // ── Downloads first: one zip with everything ──
  html += `
    <div class="downloads">
      <button class="btn btn-primary" data-zip>Download all results (ZIP)</button>
      <p class="dl-nudge">The ZIP has everything from this run: the full result tables, the quality report and the run log.</p>
    </div>`;

  // ── QC card ──
  html += `
    <section class="card">
      <div class="card-title comp-title"><span>Quality control</span><span class="comp-src">fastp</span></div>
      <div class="card-body">
        ${metrics ? `
        <div class="qc-grid">
          <div class="qc-cell"><span class="qc-val">${fmtInt(metrics.rawReads)}</span><span class="qc-key">reads in</span></div>
          <div class="qc-cell"><span class="qc-val">${(100 * metrics.retained).toFixed(1)}%</span><span class="qc-key">retained after trimming</span></div>
          <div class="qc-cell"><span class="qc-val">${(100 * (metrics.q30After ?? 0)).toFixed(1)}%</span><span class="qc-key">Q30 (after)</span></div>
          <div class="qc-cell"><span class="qc-val">${(100 * (metrics.gcBefore ?? 0)).toFixed(1)}%</span><span class="qc-key">GC content</span></div>
          ${metrics.duplication != null ? `<div class="qc-cell"><span class="qc-val">${(100 * metrics.duplication).toFixed(1)}%</span><span class="qc-key">duplication</span></div>` : ''}
        </div>
        <p class="qc-verdict qc-${verdict.tone}">${esc(verdict.text)}</p>` : `
        <p class="qc-verdict qc-warn">${esc(verdict.text)}</p>`}
      </div>
    </section>`;

  // ── One section per database, reusing the advanced-mode table ──
  results.forEach((r, i) => {
    const db = COMPREHENSIVE_DBS.find(d => d.key === r.database) || { label: r.dbLabel };
    const ok = r.exitCode === 0;
    const rows = ok && r.resTable ? parseResFile(r.resTable).rows : [];
    const countMeta = rows.length ? ` · ${rows.length} gene${rows.length === 1 ? '' : 's'} detected` : '';
    html += `
    <div class="comp-section" data-section="${i}">
      <div class="result-header ${ok ? 'ok' : 'error'}">
        <strong>${esc(db.label)}</strong>
        <span class="result-meta">${esc(r.dbLabel || '')}${countMeta}</span>
      </div>
      ${ok && r.resTable ? renderResTable(r.resTable)
        : ok ? '<p class="empty">No resistance or virulence genes detected.</p>'
        : '<p class="empty">This step failed. Open the Logs panel below for details.</p>'}
    </div>`;
  });

  // ── Downloads: everything in one zip ──
  const dlFiles = {};
  if (qc) {
    dlFiles['qc_report.json'] = { data: new TextEncoder().encode(JSON.stringify(qc, null, 2)), binary: true };
  }
  for (const r of results) {
    for (const [ext, info] of Object.entries(r.files || {})) {
      if (r.files['.log'] && ext === '.log') continue;
      dlFiles[`${r.database}${ext}`] = info;
    }
  }
  dlFiles['.log'] = { data: term.text(), binary: false };

  area.innerHTML = html;

  // The advanced-mode table bindings, scoped to each database section.
  results.forEach((r, i) => {
    const section = area.querySelector(`[data-section="${i}"]`);
    if (!section) return;
    bindSort(section);
    bindColumnResize(section);
    bindGeneViewer(section, r);
  });

  // Keep each table to its top 5 rows until the user asks for more. Sorting
  // implies interest in the full table, so a header click reveals everything.
  const TABLE_PREVIEW = 5;
  area.querySelectorAll('.comp-section').forEach(section => {
    const trs = [...section.querySelectorAll('tbody tr')];
    if (trs.length <= TABLE_PREVIEW) return;
    const extra = trs.slice(TABLE_PREVIEW);
    extra.forEach(tr => tr.classList.add('row-limited'));
    const cols = section.querySelectorAll('th').length;
    const moreRow = document.createElement('tr');
    moreRow.innerHTML = `<td colspan="${cols}"><button class="show-more-btn" type="button">Show all ${trs.length} genes</button></td>`;
    section.querySelector('tbody').appendChild(moreRow);
    const reveal = () => {
      extra.forEach(tr => tr.classList.remove('row-limited'));
      moreRow.remove();
    };
    moreRow.querySelector('button').addEventListener('click', reveal);
    section.querySelectorAll('th').forEach(th => th.addEventListener('click', reveal));
  });

  area.querySelector('[data-zip]')?.addEventListener('click', () => {
    const enc = new TextEncoder();
    const zip = {};
    for (const [fname, info] of Object.entries(dlFiles)) {
      zip[fname] = info.binary
        ? (info.data instanceof Uint8Array ? info.data : new Uint8Array(info.data))
        : enc.encode(info.data);
    }
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([window.fflate.zipSync(zip)], { type: 'application/zip' }));
    a.download = 'openpathogen_results.zip';
    a.click();
    URL.revokeObjectURL(a.href);
  });
}

function fmtInt(n) {
  if (n == null) return '—';
  if (n >= 1e6) return (n / 1e6).toFixed(1) + 'M';
  if (n >= 1e3) return (n / 1e3).toFixed(0) + 'k';
  return String(n);
}

// ── HV species profile ──
//
// Runs the hypervector profiler worker: stream FASTQ → 32 kB fingerprint →
// deconvolution against the int8 reference matrix. Same file slots and
// example data as KMA; reads never leave the worker.

function runHvProfile(files) {
  return new Promise((resolve, reject) => {
    const workerUrl = new URL('./hv-worker.js?v=2', import.meta.url);
    const worker = new Worker(workerUrl, { type: 'module' });
    let lastProgressAt = 0;
    term.push('HV profiler: encoding reads into a 32 kB hypervector, then', 'info');
    term.push('deconvolving against the pathogen reference matrix (species-level, β).', 'info');
    worker.onmessage = (e) => {
      const m = e.data;
      if (m.type === 'log') {
        term.push(m.text, 'info');
      } else if (m.type === 'progress') {
        const now = performance.now();
        if (now - lastProgressAt > 400) {
          lastProgressAt = now;
          term.push(
            `encoded ${(m.bytes / 1e6).toFixed(0)} MB · ${Math.round(m.reads / 1e3)}k reads · ` +
            `${Math.round(m.events / 1e3)}k events · ${m.mbps.toFixed(0)} MB/s`,
            'progress'
          );
        }
      } else if (m.type === 'done') {
        worker.terminate();
        renderHvResults(m);
        resolve(m);
      } else if (m.type === 'error') {
        worker.terminate();
        reject(new Error(m.message));
      }
    };
    worker.onerror = (err) => {
      worker.terminate();
      reject(new Error(err.message || 'worker failed'));
    };
    worker.postMessage({
      type: 'run',
      files,
      dbUrl: new URL('../hvprof/hv_patho_v1.hvd', location.href).href,
      tsvUrl: new URL('../hvprof/hv_patho_v1.tsv', location.href).href,
    });
  });
}

function renderHvResults(result) {
  const area = document.getElementById('results');
  const h = result.report;
  const taxa = result.taxa || [];

  const rows = taxa.map((t) => {
    const m = t.meta || { name: `row ${t.idx}`, genus: '' };
    const flags = [];
    if (t.flags & 256) flags.push('species complex');
    if (t.flags & 2) flags.push('host');
    return `
      <tr>
        <td title="${esc(m.id || '')}">${esc(m.name)}</td>
        <td>${esc(m.genus || '')}</td>
        <td class="num">${t.coverage.toFixed(2)}×</td>
        <td class="num">${(100 * t.dnaFrac).toFixed(2)}%</td>
        <td class="num">${t.z.toFixed(1)}</td>
        <td class="num">${fmtQ(t.q)}</td>
        <td>${flags.map(esc).join(', ')}</td>
      </tr>`;
  }).join('');

  area.innerHTML = `
    <div class="result-header ok">
      <strong>Species profile complete</strong>
      <span class="result-meta">HV profiler · ${result.secs.toFixed(1)}s · ${(result.totalBytes / 1e6).toFixed(0)} MB read</span>
    </div>
    <div class="downloads">
      <button class="btn btn-primary" data-hv-fp>Download fingerprint (.hvf, ${(result.fingerprint.byteLength / 1024).toFixed(0)} kB)</button>
    </div>
    <p class="hv-summary">
      ${(100 * (1 - h.unexplained)).toFixed(1)}% of the sample's k-mer mass is explained by known
      species; <strong>${(100 * h.unexplained).toFixed(1)}% is uncharacterised</strong>
      (novel organisms or low-identity relatives).
      QC spike-in z = ${h.spikeZ.toFixed(2)} (null ✓ when |z| &lt; 3).
    </p>
    ${taxa.length ? `
    <div class="table-wrap"><table class="res-table">
      <thead><tr>
        <th>Organism</th><th>Genus</th><th>Coverage</th><th>DNA fraction</th><th>z</th><th>q</th><th>Notes</th>
      </tr></thead>
      <tbody>${rows}</tbody>
    </table></div>` : '<p class="empty">No reference species detected above threshold.</p>'}`;

  area.querySelector('[data-hv-fp]')?.addEventListener('click', () => {
    const blob = new Blob([result.fingerprint], { type: 'application/octet-stream' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'sample.hvf';
    a.click();
    URL.revokeObjectURL(a.href);
  });
}

function fmtQ(q) {
  if (q < 1e-200) return '<10⁻²⁰⁰';
  if (q < 1e-3) return q.toExponential(1);
  return q.toPrecision(2);
}

// Show a running indicator (spinner in the run button + indeterminate progress
// bar) so it's clear the job is working. The analysis itself runs in a Web
// Worker, so these animations keep ticking while KMA crunches.
function setRunning(on) {
  const btn = document.getElementById('btn-run');
  const status = document.getElementById('run-status');
  document.getElementById('btn-example').disabled = on;
  document.getElementById('btn-example-full').disabled = on;
  if (on) {
    btn.innerHTML = '<span class="spinner"></span>Running…';
    status.hidden = false;
  } else {
    btn.textContent = state.mode === 'simple' ? 'Analyze sample' : 'Run analysis';
    status.hidden = true;
  }
}

// ── Results ──

function renderResults(result) {
  const area = document.getElementById('results');
  const ok = result.exitCode === 0;

  let html = `
    <div class="result-header ${ok ? 'ok' : 'error'}">
      <strong>${ok ? 'Analysis complete' : 'Analysis failed'}</strong>
      <span class="result-meta">${esc(result.sampleName)} · ${esc(result.dbLabel)} · ${result.elapsed}s</span>
    </div>`;

  html += renderDownloads(result);

  if (ok && result.resTable) html += renderResTable(result.resTable);
  else if (ok) html += '<p class="empty">No resistance or virulence genes detected.</p>';
  area.innerHTML = html;
  bindDownloads(area, result);
  bindSort(area);
  bindColumnResize(area);
  bindGeneViewer(area, result);
}

// Open the per-gene detail viewer when a result row is clicked.
function bindGeneViewer(area, result) {
  if (!result.resTable) return;
  const { rows } = parseResFile(result.resTable);
  const byTemplate = new Map(rows.map(r => [r.Template, r]));
  area.querySelectorAll('tr[data-template]').forEach(tr => {
    tr.addEventListener('click', () => {
      // Don't hijack a text selection drag.
      if (window.getSelection && String(window.getSelection()).length) return;
      const row = byTemplate.get(tr.dataset.template);
      if (row) openGeneViewer(row, result);
    });
  });
}

function renderResTable(tsv) {
  const { headers, rows } = parseResFile(tsv);
  if (!rows.length) return '<p class="empty">No resistance or virulence genes detected.</p>';

  const cols = ['Template', 'Template_Identity', 'Template_Coverage', 'Depth', 'p_value']
    .filter(c => headers.includes(c));

  const colWidths = { Template: 170, Template_Identity: 90, Template_Coverage: 90, Depth: 65, p_value: 75 };

  let html = '<div class="table-wrap"><table class="res-table"><colgroup>';
  for (const c of cols) html += `<col style="width:${colWidths[c]}px">`;
  html += '</colgroup><thead><tr>';
  for (const c of cols) html += `<th data-key="${esc(c)}">${esc(c.replace(/_/g, ' '))}<span class="arr">↕</span></th>`;
  html += '</tr></thead><tbody>';

  for (const row of rows) {
    html += `<tr data-template="${esc(row.Template || '')}">`;
    for (const c of cols) {
      let v = row[c] || '';
      if (/Identity|Coverage/.test(c) && parseFloat(v)) v = parseFloat(v).toFixed(1) + '%';
      const inner = c === 'Template' ? `<span class="gv-link">${esc(v)}</span>` : esc(v);
      html += `<td data-k="${esc(c)}" data-v="${esc(row[c] || '')}" title="${esc(row[c] || '')}">${inner}</td>`;
    }
    html += '</tr>';
  }
  return html + '</tbody></table></div>';
}

function bindSort(area) {
  const table = area.querySelector('.res-table');
  if (!table) return;
  const tbody = table.querySelector('tbody');
  let sortKey = null, sortDir = 1;
  table.querySelectorAll('th').forEach(th => {
    th.addEventListener('click', () => {
      const key = th.dataset.key;
      sortDir = (sortKey === key) ? -sortDir : 1;
      sortKey = key;
      table.querySelectorAll('th').forEach(h => {
        h.classList.toggle('sort', h === th);
        const arr = h.querySelector('.arr');
        if (arr) arr.textContent = (h === th) ? (sortDir > 0 ? '↑' : '↓') : '↕';
      });
      const rows = [...tbody.querySelectorAll('tr')];
      rows.sort((a, b) => {
        const av = a.querySelector(`td[data-k="${key}"]`)?.dataset.v || '';
        const bv = b.querySelector(`td[data-k="${key}"]`)?.dataset.v || '';
        const an = parseFloat(av), bn = parseFloat(bv);
        const numeric = !isNaN(an) && !isNaN(bn) && av.trim() !== '' && bv.trim() !== '';
        const cmp = numeric ? (an - bn) : av.localeCompare(bv);
        return cmp * sortDir;
      });
      rows.forEach(r => tbody.appendChild(r));
    });
  });
}

function bindColumnResize(area) {
  const table = area.querySelector('.res-table');
  if (!table) return;
  const colgroup = table.querySelector('colgroup');
  if (!colgroup) return;
  const cols = [...colgroup.querySelectorAll('col')];

  table.querySelectorAll('th').forEach((th, i) => {
    if (i >= cols.length) return;
    const handle = document.createElement('div');
    handle.className = 'resize-handle';
    th.appendChild(handle);

    let startX, startWidth;
    handle.addEventListener('mousedown', (e) => {
      e.preventDefault();
      e.stopPropagation();
      startX = e.pageX;
      startWidth = th.offsetWidth;
      handle.classList.add('active');

      const onMove = (e) => {
        const diff = e.pageX - startX;
        const newWidth = Math.max(40, startWidth + diff);
        cols[i].style.width = newWidth + 'px';
      };
      const onUp = () => {
        handle.classList.remove('active');
        document.removeEventListener('mousemove', onMove);
        document.removeEventListener('mouseup', onUp);
      };
      document.addEventListener('mousemove', onMove);
      document.addEventListener('mouseup', onUp);
    });
  });
}

function renderDownloads(result) {
  const entries = Object.entries(result.files);
  if (!entries.length) return '';
  let html = '<div class="downloads"><button class="btn btn-primary" data-zip>Download all (ZIP)</button>';
  for (const [ext] of entries) html += `<button class="btn btn-secondary" data-ext="${esc(ext)}">${esc(ext)}</button>`;
  return html + '</div>';
}

function bindDownloads(area, result) {
  area.querySelector('[data-zip]')?.addEventListener('click', () => {
    const enc = new TextEncoder();
    const zip = {};
    for (const [ext, info] of Object.entries(result.files)) {
      const fname = result.sampleName + '_' + result.database + ext;
      zip[fname] = info.binary
        ? (info.data instanceof Uint8Array ? info.data : new Uint8Array(info.data))
        : enc.encode(info.data);
    }
    download(new Blob([window.fflate.zipSync(zip)], { type: 'application/zip' }), result.sampleName + '_results.zip');
  });
  area.querySelectorAll('[data-ext]').forEach(btn => {
    btn.addEventListener('click', () => {
      const ext = btn.dataset.ext;
      const info = result.files[ext];
      if (!info) return;
      const fname = result.sampleName + '_' + result.database + ext;
      const d = info.binary ? (info.data instanceof Uint8Array ? info.data : new Uint8Array(info.data)) : info.data;
      download(new Blob([d], { type: info.binary ? 'application/octet-stream' : 'text/plain' }), fname);
    });
  });
}

function download(blob, name) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = name;
  a.click();
  URL.revokeObjectURL(a.href);
}

function esc(s) {
  return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
}
