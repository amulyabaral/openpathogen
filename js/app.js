import { initWasm, setLogCallback, formatBytes, isDatabaseCached, preloadDatabase } from './wasm-runtime.js';
import { loadPhenotypesDB, parseResFile } from './resfinder-db.js';
import { cacheDBFile, getCachedDBFile, clearDBCache } from './db.js';
import { openGeneViewer } from './gene-view.js';
import { fetchAssetWithProgress } from './assets.js';
import { lookupRun } from './fetch-run.js';
import { collectDetectedGenes, mountCabbageReport } from './cabbage-report.js';
import {
  runComprehensive, setFastpLog, summariseQc, qcVerdict, COMPREHENSIVE_DBS,
} from './comprehensive.js';

const state = {
  readType: 'paired',
  files: { r1: null, r2: null },
  wasmReady: false,
  running: false,
  organism: null, // set when a run's metadata names the species (ENA lookup, example)
  selectedDbs: new Set(['resfinder']),
  dbCached: new Set(), // indexes already in the IndexedDB cache
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
  setupRunOptions();
  // Optional-element wiring is null-safe: a stale cached page (or a cached
  // script against new markup) must never abort the rest of the init.
  document.getElementById('btn-run')?.addEventListener('click', doRun);
  document.getElementById('btn-example')?.addEventListener('click', () => loadExample());
  document.getElementById('btn-fetch-run')?.addEventListener('click', doFetchRun);
  document.getElementById('run-accession')?.addEventListener('keydown', (e) => {
    if (e.key === 'Enter') doFetchRun();
  });
  document.getElementById('btn-back')?.addEventListener('click', () => {
    document.getElementById('workspace').classList.remove('show-results');
  });
  document.getElementById('btn-clear-cache')?.addEventListener('click', clearCachedData);
  setupLogsToggle();
  applyReadType();

  const logToTerm = (msg, level) => term.push(msg, level);
  setLogCallback(logToTerm);

  // Deep link: ?run=<accession> fetches that public run straight away, so a
  // run can be shared or re-analyzed with a single URL.
  const runParam = new URLSearchParams(location.search).get('run');
  if (runParam && /^[ESD]RR\d+$/i.test(runParam.trim())) {
    const input = document.getElementById('run-accession');
    if (input) {
      input.value = runParam.trim();
      doFetchRun();
    }
  }

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

// ── Run options (fastp + database selection) ──
//
// No modes: fastp QC is an independent toggle and any subset of the
// databases can be selected. Defaults: ResFinder only, fastp off.

function selectedDbs() {
  return COMPREHENSIVE_DBS.filter(d => state.selectedDbs.has(d.key)).map(d => d.key);
}

function setupRunOptions() {
  document.getElementById('opt-fastp')?.addEventListener('change', () => {
    refreshFastpPanel();
    updateRunButton();
  });
  // The settings panel never inflates on its own: a "fastp settings"
  // expander appears when fastp is on, and only the user opens it.
  document.getElementById('fastp-toggle')?.addEventListener('click', () => {
    fastpExpanded = !fastpExpanded;
    refreshFastpPanel();
  });
  // Sub-toggles gate their own fields: unchecking "Adapter trimming" greys
  // out the sequence inputs, etc.
  for (const id of ['fp-trim', 'fp-quality', 'fp-length']) {
    document.getElementById(id)?.addEventListener('change', refreshFastpPanel);
  }
  setupDbSelect();
  refreshFastpPanel();
}

// ── fastp settings panel ──
//
// Hidden until the user expands it via the toggle (visible only while the
// fastp checkbox is on). Controls are replaced by a note for Nanopore
// (fastp is Illumina-oriented and runs stats-only there), and paired-only
// fields (R2 adapter, overlap correction) hide otherwise. Field values
// always carry the defaults, so what the user sees is what fastp would use.

let fastpExpanded = false;

function refreshFastpPanel() {
  const cfg = document.getElementById('fastp-config');
  if (!cfg) return;
  const on = document.getElementById('opt-fastp').checked;
  const nanopore = state.readType === 'nanopore';
  const paired = state.readType === 'paired';
  const toggle = document.getElementById('fastp-toggle');
  if (toggle) {
    toggle.hidden = !on;
    toggle.setAttribute('aria-expanded', String(on && fastpExpanded));
  }
  cfg.hidden = !(on && fastpExpanded);
  document.getElementById('fastp-panel').hidden = nanopore;
  document.getElementById('fastp-ont-note').hidden = !nanopore;
  const gate = (fieldsId, enabled) => {
    document.getElementById(fieldsId)?.querySelectorAll('input,select')
      .forEach(el => { el.disabled = !enabled; });
  };
  gate('fp-adapter-fields', document.getElementById('fp-trim').checked);
  gate('fp-quality-fields', document.getElementById('fp-quality').checked);
  gate('fp-length-fields', document.getElementById('fp-length').checked);
  document.getElementById('fp-field-adapter-r2').style.display = paired ? '' : 'none';
  document.getElementById('fp-correction-field').style.display = paired ? '' : 'none';
}

// Read the panel into the options object the fastp worker understands.
// Mirrors the defaults documented in index.html; invalid entries get a
// warning and fall back (the worker re-validates defensively).
const IUPAC_ADAPTER_RE = /^[ACGTURYSWKMBDHVN]{5,100}$/;

function collectFastpOptions() {
  const el = (id) => document.getElementById(id);
  const opts = {
    adapterTrim: el('fp-trim').checked,
    qualityFilter: el('fp-quality').checked,
    lengthFilter: el('fp-length').checked,
    polyG: el('fp-polyg').value || 'auto',
    correction: el('fp-correction').checked,
  };
  if (opts.adapterTrim) {
    const seqs = [['fp-adapter-r1', 'adapterR1', 'adapter R1'],
                  ['fp-adapter-r2', 'adapterR2', 'adapter R2']];
    for (const [id, key, label] of seqs) {
      const input = el(id);
      if (input.disabled || !input.value.trim()) continue;
      const seq = input.value.replace(/\s+/g, '').toUpperCase();
      if (IUPAC_ADAPTER_RE.test(seq)) opts[key] = seq;
      else term.push(`Ignoring ${label}: expected 5–100 IUPAC bases (A C G T U R Y S W K M B D H V N); using auto-detection.`, 'warn');
    }
  }
  if (opts.qualityFilter) {
    opts.qualifiedPhred = parseInt(el('fp-phred').value, 10);
    opts.unqualifiedPercent = parseInt(el('fp-unqual-pct').value, 10);
    opts.nBaseLimit = parseInt(el('fp-nlimit').value, 10);
  }
  if (opts.lengthFilter) {
    opts.minLength = parseInt(el('fp-minlen').value, 10);
    opts.maxLength = parseInt(el('fp-maxlen').value, 10);
  }
  return opts;
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
  document.querySelectorAll('#read-type .seg-btn').forEach(b =>
    b.setAttribute('aria-pressed', String(b.classList.contains('active'))));
  refreshFastpPanel();
  updateRunButton();
}

function setReadType(type) {
  state.readType = type;
  document.querySelectorAll('#read-type .seg-btn').forEach(b => {
    b.classList.toggle('active', b.dataset.type === type);
  });
  applyReadType();
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
  const hint = document.getElementById('run-hint');
  const haveR1 = !!state.files.r1;
  const needR2 = state.readType === 'paired';
  const haveR2 = !!state.files.r2;
  const haveDb = selectedDbs().length > 0;
  const ready = state.wasmReady && !state.running && haveDb && haveR1 && (!needR2 || haveR2);
  btn.disabled = !ready;
  btn.title = haveDb ? '' : 'Select at least one database.';
  if (hint) {
    if (ready || state.running || !state.wasmReady) {
      hint.hidden = true;
    } else if (!haveDb) {
      hint.hidden = false;
      hint.textContent = 'Select at least one database to run the analysis.';
    } else if (!haveR1) {
      hint.hidden = false;
      hint.textContent = 'Load reads to begin — choose files, fetch a public run, or load the example data.';
    } else {
      hint.hidden = false;
      hint.textContent = 'Add the reverse reads (R2) file to run paired-end analysis.';
    }
  }
}

// ── Example data ──
//
// S. aureus USA300_TCH1516 — the community-associated MRSA (ST8) reference
// strain, PVL-positive. Illumina MiSeq run SRR10341524 (~48 MB, ENA study
// PRJNA579343), fetched from ENA on demand and cached in IndexedDB, so
// nothing is served from this site.
const EXAMPLE = {
  accession: 'SRR10341524',
  organism: 'S. aureus USA300_TCH1516',
  study: 'PRJNA579343',
  note: '48 MB MiSeq run, downloaded from ENA and cached in your browser',
  files: [
    {
      slot: 'r1', name: 'SRR10341524_1.fastq.gz', bytes: 22326358,
      url: 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/024/SRR10341524/SRR10341524_1.fastq.gz',
    },
    {
      slot: 'r2', name: 'SRR10341524_2.fastq.gz', bytes: 26024204,
      url: 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/024/SRR10341524/SRR10341524_2.fastq.gz',
    },
  ],
};

async function loadExample() {
  const btn = document.getElementById('btn-example');
  if (btn.disabled) return;
  btn.disabled = true;
  const original = btn.textContent;
  try {
    setReadType('paired');
    state.organism = 'Staphylococcus aureus';
    term.push(`Example isolate: ${EXAMPLE.organism} · ENA/SRA ${EXAMPLE.accession}`, 'info');
    term.push(EXAMPLE.note, 'info');
    term.push('USA300_TCH1516: community-associated MRSA reference strain (ST8), PVL-positive.', 'info');

    for (const spec of EXAMPLE.files) {
      const cacheKey = 'example:' + spec.name;
      let data = await getCachedDBFile(cacheKey);
      if (data && data.byteLength) {
        term.push(`${spec.name} (cached, ${formatBytes(data.byteLength)})`, 'info');
        btn.textContent = 'Example (cached)';
      } else {
        term.push(`Downloading ${spec.name} from ENA…`, 'progress');
        data = await fetchAssetWithProgress(spec.url, spec.bytes, (got, total) => {
          const pct = total ? Math.min(100, Math.round(100 * got / total)) : 0;
          btn.textContent = `Example ↓ ${pct}%`;
        });
        term.push(`${spec.name} (${formatBytes(data.byteLength)})`, 'ok');
        try { await cacheDBFile(cacheKey, data); } catch (_) {
          term.push('Could not cache the FASTQ in IndexedDB; next load will re-download.', 'warn');
        }
      }
      const file = new File([data], spec.name, { type: 'application/gzip' });
      state.files[spec.slot] = file;
      setSlotName(spec.slot, `${file.name} (${formatBytes(file.size)})`);
    }
    updateRunButton();
    term.push('Ready. Choose ResFinder (mecA, blaZ, aph(3\')-III…) or VFDB (PVL, hla, ica…) and run.', 'ok');
  } catch (err) {
    term.push('Failed to load example: ' + err.message, 'error');
    expandLogs();
    term.push('ENA must be reachable from this browser (CORS). Try again later, or upload your own FASTQ.', 'info');
  } finally {
    btn.textContent = original;
    btn.disabled = false;
  }
}

// ── Database multi-select dropdown ──
//
// Selection drives downloading: ResFinder (the default) is fetched as soon
// as the page is ready, other indexes download the first time they're
// selected, and everything stays in the IndexedDB cache. Selected databases
// show a tick in the menu and appear as chips on the trigger button.

const DB_META = {
  resfinder: { short: 'ResFinder', size: '12 MB' },
  card_homolog: { short: 'CARD', size: '24 MB' },
  vfdb_core: { short: 'VFDB', size: '78 MB' },
};

const DB_DOWNLOADS = new Map(); // db key → { pct, promise } while in flight

function setupDbSelect() {
  const btn = document.getElementById('db-select-btn');
  const menu = document.getElementById('db-menu');
  const wrap = document.getElementById('db-select');
  if (!btn || !menu) return;
  const setOpen = (open) => {
    menu.hidden = !open;
    btn.setAttribute('aria-expanded', String(open));
  };
  btn.addEventListener('click', (e) => {
    e.stopPropagation();
    setOpen(menu.hidden);
  });
  document.addEventListener('click', (e) => {
    if (!menu.hidden && !wrap.contains(e.target)) setOpen(false);
  });
  document.addEventListener('keydown', (e) => {
    if (e.key === 'Escape' && !menu.hidden) setOpen(false);
  });
  document.querySelectorAll('.db-option').forEach((row) => {
    const pick = () => toggleDb(row.dataset.db);
    row.addEventListener('click', pick);
    row.addEventListener('keydown', (e) => {
      if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); pick(); }
    });
  });

  renderDbUi();
  // Cache state first, then start the default database downloading right
  // away — before any reads are loaded, so it's ready when they are.
  (async () => {
    for (const db of COMPREHENSIVE_DBS) {
      try { if (await isDatabaseCached(db.key)) state.dbCached.add(db.key); } catch (_) { /* cache check is best-effort */ }
    }
    renderDbUi();
  })();
  downloadDb('resfinder');
}

function toggleDb(key) {
  if (state.selectedDbs.has(key)) {
    state.selectedDbs.delete(key);
  } else {
    state.selectedDbs.add(key);
    downloadDb(key); // no-op when already cached or in flight
  }
  renderDbUi();
  updateRunButton();
}

// Download an index unless it's cached or already downloading. Errors are
// logged but never block: the run itself re-fetches on demand if needed.
// Zenodo serves the files without content-length, so byte-percentage is
// usually unavailable — progress then falls back to file i/n counts.
function downloadDb(key) {
  if (DB_DOWNLOADS.has(key) || state.dbCached.has(key)) return;
  const meta = DB_META[key];
  const entry = { pct: null, file: null, promise: null };
  DB_DOWNLOADS.set(key, entry);
  term.push(`Downloading ${meta.short} index (${meta.size})…`, 'info');
  entry.promise = preloadDatabase(key,
    (got, total) => {
      entry.pct = total ? Math.min(100, Math.round(100 * got / total)) : null;
      renderDbUi();
    },
    (i, n) => {
      entry.file = n > 1 ? `${i + 1}/${n}` : null;
      renderDbUi();
    })
    .then(() => {
      state.dbCached.add(key);
      term.push(`${meta.short} index cached — ready to run.`, 'ok');
    })
    .catch((err) => {
      term.push(`Could not download the ${meta.short} index: ${err.message}`, 'error');
    })
    .finally(() => {
      DB_DOWNLOADS.delete(key);
      renderDbUi();
    });
  renderDbUi();
}

const dlStatus = (dl) => dl.pct != null ? `↓ ${dl.pct}%` : dl.file != null ? `↓ ${dl.file}` : '…';

// Repaint chips (button) and ticks/progress (menu rows) from state.
function renderDbUi() {
  const chips = document.getElementById('db-chips');
  if (!chips) return;
  chips.innerHTML = '';
  for (const db of COMPREHENSIVE_DBS) {
    if (!state.selectedDbs.has(db.key)) continue;
    const meta = DB_META[db.key];
    const dl = DB_DOWNLOADS.get(db.key);
    const chip = document.createElement('span');
    chip.className = 'db-chip' + (dl ? ' loading' : '');
    const name = document.createElement('span');
    name.textContent = meta.short;
    chip.appendChild(name);
    if (dl) {
      const stat = document.createElement('span');
      stat.className = 'db-chip-stat';
      stat.textContent = dlStatus(dl);
      chip.appendChild(stat);
    }
    chips.appendChild(chip);
  }
  if (!chips.children.length) {
    const empty = document.createElement('span');
    empty.className = 'db-chip-empty';
    empty.textContent = 'No databases selected';
    chips.appendChild(empty);
  }
  document.querySelectorAll('.db-option').forEach((row) => {
    const key = row.dataset.db;
    const selected = state.selectedDbs.has(key);
    row.classList.toggle('selected', selected);
    row.setAttribute('aria-selected', String(selected));
    const sizeEl = row.querySelector('.db-opt-size');
    const dl = DB_DOWNLOADS.get(key);
    const bar = row.querySelector('.db-opt-prog > i');
    if (dl) {
      row.classList.add('downloading');
      if (sizeEl) sizeEl.textContent = dlStatus(dl);
      if (bar) bar.style.width = (dl.pct || 3) + '%';
    } else {
      row.classList.remove('downloading');
      if (bar) bar.style.width = '0';
      if (sizeEl) sizeEl.textContent = DB_META[key].size;
    }
  });
}

// ── Clear cached data ──
//
// Wipes the IndexedDB cache: database indexes, CABBAGE snapshots, the example
// isolate and any fetched public runs. Uploaded files are never persisted, so
// there is nothing else to remove. Data already loaded into memory this
// session keeps working until the page is reloaded.
async function clearCachedData() {
  if (state.running) return;
  const btn = document.getElementById('btn-clear-cache');
  if (!btn || btn.disabled) return;
  btn.disabled = true;
  try {
    const n = await clearDBCache();
    state.dbCached.clear();
    renderDbUi();
    term.push(`Cleared ${n} cached item${n === 1 ? '' : 's'} from this browser: database indexes, CABBAGE snapshots and any example or fetched-run reads.`, 'ok');
    term.push('Selected databases re-download the next time they are used. Nothing was uploaded anywhere — the cache was only ever on this device.', 'info');
  } catch (err) {
    term.push('Could not clear the cache: ' + err.message, 'error');
    expandLogs();
  } finally {
    btn.disabled = false;
  }
}

// ── Fetch a public run from ENA ──
//
// Any SRR/ERR/DRR run accession. ENA is queried first for run metadata
// (organism, platform, layout, file sizes); the read type is set from that
// metadata and the run is described to the user before the FASTQs download.
// Files are capped at 1 GB each (see fetch-run.js). Downloaded runs are
// cached in IndexedDB, so fetching the same run again is instant.

const PLATFORM_LABELS = {
  ILLUMINA: 'Illumina',
  OXFORD_NANOPORE: 'Oxford Nanopore',
  ION_TORRENT: 'Ion Torrent',
  PACIFIC_BIOSCIENCES: 'PacBio',
  LS454: '454',
};

async function doFetchRun() {
  if (state.running) return;
  const input = document.getElementById('run-accession');
  const btn = document.getElementById('btn-fetch-run');
  const metaEl = document.getElementById('run-meta');
  const acc = input.value.trim();
  if (!acc) {
    term.push('Enter a run accession, for example SRR10341524.', 'warn');
    input.focus();
    return;
  }
  if (btn.disabled) return;
  btn.disabled = true;
  const original = btn.textContent;
  metaEl.hidden = true;
  try {
    term.push(`Looking up run ${acc} at ENA…`, 'info');
    const run = await lookupRun(acc);
    state.organism = run.organism || null;
    const platform = PLATFORM_LABELS[run.platform] || run.platform || 'unknown platform';
    const kind = run.readType === 'nanopore'
      ? platform
      : `${run.layout === 'paired' ? 'paired-end' : 'single-end'} ${platform}`;
    metaEl.hidden = false;
    metaEl.innerHTML = `<b>${esc(run.organism)}</b> · ${esc([kind, run.instrument].filter(Boolean).join(', '))} · ${formatBytes(run.totalBytes)}`;
    term.push(`Detected: ${[run.organism, kind, run.instrument].filter(Boolean).join(' · ')} (${formatBytes(run.totalBytes)})`, 'ok');

    setReadType(run.readType);
    term.push(`Downloading ${run.files.length} FASTQ file${run.files.length === 1 ? '' : 's'} from ENA…`, 'info');
    for (let i = 0; i < run.files.length; i++) {
      const spec = run.files[i];
      const slot = (run.layout === 'paired' && i === 1) ? 'r2' : 'r1';
      const cacheKey = 'run:' + spec.name;
      let data = await getCachedDBFile(cacheKey);
      if (data && data.byteLength) {
        term.push(`${spec.name} (cached, ${formatBytes(data.byteLength)})`, 'info');
      } else {
        data = await fetchAssetWithProgress(spec.url, spec.bytes, (got, total) => {
          const pct = total ? Math.min(100, Math.round(100 * got / total)) : 0;
          btn.textContent = `${spec.name}  ${formatBytes(got)} / ${formatBytes(total)}  (${pct}%)`;
        });
        term.push(`${spec.name} (${formatBytes(data.byteLength)})`, 'ok');
        try { await cacheDBFile(cacheKey, data); } catch (_) {
          term.push('Could not cache this file in the browser; fetching it again will re-download it.', 'warn');
        }
      }
      const file = new File([data], spec.name, { type: 'application/gzip' });
      state.files[slot] = file;
      setSlotName(slot, `${file.name} (${formatBytes(file.size)})`);
    }
    updateRunButton();
    term.push('Files ready. Start the analysis.', 'ok');
  } catch (err) {
    term.push('Could not fetch run: ' + err.message, 'error');
    expandLogs();
    metaEl.hidden = true;
  } finally {
    btn.textContent = original;
    btn.disabled = false;
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
    const statusText = document.querySelector('#run-status .run-status-text');
    const onStep = (i, n, text) => {
      if (statusText) statusText.textContent = `Step ${i} of ${n}: ${text}`;
      term.push(`[${i}/${n}] ${text}`, 'progress');
    };
    setFastpLog((t) => term.push(t, 'info'));
    const runQc = document.getElementById('opt-fastp').checked;
    const fastpOptions = runQc ? collectFastpOptions() : undefined;
    const thresholds = {
      id_threshold: parseInt(document.getElementById('slider-id').value) / 100,
      mrc: parseInt(document.getElementById('slider-cov').value) / 100,
    };
    const result = await runComprehensive(files, state.readType, {
      onStep, runQc, dbKeys: selectedDbs(), thresholds, fastpOptions,
    });
    if (statusText) statusText.textContent = 'Analyzing… this may take a moment.';
    renderComprehensive(result, { runQc });
    document.getElementById('workspace').classList.add('show-results');
  } catch (err) {
    term.push('Error: ' + err.message, 'error');
  } finally {
    state.running = false;
    setRunning(false);
    updateRunButton();
  }
}

const DB_SHORT = { resfinder: 'ResFinder', card_homolog: 'CARD', vfdb_core: 'VFDB' };

function renderComprehensive({ qc, qcHtml, qcReads, results }, { runQc } = {}) {
  const area = document.getElementById('results');
  const metrics = summariseQc(qc);
  const verdict = qcVerdict(metrics);

  const anyOk = results.some(r => r.exitCode === 0);
  const parts = [];
  if (runQc) parts.push('fastp QC');
  results.forEach(r => parts.push(DB_SHORT[r.database] || r.dbLabel || r.database));
  let html = `
    <div class="result-header ${anyOk ? 'ok' : 'error'}">
      <strong>${anyOk ? 'Analysis complete' : 'Analysis failed'}</strong>
      <span class="result-meta">${parts.join(' + ')}</span>
    </div>`;

  // ── Downloads first: one zip with everything (except the clean reads,
  //    which can be hundreds of MB — those live on the QC card) ──
  html += `
    <div class="downloads">
      <button class="btn btn-primary" data-zip>Download all results (ZIP)</button>
      <p class="dl-nudge">The ZIP has everything from this run: the full result tables, fastp's JSON + HTML quality reports, the CABBAGE phenotype predictions and the run log. The clean reads are downloaded separately from the Quality control card below.</p>
    </div>`;

  // ── QC card (only when fastp was part of the run) ──
  if (runQc) {
    const pairedReads = (qcReads || []).length > 1;
    const readBtns = (qcReads || [])
      .map((f, i) => `<button class="qc-dl-btn" data-qc-dl="read-${i}">${pairedReads ? `Clean R${i + 1}` : 'Clean reads'} ↓ ${formatBytes(f.size)}</button>`)
      .join('');
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
          ${(qc || qcHtml) ? `
          <div class="qc-dl">
            <span class="qc-dl-label">fastp output</span>
            ${qcHtml ? '<button class="qc-dl-btn" data-qc-dl="view">View HTML report ↗</button><button class="qc-dl-btn" data-qc-dl="html">report.html ↓</button>' : ''}
            ${qc ? '<button class="qc-dl-btn" data-qc-dl="json">report.json ↓</button>' : ''}
            ${readBtns}
          </div>
          <p class="opt-note">Clean reads are the trimmed FASTQs the databases were run against — the first 400,000 reads/pairs (WebAssembly memory limit). report.json records the exact fastp command and every filtering statistic.</p>` : ''}
        </div>
      </section>`;
  }

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

  // ── Phenotype prediction (CABBAGE) mounts here ──
  html += '<div id="cabbage-report"></div>';

  // ── Downloads: everything in one zip (reports are small; the clean
  //    reads stay out of it and are served from the QC card) ──
  const dlFiles = {};
  if (qc) {
    dlFiles['fastp_report.json'] = { data: new TextEncoder().encode(JSON.stringify(qc, null, 2)), binary: true };
  }
  if (qcHtml) {
    dlFiles['fastp_report.html'] = { data: new TextEncoder().encode(qcHtml), binary: true };
  }
  for (const r of results) {
    for (const [ext, info] of Object.entries(r.files || {})) {
      if (r.files['.log'] && ext === '.log') continue;
      dlFiles[`${r.database}${ext}`] = info;
    }
  }
  dlFiles['.log'] = { data: term.text(), binary: false };

  area.innerHTML = html;

  // Per-artifact fastp downloads on the QC card. The HTML report is
  // self-contained (inline CSS/JS + data), so it can be viewed in a new
  // tab or saved and shared as-is.
  const qcDl = area.querySelector('.qc-dl');
  if (qcDl) {
    qcDl.addEventListener('click', (e) => {
      const btn = e.target.closest('[data-qc-dl]');
      if (!btn) return;
      const kind = btn.dataset.qcDl;
      if (kind === 'view' && qcHtml) {
        const url = URL.createObjectURL(new Blob([qcHtml], { type: 'text/html' }));
        window.open(url, '_blank', 'noopener');
        setTimeout(() => URL.revokeObjectURL(url), 60_000);
      } else if (kind === 'html' && qcHtml) {
        download(new Blob([qcHtml], { type: 'text/html' }), 'fastp_report.html');
      } else if (kind === 'json' && qc) {
        download(new Blob([JSON.stringify(qc, null, 2)], { type: 'application/json' }), 'fastp_report.json');
      } else if (kind.startsWith('read-')) {
        const f = (qcReads || [])[Number(kind.slice(5))];
        if (f) download(f, f.name);
      }
    });
  }

  // ── Phenotype prediction card (CABBAGE) ──
  // Mounts after the tables, auto-loads the small association table, and
  // contributes cabbage_predictions.csv to the ZIP when ready.
  const cabbageArea = area.querySelector('#cabbage-report');
  if (cabbageArea) {
    const detected = collectDetectedGenes(results, parseResFile);
    mountCabbageReport(cabbageArea, {
      organism: state.organism,
      detectedGenes: detected,
      log: (msg, level) => term.push(msg, level),
      onCsv: (name, text) => {
        dlFiles[name] = { data: new TextEncoder().encode(text), binary: true };
      },
    });
  }

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

// Show a running indicator (spinner in the run button + indeterminate progress
// bar) so it's clear the job is working. The analysis itself runs in a Web
// Worker, so these animations keep ticking while KMA crunches.
function setRunning(on) {
  const btn = document.getElementById('btn-run');
  const status = document.getElementById('run-status');
  const exampleBtn = document.getElementById('btn-example');
  if (exampleBtn) exampleBtn.disabled = on;
  const fetchBtn = document.getElementById('btn-fetch-run');
  if (fetchBtn) fetchBtn.disabled = on;
  if (on) {
    btn.innerHTML = '<span class="spinner"></span>Running…';
    status.hidden = false;
  } else {
    btn.textContent = 'Analyze sample';
    status.hidden = true;
  }
}

// ── Results ──

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
