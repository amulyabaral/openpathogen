import { initWasm, setLogCallback, formatBytes, isDatabaseCached, preloadDatabase } from './wasm-runtime.js';
import { loadPhenotypesDB, parseResFile } from './resfinder-db.js';
import { cacheDBFile, getCachedDBFile, clearDBCache } from './db.js';
import { openGeneViewer } from './gene-view.js';
import { fetchAssetWithProgress } from './assets.js';
import { lookupRun } from './fetch-run.js';
import { loadAssociations, resolveSpecies } from './cabbage.js';
import { describeHit, describeFunction, cardAroUrl } from './genes.js';
import { mountReport, DB_NAMES } from './report.js';
import {
  runComprehensive, setFastpLog, summariseQc, qcVerdict, fmtPct, COMPREHENSIVE_DBS,
} from './comprehensive.js';

const state = {
  readType: 'paired',
  files: { r1: null, r2: null },
  wasmReady: false,
  running: false,
  species: '',        // the species box; filled in for the example and ENA runs
  speciesAuto: false, // true while that text came from run metadata, not the user
  selectedDbs: new Set(['resfinder']),
  dbCached: new Set(), // indexes already in the IndexedDB cache
};

// The results on screen: files for the ZIP, the report (for species
// changes), and the sample name for the printed report's title.
let current = null;

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
  setupSpecies();
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
    window.scrollTo(0, 0); // stacked (phone) layout: back to the top of the form
  });
  document.getElementById('btn-print')?.addEventListener('click', () => window.print());
  document.getElementById('btn-zip')?.addEventListener('click', downloadZip);
  document.getElementById('btn-clear-cache')?.addEventListener('click', clearCachedData);
  setupPrint();
  setupLogsToggle();
  applyReadType();
  setupSettingsColumn();

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
    updateRunButton();
  } catch (err) {
    term.push('Init failed: ' + err.message, 'error');
    expandLogs();
  }
});

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
  // expander appears beside the checkbox when fastp is on (reserving its
  // space, so nothing moves), and only the user opens it.
  document.getElementById('fastp-toggle')?.addEventListener('click', (e) => {
    fastpExpanded = !fastpExpanded;
    refreshFastpPanel();
    // Bring the opened panel into view within the Settings column.
    const col = document.getElementById('settings-col');
    if (fastpExpanded && col && col.scrollHeight > col.clientHeight) {
      col.scrollTo({ top: e.currentTarget.offsetTop - col.offsetTop - 60, behavior: 'smooth' });
    }
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
    toggle.classList.toggle('is-off', !on);
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

// ── Settings column height ──
//
// On wide windows the whole form fits on screen. The Settings column is
// capped at the height left in the window and scrolls on its own when it
// grows (fastp settings open), so the Run button stays in view; once the
// column reaches its end, the scroll carries on to the page.

function fitSettingsColumn() {
  const col = document.getElementById('settings-col');
  const pane = document.querySelector('.pane-config');
  const cols = col?.parentElement;
  const card = document.getElementById('run');
  if (!col || !pane || !cols || !card) return;
  if (!window.matchMedia('(min-width:960px)').matches) {
    col.style.maxHeight = '';
    return;
  }
  const paneTop = pane.getBoundingClientRect().top - pane.scrollTop;
  const colTop = col.getBoundingClientRect().top - paneTop;
  // Everything below the columns: Run row, card padding and border, and the
  // pane's bottom padding. None of it depends on the column's height.
  const below = card.getBoundingClientRect().bottom - cols.getBoundingClientRect().bottom
    + parseFloat(getComputedStyle(pane).paddingBottom);
  const left = cols.firstElementChild?.offsetHeight || 0;
  col.style.maxHeight = Math.max(left, 260, Math.floor(pane.clientHeight - colTop - below)) + 'px';
}

function setupSettingsColumn() {
  fitSettingsColumn();
  window.addEventListener('resize', fitSettingsColumn);
  // The Reads column changes height as files are added or removed.
  const reads = document.getElementById('settings-col')?.previousElementSibling;
  if (reads && 'ResizeObserver' in window) new ResizeObserver(fitSettingsColumn).observe(reads);
}

// ── Read type ──
//
// Set automatically from the files (two files: paired; one file: single-end
// or Nanopore by its read lengths) or from ENA metadata; the toggle is there
// to correct it. The note under the toggle says where the setting came from.

const READ_TYPE_HINT = 'Set automatically from the files you add. Change it here if needed.';

function setupReadTypeToggle() {
  document.querySelectorAll('#read-type .seg-btn').forEach(btn => {
    btn.addEventListener('click', () => setReadType(btn.dataset.type, READ_TYPE_HINT));
  });
}

function applyReadType() {
  if (state.readType !== 'paired') state.files.r2 = null;
  document.querySelectorAll('#read-type .seg-btn').forEach(b => {
    const on = b.dataset.type === state.readType;
    b.classList.toggle('active', on);
    b.setAttribute('aria-pressed', String(on));
  });
  renderFiles();
  refreshFastpPanel();
  updateRunButton();
}

function setReadType(type, note) {
  state.readType = type;
  const el = document.getElementById('read-type-note');
  if (el && note) el.textContent = note;
  applyReadType();
}

// ── Uploads: file dialog, drag and drop, pairing ──
//
// One picker for everything: select or drop one file or a pair. Buttons
// marked data-pick="r2" fill that slot only (the "Add R2 file" row).

let pickSlot = null;

function setupUploadSlots() {
  const input = document.getElementById('file-input');
  if (!input) return;
  document.querySelectorAll('[data-pick]').forEach(btn => {
    btn.addEventListener('click', () => {
      pickSlot = btn.dataset.pick || null;
      input.click();
    });
  });
  input.addEventListener('change', () => {
    if (input.files.length) assignFiles([...input.files], pickSlot);
    input.value = '';
  });
  document.querySelectorAll('[data-clear]').forEach(btn => {
    btn.addEventListener('click', () => {
      if (state.running) return;
      state.files[btn.dataset.clear] = null;
      showUploadNote('');
      setReadType(state.readType, READ_TYPE_HINT);
    });
  });
  document.getElementById('btn-swap')?.addEventListener('click', () => {
    if (state.running) return;
    const { r1, r2 } = state.files;
    state.files.r1 = r2;
    state.files.r2 = r1;
    renderFiles();
  });
  setupDropZone();
}

// Files can be dropped anywhere on the input page (a drop elsewhere would
// make the browser open the file and leave the app). A drop on the R2 row
// goes to that slot; otherwise names decide.
function setupDropZone() {
  const grid = document.getElementById('upload-grid');
  if (!grid) return;
  let depth = 0;
  const hasFiles = (e) => [...(e.dataTransfer?.types || [])].includes('Files');
  const onForm = () => !state.running && !document.getElementById('workspace').classList.contains('show-results');
  window.addEventListener('dragenter', (e) => {
    if (!hasFiles(e)) return;
    depth++;
    if (onForm()) grid.classList.add('dragging');
  });
  window.addEventListener('dragleave', (e) => {
    if (!hasFiles(e)) return;
    depth = Math.max(0, depth - 1);
    if (!depth) grid.classList.remove('dragging');
  });
  window.addEventListener('dragover', (e) => {
    if (!hasFiles(e)) return;
    e.preventDefault();
    e.dataTransfer.dropEffect = onForm() ? 'copy' : 'none';
  });
  window.addEventListener('drop', (e) => {
    if (!hasFiles(e)) return;
    e.preventDefault();
    depth = 0;
    grid.classList.remove('dragging');
    if (!onForm()) return;
    const slot = e.target.closest?.('[data-slot]')?.dataset.slot;
    assignFiles([...e.dataTransfer.files], slot);
  });
}

// Read number from a FASTQ name: SRR1_1.fastq.gz, s_S1_L001_R2_001.fastq.gz,
// x.R1.fq → 1 or 2; 0 when the name does not say.
function readNumber(name) {
  const m = String(name).match(/[_.\-\s]R?([12])(?:_\d{3})?\.(?:fastq|fq)(?:\.gz)?$/i);
  return m ? Number(m[1]) : 0;
}

const FASTQ_NAME = /\.(fastq|fq)(\.gz)?$|\.gz$/i;

// Put chosen or dropped files into the slots. Two files make a pair (the
// read type switches to paired). A single file named as one half of a pair
// joins or starts a pair; any other single file is single-end or Nanopore,
// decided by its read lengths.
function assignFiles(list, slotHint) {
  if (state.running) return;
  const files = list.filter(f => FASTQ_NAME.test(f.name));
  const notes = [];
  if (files.length < list.length) notes.push(`Skipped ${list.length - files.length} file${list.length - files.length === 1 ? '' : 's'} that ${list.length - files.length === 1 ? 'is' : 'are'} not FASTQ (.fastq, .fq or .gz).`);
  if (!files.length) { showUploadNote(notes.join(' ')); return; }
  if (files.length >= 2) {
    let [a, b] = files;
    const na = readNumber(a.name), nb = readNumber(b.name);
    if (na === 2 || nb === 1) [a, b] = [b, a];
    else if (!na && !nb && a.name.localeCompare(b.name) > 0) [a, b] = [b, a];
    if (files.length > 2) notes.push(`Only two files can be used: ${a.name} and ${b.name}.`);
    state.files.r1 = a;
    state.files.r2 = b;
    setReadType('paired', (na || nb)
      ? 'Paired-end: R1 and R2 were matched by file name.'
      : 'Paired-end: two files. Check that R1 and R2 are the right way round.');
  } else {
    const f = files[0];
    const n = readNumber(f.name);
    const pairing = state.readType === 'paired' && (slotHint || n || (state.files.r1 && !state.files.r2));
    if (pairing) {
      const slot = slotHint || (n === 2 ? 'r2' : n === 1 ? 'r1' : (state.files.r1 && !state.files.r2 ? 'r2' : 'r1'));
      state.files[slot] = f;
      renderFiles();
    } else {
      state.files.r1 = f;
      state.files.r2 = null;
      setReadType(state.readType === 'nanopore' ? 'nanopore' : 'single', 'Checking read lengths…');
      detectLongReads(f);
    }
  }
  showUploadNote(notes.join(' '));
  // The user's own files: forget the species and metadata of a previous
  // example or ENA fetch so the report does not inherit them.
  forgetAutoSpecies();
  const meta = document.getElementById('run-meta');
  if (meta) { meta.hidden = true; meta.textContent = ''; }
  updateRunButton();
}

// Single-end Illumina or Nanopore? Illumina reads are at most 300 bp, so
// the median length of the first reads decides. Only the start of the file
// is read (and decompressed); the user can still change the type.
async function detectLongReads(file) {
  let median = 0;
  try {
    median = await sniffMedianReadLength(file);
  } catch (_) { /* unreadable start: leave the type as it is */ }
  if (state.files.r1 !== file || state.readType === 'paired') return; // replaced meanwhile
  if (!median) {
    setReadType(state.readType, READ_TYPE_HINT);
    return;
  }
  const len = median >= 1000 ? (median / 1000).toFixed(1) + ' kb' : median + ' bp';
  if (median > 400) setReadType('nanopore', `Nanopore: the reads are long (median ${len}).`);
  else setReadType('single', `Single-end: one file of short reads (median ${len}). Choose Nanopore if that is wrong.`);
}

async function sniffMedianReadLength(file) {
  let bytes = new Uint8Array(await file.slice(0, 1 << 20).arrayBuffer());
  if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
    const parts = [];
    let size = 0;
    const gz = new window.fflate.Gunzip((chunk) => { parts.push(chunk); size += chunk.length; });
    try { gz.push(bytes, false); } catch (_) { /* a truncated block ends the sample early */ }
    bytes = new Uint8Array(size);
    let at = 0;
    for (const p of parts) { bytes.set(p, at); at += p.length; }
  }
  const lines = new TextDecoder().decode(bytes).split('\n');
  const partial = lines.pop(); // the last line is probably cut off
  const lens = [];
  for (let i = 1; i < lines.length && lens.length < 500; i += 4) {
    if (lines[i - 1].startsWith('@')) lens.push(lines[i].trim().length);
  }
  // One Nanopore read can be longer than the whole sample.
  if (!lens.length) return lines.length === 1 && lines[0].startsWith('@') ? partial.length : 0;
  lens.sort((a, b) => a - b);
  return lens[lens.length >> 1];
}

// Draw the drop zone: the prompt when empty, otherwise one row per file
// (plus a row asking for R2 while a pair is incomplete).
function renderFiles() {
  const { r1, r2 } = state.files;
  const paired = state.readType === 'paired';
  const any = !!(r1 || (paired && r2));
  const set = (sel, fn) => { const el = document.querySelector(sel); if (el) fn(el); };
  set('#dz-empty', el => { el.hidden = any; });
  set('#file-list', el => { el.hidden = !any; });
  set('#file-actions', el => { el.hidden = !any; });
  set('#btn-swap', el => { el.hidden = !(paired && r1 && r2); });
  set('#upload-grid', el => el.classList.toggle('has-files', any));
  for (const slot of ['r1', 'r2']) {
    const file = state.files[slot];
    set(`.file-row[data-slot="${slot}"]`, row => {
      row.hidden = slot === 'r2' ? !paired : false;
      row.classList.toggle('missing', !file);
    });
    set(`[data-name="${slot}"]`, el => {
      el.textContent = file ? file.name : `No ${slot === 'r1' ? 'forward (R1)' : 'reverse (R2)'} file`;
      el.title = file ? file.name : '';
      el.dataset.size = file ? formatBytes(file.size) : '';
    });
    set(`[data-clear="${slot}"]`, el => { el.hidden = !file; });
  }
  for (const slot of ['r1', 'r2']) set(`[data-pick="${slot}"]`, el => { el.hidden = !!state.files[slot]; });
  set('[data-tag="r1"]', el => { el.textContent = paired ? 'R1' : state.readType === 'nanopore' ? 'ONT' : 'SE'; });
}

function setFile(slot, file) {
  state.files[slot] = file;
  renderFiles();
}

function showUploadNote(text) {
  const el = document.getElementById('upload-note');
  if (!el) return;
  el.hidden = !text;
  el.textContent = text || '';
}

// ── Species ──
//
// Optional. It adds CABBAGE tested-resistance rates to the report and can
// also be set or changed on the results page. The list of species CABBAGE
// covers loads in the background for the suggestions.

function setupSpecies() {
  const input = document.getElementById('species');
  if (!input) return;
  input.addEventListener('input', () => {
    state.species = input.value.trim();
    state.speciesAuto = false;
  });
  input.addEventListener('change', () => current?.report?.setSpecies(state.species));
  loadAssociations()
    .then((assoc) => {
      const list = document.getElementById('species-list');
      if (list) list.innerHTML = assoc.species.map(s => `<option value="${esc(s)}"></option>`).join('');
    })
    .catch((err) => term.push('CABBAGE species list unavailable: ' + err.message, 'warn'));
}

function setSpecies(text, auto) {
  state.species = text || '';
  state.speciesAuto = !!(auto && text);
  const input = document.getElementById('species');
  if (input) input.value = state.species;
}

function forgetAutoSpecies() {
  if (state.speciesAuto) setSpecies('', false);
}

// ENA names ("Staphylococcus aureus subsp. aureus USA300_TCH1516") are
// shortened to the CABBAGE species when one matches.
async function speciesFromOrganism(organism) {
  if (!organism || /^unknown/i.test(organism)) return '';
  try {
    const assoc = await loadAssociations();
    return resolveSpecies(assoc.species, organism) || organism;
  } catch (_) {
    return organism;
  }
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
      hint.textContent = 'Load reads to begin: choose or drop files, fetch a public run, or try the example.';
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
  species: 'Staphylococcus aureus',
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
  // Same thin bar as the database rows, tracking both FASTQs as one download.
  const prog = document.getElementById('example-prog');
  const bar = prog?.querySelector('i');
  const totalBytes = EXAMPLE.files.reduce((n, f) => n + f.bytes, 0);
  let doneBytes = 0;
  const setPct = (pct) => {
    btn.textContent = `Example ↓ ${pct}%`;
    if (bar) bar.style.width = Math.max(pct, 3) + '%';
  };
  try {
    setReadType('paired', 'Paired-end: the example is an Illumina MiSeq run.');
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
        if (prog) prog.hidden = false;
        setPct(Math.round(100 * doneBytes / totalBytes));
        data = await fetchAssetWithProgress(spec.url, spec.bytes, (got) => {
          setPct(Math.min(100, Math.round(100 * (doneBytes + Math.min(got, spec.bytes)) / totalBytes)));
        });
        term.push(`${spec.name} (${formatBytes(data.byteLength)})`, 'ok');
        try { await cacheDBFile(cacheKey, data); } catch (_) {
          term.push('Could not cache the FASTQ in the browser; the next load will download it again.', 'warn');
        }
      }
      doneBytes += spec.bytes;
      const file = new File([data], spec.name, { type: 'application/gzip' });
      setFile(spec.slot, file);
    }
    setSpecies(EXAMPLE.species, true);
    showUploadNote('');
    updateRunButton();
    term.push('Ready. Choose ResFinder (mecA, blaZ, aph(3\')-III…) or VFDB (PVL, hla, ica…) and run.', 'ok');
  } catch (err) {
    term.push('Failed to load example: ' + err.message, 'error');
    expandLogs();
    term.push('ENA could not be reached from this browser. Try again later, or choose your own FASTQ files.', 'info');
  } finally {
    btn.textContent = original;
    btn.disabled = false;
    if (prog) prog.hidden = true;
    if (bar) bar.style.width = '0';
  }
}

// ── Databases ──
//
// Three checkboxes. Selection drives downloading: ResFinder (the default) is
// fetched as soon as the page is ready, other indexes download the first
// time they're selected, and everything stays in the IndexedDB cache. Each
// row shows the index size, its download progress, or that it is cached.

const DB_META = {
  resfinder: { short: 'ResFinder', size: '12 MB' },
  card_homolog: { short: 'CARD', size: '24 MB' },
  vfdb_core: { short: 'VFDB', size: '78 MB' },
};

const DB_DOWNLOADS = new Map(); // db key → { pct, promise } while in flight

function setupDbSelect() {
  const boxes = document.querySelectorAll('#db-list input[type="checkbox"]');
  if (!boxes.length) return;
  boxes.forEach((cb) => {
    cb.checked = state.selectedDbs.has(cb.value);
    cb.addEventListener('change', () => toggleDb(cb.value, cb.checked));
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

function toggleDb(key, on) {
  if (on) {
    state.selectedDbs.add(key);
    downloadDb(key); // no-op when already cached or in flight
  } else {
    state.selectedDbs.delete(key);
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
      term.push(`${meta.short} index cached.`, 'ok');
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

const dlStatus = (dl) => dl.pct != null ? `↓ ${dl.pct}%` : dl.file != null ? `↓ ${dl.file}` : '↓ …';

// Repaint each database row's size / progress / cached state.
function renderDbUi() {
  document.querySelectorAll('#db-list .db-check').forEach((row) => {
    const key = row.dataset.db;
    const dl = DB_DOWNLOADS.get(key);
    const sizeEl = row.querySelector('.db-size');
    const bar = row.querySelector('.db-prog > i');
    row.classList.toggle('downloading', !!dl);
    if (dl) {
      if (sizeEl) sizeEl.textContent = dlStatus(dl);
      if (bar) bar.style.width = (dl.pct || 3) + '%';
    } else {
      if (sizeEl) sizeEl.textContent = DB_META[key].size + (state.dbCached.has(key) ? ' · cached' : '');
      if (bar) bar.style.width = '0';
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
    term.push('Selected databases download again the next time they are used. Nothing was uploaded; the cache was only on this device.', 'info');
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
// (organism, platform, layout, file sizes); the read type and species are
// set from that metadata and the run is described to the user before the
// FASTQs download. Files are capped at 1 GB each (see fetch-run.js).
// Downloaded runs are cached in IndexedDB, so fetching the same run again
// is instant.

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
    const platform = PLATFORM_LABELS[run.platform] || run.platform || 'unknown platform';
    const kind = run.readType === 'nanopore'
      ? platform
      : `${run.layout === 'paired' ? 'paired-end' : 'single-end'} ${platform}`;
    metaEl.hidden = false;
    metaEl.innerHTML = `<b>${esc(run.organism)}</b> · ${esc([kind, run.instrument].filter(Boolean).join(', '))} · ${formatBytes(run.totalBytes)}`;
    term.push(`Detected: ${[run.organism, kind, run.instrument].filter(Boolean).join(' · ')} (${formatBytes(run.totalBytes)})`, 'ok');

    setReadType(run.readType, `${kind[0].toUpperCase() + kind.slice(1)}, from the ENA run record.`);
    setSpecies(await speciesFromOrganism(run.organism), true);
    showUploadNote('');
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
          term.push('Could not cache this file in the browser; fetching it again will download it again.', 'warn');
        }
      }
      setFile(slot, new File([data], spec.name, { type: 'application/gzip' }));
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
    // A database selected moments ago may still be downloading; wait for it
    // instead of fetching the same index a second time in parallel.
    const pending = selectedDbs().map(k => DB_DOWNLOADS.get(k)).filter(Boolean);
    if (pending.length) {
      if (statusText) statusText.textContent = 'Waiting for the database download to finish…';
      term.push('Waiting for the index download to finish…', 'info');
      await Promise.all(pending.map(p => p.promise));
    }
    const result = await runComprehensive(files, state.readType, {
      onStep, runQc, dbKeys: selectedDbs(), thresholds, fastpOptions,
    });
    renderComprehensive(result, {
      runQc, thresholds, readType: state.readType, inputNames: files.map(f => f.name),
    });
    document.getElementById('workspace').classList.add('show-results');
    window.scrollTo(0, 0); // stacked (phone) layout: results start at the top
  } catch (err) {
    term.push('Error: ' + err.message, 'error');
    expandLogs();
  } finally {
    state.running = false;
    setRunning(false);
    updateRunButton();
  }
}

// Sample name from the input files (the analysis itself may have run on
// fastp's filtered_1/2 files): SRR1_1.fastq.gz → SRR1.
function sampleNameOf(names) {
  const n = names[0] || 'sample';
  const m = n.match(/[_.](?:R?[12])(?:_\d{3})?(?=\.(?:fastq|fq)(?:\.gz)?$)/i);
  return (m ? n.slice(0, m.index) : n.replace(/\.(fastq|fq)(\.gz)?$/i, '')) || n;
}

// Results page: the report first, then one collapsible gene table per
// database and the fastp card. The toolbar prints the page or downloads the
// raw output as a ZIP.
function renderComprehensive({ qc, qcHtml, qcReads, results }, { runQc, thresholds, readType, inputNames }) {
  const area = document.getElementById('results');
  const metrics = summariseQc(qc);
  const verdict = qcVerdict(metrics);
  const sampleName = sampleNameOf(inputNames);

  let html = '<div id="report"></div><div class="evidence">';

  // ── One gene table per database ──
  results.forEach((r, i) => {
    const ok = r.exitCode === 0;
    const rows = ok && r.resTable ? parseResFile(r.resTable).rows : [];
    const name = DB_NAMES[r.database] || r.dbLabel || r.database;
    const meta = !ok ? 'failed' : `${rows.length} gene${rows.length === 1 ? '' : 's'}`;
    html += `
    <details class="ev comp-section" data-section="${i}"${ok ? '' : ' open'}>
      <summary><span class="ev-name">Gene table · ${esc(name)}</span><span class="ev-meta${ok ? '' : ' ev-failed'}">${meta}</span></summary>
      <div class="ev-body">
        ${!ok ? '<p class="empty">This step failed. Open the Logs panel below for details.</p>'
          : rows.length ? renderResTable(rows, r.database)
          : '<p class="empty">No genes detected.</p>'}
      </div>
    </details>`;
  });

  // ── fastp card (only when fastp was part of the run) ──
  if (runQc) {
    const pairedReads = (qcReads || []).length > 1;
    const readBtns = (qcReads || [])
      .map((f, i) => `<button class="qc-dl-btn" type="button" data-qc-dl="read-${i}">${pairedReads ? `Clean R${i + 1}` : 'Clean reads'} ↓ ${formatBytes(f.size)}</button>`)
      .join('');
    html += `
    <details class="ev ev-qc">
      <summary><span class="ev-name">Quality control · fastp</span><span class="ev-meta qc-${esc(verdict.tone)}">${metrics ? `${fmtPct(metrics.retained)} of reads kept` : 'did not run'}</span></summary>
      <div class="ev-body">
        ${metrics ? `
        <div class="qc-grid">
          <div class="qc-cell"><span class="qc-val">${fmtCount(metrics.rawReads)}</span><span class="qc-key">reads in</span></div>
          <div class="qc-cell"><span class="qc-val">${fmtPct(metrics.retained)}</span><span class="qc-key">retained after trimming</span></div>
          <div class="qc-cell"><span class="qc-val">${fmtPct(metrics.q30After ?? 0)}</span><span class="qc-key">Q30 (after)</span></div>
          <div class="qc-cell"><span class="qc-val">${fmtPct(metrics.gcBefore ?? 0)}</span><span class="qc-key">GC content</span></div>
          ${metrics.duplication != null ? `<div class="qc-cell"><span class="qc-val">${fmtPct(metrics.duplication)}</span><span class="qc-key">duplication</span></div>` : ''}
        </div>
        <p class="qc-verdict qc-${verdict.tone}">${esc(verdict.text)}</p>` : `
        <p class="qc-verdict qc-warn">${esc(verdict.text)}</p>`}
        ${(qc || qcHtml) ? `
        <div class="qc-dl">
          <span class="qc-dl-label">fastp output</span>
          ${qcHtml ? '<button class="qc-dl-btn" type="button" data-qc-dl="view">View HTML report ↗</button><button class="qc-dl-btn" type="button" data-qc-dl="html">report.html ↓</button>' : ''}
          ${qc ? '<button class="qc-dl-btn" type="button" data-qc-dl="json">report.json ↓</button>' : ''}
          ${readBtns}
        </div>
        <p class="opt-note">Clean reads are the trimmed FASTQs the databases were run against, covering the whole sample. report.json records the exact fastp command and every filtering statistic.</p>` : ''}
      </div>
    </details>`;
  }
  html += '</div>';

  // ── ZIP contents: raw output (the clean reads, which can be hundreds of
  //    MB, stay out; they download from the fastp card) ──
  const dlFiles = {};
  if (qc) {
    dlFiles['fastp_report.json'] = { data: new TextEncoder().encode(JSON.stringify(qc, null, 2)), binary: true };
  }
  if (qcHtml) {
    dlFiles['fastp_report.html'] = { data: new TextEncoder().encode(qcHtml), binary: true };
  }
  for (const r of results) {
    for (const [ext, info] of Object.entries(r.files || {})) {
      dlFiles[`${r.database}${ext}`] = info;
    }
  }
  dlFiles['openpathogen_run.log'] = { data: term.text(), binary: false };

  area.innerHTML = html;
  current = { dlFiles, sampleName, report: null };

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

  // ── The report ──
  const openHit = (section, template) => {
    const r = results[section];
    const row = r?.resTable ? parseResFile(r.resTable).rows.find(x => x.Template === template) : null;
    if (row) openGeneViewer(row, r);
  };
  current.report = mountReport(area.querySelector('#report'), {
    results, readType, inputNames, sampleName, thresholds,
    qc: runQc ? { tone: verdict.tone, text: verdict.text } : null,
    species: state.species,
    onSpeciesChange: (text) => setSpecies(text, false),
    onCsv: (name, text) => {
      if (text == null) delete dlFiles[name];
      else dlFiles[name] = { data: new TextEncoder().encode(text), binary: true };
    },
    openHit,
    log: (msg, level) => term.push(msg, level),
  });

  // Table bindings, scoped to each database section.
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
    moreRow.className = 'show-more-row';
    moreRow.innerHTML = `<td colspan="${cols}"><button class="show-more-btn" type="button">Show all ${trs.length} genes</button></td>`;
    section.querySelector('tbody').appendChild(moreRow);
    const reveal = () => {
      extra.forEach(tr => tr.classList.remove('row-limited'));
      moreRow.remove();
    };
    moreRow.querySelector('button').addEventListener('click', reveal);
    section.querySelectorAll('th').forEach(th => th.addEventListener('click', reveal));
  });
}

function downloadZip() {
  if (!current) return;
  const enc = new TextEncoder();
  const zip = {};
  for (const [fname, info] of Object.entries(current.dlFiles)) {
    zip[fname] = info.binary
      ? (info.data instanceof Uint8Array ? info.data : new Uint8Array(info.data))
      : enc.encode(info.data);
  }
  download(new Blob([window.fflate.zipSync(zip)], { type: 'application/zip' }), 'openpathogen_results.zip');
}

// Printing (the toolbar button or the browser's own Print) lays out the
// results page as a report: collapsed sections open, every table row shows,
// and the document title becomes the default PDF file name.
function setupPrint() {
  let restore = null;
  window.addEventListener('beforeprint', () => {
    if (!current || !document.getElementById('workspace').classList.contains('show-results')) return;
    const closed = [...document.querySelectorAll('#results details:not([open])')];
    closed.forEach(d => { d.open = true; });
    const title = document.title;
    document.title = `openpathogen report - ${current.sampleName}`;
    restore = () => {
      closed.forEach(d => { d.open = false; });
      document.title = title;
    };
  });
  window.addEventListener('afterprint', () => {
    restore?.();
    restore = null;
  });
}

// Exact counts with thousands separators (515,164), never rounded to "515k".
function fmtCount(n) {
  return n == null ? '—' : Number(n).toLocaleString('en-US');
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
    btn.textContent = 'Run';
    status.hidden = true;
  }
}

// ── Gene tables ──

// Open the per-gene viewer when a row is clicked, or when its gene button
// is activated from the keyboard (the button's click bubbles to the row).
function bindGeneViewer(area, result) {
  if (!result.resTable) return;
  const { rows } = parseResFile(result.resTable);
  const byTemplate = new Map(rows.map(r => [r.Template, r]));
  area.querySelectorAll('tr[data-template]').forEach(tr => {
    tr.addEventListener('click', (e) => {
      if (e.target.closest('a')) return; // links (CARD ARO) open on their own
      // Don't hijack a text selection drag.
      if (window.getSelection && String(window.getSelection()).length) return;
      const row = byTemplate.get(tr.dataset.template);
      if (row) openGeneViewer(row, result);
    });
  });
}

const fmtNum = (v, digits, suffix) => {
  const x = parseFloat(v);
  return Number.isFinite(x) ? x.toFixed(digits) + suffix : '';
};

// Gene | what it does | identity | coverage | depth. The raw template name
// is the gene button's tooltip and the gene viewer's subtitle.
function renderResTable(rows, db) {
  const cols = [
    { key: 'gene', label: 'Gene', width: 150 },
    { key: 'fn', label: db === 'vfdb_core' ? 'Virulence factor' : 'Drug class and antibiotics' },
    { key: 'id', label: 'Identity', width: 96 },
    { key: 'cov', label: 'Coverage', width: 104 },
    { key: 'depth', label: 'Depth', width: 84 },
  ];
  let html = '<div class="table-wrap"><table class="res-table"><colgroup>';
  for (const c of cols) html += c.width ? `<col style="width:${c.width}px">` : '<col>';
  html += '</colgroup><thead><tr>';
  for (const c of cols) {
    html += `<th data-key="${c.key}" aria-sort="none"><button type="button" class="th-sort">${esc(c.label)}<span class="arr" aria-hidden="true">↕</span></button></th>`;
  }
  html += '</tr></thead><tbody>';

  for (const row of rows) {
    const hit = describeHit(row.Template, db, row);
    let fn = esc(describeFunction(hit));
    if (db === 'card_homolog') {
      if (!fn) fn = '<span class="muted">No drug class in ResFinder</span>';
      const url = cardAroUrl(hit.aro);
      if (url) fn += ` · <a href="${url}" target="_blank" rel="noopener">ARO:${esc(hit.aro)}</a>`;
    }
    const cells = {
      gene: [hit.name, `<button type="button" class="gv-link" title="${esc(row.Template)}">${esc(hit.name)}</button>`],
      fn: [describeFunction(hit), fn],
      id: [row.Template_Identity, fmtNum(row.Template_Identity, 1, '%')],
      cov: [row.Template_Coverage, fmtNum(row.Template_Coverage, 1, '%')],
      depth: [row.Depth, fmtNum(row.Depth, 1, '×')],
    };
    html += `<tr data-template="${esc(row.Template || '')}">`;
    for (const c of cols) {
      const [v, inner] = cells[c.key];
      html += `<td data-k="${c.key}" data-v="${esc(v || '')}"${c.key === 'fn' ? ' class="fn"' : ''}>${inner}</td>`;
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
    th.querySelector('.th-sort')?.addEventListener('click', () => {
      const key = th.dataset.key;
      sortDir = (sortKey === key) ? -sortDir : 1;
      sortKey = key;
      table.querySelectorAll('th').forEach(h => {
        const on = h === th;
        h.classList.toggle('sort', on);
        h.setAttribute('aria-sort', on ? (sortDir > 0 ? 'ascending' : 'descending') : 'none');
        const arr = h.querySelector('.arr');
        if (arr) arr.textContent = on ? (sortDir > 0 ? '↑' : '↓') : '↕';
      });
      const rows = [...tbody.querySelectorAll('tr[data-template]')];
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
    // A click on the handle must not sort the column.
    handle.addEventListener('click', (e) => e.stopPropagation());
  });
}

function download(blob, name) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = name;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 10_000);
}

function esc(s) {
  return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
}
