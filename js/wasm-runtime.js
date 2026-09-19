import { getCachedDBFile, cacheDBFile, deleteDBFile } from './db.js';
import { fetchAssetWithProgress } from './assets.js';
import { verifyPinned } from './integrity.js';
import { ENGINE, ENGINE_FORCED, threadCount } from './engine.js';

const DATABASES = {
  resfinder: {
    prefix: 'databases/kma_index_resfinder_2_6_0',
    shortName: 'resfinder',
    displayName: 'ResFinder 2.6.0',
  },
  card_homolog: {
    prefix: 'databases/kma_index_card_4_0_1_homolog',
    shortName: 'card_homolog',
    displayName: 'CARD 4.0.1: Protein Homolog',
  },
  vfdb_core: {
    prefix: 'databases/VFDB_setA_nt.fas.gz',
    shortName: 'vfdb_core',
    displayName: 'VFDB Set A (core virulence factors)',
  },
};
const DB_FILES = ['.comp.b', '.length.b', '.name', '.seq.b'];
const OUTPUT_EXTS = ['.res', '.frag.gz', '.aln', '.fsa', '.vcf.gz', '.mat.gz', '.mapstat', '.sam'];

let logCallback = null;

export function setLogCallback(fn) { logCallback = fn; }

function log(msg, level = 'info') {
  if (logCallback) logCallback(msg, level);
}

export function getDatabaseList() {
  return Object.entries(DATABASES).map(([key, db]) => ({ key, ...db }));
}

// ── Init / readiness ──
//
// Both builds (wasm64, and the wasm32 fallback for browsers without Memory64;
// see engine.js) use pthreads, which need cross-origin isolation
// (SharedArrayBuffer). The coi-serviceworker shim (loaded first in
// index.html) establishes it, reloading once on first visit. The KMA module
// itself is loaded per-run inside a Web Worker, so there is nothing heavy to
// pre-initialise here.

export async function initWasm() {
  if (typeof Worker === 'undefined') {
    throw new Error('Web Workers are unavailable in this browser.');
  }
  if (self.crossOriginIsolated) {
    const threads = `up to ${threadCount()} thread${threadCount() === 1 ? '' : 's'}`;
    log(ENGINE === 'wasm64'
      ? `Ready. Memory64 build, ${threads}`
      : `Ready. 32-bit build (${ENGINE_FORCED ? 'forced by ?engine=wasm32' : 'this browser has no WebAssembly Memory64'}), ${threads}`, 'ok');
  } else {
    log('Cross-origin isolation is not active. The page should reload once to enable it. If this warning stays, this browser cannot run the multithreaded engine.', 'warn');
  }
}

// ── Database bytes (IndexedDB cache or network) ──
//
// Index files are fetched same-origin first (local dev, and any deployment
// that still ships databases/), then from Zenodo — free, DOI-cited hosting
// with CORS-enabled downloads (the /api/records/…/files/<name>/content form;
// the /records/…/files/ page path lacks CORS). Verified with MD5 against the
// local indexes. To publish new indexes, upload them to a new Zenodo record
// and update this URL; the production build then stops shipping databases/.
const INDEX_BASE = 'https://zenodo.org/api/records/22102687/files/';

const DB_SHORT_NAMES = {
  resfinder: 'kma_index_resfinder_2_6_0',
  card_homolog: 'kma_index_card_4_0_1_homolog',
  vfdb_core: 'VFDB_setA_nt.fas.gz',
};

// Returns a fresh {ext: Uint8Array} map each call; the buffers are transferred
// to the run worker, so we always re-read from the IndexedDB cache rather than
// holding hundreds of MB resident on the main thread.

async function fetchDbFile(db, ext, onProgress) {
  const path = db.prefix + ext;
  let data = await getCachedDBFile(path);
  if (data && data.length) {
    try {
      // Cached bytes are re-verified once per session; a mismatch (e.g. from
      // an older, differently-built index cached under the same path) evicts
      // the entry and falls through to a fresh download.
      await verifyPinned(path, data);
      return { data, fromCache: true };
    } catch (err) {
      log(`Cached ${path} failed its integrity check. Discarding it and downloading again.`, 'warn');
      try { await deleteDBFile(path); } catch (_) { /* best-effort eviction */ }
    }
  }
  try {
    data = await fetchAssetWithProgress(path, 0, onProgress);
  } catch (err) {
    if (!INDEX_BASE) throw err;
    data = await fetchAssetWithProgress(INDEX_BASE + DB_SHORT_NAMES[db.shortName] + ext + '/content', 0, onProgress);
  }
  // Verify before use: a mismatch throws and the file is never cached, so a
  // corrupted or tampered download cannot change gene calls.
  await verifyPinned(path, data);
  if (data.length) await cacheDBFile(path, data); // never cache a failed (empty) load
  return { data, fromCache: false };
}

async function loadDatabaseBytes(db) {
  log(`Loading ${db.displayName}...`, 'info');
  const bytesByExt = {};
  for (const ext of DB_FILES) {
    const { data, fromCache } = await fetchDbFile(db, ext);
    log(`${db.prefix + ext}${fromCache ? ' (cached)' : ` (${formatBytes(data.length)})`}`, 'info');
    bytesByExt[ext] = data;
  }
  log(`${db.displayName} ready`, 'ok');
  return bytesByExt;
}

// ── Index preloading (↓ buttons) ──

export async function isDatabaseCached(key) {
  const db = DATABASES[key];
  if (!db) return false;
  for (const ext of DB_FILES) {
    const data = await getCachedDBFile(db.prefix + ext);
    if (!data || !data.length) return false;
  }
  return true;
}

// Downloads and caches every file of a database ahead of a run. onProgress is
// called with (bytesGot, bytesTotal) of the file currently transferring, and
// onFile with (index, total) as each file starts.
export async function preloadDatabase(key, onProgress, onFile) {
  const db = DATABASES[key];
  if (!db) throw new Error(`Unknown database: ${key}`);
  log(`Preloading ${db.displayName}...`, 'info');
  for (let i = 0; i < DB_FILES.length; i++) {
    onFile?.(i, DB_FILES.length);
    const { data, fromCache } = await fetchDbFile(db, DB_FILES[i], onProgress);
    log(`${db.prefix + DB_FILES[i]}${fromCache ? ' (cached)' : ` (${formatBytes(data.length)})`}`, 'info');
  }
  log(`${db.displayName} preloaded`, 'ok');
}

// ── Run single analysis ──

export async function runAnalysis(files, dbKey, config = {}) {
  if (!self.crossOriginIsolated) {
    throw new Error('Cross-origin isolation is not active. Reload the page. The multithreaded engine requires it.');
  }

  const db = DATABASES[dbKey] || DATABASES.resfinder;
  const dbBytes = await loadDatabaseBytes(db);

  const sample = detectSample(files);
  log(`Sample: ${sample.name} (${sample.type})`, 'info');

  // Plan the MEMFS layout the run worker will create + populate.
  const runDir = '/work/run';
  const inputDir = runDir + '/input';
  const mkdirs = ['/tmp', '/work', runDir, inputDir];
  // Index files write to /work/<prefix><ext>. FS.mkdir is not recursive, so
  // create every ancestor directory of the prefix in order.
  {
    const parts = db.prefix.split('/');
    parts.pop(); // drop the filename stem
    let acc = '/work';
    for (const p of parts) { acc += '/' + p; mkdirs.push(acc); }
  }

  // Database index bytes are already in memory → transfer them (zero-copy).
  const writes = [];
  const transfer = [];
  for (const [ext, data] of Object.entries(dbBytes)) {
    writes.push({ path: '/work/' + db.prefix + ext, data });
    transfer.push(data.buffer);
  }

  // Input FASTQs are passed as File objects (cheap structured clone) and streamed
  // into preallocated MEMFS files inside the worker, so a large sample is never
  // read into a full main-thread buffer and avoids MEMFS grow-and-copy spikes.
  const inputEntries = buildSafeInputEntries(files);
  const inputFiles = [];
  const inputPaths = [];
  for (const entry of inputEntries) {
    const path = inputDir + '/' + entry.safeName;
    inputFiles.push({ path, file: entry.file });
    inputPaths.push(path);
    log(`Queued ${entry.originalName || entry.safeName} (${formatBytes(entry.file.size)})`, 'info');
  }

  const threads = threadCount();
  const args = [];
  if (sample.type === 'paired') {
    args.push('-ipe', inputPaths.join(' '));
  } else {
    args.push('-i', inputPaths.join(' '));
  }
  args.push('-t_db', '/work/' + db.prefix, '-o', runDir + '/output');
  args.push('-ef');
  // Per-base depth (.mat.gz) and variant calls (.vcf.gz) for the in-page gene
  // viewer. Both are off by default in KMA; the consensus (.fsa) and alignment
  // (.aln) are already produced. Output is bounded by the length of the
  // detected templates, so this is cheap for AMR/isolate-scale runs.
  args.push('-matrix', '-vcf');
  args.push('-t', String(threads));
  args.push(...buildConfigArgs(config));
  const argsStr = args.join(' ');

  const outputs = OUTPUT_EXTS.map(ext => ({
    ext,
    path: runDir + '/output' + ext,
    binary: ext.endsWith('.gz'),
  }));

  log(`$ kma ${argsStr}`, 'cmd');
  log(`KMA alignment in progress (${threads} thread${threads === 1 ? '' : 's'})...`, 'progress');

  const t0 = performance.now();
  const result = await runInWorker({ mkdirs, writes, inputFiles, args: argsStr, outputs, engine: ENGINE, threads }, transfer);
  const elapsed = ((performance.now() - t0) / 1000).toFixed(2);

  if (result.exitCode === 0) {
    log(`Completed in ${elapsed}s`, 'ok');
  } else {
    log(`Failed (exit code ${result.exitCode}) after ${elapsed}s`, 'error');
  }

  return {
    sampleName: sample.name,
    sampleType: sample.type,
    database: db.shortName,
    dbLabel: db.displayName,
    exitCode: result.exitCode,
    elapsed,
    command: `kma ${argsStr}`,
    files: result.outputs,
    resTable: result.outputs['.res']?.data || null,
    inputFileNames: files.map(f => f.name),
  };
}

// Spawn a fresh worker for one run, stream its log lines, and tear it down when
// done — this is what frees the run's (possibly multi-GB) memory.
function runInWorker(job, transfer) {
  return new Promise((resolve, reject) => {
    let worker;
    try {
      worker = new Worker(new URL('./kma-runner.worker.js?v=wasm32a', import.meta.url));
    } catch (e) {
      reject(new Error('Failed to start KMA worker: ' + (e.message || e)));
      return;
    }
    worker.onmessage = (e) => {
      const msg = e.data;
      if (!msg) return;
      if (msg.type === 'log') {
        log(msg.text, msg.level);
      } else if (msg.type === 'done') {
        worker.terminate();
        resolve({ exitCode: msg.exitCode, outputs: msg.outputs });
      } else if (msg.type === 'error') {
        worker.terminate();
        reject(new Error(msg.message));
      }
    };
    worker.onerror = (e) => {
      worker.terminate();
      reject(new Error(e.message || 'KMA worker crashed'));
    };
    worker.postMessage({ type: 'run', ...job }, transfer);
  });
}

// ── Helpers ──

function detectSample(files) {
  const RE = /[_.](?:R?([12]))(?=\.(fastq|fq)(\.gz)?$)/i;
  if (files.length === 2) {
    const m1 = files[0].name.match(RE);
    const m2 = files[1].name.match(RE);
    if (m1 && m2 && m1[1] !== m2[1]) {
      const base = files[0].name.slice(0, m1.index);
      return { name: base || files[0].name, type: 'paired' };
    }
  }
  const name = files[0].name.replace(/\.(fastq|fq)(\.gz)?$/i, '');
  return { name, type: files.length === 2 ? 'paired' : 'single' };
}

function buildSafeInputEntries(files) {
  const used = new Set();
  return files.map((file, idx) => {
    const originalName = file.name || '';
    const safeName = makeUniqueSafeName(originalName, idx, used);
    return { file, originalName, safeName };
  });
}

function makeUniqueSafeName(name, idx, used) {
  const raw = (name || '').trim();
  const extMatch = raw.match(/(\.[A-Za-z0-9]{1,8}(?:\.[A-Za-z0-9]{1,8})?)$/);
  const ext = extMatch ? extMatch[1] : '';
  const baseRaw = ext ? raw.slice(0, -ext.length) : raw;
  const base = sanitizePathToken(baseRaw) || `input_${idx + 1}`;
  const safeExt = sanitizePathToken(ext);
  let candidate = `${base}${safeExt}`;
  if (!candidate) candidate = `input_${idx + 1}`;
  let n = 1;
  while (used.has(candidate)) {
    candidate = `${base}_${idx + 1}_${n}${safeExt}`;
    n += 1;
  }
  used.add(candidate);
  return candidate;
}

function sanitizePathToken(value) {
  return String(value).replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^_+|_+$/g, '');
}

function buildConfigArgs(config) {
  const args = [];
  if (config.nanopore) args.push('-bcNano', '-mp', '20');
  if (config.id_threshold != null) args.push('-ID', String(config.id_threshold * 100));
  if (config.mrc != null) args.push('-mrc', String(config.mrc));
  return args;
}

function formatBytes(b) {
  if (b < 1024) return b + ' B';
  if (b < 1048576) return (b / 1024).toFixed(1) + ' KB';
  if (b < 1073741824) return (b / 1048576).toFixed(1) + ' MB';
  return (b / 1073741824).toFixed(1) + ' GB';
}

export { formatBytes };
