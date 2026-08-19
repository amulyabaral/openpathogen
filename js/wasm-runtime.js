import { getCachedDBFile, cacheDBFile } from './db.js';

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

export function getDatabase(key) {
  return DATABASES[key] || DATABASES.resfinder;
}

export function getDatabaseList() {
  return Object.entries(DATABASES).map(([key, db]) => ({ key, ...db }));
}

// Number of KMA worker threads to request. Bounded by the Emscripten pthread
// pool baked into kma.js (Math.max(hardwareConcurrency||4, 4)), so this is
// always <= pool size. Threads require SharedArrayBuffer + cross-origin
// isolation; without it the (shared-memory) module cannot even load.
function threadCount() {
  if (!self.crossOriginIsolated) return 1;
  return navigator.hardwareConcurrency || 4;
}

// ── Init / readiness ──
//
// Memory64 + pthreads need cross-origin isolation (SharedArrayBuffer). The
// coi-serviceworker shim (loaded first in index.html) establishes it, reloading
// once on first visit. The KMA module itself is loaded per-run inside a Web
// Worker, so there is nothing heavy to pre-initialise here.

export async function initWasm() {
  if (typeof Worker === 'undefined') {
    throw new Error('Web Workers are unavailable in this browser.');
  }
  if (self.crossOriginIsolated) {
    log(`Ready. Memory64 build, up to ${threadCount()} threads`, 'ok');
  } else {
    log('Cross-origin isolation is not active. The page should reload once to enable it; if this warning persists, multithreading/large-sample support is unavailable in this browser.', 'warn');
  }
}

export function isReady() {
  return typeof Worker !== 'undefined' && !!self.crossOriginIsolated;
}

// ── Database bytes (IndexedDB cache or network) ──
//
// Returns a fresh {ext: Uint8Array} map each call; the buffers are transferred
// to the run worker, so we always re-read from the IndexedDB cache rather than
// holding hundreds of MB resident on the main thread.

async function loadDatabaseBytes(db) {
  log(`Loading ${db.displayName}...`, 'info');
  const bytesByExt = {};
  for (const ext of DB_FILES) {
    const cacheKey = db.prefix + ext;
    let data = await getCachedDBFile(cacheKey);
    if (data) {
      log(`${cacheKey} (cached)`, 'info');
    } else {
      const url = db.prefix + ext;
      const resp = await fetch(url);
      if (!resp.ok) throw new Error(`Failed to fetch ${url}: HTTP ${resp.status}`);
      data = new Uint8Array(await resp.arrayBuffer());
      log(`${cacheKey} (${formatBytes(data.length)})`, 'info');
      await cacheDBFile(cacheKey, data);
    }
    bytesByExt[ext] = data;
  }
  log(`${db.displayName} ready`, 'ok');
  return bytesByExt;
}

// ── Run single analysis ──

export async function runAnalysis(files, dbKey, config = {}) {
  if (!self.crossOriginIsolated) {
    throw new Error('Cross-origin isolation is not active. Reload the page. Multithreading and Memory64 require it.');
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
  const result = await runInWorker({ mkdirs, writes, inputFiles, args: argsStr, outputs }, transfer);
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
      worker = new Worker(new URL('./kma-runner.worker.js', import.meta.url));
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
