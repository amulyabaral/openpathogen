/* fastp-runner.worker.js — QC + adapter/quality trimming via fastp (WASM,
 * wasm64 + pthreads build — see fastp/WASM_BUILD.md). Loaded as a CLASSIC
 * worker so it can importScripts() the Emscripten glue, like the KMA runner.
 *
 * Job: { type: 'run', files: [File, File?], paired: bool, nanopore: bool,
 *        options?: object, engine?: 'wasm64' | 'wasm32' }
 * Posts:
 *   { type:'log', text }                     progress lines
 *   { type:'done', report, reads:[File...], html, paired }  on success
 *   { type:'error', message }                on failure
 *
 * options (all optional; the UI in index.html shows and pre-fills defaults):
 *   adapterTrim       false → --disable_adapter_trimming
 *   adapterR1/R2      custom adapter sequences (IUPAC); blank = auto
 *   qualityFilter     false → --disable_quality_filtering; otherwise filters
 *                     with qualifiedPhred (15), unqualifiedPercent (40),
 *                     nBaseLimit (5) — fastp's own defaults
 *   lengthFilter      false → --disable_length_filtering; otherwise filters
 *                     with minLength (15) and maxLength (0 = unlimited)
 *   polyG             'auto' (fastp decides) | 'on' | 'off'
 *   correction        paired-only overlap base correction (fastp: off)
 *
 * The whole sample is processed — the old 400k-read cap of the sequential
 * wasm32 build is gone. Filtered reads come back as gzip File objects (fastp
 * gzips outputs whose names end in .gz) so the (already battle-tested) KMA
 * runner can stream them straight into its own worker unchanged.
 */

// ?v= busts caches when the build changes (old glue + new wasm must not mix)
const FASTP_BUILD = 'wasm32a';
const FASTP_DIR_URL = new URL('../fastp/', self.location.href).href;

const post = (msg, transfer) => self.postMessage(msg, transfer || []);

let modulePromise = null;

// wasm32 build: a shared memory reserves its maximum up front and phones
// refuse large reservations, so take the largest maximum the device accepts
// (the module's own ceiling is 4 GB). Same rule as the KMA runner.
function createWasm32Memory() {
  const initial = 67108864 / 65536; // INITIAL_MEMORY of the fastp32 build
  for (const mb of [4096, 2048, 1024, 512]) {
    try {
      return new WebAssembly.Memory({ initial, maximum: mb * 16, shared: true });
    } catch (_) { /* try a smaller reservation */ }
  }
  return null;
}

// engine is 'wasm64' (fastp.js) or 'wasm32' (fastp32.js, for browsers
// without Memory64; see fastp/WASM_BUILD.md).
function loadModule(engine) {
  if (modulePromise) return modulePromise;
  const jsUrl = FASTP_DIR_URL + (engine === 'wasm32' ? 'fastp32.js' : 'fastp.js') + '?v=' + FASTP_BUILD;
  modulePromise = new Promise((resolve, reject) => {
    const Module = {
      noExitRuntime: true,
      locateFile(path) {
        return FASTP_DIR_URL + path + '?v=' + FASTP_BUILD;
      },
      // pthread workers must spawn from the glue itself, not this runner
      // script (self.location points at js/, like the KMA runner)
      mainScriptUrlOrBlob: jsUrl,
      print(text) { post({ type: 'log', text: String(text) }); },
      printErr(text) { post({ type: 'log', text: String(text) }); },
      onRuntimeInitialized() { resolve(Module); },
      onAbort(reason) { reject(new Error('fastp module aborted: ' + reason)); },
    };
    if (engine === 'wasm32') {
      const memory = createWasm32Memory();
      if (memory) Module.wasmMemory = memory;
    }
    self.Module = Module;
    try {
      importScripts(jsUrl);
    } catch (e) {
      reject(e);
    }
  });
  return modulePromise;
}

// Stream a File into MEMFS with preallocation (same rationale as the KMA
// runner: avoids the grow-and-copy memory spike on big inputs).
async function writeFileStreamed(FS, path, file) {
  const stream = FS.open(path, 'w');
  try {
    if (file.size > 0) FS.ftruncate(stream.fd, file.size);
    const reader = file.stream().getReader();
    let position = 0;
    for (;;) {
      const { done, value } = await reader.read();
      if (done) break;
      FS.write(stream, value, 0, value.length, position);
      position += value.length;
    }
  } finally {
    FS.close(stream);
  }
}

// Adapter sequences: 5–100 IUPAC bases. Returns the cleaned sequence, null
// for blank, or undefined when invalid (caller warns and falls back to auto).
const IUPAC_ADAPTER_RE = /^[ACGTURYSWKMBDHVN]{5,100}$/;
function cleanAdapter(seq) {
  const s = String(seq || '').replace(/\s+/g, '').toUpperCase();
  if (!s) return null;
  return IUPAC_ADAPTER_RE.test(s) ? s : undefined;
}

function clampInt(v, fallback, lo, hi) {
  const n = Math.round(Number(v));
  return Number.isFinite(n) ? Math.min(hi, Math.max(lo, n)) : fallback;
}

function buildArgs(inNames, outNames, options, { paired, nanopore }) {
  const o = options || {};
  const args = ['--in1', inNames[0]];
  if (paired) args.push('--in2', inNames[1]);
  args.push('--out1', outNames[0]);
  if (paired) args.push('--out2', outNames[1]);
  args.push(
    '--json', '/out/report.json', '--html', '/out/report.html',
    // Single-threaded on purpose: with >1 consumer thread fastp's filtered
    // output ordering (and thus its bytes) is nondeterministic, and the
    // proxied MEMFS writes serialise on the runtime thread anyway, so threads
    // buy no wall-clock win in the browser. --thread 1 keeps every run
    // byte-reproducible — the native-vs-WASM validation standard.
    '--thread', '1');

  if (nanopore) {
    // fastp is Illumina-oriented; for ONT we only want honest statistics,
    // not trimming — pass reads through untouched. User options don't apply.
    args.push('--disable_adapter_trimming', '--disable_quality_filtering',
              '--disable_trim_poly_g');
    return args;
  }

  // ── Adapter trimming ──
  if (o.adapterTrim === false) {
    args.push('--disable_adapter_trimming');
  } else {
    const a1 = cleanAdapter(o.adapterR1);
    if (a1 === undefined) {
      post({ type: 'log', text: 'Ignoring adapter R1: expected 5–100 IUPAC bases (A C G T U R Y S W K M B D H V N); using auto-detection instead.' });
    }
    const a2 = paired ? cleanAdapter(o.adapterR2) : null;
    if (a2 === undefined) {
      post({ type: 'log', text: 'Ignoring adapter R2: expected 5–100 IUPAC bases (A C G T U R Y S W K M B D H V N); using auto-detection instead.' });
    }
    if (a1) args.push('--adapter_sequence', a1);
    if (a2) args.push('--adapter_sequence_r2', a2);
    // An explicit sequence is an explicit choice — only auto-detect when the
    // user left both fields blank.
    if (paired && !a1 && !a2) args.push('--detect_adapter_for_pe');
  }

  // ── Quality filtering (on by default, with fastp's thresholds — the UI
  // pre-fills the same values; opt out to keep every trimmed read) ──
  if (o.qualityFilter === false) {
    args.push('--disable_quality_filtering');
  } else {
    args.push('--qualified_quality_phred', String(clampInt(o.qualifiedPhred, 15, 1, 40)),
              '--unqualified_percent_limit', String(clampInt(o.unqualifiedPercent, 40, 0, 100)),
              '--n_base_limit', String(clampInt(o.nBaseLimit, 5, 0, 50)));
  }

  // ── Length filtering (on by default, min 15 bp, no upper limit) ──
  if (o.lengthFilter === false) {
    args.push('--disable_length_filtering');
  } else {
    args.push('--length_required', String(clampInt(o.minLength, 15, 1, 1000)));
    const maxLen = clampInt(o.maxLength, 0, 0, 100000);
    if (maxLen > 0) args.push('--length_limit', String(maxLen));
  }

  // ── Poly-G tails (two-color chemistry: NextSeq/NovaSeq) ──
  if (o.polyG === 'on') args.push('--trim_poly_g');
  else if (o.polyG === 'off') args.push('--disable_trim_poly_g');

  // ── Overlap-based base correction (paired only) ──
  if (paired && o.correction) args.push('--correction');

  return args;
}

const sleep = ms => new Promise(r => setTimeout(r, ms));

// Run fastp via the wasm_wrapper entry point: fastp_run spawns a pthread that
// runs fastp's native threaded pipeline; the runtime thread must stay in this
// event loop to service Emscripten's proxied FS calls, so completion is
// observed by polling the wrapper's state machine.
async function runFastp(Module, args) {
  Module.ccall('fastp_run', 'number', ['string'], [args.join(' ')]);
  for (;;) {
    await sleep(100);
    const state = (() => {
      try { return Module.ccall('fastp_state', 'number'); } catch (_) { return 2; }
    })();
    if (state !== 2) continue;
    return Module.ccall('fastp_code', 'number');
  }
}

self.onmessage = async (e) => {
  const job = e.data;
  if (!job || job.type !== 'run') return;
  try {
    const Module = await loadModule(job.engine);
    const FS = Module.FS;
    for (const dir of ['/tmp', '/in', '/out']) {
      try { FS.mkdir(dir); } catch (_) { /* exists */ }
    }

    const paired = !!job.paired && job.files.length === 2;
    for (let i = 0; i < job.files.length; i++) {
      // fastp only gz-opens paths ending in .gz — match the staged name to
      // whether the incoming bytes are actually gzip
      const gz = /\.gz$/i.test(job.files[i].name || '') || /\.gz$/i.test(job.files[i].type || '');
      await writeFileStreamed(FS, `/in/r${i + 1}.fastq${gz ? '.gz' : ''}`, job.files[i]);
    }
    const inNames = job.files.map((f, i) =>
      `/in/r${i + 1}.fastq${/\.gz$/i.test(f.name || '') ? '.gz' : ''}`);
    // .gz names → fastp writes gzipped outputs, keeping full-sample output
    // memory manageable in MEMFS
    const outNames = paired ? ['/out/filtered_1.fastq.gz', '/out/filtered_2.fastq.gz']
                            : ['/out/filtered.fastq.gz'];

    const args = buildArgs(inNames, outNames, job.options, { paired, nanopore: !!job.nanopore });

    post({ type: 'log', text: `$ fastp ${args.join(' ')}` });
    const exitCode = await runFastp(Module, args);
    if (exitCode !== 0) throw new Error(`fastp exited with code ${exitCode}`);

    const reportRaw = FS.readFile('/out/report.json', { encoding: 'utf8' });
    const report = JSON.parse(reportRaw);

    const reads = [];
    const transfer = [];
    for (const name of outNames) {
      let data;
      try {
        data = FS.readFile('/out/' + name.split('/').pop());
      } catch (_) { continue; }
      if (!data || !data.length) continue;
      reads.push(new File([data], name.split('/').pop(), { type: 'application/gzip' }));
      transfer.push(data.buffer);
    }

    let html = null;
    try {
      html = FS.readFile('/out/report.html', { encoding: 'utf8' });
    } catch (_) { /* optional */ }

    post({ type: 'done', report, reads, html, paired }, transfer);
  } catch (err) {
    post({ type: 'error', message: (err && err.message) || String(err) });
  }
};
