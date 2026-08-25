/* fastp-runner.worker.js — QC + adapter/quality trimming via fastp (WASM,
 * built by biowasm, MIT). Loaded as a CLASSIC worker so it can
 * importScripts() the Emscripten glue, like the KMA runner.
 *
 * Job: { type: 'run', files: [File, File?], paired: bool, nanopore: bool }
 * Posts:
 *   { type:'log', text }                     progress lines
 *   { type:'done', report, reads:[File...], html? }  on success
 *   { type:'error', message }                on failure
 *
 * The filtered reads come back as File objects so the (already battle-tested)
 * KMA runner can stream them straight into its own MEMFS.
 */

const FASTP_JS_URL = new URL('../fastp/fastp.js', self.location.href).href;
const FASTP_DIR_URL = new URL('../fastp/', self.location.href).href;

const post = (msg, transfer) => self.postMessage(msg, transfer || []);

let modulePromise = null;

function loadModule() {
  if (modulePromise) return modulePromise;
  modulePromise = new Promise((resolve, reject) => {
    const Module = {
      noExitRuntime: true,
      locateFile(path) {
        return FASTP_DIR_URL + path;
      },
      print(text) { post({ type: 'log', text: String(text) }); },
      printErr(text) { post({ type: 'log', text: String(text) }); },
      onRuntimeInitialized() { resolve(Module); },
      onAbort(reason) { reject(new Error('fastp module aborted: ' + reason)); },
    };
    self.Module = Module;
    try {
      importScripts(FASTP_JS_URL);
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

self.onmessage = async (e) => {
  const job = e.data;
  if (!job || job.type !== 'run') return;
  try {
    const Module = await loadModule();
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

    const args = [];
    if (paired) {
      args.push('--in1', inNames[0], '--in2', inNames[1],
                '--out1', '/out/filtered_1.fastq', '--out2', '/out/filtered_2.fastq',
                '--detect_adapter_for_pe');
    } else {
      args.push('--in1', inNames[0], '--out1', '/out/filtered.fastq');
    }
    args.push('--json', '/out/report.json', '--html', '/out/report.html',
              '--thread', '1',
              // keep every read that survives trimming — KMA thresholds
              // govern gene calls, not host-side length filters
              '--disable_length_filtering', '--unqualified_percent_limit', '100',
              // cap below PACK_IN_MEM_LIMIT packs: required by the
              // sequential (no-threads) WASM build — see fastp/WASM_BUILD.md.
              // QC statistics on the first 400k reads are unchanged to
              // any practical precision.
              '--reads_to_process', '400000');
    if (job.nanopore) {
      // fastp is Illumina-oriented; for ONT we only want honest statistics,
      // not trimming — pass reads through untouched.
      args.push('--disable_adapter_trimming', '--disable_quality_filtering',
                '--disable_trim_poly_g');
    }

    post({ type: 'log', text: `$ fastp ${args.join(' ')}` });
    let exitCode = 0;
    try {
      exitCode = Module.callMain(args);
    } catch (err) {
      if (err && (err.name === 'ExitStatus' || typeof err.status === 'number')) {
        exitCode = err.status || 0;
      } else {
        throw err;
      }
    }
    if (exitCode !== 0) throw new Error(`fastp exited with code ${exitCode}`);

    const reportRaw = FS.readFile('/out/report.json', { encoding: 'utf8' });
    const report = JSON.parse(reportRaw);

    const reads = [];
    const transfer = [];
    const outNames = paired ? ['filtered_1.fastq', 'filtered_2.fastq']
                            : ['filtered.fastq'];
    for (const name of outNames) {
      let data;
      try {
        data = FS.readFile('/out/' + name);
      } catch (_) { continue; }
      if (!data || !data.length) continue;
      reads.push(new File([data], name, { type: 'text/plain' }));
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
