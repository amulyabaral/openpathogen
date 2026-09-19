/* kma-runner.worker.js — runs ONE KMA analysis in an isolated Web Worker.
 *
 * Why a worker:
 *   - KMA is multithreaded (pthreads). pthread_join() blocks the calling
 *     thread; you cannot block the main browser thread, so KMA must run here.
 *   - A fresh worker is spawned per run and terminated afterwards. That gives
 *     clean KMA globals every run and frees grown (multi-GB, Memory64) linear
 *     memory back to the OS — replacing the old whole-memory snapshot/restore.
 *
 * This is loaded as a CLASSIC worker (not type:"module") so it can
 * importScripts() the Emscripten glue, which is itself a classic script.
 */

// ?v= busts caches when the build changes (old glue + new wasm must not mix)
const KMA_BUILD = 'wasm32a';
const KMA_DIR_URL = new URL('../kma/', self.location.href).href;

let modulePromise = null;

// A shared memory reserves its maximum up front, and phones refuse large
// reservations, so the wasm32 build gets the largest maximum the device
// accepts (the module's own ceiling is 4 GB). Returns [memory, MB].
function createWasm32Memory() {
  const initial = 67108864 / 65536; // INITIAL_MEMORY in Makefile.wasm32
  for (const mb of [4096, 2048, 1024, 512]) {
    try {
      return [new WebAssembly.Memory({ initial, maximum: mb * 16, shared: true }), mb];
    } catch (_) { /* try a smaller reservation */ }
  }
  return [null, 0];
}

// Load + initialise the Emscripten module once in this worker. engine is
// 'wasm64' (kma.js) or 'wasm32' (kma32.js, for browsers without Memory64).
function loadModule(engine, threads) {
  if (modulePromise) return modulePromise;
  const jsName = engine === 'wasm32' ? 'kma32.js' : 'kma.js';
  const jsUrl = KMA_DIR_URL + jsName + '?v=' + KMA_BUILD;
  modulePromise = new Promise((resolve, reject) => {
    const Module = {
      noExitRuntime: true,
      // KMA spawns its pthread workers via `new Worker(mainScriptUrlOrBlob)`.
      // We must point that at the glue itself — importScripts() leaves
      // self.location (and Emscripten's _scriptName) pointing at THIS runner
      // script, which would spawn the wrong file as a pthread.
      mainScriptUrlOrBlob: jsUrl,
      locateFile(path) { return KMA_DIR_URL + path + '?v=' + KMA_BUILD; },
      print(text) { post({ type: 'log', level: 'stdout', text }); },
      printErr(text) { post({ type: 'log', level: 'stderr', text }); },
      onRuntimeInitialized() { resolve(Module); },
      onAbort(reason) { reject(new Error('KMA module aborted: ' + reason)); },
    };
    if (engine === 'wasm32') {
      // Same pool rule as the wasm64 build: producer and consumer stages
      // may each run -t threads, plus the pipe threads.
      Module.kmaPoolSize = 2 * Math.max(threads || 4, 4) + 2;
      const [memory, mb] = createWasm32Memory();
      if (memory) {
        Module.wasmMemory = memory;
        post({ type: 'log', level: 'info', text: `32-bit engine, memory limit ${mb >= 1024 ? mb / 1024 + ' GB' : mb + ' MB'}` });
      }
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

function post(msg, transfer) { self.postMessage(msg, transfer || []); }

// Stream a File into a MEMFS path. The file's full size is preallocated up
// front (FS.ftruncate) so MEMFS makes ONE backing allocation instead of
// repeatedly growing+copying its Uint8Array (which peaks ~2x and OOMs on large
// inputs). Reads stay fully in-memory/fast — unlike WORKERFS, whose on-demand
// FileReaderSync reads are orders of magnitude slower on big files.
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
    const Module = await loadModule(job.engine, job.threads);
    const FS = Module.FS;

    // Stage the database + input files into MEMFS at the planned paths.
    for (const dir of job.mkdirs) {
      try { FS.mkdir(dir); } catch (_) { /* already exists */ }
    }
    // Database index bytes arrive already in memory (transferred buffers).
    for (const w of job.writes) {
      FS.writeFile(w.path, w.data);
    }
    // Input FASTQs are streamed from their File objects into preallocated MEMFS
    // files (see writeFileStreamed), so a large sample avoids the grow-and-copy
    // memory spike while keeping fast in-memory reads.
    for (const inp of job.inputFiles) {
      await writeFileStreamed(FS, inp.path, inp.file);
    }

    // Run KMA. kma_main returns its exit status; on exit() Emscripten throws
    // an ExitStatus carrying the same code.
    let exitCode = 0;
    try {
      exitCode = Module.ccall('kma_run', 'number', ['string'], [job.args]);
    } catch (err) {
      if (err && (err.name === 'ExitStatus' || typeof err.status === 'number')) {
        exitCode = err.status || 0;
      } else {
        throw err;
      }
    }

    // Read outputs back out of MEMFS, transferring binary buffers zero-copy.
    const outputs = {};
    const transfer = [];
    for (const o of job.outputs) {
      try {
        if (o.binary) {
          const data = FS.readFile(o.path); // Uint8Array
          if (data && data.length) {
            outputs[o.ext] = { data, binary: true };
            transfer.push(data.buffer);
          }
        } else {
          const data = FS.readFile(o.path, { encoding: 'utf8' });
          if (data && data.trim().length) {
            outputs[o.ext] = { data, binary: false };
          }
        }
      } catch (_) { /* output not produced for this run */ }
    }

    post({ type: 'done', exitCode, outputs }, transfer);
  } catch (err) {
    post({ type: 'error', message: (err && err.message) || String(err) });
  }
};
