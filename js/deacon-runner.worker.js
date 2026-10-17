/* deacon-runner.worker.js — removes human reads with Deacon in a Web Worker
 * (deacon/, see deacon/WASM_BUILD.md). One run per worker; terminating it
 * frees the index and the module's memory.
 *
 * A MODULE worker, unlike the KMA and fastp runners: the wasm-bindgen glue
 * is an ES module. Deacon is single-threaded and 32-bit, so every browser
 * runs the same build and there is no shared memory to set up.
 *
 * Job:   { type:'run', index: Uint8Array (transferred), files: [File, File?],
 *          names: [output name per file], absThreshold, relThreshold }
 * Posts: { type:'progress', done, total }   input bytes read so far
 *        { type:'done', reads: [File…], stats, info, seconds }
 *        { type:'error', message }
 *
 * The reads stream through in chunks: decompressed, classified and the
 * kept ones gzipped as they go, so only the output is held in memory.
 * Pairs are kept or removed together (a hit in either mate counts).
 */

import init, { WasmIndex, FilterSession, PairedFilterSession } from '../deacon/deacon_wasm.js?v=0.18.0';

const WASM_URL = new URL('../deacon/deacon_wasm_bg.wasm?v=0.18.0', import.meta.url);

const post = (msg) => self.postMessage(msg);

// Input compression is read from the bytes, not the name: gzip starts 1f 8b.
async function isGzip(file) {
  const head = new Uint8Array(await file.slice(0, 2).arrayBuffer());
  return head[0] === 0x1f && head[1] === 0x8b;
}

// Reports the bytes read at most every 250 ms.
function progressReporter(total) {
  let done = 0;
  let last = 0;
  return (n) => {
    done += n;
    const now = performance.now();
    if (now - last >= 250) {
      last = now;
      post({ type: 'progress', done, total });
    }
  };
}

async function filterSingle(index, file, name, opts, onBytes) {
  const session = new FilterSession(index, true, opts.absThreshold, opts.relThreshold,
    await isGzip(file), true, false, false);
  const out = [];
  const reader = file.stream().getReader();
  try {
    for (;;) {
      const { done, value } = await reader.read();
      if (done) break;
      onBytes(value.length);
      const chunk = session.push_chunk(value);
      if (chunk.length) out.push(chunk);
    }
    out.push(session.finish());
    return { reads: [new File(out, name, { type: 'application/gzip' })], stats: session.stats() };
  } finally {
    reader.releaseLock();
    session.free();
  }
}

// R1 and R2 are read in turn, one chunk each, so neither mate's parser runs
// far ahead of the other (records pair up in order).
async function filterPaired(index, files, names, opts, onBytes) {
  const [gz1, gz2] = await Promise.all(files.map(isGzip));
  const session = new PairedFilterSession(index, true, opts.absThreshold, opts.relThreshold,
    gz1, gz2, true, true, false, false);
  const out = [[], []];
  const keep = (o) => {
    if (o.r1.length) out[0].push(o.r1);
    if (o.r2.length) out[1].push(o.r2);
  };
  const readers = files.map(f => f.stream().getReader());
  const finished = [false, false];
  try {
    while (!finished[0] || !finished[1]) {
      for (const m of [0, 1]) {
        if (finished[m]) continue;
        const { done, value } = await readers[m].read();
        if (done) {
          finished[m] = true;
          keep(m === 0 ? session.finish_r1() : session.finish_r2());
        } else {
          onBytes(value.length);
          keep(m === 0 ? session.push_r1(value) : session.push_r2(value));
        }
      }
    }
    return {
      reads: out.map((chunks, m) => new File(chunks, names[m], { type: 'application/gzip' })),
      stats: session.stats(),
    };
  } finally {
    readers.forEach(r => r.releaseLock());
    session.free();
  }
}

self.onmessage = async (e) => {
  const job = e.data;
  if (!job || job.type !== 'run') return;
  try {
    await init({ module_or_path: WASM_URL });
    const index = new WasmIndex(job.index, undefined);
    job.index = null; // the module holds its own copy now
    const info = index.info();
    const opts = { absThreshold: job.absThreshold ?? 2, relThreshold: job.relThreshold ?? 0.01 };
    const total = job.files.reduce((n, f) => n + f.size, 0);
    const onBytes = progressReporter(total);
    const t0 = performance.now();
    const result = job.files.length === 2
      ? await filterPaired(index, job.files, job.names, opts, onBytes)
      : await filterSingle(index, job.files[0], job.names[0], opts, onBytes);
    post({ type: 'done', ...result, info, seconds: (performance.now() - t0) / 1000 });
  } catch (err) {
    post({ type: 'error', message: (err && err.message) || String(err) });
  }
};
