/* hv-worker.js — streaming hypervector profiler in a Web Worker.
 *
 * Receives { type: 'run', files: [File, ...], dbUrl } and:
 *   1. instantiates hvprof.wasm (hand-rolled C-ABI module)
 *   2. fetches (or reuses from IndexedDB cache) the int8 reference matrix
 *   3. streams each FASTQ through the encoder: File.stream() →
 *      DecompressionStream('gzip') when magic says gz → 4 MiB chunks into
 *      wasm staging. The file never resides in memory in full.
 *   4. profiles and posts { type: 'done', report, fingerprint }
 *
 * The FASTQ bytes never leave the worker; only the 32 kB fingerprint,
 * the packed result doubles, and progress messages cross the boundary.
 */

import { loadHvWasm, parseHvdHeader, parseHvdMeta } from './hv-runtime.js';
import { cacheDBFile, getCachedDBFile } from './db.js';

const WASM_URL = new URL('../hvprof/hvprof.wasm', self.location.href).href;
const DB_KEY = 'hv:db:patho_v1';

const post = (msg, transfer) => self.postMessage(msg, transfer || []);

function gunzipStream(file) {
  // DecompressionStream cannot un-read the magic bytes we sniffed, so we
  // re-inject them via a small identity TransformStream.
  let sniffed = null;
  const inject = new TransformStream({
    start(controller) {
      if (sniffed) controller.enqueue(sniffed);
    },
    transform(chunk, controller) {
      if (!sniffed) {
        sniffed = chunk.subarray(0, 2);
        if (sniffed[0] === 0x1f && sniffed[1] === 0x8b) {
          controller.enqueue(new Uint8Array([0x1f, 0x8b]));
        }
        const rest = chunk.subarray(2);
        if (rest.length) controller.enqueue(rest);
        return;
      }
      controller.enqueue(chunk);
    },
  });
  return file.stream().pipeThrough(inject).pipeThrough(new DecompressionStream('gzip'));
}

async function fetchDb(dbUrl, tsvUrl, log) {
  const cached = await getCachedDBFile(DB_KEY);
  if (cached) {
    log(`reference matrix cached (${(cached.byteLength / 1e6).toFixed(2)} MB)`);
    return cached;
  }
  log(`downloading reference matrix…`);
  const resp = await fetch(dbUrl);
  if (!resp.ok) throw new Error(`db: HTTP ${resp.status}`);
  const buf = new Uint8Array(await resp.arrayBuffer());
  try {
    await cacheDBFile(DB_KEY, buf);
  } catch (_) { /* cache is best-effort */ }
  log(`reference matrix ${(buf.byteLength / 1e6).toFixed(2)} MB (cached for next run)`);
  return buf;
}

async function fetchText(url) {
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`${url}: HTTP ${resp.status}`);
  return resp.text();
}

self.onmessage = async (e) => {
  const job = e.data;
  if (!job || job.type !== 'run') return;
  const log = (t) => post({ type: 'log', text: t });
  const t0 = performance.now();
  try {
    const hv = await loadHvWasm(WASM_URL);

    const dbBytes = await fetchDb(job.dbUrl, job.tsvUrl, log);
    const tsvText = await fetchText(job.tsvUrl);
    const header = parseHvdHeader(dbBytes);
    const meta = parseHvdMeta(tsvText);
    log(`matrix: ${header.nRows} genomes × D=${1 << header.dimBits} (k=${header.k}, ρ=${(header.densNum / header.densDen).toFixed(3)})`);

    // engine first (it owns the DB), params read from the DB header so the
    // encoder and the matrix cannot disagree
    hv.reset({
      k: header.k,
      dimBits: header.dimBits,
      densNum: header.densNum,
      densDen: header.densDen,
      seed: header.seed,
      qmin: header.qmin,
    });
    const n = hv.loadDb(dbBytes, tsvText);
    if (n < 0) throw new Error(`db load failed (${n})`);

    let totalBytes = 0;
    for (const file of job.files) {
      log(`streaming ${file.name} (${(file.size / 1e6).toFixed(1)} MB)…`);
      const gzipped = /\.gz$/i.test(file.name);
      let stream = gzipped ? gunzipStream(file) : file.stream();
      const reader = stream.getReader();
      const encT0 = performance.now();
      for (;;) {
        const { done, value } = await reader.read();
        if (done) break;
        if (value && value.length) {
          hv.feed(value, false);
          totalBytes += value.length;
          const dt = (performance.now() - encT0) / 1000;
          post({
            type: 'progress',
            bytes: totalBytes,
            events: hv.stats.events,
            reads: hv.stats.reads,
            mbps: totalBytes / 1e6 / Math.max(dt, 0.001),
          });
        }
      }
    }
    hv.feed(new Uint8Array(0), true); // EOF

    const { header: rep, taxa } = hv.profile(5.0, 96);
    const fingerprint = hv.fingerprintBytes(0);
    const secs = (performance.now() - t0) / 1000;
    post({
      type: 'done',
      secs,
      report: rep,
      taxa: taxa.map((t) => ({ ...t, meta: meta[t.idx] || null })),
      fingerprint,
      totalBytes,
    }, [fingerprint.buffer]);
  } catch (err) {
    post({ type: 'error', message: (err && err.message) || String(err) });
  }
};
