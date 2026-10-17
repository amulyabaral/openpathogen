/* host.js — human read removal before gene detection.
 *
 * Deacon (deacon/, see deacon/WASM_BUILD.md) removes reads that share
 * minimizers with the human pangenome index. The index downloads once,
 * like the gene databases: fetched same-origin (in 16 MiB parts on the
 * production host; see assets.js), checked against its pinned SHA-256
 * before use and cached in IndexedDB. The filtering runs in a worker; only
 * the kept reads come back, gzipped.
 */

import { getCachedDBFile, cacheDBFile, deleteDBFile } from './db.js';
import { fetchAssetWithProgress } from './assets.js';
import { verifyPinned } from './integrity.js';

export const HOST_INDEX = {
  path: 'databases/deacon/panhuman-1.k31w61.b16.pidx',
  name: 'panhuman-1',
  bytes: 250609707,
  size: '239 MB',
};

// A remote copy (e.g. a Zenodo /content URL) tried when the site does not
// serve the index itself; empty means same-origin only.
const HOST_INDEX_URL = '';

// Deacon's defaults: a read, or a pair counted together, is human when at
// least 2 of its distinct minimizers, and at least 1% of them, are in the
// index.
export const HOST_THRESHOLDS = { absThreshold: 2, relThreshold: 0.01 };

export async function isHostIndexCached() {
  const data = await getCachedDBFile(HOST_INDEX.path);
  return !!(data && data.length);
}

// The index bytes, from the cache or the network. onProgress(got, total)
// follows a download. log(msg, level) reports what happened.
export async function loadHostIndex(onProgress, log) {
  const path = HOST_INDEX.path;
  let data = await getCachedDBFile(path);
  if (data && data.length) {
    try {
      await verifyPinned(path, data);
      return { data, fromCache: true };
    } catch (_) {
      log?.('The cached human index failed its integrity check. Discarding it and downloading again.', 'warn');
      try { await deleteDBFile(path); } catch (_) { /* best-effort eviction */ }
    }
  }
  try {
    data = await fetchAssetWithProgress(path, HOST_INDEX.bytes, onProgress);
  } catch (err) {
    if (!HOST_INDEX_URL) throw err;
    data = await fetchAssetWithProgress(HOST_INDEX_URL, HOST_INDEX.bytes, onProgress);
  }
  // A mismatch throws before anything is cached or used.
  await verifyPinned(path, data);
  try {
    await cacheDBFile(path, data);
  } catch (_) {
    log?.('Could not cache the human index in this browser (storage full or disabled); it will download again next time.', 'warn');
  }
  return { data, fromCache: false };
}

// Remove human reads from one FASTQ or a pair. names are the output file
// names. Resolves { reads: [File…], stats: {readsIn, readsOut, basesIn,
// basesOut}, info, seconds }. onProgress(done, total) counts input bytes.
export async function removeHumanReads(files, { names, onProgress, log } = {}) {
  const { data } = await loadHostIndex(null, log);
  return new Promise((resolve, reject) => {
    let worker;
    try {
      worker = new Worker(new URL('./deacon-runner.worker.js?v=0.18.0', import.meta.url), { type: 'module' });
    } catch (e) {
      reject(new Error('Could not start the Deacon worker: ' + (e.message || e)));
      return;
    }
    worker.onmessage = (e) => {
      const m = e.data;
      if (m.type === 'progress') {
        onProgress?.(m.done, m.total);
      } else if (m.type === 'done') {
        worker.terminate();
        resolve(m);
      } else if (m.type === 'error') {
        worker.terminate();
        reject(new Error(m.message));
      }
    };
    worker.onerror = (e) => {
      worker.terminate();
      reject(new Error(e.message || 'Deacon worker crashed'));
    };
    worker.postMessage({ type: 'run', index: data, files, names, ...HOST_THRESHOLDS }, [data.buffer]);
  });
}
