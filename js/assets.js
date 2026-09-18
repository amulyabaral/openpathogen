/* assets.js — fetch helper for static assets that may exceed Cloudflare
 * Pages' 25 MiB per-file cap (the VFDB index). The production build ships
 * such files as `<name>.part00`, `<name>.part01`, … ; this helper fetches
 * the parts and reassembles them in order, so nothing elsewhere in the app
 * needs to know about the cap.
 *
 * The plain URL is tried first, so local dev against the repo root
 * (unsplit files) keeps working unchanged.
 */

// Cloudflare Pages serves SPA fallback (index.html, HTTP 200, text/html)
// for unknown paths, so a 200 alone does not mean we got the real file.
function respHasAsset(resp) {
  if (!resp.ok) return false;
  const type = resp.headers.get('content-type') || '';
  return !type.includes('text/html');
}

function partPath(path, i) {
  return path + '.part' + String(i).padStart(2, '0'); // matches split -d -a 2
}

function concat(parts) {
  // tolerate both Uint8Array and ArrayBuffer inputs — ArrayBuffer has no
  // .length (only .byteLength), which silently produces an empty result
  const bufs = parts.map(p => (p instanceof Uint8Array ? p : new Uint8Array(p)));
  const total = bufs.reduce((n, p) => n + p.length, 0);
  const out = new Uint8Array(total);
  let offset = 0;
  for (const p of bufs) { out.set(p, offset); offset += p.length; }
  return out;
}

async function readStreamed(resp, total, onProgress) {
  if (!resp.body || !resp.body.getReader) {
    const buf = new Uint8Array(await resp.arrayBuffer());
    onProgress?.(buf.length, total || buf.length);
    return buf;
  }
  const reader = resp.body.getReader();
  const chunks = [];
  let received = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
    received += value.length;
    onProgress?.(received, total);
  }
  const out = new Uint8Array(received);
  let offset = 0;
  for (const c of chunks) { out.set(c, offset); offset += c.length; }
  return out;
}

export async function fetchAssetWithProgress(path, expectedBytes, onProgress) {
  const single = await fetch(path).catch(() => null);
  if (single && respHasAsset(single)) {
    const total = Number(single.headers.get('content-length')) || expectedBytes || 0;
    return readStreamed(single, total, onProgress);
  }

  // Chunked fallback: progress accumulates across parts against the
  // reassembled size (expectedBytes), not each part's own length.
  const bufs = [];
  let received = 0;
  for (let i = 0; i < 100; i++) {
    const resp = await fetch(partPath(path, i)).catch(() => null);
    if (!resp || !respHasAsset(resp)) break;
    const partTotal = Number(resp.headers.get('content-length')) || 0;
    const buf = await readStreamed(resp, partTotal,
      (got) => onProgress?.(received + got, expectedBytes));
    received += buf.length;
    bufs.push(buf);
  }
  if (!bufs.length) {
    throw new Error(`${path}: not found (HTTP ${single ? single.status : 'network error'})`);
  }
  const out = concat(bufs);
  onProgress?.(out.length, out.length);
  return out;
}
