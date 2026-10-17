/* sw.js — service worker with two jobs.
 *
 * 1. Offline use. The page, its scripts and the WebAssembly engines are
 *    cached, so a browser that was online once still loads the app and runs
 *    an analysis offline (the database indexes already live in IndexedDB).
 *    Same-origin requests are network-first: online, every file is fetched
 *    and revalidated as before; offline, the cached copy is served. Nothing
 *    cross-origin (ENA, Zenodo) is touched, and the database indexes are
 *    left out: the app keeps those in IndexedDB.
 *
 * 2. Cross-origin isolation on hosts that do not send the headers. The
 *    production host sends COOP/COEP itself (see _headers); on any other
 *    static host this worker adds them to each response and app.js reloads
 *    the page once so it is served through the worker.
 *
 * It also serves the fastp HTML report at reports/<id>.html with a policy
 * of its own, so the report's inline scripts and the vendored Plotly
 * (js/vendor/plotly-1.2.0.min.js, which app.js substitutes for the copy
 * fastp would fetch from opengene.org) can run without loosening the app's
 * Content-Security-Policy.
 */

const VERSION = '2026-10-16';
const CACHE = 'openpathogen-' + VERSION;
const REPORTS = 'openpathogen-reports';

// Fetched at install, so the page, every module and the engines are there
// before the first run, even if the first visit is the only one online.
// Every file under js/ belongs in this list (scripts/build-dist.sh checks).
const PRECACHE = [
  './', 'index.html', 'about.html', 'css/style.css', 'favicon.png', 'openpathogen_logo.png',
  'resfinder_db/phenotypes.txt', 'databases/cabbage/manifest.json',
  'js/aln-check.js', 'js/antibiotic-classes.js', 'js/app.js', 'js/assets.js', 'js/cabbage.js',
  'js/comprehensive.js', 'js/db.js', 'js/deacon-runner.worker.js', 'js/engine.js',
  'js/fastp-runner.worker.js', 'js/fetch-run.js', 'js/gene-view.js', 'js/genes.js', 'js/hamronize.js',
  'js/hits.js', 'js/host.js', 'js/integrity.js', 'js/kma-runner.worker.js', 'js/report.js',
  'js/resfinder-db.js', 'js/thresholds.js', 'js/wasm-runtime.js', 'js/vendor/fflate.min.js',
  'js/vendor/plotly-1.2.0.min.js',
  'kma/kma.js', 'kma/kma.wasm', 'kma/kma32.js', 'kma/kma32.wasm',
  'fastp/fastp.js', 'fastp/fastp.wasm', 'fastp/fastp32.js', 'fastp/fastp32.wasm',
  'deacon/deacon_wasm.js', 'deacon/deacon_wasm_bg.wasm',
];

// Same-origin paths kept for offline use (everything else same-origin is
// passed through untouched).
const CACHED = ['index.html', 'about.html', 'css/', 'js/', 'kma/', 'fastp/', 'deacon/', 'resfinder_db/',
  'databases/cabbage/manifest.json', 'favicon.png', 'openpathogen_logo.png', 'og.png'];

const REPORT_HEADERS = {
  'Content-Type': 'text/html; charset=utf-8',
  // Plotly 1.2.0 builds functions with eval at load, hence 'unsafe-eval';
  // the report can connect to nothing, so that cannot leak anything.
  'Content-Security-Policy': "default-src 'none'; script-src 'self' 'unsafe-inline' 'unsafe-eval'; "
    + "style-src 'unsafe-inline'; img-src data:; base-uri 'none'; form-action 'none'; frame-ancestors 'none'",
  'X-Content-Type-Options': 'nosniff',
  'Referrer-Policy': 'no-referrer',
  'Cache-Control': 'no-store',
};

const scopePath = new URL(self.registration.scope).pathname;

self.addEventListener('install', (event) => {
  event.waitUntil((async () => {
    const cache = await caches.open(CACHE);
    // Best effort: a file a host does not serve must not block installation.
    await Promise.allSettled(PRECACHE.map(p => cache.add(new Request(p, { cache: 'reload' }))));
    await self.skipWaiting();
  })());
});

self.addEventListener('activate', (event) => {
  event.waitUntil((async () => {
    for (const key of await caches.keys()) {
      if (key.startsWith('openpathogen-') && key !== CACHE && key !== REPORTS) await caches.delete(key);
    }
    await self.clients.claim();
  })());
});

self.addEventListener('fetch', (event) => {
  const req = event.request;
  if (req.method !== 'GET') return;
  const url = new URL(req.url);
  if (url.origin !== self.location.origin) return;
  if (req.cache === 'only-if-cached' && req.mode !== 'same-origin') return;
  if (!url.pathname.startsWith(scopePath)) return;
  const rel = url.pathname.slice(scopePath.length);
  if (rel.startsWith('reports/')) {
    event.respondWith(serveReport(req));
    return;
  }
  const cacheable = rel === '' || CACHED.some(p => rel === p || (p.endsWith('/') && rel.startsWith(p)));
  event.respondWith(networkFirst(req, cacheable));
});

// The cache key: the URL without its query, so app.js?v=66 and ?v=67 are
// one entry and the newest copy wins.
function keyOf(url) {
  const u = new URL(url);
  u.search = '';
  u.hash = '';
  return u.href;
}

async function networkFirst(req, cacheable) {
  try {
    const res = await fetch(req);
    if (cacheable && res.ok && res.status === 200 && res.type === 'basic' && !res.redirected) {
      const cache = await caches.open(CACHE);
      cache.put(keyOf(req.url), res.clone()).catch(() => {});
    }
    return isolate(res);
  } catch (err) {
    const cache = await caches.open(CACHE);
    let hit = await cache.match(keyOf(req.url));
    if (!hit && req.mode === 'navigate') hit = await cache.match(self.registration.scope);
    if (hit) return isolate(hit);
    throw err;
  }
}

// COOP/COEP on same-origin responses, for hosts that do not send them.
function isolate(res) {
  if (!res || res.status === 0 || res.type === 'opaque' || res.type === 'opaqueredirect' || res.redirected) return res;
  if (res.headers.get('Cross-Origin-Embedder-Policy') && res.headers.get('Cross-Origin-Opener-Policy')) return res;
  const headers = new Headers(res.headers);
  headers.set('Cross-Origin-Embedder-Policy', 'require-corp');
  headers.set('Cross-Origin-Opener-Policy', 'same-origin');
  return new Response(res.body, { status: res.status, statusText: res.statusText, headers });
}

// ── fastp HTML reports ──

self.addEventListener('message', (event) => {
  const data = event.data || {};
  const port = event.ports && event.ports[0];
  if (data.type !== 'store-report') return;
  event.waitUntil(storeReport(data.id, data.html)
    .then(ok => port?.postMessage({ ok }), () => port?.postMessage({ ok: false })));
});

async function storeReport(id, html) {
  if (typeof html !== 'string' || !/^[a-z0-9]{1,32}$/i.test(String(id))) return false;
  const cache = await caches.open(REPORTS);
  await cache.put(self.registration.scope + 'reports/' + id + '.html', new Response(html, { headers: REPORT_HEADERS }));
  const keys = await cache.keys(); // insertion order: drop all but the newest five
  for (const key of keys.slice(0, Math.max(0, keys.length - 5))) await cache.delete(key);
  return true;
}

async function serveReport(req) {
  const cache = await caches.open(REPORTS);
  const hit = await cache.match(keyOf(req.url));
  return hit || new Response('This report is no longer available. Open it again from the results page.',
    { status: 404, headers: { 'Content-Type': 'text/plain; charset=utf-8' } });
}
