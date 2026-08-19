// IndexedDB cache for KMA database index files only.
const DB_NAME = 'dart';
const DB_VERSION = 2;

let dbPromise = null;

function openDB() {
  if (dbPromise) return dbPromise;
  dbPromise = new Promise((resolve, reject) => {
    const req = indexedDB.open(DB_NAME, DB_VERSION);
    req.onupgradeneeded = (e) => {
      const db = e.target.result;
      // Drop legacy stores from previous versions if present.
      for (const old of ['runs', 'files']) {
        if (db.objectStoreNames.contains(old)) db.deleteObjectStore(old);
      }
      if (!db.objectStoreNames.contains('dbcache')) {
        db.createObjectStore('dbcache', { keyPath: 'key' });
      }
    };
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error);
  });
  return dbPromise;
}

function reqP(req) {
  return new Promise((resolve, reject) => {
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error);
  });
}

async function store(mode = 'readonly') {
  const db = await openDB();
  return db.transaction('dbcache', mode).objectStore('dbcache');
}

export async function cacheDBFile(key, data) {
  return reqP((await store('readwrite')).put({ key, data }));
}

export async function getCachedDBFile(key) {
  const result = await reqP((await store()).get(key));
  return result?.data || null;
}
