// hv-runtime.js — thin JS wrapper around the hypervector profiler wasm module.
//
// The wasm module exposes plain C-ABI exports over linear memory (no
// bindgen): see hv/crates/hvwasm/src/lib.rs for the contract. All views into
// wasm memory are re-created after any call that may grow the heap.

export const HV_DEFAULTS = Object.freeze({
  k: 31,
  dimBits: 13, // D = 8192 → 32 kB fingerprint
  densNum: 5,
  densDen: 100,
  seed: 0x5eed0001 >>> 0 ? 0x5eed0001n : 0n,
  qmin: 20,
});

// .hvd header layout (little-endian) — mirrors hvcore::db::HvDb::to_bytes
export function parseHvdHeader(bytes) {
  const dv = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const magic = String.fromCharCode(bytes[0], bytes[1], bytes[2], bytes[3]);
  if (magic !== 'HVP1') throw new Error('not an .hvd file');
  return {
    version: dv.getUint32(4, true),
    k: dv.getUint32(8, true),
    dimBits: dv.getUint32(12, true),
    densNum: dv.getUint32(16, true),
    densDen: dv.getUint32(20, true),
    qmin: bytes[24],
    seed: dv.getBigUint64(25, true),
    nRows: dv.getUint32(33, true),
  };
}

// sidecar TSV: idx \t id \t name \t genus \t len \t l_sel \t flags
export function parseHvdMeta(text) {
  const rows = [];
  for (const line of text.split('\n')) {
    if (!line || line.startsWith('#')) continue;
    const f = line.split('\t');
    if (f.length < 7) continue;
    rows.push({
      idx: +f[0], id: f[1], name: f[2], genus: f[3],
      len: +f[4], lSel: +f[5], flags: +f[6],
    });
  }
  return rows;
}

export class HvModule {
  constructor({ exports, memory }) {
    this.e = exports;
    this.mem = memory;
  }

  viewBytes() {
    return new Uint8Array(this.mem.buffer);
  }

  reset(params) {
    return this.e.hv_encoder_reset(
      params.k, params.dimBits, params.densNum, params.densDen,
      params.seed, params.qmin,
    );
  }

  /** Feed one chunk of FASTQ bytes. Returns cumulative event count. */
  feed(chunk, eof) {
    let off = 0;
    let events = 0;
    while (off < chunk.length || (eof && off === chunk.length)) {
      const cap = this.e.hv_stage_cap();
      const ptr = this.e.hv_stage_ptr();
      const take = Math.min(cap, chunk.length - off);
      if (take > 0) {
        new Uint8Array(this.mem.buffer, ptr, take).set(chunk.subarray(off, off + take));
        off += take;
      }
      events = this.e.hv_encoder_feed(take, eof && off >= chunk.length ? 1 : 0);
      if (eof && off >= chunk.length) break;
      if (take === 0) break; // empty chunk, non-eof: nothing to do
    }
    return events;
  }

  get stats() {
    return {
      bases: this.e.hv_stats_bases(),
      reads: this.e.hv_stats_reads(),
      events: this.e.hv_encoder_feed(0, 0),
    };
  }

  loadDb(hvd, tsvText) {
    const enc = tsvText ? new TextEncoder().encode(tsvText) : null;
    const tsvLen = enc ? enc.length : 0;
    const ptr = this.e.hv_db_stage(hvd.length + tsvLen);
    new Uint8Array(this.mem.buffer, ptr, hvd.length).set(hvd);
    if (enc) {
      new Uint8Array(this.mem.buffer, ptr + hvd.length, tsvLen).set(enc);
    }
    return this.e.hv_db_load(hvd.length, tsvLen);
  }

  profile(zmin, kmax) {
    const n = this.e.hv_profile(zmin, kmax);
    if (n < 0) throw new Error(`hv_profile failed: ${n}`);
    const len = this.e.hv_results_len();
    const ptr = this.e.hv_results_ptr();
    const dv = new DataView(this.mem.buffer, ptr, len * 8);
    const header = {
      sigma: dv.getFloat64(0, true),
      unexplained: dv.getFloat64(8, true),
      spikeZ: dv.getFloat64(16, true),
      screenMax: dv.getFloat64(24, true),
      evBound: dv.getFloat64(32, true),
      nEvents: dv.getFloat64(40, true),
      rounds: dv.getFloat64(48, true),
    };
    const taxa = [];
    for (let i = 0; i < n; i++) {
      const o = 56 + i * 48;
      taxa.push({
        idx: dv.getFloat64(o, true),
        coverage: dv.getFloat64(o + 8, true),
        z: dv.getFloat64(o + 16, true),
        q: dv.getFloat64(o + 24, true),
        dnaFrac: dv.getFloat64(o + 32, true),
        flags: dv.getFloat64(o + 40, true),
      });
    }
    return { header, taxa };
  }

  fingerprintBytes(dpEps) {
    const len = this.e.hv_fp_serialize(dpEps || 0);
    const ptr = this.e.hv_fp_ptr();
    return new Uint8Array(this.mem.buffer.slice(ptr, ptr + len));
  }

  /** Copy of the raw int32 accumulator (for debugging/parity checks). */
  accI32() {
    const ptr = this.e.hv_acc_ptr();
    const d = 1 << 13;
    return new Int32Array(this.mem.buffer.slice(ptr, ptr + d * 4));
  }
}

export async function loadHvWasm(url) {
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`${url}: HTTP ${resp.status}`);
  const bytes = await resp.arrayBuffer();
  const { instance } = await WebAssembly.instantiate(bytes, { env: {} });
  return new HvModule({ exports: instance.exports, memory: instance.exports.memory });
}
