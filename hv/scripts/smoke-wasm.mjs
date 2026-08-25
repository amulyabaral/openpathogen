#!/usr/bin/env node
// T10: native ↔ wasm concordance — the same FASTQ through the CLI encoder
// and the browser wasm module must produce bit-identical fingerprints.
//
//   node hv/scripts/smoke-wasm.mjs /tmp/hv/demo/sample.hvf /tmp/hv/demo/reads_1.fq
//                                 /tmp/hv/demo/db.hvd /tmp/hv/demo/db.tsv
//
// (use plain fastq: node has no DecompressionStream in older versions)

import { readFileSync } from 'node:fs';

const [nativeFpPath, fqPath, hvdPath, tsvPath] = process.argv.slice(2);
if (!nativeFpPath || !fqPath) {
  console.error('usage: smoke-wasm.mjs NATIVE.hvf READS.fq [DB.hvd DB.tsv]');
  process.exit(2);
}

const wasmBytes = readFileSync(new URL('../target/wasm32-unknown-unknown/release/hvwasm.wasm', import.meta.url));
const { instance } = await WebAssembly.instantiate(wasmBytes, { env: {} });
const e = instance.exports;

// --- encoder parity ---
e.hv_encoder_reset(31, 13, 5, 100, 0x5eed0001n, 20);
const fq = readFileSync(fqPath);
const chunk = 3 << 20;
let off = 0;
while (off < fq.length) {
  const take = Math.min(chunk, fq.length - off);
  const ptr = e.hv_stage_ptr();
  const cap = e.hv_stage_cap();
  const view = new Uint8Array(instance.exports.memory.buffer, ptr, take);
  view.set(fq.subarray(off, off + take));
  e.hv_encoder_feed(take, 0);
  off += take;
}
e.hv_encoder_feed(0, 1);
const fpLen = e.hv_fp_serialize(0);
const fpPtr = e.hv_fp_ptr();
const wasmFp = new Uint8Array(instance.exports.memory.buffer.slice(fpPtr, fpPtr + fpLen));

const nativeFp = readFileSync(nativeFpPath);
let same = nativeFp.length === wasmFp.length;
if (same) {
  for (let i = 0; i < nativeFp.length; i++) {
    if (nativeFp[i] !== wasmFp[i]) { same = false; break; }
  }
}
console.log(`fingerprint bytes: native=${nativeFp.length} wasm=${wasmFp.length} identical=${same}`);
if (!same) process.exit(1);

// --- profile parity (optional when DB given) ---
if (hvdPath && tsvPath) {
  const hvd = readFileSync(hvdPath);
  const tsv = readFileSync(tsvPath);
  const ptr = e.hv_db_stage(hvd.length + tsv.length);
  new Uint8Array(instance.exports.memory.buffer, ptr, hvd.length).set(hvd);
  new Uint8Array(instance.exports.memory.buffer, ptr + hvd.length, tsv.length).set(tsv);
  const n = e.hv_db_load(hvd.length, tsv.length);
  if (n < 0) { console.error(`db load failed: ${n}`); process.exit(1); }
  const taxa = e.hv_profile(5.0, 96);
  const len = e.hv_results_len();
  const rptr = e.hv_results_ptr();
  const dv = new DataView(instance.exports.memory.buffer, rptr, len * 8);
  console.log(`profile: ${taxa} taxa, unexplained=${dv.getFloat64(8, true).toFixed(4)}, ` +
    `spike_z=${dv.getFloat64(16, true).toFixed(2)}`);
  for (let i = 0; i < taxa; i++) {
    const o = 56 + i * 48;
    console.log(`  row ${dv.getFloat64(o, true)} cov=${dv.getFloat64(o + 8, true).toFixed(3)} ` +
      `z=${dv.getFloat64(o + 16, true).toFixed(1)} frac=${(100 * dv.getFloat64(o + 32, true)).toFixed(2)}%`);
  }
}
console.log('T10 OK');
