/* engine.js — which WebAssembly build of KMA and fastp this browser runs.
 *
 * wasm64 (Memory64) is the main build: Chrome 133+, Firefox 134+. Browsers
 * without Memory64 — every iOS browser (all are WebKit) and Safari as of
 * 27.0 — get the wasm32 build of the same sources (kma/kma32.*,
 * fastp/fastp32.*), which addresses at most 4 GB and runs fewer threads.
 * Both need cross-origin isolation for their threads.
 *
 * ?engine=wasm32 forces the fallback build, for testing it on desktop.
 */

function hasMemory64() {
  try {
    // A module that declares one 64-bit memory (flags 0x04), and the memory
    // descriptor the Emscripten glue itself uses.
    const mod = new Uint8Array([0, 97, 115, 109, 1, 0, 0, 0, 5, 3, 1, 4, 1]);
    if (!WebAssembly.validate(mod)) return false;
    new WebAssembly.Memory({ address: 'i64', initial: 1n, maximum: 1n });
    return true;
  } catch (_) {
    return false;
  }
}

function forced() {
  try {
    return new URLSearchParams(self.location.search).get('engine');
  } catch (_) {
    return null;
  }
}

export const ENGINE_FORCED = forced() === 'wasm32';
export const ENGINE = ENGINE_FORCED || !hasMemory64() ? 'wasm32' : 'wasm64';

// KMA threads per run. The wasm32 build caps them at 4: it is what phones
// and Safari run, and each thread is a worker holding its own engine instance.
export function threadCount() {
  if (!self.crossOriginIsolated) return 1;
  const cores = navigator.hardwareConcurrency || 4;
  return ENGINE === 'wasm32' ? Math.min(cores, 4) : cores;
}

// Above this much input the wasm32 build may run out of memory (it has at
// most 4 GB, and a phone's browser tab much less). A warning, not a limit.
export const WASM32_INPUT_WARN_BYTES = 400 * 1024 * 1024;
