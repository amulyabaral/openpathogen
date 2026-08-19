/* wasm_ringpipe.h — streaming replacement for KMA's cross-stage pipe under
 * Emscripten. See wasm_ringpipe.c for the rationale. Only compiled for WASM. */
#ifndef WASM_RINGPIPE_H
#define WASM_RINGPIPE_H
#ifdef __EMSCRIPTEN__

#include <stdio.h>

/* kmaPipe-compatible entry point backed by a shared-memory ring buffer instead
 * of the OOM-prone tmpfile used by kmaPipeWasm. Bound to the global kmaPipe
 * pointer for the WASM build (see kmapipe.c). */
FILE * kmaPipeRing(const char *cmd, const char *type, FILE *ioStream, int *status);

#endif /* __EMSCRIPTEN__ */
#endif /* WASM_RINGPIPE_H */
