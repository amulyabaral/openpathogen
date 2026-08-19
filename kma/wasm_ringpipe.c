/* wasm_ringpipe.c — streaming cross-stage pipe for KMA under Emscripten.
 *
 * KMA pipelines its stages (read-conversion -> mapping/scoring) through a child
 * connected by a pipe. Native KMA runs that child in a *thread* (kmaPipeThread
 * is the native default) sharing the address space. The WASM port previously
 * used kmaPipeWasm, which ran stages sequentially and handed data through a
 * MEMFS tmpfile() — materialising the WHOLE intermediate in linear memory, which
 * OOMs ("Array buffer allocation failed") on large samples.
 *
 * Emscripten's own pipe() (PIPEFS) can't stream here: its reads are hard-coded
 * non-blocking and its buffers are per-worker JS objects proxied to the runtime
 * thread, so a blocking read would deadlock the proxied write.
 *
 * This file keeps the proven native threading model but replaces the transport
 * with a bounded ring buffer in WASM LINEAR MEMORY. Under -pthread the linear
 * memory is a SharedArrayBuffer, so the buffer is directly readable/writable
 * from both the producer pthread and the consumer runtime thread with no FS
 * proxy — and WebAssembly atomic wait/notify (memory.atomic.wait32 /
 * memory.atomic.notify) give real blocking with backpressure. The ring is
 * wrapped as read/write FILE* ends via fopencookie(), so KMA's fread/fwrite
 * pipe sites are unchanged.
 */
#ifdef __EMSCRIPTEN__

#define _GNU_SOURCE   /* fopencookie, cookie_io_functions_t */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <pthread.h>
#include <sys/types.h>
#include <stdint.h>
#include <emscripten/threading.h>
#include "kma.h"
#include "pherror.h"
#include "threader.h"
#include "wasm_ringpipe.h"

/* 16 MiB ring — large enough to smooth producer/consumer rate mismatch, tiny
 * next to the multi-hundred-MB intermediates it replaces. */
#define RING_CAP (16u * 1024u * 1024u)

typedef struct RingPipe {
	size_t head;            /* total bytes written (producer), monotonic */
	size_t tail;            /* total bytes read (consumer), monotonic */
	int closed;             /* write end closed (producer done) */
	int reader_closed;      /* read end closed (consumer gone) */
	int refs;               /* freed when both FILE* ends close */
	int32_t wseq;           /* futex word bumped on every write / close */
	int32_t rseq;           /* futex word bumped on every read / reader-close */
	size_t cap;
	unsigned char *buf;
} RingPipe;

typedef struct RingEnd { RingPipe *ring; int is_writer; } RingEnd;

/* ── atomic helpers (all memory is shared under -pthread) ── */
#define A_LOAD(p)      __atomic_load_n((p), __ATOMIC_ACQUIRE)
#define A_STORE(p, v)  __atomic_store_n((p), (v), __ATOMIC_RELEASE)
#define A_BUMP(p)      __atomic_add_fetch((p), 1, __ATOMIC_SEQ_CST)

static void futex_wait(int32_t *addr, int32_t expected) {
	/* Blocks while *addr == expected. -1 = wait forever. Returns immediately
	   (code 1) if the value already changed, so wakeups are never lost. */
	__builtin_wasm_memory_atomic_wait32(addr, expected, -1LL);
}
static void futex_wake(int32_t *addr) {
	__builtin_wasm_memory_atomic_notify(addr, 1);
}

/* ── ring lifecycle ── */
static RingPipe * ring_create(size_t cap) {
	RingPipe *r = calloc(1, sizeof(RingPipe));
	if(!r) return NULL;
	r->buf = malloc(cap);
	if(!r->buf) { free(r); return NULL; }
	r->cap = cap;
	r->refs = 2;
	return r;
}
static void ring_end_release(RingPipe *r) {
	if(__atomic_sub_fetch(&r->refs, 1, __ATOMIC_ACQ_REL) == 0) {
		free(r->buf);
		free(r);
	}
}

/* ── producer: block on backpressure until space, then copy ── */
static ssize_t ring_write(RingPipe *r, const unsigned char *data, size_t len) {
	size_t written = 0;
	while(written < len) {
		int32_t s = A_LOAD(&r->rseq);
		size_t used = A_LOAD(&r->head) - A_LOAD(&r->tail);
		size_t freeb = r->cap - used;
		if(freeb == 0) {
			if(A_LOAD(&r->reader_closed)) return written ? (ssize_t) written : -1;
			futex_wait(&r->rseq, s);
			continue;
		}
		size_t idx = A_LOAD(&r->head) % r->cap;
		size_t n = len - written;
		if(n > freeb) n = freeb;
		if(n > r->cap - idx) n = r->cap - idx;   /* contiguous chunk */
		memcpy(r->buf + idx, data + written, n);
		A_STORE(&r->head, A_LOAD(&r->head) + n);
		A_BUMP(&r->wseq);
		futex_wake(&r->wseq);
		written += n;
	}
	return (ssize_t) written;
}

/* ── consumer: block until data or EOF, then copy one contiguous chunk ── */
static ssize_t ring_read(RingPipe *r, unsigned char *dst, size_t len) {
	if(len == 0) return 0;
	for(;;) {
		int32_t s = A_LOAD(&r->wseq);
		size_t used = A_LOAD(&r->head) - A_LOAD(&r->tail);
		if(used == 0) {
			if(A_LOAD(&r->closed)) return 0;     /* EOF */
			if(emscripten_is_main_runtime_thread()) {
				/* The producer stage reads its input through FS calls that are
				 * PROXIED to this (main runtime) thread. Hard-blocking on a futex
				 * here would stall those calls and deadlock the producer, so we
				 * service the proxy queue instead and only nap briefly. */
				emscripten_main_thread_process_queued_calls();
				__builtin_wasm_memory_atomic_wait32(&r->wseq, s, 500000LL); /* <=0.5 ms */
			} else {
				futex_wait(&r->wseq, s);
			}
			continue;
		}
		size_t idx = A_LOAD(&r->tail) % r->cap;
		size_t n = len;
		if(n > used) n = used;
		if(n > r->cap - idx) n = r->cap - idx;
		memcpy(dst, r->buf + idx, n);
		A_STORE(&r->tail, A_LOAD(&r->tail) + n);
		A_BUMP(&r->rseq);
		futex_wake(&r->rseq);
		return (ssize_t) n;
	}
}

static void ring_close_write(RingPipe *r) {
	A_STORE(&r->closed, 1);
	A_BUMP(&r->wseq);
	futex_wake(&r->wseq);              /* wake a blocked reader -> sees EOF */
}
static void ring_close_read(RingPipe *r) {
	A_STORE(&r->reader_closed, 1);
	A_BUMP(&r->rseq);
	futex_wake(&r->rseq);             /* wake a blocked writer -> aborts */
}

/* ── fopencookie glue ── */
static ssize_t rc_read(void *c, char *buf, size_t n) {
	return ring_read(((RingEnd*) c)->ring, (unsigned char*) buf, n);
}
static ssize_t rc_write(void *c, const char *buf, size_t n) {
	return ring_write(((RingEnd*) c)->ring, (const unsigned char*) buf, n);
}
static int rc_close(void *c) {
	RingEnd *e = c;
	if(e->is_writer) ring_close_write(e->ring); else ring_close_read(e->ring);
	ring_end_release(e->ring);
	free(e);
	return 0;
}

static int ring_open_pair(RingPipe *r, FILE **rd, FILE **wr) {
	RingEnd *re = malloc(sizeof(RingEnd));
	RingEnd *we = malloc(sizeof(RingEnd));
	if(!re || !we) { free(re); free(we); return -1; }
	re->ring = r; re->is_writer = 0;
	we->ring = r; we->is_writer = 1;
	cookie_io_functions_t rf = { rc_read, NULL, NULL, rc_close };
	cookie_io_functions_t wf = { NULL, rc_write, NULL, rc_close };
	*rd = fopencookie(re, "rb", rf);
	*wr = fopencookie(we, "wb", wf);
	if(!*rd || !*wr) return -1;
	return 0;
}

/* ── kmaPipe implementation (mirrors kmaPipeThread, ring transport) ── */
typedef struct RingPid {
	pthread_t id;
	FILE *fp;          /* end returned to the caller */
	FILE *ioStream;    /* end handed to the child stage */
	char *cmd;
	struct RingPid *next;
} RingPid;

static void * ringThreader(void *arg) {
	RingPid *p = arg;
	char *cmd[2];
	cmd[0] = p->cmd;
	cmd[1] = (char *) p->ioStream;     /* KMA's internal FILE*-in-argv protocol */
	kma_main(0, cmd);
	fclose(p->ioStream);               /* flush + close -> ring EOF for reader */
	return NULL;
}

FILE * kmaPipeRing(const char *cmd, const char *type, FILE *ioStream, int *status) {

	static volatile int Lock = 0;
	static RingPid *pidlist = 0;
	volatile int *lock = &Lock;
	RingPid *src, *last, *dest;
	pthread_t id;
	FILE *rd, *wr;

	if(cmd && type) {
		if((*type != 'r' && *type != 'w') || (type[1] != 0 && type[2] != 0)) {
			errno = EINVAL;
			ERROR();
		}

		dest = malloc(sizeof(RingPid));
		if(!dest) {
			ERROR();
		}
		RingPipe *ring = ring_create(RING_CAP);
		if(!ring || ring_open_pair(ring, &rd, &wr) != 0) {
			ERROR();
		}
		dest->cmd = (char *) cmd;
		if(*type == 'r') {
			dest->fp = rd;         /* caller reads the child's output */
			dest->ioStream = wr;   /* child writes */
		} else {
			dest->fp = wr;         /* caller writes to the child */
			dest->ioStream = rd;   /* child reads */
		}

		if((errno = pthread_create(&dest->id, NULL, &ringThreader, dest))) {
			ERROR();
		}

		lock(lock);
		dest->next = pidlist;
		pidlist = dest;
		unlock(lock);

		return dest->fp;
	} else {
		*status = 0;
		lock(lock);
		for(last = 0, src = pidlist; src && src->fp != ioStream; last = src, src = src->next);
		if(!src) {
			*status = 1;
			unlock(lock);
			return 0;
		}
		unlock(lock);

		id = src->id;
		if((errno = pthread_join(id, NULL))) {
			ERROR();
		}
		fclose(ioStream);          /* close caller's end -> ring freed via refcount */

		lock(lock);
		for(last = 0, src = pidlist; src && src->id != id; last = src, src = src->next);
		if(src) {
			if(!last) pidlist = src->next; else last->next = src->next;
			free(src);
		}
		unlock(lock);

		return 0;
	}
}

#endif /* __EMSCRIPTEN__ */
