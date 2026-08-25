//! hvwasm — browser bindings for hvcore.
//!
//! Hand-rolled `extern "C"` exports over linear memory (no wasm-bindgen):
//! JS stages FASTQ chunks into an input buffer, feeds the encoder, then
//! copies out the accumulator / profile results / fingerprint bytes. All
//! heavy state lives in the wasm heap; the FASTQ never does (streamed).
//!
//! Contract notes:
//!   - pointers are byte offsets into wasm memory; re-read them after any
//!     call that may grow memory (hv_db_stage, feeds that grow staging).
//!   - results are packed f64s: 7 header values then 6 per taxon
//!     [idx, coverage, z, q, dna_frac, flags].

use hvcore::db::HvDb;
use hvcore::encode::{Accumulator, EncodeParams, FastqScanner};
use hvcore::fp::Fingerprint;
use hvcore::solve::profile;
use std::cell::RefCell;

struct Engine {
    params: EncodeParams,
    acc: Accumulator,
    scanner: FastqScanner,
    staging: Vec<u8>,
    db_bytes: Vec<u8>,
    db: Option<HvDb>,
    db_meta: Vec<u8>,
    results: Vec<f64>,
    fp_bytes: Vec<u8>,
}

thread_local! {
    static ENGINE: RefCell<Option<Engine>> = const { RefCell::new(None) };
}

const STAGE_CAP: usize = 4 << 20; // 4 MiB FASTQ chunk staging

fn with_engine<T>(f: impl FnOnce(&mut Engine) -> T) -> Option<T> {
    ENGINE.with(|e| {
        let mut e = e.borrow_mut();
        match e.as_mut() {
            Some(eng) => Some(f(eng)),
            None => None,
        }
    })
}

/// (k, dim_bits, dens_num, dens_den, seed, qmin) → 0 on success.
#[no_mangle]
pub extern "C" fn hv_encoder_reset(
    k: u32,
    dim_bits: u32,
    dens_num: u32,
    dens_den: u32,
    seed: u64,
    qmin: u8,
) -> i32 {
    let params = EncodeParams {
        k,
        dim_bits,
        dens_num,
        dens_den,
        qmin,
        seed,
    };
    let acc = Accumulator::new(params);
    ENGINE.with(|e| {
        *e.borrow_mut() = Some(Engine {
            params,
            acc,
            scanner: FastqScanner::new(),
            staging: vec![0u8; STAGE_CAP],
            db_bytes: Vec::new(),
            db: None,
            db_meta: Vec::new(),
            results: Vec::new(),
            fp_bytes: Vec::new(),
        })
    });
    0
}

#[no_mangle]
pub extern "C" fn hv_stage_ptr() -> *mut u8 {
    with_engine(|e| e.staging.as_mut_ptr()).unwrap_or(std::ptr::null_mut())
}

#[no_mangle]
pub extern "C" fn hv_stage_cap() -> usize {
    with_engine(|e| e.staging.len()).unwrap_or(0)
}

/// Feed `len` staged bytes; `eof` finishes the stream. Returns cumulative
/// selected-event count (as f64 — u64-safe for JS).
#[no_mangle]
pub extern "C" fn hv_encoder_feed(len: usize, eof: i32) -> f64 {
    with_engine(|e| {
        let (staging_ptr, cap) = (e.staging.as_ptr(), e.staging.len());
        debug_assert!(len <= cap);
        let bytes = unsafe { std::slice::from_raw_parts(staging_ptr, len.min(cap)) };
        e.scanner.feed(bytes, &mut e.acc);
        if eof != 0 {
            e.scanner.finish(&mut e.acc);
        }
        e.acc.stats.n_events as f64
    })
    .unwrap_or(-1.0)
}

#[no_mangle]
pub extern "C" fn hv_stats_bases() -> f64 {
    with_engine(|e| e.acc.stats.n_bases as f64).unwrap_or(-1.0)
}

#[no_mangle]
pub extern "C" fn hv_stats_reads() -> f64 {
    with_engine(|e| e.acc.stats.n_reads as f64).unwrap_or(-1.0)
}

#[no_mangle]
pub extern "C" fn hv_acc_ptr() -> *const i32 {
    with_engine(|e| e.acc.acc.as_ptr()).unwrap_or(std::ptr::null())
}

/// Grow the DB staging buffer to `len` and return its pointer.
#[no_mangle]
pub extern "C" fn hv_db_stage(len: usize) -> *mut u8 {
    with_engine(|e| {
        if e.db_bytes.len() < len {
            e.db_bytes.resize(len, 0);
        }
        e.db_bytes.as_mut_ptr()
    })
    .unwrap_or(std::ptr::null_mut())
}

/// Load DB + sidecar metadata from staged bytes. Returns row count ≥ 0.
#[no_mangle]
pub extern "C" fn hv_db_load(hvd_len: usize, tsv_len: usize) -> i32 {
    with_engine(|e| {
        let hvd = e.db_bytes[..hvd_len.min(e.db_bytes.len())].to_vec();
        let mut db = match HvDb::from_bytes(&hvd) {
            Ok(db) => db,
            Err(_) => return -1,
        };
        if tsv_len > 0 {
            let tsv = e.db_bytes[hvd_len..(hvd_len + tsv_len).min(e.db_bytes.len())].to_vec();
            let _ = db.attach_meta_tsv(&tsv);
        }
        let n = db.n() as i32;
        e.db = Some(db);
        n
    })
    .unwrap_or(-2)
}

/// Run the profile. Returns taxon count ≥ 0 (results via hv_results_ptr).
#[no_mangle]
pub extern "C" fn hv_profile(zmin: f64, kmax: usize) -> i32 {
    with_engine(|e| {
        let db = match &e.db {
            Some(db) => db,
            None => return -1,
        };
        let params = e.params;
        let pp = hvcore::solve::ProfileParams {
            zmin,
            kmax,
            ..Default::default()
        };
        let v = e.acc.acc.clone();
        let n_events = e.acc.stats.n_events;
        let rep = profile(&v, n_events, db, pp);
        let mut out = vec![
            rep.sigma,
            rep.unexplained,
            rep.spike_z,
            rep.screen_max,
            rep.ev_bound,
            rep.n_events as f64,
            rep.rounds as f64,
        ];
        let n = rep.taxa.len();
        for t in &rep.taxa {
            out.extend_from_slice(&[t.idx as f64, t.x, t.z, t.q, t.dna_frac, t.flags as f64]);
        }
        e.results = out;
        let _ = params;
        n as i32
    })
    .unwrap_or(-2)
}

#[no_mangle]
pub extern "C" fn hv_results_ptr() -> *const f64 {
    with_engine(|e| e.results.as_ptr()).unwrap_or(std::ptr::null())
}

#[no_mangle]
pub extern "C" fn hv_results_len() -> usize {
    with_engine(|e| e.results.len()).unwrap_or(0)
}

/// Serialize the current accumulator as fingerprint bytes.
#[no_mangle]
pub extern "C" fn hv_fp_serialize(dp_eps: f64) -> usize {
    with_engine(|e| {
        let dp = if dp_eps > 0.0 {
            Some((dp_eps, 1e-6, e.params.seed ^ 0xD00D))
        } else {
            None
        };
        let fp = Fingerprint::from_accumulator(&e.acc, dp);
        e.fp_bytes = fp.to_bytes();
        e.fp_bytes.len()
    })
    .unwrap_or(0)
}

#[no_mangle]
pub extern "C" fn hv_fp_ptr() -> *const u8 {
    with_engine(|e| e.fp_bytes.as_ptr()).unwrap_or(std::ptr::null())
}

#[cfg(test)]
mod tests {
    #[test]
    fn native_builds_too() {
        // the same crate compiles natively for CLI-side parity checks
        assert!(true);
    }
}
