//! Streaming FASTQ/FASTA encoding into the fixed-length signed accumulator.
//!
//! The accumulator is the whole state of a sample: Int32Array(D). Reads are
//! independent molecules → the rolling hash resets per read; contigs of a
//! reference genome chain (no reset) but reset between records.

use crate::hash::{code_of, Projector, RollingHash};

#[derive(Clone, Copy, PartialEq, Eq)]
pub struct EncodeParams {
    pub k: u32,
    pub dim_bits: u32,
    pub dens_num: u32,
    pub dens_den: u32,
    pub qmin: u8, // Phred threshold; bases below are masked to N
    pub seed: u64,
}

impl EncodeParams {
    pub fn default_params() -> Self {
        EncodeParams {
            k: 31,
            dim_bits: 13, // D = 8192
            dens_num: 5,
            dens_den: 100,
            qmin: 20,
            seed: 0x5EED_0001,
        }
    }
    pub fn dim(&self) -> usize {
        1usize << self.dim_bits
    }
    pub fn are_compatible(&self, o: &EncodeParams) -> bool {
        self.k == o.k
            && self.dim_bits == o.dim_bits
            && self.dens_num == o.dens_num
            && self.dens_den == o.dens_den
            && self.seed == o.seed
    }
}

#[derive(Default, Clone, Copy)]
pub struct EncodeStats {
    pub n_bases: u64,   // valid ACGT bases seen
    pub n_events: u64,  // selected k-mer occurrences accumulated
    pub n_reads: u64,
    pub n_masked: u64,  // bases masked by quality / N
    pub n_records: u64, // fasta records
}

/// The sample accumulator. One i32 per bucket; saturation is impossible at
/// practical depths (expected |acc| ≲ 10^4) but we saturate defensively.
pub struct Accumulator {
    pub acc: Vec<i32>,
    pub stats: EncodeStats,
    params: EncodeParams,
    proj: Projector,
}

impl Accumulator {
    pub fn new(params: EncodeParams) -> Self {
        let proj = Projector::new(params.dim_bits, params.dens_num, params.dens_den, params.seed);
        Accumulator {
            acc: vec![0i32; proj.dim],
            stats: EncodeStats::default(),
            params,
            proj,
        }
    }

    pub fn params(&self) -> &EncodeParams {
        &self.params
    }

    pub fn reset(&mut self) {
        self.acc.iter_mut().for_each(|v| *v = 0);
        self.stats = EncodeStats::default();
    }

    /// Encode one read (fresh molecule): hash state resets first.
    pub fn encode_read(&mut self, seq: &[u8], qual: Option<&[u8]>) {
        let mut rh = RollingHash::new(self.params.k as usize, self.params.seed);
        self.stats.n_reads += 1;
        for (i, &base) in seq.iter().enumerate() {
            let code = if let Some(q) = qual {
                if i < q.len() && q[i].saturating_sub(33) < self.params.qmin {
                    self.stats.n_masked += 1;
                    4
                } else {
                    code_of(base)
                }
            } else {
                code_of(base)
            };
            if code == 4 && code_of(base) == 4 {
                self.stats.n_masked += 1;
            }
            self.stats.n_bases += (code <= 3) as u64;
            if let Some(u) = rh.push(code) {
                let (sel, bucket, sign) = self.proj.derive(u);
                if sel {
                    let v = self.acc[bucket].saturating_add(sign);
                    self.acc[bucket] = v;
                    self.stats.n_events += 1;
                }
            }
        }
    }

    /// Encode a stretch of one contig (hash state chains across calls; caller
    /// resets between records via `contig_reset`).
    pub fn encode_contig(&mut self, rh: &mut RollingHash, seq: &[u8]) {
        for &base in seq {
            let code = code_of(base);
            if let Some(u) = rh.push(code) {
                let (sel, bucket, sign) = self.proj.derive(u);
                if sel {
                    let v = self.acc[bucket].saturating_add(sign);
                    self.acc[bucket] = v;
                    self.stats.n_events += 1;
                }
            }
        }
    }

    pub fn contig_hasher(&self) -> RollingHash {
        RollingHash::new(self.params.k as usize, self.params.seed)
    }
}

// ─────────────────────────────────────────────────────────────────────────
// FASTQ byte-stream state machine (chunk-boundary safe)
// ─────────────────────────────────────────────────────────────────────────

#[derive(PartialEq, Clone, Copy)]
enum LineKind {
    Header,
    Seq,
    Plus,
    Qual,
}

pub struct FastqScanner {
    kind: LineKind,
    line: Vec<u8>,     // current line being accumulated
    seq_buf: Vec<u8>,  // completed seq line of the current read
    malformed: u64,
}

impl FastqScanner {
    pub fn new() -> Self {
        FastqScanner {
            kind: LineKind::Header,
            line: Vec::with_capacity(4096),
            seq_buf: Vec::with_capacity(4096),
            malformed: 0,
        }
    }

    /// Feed raw bytes; complete reads are dispatched to `sink`.
    pub fn feed(&mut self, bytes: &[u8], acc: &mut Accumulator) {
        for &c in bytes {
            if c == b'\n' {
                self.dispatch(acc);
                self.line.clear();
            } else if c != b'\r' {
                self.line.push(c);
            }
        }
    }

    /// Flush a trailing line at EOF (files normally end with \n; tolerate).
    pub fn finish(&mut self, acc: &mut Accumulator) {
        if !self.line.is_empty() {
            self.dispatch(acc);
            self.line.clear();
        }
    }

    fn dispatch(&mut self, acc: &mut Accumulator) {
        match self.kind {
            LineKind::Header => {
                if self.line.first() == Some(&b'@') {
                    self.kind = LineKind::Seq;
                }
                // Blank/stray lines: stay on header.
            }
            LineKind::Seq => {
                std::mem::swap(&mut self.seq_buf, &mut self.line);
                self.seq_buf.shrink_to_fit();
                self.kind = LineKind::Plus;
            }
            LineKind::Plus => {
                if self.line.first() == Some(&b'+') {
                    self.kind = LineKind::Qual;
                } else if self.line.first() == Some(&b'@') {
                    // seq line missing/empty (unusual); recover by treating
                    // this as the next header.
                    self.kind = LineKind::Seq;
                    self.seq_buf.clear();
                }
            }
            LineKind::Qual => {
                let q = &self.line;
                let s = &self.seq_buf;
                if q.len() != s.len() {
                    self.malformed += 1;
                }
                let n = q.len().min(s.len());
                acc.encode_read(&s[..n], Some(&q[..n]));
                self.seq_buf.clear();
                self.kind = LineKind::Header;
            }
        }
    }
}

impl Default for FastqScanner {
    fn default() -> Self {
        Self::new()
    }
}

// ─────────────────────────────────────────────────────────────────────────
// FASTA byte-stream state machine
// ─────────────────────────────────────────────────────────────────────────

pub struct FastaScanner {
    line: Vec<u8>,
    in_record: bool,
}

impl FastaScanner {
    pub fn new() -> Self {
        FastaScanner {
            line: Vec::with_capacity(4096),
            in_record: false,
        }
    }

    /// Feed raw bytes; sequence lines are encoded contiguously per record.
    pub fn feed(&mut self, bytes: &[u8], acc: &mut Accumulator, rh: &mut RollingHash) {
        for &c in bytes {
            if c == b'\n' {
                self.dispatch(acc, rh);
                self.line.clear();
            } else if c != b'\r' {
                self.line.push(c);
            }
        }
    }

    pub fn finish(&mut self, acc: &mut Accumulator, rh: &mut RollingHash) {
        if !self.line.is_empty() {
            self.dispatch(acc, rh);
            self.line.clear();
        }
    }

    fn dispatch(&mut self, acc: &mut Accumulator, rh: &mut RollingHash) {
        if self.line.first() == Some(&b'>') {
            rh.reset();
            acc.stats.n_records += 1;
            self.in_record = true;
        } else if self.in_record && !self.line.is_empty() {
            acc.encode_contig(rh, &self.line);
        }
    }
}

impl Default for FastaScanner {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_fastq(recs: &[(&str, &str, &str)]) -> Vec<u8> {
        let mut out = Vec::new();
        for (id, seq, q) in recs {
            out.extend_from_slice(format!("@{}\n{}\n+\n{}\n", id, seq, q).as_bytes());
        }
        out
    }

    #[test]
    fn fastq_chunk_boundary_independence() {
        // aperiodic sequence (LCG-generated) long enough for real events
        let mut rng: u64 = 0xCAFEBABE;
        let seq: Vec<u8> = (0..3000)
            .map(|_| {
                rng ^= rng << 13;
                rng ^= rng >> 7;
                rng ^= rng << 17;
                b"ACGT"[(rng % 4) as usize]
            })
            .collect();
        let seq = String::from_utf8(seq).unwrap();
        let q = "I".repeat(seq.len());
        let bytes = write_fastq(&[("r1", &seq, &q), ("r2", &seq, &q)]);

        let params = EncodeParams::default_params();
        let mut a = Accumulator::new(params);
        let mut s = FastqScanner::new();
        s.feed(&bytes, &mut a);
        s.finish(&mut a);

        let mut b = Accumulator::new(params);
        let mut t = FastqScanner::new();
        for chunk in bytes.chunks(7) {
            t.feed(chunk, &mut b);
        }
        t.finish(&mut b);

        assert_eq!(a.acc, b.acc);
        assert_eq!(a.stats.n_events, b.stats.n_events);
        assert_eq!(a.stats.n_reads, 2);
        assert!(a.stats.n_events > 0);
    }

    #[test]
    fn quality_masking() {
        let params = EncodeParams::default_params();
        let mut rng: u64 = 0x5EED;
        let seq: Vec<u8> = (0..2000)
            .map(|_| {
                rng ^= rng << 13;
                rng ^= rng >> 7;
                rng ^= rng << 17;
                b"ACGT"[(rng % 4) as usize]
            })
            .collect();
        let good = "I".repeat(seq.len()); // Q40
        let bad = "#".repeat(seq.len()); // Q2 → all masked

        let mut a = Accumulator::new(params);
        a.encode_read(&seq, Some(good.as_bytes()));
        let events_good = a.stats.n_events;

        let mut b = Accumulator::new(params);
        b.encode_read(&seq, Some(bad.as_bytes()));
        assert_eq!(b.stats.n_events, 0);
        assert!(b.stats.n_masked >= 60);
        assert!(events_good > 0);
    }

    #[test]
    fn accumulation_is_linear_in_coverage() {
        // v(2 copies of a read) ≈ 2 · v(1 copy) — trivially exact, but guards
        // against saturating/int promotion regressions.
        let params = EncodeParams::default_params();
        let seq: Vec<u8> = (0..200).map(|i| b"ACGT"[(i * 5 + 1) % 4]).collect();
        let q = "I".repeat(seq.len());
        let mut one = Accumulator::new(params);
        one.encode_read(&seq, Some(q.as_bytes()));
        let mut two = Accumulator::new(params);
        two.encode_read(&seq, Some(q.as_bytes()));
        two.encode_read(&seq, Some(q.as_bytes()));
        for (x, y) in one.acc.iter().zip(two.acc.iter()) {
            assert_eq!(2 * *x, *y);
        }
    }
}
