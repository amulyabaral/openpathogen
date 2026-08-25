//! Simulation: synthetic genomes, ANI-controlled mutants, Illumina-style
//! reads with per-base errors and quality strings. Ground truth by
//! construction — used by the calibration tests (T4–T8) and LOD sweeps.

use crate::hash::splitmix64;

pub struct Rng(pub u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(splitmix64(seed) | 1)
    }
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }
    #[inline]
    pub fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }
    pub fn unit(&mut self) -> f64 {
        self.next_u64() as f64 / u64::MAX as f64
    }
}

/// Uniform random genome (IUPAC ACGT). `gc` in [0,1] controls composition.
pub fn synth_genome(rng: &mut Rng, len: usize, gc: f64) -> Vec<u8> {
    let mut out = Vec::with_capacity(len);
    for _ in 0..len {
        let c = rng.unit();
        let b = if c < gc / 2.0 {
            b'G'
        } else if c < gc {
            b'C'
        } else if c < gc + (1.0 - gc) / 2.0 {
            b'A'
        } else {
            b'T'
        };
        out.push(b);
    }
    out
}

/// Substitution divergence at rate `rate` (per base). ANI ≈ 1 − rate.
pub fn mutate(rng: &mut Rng, seq: &[u8], rate: f64) -> Vec<u8> {
    let mut out = seq.to_vec();
    for b in out.iter_mut() {
        if rng.unit() < rate {
            *b = match *b {
                b'A' => b"CGT"[rng.below(3) as usize],
                b'C' => b"AGT"[rng.below(3) as usize],
                b'G' => b"ACT"[rng.below(3) as usize],
                _ => b"ACG"[rng.below(3) as usize],
            };
        }
    }
    out
}

#[derive(Clone)]
pub struct Read {
    pub name: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

pub struct ReadSimParams {
    pub read_len: usize,
    pub err_rate: f64,   // per-base substitution error rate
    pub fragment: usize, // mean mate gap for paired reads
    pub paired: bool,
}

impl Default for ReadSimParams {
    fn default() -> Self {
        ReadSimParams {
            read_len: 150,
            err_rate: 0.001,
            fragment: 200,
            paired: false,
        }
    }
}

/// Simulate sequencing of `genome` at k-mer coverage `coverage`.
///
/// Semantics: `coverage` is the expected number of occurrences of each
/// genomic k-mer — exactly the quantity the deconvolution fit estimates.
/// Reads start at uniform positions; n_reads = cov·glen/(read_len−k+1).
/// Error bases carry low quality (Q10–18) so Q20 masking removes most.
pub fn sim_reads(
    rng: &mut Rng,
    genome: &[u8],
    coverage: f64,
    p: &ReadSimParams,
    read_prefix: usize,
) -> Vec<Read> {
    const K: usize = 31;
    let glen = genome.len();
    let win = p.read_len.saturating_sub(K - 1).max(1);
    let n_reads = ((coverage * glen as f64) / win as f64).round() as usize;
    let mut out = Vec::with_capacity(n_reads * if p.paired { 2 } else { 1 });
    for i in 0..n_reads {
        let start = rng.below((glen - p.read_len + 1) as u64) as usize;
        out.push(mk_read(rng, genome, start, p.read_len, p.err_rate, format!("sim{read_prefix}_{i}_1")));
        if p.paired {
            let fpos = (start + p.read_len + (rng.below(60) as usize)).min(glen - p.read_len);
            out.push(mk_read(rng, genome, fpos, p.read_len, p.err_rate, format!("sim{read_prefix}_{i}_2")));
        }
    }
    out
}

fn mk_read(
    rng: &mut Rng,
    genome: &[u8],
    start: usize,
    len: usize,
    err_rate: f64,
    name: String,
) -> Read {
    let mut seq = Vec::with_capacity(len);
    let mut qual = Vec::with_capacity(len);
    for &b in &genome[start..start + len] {
        if rng.unit() < err_rate {
            seq.push(match b {
                b'A' => b"CGT"[rng.below(3) as usize],
                b'C' => b"AGT"[rng.below(3) as usize],
                b'G' => b"ACT"[rng.below(3) as usize],
                _ => b"ACG"[rng.below(3) as usize],
            });
            qual.push(33 + 10 + rng.below(9) as u8);
        } else {
            seq.push(b);
            qual.push(33 + 37); // Q37
        }
    }
    Read { name, seq, qual }
}

fn gauss(rng: &mut Rng) -> f64 {
    let u1 = rng.unit().max(1e-12);
    let u2 = rng.unit();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

/// FASTQ text for a batch of reads.
pub fn reads_to_fastq(reads: &[Read]) -> Vec<u8> {
    let mut out = Vec::with_capacity(reads.len() * 340);
    for r in reads {
        out.extend_from_slice(b"@");
        out.extend_from_slice(r.name.as_bytes());
        out.push(b'\n');
        out.extend_from_slice(&r.seq);
        out.extend_from_slice(b"\n+\n");
        out.extend_from_slice(&r.qual);
        out.push(b'\n');
    }
    out
}

/// FASTA text for a genome.
pub fn genome_to_fasta(id: &str, seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len() + 128);
    out.extend_from_slice(b">");
    out.extend_from_slice(id.as_bytes());
    out.push(b'\n');
    let mut i = 0;
    while i < seq.len() {
        let end = (i + 70).min(seq.len());
        out.extend_from_slice(&seq[i..end]);
        out.push(b'\n');
        i = end;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn synth_and_mutate() {
        let mut rng = Rng::new(42);
        let g = synth_genome(&mut rng, 100_000, 0.5);
        assert_eq!(g.len(), 100_000);
        assert!(g.iter().all(|&b| b"ACGT".contains(&b)));
        let m = mutate(&mut rng, &g, 0.01);
        let diff = g.iter().zip(m.iter()).filter(|(a, b)| a != b).count();
        assert!((diff as f64 / 100_000.0 - 0.01).abs() < 0.002, "{diff}");
    }

    #[test]
    fn read_sim_coverage_and_quality() {
        let mut rng = Rng::new(7);
        let g = synth_genome(&mut rng, 50_000, 0.5);
        let reads = sim_reads(&mut rng, &g, 1.0, &ReadSimParams::default(), 0);
        let total: usize = reads.iter().map(|r| r.seq.len()).sum();
        // ~2× fragment per frag pair; coverage ≈ 1× within Poisson noise
        let cov = total as f64 / 50_000.0;
        assert!((0.5..1.6).contains(&cov), "coverage {cov}");
        for r in &reads {
            assert_eq!(r.seq.len(), r.qual.len());
            assert!(r.qual.iter().all(|&q| q >= 33 + 9));
        }
    }
}
