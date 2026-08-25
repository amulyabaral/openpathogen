//! The sample fingerprint (.hvf): the whole sample as one bounded artifact.
//!
//! Header carries all encoder parameters so a fingerprint can only be
//! profiled against a parameter-matching DB — mismatches fail loudly.
//! The vector is i16-quantized (RMSE ≤ 0.4% of σ̂ — guarded by test T3).
//! Optional Gaussian noise implements the DP release mode (PLAN §4.7).

use crate::encode::{Accumulator, EncodeParams};

pub const FP_MAGIC: &[u8; 4] = b"HVF1";
pub const FP_FLAG_DP: u16 = 1;

#[derive(Clone)]
pub struct Fingerprint {
    pub params: EncodeParams,
    pub n_bases: u64,
    pub n_events: u64,
    pub scale: f32, // multiplier: v = q · scale
    pub dp_flags: u16,
    pub dp_eps: f32,
    pub quant: Vec<i16>,
}

/// Deterministic PRNG for the DP noise (seeded, so native/wasm agree).
struct NoiseRng(u64);
impl NoiseRng {
    fn next_u64(&mut self) -> u64 {
        // xorshift64*
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn normal(&mut self) -> f64 {
        // Box–Muller (adequate for DP noise)
        loop {
            let u1 = self.next_u64() as f64 / u64::MAX as f64;
            let u2 = self.next_u64() as f64 / u64::MAX as f64;
            if u1 > 1e-300 {
                return (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            }
        }
    }
}

impl Fingerprint {
    pub fn from_accumulator(
        acc: &Accumulator,
        dp: Option<(f64, f64, u64)>, // (eps, delta, noise_seed)
    ) -> Fingerprint {
        let maxabs = acc
            .acc
            .iter()
            .map(|&v| v.unsigned_abs() as f64)
            .fold(0.0f64, f64::max)
            .max(1.0);
        let scale = (maxabs / 32767.0) as f32;
        let mut dp_flags = 0u16;
        let mut dp_eps = 0f32;
        let mut noisy: Vec<f64> = acc.acc.iter().map(|&v| v as f64).collect();
        if let Some((eps, delta, seed)) = dp {
            // sensitivity: one read contributes ≤ ⌈ρ·readlen⌉ signed events
            let rho = acc.params().dens_num as f64 / acc.params().dens_den as f64;
            let sens = (rho * 150.0f64).sqrt();
            let sd = (2.0 * (1.25 / delta).ln()).sqrt() * sens / eps;
            let mut rng = NoiseRng(seed);
            for v in noisy.iter_mut() {
                *v += rng.normal() * sd;
            }
            dp_flags = FP_FLAG_DP;
            dp_eps = eps as f32;
        }
        let quant: Vec<i16> = noisy
            .iter()
            .map(|&v| ((v / scale as f64).round().clamp(-32767.0, 32767.0)) as i16)
            .collect();
        Fingerprint {
            params: *acc.params(),
            n_bases: acc.stats.n_bases,
            n_events: acc.stats.n_events,
            scale,
            dp_flags,
            dp_eps,
            quant,
        }
    }

    /// Dequantized accumulator-compatible vector.
    pub fn vector(&self) -> Vec<i32> {
        self.quant
            .iter()
            .map(|&q| (q as f64 * self.scale as f64).round() as i32)
            .collect()
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(80 + self.quant.len() * 2);
        out.extend_from_slice(FP_MAGIC);
        out.extend_from_slice(&1u32.to_le_bytes());
        out.extend_from_slice(&self.params.k.to_le_bytes());
        out.extend_from_slice(&self.params.dim_bits.to_le_bytes());
        out.extend_from_slice(&self.params.dens_num.to_le_bytes());
        out.extend_from_slice(&self.params.dens_den.to_le_bytes());
        out.extend_from_slice(&self.params.qmin.to_le_bytes());
        out.extend_from_slice(&self.params.seed.to_le_bytes());
        out.extend_from_slice(&self.n_bases.to_le_bytes());
        out.extend_from_slice(&self.n_events.to_le_bytes());
        out.extend_from_slice(&self.scale.to_le_bytes());
        out.extend_from_slice(&self.dp_flags.to_le_bytes());
        out.extend_from_slice(&self.dp_eps.to_le_bytes());
        for &q in &self.quant {
            out.extend_from_slice(&q.to_le_bytes());
        }
        out
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<Fingerprint, String> {
        if bytes.len() < 60 || &bytes[0..4] != FP_MAGIC {
            return Err("hvf: bad magic".into());
        }
        let rd_u32 =
            |o: usize| u32::from_le_bytes(bytes[o..o + 4].try_into().unwrap());
        let params = EncodeParams {
            k: rd_u32(8),
            dim_bits: rd_u32(12),
            dens_num: rd_u32(16),
            dens_den: rd_u32(20),
            qmin: bytes[24],
            seed: u64::from_le_bytes(bytes[25..33].try_into().unwrap()),
        };
        let d = params.dim();
        let n = bytes.len();
        if n < 59 + d * 2 {
            return Err("hvf: truncated".into());
        }
        let mut off = 33;
        let n_bases = u64::from_le_bytes(bytes[off..off + 8].try_into().unwrap());
        off += 8;
        let n_events = u64::from_le_bytes(bytes[off..off + 8].try_into().unwrap());
        off += 8;
        let scale = f32::from_le_bytes(bytes[off..off + 4].try_into().unwrap());
        off += 4;
        let dp_flags = u16::from_le_bytes(bytes[off..off + 2].try_into().unwrap());
        off += 2;
        let dp_eps = f32::from_le_bytes(bytes[off..off + 4].try_into().unwrap());
        off += 4;
        let mut quant = Vec::with_capacity(d);
        for i in 0..d {
            quant.push(i16::from_le_bytes(
                bytes[off + i * 2..off + i * 2 + 2].try_into().unwrap(),
            ));
        }
        Ok(Fingerprint {
            params,
            n_bases,
            n_events,
            scale,
            dp_flags,
            dp_eps,
            quant,
        })
    }

    /// Aggregate (federate): elementwise sum of compatible fingerprints.
    /// Population prevalence = deconvolution of the sum — free by linearity.
    pub fn aggregate(fps: &[&Fingerprint]) -> Result<Fingerprint, String> {
        if fps.is_empty() {
            return Err("aggregate: empty".into());
        }
        let first = &fps[0];
        for f in fps {
            if !first.params.are_compatible(&f.params) {
                return Err("aggregate: parameter mismatch".into());
            }
        }
        let maxabs = fps
            .iter()
            .flat_map(|f| f.quant.iter())
            .map(|&q| q.unsigned_abs() as f64)
            .fold(0.0f64, f64::max)
            .max(1.0);
        let scale = (maxabs / 32767.0) as f32;
        let quant: Vec<i16> = (0..first.quant.len())
            .map(|i| {
                let s: f64 = fps.iter().map(|f| f.quant[i] as f64 * f.scale as f64).sum();
                ((s / scale as f64).round().clamp(-32767.0, 32767.0)) as i16
            })
            .collect();
        Ok(Fingerprint {
            params: first.params,
            n_bases: fps.iter().map(|f| f.n_bases).sum(),
            n_events: fps.iter().map(|f| f.n_events).sum(),
            scale,
            dp_flags: if fps.iter().all(|f| f.dp_flags & FP_FLAG_DP != 0) {
                FP_FLAG_DP
            } else {
                0
            },
            dp_eps: fps
                .iter()
                .map(|f| f.dp_eps as f64)
                .fold(f64::INFINITY, f64::min) as f32,
            quant,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_and_quantization_error() {
        let params = EncodeParams::default_params();
        let mut acc = Accumulator::new(params);
        // many pseudo-random reads → Gaussian-ish bucket distribution, MAD > 0
        let mut rng: u64 = 0xDEADBEEF;
        let mut next = || {
            rng ^= rng << 13;
            rng ^= rng >> 7;
            rng ^= rng << 17;
            rng
        };
        for _ in 0..4000 {
            let seq: Vec<u8> = (0..250).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
            let q = vec![b'I'; seq.len()];
            acc.encode_read(&seq, Some(&q));
        }
        let fp = Fingerprint::from_accumulator(&acc, None);
        let bytes = fp.to_bytes();
        let fp2 = Fingerprint::from_bytes(&bytes).unwrap();
        assert_eq!(fp2.quant, fp.quant);
        assert_eq!(fp2.n_events, fp.n_events);

        let v = fp.vector();
        let err_rms = (0..v.len())
            .map(|i| {
                let e = (v[i] - acc.acc[i]) as f64;
                e * e
            })
            .sum::<f64>()
            .sqrt()
            / v.len() as f64;
        let sigma = crate::stats::mad_sigma(&acc.acc);
        assert!(err_rms < 0.004 * sigma, "quant err {err_rms} vs sigma {sigma}");
    }

    #[test]
    fn aggregation_is_linear() {
        let params = EncodeParams::default_params();
        let mk = |shift: u64| {
            let mut acc = Accumulator::new(params);
            let seq: Vec<u8> = (0..3000)
                .map(|i| b"ACGT"[((i + shift as usize) * 7 + 1) % 4])
                .collect();
            let q = vec![b'I'; seq.len()];
            acc.encode_read(&seq, Some(&q));
            Fingerprint::from_accumulator(&acc, None)
        };
        let a = mk(0);
        let b = mk(5);
        let agg = Fingerprint::aggregate(&[&a, &b]).unwrap();
        let va = a.vector();
        let vb = b.vector();
        let vagg = agg.vector();
        for i in 0..va.len() {
            assert!((va[i] + vb[i] - vagg[i]).abs() <= 2, "bucket {i}");
        }
    }
}
