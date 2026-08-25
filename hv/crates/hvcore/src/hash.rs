//! Strand-symmetric canonical k-mer hashing and the random-projection derive.
//!
//! Design of record (PLAN.md §4.1):
//!   forward rolling hash  h_f = Σ c_i · B^(k-1-i)          (mod 2^64)
//!   reverse rolling hash  h_r = Σ comp(c_i) · B^i           (mod 2^64)
//! Under strand swap h_f ↔ h_r exactly, so u = splitmix64(h_f ^ h_r) is a
//! strand-symmetric, avalanched 64-bit value. The splitmix finalizer is not
//! optional: low bits of raw polynomial hashes have short periods (LCG
//! structure) and would bias bucket selection.
//!
//! From one finalized u we draw selection / bucket / sign from *independent*
//! randomness (a second mix on disjoint bits), so that conditioning on
//! selection does not bias the bucket distribution.

pub const SPLITMIX_GOLDEN: u64 = 0x9E37_79B9_7F4A_7C15;

#[inline]
pub fn splitmix64(mut z: u64) -> u64 {
    z = z.wrapping_add(SPLITMIX_GOLDEN);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Modular inverse of an odd u64 mod 2^64 (Newton iteration, 6 rounds).
fn odd_inv(x: u64) -> u64 {
    let mut y = x; // correct to 3 bits since x*x ≡ 1 mod 8 for odd x
    for _ in 0..6 {
        y = y.wrapping_mul(2u64.wrapping_sub(x.wrapping_mul(y)));
    }
    y
}

pub struct RollingHash {
    k: usize,
    b: u64,
    bpow: u64,    // B^k   mod 2^64
    bpowm1: u64,  // B^(k-1) mod 2^64
    binv: u64,    // B^-1  mod 2^64
    bpow_tab: Vec<u64>, // B^0 .. B^(k-1) for the fill phase
    hf: u64,
    hr: u64,
    ring: Vec<u8>,
    ring_pos: usize,
    filled: usize, // consecutive valid codes since last reset/invalid
}

impl RollingHash {
    pub fn new(k: usize, seed: u64) -> Self {
        assert!(k >= 2 && k < 4096);
        let b = splitmix64(seed ^ 0x243F_6A88_85A3_08D3) | 1;
        let binv = odd_inv(b);
        let mut bpow = 1u64;
        for _ in 0..k {
            bpow = bpow.wrapping_mul(b);
        }
        let bpowm1 = bpow.wrapping_mul(binv);
        let mut bpow_tab = Vec::with_capacity(k);
        let mut p = 1u64;
        for _ in 0..k {
            bpow_tab.push(p);
            p = p.wrapping_mul(b);
        }
        RollingHash {
            k,
            b,
            bpow,
            bpowm1,
            binv,
            bpow_tab,
            hf: 0,
            hr: 0,
            ring: vec![0; k],
            ring_pos: 0,
            filled: 0,
        }
    }

    pub fn reset(&mut self) {
        self.hf = 0;
        self.hr = 0;
        self.ring_pos = 0;
        self.filled = 0;
    }

    /// Push the next base. `code` is 0..=3 for ACGT; anything else (4) marks
    /// the base invalid (N / low quality): `filled` resets, so the next k
    /// valid bases refill the window from scratch — exactly the k−1 spoiled
    /// windows plus recovery, with no stale ring content involved.
    ///
    /// Fill phase (no eviction, exponent = position in window):
    ///   h_f ← h_f·B + c ;  h_r ← h_r + comp(c)·B^filled
    /// Slide phase (drop outgoing o, add c):
    ///   h_f ← h_f·B − o·B^k + c ;  h_r ← (h_r − comp(o))·B⁻¹ + comp(c)·B^(k-1)
    #[inline]
    pub fn push(&mut self, code: u8) -> Option<u64> {
        if code > 3 {
            self.filled = 0;
            self.ring_pos = 0;
            return None;
        }
        if self.filled < self.k {
            self.ring[self.ring_pos] = code;
            self.ring_pos = (self.ring_pos + 1) % self.k;
            self.hf = self.hf.wrapping_mul(self.b).wrapping_add(code as u64);
            self.hr = self
                .hr
                .wrapping_add((3 - code as u64).wrapping_mul(self.bpow_tab[self.filled]));
            self.filled += 1;
            if self.filled == self.k {
                return Some(splitmix64(self.hf ^ self.hr));
            }
            None
        } else {
            let outgoing = self.ring[self.ring_pos];
            self.ring[self.ring_pos] = code;
            self.ring_pos = (self.ring_pos + 1) % self.k;
            self.hf = self
                .hf
                .wrapping_mul(self.b)
                .wrapping_sub((outgoing as u64).wrapping_mul(self.bpow))
                .wrapping_add(code as u64);
            self.hr = self
                .hr
                .wrapping_sub(3 - outgoing as u64)
                .wrapping_mul(self.binv)
                .wrapping_add((3 - code as u64).wrapping_mul(self.bpowm1));
            Some(splitmix64(self.hf ^ self.hr))
        }
    }
}

/// Selection + projection derive.
pub struct Projector {
    pub dim_bits: u32,
    threshold: u64, // selected iff u < threshold
    kappa: u64,
    pub dim: usize,
}

impl Projector {
    pub fn new(dim_bits: u32, dens_num: u32, dens_den: u32, seed: u64) -> Self {
        assert!(dens_den > 0 && dens_num <= dens_den);
        let threshold = (((dens_num as u128) << 64) / dens_den as u128) as u64;
        Projector {
            dim_bits,
            threshold,
            kappa: splitmix64(seed ^ 0x1319_8A2E_0370_7344),
            dim: 1usize << dim_bits,
        }
    }

    /// (selected, bucket, sign) with bucket/sign drawn from bits independent
    /// of the selection comparison.
    #[inline]
    pub fn derive(&self, u: u64) -> (bool, usize, i32) {
        if u >= self.threshold {
            return (false, 0, 0);
        }
        let w = splitmix64(u ^ self.kappa);
        let bucket = (w >> (64 - self.dim_bits)) as usize;
        let sign = if (w >> (64 - self.dim_bits - 1)) & 1 == 1 {
            1
        } else {
            -1
        };
        (true, bucket, sign)
    }
}

pub fn code_of(base: u8) -> u8 {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rc(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect()
    }

    fn last_window_hash(seq: &[u8], seed: u64) -> Option<u64> {
        let k = 31;
        let mut rh = RollingHash::new(k, seed);
        let mut out = None;
        for &b in seq {
            if let Some(u) = rh.push(code_of(b)) {
                out = Some(u);
            }
        }
        out
    }

    fn first_window_hash(seq: &[u8], seed: u64) -> Option<u64> {
        let k = 31;
        let mut rh = RollingHash::new(k, seed);
        let mut out = None;
        for &b in seq.iter().take(k) {
            out = rh.push(code_of(b));
        }
        out
    }

    #[test]
    fn strand_symmetry() {
        // hash(W) == hash(rc(W)): the FIRST window of rc(seq) is exactly the
        // rc of the LAST window of seq.
        let mut rng: u64 = 88172645463325252;
        for trial in 0..500 {
            let len = 31 + (trial % 60);
            let seq: Vec<u8> = (0..len)
                .map(|_| {
                    rng ^= rng << 13;
                    rng ^= rng >> 7;
                    rng ^= rng << 17;
                    b"ACGT"[(rng % 4) as usize]
                })
                .collect();
            let seed = 12345;
            assert_eq!(
                last_window_hash(&seq, seed),
                first_window_hash(&rc(&seq), seed),
                "seq={}",
                String::from_utf8_lossy(&seq)
            );
        }
    }

    #[test]
    fn rolling_matches_direct() {
        // Rolling evaluation must equal from-scratch polynomial evaluation.
        // Push i (0-indexed) completes the window starting at i−k+1.
        let k = 31usize;
        let seed = 99u64;
        let b = splitmix64(seed ^ 0x243F_6A88_85A3_08D3) | 1;
        let seq: Vec<u8> = (0..200).map(|i| b"ACGT"[(i * 7 + 3) % 4]).collect();
        let codes: Vec<u8> = seq.iter().map(|&c| code_of(c)).collect();
        for wstart in (k - 1)..seq.len() - 1 {
            let mut hf: u64 = 0;
            let mut hr: u64 = 0;
            for j in 0..k {
                let c = codes[wstart + 1 - k + j] as u64;
                hf = hf.wrapping_mul(b).wrapping_add(c);
                hr = hr.wrapping_add((3 - c).wrapping_mul(powmod(b, j)));
            }
            let mut rh2 = RollingHash::new(k, seed);
            rh2.reset();
            let mut got = None;
            for i in 0..=wstart {
                got = rh2.push(codes[i]);
            }
            let expect = splitmix64(hf ^ hr);
            assert_eq!(got, Some(expect), "window at {}", wstart + 1 - k);
        }
    }

    fn powmod(mut b: u64, mut e: usize) -> u64 {
        let mut r: u64 = 1;
        while e > 0 {
            if e & 1 == 1 {
                r = r.wrapping_mul(b);
            }
            b = b.wrapping_mul(b);
            e >>= 1;
        }
        r
    }

    #[test]
    fn bucket_sign_selection_uniformity() {
        let proj = Projector::new(13, 5, 100, 7);
        let mut counts = vec![0u64; proj.dim];
        let mut signs = [0u64; 2];
        let mut selected = 0u64;
        let n: u64 = 2_000_000;
        let mut rng: u64 = 0x9E3779B97F4A7C15;
        for i in 0..n {
            rng = splitmix64(rng.wrapping_add(i));
            let (sel, b, s) = proj.derive(rng);
            if sel {
                selected += 1;
                counts[b] += 1;
                signs[(s > 0) as usize] += 1;
            }
        }
        let dens = selected as f64 / n as f64;
        assert!((dens - 0.05).abs() < 0.002, "density {dens}");
        let expect = selected as f64 / proj.dim as f64;
        let mx = counts.iter().copied().max().unwrap() as f64;
        let mn = counts.iter().copied().min().unwrap() as f64;
        assert!(
            mx - expect < 25.0 && expect - mn < 25.0,
            "bucket spread {mn}..{mx} expect {expect}"
        );
        let s0 = signs[0] as f64;
        assert!(
            (s0 / selected as f64 - 0.5).abs() < 0.01,
            "sign balance {s0}/{}",
            selected
        );
    }

    #[test]
    fn invalid_base_spoils_window() {
        let mut rh = RollingHash::new(31, 1);
        for _ in 0..30 {
            assert!(rh.push(code_of(b'A')).is_none());
        }
        assert!(rh.push(code_of(b'A')).is_some(), "push 31 completes window 1");
        rh.push(code_of(b'N'));
        // the next 30 pushes rebuild the window; the 31st emits again
        let mut emissions = 0;
        for _ in 0..30 {
            if rh.push(code_of(b'A')).is_some() {
                emissions += 1;
            }
        }
        assert_eq!(emissions, 0, "spoiled windows must not emit");
        assert!(rh.push(code_of(b'A')).is_some(), "window must recover");
    }
}
