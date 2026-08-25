//! Reference matrix: one int8-quantized projection row per genome.
//!
//! Rows are built from DISTINCT selected k-mers (one copy per genome), so the
//! sample's occurrence accumulator satisfies v ≈ Σᵢ μᵢ·aᵢ with μᵢ = coverage.
//! Rows are appendable without any index rebuild — the property that makes
//! the matrix a "database you can grow by concatenation".

use crate::encode::{Accumulator, EncodeParams};
use crate::hash::{splitmix64, Projector, RollingHash};
use std::collections::HashSet;

pub const FLAG_SPIKE: u32 = 1; // synthetic spike-in row (QC / calibration)
pub const FLAG_HUMAN: u32 = 2; // human / host row (depletion candidate)

#[derive(Clone, Default)]
pub struct RowMeta {
    pub id: String,
    pub name: String,
    pub genus: String,
    pub len: u32,    // genome length (bp)
    pub l_sel: u32,  // distinct selected k-mers (row "mass" in events)
    pub flags: u32,
}

pub struct HvDb {
    pub params: EncodeParams,
    pub rows: Vec<i8>,   // N × D row-major
    pub scales: Vec<f32>,
    pub metas: Vec<RowMeta>,
    pub norms: Vec<f64>, // ‖aᵢ‖₂ dequantized
}

impl HvDb {
    pub fn n(&self) -> usize {
        self.metas.len()
    }
    pub fn dim(&self) -> usize {
        1usize << self.params.dim_bits
    }

    pub fn row_slice(&self, i: usize) -> &[i8] {
        let d = self.dim();
        &self.rows[i * d..(i + 1) * d]
    }

    /// Dequantized dot product ⟨v, aᵢ⟩.
    pub fn dot(&self, v: &[i32], i: usize) -> f64 {
        let scale = self.scales[i] as f64;
        let row = self.row_slice(i);
        let mut s = 0f64;
        for b in 0..v.len() {
            s += (row[b] as f64) * (v[b] as f64);
        }
        s * scale
    }    /// Dequantized row as f64 (candidates only; N×D floats would be wasteful).
    pub fn row_f64(&self, i: usize) -> Vec<f64> {
        let scale = self.scales[i] as f64;
        self.row_slice(i).iter().map(|&q| q as f64 * scale).collect()
    }

    // ── build ────────────────────────────────────────────────────────────

    /// Build a row from FASTA bytes (possibly gz-decompressed upstream).
    /// Accumulates each DISTINCT selected k-mer once (HashSet dedup).
    pub fn build_row(params: &EncodeParams, fasta: &[u8], meta: RowMeta) -> (Vec<i8>, f32, u32, u32) {
        let mut acc = Accumulator::new(*params);
        let mut seen: HashSet<u64> = HashSet::new();
        let mut rh = acc.contig_hasher();
        let proj = Projector::new(params.dim_bits, params.dens_num, params.dens_den, params.seed);

        let mut line: Vec<u8> = Vec::with_capacity(4096);
        let mut glen: u64 = 0;
        let mut flush = |line: &mut Vec<u8>, rh: &mut RollingHash, acc: &mut Accumulator, seen: &mut HashSet<u64>, glen: &mut u64, proj: &Projector| {
            if line.first() == Some(&b'>') {
                rh.reset();
            } else if !line.is_empty() {
                for &base in line.iter() {
                    let code = crate::hash::code_of(base);
                    if code <= 3 {
                        *glen += 1;
                    }
                    if let Some(u) = rh.push(code) {
                        let (sel, bucket, sign) = proj.derive(u);
                        if sel && seen.insert(u) {
                            let v = acc.acc[bucket].saturating_add(sign);
                            acc.acc[bucket] = v;
                            acc.stats.n_events += 1;
                        }
                    }
                }
            }
            line.clear();
        };

        for &c in fasta {
            if c == b'\n' {
                flush(&mut line, &mut rh, &mut acc, &mut seen, &mut glen, &proj);
            } else if c != b'\r' {
                line.push(c);
            }
        }
        flush(&mut line, &mut rh, &mut acc, &mut seen, &mut glen, &proj);

        let l_sel = acc.stats.n_events as u32;
        let maxabs = acc.acc.iter().map(|&v| v.unsigned_abs()).max().unwrap_or(1) as f64;
        let scale = (127.0 / maxabs.max(1.0)) as f32;
        let quant: Vec<i8> = acc
            .acc
            .iter()
            .map(|&v| ((v as f64 * scale as f64).round().clamp(-127.0, 127.0)) as i8)
            .collect();
        (quant, (1.0 / scale) as f32, l_sel, glen as u32)
        // note: stored scale is the MULTIPLIER (a = q · scale)
    }

    pub fn build(params: &EncodeParams, genomes: Vec<(Vec<u8>, RowMeta)>) -> HvDb {
        let d = params.dim();
        let mut rows = Vec::with_capacity(genomes.len() * d);
        let mut scales = Vec::with_capacity(genomes.len());
        let mut metas = Vec::with_capacity(genomes.len());
        for (fasta, mut meta) in genomes {
            let (q, scale, l_sel, glen) = Self::build_row(params, &fasta, meta.clone());
            meta.l_sel = l_sel;
            meta.len = if meta.len > 0 { meta.len } else { glen };
            rows.extend_from_slice(&q);
            scales.push(scale as f32);
            metas.push(meta);
        }
        let mut db = HvDb {
            params: *params,
            rows,
            scales,
            metas,
            norms: Vec::new(),
        };
        db.recompute_norms();
        db
    }

    pub fn recompute_norms(&mut self) {
        let d = self.dim();
        self.norms = (0..self.n())
            .map(|i| {
                let scale = self.scales[i] as f64;
                let row = self.row_slice(i);
                let s: f64 = row.iter().map(|&q| (q as f64) * (q as f64)).sum();
                (s.sqrt()) * scale
            })
            .collect();
    }

    // ── serialization (.hvd + sidecar tsv) ────────────────────────────────

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(64 + self.rows.len());
        out.extend_from_slice(b"HVP1");
        out.extend_from_slice(&1u32.to_le_bytes());
        out.extend_from_slice(&self.params.k.to_le_bytes());
        out.extend_from_slice(&self.params.dim_bits.to_le_bytes());
        out.extend_from_slice(&self.params.dens_num.to_le_bytes());
        out.extend_from_slice(&self.params.dens_den.to_le_bytes());
        out.extend_from_slice(&self.params.qmin.to_le_bytes());
        out.extend_from_slice(&self.params.seed.to_le_bytes());
        out.extend_from_slice(&(self.n() as u32).to_le_bytes());
        for &q in &self.rows {
            out.push(q as u8);
        }
        for &s in &self.scales {
            out.extend_from_slice(&s.to_le_bytes());
        }
        for m in &self.metas {
            out.extend_from_slice(&m.len.to_le_bytes());
            out.extend_from_slice(&m.l_sel.to_le_bytes());
            out.extend_from_slice(&m.flags.to_le_bytes());
        }
        out
    }

    pub fn from_bytes(bytes: &[u8]) -> Result<HvDb, String> {
        let bad = |m: &str| Err(format!("hvd: {m}"));
        if bytes.len() < 40 || &bytes[0..4] != b"HVP1" {
            return bad("bad magic");
        }
        let rd_u32 = |off: usize| -> u32 {
            u32::from_le_bytes([bytes[off], bytes[off + 1], bytes[off + 2], bytes[off + 3]])
        };
        let params = EncodeParams {
            k: rd_u32(8),
            dim_bits: rd_u32(12),
            dens_num: rd_u32(16),
            dens_den: rd_u32(20),
            qmin: bytes[24],
            seed: u64::from_le_bytes(
                bytes[25..33].try_into().map_err(|_| "hvd: seed".to_string())?,
            ),
        };
        let n = rd_u32(33) as usize;
        let d = params.dim();
        let need = 37 + n * d + n * 4 + n * 12;
        if bytes.len() < need {
            return bad("truncated");
        }
        let rows: Vec<i8> = bytes[37..37 + n * d].to_vec().into_iter().map(|b| b as i8).collect();
        let mut off = 37 + n * d;
        let mut scales = Vec::with_capacity(n);
        for _ in 0..n {
            scales.push(f32::from_le_bytes(bytes[off..off + 4].try_into().unwrap()));
            off += 4;
        }
        let mut metas = Vec::with_capacity(n);
        for i in 0..n {
            let len = u32::from_le_bytes(bytes[off..off + 4].try_into().unwrap());
            let l_sel = u32::from_le_bytes(bytes[off + 4..off + 8].try_into().unwrap());
            let flags = u32::from_le_bytes(bytes[off + 8..off + 12].try_into().unwrap());
            off += 12;
            metas.push(RowMeta {
                id: format!("row{}", i),
                name: format!("row{}", i),
                genus: String::new(),
                len,
                l_sel,
                flags,
            });
        }
        let mut db = HvDb {
            params,
            rows,
            scales,
            metas,
            norms: Vec::new(),
        };
        db.recompute_norms();
        Ok(db)
    }

    /// Attach sidecar metadata (idx, id, name, genus, len, l_sel, flags) —
    /// the exact format emitted by `meta_tsv` — by row order.
    pub fn attach_meta_tsv(&mut self, tsv: &[u8]) -> Result<(), String> {
        let text = String::from_utf8_lossy(tsv);
        for (i, line) in text.lines().enumerate() {
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }
            if i >= self.n() {
                break;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 4 {
                continue;
            }
            let id = f[1].trim();
            let name = f[2].trim();
            let genus = f[3].trim();
            let flags = f.get(6).copied().unwrap_or("0").trim().parse::<u32>().unwrap_or(0);
            if !id.is_empty() {
                self.metas[i].id = id.to_string();
                self.metas[i].name = if name.is_empty() { id.to_string() } else { name.to_string() };
                self.metas[i].genus = genus.to_string();
                self.metas[i].flags |= flags;
            }
        }
        Ok(())
    }

    pub fn meta_tsv(&self) -> Vec<u8> {
        let mut out = String::new();
        for (i, m) in self.metas.iter().enumerate() {
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                i, m.id, m.name, m.genus, m.len, m.l_sel, m.flags
            ));
        }
        out.into_bytes()
    }
}

/// Deterministic synthetic spike-in genome (random 100 kb) — the "standard
/// addition" QC row (PLAN.md §1.2). Two builds with the same seed must yield
/// the identical sequence, so browser and native DBs match.
pub fn spike_genome(seed: u64, len: usize) -> Vec<u8> {
    let mut rng = splitmix64(seed ^ 0x5B1CE).wrapping_mul(0x9E3779B97F4A7C15);
    let mut out = Vec::with_capacity(len + 2);
    out.extend_from_slice(b">hypv_spike_in\n");
    for _ in 0..len {
        rng ^= rng << 13;
        rng ^= rng >> 7;
        rng ^= rng << 17;
        out.push(b"ACGT"[(rng % 4) as usize]);
    }
    out.push(b'\n');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toy_fasta(seed: u64, len: usize) -> Vec<u8> {
        spike_genome(seed, len)
    }

    #[test]
    fn row_quantization_roundtrip() {
        let params = EncodeParams::default_params();
        let (_q, scale, l_sel, glen) = HvDb::build_row(&params, &toy_fasta(7, 200_000), RowMeta::default());
        assert_eq!(glen, 200_000);
        let expect_sel = (200_000.0 * 0.05) as f64;
        assert!(
            (l_sel as f64 - expect_sel).abs() / expect_sel < 0.05,
            "l_sel {l_sel} vs {expect_sel}"
        );
        assert!(scale > 0.0);
    }

    #[test]
    fn serialization_roundtrip() {
        let params = EncodeParams::default_params();
        let genomes = vec![
            (toy_fasta(1, 50_000), RowMeta { id: "g1".into(), name: "Genus one sp. A".into(), genus: "Genus".into(), ..Default::default() }),
            (toy_fasta(2, 50_000), RowMeta { id: "g2".into(), name: "Genus one sp. B".into(), genus: "Genus".into(), ..Default::default() }),
        ];
        let db = HvDb::build(&params, genomes);
        let bytes = db.to_bytes();
        let db2 = HvDb::from_bytes(&bytes).unwrap();
        assert_eq!(db2.rows, db.rows);
        assert_eq!(db2.scales, db.scales);
        assert_eq!(db2.metas.len(), 2);
        assert_eq!(db2.norms, db.norms);
    }
}
