//! Deconvolution: screen → joint NNLS → depletion → calibrated report.
//!
//! Structure (PLAN.md §4.3–§4.6), with the CLEAN / successive-interference-
//! cancellation interpretation: each round subtracts what the current support
//! explains and re-tests the residual, the standard remedy for the near–far
//! problem (a dominant genome masking a weak pathogen).
//!
//! NOTE on identity: a single occurrence-mode vector cannot distinguish
//! "exact row at coverage μ·aᵏ" from "relative at ANI a, coverage μ" — the
//! two leave the same projection. Disambiguation needs the distinct/occurrence
//! dual accumulator (PLAN §1.6, M1). v0 reports detection + mass only.

use crate::db::{HvDb, FLAG_SPIKE};
use crate::stats::{bh_qvalues, gumbel_expected_max, mad_sigma, p_one_sided};
use std::collections::HashMap;

pub const FLAG_GENUS_COMPLEX: u32 = 1 << 8;
pub const FLAG_LOW_IDENTITY: u32 = 1 << 9; // reserved for M1 dual-mode disambiguation

#[derive(Clone, Debug)]
pub struct TaxonResult {
    pub idx: usize,
    pub x: f64, // fitted coverage μ̂ (occurrences per reference copy)
    pub z: f64, // detection z: ⟨v, aᵢ⟩ / (σ̂·‖aᵢ‖)
    pub q: f64, // BH q-value over the screened family
    pub dot: f64,
    pub dna_frac: f64,
    pub flags: u32,
}

#[derive(Clone, Debug)]
pub struct ProfileReport {
    pub taxa: Vec<TaxonResult>,
    pub sigma: f64,       // final residual noise scale (per bucket)
    pub unexplained: f64, // 1 − explained events / observed events ∈ [0,1]
    pub spike_z: f64,     // null check on the spike-in row (≈ N(0,1))
    pub screen_max: f64,  // top raw z at round 0
    pub ev_bound: f64,    // expected max null z for N rows (Gumbel)
    pub n_events: u64,
    pub rounds: usize,
}

#[derive(Clone, Copy)]
pub struct ProfileParams {
    pub zmin: f64,
    pub kmax: usize,
    pub rounds: usize,
    pub genus_collapse_cos: f64,
}

impl ProfileParams {
    pub fn default_params() -> Self {
        ProfileParams {
            zmin: 5.0,
            kmax: 96,
            rounds: 3,
            genus_collapse_cos: 0.5,
        }
    }
}

impl Default for ProfileParams {
    fn default() -> Self {
        Self::default_params()
    }
}

/// profile() accepts any integer vector obeying the accumulator contract
/// (native accumulator or a dequantized fingerprint — both work).
/// `n_events` is the encoder's selected-occurrence count (stored in the
/// fingerprint header for fingerprints).
pub fn profile(v: &[i32], n_events: u64, db: &HvDb, pp: ProfileParams) -> ProfileReport {
    let n = db.n();
    let d = db.dim();
    assert_eq!(v.len(), d, "accumulator dim mismatch");

    let vf: Vec<f64> = v.iter().map(|&x| x as f64).collect();
    let sigma0 = mad_sigma(v);
    let ev_bound = gumbel_expected_max(n);

    // ── screen: matched filter against every row ─────────────────────────
    let dots: Vec<f64> = (0..n).map(|i| db.dot(v, i)).collect();
    let z_of = |i: usize, dot: f64, sigma: f64| dot / (sigma * db.norms[i]).max(1e-12);
    let screen_max = (0..n).map(|i| z_of(i, dots[i], sigma0)).fold(0.0f64, f64::max);

    // ── support evolution (depletion loop) ───────────────────────────────
    let mut support: Vec<usize> = Vec::new();
    let mut sigma = sigma0;
    let mut resid = vf.clone();
    let mut rounds_used = 0usize;

    for round in 0..pp.rounds {
        // z of every not-yet-supported row against the current residual.
        // Ordering uses the depleted residual (prioritises what the fit has
        // not yet explained); ADMISSION is gated on the σ₀-calibrated z so
        // the shrinkage of σ over rounds cannot progressively loosen the
        // entry threshold (calibration invariant across rounds).
        let mut zr: Vec<(usize, f64)> = Vec::with_capacity(n);
        for i in 0..n {
            if support.contains(&i) {
                continue;
            }
            let mut dot = 0.0;
            let row = db.row_slice(i);
            for b in 0..d {
                dot += row[b] as f64 * resid[b];
            }
            let dot = dot * db.scales[i] as f64;
            let z_rank = z_of(i, dot, sigma);
            let z_admit = z_of(i, dot, sigma0);
            if z_rank > 0.0 && z_admit > 0.0 {
                zr.push((i, z_admit));
            }
        }
        zr.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

        let admit_thr = if round == 0 { 2.5 } else { pp.zmin };
        let mut admitted = 0usize;
        for &(i, z) in &zr {
            if z < admit_thr || support.len() >= pp.kmax {
                break;
            }
            support.push(i);
            admitted += 1;
        }
        if admitted == 0 {
            break;
        }
        rounds_used = round + 1;

        // joint NNLS over the support, then rebuild the residual
        let x = nnls_support(&vf, db, &support);
        for b in 0..d {
            resid[b] = vf[b];
        }
        for (k, &i) in support.iter().enumerate() {
            let scale = db.scales[i] as f64;
            let row = db.row_slice(i);
            for b in 0..d {
                resid[b] -= x[k] * row[b] as f64 * scale;
            }
        }
        let resid_i32: Vec<i32> = resid.iter().map(|&r| r as i32).collect();
        sigma = mad_sigma(&resid_i32).max(sigma0 * 0.05);
    }

    // final fit on the settled support
    let x = nnls_support(&vf, db, &support);

    // ── significance over the full screened family (BH) ──────────────────
    // Detection z must be calibrated under the row-absent null: a raw dot's
    // variance is σ₀²‖aᵢ‖², where σ₀ is the PRE-depletion noise scale. Using
    // the depleted σ here would inflate every absent row's z by σ₀/σ_resid
    // (~2–3× on isolate-like samples) — anti-conservative exactly where FDR
    // matters. The depleted σ stays in the report as a QC number only.
    let mut p_all: Vec<f64> = Vec::with_capacity(n);
    for i in 0..n {
        let z = z_of(i, dots[i], sigma0);
        p_all.push(if z > 0.0 { p_one_sided(z) } else { 1.0 });
    }
    let q_all = bh_qvalues(&p_all);

    // ── assemble ─────────────────────────────────────────────────────────
    let mut taxa: Vec<TaxonResult> = Vec::new();
    let mut mass_total = 0.0f64;
    for (k, &i) in support.iter().enumerate() {
        if x[k] <= 0.0 {
            continue;
        }
        let z = z_of(i, dots[i], sigma0);
        if z < pp.zmin {
            continue;
        }
        mass_total += x[k] * db.metas[i].len as f64;
        taxa.push(TaxonResult {
            idx: i,
            x: x[k],
            z,
            q: q_all[i],
            dot: dots[i],
            dna_frac: 0.0,
            flags: db.metas[i].flags,
        });
    }
    for t in taxa.iter_mut() {
        if mass_total > 0.0 {
            t.dna_frac = t.x * db.metas[t.idx].len as f64 / mass_total;
        }
    }

    // genus collapse (identifiability within a coherent sub-block)
    taxa = collapse_genera(taxa, db, pp.genus_collapse_cos);

    // ── report-level numbers ────────────────────────────────────────────
    let mut explained_events = 0.0f64;
    for t in &taxa {
        explained_events += t.x * db.metas[t.idx].l_sel as f64;
    }
    let spike_idx = (0..n).find(|&i| db.metas[i].flags & FLAG_SPIKE != 0);
    let spike_z = spike_idx
        .map(|i| z_of(i, dots[i], sigma0))
        .unwrap_or(0.0);

    ProfileReport {
        taxa,
        sigma,
        unexplained: (1.0 - explained_events / n_events.max(1) as f64).clamp(0.0, 1.0),
        spike_z,
        screen_max,
        ev_bound,
        n_events,
        rounds: rounds_used,
    }
}

/// Lawson–Hanson NNLS restricted to `support`, on the Gram system.
/// Columns are unit-normalized internally; x returns in coverage units.
pub fn nnls_support(v: &[f64], db: &HvDb, support: &[usize]) -> Vec<f64> {
    let m = support.len();
    // A: m × d dequantized, columns normalized
    let mut a: Vec<Vec<f64>> = Vec::with_capacity(m);
    let mut nrm = Vec::with_capacity(m);
    for &i in support {
        let row = db.row_f64(i);
        let nn = row.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-12);
        nrm.push(nn);
        a.push(row.iter().map(|&x| x / nn).collect());
    }
    let b: Vec<f64> = a.iter().map(|r| dot(r, v)).collect();
    let mut g = vec![vec![0.0f64; m]; m];
    for j in 0..m {
        for k2 in j..m {
            let s = dot(&a[j], &a[k2]);
            g[j][k2] = s;
            g[k2][j] = s;
        }
    }
    let y = nnls(&g, &b);
    y.iter()
        .zip(nrm.iter())
        .map(|(&yy, &nn)| yy / nn)
        .collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Lawson–Hanson active-set NNLS on the Gram system G y = b, y ≥ 0.
fn nnls(g: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0f64; n];
    let mut p = vec![false; n];
    let tol = 1e-12 * (1.0 + b.iter().map(|v| v.abs()).sum::<f64>() / (n.max(1) as f64));

    for _ in 0..8 * n + 16 {
        // w = b − G x
        let w: Vec<f64> = (0..n)
            .map(|j| b[j] - (0..n).map(|l| g[j][l] * x[l]).sum::<f64>())
            .collect();
        let mut jmax = usize::MAX;
        let mut wmax = tol;
        for j in 0..n {
            if !p[j] && w[j] > wmax {
                wmax = w[j];
                jmax = j;
            }
        }
        if jmax == usize::MAX {
            break;
        }
        p[jmax] = true;

        loop {
            let idx: Vec<usize> = (0..n).filter(|&j| p[j]).collect();
            let sub_b: Vec<f64> = idx.iter().map(|&j| b[j]).collect();
            let sol = match solve_spd_sub(g, &idx, &sub_b) {
                Some(s) => s,
                None => {
                    p[jmax] = false;
                    break;
                }
            };
            if sol.iter().all(|&s| s > 0.0) {
                for (&j, &s) in idx.iter().zip(sol.iter()) {
                    x[j] = s;
                }
                break;
            }
            // interpolate to zero out the most negative coordinate
            let mut alpha = f64::INFINITY;
            let mut jout = idx[0];
            for (t, &j) in idx.iter().enumerate() {
                if sol[t] <= 0.0 {
                    let denom = x[j] - sol[t];
                    let a = if denom.abs() < 1e-300 { 0.0 } else { x[j] / denom };
                    if a < alpha {
                        alpha = a;
                        jout = j;
                    }
                }
            }
            for (t, &j) in idx.iter().enumerate() {
                x[j] += alpha * (sol[t] - x[j]);
            }
            x[jout] = 0.0;
            p[jout] = false;
        }
    }
    x
}

/// Cholesky solve of G[idx][idx] s = b_sub, with a jitter fallback.
fn solve_spd_sub(g: &[Vec<f64>], idx: &[usize], b: &[f64]) -> Option<Vec<f64>> {
    let m = idx.len();
    let mut l = vec![vec![0.0f64; m]; m];
    for jitter_pass in 0..2 {
        let jit = if jitter_pass == 0 {
            0.0
        } else {
            1e-9 * (1.0 + g[idx[0]][idx[0]])
        };
        let mut ok = true;
        for i in 0..m {
            for j in 0..=i {
                let mut s = g[idx[i]][idx[j]];
                if i == j {
                    s += jit;
                }
                for k in 0..j {
                    s -= l[i][k] * l[j][k];
                }
                if i == j {
                    if s <= 0.0 {
                        ok = false;
                        break;
                    }
                    l[i][j] = s.sqrt();
                } else {
                    l[i][j] = s / l[j][j];
                }
            }
            if !ok {
                break;
            }
        }
        if ok {
            let mut y = vec![0.0f64; m];
            for i in 0..m {
                let mut s = b[i];
                for k in 0..i {
                    s -= l[i][k] * y[k];
                }
                y[i] = s / l[i][i];
            }
            let mut s = vec![0.0f64; m];
            for i in (0..m).rev() {
                let mut t = y[i];
                for k in (i + 1)..m {
                    t -= l[k][i] * s[k];
                }
                s[i] = t / l[i][i];
            }
            return Some(s);
        }
    }
    None
}

/// Collapse same-genus species with highly correlated rows into one entry
/// (the "blend within a synthesised beam" rule — report the complex).
fn collapse_genera(taxa: Vec<TaxonResult>, db: &HvDb, cos_thr: f64) -> Vec<TaxonResult> {
    let mut out: Vec<TaxonResult> = Vec::new();
    let mut groups: HashMap<String, Vec<TaxonResult>> = HashMap::new();
    let mut singles: Vec<TaxonResult> = Vec::new();
    for t in taxa {
        let genus = db.metas[t.idx].genus.clone();
        if genus.is_empty() {
            singles.push(t);
        } else {
            groups.entry(genus).or_default().push(t);
        }
    }
    for (_, mut members) in groups {
        if members.len() == 1 {
            out.append(&mut members);
            continue;
        }
        let mut collapse = false;
        'outer: for i in 0..members.len() {
            for j in (i + 1)..members.len() {
                let ri = db.row_f64(members[i].idx);
                let rj = db.row_f64(members[j].idx);
                let cos =
                    dot(&ri, &rj) / (db.norms[members[i].idx] * db.norms[members[j].idx]);
                if cos > cos_thr {
                    collapse = true;
                    break 'outer;
                }
            }
        }
        if collapse {
            let x: f64 = members.iter().map(|t| t.x).sum();
            let idx0 = members[0].idx;
            let dot_sum: f64 = members.iter().map(|t| t.dot).sum();
            let z = members.iter().map(|t| t.z).fold(0.0f64, f64::max);
            let q = members.iter().map(|t| t.q).fold(1.0f64, f64::min);
            let dna_frac: f64 = members.iter().map(|t| t.dna_frac).sum();
            members.clear();
            out.push(TaxonResult {
                idx: idx0,
                x,
                z,
                q,
                dot: dot_sum,
                dna_frac,
                flags: FLAG_GENUS_COMPLEX,
            });
        } else {
            out.append(&mut members);
        }
    }
    out.extend(singles);
    out.sort_by(|a, b| b.x.partial_cmp(&a.x).unwrap());
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nnls_recovers_exact_nonnegative_mixture() {
        // three ~orthogonal "genomes" in D=8 dims
        let g = vec![
            vec![1.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0, 0.0, 0.05, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.02, 0.0],
        ];
        let b = vec![2.0, 0.5, 0.0]; // third absent
        let x = nnls(&g, &b);
        assert!((x[0] - 2.0).abs() < 1e-9);
        assert!((x[1] - 0.5).abs() < 1e-9);
        assert!(x[2] <= 0.0);
    }
}
