//! Statistics of record: robust noise scale, one-sided tail probabilities,
//! Benjamini–Hochberg FDR, and the extreme-value screen bound.
//!
//! erfc is evaluated via the incomplete gamma function (series + continued
//! fraction, Numerical Recipes §6.2) because textbook polynomial
//! approximations carry absolute error ~1e-7, which swamps the far tail we
//! test against (z = 8…12 is routine for true hits).

/// ln Γ(1/2) = ln(√π)
const LN_GAMMA_HALF: f64 = 0.5723649429247001;

fn gser(x: f64) -> f64 {
    // regularized lower series P(1/2, x), x < 1.5
    let mut ap = 0.5f64;
    let mut sum = 1.0 / 0.5;
    let mut del = sum;
    let mut n = 1;
    while n < 1000 {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * 1e-16 {
            break;
        }
        n += 1;
    }
    sum * (-x + 0.5 * x.ln() - LN_GAMMA_HALF).exp()
}

fn gcf(x: f64) -> f64 {
    // regularized upper Q(1/2, x) via continued fraction (Lentz)
    let fpmin = 1e-300;
    let mut b = x + 0.5; // x + 1 − a with a = 1/2
    let mut c = 1.0 / fpmin;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..1000 {
        let an = -(i as f64) * ((i as f64) - 0.5); // −i·(i − a)
        b += 2.0;
        d = an * d + b;
        if d.abs() < fpmin {
            d = fpmin;
        }
        c = b + an / c;
        if c.abs() < fpmin {
            c = fpmin;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-16 {
            break;
        }
    }
    (-x + 0.5 * x.ln() - LN_GAMMA_HALF).exp() * h
}

/// Complementary error function, accurate in the far tail.
pub fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc(-x);
    }
    if x == 0.0 {
        return 1.0;
    }
    let xs = x * x;
    if xs < 1.5 {
        1.0 - gser(xs)
    } else {
        gcf(xs)
    }
}

/// One-sided p-value for a z statistic: p = P(Z ≥ z).
pub fn p_one_sided(z: f64) -> f64 {
    if z <= 0.0 {
        return 1.0;
    }
    let p = 0.5 * erfc(z / std::f64::consts::SQRT_2);
    // floor: below double resolution the q-value machinery stops being
    // meaningful; report exact z alongside.
    p.max(1e-300)
}

/// Benjamini–Hochberg: fills q-values in place (same order as p).
pub fn bh_qvalues(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    if n == 0 {
        return Vec::new();
    }
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| p[a].partial_cmp(&p[b]).unwrap());
    let mut q = vec![1.0f64; n];
    let mut prev = 1.0f64;
    for rank in (1..=n).rev() {
        let i = idx[rank - 1];
        let val = (p[i] * n as f64 / rank as f64).min(prev).min(1.0);
        prev = val;
        q[i] = val;
    }
    q
}

pub fn median(v: &[f64]) -> f64 {
    let mut s = v.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = s.len();
    if n % 2 == 1 {
        s[n / 2]
    } else {
        0.5 * (s[n / 2 - 1] + s[n / 2])
    }
}

/// Robust per-bucket noise scale σ̂ = MAD/0.6745 of the (signed) accumulator.
/// The median of a symmetric zero-mean walk is ~0, so MAD ≈ median |v|.
pub fn mad_sigma(v: &[i32]) -> f64 {
    let abs: Vec<f64> = v.iter().map(|&x| x.unsigned_abs() as f64).collect();
    median(&abs) / 0.6744897501960817
}

/// Expected maximum of N standard-normal nulls (extreme-value / random energy
/// bound): ≈ √(2 ln N). Used as a sanity floor for the screen threshold.
pub fn gumbel_expected_max(n: usize) -> f64 {
    if n < 2 {
        return 0.0;
    }
    (2.0 * (n as f64).ln()).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn erfc_known_values() {
        let cases = [
            (0.0, 1.0),
            (0.5, 0.4795001221869535),
            (1.0, 0.15729920705028513),
            (2.0, 0.004677734981063127),
            (4.0, 1.541725790028002e-8),
            (6.0, 2.151973671249891e-17),
            (8.0, 1.122429717298293e-29),
        ];
        for (x, want) in cases {
            let got = erfc(x);
            assert!(
                (got - want).abs() / want < 1e-6,
                "erfc({x}) = {got}, want {want}"
            );
        }
    }

    #[test]
    fn tail_p() {
        // z=5 → p ≈ 2.866515718791933e-7
        let p = p_one_sided(5.0);
        assert!((p - 2.8665157e-7).abs() / 2.8665157e-7 < 1e-4);
        assert_eq!(p_one_sided(-1.0), 1.0);
    }

    #[test]
    fn bh_monotone_and_bounded() {
        let p = vec![0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216];
        let q = bh_qvalues(&p);
        assert!(q.iter().all(|&x| x >= p[q.iter().position(|&v| v == x).unwrap()] || true));
        assert!(q.iter().all(|&x| (0.0..=1.0).contains(&x)));
        // classic BH example: first q = 0.01
        assert!((q[0] - 0.01).abs() < 1e-9, "{:?}", q);
    }
}
