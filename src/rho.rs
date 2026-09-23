//! Probability that a curve (or P-1) finds a factor of a given size, from GMP-ECM's `rho.c`
//! (without Brent-Suyama's extension, which this crate does not use).
//!
//! A curve finds a prime `p` when its group order is `b1`-smooth but for at most one prime
//! factor in `(b1, b2]`. Group orders behave like random numbers around `p/e^delta` (they are
//! more often smooth, thanks to their known torsion), whose smoothness is given by Dickman's
//! rho function (in its "local" form, for numbers near `x`, from Alexander Kruppa's PhD thesis,
//! equations (5.6) and (5.10)).

use primal::Primes;
use std::sync::OnceLock;

/// Extra smoothness of the group orders of Suyama's curves (and GMP-ECM's parametrization 2,
/// which have the same torsion) with respect to random numbers (`ECM_EXTRA_SMOOTHNESS`).
const ECM_EXTRA_SMOOTHNESS: f64 = 3.134;

/// Extra smoothness of `p - 1` with respect to random numbers: `sum_q log(q)/(q - 1)^2`.
#[cfg_attr(not(any(test, feature = "bench")), allow(dead_code))]
const PM1_EXTRA_SMOOTHNESS: f64 = 1.2269688;

/// Euler-Mascheroni constant.
const EULER: f64 = 0.577_215_664_901_532_9;

/// Table steps per unit of `rho`'s argument.
const INV_H: usize = 256;
const H: f64 = 1.0 / INV_H as f64;
/// `rho` is tabulated on `[0, TABLE_MAX)`, and taken as 0 beyond.
const TABLE_MAX: usize = 10;

/// Below this bound, stage 2 primes are summed one by one.
const SUM_THRESHOLD: usize = 20_000;

/// Probability that a curve with bounds `b1` and `b2` finds a given prime factor of `digits`
/// decimal digits (taken as `10^(digits - 1/2)`).
pub fn ecm_prob(b1: f64, b2: f64, digits: f64) -> f64 {
    prob(b1, b2, 10f64.powf(digits - 0.5), ECM_EXTRA_SMOOTHNESS)
}

/// Probability that P-1 with bounds `b1` and `b2` finds a given prime factor of `digits`
/// decimal digits.
#[cfg_attr(not(any(test, feature = "bench")), allow(dead_code))]
pub fn pm1_prob(b1: f64, b2: f64, digits: f64) -> f64 {
    prob(b1, b2, 10f64.powf(digits - 0.5), PM1_EXTRA_SMOOTHNESS)
}

/// Probability that a number around `n/e^delta` is `b1`-smooth but for at most one prime factor
/// in `(b1, b2]`.
fn prob(b1: f64, b2: f64, n: f64, delta: f64) -> f64 {
    if b1 < 2.0 || n <= 1.0 {
        return 0.0;
    }
    let n = n / delta.exp();
    if n <= b1 {
        return 1.0;
    }
    let alpha = n.ln() / b1.ln();
    let stage1 = local(alpha, n);
    let mut stage2 = 0.0;
    if b2 > b1 {
        let beta = if (b1 as usize) < SUM_THRESHOLD {
            let top = b2.min(SUM_THRESHOLD as f64);
            stage2 += mu_sum(b1 as usize, top as usize, n);
            b2.ln() / top.ln()
        } else {
            b2.ln() / b1.ln()
        };
        if beta > 1.0 {
            stage2 += mu(alpha, beta, n);
        }
    }
    (stage1 + stage2).max(0.0)
}

/// `rho(i*H)` for `i < INV_H*TABLE_MAX`.
fn table() -> &'static [f64] {
    static TABLE: OnceLock<Vec<f64>> = OnceLock::new();
    TABLE.get_or_init(|| {
        let len = INV_H * TABLE_MAX;
        let mut rho = vec![0.0; len];
        for (i, r) in rho.iter_mut().enumerate().take(3 * INV_H) {
            *r = rho_exact(i as f64 * H);
        }
        // rho(i*h) = rho((i - 4)*h) - int_{(i-4)*h}^{i*h} rho(x - 1)/x dx (Boole's rule).
        for i in 3 * INV_H..len {
            let f = |k: usize| rho[i - INV_H - k] / (i - k) as f64;
            rho[i] = (rho[i - 4]
                - 2.0 / 45.0 * (7.0 * f(4) + 32.0 * f(3) + 12.0 * f(2) + 32.0 * f(1) + 7.0 * f(0)))
            .max(0.0);
        }
        rho
    })
}

/// Dilogarithm `Li2(z) = sum z^k/k^2`, for `|z| <= 1/2`.
fn dilog_series(z: f64) -> f64 {
    let (mut r, mut zk) = (0.0, z);
    for k in 1..=44 {
        r += zk / (k * k) as f64;
        zk *= z;
    }
    r
}

/// Dilogarithm, for `x <= -1`.
fn dilog(x: f64) -> f64 {
    const PI_SQR_6: f64 = std::f64::consts::PI * std::f64::consts::PI / 6.0;
    if x <= -2.0 {
        let l = (-1.0 / x).ln();
        -dilog_series(1.0 / x) - PI_SQR_6 - 0.5 * l * l
    } else {
        let log1x = (1.0 - x).ln();
        dilog_series(1.0 / (1.0 - x)) - PI_SQR_6 + log1x * (0.5 * log1x - (-x).ln())
    }
}

/// Dickman's rho, exact for `x <= 3`.
fn rho_exact(x: f64) -> f64 {
    const PI_SQR_12: f64 = std::f64::consts::PI * std::f64::consts::PI / 12.0;
    if x <= 0.0 {
        0.0
    } else if x <= 1.0 {
        1.0
    } else if x <= 2.0 {
        1.0 - x.ln()
    } else {
        1.0 - x.ln() * (1.0 - (x - 1.0).ln()) + dilog(1.0 - x) + PI_SQR_12
    }
}

/// Dickman's rho, for `alpha < TABLE_MAX` (linear interpolation of the table).
fn rho(alpha: f64) -> f64 {
    if alpha <= 3.0 {
        return rho_exact(alpha);
    }
    let t = table();
    let a = (alpha * INV_H as f64).floor() as usize;
    let rho1 = t[a];
    let rho2 = t.get(a + 1).copied().unwrap_or(0.0);
    rho1 + (rho2 - rho1) * (alpha * INV_H as f64 - a as f64)
}

/// Probability that a number near `x` is `x^(1/alpha)`-smooth.
fn local(alpha: f64, x: f64) -> f64 {
    if alpha <= 1.0 {
        rho_exact(alpha)
    } else if alpha < TABLE_MAX as f64 {
        rho(alpha) - EULER * rho(alpha - 1.0) / x.ln()
    } else {
        0.0
    }
}

/// [`local`] at `alpha = ai*H`.
fn local_i(ai: usize, x: f64) -> f64 {
    let t = table();
    if ai == 0 {
        0.0
    } else if ai <= INV_H {
        1.0
    } else if ai >= t.len() {
        0.0
    } else if ai <= 2 * INV_H {
        t[ai] - EULER / x.ln()
    } else {
        let logx = x.ln();
        t[ai] - (EULER * t[ai - INV_H] + (1.0 - EULER) * t[ai - 2 * INV_H] / logx) / logx
    }
}

/// Probability that a number near `x` is `b1`-smooth but for one prime in `(b1, b2]`, summed
/// over these primes.
fn mu_sum(b1: usize, b2: usize, x: f64) -> f64 {
    let (inv_log_b1, logx) = (1.0 / (b1 as f64).ln(), x.ln());
    Primes::all()
        .skip_while(|&p| p <= b1)
        .take_while(|&p| p <= b2)
        .map(|p| {
            let p = p as f64;
            local((logx - p.ln()) * inv_log_b1, x / p) / p
        })
        .sum()
}

/// Probability that a number near `x` has its second largest prime factor below `x^(1/alpha)`
/// and its largest one below `x^(beta/alpha)` (trapezoidal rule on the table).
fn mu(alpha: f64, beta: f64, x: f64) -> f64 {
    let max = (TABLE_MAX * INV_H) as f64;
    let ai = ((alpha - beta) * INV_H as f64).ceil().min(max).max(0.0) as usize;
    let bi = ((alpha - 1.0) * INV_H as f64).floor().min(max).max(0.0) as usize;
    let (a, b) = (ai as f64 * H, bi as f64 * H);
    let mut sum = 0.0;
    for i in ai + 1..bi {
        sum += local_i(i, x) / (alpha - i as f64 * H);
    }
    sum += 0.5 * local_i(ai, x) / (alpha - a);
    sum += 0.5 * local_i(bi, x) / (alpha - b);
    sum *= H;
    sum +=
        (a - alpha + beta) * 0.5 * (local_i(ai, x) / (alpha - a) + local(alpha - beta, x) / beta);
    sum += (alpha - 1.0 - b) * 0.5 * (local(alpha - 1.0, x) + local_i(bi, x) / (alpha - b));
    sum
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rho_values() {
        // Known values of Dickman's function.
        assert!((rho(2.5) - 0.130_319).abs() < 1e-5);
        assert!((rho(4.0) - 0.004_910_9).abs() < 1e-6);
        assert!((rho(5.0) - 0.000_354_7).abs() < 1e-7);
    }

    #[test]
    fn matches_gmp_ecm() {
        // Expected curves printed by GMP-ECM 7.0.7: `ecm -v -power 1 -param 2 B1` (default B2).
        for (b1, b2, digits, curves) in [
            (2e3, 147_396.0, 35.0, 2.2e9),
            (11e3, 1_873_422.0, 35.0, 2_924_742.0),
            (11e3, 1_873_422.0, 40.0, 1.8e8),
            (5e4, 12_746_592.0, 35.0, 71_823.0),
            (25e4, 128_992_510.0, 35.0, 5221.0),
            (25e4, 128_992_510.0, 45.0, 1_306_906.0),
            (1e6, 1_045_563_762.0, 35.0, 986.0),
            (1e6, 1_045_563_762.0, 50.0, 1_415_722.0),
            (3e6, 5_706_890_290.0, 40.0, 2547.0),
            (3e6, 5_706_890_290.0, 70.0, 6.7e9),
        ] {
            let expected = 1.0 / ecm_prob(b1, b2, digits);
            // Large values are printed with 2 significant digits.
            assert!(
                (expected / curves - 1.0).abs() < 0.03,
                "{b1} {digits}: {expected}"
            );
        }
        // `ecm -pm1 -v 1e6`: B2 = 1748900148.
        for (digits, p) in [(20.0, 0.18), (25.0, 0.039), (30.0, 0.0063), (35.0, 0.00078)] {
            let prob = pm1_prob(1e6, 1_748_900_148.0, digits);
            assert!((prob / p - 1.0).abs() < 0.05, "{digits}: {prob}");
        }
    }
}
