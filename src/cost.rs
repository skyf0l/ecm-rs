//! Cost model of stage 2, to choose between the continuations and their parameters.
//!
//! Costs are in nanoseconds, built from the measured costs of the primitive operations (on an
//! Intel i7-8750H): modular multiplications by size of `n`, and GMP products by size. Only
//! their ratios matter, which vary much less from one machine to another. The polynomial
//! costs mirror the operations of [`crate::poly`] one product at a time.

use crate::stage2_poly::PolyPlan;

/// Montgomery multiplication, by number of limbs (index 0 unused).
const MUL_NS: [f64; 17] = [
    0.0, 4.7, 9.8, 16.1, 25.8, 43.0, 58.6, 76.5, 101.6, 121.3, 148.8, 177.3, 201.8, 248.5, 275.0,
    317.5, 349.2,
];

/// Product of two `2^(10 + i)`-bit integers by GMP.
const GMP_MUL_NS: [f64; 16] = [
    179.0,
    578.0,
    1764.0,
    5430.0,
    15282.0,
    41225.0,
    109195.0,
    320103.0,
    748329.0,
    1635734.0,
    3395333.0,
    8034959.0,
    16935300.0,
    40792366.0,
    85678918.0,
    200213750.0,
];

/// Products with a factor of at most this many coefficients are schoolbook (see
/// [`crate::poly`]).
use crate::poly::SCHOOLBOOK;

/// Costs of the operations modulo a number of `bits` bits.
#[derive(Debug, Clone, Copy)]
pub struct Costs {
    bits: usize,
    /// Modular multiplication.
    mul: f64,
    /// Product on limbs, accumulated (half a modular multiplication).
    macc: f64,
    /// Reduction of a coefficient of a product.
    redc: f64,
}

impl Costs {
    /// Size of a residue in memory.
    pub fn elem_bytes(&self) -> f64 {
        (self.bits.div_ceil(64).max(1) * 8) as f64
    }

    pub fn new(bits: usize) -> Self {
        let limbs = bits.div_ceil(64).max(1);
        let mul = if limbs < MUL_NS.len() {
            MUL_NS[limbs]
        } else {
            // Plain arithmetic: roughly quadratic, and slower than Montgomery's.
            2.0 * MUL_NS[16] * (limbs as f64 / 16.0).powf(1.8)
        };
        Costs {
            bits,
            mul,
            macc: 0.6 * mul,
            redc: 0.8 * mul,
        }
    }

    /// GMP product of two `bits`-bit integers.
    fn gmp(bits: f64) -> f64 {
        let x = bits.max(1.0).log2() - 10.0;
        if x <= 0.0 {
            // Mostly call overhead below 1024 bits.
            return 20.0 + (GMP_MUL_NS[0] - 20.0) * (bits / 1024.0).powf(1.5);
        }
        let last = GMP_MUL_NS.len() - 1;
        if x >= last as f64 {
            // n log n beyond the table.
            let big = GMP_MUL_NS[last];
            let scale = 2f64.powf(x - last as f64);
            return big * scale * (x + 10.0) / (last as f64 + 10.0);
        }
        let i = x.floor() as usize;
        let f = x - i as f64;
        (GMP_MUL_NS[i].ln() * (1.0 - f) + GMP_MUL_NS[i + 1].ln() * f).exp()
    }

    /// Product of polynomials with `la` and `lb` coefficients.
    pub fn product(&self, la: usize, lb: usize) -> f64 {
        let (lo, hi) = (la.min(lb), la.max(lb));
        if lo == 0 {
            return 0.0;
        }
        let out = (la + lb - 1) as f64;
        if lo <= SCHOOLBOOK {
            return (la * lb) as f64 * self.macc + out * self.redc;
        }
        let slot = (2 * self.bits + (usize::BITS - lo.leading_zeros()) as usize) as f64;
        // GMP splits an unbalanced product into balanced ones.
        let gmp = Self::gmp(lo as f64 * slot) * hi as f64 / lo as f64;
        gmp + out * self.redc + (la + lb) as f64 * self.mul * 0.05
    }

    /// Middle product of [`crate::poly`]: `l` outputs, by a monic polynomial of degree `m`.
    fn middle(&self, l: usize, m: usize) -> f64 {
        if l.min(m) <= SCHOOLBOOK {
            (l * m) as f64 * self.macc + l as f64 * self.redc
        } else {
            self.product(m, l + m - 1)
        }
    }

    /// Product tree of `d` leaves (or the monic polynomial with these roots).
    pub fn tree(&self, d: usize) -> f64 {
        if d == 0 {
            return 0.0;
        }
        let mut cost = 0.0;
        let mut size = 1;
        while size < d {
            for start in (0..d).step_by(2 * size) {
                if start + size < d {
                    let right = size.min(d - start - size);
                    cost += self.product(size, right) + (size + right) as f64 * self.mul * 0.1;
                }
            }
            size *= 2;
        }
        cost
    }

    /// Inverse of a power series modulo `X^d`.
    pub fn inverse(&self, d: usize) -> f64 {
        let (mut cost, mut prec) = (0.0, 1);
        while prec < d {
            let next = (2 * prec).min(d);
            cost += self.product(next, prec) + self.product(next - prec, next - prec);
            prec = next;
        }
        cost
    }

    /// `h*g mod f`, all of degree `d`.
    pub fn mul_mod(&self, d: usize) -> f64 {
        self.product(d, d) + self.product(d - 1, d - 1) + self.product(d - 1, d)
    }

    /// Multipoint evaluation of a polynomial of degree `< d` on its product tree.
    pub fn evaluate(&self, d: usize) -> f64 {
        let mut cost = self.product(d, d);
        let mut size = 1;
        while size < d {
            for start in (0..d).step_by(2 * size) {
                if start + size < d {
                    let right = size.min(d - start - size);
                    cost += self.middle(size, right) + self.middle(right, size);
                }
            }
            size *= 2;
        }
        cost
    }

    /// Point `x` normalized to `z = 1` by a batch inversion of `len` points: 4 multiplications
    /// per point and one modular inversion (about 100 multiplications at these sizes).
    fn normalize(&self, len: usize) -> f64 {
        (4 * len + 100) as f64 * self.mul
    }

    /// Polynomial stage 2 of `plan`.
    pub fn poly_stage2(&self, plan: &PolyPlan) -> f64 {
        let (d1, df, giants) = plan.shape();
        // Baby steps: two chains of point additions (6M) up to d1/2, then normalized.
        let mut cost = (d1 / 6) as f64 * 6.0 * self.mul + self.normalize(df);
        if giants == 0 {
            return cost;
        }
        cost += self.tree(df) + self.inverse(df) + self.evaluate(df) + df as f64 * self.mul;
        // Each block: giant steps (6M) and their normalization, the polynomial G, then
        // H*G mod F (but for the first block). The last block may be shorter.
        let block = |len: usize| match len {
            0 => 0.0,
            _ => len as f64 * 7.0 * self.mul + self.normalize(len) + self.tree(len),
        };
        let (full, rest) = (giants / df, giants % df);
        cost += full as f64 * block(df) + block(rest);
        cost + (giants.div_ceil(df) - 1) as f64 * self.mul_mod(df)
    }

    /// Baby-step giant-step stage 2 with the giant step `d` (see [`crate::stage2`]): about
    /// 0.77 pairs (one multiplication each) per prime of `(b1, b2]`.
    pub fn pairs_stage2(&self, b1: usize, b2: usize, d: usize, baby: usize) -> f64 {
        if b2 <= b1 {
            return 0.0;
        }
        let primes = prime_count(b2) - prime_count(b1);
        let giant = (b2 - b1) / d + 2;
        let points = 6 * (d / 3) + 4 * baby + 11 * giant;
        (points as f64 + 0.77 * 1.1 * primes) * self.mul
    }
}

/// Approximation of the number of primes `<= x`.
fn prime_count(x: usize) -> f64 {
    let x = x.max(3) as f64;
    x / (x.ln() - 1.08)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gmp_interpolation() {
        assert_eq!(Costs::gmp(1024.0), GMP_MUL_NS[0]);
        assert!((Costs::gmp(2f64.powi(20)) - GMP_MUL_NS[10]).abs() < 1.0);
        let mut prev = 0.0;
        for k in 5..30 {
            let c = Costs::gmp(2f64.powi(k));
            assert!(c > prev, "{k}");
            prev = c;
        }
    }

    #[test]
    fn prime_counts() {
        for (x, pi) in [
            (1_000_000, 78_498.0),
            (100_000_000, 5_761_455.0),
            (10_000_000_000usize, 455_052_511.0),
        ] {
            assert!((prime_count(x) / pi - 1.0).abs() < 0.01, "{x}");
        }
    }
}

/// Measures the constants of this module: `cargo test --release -- --ignored primitive_costs
/// --nocapture`.
#[cfg(test)]
mod measure {
    use crate::arith::{with_arith, Arith};
    use rug::{rand::RandState, Assign, Integer};
    use std::time::Instant;

    /// Nanoseconds per modular multiplication.
    fn mul_ns<A: Arith>(a: &A) -> f64 {
        let mut rand = RandState::new();
        let x = a.residue(&a.modulus().clone().random_below(&mut rand));
        let (mut r, mut t) = (x.clone(), x.clone());
        let reps = 1_000_000u32;
        let start = Instant::now();
        for _ in 0..reps {
            a.mul(&mut t, &r, &x);
            std::mem::swap(&mut t, &mut r);
        }
        std::hint::black_box(r);
        start.elapsed().as_nanos() as f64 / reps as f64
    }

    #[test]
    #[ignore]
    fn primitive_costs() {
        let mut rand = RandState::new();
        let mut muls = Vec::new();
        for limbs in 1..=16u32 {
            let mut n = Integer::from(Integer::random_bits(64 * limbs, &mut rand));
            n.set_bit(0, true);
            n.set_bit(64 * limbs - 1, true);
            muls.push(format!("{:.1}", with_arith!(&n, |a| mul_ns(&a))));
        }
        eprintln!("MUL_NS: {}", muls.join(", "));
        let mut products = Vec::new();
        let mut p = Integer::new();
        for k in 10..=25 {
            let x = Integer::from(Integer::random_bits(1 << k, &mut rand));
            let y = Integer::from(Integer::random_bits(1 << k, &mut rand));
            let reps = (1 << 26 >> k).max(3);
            let start = Instant::now();
            for _ in 0..reps {
                p.assign(&x * &y);
            }
            products.push(format!(
                "{:.0}",
                start.elapsed().as_nanos() as f64 / reps as f64
            ));
        }
        eprintln!("GMP_MUL_NS: {}", products.join(", "));
    }
}
