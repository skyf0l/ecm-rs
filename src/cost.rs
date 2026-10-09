//! Cost model of stage 2, to choose between the continuations and their parameters.
//!
//! Costs are in nanoseconds, built from the measured costs of the primitive operations (on an
//! Intel i7-8750H): modular multiplications by size of `n`, and GMP products by size. Only
//! their ratios matter, which vary much less from one machine to another. The polynomial
//! costs mirror the operations of [`crate::poly`] one product at a time.

use crate::{
    arith::mpn,
    base2::{Base2Form, Product},
    poly,
    stage2_poly::PolyPlan,
};
use std::collections::HashMap;

/// Montgomery multiplication, by number of limbs (index 0 unused).
const MUL_NS: [f64; 17] = [
    0.0, 3.7, 9.7, 16.8, 24.1, 43.6, 56.6, 74.6, 94.8, 120.5, 143.5, 200.0, 241.0, 274.0, 311.0,
    355.0, 402.0,
];

/// Step of the stage 1 ladder (a doubling and an addition with a full base point: 6M + 4S and a
/// multiplication by a small `a24`) with Montgomery's arithmetic, by number of limbs (index 0
/// unused).
const STEP_NS: [f64; 17] = [
    0.0, 21.0, 88.0, 148.0, 235.0, 354.0, 524.0, 701.0, 872.0, 1077.0, 1517.0, 1828.0, 2180.0,
    2608.0, 2784.0, 3186.0, 3587.0,
];

/// [`STEP_NS`] with the special reduction modulo `2^k +- 1`, by limbs of `2^k` (index 0 unused).
const BASE2_STEP_NS: [f64; 25] = [
    0.0, 254.0, 233.0, 300.0, 361.0, 428.0, 495.0, 584.0, 666.0, 777.0, 939.0, 1050.0, 1143.0,
    1281.0, 1399.0, 1623.0, 1743.0, 2005.0, 2125.0, 2319.0, 2448.0, 2670.0, 2843.0, 3005.0, 3216.0,
];

/// Limb counts at which the costs above [`crate::arith::MAX_LIMBS`] limbs (the Montgomery
/// arithmetic of runtime length, GMP's products and reductions) were measured; in between and
/// beyond, they are interpolated on a log-log scale ([`large`]). Measured with the fastest of
/// repeated runs, scaled to the costs at 16 limbs of [`MUL_NS`], [`STEP_NS`] and
/// [`BASE2_STEP_NS`] measured in the same run (`large_costs` below).
pub(crate) const LARGE_LIMBS: [usize; 10] = [17, 20, 24, 32, 48, 64, 96, 128, 192, 256];

/// Montgomery multiplication, at [`LARGE_LIMBS`] limbs.
const LARGE_MUL_NS: [f64; 10] = [
    433.0, 562.0, 782.0, 1296.0, 2705.0, 4255.0, 8088.0, 12983.0, 23549.0, 38456.0,
];

/// Step of the stage 1 ladder as [`STEP_NS`], at [`LARGE_LIMBS`] limbs.
const LARGE_STEP_NS: [f64; 10] = [
    3852.0, 5142.0, 7091.0, 11699.0, 24377.0, 38778.0, 74236.0, 122712.0, 218019.0, 353342.0,
];

/// [`LARGE_STEP_NS`] with the special reduction modulo `2^k +- 1`, for `k = 64*limbs - 7`: a
/// full product and a fold (the wrap-around products of the `k` multiple of 64 are
/// [`product_ratio`] times cheaper).
const LARGE_BASE2_STEP_NS: [f64; 10] = [
    1964.0, 2446.0, 3211.0, 4975.0, 9321.0, 14696.0, 27588.0, 43284.0, 79040.0, 122516.0,
];

/// Cost at `limbs` limbs from the costs `ns` measured at [`LARGE_LIMBS`]: interpolated on a
/// log-log scale, extrapolated with the slope of the nearest two points.
fn large(ns: &[f64; 10], limbs: usize) -> f64 {
    let x = (limbs as f64).ln();
    let points = LARGE_LIMBS.map(|l| (l as f64).ln());
    let i = points
        .iter()
        .rposition(|&p| p <= x)
        .unwrap_or(0)
        .min(points.len() - 2);
    let slope = (ns[i + 1].ln() - ns[i].ln()) / (points[i + 1] - points[i]);
    (ns[i].ln() + slope * (x - points[i])).exp()
}

/// The ladder steps measured with [`BASE2_STEP_NS`] run about this much slower in the curves
/// than the ones of [`STEP_NS`]: the special reduction is used when the measured steps are
/// faster by at least this factor.
const BASE2_MARGIN: f64 = 0.9;

/// Cost of a product modulo `2^k +- 1` computed by [`Product::Wrap`] or [`Product::Fft`],
/// relative to a full product and a fold ([`Product::Fold`], `1`), by limbs of `2^k`: the
/// averages of a multiplication and a squaring, in instructions (the measures of
/// [`crate::config::BASE2_WRAP_LIMBS`] and [`crate::config::BASE2_FFT_LIMBS`]), interpolated
/// on a log scale. The gain of the wrap-around products also depends on the powers of 2 that
/// divide the number of limbs (GMP's recursion): the powers of 2 were measured.
fn product_ratio(product: Product, limbs: usize) -> f64 {
    const WRAP: [(f64, f64); 9] = [
        (8.0, 0.90),
        (16.0, 0.82),
        (32.0, 0.66),
        (64.0, 0.57),
        (128.0, 0.55),
        (256.0, 0.56),
        (512.0, 0.62),
        (1024.0, 0.58),
        (2048.0, 0.46),
    ];
    const FFT: [(f64, f64); 4] = [(448.0, 0.86), (512.0, 0.84), (1024.0, 0.68), (2048.0, 0.56)];
    let points: &[(f64, f64)] = match product {
        Product::Fold => return 1.0,
        Product::Wrap => &WRAP,
        Product::Fft { .. } => &FFT,
    };
    let x = limbs as f64;
    let last = points.len() - 1;
    if x <= points[0].0 {
        return points[0].1;
    }
    if x >= points[last].0 {
        return points[last].1;
    }
    let i = points
        .iter()
        .rposition(|p| p.0 <= x)
        .unwrap_or(0)
        .min(last - 1);
    let ((x0, y0), (x1, y1)) = (points[i], points[i + 1]);
    y0 + (y1 - y0) * (x / x0).ln() / (x1 / x0).ln()
}

/// Share of the products in the cost of a stage 1 step modulo `2^k +- 1` above 16 limbs (the
/// rest: additions, subtractions, multiplications by small constants).
const BASE2_PRODUCT_SHARE: f64 = 0.85;

/// Whether the arithmetic modulo `2^k +- 1` (the form) is faster than the arithmetic modulo
/// a number of `bits` bits: from the measured costs of the stage 1 steps with both (with the
/// wrap-around products of [`Product`] where they apply), and always above
/// [`crate::arith::MAX_LIMBS`] limbs without GMP's low-level functions (plain integers, whose
/// division costs more than a product, while `k` is at most 1.4 times `bits`).
pub(crate) fn base2_faster(bits: u32, form: Base2Form) -> bool {
    let limbs = bits.div_ceil(64) as usize;
    let base2 = form.k.div_ceil(64).max(1) as usize;
    let mont = if limbs < STEP_NS.len() {
        STEP_NS[limbs]
    } else if mpn::ENABLED {
        large(&LARGE_STEP_NS, limbs)
    } else {
        return true;
    };
    let cost = if base2 < BASE2_STEP_NS.len() {
        BASE2_STEP_NS[base2]
    } else {
        large(&LARGE_BASE2_STEP_NS, base2)
    };
    let ratio = product_ratio(Product::of(form), base2);
    let cost = cost * (1.0 - BASE2_PRODUCT_SHARE * (1.0 - ratio));
    cost < BASE2_MARGIN * mont
}

/// Packing and unpacking a coefficient of a Kronecker product, besides its share of a modular
/// multiplication.
const PACK_NS: f64 = 30.0;

/// Overhead of a pair of the baby-step giant-step continuation, besides its modular
/// multiplication.
const PAIR_NS: f64 = 2.0;

/// `b2` at which the cost of a pair was measured.
const PAIRS_B2: f64 = 12.7e6;

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

/// GMP product modulo `2^N - 1` (`mpn_mulmod_bnm1`) for `N = 2^(10 + i)` bits.
const WRAP_NS: [f64; 16] = [
    150.0, 386.0, 1514.0, 3353.0, 8727.0, 24943.0, 59166.0, 135856.0, 325854.0, 799854.0,
    1601969.0, 3253152.0, 7916731.0, 15941225.0, 38962050.0, 84787407.0,
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
    /// Cost of a modular multiplication.
    pub fn mul(&self) -> f64 {
        self.mul
    }

    /// Size of a residue in memory.
    pub fn elem_bytes(&self) -> f64 {
        (self.bits.div_ceil(64).max(1) * 8) as f64
    }

    /// Costs modulo a number of `bits` bits, with the special reduction modulo `base2` if any.
    pub fn modulo(bits: usize, base2: Option<Base2Form>) -> Self {
        let Some(form) = base2 else {
            return Self::new(bits);
        };
        // A product (GMP's), then a reduction by additions and shifts, or a wrap-around product
        // (the accumulated products of the polynomials are full ones).
        let limbs = form.k.div_ceil(64) as usize;
        let product = Self::gmp((64 * limbs) as f64);
        let wrapped = product * product_ratio(Product::of(form), limbs);
        Self {
            bits: form.k as usize + usize::from(form.plus),
            mul: 1.2 * wrapped + 25.0,
            macc: product + 10.0,
            redc: 20.0 + 2.0 * limbs as f64,
        }
    }

    pub fn new(bits: usize) -> Self {
        let limbs = bits.div_ceil(64).max(1);
        if limbs < MUL_NS.len() {
            let mul = MUL_NS[limbs];
            return Self {
                bits,
                mul,
                macc: 0.6 * mul,
                redc: 0.8 * mul,
            };
        }
        if !mpn::ENABLED {
            // Plain arithmetic: roughly quadratic, and slower than Montgomery's.
            let mul = 2.0 * MUL_NS[16] * (limbs as f64 / 16.0).powf(1.8);
            return Self {
                bits,
                mul,
                macc: 0.6 * mul,
                redc: 0.8 * mul,
            };
        }
        // GMP's product, then its Montgomery reduction: the product alone is about 45% of a
        // multiplication (50% at 17 limbs, 40% from 64 limbs), the reduction 60%.
        let mul = large(&LARGE_MUL_NS, limbs);
        Self {
            bits,
            mul,
            macc: 0.45 * mul,
            redc: 0.6 * mul,
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

    /// GMP product modulo `2^N - 1` (`mpn_mulmod_bnm1`) of `N = 64*rn` bits.
    fn gmp_wrap(bits: f64) -> f64 {
        let x = bits.max(1.0).log2() - 10.0;
        let last = WRAP_NS.len() - 1;
        if x <= 0.0 {
            return WRAP_NS[0] * bits / 1024.0;
        }
        if x >= last as f64 {
            let big = WRAP_NS[last];
            let scale = 2f64.powf(x - last as f64);
            return big * scale * (x + 10.0) / (last as f64 + 10.0);
        }
        let i = x.floor() as usize;
        let f = x - i as f64;
        (WRAP_NS[i].ln() * (1.0 - f) + WRAP_NS[i + 1].ln() * f).exp()
    }

    /// Packing `count` coefficients for a Kronecker product.
    fn pack(&self, count: usize) -> f64 {
        count as f64 * (self.mul * 0.05 + PACK_NS)
    }

    /// Product of polynomials with `la` and `lb` coefficients.
    pub fn product(&self, la: usize, lb: usize) -> f64 {
        self.part(la, lb, 0, la + lb - 1)
    }

    /// Coefficients `from..to` of the product of polynomials with `la` and `lb` coefficients
    /// ([`crate::poly::mul_part`]).
    fn part(&self, la: usize, lb: usize, from: usize, to: usize) -> f64 {
        let (lo, hi) = (la.min(lb), la.max(lb));
        if lo == 0 {
            return 0.0;
        }
        let out = (to - from) as f64;
        if lo <= SCHOOLBOOK {
            return (lo * (to - from).min(hi)) as f64 * self.macc + out * self.redc;
        }
        let unpack = out * self.redc + self.pack(la + lb);
        if let Some((_, rn)) = poly::middle_size(self.bits, la, lb, from, to) {
            return Self::gmp_wrap((64 * rn) as f64) + unpack;
        }
        let slot = poly::slot_width(self.bits, lo) as f64;
        // GMP splits an unbalanced product into balanced ones.
        Self::gmp(lo as f64 * slot) * hi as f64 / lo as f64 + unpack
    }

    /// `count` coefficients of the product of polynomials with `la` and `lb` coefficients
    /// modulo `X^L - 1`, `L >= lmin` ([`crate::poly::mul_wrap`]).
    fn wrap(&self, la: usize, lb: usize, lmin: usize, count: usize) -> f64 {
        if la.min(lb) > SCHOOLBOOK {
            if let Some(shape) = poly::wrap_size(self.bits, la, lb, lmin) {
                let unpack = count as f64 * self.redc + self.pack(la + lb);
                return Self::gmp_wrap((64 * shape.rn) as f64) + unpack;
            }
        }
        self.part(la, lb, 0, count)
    }

    /// Middle product of [`crate::poly`]: `l` outputs, by a monic polynomial of degree `m`.
    fn middle(&self, l: usize, m: usize) -> f64 {
        if l.min(m) <= SCHOOLBOOK {
            (l * m) as f64 * self.macc + l as f64 * self.redc
        } else {
            self.part(m, l + m - 1, m - 1, l + m - 1)
        }
    }

    /// Sum of `cost(left, right)` over the pairs of sibling nodes (of `left` and `right`
    /// leaves) of the product tree of `d` leaves.
    fn tree_pairs(d: usize, mut cost: impl FnMut(usize, usize) -> f64) -> f64 {
        let (mut total, mut size) = (0.0, 1);
        while size < d {
            // Nodes of `size` leaves, and maybe a last shorter one.
            let (full, rest) = (d / size, d % size);
            let nodes = full + usize::from(rest > 0);
            if rest > 0 && nodes.is_multiple_of(2) {
                total += (nodes / 2 - 1) as f64 * cost(size, size) + cost(size, rest);
            } else {
                total += (nodes / 2) as f64 * cost(size, size);
            }
            size *= 2;
        }
        total
    }

    /// Product tree of `d` leaves (or the monic polynomial with these roots).
    pub fn tree(&self, d: usize) -> f64 {
        Self::tree_pairs(d, |left, right| {
            self.product(left, right) + (left + right) as f64 * self.mul * 0.1
        })
    }

    /// Inverse of a power series modulo `X^d`.
    pub fn inverse(&self, d: usize) -> f64 {
        let (mut cost, mut prec) = (0.0, 1);
        while prec < d {
            let next = (2 * prec).min(d);
            let h = next - prec;
            cost += self.part(next, prec, prec, next) + self.part(h, h, 0, h);
            prec = next;
        }
        cost
    }

    /// `h*g mod f`, all of degree `d`.
    pub fn mul_mod(&self, d: usize) -> f64 {
        self.product(d, d) + self.part(d - 1, d - 1, 0, d - 1) + self.wrap(d + 1, d - 1, d + 1, d)
    }

    /// Multipoint evaluation of a polynomial of degree `< d` on its product tree.
    pub fn evaluate(&self, d: usize) -> f64 {
        self.part(d, d, 0, d)
            + Self::tree_pairs(d, |left, right| {
                self.middle(left, right) + self.middle(right, left)
            })
    }

    /// Point `x` normalized to `z = 1` by a batch inversion of `len` points: 4 multiplications
    /// per point and one modular inversion (about 100 multiplications at these sizes).
    fn normalize(&self, len: usize) -> f64 {
        (4 * len + 100) as f64 * self.mul
    }

    /// Polynomial stage 2 of `plan`. `cache` keeps the costs that only depend on `dF`, for
    /// the next plans with the same `dF`.
    pub fn poly_stage2(&self, plan: &PolyPlan, cache: &mut HashMap<usize, [f64; 3]>) -> f64 {
        let (d1, df, giants) = plan.shape();
        let (d2, babies) = plan.baby_shape();
        // Baby steps: two chains of point additions (6M) up to d1/2, then normalized.
        let mut cost = (d1 / 6) as f64 * 6.0 * self.mul + self.normalize(babies);
        if giants == 0 {
            return cost;
        }
        let [tree, fixed, mul_mod] = *cache.entry(df).or_insert_with(|| {
            let tree = self.tree(df);
            let fixed = tree + self.inverse(df) + self.evaluate(df);
            [tree, fixed, self.mul_mod(df)]
        });
        cost += fixed + df as f64 * self.mul;
        // Each block: giant steps (6M, and the skipped multiples of d2) and their
        // normalization, the polynomial G, then H*G mod F (but for the first block). The last
        // block may be shorter.
        let step = 6.0 * d2 as f64 / (d2 - 1).max(1) as f64 + 1.0;
        let block = |len: usize| match len {
            0 => 0.0,
            _ => {
                let g = if len == df { tree } else { self.tree(len) };
                len as f64 * step * self.mul + self.normalize(len) + g
            }
        };
        let (full, rest) = (giants / df, giants % df);
        cost += full as f64 * block(df) + block(rest);
        cost + (giants.div_ceil(df) - 1) as f64 * mul_mod
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
        // A pair also costs a subtraction and finding it in the table (a few nanoseconds), and
        // a bit more as the table outgrows the caches (measured: +1.7% when b2 doubles).
        let memory = (1.0 + 0.017 * (b2 as f64 / PAIRS_B2).log2()).max(0.9);
        let pair = (1.07 * self.mul + PAIR_NS) * memory;
        points as f64 * self.mul + 0.77 * primes * pair
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
    fn large_interpolation() {
        // Exact at the measured sizes, increasing, and continuing the costs up to 16 limbs.
        for (i, &limbs) in LARGE_LIMBS.iter().enumerate() {
            assert!((large(&LARGE_MUL_NS, limbs) / LARGE_MUL_NS[i] - 1.0).abs() < 1e-9);
        }
        let mut prev = MUL_NS[16];
        for limbs in 17..2000 {
            let mul = Costs::new(64 * limbs).mul();
            assert!(mul > prev, "{limbs}");
            prev = mul;
        }
        assert!(large(&LARGE_STEP_NS, 17) > STEP_NS[16]);
        assert!((large(&LARGE_BASE2_STEP_NS, 24) / BASE2_STEP_NS[24] - 1.0).abs() < 0.05);
    }

    #[test]
    fn base2_above_16_limbs() {
        // The special reduction stays faster above 16 limbs, even with k at 1.4 times the
        // size of n (a product 1.4 times larger, but no Montgomery reduction).
        for bits in [1025, 1100, 1536, 2048, 3072, 4096, 8192, 16384, 40000] {
            for k in [bits + 1, bits * 5 / 4, bits * 7 / 5, (bits + 64) & !63] {
                for plus in [false, true] {
                    let form = Base2Form { k, plus };
                    assert!(base2_faster(bits, form), "{bits} bits, {form}");
                }
            }
        }
    }

    #[test]
    fn wrap_around_costs() {
        // The wrap-around products cost less than a full product and a fold, and a
        // multiplication modulo 2^k +- 1 with them less than without (k one bit off).
        for limbs in [1, 8, 12, 16, 17, 100, 448, 500, 1000, 5000] {
            assert_eq!(product_ratio(Product::Fold, limbs), 1.0);
            for product in [Product::Wrap, Product::Fft { mul: 4, sqr: 4 }] {
                let ratio = product_ratio(product, limbs);
                assert!(ratio > 0.4 && ratio < 1.0, "{product:?} {limbs}");
            }
        }
        if !mpn::ENABLED {
            return;
        }
        for (k, plus) in [(1024, false), (4096, false), (28672, true), (65536, true)] {
            let bits = k as usize - 10;
            let wrap = Costs::modulo(bits, Some(Base2Form { k, plus })).mul();
            let fold = Costs::modulo(bits, Some(Base2Form { k: k - 1, plus })).mul();
            assert!(wrap < 0.95 * fold, "{k} {plus}");
            assert!(base2_faster(k - 10, Base2Form { k, plus }));
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
/// --nocapture` (and `ladder_costs` for the ladder steps, `large_costs` for the costs above 16
/// limbs).
#[cfg(test)]
mod measure {
    use crate::arith::{Arith, mpn, with_arith};
    use rug::{Assign, Integer, rand::RandState};
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
        start.elapsed().as_nanos() as f64 / f64::from(reps)
    }

    /// Nanoseconds per step of the stage 1 ladder (a doubling and an addition).
    fn step_ns<A: Arith>(a: A) -> f64 {
        use crate::{arith::Factor, curve::Curve, stop::Stop};
        let curve = Curve::new(a, &Integer::from(123_456_789));
        let x = curve
            .arith
            .factor(&(Integer::from(Integer::u_pow_u(7, 1000)) % curve.arith.modulus()));
        let k = Integer::from(Integer::u_pow_u(3, 4000));
        let reps = 3;
        let start = Instant::now();
        for _ in 0..reps {
            std::hint::black_box(curve.ladder(&x, &Factor::One, &k, Stop::NEVER));
        }
        start.elapsed().as_nanos() as f64 / f64::from(reps * k.significant_bits())
    }

    /// Nanoseconds per modular multiplication, the fastest of 5 runs (the others were slowed
    /// by other processes), with fewer repetitions for large sizes.
    fn min_mul_ns<A: Arith>(a: &A) -> f64 {
        let limbs = a.modulus().significant_bits().div_ceil(64);
        let reps = (4_000_000 / (limbs * limbs)).max(200);
        let mut rand = RandState::new();
        let x = a.residue(&a.modulus().clone().random_below(&mut rand));
        let (mut r, mut t) = (x.clone(), x.clone());
        let mut best = f64::INFINITY;
        for _ in 0..5 {
            let start = Instant::now();
            for _ in 0..reps {
                a.mul(&mut t, &r, &x);
                std::mem::swap(&mut t, &mut r);
            }
            best = best.min(start.elapsed().as_nanos() as f64 / f64::from(reps));
        }
        std::hint::black_box(r);
        best
    }

    /// Nanoseconds per step of the stage 1 ladder, the fastest of 5 ladders of 1585 bits.
    fn min_step_ns<A: Arith>(a: A) -> f64 {
        use crate::{arith::Factor, curve::Curve, stop::Stop};
        let curve = Curve::new(a, &Integer::from(123_456_789));
        let x = curve
            .arith
            .factor(&(Integer::from(Integer::u_pow_u(7, 1000)) % curve.arith.modulus()));
        let k = Integer::from(Integer::u_pow_u(3, 1000));
        let mut best = f64::INFINITY;
        for _ in 0..5 {
            let start = Instant::now();
            std::hint::black_box(curve.ladder(&x, &Factor::One, &k, Stop::NEVER));
            best = best.min(start.elapsed().as_nanos() as f64 / f64::from(k.significant_bits()));
        }
        best
    }

    #[test]
    #[ignore = "measurement (minutes), prints the constants of this module"]
    fn large_costs() {
        use crate::base2::{Base2, Base2Form};
        let mut rand = RandState::new();
        // Costs at `limbs` limbs: a multiplication and a step, and a step modulo 2^k +- 1 (both
        // forms, k not a multiple of 64: the shifted reduction).
        let mut measure = |limbs: u32| {
            let mut n = Integer::from(Integer::random_bits(64 * limbs, &mut rand));
            n.set_bit(0, true);
            n.set_bit(64 * limbs - 1, true);
            let mul = with_arith!(&n, |a| min_mul_ns(&a));
            let step = with_arith!(&n, |a| min_step_ns(a));
            let base2: f64 = [true, false]
                .map(|plus| {
                    let form = Base2Form {
                        k: 64 * limbs - 7,
                        plus,
                    };
                    min_step_ns(Base2::new(&form.value(), form))
                })
                .iter()
                .sum();
            [mul, step, base2 / 2.0]
        };
        // Scaled to the costs at 16 limbs of the tables above, measured in the same conditions
        // (the clock speed of a loaded machine varies).
        let reference = measure(16);
        let scale = [
            super::MUL_NS[16] / reference[0],
            super::STEP_NS[16] / reference[1],
            super::BASE2_STEP_NS[16] / reference[2],
        ];
        let (mut muls, mut steps, mut base2) = (Vec::new(), Vec::new(), Vec::new());
        for limbs in super::LARGE_LIMBS {
            let [mul, step, b2] = measure(limbs as u32);
            muls.push(format!("{:.0}", mul * scale[0]));
            steps.push(format!("{:.0}", step * scale[1]));
            base2.push(format!("{:.0}", b2 * scale[2]));
        }
        eprintln!("LARGE_MUL_NS: {}", muls.join(", "));
        eprintln!("LARGE_STEP_NS: {}", steps.join(", "));
        eprintln!("LARGE_BASE2_STEP_NS: {}", base2.join(", "));
    }

    #[test]
    #[ignore = "measurement (a minute), prints the constants of this module"]
    fn ladder_costs() {
        use crate::base2::{Base2, Base2Form};
        let mut rand = RandState::new();
        let mut mont = Vec::new();
        for limbs in 1..=16u32 {
            let mut n = Integer::from(Integer::random_bits(64 * limbs, &mut rand));
            n.set_bit(0, true);
            n.set_bit(64 * limbs - 1, true);
            mont.push(format!("{:.0}", with_arith!(&n, |a| step_ns(a))));
        }
        eprintln!("STEP_NS: {}", mont.join(", "));
        let mut base2 = Vec::new();
        for limbs in 1..=24u32 {
            // Both forms, k not a multiple of 64 (the shifted reduction).
            let ns: f64 = [true, false]
                .map(|plus| {
                    let form = Base2Form {
                        k: 64 * limbs - 7,
                        plus,
                    };
                    step_ns(Base2::new(&form.value(), form))
                })
                .iter()
                .sum();
            base2.push(format!("{:.0}", ns / 2.0));
        }
        eprintln!("BASE2_STEP_NS: {}", base2.join(", "));
    }

    #[test]
    #[ignore = "measurement (minutes), prints the constants of this module"]
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
                start.elapsed().as_nanos() as f64 / f64::from(reps)
            ));
        }
        eprintln!("GMP_MUL_NS: {}", products.join(", "));
        let mut wraps = Vec::new();
        let mut scratch = Vec::new();
        for k in 10..=25 {
            let rn = mpn::mulmod_bnm1_next_size(1 << (k - 6));
            let x: Vec<u64> = (0..rn as u64)
                .map(|i| i.wrapping_mul(0x9e37_79b9_7f4a_7c15))
                .collect();
            let y: Vec<u64> = (0..rn as u64)
                .map(|i| i.wrapping_mul(0xc2b2_ae3d_27d4_eb4f))
                .collect();
            let mut r = vec![0; rn];
            let reps = (1 << 26 >> k).max(3);
            let start = Instant::now();
            for _ in 0..reps {
                mpn::mulmod_bnm1(&mut r, &x, &y, &mut scratch);
            }
            let ns = start.elapsed().as_nanos() as f64 / f64::from(reps);
            // Per 2^k bits.
            wraps.push(format!("{:.0}", ns * f64::from(1 << (k - 6)) / rn as f64));
        }
        eprintln!("WRAP_NS: {}", wraps.join(", "));
    }
}
