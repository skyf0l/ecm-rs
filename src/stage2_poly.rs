//! Stage 2 of ECM with polynomial arithmetic (the "FFT continuation" of GMP-ECM).
//!
//! With the giant step `d1` (a multiple of 6 whose prime factors are `<= b1`) and `d2` (1, or
//! a small prime `<= b1` not dividing `d1`), the baby steps are `j*d2*Q` for the `phi(d1)/2`
//! numbers `j < d1/2` coprime to `d1`, and the giant steps `m*d1*Q`. A prime `l` of `(b1, b2]`
//! is `m*d1 +- j*d2` for some `j` (`l*d2^-1 mod d1` or its opposite is below `d1/2` and coprime
//! to `d1`) and `m` (`m < 0` is `|m|` with the other sign), and `l*Q = O` modulo a prime `p` of
//! `n` makes `x(m*d1*Q) = x(j*d2*Q)` modulo `p`. The giant steps with `d2 | m` are skipped:
//! `m*d1 +- j*d2` is then a multiple of `d2`, not a prime (GMP-ECM's `d2`, about `1/d2` fewer
//! giant steps). So stage 2 computes
//!
//! `g = prod_(j, m) (x(j*d2*Q) - x(m*d1*Q))`
//!
//! with polynomials: `F(X) = prod_j (X - x(j*d2*Q))` (with its product tree, of degree `dF`:
//! the baby steps, maybe with the first one repeated so that the giant steps fill whole
//! blocks), then for each block of `dF` giant steps the polynomial `G(X) = prod_m (X -
//! x(m*d1*Q))`, accumulated as `H = prod G mod F`, and at the end `g = prod_j H(x_j)` over the
//! roots of `F` by multipoint evaluation. `F` and `G` cost `O(M(dF) log dF)`, where `M(dF)` is
//! the cost of a product of polynomials of degree `dF` (Kronecker substitution: nearly
//! linear), so each prime costs much less than the one multiplication of the baby-step
//! giant-step continuation.
//!
//! As there, baby and giant steps are normalized to `z = 1` by batch inversions, and a failed
//! inversion gives a factor. This also checks the primes `l < d1/2`, which divide baby steps.

use crate::{
    arith::{Arith, PolyArith},
    curve::Xz,
    poly::{self, ProductTree, Workspace},
    stage2::{baby_steps, phi, prime_factors, Elem, Normalizer, XLine},
};
use rug::Integer;

/// Everything the polynomial stage 2 needs that only depends on the bounds.
#[derive(Debug, Clone)]
pub struct PolyPlan {
    b1: usize,
    b2: usize,
    /// Giant step, a multiple of 6.
    d1: usize,
    /// Multiplier of the baby steps, 1 or a prime `<= b1` not dividing `d1`: the giant steps
    /// `m*d1*Q` with `d2 | m` are skipped.
    d2: usize,
    /// Position of `j` among the baby steps for `j < d1/2`, or `u32::MAX` if `gcd(j, d1) > 1`.
    index: Vec<u32>,
    /// Number of baby steps, `phi(d1)/2`.
    babies: usize,
    /// Degree of `F` (at least `babies`: the first baby step is repeated) and of each `G`.
    df: usize,
    /// Giant steps `m*d1*Q` for `m` in `[m_lo, m_hi]`, `m` not a multiple of `d2`.
    m_lo: usize,
    m_hi: usize,
    /// Number of giant steps (0 if there is no prime above `d1/2` to check), in blocks of `df`
    /// (the last one possibly shorter, only if `df = babies`).
    giants: usize,
}

/// `d2` of the giant step `d1` (see [`PolyPlan`]): the smallest prime `>= 5` not dividing `d1`
/// (as GMP-ECM), if it is `<= b1` (else it may be a prime of stage 2) and at most `babies + 1`.
pub(crate) fn default_d2(b1: usize, d1: usize) -> usize {
    let d2 = [5, 7, 11, 13, 17, 19, 23]
        .into_iter()
        .find(|p| !d1.is_multiple_of(*p))
        .unwrap_or(1);
    if d2 <= b1 && d2 - 1 <= phi(d1) / 2 {
        d2
    } else {
        1
    }
}

/// Number of integers in `[lo, hi]` (`lo >= 1`) that are not multiples of `d2` (all if `d2 = 1`).
fn count_coprime(lo: usize, hi: usize, d2: usize) -> usize {
    if hi < lo {
        return 0;
    }
    if d2 == 1 {
        return hi + 1 - lo;
    }
    (hi + 1 - lo) - (hi / d2 - (lo - 1) / d2)
}

/// The `k`-th (from 1) integer `>= lo` that is not a multiple of `d2`.
fn nth_coprime(lo: usize, k: usize, d2: usize) -> usize {
    if d2 == 1 {
        return lo + k - 1;
    }
    // Rank of the wanted one among all the positive non-multiples.
    let rank = k + (lo - 1) - (lo - 1) / d2;
    (rank - 1) / (d2 - 1) * d2 + (rank - 1) % (d2 - 1) + 1
}

impl PolyPlan {
    /// Plan with the giant step `d1` (a multiple of 6 whose prime factors are `<= b1`), `d2`
    /// (1, or a prime `<= b1` not dividing `d1`), and `blocks` blocks of giant steps, all full
    /// (at most as many as blocks of `phi(d1)/2`), or 0 for blocks of `phi(d1)/2` giant steps,
    /// the last one maybe partial.
    pub fn new(b1: usize, b2: usize, d1: usize, d2: usize, blocks: usize) -> Self {
        let mut plan = Self::shape_only(b1, b2, d1, d2, blocks);
        let primes = prime_factors(d1);
        let mut babies = 0;
        plan.index = (0..d1 / 2)
            .map(|j| {
                if j == 0 || primes.iter().any(|&p| j.is_multiple_of(p)) {
                    u32::MAX
                } else {
                    babies += 1;
                    babies as u32 - 1
                }
            })
            .collect();
        debug_assert_eq!(babies, plan.babies);
        plan
    }

    /// The plan without its baby steps (`index`), cheap to compute: enough for its cost.
    pub(crate) fn shape_only(b1: usize, b2: usize, d1: usize, d2: usize, blocks: usize) -> Self {
        assert!(d1.is_multiple_of(6) && (d2 == 1 || !d1.is_multiple_of(d2)));
        assert!(d2 == 1 || d2 <= b1, "d2 must not be a prime of stage 2");
        let (half, babies) = (d1 / 2, phi(d1) / 2);
        // l = m*d1 +- j*d2 with 1 <= j < d1/2: |l - m*d1| <= (d1/2 - 1)*d2. The primes l < d1/2
        // are baby steps.
        let reach = (half - 1) * d2;
        let m_lo = (b1 + 1).saturating_sub(reach).div_ceil(d1).max(1);
        let (mut m_hi, mut df, mut giants) = (0, babies, 0);
        if b2 > b1.max(half) {
            m_hi = ((b2 + reach) / d1).max(m_lo);
            giants = count_coprime(m_lo, m_hi, d2);
            if blocks > 0 && giants > babies {
                // Whole blocks of at least babies giant steps (not more blocks than with
                // babies giant steps each; a single block stays partial).
                let blocks = blocks.min(giants.div_ceil(babies));
                df = giants.div_ceil(blocks).max(babies);
                giants = blocks * df;
                m_hi = nth_coprime(m_lo, giants, d2);
            }
        }
        PolyPlan {
            b1,
            b2,
            d1,
            d2,
            index: Vec::new(),
            babies,
            df,
            m_lo,
            m_hi,
            giants,
        }
    }

    /// Largest `b2' >= b2` such that every prime in `(b1, b2']` is checked.
    pub fn b2_covered(&self) -> usize {
        if self.giants == 0 {
            return self.b2.max(self.d1 / 2);
        }
        // A prime l is checked if its m is at most m_hi: m <= (l + (d1/2 - 1)*d2)/d1, or
        // m_hi + 1 if it is a multiple of d2 (no prime has this m).
        let mut m_hi = self.m_hi;
        if self.d2 > 1 && (m_hi + 1).is_multiple_of(self.d2) {
            m_hi += 1;
        }
        (m_hi + 1) * self.d1 - (self.d1 / 2 - 1) * self.d2 - 1
    }

    /// Giant step `d1`, degree `dF` of `F` and number of giant steps.
    pub fn shape(&self) -> (usize, usize, usize) {
        (self.d1, self.df, self.giants)
    }

    /// `d2`, and the number of baby steps (at most `dF`).
    pub fn baby_shape(&self) -> (usize, usize) {
        (self.d2, self.babies)
    }
}

/// Polynomial stage 2 on the residues of `arith`: `gcd(g, n)`, see [`crate::ecm::stage2`].
#[cfg(test)]
pub fn stage2_with<A: PolyArith>(arith: A, q: &crate::curve::Point, plan: &PolyPlan) -> Integer {
    crate::stage2::stage2_with(arith, q, &crate::stage2::Stage2Plan::Poly(plan.clone()))
}

/// Product `g` (polynomial representation), or `Err(g)` with a factor found by a failed
/// inversion.
pub(crate) fn accumulate<G: XLine>(
    curve: &G,
    q: &Xz<Elem<G>>,
    plan: &PolyPlan,
) -> Result<Elem<G>, Integer>
where
    G::A: PolyArith,
{
    let a = curve.arith();
    let (d1, d2, df) = (plan.d1, plan.d2, plan.df);
    let one = a.poly_from(&Integer::from(1));
    if plan.b2 <= plan.b1 {
        return Ok(one);
    }
    let mut normalizer = Normalizer::new(a, df);

    // Also checks the primes l < d1/2: l*Q = O modulo p makes z(l*d2*Q) = 0 mod p.
    let base = if d2 == 1 {
        q.clone()
    } else {
        curve.multiple(q, &Integer::from(d2))
    };
    let baby = baby_steps(
        curve,
        &base,
        d1 / 2,
        &plan.index,
        plan.babies,
        &mut normalizer,
    )?;
    if plan.giants == 0 {
        return Ok(one);
    }
    let zero = a.zero();
    let mut t = a.zero();
    // Leaves of a product tree: the constant coefficients -x of X - x.
    let leaf = |x: &Elem<G>, t: &mut Elem<G>| {
        let mut r = a.zero();
        a.to_poly(t, x);
        a.sub(&mut r, &zero, t);
        r
    };
    let mut ws = Workspace::new();
    let mut leaves: Vec<Elem<G>> = baby.iter().map(|x| leaf(x, &mut t)).collect();
    drop(baby);
    // F of degree dF: the first baby step again (H(x) is then just taken several times).
    leaves.resize(df, leaves[0].clone());
    let tree = ProductTree::new(a, &mut ws, leaves);
    let f = tree.root();
    let mut modulus = poly::Modulus::new(a, &mut ws, f);

    // Giant steps m*d1*Q from m_lo on, but for the multiples of d2.
    let mut scratch = curve.scratch();
    let step = curve.multiple(q, &Integer::from(d1));
    let mut r_prev = curve.multiple(q, &(Integer::from(plan.m_lo) * d1));
    let mut r = curve.multiple(q, &(Integer::from(plan.m_lo + 1) * d1));
    let mut r_next = curve.infinity();
    let mut m = plan.m_lo;
    let mut giant_x = vec![a.zero(); df];
    let mut giant_z = vec![a.zero(); df];
    let mut g = Vec::with_capacity(df);
    let mut tmp = Vec::new();
    let mut h: Vec<Elem<G>> = Vec::new();

    for first in (0..plan.giants).step_by(df) {
        let len = df.min(plan.giants - first);
        let mut i = 0;
        while i < len {
            if !m.is_multiple_of(d2) || d2 == 1 {
                giant_x[i].clone_from(&r_prev.x);
                giant_z[i].clone_from(&r_prev.z);
                i += 1;
            }
            curve.add(&mut r_next, &r, &step, &r_prev, &mut scratch);
            std::mem::swap(&mut r_prev, &mut r);
            std::mem::swap(&mut r, &mut r_next);
            m += 1;
        }
        curve.normalize(&mut normalizer, &mut giant_x[..len], &giant_z[..len])?;
        g.clear();
        g.extend(giant_x[..len].iter().map(|x| leaf(x, &mut t)));
        poly::from_roots(a, &mut ws, &mut g, &mut tmp);
        if len == df {
            // G mod F = G - F: both are monic of degree dF.
            for (g, f) in g.iter_mut().zip(f) {
                a.sub(&mut t, g, f);
                std::mem::swap(g, &mut t);
            }
        } else {
            // G mod F = G, with its leading coefficient.
            g.push(one.clone());
        }
        if first == 0 {
            g.resize(df, a.zero());
            std::mem::swap(&mut h, &mut g);
        } else {
            poly::mul_mod(a, &mut ws, &mut h, &g, &mut modulus);
        }
    }
    debug_assert!(m <= plan.m_hi + 1);

    // g = prod_j H(x_j).
    let values = poly::evaluate(a, &mut ws, &h, &tree, &mut modulus);
    let mut acc = one;
    for v in &values {
        a.poly_mul(&mut t, &acc, v);
        std::mem::swap(&mut acc, &mut t);
    }
    Ok(acc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        arith::{with_arith, Mont, Plain},
        ecm::{curve, stage1, stage1_multiplier, stage2, Param},
        stage2::{Stage2Plan, POLY_GIANT_STEPS},
    };
    use primal::Primes;

    /// Numbers `<= b2_covered` checked by `plan`: the baby steps `j` (a point at infinity makes
    /// the normalization fail) and every `|m*d1 +- j*d2|` for the giant steps `m`.
    fn covered(plan: &PolyPlan) -> Vec<bool> {
        let (d1, d2, half) = (plan.d1, plan.d2, plan.d1 / 2);
        let mut covered = vec![false; plan.b2_covered().max(d1) + 1];
        let mut cover = |l: usize| {
            if let Some(c) = covered.get_mut(l) {
                *c = true;
            }
        };
        let baby: Vec<usize> = (0..half).filter(|&j| plan.index[j] != u32::MAX).collect();
        assert_eq!(baby.len(), plan.babies);
        assert!(plan.df >= plan.babies);
        for &j in &baby {
            cover(j);
        }
        let giants: Vec<usize> = (plan.m_lo..=plan.m_hi)
            .filter(|m| d2 == 1 || !m.is_multiple_of(d2))
            .collect();
        assert_eq!(giants.len(), plan.giants);
        for m in giants {
            assert!(m >= 1);
            for &j in &baby {
                cover((m * d1).abs_diff(j * d2));
                cover(m * d1 + j * d2);
            }
        }
        covered
    }

    /// Checks that `plan` covers every prime in `(b1, b2]` and not much more.
    fn check_plan(plan: &PolyPlan, sieve: &primal::Sieve, b1: usize, b2: usize, blocks: usize) {
        let (d1, d2) = (plan.d1, plan.d2);
        let what = format!("d1 = {d1}, d2 = {d2}, ({b1}, {b2}]");
        assert!(plan.b2_covered() >= b2, "{what}");
        let covered = covered(plan);
        for l in sieve
            .primes_from(b1 + 1)
            .take_while(|&l| l <= plan.b2_covered())
        {
            assert!(covered[l], "{what}: {l} not covered");
        }
        // Not much more than needed: the last giant step is needed, but for whole blocks.
        if plan.giants > 0 && blocks == 0 {
            let reach = (d1 / 2 - 1) * d2;
            assert!(plan.m_hi * d1 <= b2 + reach + d1, "{what}");
        }
    }

    #[test]
    fn coprime_counts() {
        for d2 in [1, 5, 7, 13] {
            for lo in 1..40 {
                for k in 1..60 {
                    let m = nth_coprime(lo, k, d2);
                    assert!(d2 == 1 || !m.is_multiple_of(d2));
                    assert_eq!(count_coprime(lo, m, d2), k, "{lo} {k} {d2}");
                }
                assert_eq!(count_coprime(lo, lo - 1, d2), 0);
            }
        }
    }

    #[test]
    fn plans_cover_primes() {
        let mut rand = rug::rand::RandState::new();
        let mut random = |max: usize| {
            Integer::from(max)
                .random_below(&mut rand)
                .to_usize()
                .unwrap()
        };
        let sieve = primal::Sieve::new(4_000_000);
        for d1 in [6].into_iter().chain(POLY_GIANT_STEPS.into_iter().take(40)) {
            let b1_min = *prime_factors(d1).last().unwrap();
            let mut bounds = vec![(b1_min, 4), (b1_min, b1_min), (b1_min, b1_min + 1)];
            for b1 in [
                b1_min,
                b1_min + 1,
                23,
                d1 / 2 - 1,
                d1 / 2,
                d1 - 1,
                d1,
                d1 + 1,
                5 * d1 + 1,
                13 * d1 / 2,
            ] {
                let b1 = b1.max(b1_min);
                for b2 in [
                    b1 + 1,
                    b1 + 2,
                    d1 / 2,
                    d1 / 2 + 1,
                    d1 - 1,
                    d1 + 1,
                    2 * d1 + 1,
                    7 * d1,
                    40 * d1,
                ] {
                    bounds.push((b1, b2));
                }
                bounds.push((b1, b1 + random(20 * d1) + 1));
            }
            for (b1, b2) in bounds {
                // Without and with d2 (the default one, or any prime <= b1 not dividing d1),
                // blocks of phi(d1)/2 giant steps, or fewer whole blocks.
                let d2 = default_d2(b1, d1);
                let others: Vec<usize> = [5, 7, 11, 13, 17, 19, 23, 29]
                    .into_iter()
                    .filter(|&p| p <= b1 && !d1.is_multiple_of(p))
                    .collect();
                let other = match others.len() {
                    0 => 1,
                    len => others[random(len)],
                };
                for (d2, blocks) in [(1, 0), (d2, 0), (d2, 1 + random(4)), (other, random(4))] {
                    let plan = PolyPlan::new(b1, b2, d1, d2, blocks);
                    check_plan(&plan, &sieve, b1, b2, blocks);
                }
            }
        }
    }

    #[test]
    fn stage2_all_giant_steps() {
        // For small giant steps (several blocks, the last one partial, or whole blocks with F
        // padded; with and without d2): stage 2 finds p whenever l*Q = O modulo p for a prime l
        // in (b1, b2], whether l < d1/2 (baby step) or not.
        let n = Integer::from(4_009_823u64) * Integer::from(99_476_569u64);
        let mut found = 0;
        for d1 in [6].into_iter().chain(POLY_GIANT_STEPS.into_iter().take(12)) {
            let b1 = (*prime_factors(d1).last().unwrap()).max(30);
            let shapes = [
                (1, 0),
                (default_d2(b1, d1), 0),
                (default_d2(b1, d1), 2),
                (29, 3),
            ];
            for ((d2, blocks), b2) in shapes.into_iter().zip([
                b1 + 3 * d1 + 1000,
                b1 + 40 * d1,
                b1 + 17 * d1,
                b1 + 9 * d1 + 5000,
            ]) {
                let plan = PolyPlan::new(b1, b2, d1, d2, blocks);
                let k = stage1_multiplier(b1);
                let primes: Vec<Integer> = Primes::all()
                    .skip_while(|&l| l <= b1)
                    .take_while(|&l| l <= b2)
                    .map(Integer::from)
                    .collect();
                for sigma in 2..20 {
                    let q = stage1(
                        &curve(&n, Param::Batch2, &Integer::from(sigma)).unwrap(),
                        &k,
                    );
                    if q.z.clone().gcd(&n) != 1 {
                        continue;
                    }
                    let g = stage2_with(Mont::<1>::new(&n), &q, &plan);
                    assert_eq!(g, stage2_with(Plain::new(&n), &q, &plan));
                    let expected = primes.iter().any(|l| {
                        let g = stage1(&q, l).z.gcd(&n);
                        g != 1 && g != n
                    });
                    if expected {
                        assert_ne!(g, 1, "d1 = {d1}, {d2} {blocks}, b2 = {b2}, sigma = {sigma}");
                        found += 1;
                    }
                }
            }
        }
        assert!(found > 20, "{found}");
    }

    #[test]
    fn stage2_large_d1_limbs() {
        // Large baby-step sets (Kronecker products at every level of the tree, including the
        // middle products of the evaluation) for every kind of arithmetic: whenever the
        // baby-step giant-step continuation finds the small factor, so must the polynomial one.
        let mut rand = rug::rand::RandState::new();
        let b1 = 100;
        for (d1, blocks, b2, sizes) in [
            (2310, 0, 1_400_000, &[64, 128, 512, 1024, 1100][..]),
            (30030, 2, 1_000_000, &[64, 320][..]),
        ] {
            for &bits in sizes {
                let mut big = Integer::from(Integer::random_bits(bits, &mut rand));
                big.set_bit(bits - 1, true);
                let n = Integer::from(4_009_823) * big.next_prime();
                let plan = PolyPlan::new(b1, b2, d1, default_d2(b1, d1), blocks);
                let pairs = Stage2Plan::pairs(b1, b2);
                let k = stage1_multiplier(b1);
                let mut found = 0;
                for sigma in 2..8 {
                    let q = stage1(
                        &curve(&n, Param::Batch2, &Integer::from(sigma)).unwrap(),
                        &k,
                    );
                    if q.z.clone().gcd(&n) != 1 {
                        continue;
                    }
                    let g = with_arith!(&n, |a| stage2_with(a, &q, &plan));
                    assert_eq!(g, stage2_with(Plain::new(&n), &q, &plan));
                    if stage2(&q, &pairs) != 1 {
                        assert_ne!(g, 1, "d1 = {d1}, {bits} bits, sigma = {sigma}");
                        found += 1;
                    }
                }
                assert!(found > 0, "d1 = {d1}, {bits} bits");
            }
        }
    }
}
