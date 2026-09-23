//! Stage 2 of ECM with polynomial arithmetic (the "FFT continuation" of GMP-ECM).
//!
//! With the giant step `d1` (a multiple of 6 whose prime factors are `<= b1`), the baby steps
//! are `j*Q` for the `dF = phi(d1)/2` numbers `j < d1/2` coprime to `d1`, and the giant steps
//! `m*d1*Q`. A prime `l > d1/2` of `(b1, b2]` is `m*d1 +- j` for some `m` and `j` (`l mod d1` is
//! coprime to `d1`), and `l*Q = O` modulo a prime `p` of `n` makes `x(m*d1*Q) = x(j*Q)`
//! modulo `p`. So stage 2 computes
//!
//! `g = prod_(j, m) (x(j*Q) - x(m*d1*Q))`
//!
//! over all baby steps and `k*dF` giant steps (`k` blocks of `dF`), with polynomials:
//! `F(X) = prod_j (X - x(j*Q))` (with its product tree), then for each block the polynomial
//! `G(X) = prod_m (X - x(m*d1*Q))` over the block, accumulated as `H = prod G mod F`, and at the
//! end `g = prod_j H(x(j*Q))` by multipoint evaluation. `F` and `G` cost `O(M(dF) log dF)`,
//! where `M(dF)` is the cost of a product of polynomials of degree `dF` (Kronecker
//! substitution: nearly linear), so each prime costs much less than the one multiplication of
//! the baby-step giant-step continuation.
//!
//! As there, baby and giant steps are normalized to `z = 1` by batch inversions, and a failed
//! inversion gives a factor. This also checks the primes `l < d1/2`, which are baby steps.

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
    /// Position of `j` among the baby steps for `j < d1/2`, or `u32::MAX` if `gcd(j, d1) > 1`.
    index: Vec<u32>,
    /// Number of baby steps, `phi(d1)/2`: the degree of `F` and of each `G`.
    df: usize,
    /// First giant step `m_lo*d1*Q`.
    m_lo: usize,
    /// Number of giant steps (0 if there is no prime above `d1/2` to check), in blocks of `df`
    /// (the last one possibly shorter).
    giants: usize,
}

impl PolyPlan {
    /// Plan with the giant step `d1`: a multiple of 6 whose prime factors are `<= b1`.
    pub fn new(b1: usize, b2: usize, d1: usize) -> Self {
        let mut plan = Self::shape_only(b1, b2, d1);
        let primes = prime_factors(d1);
        let mut df = 0;
        plan.index = (0..d1 / 2)
            .map(|j| {
                if j == 0 || primes.iter().any(|&p| j.is_multiple_of(p)) {
                    u32::MAX
                } else {
                    df += 1;
                    df as u32 - 1
                }
            })
            .collect();
        debug_assert_eq!(df, plan.df);
        plan
    }

    /// The plan without its baby steps (`index`), cheap to compute: enough for its cost.
    pub(crate) fn shape_only(b1: usize, b2: usize, d1: usize) -> Self {
        assert!(d1.is_multiple_of(6));
        let half = d1 / 2;
        // l = m*d1 +- j with j < d1/2: m = (l + d1/2)/d1. The primes l <= d1/2 are baby steps.
        let m = |l: usize| (l + half) / d1;
        let m_lo = m(b1 + 1).max(1);
        let giants = if b2 > b1.max(half) {
            m(b2) + 1 - m_lo
        } else {
            0
        };
        PolyPlan {
            b1,
            b2,
            d1,
            index: Vec::new(),
            df: phi(d1) / 2,
            m_lo,
            giants,
        }
    }

    /// Largest `b2' >= b2` such that every prime in `(b1, b2']` is checked.
    #[cfg_attr(not(any(test, feature = "bench")), allow(dead_code))]
    pub fn b2_covered(&self) -> usize {
        if self.giants == 0 {
            return self.b2.max(self.d1 / 2);
        }
        // The first number not checked is (m_hi + 1)*d1 - j for the largest baby step j, at most
        // d1/2 - 1.
        let m_hi = self.m_lo + self.giants - 1;
        m_hi * self.d1 + self.d1 / 2
    }

    /// Giant step `d1`, degree `dF` of `F` and number of giant steps.
    pub fn shape(&self) -> (usize, usize, usize) {
        (self.d1, self.df, self.giants)
    }
}

/// Polynomial stage 2 on the residues of `arith`: `gcd(g, n)`, see [`crate::ecm::stage2`].
#[cfg(test)]
pub fn stage2_with<A: PolyArith>(arith: A, q: &crate::point::Point, plan: &PolyPlan) -> Integer {
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
    let (d1, df) = (plan.d1, plan.df);
    let one = a.poly_from(&Integer::from(1));
    if plan.b2 <= plan.b1 {
        return Ok(one);
    }
    let mut normalizer = Normalizer::new(a, df);

    // Also checks the primes l < d1/2: l*Q = O modulo p makes z(l*Q) = 0 mod p.
    let baby = baby_steps(curve, q, d1 / 2, &plan.index, df, &mut normalizer)?;
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
    let leaves: Vec<Elem<G>> = baby.iter().map(|x| leaf(x, &mut t)).collect();
    drop(baby);
    let tree = ProductTree::new(a, &mut ws, leaves);
    let f = tree.root();
    let inv = poly::inverse(a, &mut ws, &poly::reverse_monic(a, f, df), df);

    // Giant steps m*d1*Q from m_lo on.
    let mut scratch = curve.scratch();
    let step = curve.multiple(q, &Integer::from(d1));
    let mut r_prev = curve.multiple(q, &(Integer::from(plan.m_lo) * d1));
    let mut r = curve.multiple(q, &(Integer::from(plan.m_lo + 1) * d1));
    let mut r_next = curve.infinity();
    let mut giant_x = vec![a.zero(); df];
    let mut giant_z = vec![a.zero(); df];
    let mut g = Vec::with_capacity(df);
    let mut tmp = Vec::new();
    let mut h: Vec<Elem<G>> = Vec::new();

    for first in (0..plan.giants).step_by(df) {
        let len = df.min(plan.giants - first);
        for i in 0..len {
            giant_x[i].clone_from(&r_prev.x);
            giant_z[i].clone_from(&r_prev.z);
            curve.add(&mut r_next, &r, &step, &r_prev, &mut scratch);
            std::mem::swap(&mut r_prev, &mut r);
            std::mem::swap(&mut r, &mut r_next);
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
            g.resize(df, a.zero());
        }
        if first == 0 {
            std::mem::swap(&mut h, &mut g);
        } else {
            poly::mul_mod(a, &mut ws, &mut h, &g, f, &inv);
        }
    }

    // g = prod_j H(x_j).
    let values = poly::evaluate(a, &mut ws, &h, &tree, &inv);
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

    /// Numbers `<= b2_covered` checked by `plan`: the baby steps (a point at infinity makes the
    /// normalization fail) and every `m*d1 +- j` for the giant steps `m`.
    fn covered(plan: &PolyPlan) -> Vec<bool> {
        let (d1, half) = (plan.d1, plan.d1 / 2);
        let mut covered = vec![false; plan.b2_covered().max(d1) + 2 * d1];
        let baby: Vec<usize> = (0..half).filter(|&j| plan.index[j] != u32::MAX).collect();
        assert_eq!(baby.len(), plan.df);
        for &j in &baby {
            covered[j] = true;
        }
        for m in plan.m_lo..plan.m_lo + plan.giants {
            assert!(m >= 1);
            for &j in &baby {
                covered[m * d1 - j] = true;
                covered[m * d1 + j] = true;
            }
        }
        covered
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
        for d1 in [6].into_iter().chain(POLY_GIANT_STEPS.into_iter().take(40)) {
            let b1_min = *prime_factors(d1).last().unwrap();
            let mut bounds = vec![(b1_min, 4), (b1_min, b1_min), (b1_min, b1_min + 1)];
            for b1 in [
                b1_min,
                b1_min + 1,
                d1 / 2 - 1,
                d1 / 2,
                d1 - 1,
                d1,
                d1 + 1,
                5 * d1 + 1,
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
                ] {
                    bounds.push((b1, b2));
                }
                bounds.push((b1, b1 + random(20 * d1) + 1));
            }
            for (b1, b2) in bounds {
                let plan = PolyPlan::new(b1, b2, d1);
                assert!(plan.b2_covered() >= b2);
                let covered = covered(&plan);
                for l in Primes::all()
                    .skip_while(|&l| l <= b1)
                    .take_while(|&l| l <= plan.b2_covered())
                {
                    assert!(covered[l], "d1 = {d1}, ({b1}, {b2}]: {l} not covered");
                }
                // Not much more than needed.
                if plan.giants > 0 {
                    assert!(plan.b2_covered() <= b2.max(d1 / 2) + d1, "{d1} {b1} {b2}");
                }
            }
        }
    }

    #[test]
    fn stage2_all_giant_steps() {
        // For small giant steps (several blocks, the last one partial): stage 2 finds p whenever
        // l*Q = O modulo p for a prime l in (b1, b2], whether l < d1/2 (baby step) or not.
        let n = Integer::from(4_009_823u64) * Integer::from(99_476_569u64);
        let mut found = 0;
        for d1 in [6].into_iter().chain(POLY_GIANT_STEPS.into_iter().take(12)) {
            let b1 = (*prime_factors(d1).last().unwrap()).max(30);
            for b2 in [b1 + 3 * d1 + 1000, b1 + 40 * d1] {
                let plan = PolyPlan::new(b1, b2, d1);
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
                    if q.z_cord.clone().gcd(&n) != 1 {
                        continue;
                    }
                    let g = stage2_with(Mont::<1>::new(&n), &q, &plan);
                    assert_eq!(g, stage2_with(Plain::new(&n), &q, &plan));
                    let expected = primes.iter().any(|l| {
                        let g = q.mont_ladder(l).z_cord.gcd(&n);
                        g != 1 && g != n
                    });
                    if expected {
                        assert_ne!(g, 1, "d1 = {d1}, b2 = {b2}, sigma = {sigma}");
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
        for (d1, b2, sizes) in [
            (2310, 1_400_000, &[64, 128, 512, 1024, 1100][..]),
            (30030, 1_000_000, &[64, 320][..]),
        ] {
            for &bits in sizes {
                let mut big = Integer::from(Integer::random_bits(bits, &mut rand));
                big.set_bit(bits - 1, true);
                let n = Integer::from(4_009_823) * big.next_prime();
                let plan = PolyPlan::new(b1, b2, d1);
                let pairs = Stage2Plan::pairs(b1, b2);
                let k = stage1_multiplier(b1);
                let mut found = 0;
                for sigma in 2..8 {
                    let q = stage1(
                        &curve(&n, Param::Batch2, &Integer::from(sigma)).unwrap(),
                        &k,
                    );
                    if q.z_cord.clone().gcd(&n) != 1 {
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
