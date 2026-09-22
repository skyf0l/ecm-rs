//! Stage 2 of ECM: the baby-step giant-step standard continuation, with a wheel and prime
//! pairing.
//!
//! With the wheel modulus `D` (a multiple of small primes, all `<= b1`), the baby steps are
//! `S_j = j*Q` for `1 <= j < D` coprime to `D`, and the giant steps `R_m = m*D*Q`. Since
//! `x(P) = x(-P)`, `x(R_m) = x(S_j)` modulo a prime `p` of `n` exactly when `(m*D - j)*Q = O`
//! or `(m*D + j)*Q = O` modulo `p`: one product `x(R_m) - x(S_j)` checks the two numbers
//! `m*D - j` and `m*D + j`. Every number coprime to `D` is covered by two such pairs (with `j`
//! and with `D - j`), and the [`Pairing`] picks pairs so that most of them cover two primes.
//! Stage 2 accumulates `g = prod (x(R_m) - x(S_j))` over the chosen pairs, which only depend on
//! `b1` and `b2`: they are computed once in a [`PairPlan`], shared by all curves.
//!
//! Baby and giant steps are normalized to `z = 1` with Montgomery's batch inversion trick (one
//! modular inversion per batch, and four multiplications per point), so each pair costs one
//! subtraction and one multiplication. A non-invertible `z` means a point is at infinity
//! modulo a factor of `n`: that factor is returned right away. This also checks the primes
//! `l < D`, which are baby steps.

use crate::{
    arith::{Arith, PolyArith},
    cost::Costs,
    curve::{Curve, Xz},
    point::Point,
    stage2_poly::{self, PolyPlan},
};
use primal::Primes;
use rug::Integer;
use std::iter::Peekable;

/// Number of giant steps normalized with the same inversion.
const GIANT_BATCH: usize = 256;

/// Largest number of baby steps (they are all kept in memory).
const MAX_BABY_STEPS: usize = 1 << 17;

/// Largest pair table (in 64-bit words) kept in a [`PairPlan`]: larger ones are recomputed by
/// every curve, one batch of giant steps at a time.
const MAX_TABLE_WORDS: usize = 1 << 22;

/// Independent products accumulated in turn: their multiplications can overlap.
const ACCUMULATORS: usize = 4;

/// The primorials `P` usable as wheel base (`D = k*P`, `k < next`), with their largest prime
/// factor and the next prime `next`.
const PRIMORIALS: [(usize, usize, usize); 6] = [
    (6, 3, 5),
    (30, 5, 7),
    (210, 7, 11),
    (2310, 11, 13),
    (30030, 13, 17),
    (510_510, 17, 19),
];

/// Wheel of stage 2: the giant step `D` and its baby steps.
#[derive(Debug, Clone)]
struct Wheel {
    /// Giant step, a multiple of 6.
    d: usize,
    /// Position of `j` among the baby steps, for `j < D`, or `u32::MAX` if `gcd(j, D) > 1`.
    index: Vec<u32>,
    /// Number of baby steps: `phi(D)`.
    len: usize,
    /// 64-bit words of a row of the pair table (one bit per baby step).
    words: usize,
}

impl Wheel {
    fn new(d: usize) -> Self {
        let primes = prime_factors(d);
        let mut len = 0;
        let index = (0..d)
            .map(|j| {
                if j == 0 || primes.iter().any(|&p| j.is_multiple_of(p)) {
                    u32::MAX
                } else {
                    len += 1;
                    len as u32 - 1
                }
            })
            .collect();
        Wheel {
            d,
            index,
            len,
            words: len.div_ceil(64),
        }
    }
}

/// Distinct prime factors of `n`.
pub(crate) fn prime_factors(mut n: usize) -> Vec<usize> {
    let mut primes = Vec::new();
    let mut p = 2;
    while p * p <= n {
        if n.is_multiple_of(p) {
            primes.push(p);
            while n.is_multiple_of(p) {
                n /= p;
            }
        }
        p += 1;
    }
    if n > 1 {
        primes.push(n);
    }
    primes
}

/// Euler's totient.
pub(crate) fn phi(n: usize) -> usize {
    prime_factors(n)
        .into_iter()
        .fold(n, |phi, p| phi / p * (p - 1))
}

/// Giant step `D` minimizing the cost of stage 2 (in multiplications) outside of the pairs,
/// whose number barely depends on `D`: `D/3` point additions (6M) for the baby steps
/// `j = +-1 mod 6` and the normalization (4M) of the `phi(D)` kept ones, then one point addition
/// and one normalization per giant step.
///
/// The prime factors of `D` must be `<= b1`: the wheel skips the primes dividing `D`.
pub(crate) fn giant_step(b1: usize, b2: usize) -> usize {
    let mut best = (usize::MAX, 6);
    for (primorial, largest, next) in PRIMORIALS {
        if largest > b1 {
            break;
        }
        for k in 1..next {
            let d = k * primorial;
            let baby = phi(d);
            if baby > MAX_BABY_STEPS {
                break;
            }
            let giant = b2.saturating_sub(b1) / d + 2;
            let cost = 6 * (d / 3) + 4 * baby + 10 * giant;
            if cost < best.0 {
                best = (cost, d);
            }
        }
    }
    best.1
}

/// Chooses the pairs `(m, j)` checked by stage 2, one batch of giant steps at a time.
///
/// The numbers `+-j mod D` (for a `j < D/2`) are linked by the pairs into two chains:
/// `A_m = m*D - j` and `B_m = m*D + j` by `(m, j)`, and `B_m` and `A_(m+2)` by `(m + 1, D - j)`.
/// Going up each chain, an unchecked prime takes the pair to the next number of its chain,
/// which checks both when the next number is prime too: this greedy choice is a minimum set
/// of pairs covering the primes.
struct Pairing {
    d: usize,
    b2: usize,
    /// The primes in `(max(b1, D - 1), b2]` not yet visited.
    primes: Peekable<Primes>,
    /// Multiplier `m` of the current prime (the nearest multiple of `D` is `m*D`), and the
    /// largest number with this multiplier.
    m: usize,
    m_end: usize,
    /// By `j < D/2`: the `m` of the last `B_m` already checked by the pair `(m, j)`.
    checked_b: Vec<usize>,
    /// By `j < D/2` and parity of `m`: the `m` of the last `A_m` already checked by the pair
    /// `(m - 1, D - j)`.
    checked_a: Vec<[usize; 2]>,
    /// Pairs of the giant step after the current batch.
    next_row: Vec<u64>,
}

impl Pairing {
    fn new(wheel: &Wheel, b1: usize, b2: usize) -> Self {
        let (d, half) = (wheel.d, wheel.d / 2);
        let mut primes = Primes::all().peekable();
        // The primes l < D are baby steps.
        let lo = b1.max(d - 1);
        while primes.next_if(|&l| l <= lo).is_some() {}
        Pairing {
            d,
            b2,
            primes,
            m: 0,
            m_end: half,
            checked_b: vec![usize::MAX; half],
            checked_a: vec![[usize::MAX; 2]; half],
            next_row: vec![0; wheel.words],
        }
    }

    /// Range `[m_lo, m_hi]` of the giant steps of the pairs, empty if there is no prime to check.
    fn giant_steps(&mut self) -> (usize, usize) {
        let (d, half) = (self.d, self.d / 2);
        match self.primes.peek() {
            Some(&l) if l <= self.b2 => ((l + half) / d, (self.b2 + half) / d + 1),
            _ => (1, 0),
        }
    }

    /// Fills `rows` (`rows.len() / words` rows, zeroed) with the pairs of the giant steps from
    /// `m0` on, where `m0` follows the previous batch (or is the first giant step).
    fn fill(&mut self, wheel: &Wheel, m0: usize, rows: &mut [u64]) {
        let (d, half, words) = (self.d, self.d / 2, wheel.words);
        let m1 = m0 + rows.len() / words;
        let mut next_row = std::mem::take(&mut self.next_row);
        rows[..words].copy_from_slice(&next_row);
        next_row.fill(0);
        let mut set = |m: usize, j: usize| {
            let i = wheel.index[j] as usize;
            let word = if m < m1 {
                &mut rows[(m - m0) * words + i / 64]
            } else {
                &mut next_row[i / 64]
            };
            *word |= 1 << (i % 64);
        };
        let b2 = self.b2;
        while let Some(l) = self.primes.next_if(|&l| l <= b2 && l < m1 * d - half) {
            while l > self.m_end {
                self.m += 1;
                self.m_end += d;
            }
            let m = self.m;
            let base = m * d;
            if l < base {
                // A_m: its pair to B_m.
                let j = base - l;
                if self.checked_a[j][m % 2] != m {
                    set(m, j);
                    self.checked_b[j] = m;
                }
            } else {
                // B_m: its pair to A_(m+2).
                let j = l - base;
                if self.checked_b[j] != m {
                    set(m + 1, d - j);
                    self.checked_a[j][m % 2] = m + 2;
                }
            }
        }
        self.next_row = next_row;
    }
}

/// Everything stage 2 needs that only depends on the bounds (and the size of `n`), shared by
/// all the curves: which continuation to run, and its parameters.
#[derive(Debug, Clone)]
pub enum Stage2Plan {
    /// Baby-step giant-step continuation with prime pairing (this module): one multiplication
    /// per pair of primes, cheapest for small `b2`.
    Pairs(PairPlan),
    /// Polynomial continuation (product trees, Kronecker substitution): quasi-linear in
    /// `sqrt(b2)` per block, cheapest for large `b2`.
    Poly(PolyPlan),
}

impl Stage2Plan {
    /// Cheapest plan (according to [`Costs`]) for the primes in `(b1, b2]`, modulo numbers of the
    /// size of `n`. Requires `b1 >= 3`.
    pub fn new(n: &Integer, b1: usize, b2: usize) -> Self {
        assert!(b1 >= 3, "stage 2 requires b1 >= 3");
        let costs = Costs::new(n.significant_bits() as usize);
        let d = giant_step(b1, b2);
        let pairs = costs.pairs_stage2(b1, b2, d, phi(d));
        let (poly, poly_cost) = best_poly_plan(&costs, b1, b2);
        if poly_cost < pairs {
            Stage2Plan::Poly(poly)
        } else {
            Stage2Plan::pairs(b1, b2)
        }
    }

    /// Plan of the baby-step giant-step continuation. Requires `b1 >= 3`.
    pub fn pairs(b1: usize, b2: usize) -> Self {
        Stage2Plan::Pairs(PairPlan::new(b1, b2))
    }

    /// Plan of the polynomial continuation, for numbers of the size of `n`. Requires `b1 >= 3`.
    #[cfg_attr(not(any(test, feature = "bench")), allow(dead_code))]
    pub fn poly(n: &Integer, b1: usize, b2: usize) -> Self {
        assert!(b1 >= 3, "stage 2 requires b1 >= 3");
        let costs = Costs::new(n.significant_bits() as usize);
        Stage2Plan::Poly(best_poly_plan(&costs, b1, b2).0)
    }

    /// Largest `b2' >= b2` such that stage 2 checks every prime in `(b1, b2']`.
    #[cfg_attr(not(any(test, feature = "bench")), allow(dead_code))]
    pub fn b2(&self) -> usize {
        match self {
            Stage2Plan::Pairs(plan) => plan.b2,
            Stage2Plan::Poly(plan) => plan.b2_covered(),
        }
    }
}

/// Giant steps `d1` of the polynomial continuation: the multiples of 6 with more baby steps
/// (`phi(d1)/2`) than all the smaller ones (from GMP-ECM's `bestD`).
pub(crate) const POLY_GIANT_STEPS: [usize; 90] = [
    12, 18, 30, 42, 60, 90, 120, 150, 210, 240, 270, 330, 420, 510, 630, 840, 1050, 1260, 1470,
    1680, 1890, 2310, 2730, 3150, 3570, 3990, 4620, 5460, 6090, 6930, 8190, 9240, 10920, 12180,
    13860, 16170, 18480, 20790, 23100, 30030, 34650, 39270, 43890, 48510, 60060, 66990, 78540,
    90090, 99330, 120120, 133980, 150150, 180180, 210210, 240240, 270270, 300300, 334950, 371280,
    420420, 510510, 570570, 600600, 630630, 746130, 870870, 1021020, 1141140, 1291290, 1531530,
    1711710, 1891890, 2081310, 2312310, 2552550, 2852850, 3183180, 3573570, 3993990, 4594590,
    5105100, 5705700, 6322470, 7147140, 7987980, 8978970, 10210200, 11741730, 13123110, 14804790,
];

/// Largest memory (in bytes) used by the polynomial continuation.
const MAX_POLY_MEMORY: f64 = 256.0 * 1024.0 * 1024.0;

/// Cheapest polynomial continuation and its cost: the best giant step `d1`.
fn best_poly_plan(costs: &Costs, b1: usize, b2: usize) -> (PolyPlan, f64) {
    let mut best = (PolyPlan::shape_only(b1, b2, 6), f64::INFINITY);
    for d1 in [6].into_iter().chain(POLY_GIANT_STEPS) {
        if prime_factors(d1).iter().any(|&p| p > b1) {
            continue;
        }
        let plan = PolyPlan::shape_only(b1, b2, d1);
        let (_, df, giants) = plan.shape();
        // The product tree of F (a coefficient per leaf and level), and at the peak (measured)
        // about 44 more coefficients per leaf: the other polynomials, the Kronecker products
        // and GMP's scratch space.
        let coeffs = (usize::BITS - df.leading_zeros()) as f64 + 44.0;
        if coeffs * df as f64 * costs.elem_bytes() > MAX_POLY_MEMORY {
            break;
        }
        let cost = costs.poly_stage2(&plan);
        if cost < best.1 {
            best = (plan, cost);
        }
        if giants < df / 4 {
            // Larger giant steps only make F larger.
            break;
        }
    }
    let (_, cost) = best;
    let d1 = best.0.shape().0;
    (PolyPlan::new(b1, b2, d1), cost)
}

/// Stage 2 on the residues of `arith` with `plan`: returns `gcd(g, n)`, see
/// [`crate::ecm::stage2`].
pub fn stage2_with<A: PolyArith>(arith: A, q: &Point, plan: &Stage2Plan) -> Integer {
    match plan {
        Stage2Plan::Pairs(plan) => pairs_stage2_with(arith, q, plan),
        Stage2Plan::Poly(plan) => stage2_poly::stage2_with(arith, q, plan),
    }
}

/// Plan of the baby-step giant-step continuation: the wheel and the pairs to check.
#[derive(Debug, Clone)]
pub struct PairPlan {
    b1: usize,
    b2: usize,
    wheel: Wheel,
    /// Giant steps `m*D*Q` of the pairs, for `m` in `[m_lo, m_hi]` (empty if `m_lo > m_hi`).
    m_lo: usize,
    m_hi: usize,
    /// Pairs to check, one row of bits (by baby step) for each giant step, unless too large.
    table: Option<Vec<u64>>,
}

impl PairPlan {
    /// Plan of stage 2 for the primes in `(b1, b2]`. Requires `b1 >= 3`.
    pub fn new(b1: usize, b2: usize) -> Self {
        assert!(b1 >= 3, "stage 2 requires b1 >= 3");
        Self::with_giant_step(b1, b2, giant_step(b1, b2))
    }

    /// Plan with the giant step `d`: a multiple of 6 whose prime factors are `<= b1`.
    fn with_giant_step(b1: usize, b2: usize, d: usize) -> Self {
        let wheel = Wheel::new(d);
        let mut pairing = Pairing::new(&wheel, b1, b2);
        let (m_lo, m_hi) = pairing.giant_steps();
        let words = (m_hi + 1).saturating_sub(m_lo) * wheel.words;
        let table = (words <= MAX_TABLE_WORDS).then(|| {
            let mut table = vec![0; words];
            if words > 0 {
                pairing.fill(&wheel, m_lo, &mut table);
            }
            table
        });
        PairPlan {
            b1,
            b2,
            wheel,
            m_lo,
            m_hi,
            table,
        }
    }
}

/// Reusable buffers of the batch inversion.
pub(crate) struct Normalizer<E> {
    prefix: Vec<E>,
    inv: E,
    next: E,
    zinv: E,
    t: E,
}

impl<E: Clone> Normalizer<E> {
    pub(crate) fn new<A: Arith<Elem = E>>(a: &A, len: usize) -> Self {
        Normalizer {
            prefix: vec![a.zero(); len],
            inv: a.zero(),
            next: a.zero(),
            zinv: a.zero(),
            t: a.zero(),
        }
    }

    /// Replaces each `x[i]` by `x[i]/z[i]`, with a single modular inversion (Montgomery's trick):
    /// `3*len - 1` multiplications, plus `len` for the divisions.
    ///
    /// If some `z[i]` is not invertible, returns `Err(g)` with `g` a non-trivial factor of `n`
    /// found among the `z[i]` if possible, else `n`.
    pub(crate) fn normalize<A: Arith<Elem = E>>(
        &mut self,
        a: &A,
        x: &mut [E],
        z: &[E],
    ) -> Result<(), Integer> {
        let k = x.len();
        if k == 0 {
            return Ok(());
        }
        let prefix = &mut self.prefix[..k];
        prefix[0].clone_from(&z[0]);
        for i in 1..k {
            let (lo, hi) = prefix.split_at_mut(i);
            a.mul(&mut hi[0], &lo[i - 1], &z[i]);
        }
        let n = a.modulus();
        let Ok(inv) = a.to_integer(&prefix[k - 1]).invert(n) else {
            return Err(z
                .iter()
                .map(|z| a.gcd(z))
                .find(|g| *g != 1 && g != n)
                .unwrap_or_else(|| n.clone()));
        };
        self.inv = a.residue(&inv);
        // inv = 1/(z[0]*...*z[i]) at step i.
        for i in (1..k).rev() {
            a.mul(&mut self.zinv, &self.inv, &prefix[i - 1]);
            a.mul(&mut self.next, &self.inv, &z[i]);
            std::mem::swap(&mut self.inv, &mut self.next);
            a.mul(&mut self.t, &x[i], &self.zinv);
            std::mem::swap(&mut x[i], &mut self.t);
        }
        a.mul(&mut self.t, &x[0], &self.inv);
        std::mem::swap(&mut x[0], &mut self.t);
        Ok(())
    }
}

/// Baby-step giant-step stage 2 on the residues of `arith`: returns `gcd(g, n)`.
fn pairs_stage2_with<A: Arith>(arith: A, q: &Point, plan: &PairPlan) -> Integer {
    let curve = Curve::new(arith, &q.a_24);
    match accumulate(&curve, q, plan) {
        Ok(g) => curve.arith.gcd(&g),
        Err(g) => g,
    }
}

/// Normalized x-coordinates of the baby steps `j*Q` for `j < limit`, `j = +-1 mod 6`, with
/// `index[j] != u32::MAX` (their position in the result, of length `len`), or `Err(g)` with a
/// factor found by a failed inversion.
pub(crate) fn baby_steps<A: Arith>(
    curve: &Curve<A>,
    q: &Point,
    limit: usize,
    index: &[u32],
    len: usize,
    normalizer: &mut Normalizer<A::Elem>,
) -> Result<Vec<A::Elem>, Integer> {
    let a = &curve.arith;
    let mut scratch = curve.scratch();
    let mut xs = vec![a.zero(); len];
    let mut zs = vec![a.zero(); len];
    let mut store = |j: usize, p: &Xz<A::Elem>| {
        let i = index[j];
        if i != u32::MAX {
            xs[i as usize].clone_from(&p.x);
            zs[i as usize].clone_from(&p.z);
        }
    };

    // Q, 5Q and 6Q, then the two chains j = 1 and j = 5 mod 6: (j + 6)*Q = j*Q + 6*Q, with
    // difference (j - 6)*Q, which has the x-coordinate of (6 - j)*Q for the first step.
    let q1 = curve.point(&q.x_cord, &q.z_cord);
    let mut q2 = curve.infinity();
    curve.double(&mut q2, &q1, &mut scratch);
    let mut q3 = curve.infinity();
    curve.add(&mut q3, &q2, &q1, &q1, &mut scratch);
    let mut q5 = curve.infinity();
    curve.add(&mut q5, &q3, &q2, &q1, &mut scratch);
    let mut q6 = curve.infinity();
    curve.double(&mut q6, &q3, &mut scratch);
    let mut next = curve.infinity();
    for (first, diff, mut j) in [(q1.clone(), q5.clone(), 1), (q5, q1, 5)] {
        let (mut prev, mut cur) = (diff, first);
        while j < limit {
            store(j, &cur);
            if j + 6 >= limit {
                break;
            }
            curve.add(&mut next, &cur, &q6, &prev, &mut scratch);
            std::mem::swap(&mut prev, &mut cur);
            std::mem::swap(&mut cur, &mut next);
            j += 6;
        }
    }
    normalizer.normalize(a, &mut xs, &zs)?;
    Ok(xs)
}

/// Product `g` of stage 2, or `Err(g)` with a factor found by a failed inversion.
fn accumulate<A: Arith>(curve: &Curve<A>, q: &Point, plan: &PairPlan) -> Result<A::Elem, Integer> {
    let a = &curve.arith;
    let wheel = &plan.wheel;
    let (d, words) = (wheel.d, wheel.words);
    let mut normalizer = Normalizer::new(a, wheel.len.max(GIANT_BATCH));
    let one = a.residue(&Integer::from(1));
    if plan.b2 <= plan.b1 {
        return Ok(one);
    }

    // Also checks the primes l < D: l*Q = O modulo p makes z(l*Q) = 0 mod p.
    let baby = baby_steps(curve, q, d, &wheel.index, wheel.len, &mut normalizer)?;
    let (m_lo, m_hi) = (plan.m_lo, plan.m_hi);
    if m_lo > m_hi {
        return Ok(one);
    }

    let mut scratch = curve.scratch();
    let (xq, zq) = (a.factor(&q.x_cord), a.factor(&q.z_cord));
    let step = curve.ladder(&xq, &zq, &Integer::from(d));
    let mut r_prev = curve.ladder(&xq, &zq, &(Integer::from(m_lo) * d));
    let mut r = curve.ladder(&xq, &zq, &(Integer::from(m_lo + 1) * d));
    let mut r_next = curve.infinity();

    let batch = GIANT_BATCH.min(m_hi - m_lo + 1);
    let mut giant_x = vec![a.zero(); batch];
    let mut giant_z = vec![a.zero(); batch];
    // Without a table in the plan, the pairs are chosen batch by batch.
    let mut pairing = plan.table.is_none().then(|| {
        (
            Pairing::new(wheel, plan.b1, plan.b2),
            vec![0u64; batch * words],
        )
    });

    let mut acc: [A::Elem; ACCUMULATORS] = std::array::from_fn(|_| one.clone());
    let (mut diff, mut prod) = (a.zero(), a.zero());
    let mut turn = 0;

    for m0 in (m_lo..=m_hi).step_by(batch) {
        let len = batch.min(m_hi + 1 - m0);

        // Giant steps m0*D*Q, ..., (m0 + len - 1)*D*Q.
        for i in 0..len {
            giant_x[i].clone_from(&r_prev.x);
            giant_z[i].clone_from(&r_prev.z);
            curve.add(&mut r_next, &r, &step, &r_prev, &mut scratch);
            std::mem::swap(&mut r_prev, &mut r);
            std::mem::swap(&mut r, &mut r_next);
        }
        normalizer.normalize(a, &mut giant_x[..len], &giant_z[..len])?;

        let rows = match (&plan.table, &mut pairing) {
            (Some(table), _) => &table[(m0 - m_lo) * words..(m0 - m_lo + len) * words],
            (None, Some((pairing, rows))) => {
                let rows = &mut rows[..len * words];
                rows.fill(0);
                pairing.fill(wheel, m0, rows);
                rows
            }
            (None, None) => unreachable!(),
        };

        for (x, row) in giant_x.iter().zip(rows.chunks_exact(words)) {
            for (w, &bits) in row.iter().enumerate() {
                let mut bits = bits;
                while bits != 0 {
                    let j = w * 64 + bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    a.sub(&mut diff, x, &baby[j]);
                    a.mul(&mut prod, &acc[turn], &diff);
                    std::mem::swap(&mut acc[turn], &mut prod);
                    turn = (turn + 1) % ACCUMULATORS;
                }
            }
        }
    }

    let [mut g, rest @ ..] = acc;
    for x in rest {
        a.mul(&mut prod, &g, &x);
        std::mem::swap(&mut g, &mut prod);
    }
    Ok(g)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        arith::Mont,
        ecm::{curve, stage1, stage1_multiplier, Param},
    };

    #[test]
    fn wheel() {
        let w = Wheel::new(2310);
        assert_eq!(w.len, 480);
        assert_eq!(w.index[1], 0);
        assert_eq!(w.index[13], 1);
        assert_eq!(w.index[1155], u32::MAX);
        assert_eq!(w.index[2309], 479);
        assert_eq!(Wheel::new(6).len, 2);
    }

    #[test]
    fn auto_plan() {
        // The baby-step giant-step continuation for small b2, the polynomial one for large b2,
        // earlier for larger numbers; both cover at least (b1, b2].
        for (bits, b1, b2, poly) in [
            (128, 11_000, 1_873_422, false),
            (256, 11_000, 1_873_422, false),
            (256, 250_000, 128_992_510, true),
            (512, 50_000, 12_746_592, true),
            (1024, 3_000_000, 10_000_000_000, true),
            (2000, 50_000, 12_746_592, true),
        ] {
            let n = (Integer::from(1) << bits) - 1u32;
            let plan = Stage2Plan::new(&n, b1, b2);
            assert_eq!(
                matches!(plan, Stage2Plan::Poly(_)),
                poly,
                "{bits} {b1} {b2}"
            );
            assert!(plan.b2() >= b2);
            if let Stage2Plan::Poly(plan) = &plan {
                assert!(plan.b2_covered() < b2 + plan.shape().0);
            }
        }
    }

    #[test]
    fn giant_step_primes_below_b1() {
        for b1 in 3..100 {
            for b2 in [1, b1 + 1, 1000, 1_000_000, 100_000_000] {
                let d = giant_step(b1, b2);
                assert!(prime_factors(d).iter().all(|&p| p <= b1), "{b1} {b2} {d}");
            }
        }
        assert_eq!(giant_step(11_000, 1_873_422), 2310);
    }

    /// Every prime in `(max(b1, D - 1), b2]` is `m*D +- j` for a chosen pair `(m, j)`, whether
    /// the pairs are computed at once or batch by batch.
    #[test]
    fn pairs_cover_primes() {
        for (b1, b2) in [
            (3, 4),
            (6, 100),
            (10, 1000),
            (100, 10_001),
            (1000, 999),
            (2000, 147_396),
            (11_000, 1_873_422),
        ] {
            let plan = PairPlan::new(b1, b2);
            let (d, words) = (plan.wheel.d, plan.wheel.words);
            let table = plan.table.as_ref().unwrap();
            let mut covered = vec![false; b2 + 3 * d];
            let mut pairs = 0;
            for m in plan.m_lo..=plan.m_hi {
                for j in 1..d {
                    let i = plan.wheel.index[j];
                    if i != u32::MAX
                        && table[(m - plan.m_lo) * words + i as usize / 64] >> (i % 64) & 1 == 1
                    {
                        covered[m * d - j] = true;
                        covered[m * d + j] = true;
                        pairs += 1;
                    }
                }
            }
            let primes: Vec<usize> = Primes::all()
                .skip_while(|&l| l <= b1.max(d - 1))
                .take_while(|&l| l <= b2)
                .collect();
            for &l in &primes {
                assert!(covered[l], "{b1} {b2}: {l} not covered");
            }
            assert!(pairs <= primes.len());
            if primes.len() > 1000 {
                assert!(pairs * 10 < primes.len() * 8, "{pairs} pairs");
            }

            // Batch by batch.
            let mut pairing = Pairing::new(&plan.wheel, b1, b2);
            let mut m0 = plan.m_lo;
            while m0 <= plan.m_hi {
                let len = 7.min(plan.m_hi + 1 - m0);
                let mut rows = vec![0; len * words];
                pairing.fill(&plan.wheel, m0, &mut rows);
                let start = (m0 - plan.m_lo) * words;
                assert_eq!(rows, table[start..start + len * words]);
                m0 += len;
            }
        }
    }

    /// Every giant step the cost model can choose: `k*P` for a primorial `P` and `k` below the
    /// next prime, with at most [`MAX_BABY_STEPS`] baby steps.
    fn all_giant_steps() -> Vec<usize> {
        PRIMORIALS
            .iter()
            .flat_map(|&(primorial, _, next)| (1..next).map(move |k| k * primorial))
            .filter(|&d| phi(d) <= MAX_BABY_STEPS)
            .collect()
    }

    /// Numbers `<= b2` checked by `plan`: the baby steps `j < D` (a point at infinity makes the
    /// normalization fail) and both numbers `m*D +- j` of each pair.
    fn covered(plan: &PairPlan) -> Vec<bool> {
        let (d, words) = (plan.wheel.d, plan.wheel.words);
        let mut covered = vec![false; plan.b2.max(d) + 3 * d];
        for (j, &i) in plan.wheel.index.iter().enumerate() {
            covered[j] = i != u32::MAX;
        }
        let table = plan.table.as_ref().unwrap();
        for m in plan.m_lo..=plan.m_hi {
            assert!(m >= 1);
            for j in 1..d {
                let i = plan.wheel.index[j];
                if i != u32::MAX
                    && table[(m - plan.m_lo) * words + i as usize / 64] >> (i % 64) & 1 == 1
                {
                    covered[m * d - j] = true;
                    covered[m * d + j] = true;
                }
            }
        }
        covered
    }

    #[test]
    fn plans_cover_primes() {
        // Every giant step, with tiny bounds, b2 <= b1, b2 < D, b2 not a multiple of D, and
        // primes next to b1, b2 and multiples of D.
        let mut rand = rug::rand::RandState::new();
        let mut random = |max: usize| {
            rug::Integer::from(max)
                .random_below(&mut rand)
                .to_usize()
                .unwrap()
        };
        for d in all_giant_steps() {
            let b1_min = *prime_factors(d).last().unwrap();
            let mut bounds = vec![(b1_min, 4), (b1_min, b1_min), (b1_min, b1_min + 1)];
            for b1 in [b1_min, b1_min + 1, d - 1, d, d + 1, 3 * d / 2, 5 * d + 1] {
                for b2 in [
                    b1 + 1,
                    b1 + 2,
                    d - 2,
                    d - 1,
                    d,
                    d + 1,
                    2 * d - 1,
                    2 * d + 1,
                    7 * d,
                ] {
                    bounds.push((b1, b2));
                }
                bounds.push((b1, b1 + random(20 * d) + 1));
            }
            for (b1, b2) in bounds {
                let plan = PairPlan::with_giant_step(b1, b2, d);
                let covered = covered(&plan);
                for l in Primes::all()
                    .skip_while(|&l| l <= b1)
                    .take_while(|&l| l <= b2)
                {
                    assert!(covered[l], "D = {d}, ({b1}, {b2}]: {l} not covered");
                }
                let mut pairing = Pairing::new(&plan.wheel, b1, b2);
                let table = plan.table.as_ref().unwrap();
                let mut m0 = plan.m_lo;
                while m0 <= plan.m_hi {
                    let len = (1 + m0 % 5).min(plan.m_hi + 1 - m0);
                    let mut rows = vec![0; len * plan.wheel.words];
                    pairing.fill(&plan.wheel, m0, &mut rows);
                    let start = (m0 - plan.m_lo) * plan.wheel.words;
                    assert_eq!(
                        rows,
                        table[start..start + rows.len()],
                        "D = {d}, ({b1}, {b2}]"
                    );
                    m0 += len;
                }
            }
        }
    }

    #[test]
    fn stage2_all_giant_steps() {
        // For each giant step: stage 2 finds p whenever l*Q = O modulo p for a prime l in
        // (b1, b2], whether l < D (baby step) or not (pair). The larger wheels are only covered by
        // `plans_cover_primes`: checking l*Q = O for every l would be too slow.
        let n = Integer::from(4_009_823u64) * Integer::from(99_476_569u64);
        let mut found = 0;
        for d in all_giant_steps().into_iter().filter(|&d| d < 30030) {
            let b1 = (*prime_factors(d).last().unwrap()).max(30);
            let b2 = b1 + 3 * d + 1000;
            let plan = PairPlan::with_giant_step(b1, b2, d);
            let streamed = PairPlan {
                table: None,
                ..plan.clone()
            };
            let k = stage1_multiplier(b1);
            let primes: Vec<Integer> = Primes::all()
                .skip_while(|&l| l <= b1)
                .take_while(|&l| l <= b2)
                .map(Integer::from)
                .collect();
            for sigma in 2..30 {
                let q = stage1(
                    &curve(&n, Param::Batch2, &Integer::from(sigma)).unwrap(),
                    &k,
                );
                if q.z_cord.clone().gcd(&n) != 1 {
                    continue;
                }
                let g = pairs_stage2_with(Mont::<1>::new(&n), &q, &plan);
                assert_eq!(g, pairs_stage2_with(Mont::<1>::new(&n), &q, &streamed));
                let expected = primes.iter().any(|l| {
                    let g = q.mont_ladder(l).z_cord.gcd(&n);
                    g != 1 && g != n
                });
                if expected {
                    assert_ne!(g, 1, "D = {d}, sigma = {sigma}");
                    found += 1;
                }
            }
        }
        assert!(found > 20, "{found}");
    }

    #[test]
    fn normalize() {
        let (p, q) = (Integer::from(4_009_823u64), Integer::from(99_476_569u64));
        let n = Integer::from(&p * &q);
        let a = Mont::<1>::new(&n);
        let mut normalizer = Normalizer::new(&a, 4);
        let residues = |v: &[Integer]| v.iter().map(|v| a.residue(v)).collect::<Vec<_>>();
        let z = [3, 5, 7, 11].map(Integer::from);
        let mut x = residues(&[6, 10, 14, 22].map(Integer::from));
        normalizer.normalize(&a, &mut x, &residues(&z)).unwrap();
        assert!(x.iter().all(|x| a.to_integer(x) == 2));
        normalizer.normalize(&a, &mut [], &[]).unwrap();

        // A failed inversion returns a proper factor among the z's, else n.
        let mut x = residues(&[1, 1, 1].map(Integer::from));
        for (z, g) in [
            ([Integer::from(3), q.clone() * 2, Integer::from(5)], &q),
            ([Integer::ZERO, Integer::from(3), p.clone()], &p),
            ([p.clone(), q.clone(), Integer::from(1)], &p),
            ([Integer::from(3), Integer::ZERO, Integer::from(1)], &n),
        ] {
            assert_eq!(
                normalizer.normalize(&a, &mut x, &residues(&z)),
                Err(g.clone())
            );
        }
    }

    #[test]
    fn pairs_batch_by_batch() {
        // Without a table, the pairs are recomputed for each batch of giant steps: same product.
        let n = Integer::from(4_009_823u64) * Integer::from(99_476_569u64);
        let (b1, b2) = (1000, 1_000_000);
        let k = stage1_multiplier(b1);
        let plan = PairPlan::new(b1, b2);
        assert!(plan.m_hi - plan.m_lo > 2 * GIANT_BATCH);
        let streamed = PairPlan {
            table: None,
            ..plan.clone()
        };
        for sigma in 2..10 {
            let q = stage1(
                &curve(&n, Param::Batch2, &Integer::from(sigma)).unwrap(),
                &k,
            );
            let curve = Curve::new(Mont::<1>::new(&n), &q.a_24);
            assert_eq!(
                accumulate(&curve, &q, &plan),
                accumulate(&curve, &q, &streamed)
            );
        }
    }
}
