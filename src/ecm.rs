use crate::{
    arith::{with_arith, Arith, Factor},
    curve::{Curve, Point},
    primes::primes,
    stage2::{stage2_with, Stage2Plan},
};
#[cfg(feature = "progress-bar")]
use indicatif::ProgressBar;
use rug::{
    integer::IsPrime,
    rand::{RandGen, RandState},
    Integer,
};
use std::collections::HashMap;

/// Error occured during ecm factorization.
#[derive(thiserror::Error, Debug)]
pub enum Error {
    /// Bounds should be an even integer.
    #[error("Bounds should be an even integer")]
    BoundsNotEven,
    /// Too small bounds.
    #[error("Too small bounds")]
    BoundsTooSmall,
    /// The factorization failed.
    #[error("The factorization failed")]
    ECMFailed,
    /// The number is prime.
    #[error("The number is prime")]
    NumberIsPrime,
}

/// Number of rounds of the probabilistic primality test.
const PRIMALITY_REPS: u32 = 25;

/// Returns one factor of n using Lenstra's 2 Stage Elliptic curve Factorization,
/// with the curves of GMP-ECM's parametrization 2 (`-param 2`). Here Montgomery
/// curves and Montgomery modular arithmetic are used for fast computation of addition and
/// doubling of points.
///
/// This ECM method considers elliptic curves in Montgomery form (E : b*y^2*z = x^3 + a*x^2*z + x*z^2)
/// and involves elliptic curve operations (mod N), where the elements in Z are reduced (mod N).
/// Since N is not a prime, E over FF(N) is not really an elliptic curve but we can still do point additions
/// and doubling as if FF(N) was a field.
///
/// Stage 1: The basic algorithm involves taking a random point (P) on an elliptic curve in FF(N).
/// The compute k*P using Montgomery ladder algorithm.
/// Let q be an unknown factor of N. Then the order of the curve E, |E(FF(q))|,
/// might be a smooth number that divides k. Then we have k = l * |E(FF(q))|
/// for some l. For any point belonging to the curve E, |E(FF(q))|*P = O,
/// hence k*P = l*|E(FF(q))|*P. Thus kP.z = 0 (mod q), and the unknown factor of N (q)
/// can be recovered by taking gcd(kP.z, N).
///
/// Stage 2: This is a continuation of Stage 1 if k*P != O. The idea is to utilize
/// the fact that even if kP != 0, the value of k might miss just one large prime divisor
/// of |E(FF(q))|. In this case, we only need to compute the scalar multiplication by p
/// to get p*k*P = O. Here a second bound B2 restricts the size of possible values of p.
///
/// Parameters:
///
/// - `n`: Number to be factored.
/// - `B1`: Stage 1 Bound.
/// - `B2`: Stage 2 Bound.
/// - `max_curve`: Maximum number of curves generated.
/// - `rgen`: Random number generator.
pub fn ecm_one_factor(
    n: &Integer,
    b1: usize,
    b2: usize,
    max_curve: usize,
    rgen: &mut RandState<'_>,
    #[cfg(feature = "progress-bar")] pb: Option<&ProgressBar>,
) -> Result<Integer, Error> {
    if !b1.is_multiple_of(2) || !b2.is_multiple_of(2) {
        return Err(Error::BoundsNotEven);
    }

    // Stage 2 skips the primes of its wheel, which must be at most `b1`.
    if b1 < 6 || b2 < 4 {
        return Err(Error::BoundsTooSmall);
    }

    // BPSW only (rug runs `reps - 24` Miller-Rabin rounds on top of it): no composite is known
    // to pass it, and the caller usually already knows that `n` is composite.
    if n.is_probably_prime(PRIMALITY_REPS) != IsPrime::No {
        return Err(Error::NumberIsPrime);
    }

    #[cfg(feature = "progress-bar")]
    if let Some(pb) = pb {
        pb.set_length(max_curve as u64);
        pb.set_position(0);
    }

    let k = stage1_multiplier(b1);
    let plan = Stage2Plan::new(n, b1, b2);
    let param = Param::default();

    for _ in 0..max_curve {
        #[cfg(feature = "progress-bar")]
        if let Some(pb) = pb {
            pb.inc(1);
        }

        let sigma = random_sigma(n, param, rgen);
        match run_curve(n, param, &sigma, &k, &plan) {
            CurveOutcome::Setup(g) | CurveOutcome::Stage1(g) | CurveOutcome::Stage2(g) => {
                return Ok(g)
            }
            CurveOutcome::Failed => {}
        }
    }

    // ECM failed, Increase the bounds
    Err(Error::ECMFailed)
}

/// Random state seeded with `seed`, to draw the curves.
///
/// Seeding GMP's default generator (a Mersenne Twister) costs about 0.26 ms, a modular
/// exponentiation with a 20000-bit modulus, which is more than factoring a small number takes.
pub(crate) fn rand_state(seed: usize) -> RandState<'static> {
    RandState::new_custom_boxed(Box::new(SplitMix64(seed as u64)))
}

/// Steele, Lea and Flood's SplitMix64 generator: fast, and good enough to draw curves.
struct SplitMix64(u64);

impl RandGen for SplitMix64 {
    fn gen(&mut self) -> u32 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        ((z ^ (z >> 31)) >> 32) as u32
    }
}

/// Families of curves, named after GMP-ECM's `-param` values.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Param {
    /// Suyama's parametrization (GMP-ECM `-param 0`), see [`suyama_curve`]: `sigma` in
    /// `[6, n - 1]`.
    #[cfg_attr(not(any(test, feature = "bench")), allow(dead_code))]
    Suyama,
    /// GMP-ECM's default parametrization (`-param 1`), see [`square_curve`]: `sigma` in
    /// `[2, 2^32)`. The cheapest stage 1, but about 1.3 times more curves than with the others
    /// are needed to find a factor.
    #[cfg_attr(not(any(test, feature = "bench")), allow(dead_code))]
    Square,
    /// GMP-ECM's `-param 2`, see [`batch2_curve`]: `sigma` in `[2, 2^64)`. The default: as
    /// effective as Suyama's curves, with a stage 1 almost as cheap as with `Square`.
    #[default]
    Batch2,
}

/// Random `sigma` for the curves of `param`, in the range documented by [`Param`].
pub fn random_sigma(n: &Integer, param: Param, rgen: &mut RandState<'_>) -> Integer {
    match param {
        Param::Suyama => Integer::from(n - 6).random_below(rgen) + 6,
        Param::Square => Integer::from((1u64 << 32) - 2).random_below(rgen) + 2,
        Param::Batch2 => Integer::from(u64::MAX - 1).random_below(rgen) + 2,
    }
}

/// Starting point of the curve of `param` given by `sigma`.
///
/// Returns `Err(g)` when the curve cannot be built, where `g` is a factor of `n` found while
/// building it (`g` may be `1` or `n`: then the curve is just unusable).
pub fn curve(n: &Integer, param: Param, sigma: &Integer) -> Result<Point, Integer> {
    // Montgomery arithmetic needs an odd modulus.
    if n.is_even() {
        return Err(Integer::from(2));
    }
    match param {
        Param::Suyama => suyama_curve(n, sigma),
        Param::Square => square_curve(n, sigma),
        Param::Batch2 => batch2_curve(n, sigma),
    }
}

/// Result of running a single ECM curve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CurveOutcome {
    /// The curve could not be built, and the gcd found while building it is returned.
    Setup(Integer),
    /// Stage 1 found a non-trivial factor.
    Stage1(Integer),
    /// Stage 2 found a non-trivial factor.
    Stage2(Integer),
    /// No factor was found with this curve.
    Failed,
}

/// Runs one ECM curve (curve setup, stage 1 and stage 2) of `param` with the given `sigma`.
///
/// `k` must be the stage 1 multiplier returned by [`stage1_multiplier`] for the `b1` of `plan`.
pub fn run_curve(
    n: &Integer,
    param: Param,
    sigma: &Integer,
    k: &Integer,
    plan: &Stage2Plan,
) -> CurveOutcome {
    let q = match curve(n, param, sigma) {
        Ok(q) => q,
        // If g = 1 or n, try another curve
        Err(g) if g == 1 || &g == n => return CurveOutcome::Failed,
        Err(g) => return CurveOutcome::Setup(g),
    };

    let q = stage1(&q, k);
    let g = q.z.clone().gcd(n);

    // Stage 1 factor
    if &g != n && g != 1 {
        return CurveOutcome::Stage1(g);
    }

    // Stage 1 failure. Q.z = 0, Try another curve
    if &g == n {
        return CurveOutcome::Failed;
    }

    let g = stage2(&q, plan);

    // Stage 2 Factor found
    if &g != n && g != 1 {
        return CurveOutcome::Stage2(g);
    }

    CurveOutcome::Failed
}

/// Stage 1 multiplier: product of the largest powers of all primes `p <= b1` that are `<= b1`.
pub fn stage1_multiplier(b1: usize) -> Integer {
    prime_power_product(1, b1)
}

/// `E(hi)/E(lo)`, where `E(b)` is the product of the largest powers of the primes `p <= b` that
/// are `<= b` (so `E(b) = stage1_multiplier(b)`), for `1 <= lo <= hi`.
///
/// Computed with a product tree: quasi-linear in the size of the result.
pub(crate) fn prime_power_product(lo: usize, hi: usize) -> Integer {
    product(prime_power_words(lo, hi).map(Integer::from).collect())
}

/// `E(hi)/E(lo)` (see [`prime_power_product`]) as a sequence of 64-bit factors, from the
/// smallest primes to the largest.
pub(crate) fn prime_power_words(lo: usize, hi: usize) -> impl Iterator<Item = u64> {
    let lo = lo.max(1);
    let mut word = 1u64;
    let mut powers = primes(hi)
        // Above sqrt(hi), only the primes in (lo, hi] contribute (to the power 1).
        .filter(move |&p| p > lo || p.saturating_mul(p) <= hi)
        .flat_map(move |p| {
            let old = if p <= lo { lo.ilog(p) } else { 0 };
            std::iter::repeat_n(p as u64, (hi.ilog(p) - old) as usize)
        });
    let mut done = false;
    std::iter::from_fn(move || {
        if done {
            return None;
        }
        for q in powers.by_ref() {
            match word.checked_mul(q) {
                Some(w) => word = w,
                None => return Some(std::mem::replace(&mut word, q)),
            }
        }
        done = true;
        Some(word)
    })
}

/// Product of `values`, by a balanced product tree.
pub(crate) fn product(mut values: Vec<Integer>) -> Integer {
    while values.len() > 1 {
        let mut next = Vec::with_capacity(values.len().div_ceil(2));
        let mut it = values.into_iter();
        while let Some(a) = it.next() {
            next.push(match it.next() {
                Some(b) => a * b,
                None => a,
            });
        }
        values = next;
    }
    values.pop().unwrap_or_else(|| Integer::from(1))
}

/// Builds the starting point of a curve using Suyama's parametrization.
///
/// Returns `Err(g)` with `g = gcd(2*u^3*v, n)` when the curve cannot be built (`g` may be `n`).
pub fn suyama_curve(n: &Integer, sigma: &Integer) -> Result<Point, Integer> {
    let three = Integer::from(3);
    let u = (Integer::from(sigma * sigma) - 5u32) % n;
    let v = Integer::from(sigma * 4u32) % n;
    let u_3 = u.clone().pow_mod(&three, n).unwrap();

    // We use the elliptic curve y^2 = x^3 + a*x^2 + x
    // where a = (v - u)^3 * (3*u + v) / (4*u^3*v) - 2
    // However, we do not declare a because it is more convenient
    // to use a24 = (a + 2) / 4 in the calculation.
    let u_3_v = Integer::from(&u_3 * &v);
    let a24 = match Integer::from(&u_3_v * 16u32).invert(n) {
        Ok(inv) => Integer::from(&v - &u).pow_mod(&three, n).unwrap() * (3u32 * u + &v) * inv % n,
        // If the invert(16*u^3*v, n) doesn't exist (i.e., g != 1)
        Err(_) => return Err((u_3_v * 2u32).gcd(n)),
    };

    let v_3 = v.pow_mod(&three, n).unwrap();
    Ok(Point {
        x: u_3,
        z: v_3,
        a24,
        n: n.clone(),
    })
}

/// Builds the starting point of a curve using GMP-ECM's parametrization 1
/// (`ECM_PARAM_BATCH_SQUARE`): the curve `b*y^2 = x^3 + a*x^2 + x` with `a = 4*d - 2`,
/// `b = 16*d + 2` and `d = sigma^2/2^64 mod n`, and the point `(2 : 1)`.
///
/// Its stage 1 is cheaper than with other curves: `(a + 2)/4 = d` has the one-limb Montgomery
/// form `sigma^2` for `sigma < 2^32`, and the coordinates of the starting point are small.
///
/// `n` must be odd. Returns `Err(n)` when the curve is singular (`d = 0` or `d = 1`).
pub fn square_curve(n: &Integer, sigma: &Integer) -> Result<Point, Integer> {
    let inv = Integer::from(Integer::u_pow_u(2, 64))
        .invert(n)
        .expect("n must be odd");
    let d = Integer::from(sigma * sigma) * inv % n;
    if d == 0 || d == 1 {
        return Err(n.clone());
    }
    Ok(Point::start(d, n))
}

/// Builds the starting point of a curve using GMP-ECM's parametrization 2
/// (`ECM_PARAM_BATCH_2`): the point `(2 : 1)` on the curve `b*y^2 = x^3 + a*x^2 + x` with
/// `a = -(3*x3^4 + 6*x3^2 - 1)/(4*x3^3)` and `x3 = (3*x + y + 6)/(2*(y - 3))`, where
/// `(x, y) = sigma*(-3, 3)` on `y^2 = x^3 + 36`.
///
/// These curves have the same expected torsion (so the same probability of success) as
/// Suyama's, and a starting point with small coordinates: their stage 1 costs one full
/// multiplication per ladder step less than with Suyama's parametrization.
///
/// Requires `sigma >= 2`. Returns `Err(g)` with a factor `g` of `n` found by a failed
/// inversion (`g` may be `n`).
pub fn batch2_curve(n: &Integer, sigma: &Integer) -> Result<Point, Integer> {
    if *sigma < 2 {
        return Err(n.clone());
    }
    let md = |x: Integer| -> Integer {
        let mut x = x % n;
        if x < 0 {
            x += n;
        }
        x
    };
    let inv =
        |x: &Integer| -> Result<Integer, Integer> { x.clone().invert(n).map_err(|x| x.gcd(n)) };

    let (x, y, z) = with_arith!(n, |arith| batch2_multiple(arith, sigma));

    // Affine coordinates.
    let z_inv = inv(&z)?;
    let z_inv2 = md(z_inv.clone().square());
    let x = md(x * &z_inv2);
    let y = md(y * z_inv2 * z_inv);

    let x3 = md((3 * x + &y + 6) * inv(&md(2 * (y - 3)))?);
    let x3_2 = md(x3.clone().square());
    // a = -(3*x3^4 + 6*x3^2 - 1)/(4*x3^3), and a24 = (a + 2)/4
    let x3_4: Integer = x3_2.clone().square();
    let numerator: Integer = x3_4 * 3u32 + Integer::from(&x3_2 * 6u32) - 1u32;
    let a = md(-numerator * inv(&md(4 * x3_2 * x3))?);
    let a24 = md((a + 2) * inv(&Integer::from(4))?);
    Ok(Point::start(a24, n))
}

/// `sigma*(-3 : 3 : 1)` in Jacobian coordinates, on the curve `y^2 = x^3 + 36`.
fn batch2_multiple<A: Arith>(a: A, sigma: &Integer) -> (Integer, Integer, Integer) {
    let (px, py) = (a.residue(&Integer::from(-3)), a.residue(&Integer::from(3)));
    let (mut x, mut y, mut z) = (px.clone(), py.clone(), a.residue(&Integer::from(1)));
    let [mut t0, mut t1, mut t2, mut t3, mut t4, mut t5, mut t6, mut t7] =
        std::array::from_fn(|_| a.zero());
    for bit in (0..sigma.significant_bits() - 1).rev() {
        // Doubling, "dbl-2009-l" (a = 0).
        a.sqr(&mut t0, &x); // A = x^2
        a.sqr(&mut t1, &y); // B = y^2
        a.sqr(&mut t2, &t1); // C = B^2
        a.add(&mut t3, &x, &t1);
        a.sqr(&mut t4, &t3);
        a.sub(&mut t3, &t4, &t0);
        a.sub(&mut t4, &t3, &t2);
        a.add(&mut t3, &t4, &t4); // D = 2*((x + B)^2 - A - C)
        a.add(&mut t4, &t0, &t0);
        a.add(&mut t5, &t4, &t0); // E = 3*A
        a.sqr(&mut t0, &t5); // F = E^2
        a.sub(&mut t4, &t0, &t3);
        a.sub(&mut t0, &t4, &t3); // x3 = F - 2*D
        a.mul(&mut t4, &y, &z);
        a.add(&mut z, &t4, &t4); // z3 = 2*y*z
        a.sub(&mut t4, &t3, &t0);
        a.mul(&mut t3, &t5, &t4);
        a.add(&mut t4, &t2, &t2);
        a.add(&mut t2, &t4, &t4);
        a.add(&mut t4, &t2, &t2);
        a.sub(&mut y, &t3, &t4); // y3 = E*(D - x3) - 8*C
        std::mem::swap(&mut x, &mut t0);

        if sigma.get_bit(bit) {
            // Mixed addition of (-3, 3), "madd-2007-bl".
            a.sqr(&mut t0, &z); // Z1Z1 = z^2
            a.mul(&mut t1, &px, &t0); // U2 = px*Z1Z1
            a.mul(&mut t2, &py, &z);
            a.mul(&mut t3, &t2, &t0); // S2 = py*z*Z1Z1
            a.sub(&mut t2, &t1, &x); // H = U2 - x
            a.sqr(&mut t1, &t2); // HH = H^2
            a.add(&mut t4, &t1, &t1);
            a.add(&mut t5, &t4, &t4); // I = 4*HH
            a.mul(&mut t4, &t2, &t5); // J = H*I
            a.sub(&mut t6, &t3, &y);
            a.add(&mut t3, &t6, &t6); // r = 2*(S2 - y)
            a.mul(&mut t6, &x, &t5); // V = x*I
            a.sqr(&mut t5, &t3);
            a.sub(&mut t7, &t5, &t4);
            a.sub(&mut t5, &t7, &t6);
            a.sub(&mut x, &t5, &t6); // x3 = r^2 - J - 2*V
            a.sub(&mut t5, &t6, &x);
            a.mul(&mut t6, &t3, &t5);
            a.mul(&mut t5, &y, &t4);
            a.add(&mut t7, &t5, &t5);
            a.sub(&mut y, &t6, &t7); // y3 = r*(V - x3) - 2*y*J
            a.add(&mut t5, &z, &t2);
            a.sqr(&mut t6, &t5);
            a.sub(&mut t5, &t6, &t0);
            a.sub(&mut z, &t5, &t1); // z3 = (z + H)^2 - Z1Z1 - HH
        }
    }
    (a.to_integer(&x), a.to_integer(&y), a.to_integer(&z))
}

/// Stage 1: computes `k*P`, for `k >= 1`.
pub fn stage1(p: &Point, k: &Integer) -> Point {
    with_arith!(&p.n, |arith| stage1_with(arith, p, k))
}

fn stage1_with<A: Arith>(arith: A, p: &Point, k: &Integer) -> Point {
    let n: &Integer = &p.n;
    // Normalizing P to z = 1 saves a multiplication per ladder step. If z is not invertible, P
    // is the point at infinity modulo a factor of n, and so is k*P: P has the same gcd.
    let Ok(z_inv) = p.z.clone().invert(n) else {
        return p.clone();
    };
    let x = Integer::from(&p.x * &z_inv) % n;

    let curve = Curve::new(arith, &p.a24);
    let q = curve.ladder(&curve.arith.factor(&x), &Factor::One, k);
    Point {
        x: curve.arith.to_integer(&q.x),
        z: curve.arith.to_integer(&q.z),
        a24: p.a24.clone(),
        n: n.clone(),
    }
}

/// Stage 2: baby-step giant-step standard continuation, with baby and giant steps normalized
/// by batch inversions and prime pairing (one multiplication per pair of primes `m*D +- j`).
///
/// Returns `gcd(g, n)` where `g` is the accumulated product over the primes in `(b1, b2]` of
/// `plan` (or a factor found when normalizing the points, possibly `n`), and `1` if `b2 <= b1`.
pub fn stage2(q: &Point, plan: &Stage2Plan) -> Integer {
    with_arith!(&q.n, |arith| stage2_with(arith, q, plan))
}

/// Trial division removes the prime factors below this bound.
const TRIAL_DIVISION_BOUND: usize = 1 << 16;

/// Removes the prime factors of `n` below 2^16.
///
/// Returns the found factors with their multiplicity, and the remaining cofactor.
pub fn trial_division(n: &Integer) -> (HashMap<Integer, usize>, Integer) {
    let mut factors = HashMap::new();
    let mut n: Integer = n.clone();
    for prime in primes(TRIAL_DIVISION_BOUND - 1) {
        if n < prime * prime {
            // n is 1 or a prime.
            break;
        }
        if n.is_divisible_u(prime as u32) {
            let mut multiplicity = 0;
            while n.is_divisible_u(prime as u32) {
                n.div_exact_u_mut(prime as u32);
                multiplicity += 1;
            }
            factors.insert(Integer::from(prime), multiplicity);
        }
    }
    (factors, n)
}

/// Performs factorization using Lenstra's Elliptic curve method.
///
/// Small factors are removed by trial division, then the factors are found from the smallest
/// to the largest, as GMP-ECM recommends: curves with a first bound `B1` optimal for factors of
/// 10, 15, 20, ... digits in turn, each size for the expected number of curves to find such a
/// factor, and Pollard's P-1 method with larger bounds before each size. The time to find a
/// factor thus depends on its size much more than on the size of `n`.
///
/// Deterministic: the curves are drawn from a fixed seed.
///
/// # Parameters
///
/// - `n`: Number to be factored.
///
/// # Panics
///
/// If `n` is not positive.
pub fn ecm(
    n: &Integer,
    #[cfg(feature = "progress-bar")] pb: Option<&ProgressBar>,
) -> Result<HashMap<Integer, usize>, Error> {
    crate::driver::factor(
        n,
        1234,
        #[cfg(feature = "progress-bar")]
        pb,
    )
}

/// Performs factorization using Lenstra's Elliptic curve method, with fixed bounds.
///
/// This function repeatedly calls `ecm_one_factor` to compute the factors
/// of n. First all the small factors are taken out using trial division.
/// Then `ecm_one_factor` is used to compute one factor at a time.
///
/// # Parameters
///
/// - `n`: Number to be factored.
/// - `B1`: Stage 1 Bound.
/// - `B2`: Stage 2 Bound.
/// - `max_curve`: Maximum number of curves generated.
/// - `seed`: Initialize pseudorandom generator.
///
/// # Panics
///
/// If `n` is not positive.
pub fn ecm_with_params(
    n: &Integer,
    b1: usize,
    b2: usize,
    max_curve: usize,
    seed: usize,
    #[cfg(feature = "progress-bar")] pb: Option<&ProgressBar>,
) -> Result<HashMap<Integer, usize>, Error> {
    if !b1.is_multiple_of(2) || !b2.is_multiple_of(2) {
        return Err(Error::BoundsNotEven);
    }
    if b1 < 6 || b2 < 4 {
        return Err(Error::BoundsTooSmall);
    }
    assert!(*n > 0, "only positive numbers can be factored");

    let (mut factors, n) = trial_division(n);

    let mut rand_state = rand_state(seed);

    // Composite factors left to split, with the multiplicity they have in the original number.
    let mut queue = Vec::new();
    sort_factor(n, 1, (), &mut factors, &mut queue);

    while let Some((n, exponent, ())) = queue.pop() {
        let factor = ecm_one_factor(
            &n,
            b1,
            b2,
            max_curve,
            &mut rand_state,
            #[cfg(feature = "progress-bar")]
            pb,
        )?;

        // `factor` may itself be composite: both parts go through `sort_factor` again.
        let mut cofactor = n;
        let mut multiplicity = 0;
        while cofactor.is_divisible(&factor) {
            cofactor /= &factor;
            multiplicity += 1;
        }
        sort_factor(
            factor,
            exponent * multiplicity,
            (),
            &mut factors,
            &mut queue,
        );
        sort_factor(cofactor, exponent, (), &mut factors, &mut queue);
    }

    Ok(factors)
}

/// Records `n^exponent`: a prime goes to `factors`, a perfect power is reduced to its root, and
/// any other composite is queued for factorization, with `state`.
pub(crate) fn sort_factor<T>(
    n: Integer,
    exponent: usize,
    state: T,
    factors: &mut HashMap<Integer, usize>,
    queue: &mut Vec<(Integer, usize, T)>,
) {
    if n == 1 {
        return;
    }
    if n.is_probably_prime(PRIMALITY_REPS) != IsPrime::No {
        *factors.entry(n).or_insert(0) += exponent;
        return;
    }
    match perfect_power(&n) {
        Some((root, power)) => sort_factor(root, exponent * power as usize, state, factors, queue),
        None => queue.push((n, exponent, state)),
    }
}

/// Writes `n` as `root^power` with the smallest possible `root`, if it is a perfect power.
fn perfect_power(n: &Integer) -> Option<(Integer, u32)> {
    if !n.is_perfect_power() {
        return None;
    }
    for power in primes(n.significant_bits() as usize) {
        let (root, remainder) = n.clone().root_rem(Integer::new(), power as u32);
        if remainder == 0 {
            return Some((root, power as u32));
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use primal::Primes;
    use rug::ops::Pow;
    use std::str::FromStr;

    use super::*;
    use crate::arith::Plain;

    fn ecm(n: &Integer) -> Result<HashMap<Integer, usize>, Error> {
        super::ecm(
            n,
            #[cfg(feature = "progress-bar")]
            None,
        )
    }

    fn ecm_one_factor(
        n: &Integer,
        b1: usize,
        b2: usize,
        max_curve: usize,
    ) -> Result<Integer, Error> {
        super::ecm_one_factor(
            n,
            b1,
            b2,
            max_curve,
            &mut RandState::new(),
            #[cfg(feature = "progress-bar")]
            None,
        )
    }

    /// 4009823 * 99476569
    fn semiprime() -> Integer {
        Integer::from(4_009_823u64) * Integer::from(99_476_569u64)
    }

    #[test]
    fn suyama_curve_matches_sympy() {
        // Reference values computed with sympy's `_ecm_one_factor` curve setup.
        let n = Integer::from(398_883_434_337_287u64);
        let p = suyama_curve(&n, &Integer::from(123_456_789)).unwrap();
        assert_eq!(p.x, 397_114_098_224_516u64);
        assert_eq!(p.z, 208_271_263_140_048u64);
        assert_eq!(p.a24, 161_303_906_265_111u64);
    }

    #[test]
    fn square_curve_matches_gmp_ecm() {
        // `echo 398883434337287 | ecm -sigma 1:$sigma 300 0` finds a factor in step 1 exactly
        // for these sigma in [2, 40] (GMP-ECM 7.0.7).
        let n = semiprime();
        let k = stage1_multiplier(300);
        let expected = HashMap::from([
            (3, 4_009_823),
            (10, 4_009_823),
            (17, 99_476_569),
            (32, 4_009_823),
            (39, 4_009_823),
        ]);
        for sigma in 2..=40 {
            let p = square_curve(&n, &Integer::from(sigma)).unwrap();
            let g = stage1(&p, &k).z.gcd(&n);
            match expected.get(&sigma) {
                Some(&factor) => assert_eq!(g, factor, "sigma = {sigma}"),
                None => assert_eq!(g, 1, "sigma = {sigma}"),
            }
        }
    }

    #[test]
    fn batch2_curve_matches_gmp_ecm() {
        // `echo 398883434337287 | ecm -sigma 2:$sigma 300 0` finds a factor in step 1 (which
        // includes the curve setup) exactly for these sigma in [2, 40] (GMP-ECM 7.0.7).
        let n = semiprime();
        let k = stage1_multiplier(300);
        let (p, q) = (Integer::from(4_009_823), Integer::from(99_476_569));
        let mut expected: HashMap<u64, Integer> = [2, 3, 4, 5, 13, 15, 20, 21, 25, 26, 34, 36]
            .into_iter()
            .map(|sigma| (sigma, p.clone()))
            .collect();
        expected.extend([(9, q.clone()), (17, q), (24, n.clone())]);
        for sigma in 2..=40 {
            let g = match batch2_curve(&n, &Integer::from(sigma)) {
                Ok(p) => stage1(&p, &k).z.gcd(&n),
                Err(g) => g,
            };
            assert_eq!(
                g,
                expected
                    .get(&sigma)
                    .cloned()
                    .unwrap_or_default()
                    .max(Integer::from(1)),
                "sigma = {sigma}"
            );
        }
    }

    /// `p` and `q` are the same projective point: `x_p * z_q = x_q * z_p (mod n)`.
    fn same_point(p: &Point, q: &Point) -> bool {
        Integer::from(&p.x * &q.z) % &p.n == Integer::from(&q.x * &p.z) % &p.n
    }

    /// Checks the ladder of [`stage1`] against the plain arithmetic, and the group law
    /// `a*(b*P) = (a*b)*P` (`n` must be prime).
    fn check_stage1(n: &Integer, param: Param, sigma: u64) {
        let p = curve(n, param, &Integer::from(sigma)).unwrap();
        let plain = |k: &Integer| stage1_with(crate::arith::Plain::new(n), &p, k);
        assert!(same_point(&stage1(&p, &Integer::from(1)), &p));
        let ks = [2u32, 3, 7, 1000, 123_456_789].map(Integer::from);
        for k in ks.iter().chain([stage1_multiplier(200)].iter()) {
            assert!(same_point(&stage1(&p, k), &plain(k)), "{n} {param:?} {k}");
        }
        for a in &ks {
            for b in &ks {
                assert!(
                    same_point(
                        &stage1(&stage1(&p, b), a),
                        &stage1(&p, &Integer::from(a * b))
                    ),
                    "{n} {param:?} {a} {b}"
                );
            }
        }
    }

    #[test]
    fn multipliers() {
        let naive = |b1: usize| {
            Primes::all()
                .take_while(|&p| p <= b1)
                .fold(Integer::from(1), |k, p| k * p.pow(b1.ilog(p)))
        };
        for b1 in [1, 2, 3, 4, 5, 10, 100, 1000, 12_345, 100_000] {
            assert_eq!(stage1_multiplier(b1), naive(b1), "{b1}");
        }
        for (lo, hi) in [
            (1, 1),
            (2, 2),
            (3, 100),
            (100, 101),
            (1000, 100_000),
            (99, 12_345),
        ] {
            assert_eq!(
                prime_power_product(lo, hi) * naive(lo),
                naive(hi),
                "{lo} {hi}"
            );
        }
    }

    #[test]
    fn trial_division_bound() {
        let (factors, cofactor) = trial_division(&Integer::from(2u64 * 2 * 3 * 65_521 * 65_537));
        assert_eq!(
            factors,
            HashMap::from([
                (Integer::from(2), 2),
                (Integer::from(3), 1),
                (Integer::from(65_521), 1)
            ])
        );
        assert_eq!(cofactor, 65_537);
        let (factors, cofactor) = trial_division(&Integer::from(1));
        assert!(factors.is_empty());
        assert_eq!(cofactor, 1);
    }

    #[test]
    fn stage1_matches_reference() {
        let mut rand = RandState::new();
        // One limb, every limb count of `Mont`, and the `Plain` fallback.
        for bits in [
            40, 64, 100, 128, 192, 256, 320, 384, 448, 512, 700, 1024, 1025, 1500,
        ] {
            let mut n = Integer::from(Integer::random_bits(bits, &mut rand));
            n.set_bit(bits - 1, true);
            // Prime: building the curves never fails.
            let n = n.next_prime();
            check_stage1(&n, Param::Square, 1_234_567);
            check_stage1(&n, Param::Suyama, 1_234_567);
            check_stage1(&n, Param::Batch2, 1_234_567);
        }
    }

    #[test]
    fn stage2_matches_plain() {
        // Montgomery and plain arithmetic give the same product, so the same gcd.
        let n = semiprime();
        let (b1, b2) = (100, 10_000);
        let k = stage1_multiplier(b1);
        for plan in [Stage2Plan::pairs(b1, b2), Stage2Plan::poly(&n, b1, b2)] {
            for sigma in 2..100 {
                let q = stage1(&square_curve(&n, &Integer::from(sigma)).unwrap(), &k);
                let g = stage2(&q, &plan);
                assert_eq!(g, stage2_with(Plain::new(&n), &q, &plan));
            }
        }
    }

    #[test]
    fn stage2_finds_factor() {
        // With sigma = 9, stage 1 misses the factor and stage 2 finds it.
        let k = stage1_multiplier(2000);
        assert_eq!(
            run_curve(
                &semiprime(),
                Param::Suyama,
                &Integer::from(9),
                &k,
                &Stage2Plan::new(&semiprime(), 2000, 147_396)
            ),
            CurveOutcome::Stage2(Integer::from(4_009_823))
        );
    }

    /// Stage 2 must find a factor whenever l*Q = O modulo a factor of n for some prime
    /// b1 < l <= b2, unless it finds all factors at once (g = n).
    ///
    /// Returns the number of curves where stage 2 had a factor to find.
    fn check_stage2_primes(
        n: &Integer,
        param: Param,
        sigmas: std::ops::Range<u64>,
        b1: usize,
        b2: usize,
    ) -> usize {
        let k = stage1_multiplier(b1);
        let mut plans = vec![Stage2Plan::pairs(b1, b2), Stage2Plan::poly(n, b1, b2)];
        // And polynomial plans with d2 > 1: the default d2 with F padded to fill two whole
        // blocks, and the largest prime d2 <= b1 with partial blocks.
        let d1 = [2310, 210, 30, 6]
            .into_iter()
            .find(|&d1| crate::stage2::prime_factors(d1).iter().all(|&p| p <= b1))
            .unwrap();
        let d2 = crate::stage2_poly::default_d2(b1, d1);
        let other = [29, 23, 19, 17, 13, 11, 7, 5]
            .into_iter()
            .find(|&p| p <= b1 && !d1.is_multiple_of(p));
        for (d2, blocks) in [(d2, 2)].into_iter().chain(other.map(|p| (p, 0))) {
            let plan = crate::stage2_poly::PolyPlan::new(b1, b2, d1, d2, blocks);
            assert!(plan.b2_covered() >= b2);
            plans.push(Stage2Plan::Poly(plan));
        }
        let primes: Vec<Integer> = Primes::all()
            .skip_while(|&l| l <= b1)
            .take_while(|&l| l <= b2)
            .map(Integer::from)
            .collect();
        let mut checked = 0;
        for sigma in sigmas {
            let Ok(p) = curve(n, param, &Integer::from(sigma)) else {
                continue;
            };
            let q = stage1(&p, &k);
            if q.z.clone().gcd(n) != 1 {
                continue;
            }
            let expected = primes.iter().any(|l| {
                let g = stage1(&q, l).z.gcd(n);
                g != 1 && &g != n
            });
            for plan in &plans {
                let g = stage2(&q, plan);
                assert!(n.is_divisible(&g));
                if expected {
                    assert_ne!(
                        g, 1,
                        "stage 2 missed a factor: {plan:?} {param:?} {sigma} {b1} {b2} {n}"
                    );
                }
            }
            checked += expected as usize;
        }
        checked
    }

    #[test]
    fn stage2_checks_all_primes() {
        let n = semiprime();
        for (param, first) in [(Param::Suyama, 6), (Param::Square, 2), (Param::Batch2, 2)] {
            assert!(check_stage2_primes(&n, param, first..first + 300, 100, 10_000) > 100);
        }
    }

    #[test]
    fn stage2_checks_all_primes_bounds() {
        // Small bounds, b2 not a multiple of D, b2 <= b1, and D larger than b2.
        let n = semiprime();
        let mut checked = 0;
        for (b1, b2) in [
            (6, 4),
            (6, 6),
            (6, 8),
            (6, 100),
            (8, 50),
            (10, 1000),
            (12, 3001),
            (14, 20_000),
            (100, 10_001),
            (1000, 999),
            (2000, 147_397),
        ] {
            for param in [Param::Square, Param::Batch2] {
                checked += check_stage2_primes(&n, param, 2..80, b1, b2);
            }
        }
        assert!(checked > 100);
    }

    #[test]
    fn stage2_checks_all_primes_limbs() {
        // A small factor times a large prime: every limb count of `Mont`, and `Plain`.
        let mut rand = RandState::new();
        for bits in [64, 100, 192, 320, 512, 700, 1000, 1100] {
            let mut q = Integer::from(Integer::random_bits(bits, &mut rand));
            q.set_bit(bits - 1, true);
            let n = Integer::from(4_009_823) * q.next_prime();
            let checked = check_stage2_primes(&n, Param::Batch2, 2..60, 100, 5000)
                + check_stage2_primes(&n, Param::Square, 2..30, 30, 3000);
            assert!(checked > 5, "{bits} bits: {checked}");
        }
    }

    #[test]
    fn tiny_moduli() {
        // Degenerate curves and non-invertible setups are frequent modulo tiny numbers: a curve
        // either fails or returns a proper factor, it never panics.
        let k = stage1_multiplier(100);
        let plans = [
            Stage2Plan::pairs(100, 1000),
            Stage2Plan::poly(&Integer::from(1u64 << 40), 100, 1000),
        ];
        for n in (9u32..1500).step_by(2) {
            let n = Integer::from(n);
            if n.is_probably_prime(PRIMALITY_REPS) != IsPrime::No {
                continue;
            }
            for (param, first) in [(Param::Suyama, 6), (Param::Square, 2), (Param::Batch2, 2)] {
                for (sigma, plan) in
                    (first..first + 20).flat_map(|s| plans.iter().map(move |p| (s, p)))
                {
                    let sigma = Integer::from(sigma);
                    match run_curve(&n, param, &sigma, &k, plan) {
                        CurveOutcome::Setup(g)
                        | CurveOutcome::Stage1(g)
                        | CurveOutcome::Stage2(g) => {
                            assert!(g != 1 && g != n && n.is_divisible(&g), "{n} {param:?}")
                        }
                        CurveOutcome::Failed => {}
                    }
                }
            }
        }
    }

    #[test]
    fn small_b1() {
        // b1 < 2*sqrt(b2) used to underflow when computing the first giant step.
        let n = semiprime();
        match ecm_one_factor(&n, 100, 100_000, 50) {
            Ok(g) => assert!(g != 1 && g != n && n.is_divisible(&g)),
            Err(e) => assert!(matches!(e, Error::ECMFailed)),
        }
    }

    fn ecm_with_params(
        n: &Integer,
        b1: usize,
        b2: usize,
        max_curve: usize,
    ) -> Result<HashMap<Integer, usize>, Error> {
        super::ecm_with_params(
            n,
            b1,
            b2,
            max_curve,
            1234,
            #[cfg(feature = "progress-bar")]
            None,
        )
    }

    #[test]
    fn failure_is_an_error() {
        // Bounds far too small for a 8-digit factor: no wrong factorization, an error.
        assert!(matches!(
            ecm_with_params(&semiprime(), 6, 100, 2),
            Err(Error::ECMFailed)
        ));
    }

    #[test]
    fn driver_checks_bounds() {
        assert!(matches!(
            ecm_with_params(&semiprime(), 3, 100, 10),
            Err(Error::BoundsNotEven)
        ));
        assert!(matches!(
            ecm_with_params(&semiprime(), 4, 100, 10),
            Err(Error::BoundsTooSmall)
        ));
    }

    #[test]
    fn small_numbers() {
        assert!(ecm(&Integer::from(1)).unwrap().is_empty());
        for n in 2u32..3000 {
            let factors = ecm(&Integer::from(n)).unwrap();
            let mut product = Integer::from(1);
            for (p, &e) in &factors {
                assert_ne!(p.is_probably_prime(30), IsPrime::No, "{n}: {p}");
                product *= p.clone().pow(e as u32);
            }
            assert_eq!(product, n);
        }
    }

    #[test]
    #[should_panic(expected = "only positive numbers")]
    fn zero() {
        let _ = ecm(&Integer::ZERO);
    }

    #[test]
    fn perfect_power_of_prime() {
        let p = Integer::from(2_802_377);
        let n = p.clone().pow(3);
        assert_eq!(
            ecm_with_params(&n, 2000, 147_396, 10).unwrap(),
            HashMap::from([(p, 3)])
        );
    }

    #[test]
    fn perfect_power_of_composite() {
        // (4009823 * 99476569)^2: the root is composite and still has to be split.
        let n = semiprime().pow(2);
        assert_eq!(
            ecm_with_params(&n, 2000, 147_396, 100).unwrap(),
            HashMap::from([
                (Integer::from(4_009_823), 2),
                (Integer::from(99_476_569), 2)
            ])
        );
    }

    #[test]
    fn too_small_bounds() {
        assert!(matches!(
            ecm_one_factor(&semiprime(), 4, 100, 10),
            Err(Error::BoundsTooSmall)
        ));
    }

    #[test]
    fn no_curve() {
        assert!(matches!(
            ecm_one_factor(&semiprime(), 2000, 147_396, 0),
            Err(Error::ECMFailed)
        ));
    }

    /// Factorizations from sympy's ECM tests.
    #[test]
    fn sympy() {
        let cases: [(&str, &[(&str, usize)]); 9] = [
            ("398883434337287", &[("99476569", 1), ("4009823", 1)]),
            (
                "46167045131415113",
                &[("43", 1), ("2634823", 1), ("407485517", 1)],
            ),
            (
                "64211816600515193",
                &[("281719", 1), ("359641", 1), ("633767", 1)],
            ),
            (
                "168541512131094651323",
                &[("79", 1), ("113", 1), ("11011069", 1), ("1714635721", 1)],
            ),
            (
                "631211032315670776841",
                &[("9312934919", 1), ("67777885039", 1)],
            ),
            (
                "4132846513818654136451",
                &[("47", 1), ("160343", 1), ("2802377", 1), ("195692803", 1)],
            ),
            (
                "4516511326451341281684513",
                &[("3", 2), ("39869", 1), ("131743543", 1), ("95542348571", 1)],
            ),
            (
                "3146531246531241245132451321",
                &[("3", 1), ("100327907731", 1), ("10454157497791297", 1)],
            ),
            (
                "4269021180054189416198169786894227",
                &[
                    ("184039", 1),
                    ("241603", 1),
                    ("333331", 1),
                    ("477973", 1),
                    ("618619", 1),
                    ("974123", 1),
                ],
            ),
        ];
        for (n, factors) in cases {
            let expected = factors
                .iter()
                .map(|&(p, e)| (Integer::from_str(p).unwrap(), e))
                .collect();
            assert_eq!(
                ecm(&Integer::from_str(n).unwrap()).unwrap(),
                expected,
                "{n}"
            );
        }
    }

    #[test]
    fn same_factors() {
        assert_eq!(
            ecm(&Integer::from_str("7853316850129").unwrap()).unwrap(),
            HashMap::from([(Integer::from_str("2802377").unwrap(), 2)])
        );
    }

    #[test]
    fn small_prime() {
        assert_eq!(
            ecm(&Integer::from(17)).unwrap(),
            HashMap::from([(Integer::from(17), 1)])
        );
    }

    #[test]
    fn big_prime() {
        assert_eq!(
            ecm(&Integer::from_str("21472883178031195225853317139").unwrap()).unwrap(),
            HashMap::from([(
                Integer::from_str("21472883178031195225853317139").unwrap(),
                1
            )])
        );
    }
}
