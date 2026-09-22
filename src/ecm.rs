use crate::{
    arith::{with_arith, Arith, Factor},
    curve::Curve,
    point::Point,
};
#[cfg(feature = "progress-bar")]
use indicatif::ProgressBar;
use primal::Primes;
use rug::{integer::IsPrime, rand::RandState, Integer};
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
/// with the curves of GMP-ECM's parametrization 2 (see [`Param::Batch2`]). Here Montgomery
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
/// hence k*P = l*|E(FF(q))|*P. Thus kP.z_cord = 0 (mod q), and the unknown factor of N (q)
/// can be recovered by taking gcd(kP.z_cord, N).
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

    // Stage 2 needs at least 2 baby steps: `stage2_d(b1, b2) >= 2`.
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
    let param = Param::default();

    for _ in 0..max_curve {
        #[cfg(feature = "progress-bar")]
        if let Some(pb) = pb {
            pb.inc(1);
        }

        let sigma = random_sigma(n, param, rgen);
        match run_curve(n, param, &sigma, &k, b1, b2) {
            CurveOutcome::Setup(g) | CurveOutcome::Stage1(g) | CurveOutcome::Stage2(g) => {
                return Ok(g)
            }
            CurveOutcome::Failed => {}
        }
    }

    // ECM failed, Increase the bounds
    Err(Error::ECMFailed)
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
/// `k` must be the stage 1 multiplier returned by [`stage1_multiplier`] for `b1`.
pub fn run_curve(
    n: &Integer,
    param: Param,
    sigma: &Integer,
    k: &Integer,
    b1: usize,
    b2: usize,
) -> CurveOutcome {
    let q = match curve(n, param, sigma) {
        Ok(q) => q,
        // If g = 1 or n, try another curve
        Err(g) if g == 1 || &g == n => return CurveOutcome::Failed,
        Err(g) => return CurveOutcome::Setup(g),
    };

    let q = stage1(&q, k);
    let g = q.z_cord.clone().gcd(n);

    // Stage 1 factor
    if &g != n && g != 1 {
        return CurveOutcome::Stage1(g);
    }

    // Stage 1 failure. Q.z = 0, Try another curve
    if &g == n {
        return CurveOutcome::Failed;
    }

    let g = stage2(&q, b1, b2);

    // Stage 2 Factor found
    if &g != n && g != 1 {
        return CurveOutcome::Stage2(g);
    }

    CurveOutcome::Failed
}

/// Stage 1 multiplier: product of the largest powers of all primes `p <= b1` that are `<= b1`.
pub fn stage1_multiplier(b1: usize) -> Integer {
    let mut k = Integer::from(1);
    for p in Primes::all().take_while(|&p| p <= b1) {
        k *= p.pow(b1.ilog(p));
    }
    k
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
    Ok(Point::new(u_3, v_3, a24, n.clone()))
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
    Ok(Point::new(2.into(), 1.into(), d, n.clone()))
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
    Ok(Point::new(2.into(), 1.into(), a24, n.clone()))
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
    with_arith!(&p.modulus, |arith| stage1_with(arith, p, k))
}

fn stage1_with<A: Arith>(arith: A, p: &Point, k: &Integer) -> Point {
    let n: &Integer = &p.modulus;
    // Normalizing P to z = 1 saves a multiplication per ladder step. If z is not invertible, P
    // is the point at infinity modulo a factor of n, and so is k*P: P has the same gcd.
    let Ok(z_inv) = p.z_cord.clone().invert(n) else {
        return p.clone();
    };
    let x = Integer::from(&p.x_cord * &z_inv) % n;

    let curve = Curve::new(arith, &p.a_24);
    let q = curve.ladder(&curve.arith.factor(&x), &Factor::One, k);
    p.on_same_curve(curve.arith.to_integer(&q.x), curve.arith.to_integer(&q.z))
}

/// Number of baby steps `D` of stage 2.
///
/// `D <= b1 / 2 - 1` keeps `b1 - 2*D` positive: it is the multiplier of the first giant step.
fn stage2_d(b1: usize, b2: usize) -> usize {
    b2.isqrt().min(b1 / 2 - 1)
}

/// Stage 2 - Improved Standard Continuation.
///
/// Returns `gcd(g, n)` where `g` is the accumulated product over the primes in `(b1, b2]`.
/// Requires `b1 >= 6` and `b2 >= 4`, so that `D >= 2`.
///
/// The x-coordinates of a point and its inverse are equal, so the primes `r - (2i + 1)` and
/// `r + (2i + 1)` are both checked by comparing `r*Q` with `S[i] = (2i + 1)*Q`: each giant step
/// covers `4*D` instead of `2*D`.
pub fn stage2(q: &Point, b1: usize, b2: usize) -> Integer {
    with_arith!(&q.modulus, |arith| stage2_with(arith, q, b1, b2))
}

fn stage2_with<A: Arith>(arith: A, q: &Point, b1: usize, b2: usize) -> Integer {
    let d = stage2_d(b1, b2);
    let two_d = 2 * d;
    let curve = Curve::new(arith, &q.a_24);
    let a = &curve.arith;
    let mut scratch = curve.scratch();

    // S[i] = (2*i + 1)*Q
    let q_xz = curve.point(&q.x_cord, &q.z_cord);
    let mut q2 = curve.infinity();
    curve.double(&mut q2, &q_xz, &mut scratch);
    let mut s = Vec::with_capacity(d);
    s.push(q_xz.clone());
    let mut q3 = curve.infinity();
    curve.add(&mut q3, &q2, &q_xz, &q_xz, &mut scratch);
    s.push(q3);
    for i in 2..d {
        let mut next = curve.infinity();
        curve.add(&mut next, &s[i - 1], &q2, &s[i - 2], &mut scratch);
        s.push(next);
    }
    let beta: Vec<A::Elem> = s
        .iter()
        .map(|s| {
            let mut b = a.zero();
            a.mul(&mut b, &s.x, &s.z);
            b
        })
        .collect();

    let mut g = a.residue(&Integer::from(1));
    let (xq, zq) = (a.factor(&q.x_cord), a.factor(&q.z_cord));
    let w = curve.ladder(&xq, &zq, &Integer::from(2 * two_d));
    let mut t = curve.ladder(&xq, &zq, &Integer::from(b1 - two_d));
    let mut r = curve.ladder(&xq, &zq, &Integer::from(b1 + two_d));
    let mut next = curve.infinity();

    // Each prime is checked once, even when `rr - delta` and `rr + delta` are both prime.
    let mut seen = vec![false; d];
    let mut deltas: Vec<usize> = Vec::with_capacity(two_d);
    let mut primes = Primes::all().skip_while(|&p| p <= b1).peekable();
    // Reused by every prime, to avoid allocating in the inner loop.
    let (mut f, mut diff, mut sum, mut alpha, mut acc) =
        (a.zero(), a.zero(), a.zero(), a.zero(), a.zero());

    for rr in (b1 + two_d..b2 + two_d).step_by(2 * two_d) {
        // R = rr*Q, and the primes of this giant step are rr +/- (2*delta + 1)
        deltas.clear();
        while let Some(p) = primes.next_if(|&p| p < rr + two_d) {
            let delta = p.abs_diff(rr) >> 1;
            if !seen[delta] {
                seen[delta] = true;
                deltas.push(delta);
            }
        }

        a.mul(&mut alpha, &r.x, &r.z);
        for &delta in &deltas {
            seen[delta] = false;
            // We want to calculate
            // f = R.x_cord * S[delta].z_cord - S[delta].x_cord * R.z_cord
            //   = (R.x - S.x) * (R.z + S.z) - alpha + beta[delta]
            a.sub(&mut diff, &r.x, &s[delta].x);
            a.add(&mut sum, &r.z, &s[delta].z);
            a.mul(&mut f, &diff, &sum);
            a.sub(&mut diff, &f, &alpha);
            a.add(&mut f, &diff, &beta[delta]);
            a.mul(&mut acc, &g, &f);
            std::mem::swap(&mut g, &mut acc);
        }

        // T, R = R, R + W: R + W is computed from the old R, with difference T = R - W
        curve.add(&mut next, &r, &w, &t, &mut scratch);
        std::mem::swap(&mut t, &mut r);
        std::mem::swap(&mut r, &mut next);
    }
    a.gcd(&g)
}

/// Removes the factors of `n` among the first 100 000 primes.
///
/// Returns the found factors with their multiplicity, and the remaining cofactor.
pub fn trial_division(n: &Integer) -> (HashMap<Integer, usize>, Integer) {
    let mut factors = HashMap::new();
    let mut n: Integer = n.clone();
    for prime in Primes::all().take(100_000) {
        if n.is_divisible_u(prime as u32) {
            let prime = Integer::from(prime);
            while n.is_divisible(&prime) {
                n /= &prime;
                *factors.entry(prime.clone()).or_insert(0) += 1;
            }
        }
    }
    (factors, n)
}

/// Default `(B1, B2, max_curve)` for a number of `digits` decimal digits.
///
/// Optimal params retrieved from <https://gitlab.inria.fr/zimmerma/ecm>
pub fn optimal_params(digits: usize) -> (usize, usize, usize) {
    match digits {
        1..=10 => (2_000, 160_000, 35),
        11..=15 => (5_000, 500_000, 500),
        16..=20 => (11_000, 1_900_000, 74),
        21..=25 => (50_000, 13_000_000, 214),
        26..=30 => (250_000, 130_000_000, 430),
        31..=35 => (1_000_000, 1_000_000_000, 904),
        36..=40 => (3_000_000, 5_700_000_000, 2350),
        41..=45 => (11_000_000, 35_000_000_000, 4480),
        46..=50 => (44_000_000, 240_000_000_000, 7553),
        51..=55 => (110_000_000, 780_000_000_000, 17769),
        56..=60 => (260_000_000, 3_200_000_000_000, 42017),
        _ => (850_000_000, 16_000_000_000_000, 69408),
    }
}

/// Performs factorization using Lenstra's Elliptic curve method.
///
/// This function repeatedly calls `ecm_one_factor` to compute the factors
/// of n. First all the small factors are taken out using trial division.
/// Then `ecm_one_factor` is used to compute one factor at a time.
///
/// # Parameters
///
/// - `n`: Number to be factored.
pub fn ecm(
    n: &Integer,
    #[cfg(feature = "progress-bar")] pb: Option<&ProgressBar>,
) -> Result<HashMap<Integer, usize>, Error> {
    let optimal_params = optimal_params(n.to_string().len());

    ecm_with_params(
        n,
        optimal_params.0,
        optimal_params.1,
        optimal_params.2,
        1234,
        #[cfg(feature = "progress-bar")]
        pb,
    )
}

/// Performs factorization using Lenstra's Elliptic curve method.
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

    let (mut factors, n) = trial_division(n);

    let mut rand_state = RandState::new();
    rand_state.seed(&seed.into());

    // Composite factors left to split, with the multiplicity they have in the original number.
    let mut queue = Vec::new();
    sort_factor(n, 1, &mut factors, &mut queue);

    while let Some((n, exponent)) = queue.pop() {
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
        sort_factor(factor, exponent * multiplicity, &mut factors, &mut queue);
        sort_factor(cofactor, exponent, &mut factors, &mut queue);
    }

    Ok(factors)
}

/// Records `n^exponent`: a prime goes to `factors`, a perfect power is reduced to its root, and
/// any other composite is queued for [`ecm_one_factor`].
fn sort_factor(
    n: Integer,
    exponent: usize,
    factors: &mut HashMap<Integer, usize>,
    queue: &mut Vec<(Integer, usize)>,
) {
    if n == 1 {
        return;
    }
    if n.is_probably_prime(PRIMALITY_REPS) != IsPrime::No {
        *factors.entry(n).or_insert(0) += exponent;
        return;
    }
    match perfect_power(&n) {
        Some((root, power)) => sort_factor(root, exponent * power as usize, factors, queue),
        None => queue.push((n, exponent)),
    }
}

/// Writes `n` as `root^power` with the smallest possible `root`, if it is a perfect power.
fn perfect_power(n: &Integer) -> Option<(Integer, u32)> {
    if !n.is_perfect_power() {
        return None;
    }
    for power in Primes::all().take_while(|&p| p as u32 <= n.significant_bits()) {
        let (root, remainder) = n.clone().root_rem(Integer::new(), power as u32);
        if remainder == 0 {
            return Some((root, power as u32));
        }
    }
    None
}

#[cfg(test)]
mod tests {
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
        assert_eq!(p.x_cord, 397_114_098_224_516u64);
        assert_eq!(p.z_cord, 208_271_263_140_048u64);
        assert_eq!(*p.a_24, 161_303_906_265_111u64);
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
            let g = stage1(&p, &k).z_cord.gcd(&n);
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
                Ok(p) => stage1(&p, &k).z_cord.gcd(&n),
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

    /// Checks the Montgomery ladder of [`stage1`] against the reference [`Point::mont_ladder`].
    fn check_stage1(n: &Integer, param: Param, sigma: u64) {
        let p = curve(n, param, &Integer::from(sigma)).unwrap();
        for k in [1u32, 2, 3, 7, 1000, 123_456_789] {
            let k = Integer::from(k);
            assert_eq!(
                stage1(&p, &k),
                p.mont_ladder(&k),
                "{n} {param:?} {sigma} {k}"
            );
        }
        let k = stage1_multiplier(200);
        assert_eq!(stage1(&p, &k), p.mont_ladder(&k));
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
        for sigma in 2..100 {
            let q = stage1(&square_curve(&n, &Integer::from(sigma)).unwrap(), &k);
            let g = stage2(&q, b1, b2);
            assert_eq!(g, stage2_with(Plain::new(&n), &q, b1, b2));
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
                2000,
                147_396
            ),
            CurveOutcome::Stage2(Integer::from(4_009_823))
        );
    }

    /// Stage 2 must find a factor whenever l*Q = O modulo a factor of n for some prime
    /// b1 < l <= b2, unless it finds all factors at once (g = n).
    fn check_stage2_primes(param: Param, sigmas: std::ops::Range<u64>) {
        let n = semiprime();
        let (b1, b2) = (100, 10_000);
        let k = stage1_multiplier(b1);
        let mut checked = 0;
        for sigma in sigmas {
            let q = stage1(&curve(&n, param, &Integer::from(sigma)).unwrap(), &k);
            if q.z_cord.clone().gcd(&n) != 1 {
                continue;
            }
            let expected = Primes::all()
                .skip_while(|&l| l <= b1)
                .take_while(|&l| l <= b2)
                .any(|l| {
                    let g = q.mont_ladder(&Integer::from(l)).z_cord.gcd(&n);
                    g != 1 && g != n
                });
            if expected {
                let g = stage2(&q, b1, b2);
                assert_ne!(g, 1, "stage 2 missed a factor with sigma = {sigma}");
                checked += 1;
            }
        }
        assert!(checked > 100);
    }

    #[test]
    fn stage2_checks_all_primes() {
        check_stage2_primes(Param::Suyama, 6..300);
        check_stage2_primes(Param::Square, 2..300);
        check_stage2_primes(Param::Batch2, 2..300);
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

    #[test]
    fn sympy_1() {
        assert_eq!(
            ecm(&Integer::from_str("398883434337287").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("99476569").unwrap(), 1),
                (Integer::from_str("4009823").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_2() {
        assert_eq!(
            ecm(&Integer::from_str("46167045131415113").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("43").unwrap(), 1),
                (Integer::from_str("2634823").unwrap(), 1),
                (Integer::from_str("407485517").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_3() {
        assert_eq!(
            ecm(&Integer::from_str("64211816600515193").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("281719").unwrap(), 1),
                (Integer::from_str("359641").unwrap(), 1),
                (Integer::from_str("633767").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_4() {
        assert_eq!(
            ecm(&Integer::from_str("168541512131094651323").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("79").unwrap(), 1),
                (Integer::from_str("113").unwrap(), 1),
                (Integer::from_str("11011069").unwrap(), 1),
                (Integer::from_str("1714635721").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_5() {
        assert_eq!(
            ecm(&Integer::from_str("631211032315670776841").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("9312934919").unwrap(), 1),
                (Integer::from_str("67777885039").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_6() {
        assert_eq!(
            ecm(&Integer::from_str("4132846513818654136451").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("47").unwrap(), 1),
                (Integer::from_str("160343").unwrap(), 1),
                (Integer::from_str("2802377").unwrap(), 1),
                (Integer::from_str("195692803").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_7() {
        assert_eq!(
            ecm(&Integer::from_str("4516511326451341281684513").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("3").unwrap(), 2),
                (Integer::from_str("39869").unwrap(), 1),
                (Integer::from_str("131743543").unwrap(), 1),
                (Integer::from_str("95542348571").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_8() {
        assert_eq!(
            ecm(&Integer::from_str("3146531246531241245132451321").unwrap(),).unwrap(),
            HashMap::from([
                (Integer::from_str("3").unwrap(), 1),
                (Integer::from_str("100327907731").unwrap(), 1),
                (Integer::from_str("10454157497791297").unwrap(), 1),
            ])
        );
    }

    #[test]
    fn sympy_9() {
        assert_eq!(
            ecm(&Integer::from_str("4269021180054189416198169786894227").unwrap()).unwrap(),
            HashMap::from([
                (Integer::from_str("184039").unwrap(), 1),
                (Integer::from_str("241603").unwrap(), 1),
                (Integer::from_str("333331").unwrap(), 1),
                (Integer::from_str("477973").unwrap(), 1),
                (Integer::from_str("618619").unwrap(), 1),
                (Integer::from_str("974123").unwrap(), 1),
            ])
        );
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
