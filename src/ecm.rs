use crate::point::Point;
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

/// Returns one factor of n using Lenstra's 2 Stage Elliptic curve Factorization
/// with Suyama's Parameterization. Here Montgomery arithmetic is used for fast
/// computation of addition and doubling of points in elliptic curve.
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

    if n.is_probably_prime(1000) != IsPrime::No {
        return Err(Error::NumberIsPrime);
    }

    #[cfg(feature = "progress-bar")]
    if let Some(pb) = pb {
        pb.set_length(max_curve as u64);
        pb.set_position(0);
    }

    let k = stage1_multiplier(b1);
    let sigma_range = Integer::from(n - 6);

    for _ in 0..max_curve {
        #[cfg(feature = "progress-bar")]
        if let Some(pb) = pb {
            pb.inc(1);
        }

        // Suyama's parametrization: sigma in [6, n - 1]
        let sigma = sigma_range.clone().random_below(rgen) + 6;
        match run_curve(n, &sigma, &k, b1, b2) {
            CurveOutcome::Setup(g) | CurveOutcome::Stage1(g) | CurveOutcome::Stage2(g) => {
                return Ok(g)
            }
            CurveOutcome::Failed => {}
        }
    }

    // ECM failed, Increase the bounds
    Err(Error::ECMFailed)
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

/// Runs one ECM curve (curve setup, stage 1 and stage 2) with the given `sigma`.
///
/// `k` must be the stage 1 multiplier returned by [`stage1_multiplier`] for `b1`.
pub fn run_curve(n: &Integer, sigma: &Integer, k: &Integer, b1: usize, b2: usize) -> CurveOutcome {
    let q = match suyama_curve(n, sigma) {
        Ok(q) => q,
        // If g = n, try another curve
        Err(g) if &g == n => return CurveOutcome::Failed,
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

/// Stage 1: computes `k*P`.
pub fn stage1(p: &Point, k: &Integer) -> Point {
    p.mont_ladder(k)
}

/// Number of baby steps `D` of stage 2.
///
/// `D <= b1 / 2 - 1` keeps `b1 - 1 - 2*D` positive: it is the multiplier of the first giant step.
fn stage2_d(b1: usize, b2: usize) -> usize {
    b2.isqrt().min(b1 / 2 - 1)
}

/// Stage 2 - Improved Standard Continuation.
///
/// Returns `gcd(g, n)` where `g` is the accumulated product over the primes in `(b1, b2]`.
/// Requires `b1 >= 6` and `b2 >= 4`, so that `D >= 2`.
pub fn stage2(q: &Point, b1: usize, b2: usize) -> Integer {
    let n = &q.modulus;
    let d = stage2_d(b1, b2);
    let two_d = 2 * d;
    let mut beta: Vec<Integer> = vec![Integer::default(); d + 1];
    let mut s: Vec<Point> = vec![Point::default(); d + 1];

    s[1] = q.double();
    s[2] = s[1].double();
    beta[1] = Integer::from(&s[1].x_cord * &s[1].z_cord) % n;
    beta[2] = Integer::from(&s[2].x_cord * &s[2].z_cord) % n;

    for d in 3..=(d) {
        s[d] = s[d - 1].add(&s[1], &s[d - 2]);
        beta[d] = Integer::from(&s[d].x_cord * &s[d].z_cord) % n;
    }

    let mut g = Integer::from(1);
    let b = b1 - 1;
    let mut t = q.mont_ladder(&Integer::from(b - two_d));
    let mut r = q.mont_ladder(&Integer::from(b));

    // R = rr*Q: primes q in [rr + 2, rr + 2*D] are checked with S[delta] = 2*delta*Q,
    // where q = rr + 2*delta.
    let mut primes = Primes::all().skip_while(|&q| q < b + 2).peekable();
    for rr in (b..b2).step_by(two_d) {
        let alpha = Integer::from(&r.x_cord * &r.z_cord) % n;
        while let Some(q) = primes.next_if(|&q| q <= rr + two_d) {
            let delta = (q - rr) / 2;
            // We want to calculate
            // f = R.x_cord * S[delta].z_cord - S[delta].x_cord * R.z_cord
            let f = Integer::from(&r.x_cord - &s[delta].x_cord)
                * Integer::from(&r.z_cord + &s[delta].z_cord)
                - &alpha
                + &beta[delta];
            g = (g * f) % n;
        }
        // T, R = R, R + S[D]: R + S[D] is computed from the old R, with difference T = R - S[D]
        let next = r.add(&s[d], &t);
        t = std::mem::replace(&mut r, next);
    }
    g.gcd(n)
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
    let (mut factors, mut n) = trial_division(n);

    let mut rand_state = RandState::new();
    rand_state.seed(&seed.into());

    while n != 1 {
        let factor = ecm_one_factor(
            &n,
            b1,
            b2,
            max_curve,
            &mut rand_state,
            #[cfg(feature = "progress-bar")]
            pb,
        )
        .unwrap_or(n.clone());

        while n.is_divisible(&factor) {
            n /= &factor;
            *factors.entry(factor.clone()).or_insert(0) += 1;
        }
    }

    Ok(factors)
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

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
        assert_eq!(p.a_24, 161_303_906_265_111u64);
    }

    #[test]
    fn stage2_finds_factor() {
        // With sigma = 9, stage 1 misses the factor and stage 2 finds it.
        let k = stage1_multiplier(2000);
        assert_eq!(
            run_curve(&semiprime(), &Integer::from(9), &k, 2000, 147_396),
            CurveOutcome::Stage2(Integer::from(4_009_823))
        );
    }

    #[test]
    fn stage2_checks_all_primes() {
        // Stage 2 must find a factor whenever l*Q = O modulo a factor of n for some prime
        // b1 < l <= b2, unless it finds all factors at once (g = n).
        let n = semiprime();
        let (b1, b2) = (100, 10_000);
        let k = stage1_multiplier(b1);
        let mut checked = 0;
        for sigma in 6..300 {
            let q = stage1(&suyama_curve(&n, &Integer::from(sigma)).unwrap(), &k);
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
    fn small_b1() {
        // b1 < 2*sqrt(b2) used to underflow when computing the first giant step.
        let n = semiprime();
        match ecm_one_factor(&n, 100, 100_000, 50) {
            Ok(g) => assert!(g != 1 && g != n && n.is_divisible(&g)),
            Err(e) => assert!(matches!(e, Error::ECMFailed)),
        }
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
