//! Driver of [`crate::ecm()`]: finds the factors from the smallest to the largest, like
//! GMP-ECM's recommended usage.
//!
//! After trial division, each composite goes through levels of increasing factor size (GMP-ECM's
//! table of optimal `B1` for factors of 15, 20, ..., 65 digits, with its default `B2`). A level
//! runs its expected number of curves (from GMP-ECM's probability model, for the `B2` stage 2
//! really covers), then the next level starts: the factors are found about in the time it
//! takes to find a factor of their size, whatever the size of the number. The levels stop at
//! half the digits of the number, as its smallest factor is below its square root. Before the
//! curves of a level, P-1 (see [`crate::pm1`]) is extended to a bound proportional to the `B1`
//! of the level: it costs about as much as a few curves, and finds the factors `p` with a
//! smooth `p - 1` of much larger sizes.
//!
//! When a factor is found, the cofactor keeps the progress made: its factors are larger than
//! the ones already searched for.

use crate::{
    cost::Costs,
    ecm::{
        CurveOutcome, Error, Param, rand_state, random_sigma, run_curve, sort_factor,
        stage1_multiplier, trial_division,
    },
    pm1::Pm1,
    rho::ecm_prob,
    stage2::Stage2Plan,
};
#[cfg(feature = "progress-bar")]
use indicatif::ProgressBar;
use rug::{Integer, rand::RandState};
use std::{collections::HashMap, rc::Rc};

/// A level of the search: factors of up to `digits` digits, with the bounds of the curves.
struct Level {
    digits: u32,
    b1: usize,
    b2: usize,
}

/// Levels of the search: GMP-ECM's optimal `B1` for each factor size (README, table 1) with its
/// default `B2` (`ecm -v`), below a level for factors up to 10 digits.
const LEVELS: [Level; 12] = [
    Level {
        digits: 10,
        b1: 300,
        b2: 9_846,
    },
    Level {
        digits: 15,
        b1: 2_000,
        b2: 147_396,
    },
    Level {
        digits: 20,
        b1: 11_000,
        b2: 1_873_422,
    },
    Level {
        digits: 25,
        b1: 50_000,
        b2: 12_746_592,
    },
    Level {
        digits: 30,
        b1: 250_000,
        b2: 128_992_510,
    },
    Level {
        digits: 35,
        b1: 1_000_000,
        b2: 1_045_563_762,
    },
    Level {
        digits: 40,
        b1: 3_000_000,
        b2: 5_706_890_290,
    },
    Level {
        digits: 45,
        b1: 11_000_000,
        b2: 35_133_391_030,
    },
    Level {
        digits: 50,
        b1: 43_000_000,
        b2: 240_490_660_426,
    },
    Level {
        digits: 55,
        b1: 110_000_000,
        b2: 776_278_396_540,
    },
    Level {
        digits: 60,
        b1: 260_000_000,
        b2: 3_178_559_884_516,
    },
    Level {
        digits: 65,
        b1: 850_000_000,
        b2: 15_892_628_251_516,
    },
];

/// Stage 1 bound of P-1 at a level, relative to the `B1` of the curves.
const PM1_B1_RATIO: usize = 20;

/// Stage 2 of P-1 may cost this many times as much as its stage 1.
const PM1_STAGE2_RATIO: f64 = 1.0;

/// At the last level, rounds of the expected number of curves before giving up: the probability
/// to miss a factor of the size of the level is then `e^-TOP_LEVEL_ROUNDS`.
const TOP_LEVEL_ROUNDS: usize = 10;

/// Progress of the search on a number, inherited by its cofactors.
#[derive(Clone)]
struct Progress {
    /// Current level, and number of curves already run at this level.
    level: usize,
    curves: usize,
    /// Rounds of curves completed at the last level.
    rounds: usize,
    /// P-1 state (`None` once P-1 is useless: it found all the factors at once), and the last
    /// level it ran at.
    pm1: Option<Pm1>,
    pm1_level: Option<usize>,
}

impl Progress {
    fn new(pm1: Option<Pm1>) -> Self {
        Progress {
            level: 0,
            curves: 0,
            rounds: 0,
            pm1,
            pm1_level: None,
        }
    }
}

/// What the curves of a level share.
struct Context<'a> {
    rand: RandState<'a>,
    /// Stage 1 multipliers by `b1`, and stage 2 plans by bounds and size of the number.
    multipliers: HashMap<usize, Rc<Integer>>,
    plans: HashMap<(usize, usize, u32), Rc<Stage2Plan>>,
    #[cfg(feature = "progress-bar")]
    pb: Option<&'a ProgressBar>,
}

impl Context<'_> {
    fn multiplier(&mut self, b1: usize) -> Rc<Integer> {
        self.multipliers
            .entry(b1)
            .or_insert_with(|| Rc::new(stage1_multiplier(b1)))
            .clone()
    }

    fn plan(&mut self, n: &Integer, b1: usize, b2: usize) -> Rc<Stage2Plan> {
        self.plans
            .entry((b1, b2, n.significant_bits()))
            .or_insert_with(|| Rc::new(Stage2Plan::new(n, b1, b2)))
            .clone()
    }
}

/// Factors `n` (see [`crate::ecm()`]), with the random state seeded by `seed`.
///
/// # Panics
///
/// If `n` is not positive.
pub fn factor(
    n: &Integer,
    seed: usize,
    #[cfg(feature = "progress-bar")] pb: Option<&ProgressBar>,
) -> Result<HashMap<Integer, usize>, Error> {
    assert!(*n > 0, "only positive numbers can be factored");
    let (mut factors, n) = trial_division(n);
    let mut ctx = Context {
        rand: rand_state(seed),
        multipliers: HashMap::new(),
        plans: HashMap::new(),
        #[cfg(feature = "progress-bar")]
        pb,
    };

    let mut queue = Vec::new();
    sort_factor(
        n,
        1,
        Progress::new(Some(Pm1::new())),
        &mut factors,
        &mut queue,
    );
    while let Some((n, exponent, mut progress)) = queue.pop() {
        let factor = find_factor(&n, &mut progress, &mut ctx)?;
        let mut cofactor = n;
        let mut multiplicity = 0;
        while cofactor.is_divisible(&factor) {
            cofactor /= &factor;
            multiplicity += 1;
        }
        // A composite factor was found by a curve or P-1 that found all its factors at once:
        // other curves split it, from the first level.
        sort_factor(
            factor,
            exponent * multiplicity,
            Progress::new(None),
            &mut factors,
            &mut queue,
        );
        sort_factor(cofactor, exponent, progress, &mut factors, &mut queue);
    }
    Ok(factors)
}

/// Index of the last level for `n`: its smallest factor has at most half its digits.
fn top_level(n: &Integer) -> usize {
    let digits = (n.significant_bits() as f64 * std::f64::consts::LOG10_2).ceil() as u32;
    let half = digits.div_ceil(2);
    LEVELS
        .iter()
        .position(|level| level.digits >= half)
        .unwrap_or(LEVELS.len() - 1)
}

/// A proper factor of the composite `n`, resuming the search from `progress`.
fn find_factor(
    n: &Integer,
    progress: &mut Progress,
    ctx: &mut Context<'_>,
) -> Result<Integer, Error> {
    let top = top_level(n);
    let param = Param::default();
    loop {
        let index = progress.level.min(top);
        let level = &LEVELS[index];
        if let Some(g) = pm1(n, index, progress) {
            return Ok(g);
        }

        let k = ctx.multiplier(level.b1);
        let plan = ctx.plan(n, level.b1, level.b2);
        let prob = ecm_prob(level.b1 as f64, plan.b2() as f64, level.digits as f64);
        let curves = (1.0 / prob).ceil().max(1.0) as usize;

        #[cfg(feature = "progress-bar")]
        if let Some(pb) = ctx.pb {
            pb.set_length(curves as u64);
            pb.set_position(progress.curves as u64);
        }

        while progress.curves < curves {
            progress.curves += 1;
            #[cfg(feature = "progress-bar")]
            if let Some(pb) = ctx.pb {
                pb.inc(1);
            }
            let sigma = random_sigma(n, param, &mut ctx.rand);
            match run_curve(n, param, &sigma, &k, &plan) {
                CurveOutcome::Setup(g) | CurveOutcome::Stage1(g) | CurveOutcome::Stage2(g) => {
                    return Ok(g);
                }
                CurveOutcome::Failed => {}
            }
        }

        progress.curves = 0;
        if index < top {
            progress.level = index + 1;
        } else {
            progress.rounds += 1;
            if progress.rounds >= TOP_LEVEL_ROUNDS {
                return Err(Error::ECMFailed);
            }
        }
    }
}

/// Extends P-1 for the level `index` (once per level): returns a proper factor of `n` if it
/// finds one.
fn pm1(n: &Integer, index: usize, progress: &mut Progress) -> Option<Integer> {
    if progress.pm1_level >= Some(index) {
        return None;
    }
    progress.pm1_level = Some(index);
    let pm1 = progress.pm1.as_mut()?;
    let b1 = LEVELS[index].b1 * PM1_B1_RATIO;
    if b1 <= pm1.b1() {
        return None;
    }
    let g = pm1.stage1(n, b1);
    if &g == n {
        // Every p - 1 is smooth: P-1 cannot separate the factors.
        progress.pm1 = None;
        return None;
    }
    if g != 1 {
        return Some(g);
    }
    let plan = Stage2Plan::new(n, b1, pm1_b2(n, b1));
    let g = pm1.stage2(n, &plan);
    (g != 1 && &g != n).then_some(g)
}

/// Stage 2 bound of P-1 after a stage 1 up to `b1`: the largest `b1*2^i` whose stage 2 costs at
/// most [`PM1_STAGE2_RATIO`] times stage 1 (`1.44*b1` modular squarings).
fn pm1_b2(n: &Integer, b1: usize) -> usize {
    let bits = n.significant_bits() as usize;
    let budget = PM1_STAGE2_RATIO * 1.44 * b1 as f64 * 1.3 * Costs::new(bits).mul();
    let mut b2 = b1;
    while b2 < usize::MAX / 4 && Stage2Plan::cost_at_most(bits, b1, 2 * b2, budget) {
        b2 *= 2;
    }
    b2
}

#[cfg(test)]
mod tests {
    use super::*;
    use rug::ops::Pow;
    use std::str::FromStr;

    fn factor(n: &Integer, seed: usize) -> HashMap<Integer, usize> {
        super::factor(
            n,
            seed,
            #[cfg(feature = "progress-bar")]
            None,
        )
        .unwrap()
    }

    #[test]
    fn levels() {
        for window in LEVELS.windows(2) {
            assert!(window[0].digits < window[1].digits);
            assert!(window[0].b1 < window[1].b1 && window[0].b2 < window[1].b2);
        }
        let digits = |d: u32| Integer::from(10).pow(d - 1) + 1u32;
        assert_eq!(top_level(&digits(10)), 0);
        assert_eq!(top_level(&digits(20)), 0);
        assert_eq!(top_level(&digits(21)), 1);
        assert_eq!(top_level(&digits(60)), 4);
        assert_eq!(top_level(&digits(61)), 5);
        assert_eq!(top_level(&digits(1000)), LEVELS.len() - 1);
    }

    #[test]
    fn expected_curves() {
        // The expected curves of GMP-ECM's table 1 (with a stage 2 without Brent-Suyama's
        // extension): the bounds are the right ones for each size.
        for (level, curves) in LEVELS[2..6].iter().zip([86.0, 221.0, 454.0, 986.0]) {
            let expected = 1.0 / ecm_prob(level.b1 as f64, level.b2 as f64, level.digits as f64);
            assert!((expected / curves - 1.0).abs() < 0.01, "{}", level.b1);
        }
    }

    #[test]
    fn factors_by_size() {
        // Factors of 6 to 16 digits (2 of them with multiplicity), found in the order of their
        // size, whatever the seed.
        let primes = [
            "100003",
            "1000000007",
            "2147483647",
            "99999999977",
            "1000000000039",
            "9007199254740881",
        ]
        .map(|p| Integer::from_str(p).unwrap());
        let n = primes.iter().product::<Integer>() * &primes[0] * &primes[3];
        let mut expected: HashMap<Integer, usize> = primes.iter().map(|p| (p.clone(), 1)).collect();
        expected.insert(primes[0].clone(), 2);
        expected.insert(primes[3].clone(), 2);
        for seed in 0..3 {
            assert_eq!(factor(&n, seed), expected);
        }
    }

    #[test]
    fn pm1_finds_large_factor() {
        // p - 1 = 2^2 * 5743 * 6037 * 30259 * 30949 * 34039 * 300007: P-1 finds this 28-digit
        // factor at the second level (B1 = 40000), in stage 2, long before the curves could.
        let p = Integer::from_str("1326262092842391910564284053").unwrap();
        let q = Integer::from(10).pow(100) + 267u32;
        let n = Integer::from(&p * &q);
        let mut progress = Progress::new(Some(Pm1::new()));
        let found = (0..3).find_map(|level| pm1(&n, level, &mut progress).map(|g| (level, g)));
        assert_eq!(found, Some((1, p)));
        assert_eq!(progress.pm1.unwrap().b1(), LEVELS[1].b1 * PM1_B1_RATIO);
    }

    #[test]
    fn deterministic() {
        let n = Integer::from_str("7060005655815754299976961394452809").unwrap();
        assert_eq!(factor(&n, 7), factor(&n, 7));
    }
}
