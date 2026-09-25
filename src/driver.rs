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
//!
//! With fixed bounds ([`crate::ecm_with_params`]), each composite runs curves with these bounds
//! only. With [`crate::Algorithm::Pm1`] or [`crate::Algorithm::Pp1`], P-1 or P+1 runs alone
//! (with the P-1 bounds of the levels, or once with fixed bounds).
//!
//! P+1 is not part of the default search: after P-1 (which finds the factors with a smooth
//! `p - 1`, P+1 with the seed `2/7` finds half of them again), it only finds the factors `p = 2
//! mod 3` with a smooth `p + 1`. At its best `B1` (1 to 2 times the `B1` of the curves), it finds
//! them at the rate of the curves, and with more, below it: with the costs measured on 40, 60 and
//! 100-digit numbers, the expected time to find factors of 15 to 30 digits changes by less than
//! 0.1% (and grows by 0.3-4% with 10 times the `B1` of the curves).
//!
//! The searches report [`Event`]s to the handler of the [`crate::Factorizer`], which may
//! interrupt them.

use crate::{
    base2::{Base2Form, Base2Mode},
    cost::Costs,
    ecm::{
        CurveOutcome, Error, Param, is_prime, perfect_power, random_sigma, run_curve_timed,
        small_factor, stage1_multiplier, trial_division,
    },
    events::{Event, EventHandler, Method},
    factorizer::{Algorithm, Factorization},
    parallel::Pool,
    pm1::Pm1,
    pp1::Pp1,
    rho::ecm_prob,
    stage2::Stage2Plan,
    stop::Stop,
};
use rug::{Integer, rand::RandState};

mod pipeline;
#[cfg(test)]
pub(crate) use pipeline::PANIC_SIGMA;
use std::{
    collections::HashMap,
    sync::Arc,
    time::{Duration, Instant},
};

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

/// Stage 1 bound of P-1 at a level (and of P+1, when it runs alone), relative to the `B1` of the
/// curves.
const PM1_B1_RATIO: usize = 20;

/// Stage 2 of P-1 may cost this many times as much as its stage 1.
const PM1_STAGE2_RATIO: f64 = 1.0;

/// At the last level, rounds of the expected number of curves before giving up: the probability
/// to miss a factor of the size of the level is then `e^-TOP_LEVEL_ROUNDS`.
const TOP_LEVEL_ROUNDS: usize = 10;

/// Stage 2 bound for curves with the stage 1 bound `b1` when none is given: GMP-ECM's default
/// `B2` for the `B1` of the levels, interpolated (and extrapolated) linearly in `log B1`,
/// `log B2`, rounded to an even number.
#[must_use]
pub(crate) fn default_b2(b1: usize) -> usize {
    let i = LEVELS
        .iter()
        .position(|level| level.b1 >= b1)
        .unwrap_or(LEVELS.len() - 1)
        .clamp(1, LEVELS.len() - 1);
    let (lo, hi) = (&LEVELS[i - 1], &LEVELS[i]);
    let slope = (hi.b2 as f64 / lo.b2 as f64).ln() / (hi.b1 as f64 / lo.b1 as f64).ln();
    // Below the first level, the ratio B2/B1 of the first level.
    let slope = if b1 < lo.b1 { 1.0 } else { slope };
    let b2 = lo.b2 as f64 * (b1 as f64 / lo.b1 as f64).powf(slope);
    (b2.min(usize::MAX as f64 / 2.0) as usize & !1).max(4)
}

/// How the composite parts are searched.
#[derive(Debug, Clone, Copy)]
pub(crate) enum Mode {
    /// Levels of increasing factor size, with P-1 before each ([`crate::ecm()`]).
    Levels,
    /// Fixed bounds, with at most `curves` curves per composite part ([`crate::ecm_with_params`]).
    Fixed {
        b1: usize,
        b2: usize,
        curves: Option<usize>,
    },
}

impl Mode {
    /// Fixed bounds, checked.
    pub(crate) fn fixed(b1: usize, b2: usize, curves: Option<usize>) -> Result<Self, Error> {
        if !b1.is_multiple_of(2) || !b2.is_multiple_of(2) {
            return Err(Error::BoundsNotEven);
        }
        // Stage 2 skips the primes of its wheel, which must be at most `b1`.
        if b1 < 6 || b2 < 4 {
            return Err(Error::BoundsTooSmall);
        }
        Ok(Self::Fixed { b1, b2, curves })
    }
}

/// P-1 or P+1 on a number, resumable.
#[derive(Clone)]
enum PlusMinus {
    Pm1(Pm1),
    Pp1(Pp1),
}

impl PlusMinus {
    fn b1(&self) -> usize {
        match self {
            Self::Pm1(pm1) => pm1.b1(),
            Self::Pp1(pp1) => pp1.b1(),
        }
    }

    fn stage1(
        &mut self,
        n: &Integer,
        b1: usize,
        base2: Option<Base2Form>,
        stop: Stop<'_>,
    ) -> Integer {
        match self {
            Self::Pm1(pm1) => pm1.stage1_until(n, b1, base2, stop),
            Self::Pp1(pp1) => pp1.stage1_until(n, b1, base2, stop),
        }
    }

    fn stage2(
        &self,
        n: &Integer,
        plan: &Stage2Plan,
        base2: Option<Base2Form>,
        stop: Stop<'_>,
    ) -> Integer {
        match self {
            Self::Pm1(pm1) => pm1.stage2_until(n, plan, base2, stop),
            Self::Pp1(pp1) => pp1.stage2_until(n, plan, base2, stop),
        }
    }

    /// How stage 1 and stage 2 find factors.
    fn methods(&self) -> [Method; 2] {
        match self {
            Self::Pm1(_) => [Method::Pm1Stage1, Method::Pm1Stage2],
            Self::Pp1(_) => [Method::Pp1Stage1, Method::Pp1Stage2],
        }
    }

    /// The event of a run.
    fn event<'a>(
        &self,
        n: &'a Integer,
        b1: usize,
        b2: Option<usize>,
        [stage1, stage2]: [Duration; 2],
    ) -> Event<'a> {
        match self {
            Self::Pm1(_) => Event::Pm1 {
                n,
                b1,
                b2,
                stage1,
                stage2,
            },
            Self::Pp1(_) => Event::Pp1 {
                n,
                b1,
                b2,
                stage1,
                stage2,
            },
        }
    }
}

/// Progress of the search on a number, inherited by its cofactors.
#[derive(Clone)]
struct Progress {
    /// Current level, and number of curves already run at this level.
    level: usize,
    curves: usize,
    /// Rounds of curves completed at the last level.
    rounds: usize,
    /// P-1 (or P+1) state (`None` once it is useless: it found all the factors at once), and
    /// the last level it ran at.
    pm: Option<PlusMinus>,
    pm_level: Option<usize>,
}

impl Progress {
    fn new(pm: Option<PlusMinus>) -> Self {
        Self {
            level: 0,
            curves: 0,
            rounds: 0,
            pm,
            pm_level: None,
        }
    }
}

/// The event handler, and whether it (or the [`Stop`]) interrupted the factorization.
struct Events<'h, H> {
    handler: &'h mut H,
    stop: Stop<'h>,
    interrupted: bool,
}

impl<H: EventHandler> Events<'_, H> {
    /// Whether the factorization is interrupted: by the handler, or now by the [`Stop`].
    fn interrupted(&mut self) -> bool {
        if !self.interrupted && self.stop.requested() {
            self.interrupted = true;
        }
        self.interrupted
    }

    /// [`Error::Interrupted`] if the factorization is interrupted (see [`Events::interrupted`]).
    fn check(&mut self) -> Result<(), Error> {
        if self.interrupted() {
            Err(Error::Interrupted)
        } else {
            Ok(())
        }
    }

    /// Reports `event`, unless the handler already interrupted the factorization: then (or if it
    /// does now) returns [`Error::Interrupted`].
    #[inline]
    fn emit(&mut self, event: Event<'_>) -> Result<(), Error> {
        if !H::ENABLED {
            return Ok(());
        }
        if !self.interrupted && self.handler.handle(&event).is_continue() {
            return Ok(());
        }
        self.interrupted = true;
        Err(Error::Interrupted)
    }
}

/// Time elapsed since `start` (zero without events).
fn since(start: Option<Instant>) -> Duration {
    start.map_or(Duration::ZERO, |start| start.elapsed())
}

/// Stage 2 plans by bounds `(b1, b2)`, size of the number (bits) and form of its special
/// reduction in stage 1.
type PlanKey = (usize, usize, u32, Option<Base2Form>);

/// A factorization: the options, and what the searches share (random state, stage 1
/// multipliers and stage 2 plans).
pub(crate) struct Engine<'a, 'r, H> {
    mode: Mode,
    param: Param,
    /// Parameter of the next curve, if fixed.
    sigma: Option<Integer>,
    algorithm: Algorithm,
    /// Seed of P-1 or P+1, if not the default one.
    x0: Option<(Integer, Integer)>,
    max_memory: usize,
    base2: Base2Mode,
    rand: &'a mut RandState<'r>,
    events: Events<'a, H>,
    /// Stage 1 multipliers by `b1`, and stage 2 plans by bounds, size of the number and form
    /// of its special reduction, with the form stage 2 uses (see [`Engine::plan`]).
    multipliers: HashMap<usize, Arc<Integer>>,
    plans: HashMap<PlanKey, (Arc<Stage2Plan>, Option<Base2Form>)>,
    /// Threads running curves (see [`crate::Factorizer::threads`]), and their pool, started
    /// at the first level that runs curves in parallel.
    threads: usize,
    pool: Option<Pool<pipeline::Output>>,
}

impl<'a, 'r, H: EventHandler> Engine<'a, 'r, H> {
    #[allow(clippy::too_many_arguments, reason = "the options of a Factorizer")]
    pub(crate) fn new(
        mode: Mode,
        param: Param,
        sigma: Option<Integer>,
        (algorithm, x0): (Algorithm, Option<(Integer, Integer)>),
        max_memory: usize,
        base2: Base2Mode,
        rand: &'a mut RandState<'r>,
        (handler, stop): (&'a mut H, Stop<'a>),
        threads: usize,
    ) -> Self {
        Self {
            mode,
            param,
            sigma,
            algorithm,
            x0,
            max_memory,
            base2,
            rand,
            events: Events {
                handler,
                stop,
                interrupted: false,
            },
            multipliers: HashMap::new(),
            plans: HashMap::new(),
            threads,
            pool: None,
        }
    }

    fn multiplier(&mut self, b1: usize) -> Arc<Integer> {
        self.multipliers
            .entry(b1)
            .or_insert_with(|| Arc::new(stage1_multiplier(b1)))
            .clone()
    }

    /// The stage 2 plan for `n` with the special reduction modulo `base2` (if any) in stage 1,
    /// and the form of the special reduction in stage 2: `base2`, or none if computing modulo
    /// `n` is cheaper there (with [`Base2Mode::Auto`], see [`Stage2Plan::cheapest`]).
    fn plan(
        &mut self,
        n: &Integer,
        b1: usize,
        b2: usize,
        base2: Option<Base2Form>,
    ) -> (Arc<Stage2Plan>, Option<Base2Form>) {
        let max_memory = self.max_memory;
        let either = self.base2 == Base2Mode::Auto;
        self.plans
            .entry((b1, b2, n.significant_bits(), base2))
            .or_insert_with(|| {
                let (plan, form) = Stage2Plan::cheapest(n, (b1, b2), max_memory, base2, either);
                (Arc::new(plan), form)
            })
            .clone()
    }

    /// Factors `n` completely (see [`crate::Factorizer::factor_partial`]).
    ///
    /// # Panics
    ///
    /// If `n` is not positive.
    pub(crate) fn factor(&mut self, n: &Integer) -> Factorization {
        assert!(*n > 0, "only positive numbers can be factored");
        let (factors, cofactor) = trial_division(n);
        let mut done = Factorization {
            primes: HashMap::new(),
            unfactored: Vec::new(),
            error: None,
        };
        let mut queue = Vec::new();
        let _ = self.events.emit(Event::TrialDivision {
            factors: &factors,
            cofactor: &cofactor,
        });
        let mut small: Vec<_> = factors.into_iter().collect();
        small.sort_unstable();
        for (p, exponent) in small {
            let _ = self.events.emit(Event::Prime { p: &p, exponent });
            done.primes.insert(p, exponent);
        }
        let pm = matches!(self.mode, Mode::Levels).then(|| self.plus_minus());
        self.sort(cofactor, 1, Progress::new(pm), &mut done.primes, &mut queue);

        while let Some((n, exponent, mut progress)) = queue.pop() {
            if self.events.interrupted() {
                done.unfactored.push((n, exponent));
                continue;
            }
            let (factor, method) = match self.find(&n, &mut progress) {
                Ok(found) => found,
                Err(error) => {
                    done.unfactored.push((n, exponent));
                    if error == Error::Interrupted || done.error.is_none() {
                        done.error = Some(error);
                    }
                    continue;
                }
            };
            let _ = self.events.emit(Event::Factor {
                n: &n,
                factor: &factor,
                method,
            });
            let mut cofactor = n;
            let mut multiplicity = 0;
            while cofactor.is_divisible(&factor) {
                cofactor /= &factor;
                multiplicity += 1;
            }
            // A composite factor was found by a curve or P-1 that found all its factors at
            // once: other curves split it, from the first level.
            self.sort(
                factor,
                exponent * multiplicity,
                Progress::new(None),
                &mut done.primes,
                &mut queue,
            );
            self.sort(cofactor, exponent, progress, &mut done.primes, &mut queue);
        }
        if self.events.interrupted {
            done.error = Some(Error::Interrupted);
        }
        done
    }

    /// Records `n^exponent`: a prime goes to `primes`, a perfect power is reduced to its root,
    /// and any other composite is queued for factorization, with `progress`.
    fn sort(
        &mut self,
        n: Integer,
        exponent: usize,
        progress: Progress,
        primes: &mut HashMap<Integer, usize>,
        queue: &mut Vec<(Integer, usize, Progress)>,
    ) {
        if n == 1 {
            return;
        }
        if is_prime(&n) {
            let _ = self.events.emit(Event::Prime { p: &n, exponent });
            *primes.entry(n).or_insert(0) += exponent;
            return;
        }
        match perfect_power(&n) {
            Some((root, power)) => {
                let _ = self.events.emit(Event::Factor {
                    n: &n,
                    factor: &root,
                    method: Method::PerfectPower,
                });
                self.sort(root, exponent * power as usize, progress, primes, queue);
            }
            None => queue.push((n, exponent, progress)),
        }
    }

    /// A proper factor of `n > 1` (see [`crate::Factorizer::find_factor`]), after trial
    /// division if `trial`.
    ///
    /// # Panics
    ///
    /// If `n <= 1`.
    pub(crate) fn find_one(&mut self, n: &Integer, trial: bool) -> Result<Integer, Error> {
        assert!(*n > 1, "only numbers greater than 1 have a proper factor");
        let (factor, method) = self.find_one_method(n, trial)?;
        let _ = self.events.emit(Event::Factor {
            n,
            factor: &factor,
            method,
        });
        Ok(factor)
    }

    fn find_one_method(&mut self, n: &Integer, trial: bool) -> Result<(Integer, Method), Error> {
        if let Some(p) = small_factor(n).filter(|_| trial) {
            return Ok((p, Method::TrialDivision));
        }
        if is_prime(n) {
            return Err(Error::NumberIsPrime);
        }
        // A perfect power is split by its root at once (modulo a power of a small prime, the
        // curves may all find its whole power).
        if let Some((root, _)) = perfect_power(n) {
            return Ok((root, Method::PerfectPower));
        }
        let pm = matches!(self.mode, Mode::Levels).then(|| self.plus_minus());
        self.find(n, &mut Progress::new(pm))
    }

    /// P-1, or P+1 with [`Algorithm::Pp1`], from the start.
    fn plus_minus(&self) -> PlusMinus {
        let x0 = self.x0.clone();
        match self.algorithm {
            Algorithm::Pp1 => PlusMinus::Pp1(x0.map_or_else(Pp1::default, |(a, b)| Pp1::new(a, b))),
            _ => PlusMinus::Pm1(x0.map_or_else(Pm1::new, |(a, b)| Pm1::with_x0(a, b))),
        }
    }

    /// A proper factor of the composite `n`, resuming the search from `progress`.
    fn find(&mut self, n: &Integer, progress: &mut Progress) -> Result<(Integer, Method), Error> {
        let base2 = self.base2.form(n).map_err(Error::InvalidOption)?;
        if let Some(form) = base2 {
            self.events.emit(Event::Base2 {
                n,
                k: form.signed(),
            })?;
        }
        match self.mode {
            Mode::Levels if self.parallel() && self.algorithm == Algorithm::Ecm => {
                // The first level takes less time than starting threads.
                match self.find_by_levels(n, progress, base2, 1)? {
                    Some(found) => Ok(found),
                    None => self.find_parallel(n, progress, base2, None),
                }
            }
            Mode::Levels => self
                .find_by_levels(n, progress, base2, usize::MAX)
                .map(|found| found.expect("no last level")),
            Mode::Fixed { b1, b2, .. } if self.algorithm != Algorithm::Ecm => {
                self.pm_fixed(n, b1, b2, base2)
            }
            Mode::Fixed { b1, b2, curves } if self.parallel() => {
                self.find_parallel(n, progress, base2, Some((b1, b2, curves)))
            }
            Mode::Fixed { b1, b2, curves } => self.curves(n, (b1, b2), None, curves, &mut 0, base2),
        }
    }

    /// Searches `n` level by level, from `progress` on (with the special reduction modulo
    /// `base2` if any), up to the level `until` (excluded): `None` when it is reached.
    fn find_by_levels(
        &mut self,
        n: &Integer,
        progress: &mut Progress,
        base2: Option<Base2Form>,
        until: usize,
    ) -> Result<Option<(Integer, Method)>, Error> {
        let top = top_level(n);
        loop {
            if progress.level >= until {
                return Ok(None);
            }
            let index = progress.level.min(top);
            let level = &LEVELS[index];
            if let Some(found) = self.pm_level(n, index, progress, base2)? {
                return Ok(Some(found));
            }

            if self.algorithm != Algorithm::Ecm {
                if index == top || progress.pm.is_none() {
                    return Err(Error::ECMFailed);
                }
                progress.level = index + 1;
                continue;
            }

            let prob = {
                let (plan, _) = self.plan(n, level.b1, level.b2, base2);
                ecm_prob(level.b1 as f64, plan.b2() as f64, f64::from(level.digits))
            };
            let curves = (1.0 / prob).ceil().max(1.0) as usize;
            match self.curves(
                n,
                (level.b1, level.b2),
                Some(level.digits),
                Some(curves),
                &mut progress.curves,
                base2,
            ) {
                Err(Error::ECMFailed) => {}
                result => return result.map(Some),
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

    /// Runs curves on `n` with the bounds `(b1, b2)`, `done` of them already run, until
    /// `curves` in all (without limit if `None`) or a factor is found.
    fn curves(
        &mut self,
        n: &Integer,
        (b1, b2): (usize, usize),
        digits: Option<u32>,
        curves: Option<usize>,
        done: &mut usize,
        base2: Option<Base2Form>,
    ) -> Result<(Integer, Method), Error> {
        if curves.is_some_and(|curves| *done >= curves) {
            // Resumed on a cofactor after all the curves of the level.
            return Err(Error::ECMFailed);
        }
        let k = self.multiplier(b1);
        let (plan, base2_stage2) = self.plan(n, b1, b2, base2);
        let base2 = [base2, base2_stage2];
        self.events.emit(Event::Level {
            n,
            digits,
            b1,
            b2: plan.b2(),
            curves,
            done: *done,
        })?;

        let param = self.param;
        while curves.is_none_or(|curves| *done < curves) {
            self.events.check()?;
            *done += 1;
            let sigma = match &mut self.sigma {
                Some(sigma) => {
                    let next = Integer::from(&*sigma + 1u32);
                    std::mem::replace(sigma, next)
                }
                None => random_sigma(n, param, self.rand),
            };
            let stop = self.events.stop;
            let (outcome, [stage1, stage2]) = if H::ENABLED {
                run_curve_timed::<true>(n, param, &sigma, &k, &plan, base2, stop)
            } else {
                run_curve_timed::<false>(n, param, &sigma, &k, &plan, base2, stop)
            };
            if outcome == CurveOutcome::Failed {
                // The curve may have stopped early: not reported.
                self.events.check()?;
            }
            let event = self.events.emit(Event::Curve {
                n,
                param,
                sigma: &sigma,
                index: *done,
                stage1,
                stage2,
            });
            match outcome {
                CurveOutcome::Setup(g) => return Ok((g, Method::EcmSetup)),
                CurveOutcome::Stage1(g) => return Ok((g, Method::EcmStage1)),
                CurveOutcome::Stage2(g) => return Ok((g, Method::EcmStage2)),
                CurveOutcome::Failed => event?,
            }
        }
        Err(Error::ECMFailed)
    }

    /// Extends P-1 (or P+1) for the level `index` (once per level): returns a proper factor of
    /// `n` if it finds one.
    fn pm_level(
        &mut self,
        n: &Integer,
        index: usize,
        progress: &mut Progress,
        base2: Option<Base2Form>,
    ) -> Result<Option<(Integer, Method)>, Error> {
        if progress.pm_level >= Some(index) {
            return Ok(None);
        }
        progress.pm_level = Some(index);
        let Some(pm) = progress.pm.as_mut() else {
            return Ok(None);
        };
        let b1 = LEVELS[index].b1 * PM1_B1_RATIO;
        if b1 <= pm.b1() {
            return Ok(None);
        }
        let max_memory = self.max_memory;
        let either = self.base2 == Base2Mode::Auto;
        let plan = || {
            let b2 = pm1_b2(n, b1, max_memory, (base2, either));
            Stage2Plan::cheapest(n, (b1, b2), max_memory, base2, either)
        };
        let (found, all) = run_plus_minus(pm, n, (b1, plan), base2, &mut self.events)?;
        if all {
            // Every p - 1 (or p + 1) is smooth: P-1 (or P+1) cannot separate the factors.
            progress.pm = None;
        }
        Ok(found)
    }

    /// P-1 (or P+1) with fixed bounds.
    fn pm_fixed(
        &mut self,
        n: &Integer,
        b1: usize,
        b2: usize,
        base2: Option<Base2Form>,
    ) -> Result<(Integer, Method), Error> {
        let mut pm = self.plus_minus();
        let max_memory = self.max_memory;
        let either = self.base2 == Base2Mode::Auto;
        let plan = || Stage2Plan::cheapest(n, (b1, b2), max_memory, base2, either);
        let (found, _) = run_plus_minus(&mut pm, n, (b1, plan), base2, &mut self.events)?;
        found.ok_or(Error::ECMFailed)
    }
}

/// Runs P-1 (or P+1) on `n` up to `b1` (resuming `pm`, with the special reduction modulo
/// `base2` if any), then stage 2 with `plan()` (the plan and the form of its special reduction)
/// if stage 1 finds nothing: returns the proper factor found, if any, and whether stage 1 found
/// all the factors of `n` at once.
fn run_plus_minus<H: EventHandler>(
    pm: &mut PlusMinus,
    n: &Integer,
    (b1, plan): (usize, impl FnOnce() -> (Stage2Plan, Option<Base2Form>)),
    base2: Option<Base2Form>,
    events: &mut Events<'_, H>,
) -> Result<(Option<(Integer, Method)>, bool), Error> {
    let run = PlusMinusRun::new(pm, n, (b1, plan), base2, H::ENABLED, events.stop);
    run.report(pm, n, b1, events)
}

/// A run of P-1 (or P+1), to report (see [`run_plus_minus`]).
pub(super) struct PlusMinusRun {
    /// The gcd after stage 1.
    stage1: Integer,
    /// The gcd after stage 2 and the bound it covered, if it ran.
    stage2: Option<(Integer, usize)>,
    /// Durations of the stages (zero if not measured).
    times: [Duration; 2],
    /// Whether it stopped early, without finding a factor: not reported.
    stopped: bool,
}

impl PlusMinusRun {
    /// Runs `pm` (see [`run_plus_minus`]), measuring the durations if `timed`.
    fn new(
        pm: &mut PlusMinus,
        n: &Integer,
        (b1, plan): (usize, impl FnOnce() -> (Stage2Plan, Option<Base2Form>)),
        base2: Option<Base2Form>,
        timed: bool,
        stop: Stop<'_>,
    ) -> Self {
        let start = timed.then(Instant::now);
        let g = pm.stage1(n, b1, base2, stop);
        let stage1 = since(start);
        if g != 1 {
            return Self {
                stage1: g,
                stage2: None,
                times: [stage1, Duration::ZERO],
                stopped: false,
            };
        }
        if stop.requested() {
            // Stage 1 may have stopped early.
            return Self {
                stage1: g,
                stage2: None,
                times: [stage1, Duration::ZERO],
                stopped: true,
            };
        }
        let start = timed.then(Instant::now);
        let (plan, base2) = plan();
        let g = pm.stage2(n, &plan, base2, stop);
        let stage2 = since(start);
        let stopped = (g == 1 || &g == n) && stop.requested();
        Self {
            stage1: Integer::from(1),
            stage2: Some((g, plan.b2())),
            times: [stage1, stage2],
            stopped,
        }
    }

    /// Whether it found a proper factor of `n`.
    fn found(&self, n: &Integer) -> bool {
        let proper = |g: &Integer| *g != 1 && g != n;
        proper(&self.stage1) || self.stage2.as_ref().is_some_and(|(g, _)| proper(g))
    }

    /// Whether stage 1 found all the factors of `n` at once.
    fn all(&self, n: &Integer) -> bool {
        &self.stage1 == n
    }

    /// Reports the run of `pm` up to `b1` on `n`: returns the proper factor found, if any, and
    /// whether stage 1 found all the factors of `n` at once, or [`Error::Interrupted`].
    fn report<H: EventHandler>(
        self,
        pm: &PlusMinus,
        n: &Integer,
        b1: usize,
        events: &mut Events<'_, H>,
    ) -> Result<(Option<(Integer, Method)>, bool), Error> {
        if self.stopped {
            events.interrupted = true;
            return Err(Error::Interrupted);
        }
        let [method1, method2] = pm.methods();
        let all = self.all(n);
        let (found, b2) = match self.stage2 {
            None => ((!all).then_some((self.stage1, method1)), None),
            Some((g, b2)) => ((g != 1 && &g != n).then_some((g, method2)), Some(b2)),
        };
        let event = events.emit(pm.event(n, b1, b2, self.times));
        // A factor found is kept, even if the handler interrupts.
        if found.is_none() {
            event?;
        }
        Ok((found, all))
    }
}

/// Index of the last level for `n`: its smallest factor has at most half its digits.
fn top_level(n: &Integer) -> usize {
    let digits = (f64::from(n.significant_bits()) * std::f64::consts::LOG10_2).ceil() as u32;
    let half = digits.div_ceil(2);
    LEVELS
        .iter()
        .position(|level| level.digits >= half)
        .unwrap_or(LEVELS.len() - 1)
}

/// Stage 2 bound of P-1 after a stage 1 up to `b1`: the largest `b1*2^i` whose stage 2 costs at
/// most [`PM1_STAGE2_RATIO`] times stage 1 (`1.44*b1` modular squarings), with the special
/// reduction modulo `base2` if any (in stage 2 too, unless `either` and computing modulo `n`
/// is cheaper there, see [`Stage2Plan::cheapest`]).
fn pm1_b2(
    n: &Integer,
    b1: usize,
    max_memory: usize,
    (base2, either): (Option<Base2Form>, bool),
) -> usize {
    let bits = n.significant_bits() as usize;
    let costs = Costs::modulo(bits, base2);
    let budget = PM1_STAGE2_RATIO * 1.44 * b1 as f64 * 1.3 * costs.mul();
    let fits = |b2| {
        Stage2Plan::cost_at_most((bits, base2), b1, b2, budget, max_memory)
            || (either
                && base2.is_some()
                && Stage2Plan::cost_at_most((bits, None), b1, b2, budget, max_memory))
    };
    let mut b2 = b1;
    while b2 < usize::MAX / 4 && fits(2 * b2) {
        b2 *= 2;
    }
    b2
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ecm::rand_state, stage2::MAX_POLY_MEMORY};
    use rug::ops::Pow;
    use std::str::FromStr;

    fn factor(n: &Integer, seed: u64) -> HashMap<Integer, usize> {
        crate::Factorizer::new().seed(seed).factor(n).unwrap()
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
            let expected =
                1.0 / ecm_prob(level.b1 as f64, level.b2 as f64, f64::from(level.digits));
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
        let mut progress = Progress::new(Some(PlusMinus::Pm1(Pm1::new())));
        let (mut rand, mut handler) = (rand_state(0), crate::events::NoEvents);
        let mut engine = Engine::new(
            Mode::Levels,
            Param::default(),
            None,
            (Algorithm::Ecm, None),
            MAX_POLY_MEMORY,
            Base2Mode::Auto,
            &mut rand,
            (&mut handler, Stop::NEVER),
            1,
        );
        let found = (0..3).find_map(|level| {
            let found = engine.pm_level(&n, level, &mut progress, None).unwrap();
            found.map(|(g, method)| (level, g, method))
        });
        assert_eq!(found, Some((1, p, Method::Pm1Stage2)));
        assert_eq!(progress.pm.unwrap().b1(), LEVELS[1].b1 * PM1_B1_RATIO);
    }

    #[test]
    fn deterministic() {
        let n = Integer::from_str("7060005655815754299976961394452809").unwrap();
        assert_eq!(factor(&n, 7), factor(&n, 7));
    }
}
