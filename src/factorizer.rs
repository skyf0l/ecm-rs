//! [`Factorizer`]: the configurable entry point, with events and cancellation.

use crate::{
    driver::{Engine, Mode},
    ecm::{Error, Param, rand_state},
    events::{Event, EventHandler, NoEvents},
    stage2::MAX_POLY_MEMORY,
    stop::Stop,
};
use rug::Integer;
use std::{
    collections::HashMap,
    fmt,
    ops::ControlFlow,
    sync::{Arc, atomic::AtomicBool},
    time::{Duration, Instant},
};

/// Seed of [`crate::ecm()`], and of a [`Factorizer`] by default.
pub(crate) const DEFAULT_SEED: u64 = 1234;

/// The factoring method of a [`Factorizer`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub enum Algorithm {
    /// Elliptic curves, after Pollard's P-1 method with larger bounds (except with fixed
    /// bounds).
    #[default]
    Ecm,
    /// Pollard's P-1 method only: finds the factors `p` with a smooth `p - 1`.
    Pm1,
    /// Williams' P+1 method only: finds the factors `p` with a smooth `p + 1` (or `p - 1`,
    /// depending on the seed and `p`, see [`Factorizer::x0`]).
    ///
    /// It is not part of [`Algorithm::Ecm`]: after P-1, it finds too few factors for its cost
    /// (the curves find random factors of 15 to 30 digits as fast without it).
    Pp1,
}

impl fmt::Display for Algorithm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Ecm => "ECM",
            Self::Pm1 => "P-1",
            Self::Pp1 => "P+1",
        })
    }
}

/// Factorization with options: seed, fixed bounds, curves, P-1 or P+1 only, stage 2 memory, a
/// callback receiving the [`Event`]s of the factorization, which can interrupt it, and an
/// interruption flag or a timeout, which interrupt it even during a curve.
///
/// By default, the factors are found from the smallest to the largest, as with [`crate::ecm()`]
/// (which is `Factorizer::new().factor(n)`). With fixed bounds ([`Factorizer::b1`]), each
/// composite part of the number runs curves with these bounds instead, as with
/// [`crate::ecm_with_params`].
///
/// The results only depend on the number and the options, not on the callback: the same seed
/// gives the same factors, found by the same curves.
///
/// ```
/// use ecm::{Event, Factorizer};
/// use rug::Integer;
/// use std::ops::ControlFlow;
///
/// let n: Integer = "398883434337287".parse().unwrap();
/// let mut curves = 0;
/// let factors = Factorizer::new()
///     .seed(42)
///     .on_event(|event| {
///         if let Event::Curve { .. } = event {
///             curves += 1;
///         }
///         ControlFlow::Continue(())
///     })
///     .factor(&n)
///     .unwrap();
/// assert_eq!(factors.len(), 2);
/// ```
#[derive(Debug, Clone)]
pub struct Factorizer<H = NoEvents> {
    seed: u64,
    param: Param,
    b1: Option<usize>,
    b2: Option<usize>,
    curves: Option<usize>,
    sigma: Option<Integer>,
    algorithm: Algorithm,
    x0: Option<(Integer, Integer)>,
    max_memory: usize,
    interrupt: Option<Arc<AtomicBool>>,
    timeout: Option<Duration>,
    handler: H,
}

impl Default for Factorizer {
    fn default() -> Self {
        Self::new()
    }
}

impl Factorizer {
    /// Default options: those of [`crate::ecm()`].
    #[must_use]
    pub fn new() -> Self {
        Self {
            seed: DEFAULT_SEED,
            param: Param::default(),
            b1: None,
            b2: None,
            curves: None,
            sigma: None,
            algorithm: Algorithm::Ecm,
            x0: None,
            max_memory: MAX_POLY_MEMORY,
            interrupt: None,
            timeout: None,
            handler: NoEvents,
        }
    }
}

impl<H: EventHandler> Factorizer<H> {
    /// Seed of the random curves (default: the seed of [`crate::ecm()`]).
    #[must_use]
    pub const fn seed(mut self, seed: u64) -> Self {
        self.seed = seed;
        self
    }

    /// Parametrization of the curves (default: [`Param::Batch2`]).
    #[must_use]
    pub const fn param(mut self, param: Param) -> Self {
        self.param = param;
        self
    }

    /// Fixed bounds: each composite part runs curves with the stage 1 bound `b1` (even, at
    /// least 6) instead of the bounds for factors of increasing size.
    #[must_use]
    pub const fn b1(mut self, b1: usize) -> Self {
        self.b1 = Some(b1);
        self
    }

    /// Stage 2 bound with fixed bounds (even, at least 4; no stage 2 if at most `b1`). By
    /// default, close to GMP-ECM's default for `b1`. Requires [`Factorizer::b1`].
    #[must_use]
    pub const fn b2(mut self, b2: usize) -> Self {
        self.b2 = Some(b2);
        self
    }

    /// Largest number of curves per composite part with fixed bounds (by default, unlimited:
    /// until a factor is found, or the callback interrupts). Requires [`Factorizer::b1`].
    #[must_use]
    pub const fn curves(mut self, curves: usize) -> Self {
        self.curves = Some(curves);
        self
    }

    /// Parameter of the first curve with fixed bounds, instead of a random one (the next curves
    /// take `sigma + 1`, `sigma + 2`, ...), in the range of the parametrization (see
    /// [`Param`]; for [`Param::Suyama`], at least 6, and taken modulo the number). Requires
    /// [`Factorizer::b1`].
    #[must_use]
    pub fn sigma(mut self, sigma: Integer) -> Self {
        self.sigma = Some(sigma);
        self
    }

    /// The factoring method (default: [`Algorithm::Ecm`]). With [`Algorithm::Pm1`] or
    /// [`Algorithm::Pp1`], no curves: with fixed bounds, the method runs once with `b1` and
    /// `b2` per composite part; otherwise it runs with the bounds of P-1 at each level.
    #[must_use]
    pub const fn algorithm(mut self, algorithm: Algorithm) -> Self {
        self.algorithm = algorithm;
        self
    }

    /// Starting value `x0 = numerator/denominator` (modulo the number) of P-1 or P+1, as
    /// GMP-ECM's `-x0` (default: 3 for P-1, `2/7` for P+1). Requires [`Algorithm::Pm1`] or
    /// [`Algorithm::Pp1`].
    ///
    /// P+1 works modulo a prime `p` in a group of order `p + 1` if `x0^2 - 4` is not a square
    /// modulo `p`, else of order `p - 1`: about half of the seeds find a factor with a smooth
    /// `p + 1`. `2/7` (orders multiple of 6: `p + 1` for `p = 2 mod 3`) and `6/5` (multiple of
    /// 4: `p + 1` for `p = 3 mod 4`) do slightly better than random seeds.
    #[must_use]
    pub fn x0(mut self, numerator: Integer, denominator: Integer) -> Self {
        self.x0 = Some((numerator, denominator));
        self
    }

    /// Memory (in bytes) the polynomial stage 2 may use (default: 256 MiB): a smaller limit
    /// makes it slower. The baby-step giant-step stage 2, chosen when no polynomial one fits,
    /// uses at most about 32 MiB plus its baby steps.
    #[must_use]
    pub const fn max_memory(mut self, bytes: usize) -> Self {
        self.max_memory = bytes;
        self
    }

    /// Interrupts the factorizations when `flag` is set (from another thread, or a signal
    /// handler): they return [`Error::Interrupted`] soon after, even in the middle of a curve
    /// (see [`Factorizer::factor_partial`] for what was found so far). The flag is checked
    /// every few milliseconds, but for the largest polynomial products of stage 2, which take
    /// up to a few tenths of a second for a 1024-bit number with the bounds of 35-digit factors.
    /// The flag is never reset: a factorization started while it is set is interrupted at
    /// once.
    #[must_use]
    pub fn interrupt_flag(mut self, flag: Arc<AtomicBool>) -> Self {
        self.interrupt = Some(flag);
        self
    }

    /// Interrupts each factorization after `timeout` (from the call to [`Factorizer::factor`],
    /// [`Factorizer::factor_partial`] or [`Factorizer::find_factor`]), as
    /// [`Factorizer::interrupt_flag`] does.
    #[must_use]
    pub const fn timeout(mut self, timeout: Duration) -> Self {
        self.timeout = Some(timeout);
        self
    }

    /// Calls `f` with each [`Event`] of the factorizations: returning [`ControlFlow::Break`]
    /// interrupts the factorization, which returns [`Error::Interrupted`] (see
    /// [`Factorizer::factor_partial`] for what was found so far). The events come between
    /// curves and stages: an interruption takes effect at the next event, at most a curve or
    /// a stage of P-1 later (use [`Factorizer::interrupt_flag`] or [`Factorizer::timeout`] to
    /// interrupt the factorization during a curve). No event comes after an interruption.
    ///
    /// Without it, the events are not even built.
    #[must_use]
    pub fn on_event<F: FnMut(&Event<'_>) -> ControlFlow<()>>(self, f: F) -> Factorizer<F> {
        Factorizer {
            seed: self.seed,
            param: self.param,
            b1: self.b1,
            b2: self.b2,
            curves: self.curves,
            sigma: self.sigma,
            algorithm: self.algorithm,
            x0: self.x0,
            max_memory: self.max_memory,
            interrupt: self.interrupt,
            timeout: self.timeout,
            handler: f,
        }
    }

    /// The complete factorization of `n`: its prime factors with their multiplicity (empty for
    /// `n = 1`).
    ///
    /// # Errors
    ///
    /// [`Error::BoundsNotEven`], [`Error::BoundsTooSmall`] and [`Error::InvalidOption`] for
    /// invalid options, [`Error::Interrupted`] if the callback, the interruption flag or the
    /// timeout interrupts the factorization, and [`Error::ECMFailed`] if a composite part is not split (with fixed bounds, after the
    /// maximum number of curves; otherwise, as [`crate::ecm()`]).
    ///
    /// # Panics
    ///
    /// If `n` is not positive.
    pub fn factor(&mut self, n: &Integer) -> Result<HashMap<Integer, usize>, Error> {
        self.factor_partial(n).into_result()
    }

    /// As [`Factorizer::factor`], but on failure or interruption, also returns what was found
    /// so far: the prime factors, and the parts still to factor.
    ///
    /// # Panics
    ///
    /// If `n` is not positive (and the options are valid).
    pub fn factor_partial(&mut self, n: &Integer) -> Factorization {
        let mode = match self.mode() {
            Ok(mode) => mode,
            Err(error) => {
                return Factorization {
                    primes: HashMap::new(),
                    unfactored: vec![(n.clone(), 1)],
                    error: Some(error),
                };
            }
        };
        let mut rand = rand_state(self.seed);
        self.engine(mode, &mut rand).factor(n)
    }

    /// A proper factor of `n` (maybe composite): by trial division, then as
    /// [`Factorizer::factor`] searches, stopping at the first factor found.
    ///
    /// # Errors
    ///
    /// [`Error::NumberIsPrime`] if `n` is prime, and otherwise as [`Factorizer::factor`].
    ///
    /// # Panics
    ///
    /// If `n <= 1`: it has no proper factor.
    pub fn find_factor(&mut self, n: &Integer) -> Result<Integer, Error> {
        let mode = self.mode()?;
        let mut rand = rand_state(self.seed);
        self.engine(mode, &mut rand).find_one(n, true)
    }

    fn engine<'a, 'r>(
        &'a mut self,
        mode: Mode,
        rand: &'a mut rug::rand::RandState<'r>,
    ) -> Engine<'a, 'r, H> {
        let deadline = self
            .timeout
            .and_then(|timeout| Instant::now().checked_add(timeout));
        let stop = Stop::new(self.interrupt.as_deref(), deadline);
        Engine::new(
            mode,
            self.param,
            self.sigma.clone(),
            (self.algorithm, self.x0.clone()),
            self.max_memory,
            rand,
            &mut self.handler,
            stop,
        )
    }

    /// Checks the options.
    fn mode(&self) -> Result<Mode, Error> {
        if let Some((_, denominator)) = &self.x0 {
            if self.algorithm == Algorithm::Ecm {
                return Err(Error::InvalidOption("x0 is for P-1 and P+1, not curves"));
            }
            if *denominator == 0 {
                return Err(Error::InvalidOption("x0 with a zero denominator"));
            }
        }
        let Some(b1) = self.b1 else {
            if self.b2.is_some() || self.curves.is_some() || self.sigma.is_some() {
                return Err(Error::InvalidOption("b2, curves and sigma require b1"));
            }
            return Ok(Mode::Levels);
        };
        if let Some(sigma) = &self.sigma {
            if self.algorithm != Algorithm::Ecm {
                return Err(Error::InvalidOption("sigma is for curves, not P-1 or P+1"));
            }
            let valid = match self.param {
                Param::Suyama => *sigma >= 6,
                Param::Square => *sigma >= 2 && *sigma < 1u64 << 32,
                Param::Batch2 => *sigma >= 2 && sigma.significant_bits() <= 64,
            };
            if !valid {
                return Err(Error::InvalidOption(
                    "sigma out of the range of the parametrization",
                ));
            }
        }
        let b2 = self.b2.unwrap_or_else(|| crate::driver::default_b2(b1));
        Mode::fixed(b1, b2, self.curves)
    }
}

/// Result of [`Factorizer::factor_partial`]: the factors found, complete unless there is an
/// error.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Factorization {
    /// Prime factors found, with their multiplicity.
    pub primes: HashMap<Integer, usize>,
    /// Parts of the number not factored (composite, or not tested yet), with their
    /// multiplicity: the number is the product of these and of the primes.
    pub unfactored: Vec<(Integer, usize)>,
    /// Why the factorization is incomplete, if it is.
    pub error: Option<Error>,
}

impl Factorization {
    /// The prime factors if the factorization is complete, the error otherwise.
    ///
    /// # Errors
    ///
    /// [`Factorization::error`], if any.
    pub fn into_result(self) -> Result<HashMap<Integer, usize>, Error> {
        match self.error {
            None => Ok(self.primes),
            Some(error) => Err(error),
        }
    }
}
