//! Events of a factorization, reported to the callback of [`crate::Factorizer::on_event`].

use crate::ecm::Param;
use rug::Integer;
use std::{collections::HashMap, fmt, ops::ControlFlow, time::Duration};

/// Something that happened during a factorization, see [`crate::Factorizer::on_event`].
///
/// The data is borrowed from the factorization: copy what must outlive the callback.
///
/// Order of the events of [`crate::Factorizer::factor`]: [`Event::TrialDivision`] first, with
/// the [`Event::Prime`] events of its factors, then for each composite part, the searches
/// ([`Event::Pm1`] or [`Event::Pp1`], [`Event::Level`] followed by its [`Event::Curve`] events),
/// until a
/// [`Event::Factor`] splits it, followed by the [`Event::Prime`] (and [`Event::Factor`] for a
/// perfect power) events of the parts. The [`Event::Prime`] events together are the complete
/// factorization: the product of `p^exponent` over them is the number.
#[non_exhaustive]
#[derive(Debug, Clone, Copy)]
pub enum Event<'a> {
    /// Trial division removed the prime factors below 2^16 (each is also reported by an
    /// [`Event::Prime`] event, right after this one).
    #[non_exhaustive]
    TrialDivision {
        /// The prime factors found, with their multiplicity.
        factors: &'a HashMap<Integer, usize>,
        /// What remains of the number (`1` if trial division factored it completely).
        cofactor: &'a Integer,
    },
    /// The searches on the composite `n` compute modulo `2^k + 1` (if `k > 0`) or `2^-k - 1`
    /// (if `k < 0`), a multiple of `n`, with a special reduction (see [`crate::Base2Mode`]):
    /// before the other events of `n`.
    #[non_exhaustive]
    Base2 {
        /// The composite number searched.
        n: &'a Integer,
        /// The exponent, signed as GMP-ECM's `-base2`.
        k: i64,
    },
    /// Pollard's P-1 method ran on `n`.
    #[non_exhaustive]
    Pm1 {
        /// The composite number searched.
        n: &'a Integer,
        /// Stage 1 bound (stage 1 resumes from the previous bound on the same number).
        b1: usize,
        /// Stage 2 bound (every prime up to it is covered), `None` if stage 2 did not run
        /// (stage 1 found a factor, or all the factors at once).
        b2: Option<usize>,
        /// Duration of stage 1.
        stage1: Duration,
        /// Duration of stage 2.
        stage2: Duration,
    },
    /// Williams' P+1 method ran on `n` (see [`crate::Algorithm::Pp1`]).
    #[non_exhaustive]
    Pp1 {
        /// The composite number searched.
        n: &'a Integer,
        /// Stage 1 bound (stage 1 resumes from the previous bound on the same number).
        b1: usize,
        /// Stage 2 bound (every prime up to it is covered), `None` if stage 2 did not run
        /// (stage 1 found a factor, or all the factors at once).
        b2: Option<usize>,
        /// Duration of stage 1.
        stage1: Duration,
        /// Duration of stage 2.
        stage2: Duration,
    },
    /// Curves with the bounds `b1` and `b2` start (or resume) on `n`.
    #[non_exhaustive]
    Level {
        /// The composite number searched.
        n: &'a Integer,
        /// Size in digits of the factors this level targets, `None` with fixed bounds.
        digits: Option<u32>,
        /// Stage 1 bound.
        b1: usize,
        /// Stage 2 bound actually covered (at least the one asked for).
        b2: usize,
        /// Number of curves of the level: the expected number to find a factor of `digits`
        /// digits, or with fixed bounds the maximum number of curves (`None` if unlimited).
        curves: Option<usize>,
        /// Curves of the level already run (on a multiple of `n`, before it was split).
        done: usize,
    },
    /// A curve was run, without finding a factor or before the [`Event::Factor`] event of the
    /// factor it found.
    #[non_exhaustive]
    Curve {
        /// The composite number searched.
        n: &'a Integer,
        /// Parametrization of the curve.
        param: Param,
        /// Parameter of the curve: GMP-ECM finds the same factor with `-param {param} -sigma
        /// {sigma}` and the same bounds.
        sigma: &'a Integer,
        /// Index of the curve in its level, from 1 (see [`Event::Level`]).
        index: usize,
        /// Duration of the curve setup and stage 1.
        stage1: Duration,
        /// Duration of stage 2 (zero if it did not run).
        stage2: Duration,
    },
    /// `n` was split: `factor` is a proper factor of `n` (maybe composite).
    #[non_exhaustive]
    Factor {
        /// The number split.
        n: &'a Integer,
        /// A proper factor of `n`.
        factor: &'a Integer,
        /// How it was found.
        method: Method,
    },
    /// `p^exponent` divides the number: `p` is a prime factor, with this contribution to its
    /// multiplicity (the same prime may be reported more than once).
    #[non_exhaustive]
    Prime {
        /// A prime factor (probably prime: BPSW, no counterexample is known).
        p: &'a Integer,
        /// Multiplicity found with this event.
        exponent: usize,
    },
}

/// How a factor was found.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Method {
    /// Trial division by the primes below 2^16.
    TrialDivision,
    /// The number is a perfect power, split by its root.
    PerfectPower,
    /// Stage 1 of Pollard's P-1 method.
    Pm1Stage1,
    /// Stage 2 of Pollard's P-1 method.
    Pm1Stage2,
    /// Stage 1 of Williams' P+1 method (or the seed, not defined modulo a factor).
    Pp1Stage1,
    /// Stage 2 of Williams' P+1 method.
    Pp1Stage2,
    /// Building a curve (a failed modular inversion).
    EcmSetup,
    /// Stage 1 of a curve.
    EcmStage1,
    /// Stage 2 of a curve.
    EcmStage2,
}

impl fmt::Display for Method {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::TrialDivision => "trial division",
            Self::PerfectPower => "perfect power",
            Self::Pm1Stage1 => "P-1 stage 1",
            Self::Pm1Stage2 => "P-1 stage 2",
            Self::Pp1Stage1 => "P+1 stage 1",
            Self::Pp1Stage2 => "P+1 stage 2",
            Self::EcmSetup => "ECM curve setup",
            Self::EcmStage1 => "ECM stage 1",
            Self::EcmStage2 => "ECM stage 2",
        })
    }
}

/// Receives the [`Event`]s of a factorization: implemented by the closures
/// `FnMut(&Event<'_>) -> ControlFlow<()>` (see [`crate::Factorizer::on_event`]).
pub trait EventHandler {
    /// Whether the events are wanted: when `false`, the factorization neither builds them nor
    /// measures the durations they report.
    const ENABLED: bool = true;

    /// Handles `event`: [`ControlFlow::Break`] interrupts the factorization (see
    /// [`crate::Error::Interrupted`]).
    fn handle(&mut self, event: &Event<'_>) -> ControlFlow<()>;
}

impl<F: FnMut(&Event<'_>) -> ControlFlow<()>> EventHandler for F {
    fn handle(&mut self, event: &Event<'_>) -> ControlFlow<()> {
        self(event)
    }
}

/// The default [`EventHandler`]: ignores the events, at no cost.
#[derive(Debug, Clone, Copy, Default)]
pub struct NoEvents;

impl EventHandler for NoEvents {
    const ENABLED: bool = false;

    fn handle(&mut self, _: &Event<'_>) -> ControlFlow<()> {
        ControlFlow::Continue(())
    }
}
