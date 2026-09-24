//! Interruption of the long computations (curves, P-1) from inside, see
//! [`crate::Factorizer::interrupt_flag`] and [`crate::Factorizer::timeout`].

use std::{
    sync::atomic::{AtomicBool, Ordering},
    time::Instant,
};

/// Ladder steps (or equivalent work) between two checks of a [`Stop`]: a few milliseconds at
/// 1024 bits, while a check (an atomic load, and reading the clock with a deadline) costs less
/// than a ladder step even at 64 bits.
pub(crate) const STOP_INTERVAL: u32 = 1 << 10;

/// When to stop a computation: when a flag is set, or at a deadline. The default never stops.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Stop<'a> {
    flag: Option<&'a AtomicBool>,
    deadline: Option<Instant>,
}

impl<'a> Stop<'a> {
    /// Never stops.
    pub(crate) const NEVER: Stop<'static> = Stop {
        flag: None,
        deadline: None,
    };

    pub(crate) fn new(flag: Option<&'a AtomicBool>, deadline: Option<Instant>) -> Self {
        Self { flag, deadline }
    }

    /// Whether it never stops.
    pub(crate) fn is_never(self) -> bool {
        self.flag.is_none() && self.deadline.is_none()
    }

    /// Whether the computation must stop now.
    #[inline]
    pub(crate) fn requested(self) -> bool {
        self.flag.is_some_and(|flag| flag.load(Ordering::Relaxed))
            || self
                .deadline
                .is_some_and(|deadline| Instant::now() >= deadline)
    }
}
