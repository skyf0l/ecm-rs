//! Prime generation for the whole crate.
//!
//! `primal::Primes::all()` sieves a large first segment when it is created (about 0.2 ms), which
//! is more than factoring a small number takes: for small bounds, a sieve of exactly the needed
//! range is much cheaper (23 us up to 2^16).

use primal::{Primes, Sieve};
use std::sync::{Arc, Mutex, PoisonError};

/// Up to this bound, the primes come from a sieve of exactly the range, else from
/// `primal::Primes::all()` (streamed, in constant memory).
const SMALL_SIEVE: usize = 1 << 22;

/// The primes up to the largest bound `<= SMALL_SIEVE` asked for so far, and this bound: the
/// stage 2 plans and the probabilities of the levels of [`crate::ecm()`] all need the primes
/// up to (different) bounds, and collecting them costs more than the plans themselves for small
/// numbers. At most 1.2 MB (the 295947 primes below `SMALL_SIEVE`).
static SIEVED: Mutex<Option<(usize, Arc<[u32]>)>> = Mutex::new(None);

/// The primes `<= hi`, in increasing order.
pub(crate) fn primes(hi: usize) -> PrimesUpTo {
    if hi <= SMALL_SIEVE {
        let primes = sieved(hi);
        let end = primes.partition_point(|&p| p as usize <= hi);
        PrimesUpTo::Small(primes, 0..end)
    } else {
        PrimesUpTo::Large(Box::new(Primes::all()), hi)
    }
}

/// The primes up to at least `hi <= SMALL_SIEVE` (from [`SIEVED`], or sieved now).
fn sieved(hi: usize) -> Arc<[u32]> {
    let mut cache = SIEVED.lock().unwrap_or_else(PoisonError::into_inner);
    if let Some((bound, primes)) = &*cache {
        if *bound >= hi {
            return Arc::clone(primes);
        }
    }
    // pi(x) < 1.26*x/ln(x).
    let x = hi.max(2) as f64;
    let mut list = Vec::with_capacity((1.26 * x / x.ln()) as usize + 1);
    list.extend(
        Sieve::new(hi)
            .primes_from(0)
            .take_while(|&p| p <= hi)
            .map(|p| p as u32),
    );
    let primes: Arc<[u32]> = list.into();
    *cache = Some((hi, Arc::clone(&primes)));
    primes
}

/// Iterator of [`primes`].
pub(crate) enum PrimesUpTo {
    /// Sieved at once: the primes at these indices.
    Small(Arc<[u32]>, std::ops::Range<usize>),
    /// Streamed, up to the bound.
    Large(Box<Primes>, usize),
}

impl Iterator for PrimesUpTo {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<usize> {
        match self {
            Self::Small(primes, range) => range.next().map(|i| primes[i] as usize),
            Self::Large(primes, hi) => primes.next().filter(|p| p <= hi),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_primal() {
        for hi in [
            0,
            1,
            2,
            3,
            4,
            5,
            30,
            31,
            65535,
            65536,
            65537,
            1 << 20,
            SMALL_SIEVE,
        ] {
            let expected: Vec<usize> = Primes::all().take_while(|&p| p <= hi).collect();
            assert_eq!(primes(hi).collect::<Vec<_>>(), expected, "{hi}");
        }
        // Streamed: the same primes near the threshold.
        let hi = SMALL_SIEVE + 1000;
        let tail: Vec<usize> = primes(hi).skip_while(|&p| p < SMALL_SIEVE - 1000).collect();
        let expected: Vec<usize> = Primes::all()
            .skip_while(|&p| p < SMALL_SIEVE - 1000)
            .take_while(|&p| p <= hi)
            .collect();
        assert_eq!(tail, expected);
    }
}
