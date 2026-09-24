//! Prime generation for the whole crate.
//!
//! `primal::Primes::all()` sieves a large first segment when it is created (about 0.2 ms), which
//! is more than factoring a small number takes: for small bounds, a sieve of exactly the needed
//! range is much cheaper (23 us up to 2^16).

use primal::{Primes, Sieve};

/// Up to this bound, the primes come from a sieve of exactly the range, else from
/// `primal::Primes::all()` (streamed, in constant memory).
const SMALL_SIEVE: usize = 1 << 22;

/// The primes `<= hi`, in increasing order.
pub(crate) fn primes(hi: usize) -> PrimesUpTo {
    if hi <= SMALL_SIEVE {
        let primes: Vec<u32> = Sieve::new(hi)
            .primes_from(0)
            .take_while(|&p| p <= hi)
            .map(|p| p as u32)
            .collect();
        PrimesUpTo::Small(primes.into_iter())
    } else {
        PrimesUpTo::Large(Box::new(Primes::all()), hi)
    }
}

/// Iterator of [`primes`].
pub(crate) enum PrimesUpTo {
    /// Sieved at once.
    Small(std::vec::IntoIter<u32>),
    /// Streamed, up to the bound.
    Large(Box<Primes>, usize),
}

impl Iterator for PrimesUpTo {
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<usize> {
        match self {
            Self::Small(primes) => primes.next().map(|p| p as usize),
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
