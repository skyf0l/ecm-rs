//! Pollard's P-1 method, the essentials of GMP-ECM's `pm1.c`.
//!
//! Stage 1 computes `x = x0^E mod n` with `E` the product of the largest powers of the primes
//! `<= b1` that are `<= b1` (the same multiplier as ECM's stage 1): a prime `p` of `n` with
//! `b1`-smooth `p - 1` divides `x - 1`. Stage 1 can be resumed with a larger `b1`, and `x` stays
//! valid modulo any divisor of `n`.
//!
//! Stage 2 checks the primes `l` in `(b1, b2]` with the same continuations as ECM (see
//! [`crate::stage2`]), on the Lucas sequence `V_k = x^k + x^-k` (see [`crate::lucas`]), with
//! `x^(m*D) * (V_(m*D) - V_j) = (x^(m*D - j) - 1) * (x^(m*D + j) - 1)`, so `V_(m*D) = V_j`
//! modulo `p` checks `m*D +- j` as `x(m*D*Q) = x(j*Q)` does for a curve. The identity element
//! is `V_k = 2`, that is `x^k = 1`.

use crate::{
    ecm::{prime_power_words, product},
    lucas,
    stage2::Stage2Plan,
    stop::Stop,
};
use rug::Integer;

/// Starting value `x0` of stage 1 (the choice barely matters: `E` has a large power of 2).
const X0: u32 = 3;

/// 64-bit factors of the exponent multiplied together before each modular exponentiation
/// (about 2^20 bits): bounds the memory used.
const CHUNK_WORDS: usize = 1 << 14;

/// [`CHUNK_WORDS`] when stage 1 may have to stop (see [`Stop`]), divided by the number of limbs
/// of the number: an exponentiation takes a few milliseconds. GMP then uses smaller windows,
/// which costs about 2% more instructions at 1024 bits.
const STOP_CHUNK_WORDS: usize = 1 << 13;

/// State of P-1 on a number: stage 1 done up to `b1`.
#[derive(Debug, Clone)]
pub struct Pm1 {
    /// `x0^E(b1)` modulo the number (or a multiple of it).
    x: Integer,
    b1: usize,
}

impl Default for Pm1 {
    fn default() -> Self {
        Self {
            x: Integer::from(X0),
            b1: 1,
        }
    }
}

impl Pm1 {
    /// Nothing done yet.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Stage 1 bound reached so far.
    #[must_use]
    pub fn b1(&self) -> usize {
        self.b1
    }

    /// Extends stage 1 to `b1` modulo `n` (a divisor of the previous numbers, if any): returns
    /// `gcd(x - 1, n)`.
    #[cfg_attr(not(any(test, feature = "bench")), allow(dead_code))]
    pub fn stage1(&mut self, n: &Integer, b1: usize) -> Integer {
        self.stage1_until(n, b1, Stop::NEVER)
    }

    /// [`Pm1::stage1`], stopping early if `stop` is requested: the bound reached is then not
    /// updated, and `x` stays valid (a power of the previous one).
    pub(crate) fn stage1_until(&mut self, n: &Integer, b1: usize, stop: Stop<'_>) -> Integer {
        self.x %= n;
        if self.b1 < b1 {
            let chunk_words = if stop.is_never() {
                CHUNK_WORDS
            } else {
                (STOP_CHUNK_WORDS / n.significant_digits::<u64>().max(1)).max(16)
            };
            // The exponent by chunks of chunk_words words.
            let mut words = prime_power_words(self.b1, b1).peekable();
            let mut stopped = false;
            while words.peek().is_some() {
                if stop.requested() {
                    stopped = true;
                    break;
                }
                let chunk: Vec<Integer> = words
                    .by_ref()
                    .take(chunk_words)
                    .map(Integer::from)
                    .collect();
                self.x
                    .pow_mod_mut(&product(chunk), n)
                    .expect("positive exponent");
            }
            if !stopped {
                self.b1 = b1;
            }
        }
        Integer::from(&self.x - 1u32).gcd(n)
    }

    /// Stage 2 modulo `n` with `plan`, whose `b1` must be at most the stage 1 bound: returns
    /// `gcd(g, n)` (see [`crate::ecm::stage2`]).
    #[cfg_attr(not(any(test, feature = "bench")), allow(dead_code))]
    #[must_use]
    pub fn stage2(&self, n: &Integer, plan: &Stage2Plan) -> Integer {
        self.stage2_until(n, plan, Stop::NEVER)
    }

    /// [`Pm1::stage2`], stopping early (with the gcd of a partial product) if `stop` is
    /// requested.
    pub(crate) fn stage2_until(&self, n: &Integer, plan: &Stage2Plan, stop: Stop<'_>) -> Integer {
        let x = Integer::from(&self.x % n);
        let Ok(inv) = x.clone().invert(n) else {
            return x.gcd(n);
        };
        lucas::stage2(n, &((inv + x) % n), plan, stop)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use primal::Primes;
    use rug::{ops::Pow, rand::RandState};

    #[test]
    fn stage1_resumes() {
        let n = Integer::from(1_000_000_007u64) * Integer::from(998_244_353u64);
        let mut direct = Pm1::new();
        direct.stage1(&n, 100_000);
        let mut resumed = Pm1::new();
        for b1 in [10, 11, 100, 2_000, 99_999, 100_000] {
            resumed.stage1(&n, b1);
        }
        assert_eq!(resumed.b1(), 100_000);
        assert_eq!(resumed.x, direct.x);
        // Larger than a chunk of the exponent.
        let b1 = 2_000_000;
        let mut pm1 = Pm1::new();
        pm1.stage1(&n, b1);
        let e = crate::ecm::stage1_multiplier(b1);
        assert_eq!(pm1.x, Integer::from(X0).pow_mod(&e, &n).unwrap());
    }

    #[test]
    fn stage1_finds_smooth() {
        // 1000000007 - 1 = 2 * 500000003, 998244353 - 1 = 2^23 * 7 * 17.
        let n = Integer::from(1_000_000_007u64) * Integer::from(998_244_353u64);
        assert_eq!(Pm1::new().stage1(&n, 10), 1);
        assert_eq!(Pm1::new().stage1(&n, 1 << 23), 998_244_353);
    }

    /// Stage 2 must find `p` whenever `x^l = 1` modulo `p` for some prime `b1 < l <= b2` after
    /// stage 1 (`x != 1`), for each small prime `p` times the large prime `q`.
    ///
    /// Returns the number of primes `p` where stage 2 had a factor to find.
    fn check_stage2_primes(q: &Integer, b1: usize, b2: usize, primes: &[u64]) -> usize {
        let ls: Vec<u64> = Primes::all()
            .skip_while(|&l| l <= b1)
            .take_while(|&l| l <= b2)
            .map(|l| l as u64)
            .collect();
        let mut checked = 0;
        for &p in primes {
            let n = Integer::from(q * p);
            let plans = [Stage2Plan::pairs(b1, b2), Stage2Plan::poly(&n, b1, b2)];
            let mut pm1 = Pm1::new();
            if pm1.stage1(&n, b1) != 1 {
                continue;
            }
            let y = Integer::from(&pm1.x % p);
            let expected = ls.iter().any(|&l| {
                y.clone()
                    .pow_mod(&Integer::from(l), &Integer::from(p))
                    .unwrap()
                    == 1
            });
            for plan in &plans {
                let g = pm1.stage2(&n, plan);
                assert!(n.is_divisible(&g));
                if expected {
                    assert_eq!(g, p, "stage 2 missed p = {p}: {plan:?} {b1} {b2}");
                }
            }
            checked += usize::from(expected);
        }
        checked
    }

    #[test]
    fn stage2_checks_all_primes() {
        let q = Integer::from(10).pow(30) + 57u32;
        assert_ne!(q.is_probably_prime(30), rug::integer::IsPrime::No);
        let primes: Vec<u64> = Primes::all()
            .skip_while(|&p| p < 100_000)
            .take(800)
            .map(|p| p as u64)
            .collect();
        let mut checked = 0;
        for (b1, b2) in [
            (30, 1000),
            (100, 5000),
            (100, 10_001),
            (200, 20_000),
            (1000, 999),
        ] {
            checked += check_stage2_primes(&q, b1, b2, &primes);
        }
        assert!(checked > 100, "{checked}");
    }

    #[test]
    fn stage2_checks_all_primes_limbs() {
        // Every limb count of `Mont`, and `Plain`.
        let mut rand = RandState::new();
        let primes: Vec<u64> = Primes::all()
            .skip_while(|&p| p < 10_000)
            .take(100)
            .map(|p| p as u64)
            .collect();
        for bits in [64, 100, 192, 320, 512, 700, 1000, 1100] {
            let mut q = Integer::from(Integer::random_bits(bits, &mut rand));
            q.set_bit(bits - 1, true);
            let checked = check_stage2_primes(&q.next_prime(), 50, 3000, &primes);
            assert!(checked > 5, "{bits} bits: {checked}");
        }
    }
}
