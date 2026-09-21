//! Deterministic inputs shared by the benchmarks and the `success_rate` example.

#![allow(dead_code)]

use rug::{rand::RandState, Integer};

/// Seed used for every generated input, so all runs measure exactly the same work.
pub const SEED: u64 = 0xEC0;

/// Fixed Suyama parameter for single-curve benchmarks.
pub const SIGMA: u64 = 1_234_567;

/// Stage bounds recommended by GMP-ECM for a factor of `digits` decimal digits:
/// `(digits, B1, B2)`, with `B2` rounded to an even number.
pub const GMP_ECM_BOUNDS: [(u32, usize, usize); 3] = [
    (15, 2_000, 147_396),
    (20, 11_000, 1_873_422),
    (25, 50_000, 12_746_592),
];

/// Random prime of exactly `bits` bits, derived from `seed`.
pub fn prime_bits(bits: u32, seed: u64) -> Integer {
    let mut rand = RandState::new();
    rand.seed(&Integer::from(seed));
    loop {
        let mut x = Integer::from(Integer::random_bits(bits, &mut rand));
        x.set_bit(bits - 1, true);
        let p = x.next_prime();
        if p.significant_bits() == bits {
            return p;
        }
    }
}

/// Random prime of exactly `digits` decimal digits, derived from `seed`.
pub fn prime_digits(digits: u32, seed: u64) -> Integer {
    let mut rand = RandState::new();
    rand.seed(&Integer::from(seed));
    let low = Integer::from(Integer::u_pow_u(10, digits - 1));
    let high = Integer::from(Integer::u_pow_u(10, digits));
    loop {
        let range = Integer::from(&high - &low);
        let p = (low.clone() + range.random_below(&mut rand)).next_prime();
        if p < high {
            return p;
        }
    }
}

/// Product of two random primes of `bits / 2` bits: no small factors, so every
/// operation runs on a full-size modulus.
pub fn semiprime_bits(bits: u32, seed: u64) -> Integer {
    prime_bits(bits / 2, seed) * prime_bits(bits - bits / 2, seed.wrapping_add(1))
}
