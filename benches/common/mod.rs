//! Deterministic inputs shared by the benchmarks and the examples.

#![allow(dead_code)]

use ecm::bench::{curve, stage1, stage1_multiplier, Param, Point};
use rug::{rand::RandState, Integer};

/// Seed used for every generated input, so all runs measure exactly the same work.
pub const SEED: u64 = 0xEC0;

/// Fixed curve parameter `sigma` for single-curve benchmarks (valid for every `Param`).
pub const SIGMA: u64 = 1_234_567;

/// Stage bounds recommended by GMP-ECM for a factor of `digits` decimal digits:
/// `(digits, B1, B2)`, with `B2` rounded to an even number.
pub const GMP_ECM_BOUNDS: [(u32, usize, usize); 5] = [
    (15, 2_000, 147_396),
    (20, 11_000, 1_873_422),
    (25, 50_000, 12_746_592),
    (30, 250_000, 128_992_510),
    (35, 1_000_000, 1_045_563_762),
];

/// `(B1, B2)` of [`GMP_ECM_BOUNDS`] for a factor of `digits` digits.
pub fn bounds(digits: u32) -> (usize, usize) {
    let &(_, b1, b2) = GMP_ECM_BOUNDS
        .iter()
        .find(|(d, _, _)| *d == digits)
        .unwrap_or_else(|| panic!("no bounds for {digits}-digit factors"));
    (b1, b2)
}

/// Per-curve wall-clock rows `(bits, B1, B2)` (`benches/walltime.rs`, `examples/per_curve.rs`):
/// the 256 and 1024-bit rows of the comparison with GMP-ECM.
pub const PER_CURVE_ROWS: [(u32, usize, usize); 6] = [
    (256, 50_000, 12_746_592),
    (256, 250_000, 128_992_510),
    (256, 1_000_000, 1_045_563_762),
    (1024, 50_000, 12_746_592),
    (1024, 250_000, 128_992_510),
    (1024, 1_000_000, 1_045_563_762),
];

/// Decimal digits of the cofactor of the `success_rate` numbers: larger than every target
/// factor, so the factor found is (almost always) the target one.
pub const COFACTOR_DIGITS: u32 = 40;

/// The `i`-th number of the `success_rate` example for `digits`-digit factors: a random
/// `digits`-digit prime times a random [`COFACTOR_DIGITS`]-digit prime.
pub fn success_rate_number(digits: u32, i: u64) -> Integer {
    let seed = SEED + u64::from(digits) * 1_000_000 + i;
    prime_digits(digits, seed) * prime_digits(COFACTOR_DIGITS, seed + 500_000)
}

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

/// `count` random residues modulo `n`, derived from `seed`.
pub fn residues(n: &Integer, count: usize, seed: u64) -> Vec<Integer> {
    let mut rand = RandState::new();
    rand.seed(&Integer::from(seed));
    (0..count)
        .map(|_| n.clone().random_below(&mut rand))
        .collect()
}

/// Starting point of the curve of `param` given by [`SIGMA`].
pub fn curve_point(n: &Integer, param: Param) -> Point {
    curve(n, param, &Integer::from(SIGMA)).expect("sigma must give a valid curve")
}

/// `B1` of the cheap stage 1 run by [`stage2_point`].
const STAGE2_POINT_B1: usize = 1_000;

/// A point to run stage 2 from: the first curve (default `Param`) from [`SIGMA`] on after a
/// small stage 1 (`B1 = 1000`) that finds no factor, so stage 2 always runs completely.
///
/// Stage 2 costs the same from any point: a real stage 1 (with the plan's `B1`) would only make
/// the benchmark setup much longer (it runs under Valgrind too).
pub fn stage2_point(n: &Integer) -> Point {
    let k = stage1_multiplier(STAGE2_POINT_B1);
    (SIGMA..)
        .filter_map(|sigma| curve(n, Param::default(), &Integer::from(sigma)).ok())
        .map(|p| stage1(&p, &k))
        .find(|q| q.z.clone().gcd(n) == 1)
        .unwrap()
}
