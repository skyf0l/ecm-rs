//! Instruction-count benchmarks (Valgrind/Callgrind) of the building blocks of ECM.
//!
//! Every input is deterministic (fixed seeds, fixed `sigma`), and the setup work
//! (building curves, running stage 1 before stage 2, ...) is not measured.
//!
//! Run with `cargo bench --features bench --bench ops` (requires Valgrind and `gungraun-runner`).

mod common;

use common::{prime_bits, semiprime_bits, GMP_ECM_BOUNDS, SEED, SIGMA};
use ecm::bench::{stage1, stage1_multiplier, stage2, suyama_curve, trial_division, Point};
use gungraun::{library_benchmark, library_benchmark_group, main};
use rug::{integer::IsPrime, Integer};
use std::{collections::HashMap, hint::black_box};

/// Starting point of the curve given by `SIGMA`, on a `bits`-bit modulus.
fn curve_point(bits: u32) -> Point {
    let n = semiprime_bits(bits, SEED);
    suyama_curve(&n, &Integer::from(SIGMA)).expect("sigma must give a valid curve")
}

/// `(2P, P, P)`: arguments of `2P.add(P, P)`.
fn add_input(bits: u32) -> (Point, Point, Point) {
    let p = curve_point(bits);
    (p.double(), p.clone(), p)
}

/// Starting point and stage 1 multiplier for `b1`.
fn stage1_input(bits: u32, b1: usize) -> (Point, Integer) {
    (curve_point(bits), stage1_multiplier(b1))
}

/// Stage 1 output and the bounds, for the first curve from `SIGMA` on where stage 1 finds
/// no factor, so stage 2 always runs completely.
fn stage2_input(bits: u32, b1: usize, b2: usize) -> (Point, usize, usize) {
    let n = semiprime_bits(bits, SEED);
    let k = stage1_multiplier(b1);
    (SIGMA..)
        .filter_map(|sigma| suyama_curve(&n, &Integer::from(sigma)).ok())
        .map(|p| stage1(&p, &k))
        .find(|q| q.z_cord.clone().gcd(&n) == 1)
        .map(|q| (q, b1, b2))
        .unwrap()
}

/// Rounds used by `ecm_one_factor` and the driver.
const PRIMALITY_REPS: u32 = 25;

const B1_15: usize = GMP_ECM_BOUNDS[0].1;
const B2_15: usize = GMP_ECM_BOUNDS[0].2;
const B1_20: usize = GMP_ECM_BOUNDS[1].1;
const B2_20: usize = GMP_ECM_BOUNDS[1].2;
const B1_25: usize = GMP_ECM_BOUNDS[2].1;

#[library_benchmark]
#[bench::bits_64(curve_point(64))]
#[bench::bits_128(curve_point(128))]
#[bench::bits_256(curve_point(256))]
#[bench::bits_512(curve_point(512))]
fn point_double(p: Point) -> Point {
    black_box(black_box(&p).double())
}

#[library_benchmark]
#[bench::bits_64(add_input(64))]
#[bench::bits_128(add_input(128))]
#[bench::bits_256(add_input(256))]
#[bench::bits_512(add_input(512))]
fn point_add(input: (Point, Point, Point)) -> Point {
    let (q, p, diff) = black_box(&input);
    black_box(q.add(p, diff))
}

library_benchmark_group!(name = point, benchmarks = [point_double, point_add]);

#[library_benchmark]
#[bench::bits_64_b1_11k(stage1_input(64, B1_20))]
#[bench::bits_128_b1_11k(stage1_input(128, B1_20))]
#[bench::bits_256_b1_11k(stage1_input(256, B1_20))]
#[bench::bits_512_b1_11k(stage1_input(512, B1_20))]
#[bench::bits_256_b1_2k(stage1_input(256, B1_15))]
fn curve_stage1(input: (Point, Integer)) -> Point {
    let (p, k) = black_box(&input);
    black_box(stage1(p, k))
}

#[library_benchmark]
#[bench::bits_64_b2_1_9m(stage2_input(64, B1_20, B2_20))]
#[bench::bits_128_b2_1_9m(stage2_input(128, B1_20, B2_20))]
#[bench::bits_256_b2_1_9m(stage2_input(256, B1_20, B2_20))]
#[bench::bits_512_b2_1_9m(stage2_input(512, B1_20, B2_20))]
#[bench::bits_256_b2_147k(stage2_input(256, B1_15, B2_15))]
fn curve_stage2(input: (Point, usize, usize)) -> Integer {
    let (q, b1, b2) = black_box(&input);
    black_box(stage2(q, *b1, *b2))
}

library_benchmark_group!(name = curve, benchmarks = [curve_stage1, curve_stage2]);

#[library_benchmark]
#[bench::bits_256(semiprime_bits(256, SEED))]
fn setup_trial_division(n: Integer) -> (HashMap<Integer, usize>, Integer) {
    black_box(trial_division(black_box(&n)))
}

// Same primality test as `ecm_one_factor` and the driver, on a prime (worst case: all rounds run).
#[library_benchmark]
#[bench::bits_128(prime_bits(128, SEED))]
#[bench::bits_256(prime_bits(256, SEED))]
#[bench::bits_512(prime_bits(512, SEED))]
fn setup_primality(p: Integer) -> IsPrime {
    black_box(black_box(&p).is_probably_prime(PRIMALITY_REPS))
}

#[library_benchmark]
#[bench::b1_2k(B1_15)]
#[bench::b1_11k(B1_20)]
#[bench::b1_50k(B1_25)]
fn setup_stage1_multiplier(b1: usize) -> Integer {
    black_box(stage1_multiplier(black_box(b1)))
}

library_benchmark_group!(
    name = setup,
    benchmarks = [
        setup_trial_division,
        setup_primality,
        setup_stage1_multiplier
    ]
);

main!(library_benchmark_groups = point, curve, setup);
