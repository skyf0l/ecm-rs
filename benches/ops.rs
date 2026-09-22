//! Instruction-count benchmarks (Valgrind/Callgrind) of the building blocks of ECM.
//!
//! Every input is deterministic (fixed seeds, fixed `sigma`), and the setup work
//! (building curves, running stage 1 before stage 2, ...) is not measured.
//!
//! Run with `cargo bench --features bench --bench ops` (requires Valgrind and `gungraun-runner`).

mod common;

use common::{prime_bits, semiprime_bits, GMP_ECM_BOUNDS, SEED, SIGMA};
use ecm::bench::{
    curve, stage1, stage1_multiplier, stage2, trial_division, Param, Point, Stage2Plan,
};
use gungraun::{library_benchmark, library_benchmark_group, main};
use rug::{integer::IsPrime, Integer};
use std::{collections::HashMap, hint::black_box};

/// Starting point of the curve of `param` given by `SIGMA`, on a `bits`-bit modulus.
fn curve_point(bits: u32, param: Param) -> Point {
    let n = semiprime_bits(bits, SEED);
    curve(&n, param, &Integer::from(SIGMA)).expect("sigma must give a valid curve")
}

/// Starting point and stage 1 multiplier for `b1`.
fn stage1_input(bits: u32, b1: usize, param: Param) -> (Point, Integer) {
    (curve_point(bits, param), stage1_multiplier(b1))
}

/// Stage 1 output and the plan of stage 2 (the cheapest one, or the polynomial one if `poly`),
/// for the first curve from `SIGMA` on where stage 1 finds no factor, so stage 2 always runs
/// completely.
fn stage2_input(bits: u32, b1: usize, b2: usize, poly: bool) -> (Point, Stage2Plan) {
    let n = semiprime_bits(bits, SEED);
    let k = stage1_multiplier(b1);
    let plan = if poly {
        Stage2Plan::poly(&n, b1, b2)
    } else {
        Stage2Plan::new(&n, b1, b2)
    };
    (SIGMA..)
        .filter_map(|sigma| curve(&n, Param::default(), &Integer::from(sigma)).ok())
        .map(|p| stage1(&p, &k))
        .find(|q| q.z_cord.clone().gcd(&n) == 1)
        .map(|q| (q, plan))
        .unwrap()
}

/// Rounds used by `ecm_one_factor` and the driver.
const PRIMALITY_REPS: u32 = 25;

const B1_15: usize = GMP_ECM_BOUNDS[0].1;
const B2_15: usize = GMP_ECM_BOUNDS[0].2;
const B1_20: usize = GMP_ECM_BOUNDS[1].1;
const B2_20: usize = GMP_ECM_BOUNDS[1].2;
const B1_25: usize = GMP_ECM_BOUNDS[2].1;
const B2_25: usize = GMP_ECM_BOUNDS[2].2;

#[library_benchmark]
#[bench::bits_64_b1_11k(stage1_input(64, B1_20, Param::Square))]
#[bench::bits_128_b1_11k(stage1_input(128, B1_20, Param::Square))]
#[bench::bits_256_b1_11k(stage1_input(256, B1_20, Param::Square))]
#[bench::bits_512_b1_11k(stage1_input(512, B1_20, Param::Square))]
#[bench::bits_1024_b1_11k(stage1_input(1024, B1_20, Param::Square))]
#[bench::bits_256_b1_2k(stage1_input(256, B1_15, Param::Square))]
#[bench::suyama_bits_256_b1_11k(stage1_input(256, B1_20, Param::Suyama))]
fn curve_stage1(input: (Point, Integer)) -> Point {
    let (p, k) = black_box(&input);
    black_box(stage1(p, k))
}

#[library_benchmark]
#[bench::bits_64_b2_1_9m(stage2_input(64, B1_20, B2_20, false))]
#[bench::bits_128_b2_1_9m(stage2_input(128, B1_20, B2_20, false))]
#[bench::bits_256_b2_1_9m(stage2_input(256, B1_20, B2_20, false))]
#[bench::bits_512_b2_1_9m(stage2_input(512, B1_20, B2_20, false))]
#[bench::bits_256_b2_147k(stage2_input(256, B1_15, B2_15, false))]
#[bench::poly_bits_256_b2_1_9m(stage2_input(256, B1_20, B2_20, true))]
#[bench::poly_bits_512_b2_12_7m(stage2_input(512, B1_25, B2_25, true))]
fn curve_stage2(input: (Point, Stage2Plan)) -> Integer {
    let (q, plan) = black_box(&input);
    black_box(stage2(q, plan))
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

#[library_benchmark]
#[bench::b2_147k((semiprime_bits(256, SEED), B1_15, B2_15))]
#[bench::b2_1_9m((semiprime_bits(256, SEED), B1_20, B2_20))]
#[bench::b2_12_7m((semiprime_bits(256, SEED), B1_25, B2_25))]
fn setup_stage2_plan(input: (Integer, usize, usize)) -> Stage2Plan {
    let (n, b1, b2) = black_box(&input);
    black_box(Stage2Plan::new(n, *b1, *b2))
}

library_benchmark_group!(
    name = setup,
    benchmarks = [
        setup_trial_division,
        setup_primality,
        setup_stage1_multiplier,
        setup_stage2_plan
    ]
);

main!(library_benchmark_groups = curve, setup);
