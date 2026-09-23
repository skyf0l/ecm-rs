//! Instruction-count benchmarks (Valgrind/Callgrind) of the building blocks of ECM.
//!
//! Every input is deterministic (fixed seeds, fixed `sigma`), and the setup work
//! (building curves, stage 2 plans, stage 1 multipliers, ...) is not measured.
//!
//! - `arith`: batches of modular multiplications and squarings, per number of limbs.
//! - `curve`: stage 1 and stage 2 of one curve, and one complete curve at the bounds of the
//!   `success_rate` example (instructions per curve; the CI report multiplies them by the
//!   expected number of curves).
//! - `setup`: work shared by the curves (stage 2 plan, stage 1 multiplier, P-1 stage 1).
//!
//! Run with `cargo bench --features bench --bench ops` (requires Valgrind and `gungraun-runner`).

mod common;

use common::{
    bounds, curve_point, residues, semiprime_bits, stage2_point, success_rate_number,
    GMP_ECM_BOUNDS, SEED, SIGMA,
};
use ecm::bench::{
    run_curve, stage1, stage1_multiplier, stage2, ArithBatch, CurveOutcome, Param, Pm1, Point,
    Stage2Plan,
};
use gungraun::{library_benchmark, library_benchmark_group, main};
use rug::Integer;
use std::hint::black_box;

/// Operations per `arith` benchmark.
const ARITH_OPS: usize = 10_000;

/// Distinct residues multiplied by `arith_mul`.
const ARITH_VALUES: usize = 16;

/// Modular arithmetic on a modulus of `limbs` 64-bit limbs, as the curve code dispatches it
/// (Montgomery up to 10 limbs, with GMP's `mpn` functions from 11).
fn arith_input(limbs: u32) -> ArithBatch {
    let n = semiprime_bits(64 * limbs, SEED);
    assert_eq!(ArithBatch::limbs(&n), limbs as usize);
    ArithBatch::new(&n, &residues(&n, ARITH_VALUES, SEED))
}

#[library_benchmark]
#[bench::limbs_1(arith_input(1))]
#[bench::limbs_2(arith_input(2))]
#[bench::limbs_4(arith_input(4))]
#[bench::limbs_8(arith_input(8))]
#[bench::limbs_11(arith_input(11))]
#[bench::limbs_16(arith_input(16))]
fn arith_mul(batch: ArithBatch) -> Integer {
    black_box(black_box(&batch).run(ARITH_OPS, false))
}

#[library_benchmark]
#[bench::limbs_1(arith_input(1))]
#[bench::limbs_2(arith_input(2))]
#[bench::limbs_4(arith_input(4))]
#[bench::limbs_8(arith_input(8))]
#[bench::limbs_11(arith_input(11))]
#[bench::limbs_16(arith_input(16))]
fn arith_sqr(batch: ArithBatch) -> Integer {
    black_box(black_box(&batch).run(ARITH_OPS, true))
}

library_benchmark_group!(name = arith, benchmarks = [arith_mul, arith_sqr]);

/// Starting point of the curve of `param` given by `SIGMA` on a `bits`-bit modulus, and the
/// stage 1 multiplier for `b1`.
fn stage1_input(bits: u32, b1: usize, param: Param) -> (Point, Integer) {
    let n = semiprime_bits(bits, SEED);
    (curve_point(&n, param), stage1_multiplier(b1))
}

/// A point to run stage 2 from on a `bits`-bit modulus, and the plan of stage 2: the
/// baby-step giant-step one, or the polynomial one if `poly` (not the cost model's choice, so
/// that each benchmark keeps measuring the same continuation).
fn stage2_input(bits: u32, b1: usize, b2: usize, poly: bool) -> (Point, Stage2Plan) {
    let n = semiprime_bits(bits, SEED);
    let plan = if poly {
        Stage2Plan::poly(&n, b1, b2)
    } else {
        Stage2Plan::pairs(b1, b2)
    };
    (stage2_point(&n), plan)
}

const B1_20: usize = GMP_ECM_BOUNDS[1].1;
const B2_20: usize = GMP_ECM_BOUNDS[1].2;
const B1_25: usize = GMP_ECM_BOUNDS[2].1;
const B2_25: usize = GMP_ECM_BOUNDS[2].2;
const B1_35: usize = GMP_ECM_BOUNDS[4].1;
const B2_35: usize = GMP_ECM_BOUNDS[4].2;

#[library_benchmark]
#[bench::bits_128(stage1_input(128, B1_20, Param::Batch2))]
#[bench::bits_256(stage1_input(256, B1_20, Param::Batch2))]
#[bench::bits_512(stage1_input(512, B1_20, Param::Batch2))]
#[bench::bits_1024(stage1_input(1024, B1_20, Param::Batch2))]
#[bench::square_bits_256(stage1_input(256, B1_20, Param::Square))]
#[bench::suyama_bits_256(stage1_input(256, B1_20, Param::Suyama))]
fn stage1_b1_11k(input: (Point, Integer)) -> Point {
    let (p, k) = black_box(&input);
    black_box(stage1(p, k))
}

#[library_benchmark]
#[bench::bits_128_b2_1_9m(stage2_input(128, B1_20, B2_20, false))]
#[bench::bits_256_b2_12_7m(stage2_input(256, B1_25, B2_25, false))]
fn stage2_pairs(input: (Point, Stage2Plan)) -> Integer {
    let (q, plan) = black_box(&input);
    let g = stage2(q, plan);
    assert_eq!(g, 1, "stage 2 must run completely");
    black_box(g)
}

#[library_benchmark]
#[bench::bits_512_b2_12_7m(stage2_input(512, B1_25, B2_25, true))]
#[bench::bits_1024_b2_12_7m(stage2_input(1024, B1_25, B2_25, true))]
fn stage2_poly(input: (Point, Stage2Plan)) -> Integer {
    let (q, plan) = black_box(&input);
    let g = stage2(q, plan);
    assert_eq!(g, 1, "stage 2 must run completely");
    black_box(g)
}

/// Everything one curve of the `success_rate` example needs for `digits`-digit factors: its
/// first number (a `digits + 40`-digit `n`), `sigma`, the stage 1 multiplier and the stage 2
/// plan of the cost model.
fn one_curve_input(digits: u32) -> (Integer, Integer, Integer, Stage2Plan) {
    let n = success_rate_number(digits, 0);
    let (b1, b2) = bounds(digits);
    let plan = Stage2Plan::new(&n, b1, b2);
    (n, Integer::from(SIGMA), stage1_multiplier(b1), plan)
}

// One complete curve (curve setup, stage 1, stage 2) at the bounds of the `success_rate`
// example (GMP-ECM's for `digits`-digit factors): the cost side of `cost per curve x expected
// curves`. The curve finds no factor, so both stages run completely.
#[library_benchmark]
#[bench::p15(one_curve_input(15))]
#[bench::p20(one_curve_input(20))]
#[bench::p25(one_curve_input(25))]
fn one_curve(input: (Integer, Integer, Integer, Stage2Plan)) -> CurveOutcome {
    let (n, sigma, k, plan) = black_box(&input);
    let outcome = run_curve(n, Param::default(), sigma, k, plan);
    assert_eq!(
        outcome,
        CurveOutcome::Failed,
        "the curve must run both stages"
    );
    black_box(outcome)
}

library_benchmark_group!(
    name = curve,
    benchmarks = [stage1_b1_11k, stage2_pairs, stage2_poly, one_curve]
);

// Stage 2 plan search (sieve, cost model): with `B2 = 1.9M` the baby-step giant-step plan,
// with `B2 = 1e9` the polynomial one.
#[library_benchmark]
#[bench::bits_256_b2_1_9m((semiprime_bits(256, SEED), B1_20, B2_20))]
#[bench::bits_256_b2_1e9((semiprime_bits(256, SEED), B1_35, B2_35))]
fn setup_stage2_plan(input: (Integer, usize, usize)) -> Stage2Plan {
    let (n, b1, b2) = black_box(&input);
    black_box(Stage2Plan::new(n, *b1, *b2))
}

#[library_benchmark]
#[bench::b1_1m(B1_35)]
fn setup_stage1_multiplier(b1: usize) -> Integer {
    black_box(stage1_multiplier(black_box(b1)))
}

// P-1 stage 1 as run by the driver before the curves of a level (B1 = 20 x the curves' B1).
#[library_benchmark]
#[bench::b1_40k((semiprime_bits(256, SEED), 20 * GMP_ECM_BOUNDS[0].1))]
fn setup_pm1_stage1(input: (Integer, usize)) -> Integer {
    let (n, b1) = black_box(&input);
    black_box(Pm1::new().stage1(n, *b1))
}

library_benchmark_group!(
    name = setup,
    benchmarks = [setup_stage2_plan, setup_stage1_multiplier, setup_pm1_stage1]
);

main!(library_benchmark_groups = arith, curve, setup);
