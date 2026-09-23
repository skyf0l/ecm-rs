//! Wall-clock benchmarks (criterion) of stage 1 and stage 2 of one curve, for local use.
//!
//! Complements the instruction-count benchmarks with what Valgrind does not see: cache and
//! memory effects, CPU-specific GMP code, and the build profile (inlining across crates with
//! `lto = "fat"`, register spills of the carry chains with `codegen-units = 1`), on the sizes
//! and bounds of the comparison with GMP-ECM (256 and 1024 bits, `B1`/`B2` from 50k/12.7M to
//! 1M/1e9). Too noisy for shared CI runners.
//!
//! ```text
//! # Default bench profile
//! cargo bench --features bench --bench walltime -- --save-baseline before
//! cargo bench --features bench --bench walltime -- --baseline before
//! # With the optimizations of a final build (`lto = "fat"`, `codegen-units = 1`)
//! cargo bench --profile bench-lto --features bench --bench walltime
//! # Only some rows
//! cargo bench --features bench --bench walltime -- 'stage2/bits_256'
//! ```
//!
//! A full run takes about 3 minutes (the 1M/1e9 rows cost seconds per iteration). Stage 2 runs
//! with the plan of the cost model, built beforehand (not measured). See also
//! `examples/per_curve.rs` and `scripts/compare_gmp_ecm.sh` for a comparison with GMP-ECM.

mod common;

use common::{curve_point, semiprime_bits, stage2_point, PER_CURVE_ROWS, SEED};
use criterion::{criterion_group, criterion_main, Criterion, SamplingMode};
use ecm::bench::{stage1, stage1_multiplier, stage2, Param, Stage2Plan};
use std::{hint::black_box, time::Duration};

/// Benchmark id of a row: `bits_256_b1_50k`.
fn id(bits: u32, b1: usize) -> String {
    format!("bits_{bits}_b1_{}k", b1 / 1000)
}

fn per_curve(c: &mut Criterion) {
    let mut group = c.benchmark_group("stage1");
    group
        .sample_size(10)
        .sampling_mode(SamplingMode::Flat)
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(5));
    for (bits, b1, _) in PER_CURVE_ROWS {
        let n = semiprime_bits(bits, SEED);
        let p = curve_point(&n, Param::default());
        let k = stage1_multiplier(b1);
        group.bench_function(id(bits, b1), |b| b.iter(|| stage1(black_box(&p), &k)));
    }
    group.finish();

    let mut group = c.benchmark_group("stage2");
    group
        .sample_size(10)
        .sampling_mode(SamplingMode::Flat)
        .warm_up_time(Duration::from_secs(1))
        .measurement_time(Duration::from_secs(5));
    for (bits, b1, b2) in PER_CURVE_ROWS {
        let n = semiprime_bits(bits, SEED);
        let q = stage2_point(&n);
        let plan = Stage2Plan::new(&n, b1, b2);
        group.bench_function(id(bits, b1), |b| b.iter(|| stage2(black_box(&q), &plan)));
    }
    group.finish();
}

criterion_group!(benches, per_curve);
criterion_main!(benches);
