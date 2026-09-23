//! Wall-clock benchmarks (criterion) of complete factorizations, for local use.
//!
//! Complements the instruction-count benchmarks with effects Valgrind does not measure
//! (cache, allocator, CPU-specific GMP code), on numbers too large to run under Valgrind.
//! Too noisy for shared CI runners: compare locally with
//! `cargo bench --features bench --bench walltime -- --save-baseline before` then `-- --baseline before`.
//!
//! Each iteration factors the number with several seeds, to average out the luck of the
//! curves drawn.

use criterion::{criterion_group, criterion_main, Criterion};
use ecm::bench::factor;
use rug::Integer;
use std::{hint::black_box, str::FromStr};

const SEEDS: usize = 3;

/// `(name, number)`, factored as by `ecm::ecm` (with other seeds).
const NUMBERS: [(&str, &str); 3] = [
    ("digits_21", "631211032315670776841"),
    ("digits_22", "4132846513818654136451"),
    ("digits_25", "4516511326451341281684513"),
];

fn factorize(c: &mut Criterion) {
    let mut group = c.benchmark_group("factorize");
    group.sample_size(10);

    for (name, n) in NUMBERS {
        let n = Integer::from_str(n).unwrap();
        group.bench_function(name, |b| {
            b.iter(|| {
                for seed in 0..SEEDS {
                    black_box(factor(
                        black_box(&n),
                        seed,
                        #[cfg(feature = "progress-bar")]
                        None,
                    ))
                    .unwrap();
                }
            })
        });
    }

    group.finish();
}

criterion_group!(benches, factorize);
criterion_main!(benches);
