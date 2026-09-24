//! Instruction-count benchmarks (Valgrind/Callgrind) of the driver overhead: complete
//! factorizations (`ecm::ecm`) of small numbers of the README.
//!
//! These numbers are factored in a few curves (or by trial division, P-1, ...), so the cost is
//! mostly the driver's: trial division, primality tests, stage 1 multipliers, stage 2 plans,
//! choice of the bounds. The cost of the curves themselves is measured by `ops.rs`.
//!
//! The run time depends on luck (which curves are drawn), so each benchmark factors the same
//! number with several seeds and measures the total. Results are checked: a wrong factorization
//! makes the benchmark fail instead of reporting a fake speedup.
//!
//! Run with `cargo bench --features bench --bench e2e` (requires Valgrind and `gungraun-runner`).

use ecm::bench::factor;
use gungraun::{library_benchmark, library_benchmark_group, main};
use rug::{Integer, integer::IsPrime, ops::Pow};
use std::{collections::HashMap, hint::black_box, str::FromStr};

/// Factors `n` as `ecm::ecm` does, once per seed, and checks the results.
fn factor_with_seeds(n: &str, seeds: u64) -> Vec<HashMap<Integer, usize>> {
    let n = Integer::from_str(n).unwrap();
    (0..seeds as usize)
        .map(|seed| {
            let factors = factor(
                black_box(&n),
                seed,
                #[cfg(feature = "progress-bar")]
                None,
            )
            .unwrap();
            check(&n, &factors);
            factors
        })
        .collect()
}

fn check(n: &Integer, factors: &HashMap<Integer, usize>) {
    let product = factors.iter().fold(Integer::from(1), |acc, (p, e)| {
        acc * p.clone().pow(*e as u32)
    });
    assert_eq!(&product, n, "wrong factorization of {n}: {factors:?}");
    for p in factors.keys() {
        assert_ne!(
            p.is_probably_prime(30),
            IsPrime::No,
            "composite factor {p} of {n}"
        );
    }
}

#[library_benchmark]
#[bench::digits_15(("398883434337287", 5))]
#[bench::digits_17(("46167045131415113", 5))]
#[bench::digits_25(("4516511326451341281684513", 3))]
fn factorize(input: (&str, u64)) -> Vec<HashMap<Integer, usize>> {
    black_box(factor_with_seeds(input.0, input.1))
}

library_benchmark_group!(name = e2e, benchmarks = [factorize]);

main!(library_benchmark_groups = e2e);
