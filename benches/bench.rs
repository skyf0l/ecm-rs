#[macro_use]
extern crate criterion;
use criterion::Criterion;
use ecm::Error;
use rug::Integer;
use std::{collections::HashMap, str::FromStr};

fn ecm(n: &Integer) -> Result<HashMap<Integer, usize>, Error> {
    ecm::ecm(
        n,
        #[cfg(feature = "progress-bar")]
        None,
    )
}

fn bench_factorization(c: &mut Criterion) {
    let mut group = c.benchmark_group("ecm");
    group.sample_size(10);

    group.bench_function("sympy_1", |b| {
        let n = Integer::from_str("398883434337287").unwrap();
        b.iter(|| ecm(&n))
    });

    group.bench_function("sympy_2", |b| {
        let n = Integer::from_str("46167045131415113").unwrap();
        b.iter(|| ecm(&n))
    });

    group.bench_function("sympy_3", |b| {
        let n = Integer::from_str("64211816600515193").unwrap();
        b.iter(|| ecm(&n))
    });

    group.bench_function("sympy_4", |b| {
        let n = Integer::from_str("168541512131094651323").unwrap();
        b.iter(|| ecm(&n))
    });

    group.bench_function("sympy_5", |b| {
        let n = Integer::from_str("631211032315670776841").unwrap();
        b.iter(|| ecm(&n))
    });

    group.bench_function("sympy_6", |b| {
        let n = Integer::from_str("4132846513818654136451").unwrap();
        b.iter(|| ecm(&n))
    });

    group.bench_function("sympy_7", |b| {
        let n = Integer::from_str("4516511326451341281684513").unwrap();
        b.iter(|| ecm(&n))
    });

    group.bench_function("sympy_8", |b| {
        let n = Integer::from_str("3146531246531241245132451321").unwrap();
        b.iter(|| ecm(&n))
    });

    group.bench_function("sympy_9", |b| {
        let n = Integer::from_str("4269021180054189416198169786894227").unwrap();
        b.iter(|| ecm(&n))
    });

    // This one takes too long to be included in regular benchmarks
    // group.bench_function("sympy_10", |b| {
    //     let n = Integer::from_str("7060005655815754299976961394452809").unwrap();
    //     b.iter(|| ecm(&n))
    // });

    group.finish();
}

criterion_group!(benches, bench_factorization);
criterion_main!(benches);
