#[macro_use]
extern crate criterion;
use criterion::Criterion;
use rug::Integer;
use std::str::FromStr;

extern crate ecm;
use ecm::ecm;

fn bench_factorization(c: &mut Criterion) {
    let mut group = c.benchmark_group("ecm");
    group.sample_size(10);

    group.bench_function("small_semiprime", |b| {
        let number = 9_804_659_461_513_846_513u64.into();
        b.iter(|| ecm(&number))
    });

    group.bench_function("medium_number", |b| {
        let number = 340_282_366_920_938_463_463_374_607_431_768_211_455u128.into();
        b.iter(|| ecm(&number))
    });

    group.bench_function("huge_number_with_many_factors", |b| {
        let n = Integer::from_str("53023048477898510482458653811989372419218837371063796289274915809338617743188682925802500284507148707878032855317656447646692819864000830091917098096680351826491670941475242192145621534512908814094787326864206640218221584804721345600283196760492110735337736483640934843900303550921096343394960797997930911828100786536812962014994852325098611681").unwrap();
        b.iter(|| ecm(&n))
    });

    group.finish();
}

criterion_group!(benches, bench_factorization);
criterion_main!(benches);
