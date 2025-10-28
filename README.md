# ecm-rs

[![CI](https://github.com/skyf0l/ecm-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/skyf0l/ecm-rs/actions/workflows/ci.yml)
[![Benchmark](https://github.com/skyf0l/ecm-rs/actions/workflows/benchmark.yml/badge.svg)](https://github.com/skyf0l/ecm-rs/actions/workflows/benchmark.yml)
[![Crate.io](https://img.shields.io/crates/v/ecm.svg)](https://crates.io/crates/ecm)
[![Codecov](https://codecov.io/gh/skyf0l/ecm-rs/branch/main/graph/badge.svg)](https://codecov.io/gh/skyf0l/ecm-rs)

Lenstra's Elliptic Curve Factorization Implementation with Big Integers.

The code is based on the [sympy](https://github.com/sympy/sympy) implementation and translated to Rust.

Based on [rug](https://crates.io/crates/rug), it can use [arbitrary-precision numbers (aka BigNum)](https://en.wikipedia.org/wiki/Arbitrary-precision_arithmetic).

## Performance

[View continuous benchmark results](https://skyf0l.github.io/ecm-rs/dev/bench/)

Using a `Intel(R) Core(TM) i7-8750H CPU @ 2.20GHz` CPU, the following results were obtained:

| Number                             | sympy   | ecm-rs | sympy / ecm-rs |
| ---------------------------------- | ------- | ------ | -------------- |
| 398883434337287                    | 0.074s  | 0.057s | 1.23x faster   |
| 46167045131415113                  | 0.148s  | 0.039s | 3.8x faster    |
| 64211816600515193                  | 0.552s  | 0.017s | 32.47x faster  |
| 168541512131094651323              | 0.071s  | 0.038s | 1.87x faster   |
| 631211032315670776841              | 0.081s  | 0.128s | 0.63x faster   |
| 4132846513818654136451             | 0.266s  | 0.038s | 7.0x faster    |
| 4516511326451341281684513          | 0.495s  | 0.038s | 13.03x faster  |
| 3146531246531241245132451321       | 1.22s   | 0.22s  | 5.55x faster   |
| 4269021180054189416198169786894227 | 1.916s  | 0.018s | 106.44x faster |
| 7060005655815754299976961394452809 | 13.555s | 3.467s | 3.91x faster   |

## License

Licensed under either of

- Apache License, Version 2.0
  ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
- MIT license
  ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.
