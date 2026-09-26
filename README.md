# ecm-rs

[![CI](https://github.com/skyf0l/ecm-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/skyf0l/ecm-rs/actions/workflows/ci.yml)
[![Crate.io](https://img.shields.io/crates/v/ecm.svg)](https://crates.io/crates/ecm)
[![Codecov](https://codecov.io/gh/skyf0l/ecm-rs/branch/main/graph/badge.svg)](https://codecov.io/gh/skyf0l/ecm-rs)

Lenstra's Elliptic Curve Factorization Implementation with Big Integers.

The code is based on the [sympy](https://github.com/sympy/sympy) implementation and translated to Rust.

Based on [rug](https://crates.io/crates/rug), it can use [arbitrary-precision numbers (aka BigNum)](https://en.wikipedia.org/wiki/Arbitrary-precision_arithmetic).

## Performance

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

## Credits

- [GMP-ECM](https://gitlab.inria.fr/zimmerma/ecm), by Paul Zimmermann, Alexander Kruppa and
  the other authors listed in its
  [AUTHORS](https://gitlab.inria.fr/zimmerma/ecm/-/blob/master/AUTHORS) file, under the GNU
  LGPL version 3 or later (its library): ecm-rs translates its algorithms and code (modular
  arithmetic, stages 1 and 2, curve parametrizations, P-1, P+1, bounds and probability model).
- [SymPy](https://github.com/sympy/sympy)'s ECM (`sympy.ntheory.ecm`), from which the first
  versions were translated, Copyright (c) 2006-2023 SymPy Development Team, under the BSD
  3-Clause license (see [LICENSE-SYMPY](LICENSE-SYMPY)).

## License

Copyright (C) 2023-2026 skyf0l. Portions translated from GMP-ECM, Copyright (C) the GMP-ECM
authors.

This library is free software: you can redistribute it and/or modify it under the terms of the
GNU Lesser General Public License as published by the Free Software Foundation, either version 3
of the License, or (at your option) any later version. See [COPYING.LESSER](COPYING.LESSER) and
[COPYING](COPYING).

Versions up to 1.0.2 were released under MIT OR Apache-2.0 and remain available under those
terms.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion
in the work by you shall be licensed under the LGPL-3.0-or-later, without any additional terms
or conditions.
