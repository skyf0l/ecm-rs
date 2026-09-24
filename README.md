# ecm-rs

[![CI](https://github.com/skyf0l/ecm-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/skyf0l/ecm-rs/actions/workflows/ci.yml)
[![Crate.io](https://img.shields.io/crates/v/ecm.svg)](https://crates.io/crates/ecm)
[![Codecov](https://codecov.io/gh/skyf0l/ecm-rs/branch/main/graph/badge.svg)](https://codecov.io/gh/skyf0l/ecm-rs)

Lenstra's Elliptic Curve Factorization Implementation with Big Integers.

Based on [rug](https://crates.io/crates/rug), it can use [arbitrary-precision numbers (aka BigNum)](https://en.wikipedia.org/wiki/Arbitrary-precision_arithmetic).

## Algorithm

The implementation started as a translation of sympy's, and now uses the techniques of
[GMP-ECM](https://gitlab.inria.fr/zimmerma/ecm) (whose code and papers it follows closely):

- Montgomery modular arithmetic on fixed-size limb arrays (up to 1024 bits, with GMP's
  low-level functions from 641 bits; GMP integers above).
- GMP-ECM's curves with parametrization 2 (`-param 2`): small starting point, same torsion as
  Suyama's curves.
- Stage 2: baby-step giant-step continuation with prime pairing for small `B2`, and the
  polynomial ("FFT") continuation (product trees, multipoint evaluation, Kronecker substitution)
  for large `B2`, chosen by a cost model.
- `ecm` finds the factors from the smallest to the largest: after trial division by the primes
  below 2^16, curves are run with GMP-ECM's optimal `B1` (and default `B2`) for factors of 10,
  15, 20, ... digits in turn, each for the expected number of curves given by GMP-ECM's
  probability model, and Pollard's P-1 method (with a `B1` 20 times larger) runs before each
  size. The time to find a factor depends on its size much more than on the size of the number.
- `ecm_with_params` and `ecm_one_factor` run curves with fixed bounds.

## Performance

Using a `Intel(R) Core(TM) i7-8750H CPU @ 2.20GHz` CPU (one thread), time to factor completely
with `ecm` (fixed seed, deterministic), and with GMP-ECM 7.0.7 (`ecm -c 100000 B1`, `B1` for the
size of the second largest prime factor, mean of 5 runs, including about 2.4 ms of process
startup):

| Number                             | sympy   | ecm-rs 1.0.2 | ecm-rs  | GMP-ECM |
| ---------------------------------- | ------- | ------------ | ------- | ------- |
| 398883434337287                    | 0.074s  | 0.057s       | 0.0002s | 0.0036s |
| 46167045131415113                  | 0.148s  | 0.039s       | 0.0002s | 0.0047s |
| 64211816600515193                  | 0.552s  | 0.017s       | 0.0005s | 0.0039s |
| 168541512131094651323              | 0.071s  | 0.038s       | 0.0002s | 0.0039s |
| 631211032315670776841              | 0.081s  | 0.128s       | 0.0038s | 0.0085s |
| 4132846513818654136451             | 0.266s  | 0.038s       | 0.0004s | 0.0047s |
| 4516511326451341281684513          | 0.495s  | 0.038s       | 0.0002s | 0.0057s |
| 3146531246531241245132451321       | 1.22s   | 0.22s        | 0.0073s | 0.0105s |
| 4269021180054189416198169786894227 | 1.916s  | 0.018s       | 0.0009s | 0.0041s |
| 7060005655815754299976961394452809 | 13.555s | 3.467s       | 0.0057s | 0.056s  |

Numbers of 60 and 80 digits with one small prime factor (and a prime cofactor): mean time to
factor completely with `ecm` (3 numbers, 5 seeds each for 60 digits, 3 for 80 digits), and
to find the factor with GMP-ECM 7.0.7 given the optimal `B1` for its size (`ecm -c 100000 -one
B1`, default parametrization, as many runs). The variance is large: single runs range from
below 0.1x to 5x the mean. ecm-rs 1.0.2 chose its bounds from the size of the number
(`B1 = 26e7` for 60 digits): more than an hour per curve.

| Factor    | Number    | ecm-rs | GMP-ECM (`B1`)  |
| --------- | --------- | ------ | --------------- |
| 15 digits | 60 digits | 0.062s | 0.153s (2000)   |
| 20 digits | 60 digits | 0.55s  | 1.10s (11000)   |
| 25 digits | 60 digits | 7.3s   | 14.8s (50000)   |
| 30 digits | 80 digits | 72s    | 138s (250000)   |

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
