# ecm-rs

[![crates.io](https://img.shields.io/crates/v/ecm.svg)](https://crates.io/crates/ecm)
[![docs.rs](https://img.shields.io/docsrs/ecm)](https://docs.rs/ecm)
[![MSRV](https://img.shields.io/crates/msrv/ecm)](https://crates.io/crates/ecm)
[![License](https://img.shields.io/crates/l/ecm)](#license)

Fast integer factorization with Lenstra's elliptic curve method (ECM), in Rust, following [GMP-ECM](https://gitlab.inria.fr/zimmerma/ecm) and built on [rug](https://crates.io/crates/rug) (GMP) for arbitrary-precision integers.

## Usage

```rust
use ecm::{Algorithm, Event, Factorizer, ecm};
use rug::Integer;
use std::{ops::ControlFlow, time::Duration};

let n: Integer = "4516511326451341281684513".parse().unwrap();

// Complete factorization: 3^2 * 39869 * 131743543 * 95542348571.
let factors = ecm(&n).unwrap();
assert_eq!(factors[&Integer::from(3)], 2);
assert_eq!(factors.len(), 4);

// With a seed, a timeout and progress events (returning `ControlFlow::Break` interrupts the
// factorization too, as `interrupt_flag` does): `factor_partial` also returns what was found
// when the factorization is interrupted (or fails).
let result = Factorizer::new()
    .seed(42)
    .timeout(Duration::from_secs(60))
    .on_event(|event| {
        match event {
            Event::Level { digits: Some(digits), curves: Some(curves), .. } => {
                println!("factors of {digits} digits: {curves} curves");
            }
            Event::Factor { factor, method, .. } => println!("{factor} found by {method}"),
            _ => {}
        }
        ControlFlow::Continue(())
    })
    .factor_partial(&n);
println!("primes: {:?}, unfactored: {:?}", result.primes, result.unfactored);

// On 4 threads (by default, 1): the same factors, curves and events (but their durations)
// as with one; the callback still runs on this thread.
let factors = Factorizer::new().threads(4).factor(&n).unwrap();
assert_eq!(factors.len(), 4);

// Fixed bounds (as GMP-ECM's `ecm -c 100 11000 1873422`), or only P-1 or P+1 (as
// `ecm -pp1 -x0 2/7 100000`).
let factors = Factorizer::new()
    .b1(11_000)
    .b2(1_873_422)
    .curves(100)
    .factor(&n)
    .unwrap();
assert_eq!(factors.len(), 4);
let factor = Factorizer::new().algorithm(Algorithm::Pm1).b1(100_000).find_factor(&n);
let factor = Factorizer::new()
    .algorithm(Algorithm::Pp1)
    .x0(2.into(), 7.into())
    .b1(100_000)
    .find_factor(&n);
```

## Command line tool

`ecm-rs`, with GMP-ECM-like options (the library alone compiles none of its dependencies):

```sh
cargo install ecm --features cli
```

Each number (an argument, or a line of the standard input) is factored completely, and printed
with its prime factors in increasing order (parts marked `(composite)` when the factorization is
incomplete: `--one`, `--timeout`, Ctrl-C, `--curves`, `--pm1` or `--pp1` exhausted). The numbers
can be expressions: `+ - * /`, `^`, parentheses, `n!` (factorial), `n#` (primorial), as
GMP-ECM's.

```text
$ ecm-rs '2^67-1' '10!+1' 17
147573952589676412927 = 193707721 * 761838257287
3628801 = 11 * 329891
17 = 17

$ ecm-rs --json '2^67-1'
{"input":"2^67-1","n":"147573952589676412927","factors":[{"p":"193707721","exponent":1},{"p":"761838257287","exponent":1}],"unfactored":[],"complete":true,"error":null,"time":0.000326099}

$ echo 15658598057181786459081452046251445462002559800474409088109 | ecm-rs -v --b1 11000 --sigma 1:1176292814
Input number is 15658598057181786459081452046251445462002559800474409088109 (59 digits)
Trial division below 2^16: no factor, cofactor has 59 digits
Using B1=11000, B2=1873420 on C59: curves 1
Curve 1/1: sigma=1:1176292814, Step 1 took 5.2ms, Step 2 took 5.8ms
********** Factor found by ECM stage 2 (B1=11000, B2=1873420, sigma=1:1176292814): 13507140964289979319
Found prime factor of 20 digits: 13507140964289979319
Prime cofactor 1159282937712710938601499347662537052411 has 40 digits
15658598057181786459081452046251445462002559800474409088109 = 13507140964289979319 * 1159282937712710938601499347662537052411
```

The curves run on all the available threads by default (`-t N` to choose): the output is the
same with any number of threads (but the times of `-v`).

On a terminal, a progress bar (on stderr) shows the level being run, with the curves done out of
the expected number and the expected time of the level. Results go to stdout, diagnostics to
stderr. Ctrl-C, `--timeout SECS` (for each number) and `--total-timeout SECS` (for the whole
run) print what was found so far. With `--json`, an invalid input gives `"n": null`,
`"error": "invalid"` and a `"message"`. See `ecm-rs --help` for all the options.

| GMP-ECM              | ecm-rs                                                   |
| -------------------- | -------------------------------------------------------- |
| `ecm B1`             | `ecm-rs --b1 B1 -c 1` (GMP-ECM runs one curve by default) |
| `ecm B1 B2`          | `ecm-rs --b1 B1 --b2 B2 -c 1`                            |
| `ecm -c N B1`        | `ecm-rs --b1 B1 -c N` (`N` curves per composite part)    |
| `ecm -sigma 1:S B1`  | `ecm-rs --b1 B1 --sigma 1:S` (one curve, as GMP-ECM)     |
| `ecm -param 0 B1`    | `ecm-rs --b1 B1 -c 1 --param 0`                          |
| `ecm -one ...`       | `ecm-rs --one ...`                                       |
| `ecm -pm1 B1 B2`     | `ecm-rs --pm1 --b1 B1 --b2 B2`                           |
| `ecm -pp1 -x0 2/7 B1 B2` | `ecm-rs --pp1 --b1 B1 --b2 B2` (seed `2/7` by default; GMP-ECM: random) |
| `ecm -pm1/-pp1 -x0 X`| `ecm-rs --pm1/--pp1 --x0 X`                              |
| `ecm -maxmem MB`     | `ecm-rs --maxmem MB`                                     |
| `ecm -base2 K`, `ecm -nobase2` | `ecm-rs --base2 K`, `ecm-rs --nobase2` (automatic by default, as GMP-ECM) |
| `ecm -primetest`     | `ecm-rs --primetest` (prints `N: prime` or `N: composite`) |
| `ecm -printconfig`   | `ecm-rs --printconfig`                                   |
| `ecm -q`, `ecm -v`   | `ecm-rs -q`, `ecm-rs -v`                                 |
| (none)               | `ecm-rs --b1 B1`: curves until a factor is found         |
| (none)               | `ecm-rs N`: complete factorization, bounds by factor size |

Exit status: bits as GMP-ECM's, for the last number (bits 1 and 16 for any number).

| Status | Meaning                                                                 |
| ------ | ----------------------------------------------------------------------- |
| 0      | No factor found (curves exhausted; `--primetest`: composite)            |
| 1      | Error: an invalid number (the others are still processed)               |
| 2      | A composite factor found, the cofactor is composite (`--one`)           |
| 6      | A prime factor found, the cofactor is composite (incomplete)            |
| 8      | The number is prime (or 1)                                              |
| 10     | A composite factor found, the cofactor is prime (`--one`)               |
| 14     | Factored completely (`--one`: a prime factor, a prime cofactor)         |
| +16    | `--timeout` or `--total-timeout` interrupted a factorization            |
| 64     | Invalid command line                                                    |
| 130    | Interrupted by Ctrl-C                                                   |

## Algorithm

The implementation started as a translation of sympy's, and now uses the techniques of
[GMP-ECM](https://gitlab.inria.fr/zimmerma/ecm) (whose code and papers it follows closely):

- Montgomery modular arithmetic on fixed-size limb arrays (up to 1024 bits, with GMP's
  low-level functions from 641 bits; GMP integers above).
- GMP-ECM's "special division" for the divisors of `2^k +- 1` (Mersenne, Fermat and
  Cunningham numbers and their cofactors): with `k` at most 1.4 times their size, the curves,
  P-1 and P+1 compute modulo `2^k +- 1`, where a product is reduced by a shift and an addition
  instead of a Montgomery reduction, when the measured costs say it is faster (from about 400
  bits with `k` close to their size, 770 bits with `k` up to 1.4 times). Stage 1 is 1.1 to 2.8
  times faster (2.4 times at 1024 bits, 2.6 at 2048 and 4096 bits: as fast as GMP-ECM's),
  stage 2 up to 1.5 times (it stays modulo the number when that is cheaper).
  `Factorizer::base2` (`Base2Mode`) and `--base2 K`/`--nobase2` force it or turn it off.
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
- Williams' P+1 method (`Algorithm::Pp1`, `--pp1`) runs alone when asked, as GMP-ECM's
  `-pp1`: stage 1 with Montgomery's PRAC Lucas chains, stage 2 shared with P-1 and the curves.
  It is not part of `ecm`: after P-1, the factors it adds (`p = 2 mod 3` with a smooth `p + 1`)
  do not make the search faster (expected time within 0.1% at best, by GMP-ECM's model with
  measured costs).
- `ecm_with_params` and `ecm_one_factor` run curves with fixed bounds.
- Threads (`Factorizer::threads`, `-t`): the curves run in parallel, with P-1 and the curves
  of the next level alongside them, in the order of a search on one thread, whose results they
  give: the factor found is the one of the first curve (or P-1 run) in this order that finds
  one, once all the ones before it ran. The first level runs on the calling thread.
- `Factorizer` has all the options (seed, fixed bounds, curves, `sigma`, parametrization,
  P-1 or P+1 only, stage 2 memory, special division, threads), and reports events (levels,
  curves with their `sigma` and stage durations, P-1 and P+1 runs, factors and primes) to a
  callback, which can interrupt the factorization. An interruption flag or a timeout interrupt it even during a
  curve.

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

With threads (`Factorizer::threads`, `ecm-rs -t N`), same machine (6 cores, 12 hardware
threads), mean time to factor completely the numbers above (60 digits with a 20 or 25-digit
factor: 3 numbers x 5 or 3 seeds; 80 digits with a 30-digit factor: 3 numbers), and curves per
second with fixed bounds on a 77-digit number. The results are the same with any number of
threads; the second hardware thread of a core adds little (the arithmetic saturates it).

| Threads                    | 1      | 2      | 4      | 6      | 12     |
| -------------------------- | ------ | ------ | ------ | ------ | ------ |
| 20 digits, 60-digit number | 0.35s  | 0.19s  | 0.11s  | 0.11s  | 0.10s  |
| 25 digits, 60-digit number | 6.2s   | 3.2s   | 1.7s   | 1.4s   | 1.4s   |
| 30 digits, 80-digit number | 99s    | 50s    | 28s    | 22s    | 22s    |
| curves/s, `B1` = 11000     | 155    | 302    | 581    | 758    | 818    |
| curves/s, `B1` = 50000     | 29.8   | 58.9   | 106    | 126    | 134    |
| curves/s, `B1` = 250000    | 6.1    | 12.0   | 22.9   | 26.1   | 28.8   |

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
[COPYING](./COPYING).

Versions up to 1.0.2 were released under MIT OR Apache-2.0 and remain available under those
terms.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion
in the work by you shall be licensed under the LGPL-3.0-or-later, without any additional terms
or conditions.
