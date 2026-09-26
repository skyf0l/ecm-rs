# Changelog

## 2.0.0

A port of GMP-ECM's techniques: much faster curves, a driver that finds the factors from the
smallest to the largest, and a builder with progress events and cancellation.

### Breaking changes

- The `progress-bar` feature (and the `indicatif` dependency) is removed: `ecm`,
  `ecm_with_params` and `ecm_one_factor` no longer take a progress bar. For progress, use
  `Factorizer::on_event`.
- `Error` is `#[non_exhaustive]`, with the new variants `Interrupted` (the event callback
  stopped the factorization) and `InvalidOption` (incompatible `Factorizer` options). It now
  derives `Clone`, `PartialEq` and `Eq`.
- `ecm` and `ecm_with_params` return `Error::ECMFailed` when a composite part is not split
  (after `max_curve` curves for `ecm_with_params`), instead of returning that composite as if
  it were a prime factor. The other parts are still factored (`Factorizer::factor_partial`
  returns them). A composite factor found by a curve is factored too, instead of being
  returned as a prime.
- `ecm_one_factor` runs at most `max_curve` curves (1.x ran one more).
- `ecm_one_factor` and `ecm_with_params` return `Error::BoundsTooSmall` if `b1 < 6` or
  `b2 < 4`.
- `ecm` and `ecm_with_params` panic if `n <= 0` (they used to loop forever on 0), and
  `ecm_one_factor` panics if `n <= 1`: these numbers have no (proper) factorization.
- Trial division removes the primes below 2^16 (it was the first 100,000 primes, up to
  1,299,709): larger factors are found by P-1 and the curves.
- The curves are GMP-ECM's parametrization 2 (Suyama's before), drawn from a different
  generator: the same seed gives other curves than 1.x (the factors found are the same).
- `ecm` chooses its bounds by factor size (GMP-ECM's table of optimal `B1` for 10 to 65
  digits, with P-1 before each size), not from the size of the number.
- Edition 2024, minimum supported Rust version 1.87.

### Added

- `Factorizer`: seed, parametrization (`Param`, GMP-ECM's `-param 0/1/2`), fixed bounds
  (`b1`, `b2`, `curves`, `sigma`), P-1 or P+1 only (`algorithm`, with the starting value
  `x0`), memory of the polynomial stage 2
  (`max_memory`), `factor`, `factor_partial` (the primes and the unfactored parts found so far
  on failure or interruption) and `find_factor` (one proper factor).
- Williams' P+1 method (`Algorithm::Pp1`, GMP-ECM's `-pp1`): stage 1 with PRAC Lucas chains
  on the Montgomery arithmetic, stage 2 with the continuations of the curves, seed `2/7` by
  default. It only runs when asked: after P-1, it does not make the default search faster.
- GMP-ECM's special division (`-base2`): the numbers dividing `2^k +- 1` with `k` at most 1.4
  times their size compute modulo `2^k +- 1` (reduction by a shift and an addition) in the
  curves, P-1 and P+1, when it is faster than the Montgomery arithmetic (by measured costs).
  `Factorizer::base2` with `Base2Mode::{Auto, Off, Force(k)}`, `Event::Base2`, and in `ecm-rs`
  `--base2 K`, `--nobase2` and `-v`'s "Using special division for factor of 2^k+1".
- `Event` (trial division, P-1 and P+1 runs, levels with their expected number of curves,
  curves with their `sigma` and stage durations, factors with their `Method`, primes), reported
  to the `Factorizer::on_event` callback, whose `ControlFlow::Break` interrupts the factorization
  between two curves. Without a callback, the events cost nothing.
- `Factorizer::interrupt_flag` (an `Arc<AtomicBool>`, for a Ctrl-C handler or another thread)
  and `Factorizer::timeout` interrupt the factorization even during a curve or P-1: usually
  within a few milliseconds.
- `Factorizer::threads` (default 1; 0 for the available parallelism): the curves run in
  parallel, with P-1 and the next levels alongside them, with the same results (factors,
  curves, events but their durations) whatever the number of threads. The callback stays on
  the calling thread (no `Send` bound). `ecm`, `ecm_with_params` and `ecm_one_factor` use one
  thread.
- The `ecm-rs` command line tool (`cargo install ecm --features cli`): complete factorization
  by default, GMP-ECM-like options (`--b1`, `--b2`, `-c`, `--sigma`, `--param`, `--one`,
  `--pm1`, `--pp1`, `--x0`, `--maxmem`, `--base2`, `--nobase2`, `--primetest`, `--printconfig`,
  `-q`, `-v`), `-t`/`--threads` (default: the available parallelism), GMP-ECM's input
  expressions, `--json`, `--seed`, `--timeout` (for each number)
  and `--total-timeout`, a progress bar, partial results on Ctrl-C, and GMP-ECM's exit status
  bits. The `cli` feature only adds its dependencies (`clap`, `ctrlc`, `indicatif`,
  `serde_json`): the library compiles none of them.

### Performance

- Montgomery arithmetic on fixed-size limb arrays (up to 1024 bits; GMP's low-level functions
  from 641 bits), with BMI2/ADX copies of the hot loops chosen at run time.
- Special division for the divisors of `2^k +- 1`: stage 1 of the curves and of P-1 2.4 to 2.8
  times faster from 1024 bits (as fast as GMP-ECM's), stage 2 up to 1.5 times (modulo the
  number when that is cheaper).
- Stage 1: GMP-ECM's parametrization 2 curves (small starting point), a product-tree
  multiplier.
- Stage 2: baby-step giant-step continuation with prime pairing for small `B2`, and GMP-ECM's
  polynomial continuation (product trees, multipoint evaluation, Kronecker substitution) for
  large `B2`, chosen by a cost model; GMP-ECM's default `B2`.
- Driver: factors found from the smallest to the largest (expected curves from GMP-ECM's
  probability model), P-1 before each size, cofactors resuming the search.
- The numbers of the README are factored in 0.2 to 7.3 ms (1.0.2: 17 ms to 3.5 s); a 25-digit
  factor of a 60-digit number is found in about 7 s (1.0.2: over an hour per curve at its
  bounds), about twice as fast as GMP-ECM 7.0.7 with the optimal `B1` for the factor size.
