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
- `ecm_with_params` returns `Error::ECMFailed` when a composite part is not split after
  `max_curve` curves, instead of returning that composite as if it were a prime factor. The
  other parts are still factored (`Factorizer::factor_partial` returns them).
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
  (`b1`, `b2`, `curves`, `sigma`), P-1 only (`pm1`), memory of the polynomial stage 2
  (`max_memory`), `factor`, `factor_partial` (the primes and the unfactored parts found so far
  on failure or interruption) and `find_factor` (one proper factor).
- `Event` (trial division, P-1 runs, levels with their expected number of curves, curves with
  their `sigma` and stage durations, factors with their `Method`, primes), reported to the
  `Factorizer::on_event` callback, whose `ControlFlow::Break` interrupts the factorization.
  Without a callback, the events cost nothing.

### Performance

- Montgomery arithmetic on fixed-size limb arrays (up to 1024 bits; GMP's low-level functions
  from 641 bits), with BMI2/ADX copies of the hot loops chosen at run time.
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
