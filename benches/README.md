# Benchmarks

ECM run time is `cost per curve x curves until a factor is found`, and the second term is
random. So these three measures are kept separate:

| What | Where | Measure | CI |
|---|---|---|---|
| Cost: modular multiplication/squaring per limb count, stage 1 / stage 2 of one curve, one complete curve at the `success_rate` bounds, setup (stage 2 plan, stage 1 multiplier, P-1 stage 1) | `ops.rs` | Instruction counts (Valgrind), deterministic | yes |
| Driver overhead: complete factorizations of small README numbers, over several seeds | `e2e.rs` | Instruction counts (Valgrind), deterministic | yes |
| Effectiveness: expected curves to find a 15/20/25-digit factor | `examples/success_rate.rs` | Curves tried / factors found, deterministic | yes |
| Expected cost to find a factor: expected curves x instructions per curve (`ops::curve::one_curve`) | CI report | Derived | yes |
| Stage 1 / stage 2 of one curve at 256 and 1024 bits, `B1` up to 1M | `walltime.rs` | Wall-clock time (criterion) | no, too noisy on shared runners |
| Same rows, side by side with GMP-ECM | `scripts/compare_gmp_ecm.sh` (`examples/per_curve.rs`) | Wall-clock time | no |

`ops.rs` groups:

- `arith`: 10000 chained multiplications (`arith_mul`) or squarings (`arith_sqr`) modulo 1, 2,
  4, 8, 11 and 16-limb numbers, dispatched as in the curve code (our Montgomery code up to 10
  limbs, GMP's `mpn` functions from 11): shows inlining regressions and the 10/11 limb switch.
- `curve`: stage 1 (`B1 = 11000`) with the default curves (`-param 2`) at 128 to 1024 bits, and
  with `-param 1` and Suyama's at 256 bits; stage 2 with the baby-step giant-step continuation
  (128 bits / `B2 = 1.9M`, 256 bits / 12.7M) and the polynomial one (512 and 1024 bits / 12.7M);
  `one_curve`: one complete curve at the bounds of `success_rate` for 15, 20 and 25-digit
  factors, on its first number (a factor times a 40-digit prime).
- `setup`: stage 2 plan search at 256 bits (`B2 = 1.9M` and 1e9), stage 1 multiplier for
  `B1 = 1M`, P-1 stage 1.

## Instruction counts

Requires [Valgrind](https://valgrind.org/) and `gungraun-runner` with the same version as the
`gungraun` dev-dependency:

```sh
cargo install gungraun-runner --version 0.19.4 --locked

cargo bench --features bench --bench ops --bench e2e -- --parallel=auto
# Only some benchmarks: FILE::GROUP::FUNCTION::ID wildcard
cargo bench --features bench --bench ops -- 'ops::curve::*'
```

Compare two versions of the code:

```sh
git checkout main && cargo bench --features bench --bench ops --bench e2e -- --save-baseline=main
git checkout -    && cargo bench --features bench --bench ops --bench e2e -- --baseline=main
```

## Success rate

```sh
cargo run --release --features bench --example success_rate -- --sizes 15,20,25 --json out.json
cargo run --release --features bench --example success_rate -- --compare out.json
```

## Wall-clock time

Instruction counts don't see cache effects nor the build profile: inlining across crates with
`lto = "fat"`, register spills of the carry chains with `codegen-units = 1`. `walltime.rs`
measures stage 1 and stage 2 of one curve for the rows of the comparison with GMP-ECM (about 3
minutes), with the default profile or with the `bench-lto` profile (`bench` + `lto = "fat"`,
`codegen-units = 1`):

```sh
cargo bench --features bench --bench walltime -- --save-baseline before
cargo bench --features bench --bench walltime -- --baseline before
cargo bench --profile bench-lto --features bench --bench walltime
```

Side by side with a local GMP-ECM binary (same numbers, same curves: `ecm -c 1 -param 2 -sigma
2:1234567 B1 B2`), from `examples/per_curve.rs`:

```sh
taskset -c 2 scripts/compare_gmp_ecm.sh path/to/ecm [--reps 3] [--bits 256,1024]
PROFILE=bench-lto scripts/compare_gmp_ecm.sh path/to/ecm
```

## CI

`.github/workflows/bench.yml`:

- On pull requests, the base branch and the PR run in the same job and are compared there. The
  job fails if an instruction count grows by more than 2%.
- On pull requests, one comment compares the PR with its base branch, worst changes first
  (updated on each push), and adds the expected cost to find a factor: expected curves
  (`success_rate`) x instructions per curve (`one_curve`), for the PR and its base branch.
- The benchmarks are built without debug info in CI (`CARGO_PROFILE_BENCH_DEBUG=false`): half
  the build time, same instruction counts.
- On pushes to `main`, results are stored on the `gh-pages` branch by
  [github-action-benchmark](https://github.com/benchmark-action/github-action-benchmark).
