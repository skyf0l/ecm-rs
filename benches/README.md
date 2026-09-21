# Benchmarks

ECM run time is `cost per curve × curves until a factor is found`, and the second term is
random. So these three measures are kept separate:

| What | Where | Measure | CI |
|---|---|---|---|
| Cost: point ops, stage 1 / stage 2 of one curve, setup (trial division, primality test, stage 1 multiplier) | `ops.rs` | Instruction counts (Valgrind), deterministic | ✅ |
| Complete factorizations of small numbers, over several seeds | `e2e.rs` | Instruction counts (Valgrind), deterministic | ✅ |
| Effectiveness: expected curves to find a 15/20/25-digit factor | `examples/success_rate.rs` | Curves tried / factors found, deterministic | ✅ |
| Complete factorizations of larger numbers | `walltime.rs` | Wall-clock time (criterion) | ❌ too noisy on shared runners |

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

```sh
cargo bench --features bench --bench walltime -- --save-baseline before
cargo bench --features bench --bench walltime -- --baseline before
```

## CI

`.github/workflows/bench.yml`:

- On pull requests, the base branch and the PR run in the same job and are compared there. The
  job fails if an instruction count grows by more than 2%.
- On pull requests, github-action-benchmark also comments a comparison with the last `main` results
  on the PR's commit.
- On pushes to `main`, results are stored on the `gh-pages` branch by
  [github-action-benchmark](https://github.com/benchmark-action/github-action-benchmark).
