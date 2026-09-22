//! Measures how effective each ECM curve is: the fraction of curves that find a factor of a
//! given size, and the expected number of curves to find it.
//!
//! This complements the instruction-count benchmarks: expected time to find a factor is
//! roughly `cost per curve × expected curves`. A change can make curves cheaper but less
//! effective (or the opposite), and only both numbers together tell if it is an improvement.
//!
//! Everything is deterministic (fixed numbers and `sigma` values), so results only change
//! when the algorithm changes.
//!
//! ```text
//! cargo run --release --features bench --example success_rate -- \
//!     [--sizes 15,20,25] [--numbers 20] [--curves 50] [--json out.json] [--compare base.json] \
//!     [--markdown comparison.md]
//! ```
//!
//! - `--json`: writes results in the `customSmallerIsBetter` format of
//!   github-action-benchmark.
//! - `--compare`: prints a markdown table comparing the results with a previous `--json` file.
//!   A change is only marked better or worse when the 95% confidence intervals don't overlap.
//! - `--markdown`: also writes that table to a file.

#[path = "../benches/common/mod.rs"]
mod common;

use common::{prime_digits, GMP_ECM_BOUNDS, SEED};
use ecm::bench::{run_curve, stage1_multiplier, CurveOutcome};
use rug::{rand::RandState, Integer};
use serde_json::{json, Value};
use std::{
    sync::atomic::{AtomicUsize, Ordering},
    thread,
    time::Instant,
};

/// Decimal digits of the cofactor: larger than every target factor, so the factor found is
/// (almost always) the target one.
const COFACTOR_DIGITS: u32 = 40;

/// Approximate expected curves from GMP-ECM's documentation for the same bounds. For reference
/// only: GMP-ECM uses a much stronger stage 2.
const GMP_ECM_EXPECTED_CURVES: [(u32, u32); 3] = [(15, 25), (20, 90), (25, 300)];

struct Args {
    sizes: Vec<u32>,
    numbers: u64,
    curves: u64,
    json: Option<String>,
    compare: Option<String>,
    markdown: Option<String>,
}

fn parse_args() -> Args {
    let mut args = Args {
        sizes: vec![15, 20],
        numbers: 20,
        curves: 50,
        json: None,
        compare: None,
        markdown: None,
    };
    let mut it = std::env::args().skip(1);
    while let Some(arg) = it.next() {
        let mut value = || {
            it.next()
                .unwrap_or_else(|| panic!("missing value for {arg}"))
        };
        match arg.as_str() {
            "--sizes" => {
                args.sizes = value().split(',').map(|s| s.parse().unwrap()).collect();
            }
            "--numbers" => args.numbers = value().parse().unwrap(),
            "--curves" => args.curves = value().parse().unwrap(),
            "--json" => args.json = Some(value()),
            "--compare" => args.compare = Some(value()),
            "--markdown" => args.markdown = Some(value()),
            _ => panic!("unknown argument {arg}"),
        }
    }
    args
}

#[derive(Default)]
struct Counts {
    curves: u64,
    stage1: u64,
    stage2: u64,
    setup: u64,
}

impl Counts {
    fn found(&self) -> u64 {
        self.stage1 + self.stage2 + self.setup
    }
}

/// Runs `curves` curves on each of `numbers` composites `p * q`, with `p` of `digits` digits.
fn measure(digits: u32, b1: usize, b2: usize, numbers: u64, curves: u64) -> Counts {
    let k = stage1_multiplier(b1);

    // (n, sigma) of every curve, drawn up front so results don't depend on thread scheduling.
    let mut tasks = Vec::new();
    for i in 0..numbers {
        let seed = SEED + u64::from(digits) * 1_000_000 + i;
        let n = prime_digits(digits, seed) * prime_digits(COFACTOR_DIGITS, seed + 500_000);
        let mut rand = RandState::new();
        rand.seed(&Integer::from(seed));
        let range = Integer::from(&n - 7);
        for _ in 0..curves {
            let sigma = Integer::from(6) + range.clone().random_below(&mut rand);
            tasks.push((n.clone(), sigma));
        }
    }

    let next = AtomicUsize::new(0);
    let threads = thread::available_parallelism().map_or(1, |n| n.get());
    let results: Vec<Counts> = thread::scope(|s| {
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                s.spawn(|| {
                    let mut counts = Counts::default();
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        let Some((n, sigma)) = tasks.get(i) else {
                            return counts;
                        };
                        counts.curves += 1;
                        match run_curve(n, sigma, &k, b1, b2) {
                            CurveOutcome::Stage1(_) => counts.stage1 += 1,
                            CurveOutcome::Stage2(_) => counts.stage2 += 1,
                            // `Setup` may return `n` itself: only count proper factors.
                            CurveOutcome::Setup(g) if &g != n => counts.setup += 1,
                            CurveOutcome::Setup(_) | CurveOutcome::Failed => {}
                        }
                    }
                })
            })
            .collect();
        workers.into_iter().map(|w| w.join().unwrap()).collect()
    });

    results
        .into_iter()
        .fold(Counts::default(), |acc, c| Counts {
            curves: acc.curves + c.curves,
            stage1: acc.stage1 + c.stage1,
            stage2: acc.stage2 + c.stage2,
            setup: acc.setup + c.setup,
        })
}

/// 95% Wilson score interval of a success rate.
fn wilson(successes: u64, trials: u64) -> (f64, f64) {
    let (s, n, z) = (successes as f64, trials as f64, 1.96_f64);
    let p = s / n;
    let denominator = 1.0 + z * z / n;
    let center = (p + z * z / (2.0 * n)) / denominator;
    let half = z * (p * (1.0 - p) / n + z * z / (4.0 * n * n)).sqrt() / denominator;
    ((center - half).max(0.0), (center + half).min(1.0))
}

fn main() {
    let args = parse_args();
    let mut results = Vec::new();

    for &digits in &args.sizes {
        let &(_, b1, b2) = GMP_ECM_BOUNDS
            .iter()
            .find(|(d, _, _)| *d == digits)
            .unwrap_or_else(|| panic!("no bounds for {digits}-digit factors"));

        let start = Instant::now();
        let counts = measure(digits, b1, b2, args.numbers, args.curves);
        let elapsed = start.elapsed().as_secs_f64();

        let found = counts.found();
        let (low, high) = wilson(found, counts.curves);
        // With no success, the expected number of curves is unknown: report the number of
        // curves tried, which is a lower bound.
        let expected = if found == 0 {
            counts.curves as f64
        } else {
            counts.curves as f64 / found as f64
        };
        let interval = if low == 0.0 {
            format!("{:.0}..inf", 1.0 / high)
        } else {
            format!("{:.0}..{:.0}", 1.0 / high, 1.0 / low)
        };
        let reference = GMP_ECM_EXPECTED_CURVES
            .iter()
            .find(|(d, _)| *d == digits)
            .map_or(0, |(_, c)| *c);

        println!(
            "{digits}-digit factor (B1={b1}, B2={b2}): found {found}/{} curves \
             (stage 1: {}, stage 2: {}, setup: {}), expected curves {expected:.1} \
             (95% CI {interval}), GMP-ECM ~{reference}, {elapsed:.1}s",
            counts.curves, counts.stage1, counts.stage2, counts.setup
        );

        results.push(json!({
            "name": format!("expected curves, {digits}-digit factor (B1={b1}, B2={b2})"),
            "unit": "curves",
            "value": expected,
            "extra": format!(
                "found {found}/{} curves (stage 1: {}, stage 2: {}, setup: {})\n\
                 95% CI: {interval} curves{}\nGMP-ECM reference: ~{reference} curves",
                counts.curves,
                counts.stage1,
                counts.stage2,
                counts.setup,
                if found == 0 { "\nno factor found: value is a lower bound" } else { "" },
            ),
        }));
    }

    if let Some(path) = &args.json {
        std::fs::write(path, serde_json::to_string_pretty(&results).unwrap()).unwrap();
    }

    if let Some(path) = &args.compare {
        let table = comparison(&std::fs::read_to_string(path).unwrap(), &results);
        println!("\n{table}");
        if let Some(markdown) = &args.markdown {
            std::fs::write(markdown, table).unwrap();
        }
    }
}

/// 95% confidence interval of the expected curves, from the `extra` field of a result.
fn interval(result: &Value) -> (f64, f64) {
    let extra = result["extra"].as_str().unwrap();
    let start = extra.find("95% CI: ").unwrap() + "95% CI: ".len();
    let end = start + extra[start..].find(" curves").unwrap();
    let (low, high) = extra[start..end].split_once("..").unwrap();
    (low.parse().unwrap(), high.parse().unwrap_or(f64::INFINITY))
}

/// Markdown table comparing `results` with a previous JSON output (fewer curves is better).
fn comparison(base: &str, results: &[Value]) -> String {
    let base: Vec<Value> = serde_json::from_str(base).unwrap();
    let mut table = String::from(
        "### Curve success rate\n\n| | Factor | Base | PR | Change |\n|---|---|---:|---:|---:|\n",
    );
    for result in results {
        let name = result["name"].as_str().unwrap();
        let name = name.strip_prefix("expected curves, ").unwrap_or(name);
        let value = result["value"].as_f64().unwrap();
        let (low, high) = interval(result);
        let row = match base.iter().find(|b| b["name"] == result["name"]) {
            Some(b) => {
                let old = b["value"].as_f64().unwrap();
                let (old_low, old_high) = interval(b);
                let mark = if high < old_low {
                    "🟢"
                } else if low > old_high {
                    "🔴"
                } else {
                    "⚪"
                };
                let change = (value - old) / old * 100.0;
                format!(
                    "| {mark} | {name} | {old:.1} ({old_low}..{old_high}) | {value:.1} ({low}..{high}) \
                     | {change:+.1}% |"
                )
            }
            None => format!("| | {name} | - | {value:.1} ({low}..{high}) | new |"),
        };
        table.push_str(&row);
        table.push('\n');
    }
    table
}
