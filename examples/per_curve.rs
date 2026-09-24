//! Wall-clock time of stage 1 and stage 2 of one curve, on the rows of the comparison with
//! GMP-ECM (`scripts/compare_gmp_ecm.sh`): balanced semiprimes of 256 and 1024 bits, `B1`/`B2`
//! from 50k/12.7M to 1M/1e9, curve `-param 2` with `sigma = 1234567`.
//!
//! ```text
//! cargo run --release --features bench --example per_curve -- [--reps 3] [--bits 256,1024]
//! cargo run --release --features bench --example per_curve -- --numbers
//! ```
//!
//! - Default: one line per row, `bits B1 B2 plan stage1_ms stage2_ms` (medians over `--reps`
//!   runs; the stage 2 plan, chosen by the cost model, is built beforehand and not measured).
//! - `--numbers`: prints `bits B1 B2 sigma n` per row instead, to run the same curves elsewhere
//!   (`echo n | ecm -c 1 -param 2 -sigma 2:sigma B1 B2`).
//!
//! Build with `--profile bench-lto` instead of `--release` for the optimizations of a final
//! build (`lto = "fat"`, `codegen-units = 1`).

#[path = "../benches/common/mod.rs"]
mod common;

use common::{PER_CURVE_ROWS, SEED, SIGMA, curve_point, semiprime_bits};
use ecm::bench::{Param, Stage2Plan, stage1, stage1_multiplier, stage2};
use std::time::Instant;

fn main() {
    let mut reps = 3;
    let mut bits: Option<Vec<u32>> = None;
    let mut numbers = false;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--reps" => reps = args.next().expect("missing --reps value").parse().unwrap(),
            "--bits" => {
                let list = args.next().expect("missing --bits value");
                bits = Some(list.split(',').map(|b| b.parse().unwrap()).collect());
            }
            "--numbers" => numbers = true,
            _ => panic!("unknown argument {arg}"),
        }
    }

    let rows = PER_CURVE_ROWS
        .into_iter()
        .filter(|(b, _, _)| bits.as_ref().is_none_or(|bits| bits.contains(b)));

    if !numbers {
        println!("# bits B1 B2 plan stage1_ms stage2_ms");
    }
    for (bits, b1, b2) in rows {
        let n = semiprime_bits(bits, SEED);
        if numbers {
            println!("{bits} {b1} {b2} {SIGMA} {n}");
            continue;
        }
        let p = curve_point(&n, Param::Batch2);
        let k = stage1_multiplier(b1);
        let plan = Stage2Plan::new(&n, b1, b2);
        let (mut t1, mut t2) = (Vec::new(), Vec::new());
        for _ in 0..reps {
            let start = Instant::now();
            let q = stage1(&p, &k);
            t1.push(start.elapsed().as_secs_f64() * 1e3);
            assert_eq!(q.z.clone().gcd(&n), 1, "stage 1 found a factor");
            let start = Instant::now();
            let g = stage2(&q, &plan);
            t2.push(start.elapsed().as_secs_f64() * 1e3);
            assert_eq!(g, 1, "stage 2 found a factor");
        }
        let kind = match plan {
            Stage2Plan::Pairs(_) => "pairs",
            Stage2Plan::Poly(_) => "poly",
        };
        println!(
            "{bits} {b1} {b2} {kind} {:.1} {:.1}",
            median(&mut t1),
            median(&mut t2)
        );
    }
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}
