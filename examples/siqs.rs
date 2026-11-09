//! Wall-clock time of the self-initializing quadratic sieve on one thread, on balanced
//! semiprimes (two primes of half the size), with its statistics: for the tuning of its
//! parameters.
//!
//! ```text
//! cargo run --release --features bench --example siqs -- [--sizes 50,60,70] [--count 3]
//! ```

use ecm::bench::siqs_factor_stats;
use rug::{Integer, rand::RandState};
use std::time::Instant;

/// A semiprime of `digits` digits, product of two primes of about half its size.
fn semiprime(digits: u32, index: u64) -> Integer {
    let mut rand = RandState::new();
    rand.seed(&Integer::from(1000 * u64::from(digits) + index));
    loop {
        let half = digits / 2;
        let low = Integer::from(Integer::u_pow_u(10, half - 1));
        let p = (Integer::from(low.random_below_ref(&mut rand)) * 9u32 + &low).next_prime();
        let low = Integer::from(Integer::u_pow_u(10, digits - half - 1));
        let q = (Integer::from(low.random_below_ref(&mut rand)) * 9u32 + &low).next_prime();
        let n = Integer::from(&p * &q);
        if n.to_string().len() == digits as usize {
            return n;
        }
    }
}

fn main() {
    let mut sizes = vec![40, 50, 60];
    let mut count = 1;
    let mut print = false;
    let mut skip = 0;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--sizes" => {
                let list = args.next().expect("missing --sizes value");
                sizes = list.split(',').map(|d| d.parse().unwrap()).collect();
            }
            "--count" => count = args.next().expect("missing --count value").parse().unwrap(),
            "--numbers" => print = true,
            "--skip" => skip = args.next().expect("missing --skip value").parse().unwrap(),
            _ => panic!("unknown argument {arg}"),
        }
    }
    if print {
        for &d in &sizes {
            for index in 0..count {
                println!("{}", semiprime(d, index));
            }
        }
        return;
    }
    println!(
        "digits k fb interval s fulls combined partials polys matrix deps sieve_s la_s total_s"
    );
    for &d in &sizes {
        for index in skip..count {
            let n = semiprime(d, index);
            let start = Instant::now();
            let stats = siqs_factor_stats(&n, 0).expect("SIQS applies");
            let total = start.elapsed();
            assert_eq!(stats.parts.iter().product::<Integer>(), n);
            println!(
                "{d} {} {} {} {} {} {} {} {} {}x{} {} {:.3} {:.3} {:.3}",
                stats.multiplier,
                stats.factor_base,
                stats.interval,
                stats.a_primes,
                stats.fulls,
                stats.combined,
                stats.partials,
                stats.polynomials,
                stats.matrix.0,
                stats.matrix.1,
                stats.matrix.2,
                stats.sieve.as_secs_f64(),
                stats.linear_algebra.as_secs_f64(),
                total.as_secs_f64()
            );
        }
    }
}
