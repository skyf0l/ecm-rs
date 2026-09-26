//! Curves in parallel (`Factorizer::threads`): same results and events as with one thread,
//! interruptions, no deadlock.

use ecm::{Base2Mode, Error, Event, Factorization, Factorizer, Method, Param};
use rug::{Integer, ops::Pow, rand::RandState};
use std::{
    collections::HashMap,
    ops::ControlFlow,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
    thread,
    time::{Duration, Instant},
};

/// The data of an event, but the durations.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Owned {
    TrialDivision(Integer),
    Base2(Integer, i64),
    Pm1(Integer, usize, Option<usize>),
    Level(Integer, Option<u32>, usize, usize, Option<usize>, usize),
    Curve(Integer, Param, Integer, usize),
    Factor(Integer, Integer, Method),
    Prime(Integer, usize),
    Other,
}

fn owned(event: &Event<'_>) -> Owned {
    match *event {
        Event::TrialDivision { cofactor, .. } => Owned::TrialDivision(cofactor.clone()),
        Event::Base2 { n, k, .. } => Owned::Base2(n.clone(), k),
        Event::Pm1 { n, b1, b2, .. } => Owned::Pm1(n.clone(), b1, b2),
        Event::Level {
            n,
            digits,
            b1,
            b2,
            curves,
            done,
            ..
        } => Owned::Level(n.clone(), digits, b1, b2, curves, done),
        Event::Curve {
            n,
            param,
            sigma,
            index,
            stage1,
            ..
        } => {
            assert!(stage1 > Duration::ZERO);
            Owned::Curve(n.clone(), param, sigma.clone(), index)
        }
        Event::Factor {
            n, factor, method, ..
        } => Owned::Factor(n.clone(), factor.clone(), method),
        Event::Prime { p, exponent, .. } => Owned::Prime(p.clone(), exponent),
        _ => Owned::Other,
    }
}

/// Factors `n` with `factorizer`, recording the events, until `stop` returns `true`.
fn record(
    factorizer: Factorizer,
    n: &Integer,
    mut stop: impl FnMut(&Owned) -> bool,
) -> (Factorization, Vec<Owned>) {
    let mut events = Vec::new();
    let mut stopped = false;
    let result = factorizer
        .on_event(|event| {
            assert!(!stopped, "event after the interruption");
            events.push(owned(event));
            stopped = stop(events.last().unwrap());
            if stopped {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        })
        .factor_partial(n);
    (result, events)
}

fn product(factors: &HashMap<Integer, usize>) -> Integer {
    factors
        .iter()
        .map(|(p, &e)| p.clone().pow(e as u32))
        .product()
}

/// The primes and the unfactored parts multiply back to `n`.
fn check_partial(n: &Integer, result: &Factorization) {
    let unfactored: Integer = result
        .unfactored
        .iter()
        .map(|(m, e)| m.clone().pow(*e as u32))
        .product();
    assert_eq!(product(&result.primes) * unfactored, *n);
}

/// Products of 2 or 3 random primes of 7 to 20 digits.
fn numbers(count: usize) -> Vec<Integer> {
    let mut rand = RandState::new();
    rand.seed(&Integer::from(2024));
    (0..count)
        .map(|i| {
            let primes = 2 + i % 2;
            (0..primes)
                .map(|_| {
                    let bits = 24 + rand.bits(6) % 43;
                    Integer::from(Integer::random_bits(bits, &mut rand)).next_prime()
                })
                .product()
        })
        .collect()
}

const THREADS: [usize; 3] = [2, 4, 8];

#[test]
fn same_results_and_events() {
    let mut many_curves = 0;
    for (i, n) in numbers(12).iter().enumerate() {
        for seed in 0..3 {
            let factorizer = Factorizer::new().seed(seed + 10 * i as u64);
            let (single, events) = record(factorizer.clone(), n, |_| false);
            let single = single.into_result().unwrap();
            assert_eq!(product(&single), *n);
            let curves = events
                .iter()
                .filter(|e| matches!(e, Owned::Curve(..)))
                .count();
            many_curves += usize::from(curves >= 10);
            for threads in THREADS {
                let (result, parallel) = record(factorizer.clone().threads(threads), n, |_| false);
                assert_eq!(result.into_result().unwrap(), single, "{n} {threads}");
                assert_eq!(parallel, events, "{n} seed {seed} threads {threads}");
                // Without events.
                let plain = factorizer.clone().threads(threads).factor(n).unwrap();
                assert_eq!(plain, single);
            }
        }
    }
    assert!(many_curves >= 12, "{many_curves}");
}

/// A prime `p` of at least `bits` bits with `p - 1` a product of distinct primes below 250.
fn smooth_prime(bits: u32, seed: u32) -> Integer {
    let mut rand = RandState::new();
    rand.seed(&Integer::from(seed));
    loop {
        let mut small: Vec<u32> = (3..250u32)
            .filter(|&q| (2..q).all(|d| q % d != 0))
            .collect();
        let mut p = Integer::from(2);
        while p.significant_bits() < bits {
            p *= small.swap_remove(rand.below(small.len() as u32) as usize);
        }
        p += 1u32;
        if p.is_probably_prime(30) != rug::integer::IsPrime::No {
            return p;
        }
    }
}

#[test]
fn same_results_special_cases() {
    let q = Integer::from(10u64.pow(19)).next_prime();
    let (p1, p2) = (smooth_prime(40, 3), smooth_prime(45, 4));
    let cases = [
        // P-1 finds the product of two factors (composite, searched again), or all the factors
        // at once (not run again).
        (Factorizer::new(), Integer::from(&p1 * &p2) * &q),
        (Factorizer::new(), Integer::from(&p1 * &p2)),
        // P-1 finds a 28-digit factor at the second level (stage 2), with curves running.
        (
            Factorizer::new(),
            Integer::from_str_radix("1326262092842391910564284053", 10).unwrap()
                * (Integer::from(10).pow(60) + 7u32).next_prime(),
        ),
        // Special reduction modulo 2^256 + 1 (= p16 * p62).
        (
            Factorizer::new().base2(Base2Mode::Force(256)),
            (Integer::from(1) << 256u32) + 1u32,
        ),
        // Several factors of each size, with cofactors resuming the levels.
        (
            Factorizer::new().seed(5),
            numbers(6).into_iter().product::<Integer>(),
        ),
    ];
    for (factorizer, n) in cases {
        let (single, events) = record(factorizer.clone(), &n, |_| false);
        assert_eq!(product(&single.clone().into_result().unwrap()), n);
        for threads in THREADS {
            let (result, parallel) = record(factorizer.clone().threads(threads), &n, |_| false);
            assert_eq!(result, single, "{n} {threads}");
            assert_eq!(parallel, events, "{n} threads {threads}");
        }
    }
}

#[test]
fn fixed_bounds_and_find_factor() {
    let p20 = Integer::from(10u64.pow(19)).next_prime();
    let n = Integer::from(&p20 * 1_000_000_007u32) * Integer::from(10u64.pow(18)).next_prime();
    let cases = [
        Factorizer::new().b1(2_000),
        Factorizer::new().b1(11_000).param(Param::Suyama).seed(3),
        Factorizer::new().b1(11_000).param(Param::Square),
        Factorizer::new()
            .b1(11_000)
            .sigma(Integer::from(1_000_000))
            .curves(40),
        // Fails: fewer curves than needed.
        Factorizer::new().b1(1_000).b2(10_000).curves(3),
    ];
    for factorizer in cases {
        let (single, events) = record(factorizer.clone(), &n, |_| false);
        let one = factorizer.clone().find_factor(&n);
        for threads in THREADS {
            let (result, parallel) = record(factorizer.clone().threads(threads), &n, |_| false);
            assert_eq!(result, single);
            assert_eq!(parallel, events);
            assert_eq!(factorizer.clone().threads(threads).find_factor(&n), one);
        }
    }
    // All the available threads.
    let factors = Factorizer::new().threads(0).factor(&n).unwrap();
    assert_eq!(product(&factors), n);
}

/// 60 digits: a 25-digit factor, long to find.
fn hard() -> Integer {
    Integer::from_str_radix("1000000000000000000000007", 10).unwrap()
        * Integer::from_str_radix("100000000000000000000000000000000067", 10).unwrap()
}

#[test]
fn callback_interruptions() {
    // Interrupted at a curve (with the curves after it running): no deadlock, no event after,
    // the events of one thread up to there.
    let n = hard() * 12u32;
    let (_, all) = record(Factorizer::new(), &n, |e| {
        matches!(e, Owned::Curve(_, _, _, 60))
    });
    for threads in THREADS {
        for stop_at in [1, 2, 7, 30, 60] {
            let start = Instant::now();
            let (result, events) = record(
                Factorizer::new().threads(threads),
                &n,
                |e| matches!(e, Owned::Curve(_, _, _, i) if *i == stop_at),
            );
            assert!(start.elapsed() < Duration::from_secs(2));
            assert_eq!(result.error, Some(Error::Interrupted));
            check_partial(&n, &result);
            assert_eq!(events, all[..events.len()]);
        }
    }

    // Interrupted at the factor found by a curve: the factor is kept.
    let p = Integer::from(10u64.pow(11)).next_prime();
    let n = Integer::from(&p * 1_000_000_007u32) * Integer::from(10u64.pow(19)).next_prime();
    let stop = |e: &Owned| {
        matches!(
            e,
            Owned::Factor(_, _, Method::EcmStage1 | Method::EcmStage2)
        )
    };
    let (single, events) = record(Factorizer::new(), &n, stop);
    for threads in THREADS {
        let (result, parallel) = record(Factorizer::new().threads(threads), &n, stop);
        assert_eq!(result, single);
        assert_eq!(parallel, events);
    }
}

fn rsa1024() -> Integer {
    let p = (Integer::from(1) << 511u32) + 1_000u32;
    let q = (Integer::from(3) << 510u32) + 1_000u32;
    p.next_prime() * q.next_prime()
}

#[test]
fn interrupt_flag_and_timeout() {
    let n = rsa1024();
    for threads in [2, 4] {
        for (factorizer, delay) in [
            // Stage 1 of the curves (about 5 s each).
            (Factorizer::new().b1(1_000_000), 300),
            // Stage 2 of the curves.
            (Factorizer::new().b1(11_000).b2(1_045_563_762), 600),
            // The first level: P-1, then curves.
            (Factorizer::new(), 400),
        ] {
            let flag = Arc::new(AtomicBool::new(false));
            let setter = {
                let flag = Arc::clone(&flag);
                thread::spawn(move || {
                    thread::sleep(Duration::from_millis(delay));
                    flag.store(true, Ordering::Relaxed);
                    Instant::now()
                })
            };
            let result = factorizer
                .threads(threads)
                .interrupt_flag(flag)
                .factor_partial(&n);
            let latency = setter.join().unwrap().elapsed();
            assert!(latency < Duration::from_millis(100), "{latency:?}");
            assert_eq!(result.error, Some(Error::Interrupted));
            assert_eq!(result.unfactored, [(n.clone(), 1)]);
        }

        let n = Integer::from(&n * 7u32);
        let start = Instant::now();
        let result = Factorizer::new()
            .threads(threads)
            .timeout(Duration::from_millis(200))
            .factor_partial(&n);
        let elapsed = start.elapsed();
        assert!(elapsed < Duration::from_millis(300), "{elapsed:?}");
        assert_eq!(result.error, Some(Error::Interrupted));
        assert_eq!(result.primes, HashMap::from([(Integer::from(7), 1)]));
    }
}
