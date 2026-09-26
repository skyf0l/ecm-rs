//! The public `Factorizer` API: events, cancellation, fixed bounds, determinism.

use ecm::{
    Algorithm, Base2Mode, Error, Event, Factorization, Factorizer, Method, Param, ecm,
    ecm_with_params,
};
use rug::{Integer, ops::Pow};
use std::{
    collections::HashMap,
    ops::ControlFlow,
    str::FromStr,
    sync::{
        Arc,
        atomic::{AtomicBool, Ordering},
    },
    thread,
    time::{Duration, Instant},
};

fn int(s: &str) -> Integer {
    Integer::from_str(s).unwrap()
}

/// An owned copy of the data of an event.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Owned {
    TrialDivision(HashMap<Integer, usize>, Integer),
    Base2(Integer, i64),
    Pm1(Integer, usize, Option<usize>),
    Pp1(Integer, usize, Option<usize>),
    Level(Integer, Option<u32>, usize, usize, Option<usize>, usize),
    Curve(Integer, Param, Integer, usize),
    Factor(Integer, Integer, Method),
    Prime(Integer, usize),
}

fn owned(event: &Event<'_>) -> Owned {
    match *event {
        Event::TrialDivision {
            factors, cofactor, ..
        } => Owned::TrialDivision(factors.clone(), cofactor.clone()),
        Event::Base2 { n, k, .. } => Owned::Base2(n.clone(), k),
        Event::Pm1 {
            n,
            b1,
            b2,
            stage1,
            stage2,
            ..
        } => {
            assert!(b2.is_some() || stage2 == Duration::ZERO);
            assert!(stage1 > Duration::ZERO);
            Owned::Pm1(n.clone(), b1, b2)
        }
        Event::Pp1 {
            n,
            b1,
            b2,
            stage1,
            stage2,
            ..
        } => {
            assert!(b2.is_some() || stage2 == Duration::ZERO);
            assert!(stage1 > Duration::ZERO);
            Owned::Pp1(n.clone(), b1, b2)
        }
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
        _ => panic!("unknown event"),
    }
}

/// Factors `n` with `factorizer`, recording the events.
fn record<F>(factorizer: Factorizer, n: &Integer, stop: F) -> (Factorization, Vec<Owned>)
where
    F: FnMut(&Owned) -> bool,
{
    let (result, events, _) = record_timed(factorizer, n, stop);
    (result, events)
}

/// [`record`], with the time from the interruption (the callback returning `Break`) to the
/// return of the factorization.
fn record_timed<F>(
    factorizer: Factorizer,
    n: &Integer,
    mut stop: F,
) -> (Factorization, Vec<Owned>, Option<Duration>)
where
    F: FnMut(&Owned) -> bool,
{
    let mut events = Vec::new();
    let mut stopped = None;
    let result = factorizer
        .on_event(|event| {
            assert!(stopped.is_none(), "event after the interruption");
            events.push(owned(event));
            if stop(events.last().unwrap()) {
                stopped = Some(Instant::now());
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        })
        .factor_partial(n);
    (result, events, stopped.map(|at| at.elapsed()))
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

/// Checks the consistency of the events of the complete factorization of `n`.
fn check_events(n: &Integer, factors: &HashMap<Integer, usize>, events: &[Owned]) {
    assert!(matches!(events[0], Owned::TrialDivision(..)));
    let mut primes = HashMap::new();
    let mut level = None;
    let mut last_curve = 0;
    for event in events {
        match event {
            Owned::Prime(p, e) => *primes.entry(p.clone()).or_insert(0) += e,
            Owned::Factor(m, g, _) => {
                assert!(*g > 1 && g < m && m.is_divisible(g), "{m} {g}");
                assert!(n.is_divisible(m));
            }
            Owned::Level(m, digits, b1, b2, curves, done) => {
                assert!(n.is_divisible(m));
                assert!(digits.is_some() && b1 < b2 && done < &curves.unwrap());
                level = Some(m.clone());
                last_curve = *done;
            }
            Owned::Curve(m, param, sigma, index) => {
                assert_eq!(Some(m), level.as_ref());
                assert_eq!(*param, Param::Batch2);
                assert!(*sigma >= 2 && *sigma < u64::MAX);
                assert_eq!(*index, last_curve + 1);
                last_curve = *index;
            }
            Owned::Pm1(m, b1, b2) | Owned::Pp1(m, b1, b2) => {
                assert!(n.is_divisible(m));
                assert!(b2.is_none_or(|b2| b2 > *b1));
                level = None;
            }
            Owned::Base2(m, k) => {
                assert!(n.is_divisible(m));
                let power = Integer::from(1) << (k.unsigned_abs() as u32);
                let form = if *k > 0 { power + 1u32 } else { power - 1u32 };
                assert!(form.is_divisible(m));
            }
            Owned::TrialDivision(..) => {}
        }
    }
    assert_eq!(&primes, factors);
    assert_eq!(product(factors), *n);
}

/// 43 * 2634823 * 407485517 * 1000000007^2 * 99999999977 * 1000000000039
fn sample() -> Integer {
    int("46167045131415113") * int("1000000007").pow(2) * int("99999999977") * int("1000000000039")
}

#[test]
fn events_are_consistent() {
    let n = sample();
    let (result, events) = record(Factorizer::new(), &n, |_| false);
    let factors = result.into_result().unwrap();
    assert_eq!(factors, ecm(&n).unwrap());
    check_events(&n, &factors, &events);
    assert!(events.iter().any(|e| matches!(e, Owned::Curve(..))));
    assert!(events.iter().any(|e| matches!(e, Owned::Pm1(..))));

    // The expected curves grow with the size of the factors.
    let mut curves: Vec<_> = events
        .iter()
        .filter_map(|e| match e {
            Owned::Level(_, Some(digits), _, _, Some(curves), _) => Some((*digits, *curves)),
            _ => None,
        })
        .collect();
    curves.sort_unstable();
    assert!(curves.windows(2).all(|w| w[0].1 <= w[1].1), "{curves:?}");
    // GMP-ECM's expected curves for 15-digit factors with B1 = 2000 (up to the B2 covered).
    if let Some(&(_, c)) = curves.iter().find(|(d, _)| *d == 15) {
        assert!((20..40).contains(&c), "{c}");
    }
}

#[test]
fn events_of_small_numbers() {
    for n in [1u64, 2, 17, 1024, 3 * 5 * 7, 65_537 * 65_539] {
        let n = Integer::from(n);
        let (result, events) = record(Factorizer::new(), &n, |_| false);
        check_events(&n, &result.into_result().unwrap(), &events);
    }
    // A perfect power of a composite.
    let n = (Integer::from(4_009_823u32) * 99_476_569u32).pow(3);
    let (result, events) = record(Factorizer::new(), &n, |_| false);
    check_events(&n, &result.into_result().unwrap(), &events);
    assert!(
        events
            .iter()
            .any(|e| matches!(e, Owned::Factor(_, _, Method::PerfectPower)))
    );
}

#[test]
fn callback_does_not_change_results() {
    let n = sample();
    for seed in 0..3 {
        let plain = Factorizer::new().seed(seed).factor(&n).unwrap();
        let (with_events, events) = record(Factorizer::new().seed(seed), &n, |_| false);
        assert_eq!(with_events.into_result().unwrap(), plain);
        // Same seed, same curves.
        let (_, again) = record(Factorizer::new().seed(seed), &n, |_| false);
        let sigmas = |events: &[Owned]| -> Vec<Integer> {
            events
                .iter()
                .filter_map(|e| match e {
                    Owned::Curve(_, _, sigma, _) => Some(sigma.clone()),
                    _ => None,
                })
                .collect()
        };
        assert_eq!(sigmas(&events), sigmas(&again));
    }
    assert_eq!(
        Factorizer::new().factor(&n).unwrap(),
        ecm(&n).unwrap(),
        "the default seed is the one of ecm()"
    );
}

/// 60 digits: a 25-digit factor, long to find.
fn hard() -> Integer {
    int("1000000000000000000000007") * int("100000000000000000000000000000000067")
}

/// Interrupts at the first event matching `stop`: the factorization returns at once, with the
/// parts found so far.
fn interrupt(n: &Integer, factorizer: Factorizer, stop: impl Fn(&Owned) -> bool) -> Vec<Owned> {
    let (result, events, latency) = record_timed(factorizer.clone(), n, stop);
    let latency = latency.unwrap();
    assert!(latency < Duration::from_secs(1), "{latency:?}");
    assert_eq!(result.error, Some(Error::Interrupted));
    check_partial(n, &result);
    // factor() returns the error.
    let (mut calls, stop_at) = (0, events.len());
    let error = factorizer
        .on_event(|_| {
            calls += 1;
            if calls == stop_at {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        })
        .factor(n);
    assert_eq!(error, Err(Error::Interrupted));
    events
}

#[test]
fn interruptions() {
    let n = hard() * 12u32;
    // Trial division.
    let events = interrupt(&n, Factorizer::new(), |e| {
        matches!(e, Owned::TrialDivision(..))
    });
    assert_eq!(events.len(), 1);
    // P-1.
    let events = interrupt(&n, Factorizer::new(), |e| matches!(e, Owned::Pm1(..)));
    assert!(matches!(events.last(), Some(Owned::Pm1(..))));
    // A level, a curve, the 20th curve of the second level.
    interrupt(&n, Factorizer::new(), |e| matches!(e, Owned::Level(..)));
    interrupt(&n, Factorizer::new(), |e| matches!(e, Owned::Curve(..)));
    interrupt(&n, Factorizer::new(), |e| {
        matches!(e, Owned::Curve(_, _, _, 20))
    });
    // Fixed bounds, P-1 only.
    interrupt(&n, Factorizer::new().b1(11_000), |e| {
        matches!(e, Owned::Curve(_, _, _, 5))
    });
    interrupt(
        &n,
        Factorizer::new().b1(11_000).algorithm(Algorithm::Pm1),
        |e| matches!(e, Owned::Pm1(..)),
    );

    // Interrupted while a factor was found: the factor is kept.
    let n = sample();
    let (result, events) = record(Factorizer::new(), &n, |e| {
        matches!(
            e,
            Owned::Factor(_, _, Method::EcmStage1 | Method::EcmStage2)
        )
    });
    assert_eq!(result.error, Some(Error::Interrupted));
    check_partial(&n, &result);
    let Some(Owned::Factor(_, factor, _)) = events.last() else {
        panic!()
    };
    assert!(
        result.primes.contains_key(factor) || result.unfactored.iter().any(|(m, _)| m == factor)
    );
}

#[test]
fn find_factor() {
    let mut events = Vec::new();
    let mut factorizer = Factorizer::new().on_event(|e| {
        events.push(owned(e));
        ControlFlow::Continue(())
    });
    assert_eq!(factorizer.find_factor(&(hard() * 7u32)).unwrap(), 7);
    assert_eq!(factorizer.find_factor(&Integer::from(3 * 3)).unwrap(), 3);
    assert_eq!(
        factorizer
            .find_factor(&Integer::from(65_537u64 * 65_537))
            .unwrap(),
        65_537
    );
    assert_eq!(
        factorizer.find_factor(&int("1000000007")),
        Err(Error::NumberIsPrime)
    );
    let n = int("398883434337287");
    let g = factorizer.find_factor(&n).unwrap();
    assert!(g == 4_009_823 || g == 99_476_569);
    drop(factorizer);
    let methods: Vec<_> = events
        .iter()
        .filter_map(|e| match e {
            Owned::Factor(_, _, method) => Some(*method),
            _ => None,
        })
        .collect();
    // 65537 is above the bound of trial division.
    assert_eq!(
        methods[..3],
        [
            Method::TrialDivision,
            Method::TrialDivision,
            Method::PerfectPower
        ]
    );
    assert!(matches!(
        methods[3],
        Method::EcmStage1 | Method::EcmStage2 | Method::Pm1Stage1 | Method::Pm1Stage2
    ));
}

#[test]
fn fixed_bounds_as_ecm_with_params() {
    for (n, b1, b2, curves) in [
        (sample(), 2_000, 147_396, 1000),
        (int("398883434337287"), 300, 9_846, 100),
        (
            int("4269021180054189416198169786894227"),
            11_000,
            1_873_422,
            500,
        ),
        (int("398883434337287"), 6, 8, 3),
    ] {
        for seed in 0..3 {
            let expected = ecm_with_params(&n, b1, b2, curves, seed as usize);
            let factorizer = Factorizer::new().seed(seed).b1(b1).b2(b2).curves(curves);
            let (result, events) = record(factorizer, &n, |_| false);
            assert_eq!(result.clone().into_result(), expected);
            check_partial(&n, &result);
            for event in &events {
                match event {
                    Owned::Level(_, digits, lb1, lb2, lcurves, done) => {
                        assert_eq!(
                            (*digits, *lb1, *lcurves, *done),
                            (None, b1, Some(curves), 0)
                        );
                        assert!(*lb2 >= b2);
                    }
                    Owned::Pm1(..) => panic!("no P-1 with fixed bounds"),
                    _ => {}
                }
            }
        }
    }
}

#[test]
fn fixed_sigma() {
    // `echo 398883434337287 | ecm -param 1 -sigma 1:3 300 0` finds 4009823 in stage 1, and
    // nothing with sigma = 2.
    let n = int("398883434337287");
    let fixed = || Factorizer::new().param(Param::Square).b1(300).b2(300);
    let (result, events) = record(fixed().sigma(Integer::from(2)).curves(1), &n, |_| false);
    assert_eq!(result.error, Some(Error::ECMFailed));
    assert_eq!(result.unfactored, [(n.clone(), 1)]);
    assert_eq!(events.len(), 3);
    // The next curve takes sigma = 3.
    let (result, events) = record(fixed().sigma(Integer::from(2)).curves(2), &n, |_| false);
    assert_eq!(result.error, None);
    let curves: Vec<_> = events
        .iter()
        .filter_map(|e| match e {
            Owned::Curve(_, Param::Square, sigma, index) => Some((sigma.clone(), *index)),
            _ => None,
        })
        .collect();
    assert_eq!(curves, [(Integer::from(2), 1), (Integer::from(3), 2)]);
    assert!(events.contains(&Owned::Factor(
        n.clone(),
        Integer::from(4_009_823),
        Method::EcmStage1
    )));
    assert_eq!(
        fixed().sigma(Integer::from(3)).find_factor(&n),
        Ok(Integer::from(4_009_823))
    );
}

#[test]
fn pm1_only() {
    // p - 1 = 2^2 * 5743 * 6037 * 30259 * 30949 * 34039 * 300007.
    let p = int("1326262092842391910564284053");
    let n = Integer::from(10).pow(100) + 267u32;
    let n = Integer::from(&p * &n);
    let pm1 = |b1, b2| Factorizer::new().algorithm(Algorithm::Pm1).b1(b1).b2(b2);
    let (result, events) = record(pm1(40_000, 400_000), &n, |_| false);
    assert_eq!(result.error, None);
    assert!(events.contains(&Owned::Factor(n.clone(), p.clone(), Method::Pm1Stage2)));
    assert!(!events.iter().any(|e| matches!(e, Owned::Curve(..))));
    assert_eq!(pm1(400_000, 400_000).find_factor(&n), Ok(p.clone()));
    assert_eq!(pm1(40_000, 200_000).find_factor(&n), Err(Error::ECMFailed));
    // Without fixed bounds: P-1 of the levels only.
    let (result, events) = record(Factorizer::new().algorithm(Algorithm::Pm1), &n, |_| false);
    assert_eq!(result.error, None);
    assert!(
        !events
            .iter()
            .any(|e| matches!(e, Owned::Curve(..) | Owned::Level(..)))
    );
}

#[test]
fn pp1_only() {
    // p + 1 = 2^2 * 3 * 5813 * 14683 * 18691 * 35089 * 39227 * 300017, p - 1 has a 23-digit
    // prime factor.
    let p = int("7905527545387271831388442067");
    let n = Integer::from(10).pow(100) + 267u32;
    let n = Integer::from(&p * &n);
    let pp1 = |b1, b2| Factorizer::new().algorithm(Algorithm::Pp1).b1(b1).b2(b2);
    let (result, events) = record(pp1(40_000, 400_000), &n, |_| false);
    assert_eq!(result.error, None);
    assert!(events.contains(&Owned::Pp1(n.clone(), 40_000, Some(400_000))));
    assert!(events.contains(&Owned::Factor(n.clone(), p.clone(), Method::Pp1Stage2)));
    assert!(
        !events
            .iter()
            .any(|e| matches!(e, Owned::Curve(..) | Owned::Pm1(..)))
    );
    let (_, events) = record(pp1(400_000, 400_000), &n, |_| false);
    assert!(events.contains(&Owned::Factor(n.clone(), p.clone(), Method::Pp1Stage1)));
    assert_eq!(pp1(40_000, 200_000).find_factor(&n), Err(Error::ECMFailed));
    // P-1 does not find it.
    let pm1 = Factorizer::new()
        .algorithm(Algorithm::Pm1)
        .b1(400_000)
        .b2(4_000_000);
    assert_eq!(pm1.clone().find_factor(&n), Err(Error::ECMFailed));
    // The seeds: 6/5 and 3 work in the group of order p + 1 too, 4 in the one of order p - 1.
    let seed = |num: u32, den: u32| pp1(40_000, 400_000).x0(num.into(), den.into());
    assert_eq!(seed(6, 5).find_factor(&n), Ok(p.clone()));
    assert_eq!(seed(3, 1).find_factor(&n), Ok(p.clone()));
    assert_eq!(seed(4, 1).find_factor(&n), Err(Error::ECMFailed));
    // A seed not defined modulo p: its denominator is a factor.
    let q = Integer::from(1_000_003);
    let m = Integer::from(&q * &p);
    assert_eq!(
        pp1(40_000, 400_000).x0(1.into(), q.clone()).find_factor(&m),
        Ok(q.clone())
    );
    // Without fixed bounds: P+1 with the P-1 bounds of the levels only (B1 = 40000 at the
    // second one).
    let (result, events) = record(Factorizer::new().algorithm(Algorithm::Pp1), &n, |_| false);
    assert_eq!(result.error, None);
    assert!(
        !events
            .iter()
            .any(|e| matches!(e, Owned::Curve(..) | Owned::Level(..) | Owned::Pm1(..)))
    );
    assert!(events.contains(&Owned::Factor(n.clone(), p.clone(), Method::Pp1Stage2)));
    // The starting value of P-1 (x0 = 1 would find all the factors at once: none).
    let pm1 = |num: i32| pm1.clone().x0(num.into(), 1.into());
    assert!(matches!(
        pm1(1).find_factor(&n),
        Err(Error::InvalidOption(_))
    ));
    assert_eq!(pm1(5).find_factor(&n), Err(Error::ECMFailed));
}

#[test]
fn max_memory() {
    // Stage 2 with little memory: same factors.
    let n = int("4269021180054189416198169786894227");
    let factors = ecm(&n).unwrap();
    for bytes in [0, 1 << 20] {
        assert_eq!(
            Factorizer::new().max_memory(bytes).factor(&n).unwrap(),
            factors
        );
        assert_eq!(
            Factorizer::new()
                .max_memory(bytes)
                .b1(50_000)
                .b2(20_000_000)
                .curves(100)
                .factor(&n)
                .unwrap(),
            factors
        );
    }
}

#[test]
fn invalid_options() {
    let n = sample();
    for mut factorizer in [
        Factorizer::new().b2(1000),
        Factorizer::new().curves(10),
        Factorizer::new().sigma(Integer::from(7)),
        Factorizer::new()
            .b1(2000)
            .algorithm(Algorithm::Pm1)
            .sigma(Integer::from(7)),
        Factorizer::new()
            .b1(2000)
            .algorithm(Algorithm::Pp1)
            .sigma(Integer::from(7)),
        Factorizer::new().x0(2.into(), 7.into()),
        Factorizer::new().b1(2000).x0(2.into(), 7.into()),
        Factorizer::new()
            .algorithm(Algorithm::Pp1)
            .x0(2.into(), 0.into()),
        Factorizer::new()
            .algorithm(Algorithm::Pm1)
            .b1(2000)
            .x0(2.into(), 0.into()),
        Factorizer::new()
            .algorithm(Algorithm::Pp1)
            .x0((-4).into(), 2.into()),
        Factorizer::new()
            .algorithm(Algorithm::Pp1)
            .x0(0.into(), 5.into()),
        Factorizer::new()
            .algorithm(Algorithm::Pm1)
            .b1(2000)
            .x0(3.into(), (-3).into()),
        Factorizer::new().b1(2000).sigma(Integer::from(1)),
        Factorizer::new().b1(2000).sigma(Integer::from(1) << 64),
        Factorizer::new()
            .b1(2000)
            .param(Param::Square)
            .sigma(Integer::from(1) << 32),
        Factorizer::new()
            .b1(2000)
            .param(Param::Suyama)
            .sigma(Integer::from(5)),
    ] {
        assert!(matches!(
            factorizer.factor(&n),
            Err(Error::InvalidOption(_))
        ));
        assert!(matches!(
            factorizer.find_factor(&n),
            Err(Error::InvalidOption(_))
        ));
    }
    assert_eq!(
        Factorizer::new().b1(2001).factor(&n),
        Err(Error::BoundsNotEven)
    );
    assert_eq!(
        Factorizer::new().b1(4).factor(&n),
        Err(Error::BoundsTooSmall)
    );
    // The default b2 is valid.
    for b1 in [6, 100, 2000, 3000, 11_000, 1_000_000, 2_000_000_000] {
        let mut events = Vec::new();
        let _ = Factorizer::new()
            .b1(b1)
            .on_event(|e| {
                events.push(owned(e));
                ControlFlow::Break(())
            })
            .factor(&hard());
        assert!(
            matches!(events.last(), Some(Owned::TrialDivision(..))),
            "{b1}"
        );
    }
    for param in [Param::Suyama, Param::Square, Param::Batch2] {
        let value = u8::from(param);
        assert_eq!(Param::try_from(value), Ok(param));
        assert_eq!(param.to_string(), value.to_string());
    }
    assert!(Param::try_from(3).is_err());
}

#[test]
#[should_panic(expected = "only positive numbers")]
fn factor_zero() {
    let _ = Factorizer::new().factor(&Integer::ZERO);
}

#[test]
#[should_panic(expected = "greater than 1")]
fn find_factor_of_one() {
    let _ = Factorizer::new().find_factor(&Integer::from(1));
}

/// A 1024-bit semiprime, far too hard for ECM.
fn rsa1024() -> Integer {
    let p = (Integer::from(1) << 511u32) + 1_000u32;
    let q = (Integer::from(3) << 510u32) + 1_000u32;
    p.next_prime() * q.next_prime()
}

/// Factors `n` with `factorizer`, setting its interruption flag after `delay`: returns the
/// result, and the time the factorization took to return after the flag was set.
fn interrupt_after(
    n: &Integer,
    factorizer: Factorizer,
    delay: Duration,
) -> (Factorization, Duration) {
    let flag = Arc::new(AtomicBool::new(false));
    let setter = {
        let flag = Arc::clone(&flag);
        thread::spawn(move || {
            thread::sleep(delay);
            flag.store(true, Ordering::Relaxed);
            Instant::now()
        })
    };
    let result = factorizer.interrupt_flag(flag).factor_partial(n);
    let returned = Instant::now();
    let set = setter.join().unwrap();
    (result, returned.saturating_duration_since(set))
}

/// The interruption flag stops the factorization during a curve or P-1, a few milliseconds
/// after it is set (much less than the bound checked here), with the bounds of 35-digit factors
/// (larger bounds make the largest polynomial products of stage 2 take longer).
#[test]
fn interrupt_flag_is_prompt() {
    let n = rsa1024();
    let prompt = Duration::from_millis(300);
    let b2 = 1_045_563_762;
    for (factorizer, delay) in [
        // Stage 1 of the first curve (about 5 s).
        (Factorizer::new().b1(1_000_000), 300),
        // The polynomial stage 2 of the first curve (about 1.5 s), after a short stage 1.
        (Factorizer::new().b1(11_000).b2(b2), 400),
        (Factorizer::new().b1(11_000).b2(b2), 900),
        // Stage 1 of P-1, then stage 2.
        (
            Factorizer::new().algorithm(Algorithm::Pm1).b1(20_000_000),
            300,
        ),
        (
            Factorizer::new()
                .algorithm(Algorithm::Pm1)
                .b1(100_000)
                .b2(b2),
            300,
        ),
        // Stage 1 of P+1, then stage 2.
        (
            Factorizer::new().algorithm(Algorithm::Pp1).b1(20_000_000),
            300,
        ),
        (
            Factorizer::new()
                .algorithm(Algorithm::Pp1)
                .b1(100_000)
                .b2(b2),
            300,
        ),
        // The first level: P-1, then curves.
        (Factorizer::new(), 300),
    ] {
        let case = format!("{factorizer:?}, {delay} ms");
        let (result, latency) = interrupt_after(&n, factorizer, Duration::from_millis(delay));
        assert!(latency < prompt, "{case}: {latency:?}");
        assert_eq!(result.error, Some(Error::Interrupted));
        assert_eq!(result.unfactored, [(n.clone(), 1)]);
    }
}

#[test]
fn interrupt_flag_and_timeout() {
    // Set before: interrupted after trial division, with its factors.
    let n = hard() * 12u32;
    let flag = Arc::new(AtomicBool::new(true));
    let mut events = Vec::new();
    let result = Factorizer::new()
        .interrupt_flag(flag)
        .on_event(|e| {
            events.push(owned(e));
            ControlFlow::Continue(())
        })
        .factor_partial(&n);
    assert_eq!(result.error, Some(Error::Interrupted));
    check_partial(&n, &result);
    assert_eq!(result.primes[&Integer::from(2)], 2);
    assert_eq!(result.unfactored.len(), 1);
    assert!(
        !events
            .iter()
            .any(|e| matches!(e, Owned::Curve(..) | Owned::Pm1(..)))
    );

    // Timeout, for each call.
    let n = rsa1024() * 7u32;
    let mut factorizer = Factorizer::new().timeout(Duration::from_millis(200));
    for _ in 0..2 {
        let start = Instant::now();
        let result = factorizer.factor_partial(&n);
        let elapsed = start.elapsed();
        assert!(elapsed >= Duration::from_millis(200), "{elapsed:?}");
        assert!(elapsed < Duration::from_millis(500), "{elapsed:?}");
        assert_eq!(result.error, Some(Error::Interrupted));
        check_partial(&n, &result);
        assert_eq!(result.primes, HashMap::from([(Integer::from(7), 1)]));
    }
    assert_eq!(factorizer.find_factor(&rsa1024()), Err(Error::Interrupted));

    // Neither changes the results when it does not interrupt.
    let n = sample();
    let flag = Arc::new(AtomicBool::new(false));
    assert_eq!(
        Factorizer::new()
            .interrupt_flag(flag)
            .timeout(Duration::from_secs(3600))
            .factor(&n),
        ecm(&n)
    );
}

#[test]
fn interrupt_at_every_event() {
    // Wherever the callback interrupts, the primes and the unfactored parts multiply to n, and
    // they are the ones of the events so far.
    let n = sample() * Integer::from(4_009_823u32).pow(3);
    for factorizer in [Factorizer::new(), Factorizer::new().b1(2_000).curves(500)] {
        let (_, all) = record(factorizer.clone(), &n, |_| false);
        for stop in 1..=all.len() {
            let mut count = 0;
            let (result, events) = record(factorizer.clone(), &n, |_| {
                count += 1;
                count == stop
            });
            assert_eq!(events[..], all[..stop]);
            assert_eq!(result.error, Some(Error::Interrupted));
            check_partial(&n, &result);
            let mut primes = HashMap::new();
            for event in &events {
                if let Owned::Prime(p, e) = event {
                    *primes.entry(p.clone()).or_insert(0) += e;
                }
            }
            for (p, e) in &primes {
                assert!(result.primes[p] >= *e);
            }
            for (m, _) in &result.unfactored {
                assert!(*m > 1 && !result.primes.contains_key(m));
            }
        }
    }
}

#[test]
fn base2() {
    // 2^423 + 1 and 2^425 + 1: composite parts of 400 bits and less, those above about 320 bits
    // with the special reduction.
    for (n, k) in [
        (int("2").pow(423) + 1u32, 423),
        (int("2").pow(425) + 1u32, 425),
    ] {
        let (auto, events) = record(Factorizer::new(), &n, |_| false);
        let factors = auto.into_result().unwrap();
        check_events(&n, &factors, &events);
        let base2: Vec<_> = events
            .iter()
            .filter_map(|event| match event {
                Owned::Base2(m, k) => Some((m.significant_bits(), *k)),
                _ => None,
            })
            .collect();
        assert!(!base2.is_empty() && base2.iter().all(|&(bits, kk)| kk == k && bits > 300));

        let (off, events) = record(Factorizer::new().base2(Base2Mode::Off), &n, |_| false);
        assert_eq!(off.into_result().unwrap(), factors);
        assert!(!events.iter().any(|event| matches!(event, Owned::Base2(..))));

        // Forced, also modulo 2^(2k) - 1, a multiple of 2^k + 1, and for every part.
        for force in [k, -2 * k] {
            let (forced, events) =
                record(Factorizer::new().base2(Base2Mode::Force(force)), &n, |_| {
                    false
                });
            assert_eq!(forced.into_result().unwrap(), factors, "{force}");
            let parts = events
                .iter()
                .filter(|e| matches!(e, Owned::Level(..) | Owned::Pm1(..)))
                .count();
            assert!(
                parts > 0
                    && events
                        .iter()
                        .filter(|e| matches!(e, Owned::Base2(_, kk) if *kk == force))
                        .count()
                        > 0
            );
        }
    }
    // Forced modulo a number the composite parts do not divide.
    let n = int("2").pow(425) + 1u32;
    let result = Factorizer::new()
        .base2(Base2Mode::Force(-425))
        .factor_partial(&n);
    assert!(matches!(result.error, Some(Error::InvalidOption(_))));
    check_partial(&n, &result);
    // 3 * (2^67 - 1): trial division removes the 3, the rest divides 2^67 - 1.
    let n = int("3") * (int("2").pow(67) - 1u32);
    let expected: HashMap<Integer, usize> = [("3", 1), ("193707721", 1), ("761838257287", 1)]
        .map(|(p, e)| (int(p), e))
        .into();
    for mode in [
        Base2Mode::Force(-67),
        Base2Mode::Force(-134),
        Base2Mode::Auto,
    ] {
        assert_eq!(Factorizer::new().base2(mode).factor(&n).unwrap(), expected);
    }
    assert_eq!(
        Factorizer::new()
            .base2(Base2Mode::Force(-67))
            .find_factor(&(int("1000003") * int("1000033"))),
        Err(Error::InvalidOption(
            "base2: the number does not divide 2^k+1 (2^-k-1 if k < 0)"
        ))
    );
    // P-1 and P+1 with fixed bounds, forced.
    let n = int("2").pow(128) + 1u32;
    for algorithm in [Algorithm::Pm1, Algorithm::Pp1] {
        let run = |mode| {
            Factorizer::new()
                .algorithm(algorithm)
                .b1(100_000)
                .b2(10_000_000)
                .base2(mode)
                .find_factor(&n)
        };
        assert_eq!(
            run(Base2Mode::Force(128)),
            run(Base2Mode::Off),
            "{algorithm}"
        );
    }
}
