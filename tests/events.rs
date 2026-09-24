//! The public `Factorizer` API: events, cancellation, fixed bounds, determinism.

use ecm::{Error, Event, Factorization, Factorizer, Method, Param, ecm, ecm_with_params};
use rug::{Integer, ops::Pow};
use std::{
    collections::HashMap,
    ops::ControlFlow,
    str::FromStr,
    time::{Duration, Instant},
};

fn int(s: &str) -> Integer {
    Integer::from_str(s).unwrap()
}

/// An owned copy of the data of an event.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Owned {
    TrialDivision(HashMap<Integer, usize>, Integer),
    Pm1(Integer, usize, Option<usize>),
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
fn record<F>(factorizer: Factorizer, n: &Integer, mut stop: F) -> (Factorization, Vec<Owned>)
where
    F: FnMut(&Owned) -> bool,
{
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
            Owned::Pm1(m, b1, b2) => {
                assert!(n.is_divisible(m));
                assert!(b2.is_none_or(|b2| b2 > *b1));
                level = None;
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
    let start = Instant::now();
    let (result, events) = record(factorizer.clone(), n, stop);
    assert!(
        start.elapsed() < Duration::from_secs(2),
        "{:?}",
        start.elapsed()
    );
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
    interrupt(&n, Factorizer::new().b1(11_000).pm1(true), |e| {
        matches!(e, Owned::Pm1(..))
    });

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
    let pm1 = |b1, b2| Factorizer::new().pm1(true).b1(b1).b2(b2);
    let (result, events) = record(pm1(40_000, 400_000), &n, |_| false);
    assert_eq!(result.error, None);
    assert!(events.contains(&Owned::Factor(n.clone(), p.clone(), Method::Pm1Stage2)));
    assert!(!events.iter().any(|e| matches!(e, Owned::Curve(..))));
    assert_eq!(pm1(400_000, 400_000).find_factor(&n), Ok(p.clone()));
    assert_eq!(pm1(40_000, 200_000).find_factor(&n), Err(Error::ECMFailed));
    // Without fixed bounds: P-1 of the levels only.
    let (result, events) = record(Factorizer::new().pm1(true), &n, |_| false);
    assert_eq!(result.error, None);
    assert!(
        !events
            .iter()
            .any(|e| matches!(e, Owned::Curve(..) | Owned::Level(..)))
    );
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
        Factorizer::new().b1(2000).pm1(true).sigma(Integer::from(7)),
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
