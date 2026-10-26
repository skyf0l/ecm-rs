use super::*;
use rug::ops::{Pow, RemRounding};

/// SIQS on one thread until a split: the parts found.
pub(crate) fn factor(n: &Integer, seed: u64) -> Vec<Integer> {
    let mut siqs = Siqs::new(n, seed).unwrap_or_else(|e| panic!("{n}: {e:?}"));
    loop {
        while !siqs.enough() {
            let a = siqs.next_a();
            let batch = sieve_a(siqs.shared(), &a, Stop::NEVER).unwrap();
            siqs.merge(batch);
        }
        match siqs.finish(Stop::NEVER).unwrap().0 {
            Some(parts) => return parts,
            None => siqs.more(),
        }
    }
}

/// The product of the `k`-th primes after `10^(digits-1)` and after `2 * 10^(digits-1)`... a
/// semiprime with two factors of `digits` digits.
fn semiprime(digits: u32, k: u32) -> (Integer, Integer, Integer) {
    let base = Integer::from(10).pow(digits - 1);
    let mut p = Integer::from(&base * (3 + k)) + 12345u32 * k;
    p.next_prime_mut();
    let mut q = Integer::from(&base * (7 + 2 * k)) + 54321u32 * k;
    q.next_prime_mut();
    let n = Integer::from(&p * &q);
    (n, p, q)
}

#[test]
fn relations_are_valid() {
    let (n, _, _) = semiprime(20, 1);
    let mut siqs = Siqs::new(&n, 1).unwrap();
    let sh = siqs.shared().clone();
    let mut fulls = 0;
    let mut partials = 0;
    for _ in 0..4 {
        let a = siqs.next_a();
        let batch = sieve_a(&sh, &a, Stop::NEVER).unwrap();
        assert_eq!(batch.polynomials, 1 << (sh.s - 1));
        for r in batch.fulls.iter().chain(&batch.partials) {
            let mut value = Integer::from(r.large);
            for &f in &r.factors {
                if f == 0 {
                    value = -value;
                } else {
                    value *= sh.fb.primes[f as usize - 1];
                }
            }
            let y2 = Integer::from(r.y.square_ref()) % &n;
            let value = value.rem_euc(&n);
            assert_eq!(y2, value, "{r:?}");
            assert!(r.factors.windows(2).all(|w| w[0] <= w[1]));
        }
        fulls += batch.fulls.len();
        partials += batch.partials.len();
    }
    assert!(fulls > 10 && partials > 10, "{fulls} {partials}");
}

#[test]
fn a_values() {
    let (n, _, _) = semiprime(30, 2);
    let mut siqs = Siqs::new(&n, 5).unwrap();
    let sh = siqs.shared().clone();
    let mut seen = HashSet::new();
    for i in 0..200 {
        let a = siqs.next_a();
        assert_eq!(a.index, i);
        assert_eq!(a.q.len(), sh.s);
        assert!(seen.insert(a.q.clone()));
        let bits: f64 = a.q.iter().map(|&i| f64::from(sh.fb.primes[i]).log2()).sum();
        assert!((bits - sh.a_bits).abs() < 1.0, "{bits} vs {}", sh.a_bits);
    }
    // The same sequence from the same seed.
    let mut again = Siqs::new(&n, 5).unwrap();
    let mut first = Siqs::new(&n, 5).unwrap();
    for _ in 0..20 {
        assert_eq!(again.next_a().q, first.next_a().q);
    }
}

#[test]
fn semiprimes() {
    for (digits, k) in [(10, 0), (12, 1), (15, 2), (18, 3), (20, 4), (25, 5)] {
        let (n, p, q) = semiprime(digits, k);
        let parts = factor(&n, u64::from(k));
        let mut expected = vec![p, q];
        expected.sort_unstable();
        assert_eq!(parts, expected, "{n}");
    }
}

#[test]
fn three_factors() {
    let p = Integer::from(10u64.pow(14)).next_prime();
    let q = Integer::from(3 * 10u64.pow(14)).next_prime();
    let r = Integer::from(7 * 10u64.pow(14)).next_prime();
    let n = Integer::from(&p * &q) * &r;
    let parts = factor(&n, 0);
    assert_eq!(parts.iter().product::<Integer>(), n);
    assert!(parts.len() >= 2);
    for part in &parts {
        assert!(n.is_divisible(part) && *part > 1 && *part < n);
    }
}

#[test]
fn guards() {
    let (n, p, _) = semiprime(15, 1);
    // Perfect powers.
    let square = Integer::from(p.square_ref());
    assert_eq!(Siqs::new(&square, 0).err(), Some(Setup::Factor(p.clone())));
    let fourth = Integer::from(square.square_ref());
    assert_eq!(
        Siqs::new(&fourth, 0).err(),
        Some(Setup::Factor(square.clone()))
    );
    let cube = Integer::from(&n * &n) * &n;
    assert_eq!(Siqs::new(&cube, 0).err(), Some(Setup::Factor(n.clone())));
    // Even, prime, too small, too large.
    let even = Integer::from(&n * 2u32);
    assert_eq!(
        Siqs::new(&even, 0).err(),
        Some(Setup::Factor(Integer::from(2)))
    );
    assert_eq!(Siqs::new(&p, 0).err(), Some(Setup::Unsupported));
    assert_eq!(
        Siqs::new(&Integer::from(1_000_003u64 * 1_000_033), 0).err(),
        Some(Setup::Unsupported)
    );
    assert_eq!(
        Siqs::new(&Integer::from(10).pow(120), 0).err(),
        Some(Setup::Unsupported)
    );
    // A small factor, in the factor base (or dividing the multiplier).
    let big = Integer::from(10).pow(60).next_prime();
    for small in [3u32, 5, 7, 101, 7919, 65537] {
        let m = Integer::from(&big * small);
        assert_eq!(
            Siqs::new(&m, 0).err(),
            Some(Setup::Factor(Integer::from(small))),
            "{small}"
        );
    }
}

#[test]
fn digit_count() {
    for d in 1..60 {
        assert_eq!(digits(&Integer::from(10).pow(d - 1)), d);
        assert_eq!(digits(&(Integer::from(10).pow(d) - 1u32)), d);
    }
}

#[test]
fn params() {
    let mut last = 0;
    for d in 20..=110 {
        let p = Params::for_digits(f64::from(d));
        assert!(p.fb >= last && p.blocks >= 1 && p.lp_mult > 1);
        last = p.fb;
    }
}
