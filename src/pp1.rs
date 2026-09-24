//! Williams' P+1 method, the essentials of GMP-ECM's `pp1.c`.
//!
//! Stage 1 computes `V_E(x0)` modulo `n` (see [`crate::lucas`]), with `E` the stage 1 multiplier
//! of ECM and P-1: modulo a prime `p` of `n`, it is `2` when `E` is a multiple of the order of
//! the root `a` of `X^2 - x0*X + 1`, which divides `p + 1` when `x0^2 - 4` is not a square
//! modulo `p`, else `p - 1`. Half of the seeds `x0` thus find the factors `p` with a smooth
//! `p + 1`, the other half behave as P-1. The seed `2/7` (the default, as GMP-ECM recommends)
//! makes the order a multiple of 6: `x0^2 - 4 = -3*(8/7)^2`, so P+1 works for the primes
//! `p = 2 mod 3` (the order is `p + 1`), and P-1 for the others; `6/5` makes it a multiple of 4.
//! Stage 1 can be resumed with a larger `b1`, and `V` stays valid modulo any divisor of `n`.
//!
//! Stage 2 checks the primes `l` in `(b1, b2]` with the same continuations as ECM and P-1 (see
//! [`crate::stage2`]), in the same Lucas sequence.

use crate::{base2::Base2Form, lucas, stage2::Stage2Plan, stop::Stop};
use rug::Integer;

/// Default seed `x0 = 2/7` (numerator, denominator).
const DEFAULT_X0: (u32, u32) = (2, 7);

/// State of P+1 on a number: stage 1 done up to `b1` from the seed `x0`.
#[derive(Debug, Clone)]
pub struct Pp1 {
    /// Numerator and denominator of the seed.
    x0: (Integer, Integer),
    /// `V_E(b1)(x0)` modulo the number (or a multiple of it), `None` before stage 1.
    v: Option<Integer>,
    b1: usize,
}

impl Default for Pp1 {
    fn default() -> Self {
        Self::new(DEFAULT_X0.0.into(), DEFAULT_X0.1.into())
    }
}

impl Pp1 {
    /// Nothing done yet, with the seed `x0 = numerator/denominator` (`denominator != 0`).
    #[must_use]
    pub fn new(numerator: Integer, denominator: Integer) -> Self {
        Self {
            x0: (numerator, denominator),
            v: None,
            b1: 1,
        }
    }

    /// Stage 1 bound reached so far.
    #[must_use]
    pub fn b1(&self) -> usize {
        self.b1
    }

    /// `V_E(x0)` after stage 1 (GMP-ECM's residue `x` at the end of stage 1, printed by
    /// `ecm -pp1 -v -v -v`), modulo the last number.
    #[cfg_attr(not(any(test, feature = "bench")), allow(dead_code))]
    #[must_use]
    pub fn residue(&self) -> Option<&Integer> {
        self.v.as_ref()
    }

    /// Extends stage 1 to `b1` modulo `n` (a divisor of the previous numbers, if any): returns
    /// `gcd(V - 2, n)` (or `gcd(denominator, n)` if the seed is not defined modulo `n`).
    #[cfg_attr(not(any(test, feature = "bench")), allow(dead_code))]
    pub fn stage1(&mut self, n: &Integer, b1: usize) -> Integer {
        self.stage1_until(n, b1, Base2Form::detect(n), Stop::NEVER)
    }

    /// [`Pp1::stage1`] with the special reduction modulo `base2` (a multiple of `n`) if any,
    /// stopping early if `stop` is requested: the bound reached is then not updated, and `V` is
    /// kept.
    pub(crate) fn stage1_until(
        &mut self,
        n: &Integer,
        b1: usize,
        base2: Option<Base2Form>,
        stop: Stop<'_>,
    ) -> Integer {
        let mut v = match self.v.take() {
            Some(v) => v % n,
            None => match lucas::rational(&self.x0.0, &self.x0.1, n) {
                Ok(v) => v,
                Err(g) => return g,
            },
        };
        if self.b1 < b1 {
            if let Some(w) = lucas::stage1(n, &v, (self.b1, b1), base2, stop) {
                v = w;
                self.b1 = b1;
            }
        }
        let g = Integer::from(&v - 2u32).gcd(n);
        self.v = Some(v);
        g
    }

    /// Stage 2 modulo `n` with `plan`, whose `b1` must be at most the stage 1 bound: returns
    /// `gcd(g, n)` (see [`crate::ecm::stage2`]).
    ///
    /// # Panics
    ///
    /// If stage 1 did not run.
    #[cfg_attr(not(any(test, feature = "bench")), allow(dead_code))]
    #[must_use]
    pub fn stage2(&self, n: &Integer, plan: &Stage2Plan) -> Integer {
        self.stage2_until(n, plan, Base2Form::detect(n), Stop::NEVER)
    }

    /// [`Pp1::stage2`] with the special reduction modulo `base2` (a multiple of `n`) if any,
    /// stopping early (with the gcd of a partial product) if `stop` is requested.
    pub(crate) fn stage2_until(
        &self,
        n: &Integer,
        plan: &Stage2Plan,
        base2: Option<Base2Form>,
        stop: Stop<'_>,
    ) -> Integer {
        let v = self.v.as_ref().expect("stage 1 runs first");
        lucas::stage2(n, &Integer::from(v % n), plan, base2, stop)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        arith::{Arith, Plain},
        lucas::Lucas,
        stage2::XLine,
    };
    use primal::Primes;
    use rug::{ops::Pow, rand::RandState};

    /// Order of `a` (root of `X^2 - x0*X + 1`) modulo the prime `p`: a divisor of `p - 1` if
    /// `x0^2 - 4` is a square modulo `p`, else of `p + 1`. `None` if `x0^2 = 4` modulo `p`.
    fn group_order(x0: &Integer, p: u64) -> Option<u64> {
        let p_ = Integer::from(p);
        let d = Integer::from(x0 * x0) - 4u32;
        match d.jacobi(&p_) {
            0 => None,
            1 => Some(p - 1),
            _ => Some(p + 1),
        }
    }

    /// `x0 = num/den` modulo `m`.
    fn seed(m: &Integer) -> Integer {
        let (num, den) = DEFAULT_X0;
        (Integer::from(den).invert(m).unwrap() * num).modulo(m)
    }

    /// `V_k(x0)` modulo the odd `m`.
    fn v(x0: &Integer, k: u64, m: &Integer) -> Integer {
        let lucas = Lucas::new(Plain::new(m));
        let r = lucas.multiple(&lucas.element(x0), &Integer::from(k));
        lucas.arith.to_integer(&r.x)
    }

    /// Smallest `k` dividing `order` with `V_k(x0) = 2` modulo `p`.
    fn element_order(x0: &Integer, p: u64, order: u64) -> u64 {
        let p_ = Integer::from(p);
        let (mut k, mut rest) = (order, order);
        let mut q = 2;
        while rest > 1 {
            if q * q > rest {
                q = rest;
            }
            if rest.is_multiple_of(q) {
                while rest.is_multiple_of(q) {
                    rest /= q;
                }
                while k.is_multiple_of(q) && v(x0, k / q, &p_) == 2 {
                    k /= q;
                }
            }
            q += 1;
        }
        assert_eq!(v(x0, k, &p_), 2);
        k
    }

    #[test]
    fn stage1_finds_smooth_p_plus_1() {
        // p + 1 = 2^5 * 3 * 5 * 7 * 11 * 13 * 17 * 23 (p = 2 mod 3: found by the seed 2/7);
        // q - 1 and q + 1 have a large prime factor.
        let p = Integer::from(187_867_679u64);
        let q = Integer::from(10).pow(20) + 39u32;
        let n = Integer::from(&p * &q);
        assert_eq!(Pp1::default().stage1(&n, 22), 1);
        assert_eq!(Pp1::default().stage1(&n, 32), p);
        let mut resumed = Pp1::default();
        assert_eq!(resumed.stage1(&n, 10), 1);
        assert_eq!(resumed.stage1(&n, 32), p);
        // A seed not defined modulo n.
        assert_eq!(Pp1::new(1.into(), p.clone()).stage1(&n, 10), p);
    }

    #[test]
    fn stage1_resumes() {
        let n = Integer::from(1_000_000_007u64) * Integer::from(998_244_353u64);
        let mut direct = Pp1::default();
        direct.stage1(&n, 100_000);
        let mut resumed = Pp1::default();
        for b1 in [10, 11, 100, 2_000, 99_999, 100_000] {
            resumed.stage1(&n, b1);
        }
        assert_eq!(resumed.b1(), 100_000);
        assert_eq!(resumed.residue(), direct.residue());
        // Modulo a divisor.
        let p = Integer::from(1_000_000_007u64);
        let mut divisor = direct.clone();
        divisor.stage1(&p, 200_000);
        let mut fresh = Pp1::default();
        fresh.stage1(&p, 200_000);
        assert_eq!(divisor.v, fresh.v);
    }

    /// Stage 2 must find `p` whenever the order of `a` modulo `p` after stage 1 is a prime `l`
    /// in `(b1, b2]`, for each small prime `p` times the large prime `q`.
    ///
    /// Returns the number of primes `p` where stage 2 had a factor to find, with an order
    /// dividing `p + 1` and `p - 1`.
    fn check_stage2_primes(q: &Integer, b1: usize, b2: usize, primes: &[u64]) -> [usize; 2] {
        let e = crate::ecm::stage1_multiplier(b1);
        let mut checked = [0; 2];
        for &p in primes {
            let n = Integer::from(q * p);
            let plans = [Stage2Plan::pairs(b1, b2), Stage2Plan::poly(&n, b1, b2)];
            let mut pp1 = Pp1::default();
            if pp1.stage1(&n, b1) != 1 {
                continue;
            }
            let x0 = seed(&Integer::from(p));
            let order = group_order(&x0, p).unwrap();
            // Order of V_E: the order of a divided by its gcd with E.
            let k = element_order(&x0, p, order);
            let k = k / Integer::from(k).gcd(&e).to_u64().unwrap();
            let expected = k > b1 as u64
                && k <= b2 as u64
                && Integer::from(k).is_probably_prime(20) != rug::integer::IsPrime::No;
            for plan in &plans {
                let g = pp1.stage2(&n, plan);
                assert!(n.is_divisible(&g));
                if expected {
                    assert_eq!(g, p, "stage 2 missed p = {p}: {plan:?} {b1} {b2}");
                }
            }
            if expected {
                checked[usize::from(order == p - 1)] += 1;
            }
        }
        checked
    }

    #[test]
    fn stage2_checks_all_primes() {
        let q = Integer::from(10).pow(30) + 57u32;
        let primes: Vec<u64> = Primes::all()
            .skip_while(|&p| p < 100_000)
            .take(800)
            .map(|p| p as u64)
            .collect();
        let mut checked = [0; 2];
        for (b1, b2) in [(30, 1000), (100, 5000), (100, 10_001), (200, 20_000)] {
            let [plus, minus] = check_stage2_primes(&q, b1, b2, &primes);
            checked = [checked[0] + plus, checked[1] + minus];
        }
        assert!(checked[0] > 50 && checked[1] > 50, "{checked:?}");
    }

    #[test]
    fn stage2_checks_all_primes_limbs() {
        // Every limb count of `Mont`, and `Plain`.
        let mut rand = RandState::new();
        let primes: Vec<u64> = Primes::all()
            .skip_while(|&p| p < 10_000)
            .take(100)
            .map(|p| p as u64)
            .collect();
        for bits in [64, 100, 192, 320, 512, 700, 1000, 1100] {
            let mut q = Integer::from(Integer::random_bits(bits, &mut rand));
            q.set_bit(bits - 1, true);
            let [plus, minus] = check_stage2_primes(&q.next_prime(), 50, 3000, &primes);
            assert!(plus + minus > 5, "{bits} bits: {plus} {minus}");
        }
    }

    #[test]
    fn base2_matches_generic() {
        for k in [-101, 128, -263, 512, -1061] {
            let (n, form) = crate::base2::cofactor_of(k);
            let (mut on, mut off) = (Pp1::default(), Pp1::default());
            for b1 in [1_000, 30_000] {
                let g = on.stage1_until(&n, b1, Some(form), Stop::NEVER);
                assert_eq!(g, off.stage1_until(&n, b1, None, Stop::NEVER), "{form}");
                assert_eq!(on.v, off.v, "{form}");
            }
            for plan in [
                Stage2Plan::pairs(30_000, 3_000_000),
                Stage2Plan::poly(&n, 30_000, 3_000_000),
            ] {
                let g = on.stage2_until(&n, &plan, Some(form), Stop::NEVER);
                assert_eq!(g, off.stage2_until(&n, &plan, None, Stop::NEVER), "{form}");
            }
        }
    }
}
