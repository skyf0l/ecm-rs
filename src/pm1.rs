//! Pollard's P-1 method, the essentials of GMP-ECM's `pm1.c`.
//!
//! Stage 1 computes `x = x0^E mod n` with `E` the product of the largest powers of the primes
//! `<= b1` that are `<= b1` (the same multiplier as ECM's stage 1): a prime `p` of `n` with
//! `b1`-smooth `p - 1` divides `x - 1`. Stage 1 can be resumed with a larger `b1`, and `x` stays
//! valid modulo any divisor of `n`.
//!
//! Stage 2 checks the primes `l` in `(b1, b2]` with the same continuations as ECM (see
//! [`crate::stage2`]), on the Lucas sequence `V_k = x^k + x^-k`: `V_(k+l) = V_k*V_l - V_(k-l)`
//! and `V_2k = V_k^2 - 2` are differential additions, and
//! `x^(m*D) * (V_(m*D) - V_j) = (x^(m*D - j) - 1) * (x^(m*D + j) - 1)`, so `V_(m*D) = V_j`
//! modulo `p` checks `m*D +- j` as `x(m*D*Q) = x(j*Q)` does for a curve. The identity element
//! is `V_k = 2`, that is `x^k = 1`.

use crate::{
    arith::{with_arith, Arith, PolyArith},
    curve::{Scratch, Xz},
    ecm::{prime_power_words, product},
    stage2::{stage2_group, Normalizer, Stage2Plan, XLine},
};
use rug::Integer;

/// Starting value `x0` of stage 1 (the choice barely matters: `E` has a large power of 2).
const X0: u32 = 3;

/// 64-bit factors of the exponent multiplied together before each modular exponentiation (about
/// 2^20 bits): bounds the memory used.
const CHUNK_WORDS: usize = 1 << 14;

/// State of P-1 on a number: stage 1 done up to `b1`.
#[derive(Debug, Clone)]
pub struct Pm1 {
    /// `x0^E(b1)` modulo the number (or a multiple of it).
    x: Integer,
    b1: usize,
}

impl Default for Pm1 {
    fn default() -> Self {
        Pm1 {
            x: Integer::from(X0),
            b1: 1,
        }
    }
}

impl Pm1 {
    /// Nothing done yet.
    pub fn new() -> Self {
        Self::default()
    }

    /// Stage 1 bound reached so far.
    pub fn b1(&self) -> usize {
        self.b1
    }

    /// Extends stage 1 to `b1` modulo `n` (a divisor of the previous numbers, if any): returns
    /// `gcd(x - 1, n)`.
    pub fn stage1(&mut self, n: &Integer, b1: usize) -> Integer {
        self.x %= n;
        if self.b1 < b1 {
            // The exponent by chunks of about CHUNK_BITS bits.
            let mut words = prime_power_words(self.b1, b1).peekable();
            while words.peek().is_some() {
                let chunk: Vec<Integer> = words
                    .by_ref()
                    .take(CHUNK_WORDS)
                    .map(Integer::from)
                    .collect();
                self.x
                    .pow_mod_mut(&product(chunk), n)
                    .expect("positive exponent");
            }
            self.b1 = b1;
        }
        Integer::from(&self.x - 1u32).gcd(n)
    }

    /// Stage 2 modulo `n` with `plan`, whose `b1` must be at most the stage 1 bound: returns
    /// `gcd(g, n)` (see [`crate::ecm::stage2`]).
    pub fn stage2(&self, n: &Integer, plan: &Stage2Plan) -> Integer {
        let x = Integer::from(&self.x % n);
        with_arith!(n, |arith| stage2_with(arith, &x, plan))
    }
}

fn stage2_with<A: PolyArith>(arith: A, x: &Integer, plan: &Stage2Plan) -> Integer {
    let n = arith.modulus();
    let Ok(inv) = x.clone().invert(n) else {
        return x.clone().gcd(n);
    };
    let v1 = (inv + x) % n;
    let lucas = Lucas::new(arith);
    let start = lucas.element(&v1);
    stage2_group(&lucas, &start, plan)
}

/// Lucas sequences `V_k = x^k + x^-k` modulo `n`, as elements `(V_k : V_k - 2)`: the second
/// coordinate vanishes modulo `p` exactly when `x^k = 1` modulo `p`.
struct Lucas<A: Arith> {
    arith: A,
    two: A::Elem,
}

impl<A: Arith> Lucas<A> {
    fn new(arith: A) -> Self {
        let two = arith.residue(&Integer::from(2));
        Lucas { arith, two }
    }

    /// The element `V`.
    fn element(&self, v: &Integer) -> Xz<A::Elem> {
        let x = self.arith.residue(v);
        let mut z = self.arith.zero();
        self.arith.sub(&mut z, &x, &self.two);
        Xz { x, z }
    }
}

impl<A: Arith> XLine for Lucas<A> {
    type A = A;

    fn arith(&self) -> &A {
        &self.arith
    }

    fn infinity(&self) -> Xz<A::Elem> {
        Xz {
            x: self.two.clone(),
            z: self.arith.zero(),
        }
    }

    fn scratch(&self) -> Scratch<A::Elem> {
        Scratch(std::array::from_fn(|_| self.arith.zero()))
    }

    fn double(&self, r: &mut Xz<A::Elem>, p: &Xz<A::Elem>, scratch: &mut Scratch<A::Elem>) {
        let a = &self.arith;
        let t = &mut scratch.0[0];
        a.sqr(t, &p.x);
        a.sub(&mut r.x, t, &self.two);
        a.sub(&mut r.z, &r.x, &self.two);
    }

    fn add(
        &self,
        r: &mut Xz<A::Elem>,
        p: &Xz<A::Elem>,
        q: &Xz<A::Elem>,
        diff: &Xz<A::Elem>,
        scratch: &mut Scratch<A::Elem>,
    ) {
        let a = &self.arith;
        let t = &mut scratch.0[0];
        a.mul(t, &p.x, &q.x);
        a.sub(&mut r.x, t, &diff.x);
        a.sub(&mut r.z, &r.x, &self.two);
    }

    fn multiple(&self, p: &Xz<A::Elem>, k: &Integer) -> Xz<A::Elem> {
        let a = &self.arith;
        // (u, w) = (V_j, V_(j + 1)), with j the bits of k above the current one.
        let mut u = p.x.clone();
        let mut w = a.zero();
        let mut t = a.zero();
        a.sqr(&mut t, &u);
        a.sub(&mut w, &t, &self.two);
        for bit in (0..k.significant_bits().saturating_sub(1)).rev() {
            a.mul(&mut t, &u, &w);
            if k.get_bit(bit) {
                a.sub(&mut u, &t, &p.x);
                a.sqr(&mut t, &w);
                a.sub(&mut w, &t, &self.two);
            } else {
                a.sub(&mut w, &t, &p.x);
                a.sqr(&mut t, &u);
                a.sub(&mut u, &t, &self.two);
            }
        }
        let mut z = a.zero();
        a.sub(&mut z, &u, &self.two);
        Xz { x: u, z }
    }

    /// The `x` are already normalized: only checks that no `z` shares a factor with `n`, with a
    /// product and a single gcd.
    fn normalize(
        &self,
        _normalizer: &mut Normalizer<A::Elem>,
        _x: &mut [A::Elem],
        z: &[A::Elem],
    ) -> Result<(), Integer> {
        let a = &self.arith;
        let (mut acc, mut t) = (a.residue(&Integer::from(1)), a.zero());
        for z in z {
            a.mul(&mut t, &acc, z);
            std::mem::swap(&mut acc, &mut t);
        }
        let g = a.gcd(&acc);
        let n = a.modulus();
        if g == 1 {
            Ok(())
        } else if &g != n {
            Err(g)
        } else {
            Err(z
                .iter()
                .map(|z| a.gcd(z))
                .find(|g| *g != 1 && g != n)
                .unwrap_or(g))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arith::{Mont, Plain};
    use primal::Primes;
    use rug::{ops::Pow, rand::RandState};

    #[test]
    fn lucas_multiples() {
        let n = Integer::from(1_000_000_007u64) * Integer::from(998_244_353u64);
        let x = Integer::from(123_456_789);
        let v = |k: u32| {
            let a = x.clone().pow_mod(&Integer::from(k), &n).unwrap();
            (Integer::from(a.invert_ref(&n).unwrap()) + a) % &n
        };
        fn check<A: Arith>(lucas: Lucas<A>, n: &Integer, v: &dyn Fn(u32) -> Integer) {
            let a = &lucas.arith;
            let p = lucas.element(&v(1));
            let mut scratch = lucas.scratch();
            let (mut r, mut s) = (lucas.infinity(), lucas.infinity());
            lucas.double(&mut r, &p, &mut scratch);
            assert_eq!(a.to_integer(&r.x), v(2));
            lucas.add(&mut s, &r, &p, &p, &mut scratch);
            assert_eq!(a.to_integer(&s.x), v(3));
            assert_eq!(a.to_integer(&s.z), (v(3) + n - 2u32) % n);
            for k in [1, 2, 3, 4, 5, 6, 7, 100, 1000, 65_537, 1_234_567] {
                let m = lucas.multiple(&p, &Integer::from(k));
                assert_eq!(a.to_integer(&m.x), v(k), "{k}");
            }
        }
        check(Lucas::new(Mont::<2>::new(&n)), &n, &v);
        check(Lucas::new(Plain::new(&n)), &n, &v);
    }

    #[test]
    fn stage1_resumes() {
        let n = Integer::from(1_000_000_007u64) * Integer::from(998_244_353u64);
        let mut direct = Pm1::new();
        direct.stage1(&n, 100_000);
        let mut resumed = Pm1::new();
        for b1 in [10, 11, 100, 2_000, 99_999, 100_000] {
            resumed.stage1(&n, b1);
        }
        assert_eq!(resumed.b1(), 100_000);
        assert_eq!(resumed.x, direct.x);
        // Larger than a chunk of the exponent.
        let b1 = 2_000_000;
        let mut pm1 = Pm1::new();
        pm1.stage1(&n, b1);
        let e = crate::ecm::stage1_multiplier(b1);
        assert_eq!(pm1.x, Integer::from(X0).pow_mod(&e, &n).unwrap());
    }

    #[test]
    fn stage1_finds_smooth() {
        // 1000000007 - 1 = 2 * 500000003, 998244353 - 1 = 2^23 * 7 * 17.
        let n = Integer::from(1_000_000_007u64) * Integer::from(998_244_353u64);
        assert_eq!(Pm1::new().stage1(&n, 10), 1);
        assert_eq!(Pm1::new().stage1(&n, 1 << 23), 998_244_353);
    }

    /// Stage 2 must find `p` whenever `x^l = 1` modulo `p` for some prime `b1 < l <= b2` after
    /// stage 1 (`x != 1`), for each small prime `p` times the large prime `q`.
    ///
    /// Returns the number of primes `p` where stage 2 had a factor to find.
    fn check_stage2_primes(q: &Integer, b1: usize, b2: usize, primes: &[u64]) -> usize {
        let ls: Vec<u64> = Primes::all()
            .skip_while(|&l| l <= b1)
            .take_while(|&l| l <= b2)
            .map(|l| l as u64)
            .collect();
        let mut checked = 0;
        for &p in primes {
            let n = Integer::from(q * p);
            let plans = [Stage2Plan::pairs(b1, b2), Stage2Plan::poly(&n, b1, b2)];
            let mut pm1 = Pm1::new();
            if pm1.stage1(&n, b1) != 1 {
                continue;
            }
            let y = Integer::from(&pm1.x % p);
            let expected = ls.iter().any(|&l| {
                y.clone()
                    .pow_mod(&Integer::from(l), &Integer::from(p))
                    .unwrap()
                    == 1
            });
            for plan in &plans {
                let g = pm1.stage2(&n, plan);
                assert!(n.is_divisible(&g));
                if expected {
                    assert_eq!(g, p, "stage 2 missed p = {p}: {plan:?} {b1} {b2}");
                }
            }
            checked += expected as usize;
        }
        checked
    }

    #[test]
    fn stage2_checks_all_primes() {
        let q = Integer::from(10).pow(30) + 57u32;
        assert_ne!(q.is_probably_prime(30), rug::integer::IsPrime::No);
        let primes: Vec<u64> = Primes::all()
            .skip_while(|&p| p < 100_000)
            .take(800)
            .map(|p| p as u64)
            .collect();
        let mut checked = 0;
        for (b1, b2) in [
            (30, 1000),
            (100, 5000),
            (100, 10_001),
            (200, 20_000),
            (1000, 999),
        ] {
            checked += check_stage2_primes(&q, b1, b2, &primes);
        }
        assert!(checked > 100, "{checked}");
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
            let checked = check_stage2_primes(&q.next_prime(), 50, 3000, &primes);
            assert!(checked > 5, "{bits} bits: {checked}");
        }
    }
}
