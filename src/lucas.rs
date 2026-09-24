//! Lucas sequences `V_k = a^k + a^-k` modulo `n`, the group of P-1's stage 2 and of P+1.
//!
//! `V_k` only depends on `V_1 = a + 1/a`, through `V_(j+k) = V_j*V_k - V_(j-k)` and
//! `V_2k = V_k^2 - 2`: differential additions, as for the x-coordinates of a Montgomery curve.
//! `V_jk(V_1) = V_j(V_k(V_1))`, so `V_E` for a product `E` is computed factor by factor. For
//! P-1, `a` is an integer modulo `n` (`V_1 = x + 1/x`). For P+1, `V_1 = x0` and `a` is a root of
//! `X^2 - x0*X + 1`: modulo a prime `p`, in the group of order `p - 1` of the integers when
//! `x0^2 - 4` is a square modulo `p`, else in the group of order `p + 1` of the norm-1 elements
//! of the field with `p^2` elements.
//!
//! As elements of [`XLine`], `(V_k : V_k - 2)`: the second coordinate vanishes modulo `p`
//! exactly when `a^k = 1` modulo `p` (`V_k = 2`).

use crate::{
    arith::{Arith, with_arith},
    base2::Base2Form,
    curve::{Scratch, Xz},
    primes::primes,
    stage2::{Normalizer, Stage2Plan, XLine, stage2_group},
    stop::{STOP_INTERVAL, Stop},
};
use rug::Integer;

/// `num/den` modulo `n`, or `Err(gcd(den, n))` if `den` is not invertible.
pub(crate) fn rational(num: &Integer, den: &Integer, n: &Integer) -> Result<Integer, Integer> {
    match den.invert_ref(n) {
        Some(inv) => Ok((Integer::from(inv) * num).modulo(n)),
        None => Err(den.gcd_ref(n).into()),
    }
}

/// Stage 2 from `V_1 = v` modulo `n` (with the special reduction modulo `base2`, a multiple of
/// `n`, if any): checks the primes `l` in `(b1, b2]` of `plan` (see [`crate::stage2`]), and
/// returns `gcd(g, n)`, partial if `stop` is requested.
pub(crate) fn stage2(
    n: &Integer,
    v: &Integer,
    plan: &Stage2Plan,
    base2: Option<Base2Form>,
    stop: Stop<'_>,
) -> Integer {
    with_arith!(n, base2, |arith| {
        let lucas = Lucas::new(arith);
        let start = lucas.element(v);
        stage2_group(&lucas, &start, plan, stop)
    })
}

/// `V_E(v)` modulo `n`, where `E` is the product of the largest powers of the primes `<= hi`
/// that are `<= hi`, divided by the same product for `lo` (see [`crate::ecm::stage1_multiplier`]),
/// with the special reduction modulo `base2` (a multiple of `n`) if any.
///
/// Returns `None` if `stop` was requested (the value is then lost: `v` is kept).
pub(crate) fn stage1(
    n: &Integer,
    v: &Integer,
    (lo, hi): (usize, usize),
    base2: Option<Base2Form>,
    stop: Stop<'_>,
) -> Option<Integer> {
    with_arith!(n, base2, |arith| {
        let lucas = Lucas::new(arith);
        let mut x = lucas.arith.residue(v);
        lucas
            .stage1(&mut x, lo, hi, stop)
            .then(|| lucas.arith.to_integer(&x))
    })
}

/// Lucas sequences modulo `n`, on the residues of `A`.
pub(crate) struct Lucas<A: Arith> {
    pub(crate) arith: A,
    two: A::Elem,
}

/// 1 over the golden ratio: [`Lucas::prac`] starts its chains with `r = k*VAL`.
const VAL: f64 = 0.618_033_988_749_894_9;

/// Temporaries of [`Lucas::prac`].
struct Prac<E> {
    b: E,
    c: E,
    t: E,
    u: E,
    prod: E,
    res: E,
}

impl<A: Arith> Lucas<A> {
    pub(crate) fn new(arith: A) -> Self {
        let two = arith.residue(&Integer::from(2));
        Self { arith, two }
    }

    /// The element `V`.
    pub(crate) fn element(&self, v: &Integer) -> Xz<A::Elem> {
        let x = self.arith.residue(v);
        let mut z = self.arith.zero();
        self.arith.sub(&mut z, &x, &self.two);
        Xz { x, z }
    }

    /// `x = V_k(x)` with the binary ladder (1M + 1S per bit), for `k >= 1`: returns `false`,
    /// with a meaningless `x`, if `stop` is requested.
    fn ladder(&self, x: &mut A::Elem, k: &Integer, stop: Stop<'_>) -> bool {
        let a = &self.arith;
        let (mut u, mut w, mut t) = (x.clone(), a.zero(), a.zero());
        // (u, w) = (V_j, V_(j + 1)), with j the bits of k above the current one.
        a.sqr(&mut t, &u);
        a.sub(&mut w, &t, &self.two);
        let mut end = k.significant_bits().saturating_sub(1);
        while end > 0 {
            let start = end.saturating_sub(STOP_INTERVAL);
            for bit in (start..end).rev() {
                a.mul(&mut t, &u, &w);
                if k.get_bit(bit) {
                    a.sub(&mut u, &t, x);
                    a.sqr(&mut t, &w);
                    a.sub(&mut w, &t, &self.two);
                } else {
                    a.sub(&mut w, &t, x);
                    a.sqr(&mut t, &u);
                    a.sub(&mut u, &t, &self.two);
                }
            }
            end = start;
            if end > 0 && stop.requested() {
                return false;
            }
        }
        *x = u;
        true
    }

    /// `r = p*q - s`: `V_(j+k)` from `V_j`, `V_k` and `V_(j-k)`.
    #[inline(always)]
    fn add3(&self, r: &mut A::Elem, p: &A::Elem, q: &A::Elem, s: &A::Elem, prod: &mut A::Elem) {
        self.arith.mul(prod, p, q);
        self.arith.sub(r, prod, s);
    }

    /// `r = p^2 - 2`: `V_2k` from `V_k`.
    #[inline(always)]
    fn dup(&self, r: &mut A::Elem, p: &A::Elem, prod: &mut A::Elem) {
        self.arith.sqr(prod, p);
        self.arith.sub(r, prod, &self.two);
    }

    /// `x = V_k(x)` for a prime `k >= 3` with Montgomery's PRAC Lucas chains (GMP-ECM's
    /// `pp1_mul_prac`): about 1.6 multiplications per bit of `k`, against 2 for the ladder.
    /// Some composite `k` (4, for one) never end.
    #[inline(always)]
    fn prac(&self, x: &mut A::Elem, k: u64, s: &mut Prac<A::Elem>) {
        debug_assert!((3..1 << 52).contains(&k));
        let Prac {
            b,
            c,
            t,
            u,
            prod,
            res,
        } = s;
        let mut a = std::mem::replace(x, self.arith.zero());
        // The chain keeps (A, B, C) = (V_i, V_j, V_(i-j)) with d*i + e*j = k (roughly).
        let r = (k as f64 * VAL + 0.5) as u64;
        let (mut d, mut e) = (k - r, 2 * r - k);
        b.clone_from(&a);
        c.clone_from(&a);
        self.dup(res, &a, prod);
        std::mem::swap(&mut a, res);
        while d != e {
            if d < e {
                std::mem::swap(&mut d, &mut e);
                std::mem::swap(&mut a, b);
            }
            if d - e <= e / 4 && (d + e).is_multiple_of(3) {
                d = (2 * d - e) / 3;
                e = (e - d) / 2;
                self.add3(t, &a, b, c, prod); // T = f(A, B, C)
                self.add3(u, t, &a, b, prod); // T2 = f(T, A, B)
                self.add3(res, b, t, &a, prod); // B = f(B, T, A)
                std::mem::swap(b, res);
                std::mem::swap(&mut a, u);
            } else if d - e <= e / 4 && (d - e).is_multiple_of(6) {
                d = (d - e) / 2;
                self.add3(res, &a, b, c, prod); // B = f(A, B, C)
                std::mem::swap(b, res);
                self.dup(res, &a, prod); // A = 2A
                std::mem::swap(&mut a, res);
            } else if d <= 4 * e {
                d -= e;
                self.add3(res, b, &a, c, prod); // C = f(B, A, C), then swap B and C
                std::mem::swap(c, b);
                std::mem::swap(b, res);
            } else if (d + e).is_multiple_of(2) {
                d = (d - e) / 2;
                self.add3(res, b, &a, c, prod); // B = f(B, A, C)
                std::mem::swap(b, res);
                self.dup(res, &a, prod); // A = 2A
                std::mem::swap(&mut a, res);
            } else if d.is_multiple_of(2) {
                // d + e is odd.
                d /= 2;
                self.add3(res, c, &a, b, prod); // C = f(C, A, B)
                std::mem::swap(c, res);
                self.dup(res, &a, prod); // A = 2A
                std::mem::swap(&mut a, res);
            } else if d.is_multiple_of(3) {
                // d is odd, e is even.
                d = d / 3 - e;
                self.dup(t, &a, prod); // T = 2A
                self.add3(u, &a, b, c, prod); // T2 = f(A, B, C)
                self.add3(res, t, &a, &a, prod); // A = f(T, A, A)
                std::mem::swap(&mut a, res);
                self.add3(res, t, u, c, prod); // C = f(T, T2, C), then swap B and C
                std::mem::swap(c, b);
                std::mem::swap(b, res);
            } else if (d + e).is_multiple_of(3) {
                d = (d - 2 * e) / 3;
                self.add3(t, &a, b, c, prod); // T = f(A, B, C)
                self.add3(res, t, &a, b, prod); // B = f(T, A, B)
                std::mem::swap(b, res);
                self.dup(t, &a, prod); // A = 3A
                self.add3(res, &a, t, &a, prod);
                std::mem::swap(&mut a, res);
            } else if (d - e).is_multiple_of(3) {
                d = (d - e) / 3;
                self.add3(t, &a, b, c, prod); // T = f(A, B, C)
                self.add3(res, c, &a, b, prod); // C = f(C, A, B)
                std::mem::swap(c, res);
                std::mem::swap(b, t); // swap B and T
                self.dup(t, &a, prod); // A = 3A
                self.add3(res, &a, t, &a, prod);
                std::mem::swap(&mut a, res);
            } else {
                // e is even.
                e /= 2;
                self.add3(res, c, b, &a, prod); // C = f(C, B, A)
                std::mem::swap(c, res);
                self.dup(res, b, prod); // B = 2B
                std::mem::swap(b, res);
            }
        }
        self.add3(x, &a, b, c, prod);
    }

    /// Stage 1 of P+1 from `lo` to `hi` (see [`stage1`]): the powers of the primes up to
    /// `sqrt(hi)` multiplied together, with the ladder, then each larger prime with
    /// [`Lucas::prac`], as GMP-ECM does. Returns `false` if `stop` was requested.
    fn stage1(&self, x: &mut A::Elem, lo: usize, hi: usize, stop: Stop<'_>) -> bool {
        #[cfg(target_arch = "x86_64")]
        if crate::arith::has_bmi2_adx() {
            // SAFETY: the CPU has the features `stage1_bmi2` is compiled for.
            return unsafe { self.stage1_bmi2(x, lo, hi, stop) };
        }
        self.stage1_generic(x, lo, hi, stop)
    }

    /// [`Lucas::stage1`] compiled with BMI2 and ADX (see [`crate::arith::has_bmi2_adx`]).
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "bmi2,adx")]
    fn stage1_bmi2(&self, x: &mut A::Elem, lo: usize, hi: usize, stop: Stop<'_>) -> bool {
        self.stage1_generic(x, lo, hi, stop)
    }

    #[inline(always)]
    fn stage1_generic(&self, x: &mut A::Elem, lo: usize, hi: usize, stop: Stop<'_>) -> bool {
        let lo = lo.max(1);
        if lo >= hi {
            return true;
        }
        // 2 is not a valid input of `prac`.
        let root = hi.isqrt().max(2);
        let small = crate::ecm::prime_power_product_below(lo, hi, root);
        if !self.ladder(x, &small, stop) {
            return false;
        }
        let mut s = Prac {
            b: self.arith.zero(),
            c: self.arith.zero(),
            t: self.arith.zero(),
            u: self.arith.zero(),
            prod: self.arith.zero(),
            res: self.arith.zero(),
        };
        // About STOP_INTERVAL chain steps between two checks of `stop`.
        let interval = (STOP_INTERVAL as usize / hi.ilog2().max(1) as usize).max(1);
        for (i, p) in primes(hi).skip_while(|&p| p <= root.max(lo)).enumerate() {
            self.prac(x, p as u64, &mut s);
            if i % interval == interval - 1 && stop.requested() {
                return false;
            }
        }
        true
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
        self.dup(&mut r.x, &p.x, &mut scratch.0[0]);
        self.arith.sub(&mut r.z, &r.x, &self.two);
    }

    fn add(
        &self,
        r: &mut Xz<A::Elem>,
        p: &Xz<A::Elem>,
        q: &Xz<A::Elem>,
        diff: &Xz<A::Elem>,
        scratch: &mut Scratch<A::Elem>,
    ) {
        self.add3(&mut r.x, &p.x, &q.x, &diff.x, &mut scratch.0[0]);
        self.arith.sub(&mut r.z, &r.x, &self.two);
    }

    fn multiple(&self, p: &Xz<A::Elem>, k: &Integer) -> Xz<A::Elem> {
        let mut x = p.x.clone();
        self.ladder(&mut x, k, Stop::NEVER);
        let mut z = self.arith.zero();
        self.arith.sub(&mut z, &x, &self.two);
        Xz { x, z }
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

    fn n() -> Integer {
        Integer::from(1_000_000_007u64) * Integer::from(998_244_353u64)
    }

    /// `V_k(v)` modulo `n` from the definition, with `v = x + 1/x`.
    fn reference(n: &Integer) -> (Integer, impl Fn(u64) -> Integer) {
        let x = Integer::from(123_456_789);
        let v1 = (Integer::from(x.invert_ref(n).unwrap()) + &x) % n;
        let n = n.clone();
        (v1, move |k: u64| {
            let a = x.clone().pow_mod(&Integer::from(k), &n).unwrap();
            (Integer::from(a.invert_ref(&n).unwrap()) + a) % &n
        })
    }

    #[test]
    fn multiples() {
        fn check<A: Arith>(lucas: &Lucas<A>, n: &Integer) {
            let (v1, v) = reference(n);
            let a = &lucas.arith;
            let p = lucas.element(&v1);
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
        let n = n();
        check(&Lucas::new(Mont::<2>::new(&n)), &n);
        check(&Lucas::new(Plain::new(&n)), &n);
    }

    #[test]
    fn prac_matches_ladder() {
        fn check<A: Arith>(lucas: &Lucas<A>, n: &Integer) {
            let (v1, v) = reference(n);
            let a = &lucas.arith;
            let x = a.residue(&v1);
            let mut s = Prac {
                b: a.zero(),
                c: a.zero(),
                t: a.zero(),
                u: a.zero(),
                prod: a.zero(),
                res: a.zero(),
            };
            let large = [
                999_983,
                4_294_967_291,
                1_000_000_000_039,
                4_503_599_627_370_449,
            ];
            let small = primal::Primes::all().skip(1).take_while(|&p| p < 30_000);
            for k in small.map(|p| p as u64).chain(large) {
                let mut y = x.clone();
                lucas.prac(&mut y, k, &mut s);
                assert_eq!(a.to_integer(&y), v(k), "{k}");
            }
        }
        let n = n();
        check(&Lucas::new(Mont::<2>::new(&n)), &n);
        check(&Lucas::new(Plain::new(&n)), &n);
    }

    #[test]
    fn stage1_is_the_multiplier() {
        let n = n();
        let (v1, _) = reference(&n);
        for (lo, hi) in [
            (1, 2),
            (1, 3),
            (1, 10),
            (1, 1000),
            (10, 1000),
            (500, 1000),
            (1000, 1000),
        ] {
            let mut expected = Lucas::new(Plain::new(&n)).element(&v1);
            let k = crate::ecm::prime_power_product(lo, hi);
            expected = Lucas::new(Plain::new(&n)).multiple(&expected, &k);
            let got = stage1(&n, &v1, (lo, hi), None, Stop::NEVER).unwrap();
            assert_eq!(got, expected.x, "{lo} {hi}");
        }
        // Resumed.
        let direct = stage1(&n, &v1, (1, 100_000), None, Stop::NEVER).unwrap();
        let mut v = v1;
        for (lo, hi) in [
            (1, 10),
            (10, 11),
            (11, 2000),
            (2000, 99_999),
            (99_999, 100_000),
        ] {
            v = stage1(&n, &v, (lo, hi), None, Stop::NEVER).unwrap();
        }
        assert_eq!(v, direct);
    }

    #[test]
    #[cfg(target_arch = "x86_64")]
    fn stage1_generic_matches_dispatched() {
        let (v1, _) = reference(&n());
        for bits in [60, 128, 300, 1000] {
            let n = Integer::from(Integer::u_pow_u(2, bits)) - 1u32;
            let run = || stage1(&n, &v1, (1, 20_000), None, Stop::NEVER).unwrap();
            crate::arith::GENERIC_ONLY.set(true);
            let generic = run();
            crate::arith::GENERIC_ONLY.set(false);
            assert_eq!(generic, run(), "{bits} bits");
        }
    }
}
