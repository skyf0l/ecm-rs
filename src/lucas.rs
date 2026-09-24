//! Lucas sequences `V_k = a^k + a^-k` modulo `n`, the group of P-1's stage 2.
//!
//! `V_k` only depends on `V_1 = a + 1/a`, through `V_(j+k) = V_j*V_k - V_(j-k)` and
//! `V_2k = V_k^2 - 2`: differential additions, as for the x-coordinates of a Montgomery curve.
//! For P-1, `a` is an integer modulo `n` (`V_1 = x + 1/x`).
//!
//! As elements of [`XLine`], `(V_k : V_k - 2)`: the second coordinate vanishes modulo `p`
//! exactly when `a^k = 1` modulo `p` (`V_k = 2`).

use crate::{
    arith::{Arith, with_arith},
    curve::{Scratch, Xz},
    stage2::{Normalizer, Stage2Plan, XLine, stage2_group},
    stop::{STOP_INTERVAL, Stop},
};
use rug::Integer;

/// Stage 2 from `V_1 = v` modulo `n`: checks the primes `l` in `(b1, b2]` of `plan` (see
/// [`crate::stage2`]), and returns `gcd(g, n)`, partial if `stop` is requested.
pub(crate) fn stage2(n: &Integer, v: &Integer, plan: &Stage2Plan, stop: Stop<'_>) -> Integer {
    with_arith!(n, |arith| {
        let lucas = Lucas::new(arith);
        let start = lucas.element(v);
        stage2_group(&lucas, &start, plan, stop)
    })
}

/// Lucas sequences modulo `n`, on the residues of `A`.
pub(crate) struct Lucas<A: Arith> {
    pub(crate) arith: A,
    two: A::Elem,
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
}
