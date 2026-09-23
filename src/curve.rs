//! Fast `(X : Z)` arithmetic on Montgomery curves `b*y^2 = x^3 + a*x^2 + x`, on residues of any
//! [`Arith`] implementation.
//!
//! This is the hot path of both stages: [`Point`](crate::point::Point) is only the interface
//! type, converted to and from this representation once per stage.

use crate::arith::{Arith, Factor};
use rug::Integer;

/// Projective point `(x : z)` (the `y` coordinate is never needed).
#[derive(Debug, Clone)]
pub struct Xz<E> {
    pub x: E,
    pub z: E,
}

/// Temporaries of the point operations, allocated once per ladder (or stage).
pub struct Scratch<E>(pub [E; 4]);

/// A Montgomery curve, given by `a24 = (a + 2)/4`, over the residues of `A`.
pub struct Curve<A: Arith> {
    pub arith: A,
    a24: Factor<A::Elem>,
}

impl<A: Arith> Curve<A> {
    /// Curve with parameter `a24 = (a + 2)/4`.
    pub fn new(arith: A, a24: &Integer) -> Self {
        let a24 = arith.factor(a24);
        Curve { arith, a24 }
    }

    /// Point with the given coordinates.
    pub fn point(&self, x: &Integer, z: &Integer) -> Xz<A::Elem> {
        Xz {
            x: self.arith.residue(x),
            z: self.arith.residue(z),
        }
    }

    /// Point at infinity, used as a placeholder.
    pub fn infinity(&self) -> Xz<A::Elem> {
        Xz {
            x: self.arith.zero(),
            z: self.arith.zero(),
        }
    }

    /// Temporaries for the point operations.
    pub fn scratch(&self) -> Scratch<A::Elem> {
        Scratch(std::array::from_fn(|_| self.arith.zero()))
    }

    /// `r = 2*p`: 2M + 2S, plus a multiplication by `a24` (cheap if it is small).
    pub fn double(&self, r: &mut Xz<A::Elem>, p: &Xz<A::Elem>, scratch: &mut Scratch<A::Elem>) {
        let a = &self.arith;
        let [s, d, t, w] = &mut scratch.0;
        a.add(s, &p.x, &p.z);
        a.sqr(t, s); // (x + z)^2
        a.sub(d, &p.x, &p.z);
        a.sqr(w, d); // (x - z)^2
        a.mul(&mut r.x, t, w);
        a.sub(s, t, w); // 4xz
        a.mul_factor(d, s, &self.a24);
        a.add(t, w, d);
        a.mul(&mut r.z, s, t);
    }

    /// Computes `(u + v)^2` in `scratch[2]` and `(u - v)^2` in `scratch[3]`, with
    /// `u = (xp - zp)(xq + zq)` and `v = (xp + zp)(xq - zq)`: 2M + 2S.
    ///
    /// Then `p + q = (zd * (u + v)^2 : xd * (u - v)^2)`, where `(xd : zd) = p - q`.
    #[inline]
    fn add_squares(&self, p: &Xz<A::Elem>, q: &Xz<A::Elem>, scratch: &mut Scratch<A::Elem>) {
        let a = &self.arith;
        let [s, d, t, w] = &mut scratch.0;
        a.sub(s, &p.x, &p.z);
        a.add(d, &q.x, &q.z);
        a.mul(t, s, d); // u
        a.add(s, &p.x, &p.z);
        a.sub(d, &q.x, &q.z);
        a.mul(w, s, d); // v
        a.add(s, t, w);
        a.sub(d, t, w);
        a.sqr(t, s);
        a.sqr(w, d);
    }

    /// `r = p + q`, where `diff = p - q`: 4M + 2S.
    pub fn add(
        &self,
        r: &mut Xz<A::Elem>,
        p: &Xz<A::Elem>,
        q: &Xz<A::Elem>,
        diff: &Xz<A::Elem>,
        scratch: &mut Scratch<A::Elem>,
    ) {
        self.add_squares(p, q, scratch);
        let [_, _, t, w] = &scratch.0;
        self.arith.mul(&mut r.x, t, &diff.z);
        self.arith.mul(&mut r.z, w, &diff.x);
    }

    /// One Montgomery ladder step: `(p, q) = (2*p, p + q)`, where `q - p = (xd : zd)`.
    ///
    /// 4M + 4S, plus the multiplications by `xd`, `zd` and `a24`: with GMP-ECM's parametrization
    /// 1 (`(xd : zd) = (2 : 1)` and a small `a24`), these are two additions and a one-limb
    /// multiplication.
    #[inline(always)]
    pub fn dup_add(
        &self,
        p: &mut Xz<A::Elem>,
        q: &mut Xz<A::Elem>,
        xd: &Factor<A::Elem>,
        zd: &Factor<A::Elem>,
        scratch: &mut Scratch<A::Elem>,
    ) {
        let a = &self.arith;
        let [s, d, t, w] = &mut scratch.0;
        a.add(s, &p.x, &p.z);
        a.sub(d, &p.x, &p.z);
        a.add(t, &q.x, &q.z);
        a.sub(w, &q.x, &q.z);
        a.mul(&mut q.x, d, t); // u = (xp - zp)(xq + zq)
        a.mul(&mut q.z, s, w); // v = (xp + zp)(xq - zq)
        a.sqr(t, s); // (xp + zp)^2
        a.sqr(w, d); // (xp - zp)^2
        a.mul(&mut p.x, t, w);
        a.sub(s, t, w);
        a.mul_factor(d, s, &self.a24);
        a.add(t, w, d);
        a.mul(&mut p.z, s, t);
        a.add(s, &q.x, &q.z);
        a.sub(d, &q.x, &q.z);
        a.sqr(t, s);
        a.sqr(w, d);
        a.mul_factor(&mut q.x, t, zd);
        a.mul_factor(&mut q.z, w, xd);
    }

    /// `k*P` with the Montgomery ladder, where `P = (xp : zp)`, for `k >= 1`.
    pub fn ladder(&self, xp: &Factor<A::Elem>, zp: &Factor<A::Elem>, k: &Integer) -> Xz<A::Elem> {
        let a = &self.arith;
        let mut scratch = self.scratch();
        // (p, q) = (j*P, (j + 1)*P), with j the bits of k above the current one.
        let mut p = Xz {
            x: a.factor_elem(xp),
            z: a.factor_elem(zp),
        };
        let mut q = self.infinity();
        self.double(&mut q, &p, &mut scratch);
        for bit in (0..k.significant_bits().saturating_sub(1)).rev() {
            if k.get_bit(bit) {
                self.dup_add(&mut q, &mut p, xp, zp, &mut scratch);
            } else {
                self.dup_add(&mut p, &mut q, xp, zp, &mut scratch);
            }
        }
        p
    }
}
