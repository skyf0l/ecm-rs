//! Modular arithmetic for the curve computations.
//!
//! Residues modulo an odd `n` of at most [`MAX_LIMBS`] 64-bit limbs are kept in Montgomery
//! representation (`x` is stored as `x*R mod n` with `R = 2^(64*limbs)`), in fixed-size limb
//! arrays: a modular multiplication is one interleaved multiply-and-reduce (CIOS) pass, with no
//! division and no allocation. Larger (or even) moduli fall back to [`Plain`] big integer
//! arithmetic.
//!
//! The representation never leaks: values go in with [`Arith::residue`] and out with
//! [`Arith::to_integer`], and [`Arith::gcd`] gives `gcd(x, n)` directly (`R` is coprime to an
//! odd `n`).

use rug::{integer::Order, Assign, Integer};

/// Largest number of 64-bit limbs handled by [`Mont`]; larger moduli use [`Plain`].
pub const MAX_LIMBS: usize = 16;

/// Arithmetic modulo a fixed `n`, on residues of type [`Arith::Elem`].
///
/// Every operation writes its result to `r`, which never aliases an operand.
pub trait Arith {
    /// Residue modulo `n`, in the internal representation.
    type Elem: Clone;

    /// The modulus `n`.
    fn modulus(&self) -> &Integer;
    /// The residue `0`.
    fn zero(&self) -> Self::Elem;
    /// Residue of `x` (any integer).
    fn residue(&self, x: &Integer) -> Self::Elem;
    /// Value in `[0, n)` of the residue `x`.
    fn to_integer(&self, x: &Self::Elem) -> Integer;
    /// `gcd(x, n)`.
    fn gcd(&self, x: &Self::Elem) -> Integer;
    /// `r = a * b`.
    fn mul(&self, r: &mut Self::Elem, a: &Self::Elem, b: &Self::Elem);
    /// `r = a^2`.
    fn sqr(&self, r: &mut Self::Elem, a: &Self::Elem);
    /// `r = a + b`.
    fn add(&self, r: &mut Self::Elem, a: &Self::Elem, b: &Self::Elem);
    /// `r = a - b`.
    fn sub(&self, r: &mut Self::Elem, a: &Self::Elem, b: &Self::Elem);
    /// One-limb form of the residue `x`, for [`Arith::mul_small`], if it has one.
    fn small(&self, x: &Integer) -> Option<u64>;
    /// `r = a * x`, where `c = small(x)`: much cheaper than a full multiplication.
    fn mul_small(&self, r: &mut Self::Elem, a: &Self::Elem, c: u64);

    /// Multiplier `x`, in the cheapest form to multiply by.
    fn factor(&self, x: &Integer) -> Factor<Self::Elem> {
        let x = reduce(x, self.modulus());
        if x == 1 {
            Factor::One
        } else if x == 2 {
            Factor::Two
        } else if let Some(c) = self.small(&x) {
            Factor::Small(c)
        } else {
            Factor::Full(self.residue(&x))
        }
    }

    /// `r = a * f`.
    #[inline]
    fn mul_factor(&self, r: &mut Self::Elem, a: &Self::Elem, f: &Factor<Self::Elem>) {
        match f {
            Factor::One => r.clone_from(a),
            Factor::Two => self.add(r, a, a),
            Factor::Small(c) => self.mul_small(r, a, *c),
            Factor::Full(b) => self.mul(r, a, b),
        }
    }

    /// Residue of the multiplier `f`.
    fn factor_elem(&self, f: &Factor<Self::Elem>) -> Self::Elem {
        let mut r = self.zero();
        let one = self.residue(&Integer::from(1));
        self.mul_factor(&mut r, &one, f);
        r
    }
}

/// A constant multiplier, stored in the form that is cheapest to multiply by.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Factor<E> {
    /// `1`: a copy.
    One,
    /// `2`: an addition.
    Two,
    /// A residue with a one-limb form (see [`Arith::small`]).
    Small(u64),
    /// Any other residue: a full multiplication.
    Full(E),
}

/// Calls `$body` with `$a` bound to the fastest [`Arith`] implementation for the modulus `$n`.
///
/// `$body` is compiled once per implementation, so it should be a call to a generic function.
macro_rules! with_arith {
    ($n:expr, |$a:ident| $body:expr) => {{
        let n: &rug::Integer = $n;
        match $crate::arith::mont_limbs(n) {
            1 => {
                let $a = $crate::arith::Mont::<1>::new(n);
                $body
            }
            2 => {
                let $a = $crate::arith::Mont::<2>::new(n);
                $body
            }
            3 => {
                let $a = $crate::arith::Mont::<3>::new(n);
                $body
            }
            4 => {
                let $a = $crate::arith::Mont::<4>::new(n);
                $body
            }
            5 => {
                let $a = $crate::arith::Mont::<5>::new(n);
                $body
            }
            6 => {
                let $a = $crate::arith::Mont::<6>::new(n);
                $body
            }
            7 => {
                let $a = $crate::arith::Mont::<7>::new(n);
                $body
            }
            8 => {
                let $a = $crate::arith::Mont::<8>::new(n);
                $body
            }
            9 => {
                let $a = $crate::arith::Mont::<9>::new(n);
                $body
            }
            10 => {
                let $a = $crate::arith::Mont::<10>::new(n);
                $body
            }
            11 => {
                let $a = $crate::arith::Mont::<11>::new(n);
                $body
            }
            12 => {
                let $a = $crate::arith::Mont::<12>::new(n);
                $body
            }
            13 => {
                let $a = $crate::arith::Mont::<13>::new(n);
                $body
            }
            14 => {
                let $a = $crate::arith::Mont::<14>::new(n);
                $body
            }
            15 => {
                let $a = $crate::arith::Mont::<15>::new(n);
                $body
            }
            16 => {
                let $a = $crate::arith::Mont::<16>::new(n);
                $body
            }
            _ => {
                let $a = $crate::arith::Plain::new(n);
                $body
            }
        }
    }};
}
pub(crate) use with_arith;

/// Number of limbs of the [`Mont`] implementation for `n`, or `0` if `n` needs [`Plain`].
pub fn mont_limbs(n: &Integer) -> usize {
    let limbs = n.significant_bits().div_ceil(64) as usize;
    if n.is_odd() && *n > 1 && limbs <= MAX_LIMBS {
        limbs
    } else {
        0
    }
}

/// `x mod n`, in `[0, n)`.
fn reduce(x: &Integer, n: &Integer) -> Integer {
    let mut r = Integer::from(x % n);
    if r < 0 {
        r += n;
    }
    r
}

/// `lo + a * b + carry`, as `(low limb, high limb)`: never overflows.
#[inline(always)]
fn mac(lo: u64, a: u64, b: u64, carry: u64) -> (u64, u64) {
    let t = lo as u128 + a as u128 * b as u128 + carry as u128;
    (t as u64, (t >> 64) as u64)
}

/// `a + b + carry` (carry is 0 or 1), as `(sum, carry out)`.
#[inline(always)]
fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let (s, c1) = a.overflowing_add(b);
    let (s, c2) = s.overflowing_add(carry);
    (s, (c1 | c2) as u64)
}

/// `a - b - borrow` (borrow is 0 or 1), as `(difference, borrow out)`.
#[inline(always)]
fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let (d, b1) = a.overflowing_sub(b);
    let (d, b2) = d.overflowing_sub(borrow);
    (d, (b1 | b2) as u64)
}

/// Montgomery arithmetic modulo an odd `n` of at most `N` limbs.
#[derive(Debug, Clone)]
pub struct Mont<const N: usize> {
    /// Limbs of `n`, least significant first.
    m: [u64; N],
    /// `-1/n mod 2^64`.
    ninv: u64,
    n: Integer,
}

impl<const N: usize> Mont<N> {
    /// Arithmetic modulo `n`, which must be odd and fit in `N` limbs.
    pub fn new(n: &Integer) -> Self {
        assert!(n.is_odd() && n.significant_bits() as usize <= 64 * N);
        let m = Self::limbs(n);
        // Newton iteration: each step doubles the number of correct low bits of 1/m[0].
        let mut inv: u64 = 1;
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(m[0].wrapping_mul(inv)));
        }
        Mont {
            m,
            ninv: inv.wrapping_neg(),
            n: n.clone(),
        }
    }

    /// Limbs of `0 <= x < 2^(64*N)`.
    fn limbs(x: &Integer) -> [u64; N] {
        let mut limbs = [0; N];
        x.write_digits(&mut limbs, Order::Lsf);
        limbs
    }

    /// `t - n` if `t >= n` (with `t = t + carry*2^(64*N) < 2n`), else `t`.
    #[inline(always)]
    fn reduce_once(&self, r: &mut [u64; N], t: &[u64; N], carry: u64) {
        let mut d = [0; N];
        let mut borrow = 0;
        for j in 0..N {
            (d[j], borrow) = sbb(t[j], self.m[j], borrow);
        }
        // Keep `t` only if `t < n`: no carry, and the subtraction borrowed.
        let keep = (borrow & !carry).wrapping_neg();
        for j in 0..N {
            r[j] = (t[j] & keep) | (d[j] & !keep);
        }
    }
}

impl<const N: usize> Arith for Mont<N> {
    type Elem = [u64; N];

    fn modulus(&self) -> &Integer {
        &self.n
    }

    fn zero(&self) -> [u64; N] {
        [0; N]
    }

    fn residue(&self, x: &Integer) -> [u64; N] {
        let x = reduce(x, &self.n) << (64 * N as u32);
        Self::limbs(&(x % &self.n))
    }

    fn to_integer(&self, x: &[u64; N]) -> Integer {
        let mut one = [0; N];
        one[0] = 1;
        let mut r = [0; N];
        self.mul(&mut r, x, &one);
        Integer::from_digits(&r, Order::Lsf)
    }

    fn gcd(&self, x: &[u64; N]) -> Integer {
        Integer::from_digits(x, Order::Lsf).gcd(&self.n)
    }

    /// Montgomery multiplication, CIOS method: `r = a*b/R mod n`.
    #[inline]
    fn mul(&self, r: &mut [u64; N], a: &[u64; N], b: &[u64; N]) {
        let m = &self.m;
        // t = (t_hi, t_n, t[N-1..0]) < 2n at the end of every iteration.
        let mut t = [0u64; N];
        let mut t_n = 0u64;
        for &bi in b.iter() {
            let mut c = 0;
            for j in 0..N {
                (t[j], c) = mac(t[j], a[j], bi, c);
            }
            let (s, t_hi) = adc(t_n, c, 0);
            // Add q*n, with q such that the low limb becomes 0, and shift down by one limb.
            let q = t[0].wrapping_mul(self.ninv);
            let (_, mut c) = mac(t[0], q, m[0], 0);
            for j in 1..N {
                (t[j - 1], c) = mac(t[j], q, m[j], c);
            }
            let carry;
            (t[N - 1], carry) = adc(s, c, 0);
            t_n = t_hi + carry;
        }
        self.reduce_once(r, &t, t_n);
    }

    /// For mid sizes, separate squaring (half the products of a multiplication) and reduction
    /// are faster than CIOS.
    #[inline]
    fn sqr(&self, r: &mut [u64; N], a: &[u64; N]) {
        if !(3..=12).contains(&N) {
            return self.mul(r, a, a);
        }
        let m = &self.m;
        let mut wide = [[0u64; N]; 2];
        let t = wide.as_flattened_mut();
        // Products a[i]*a[j] for i < j, then doubled.
        for i in 0..N {
            let mut c = 0;
            for j in i + 1..N {
                (t[i + j], c) = mac(t[i + j], a[j], a[i], c);
            }
            t[i + N] = c;
        }
        let mut c = 0;
        for x in t.iter_mut() {
            (*x, c) = ((*x << 1) | c, *x >> 63);
        }
        // Squares a[i]^2.
        let mut carry = 0;
        for i in 0..N {
            let (lo, hi) = mac(0, a[i], a[i], 0);
            (t[2 * i], carry) = adc(t[2 * i], lo, carry);
            (t[2 * i + 1], carry) = adc(t[2 * i + 1], hi, carry);
        }
        // Montgomery reduction, one limb at a time.
        let mut hi = 0;
        for i in 0..N {
            let q = t[i].wrapping_mul(self.ninv);
            let mut c = 0;
            for j in 0..N {
                (t[i + j], c) = mac(t[i + j], q, m[j], c);
            }
            (t[i + N], c) = adc(t[i + N], c, 0);
            let carry;
            (t[i + N], carry) = adc(t[i + N], hi, 0);
            hi = c | carry;
        }
        self.reduce_once(r, &wide[1], hi);
    }

    #[inline]
    fn add(&self, r: &mut [u64; N], a: &[u64; N], b: &[u64; N]) {
        let mut s = [0; N];
        let mut carry = 0;
        for j in 0..N {
            (s[j], carry) = adc(a[j], b[j], carry);
        }
        self.reduce_once(r, &s, carry);
    }

    #[inline]
    fn sub(&self, r: &mut [u64; N], a: &[u64; N], b: &[u64; N]) {
        let mut borrow = 0;
        for j in 0..N {
            (r[j], borrow) = sbb(a[j], b[j], borrow);
        }
        // Add n back if the subtraction borrowed.
        let mask = borrow.wrapping_neg();
        let mut carry = 0;
        for (r, &m) in r.iter_mut().zip(&self.m) {
            (*r, carry) = adc(*r, m & mask, carry);
        }
    }

    /// The one-limb form of `x` is `c = x*2^64 mod n`, if `c < 2^64`: then `a*x*R` is the
    /// Montgomery reduction by one limb `a*R*c/2^64 mod n`.
    fn small(&self, x: &Integer) -> Option<u64> {
        (Integer::from(x << 64) % &self.n).to_u64()
    }

    #[inline]
    fn mul_small(&self, r: &mut [u64; N], a: &[u64; N], c: u64) {
        let m = &self.m;
        // t = a*c < n*2^64
        let mut t = [0u64; N];
        let mut carry = 0;
        for j in 0..N {
            (t[j], carry) = mac(0, a[j], c, carry);
        }
        // (t + q*n) / 2^64 < 2n
        let q = t[0].wrapping_mul(self.ninv);
        let (_, mut c) = mac(t[0], q, m[0], 0);
        let mut s = [0u64; N];
        for j in 1..N {
            (s[j - 1], c) = mac(t[j], q, m[j], c);
        }
        let hi;
        (s[N - 1], hi) = adc(carry, c, 0);
        self.reduce_once(r, &s, hi);
    }
}

/// Plain big integer arithmetic modulo any `n`: residues are integers in `[0, n)`.
///
/// Fallback for the moduli [`Mont`] doesn't handle.
#[derive(Debug, Clone)]
pub struct Plain {
    n: Integer,
}

impl Plain {
    /// Arithmetic modulo `n > 0`.
    pub fn new(n: &Integer) -> Self {
        Plain { n: n.clone() }
    }
}

impl Arith for Plain {
    type Elem = Integer;

    fn modulus(&self) -> &Integer {
        &self.n
    }

    fn zero(&self) -> Integer {
        // Large enough for a product: computing in place never reallocates.
        Integer::with_capacity(2 * self.n.significant_bits() as usize + 64)
    }

    fn residue(&self, x: &Integer) -> Integer {
        let mut r = self.zero();
        r.assign(reduce(x, &self.n));
        r
    }

    fn to_integer(&self, x: &Integer) -> Integer {
        x.clone()
    }

    fn gcd(&self, x: &Integer) -> Integer {
        x.clone().gcd(&self.n)
    }

    fn mul(&self, r: &mut Integer, a: &Integer, b: &Integer) {
        r.assign(a * b);
        *r %= &self.n;
    }

    fn sqr(&self, r: &mut Integer, a: &Integer) {
        r.assign(a.square_ref());
        *r %= &self.n;
    }

    fn add(&self, r: &mut Integer, a: &Integer, b: &Integer) {
        r.assign(a + b);
        if *r >= self.n {
            *r -= &self.n;
        }
    }

    fn sub(&self, r: &mut Integer, a: &Integer, b: &Integer) {
        r.assign(a - b);
        if *r < 0 {
            *r += &self.n;
        }
    }

    fn small(&self, x: &Integer) -> Option<u64> {
        reduce(x, &self.n).to_u64()
    }

    fn mul_small(&self, r: &mut Integer, a: &Integer, c: u64) {
        r.assign(a * c);
        *r %= &self.n;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rug::rand::RandState;

    /// Checks every operation of `A` against plain integer arithmetic, on random residues.
    fn check<A: Arith>(a: &A, rand: &mut RandState<'_>) {
        let n = a.modulus().clone();
        let mut values: Vec<Integer> = (0..20).map(|_| n.clone().random_below(rand)).collect();
        values.extend([Integer::ZERO, Integer::from(1), Integer::from(&n - 1)]);
        let mut r = a.zero();
        for x in &values {
            let ex = a.residue(x);
            assert_eq!(a.to_integer(&ex), *x);
            assert_eq!(a.gcd(&ex), x.clone().gcd(&n));
            a.sqr(&mut r, &ex);
            assert_eq!(a.to_integer(&r), Integer::from(x * x) % &n);
            for y in &values {
                let ey = a.residue(y);
                a.mul(&mut r, &ex, &ey);
                assert_eq!(a.to_integer(&r), Integer::from(x * y) % &n);
                a.add(&mut r, &ex, &ey);
                assert_eq!(a.to_integer(&r), Integer::from(x + y) % &n);
                a.sub(&mut r, &ex, &ey);
                assert_eq!(a.to_integer(&r), reduce(&Integer::from(x - y), &n));
                let f = a.factor(y);
                a.mul_factor(&mut r, &ex, &f);
                assert_eq!(a.to_integer(&r), Integer::from(x * y) % &n);
                assert_eq!(a.to_integer(&a.factor_elem(&f)), *y);
            }
        }
    }

    /// A value with a one-limb form: `c/2^64 mod n` for [`Mont`], `c` for [`Plain`].
    fn check_small<A: Arith>(a: &A, x: &Integer, rand: &mut RandState<'_>) {
        let n = a.modulus();
        let c = a.small(x).expect("x must have a one-limb form");
        let y = n.clone().random_below(rand);
        let mut r = a.zero();
        a.mul_small(&mut r, &a.residue(&y), c);
        assert_eq!(a.to_integer(&r), Integer::from(&y * x) % n);
        assert!(matches!(a.factor(x), Factor::Small(_)) || *x <= 2);
    }

    #[test]
    fn mont_matches_integers() {
        let mut rand = RandState::new();
        for bits in [2, 5, 63, 64, 65, 127, 128, 200, 256, 511, 512, 1000, 1024] {
            for _ in 0..3 {
                let mut n = Integer::from(Integer::random_bits(bits, &mut rand));
                n.set_bit(0, true);
                n.set_bit(bits - 1, true);
                let limbs = mont_limbs(&n);
                assert_eq!(limbs, bits.div_ceil(64) as usize);
                with_arith!(&n, |a| check(&a, &mut rand));
                // sigma^2 / 2^64 with sigma < 2^32 has the one-limb form sigma^2 (mod n).
                let sigma = Integer::from(Integer::random_bits(32, &mut rand));
                let inv = Integer::from(Integer::u_pow_u(2, 64)).invert(&n).unwrap();
                let d = Integer::from(&sigma * &sigma) * inv % &n;
                with_arith!(&n, |a| check_small(&a, &d, &mut rand));
            }
        }
    }

    #[test]
    fn plain_matches_integers() {
        let mut rand = RandState::new();
        for n in [
            Integer::from(1_000_000u32),
            Integer::from(Integer::u_pow_u(2, 1100)) + 15u32,
            Integer::from(Integer::u_pow_u(3, 1000)),
        ] {
            assert_eq!(mont_limbs(&n), 0);
            let a = Plain::new(&n);
            check(&a, &mut rand);
            check_small(&a, &Integer::from(12345), &mut rand);
        }
    }
}
