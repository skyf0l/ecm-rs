//! Modular arithmetic for the curve computations.
//!
//! Residues modulo an odd `n` are kept in Montgomery representation (`x` is stored as `x*R mod
//! n` with `R = 2^(64*limbs)`), with no division and no allocation in the operations:
//!
//! - up to [`MAX_LIMBS`] 64-bit limbs, in fixed-size limb arrays ([`Mont`]): a modular
//!   multiplication is one interleaved multiply-and-reduce (CIOS) pass, or from [`GMP_LIMBS`]
//!   limbs GMP's product and Montgomery reduction;
//! - above, in vectors of runtime length ([`MontLarge`]): GMP's product, then GMP's Montgomery
//!   reduction, one limb at a time, two from [`REDC_2_LIMBS`] limbs, and subquadratic from
//!   [`REDC_N_LIMBS`] limbs (GMP-ECM's `MOD_MODMULN`, and GMP's reduction for large moduli).
//!
//! Even moduli (and the large ones if GMP's limbs are not 64-bit, see [`mpn::ENABLED`]) fall
//! back to [`Plain`] big integer arithmetic, with a division per product.
//!
//! The representation never leaks: values go in with [`Arith::residue`] and out with
//! [`Arith::to_integer`], and [`Arith::gcd`] gives `gcd(x, n)` directly (`R` is coprime to an
//! odd `n`).

use std::cell::RefCell;

use rug::{Assign, Integer, integer::Order};

pub use crate::config::MAX_LIMBS;
use crate::config::{GMP_LIMBS, REDC_2_LIMBS, REDC_N_LIMBS};

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
    #[inline(always)]
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
///
/// With a form `$base2` (an `Option<Base2Form>`, see [`crate::base2`]), the special reduction
/// modulo `2^k +- 1` of [`crate::base2::Base2`] if it is `Some`.
macro_rules! with_arith {
    ($n:expr, $base2:expr, |$a:ident| $body:expr) => {{
        match $base2 {
            Some(form) => {
                let $a = $crate::base2::Base2::new($n, form);
                $body
            }
            None => $crate::arith::with_arith!($n, |$a| $body),
        }
    }};
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
            0 => {
                let $a = $crate::arith::Plain::new(n);
                $body
            }
            _ => {
                let $a = $crate::arith::MontLarge::new(n);
                $body
            }
        }
    }};
}
pub(crate) use with_arith;

/// Number of limbs of the Montgomery implementation for `n` ([`Mont`] up to [`MAX_LIMBS`],
/// [`MontLarge`] above), or `0` if `n` needs [`Plain`] (even, or large without
/// [`mpn::ENABLED`]).
pub fn mont_limbs(n: &Integer) -> usize {
    let limbs = n.significant_bits().div_ceil(64) as usize;
    if n.is_odd() && *n > 1 && (limbs <= MAX_LIMBS || mpn::ENABLED) {
        limbs
    } else {
        0
    }
}

/// `x mod n`, in `[0, n)`.
pub(crate) fn reduce(x: &Integer, n: &Integer) -> Integer {
    let mut r = Integer::from(x % n);
    if r < 0 {
        r += n;
    }
    r
}

/// Whether the CPU has BMI2 (`mulx`) and ADX (`adcx`, `adox`).
///
/// The limb arithmetic is generic Rust, compiled for the baseline `x86_64` unless the build
/// enables more (`-C target-cpu`): the hottest loops (the stage 1 ladder, the pairs of stage 2)
/// also have a copy compiled with these extensions, called when this returns `true`. There,
/// `mulx` leaves the carry flag alone, so the carry chains of the Montgomery multiplications
/// stay in the flags: 10-13% faster at 2-8 limbs. That gain needs the default build profile:
/// with `codegen-units = 1` (or fat LTO) LLVM rewrites the chains with `setb` spills, and both
/// copies run at the same speed.
///
/// The standard library caches the detection: this is two relaxed atomic loads.
#[cfg(target_arch = "x86_64")]
#[inline]
pub fn has_bmi2_adx() -> bool {
    #[cfg(test)]
    if GENERIC_ONLY.get() {
        return false;
    }
    crate::config::bmi2_adx_detected()
}

#[cfg(all(test, target_arch = "x86_64"))]
thread_local! {
    /// Tests only: makes [`has_bmi2_adx`] return `false` on this thread, to run the generic
    /// copies on CPUs that have the extensions.
    pub static GENERIC_ONLY: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// `lo + a * b + carry`, as `(low limb, high limb)`: never overflows.
#[inline(always)]
pub(crate) fn mac(lo: u64, a: u64, b: u64, carry: u64) -> (u64, u64) {
    let t = u128::from(lo) + u128::from(a) * u128::from(b) + u128::from(carry);
    (t as u64, (t >> 64) as u64)
}

/// `a + b + carry` (carry is 0 or 1), as `(sum, carry out)`.
#[inline(always)]
pub(crate) fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let (s, c1) = a.overflowing_add(b);
    let (s, c2) = s.overflowing_add(carry);
    (s, u64::from(c1 | c2))
}

/// `a - b - borrow` (borrow is 0 or 1), as `(difference, borrow out)`.
#[inline(always)]
pub(crate) fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
    let (d, b1) = a.overflowing_sub(b);
    let (d, b2) = d.overflowing_sub(borrow);
    (d, u64::from(b1 | b2))
}

/// Montgomery arithmetic modulo an odd `n` of at most `N` limbs.
#[derive(Debug, Clone)]
pub struct Mont<const N: usize> {
    /// Limbs of `n`, least significant first.
    m: [u64; N],
    /// `-1/n mod 2^64`.
    ninv: u64,
    /// `2^64*R mod n`: multiplying by it maps a residue to the polynomial representation.
    to_poly: [u64; N],
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
        let mut to_poly = Integer::from(1) << (64 * N as u32 + 64);
        to_poly %= n;
        Self {
            m,
            ninv: inv.wrapping_neg(),
            to_poly: Self::limbs(&to_poly),
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

    /// Montgomery multiplication, CIOS method: `r = a*b/R mod n`.
    #[inline(always)]
    fn cios(&self, r: &mut [u64; N], a: &[u64; N], b: &[u64; N]) {
        let m = &self.m;
        // t = (t_hi, t_n, t[N-1..0]) < 2n at the end of every iteration.
        let mut t = [0u64; N];
        let mut t_n = 0u64;
        for &bi in b {
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

    /// `r = a*b/R mod n` (`r = a^2/R mod n` if `b` is `None`) with GMP: a product, then a
    /// Montgomery reduction.
    ///
    /// Out of line: at these sizes, a call costs nothing next to the operation, and keeps the
    /// (inlined) callers small.
    #[inline(never)]
    fn gmp_mul(&self, r: &mut [u64; N], a: &[u64; N], b: Option<&[u64; N]>) {
        let mut wide = [[0u64; N]; 2];
        let wide = wide.as_flattened_mut();
        match b {
            Some(b) => mpn::mul(wide, a, b),
            None => mpn::sqr(wide, a),
        }
        let mut t = [0u64; N];
        // a*b < n^2 < n*R, so t + carry*R = (a*b + q*n)/R < 2n.
        let carry = mpn::redc(&mut t, wide, &self.m, self.ninv);
        self.reduce_once(r, &t, carry);
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

    /// Montgomery multiplication: `r = a*b/R mod n`, with GMP from [`GMP_LIMBS`] limbs.
    #[inline(always)]
    fn mul(&self, r: &mut [u64; N], a: &[u64; N], b: &[u64; N]) {
        if mpn::ENABLED && N >= GMP_LIMBS {
            self.gmp_mul(r, a, Some(b));
        } else {
            self.cios(r, a, b);
        }
    }

    /// For mid sizes, separate squaring (half the products of a multiplication) and reduction
    /// are faster than CIOS; from [`GMP_LIMBS`], GMP's squaring and reduction are faster.
    #[inline(always)]
    fn sqr(&self, r: &mut [u64; N], a: &[u64; N]) {
        if mpn::ENABLED && N >= GMP_LIMBS {
            return self.gmp_mul(r, a, None);
        }
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

    #[inline(always)]
    fn add(&self, r: &mut [u64; N], a: &[u64; N], b: &[u64; N]) {
        let mut s = [0; N];
        let mut carry = 0;
        for j in 0..N {
            (s[j], carry) = adc(a[j], b[j], carry);
        }
        self.reduce_once(r, &s, carry);
    }

    #[inline(always)]
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

    #[inline(always)]
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

/// `x^e` for `e > 0`, by left-to-right sliding windows: about one squaring per bit of `e`,
/// and one multiplication per window.
pub(crate) fn pow<A: Arith>(a: &A, x: &A::Elem, e: &Integer) -> A::Elem {
    assert!(*e > 0, "positive exponent");
    let bits = e.significant_bits();
    // Window width: the odd powers below 2^width cost 2^(width - 1) multiplications, a window
    // saves bits/(width + 1) - bits/(width + 2) of them.
    let width = match bits {
        0..64 => 2,
        64..512 => 4,
        512..8192 => 5,
        8192..65536 => 6,
        _ => 7,
    };
    // Odd powers x, x^3, ..., x^(2^width - 1).
    let mut x2 = a.zero();
    a.sqr(&mut x2, x);
    let mut odd = vec![x.clone()];
    for i in 1..1usize << (width - 1) {
        let mut t = a.zero();
        a.mul(&mut t, &odd[i - 1], &x2);
        odd.push(t);
    }
    let (mut r, mut t) = (a.zero(), a.zero());
    let mut first = true;
    // Bits i and below remain.
    let mut i = i64::from(bits) - 1;
    while i >= 0 {
        if !e.get_bit(i as u32) {
            a.sqr(&mut t, &r);
            std::mem::swap(&mut r, &mut t);
            i -= 1;
            continue;
        }
        // The window i..=j, ending with a one.
        let mut j = (i - i64::from(width) + 1).max(0);
        while !e.get_bit(j as u32) {
            j += 1;
        }
        let mut value = 0;
        for b in (j..=i).rev() {
            value = (value << 1) | usize::from(e.get_bit(b as u32));
        }
        if first {
            r.clone_from(&odd[value >> 1]);
            first = false;
        } else {
            for _ in j..=i {
                a.sqr(&mut t, &r);
                std::mem::swap(&mut r, &mut t);
            }
            a.mul(&mut t, &r, &odd[value >> 1]);
            std::mem::swap(&mut r, &mut t);
        }
        i = j - 1;
    }
    r
}

/// GMP's low-level (`mpn`) functions on 64-bit limbs, for the large sizes of [`Mont`] and the
/// Kronecker products of [`crate::poly`].
pub(crate) mod mpn {
    use gmp_mpfr_sys::gmp;

    /// Whether GMP's limbs are our 64-bit limbs: if not, these functions must not be called.
    pub const ENABLED: bool = crate::config::MPN_ENABLED;

    unsafe extern "C" {
        /// `mpn_redc_1(rp, up, mp, n, invm)`: Montgomery reduction by `n` limbs of the `2n`
        /// limbs `up` (clobbered) modulo the `n` limbs `mp`, with `invm = -1/mp[0] mod 2^64`.
        /// Writes `n` limbs to `rp` and returns the carry out of them.
        ///
        /// Internal to GMP (not in `gmp.h`, declared `__GMP_DECLSPEC` in `gmp-impl.h`), but
        /// exported with this signature by every build of GMP since 5.1 (before, it returned
        /// nothing), fat builds included. `gmp-mpfr-sys` builds GMP 6.3, or links a system GMP 6
        /// of minor version at least 3. A public equivalent (a loop of `mpn_addmul_1`, then
        /// `mpn_add_n`) makes stage 1 6-10% slower from 11 limbs.
        #[link_name = "__gmpn_redc_1"]
        fn mpn_redc_1(
            rp: *mut gmp::limb_t,
            up: *mut gmp::limb_t,
            mp: *const gmp::limb_t,
            n: gmp::size_t,
            invm: gmp::limb_t,
        ) -> gmp::limb_t;

        /// `mpn_redc_2(rp, up, mp, n, mip)`: as `mpn_redc_1`, two limbs at a time, with the two
        /// limbs `mip = -1/mp mod 2^128`. Returns the carry out of the `n` limbs of `rp`.
        ///
        /// Internal to GMP like `mpn_redc_1` (`__GMP_DECLSPEC` in `gmp-impl.h`), with this
        /// signature since GMP 5.1 (as `mpn_redc_1`); assembly on x86_64, C elsewhere.
        #[link_name = "__gmpn_redc_2"]
        fn mpn_redc_2(
            rp: *mut gmp::limb_t,
            up: *mut gmp::limb_t,
            mp: *const gmp::limb_t,
            n: gmp::size_t,
            mip: *const gmp::limb_t,
        ) -> gmp::limb_t;

        /// `mpn_redc_n(rp, up, mp, n, ip)`: subquadratic Montgomery reduction of the `2n` limbs
        /// `up` modulo the `n > 8` limbs `mp`, with the `n` limbs `ip = 1/mp mod B^n` (`B =
        /// 2^64`; the opposite sign of `mpn_redc_1`'s inverse): `q = up*ip mod B^n` by a low
        /// half product (`mpn_mullo_n`), then the high half of `q*mp` by a wrap-around product
        /// (`mpn_mulmod_bnm1`, whose wrapped part is known from the low half of `up`), and `rp
        /// = up/B^n - (q*mp)/B^n`, plus `mp` if that is negative: `rp` is `(up - q*mp)/B^n`
        /// modulo `mp`, below `mp` if `up < mp*B^n`. `up` is left unchanged.
        ///
        /// Internal to GMP like `mpn_redc_1` (`__GMP_DECLSPEC` in `gmp-impl.h`), with this
        /// signature since GMP 5.0; generic C code (the one of `mpz_powm` for large moduli), so
        /// in every build. Its temporary space is on the stack (`TMP_ALLOC`) at these sizes.
        #[link_name = "__gmpn_redc_n"]
        fn mpn_redc_n(
            rp: *mut gmp::limb_t,
            up: *mut gmp::limb_t,
            mp: *const gmp::limb_t,
            n: gmp::size_t,
            ip: *const gmp::limb_t,
        );

        /// `mpn_mulmod_bnm1(rp, rn, ap, an, bp, bn, tp)`: `{ap, an} * {bp, bn} mod (B^rn - 1)`
        /// (`B = 2^64`) to the `min(rn, an + bn)` limbs `rp`, for `0 < bn <= an <= rn` and
        /// `an + bn > rn/2`, with the scratch space `tp` of `2*rn + 4` limbs. A non-zero
        /// multiple of `B^rn - 1` gives `B^rn - 1`. A wrap-around product: the half of the
        /// FFT sizes of GMP's own multiplication (`mpn_mul` of large numbers is this with
        /// `rn >= an + bn`), so about half the cost of a full product.
        ///
        /// Internal to GMP like `mpn_redc_1` (`__GMP_DECLSPEC` in `gmp-impl.h`), with this
        /// signature since GMP 5.0; generic C code, so in every build.
        #[link_name = "__gmpn_mulmod_bnm1"]
        fn mpn_mulmod_bnm1(
            rp: *mut gmp::limb_t,
            rn: gmp::size_t,
            ap: *const gmp::limb_t,
            an: gmp::size_t,
            bp: *const gmp::limb_t,
            bn: gmp::size_t,
            tp: *mut gmp::limb_t,
        );

        /// `mpn_mulmod_bnm1_next_size(n)`: the smallest `rn >= n` efficient for
        /// `mpn_mulmod_bnm1` (internal, as `mpn_mulmod_bnm1`).
        #[link_name = "__gmpn_mulmod_bnm1_next_size"]
        fn mpn_mulmod_bnm1_next_size(n: gmp::size_t) -> gmp::size_t;

        /// `mpn_sqrmod_bnm1(rp, rn, ap, an, tp)`: `{ap, an}^2 mod (B^rn - 1)` to the
        /// `min(rn, 2*an)` limbs `rp`, for `rn/4 < an <= rn`, with the scratch space `tp` of
        /// `2*rn + 3` limbs at most (`rn + 3`, plus `an` if `an > rn/2`, as `gmp-impl.h`'s
        /// `mpn_sqrmod_bnm1_itch`): the squaring of `mpn_mulmod_bnm1`, with the same
        /// representation of `0`.
        ///
        /// Internal to GMP like `mpn_mulmod_bnm1` (`__GMP_DECLSPEC` in `gmp-impl.h`), with this
        /// signature since GMP 5.0; generic C code, so in every build.
        #[link_name = "__gmpn_sqrmod_bnm1"]
        fn mpn_sqrmod_bnm1(
            rp: *mut gmp::limb_t,
            rn: gmp::size_t,
            ap: *const gmp::limb_t,
            an: gmp::size_t,
            tp: *mut gmp::limb_t,
        );

        /// `mpn_mul_fft(op, pl, n, nl, m, ml, k)`: `{n, nl} * {m, ml} mod (B^pl + 1)` by
        /// Schoenhage-Strassen's FFT of `2^k` pieces, for `pl` a multiple of `2^k` (then
        /// `mpn_fft_next_size(pl, k) = pl`, which it asserts) and `2*pl/2^k` well below `pl`
        /// (`k >= 4`). Writes the `pl` limbs `op` and returns the limb `op[pl]` (`1` only for
        /// the result `B^pl`): the result is below `B^pl + 1`. A squaring if `n == m` and
        /// `nl == ml`. Its temporary space is allocated by GMP (`TMP_ALLOC`).
        ///
        /// Internal to GMP like `mpn_redc_1` (`__GMP_DECLSPEC` in `gmp-impl.h`), with this
        /// signature since GMP 5 (GMP-ECM calls it for the Fermat numbers from `2^32768 +
        /// 1`, `mpmod.c`); generic C code, so in every build.
        #[link_name = "__gmpn_mul_fft"]
        fn mpn_mul_fft(
            op: *mut gmp::limb_t,
            pl: gmp::size_t,
            n: *const gmp::limb_t,
            nl: gmp::size_t,
            m: *const gmp::limb_t,
            ml: gmp::size_t,
            k: std::ffi::c_int,
        ) -> gmp::limb_t;

        /// `mpn_fft_best_k(n, sqr)`: the best `k` of `mpn_mul_fft` for `n` limbs, from GMP's
        /// tuned tables (internal, as `mpn_mul_fft`).
        #[link_name = "__gmpn_fft_best_k"]
        fn mpn_fft_best_k(n: gmp::size_t, sqr: std::ffi::c_int) -> std::ffi::c_int;
    }

    /// `r = a^2 mod (2^(64*rn) - 1)`, `rn = r.len()`, for `rn >= a.len() >= 1`, with
    /// [`mulmod_bnm1`]'s representation of `0`. `scratch` must have `2*rn + 4` limbs.
    #[inline(always)]
    pub fn sqrmod_bnm1(r: &mut [u64], a: &[u64], scratch: &mut [u64]) {
        let rn = r.len();
        assert!(ENABLED && rn >= a.len() && 4 * a.len() > rn && scratch.len() >= 2 * rn + 4);
        // SAFETY: the limbs are 64 bits (ENABLED), the sizes are those mpn_sqrmod_bnm1 requires
        // (checked above), the scratch space has at least the rn + 3 + an limbs it needs, and
        // the result cannot overlap the operand or the scratch space (all borrowed).
        unsafe {
            mpn_sqrmod_bnm1(
                r.as_mut_ptr().cast(),
                rn as gmp::size_t,
                a.as_ptr().cast(),
                a.len() as gmp::size_t,
                scratch.as_mut_ptr().cast(),
            );
        }
        // Only min(rn, 2*an) limbs are written.
        let written = rn.min(2 * a.len());
        r[written..].fill(0);
    }

    /// [`mulmod_bnm1`] with operands of `rn` limbs and the scratch space of `2*rn + 4` limbs
    /// given (no allocation).
    #[inline(always)]
    pub fn mulmod_bnm1_n(r: &mut [u64], a: &[u64], b: &[u64], scratch: &mut [u64]) {
        let rn = r.len();
        assert!(ENABLED && a.len() == rn && b.len() == rn && rn > 0 && scratch.len() >= 2 * rn + 4);
        // SAFETY: as in `mulmod_bnm1` (an = bn = rn), with the 2rn + 4 limbs of scratch space.
        unsafe {
            mpn_mulmod_bnm1(
                r.as_mut_ptr().cast(),
                rn as gmp::size_t,
                a.as_ptr().cast(),
                rn as gmp::size_t,
                b.as_ptr().cast(),
                rn as gmp::size_t,
                scratch.as_mut_ptr().cast(),
            );
        }
    }

    /// `r = a*b mod (2^(64*pl) + 1)` (`a^2` if `b` is `None`), `pl = r.len() - 1`, the limb
    /// `r[pl]` set only for the result `2^(64*pl)`, by an FFT of `2^k` pieces: `pl` must be a
    /// multiple of `2^k`, `k >= 4`.
    pub fn mul_fft(r: &mut [u64], a: &[u64], b: Option<&[u64]>, k: u32) {
        let pl = r.len() - 1;
        let b = b.unwrap_or(a);
        assert!(
            ENABLED
                && (4..=30).contains(&k)
                && pl.is_multiple_of(1 << k)
                && pl >> k >= 1
                && !a.is_empty()
                && !b.is_empty()
        );
        // SAFETY: the limbs are 64 bits (ENABLED), pl is a multiple of 2^k with k >= 4 (checked
        // above: mpn_mul_fft asserts it, and 2*pl/2^k + O(1) < pl), the operands have their
        // given lengths, and `r` cannot overlap them (borrowed mutably). A squaring passes the
        // same pointer and length twice, which mpn_mul_fft detects.
        unsafe {
            r[pl] = mpn_mul_fft(
                r.as_mut_ptr().cast(),
                pl as gmp::size_t,
                a.as_ptr().cast(),
                a.len() as gmp::size_t,
                b.as_ptr().cast(),
                b.len() as gmp::size_t,
                k as std::ffi::c_int,
            ) as u64;
        }
    }

    /// The `k` of [`mul_fft`] for `pl` limbs: GMP's best one (`mpn_fft_best_k`), lowered until
    /// `2^k` divides `pl` (as `mpn_mulmod_bnm1` does), if still at least 4.
    pub fn fft_k(pl: usize, sqr: bool) -> Option<u32> {
        assert!(ENABLED && pl > 0);
        // SAFETY: a pure function of its arguments (ENABLED: the library has it with this
        // signature).
        let best = unsafe { mpn_fft_best_k(pl as gmp::size_t, std::ffi::c_int::from(sqr)) };
        let k = (4..=u32::try_from(best).ok()?.min(30))
            .rev()
            .find(|&k| pl.is_multiple_of(1 << k))?;
        Some(k)
    }

    /// `r = a*b`, with `a.len() >= b.len() >= 1` and `r` of `a.len() + b.len()` limbs.
    pub fn mul_long(r: &mut [u64], a: &[u64], b: &[u64]) {
        assert!(ENABLED && a.len() >= b.len() && !b.is_empty() && r.len() == a.len() + b.len());
        // SAFETY: the limbs are 64 bits (ENABLED), the sizes are those mpn_mul requires
        // (checked above), `r` cannot overlap the operands (it is borrowed mutably).
        unsafe {
            gmp::mpn_mul(
                r.as_mut_ptr().cast(),
                a.as_ptr().cast(),
                a.len() as gmp::size_t,
                b.as_ptr().cast(),
                b.len() as gmp::size_t,
            );
        }
    }

    /// `r = a*b mod (2^(64*rn) - 1)`, `rn = r.len()`, for `rn >= a.len() >= b.len() >= 1` and
    /// `a.len() + b.len() > rn/2`: a non-zero multiple of `2^(64*rn) - 1` gives
    /// `2^(64*rn) - 1`, so the result is exact if the product is known to be below it.
    /// `scratch` is resized as needed.
    pub fn mulmod_bnm1(r: &mut [u64], a: &[u64], b: &[u64], scratch: &mut Vec<u64>) {
        let rn = r.len();
        assert!(
            ENABLED
                && rn >= a.len()
                && a.len() >= b.len()
                && !b.is_empty()
                && a.len() + b.len() > rn / 2
        );
        scratch.resize(2 * rn + 4, 0);
        // SAFETY: the limbs are 64 bits (ENABLED), the sizes are those mpn_mulmod_bnm1
        // requires (checked above), the scratch space has the 2rn + 4 limbs it needs, and the
        // result cannot overlap the operands or the scratch space (all borrowed).
        unsafe {
            mpn_mulmod_bnm1(
                r.as_mut_ptr().cast(),
                rn as gmp::size_t,
                a.as_ptr().cast(),
                a.len() as gmp::size_t,
                b.as_ptr().cast(),
                b.len() as gmp::size_t,
                scratch.as_mut_ptr().cast(),
            );
        }
        // Only min(rn, an + bn) limbs are written.
        let written = rn.min(a.len() + b.len());
        r[written..].fill(0);
    }

    /// The smallest size `>= n` efficient for [`mulmod_bnm1`].
    pub fn mulmod_bnm1_next_size(n: usize) -> usize {
        assert!(ENABLED && n > 0);
        // SAFETY: a pure function of n (ENABLED: the library has it with this signature).
        unsafe { mpn_mulmod_bnm1_next_size(n as gmp::size_t) as usize }
    }

    /// `r = a*b`, with `r` of `2*a.len()` limbs.
    #[inline(always)]
    pub fn mul(r: &mut [u64], a: &[u64], b: &[u64]) {
        assert!(ENABLED && a.len() == b.len() && r.len() == 2 * a.len() && !a.is_empty());
        // SAFETY: the limbs are 64 bits (ENABLED), `r` has room for the 2n limbs of the product
        // and cannot overlap the operands (it is borrowed mutably), and n > 0.
        unsafe {
            gmp::mpn_mul_n(
                r.as_mut_ptr().cast(),
                a.as_ptr().cast(),
                b.as_ptr().cast(),
                a.len() as gmp::size_t,
            );
        }
    }

    /// `r = a^2`, with `r` of `2*a.len()` limbs.
    #[inline(always)]
    pub fn sqr(r: &mut [u64], a: &[u64]) {
        assert!(ENABLED && r.len() == 2 * a.len() && !a.is_empty());
        // SAFETY: as in `mul`.
        unsafe {
            gmp::mpn_sqr(
                r.as_mut_ptr().cast(),
                a.as_ptr().cast(),
                a.len() as gmp::size_t,
            );
        }
    }

    /// `r + carry*2^(64*n) = (t + q*m)/2^(64*n)` for the `q < 2^(64*n)` that makes it exact:
    /// Montgomery reduction of `t` (`2n` limbs, clobbered) modulo the odd `m` (`n` limbs), with
    /// `ninv = -1/m mod 2^64`. Returns `carry`.
    #[inline(always)]
    pub fn redc(r: &mut [u64], t: &mut [u64], m: &[u64], ninv: u64) -> u64 {
        let n = m.len();
        assert!(ENABLED && r.len() == n && t.len() == 2 * n && n > 0 && m[0] & 1 == 1);
        // SAFETY: the limbs are 64 bits (ENABLED), the buffers have the sizes mpn_redc_1
        // expects (checked above) and are distinct (`r` and `t` are borrowed mutably).
        unsafe {
            mpn_redc_1(
                r.as_mut_ptr().cast(),
                t.as_mut_ptr().cast(),
                m.as_ptr().cast(),
                n as gmp::size_t,
                ninv as gmp::limb_t,
            ) as u64
        }
    }

    /// [`redc`] two limbs at a time, with `ninv = -1/m mod 2^128` (least significant limb
    /// first).
    #[inline(always)]
    pub fn redc_2(r: &mut [u64], t: &mut [u64], m: &[u64], ninv: &[u64; 2]) -> u64 {
        let n = m.len();
        assert!(ENABLED && r.len() == n && t.len() == 2 * n && n >= 2 && m[0] & 1 == 1);
        // SAFETY: as in `redc`, and `ninv` has the two limbs mpn_redc_2 reads.
        unsafe {
            mpn_redc_2(
                r.as_mut_ptr().cast(),
                t.as_mut_ptr().cast(),
                m.as_ptr().cast(),
                n as gmp::size_t,
                ninv.as_ptr().cast(),
            ) as u64
        }
    }

    /// `r = (t - q*m)/2^(64*n)` modulo `m` in `[0, 2^(64*n))`, for the `q < 2^(64*n)` that makes
    /// it exact: subquadratic Montgomery reduction of `t` (`2n` limbs, `n > 8`) modulo the odd
    /// `m` (`n` limbs), with `inv = 1/m mod 2^(64*n)` (`n` limbs). Below `m` if `t < m*2^(64*n)`,
    /// at most `t/2^(64*n)` in any case.
    #[inline(always)]
    pub fn redc_n(r: &mut [u64], t: &mut [u64], m: &[u64], inv: &[u64]) {
        let n = m.len();
        assert!(
            ENABLED && r.len() == n && t.len() == 2 * n && inv.len() == n && n > 8 && m[0] & 1 == 1
        );
        // SAFETY: the limbs are 64 bits (ENABLED), the buffers have the sizes mpn_redc_n
        // expects (checked above; it asserts n > 8), and `r` and `t` are distinct from each
        // other and from `m` and `inv` (borrowed mutably).
        unsafe {
            mpn_redc_n(
                r.as_mut_ptr().cast(),
                t.as_mut_ptr().cast(),
                m.as_ptr().cast(),
                n as gmp::size_t,
                inv.as_ptr().cast(),
            );
        }
    }

    /// `r = a + b` (all of the same length), returns the carry out.
    #[inline(always)]
    pub fn add(r: &mut [u64], a: &[u64], b: &[u64]) -> u64 {
        let n = r.len();
        assert!(ENABLED && a.len() == n && b.len() == n && n > 0);
        // SAFETY: the limbs are 64 bits (ENABLED), the three operands have n > 0 limbs, and
        // `r` cannot overlap the operands (it is borrowed mutably).
        unsafe {
            gmp::mpn_add_n(
                r.as_mut_ptr().cast(),
                a.as_ptr().cast(),
                b.as_ptr().cast(),
                n as gmp::size_t,
            ) as u64
        }
    }

    /// `r = a - b` (all of the same length), returns the borrow out.
    #[inline(always)]
    pub fn sub(r: &mut [u64], a: &[u64], b: &[u64]) -> u64 {
        let n = r.len();
        assert!(ENABLED && a.len() == n && b.len() == n && n > 0);
        // SAFETY: as in `add`.
        unsafe {
            gmp::mpn_sub_n(
                r.as_mut_ptr().cast(),
                a.as_ptr().cast(),
                b.as_ptr().cast(),
                n as gmp::size_t,
            ) as u64
        }
    }

    /// `r += b`, for `r.len() >= b.len() >= 1`: returns the carry out of `r`.
    #[inline(always)]
    pub fn add_assign(r: &mut [u64], b: &[u64]) -> u64 {
        assert!(ENABLED && r.len() >= b.len() && !b.is_empty());
        let rp = r.as_mut_ptr().cast();
        // SAFETY: the limbs are 64 bits (ENABLED), the sizes are those mpn_add requires
        // (checked above), the destination is the first source (in place, which GMP allows),
        // and `b` cannot overlap `r` (borrowed mutably).
        unsafe {
            gmp::mpn_add(
                rp,
                rp,
                r.len() as gmp::size_t,
                b.as_ptr().cast(),
                b.len() as gmp::size_t,
            ) as u64
        }
    }

    /// `r -= b` (of the same length), returns the borrow out.
    #[inline(always)]
    pub fn sub_assign(r: &mut [u64], b: &[u64]) -> u64 {
        let n = r.len();
        assert!(ENABLED && b.len() == n && n > 0);
        let rp = r.as_mut_ptr().cast();
        // SAFETY: as in `add_assign`.
        unsafe { gmp::mpn_sub_n(rp, rp, b.as_ptr().cast(), n as gmp::size_t) as u64 }
    }

    /// Whether `a >= b` (of the same length).
    #[inline(always)]
    pub fn ge(a: &[u64], b: &[u64]) -> bool {
        assert!(ENABLED && a.len() == b.len() && !a.is_empty());
        // SAFETY: the limbs are 64 bits (ENABLED), both operands have the n > 0 limbs read.
        unsafe { gmp::mpn_cmp(a.as_ptr().cast(), b.as_ptr().cast(), a.len() as gmp::size_t) >= 0 }
    }

    /// `r = a*c` (of the same length), returns the high limb.
    #[inline(always)]
    pub fn mul_1(r: &mut [u64], a: &[u64], c: u64) -> u64 {
        let n = r.len();
        assert!(ENABLED && a.len() == n && n > 0);
        // SAFETY: as in `add`.
        unsafe {
            gmp::mpn_mul_1(
                r.as_mut_ptr().cast(),
                a.as_ptr().cast(),
                n as gmp::size_t,
                c as gmp::limb_t,
            ) as u64
        }
    }

    /// `r += a*c` (of the same length), returns the carry limb out.
    #[inline(always)]
    pub fn addmul_1(r: &mut [u64], a: &[u64], c: u64) -> u64 {
        let n = r.len();
        assert!(ENABLED && a.len() == n && n > 0);
        // SAFETY: as in `add`.
        unsafe {
            gmp::mpn_addmul_1(
                r.as_mut_ptr().cast(),
                a.as_ptr().cast(),
                n as gmp::size_t,
                c as gmp::limb_t,
            ) as u64
        }
    }
}

/// Montgomery reduction of [`MontLarge`], by size of the modulus.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Redc {
    /// GMP's `mpn_redc_1`, one limb at a time (quadratic).
    One,
    /// GMP's `mpn_redc_2`, two limbs at a time (quadratic), from [`REDC_2_LIMBS`] limbs.
    Two,
    /// GMP's subquadratic `mpn_redc_n`, from [`REDC_N_LIMBS`] limbs.
    Subquadratic,
}

/// Montgomery arithmetic modulo an odd `n` of more than [`MAX_LIMBS`] limbs, on residues of
/// runtime length: vectors of the `limbs` limbs of `x*R mod n` (`R = 2^(64*limbs)`), in
/// `[0, n)`.
///
/// A modular multiplication is GMP's product (`mpn_mul_n` or `mpn_sqr`), then GMP's Montgomery
/// reduction: `mpn_redc_1` or `mpn_redc_2` (quadratic, GMP-ECM's `MOD_MODMULN`), then the
/// subquadratic `mpn_redc_n` from [`REDC_N_LIMBS`] limbs; no division and no allocation (the
/// buffers are the arithmetic's). Only with GMP's low-level functions ([`mpn::ENABLED`]):
/// without them, these moduli use [`Plain`].
#[derive(Debug)]
pub struct MontLarge {
    /// Limbs of `n`, least significant first.
    m: Vec<u64>,
    /// `-1/n mod 2^128`, least significant limb first (`ninv[0] = -1/n mod 2^64`).
    ninv: [u64; 2],
    /// `1/n mod R`, for [`Redc::Subquadratic`] (empty otherwise).
    inv: Vec<u64>,
    redc: Redc,
    /// `2^64*R mod n`: multiplying by it maps a residue to the polynomial representation.
    to_poly: Vec<u64>,
    /// Products and wide values, of `2*limbs + 2` limbs.
    scratch: RefCell<Vec<u64>>,
    n: Integer,
}

impl MontLarge {
    /// Arithmetic modulo `n`, which must be odd and above `1`, with [`mpn::ENABLED`].
    pub fn new(n: &Integer) -> Self {
        assert!(mpn::ENABLED && n.is_odd() && *n > 1);
        let limbs = n.significant_bits().div_ceil(64) as usize;
        let r = Integer::from(1) << (64 * limbs as u32);
        let inv = n.clone().invert(&r).expect("odd modulus");
        let two128 = Integer::from(1) << 128;
        let ninv: Integer = &two128 - Integer::from(&inv % &two128);
        let mut ninv2 = [0; 2];
        ninv.write_digits(&mut ninv2, Order::Lsf);
        let redc = if limbs >= REDC_N_LIMBS {
            Redc::Subquadratic
        } else if limbs >= REDC_2_LIMBS {
            Redc::Two
        } else {
            Redc::One
        };
        let to_poly = Integer::from(1) << (64 * limbs as u32 + 64);
        let mut a = Self {
            m: Vec::new(),
            ninv: ninv2,
            inv: Vec::new(),
            redc,
            to_poly: Vec::new(),
            scratch: RefCell::new(vec![0; 2 * limbs + 2]),
            n: n.clone(),
        };
        a.m = a.limbs_of(n);
        if redc == Redc::Subquadratic {
            a.inv = a.limbs_of(&inv);
        }
        a.to_poly = a.limbs_of(&(to_poly % n));
        a
    }

    /// Limbs of `0 <= x < R`.
    fn limbs_of(&self, x: &Integer) -> Vec<u64> {
        let mut v = vec![0; self.n.significant_bits().div_ceil(64) as usize];
        x.write_digits(&mut v, Order::Lsf);
        v
    }

    /// `r = r - n` if `r >= n` (with `r = r + carry*R < 2n`).
    #[inline(always)]
    fn reduce_once(&self, r: &mut [u64], carry: u64) {
        if carry != 0 || mpn::ge(r, &self.m) {
            mpn::sub_assign(r, &self.m);
        }
    }

    /// `r = t/R mod n` in `[0, n)`, for `t` of `2*limbs` limbs (clobbered) below `n*(R + 1)`.
    ///
    /// Then `(t + q*n)/R < 2n` for the quadratic reductions (`q < R`), and `(t - q*n)/R` is in
    /// `(-n, n]` for the subquadratic one, which adds `n` if it is negative: one subtraction
    /// of `n` at most.
    #[inline(always)]
    fn redc(&self, r: &mut [u64], t: &mut [u64]) {
        match self.redc {
            Redc::One => {
                let carry = mpn::redc(r, t, &self.m, self.ninv[0]);
                self.reduce_once(r, carry);
            }
            Redc::Two => {
                let carry = mpn::redc_2(r, t, &self.m, &self.ninv);
                self.reduce_once(r, carry);
            }
            Redc::Subquadratic => {
                mpn::redc_n(r, t, &self.m, &self.inv);
                self.reduce_once(r, 0);
            }
        }
    }

    /// `r = a*b/R mod n` (`r = a^2/R mod n` if `b` is `None`).
    #[inline(always)]
    fn mul_or_sqr(&self, r: &mut [u64], a: &[u64], b: Option<&[u64]>) {
        let mut scratch = self.scratch.borrow_mut();
        let wide = &mut scratch[..2 * self.m.len()];
        match b {
            Some(b) => mpn::mul(wide, a, b),
            None => mpn::sqr(wide, a),
        }
        self.redc(r, wide);
    }
}

impl Arith for MontLarge {
    type Elem = Vec<u64>;

    fn modulus(&self) -> &Integer {
        &self.n
    }

    fn zero(&self) -> Vec<u64> {
        vec![0; self.m.len()]
    }

    fn residue(&self, x: &Integer) -> Vec<u64> {
        let x = reduce(x, &self.n) << (64 * self.m.len() as u32);
        self.limbs_of(&(x % &self.n))
    }

    fn to_integer(&self, x: &Vec<u64>) -> Integer {
        let mut r = self.zero();
        {
            let mut scratch = self.scratch.borrow_mut();
            let wide = &mut scratch[..2 * self.m.len()];
            wide.fill(0);
            wide[..x.len()].copy_from_slice(x);
            self.redc(&mut r, wide);
        }
        Integer::from_digits(&r, Order::Lsf)
    }

    fn gcd(&self, x: &Vec<u64>) -> Integer {
        Integer::from_digits(x, Order::Lsf).gcd(&self.n)
    }

    #[inline(always)]
    fn mul(&self, r: &mut Vec<u64>, a: &Vec<u64>, b: &Vec<u64>) {
        self.mul_or_sqr(r, a, Some(b));
    }

    #[inline(always)]
    fn sqr(&self, r: &mut Vec<u64>, a: &Vec<u64>) {
        self.mul_or_sqr(r, a, None);
    }

    #[inline(always)]
    fn add(&self, r: &mut Vec<u64>, a: &Vec<u64>, b: &Vec<u64>) {
        let carry = mpn::add(r, a, b);
        self.reduce_once(r, carry);
    }

    #[inline(always)]
    fn sub(&self, r: &mut Vec<u64>, a: &Vec<u64>, b: &Vec<u64>) {
        if mpn::sub(r, a, b) != 0 {
            // Add n back: the carry out cancels the borrow.
            mpn::add_assign(r, &self.m);
        }
    }

    /// As [`Mont`]: the one-limb form of `x` is `c = x*2^64 mod n`, if `c < 2^64`.
    fn small(&self, x: &Integer) -> Option<u64> {
        (Integer::from(x << 64) % &self.n).to_u64()
    }

    /// Montgomery reduction by one limb of `a*c`: `(a*c + q*n)/2^64 < 2n`.
    #[inline(always)]
    fn mul_small(&self, r: &mut Vec<u64>, a: &Vec<u64>, c: u64) {
        let limbs = self.m.len();
        let mut scratch = self.scratch.borrow_mut();
        let t = &mut scratch[..limbs + 2];
        t[limbs] = mpn::mul_1(&mut t[..limbs], a, c);
        t[limbs + 1] = 0;
        let q = t[0].wrapping_mul(self.ninv[0]);
        let carry = mpn::addmul_1(&mut t[..limbs], &self.m, q);
        mpn::add_assign(&mut t[limbs..], &[carry]);
        r.copy_from_slice(&t[1..=limbs]);
        self.reduce_once(r, t[limbs + 1]);
    }
}

/// Plain big integer arithmetic modulo any `n`: residues are integers in `[0, n)`.
///
/// Fallback for the moduli [`Mont`] and [`MontLarge`] don't handle.
#[derive(Debug, Clone)]
pub struct Plain {
    n: Integer,
}

impl Plain {
    /// Arithmetic modulo `n > 0`.
    pub fn new(n: &Integer) -> Self {
        Self { n: n.clone() }
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

/// Arithmetic on the coefficients of polynomials (see [`crate::poly`]), in a representation of
/// their own where products are formed on plain limbs and reduced later: Kronecker substitution
/// packs whole polynomials into big integers, so a coefficient of a product is a sum of many
/// products of residues.
///
/// The polynomial representation of `x` is `x*R' mod n` (`R' = 2^(64*(limbs + 1))` for [`Mont`]
/// and [`MontLarge`], `R' = 1` for [`Plain`]), a value in `[0, n)`: [`PolyArith::redc_wide`] divides by `R'`, so it
/// maps the product of two such values back to the representation of the product. For
/// [`crate::base2::Base2`], `R' = 1` and the value is the residue modulo `M = 2^k +- 1`, below
/// [`PolyArith::value_bits`] (not `n`).
pub trait PolyArith: Arith {
    /// Number of limbs of the values: they are below `2^(64*limbs)`.
    fn limbs(&self) -> usize;
    /// Bits of the values: they are below `2^value_bits` (those of `n` by default).
    fn value_bits(&self) -> usize {
        self.modulus().significant_bits() as usize
    }
    /// `r` = polynomial representation of the residue `x`.
    fn to_poly(&self, r: &mut Self::Elem, x: &Self::Elem);
    /// Writes the limbs of the value of `x` (polynomial representation) to `out`, of length
    /// [`PolyArith::limbs`].
    fn write_limbs(&self, x: &Self::Elem, out: &mut [u64]);
    /// `acc += x*y` (values of the polynomial representation), `acc` having at least
    /// `2*limbs + 1` limbs (the carry out of `acc` is lost).
    fn mul_acc(&self, acc: &mut [u64], x: &Self::Elem, y: &Self::Elem);
    /// `r = t/R' mod n`, for `t < n*R'` (any `t` for `Base2`) given by at most `2*limbs + 2`
    /// limbs.
    fn redc_wide(&self, r: &mut Self::Elem, t: &[u64]);

    /// `r = x*y` in the polynomial representation.
    fn poly_mul(&self, r: &mut Self::Elem, x: &Self::Elem, y: &Self::Elem) {
        let mut acc = vec![0; 2 * self.limbs() + 2];
        self.mul_acc(&mut acc, x, y);
        self.redc_wide(r, &acc);
    }

    /// Polynomial representation of `x` (any integer).
    fn poly_from(&self, x: &Integer) -> Self::Elem {
        let mut r = self.zero();
        self.to_poly(&mut r, &self.residue(x));
        r
    }

    /// Value in `[0, n)` of `x`, in the polynomial representation.
    #[cfg(test)]
    fn poly_value(&self, x: &Self::Elem) -> Integer {
        let mut limbs = vec![0; self.limbs()];
        self.write_limbs(x, &mut limbs);
        let mut r = self.zero();
        self.redc_wide(&mut r, &limbs);
        self.write_limbs(&r, &mut limbs);
        Integer::from_digits(&limbs, Order::Lsf) % self.modulus()
    }
}

/// `acc += x*y` on limbs, `acc` having at least `x.len() + y.len()` limbs; returns the carry
/// out of `acc`.
pub(crate) fn mul_acc_limbs(acc: &mut [u64], x: &[u64], y: &[u64]) -> u64 {
    let mut top = 0;
    for (i, &yi) in y.iter().enumerate() {
        let mut c = 0;
        for (j, &xj) in x.iter().enumerate() {
            (acc[i + j], c) = mac(acc[i + j], xj, yi, c);
        }
        // Propagate the carry to the end of acc.
        for a in &mut acc[i + x.len()..] {
            if c == 0 {
                break;
            }
            (*a, c) = adc(*a, c, 0);
        }
        top += c;
    }
    top
}

impl<const N: usize> PolyArith for Mont<N> {
    fn limbs(&self) -> usize {
        N
    }

    fn to_poly(&self, r: &mut [u64; N], x: &[u64; N]) {
        self.mul(r, x, &self.to_poly);
    }

    fn write_limbs(&self, x: &[u64; N], out: &mut [u64]) {
        out.copy_from_slice(x);
    }

    #[inline]
    fn mul_acc(&self, acc: &mut [u64], x: &[u64; N], y: &[u64; N]) {
        // Product on 2N limbs, then added to acc.
        let mut p = [0u64; 2 * MAX_LIMBS];
        let p = &mut p[..2 * N];
        for i in 0..N {
            let mut c = 0;
            for j in 0..N {
                (p[i + j], c) = mac(p[i + j], x[j], y[i], c);
            }
            p[i + N] = c;
        }
        let mut carry = 0;
        for (a, &p) in acc.iter_mut().zip(p.iter()) {
            (*a, carry) = adc(*a, p, carry);
        }
        for a in &mut acc[2 * N..] {
            (*a, carry) = adc(*a, 0, carry);
        }
    }

    /// Montgomery reduction by `N + 1` limbs: `(t + q*n)/R' < t/R' + n < 2n`.
    fn redc_wide(&self, r: &mut [u64; N], t: &[u64]) {
        let m = &self.m;
        let mut buf = [0u64; 2 * MAX_LIMBS + 2];
        let buf = &mut buf[..2 * N + 2];
        buf[..t.len()].copy_from_slice(t);
        // Step i clears limb i; hi is the carry into limb i + N + 1.
        let mut hi = 0;
        for i in 0..=N {
            let q = buf[i].wrapping_mul(self.ninv);
            let mut c = 0;
            for j in 0..N {
                (buf[i + j], c) = mac(buf[i + j], q, m[j], c);
            }
            let (s, c1) = adc(buf[i + N], c, 0);
            let (s, c2) = adc(s, hi, 0);
            buf[i + N] = s;
            hi = c1 | c2;
        }
        let top = buf[2 * N + 1] + hi;
        let mut v = [0; N];
        v.copy_from_slice(&buf[N + 1..2 * N + 1]);
        self.reduce_once(r, &v, top);
    }
}

impl PolyArith for MontLarge {
    fn limbs(&self) -> usize {
        self.m.len()
    }

    fn to_poly(&self, r: &mut Vec<u64>, x: &Vec<u64>) {
        self.mul(r, x, &self.to_poly);
    }

    fn write_limbs(&self, x: &Vec<u64>, out: &mut [u64]) {
        out.copy_from_slice(x);
    }

    fn mul_acc(&self, acc: &mut [u64], x: &Vec<u64>, y: &Vec<u64>) {
        let mut scratch = self.scratch.borrow_mut();
        let p = &mut scratch[..2 * self.m.len()];
        mpn::mul(p, x, y);
        mpn::add_assign(acc, p);
    }

    /// Montgomery reduction by one limb, then by `limbs` limbs: for `t < n*R'`, `(t +
    /// q*n)/2^64 < n*R + n` (`q < 2^64`) fits in `2*limbs` limbs, as [`MontLarge::redc`]
    /// needs.
    fn redc_wide(&self, r: &mut Vec<u64>, t: &[u64]) {
        let limbs = self.m.len();
        let mut scratch = self.scratch.borrow_mut();
        let buf = &mut scratch[..2 * limbs + 2];
        buf[..t.len()].copy_from_slice(t);
        buf[t.len()..].fill(0);
        let q = buf[0].wrapping_mul(self.ninv[0]);
        let carry = mpn::addmul_1(&mut buf[..limbs], &self.m, q);
        mpn::add_assign(&mut buf[limbs..], &[carry]);
        debug_assert_eq!(buf[2 * limbs + 1], 0);
        self.redc(r, &mut buf[1..2 * limbs + 1]);
    }
}

impl PolyArith for Plain {
    fn limbs(&self) -> usize {
        self.n.significant_bits().div_ceil(64) as usize
    }

    fn to_poly(&self, r: &mut Integer, x: &Integer) {
        r.assign(x);
    }

    fn write_limbs(&self, x: &Integer, out: &mut [u64]) {
        x.write_digits(out, Order::Lsf);
    }

    fn mul_acc(&self, acc: &mut [u64], x: &Integer, y: &Integer) {
        let limbs = |v: &Integer| {
            let mut digits = vec![0; v.significant_digits::<u64>()];
            v.write_digits(&mut digits, Order::Lsf);
            digits
        };
        mul_acc_limbs(acc, &limbs(x), &limbs(y));
    }

    fn redc_wide(&self, r: &mut Integer, t: &[u64]) {
        r.assign_digits(t, Order::Lsf);
        *r %= &self.n;
    }
}

/// A chain of modular multiplications or squarings on fixed residues, for benchmarks.
///
/// Dispatched like the curve code (Montgomery on fixed-size arrays up to 16 limbs, on vectors
/// above, plain integers for even moduli, the special reduction for the divisors of `2^k +- 1`
/// it applies to): the setup (conversion to the
/// internal representation) is done by [`ArithBatch::new`], and [`ArithBatch::run`] only
/// computes.
#[cfg(feature = "bench")]
pub struct ArithBatch(Box<dyn BatchRun>);

#[cfg(feature = "bench")]
trait BatchRun {
    fn run(&self, ops: usize, square: bool) -> Integer;
}

#[cfg(feature = "bench")]
struct Batch<A: Arith> {
    arith: A,
    values: Vec<A::Elem>,
}

#[cfg(feature = "bench")]
impl<A: Arith> BatchRun for Batch<A> {
    fn run(&self, ops: usize, square: bool) -> Integer {
        let a = &self.arith;
        let mut acc = self.values[0].clone();
        let mut t = a.zero();
        for i in 0..ops {
            if square {
                a.sqr(&mut t, &acc);
            } else {
                a.mul(&mut t, &acc, &self.values[i % self.values.len()]);
            }
            std::mem::swap(&mut acc, &mut t);
        }
        a.to_integer(&acc)
    }
}

#[cfg(feature = "bench")]
impl ArithBatch {
    /// Arithmetic modulo `n` (the implementation the curve code would use), on the residues of
    /// `values`.
    ///
    /// # Panics
    ///
    /// If `values` is empty.
    #[must_use]
    pub fn new(n: &Integer, values: &[Integer]) -> Self {
        assert!(!values.is_empty());
        with_arith!(n, crate::base2::Base2Form::detect(n), |arith| {
            let values = values.iter().map(|x| arith.residue(x)).collect();
            Self(Box::new(Batch { arith, values }))
        })
    }

    /// Number of limbs of the Montgomery implementation for `n` (on fixed-size arrays up to 16
    /// limbs, on vectors above), or `0` for plain integers (the special reduction may be used
    /// instead, see [`ArithBatch::new`]).
    #[must_use]
    pub fn limbs(n: &Integer) -> usize {
        mont_limbs(n)
    }

    /// `ops` chained multiplications `acc = acc * values[i % len]` (or squarings `acc = acc^2`)
    /// from `acc = values[0]`: returns the final `acc`.
    #[must_use]
    pub fn run(&self, ops: usize, square: bool) -> Integer {
        self.0.run(ops, square)
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
        for bits in [
            2, 5, 63, 64, 65, 127, 128, 200, 256, 511, 512, 640, 641, 703, 704, 1000, 1024,
        ] {
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

    /// Random value below `2^bits` whose limbs are mostly `0`, `1` or all ones: carry and
    /// borrow edge cases that uniform random values almost never hit.
    fn adversarial(bits: u32, rand: &mut RandState<'_>) -> Integer {
        let mut x = Integer::new();
        for i in 0..bits.div_ceil(64) {
            let limb = match rand.bits(2) {
                0 => Integer::ZERO,
                1 => Integer::from(1),
                2 => Integer::from(u64::MAX),
                _ => Integer::from(Integer::random_bits(64, rand)),
            };
            x += limb << (64 * i);
        }
        x.keep_bits(bits)
    }

    /// Limb counts of the edge case tests: every one of [`Mont`], then [`MontLarge`] around
    /// its thresholds.
    fn edge_limbs() -> Vec<u32> {
        let mut limbs: Vec<u32> = (1..=MAX_LIMBS as u32 + 2).collect();
        for t in [REDC_2_LIMBS, REDC_N_LIMBS] {
            limbs.extend([t as u32 - 1, t as u32]);
        }
        limbs
    }

    #[test]
    fn mont_edge_cases() {
        let mut rand = RandState::new();
        for limbs in edge_limbs() {
            let bits = 64 * limbs;
            let r = Integer::from(Integer::u_pow_u(2, bits));
            let mut moduli = vec![
                // Largest modulus, all limbs ones.
                Integer::from(&r - 1),
                Integer::from(&r - 3),
                // Smallest modulus of this size: top limb 1.
                Integer::from(Integer::u_pow_u(2, bits - 64)) + 1,
                // Top limb 1, low limbs all ones.
                Integer::from(Integer::u_pow_u(2, bits - 63)) - 1,
                // Just above R/2.
                Integer::from(&r >> 1) + 1,
            ];
            moduli.extend((0..8).map(|_| adversarial(bits, &mut rand) | 1u32));
            for n in moduli {
                if n <= 1 {
                    continue;
                }
                let size = n.significant_bits().div_ceil(64) as usize;
                let expected = if n.is_odd() && (size <= MAX_LIMBS || mpn::ENABLED) {
                    size
                } else {
                    0
                };
                assert_eq!(mont_limbs(&n), expected);
                with_arith!(&n, |a| check_edges(&a, &mut rand));
            }
        }
    }

    /// Checks [`Arith::mul`], [`Arith::sqr`] and [`Arith::mul_small`] on residues near `0`,
    /// near `n` and with limbs all ones.
    fn check_edges<A: Arith>(a: &A, rand: &mut RandState<'_>) {
        let n = a.modulus().clone();
        let bits = n.significant_bits();
        let mut values = vec![
            Integer::ZERO,
            Integer::from(1),
            Integer::from(2),
            Integer::from(&n - 1),
            Integer::from(&n - 2),
            Integer::from(&n >> 1),
            Integer::from(&n >> 1) + 1u32,
        ];
        values.extend((0..10).map(|_| adversarial(bits, rand) % &n));
        values.extend((0..5).map(|_| n.clone().random_below(rand)));
        for x in &mut values {
            *x = reduce(x, &n);
        }
        let mut r = a.zero();
        for x in &values {
            let ex = a.residue(x);
            assert_eq!(a.to_integer(&ex), *x, "residue {x} mod {n}");
            a.sqr(&mut r, &ex);
            assert_eq!(a.to_integer(&r), Integer::from(x * x) % &n, "{x}^2 mod {n}");
            for c in [1, 2, u64::MAX - 1, u64::MAX, u64::from(rand.bits(32))] {
                // One-limb form `c` of `c/2^64` (Mont) or of `c` (Plain).
                let c_int = Integer::from(c);
                for y in [&c_int * inverse_r(&n) % &n, c_int % &n] {
                    if a.small(&y) == Some(c) {
                        a.mul_small(&mut r, &ex, c);
                        let expected = Integer::from(x * &y) % &n;
                        assert_eq!(a.to_integer(&r), expected, "{x}*{y} mod {n}");
                    }
                }
            }
            for y in &values {
                let ey = a.residue(y);
                a.mul(&mut r, &ex, &ey);
                assert_eq!(
                    a.to_integer(&r),
                    Integer::from(x * y) % &n,
                    "{x}*{y} mod {n}"
                );
                a.add(&mut r, &ex, &ey);
                assert_eq!(a.to_integer(&r), Integer::from(x + y) % &n);
                a.sub(&mut r, &ex, &ey);
                assert_eq!(a.to_integer(&r), reduce(&Integer::from(x - y), &n));
            }
        }
    }

    /// The GMP path of [`Mont`] against its own CIOS, including the carry out of the reduction.
    #[test]
    fn gmp_matches_cios() {
        fn check<const N: usize>(rand: &mut RandState<'_>) {
            let r = Integer::from(Integer::u_pow_u(2, 64 * N as u32));
            let mut moduli = vec![Integer::from(&r - 1), Integer::from(&r >> 1) + 1u32];
            moduli.extend((0..6).map(|_| adversarial(64 * N as u32, rand) | 1u32));
            for n in moduli
                .into_iter()
                .filter(|n| n.significant_bits() > 64 * N as u32 - 64)
            {
                let a = Mont::<N>::new(&n);
                let mut values: Vec<[u64; N]> = vec![a.residue(&Integer::from(&n - 1))];
                values.extend((0..6).map(|_| a.residue(&(adversarial(64 * N as u32, rand) % &n))));
                let (mut r1, mut r2) = ([0; N], [0; N]);
                for x in &values {
                    a.cios(&mut r1, x, x);
                    a.gmp_mul(&mut r2, x, None);
                    assert_eq!(r1, r2, "square mod {n}");
                    for y in &values {
                        a.cios(&mut r1, x, y);
                        a.gmp_mul(&mut r2, x, Some(y));
                        assert_eq!(r1, r2, "product mod {n}");
                    }
                }
            }
        }
        if !mpn::ENABLED {
            return;
        }
        let mut rand = RandState::new();
        check::<1>(&mut rand);
        check::<4>(&mut rand);
        check::<11>(&mut rand);
        check::<12>(&mut rand);
        check::<13>(&mut rand);
        check::<16>(&mut rand);
    }

    /// `1/2^64 mod n`, or `1` if `n` is even.
    fn inverse_r(n: &Integer) -> Integer {
        Integer::from(Integer::u_pow_u(2, 64))
            .invert(n)
            .unwrap_or_else(|_| Integer::from(1))
    }

    #[test]
    fn pow_matches_integers() {
        let mut rand = RandState::new();
        let n = Integer::from(Integer::u_pow_u(2, 300)) - 153u32;
        for bits in [1, 2, 5, 63, 64, 100, 511, 512, 3000, 9000, 70000] {
            for _ in 0..2 {
                let mut e = Integer::from(Integer::random_bits(bits, &mut rand));
                e.set_bit(bits - 1, true);
                let x = n.clone().random_below(&mut rand);
                let expected = x.clone().pow_mod(&e, &n).unwrap();
                with_arith!(&n, |a| assert_eq!(
                    a.to_integer(&pow(&a, &a.residue(&x), &e)),
                    expected
                ));
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
            assert_eq!(mont_limbs(&n) == 0, n.is_even() || !mpn::ENABLED);
            let a = Plain::new(&n);
            check(&a, &mut rand);
            check_small(&a, &Integer::from(12345), &mut rand);
        }
    }

    /// Random odd modulus of exactly `limbs` limbs.
    fn random_modulus(limbs: u32, rand: &mut RandState<'_>) -> Integer {
        let mut n = Integer::from(Integer::random_bits(64 * limbs, rand));
        n.set_bit(0, true);
        n.set_bit(64 * limbs - 1, true);
        n
    }

    /// Limb counts of the [`MontLarge`] tests: from the first one, around the thresholds of
    /// the reductions, and large.
    fn large_limbs() -> Vec<u32> {
        let (two, sub) = (REDC_2_LIMBS as u32, REDC_N_LIMBS as u32);
        let mut limbs = vec![MAX_LIMBS as u32 + 1, MAX_LIMBS as u32 + 2];
        limbs.extend([
            two - 1,
            two,
            two + 1,
            32,
            sub - 1,
            sub,
            sub + 1,
            64,
            100,
            128,
            256,
            300,
        ]);
        limbs
    }

    /// [`MontLarge`] against [`Plain`], on random operands and chains of operations, including
    /// the polynomial representation.
    #[test]
    fn mont_large_matches_plain() {
        if !mpn::ENABLED {
            return;
        }
        let mut rand = RandState::new();
        for limbs in large_limbs() {
            let r = Integer::from(1) << (64 * limbs);
            let mut moduli = vec![
                random_modulus(limbs, &mut rand),
                Integer::from(&r - 1),
                Integer::from(&r >> 1) + 1u32,
                (Integer::from(1) << (64 * limbs - 64)) + 1u32,
            ];
            moduli.push(adversarial(64 * limbs, &mut rand) | 1u32);
            // A known factor f, for the gcds.
            let f = Integer::from(u64::MAX - 58);
            moduli.push(&f * random_modulus(limbs - 1, &mut rand));
            for n in moduli {
                let size = n.significant_bits().div_ceil(64);
                assert_eq!(mont_limbs(&n), size as usize);
                if size as usize <= MAX_LIMBS {
                    continue;
                }
                let (a, p) = (MontLarge::new(&n), Plain::new(&n));
                let multiple = Integer::from(&f * 12345u32) % &n;
                compare_with_plain(&a, &p, multiple, &mut rand);
            }
        }
    }

    /// Every operation of `a` against `p`, on random values, the extremes, and `extra`.
    fn compare_with_plain(a: &MontLarge, p: &Plain, extra: Integer, rand: &mut RandState<'_>) {
        let n = a.modulus().clone();
        let mut values: Vec<Integer> = (0..6).map(|_| n.clone().random_below(rand)).collect();
        values.extend([
            Integer::ZERO,
            Integer::from(1),
            Integer::from(&n - 1),
            adversarial(n.significant_bits(), rand) % &n,
            extra,
        ]);
        let (mut r, mut rp) = (a.zero(), p.zero());
        for x in &values {
            let (ex, px) = (a.residue(x), p.residue(x));
            assert_eq!(a.to_integer(&ex), *x);
            assert_eq!(a.gcd(&ex), p.gcd(&px));
            a.sqr(&mut r, &ex);
            p.sqr(&mut rp, &px);
            assert_eq!(a.to_integer(&r), rp, "square mod {n}");
            for c in [1, 2, 3, u64::MAX, rand.bits(32).into()] {
                // The residue with the one-limb form c: c/2^64 mod n.
                let y = Integer::from(c) * inverse_r(&n) % &n;
                assert_eq!(a.small(&y), Some(c));
                a.mul_small(&mut r, &ex, c);
                p.mul(&mut rp, &px, &p.residue(&y));
                assert_eq!(a.to_integer(&r), rp, "{x}*{c}/2^64 mod {n}");
            }
            for y in &values {
                let (ey, py) = (a.residue(y), p.residue(y));
                a.mul(&mut r, &ex, &ey);
                p.mul(&mut rp, &px, &py);
                assert_eq!(a.to_integer(&r), rp, "product mod {n}");
                a.add(&mut r, &ex, &ey);
                p.add(&mut rp, &px, &py);
                assert_eq!(a.to_integer(&r), rp, "sum mod {n}");
                a.sub(&mut r, &ex, &ey);
                p.sub(&mut rp, &px, &py);
                assert_eq!(a.to_integer(&r), rp, "difference mod {n}");
                let (f, fp) = (a.factor(y), p.factor(y));
                a.mul_factor(&mut r, &ex, &f);
                p.mul_factor(&mut rp, &px, &fp);
                assert_eq!(a.to_integer(&r), rp, "factor mod {n}");
                // Polynomial representation.
                a.poly_mul(&mut r, &a.poly_from(x), &a.poly_from(y));
                assert_eq!(a.poly_value(&r), Integer::from(x * y) % &n);
            }
        }
        // A chain mixing every operation, from random values: the residues stay reduced.
        let (mut x, mut y) = (a.residue(&values[0]), a.residue(&values[1]));
        let (mut px, mut py) = (p.residue(&values[0]), p.residue(&values[1]));
        let f = a.factor(&Integer::from(12345));
        let fp = p.factor(&Integer::from(12345));
        for i in 0..50 {
            match i % 5 {
                0 => {
                    a.mul(&mut r, &x, &y);
                    p.mul(&mut rp, &px, &py);
                }
                1 => {
                    a.sqr(&mut r, &x);
                    p.sqr(&mut rp, &px);
                }
                2 => {
                    a.sub(&mut r, &y, &x);
                    p.sub(&mut rp, &py, &px);
                }
                3 => {
                    a.add(&mut r, &x, &y);
                    p.add(&mut rp, &px, &py);
                }
                _ => {
                    a.mul_factor(&mut r, &x, &f);
                    p.mul_factor(&mut rp, &px, &fp);
                }
            }
            std::mem::swap(&mut x, &mut y);
            std::mem::swap(&mut px, &mut py);
            std::mem::swap(&mut y, &mut r);
            std::mem::swap(&mut py, &mut rp);
            assert!(Integer::from_digits(&y, Order::Lsf) < n, "step {i} mod {n}");
            assert_eq!(a.to_integer(&y), py, "step {i} mod {n}");
        }
    }

    /// The three reductions of [`MontLarge`] agree, whatever the size (the subquadratic one
    /// needs more than 8 limbs).
    #[test]
    fn mont_large_reductions_agree() {
        if !mpn::ENABLED {
            return;
        }
        let mut rand = RandState::new();
        for limbs in [17, 20, 33, 47, 48, 64, 130] {
            let n = random_modulus(limbs, &mut rand);
            let mut arith = MontLarge::new(&n);
            let inv = n
                .clone()
                .invert(&(Integer::from(1) << (64 * limbs)))
                .unwrap();
            arith.inv = arith.limbs_of(&inv);
            let values: Vec<Vec<u64>> = (0..5)
                .map(|_| arith.residue(&n.clone().random_below(&mut rand)))
                .chain([arith.residue(&Integer::from(&n - 1))])
                .collect();
            for x in &values {
                for y in &values {
                    let results = [Redc::One, Redc::Two, Redc::Subquadratic].map(|redc| {
                        arith.redc = redc;
                        let mut r = arith.zero();
                        arith.mul(&mut r, x, y);
                        let mut s = arith.zero();
                        arith.sqr(&mut s, x);
                        (r, s)
                    });
                    assert_eq!(results[0], results[1], "{limbs} limbs");
                    assert_eq!(results[0], results[2], "{limbs} limbs");
                }
            }
        }
    }
}
