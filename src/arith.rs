//! Modular arithmetic for the curve computations.
//!
//! Residues modulo an odd `n` of at most [`MAX_LIMBS`] 64-bit limbs are kept in Montgomery
//! representation (`x` is stored as `x*R mod n` with `R = 2^(64*limbs)`), in fixed-size limb
//! arrays: a modular multiplication is one interleaved multiply-and-reduce (CIOS) pass (from
//! [`GMP_LIMBS`] limbs, GMP's product and Montgomery reduction), with no division and no
//! allocation. Larger (or even) moduli fall back to [`Plain`] big integer arithmetic.
//!
//! The representation never leaks: values go in with [`Arith::residue`] and out with
//! [`Arith::to_integer`], and [`Arith::gcd`] gives `gcd(x, n)` directly (`R` is coprime to an
//! odd `n`).

use rug::{Assign, Integer, integer::Order};

use crate::config::GMP_LIMBS;
pub use crate::config::MAX_LIMBS;

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
fn mac(lo: u64, a: u64, b: u64, carry: u64) -> (u64, u64) {
    let t = u128::from(lo) + u128::from(a) * u128::from(b) + u128::from(carry);
    (t as u64, (t >> 64) as u64)
}

/// `a + b + carry` (carry is 0 or 1), as `(sum, carry out)`.
#[inline(always)]
fn adc(a: u64, b: u64, carry: u64) -> (u64, u64) {
    let (s, c1) = a.overflowing_add(b);
    let (s, c2) = s.overflowing_add(carry);
    (s, u64::from(c1 | c2))
}

/// `a - b - borrow` (borrow is 0 or 1), as `(difference, borrow out)`.
#[inline(always)]
fn sbb(a: u64, b: u64, borrow: u64) -> (u64, u64) {
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
/// The polynomial representation of `x` is `x*R' mod n` (`R' = 2^(64*(N + 1))` for [`Mont`],
/// `R' = 1` for [`Plain`]), a value in `[0, n)`: [`PolyArith::redc_wide`] divides by `R'`, so it
/// maps the product of two such values back to the representation of the product.
pub trait PolyArith: Arith {
    /// Number of limbs of the values: `n < 2^(64*limbs)`.
    fn limbs(&self) -> usize;
    /// `r` = polynomial representation of the residue `x`.
    fn to_poly(&self, r: &mut Self::Elem, x: &Self::Elem);
    /// Writes the limbs of the value of `x` (polynomial representation) to `out`, of length
    /// [`PolyArith::limbs`].
    fn write_limbs(&self, x: &Self::Elem, out: &mut [u64]);
    /// `acc += x*y` (values of the polynomial representation), `acc` having at least
    /// `2*limbs + 1` limbs (the carry out of `acc` is lost).
    fn mul_acc(&self, acc: &mut [u64], x: &Self::Elem, y: &Self::Elem);
    /// `r = t/R' mod n`, for `t < n*R'` given by at most `2*limbs + 2` limbs.
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
        Integer::from_digits(&limbs, Order::Lsf)
    }
}

/// `acc += x*y` on limbs, `acc` having at least `x.len() + y.len()` limbs; returns the carry
/// out of `acc`.
fn mul_acc_limbs(acc: &mut [u64], x: &[u64], y: &[u64]) -> u64 {
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
/// Dispatched like the curve code (Montgomery up to 16 limbs, plain integers above): the setup
/// (conversion to the internal representation) is done by [`ArithBatch::new`], and
/// [`ArithBatch::run`] only computes.
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
        with_arith!(n, |arith| {
            let values = values.iter().map(|x| arith.residue(x)).collect();
            Self(Box::new(Batch { arith, values }))
        })
    }

    /// Number of limbs of the Montgomery implementation used, or `0` for plain integers.
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

    #[test]
    fn mont_edge_cases() {
        let mut rand = RandState::new();
        for limbs in 1..=MAX_LIMBS as u32 + 1 {
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
                let expected = if n.is_odd() && size <= MAX_LIMBS {
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
