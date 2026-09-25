//! Special reduction modulo `M = 2^k +- 1` for the divisors `n` of such numbers (Mersenne,
//! Fermat and Cunningham numbers and their cofactors), GMP-ECM's "special division" (`-base2`).
//!
//! The arithmetic ([`Base2`]) works modulo `M` instead of `n`: `n` divides `M`, so the residues
//! modulo `M` map to the residues modulo `n`, and every result reduced modulo `n` (the gcds, the
//! values out) is the one the arithmetic modulo `n` would give. A product `x = hi*2^k + lo` is
//! reduced to `lo + hi` (`2^k = 1` modulo `2^k - 1`) or `lo - hi` (`2^k = -1` modulo `2^k + 1`),
//! a shift and an addition: no Montgomery reduction, whose cost is that of a second product.
//!
//! [`Base2Mode`] chooses when: [`Base2Form::detect`] recognizes the numbers that divide
//! `2^k +- 1` for a `k` not much larger than their size, with GMP-ECM's rule.

use std::{cell::RefCell, fmt};

use rug::{Integer, integer::Order};

use crate::{
    arith::{Arith, PolyArith, adc, mac, mpn, mul_acc_limbs, reduce, sbb},
    config::{BASE2_MIN_EXPONENT, BASE2_THRESHOLD},
    cost::base2_faster,
};

/// When to use the special reduction modulo `2^k +- 1` (see [`crate::Factorizer::base2`]).
#[non_exhaustive]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub enum Base2Mode {
    /// For the numbers that divide `2^k +- 1` with `k` at most 1.4 times their size in bits
    /// (as GMP-ECM), when it is faster: from about 400 bits with `k` close to their size, from
    /// 770 bits with `k` up to 1.4 times their size (stage 2 computes modulo the number when
    /// that is cheaper, for `k` well above its size).
    #[default]
    Auto,
    /// Never (GMP-ECM's `-nobase2`).
    Off,
    /// Modulo `2^k + 1` if `k > 0`, `2^-k - 1` if `k < 0` (GMP-ECM's `-base2 k`), which the
    /// number must divide ([`crate::Error::InvalidOption`] otherwise), whatever its size.
    Force(i64),
}

impl Base2Mode {
    /// The form of `n` to reduce modulo, if any.
    ///
    /// # Errors
    ///
    /// With [`Base2Mode::Force`], if `n` does not divide the number.
    pub(crate) fn form(self, n: &Integer) -> Result<Option<Base2Form>, &'static str> {
        match self {
            Self::Auto => Ok(Base2Form::detect(n)),
            Self::Off => Ok(None),
            Self::Force(k) => Base2Form::from_signed(k)
                .filter(|form| form.divides(n))
                .map(Some)
                .ok_or(BASE2_ERROR),
        }
    }
}

/// The error of a [`Base2Mode::Force`] that does not apply to the number.
pub(crate) const BASE2_ERROR: &str = "base2: the number does not divide 2^k+1 (2^-k-1 if k < 0)";

/// `2^k + 1` (`plus`) or `2^k - 1`, a multiple of the number to factor.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct Base2Form {
    pub k: u32,
    pub plus: bool,
}

impl Base2Form {
    /// The form of `k` with GMP-ECM's convention: `2^k + 1` if `k > 0`, `2^-k - 1` if `k < 0`.
    pub fn from_signed(k: i64) -> Option<Self> {
        let plus = k > 0;
        let k = u32::try_from(k.unsigned_abs()).ok().filter(|&k| k > 0)?;
        Some(Self { k, plus })
    }

    /// GMP-ECM's convention: `k` for `2^k + 1`, `-k` for `2^k - 1`.
    pub fn signed(self) -> i64 {
        if self.plus {
            i64::from(self.k)
        } else {
            -i64::from(self.k)
        }
    }

    /// `2^k +- 1`.
    pub fn value(self) -> Integer {
        let power = Integer::from(1) << self.k;
        if self.plus {
            power + 1u32
        } else {
            power - 1u32
        }
    }

    /// Whether `n > 1` divides `2^k +- 1`.
    pub fn divides(self, n: &Integer) -> bool {
        if *n <= 1 || n.is_even() {
            return false;
        }
        let power = Integer::from(2)
            .pow_mod(&Integer::from(self.k), n)
            .expect("positive exponent");
        if self.plus {
            power == Integer::from(n - 1u32)
        } else {
            power == 1
        }
    }

    /// The form `2^k +- 1` of `n` worth using: GMP-ECM's `isbase2` (`k` at least
    /// [`BASE2_MIN_EXPONENT`] and at most [`BASE2_THRESHOLD`] times the size of `n`), if the
    /// arithmetic modulo `2^k +- 1` is faster than modulo `n` (from about 400 bits with `k`
    /// close to the size of `n`, from 770 bits with `k` up to 1.4 times, see
    /// [`crate::cost::base2_faster`]).
    pub fn detect(n: &Integer) -> Option<Self> {
        let bits = n.significant_bits();
        let form = Self::find(n)?;
        let lo = f64::from(bits - 1);
        (form.k >= BASE2_MIN_EXPONENT
            && f64::from(form.k) <= BASE2_THRESHOLD * lo
            && base2_faster(bits, form.k))
        .then_some(form)
    }

    /// The `2^k +- 1` with `lo < k <= 2*lo` that `n > 2` divides, if any (then unique), or `n`
    /// itself if it is `2^lo + 1`, with `2^lo <= n < 2^(lo + 1)`: `2^(2*lo) mod n` is
    /// `2^(2*lo - k)` if `n` divides `2^k - 1`, `n - 2^(2*lo - k)` if `n` divides `2^k + 1` (one
    /// test, as GMP-ECM's `isbase2`).
    fn find(n: &Integer) -> Option<Self> {
        if *n <= 2 || n.is_even() {
            return None;
        }
        let lo = n.significant_bits() - 1;
        let w = Integer::from(Integer::u_pow_u(2, 2 * lo)) % n;
        let power_of_two = |x: &Integer| x.is_power_of_two().then(|| x.significant_bits() - 1);
        if w == 1 {
            // n divides 2^(2*lo) - 1 = (2^lo - 1)(2^lo + 1), and is above 2^lo - 1: 2^lo + 1.
            let form = Self { k: lo, plus: true };
            return (*n == form.value()).then_some(form);
        }
        if let Some(j) = power_of_two(&w) {
            return Some(Self {
                k: 2 * lo - j,
                plus: false,
            });
        }
        let j = power_of_two(&Integer::from(n - &w))?;
        Some(Self {
            k: 2 * lo - j,
            plus: true,
        })
    }
}

impl fmt::Display for Base2Form {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "2^{}{}1", self.k, if self.plus { '+' } else { '-' })
    }
}

/// Arithmetic modulo `n` computed modulo `M = 2^k +- 1`, a multiple of `n`: residues are values
/// in `[0, M]` for `2^k - 1`, in `[0, M)` for `2^k + 1`, on `limbs` limbs.
///
/// A product (GMP's `mpn_mul_n` or `mpn_sqr`, on the limbs below `2^k`) is reduced by one
/// addition or subtraction of its halves (see the module documentation), then at most one of
/// `M`: about the cost of the product alone.
#[derive(Debug)]
pub struct Base2 {
    n: Integer,
    k: usize,
    plus: bool,
    /// Limbs of the residues: `M` and the sums of two residues fit (with a carry).
    limbs: usize,
    /// Limbs of the factors of a product: `ceil(k/64)`. One less than `limbs` for `2^k + 1`
    /// when `64` divides `k`: the top limb is then `0`, but for the residue `2^k = -1`.
    mul_limbs: usize,
    /// The residue `0`, for negations.
    zero: Vec<u64>,
    /// Buffers: products, the chunks of the wide reductions, and the products of
    /// [`PolyArith::mul_acc`].
    scratch: RefCell<Scratch>,
}

impl Base2 {
    /// Arithmetic modulo `n`, a divisor of `2^k +- 1` (the form).
    pub(crate) fn new(n: &Integer, form: Base2Form) -> Self {
        debug_assert!(form.divides(n));
        let k = form.k as usize;
        let limbs = (k + usize::from(form.plus)).div_ceil(64);
        Self {
            n: n.clone(),
            k,
            plus: form.plus,
            limbs,
            mul_limbs: k.div_ceil(64),
            zero: vec![0; limbs],
            scratch: RefCell::new(Scratch {
                // Room for a product and a few more (zero) limbs to read.
                product: vec![0; 2 * limbs + 2],
                chunk: vec![0; limbs],
                wide: vec![0; 2 * limbs],
            }),
        }
    }

    /// `r = a + b`.
    #[inline(always)]
    fn add_to(&self, r: &mut [u64], a: &[u64], b: &[u64]) {
        let mut carry = 0;
        for ((r, &a), &b) in r.iter_mut().zip(a).zip(b) {
            (*r, carry) = adc(a, b, carry);
        }
        self.wrap_sum(r, carry);
    }

    /// `r = r mod M` in the range of the residues, for `r + carry*2^(64*limbs)` the sum of two
    /// residues.
    #[inline(always)]
    fn wrap_sum(&self, r: &mut [u64], carry: u64) {
        if !self.plus {
            return self.wrap_up(r, carry);
        }
        // r = h*2^k + l with h <= 2 (the sum is at most 2^(k + 1)), and 2^k = -1: l - h.
        let (w, s) = (self.k / 64, self.k % 64);
        let h = if s == 0 {
            r[w]
        } else {
            // A carry only if s = 63.
            (r[w] >> s) | (carry << (63 - s) << 1)
        };
        r[w] &= (1 << s) - 1;
        let mut borrow;
        (r[0], borrow) = sbb(r[0], h, 0);
        let mut j = 1;
        while borrow != 0 && j < r.len() {
            (r[j], borrow) = sbb(r[j], 0, borrow);
            j += 1;
        }
        self.wrap_down(r, borrow);
    }

    /// `r = a - b`.
    #[inline(always)]
    fn sub_to(&self, r: &mut [u64], a: &[u64], b: &[u64]) {
        let mut borrow = 0;
        for ((r, &a), &b) in r.iter_mut().zip(a).zip(b) {
            (*r, borrow) = sbb(a, b, borrow);
        }
        self.wrap_down(r, borrow);
    }

    /// `r += b`.
    fn add_assign(&self, r: &mut [u64], b: &[u64]) {
        let mut carry = 0;
        for (r, &b) in r.iter_mut().zip(b) {
            (*r, carry) = adc(*r, b, carry);
        }
        self.wrap_sum(r, carry);
    }

    /// `r -= b`.
    fn sub_assign(&self, r: &mut [u64], b: &[u64]) {
        let mut borrow = 0;
        for (r, &b) in r.iter_mut().zip(b) {
            (*r, borrow) = sbb(*r, b, borrow);
        }
        self.wrap_down(r, borrow);
    }

    /// `p = a*b` (`a^2` if `b` is `None`) on `mul_limbs` limbs, `p` of `2*mul_limbs` limbs.
    #[inline(always)]
    fn product(&self, p: &mut [u64], a: &[u64], b: Option<&[u64]>) {
        let ml = self.mul_limbs;
        let a = &a[..ml];
        if mpn::ENABLED {
            match b {
                Some(b) => mpn::mul(p, a, &b[..ml]),
                None => mpn::sqr(p, a),
            }
        } else {
            p.fill(0);
            mul_acc_limbs(p, a, &b.unwrap_or(a)[..ml]);
        }
    }

    /// `r = p mod M`, for a product `p <= 2^(2*k)` of two residues (its limbs from `2*mul_limbs`
    /// on are zero): `lo +- hi`, with `p = hi*2^k + lo`.
    #[inline(always)]
    fn fold(&self, r: &mut [u64], p: &[u64]) {
        let (w, s) = (self.k / 64, self.k % 64);
        let limbs = self.limbs;
        // lo: the limbs below w, then the low s bits of limb w (none if s = 0).
        let top = p[w] & ((1u64 << s) - 1);
        if s == 0 {
            let hi = &p[w..w + limbs];
            if self.plus {
                let mut borrow = 0;
                for j in 0..w {
                    (r[j], borrow) = sbb(p[j], hi[j], borrow);
                }
                (r[w], borrow) = sbb(0, hi[w], borrow);
                self.wrap_down(r, borrow);
            } else {
                let mut carry = 0;
                for j in 0..w {
                    (r[j], carry) = adc(p[j], hi[j], carry);
                }
                self.wrap_up(r, carry);
            }
        } else {
            let hi = |j: usize| (p[w + j] >> s) | (p[w + j + 1] << (64 - s));
            // limbs = w + 1.
            if self.plus {
                let mut borrow = 0;
                for j in 0..w {
                    (r[j], borrow) = sbb(p[j], hi(j), borrow);
                }
                (r[w], borrow) = sbb(top, hi(w), borrow);
                self.wrap_down(r, borrow);
            } else {
                let mut carry = 0;
                for j in 0..w {
                    (r[j], carry) = adc(p[j], hi(j), carry);
                }
                (r[w], carry) = adc(top, hi(w), carry);
                self.wrap_up(r, carry);
            }
        }
    }

    /// For `2^k - 1`: `r = r mod M` in `[0, M]`, for `r + carry*2^(64*limbs) <= 2*M` (a sum of
    /// two residues): `2^k = 1`, so the bit `k` (`carry` if `64` divides `k`) moves to the bottom.
    #[inline(always)]
    fn wrap_up(&self, r: &mut [u64], carry: u64) {
        let (w, s) = (self.k / 64, self.k % 64);
        let bit = if s == 0 {
            carry
        } else {
            let bit = r[w] >> s;
            r[w] &= (1 << s) - 1;
            bit
        };
        increment(r, bit);
    }

    /// `r = r mod M` in the range of the residues, for `r` the difference `d` of a residue and
    /// a value at most `2^k` (for `2^k + 1`) or `M` (for `2^k - 1`), `r = d + 2^(64*limbs)` if
    /// `borrow`: `d + M` is then `d + 2^k +- 1`, the `k` low bits of `r`, plus or minus `1`.
    #[inline(always)]
    fn wrap_down(&self, r: &mut [u64], borrow: u64) {
        let (w, s) = (self.k / 64, self.k % 64);
        if w < self.limbs {
            // Keep the k low bits if borrow.
            r[w] &= ((1 << s) - 1) | borrow.wrapping_sub(1);
        }
        if self.plus {
            increment(r, borrow);
        } else {
            decrement(r, borrow);
        }
    }

    /// `r = a*b` (`a^2` if `b` is `None`).
    #[inline(always)]
    fn mul_or_sqr(&self, r: &mut [u64], a: &[u64], b: Option<&[u64]>) {
        let ml = self.mul_limbs;
        if ml < self.limbs {
            // 2^k + 1 with 64 | k: the residue 2^k is -1.
            if a[ml] != 0 {
                return self.sub_to(r, &self.zero, b.unwrap_or(a));
            }
            if b.is_some_and(|b| b[ml] != 0) {
                return self.sub_to(r, &self.zero, a);
            }
        }
        let p = &mut self.scratch.borrow_mut().product;
        self.product(&mut p[..2 * ml], a, b);
        self.fold(r, p);
    }

    /// `r = t mod M`, for any `t`: the sum of its chunks of `k` bits, with alternating signs
    /// for `2^k + 1` (each chunk is below `2^k`, a residue).
    fn reduce_wide(&self, r: &mut [u64], t: &[u64]) {
        let len = t.iter().rposition(|&x| x != 0).map_or(0, |i| i + 1);
        let bits = 64 * len;
        let (w, s) = (self.k / 64, self.k % 64);
        let chunk = &mut self.scratch.borrow_mut().chunk;
        r.fill(0);
        let limb = |i: usize| t.get(i).copied().unwrap_or(0);
        for (i, start) in (0..bits).step_by(self.k).enumerate() {
            let (word, shift) = (start / 64, start % 64);
            for (j, c) in chunk.iter_mut().enumerate() {
                *c = if shift == 0 {
                    limb(word + j)
                } else {
                    (limb(word + j) >> shift) | (limb(word + j + 1) << (64 - shift))
                };
            }
            // Keep the k low bits.
            if w < self.limbs {
                chunk[w] &= (1u64 << s) - 1;
                chunk[w + 1..].fill(0);
            }
            if self.plus && i % 2 == 1 {
                self.sub_assign(r, chunk);
            } else {
                self.add_assign(r, chunk);
            }
        }
    }
}

/// Buffers of [`Base2`].
#[derive(Debug)]
struct Scratch {
    product: Vec<u64>,
    chunk: Vec<u64>,
    wide: Vec<u64>,
}

/// `r += c` (`c` is 0 or 1; the carry out of `r` is lost).
#[inline(always)]
fn increment(r: &mut [u64], c: u64) {
    let mut carry;
    (r[0], carry) = adc(r[0], c, 0);
    let mut j = 1;
    while carry != 0 && j < r.len() {
        (r[j], carry) = adc(r[j], 0, carry);
        j += 1;
    }
}

/// `r -= c` (`c` is 0 or 1, at most `r`).
#[inline(always)]
fn decrement(r: &mut [u64], c: u64) {
    let mut borrow;
    (r[0], borrow) = sbb(r[0], c, 0);
    let mut j = 1;
    while borrow != 0 && j < r.len() {
        (r[j], borrow) = sbb(r[j], 0, borrow);
        j += 1;
    }
}

impl Arith for Base2 {
    type Elem = Vec<u64>;

    fn modulus(&self) -> &Integer {
        &self.n
    }

    fn zero(&self) -> Vec<u64> {
        self.zero.clone()
    }

    fn residue(&self, x: &Integer) -> Vec<u64> {
        // Below n, a divisor of M: a residue.
        let mut r = self.zero();
        reduce(x, &self.n).write_digits(&mut r, Order::Lsf);
        r
    }

    fn to_integer(&self, x: &Vec<u64>) -> Integer {
        Integer::from_digits(x, Order::Lsf) % &self.n
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
        self.add_to(r, a, b);
    }

    #[inline(always)]
    fn sub(&self, r: &mut Vec<u64>, a: &Vec<u64>, b: &Vec<u64>) {
        self.sub_to(r, a, b);
    }

    /// As [`crate::arith::Plain`]: the one-limb form of `x` is its value modulo `n`.
    fn small(&self, x: &Integer) -> Option<u64> {
        reduce(x, &self.n).to_u64()
    }

    fn mul_small(&self, r: &mut Vec<u64>, a: &Vec<u64>, c: u64) {
        let limbs = self.limbs;
        if self.mul_limbs == 1 {
            // k <= 64: a*c on at most 3 limbs.
            let mut t = [0; 3];
            let mut carry = 0;
            for (t, &a) in t.iter_mut().zip(a) {
                (*t, carry) = mac(0, a, c, carry);
            }
            t[limbs] = carry;
            return self.reduce_wide(r, &t);
        }
        // a*c < 2^(k + 65) <= 2^(2*k): hi < 2^64, folded as a product.
        let p = &mut self.scratch.borrow_mut().product;
        let mut carry = 0;
        for (p, &a) in p.iter_mut().zip(a) {
            (*p, carry) = mac(0, a, c, carry);
        }
        p[limbs] = carry;
        p[limbs + 1..2 * self.mul_limbs].fill(0);
        self.fold(r, p);
    }
}

/// The polynomial representation is the residue itself (`R' = 1`).
impl PolyArith for Base2 {
    fn limbs(&self) -> usize {
        self.limbs
    }

    fn value_bits(&self) -> usize {
        // M, or 2^k for 2^k + 1.
        self.k + usize::from(self.plus)
    }

    fn to_poly(&self, r: &mut Vec<u64>, x: &Vec<u64>) {
        r.clone_from(x);
    }

    fn write_limbs(&self, x: &Vec<u64>, out: &mut [u64]) {
        out.copy_from_slice(x);
    }

    fn mul_acc(&self, acc: &mut [u64], x: &Vec<u64>, y: &Vec<u64>) {
        if !mpn::ENABLED {
            mul_acc_limbs(acc, x, y);
            return;
        }
        let p = &mut self.scratch.borrow_mut().wide;
        mpn::mul(p, x, y);
        let mut carry = 0;
        for (a, &p) in acc.iter_mut().zip(p.iter()) {
            (*a, carry) = adc(*a, p, carry);
        }
        increment(&mut acc[p.len()..], carry);
    }

    fn redc_wide(&self, r: &mut Vec<u64>, t: &[u64]) {
        self.reduce_wide(r, t);
    }
}

/// Tests: `2^k + 1` (`k > 0`) or `2^-k - 1` (`k < 0`) divided by its prime factors below 2^16,
/// and its form.
#[cfg(test)]
pub(crate) fn cofactor_of(k: i64) -> (Integer, Base2Form) {
    let form = Base2Form::from_signed(k).unwrap();
    let mut n = form.value();
    for p in primal::Primes::all().take_while(|&p| p < 1 << 16) {
        while n.is_divisible_u(p as u32) && n > p {
            n /= p as u32;
        }
    }
    (n, form)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rug::rand::RandState;

    fn cofactor(k: i64) -> Integer {
        cofactor_of(k).0
    }

    /// Numbers `2^k +- 1` and their cofactors, for `k` around the word boundaries and some
    /// larger ones.
    fn moduli() -> Vec<(Integer, Base2Form)> {
        // Numbers much smaller than their 2^k +- 1 too.
        let mut out = vec![
            (
                Base2Form::from_signed(64).unwrap().value(),
                Base2Form {
                    k: 1088,
                    plus: true,
                },
            ),
            (
                Base2Form::from_signed(-61).unwrap().value(),
                Base2Form {
                    k: 610,
                    plus: false,
                },
            ),
        ];
        for k in [
            17, 31, 61, 63, 64, 65, 67, 101, 127, 128, 129, 191, 192, 193, 255, 256, 257, 511, 512,
            513, 521, 607, 1023, 1024, 1025, 1061, 1279, 2048, 2203,
        ] {
            for plus in [false, true] {
                let form = Base2Form { k, plus };
                out.push((form.value(), form));
                let n = cofactor(form.signed());
                if n != form.value() && n > 1u32 << 16 {
                    out.push((n, form));
                }
            }
        }
        out
    }

    /// A residue is in `[0, M]` (`2^k - 1`) or `[0, M)` (`2^k + 1`).
    fn assert_in_range(r: &[u64], m: &Integer, form: Base2Form) {
        let v = Integer::from_digits(r, Order::Lsf);
        assert!(
            if form.plus { v < *m } else { v <= *m },
            "{v} out of range for {form}"
        );
    }

    /// Checks every operation of [`Base2`] modulo `n` against integers, on random residues and
    /// on the extremes of the representation.
    #[test]
    fn matches_integers() {
        let mut rand = RandState::new();
        for (n, form) in moduli() {
            let a = Base2::new(&n, form);
            let m = form.value();
            let mut values: Vec<Integer> =
                (0..10).map(|_| n.clone().random_below(&mut rand)).collect();
            values.extend([Integer::ZERO, Integer::from(1), Integer::from(&n - 1)]);
            let mut all: Vec<(Integer, Vec<u64>)> = values
                .into_iter()
                .map(|x| (x.clone(), a.residue(&x)))
                .collect();
            // The largest residues, M (0) for 2^k - 1 and 2^k (-1) for 2^k + 1, and the one below:
            // products and sums reach them.
            let top = if form.plus {
                Integer::from(&m - 1)
            } else {
                m.clone()
            };
            for x in [Integer::from(&top - 1), top] {
                let mut e = a.zero();
                x.write_digits(&mut e, Order::Lsf);
                all.push((x % &n, e));
            }
            let mut r = a.zero();
            for (x, ex) in &all {
                assert_eq!(a.to_integer(ex), *x);
                assert_eq!(a.gcd(ex), x.clone().gcd(&n));
                a.sqr(&mut r, ex);
                assert_in_range(&r, &m, form);
                assert_eq!(a.to_integer(&r), Integer::from(x * x) % &n, "{x}^2 mod {n}");
                for c in [0, 1, 2, 3, u64::MAX] {
                    a.mul_small(&mut r, ex, c);
                    assert_in_range(&r, &m, form);
                    assert_eq!(a.to_integer(&r), Integer::from(x * c) % &n);
                }
                for (y, ey) in &all {
                    a.mul(&mut r, ex, ey);
                    assert_in_range(&r, &m, form);
                    assert_eq!(
                        a.to_integer(&r),
                        Integer::from(x * y) % &n,
                        "{x}*{y} mod {n}"
                    );
                    a.add(&mut r, ex, ey);
                    assert_in_range(&r, &m, form);
                    assert_eq!(a.to_integer(&r), Integer::from(x + y) % &n);
                    a.sub(&mut r, ex, ey);
                    assert_in_range(&r, &m, form);
                    assert_eq!(a.to_integer(&r), reduce(&Integer::from(x - y), &n));
                    let f = a.factor(y);
                    a.mul_factor(&mut r, ex, &f);
                    assert_eq!(a.to_integer(&r), Integer::from(x * y) % &n);
                }
            }
        }
    }

    /// Long chains of operations stay in range and match the integers.
    #[test]
    fn chains() {
        let mut rand = RandState::new();
        for (n, form) in moduli() {
            let a = Base2::new(&n, form);
            let m = form.value();
            let x = n.clone().random_below(&mut rand);
            let (mut e, mut t) = (a.residue(&x), a.zero());
            let mut v = x;
            for i in 0..300u32 {
                if i % 3 == 0 {
                    a.sqr(&mut t, &e);
                    v.square_mut();
                } else {
                    let y = a.residue(&Integer::from(i));
                    a.mul(&mut t, &e, &y);
                    std::mem::swap(&mut e, &mut t);
                    a.sub(&mut t, &e, &y);
                    v = v * i - i;
                }
                v = reduce(&v, &n);
                std::mem::swap(&mut e, &mut t);
                assert_in_range(&e, &m, form);
            }
            assert_eq!(a.to_integer(&e), v);
        }
    }

    #[test]
    fn wide_reduction() {
        let mut rand = RandState::new();
        for (n, form) in moduli() {
            let a = Base2::new(&n, form);
            let m = form.value();
            for limbs in [1, a.limbs, 2 * a.limbs + 2] {
                let t = Integer::from(Integer::random_bits(64 * limbs as u32, &mut rand));
                let mut digits = vec![0; limbs];
                t.write_digits(&mut digits, Order::Lsf);
                let mut r = a.zero();
                a.redc_wide(&mut r, &digits);
                assert_in_range(&r, &m, form);
                assert_eq!(
                    a.to_integer(&r),
                    Integer::from(&t % &n),
                    "{t} mod {n} ({form})"
                );
            }
            let (x, y) = (
                n.clone().random_below(&mut rand),
                n.clone().random_below(&mut rand),
            );
            let mut r = a.zero();
            a.poly_mul(&mut r, &a.poly_from(&x), &a.poly_from(&y));
            assert_eq!(a.poly_value(&r), Integer::from(&x * &y) % &n);
        }
    }

    #[test]
    fn detection() {
        let form = |k: i64| Base2Form::from_signed(k).unwrap();
        // M1061 (prime), Fermat numbers and cofactors.
        for k in [-1061, 1024, 2048, -1279, 1000, -2000] {
            assert_eq!(Base2Form::detect(&form(k).value()), Some(form(k)));
            let n = cofactor(k);
            assert_eq!(Base2Form::detect(&n), Some(form(k)), "{k}");
            assert!(form(k).divides(&n) && form(-2 * k.abs()).divides(&n));
        }
        // The primitive part of 2^750 + 1 (500 bits): k above 1.4 times its size.
        let n = form(750).value() / form(250).value();
        assert_eq!(Base2Form::find(&n), Some(form(750)));
        assert_eq!(Base2Form::detect(&n), None);
        assert_eq!(Base2Mode::Force(750).form(&n), Ok(Some(form(750))));
        assert_eq!(Base2Mode::Force(1500).form(&n), Err(BASE2_ERROR));
        assert_eq!(Base2Mode::Force(-1500).form(&n), Ok(Some(form(-1500))));
        assert_eq!(Base2Mode::Force(0).form(&n), Err(BASE2_ERROR));
        assert_eq!(Base2Mode::Off.form(&n), Ok(None));
        // Below BASE2_MIN_BITS.
        assert_eq!(Base2Form::detect(&form(-127).value()), None);
        assert_eq!(Base2Form::find(&form(-127).value()), Some(form(-127)));
        assert_eq!(Base2Form::find(&form(64).value()), Some(form(64)));
        assert_eq!(Base2Form::find(&Integer::from(3)), Some(form(1)));
        assert_eq!(Base2Form::find(&Integer::from(5)), Some(form(2)));
        // Other numbers.
        let mut rand = RandState::new();
        for bits in [600, 1000, 2000] {
            let n = Integer::from(Integer::random_bits(bits, &mut rand)) | 1u32;
            assert_eq!(Base2Form::detect(&n), None);
        }
        assert_eq!(Base2Form::detect(&(form(-1061).value() + 4u32)), None);
        assert_eq!(form(-5).to_string(), "2^5-1");
        assert_eq!(form(7).to_string(), "2^7+1");
        assert_eq!((form(7).signed(), form(-7).signed()), (7, -7));
        assert_eq!(Base2Form::from_signed(1 << 40), None);
    }

    /// Polynomial products (schoolbook, Kronecker, middle and wrap-around) on residues with the
    /// largest values of the representation: the slots must hold their sums of products.
    #[test]
    fn polynomial_products() {
        use crate::poly::{Operand, Workspace, mul, mul_part, mul_wrap};
        let mut rand = RandState::new();
        let naive = |x: &[Integer], y: &[Integer], n: &Integer| {
            let mut r = vec![Integer::new(); x.len() + y.len() - 1];
            for (i, x) in x.iter().enumerate() {
                for (j, y) in y.iter().enumerate() {
                    r[i + j] += Integer::from(x * y);
                }
            }
            r.into_iter().map(|r| r % n).collect::<Vec<_>>()
        };
        // With n much smaller than 2^k +- 1 too: the values are those of the representation.
        let small = |k: i64, of: i64| {
            (
                Base2Form::from_signed(k).unwrap().value(),
                cofactor_of(of).1,
            )
        };
        let moduli = [-127, 128, 192, -521, 1024, -1061].map(cofactor_of);
        for (n, form) in moduli
            .into_iter()
            .chain([small(64, 1088), small(-127, -1016)])
        {
            let a = Base2::new(&n, form);
            let top = if form.plus {
                form.value() - 1u32
            } else {
                form.value()
            };
            let mut largest = a.zero();
            top.write_digits(&mut largest, Order::Lsf);
            let mut ws = Workspace::new();
            for (lx, ly) in [
                (1, 1),
                (3, 4),
                (12, 13),
                (13, 13),
                (20, 7),
                (40, 41),
                (100, 3),
                (64, 64),
                (300, 250),
            ] {
                let poly = |len: usize, rand: &mut RandState<'_>| -> Vec<Vec<u64>> {
                    (0..len)
                        .map(|i| {
                            if i % 3 == 0 {
                                largest.clone()
                            } else {
                                a.residue(&n.clone().random_below(rand))
                            }
                        })
                        .collect()
                };
                let (x, y) = (poly(lx, &mut rand), poly(ly, &mut rand));
                let values = |x: &[Vec<u64>]| x.iter().map(|x| a.to_integer(x)).collect::<Vec<_>>();
                let expected = naive(&values(&x), &values(&y), &n);
                let mut r = vec![a.zero(); lx + ly - 1];
                mul(&a, &mut ws, &mut r, &x, &y);
                assert_eq!(values(&r), expected, "{form} {lx} {ly}");
                // A middle product.
                let from = (lx + ly) / 3;
                let mut part = vec![a.zero(); lx + ly - 1 - from];
                mul_part(
                    &a,
                    &mut ws,
                    &mut part,
                    from,
                    Operand::new(&x),
                    Operand::new(&y),
                );
                assert_eq!(
                    values(&part),
                    expected[from..],
                    "{form} {lx} {ly} from {from}"
                );
                // Wrap-around: coefficient t is c_t + c_(t + L).
                let lmin = lx.max(ly) + 5;
                let mut wrap = vec![a.zero(); lmin];
                let l = mul_wrap(
                    &a,
                    &mut ws,
                    &mut wrap,
                    0,
                    Operand::new(&x),
                    Operand::new(&y),
                    lmin,
                );
                for (t, w) in values(&wrap).iter().enumerate().take(l.min(lmin)) {
                    let mut c = expected.get(t).cloned().unwrap_or_default();
                    c += expected.get(t + l).cloned().unwrap_or_default();
                    assert_eq!(*w, c % &n, "{form} {lx} {ly} wrap {t}");
                }
            }
        }
    }

    /// Values of the representation that make carries and borrows ripple through every limb
    /// (all-ones limbs, `2^k - 2^j`) and the extremes, for [`adversarial_values`].
    fn edge_values(form: Base2Form, rand: &mut RandState<'_>) -> Vec<Integer> {
        let k = form.k;
        let m = form.value();
        let top = if form.plus { Integer::from(&m - 1) } else { m };
        let pow = |j: u32| Integer::from(1) << j;
        let mut v = vec![
            Integer::ZERO,
            Integer::from(1),
            Integer::from(2),
            Integer::from(&top - 1),
            Integer::from(&top - 2),
            pow(k - 1),
            pow(k - 1) - 1u32,
            pow(k - 1) + 1u32,
            top.clone(),
        ];
        for j in [64, 128, 64 * (k / 64), k - 1, k.saturating_sub(64)] {
            if j > 0 && j <= k {
                v.push(pow(j) - 1u32);
            }
        }
        for j in [0, 1, 63, 64, 65, k / 2] {
            if j < k {
                v.push(pow(k) - pow(j));
            }
        }
        v.extend((0..4).map(|_| Integer::from(top.random_below_ref(rand)) + 1u32));
        v.retain(|x| *x <= top);
        v
    }

    /// Every operation on values anywhere in the representation (not only the residues of
    /// numbers below `n`), checked modulo `M` itself, for `k` from 2 (below
    /// [`BASE2_MIN_EXPONENT`], as `Force` allows) to 4097, `n` = `M`, a cofactor or `3`.
    #[test]
    fn adversarial_values() {
        let mut rand = RandState::new();
        for k in [
            2, 3, 5, 8, 15, 16, 17, 31, 32, 33, 63, 64, 65, 127, 128, 129, 191, 192, 1024, 1025,
            2047, 4096, 4097,
        ] {
            for plus in [false, true] {
                let form = Base2Form { k, plus };
                let m = form.value();
                let mut moduli = vec![m.clone(), cofactor(form.signed())];
                // 3 divides 2^k - 1 for k even, 2^k + 1 for k odd.
                if (k % 2 == 0) != plus {
                    moduli.push(Integer::from(3));
                }
                for n in moduli.iter().filter(|n| **n > 1) {
                    check_edges(form, n, &mut rand);
                }
            }
        }
    }

    fn check_edges(form: Base2Form, n: &Integer, rand: &mut RandState<'_>) {
        let m = form.value();
        let a = Base2::new(n, form);
        let values = edge_values(form, rand);
        let elem = |x: &Integer| {
            let mut e = a.zero();
            x.write_digits(&mut e, Order::Lsf);
            e
        };
        // The value modulo M of a result, checked in range.
        let value = |r: &[u64]| {
            assert_in_range(r, &m, form);
            Integer::from_digits(r, Order::Lsf) % &m
        };
        let mut r = a.zero();
        let mut acc = vec![0; 2 * a.limbs + 2];
        for x in &values {
            let ex = elem(x);
            assert_eq!(a.to_integer(&ex), Integer::from(x % n), "{x} mod {n}");
            assert_eq!(a.gcd(&ex), Integer::from(x % n).gcd(n));
            a.sqr(&mut r, &ex);
            assert_eq!(value(&r), Integer::from(x * x) % &m, "{x}^2 ({form})");
            for c in [0, 1, 3, 1 << 63, u64::MAX] {
                a.mul_small(&mut r, &ex, c);
                assert_eq!(value(&r), Integer::from(x * c) % &m, "{x}*{c} ({form})");
            }
            for y in &values {
                let ey = elem(y);
                a.mul(&mut r, &ex, &ey);
                assert_eq!(value(&r), Integer::from(x * y) % &m, "{x}*{y} ({form})");
                a.add(&mut r, &ex, &ey);
                assert_eq!(value(&r), Integer::from(x + y) % &m, "{x}+{y} ({form})");
                a.sub(&mut r, &ex, &ey);
                assert_eq!(
                    value(&r),
                    reduce(&Integer::from(x - y), &m),
                    "{x}-{y} ({form})"
                );
                acc.fill(0);
                a.mul_acc(&mut acc, &ex, &ey);
                a.mul_acc(&mut acc, &ex, &ey);
                a.redc_wide(&mut r, &acc);
                assert_eq!(
                    value(&r),
                    Integer::from(x * y) * 2u32 % &m,
                    "2*{x}*{y} ({form})"
                );
            }
        }
        for len in 1..=2 * a.limbs + 2 {
            let bits = 64 * len as u32;
            for t in [
                (Integer::from(1) << bits) - 1u32,
                Integer::from(1) << (bits - 1),
                Integer::from(Integer::random_bits(bits, rand)),
            ] {
                let mut digits = vec![0; len];
                t.write_digits(&mut digits, Order::Lsf);
                a.redc_wide(&mut r, &digits);
                assert_eq!(value(&r), Integer::from(&t % &m), "{t} ({form})");
            }
        }
        for x in [
            Integer::from(n * 3u32) + 5u32,
            Integer::from(-7),
            Integer::from(-n),
        ] {
            let e = a.residue(&x);
            assert_in_range(&e, &m, form);
            assert_eq!(a.to_integer(&e), reduce(&x, n));
        }
        // The sliding windows exponentiation of P-1.
        let x = Integer::from(n.random_below_ref(rand));
        for bits in [1, 7, 300] {
            let e = Integer::from(Integer::random_bits(bits, rand)) | 1u32;
            let expected = x.clone().pow_mod(&e, n).unwrap();
            assert_eq!(
                a.to_integer(&crate::arith::pow(&a, &a.residue(&x), &e)),
                expected
            );
        }
    }

    /// Polynomial products of polynomials whose coefficients all are the largest value of the
    /// representation: the Kronecker slots hold the largest sums of products.
    #[test]
    fn worst_case_slots() {
        use crate::poly::{Operand, Workspace, mul, mul_part, mul_wrap};
        let small = |k: i64, of: i64| {
            (
                Base2Form::from_signed(k).unwrap().value(),
                cofactor_of(of).1,
            )
        };
        for (n, form) in [64, -64, -127, 1024, -1061, 1025]
            .map(cofactor_of)
            .into_iter()
            .chain([small(64, 1088), small(-61, -610)])
        {
            let a = Base2::new(&n, form);
            let top = if form.plus {
                form.value() - 1u32
            } else {
                form.value()
            };
            let mut largest = a.zero();
            top.write_digits(&mut largest, Order::Lsf);
            let square = Integer::from(a.to_integer(&largest).square_ref());
            let mut ws = Workspace::new();
            for (lx, ly) in [(13, 13), (14, 13), (100, 100), (255, 256), (700, 13)] {
                let (x, y) = (vec![largest.clone(); lx], vec![largest.clone(); ly]);
                // Coefficient t is the sum of min(t + 1, lx, ly, lx + ly - 1 - t) squares.
                let expected: Vec<Integer> = (0..lx + ly - 1)
                    .map(|t| {
                        let terms = (t + 1).min(lx).min(ly).min(lx + ly - 1 - t);
                        Integer::from(&square * terms as u32) % &n
                    })
                    .collect();
                let values = |x: &[Vec<u64>]| x.iter().map(|x| a.to_integer(x)).collect::<Vec<_>>();
                let mut r = vec![a.zero(); lx + ly - 1];
                mul(&a, &mut ws, &mut r, &x, &y);
                assert_eq!(values(&r), expected, "{form} {lx} {ly}");
                for from in [1, (lx + ly) / 3, lx.max(ly) - 1] {
                    let mut part = vec![a.zero(); lx + ly - 1 - from];
                    let (ox, oy) = (Operand::new(&x), Operand::new(&y));
                    mul_part(&a, &mut ws, &mut part, from, ox, oy);
                    assert_eq!(values(&part), expected[from..], "{form} {lx} {ly} {from}");
                }
                for lmin in [lx.max(ly), lx.max(ly) + 5] {
                    let mut wrap = vec![a.zero(); lmin];
                    let (ox, oy) = (Operand::new(&x), Operand::new(&y));
                    let l = mul_wrap(&a, &mut ws, &mut wrap, 0, ox, oy, lmin);
                    for (t, w) in values(&wrap).iter().enumerate().take(l.min(lmin)) {
                        let mut c = expected.get(t).cloned().unwrap_or_default();
                        c += expected.get(t + l).cloned().unwrap_or_default();
                        assert_eq!(*w, c % &n, "{form} {lx} {ly} wrap {t}");
                    }
                }
            }
        }
    }
}
