//! Stage 1 on Edwards curves `x^2 + y^2 = 1 + d*x^2*y^2` with torsion group `Z/12`
//! ([`Param::Edwards12`](crate::Param::Edwards12)).
//!
//! The curves are the ones of Bernstein, Birkner, Lange and Peters, "ECM using Edwards curves"
//! (EECM-MPFQ, <https://eprint.iacr.org/2008/016>, theorem 7.8, after Montgomery): for a point
//! `(s, t)` of `t^2 = s^3 - 12*s`, the curve with
//! `d = -(s - 2)^3*(s + 6)^3*(s^2 - 12*s - 12)/(1024*s^2*t^2)` has torsion group `Z/12` over
//! `Q` and the non-torsion point
//! `x1 = 8*t*(s^2 + 12)/((s - 2)*(s + 6)*(s^2 + 12*s - 12))`,
//! `y1 = -4*s*(s^2 - 12*s - 12)/((s - 2)*(s + 6)*(s^2 - 12))`. Here `(s, t) = sigma*(-2, -4)`.
//!
//! Stage 1 computes `k*P` with a signed sliding window (width-`w` NAF) in extended coordinates
//! `(X : Y : Z : T)` (Hisil, Wong, Carter and Dawson, "Twisted Edwards curves revisited"): a
//! doubling costs 3M + 4S (4M + 4S when the next operation is an addition, which needs `T`),
//! and the addition of a precomputed odd multiple normalized to `Z = 1` 7M, once every
//! `w + 2` bits on average; the Montgomery ladder costs 5M + 4S per bit with the curves of
//! parametrization 2. The result is mapped to the birationally equivalent Montgomery curve
//! `B*v^2 = u^3 + A*u^2 + u` (`u = (1 + y)/(1 - y)`, `(A + 2)/4 = 1/(1 - d)`), where stage 2
//! runs unchanged.
//!
//! These are `a = 1` curves: twisted Edwards curves with `a = -1`, whose additions cost 1M less,
//! cannot have torsion group `Z/12` or `Z/2 x Z/8` over `Q` (EECM-MPFQ, theorem 6.11), and the
//! additions are only one bit in `w + 2`.

use crate::{
    arith::Arith,
    stop::{STOP_INTERVAL, Stop},
};
use rug::Integer;

/// Starting point `(x, y)` of an Edwards curve with parameter `d`, in affine coordinates
/// modulo `n`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Edwards {
    /// `x` coordinate.
    pub x: Integer,
    /// `y` coordinate.
    pub y: Integer,
    /// Parameter `d` of the curve.
    pub d: Integer,
}

/// `r = c*x` for a small constant `c >= 1`, by additions.
fn mul_const<A: Arith>(a: &A, x: &A::Elem, c: u32) -> A::Elem {
    debug_assert!(c >= 1);
    let mut r = x.clone();
    let mut t = a.zero();
    for bit in (0..c.ilog2()).rev() {
        a.add(&mut t, &r, &r);
        if c >> bit & 1 == 1 {
            a.add(&mut r, &t, x);
        } else {
            std::mem::swap(&mut r, &mut t);
        }
    }
    r
}

/// `sigma*(-2 : -4 : 1)` in Jacobian coordinates on `t^2 = s^3 - 12*s`, for `sigma >= 1`.
fn multiple<A: Arith>(a: &A, sigma: &Integer) -> [A::Elem; 3] {
    let (px, py) = (a.residue(&Integer::from(-2)), a.residue(&Integer::from(-4)));
    let (mut x, mut y, mut z) = (px.clone(), py.clone(), a.residue(&Integer::from(1)));
    let [
        mut t0,
        mut t1,
        mut t2,
        mut t3,
        mut t4,
        mut t5,
        mut t6,
        mut t7,
    ] = std::array::from_fn(|_| a.zero());
    for bit in (0..sigma.significant_bits() - 1).rev() {
        // Doubling, "dbl-2007-bl" with a = -12.
        a.sqr(&mut t0, &x); // XX
        a.sqr(&mut t1, &y); // YY
        a.sqr(&mut t2, &t1); // YYYY
        a.sqr(&mut t3, &z); // ZZ
        a.add(&mut t4, &x, &t1);
        a.sqr(&mut t5, &t4);
        a.sub(&mut t4, &t5, &t0);
        a.sub(&mut t5, &t4, &t2);
        a.add(&mut t4, &t5, &t5); // S = 2*((X + YY)^2 - XX - YYYY)
        a.sqr(&mut t5, &t3);
        let t12 = mul_const(a, &t5, 12);
        let t0x3 = mul_const(a, &t0, 3);
        a.sub(&mut t5, &t0x3, &t12); // M = 3*XX - 12*ZZ^2
        a.sqr(&mut t6, &t5);
        a.sub(&mut t7, &t6, &t4);
        a.sub(&mut t6, &t7, &t4); // T = M^2 - 2*S
        a.add(&mut t7, &y, &z);
        a.sqr(&mut z, &t7);
        a.sub(&mut t7, &z, &t1);
        a.sub(&mut z, &t7, &t3); // Z3 = (Y + Z)^2 - YY - ZZ
        a.sub(&mut t7, &t4, &t6);
        a.mul(&mut t0, &t5, &t7);
        let t2x8 = mul_const(a, &t2, 8);
        a.sub(&mut y, &t0, &t2x8); // Y3 = M*(S - T) - 8*YYYY
        std::mem::swap(&mut x, &mut t6); // X3 = T

        if sigma.get_bit(bit) {
            // Mixed addition of (-2, -4), "madd-2007-bl".
            a.sqr(&mut t0, &z); // Z1Z1 = z^2
            a.mul(&mut t1, &px, &t0); // U2 = px*Z1Z1
            a.mul(&mut t2, &py, &z);
            a.mul(&mut t3, &t2, &t0); // S2 = py*z*Z1Z1
            a.sub(&mut t2, &t1, &x); // H = U2 - x
            a.sqr(&mut t1, &t2); // HH = H^2
            a.add(&mut t4, &t1, &t1);
            a.add(&mut t5, &t4, &t4); // I = 4*HH
            a.mul(&mut t4, &t2, &t5); // J = H*I
            a.sub(&mut t6, &t3, &y);
            a.add(&mut t3, &t6, &t6); // r = 2*(S2 - y)
            a.mul(&mut t6, &x, &t5); // V = x*I
            a.sqr(&mut t5, &t3);
            a.sub(&mut t7, &t5, &t4);
            a.sub(&mut t5, &t7, &t6);
            a.sub(&mut x, &t5, &t6); // x3 = r^2 - J - 2*V
            a.sub(&mut t5, &t6, &x);
            a.mul(&mut t6, &t3, &t5);
            a.mul(&mut t5, &y, &t4);
            a.add(&mut t7, &t5, &t5);
            a.sub(&mut y, &t6, &t7); // y3 = r*(V - x3) - 2*y*J
            a.add(&mut t5, &z, &t2);
            a.sqr(&mut t6, &t5);
            a.sub(&mut t5, &t6, &t0);
            a.sub(&mut z, &t5, &t1); // z3 = (z + H)^2 - Z1Z1 - HH
        }
    }
    [x, y, z]
}

/// The curve of [`Param::Edwards12`](crate::Param::Edwards12) given by `sigma >= 2`: its
/// starting point and `d`, and `a24 = (A + 2)/4 = 1/(1 - d)` of the equivalent Montgomery
/// curve. One modular inversion.
///
/// # Errors
///
/// `Err(g)` with `g = gcd(D, n)` (possibly `n`), where `D` is the product of the denominators
/// and of the numerator of `d`: `sigma*(-2, -4)` is the point at infinity or one of the
/// excluded points modulo a factor of `n`, or the curve is singular modulo it.
pub(crate) fn curve<A: Arith>(a: &A, sigma: &Integer) -> Result<(Edwards, Integer), Integer> {
    let [x, y, z] = multiple(a, sigma);
    let [mut t0, mut t1] = std::array::from_fn(|_| a.zero());
    // With (s, t) = (X/Z^2, Y/Z^3) and Z2 = Z^2, everything times a power of Z:
    // d = dn/dd with dn = -(X - 2*Z2)^3*(X + 6*Z2)^3*c, c = X^2 - 12*X*Z2 - 12*Z2^2 and
    // dd = 1024*X^2*Y^2*Z^6, x1 = nx/dx with nx = 8*Y*Z*(X^2 + 12*Z2^2) and
    // dx = (X - 2*Z2)*(X + 6*Z2)*(X^2 + 12*X*Z2 - 12*Z2^2), y1 = ny/dy with
    // ny = -4*X*Z2*c and dy = (X - 2*Z2)*(X + 6*Z2)*(X^2 - 12*Z2^2).
    let mut z2 = a.zero();
    a.sqr(&mut z2, &z);
    let mut xx = a.zero();
    a.sqr(&mut xx, &x);
    let mut z4 = a.zero();
    a.sqr(&mut z4, &z2);
    let z4x12 = mul_const(a, &z4, 12);
    a.mul(&mut t0, &x, &z2);
    let xz2x12 = mul_const(a, &t0, 12);
    let z2x2 = mul_const(a, &z2, 2);
    let z2x6 = mul_const(a, &z2, 6);
    let mut am = a.zero(); // X - 2*Z2
    a.sub(&mut am, &x, &z2x2);
    let mut bp = a.zero(); // X + 6*Z2
    a.add(&mut bp, &x, &z2x6);
    let mut ab = a.zero();
    a.mul(&mut ab, &am, &bp);
    let mut c = a.zero(); // X^2 - 12*X*Z2 - 12*Z2^2
    a.sub(&mut t0, &xx, &xz2x12);
    a.sub(&mut c, &t0, &z4x12);
    // dn = -(ab)^3 * c
    let mut dn = a.zero();
    a.sqr(&mut t0, &ab);
    a.mul(&mut t1, &t0, &ab);
    a.mul(&mut t0, &t1, &c);
    a.sub(&mut dn, &a.zero(), &t0);
    // dd = 1024*(X*Y*Z^3)^2
    let mut dd = a.zero();
    a.mul(&mut t0, &x, &y);
    a.mul(&mut t1, &t0, &z2);
    a.mul(&mut t0, &t1, &z);
    a.sqr(&mut t1, &t0);
    for _ in 0..10 {
        a.add(&mut t0, &t1, &t1);
        std::mem::swap(&mut t0, &mut t1);
    }
    std::mem::swap(&mut dd, &mut t1);
    // 1 - d = (dd - dn)/dd
    let mut de = a.zero();
    a.sub(&mut de, &dd, &dn);
    // nx, dx
    let mut nx = a.zero();
    a.add(&mut t0, &xx, &z4x12);
    a.mul(&mut t1, &t0, &z);
    a.mul(&mut t0, &t1, &y);
    t1 = mul_const(a, &t0, 8);
    std::mem::swap(&mut nx, &mut t1);
    let mut dx = a.zero();
    a.add(&mut t0, &xx, &xz2x12);
    a.sub(&mut t1, &t0, &z4x12);
    a.mul(&mut dx, &ab, &t1);
    // ny, dy
    let mut ny = a.zero();
    a.mul(&mut t0, &x, &z2);
    a.mul(&mut t1, &t0, &c);
    t0 = mul_const(a, &t1, 4);
    a.sub(&mut ny, &a.zero(), &t0);
    let mut dy = a.zero();
    a.sub(&mut t0, &xx, &z4x12);
    a.mul(&mut dy, &ab, &t0);

    // One inversion for dd, de, dx and dy (and the check of dn).
    let values = [&dd, &de, &dx, &dy, &dn];
    let mut prefix = Vec::with_capacity(values.len());
    let mut acc = values[0].clone();
    prefix.push(acc.clone());
    for v in &values[1..] {
        a.mul(&mut t0, &acc, v);
        std::mem::swap(&mut acc, &mut t0);
        prefix.push(acc.clone());
    }
    let product = a.to_integer(&acc);
    let Ok(inv) = product.invert(a.modulus()) else {
        return Err(a.gcd(&acc));
    };
    // inverses[i] = 1/values[i]
    let mut run = a.residue(&inv);
    let mut inverses: Vec<A::Elem> = vec![a.zero(); values.len()];
    for i in (1..values.len()).rev() {
        a.mul(&mut inverses[i], &run, &prefix[i - 1]);
        a.mul(&mut t0, &run, values[i]);
        std::mem::swap(&mut run, &mut t0);
    }
    inverses[0] = run;
    a.mul(&mut t0, &dn, &inverses[0]);
    let d = a.to_integer(&t0);
    a.mul(&mut t0, &dd, &inverses[1]);
    let a24 = a.to_integer(&t0);
    a.mul(&mut t0, &nx, &inverses[2]);
    let x1 = a.to_integer(&t0);
    a.mul(&mut t0, &ny, &inverses[3]);
    let y1 = a.to_integer(&t0);
    Ok((Edwards { x: x1, y: y1, d }, a24))
}

/// Point in extended coordinates `(X : Y : Z : T)`, `x = X/Z`, `y = Y/Z`, `x*y = T/Z`.
struct Ext<E> {
    x: E,
    y: E,
    z: E,
    t: E,
}

/// Precomputed affine point: `x`, `y`, `y + x`, `y - x` and `d*x*y`.
struct Affine<E> {
    x: E,
    y: E,
    ypx: E,
    ymx: E,
    dt: E,
}

/// Temporaries of the point operations.
struct Scratch<E>([E; 5]);

/// `p = 2*p` ("dbl-2008-hwcd" with `a = 1`): 3M + 4S, and one more M for `T` if `with_t`.
#[inline(always)]
fn double<A: Arith>(a: &A, p: &mut Ext<A::Elem>, s: &mut Scratch<A::Elem>, with_t: bool) {
    let [s0, s1, s2, s3, s4] = &mut s.0;
    a.sqr(s0, &p.x); // A = X^2
    a.sqr(s1, &p.y); // B = Y^2
    a.add(s2, &p.x, &p.y);
    a.sqr(s3, s2);
    a.sub(s2, s3, s0);
    a.sub(s3, s2, s1); // E = (X + Y)^2 - A - B
    a.add(s2, s0, s1); // G = A + B
    a.sub(s4, s0, s1); // H = A - B
    a.sqr(s0, &p.z);
    a.add(s1, s0, s0); // C = 2*Z^2
    a.sub(s0, s2, s1); // F = G - C
    a.mul(&mut p.x, s3, s0); // E*F
    a.mul(&mut p.y, s2, s4); // G*H
    a.mul(&mut p.z, s0, s2); // F*G
    if with_t {
        a.mul(&mut p.t, s3, s4); // E*H
    }
}

/// `p = p + q` or `p = p - q` if `neg` ("madd-2008-hwcd" with `a = 1`, `q` normalized to
/// `Z = 1`): 7M, and one more M for `T` if `with_t`. `p` must have `T`.
#[inline(always)]
fn add<A: Arith>(
    a: &A,
    p: &mut Ext<A::Elem>,
    q: &Affine<A::Elem>,
    neg: bool,
    s: &mut Scratch<A::Elem>,
    with_t: bool,
) {
    let [s0, s1, s2, s3, s4] = &mut s.0;
    // With -q = (-x, y): A and C change sign, and x + y becomes y - x.
    a.mul(s0, &p.x, &q.x); // +-A = X1*x2
    a.mul(s1, &p.y, &q.y); // B = Y1*y2
    a.mul(s2, &p.t, &q.dt); // +-C = T1*d*x2*y2
    a.add(s3, &p.x, &p.y);
    a.mul(s4, s3, if neg { &q.ymx } else { &q.ypx });
    if neg {
        a.add(s3, s4, s0);
        a.sub(s4, s3, s1); // E = (X1 + Y1)*(y2 - x2) + A - B
        a.add(s3, &p.z, s2); // F = Z1 + C
        a.add(&mut p.t, s1, s0); // H = B + A
        a.sub(s0, &p.z, s2); // G = Z1 - C
    } else {
        a.sub(s3, s4, s0);
        a.sub(s4, s3, s1); // E = (X1 + Y1)*(x2 + y2) - A - B
        a.sub(s3, &p.z, s2); // F = Z1 - C
        a.sub(&mut p.t, s1, s0); // H = B - A
        a.add(s0, &p.z, s2); // G = Z1 + C
    }
    a.mul(&mut p.x, s4, s3); // E*F
    a.mul(&mut p.y, s0, &p.t); // G*H
    a.mul(&mut p.z, s3, s0); // F*G
    if with_t {
        a.mul(s1, s4, &p.t); // E*H
        std::mem::swap(&mut p.t, s1);
    }
}

/// `p = p + q` for two points with `T` ("add-2008-hwcd" with `a = 1`, `d2 = d*T2`
/// precomputed): 8M, and one more for `T`.
fn add_full<A: Arith>(
    a: &A,
    p: &Ext<A::Elem>,
    q: &Ext<A::Elem>,
    dt2: &A::Elem,
    s: &mut Scratch<A::Elem>,
) -> Ext<A::Elem> {
    let [s0, s1, s2, s3, s4] = &mut s.0;
    let mut r = Ext {
        x: a.zero(),
        y: a.zero(),
        z: a.zero(),
        t: a.zero(),
    };
    a.mul(s0, &p.x, &q.x); // A
    a.mul(s1, &p.y, &q.y); // B
    a.mul(s2, &p.t, dt2); // C
    a.mul(&mut r.z, &p.z, &q.z); // D
    a.add(s3, &p.x, &p.y);
    a.add(s4, &q.x, &q.y);
    a.mul(&mut r.t, s3, s4);
    a.sub(s3, &r.t, s0);
    a.sub(s4, s3, s1); // E
    a.sub(s3, &r.z, s2); // F = D - C
    a.add(&mut r.x, &r.z, s2); // G = D + C
    a.sub(s2, s1, s0); // H = B - A
    a.mul(&mut r.y, &r.x, s2); // G*H
    a.mul(&mut r.z, s3, &r.x); // F*G
    a.mul(&mut r.t, s4, s2); // E*H
    a.mul(&mut r.x, s4, s3); // E*F
    r
}

/// Width of the windows for a multiplier of `bits` bits modulo a number of `limbs` limbs: the
/// table of the `2^(w - 1)` odd multiples up to `(2^w - 1)*P` costs about 16M per point, and an
/// addition 8M every `w + 2` bits on average. The table stays within [`MAX_TABLE_BYTES`]:
/// larger ones fall out of the L2 cache, and stage 1 runs slower (measured at 1024 bits with
/// `B1 = 1M`: 20% slower with 2048 points than with 256).
fn window(bits: u32, limbs: usize) -> u32 {
    let entry = 5 * 8 * limbs.max(1);
    (1..=MAX_WINDOW)
        .filter(|&w| w == 1 || (entry << (w - 1)) <= MAX_TABLE_BYTES)
        .min_by_key(|&w| u64::from(bits) * 8 / u64::from(w + 2) + (16u64 << (w - 1)))
        .unwrap()
}

/// Largest window width (a table of 2^(MAX_WINDOW - 1) points).
const MAX_WINDOW: u32 = 12;

/// Largest size of the table of precomputed points (5 residues per point).
const MAX_TABLE_BYTES: usize = 256 << 10;

/// Signed digits of `k >= 1` with windows of width `w + 1`: `k = sum digit*2^position` with
/// odd digits in `(-2^w, 2^w)`, from the most significant one.
fn wnaf(k: &Integer, w: u32) -> Vec<(u32, i32)> {
    let len = k.significant_bits();
    let width = w + 1;
    let mut digits = Vec::new();
    let mut carry = 0u32;
    let mut bit = 0;
    while bit < len {
        if u32::from(k.get_bit(bit)) == carry {
            bit += 1;
            continue;
        }
        let now = width.min(len - bit);
        let mut word = 0i64;
        for i in (0..now).rev() {
            word = word << 1 | i64::from(k.get_bit(bit + i));
        }
        word += i64::from(carry);
        carry = u32::from(word >> (width - 1) & 1 == 1);
        word -= i64::from(carry) << width;
        digits.push((bit, word as i32));
        bit += now;
    }
    if carry == 1 {
        digits.push((len, 1));
    }
    digits.reverse();
    digits
}

/// Stage 1 on the Edwards curve of `e`: `k*P` for `k >= 1`, as the point `(u : w)` of the
/// equivalent Montgomery curve. Stops early (with a meaningless point) if `stop` is requested.
///
/// If a factor of `n` shows up when normalizing the precomputed points, returns `(1 : g)` with
/// `gcd(g, n)` that factor (or `n`).
pub(crate) fn stage1<A: Arith>(a: &A, e: &Edwards, k: &Integer, stop: Stop<'_>) -> [A::Elem; 2] {
    #[cfg(target_arch = "x86_64")]
    if crate::arith::has_bmi2_adx() {
        // SAFETY: the CPU has the features `stage1_bmi2` is compiled for.
        return unsafe { stage1_bmi2(a, e, k, stop) };
    }
    stage1_generic(a, e, k, stop)
}

/// [`stage1`] compiled with BMI2 and ADX (see [`crate::arith::has_bmi2_adx`]).
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "bmi2,adx")]
fn stage1_bmi2<A: Arith>(a: &A, e: &Edwards, k: &Integer, stop: Stop<'_>) -> [A::Elem; 2] {
    stage1_generic(a, e, k, stop)
}

#[inline(always)]
fn stage1_generic<A: Arith>(a: &A, e: &Edwards, k: &Integer, stop: Stop<'_>) -> [A::Elem; 2] {
    let limbs = a.modulus().significant_bits().div_ceil(64) as usize;
    let w = window(k.significant_bits(), limbs);
    let digits = wnaf(k, w);
    let table = match table(a, e, 1 << (w - 1)) {
        Ok(table) => table,
        Err(g) => return [a.residue(&Integer::from(1)), g],
    };
    let mut s = Scratch(std::array::from_fn(|_| a.zero()));
    let (top, first) = digits[0];
    let q = &table[(first.unsigned_abs() / 2) as usize];
    let mut p = Ext {
        x: q.x.clone(),
        y: q.y.clone(),
        z: a.residue(&Integer::from(1)),
        t: a.zero(),
    };
    if first < 0 {
        let x = std::mem::replace(&mut p.x, a.zero());
        a.sub(&mut p.x, &a.zero(), &x);
    }
    let mut position = top;
    let mut since_check = 0;
    for &(next, digit) in digits[1..].iter().chain(std::iter::once(&(0, 0))) {
        let doublings = position - next;
        for i in 0..doublings {
            double(a, &mut p, &mut s, digit != 0 && i + 1 == doublings);
        }
        if digit != 0 {
            let q = &table[(digit.unsigned_abs() / 2) as usize];
            add(a, &mut p, q, digit < 0, &mut s, false);
        }
        position = next;
        since_check += doublings;
        if since_check >= STOP_INTERVAL {
            since_check = 0;
            if stop.requested() {
                break;
            }
        }
    }
    // u = (1 + y)/(1 - y) = (Z + Y)/(Z - Y)
    let mut u = a.zero();
    let mut w = a.zero();
    a.add(&mut u, &p.z, &p.y);
    a.sub(&mut w, &p.z, &p.y);
    [u, w]
}

/// The odd multiples `P, 3P, ..., (2*count - 1)*P`, normalized with one inversion.
///
/// # Errors
///
/// `Err(g)` where `gcd(g, n) != 1` if the normalization fails.
fn table<A: Arith>(a: &A, e: &Edwards, count: usize) -> Result<Vec<Affine<A::Elem>>, A::Elem> {
    let mut s = Scratch(std::array::from_fn(|_| a.zero()));
    let one = a.residue(&Integer::from(1));
    let d = a.residue(&e.d);
    let (x, y) = (a.residue(&e.x), a.residue(&e.y));
    let mut t = a.zero();
    a.mul(&mut t, &x, &y);
    let p = Ext {
        x,
        y,
        z: one.clone(),
        t,
    };
    let mut points = Vec::with_capacity(count);
    if count > 1 {
        let mut p2 = Ext {
            x: p.x.clone(),
            y: p.y.clone(),
            z: p.z.clone(),
            t: a.zero(),
        };
        double(a, &mut p2, &mut s, true);
        let mut dt2 = a.zero();
        a.mul(&mut dt2, &p2.t, &d);
        points.push(p);
        for i in 1..count {
            let next = add_full(a, &points[i - 1], &p2, &dt2, &mut s);
            points.push(next);
        }
    } else {
        points.push(p);
    }
    // Batch inversion of the Z (Montgomery's trick).
    let mut prefix = Vec::with_capacity(count);
    let mut acc = one.clone();
    let mut tmp = a.zero();
    for point in &points {
        a.mul(&mut tmp, &acc, &point.z);
        std::mem::swap(&mut acc, &mut tmp);
        prefix.push(acc.clone());
    }
    let Ok(inv) = a.to_integer(&acc).invert(a.modulus()) else {
        return Err(acc);
    };
    let mut run = a.residue(&inv);
    let mut table: Vec<Affine<A::Elem>> = Vec::with_capacity(count);
    for i in (0..count).rev() {
        let mut zi = a.zero();
        if i > 0 {
            a.mul(&mut zi, &run, &prefix[i - 1]);
            a.mul(&mut tmp, &run, &points[i].z);
            std::mem::swap(&mut run, &mut tmp);
        } else {
            zi = run.clone();
        }
        let point = &points[i];
        let (mut x, mut y, mut xy, mut dt) = (a.zero(), a.zero(), a.zero(), a.zero());
        a.mul(&mut x, &point.x, &zi);
        a.mul(&mut y, &point.y, &zi);
        a.mul(&mut xy, &x, &y);
        a.mul(&mut dt, &xy, &d);
        let (mut ypx, mut ymx) = (a.zero(), a.zero());
        a.add(&mut ypx, &y, &x);
        a.sub(&mut ymx, &y, &x);
        table.push(Affine { x, y, ypx, ymx, dt });
    }
    table.reverse();
    Ok(table)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wnaf_digits() {
        let mut rand = rug::rand::RandState::new();
        for w in 1..=MAX_WINDOW {
            for bits in [1, 2, 3, 10, 64, 65, 200, 1000] {
                for _ in 0..20 {
                    let k = Integer::from(Integer::random_bits(bits, &mut rand)) + 1u32;
                    let digits = wnaf(&k, w);
                    let mut sum = Integer::new();
                    let mut last = u32::MAX;
                    for &(position, digit) in &digits {
                        assert!(digit % 2 != 0 && digit.unsigned_abs() < 1 << w, "{k} {w}");
                        // Nonzero digits are at least w + 1 positions apart.
                        assert!(last == u32::MAX || last > position + w);
                        last = position;
                        sum += Integer::from(digit) << position;
                    }
                    assert_eq!(sum, k, "{w}");
                }
            }
        }
    }
}
