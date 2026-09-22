//! Polynomial arithmetic modulo `n`, for the polynomial stage 2 (see [`crate::stage2_poly`]).
//!
//! Coefficients are residues in the polynomial representation of [`PolyArith`]. Products of
//! large polynomials use Kronecker substitution: each polynomial is packed into one big integer
//! (a coefficient every `s` bits, with `s` large enough that the coefficients of the product
//! don't overlap), the integers are multiplied by GMP (FFT multiplication for large sizes) and
//! the coefficients of the product are unpacked and reduced. Small products are schoolbook,
//! with one reduction per coefficient too.
//!
//! A monic polynomial of degree `d` is stored as its `d` low coefficients, the leading `1` is
//! implicit.

use crate::arith::PolyArith;
use rug::{integer::Order, Assign, Integer};

/// Products with a factor of at most this many coefficients are schoolbook.
pub const SCHOOLBOOK: usize = 12;

/// Reusable buffers of the polynomial products.
pub struct Workspace {
    x: Integer,
    y: Integer,
    product: Integer,
    limbs: Vec<u64>,
    acc: Vec<u64>,
}

impl Workspace {
    pub fn new() -> Self {
        Workspace {
            x: Integer::new(),
            y: Integer::new(),
            product: Integer::new(),
            limbs: Vec::new(),
            acc: Vec::new(),
        }
    }
}

/// `r = x*y`, with `r.len() = x.len() + y.len() - 1` (`x` and `y` not empty).
pub fn mul<A: PolyArith>(
    a: &A,
    ws: &mut Workspace,
    r: &mut [A::Elem],
    x: &[A::Elem],
    y: &[A::Elem],
) {
    debug_assert_eq!(r.len() + 1, x.len() + y.len());
    if x.len().min(y.len()) <= SCHOOLBOOK {
        schoolbook(a, ws, r, x, y);
    } else {
        kronecker(a, ws, r, x, y);
    }
}

/// Schoolbook product: each coefficient of `r` is accumulated on limbs, then reduced once.
fn schoolbook<A: PolyArith>(
    a: &A,
    ws: &mut Workspace,
    r: &mut [A::Elem],
    x: &[A::Elem],
    y: &[A::Elem],
) {
    let width = 2 * a.limbs() + 2;
    ws.acc.resize(width, 0);
    for (k, r) in r.iter_mut().enumerate() {
        ws.acc.fill(0);
        let lo = k.saturating_sub(y.len() - 1);
        let hi = k.min(x.len() - 1);
        for i in lo..=hi {
            a.mul_acc(&mut ws.acc, &x[i], &y[k - i]);
        }
        a.redc_wide(r, &ws.acc);
    }
}

/// Bits per coefficient in a Kronecker product where the shorter factor has `len`
/// coefficients: each coefficient of the product is a sum of at most `len` products of values
/// `< n`.
fn slot_bits(n: &Integer, len: usize) -> usize {
    2 * n.significant_bits() as usize + (usize::BITS - len.leading_zeros()) as usize
}

/// Packs the values of `x` into `buf`, one every `s` bits.
fn pack<A: PolyArith>(a: &A, buf: &mut Vec<u64>, x: &[A::Elem], s: usize) {
    let limbs = a.limbs();
    buf.clear();
    buf.resize((x.len() * s).div_ceil(64) + limbs + 1, 0);
    let mut v = vec![0; limbs];
    for (i, x) in x.iter().enumerate() {
        a.write_limbs(x, &mut v);
        let (word, shift) = ((i * s) / 64, (i * s) % 64);
        if shift == 0 {
            for (j, &v) in v.iter().enumerate() {
                buf[word + j] |= v;
            }
        } else {
            for (j, &v) in v.iter().enumerate() {
                buf[word + j] |= v << shift;
                buf[word + j + 1] |= v >> (64 - shift);
            }
        }
    }
}

/// Kronecker substitution: `r = x*y` with one big integer product.
fn kronecker<A: PolyArith>(
    a: &A,
    ws: &mut Workspace,
    r: &mut [A::Elem],
    x: &[A::Elem],
    y: &[A::Elem],
) {
    let s = slot_bits(a.modulus(), x.len().min(y.len()));
    pack(a, &mut ws.limbs, x, s);
    ws.x.assign_digits(&ws.limbs, Order::Lsf);
    pack(a, &mut ws.limbs, y, s);
    ws.y.assign_digits(&ws.limbs, Order::Lsf);
    ws.product.assign(&ws.x * &ws.y);

    // Unpack: coefficient k is the s bits at k*s.
    let width = s.div_ceil(64);
    let buf = &mut ws.limbs;
    buf.clear();
    buf.resize((r.len() * s).div_ceil(64) + width + 1, 0);
    let len = ws.product.significant_digits::<u64>();
    ws.product.write_digits(&mut buf[..len], Order::Lsf);
    let top_mask = if s.is_multiple_of(64) {
        u64::MAX
    } else {
        (1 << (s % 64)) - 1
    };
    ws.acc.resize(width, 0);
    for (k, r) in r.iter_mut().enumerate() {
        let (word, shift) = ((k * s) / 64, (k * s) % 64);
        let t = &mut ws.acc;
        if shift == 0 {
            t.copy_from_slice(&buf[word..word + width]);
        } else {
            for j in 0..width {
                t[j] = (buf[word + j] >> shift) | (buf[word + j + 1] << (64 - shift));
            }
        }
        t[width - 1] &= top_mask;
        a.redc_wide(r, t);
    }
}

/// `r = x*y` for monic `x` and `y` (implicit leading coefficients): `r.len() = x.len() + y.len()`.
pub fn mul_monic<A: PolyArith>(
    a: &A,
    ws: &mut Workspace,
    r: &mut [A::Elem],
    x: &[A::Elem],
    y: &[A::Elem],
) {
    let (lx, ly) = (x.len(), y.len());
    debug_assert_eq!(r.len(), lx + ly);
    mul(a, ws, &mut r[..lx + ly - 1], x, y);
    r[lx + ly - 1] = a.zero();
    let mut t = a.zero();
    for (r, y) in r[lx..].iter_mut().zip(y) {
        a.add(&mut t, r, y);
        std::mem::swap(r, &mut t);
    }
    for (r, x) in r[ly..].iter_mut().zip(x) {
        a.add(&mut t, r, x);
        std::mem::swap(r, &mut t);
    }
}

/// Product tree of the polynomials `X - root`: level `i` holds the monic products of `2^i`
/// consecutive leaves (the last one possibly fewer), each stored at the position of its first
/// leaf, so that every level has one coefficient per leaf. Level 0 holds the `-root`.
pub struct ProductTree<E> {
    pub levels: Vec<Vec<E>>,
}

/// Degree of the node of level `level` starting at leaf `start`, in a tree of `len` leaves.
fn node_len(level: usize, start: usize, len: usize) -> usize {
    (1 << level).min(len - start)
}

/// Monic polynomial of the next level: products of adjacent pairs of nodes of `level`.
fn next_level<A: PolyArith>(
    a: &A,
    ws: &mut Workspace,
    level: usize,
    cur: &[A::Elem],
    next: &mut [A::Elem],
) {
    let len = cur.len();
    let size = 1 << level;
    for start in (0..len).step_by(2 * size) {
        let left = node_len(level, start, len);
        if start + size >= len {
            // Single child.
            next[start..start + left].clone_from_slice(&cur[start..start + left]);
            continue;
        }
        let right = node_len(level, start + size, len);
        mul_monic(
            a,
            ws,
            &mut next[start..start + left + right],
            &cur[start..start + left],
            &cur[start + size..start + size + right],
        );
    }
}

/// Number of levels above the leaves: the root is at level `levels(len)`.
fn height(len: usize) -> usize {
    (usize::BITS - (len - 1).leading_zeros()) as usize
}

impl<E: Clone> ProductTree<E> {
    /// Product tree of the leaves `leaves` (the constant coefficients `-root`, not empty).
    pub fn new<A: PolyArith<Elem = E>>(a: &A, ws: &mut Workspace, leaves: Vec<E>) -> Self {
        let mut levels = vec![leaves];
        for level in 0..height(levels[0].len()) {
            let mut next = vec![a.zero(); levels[0].len()];
            next_level(a, ws, level, &levels[level], &mut next);
            levels.push(next);
        }
        ProductTree { levels }
    }

    /// The monic polynomial at the root, `prod (X - root)`.
    pub fn root(&self) -> &[E] {
        self.levels.last().unwrap()
    }
}

/// Monic `prod (X - root)` from the leaves `-root` (not empty), without keeping the tree:
/// `leaves` is overwritten.
pub fn from_roots<A: PolyArith>(
    a: &A,
    ws: &mut Workspace,
    leaves: &mut Vec<A::Elem>,
    tmp: &mut Vec<A::Elem>,
) {
    let len = leaves.len();
    tmp.resize(len, a.zero());
    for level in 0..height(len) {
        next_level(a, ws, level, leaves, tmp);
        std::mem::swap(leaves, tmp);
    }
}

/// Inverse of the power series `f` (`f[0] = 1`, as for the reverse of a monic polynomial)
/// modulo `X^len`, by Newton iteration: `g = g - g*(f*g - 1)`, doubling the precision.
pub fn inverse<A: PolyArith>(a: &A, ws: &mut Workspace, f: &[A::Elem], len: usize) -> Vec<A::Elem> {
    let mut g = vec![a.zero(); len];
    if len == 0 {
        return g;
    }
    g[0].clone_from(&f[0]);
    let mut prec = 1;
    let mut e = Vec::new();
    let mut d = Vec::new();
    while prec < len {
        let next = (2 * prec).min(len);
        let fl = next.min(f.len());
        // f*g = 1 + X^prec * e (mod X^next).
        e.resize(fl + prec - 1, a.zero());
        mul(a, ws, &mut e, &f[..fl], &g[..prec]);
        let h = next - prec;
        if e.len() < next {
            e.resize(next, a.zero());
        }
        // g[prec..next] = -(g*e)[..h].
        d.resize(2 * h - 1, a.zero());
        mul(a, ws, &mut d, &g[..h], &e[prec..next]);
        let zero = a.zero();
        for (g, d) in g[prec..next].iter_mut().zip(&d) {
            a.sub(g, &zero, d);
        }
        prec = next;
    }
    g
}

/// Reverse of the monic `f` of degree `f.len()`, `X^deg f(1/X)`, as a power series of `len`
/// terms.
pub fn reverse_monic<A: PolyArith>(a: &A, f: &[A::Elem], len: usize) -> Vec<A::Elem> {
    let one = a.poly_from(&Integer::from(1));
    let mut r = vec![a.zero(); len.min(f.len() + 1)];
    r[0] = one;
    for (i, r) in r.iter_mut().enumerate().skip(1) {
        r.clone_from(&f[f.len() - i]);
    }
    r
}

/// `h = h*g mod f`, with `f` monic of degree `d = h.len() = g.len()` (so `h`, `g` of degree
/// `< d`), and `inv` the inverse of the reverse of `f` modulo `X^(d - 1)`.
pub fn mul_mod<A: PolyArith>(
    a: &A,
    ws: &mut Workspace,
    h: &mut [A::Elem],
    g: &[A::Elem],
    f: &[A::Elem],
    inv: &[A::Elem],
) {
    let d = f.len();
    let mut p = vec![a.zero(); 2 * d - 1];
    mul(a, ws, &mut p, h, g);
    if d == 1 {
        h[0].clone_from(&p[0]);
        return;
    }
    // Quotient: reverse(q) = reverse(p) / reverse(f) mod X^(d - 1).
    let rev_p: Vec<A::Elem> = p[d..].iter().rev().cloned().collect();
    let mut rev_q = vec![a.zero(); 2 * (d - 1) - 1];
    mul(a, ws, &mut rev_q, &rev_p, &inv[..d - 1]);
    let q: Vec<A::Elem> = rev_q[..d - 1].iter().rev().cloned().collect();
    // Remainder: p - q*f, whose degree is < d (q*X^d only has terms of degree >= d).
    let mut qf = vec![a.zero(); 2 * d - 2];
    mul(a, ws, &mut qf, &q, f);
    for i in 0..d {
        a.sub(&mut h[i], &p[i], &qf[i]);
    }
}

/// Values of `h` (degree `< d`) at the `d` roots of the product tree `tree` (in the order of the
/// leaves), with `inv` the inverse of the reverse of the root of the tree modulo `X^d`.
///
/// Transposed algorithm (Bostan, Lecerf and Schost, "Tellegen's principle into practice"): at a
/// node `P`, keep the first `deg P` coefficients `c_1, c_2, ...` of `(h mod P)/P` as a series in
/// `1/X`. For a child `Q` of `P = Q*R`, `(h mod Q)/Q` is the part of `R * (h mod P)/P` with
/// negative powers: a middle product by the sibling `R`. At a leaf `X - x`, `c_1 = h(x)`.
pub fn evaluate<A: PolyArith>(
    a: &A,
    ws: &mut Workspace,
    h: &[A::Elem],
    tree: &ProductTree<A::Elem>,
    inv: &[A::Elem],
) -> Vec<A::Elem> {
    let d = h.len();
    // Root: h/f = X^-1 * rev(h)(1/X) / rev(f)(1/X), so c_(t+1) = (rev(h) * inv)[t].
    let rev_h: Vec<A::Elem> = h.iter().rev().cloned().collect();
    let mut prod = vec![a.zero(); 2 * d - 1];
    mul(a, ws, &mut prod, &rev_h, &inv[..d]);
    let mut series = prod;
    series.truncate(d);
    let mut next = vec![a.zero(); d];
    let mut t = a.zero();
    let (mut rev, mut buf) = (Vec::new(), Vec::new());
    for level in (0..tree.levels.len() - 1).rev() {
        let size = 1 << level;
        let children = &tree.levels[level];
        for start in (0..d).step_by(2 * size) {
            let left = node_len(level, start, d);
            if start + size >= d {
                next[start..start + left].clone_from_slice(&series[start..start + left]);
                continue;
            }
            let right = node_len(level, start + size, d);
            let s = &series[start..start + left + right];
            let (q, r) = (
                &children[start..start + left],
                &children[start + size..start + size + right],
            );
            // Left child: sum_(i <= right) r_i * s[t + i] for t < left (r_right = 1).
            let out = &mut next[start..start + left];
            middle(a, ws, out, r, s, &mut rev, &mut buf, &mut t);
            // Right child: sum_(i <= left) q_i * s[t + i] for t < right.
            let out = &mut next[start + left..start + left + right];
            middle(a, ws, out, q, s, &mut rev, &mut buf, &mut t);
        }
        std::mem::swap(&mut series, &mut next);
    }
    series
}

/// `out[t] = sum_(i <= m) r_i * s[t + i]` for `t < out.len()`, with `r` monic of degree `m`
/// (`r_m = 1`) and `s.len() >= out.len() + m`. `rev` and `prod` are scratch space.
#[allow(clippy::too_many_arguments)]
fn middle<A: PolyArith>(
    a: &A,
    ws: &mut Workspace,
    out: &mut [A::Elem],
    r: &[A::Elem],
    s: &[A::Elem],
    rev: &mut Vec<A::Elem>,
    prod: &mut Vec<A::Elem>,
    t: &mut A::Elem,
) {
    let (l, m) = (out.len(), r.len());
    if l.min(m) <= SCHOOLBOOK {
        ws.acc.resize(2 * a.limbs() + 2, 0);
        for (k, out) in out.iter_mut().enumerate() {
            ws.acc.fill(0);
            for (r, s) in r.iter().zip(&s[k..]) {
                a.mul_acc(&mut ws.acc, r, s);
            }
            a.redc_wide(t, &ws.acc);
            a.add(out, t, &s[k + m]);
        }
        return;
    }
    // (reverse(r) * s)[t + m - 1] = sum_i r_i * s[t + i].
    rev.clear();
    rev.extend(r.iter().rev().cloned());
    let s_used = &s[..l + m - 1];
    prod.resize(m + s_used.len() - 1, a.zero());
    mul(a, ws, prod, rev, s_used);
    for (i, out) in out.iter_mut().enumerate() {
        a.add(out, &prod[i + m - 1], &s[i + m]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arith::{with_arith, Arith};
    use rug::rand::RandState;

    /// Plain integer polynomial helpers modulo `n`, on values.
    fn values<A: PolyArith>(a: &A, x: &[A::Elem]) -> Vec<Integer> {
        // Values of the polynomial representation: multiply by 1 (representation of 1).
        x.iter().map(|x| a.poly_value(x)).collect()
    }

    fn random_poly<A: PolyArith>(a: &A, len: usize, rand: &mut RandState<'_>) -> Vec<A::Elem> {
        (0..len)
            .map(|_| a.poly_from(&a.modulus().clone().random_below(rand)))
            .collect()
    }

    fn naive_mul(x: &[Integer], y: &[Integer], n: &Integer) -> Vec<Integer> {
        let mut r = vec![Integer::new(); x.len() + y.len() - 1];
        for (i, x) in x.iter().enumerate() {
            for (j, y) in y.iter().enumerate() {
                r[i + j] += Integer::from(x * y);
            }
        }
        r.iter().map(|r| Integer::from(r % n)).collect()
    }

    fn check_mul<A: PolyArith>(a: &A, rand: &mut RandState<'_>) {
        let n = a.modulus().clone();
        let mut ws = Workspace::new();
        for (lx, ly) in [
            (1, 1),
            (1, 5),
            (3, 4),
            (12, 13),
            (13, 13),
            (20, 7),
            (40, 41),
            (100, 3),
            (64, 64),
        ] {
            let x = random_poly(a, lx, rand);
            let mut y = random_poly(a, ly, rand);
            // Largest values: the slots must hold the sums of products.
            if ly > 2 {
                for y in &mut y[..2] {
                    *y = a.poly_from(&Integer::from(&n - 1));
                }
            }
            let mut r = vec![a.zero(); lx + ly - 1];
            mul(a, &mut ws, &mut r, &x, &y);
            assert_eq!(
                values(a, &r),
                naive_mul(&values(a, &x), &values(a, &y), &n),
                "{lx} {ly}"
            );
            // All n - 1.
            let big = vec![a.poly_from(&Integer::from(&n - 1)); lx];
            mul(
                a,
                &mut ws,
                &mut r,
                &big,
                &big[..].iter().cycle().take(ly).cloned().collect::<Vec<_>>(),
            );
            let vb = values(a, &big);
            let vy: Vec<Integer> = vb.iter().cycle().take(ly).cloned().collect();
            assert_eq!(values(a, &r), naive_mul(&vb, &vy, &n));
        }
    }

    #[test]
    fn products() {
        let mut rand = RandState::new();
        for bits in [3, 20, 64, 65, 128, 200, 256, 512, 1000, 1024, 1100] {
            for _ in 0..2 {
                let mut n = Integer::from(Integer::random_bits(bits, &mut rand));
                n.set_bit(0, true);
                n.set_bit(bits - 1, true);
                with_arith!(&n, |a| check_mul(&a, &mut rand));
            }
            // Largest modulus of the size.
            let n = (Integer::from(1) << bits) - 1u32;
            with_arith!(&n, |a| check_mul(&a, &mut rand));
        }
    }

    /// Polynomial of `len` coefficients whose values in the polynomial representation are all
    /// `n - 1`: the largest sums of products in a Kronecker slot or a schoolbook accumulator.
    fn largest_poly<A: PolyArith>(a: &A, len: usize) -> Vec<A::Elem> {
        let n = a.modulus();
        let r_prime = if crate::arith::mont_limbs(n) == 0 {
            Integer::from(1)
        } else {
            Integer::from(1) << (64 * (a.limbs() as u32 + 1))
        };
        let x = Integer::from(n - 1u32) * r_prime.invert(n).unwrap() % n;
        let x = a.poly_from(&x);
        let mut v = vec![0; a.limbs()];
        a.write_limbs(&x, &mut v);
        assert_eq!(
            Integer::from_digits(&v, Order::Lsf),
            Integer::from(n - 1u32)
        );
        vec![x; len]
    }

    fn check_random_products<A: PolyArith>(a: &A, rand: &mut RandState<'_>) {
        let n = a.modulus().clone();
        let mut ws = Workspace::new();
        let below = |max: usize, rand: &mut RandState<'_>| {
            Integer::from(max).random_below(rand).to_usize().unwrap()
        };
        for i in 0..12 {
            let lx = 1 + below(if i % 3 == 0 { 200 } else { 30 }, rand);
            let ly = 1 + below(if i % 4 == 0 { 200 } else { 30 }, rand);
            let (x, y) = if i % 2 == 0 {
                (largest_poly(a, lx), largest_poly(a, ly))
            } else {
                (random_poly(a, lx, rand), random_poly(a, ly, rand))
            };
            let (xv, yv) = (values(a, &x), values(a, &y));
            let mut r = vec![a.zero(); lx + ly - 1];
            mul(a, &mut ws, &mut r, &x, &y);
            assert_eq!(values(a, &r), naive_mul(&xv, &yv, &n), "{lx} {ly}");

            // Monic factors.
            let one = || std::iter::once(Integer::from(1));
            let (xm, ym): (Vec<_>, Vec<_>) = (
                xv.iter().cloned().chain(one()).collect(),
                yv.iter().cloned().chain(one()).collect(),
            );
            let mut r = vec![a.zero(); lx + ly];
            mul_monic(a, &mut ws, &mut r, &x, &y);
            assert_eq!(
                values(a, &r),
                naive_mul(&xm, &ym, &n)[..lx + ly],
                "{lx} {ly}"
            );

            // Middle product by the monic y, lx outputs.
            let s = if i % 2 == 0 {
                largest_poly(a, lx + ly)
            } else {
                random_poly(a, lx + ly, rand)
            };
            let sv = values(a, &s);
            let mut out = vec![a.zero(); lx];
            let (mut rev, mut buf, mut t) = (Vec::new(), Vec::new(), a.zero());
            middle(a, &mut ws, &mut out, &y, &s, &mut rev, &mut buf, &mut t);
            for (k, out) in values(a, &out).iter().enumerate() {
                let sum = ym
                    .iter()
                    .zip(&sv[k..])
                    .fold(Integer::new(), |acc, (y, s)| acc + y * s);
                assert_eq!(*out, sum % &n, "{lx} {ly}");
            }

            // Inverse of a series given by fewer terms than the precision.
            let mut f = x;
            f[0] = a.poly_from(&Integer::from(1));
            let len = lx + below(40, rand);
            let g = inverse(a, &mut ws, &f, len);
            let p = naive_mul(&values(a, &f), &values(a, &g), &n);
            assert!(
                p[..len]
                    .iter()
                    .enumerate()
                    .all(|(k, p)| *p == (k == 0) as u32),
                "{lx} {len}"
            );
        }
    }

    #[test]
    fn random_products() {
        // Random lengths around the schoolbook threshold and largest coefficients, for every
        // limb count and Plain (even or too large n).
        let mut rand = RandState::new();
        for limbs in 1..=17u32 {
            for n in [
                (Integer::from(1) << (64 * limbs)) - 1u32,
                (Integer::from(1) << (64 * limbs - 37)) + 1u32,
                (Integer::from(1) << (64 * limbs)) - 2u32,
            ] {
                with_arith!(&n, |a| check_random_products(&a, &mut rand));
            }
        }
    }

    /// Evaluates `x` at `v` modulo `n`.
    fn eval(x: &[Integer], v: &Integer, n: &Integer) -> Integer {
        x.iter()
            .rev()
            .fold(Integer::new(), |acc, c| (acc * v + c) % n)
    }

    fn check_tree<A: PolyArith>(a: &A, rand: &mut RandState<'_>) {
        let n = a.modulus().clone();
        let mut ws = Workspace::new();
        for d in [1, 2, 3, 5, 8, 13, 16, 31, 33, 100, 129, 300] {
            let roots: Vec<Integer> = (0..d).map(|_| n.clone().random_below(rand)).collect();
            let leaves: Vec<A::Elem> = roots
                .iter()
                .map(|r| a.poly_from(&Integer::from(&n - r)))
                .collect();
            let tree = ProductTree::new(a, &mut ws, leaves.clone());
            let mut f = values(a, tree.root());
            f.push(Integer::from(1));
            for r in &roots {
                assert_eq!(eval(&f, r, &n), 0, "d = {d}");
            }
            let mut leaves = leaves;
            let mut tmp = Vec::new();
            from_roots(a, &mut ws, &mut leaves, &mut tmp);
            assert_eq!(values(a, &leaves), values(a, tree.root()));

            // Inverse of the reverse.
            let rev = reverse_monic(a, tree.root(), d);
            let inv = inverse(a, &mut ws, &rev, d);
            let mut p = vec![a.zero(); 2 * d - 1];
            mul(
                a,
                &mut ws,
                &mut p,
                &rev.iter()
                    .cloned()
                    .chain(std::iter::repeat(a.zero()))
                    .take(d)
                    .collect::<Vec<_>>(),
                &inv,
            );
            let pv = values(a, &p[..d]);
            assert_eq!(pv[0], 1);
            assert!(pv[1..].iter().all(|v| *v == 0), "d = {d}");

            // h*g mod f, then evaluation at the roots.
            let mut h = random_poly(a, d, rand);
            let g = random_poly(a, d, rand);
            let (hv, gv) = (values(a, &h), values(a, &g));
            mul_mod(a, &mut ws, &mut h, &g, tree.root(), &inv);
            let got = evaluate(a, &mut ws, &h, &tree, &inv);
            for (r, got) in roots.iter().zip(values(a, &got)) {
                let expected = eval(&hv, r, &n) * eval(&gv, r, &n) % &n;
                assert_eq!(got, expected, "d = {d}");
            }
        }
    }

    #[test]
    fn tree_inverse_mul_mod_evaluate() {
        let mut rand = RandState::new();
        for bits in [20, 64, 100, 256, 1024, 1100] {
            let mut n = Integer::from(Integer::random_bits(bits, &mut rand));
            n.set_bit(0, true);
            n.set_bit(bits - 1, true);
            with_arith!(&n, |a| check_tree(&a, &mut rand));
        }
    }

    #[test]
    fn redc_wide_bounds() {
        // Largest input t < n*R': (2^(64(N+1)) - 1) * (n - 1) / 2^0 ... check against integers.
        let mut rand = RandState::new();
        for bits in [64, 128, 192, 1024] {
            for n in [
                (Integer::from(1) << bits) - 1u32,
                (Integer::from(1) << (bits - 1)) + 1u32,
            ] {
                with_arith!(&n, |a| {
                    let limbs = a.limbs();
                    let r_prime = Integer::from(1) << (64 * (limbs as u32 + 1));
                    let r_prime = if crate::arith::mont_limbs(&n) == 0 {
                        Integer::from(1)
                    } else {
                        r_prime
                    };
                    let max = Integer::from(&n * &r_prime) - 1u32;
                    for t in [
                        max.clone(),
                        Integer::from(&max >> 1),
                        n.clone().random_below(&mut rand),
                    ] {
                        let mut digits = vec![0; 2 * limbs + 2];
                        t.write_digits(&mut digits[..t.significant_digits::<u64>()], Order::Lsf);
                        let mut r = a.zero();
                        a.redc_wide(&mut r, &digits);
                        let inv = r_prime.clone().invert(&n).unwrap();
                        let expected = Integer::from(&t * &inv) % &n;
                        let mut v = vec![0; limbs];
                        a.write_limbs(&r, &mut v);
                        assert_eq!(Integer::from_digits(&v, Order::Lsf), expected);
                    }
                });
            }
        }
    }
}
