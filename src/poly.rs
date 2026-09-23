//! Polynomial arithmetic modulo `n`, for the polynomial stage 2 (see [`crate::stage2_poly`]).
//!
//! Coefficients are residues in the polynomial representation of [`PolyArith`]. Products of
//! large polynomials use Kronecker substitution: each polynomial is packed into one big integer
//! (a coefficient every `s` bits, with `s` large enough that the coefficients of the product
//! don't overlap), the integers are multiplied by GMP (FFT multiplication for large sizes) and
//! the coefficients of the product are unpacked and reduced. Small products are schoolbook,
//! with one reduction per coefficient too.
//!
//! When only some coefficients of a product are needed, a wrap-around product (modulo
//! `X^L - 1`, a product of integers modulo `2^(L*s) - 1`: GMP's `mpn_mulmod_bnm1`, about half
//! the cost of a full product) is often enough: for the middle products of the multipoint
//! evaluation and of the Newton iteration, and for `q*f` in a reduction modulo `f`, whose high
//! coefficients are known.
//!
//! A monic polynomial of degree `d` is stored as its `d` low coefficients, the leading `1` is
//! implicit.

use crate::arith::{mpn, PolyArith};
use rug::{integer::Order, Assign, Integer};
use std::collections::HashMap;

/// Products with a factor of at most this many coefficients are schoolbook.
pub const SCHOOLBOOK: usize = 12;

/// Whether the Kronecker products use GMP's low-level functions ([`mpn::ENABLED`]), without
/// which there are only full products (of `Integer`s). Tests also check them without.
fn low_level() -> bool {
    #[cfg(test)]
    if tests::NO_LOW_LEVEL.get() {
        return false;
    }
    mpn::ENABLED
}

/// Reusable buffers of the polynomial products.
pub struct Workspace {
    /// The packed operands.
    x: Vec<u64>,
    y: Vec<u64>,
    bufs: Buffers,
    /// Shapes of the wrap-around products, by `(bits, lx, ly, lmin)`: see [`wrap_size`].
    shapes: HashMap<(usize, usize, usize, usize), Option<Shape>>,
}

/// Buffers of a Kronecker product.
struct Buffers {
    product: Vec<u64>,
    scratch: Vec<u64>,
    acc: Vec<u64>,
    /// Without GMP's low-level functions: the operands and the product as integers.
    ints: [Integer; 3],
}

impl Workspace {
    pub fn new() -> Self {
        Workspace {
            x: Vec::new(),
            y: Vec::new(),
            bufs: Buffers {
                product: Vec::new(),
                scratch: Vec::new(),
                acc: Vec::new(),
                ints: [Integer::new(), Integer::new(), Integer::new()],
            },
            shapes: HashMap::new(),
        }
    }
}

/// A polynomial packed for Kronecker products, kept to be reused while the slot width and the
/// polynomial don't change (the polynomial is the caller's responsibility).
#[derive(Default)]
pub struct Packed {
    s: usize,
    len: usize,
    limbs: Vec<u64>,
    size: usize,
}

impl Packed {
    /// The limbs of `x` packed with slots of `s` bits, packed again if `s` or the length
    /// changed.
    fn get<A: PolyArith>(&mut self, a: &A, x: &[A::Elem], s: usize) -> &[u64] {
        if self.s != s || self.len != x.len() || self.limbs.is_empty() {
            self.size = pack(a, &mut self.limbs, x, s);
            (self.s, self.len) = (s, x.len());
        }
        &self.limbs[..self.size]
    }
}

/// An operand of a product: its coefficients, and maybe where to keep its packed form.
pub struct Operand<'a, 'b, E> {
    coeffs: &'a [E],
    cache: Option<&'b mut Packed>,
}

impl<'a, 'b, E> Operand<'a, 'b, E> {
    /// An operand packed again for each product.
    pub fn new(coeffs: &'a [E]) -> Self {
        Operand {
            coeffs,
            cache: None,
        }
    }

    /// An operand whose packed form is kept in `cache` (packed again if the slot width or
    /// the length change, but not if only the coefficients do).
    pub fn cached(coeffs: &'a [E], cache: &'b mut Packed) -> Self {
        Operand {
            coeffs,
            cache: Some(cache),
        }
    }

    fn len(&self) -> usize {
        self.coeffs.len()
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
    mul_part(a, ws, r, 0, Operand::new(x), Operand::new(y));
}

/// `out[t]` = coefficient `from + t` of `x*y` (zero above its degree), unpacking and reducing
/// only these: a short product when `from = 0`, a middle product otherwise.
///
/// For a middle product, the integer product modulo `2^N - 1` (GMP's `mpn_mulmod_bnm1`, about
/// half the cost of a full product for `N` about half its size) is enough when `N` covers the
/// needed coefficients and the part above `N`, which wraps around to the low bits, stays below
/// the needed ones: with a spare bit per slot, the low coefficients and the wrapped part both
/// are below `2^(from*s - 1)`, so their sum doesn't carry into the coefficient `from`.
pub fn mul_part<A: PolyArith>(
    a: &A,
    ws: &mut Workspace,
    out: &mut [A::Elem],
    from: usize,
    x: Operand<'_, '_, A::Elem>,
    y: Operand<'_, '_, A::Elem>,
) {
    let (lx, ly) = (x.len(), y.len());
    debug_assert!(lx > 0 && ly > 0);
    if lx.min(ly) <= SCHOOLBOOK {
        schoolbook(a, &mut ws.bufs.acc, out, from, x.coeffs, y.coeffs, None);
        return;
    }
    let bits = a.modulus().significant_bits() as usize;
    if let Some((s, rn)) = middle_size(bits, lx, ly, from, from + out.len()) {
        kronecker(a, ws, out, from, x, y, s, Some(rn));
    } else {
        kronecker(a, ws, out, from, x, y, slot_width(bits, lx.min(ly)), None);
    }
}

/// The slot width `s` and the size `rn` of a middle product (see [`mul_part`]) of polynomials
/// of `lx` and `ly > SCHOOLBOOK` coefficients modulo a number of `bits` bits, for the
/// coefficients `from..to`, if cheaper than a full product.
pub(crate) fn middle_size(
    bits: usize,
    lx: usize,
    ly: usize,
    from: usize,
    to: usize,
) -> Option<(usize, usize)> {
    if from == 0 || !low_level() {
        return None;
    }
    let (s, n) = middle_bits(bits, lx, ly, from, to);
    let rn = mpn::mulmod_bnm1_next_size(n.div_ceil(64));
    let full = ((lx + ly - 1) * s).div_ceil(64);
    // Cheaper (a product modulo 2^N - 1 costs about that of a full product of N/2 bits), and
    // valid for GMP: an + bn > rn/2.
    (rn * 5 < full * 4 && ((lx + ly) * s) / 64 > rn / 2).then_some((s, rn))
}

/// The slot width `s` of a middle product (see [`mul_part`]) and the fewest bits `N` of its
/// product modulo `2^N - 1`: `N` covers the operands and the coefficients `from..to`, and the
/// part above `N` (below `2^(len*s - 1 - N)`, `len = lx + ly - 1`) wraps around below
/// `2^(from*s - 2)`, so that adding it to the coefficients below `from` (below
/// `2^(from*s - 1)`) doesn't carry into the coefficient `from`.
fn middle_bits(bits: usize, lx: usize, ly: usize, from: usize, to: usize) -> (usize, usize) {
    // With a spare bit: every coefficient is below 2^(s - 1).
    let s = slot_width(bits, lx.min(ly)) + 1;
    let len = lx + ly - 1;
    let n = (to.max(lx).max(ly) * s).max((len.saturating_sub(from)) * s + 1);
    (s, n)
}

/// Wrap-around product: `out[t]` = coefficient `from + t` of `x*y mod (X^L - 1)`, for some
/// `L >= lmin` (returned). Requires `x.len()`, `y.len()` and `from + out.len()` at most
/// `lmin`. Coefficient `k < L` is `c_k + c_(k + L)` for the coefficients `c` of `x*y` (which
/// has `x.len() + y.len() - 1` of them): exact when `k + L` is above the degree.
pub fn mul_wrap<A: PolyArith>(
    a: &A,
    ws: &mut Workspace,
    out: &mut [A::Elem],
    from: usize,
    x: Operand<'_, '_, A::Elem>,
    y: Operand<'_, '_, A::Elem>,
    lmin: usize,
) -> usize {
    let (lx, ly) = (x.len(), y.len());
    debug_assert!(lx > 0 && ly > 0 && lx.max(ly) <= lmin && from + out.len() <= lmin);
    if lx.min(ly) <= SCHOOLBOOK {
        schoolbook(
            a,
            &mut ws.bufs.acc,
            out,
            from,
            x.coeffs,
            y.coeffs,
            Some(lmin),
        );
        return lmin;
    }
    let bits = a.modulus().significant_bits() as usize;
    let key = (bits, lx, ly, lmin);
    match *ws
        .shapes
        .entry(key)
        .or_insert_with(|| wrap_size(bits, lx, ly, lmin))
    {
        Some(shape) => {
            kronecker(a, ws, out, from, x, y, shape.s, Some(shape.rn));
            shape.l
        }
        None => {
            // Full product: exact, the same as a wrap-around at L >= its length.
            mul_part(a, ws, out, from, x, y);
            lmin.max(lx + ly - 1)
        }
    }
}

/// The wrap-around product of [`mul_wrap`] for polynomials of `lx` and `ly > SCHOOLBOOK`
/// coefficients modulo a number of `bits` bits, if cheaper than a full product.
pub(crate) fn wrap_size(bits: usize, lx: usize, ly: usize, lmin: usize) -> Option<Shape> {
    // A wrapped coefficient is a sum of at most min(lx, ly) products (both are <= L), as in
    // a full product.
    let smin = slot_width(bits, lx.min(ly));
    let shape = wrap_shape(lmin, smin)?;
    let full = ((lx + ly) * smin).div_ceil(64);
    // Cheaper, and valid for GMP: an + bn > rn/2.
    (shape.rn * 5 < full * 4 && ((lx + ly) * shape.s) / 64 > shape.rn / 2).then_some(shape)
}

/// Size of a wrap-around product modulo `X^L - 1` and `2^(64*rn) - 1`, with slots of `s`
/// bits: `L*s = 64*rn`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Shape {
    pub l: usize,
    pub s: usize,
    pub rn: usize,
}

/// The smallest wrap-around product with `L >= lmin` and slots of `s >= smin` bits, among the
/// sizes `rn` efficient for GMP (`None` without GMP's low-level functions).
fn wrap_shape(lmin: usize, smin: usize) -> Option<Shape> {
    if !low_level() {
        return None;
    }
    let mut rn = mpn::mulmod_bnm1_next_size((lmin * smin).div_ceil(64));
    for _ in 0..200 {
        let bits = 64 * rn;
        // L = bits/s >= lmin: s <= bits/lmin.
        if let Some(s) = (smin..=bits / lmin).find(|s| bits.is_multiple_of(*s)) {
            return Some(Shape { l: bits / s, s, rn });
        }
        rn = mpn::mulmod_bnm1_next_size(rn + 1);
    }
    None
}

/// Schoolbook product: `out[t]` = coefficient `from + t` of `x*y`, or of `x*y mod (X^L - 1)`
/// with `wrap = Some(L)` (`x.len()`, `y.len() <= L`), accumulated on limbs, then reduced once.
fn schoolbook<A: PolyArith>(
    a: &A,
    acc: &mut Vec<u64>,
    out: &mut [A::Elem],
    from: usize,
    x: &[A::Elem],
    y: &[A::Elem],
    wrap: Option<usize>,
) {
    let width = 2 * a.limbs() + 2;
    acc.resize(width, 0);
    let (lx, ly) = (x.len(), y.len());
    for (t, r) in out.iter_mut().enumerate() {
        acc.fill(0);
        // At most one term per i: at most min(lx, ly) terms, wrapped or not.
        for k in std::iter::once(from + t).chain(wrap.map(|l| from + t + l)) {
            if k + 1 >= lx + ly {
                continue;
            }
            let lo = k.saturating_sub(ly - 1);
            let hi = k.min(lx - 1);
            for i in lo..=hi {
                a.mul_acc(acc, &x[i], &y[k - i]);
            }
        }
        a.redc_wide(r, acc);
    }
}

/// Bits per coefficient in a Kronecker product where the shorter factor has `len`
/// coefficients: each coefficient of the product is a sum of at most `len` products of values
/// `< n`
/// (`n` of `bits` bits).
pub(crate) fn slot_width(bits: usize, len: usize) -> usize {
    2 * bits + (usize::BITS - len.leading_zeros()) as usize
}

/// Packs the values of `x` into `buf`, one every `s` bits. Returns the number of limbs they
/// take (`buf` has a few more, all zero).
fn pack<A: PolyArith>(a: &A, buf: &mut Vec<u64>, x: &[A::Elem], s: usize) -> usize {
    let limbs = a.limbs();
    let size = (x.len() * s).div_ceil(64);
    buf.clear();
    buf.resize(size + limbs + 1, 0);
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
    size
}

/// Kronecker substitution: coefficients `from..from + out.len()` of `x*y` with slots of `s`
/// bits, with one big integer product, or modulo `2^(64*rn) - 1` with `rn = Some(rn)` (a
/// wrap-around product modulo `X^L - 1` if `s*L = 64*rn`, see [`mul_wrap`], or a middle product,
/// see [`mul_part`]).
#[allow(clippy::too_many_arguments)]
fn kronecker<A: PolyArith>(
    a: &A,
    ws: &mut Workspace,
    out: &mut [A::Elem],
    from: usize,
    x: Operand<'_, '_, A::Elem>,
    y: Operand<'_, '_, A::Elem>,
    s: usize,
    rn: Option<usize>,
) {
    let Workspace {
        x: bx, y: by, bufs, ..
    } = ws;
    let xl = match x.cache {
        Some(cache) => cache.get(a, x.coeffs, s),
        None => {
            let size = pack(a, bx, x.coeffs, s);
            &bx[..size]
        }
    };
    let yl = match y.cache {
        Some(cache) => cache.get(a, y.coeffs, s),
        None => {
            let size = pack(a, by, y.coeffs, s);
            &by[..size]
        }
    };
    let (xl, yl) = if xl.len() >= yl.len() {
        (xl, yl)
    } else {
        (yl, xl)
    };

    // The product, with room to read every requested coefficient.
    let width = s.div_ceil(64);
    let len = rn.unwrap_or(xl.len() + yl.len());
    let room = ((from + out.len()) * s).div_ceil(64).max(len) + width + 1;
    let product = &mut bufs.product;
    product.clear();
    product.resize(room, 0);
    match rn {
        Some(rn) => mpn::mulmod_bnm1(&mut product[..rn], xl, yl, &mut bufs.scratch),
        None if low_level() => mpn::mul_long(&mut product[..len], xl, yl),
        None => {
            let [ix, iy, ip] = &mut bufs.ints;
            ix.assign_digits(xl, Order::Lsf);
            iy.assign_digits(yl, Order::Lsf);
            ip.assign(&*ix * &*iy);
            let size = ip.significant_digits::<u64>();
            ip.write_digits(&mut product[..size], Order::Lsf);
        }
    }

    // Unpack: coefficient k is the s bits at k*s.
    let top_mask = if s.is_multiple_of(64) {
        u64::MAX
    } else {
        (1 << (s % 64)) - 1
    };
    let t = &mut bufs.acc;
    t.resize(width, 0);
    for (k, r) in (from..).zip(out.iter_mut()) {
        let (word, shift) = ((k * s) / 64, (k * s) % 64);
        if shift == 0 {
            t.copy_from_slice(&product[word..word + width]);
        } else {
            for j in 0..width {
                t[j] = (product[word + j] >> shift) | (product[word + j + 1] << (64 - shift));
            }
        }
        t[width - 1] &= top_mask;
        // The value is below 2^smin: at most 2*limbs + 1 limbs, even if the slot is wider.
        a.redc_wide(r, &t[..width.min(2 * a.limbs() + 2)]);
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
        let h = next - prec;
        // f*g = 1 + X^prec * e (mod X^next): a middle product.
        e.resize(h, a.zero());
        let x = Operand::new(&f[..fl]);
        mul_part(a, ws, &mut e, prec, x, Operand::new(&g[..prec]));
        // g[prec..next] = -(g*e)[..h].
        d.resize(h, a.zero());
        mul_part(a, ws, &mut d, 0, Operand::new(&g[..h]), Operand::new(&e));
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

/// A monic modulus `f` of degree `d`, with what the reductions modulo `f` need: the inverse of
/// its reverse modulo `X^d`, and both packed for the products that reuse them.
pub struct Modulus<E> {
    /// The `d + 1` coefficients of `f`, the leading 1 included.
    f: Vec<E>,
    /// Inverse of the reverse of `f` modulo `X^d`.
    inv: Vec<E>,
    packed_f: Packed,
    packed_inv: Packed,
}

impl<E: Clone> Modulus<E> {
    /// The modulus `f` (monic, the leading 1 implicit, degree `f.len() >= 1`).
    pub fn new<A: PolyArith<Elem = E>>(a: &A, ws: &mut Workspace, f: &[E]) -> Self {
        let d = f.len();
        let inv = inverse(a, ws, &reverse_monic(a, f, d), d);
        let mut full = f.to_vec();
        full.push(a.poly_from(&Integer::from(1)));
        Modulus {
            f: full,
            inv,
            packed_f: Packed::default(),
            packed_inv: Packed::default(),
        }
    }

    /// Degree of `f`.
    pub fn degree(&self) -> usize {
        self.f.len() - 1
    }
}

/// `h = h*g mod f`, with `f` monic of degree `d = h.len()` (so `h` of degree `< d`) and `g` of
/// degree `< d` (`g.len() <= d`).
///
/// Quotient `q` from the high coefficients of `p = h*g` and the inverse of the reverse of `f`
/// (a short product), then `p - q*f`: its degree is `< d` and the coefficients of `q*f` of
/// degree `>= d` are those of `p`, so a wrap-around product modulo `X^L - 1` (`L > d`) is
/// enough.
pub fn mul_mod<A: PolyArith>(
    a: &A,
    ws: &mut Workspace,
    h: &mut [A::Elem],
    g: &[A::Elem],
    m: &mut Modulus<A::Elem>,
) {
    let d = m.degree();
    debug_assert!(h.len() == d && !g.is_empty() && g.len() <= d);
    let len = d + g.len() - 1;
    let mut p = vec![a.zero(); len];
    mul(a, ws, &mut p, h, g);
    if len <= d {
        h[..len].clone_from_slice(&p);
        h[len..].fill(a.zero());
        return;
    }
    // Quotient, t coefficients: reverse(q) = reverse(p) / reverse(f) mod X^t.
    let t = len - d;
    let rev_p: Vec<A::Elem> = p[d..].iter().rev().cloned().collect();
    let mut rev_q = vec![a.zero(); t];
    let inv = Operand::cached(&m.inv[..t], &mut m.packed_inv);
    mul_part(a, ws, &mut rev_q, 0, Operand::new(&rev_p), inv);
    let q: Vec<A::Elem> = rev_q.into_iter().rev().collect();
    // Remainder: coefficient k < d of p - q*f, with (q*f mod (X^L - 1))_k = c_k + c_(k + L)
    // and c_(k + L) = p_(k + L) (degree >= d).
    let mut w = vec![a.zero(); d];
    let f = Operand::cached(&m.f, &mut m.packed_f);
    let l = mul_wrap(a, ws, &mut w, 0, f, Operand::new(&q), d + 1);
    let mut c = a.zero();
    for (k, h) in h.iter_mut().enumerate() {
        match p.get(k + l) {
            Some(high) => a.sub(&mut c, &w[k], high),
            None => c.clone_from(&w[k]),
        }
        a.sub(h, &p[k], &c);
    }
}

/// Values of `h` (degree `< d`) at the `d` roots of the product tree `tree` (in the order of the
/// leaves), whose root is the modulus `m`.
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
    m: &mut Modulus<A::Elem>,
) -> Vec<A::Elem> {
    let d = h.len();
    debug_assert_eq!(d, m.degree());
    // Root: h/f = X^-1 * rev(h)(1/X) / rev(f)(1/X), so c_(t+1) = (rev(h) * inv)[t].
    let rev_h: Vec<A::Elem> = h.iter().rev().cloned().collect();
    let mut series = vec![a.zero(); d];
    let inv = Operand::cached(&m.inv[..d], &mut m.packed_inv);
    mul_part(a, ws, &mut series, 0, Operand::new(&rev_h), inv);
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
        let acc = &mut ws.bufs.acc;
        acc.resize(2 * a.limbs() + 2, 0);
        for (k, out) in out.iter_mut().enumerate() {
            acc.fill(0);
            for (r, s) in r.iter().zip(&s[k..]) {
                a.mul_acc(acc, r, s);
            }
            a.redc_wide(t, acc);
            a.add(out, t, &s[k + m]);
        }
        return;
    }
    // (reverse(r) * s)[t + m - 1] = sum_i r_i * s[t + i]: a middle product.
    rev.clear();
    rev.extend(r.iter().rev().cloned());
    let s_used = &s[..l + m - 1];
    prod.resize(l, a.zero());
    let (x, y) = (Operand::new(&rev[..]), Operand::new(s_used));
    mul_part(a, ws, prod, m - 1, x, y);
    for (i, out) in out.iter_mut().enumerate() {
        a.add(out, &prod[i], &s[i + m]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arith::{with_arith, Arith};
    use rug::rand::RandState;
    use std::cell::Cell;

    thread_local! {
        /// Products without GMP's low-level functions in this thread (see [`low_level`]).
        pub(super) static NO_LOW_LEVEL: Cell<bool> = const { Cell::new(false) };
    }

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

            // Parts of the product, and wrap-around products (the same operands, packed once).
            let full = naive_mul(&xv, &yv, &n);
            let (mut px, mut py) = (Packed::default(), Packed::default());
            for _ in 0..3 {
                let from = below(lx + ly, rand);
                let count = 1 + below(lx + ly + 3 - from, rand);
                let mut out = vec![a.zero(); count];
                let (ox, oy) = (Operand::cached(&x, &mut px), Operand::cached(&y, &mut py));
                mul_part(a, &mut ws, &mut out, from, ox, oy);
                for (k, out) in (from..).zip(values(a, &out)) {
                    let expected = full.get(k).cloned().unwrap_or_default();
                    assert_eq!(out, expected, "{lx} {ly} part {k}");
                }
                let lmin = lx.max(ly) + below(lx + ly, rand);
                let from = below(lmin, rand);
                let mut out = vec![a.zero(); 1 + below(lmin - from, rand)];
                let (ox, oy) = (Operand::cached(&x, &mut px), Operand::new(&y));
                let l = mul_wrap(a, &mut ws, &mut out, from, ox, oy, lmin);
                assert!(l >= lmin);
                for (k, out) in (from..).zip(values(a, &out)) {
                    let expected = full.iter().skip(k).step_by(l).sum::<Integer>() % &n;
                    assert_eq!(out, expected, "{lx} {ly} wrap {k} mod X^{l} - 1");
                }
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
    fn wrap_shapes() {
        if !mpn::ENABLED {
            return;
        }
        for lmin in [
            13, 14, 100, 129, 257, 481, 1000, 1441, 1921, 2048, 2881, 6000, 11521,
        ] {
            for smin in [67, 130, 523, 1035, 2059, 2200] {
                let shape = wrap_shape(lmin, smin).unwrap();
                assert!(shape.l >= lmin && shape.s >= smin);
                assert_eq!(shape.l * shape.s, 64 * shape.rn);
                // Not much larger than needed, but for small products (rarely wrapped).
                let waste = (shape.rn * 64) as f64 / (lmin * smin) as f64;
                if lmin >= 1000 {
                    assert!(waste < 1.2, "{lmin} {smin} {shape:?}");
                }
            }
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

    #[test]
    fn products_without_low_level() {
        // The products (full ones only) and reductions without GMP's low-level functions.
        NO_LOW_LEVEL.set(true);
        assert!(middle_size(1024, 100, 100, 50, 150).is_none());
        assert!(wrap_size(1024, 100, 100, 150).is_none());
        let mut rand = RandState::new();
        for limbs in [1, 4, 11, 16u32] {
            let n = (Integer::from(1) << (64 * limbs)) - 1u32;
            with_arith!(&n, |a| {
                check_random_products(&a, &mut rand);
                check_tree(&a, &mut rand);
            });
        }
        NO_LOW_LEVEL.set(false);
    }

    /// Coefficients of `x*y` (values), or of `x*y mod (X^l - 1)` with `wrap = Some(l)`.
    fn naive_product(
        x: &[Integer],
        y: &[Integer],
        n: &Integer,
        wrap: Option<usize>,
    ) -> Vec<Integer> {
        let full = naive_mul(x, y, n);
        match wrap {
            None => full,
            Some(l) => (0..l)
                .map(|k| full.iter().skip(k).step_by(l).sum::<Integer>() % n)
                .collect(),
        }
    }

    /// Middle and wrap-around products at the tightest sizes: the fewest bits for a middle
    /// product (not rounded up by `mulmod_bnm1_next_size`), slots of exactly `slot_width` bits
    /// for a wrap-around product with `lx = ly = L`, the largest coefficients (all `n - 1`,
    /// for `n = 2^(64*limbs) - 1`: the largest sums in a slot), then the same by `mul_part` and
    /// `mul_wrap`.
    fn check_tight<A: PolyArith>(a: &A, rand: &mut RandState<'_>) {
        if !mpn::ENABLED {
            return;
        }
        let n = a.modulus().clone();
        let bits = n.significant_bits() as usize;
        let mut ws = Workspace::new();
        let limbs = |len: usize, s: usize| (len * s).div_ceil(64);
        for (lx, ly) in [
            (13, 13),
            (13, 14),
            (15, 13),
            (16, 16),
            (17, 40),
            (31, 32),
            (33, 31),
            (64, 64),
            (63, 65),
            (100, 29),
            (127, 128),
        ] {
            let len = lx + ly - 1;
            for kind in 0..3 {
                let (x, y) = match kind {
                    0 => (largest_poly(a, lx), largest_poly(a, ly)),
                    1 => (random_poly(a, lx, rand), largest_poly(a, ly)),
                    _ => (random_poly(a, lx, rand), random_poly(a, ly, rand)),
                };
                let (xv, yv) = (values(a, &x), values(a, &y));
                let full = naive_product(&xv, &yv, &n, None);
                let min = lx.min(ly);
                for from in [
                    1,
                    2,
                    min - 1,
                    min,
                    min + 1,
                    len / 2,
                    len - min,
                    len - 2,
                    len - 1,
                ] {
                    for to in [from + 1, from + min, len - 1, len, len + 2] {
                        if to <= from {
                            continue;
                        }
                        let (s, nbits) = middle_bits(bits, lx, ly, from, to);
                        let rn = nbits.div_ceil(64);
                        // GMP's requirements (an + bn > rn/2, larger operand at most rn).
                        if limbs(lx, s) + limbs(ly, s) <= rn / 2 || limbs(lx.max(ly), s) > rn {
                            continue;
                        }
                        let mut out = vec![a.zero(); to - from];
                        let (ox, oy) = (Operand::new(&x[..]), Operand::new(&y[..]));
                        kronecker(a, &mut ws, &mut out, from, ox, oy, s, Some(rn));
                        for (k, out) in (from..).zip(values(a, &out)) {
                            let expected = full.get(k).cloned().unwrap_or_default();
                            assert_eq!(out, expected, "{bits}: {lx} {ly} {from}..{to} rn = {rn}");
                        }
                        let (ox, oy) = (Operand::new(&x[..]), Operand::new(&y[..]));
                        mul_part(a, &mut ws, &mut out, from, ox, oy);
                        for (k, out) in (from..).zip(values(a, &out)) {
                            let expected = full.get(k).cloned().unwrap_or_default();
                            assert_eq!(out, expected, "{bits}: {lx} {ly} part {from}..{to}");
                        }
                    }
                }

                // Wrap-around, L = lx = ly: min(lx, ly) terms in every coefficient.
                let l = lx.max(ly);
                let (x, y, xv, yv) = if lx == ly {
                    (x, y, xv, yv)
                } else {
                    let (x, y) = (largest_poly(a, l), largest_poly(a, l));
                    let (xv, yv) = (values(a, &x), values(a, &y));
                    (x, y, xv, yv)
                };
                let wrapped = naive_product(&xv, &yv, &n, Some(l));
                let s = slot_width(bits, l);
                // L*s bits must be whole limbs: L a multiple of 64/gcd(s, 64).
                if (l * s).is_multiple_of(64) {
                    let mut out = vec![a.zero(); l];
                    let (ox, oy) = (Operand::new(&x[..]), Operand::new(&y[..]));
                    kronecker(a, &mut ws, &mut out, 0, ox, oy, s, Some(l * s / 64));
                    assert_eq!(values(a, &out), wrapped, "{bits}: wrap {l} slots of {s}");
                }
                let mut out = vec![a.zero(); l];
                let (ox, oy) = (Operand::new(&x[..]), Operand::new(&y[..]));
                let lw = mul_wrap(a, &mut ws, &mut out, 0, ox, oy, l);
                let expected = naive_product(&xv, &yv, &n, Some(lw));
                assert_eq!(values(a, &out)[..], expected[..l], "{bits}: {l}");
            }
        }
        // Wrap-around products with slots of exactly slot_width bits: the smallest L >= 13
        // with L*s a multiple of 64.
        for m in [13, 15, 31, 63, 127] {
            let s = slot_width(bits, m);
            let l = (m..).find(|l| (l * s).is_multiple_of(64)).unwrap();
            if l > 400 {
                continue;
            }
            let (x, y) = (largest_poly(a, l), largest_poly(a, m));
            let wrapped = naive_product(&values(a, &x), &values(a, &y), &n, Some(l));
            let mut out = vec![a.zero(); l];
            let (ox, oy) = (Operand::new(&x[..]), Operand::new(&y[..]));
            kronecker(a, &mut ws, &mut out, 0, ox, oy, s, Some(l * s / 64));
            assert_eq!(
                values(a, &out),
                wrapped,
                "{bits}: wrap {l} x {m}, slots of {s}"
            );
        }
    }

    #[test]
    fn tight_products() {
        let mut rand = RandState::new();
        for limbs in 1..=16u32 {
            let n = (Integer::from(1) << (64 * limbs)) - 1u32;
            with_arith!(&n, |a| check_tight(&a, &mut rand));
        }
        // Odd sizes, and Plain.
        for n in [
            (Integer::from(1) << 100) - 3u32,
            (Integer::from(1) << 700) - 1u32,
            (Integer::from(1) << 1100) - 1u32,
            (Integer::from(1) << 256) - 2u32,
        ] {
            with_arith!(&n, |a| check_tight(&a, &mut rand));
        }
    }

    #[test]
    fn product_sizes() {
        // The sizes of the middle and wrap-around products satisfy the bounds and GMP's
        // requirements, including where mulmod_bnm1_next_size rounds up (FFT sizes).
        if !mpn::ENABLED {
            return;
        }
        let limbs = |len: usize, s: usize| (len * s).div_ceil(64);
        let mut lens: Vec<usize> = (13..80).collect();
        for p in 6..16 {
            lens.extend([
                (1 << p) - 1,
                1 << p,
                (1 << p) + 1,
                3 << (p - 1),
                5 << (p - 2),
            ]);
        }
        lens.extend([1440, 1920, 2880, 5760, 11520, 23040]);
        for bits in [2, 64, 65, 128, 300, 640, 641, 1024, 1100] {
            for &lx in &lens {
                for &ly in &[13, 14, 64, 100, lx / 2 + 13, lx, lx + 1] {
                    let (min, len) = (lx.min(ly), lx + ly - 1);
                    let smin = slot_width(bits, min);
                    // A coefficient is a sum of at most min products of values < 2^bits.
                    let max = (Integer::from(1) << bits as u32) - 1u32;
                    assert!((max.square() * min) >> smin as u32 == 0);
                    for from in [1, min - 1, len / 3, len - min, len - 1] {
                        for to in [from + 1, len.min(from + lx), len] {
                            if let Some((s, rn)) = middle_size(bits, lx, ly, from, to) {
                                assert!(s > smin && 64 * rn >= to.max(lx).max(ly) * s);
                                assert!(64 * rn > (len - from) * s);
                                assert!(limbs(lx, s) + limbs(ly, s) > rn / 2);
                                assert!(limbs(lx.max(ly), s) <= rn);
                            }
                        }
                    }
                    for lmin in [lx.max(ly), lx.max(ly) + 1, len - 1, len] {
                        if let Some(shape) = wrap_size(bits, lx, ly, lmin) {
                            assert!(shape.s >= smin && shape.l >= lmin);
                            assert_eq!(shape.l * shape.s, 64 * shape.rn);
                            assert!(limbs(lx, shape.s) + limbs(ly, shape.s) > shape.rn / 2);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn large_middle_and_wrap_products() {
        // Sizes where GMP multiplies by FFT and mulmod_bnm1_next_size rounds up: middle and
        // wrap-around products against the full product.
        let mut rand = RandState::new();
        for (bits, lens) in [(64, &[1000, 2100, 4099][..]), (1024, &[700, 1921][..])] {
            let n = (Integer::from(1) << bits) - 1u32;
            with_arith!(&n, |a| {
                let mut ws = Workspace::new();
                for &l in lens {
                    for (lx, ly) in [(l, l), (l, l / 2 + 7), (l / 3 + 1, l)] {
                        let (x, y) = if l % 2 == 0 {
                            (largest_poly(&a, lx), largest_poly(&a, ly))
                        } else {
                            (random_poly(&a, lx, &mut rand), largest_poly(&a, ly))
                        };
                        let len = lx + ly - 1;
                        let mut full = vec![a.zero(); len];
                        mul(&a, &mut ws, &mut full, &x, &y);
                        let full = values(&a, &full);
                        for from in [1, lx.min(ly) - 1, len / 2, len - lx.min(ly)] {
                            let mut out = vec![a.zero(); lx.max(ly).min(len - from)];
                            let (ox, oy) = (Operand::new(&x[..]), Operand::new(&y[..]));
                            mul_part(&a, &mut ws, &mut out, from, ox, oy);
                            assert_eq!(values(&a, &out)[..], full[from..from + out.len()]);
                        }
                        let lmin = lx.max(ly) + 1;
                        let mut out = vec![a.zero(); lmin];
                        let (ox, oy) = (Operand::new(&x[..]), Operand::new(&y[..]));
                        let lw = mul_wrap(&a, &mut ws, &mut out, 0, ox, oy, lmin);
                        for (k, out) in values(&a, &out).into_iter().enumerate() {
                            let expected = full.iter().skip(k).step_by(lw).sum::<Integer>() % &n;
                            assert_eq!(out, expected, "{bits}: {lx} {ly} wrap {k} mod X^{lw} - 1");
                        }
                    }
                }
            });
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

            // h*g mod f (several times, reusing the packed f and inverse; with g shorter, and
            // with the largest values), then evaluation at the roots.
            let mut modulus = Modulus::new(a, &mut ws, tree.root());
            assert_eq!(values(a, &modulus.inv), values(a, &inv));
            let mut h = random_poly(a, d, rand);
            let mut expected_h = values(a, &h);
            for (i, gl) in [d, d, 1 + d / 3, d, d].into_iter().enumerate() {
                let g = if i == 3 {
                    largest_poly(a, gl)
                } else {
                    random_poly(a, gl, rand)
                };
                if i == 3 {
                    h = largest_poly(a, d);
                    expected_h = values(a, &h);
                }
                mul_mod(a, &mut ws, &mut h, &g, &mut modulus);
                expected_h = naive_rem(&naive_mul(&expected_h, &values(a, &g), &n), &f, &n);
                assert_eq!(values(a, &h), expected_h, "d = {d}, {i}");
            }
            let got = evaluate(a, &mut ws, &h, &tree, &mut modulus);
            for (r, got) in roots.iter().zip(values(a, &got)) {
                assert_eq!(got, eval(&expected_h, r, &n), "d = {d}");
            }
        }
    }

    /// Remainder of `x` modulo the monic `f` (all its coefficients, leading 1 included).
    fn naive_rem(x: &[Integer], f: &[Integer], n: &Integer) -> Vec<Integer> {
        let d = f.len() - 1;
        let mut x = x.to_vec();
        for k in (d..x.len()).rev() {
            let q = x[k].clone();
            for (i, f) in f.iter().enumerate() {
                x[k - d + i] = Integer::from(&x[k - d + i] - &q * f) % n;
            }
        }
        x.resize(d, Integer::new());
        x.iter().map(|x| (Integer::from(x % n) + n) % n).collect()
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
