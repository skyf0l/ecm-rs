//! Montgomery's block Lanczos algorithm over GF(2) ("A block Lanczos algorithm for finding
//! dependencies over GF(2)", EUROCRYPT 1995), with blocks of 64 vectors: vectors in the null
//! space of a sparse matrix `B` with more columns than rows, from the iteration on the
//! symmetric `A = B^T B`.

use crate::stop::Stop;

/// A 64x64 matrix over GF(2): row `i` is `m[i]`, column `j` its bit `j`.
type M64 = [u64; 64];

/// A sparse matrix over GF(2) by columns: the rows of the nonzero entries of each column.
pub(crate) struct Sparse<'a> {
    pub rows: usize,
    pub cols: &'a [Vec<u32>],
}

impl Sparse<'_> {
    /// `B v` (one vector of 64 bits per row).
    fn mul(&self, v: &[u64], out: &mut [u64]) {
        out.fill(0);
        for (col, &x) in self.cols.iter().zip(v) {
            if x != 0 {
                for &r in col {
                    out[r as usize] ^= x;
                }
            }
        }
    }

    /// `B^T w`.
    fn mul_transpose(&self, w: &[u64], out: &mut [u64]) {
        for (col, o) in self.cols.iter().zip(out.iter_mut()) {
            *o = col.iter().fold(0, |acc, &r| acc ^ w[r as usize]);
        }
    }

    /// `A v = B^T B v`, with `tmp` of one entry per row.
    fn mul_a(&self, v: &[u64], tmp: &mut [u64], out: &mut [u64]) {
        self.mul(v, tmp);
        self.mul_transpose(tmp, out);
    }
}

/// `V^T W` for two blocks of vectors.
fn inner(v: &[u64], w: &[u64]) -> M64 {
    // By byte of v: table[byte position][byte value] accumulates the w with this byte.
    let mut table = vec![[0u64; 256]; 8];
    for (&x, &y) in v.iter().zip(w) {
        for (b, t) in table.iter_mut().enumerate() {
            t[((x >> (8 * b)) & 0xff) as usize] ^= y;
        }
    }
    let mut out = [0u64; 64];
    for (b, t) in table.iter().enumerate() {
        for bit in 0..8 {
            let mut acc = 0;
            for (value, &entry) in t.iter().enumerate() {
                if value >> bit & 1 == 1 {
                    acc ^= entry;
                }
            }
            out[8 * b + bit] = acc;
        }
    }
    out
}

/// Tables of the products of the bytes of a vector by `m`: `t[b][x]` for the byte `b` of value
/// `x`.
fn byte_tables(m: &M64) -> Vec<[u64; 256]> {
    let mut tables = vec![[0u64; 256]; 8];
    for (b, t) in tables.iter_mut().enumerate() {
        for x in 1..256usize {
            let low = x & (x - 1);
            let bit = x.trailing_zeros() as usize;
            t[x] = t[low] ^ m[8 * b + bit];
        }
    }
    tables
}

/// `out ^= V M`.
fn mul_add(v: &[u64], m: &M64, out: &mut [u64]) {
    let tables = byte_tables(m);
    for (o, &x) in out.iter_mut().zip(v) {
        let mut acc = 0;
        for (b, t) in tables.iter().enumerate() {
            acc ^= t[((x >> (8 * b)) & 0xff) as usize];
        }
        *o ^= acc;
    }
}

/// `A B` of 64x64 matrices.
fn mul_mm(a: &M64, b: &M64) -> M64 {
    let mut out = [0u64; 64];
    for (o, &row) in out.iter_mut().zip(a) {
        let mut acc = 0;
        let mut bits = row;
        while bits != 0 {
            acc ^= b[bits.trailing_zeros() as usize];
            bits &= bits - 1;
        }
        *o = acc;
    }
    out
}

/// `M S S^T`: the columns of `m` out of the mask `s` cleared.
fn mask_cols(m: &M64, s: u64) -> M64 {
    m.map(|row| row & s)
}

fn identity() -> M64 {
    std::array::from_fn(|i| 1u64 << i)
}

fn add(a: &M64, b: &M64) -> M64 {
    std::array::from_fn(|i| a[i] ^ b[i])
}

/// Montgomery's choice of the columns `S_i` and of `W_i^-1 = S (S^T T S)^-1 S^T` for
/// `T = V_i^T A V_i`, with priority to the columns not in `last` (`S_{i-1}`).
pub(crate) fn choose_s(t: &M64, last: u64) -> (M64, u64) {
    let mut m: [[u64; 2]; 64] = std::array::from_fn(|i| [t[i], 1u64 << i]);
    let order: Vec<usize> = (0..64)
        .filter(|&i| last >> i & 1 == 0)
        .chain((0..64).filter(|&i| last >> i & 1 == 1))
        .collect();
    let mut s = 0u64;
    for j in 0..64 {
        let c = order[j];
        let bit = 1u64 << c;
        if let Some(k) = (j..64).find(|&k| m[order[k]][0] & bit != 0) {
            m.swap(order[j], order[k]);
            s |= bit;
            let pivot = m[c];
            for (r, row) in m.iter_mut().enumerate() {
                if r != c && row[0] & bit != 0 {
                    row[0] ^= pivot[0];
                    row[1] ^= pivot[1];
                }
            }
        } else {
            let k = (j..64)
                .find(|&k| m[order[k]][1] & bit != 0)
                .expect("an invertible matrix");
            m.swap(order[j], order[k]);
            let pivot = m[c];
            for (r, row) in m.iter_mut().enumerate() {
                if r != c && row[1] & bit != 0 {
                    row[0] ^= pivot[0];
                    row[1] ^= pivot[1];
                }
            }
            m[c] = [0, 0];
        }
    }
    (std::array::from_fn(|i| m[i][1]), s)
}

/// Vectors of the null space of `b` (bit `k` of entry `j`: coordinate `j` of the vector `k`),
/// from the random start `seed`; `Err` if `stop` was requested. Some vectors may be zero.
pub(crate) fn null_space(b: &Sparse<'_>, seed: u64, stop: Stop<'_>) -> Result<Vec<u64>, ()> {
    let n = b.cols.len();
    let mut state = seed;
    let y: Vec<u64> = (0..n).map(|_| super::next_random(&mut state)).collect();
    let mut tmp = vec![0u64; b.rows];
    let mut v0 = vec![0u64; n];
    b.mul_a(&y, &mut tmp, &mut v0);
    let rhs = v0.clone();

    let mut x = vec![0u64; n];
    let mut v1 = vec![0u64; n];
    let mut v2 = vec![0u64; n];
    let mut av = vec![0u64; n];
    let (mut winv1, mut winv2) = ([0u64; 64], [0u64; 64]);
    let (mut vtav1, mut vta2v1) = ([0u64; 64], [0u64; 64]);
    let mut s1 = u64::MAX;
    // At most n / 63 iterations are expected: more is a failure.
    let max_iterations = n / 32 + 100;
    let mut iterations = 0;
    loop {
        if stop.requested() {
            return Err(());
        }
        iterations += 1;
        if iterations > max_iterations {
            break;
        }
        b.mul_a(&v0, &mut tmp, &mut av);
        let vtav = inner(&v0, &av);
        if vtav.iter().all(|&r| r == 0) {
            break;
        }
        let vta2v = inner(&av, &av);
        let (winv, s) = choose_s(&vtav, s1);
        if s == 0 {
            break;
        }
        // x += V Winv V^T b.
        let vtb = inner(&v0, &rhs);
        mul_add(&v0, &mul_mm(&winv, &vtb), &mut x);

        // D = I - Winv (V^T A^2 V S S^T + V^T A V).
        let d = add(
            &identity(),
            &mul_mm(&winv, &add(&mask_cols(&vta2v, s), &vtav)),
        );
        // E = -Winv_{i-1} V^T A V S S^T.
        let e = mul_mm(&winv1, &mask_cols(&vtav, s));
        // F = -Winv_{i-2} (I - V_{i-1}^T A V_{i-1} Winv_{i-1})
        //     (V_{i-1}^T A^2 V_{i-1} S_{i-1} S_{i-1}^T + V_{i-1}^T A V_{i-1}) S S^T.
        let f = mask_cols(
            &mul_mm(
                &winv2,
                &mul_mm(
                    &add(&identity(), &mul_mm(&vtav1, &winv1)),
                    &add(&mask_cols(&vta2v1, s1), &vtav1),
                ),
            ),
            s,
        );
        // V_{i+1} = A V S S^T + V D + V_{i-1} E + V_{i-2} F.
        let mut next: Vec<u64> = av.iter().map(|&a| a & s).collect();
        mul_add(&v0, &d, &mut next);
        mul_add(&v1, &e, &mut next);
        mul_add(&v2, &f, &mut next);

        v2 = std::mem::replace(&mut v1, std::mem::replace(&mut v0, next));
        (winv2, winv1) = (winv1, winv);
        (vtav1, vta2v1, s1) = (vtav, vta2v, s);
    }

    // Combinations of the columns of [x - y, V_m] in the null space of B.
    let z: Vec<(u64, u64)> = x
        .iter()
        .zip(&y)
        .zip(&v0)
        .map(|((&x, &y), &v)| (x ^ y, v))
        .collect();
    let mut bz = vec![0u128; b.rows];
    for (col, &(z0, z1)) in b.cols.iter().zip(&z) {
        let value = u128::from(z0) | u128::from(z1) << 64;
        if value != 0 {
            for &r in col {
                bz[r as usize] ^= value;
            }
        }
    }
    let kernel = kernel_128(&bz);
    let mut out = vec![0u64; n];
    let mut count = 0;
    for c in kernel {
        if count == 64 {
            break;
        }
        let (c0, c1) = (c as u64, (c >> 64) as u64);
        let mut any = false;
        let bit = 1u64 << count;
        for (o, &(z0, z1)) in out.iter_mut().zip(&z) {
            if ((z0 & c0).count_ones() + (z1 & c1).count_ones()) % 2 == 1 {
                *o |= bit;
                any = true;
            }
        }
        if any {
            count += 1;
        }
    }
    Ok(out)
}

/// A basis of the kernel of the map `c -> (parity(row & c))` for the rows given.
fn kernel_128(rows: &[u128]) -> Vec<u128> {
    // Echelon basis of the row space, by pivot bit.
    let mut basis: Vec<u128> = Vec::new();
    for &row in rows {
        let mut r = row;
        for &b in &basis {
            let pivot = 1u128 << b.trailing_zeros();
            if r & pivot != 0 {
                r ^= b;
            }
        }
        if r != 0 {
            // Keep the basis reduced: clear the new pivot from the others.
            let pivot = 1u128 << r.trailing_zeros();
            for b in &mut basis {
                if *b & pivot != 0 {
                    *b ^= r;
                }
            }
            basis.push(r);
        }
    }
    let pivots: u128 = basis
        .iter()
        .fold(0, |acc, b| acc | 1u128 << b.trailing_zeros());
    (0..128)
        .filter(|&f| pivots >> f & 1 == 0)
        .map(|f| {
            let mut c = 1u128 << f;
            for &b in &basis {
                if b >> f & 1 == 1 {
                    c |= 1u128 << b.trailing_zeros();
                }
            }
            c
        })
        .collect()
}

#[cfg(test)]
#[allow(clippy::needless_range_loop, reason = "matrix indices")]
mod tests {
    use super::*;

    fn random_m64(state: &mut u64) -> M64 {
        std::array::from_fn(|_| super::super::next_random(state))
    }

    #[test]
    fn dense_products() {
        let mut state = 1;
        let a = random_m64(&mut state);
        let b = random_m64(&mut state);
        let v: Vec<u64> = (0..200)
            .map(|_| super::super::next_random(&mut state))
            .collect();
        let w: Vec<u64> = (0..200)
            .map(|_| super::super::next_random(&mut state))
            .collect();
        // (V A) B = V (A B).
        let mut va = vec![0; 200];
        mul_add(&v, &a, &mut va);
        let mut vab = vec![0; 200];
        mul_add(&va, &b, &mut vab);
        let mut v_ab = vec![0; 200];
        mul_add(&v, &mul_mm(&a, &b), &mut v_ab);
        assert_eq!(vab, v_ab);
        // (V A)^T W = A^T (V^T W): checked entry by entry.
        let t = inner(&v, &w);
        for i in 0..64 {
            for j in 0..64 {
                let expected = v
                    .iter()
                    .zip(&w)
                    .filter(|&(x, y)| x >> i & 1 == 1 && y >> j & 1 == 1)
                    .count()
                    % 2;
                assert_eq!((t[i] >> j & 1) as usize, expected);
            }
        }
    }

    #[test]
    fn choose_s_inverts() {
        let mut state = 7;
        for rank_mask in [u64::MAX, 0x0000_ffff_ffff_ffff, 0x5555_5555_5555_5555] {
            // A symmetric T = U^T U of rank at most the popcount of the mask.
            let u = random_m64(&mut state).map(|r| r & rank_mask);
            let mut ut = [0u64; 64];
            for i in 0..64 {
                for j in 0..64 {
                    if u[j] >> i & 1 == 1 {
                        ut[i] |= 1 << j;
                    }
                }
            }
            let t = mul_mm(&ut, &u);
            let (winv, s) = choose_s(&t, 0);
            // Winv T restricted to S is the identity on S.
            let prod = mul_mm(&winv, &t);
            for i in 0..64 {
                if s >> i & 1 == 1 {
                    assert_eq!(prod[i] & s, 1 << i, "row {i}");
                } else {
                    assert_eq!(winv[i], 0);
                }
            }
            assert!(s.count_ones() >= rank_mask.count_ones().min(40) / 2);
        }
    }

    #[test]
    fn kernel() {
        let mut state = 3;
        let rows: Vec<u128> = (0..300)
            .map(|_| {
                let r = u128::from(super::super::next_random(&mut state))
                    | u128::from(super::super::next_random(&mut state)) << 64;
                // Rank at most 100.
                r & ((1u128 << 100) - 1)
            })
            .collect();
        let kernel = kernel_128(&rows);
        assert_eq!(kernel.len(), 28);
        for c in kernel {
            assert_ne!(c, 0);
            for &r in &rows {
                assert_eq!((r & c).count_ones() % 2, 0);
            }
        }
    }
}
