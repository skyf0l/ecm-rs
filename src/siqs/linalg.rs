//! Dependencies between the relations: the matrix over GF(2) of their exponents (a column per
//! relation), without its singletons, then block Lanczos (see [`super::lanczos`]), or a dense
//! Gaussian elimination for the small matrices.

use super::lanczos::{Sparse, null_space};
use crate::stop::Stop;

/// Matrices with at most this many columns are solved by dense Gaussian elimination.
const DENSE_MAX: usize = 1500;

/// Attempts of block Lanczos (with other random starts) before the dense elimination.
const LANCZOS_ATTEMPTS: u64 = 3;

/// The dependencies between the columns `cols` (each the sorted rows of its nonzero entries,
/// below `rows`): bit `k` of entry `j` says whether column `j` is in the dependency `k` (up to
/// 64), with the size (rows, columns) of the matrix solved; `Err` if `stop` was requested.
pub(crate) fn dependencies(
    cols: &[Vec<u32>],
    rows: usize,
    seed: u64,
    stop: Stop<'_>,
) -> Result<(Vec<u64>, (usize, usize)), ()> {
    // Remove the singletons: a row with one nonzero entry, and its column, until none.
    let mut weight = vec![0u32; rows];
    for col in cols {
        for &r in col {
            weight[r as usize] += 1;
        }
    }
    let mut active: Vec<bool> = cols.iter().map(|col| !col.is_empty()).collect();
    loop {
        let mut changed = false;
        for (j, col) in cols.iter().enumerate() {
            if active[j] && col.iter().any(|&r| weight[r as usize] == 1) {
                active[j] = false;
                changed = true;
                for &r in col {
                    weight[r as usize] -= 1;
                }
            }
        }
        if !changed {
            break;
        }
    }
    // Renumber the rows in use, keep the active columns.
    let mut row_index = vec![u32::MAX; rows];
    let mut used_rows = 0u32;
    for (r, &w) in weight.iter().enumerate() {
        if w > 0 {
            row_index[r] = used_rows;
            used_rows += 1;
        }
    }
    let kept: Vec<usize> = (0..cols.len()).filter(|&j| active[j]).collect();
    let matrix: Vec<Vec<u32>> = kept
        .iter()
        .map(|&j| cols[j].iter().map(|&r| row_index[r as usize]).collect())
        .collect();
    let used_rows = used_rows as usize;
    let size = (used_rows, matrix.len());
    let mut out = vec![0u64; cols.len()];
    if matrix.is_empty() {
        return Ok((out, size));
    }

    let deps = if matrix.len() <= DENSE_MAX {
        dense(&matrix, used_rows)
    } else {
        let sparse = Sparse {
            rows: used_rows,
            cols: &matrix,
        };
        let mut found = None;
        for attempt in 0..LANCZOS_ATTEMPTS {
            let deps = null_space(&sparse, seed.wrapping_add(attempt), stop)?;
            let count = deps.iter().fold(0u64, |acc, &d| acc | d).count_ones();
            if count >= 8 {
                found = Some(deps);
                break;
            }
        }
        match found {
            Some(deps) => deps,
            None => dense(&matrix, used_rows),
        }
    };
    for (&j, &d) in kept.iter().zip(&deps) {
        out[j] = d;
    }
    Ok((out, size))
}

/// Dependencies between the columns of a small matrix, by Gaussian elimination on the columns
/// (each with the identity of the combination it is).
fn dense(cols: &[Vec<u32>], rows: usize) -> Vec<u64> {
    let n = cols.len();
    let row_words = rows.div_ceil(64);
    let id_words = n.div_ceil(64);
    let width = row_words + id_words;
    let mut m = vec![0u64; n * width];
    for (j, col) in cols.iter().enumerate() {
        let line = &mut m[j * width..(j + 1) * width];
        for &r in col {
            line[r as usize / 64] ^= 1 << (r % 64);
        }
        line[row_words + j / 64] |= 1 << (j % 64);
    }
    let mut pivoted = vec![false; n];
    for bit in 0..rows {
        let (word, mask) = (bit / 64, 1u64 << (bit % 64));
        let Some(p) = (0..n).find(|&j| !pivoted[j] && m[j * width + word] & mask != 0) else {
            continue;
        };
        pivoted[p] = true;
        let pivot: Vec<u64> = m[p * width..(p + 1) * width].to_vec();
        for j in 0..n {
            if j != p && m[j * width + word] & mask != 0 {
                for (x, &y) in m[j * width..(j + 1) * width].iter_mut().zip(&pivot) {
                    *x ^= y;
                }
            }
        }
    }
    let mut out = vec![0u64; n];
    let mut count = 0;
    for j in 0..n {
        if count == 64 {
            break;
        }
        let line = &m[j * width..(j + 1) * width];
        if line[..row_words].iter().all(|&w| w == 0) {
            for (k, o) in out.iter_mut().enumerate() {
                if line[row_words + k / 64] >> (k % 64) & 1 == 1 {
                    *o |= 1 << count;
                }
            }
            count += 1;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A random sparse matrix: `cols` columns of about `weight` entries below `rows`, with
    /// denser low rows.
    fn random_matrix(rows: usize, cols: usize, weight: usize, seed: u64) -> Vec<Vec<u32>> {
        let mut state = seed;
        (0..cols)
            .map(|_| {
                let mut col: Vec<u32> = (0..weight)
                    .map(|i| {
                        let r = super::super::next_random(&mut state);
                        // Half the entries in the first rows.
                        let range = if i % 2 == 0 { rows.min(64) } else { rows };
                        (r % range as u64) as u32
                    })
                    .collect();
                col.sort_unstable();
                // Pairs cancel over GF(2).
                let mut odd = Vec::new();
                for r in col {
                    if odd.last() == Some(&r) {
                        odd.pop();
                    } else {
                        odd.push(r);
                    }
                }
                odd
            })
            .collect()
    }

    fn check(cols: &[Vec<u32>], rows: usize, deps: &[u64]) -> u32 {
        let mut found = 0;
        for bit in 0..64 {
            let mut sum = vec![0u8; rows];
            let mut any = false;
            for (col, &d) in cols.iter().zip(deps) {
                if d >> bit & 1 == 1 {
                    any = true;
                    for &r in col {
                        sum[r as usize] ^= 1;
                    }
                }
            }
            if any {
                assert!(sum.iter().all(|&s| s == 0), "dependency {bit}");
                found += 1;
            }
        }
        found
    }

    #[test]
    fn dense_dependencies() {
        let cols = random_matrix(300, 330, 12, 1);
        let (deps, _) = dependencies(&cols, 300, 1, Stop::NEVER).unwrap();
        assert!(check(&cols, 300, &deps) >= 30);
    }

    #[test]
    fn lanczos_dependencies() {
        for (rows, seed) in [(3000, 2), (5000, 3)] {
            let cols = random_matrix(rows, rows + 80, 20, seed);
            let (deps, (r, c)) = dependencies(&cols, rows, seed, Stop::NEVER).unwrap();
            assert!(c > DENSE_MAX && r <= rows);
            let found = check(&cols, rows, &deps);
            assert!(found >= 40, "{found} dependencies");
            // The same from Lanczos only.
            let (deps, _) = dependencies(&cols, rows, seed + 100, Stop::NEVER).unwrap();
            assert!(check(&cols, rows, &deps) >= 40);
        }
    }
}
