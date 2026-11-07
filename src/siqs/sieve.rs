//! Sieving the polynomials of one `A` (see [`sieve_a`]): self-initialization of the `B` values
//! in Gray code order, a block sieve for the primes below the block size, a bucket sieve for
//! the larger ones, and the trial division of the candidates.

use super::{ASpec, BLOCK, BLOCK_BITS, Relation, Shared, fb::inverse};
use crate::stop::Stop;
use rug::Integer;

/// The relations found with the polynomials of one `A`.
#[derive(Debug, Default)]
pub(crate) struct Batch {
    /// Relations `y^2 = product of factor base primes` (modulo `n`).
    pub fulls: Vec<Relation>,
    /// Relations with one prime above the factor base.
    pub partials: Vec<Relation>,
    /// Polynomials sieved.
    pub polynomials: usize,
}

/// The `B` value of a polynomial and what is needed to evaluate it.
struct Poly {
    a: Integer,
    b: Integer,
    c: Integer,
}

/// Adds `d` to the root `r` modulo `p` (`r, d < p`), without branches (vectorizable).
#[inline(always)]
fn add_mod(r: u32, d: u32, p: u32) -> u32 {
    let t = r + d;
    t.min(t.wrapping_sub(p))
}

/// Subtracts `d` from the root `r` modulo `p` (`r, d < p`), without branches.
#[inline(always)]
fn sub_mod(r: u32, d: u32, p: u32) -> u32 {
    let t = r.wrapping_sub(d);
    t.min(t.wrapping_add(p))
}

/// `a b mod p` for `a, b < p < 2^26` (the quotient from floating point: `a b` is exact there).
#[inline(always)]
fn mul_mod(a: u32, b: u32, p: u32, inv_p: f64) -> u32 {
    let x = u64::from(a) * u64::from(b);
    let q = (f64::from(a) * f64::from(b) * inv_p) as u64;
    // The quotient is off by at most one.
    let r = x.wrapping_sub(q.wrapping_mul(u64::from(p))) as i64;
    let p = i64::from(p);
    let r = if r < 0 {
        r + p
    } else if r >= p {
        r - p
    } else {
        r
    };
    r as u32
}

/// The large prime buckets of a polynomial: for each block and slice of primes of the same
/// logarithm, the entries `(index of the prime) << BLOCK_BITS | offset in the block`. A prime
/// above the block size hits a block at most once per root: a bucket has room for twice the
/// primes of its slice, and the entries past the interval go to a last, discarded block.
struct Buckets {
    data: Vec<u32>,
    /// Start of the bucket (block, slice) in `data`, at `block * slices + slice`.
    start: Vec<usize>,
    count: Vec<usize>,
    slices: usize,
}

impl Buckets {
    fn new(sh: &Shared) -> Self {
        let slices = sh.slices.len();
        let mut start = Vec::with_capacity((sh.blocks + 1) * slices);
        let mut size = 0;
        for _ in 0..=sh.blocks {
            for slice in &sh.slices {
                start.push(size);
                size += 2 * (slice.end - slice.start);
            }
        }
        Self {
            data: vec![0; size],
            count: vec![0; start.len()],
            start,
            slices,
        }
    }

    fn clear(&mut self) {
        self.count.fill(0);
    }

    /// The entries of the bucket `(block, slice)`.
    fn get(&self, block: usize, slice: usize) -> &[u32] {
        let b = block * self.slices + slice;
        &self.data[self.start[b]..self.start[b] + self.count[b]]
    }
}

/// Sieves the `2^(s-1)` polynomials `((A x + B)^2 - kn) / A` of the `A` of `spec` over the
/// interval of `shared`, returning their relations, or `None` if `stop` is requested (checked
/// after each polynomial).
pub(crate) fn sieve_a(sh: &Shared, spec: &ASpec, stop: Stop<'_>) -> Option<Batch> {
    #[cfg(target_arch = "x86_64")]
    if std::is_x86_feature_detected!("avx2") {
        // SAFETY: the CPU has the features `sieve_a_avx2` is compiled for.
        return unsafe { sieve_a_avx2(sh, spec, stop) };
    }
    sieve_a_generic::<false>(sh, spec, stop)
}

/// [`sieve_a`] compiled with AVX2: the root updates and the trial division by 16 primes at a
/// time are vectorized with 32-bit lanes (`vpminud`, `vpmulld`), which SSE2 lacks.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn sieve_a_avx2(sh: &Shared, spec: &ASpec, stop: Stop<'_>) -> Option<Batch> {
    sieve_a_generic::<true>(sh, spec, stop)
}

/// The roots of the first polynomial of `A` (as indices in the interval: `x + m`) and their
/// updates for each `B_j`, computed from the primes of `A` and the `gamma_j` of the `B_j`
/// (`B_j = A / q_j * gamma_j`) without multiprecision arithmetic.
fn first_roots(
    sh: &Shared,
    spec: &ASpec,
    gammas: &[u32],
    logs: &[u8],
) -> (Vec<u32>, Vec<u32>, Vec<Vec<u32>>) {
    let fb = &sh.fb;
    let len = fb.len();
    let s = spec.q.len();
    let mut r1 = vec![0u32; len];
    let mut r2 = vec![0u32; len];
    let mut deltas = vec![vec![0u32; len]; s];
    let qs: Vec<u32> = spec.q.iter().map(|&i| fb.primes[i]).collect();
    let mut prefix = vec![1u32; s + 1];
    let mut residues = vec![0u32; s];
    for i in sh.first_sieved..len {
        if logs[i] == 0 {
            continue;
        }
        let p = fb.primes[i];
        let inv_p = 1.0 / f64::from(p);
        for j in 0..s {
            residues[j] = qs[j] % p;
            prefix[j + 1] = mul_mod(prefix[j], residues[j], p, inv_p);
        }
        let ainv = inverse(prefix[s], p);
        // From the last q_j down: others = A / q_j mod p = prefix[j] * (product after j).
        let mut suffix = 1;
        let mut b = 0;
        for j in (0..s).rev() {
            let others = mul_mod(prefix[j], suffix, p, inv_p);
            let gamma = gammas[j] % p;
            let bj = mul_mod(gamma, others, p, inv_p);
            b = add_mod(b, bj, p);
            // 2 B_j / A = 2 gamma_j / q_j.
            let inv_q = mul_mod(ainv, others, p, inv_p);
            deltas[j][i] = mul_mod(add_mod(gamma, gamma, p), inv_q, p, inv_p);
            suffix = mul_mod(suffix, residues[j], p, inv_p);
        }
        let t = fb.sqrt[i];
        let m = sh.m % p;
        let x1 = mul_mod(ainv, sub_mod(t, b, p), p, inv_p);
        let x2 = mul_mod(ainv, sub_mod(p - t, b, p), p, inv_p);
        r1[i] = add_mod(x1, m, p);
        r2[i] = add_mod(x2, m, p);
    }
    (r1, r2, deltas)
}

/// [`sieve_a`], with the AVX2 intrinsics if `AVX2` (only called from [`sieve_a_avx2`]).
#[inline(always)]
fn sieve_a_generic<const AVX2: bool>(sh: &Shared, spec: &ASpec, stop: Stop<'_>) -> Option<Batch> {
    let fb = &sh.fb;
    let s = spec.q.len();
    let interval = sh.interval() as u32;

    // A, and the B_j: B_j = (A/q_j) * gamma_j with B_j^2 = kn mod q_j, 0 mod the other q_i.
    let a = spec
        .q
        .iter()
        .fold(Integer::from(1), |a, &i| a * fb.primes[i]);
    let gammas: Vec<u32> = spec
        .q
        .iter()
        .map(|&i| {
            let q = fb.primes[i];
            let a_q = Integer::from(&a / q);
            let gamma = super::fb::mul_mod(fb.sqrt[i], inverse(a_q.mod_u(q), q), q);
            gamma.min(q - gamma)
        })
        .collect();
    let bj: Vec<Integer> = spec
        .q
        .iter()
        .zip(&gammas)
        .map(|(&i, &gamma)| Integer::from(&a / fb.primes[i]) * gamma)
        .collect();
    let mut b = bj.iter().fold(Integer::new(), |b, bj| b + bj);
    // e_j = +1 (false) or -1 (true), B = sum of e_j B_j.
    let mut negative = vec![false; s];

    // Not sieved: the small primes, the primes of A and those dividing the multiplier.
    let mut logs = fb.logs.clone();
    logs[..sh.first_sieved].fill(0);
    for &i in spec.q.iter().chain(&sh.k_primes) {
        logs[i] = 0;
    }
    let (mut r1, mut r2, deltas) = first_roots(sh, spec, &gammas, &logs);

    let blocks = sh.blocks;
    let slices = &sh.slices;
    let nslices = slices.len();
    let mut buckets = Buckets::new(sh);
    // Padded for the trial division by 16 primes at a time.
    let mut pos1 = vec![0u32; sh.first_sieved + sh.td_inv.len()];
    let mut pos2 = vec![0u32; sh.first_sieved + sh.td_inv.len()];
    let mut sieve = Box::new([0u8; BLOCK]);
    let mut batch = Batch::default();
    let mut candidates: Vec<u32> = Vec::new();
    let mut mark = vec![0u16; BLOCK];
    let mut large_hits: Vec<(u16, u32)> = Vec::new();
    let two_b_count = 1usize << (s - 1);

    for g in 0..two_b_count {
        if stop.requested() {
            return None;
        }
        // The next B in Gray code order: flip the sign of B_v.
        if g > 0 {
            let v = g.trailing_zeros() as usize;
            negative[v] = !negative[v];
            let range = sh.first_sieved..fb.len();
            let (p, d) = (&fb.primes[range.clone()], &deltas[v][range.clone()]);
            let (x1, x2) = (&mut r1[range.clone()], &mut r2[range]);
            if negative[v] {
                // B decreases by 2 B_v: the roots A^-1 (+-t - B) increase by 2 B_v A^-1.
                b -= &bj[v];
                b -= &bj[v];
                for (((x1, x2), &p), &d) in x1.iter_mut().zip(x2.iter_mut()).zip(p).zip(d) {
                    *x1 = add_mod(*x1, d, p);
                    *x2 = add_mod(*x2, d, p);
                }
            } else {
                b += &bj[v];
                b += &bj[v];
                for (((x1, x2), &p), &d) in x1.iter_mut().zip(x2.iter_mut()).zip(p).zip(d) {
                    *x1 = sub_mod(*x1, d, p);
                    *x2 = sub_mod(*x2, d, p);
                }
            }
        }
        let c = {
            let mut c = Integer::from(b.square_ref());
            c -= &sh.kn;
            c.div_exact_mut(&a);
            c
        };
        let poly = Poly {
            a: a.clone(),
            b: b.clone(),
            c,
        };

        // Bucket sieve of the large primes.
        buckets.clear();
        for (si, slice) in slices.iter().enumerate() {
            let (data, count, start) = (&mut buckets.data, &mut buckets.count, &buckets.start);
            if fb.primes[slice.start] >= interval {
                // At most one hit per root in the interval: no branch, the misses go to the
                // last block.
                for i in slice.start..slice.end {
                    let tag = (i as u32) << BLOCK_BITS;
                    for x in [r1[i], r2[i]] {
                        let block = ((x >> BLOCK_BITS) as usize).min(blocks);
                        let bucket = block * nslices + si;
                        data[start[bucket] + count[bucket]] = tag | (x & (BLOCK as u32 - 1));
                        count[bucket] += 1;
                    }
                }
            } else {
                for i in slice.start..slice.end {
                    let p = fb.primes[i];
                    let tag = (i as u32) << BLOCK_BITS;
                    for root in [r1[i], r2[i]] {
                        let mut x = root;
                        while x < interval {
                            let bucket = (x >> BLOCK_BITS) as usize * nslices + si;
                            data[start[bucket] + count[bucket]] = tag | (x & (BLOCK as u32 - 1));
                            count[bucket] += 1;
                            x += p;
                        }
                    }
                }
            }
        }

        let medium = sh.first_sieved..sh.first_large;
        pos1[medium.clone()].copy_from_slice(&r1[medium.clone()]);
        pos2[medium.clone()].copy_from_slice(&r2[medium.clone()]);
        for block in 0..blocks {
            sieve.fill(0);
            sieve_primes(
                &mut sieve,
                &fb.primes[medium.clone()],
                &logs[medium.clone()],
                &mut pos1[medium.clone()],
                &mut pos2[medium.clone()],
            );
            // Larger primes, from their buckets.
            for (si, slice) in slices.iter().enumerate() {
                apply_bucket(&mut sieve, buckets.get(block, si), slice.log);
            }
            // Candidates.
            candidates.clear();
            for (chunk_index, chunk) in sieve.chunks_exact(64).enumerate() {
                let max = chunk.iter().fold(0u8, |acc, &v| acc.max(v));
                if max < sh.threshold {
                    continue;
                }
                for (j, &v) in chunk.iter().enumerate() {
                    if v >= sh.threshold {
                        candidates.push((chunk_index * 64 + j) as u32);
                    }
                }
            }
            if candidates.is_empty() {
                continue;
            }
            // The large primes hitting the candidates, from the buckets of the block.
            for (k, &offset) in candidates.iter().enumerate() {
                mark[offset as usize] = k as u16 + 1;
            }
            large_hits.clear();
            for si in 0..nslices {
                for &e in buckets.get(block, si) {
                    let k = mark[(e & (BLOCK as u32 - 1)) as usize];
                    if k != 0 {
                        large_hits.push((k - 1, e >> BLOCK_BITS));
                    }
                }
            }
            for &offset in &candidates {
                mark[offset as usize] = 0;
            }
            for (k, &offset) in candidates.iter().enumerate() {
                let large = large_hits
                    .iter()
                    .filter(|&&(c, _)| usize::from(c) == k)
                    .map(|&(_, i)| i as usize);
                let index = (block * BLOCK) as u32 + offset;
                let positions = (&pos1[..], &pos2[..]);
                if let Some(relation) =
                    trial_divide::<AVX2>(sh, spec, &poly, (index, offset), positions, large)
                {
                    if relation.large == 1 {
                        batch.fulls.push(relation);
                    } else {
                        batch.partials.push(relation);
                    }
                }
            }
        }
        batch.polynomials += 1;
    }
    Some(batch)
}

/// Adds `log` at the positions of the primes in the block (`pos1` and `pos2`, below the prime,
/// which is below the block size), and moves the positions to the next block.
#[inline(always)]
fn sieve_primes(
    sieve: &mut [u8; BLOCK],
    primes: &[u32],
    logs: &[u8],
    pos1: &mut [u32],
    pos2: &mut [u32],
) {
    for (((&p, &lg), x1), x2) in primes.iter().zip(logs).zip(pos1).zip(pos2) {
        let step = p as usize;
        let (mut lo, mut hi) = (*x1 as usize, *x2 as usize);
        if lo > hi {
            std::mem::swap(&mut lo, &mut hi);
        }
        while hi < BLOCK {
            // SAFETY: lo <= hi < BLOCK.
            unsafe {
                *sieve.get_unchecked_mut(lo) = sieve.get_unchecked(lo).wrapping_add(lg);
                *sieve.get_unchecked_mut(hi) = sieve.get_unchecked(hi).wrapping_add(lg);
            }
            lo += step;
            hi += step;
        }
        if lo < BLOCK {
            sieve[lo] = sieve[lo].wrapping_add(lg);
            lo += step;
        }
        // Both are now in [BLOCK, BLOCK + p).
        *x1 = (lo - BLOCK) as u32;
        *x2 = (hi - BLOCK) as u32;
    }
}

/// For 16 primes (their `inv`, `lim` and positions `x1`, `x2`, see [`trial_divide`]): bit `j`
/// set if the prime `j` divides `x1 + t` or `x2 + t`.
#[inline(always)]
fn divisible_mask<const AVX2: bool>(t: u32, [inv, lim, x1, x2]: [&[u32]; 4]) -> u32 {
    #[cfg(target_arch = "x86_64")]
    if AVX2 {
        // SAFETY: only called with AVX2 from `sieve_a_avx2`, on 16 lanes.
        return unsafe { divisible_mask_avx2(t, [inv, lim, x1, x2]) };
    }
    let mut mask = 0;
    for j in 0..16 {
        let d1 = x1[j].wrapping_add(t).wrapping_mul(inv[j]);
        let d2 = x2[j].wrapping_add(t).wrapping_mul(inv[j]);
        mask |= u32::from((d1 <= lim[j]) | (d2 <= lim[j])) << j;
    }
    mask
}

/// [`divisible_mask`] with AVX2 intrinsics: two times 8 lanes.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn divisible_mask_avx2(t: u32, [inv, lim, x1, x2]: [&[u32]; 4]) -> u32 {
    use std::arch::x86_64::{
        __m256i, _mm256_add_epi32, _mm256_castsi256_ps, _mm256_cmpeq_epi32, _mm256_loadu_si256,
        _mm256_min_epu32, _mm256_movemask_ps, _mm256_mullo_epi32, _mm256_or_si256,
        _mm256_set1_epi32,
    };
    assert!(inv.len() >= 16 && lim.len() >= 16 && x1.len() >= 16 && x2.len() >= 16);
    let tv = _mm256_set1_epi32(t as i32);
    let mut mask = 0;
    for half in 0..2 {
        let at = |v: &[u32]| v[8 * half..].as_ptr().cast::<__m256i>();
        // SAFETY: 8 lanes from 8 * half, within the 16 checked above (unaligned loads).
        let (inv, lim, x1, x2) = unsafe {
            (
                _mm256_loadu_si256(at(inv)),
                _mm256_loadu_si256(at(lim)),
                _mm256_loadu_si256(at(x1)),
                _mm256_loadu_si256(at(x2)),
            )
        };
        let d1 = _mm256_mullo_epi32(_mm256_add_epi32(x1, tv), inv);
        let d2 = _mm256_mullo_epi32(_mm256_add_epi32(x2, tv), inv);
        // d <= lim (unsigned) if min(d, lim) == d.
        let c1 = _mm256_cmpeq_epi32(_mm256_min_epu32(d1, lim), d1);
        let c2 = _mm256_cmpeq_epi32(_mm256_min_epu32(d2, lim), d2);
        let bits = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_or_si256(c1, c2)));
        mask |= (bits as u32) << (8 * half);
    }
    mask
}

/// Adds `log` at the offsets of the bucket entries.
#[inline(always)]
fn apply_bucket(sieve: &mut [u8; BLOCK], bucket: &[u32], log: u8) {
    for &e in bucket {
        let x = (e & (BLOCK as u32 - 1)) as usize;
        sieve[x] = sieve[x].wrapping_add(log);
    }
}

/// The relation of the candidate at `index` in the interval (`offset` in its block), if its
/// value is smooth but for at most one prime below the large prime bound.
#[inline(always)]
fn trial_divide<const AVX2: bool>(
    sh: &Shared,
    spec: &ASpec,
    poly: &Poly,
    (index, offset): (u32, u32),
    (pos1, pos2): (&[u32], &[u32]),
    large: impl Iterator<Item = usize>,
) -> Option<Relation> {
    let fb = &sh.fb;
    let x = i64::from(index) - i64::from(sh.m);
    // y = A x + B, and Q(x) = (y^2 - kn) / A = (A x + 2B) x + C.
    let mut y = Integer::from(&poly.a * x);
    y += &poly.b;
    let mut q = Integer::from(&y + &poly.b);
    q *= x;
    q += &poly.c;
    let mut factors: Vec<u32> = Vec::with_capacity(32);
    if q < 0 {
        factors.push(0);
        q = -q;
    }
    if q == 0 {
        return None;
    }
    let twos = q.find_one(0).expect("nonzero");
    q >>= twos;
    factors.extend(std::iter::repeat_n(1, twos as usize));
    let divide = |q: &mut Integer, i: usize, factors: &mut Vec<u32>| {
        let p = fb.primes[i];
        while q.is_divisible_u(p) {
            q.div_exact_u_mut(p);
            factors.push(i as u32 + 1);
        }
    };
    // Primes not sieved (small, or dividing the multiplier).
    for i in 1..sh.first_sieved {
        divide(&mut q, i, &mut factors);
    }
    for &i in &sh.k_primes {
        if i >= sh.first_sieved {
            divide(&mut q, i, &mut factors);
        }
    }
    // The primes of A: once for A, and in Q(x).
    for &i in &spec.q {
        factors.push(i as u32 + 1);
        divide(&mut q, i, &mut factors);
    }
    // Primes below the block size: the next position of a root (after the block) is at a
    // multiple of p from the candidate. p divides d if d / p (by the inverse of p modulo 2^32)
    // is at most (2^32 - 1) / p.
    let t = BLOCK as u32 - offset;
    let range = sh.first_sieved..sh.first_sieved + sh.td_inv.len();
    let (inv, lim) = (&sh.td_inv[..], &sh.td_lim[..]);
    let (x1, x2) = (&pos1[range.clone()], &pos2[range]);
    let chunks = inv
        .chunks_exact(16)
        .zip(lim.chunks_exact(16))
        .zip(x1.chunks_exact(16).zip(x2.chunks_exact(16)));
    for (c, ((inv, lim), (x1, x2))) in chunks.enumerate() {
        let mut flags = divisible_mask::<AVX2>(t, [inv, lim, x1, x2]);
        while flags != 0 {
            let j = flags.trailing_zeros() as usize;
            flags &= flags - 1;
            divide(&mut q, sh.first_sieved + 16 * c + j, &mut factors);
        }
    }
    // Larger primes: the entries of the block's buckets at the candidate.
    for i in large {
        divide(&mut q, i, &mut factors);
    }
    // Tests: the roots found every prime of the factor base dividing Q(x).
    #[cfg(test)]
    for &p in &fb.primes {
        assert!(!q.is_divisible_u(p), "{p} missed at {x}");
    }
    let large = if q == 1 {
        1
    } else {
        match q.to_u64() {
            Some(large) if large < sh.lp_bound => large,
            _ => return None,
        }
    };
    factors.sort_unstable();
    if y < 0 {
        y = -y;
    }
    y %= &sh.n;
    Some(Relation { y, factors, large })
}
