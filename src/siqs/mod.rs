//! The self-initializing quadratic sieve (SIQS), for the composites without small factors
//! (after trial division and the first levels of curves, see [`crate::driver`]).
//!
//! Following Contini's thesis and the design of msieve and YAFU (studied, not copied):
//!
//! - Knuth-Schroeppel multiplier `k`, factor base of the primes modulo which `kn` is a square.
//! - Polynomials `Q(x) = ((A x + B)^2 - kn) / A` with `A` a product of `s` primes of the factor
//!   base close to `sqrt(2 kn) / M` (for `x` in `[-M, M)`), and the `2^(s-1)` values of `B`
//!   enumerated in Gray code order: the roots of `Q` modulo each prime are updated by one
//!   addition (self-initialization).
//! - The interval is sieved in blocks of 32 KiB (the L1 data cache) with rounded logarithms:
//!   the smallest primes are not sieved (small prime variation, the threshold accounts for
//!   them), the primes below the block size are sieved block by block, and the larger ones
//!   are pushed to per-block buckets once per polynomial (bucket sieve).
//! - Candidates are trial divided with the roots (and the bucket entries of their block), and
//!   kept if smooth but for one prime below a large prime bound (single large prime
//!   variation): two such relations with the same large prime make one.
//! - Linear algebra over GF(2) (see [`linalg`]): singletons removed, then block Lanczos (dense
//!   Gaussian elimination for the small matrices), and the square root step on each
//!   dependency.
//!
//! The polynomials are sieved by batches (all the `B` of one `A`); the batches are drawn in
//! a fixed order from the seed, and merged in this order (on several threads too), until
//! enough relations: the result only depends on the number and the seed.

mod fb;
mod lanczos;
mod linalg;
mod sieve;

pub(crate) use sieve::{Batch, sieve_a};

use crate::stop::Stop;
use rug::{Integer, ops::Pow};
use std::collections::{HashMap, HashSet};

/// Sieve block, in bytes (and entries): the L1 data cache of most x86 CPUs.
const BLOCK_BITS: u32 = 15;
const BLOCK: usize = 1 << BLOCK_BITS;

/// The primes below this bound are not sieved, but for their contribution to the threshold.
const SMALL_PRIME: u32 = 40;

/// Relations beyond the columns of the matrix.
const EXTRA_RELATIONS: usize = 64;

/// Smallest and largest sizes (in decimal digits) SIQS accepts.
pub(crate) const MIN_DIGITS: u32 = 20;
pub(crate) const MAX_DIGITS: u32 = 110;

/// Parameters by size of the number (decimal digits), interpolated between the rows.
struct Row {
    digits: u32,
    /// Primes in the factor base.
    fb: usize,
    /// Sieve interval `[-M, M)`, in blocks.
    blocks: usize,
    /// Large prime bound, as a multiple of the largest prime of the factor base.
    lp_mult: u32,
    /// Bits below `log2 |Q(x)|` (besides the large prime) where the candidates start.
    fudge: f64,
}

const PARAMS: [Row; 10] = [
    Row {
        digits: 20,
        fb: 60,
        blocks: 1,
        lp_mult: 20,
        fudge: 4.0,
    },
    Row {
        digits: 30,
        fb: 150,
        blocks: 1,
        lp_mult: 30,
        fudge: 4.0,
    },
    Row {
        digits: 40,
        fb: 350,
        blocks: 1,
        lp_mult: 40,
        fudge: 5.0,
    },
    Row {
        digits: 50,
        fb: 1000,
        blocks: 2,
        lp_mult: 50,
        fudge: 6.0,
    },
    Row {
        digits: 60,
        fb: 2500,
        blocks: 2,
        lp_mult: 60,
        fudge: 7.0,
    },
    Row {
        digits: 70,
        fb: 6000,
        blocks: 4,
        lp_mult: 70,
        fudge: 8.0,
    },
    Row {
        digits: 80,
        fb: 13000,
        blocks: 6,
        lp_mult: 80,
        fudge: 9.0,
    },
    Row {
        digits: 90,
        fb: 28000,
        blocks: 8,
        lp_mult: 90,
        fudge: 10.0,
    },
    Row {
        digits: 100,
        fb: 55000,
        blocks: 12,
        lp_mult: 100,
        fudge: 11.0,
    },
    Row {
        digits: 110,
        fb: 90000,
        blocks: 16,
        lp_mult: 110,
        fudge: 12.0,
    },
];

/// The parameters for `digits` digits.
struct Params {
    fb: usize,
    blocks: usize,
    lp_mult: u32,
    fudge: f64,
}

impl Params {
    fn for_digits(digits: f64) -> Self {
        let last = PARAMS.len() - 1;
        let i = PARAMS
            .iter()
            .position(|row| f64::from(row.digits) >= digits)
            .unwrap_or(last)
            .clamp(1, last);
        let (lo, hi) = (&PARAMS[i - 1], &PARAMS[i]);
        let t =
            ((digits - f64::from(lo.digits)) / f64::from(hi.digits - lo.digits)).clamp(0.0, 1.0);
        let lerp = |a: f64, b: f64| a + (b - a) * t;
        let mut params = Self {
            // Geometric interpolation of the factor base size.
            fb: ((lo.fb as f64).ln() * (1.0 - t) + (hi.fb as f64).ln() * t).exp() as usize,
            blocks: lerp(lo.blocks as f64, hi.blocks as f64).round() as usize,
            lp_mult: lerp(f64::from(lo.lp_mult), f64::from(hi.lp_mult)).round() as u32,
            fudge: lerp(lo.fudge, hi.fudge),
        };
        // Tuning only (the `siqs` example): `SIQS_TUNE=fb,blocks,lp_mult,fudge` (0: default).
        #[cfg(feature = "bench")]
        if let Ok(v) = std::env::var("SIQS_TUNE") {
            let v: Vec<f64> = v.split(',').map(|x| x.parse().unwrap()).collect();
            if v[0] > 0.0 {
                params.fb = v[0] as usize;
            }
            if v[1] > 0.0 {
                params.blocks = v[1] as usize;
            }
            if v[2] > 0.0 {
                params.lp_mult = v[2] as u32;
            }
            if v[3] > 0.0 {
                params.fudge = v[3];
            }
        }
        params
    }
}

/// A relation: `y^2 = (-1)^e0 * product of the factor base primes * large^2`... modulo `n`,
/// more precisely `y^2 = product of the factors * large (mod n)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Relation {
    /// `|A x + B| mod n`.
    pub y: Integer,
    /// Columns of the factors, with multiplicity, sorted: 0 for `-1`, `i + 1` for the prime
    /// `i` of the factor base.
    pub factors: Vec<u32>,
    /// The prime above the factor base, 1 if none.
    pub large: u64,
}

/// Primes of the factor base with the same rounded logarithm, sieved from buckets.
#[derive(Debug, Clone)]
pub(crate) struct Slice {
    start: usize,
    end: usize,
    log: u8,
}

/// What the sieving of all the polynomials shares.
pub(crate) struct Shared {
    n: Integer,
    kn: Integer,
    k: u32,
    fb: fb::FactorBase,
    /// Half the sieve interval: `x` in `[-m, m)`.
    m: u32,
    blocks: usize,
    threshold: u8,
    lp_bound: u64,
    /// Index of the first prime sieved, and of the first one sieved from buckets.
    first_sieved: usize,
    first_large: usize,
    /// Indices of the primes dividing the multiplier.
    k_primes: Vec<usize>,
    slices: Vec<Slice>,
    /// Inverses modulo 2^32 and bounds (see [`fb::FactorBase::inv32`]) of the primes from
    /// `first_sieved` to `first_large`, padded to a multiple of 16 with values never dividing.
    td_inv: Vec<u32>,
    td_lim: Vec<u32>,
    /// Primes in each `A`, range of their indices, and `log2` of the ideal `A`.
    s: usize,
    a_range: (usize, usize),
    a_bits: f64,
}

impl Shared {
    fn interval(&self) -> usize {
        self.blocks * BLOCK
    }
}

/// The primes of an `A` (indices in the factor base), and its position in the sequence.
#[derive(Debug, Clone)]
pub(crate) struct ASpec {
    pub index: usize,
    q: Vec<usize>,
}

/// Why SIQS did not start.
#[derive(Debug, PartialEq, Eq)]
pub(crate) enum Setup {
    /// A factor found while setting up (a prime of the factor base, or the multiplier, divides
    /// `n`; or `n` is a perfect power, even or too small).
    Factor(Integer),
    /// Out of the range of the parameters, or not composite.
    Unsupported,
}

/// A quadratic sieve on a number: the setup, the sequence of the `A` values, and the relations
/// collected.
pub(crate) struct Siqs {
    shared: std::sync::Arc<Shared>,
    rng: u64,
    used: HashSet<Vec<usize>>,
    next: usize,
    /// `A` values drawn but not sieved (see [`Siqs::give_back`]), next in the sequence.
    given_back: std::collections::VecDeque<ASpec>,
    /// Relations with all their factors in the factor base.
    fulls: Vec<Relation>,
    /// Relations with a large prime, grouped by it (in the order found).
    partials: HashMap<u64, Vec<Relation>>,
    /// Pairs of relations with the same large prime (sum over the groups of their size - 1).
    combined: usize,
    seen: HashSet<Integer>,
    needed: usize,
    polynomials: usize,
}

/// Summary of the setup, for the events.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Info {
    pub multiplier: u32,
    pub factor_base: usize,
    pub interval: usize,
    pub large_prime_bound: u64,
}

/// SplitMix64: a small generator for the choice of the `A` values.
fn next_random(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// Number of decimal digits of `n > 0`.
pub(crate) fn digits(n: &Integer) -> u32 {
    let d = (f64::from(n.significant_bits()) * std::f64::consts::LOG10_2).ceil() as u32;
    if *n < Integer::from(10).pow(d.saturating_sub(1)) {
        d - 1
    } else {
        d
    }
}

impl Siqs {
    /// Sets up SIQS on the odd composite `n`, not a perfect power, with the choices of the
    /// polynomials drawn from `seed`.
    pub(crate) fn new(n: &Integer, seed: u64) -> Result<Self, Setup> {
        let d = digits(n);
        if !(MIN_DIGITS..=MAX_DIGITS).contains(&d) || fb::is_prime(n) {
            return Err(Setup::Unsupported);
        }
        if n.is_even() {
            return Err(Setup::Factor(Integer::from(2)));
        }
        if n.is_perfect_power() {
            let root = (2..n.significant_bits())
                .find_map(|e| {
                    let (root, rem) = n.root_rem_ref(e).into();
                    let (root, rem): (Integer, Integer) = (root, rem);
                    (rem == 0).then_some(root)
                })
                .expect("a perfect power");
            return Err(Setup::Factor(root));
        }
        let k = fb::multiplier(n);
        let g = n.gcd_u_ref(k).into();
        let g: Option<u32> = g;
        if let Some(g) = g.filter(|&g| g > 1) {
            return Err(Setup::Factor(Integer::from(g)));
        }
        let kn = Integer::from(n * k);
        let params = Params::for_digits(f64::from(d));
        let blocks = params.blocks.max(1);
        let m = (blocks * BLOCK / 2) as u32;
        let log_q = f64::from(m).log2() + f64::from(kn.significant_bits()) / 2.0 - 0.5;
        // Logarithms scaled for the sums to fit in a byte.
        let scale = (240.0 / log_q).min(1.0);
        let fb = fb::FactorBase::new(n, k, params.fb.max(30), scale)
            .map_err(|fb::Divides(p)| Setup::Factor(Integer::from(p)))?;
        let pmax = fb.max();
        let lp_bound = (u64::from(pmax) * u64::from(params.lp_mult))
            .min(u64::from(pmax) * u64::from(pmax))
            .max(u64::from(pmax) + 1);
        let threshold = ((log_q - (lp_bound as f64).log2() - params.fudge) * scale)
            .clamp(8.0, 250.0)
            .round() as u8;
        let first_sieved = fb
            .primes
            .iter()
            .position(|&p| p >= SMALL_PRIME)
            .unwrap_or(fb.len());
        let first_large = fb
            .primes
            .iter()
            .position(|&p| p as usize >= BLOCK)
            .unwrap_or(fb.len())
            .max(first_sieved);
        let k_primes = (0..fb.len()).filter(|&i| fb.divides_k[i]).collect();
        // Slices of the same logarithm, and on the same side of the interval length (the
        // primes above it hit it at most once per root).
        let interval = (blocks * BLOCK) as u32;
        let mut slices: Vec<Slice> = Vec::new();
        for i in first_large..fb.len() {
            let above = |i: usize| fb.primes[i] >= interval;
            match slices.last_mut() {
                Some(slice) if slice.log == fb.logs[i] && above(slice.start) == above(i) => {
                    slice.end = i + 1;
                }
                _ => slices.push(Slice {
                    start: i,
                    end: i + 1,
                    log: fb.logs[i],
                }),
            }
        }
        // A: s primes near 2^(a_bits / s), of a size depending on the factor base.
        let a_bits = (f64::from(kn.significant_bits()) + 1.0) / 2.0 - f64::from(m).log2();
        let preferred = (f64::from(pmax) / 8.0).clamp(300.0, 4000.0).log2();
        let s = ((a_bits / preferred).round() as usize).max(1);
        let q_bits = a_bits / s as f64;
        let range_of = |lo: f64, hi: f64| {
            let lo = fb.primes.partition_point(|&p| f64::from(p) < lo);
            let hi = fb.primes.partition_point(|&p| f64::from(p) <= hi);
            (lo.max(first_sieved), hi.min(first_large))
        };
        let mut a_range = range_of(2f64.powf(q_bits) / 1.5, 2f64.powf(q_bits) * 1.5);
        let mut widen = 1.5;
        while a_range.1 < a_range.0 + 2 * s + 8 && widen < 64.0 {
            widen *= 1.5;
            a_range = range_of(2f64.powf(q_bits) / widen, 2f64.powf(q_bits) * widen);
        }
        if a_range.1 < a_range.0 + s + 2 {
            return Err(Setup::Unsupported);
        }
        let needed = fb.len() + 1 + EXTRA_RELATIONS;
        let padded = (first_large - first_sieved).next_multiple_of(16);
        let mut td_inv = fb.inv32[first_sieved..first_large].to_vec();
        let mut td_lim = fb.lim32[first_sieved..first_large].to_vec();
        td_inv.resize(padded, 1);
        td_lim.resize(padded, 0);
        let shared = Shared {
            n: n.clone(),
            kn,
            k,
            fb,
            m,
            blocks,
            threshold,
            lp_bound,
            first_sieved,
            first_large,
            k_primes,
            slices,
            td_inv,
            td_lim,
            s,
            a_range,
            a_bits,
        };
        Ok(Self {
            shared: std::sync::Arc::new(shared),
            rng: seed,
            used: HashSet::new(),
            next: 0,
            given_back: std::collections::VecDeque::new(),
            fulls: Vec::new(),
            partials: HashMap::new(),
            combined: 0,
            seen: HashSet::new(),
            needed,
            polynomials: 0,
        })
    }

    pub(crate) fn shared(&self) -> &std::sync::Arc<Shared> {
        &self.shared
    }

    pub(crate) fn info(&self) -> Info {
        let sh = &self.shared;
        Info {
            multiplier: sh.k,
            factor_base: sh.fb.len(),
            interval: sh.interval(),
            large_prime_bound: sh.lp_bound,
        }
    }

    /// Relations found (full ones and pairs of partial ones), relations needed, and
    /// polynomials sieved so far.
    pub(crate) fn progress(&self) -> (usize, usize, usize, usize) {
        (
            self.fulls.len(),
            self.combined,
            self.needed,
            self.polynomials,
        )
    }

    /// Whether there are enough relations for the linear algebra.
    pub(crate) fn enough(&self) -> bool {
        self.fulls.len() + self.combined >= self.needed
    }

    /// Gives back the `A` values drawn but not merged, in the order of the sequence: they are
    /// the next ones again.
    pub(crate) fn give_back(&mut self, specs: Vec<ASpec>) {
        for spec in specs.into_iter().rev() {
            self.given_back.push_front(spec);
        }
    }

    /// The next `A` of the sequence.
    pub(crate) fn next_a(&mut self) -> ASpec {
        if let Some(spec) = self.given_back.pop_front() {
            return spec;
        }
        let sh = std::sync::Arc::clone(&self.shared);
        let fb = &sh.fb;
        let (lo, hi) = sh.a_range;
        let s = sh.s;
        let usable = |i: usize| i >= sh.first_sieved && i < sh.first_large && !fb.divides_k[i];
        let log = |i: usize| f64::from(fb.primes[i]).log2();
        let mut best: Option<(f64, Vec<usize>)> = None;
        for attempt in 0..2000 {
            // s - 1 random primes of the range (all of them at random after many attempts:
            // the small numbers have few values of A).
            let (from, to) = if attempt < 1000 {
                (lo, hi)
            } else {
                (sh.first_sieved, sh.first_large)
            };
            let mut q: Vec<usize> = Vec::with_capacity(s);
            let mut tries = 0;
            while q.len() < s - 1 && tries < 100 * s {
                tries += 1;
                let i = from + (next_random(&mut self.rng) % (to - from) as u64) as usize;
                if !q.contains(&i) && usable(i) {
                    q.push(i);
                }
            }
            if q.len() < s - 1 {
                continue;
            }
            let rest = sh.a_bits - q.iter().map(|&i| log(i)).sum::<f64>();
            // The last prime: the closest to the remaining size, among the sieved primes
            // below the block size, that makes a new A.
            let target = 2f64.powf(rest);
            let j = fb.primes[..sh.first_large].partition_point(|&p| f64::from(p) < target);
            let found = (0..16).find_map(|d| {
                [j + d, j.wrapping_sub(d + 1)].into_iter().find_map(|i| {
                    if i >= sh.first_large || !usable(i) || q.contains(&i) {
                        return None;
                    }
                    let mut candidate = q.clone();
                    candidate.push(i);
                    candidate.sort_unstable();
                    (!self.used.contains(&candidate)).then(|| ((log(i) - rest).abs(), candidate))
                })
            });
            let Some((error, candidate)) = found else {
                continue;
            };
            if best.as_ref().is_none_or(|(e, _)| error < *e) {
                best = Some((error, candidate));
            }
            if error < 0.25 || (attempt > 100 && best.is_some()) {
                break;
            }
        }
        let (_, q) = best.expect("an unused A value");
        self.used.insert(q.clone());
        let index = self.next;
        self.next += 1;
        ASpec { index, q }
    }

    /// Adds the relations of a batch (the batches must come in the order of their `A`).
    pub(crate) fn merge(&mut self, batch: Batch) {
        self.polynomials += batch.polynomials;
        for relation in batch.fulls {
            if self.seen.insert(relation.y.clone()) {
                self.fulls.push(relation);
            }
        }
        for relation in batch.partials {
            if self.seen.insert(relation.y.clone()) {
                let group = self.partials.entry(relation.large).or_default();
                if !group.is_empty() {
                    self.combined += 1;
                }
                group.push(relation);
            }
        }
    }

    /// Asks for more relations (after a linear algebra that found no factor).
    pub(crate) fn more(&mut self) {
        self.needed += self.needed / 20 + EXTRA_RELATIONS;
    }

    /// The relations of the matrix: the full ones, then each partial one with the first of its
    /// group (in a fixed order: by large prime).
    fn relations(&self) -> Vec<(&Relation, Option<&Relation>)> {
        let mut relations: Vec<_> = self.fulls.iter().map(|r| (r, None)).collect();
        let mut larges: Vec<_> = self.partials.keys().copied().collect();
        larges.sort_unstable();
        for large in larges {
            let group = &self.partials[&large];
            for other in &group[1..] {
                relations.push((&group[0], Some(other)));
            }
        }
        relations
    }

    /// The linear algebra and the square roots: the prime factors of `n` found (at least two
    /// coprime factors whose product is `n`), `None` if no dependency split `n`, or `Err` if
    /// `stop` was requested. Also returns the size of the matrix and the number of
    /// dependencies.
    pub(crate) fn finish(
        &self,
        stop: Stop<'_>,
    ) -> Result<(Option<Vec<Integer>>, (usize, usize, usize)), ()> {
        let sh = &self.shared;
        let relations = self.relations();
        let columns: Vec<Vec<u32>> = relations
            .iter()
            .map(|(r, other)| {
                let mut odd: Vec<u32> = Vec::new();
                let mut push_all = |factors: &[u32]| {
                    for &f in factors {
                        match odd.iter().position(|&g| g == f) {
                            Some(i) => {
                                odd.swap_remove(i);
                            }
                            None => odd.push(f),
                        }
                    }
                };
                push_all(&r.factors);
                if let Some(other) = other {
                    push_all(&other.factors);
                }
                odd.sort_unstable();
                odd
            })
            .collect();
        let rows = sh.fb.len() + 1;
        let (deps, size) = linalg::dependencies(&columns, rows, next_random_seed(sh), stop)?;
        let dependencies = deps.iter().fold(0u64, |acc, &d| acc | d).count_ones() as usize;
        let size = (size.0, size.1, dependencies);
        let n = &sh.n;
        let mut parts = vec![n.clone()];
        for bit in 0..64 {
            if stop.requested() {
                return Err(());
            }
            let chosen: Vec<usize> = (0..relations.len())
                .filter(|&i| deps[i] >> bit & 1 == 1)
                .collect();
            if chosen.is_empty() {
                continue;
            }
            let Some(g) = self.square_root(&relations, &chosen) else {
                continue;
            };
            parts = parts
                .into_iter()
                .flat_map(|part| {
                    let d = Integer::from(part.gcd_ref(&g));
                    if d == 1 || d == part {
                        vec![part]
                    } else {
                        let other = Integer::from(&part / &d);
                        vec![d, other]
                    }
                })
                .collect();
            if parts.iter().all(fb::is_prime) {
                break;
            }
        }
        if parts.len() < 2 {
            return Ok((None, size));
        }
        parts.sort_unstable();
        Ok((Some(parts), size))
    }

    /// `x - y` for the dependency `chosen` (`x^2 = y^2 mod n`), whose gcd with `n` may split it;
    /// `None` if the exponents are not all even (never, for a true dependency).
    fn square_root(
        &self,
        relations: &[(&Relation, Option<&Relation>)],
        chosen: &[usize],
    ) -> Option<Integer> {
        let sh = &self.shared;
        let n = &sh.n;
        let mut exponents = vec![0u32; sh.fb.len() + 1];
        let mut x = Integer::from(1);
        let mut y = Integer::from(1);
        for &i in chosen {
            let (r, other) = relations[i];
            for relation in std::iter::once(r).chain(other) {
                x *= &relation.y;
                x %= n;
                for &f in &relation.factors {
                    exponents[f as usize] += 1;
                }
            }
            if other.is_some() {
                y *= r.large;
                y %= n;
            }
        }
        if exponents.iter().any(|e| e % 2 == 1) {
            return None;
        }
        for (i, &e) in exponents.iter().enumerate().skip(1) {
            if e > 0 {
                let p = Integer::from(sh.fb.primes[i - 1]);
                let power = p
                    .pow_mod(&Integer::from(e / 2), n)
                    .expect("a nonnegative exponent");
                y *= power;
                y %= n;
            }
        }
        Some(x - y)
    }
}

/// The seed of the random blocks of the linear algebra: fixed for a number.
fn next_random_seed(sh: &Shared) -> u64 {
    let mut state = sh.n.mod_u(u32::MAX) as u64 ^ 0x5eed;
    next_random(&mut state)
}

/// Statistics of a SIQS run on one thread, for the tuning example.
#[cfg(feature = "bench")]
#[derive(Debug, Clone)]
#[allow(missing_docs, reason = "internals for the benchmarks")]
pub struct Stats {
    /// Parts of the number found (coprime, sorted).
    pub parts: Vec<Integer>,
    pub multiplier: u32,
    pub factor_base: usize,
    pub interval: usize,
    pub large_prime_bound: u64,
    /// Primes in each `A`.
    pub a_primes: usize,
    /// Full relations, pairs of partial ones, partial ones, polynomials sieved.
    pub fulls: usize,
    pub combined: usize,
    pub partials: usize,
    pub polynomials: usize,
    /// Rows and columns of the matrix after the removal of the singletons, dependencies.
    pub matrix: (usize, usize, usize),
    /// Durations of the sieving and of the linear algebra with the square roots.
    pub sieve: std::time::Duration,
    pub linear_algebra: std::time::Duration,
}

/// Runs SIQS on `n` on one thread (see [`Stats`]); `None` if it does not apply to `n`.
#[cfg(feature = "bench")]
pub fn factor_stats(n: &Integer, seed: u64) -> Option<Stats> {
    let mut siqs = Siqs::new(n, seed).ok()?;
    let mut sieve = std::time::Duration::ZERO;
    let mut linear_algebra = std::time::Duration::ZERO;
    loop {
        let start = std::time::Instant::now();
        while !siqs.enough() {
            let a = siqs.next_a();
            let batch = sieve_a(&siqs.shared, &a, Stop::NEVER).expect("never stopped");
            siqs.merge(batch);
        }
        sieve += start.elapsed();
        let start = std::time::Instant::now();
        let (parts, matrix) = siqs.finish(Stop::NEVER).expect("never stopped");
        linear_algebra += start.elapsed();
        if let Some(parts) = parts {
            let info = siqs.info();
            return Some(Stats {
                parts,
                multiplier: info.multiplier,
                factor_base: info.factor_base,
                interval: info.interval,
                large_prime_bound: info.large_prime_bound,
                a_primes: siqs.shared.s,
                fulls: siqs.fulls.len(),
                combined: siqs.combined,
                partials: siqs.partials.values().map(Vec::len).sum(),
                polynomials: siqs.polynomials,
                matrix,
                sieve,
                linear_algebra,
            });
        }
        siqs.more();
    }
}

#[cfg(test)]
mod tests;
