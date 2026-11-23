//! The multiplier and the factor base of the quadratic sieve, with the arithmetic modulo its
//! primes.

use rug::{Integer, integer::IsPrime};

/// Squarefree odd multipliers tried by [`multiplier`].
const MULTIPLIERS: [u32; 23] = [
    1, 3, 5, 7, 11, 13, 15, 17, 19, 21, 23, 29, 31, 33, 35, 37, 39, 41, 43, 47, 51, 53, 55,
];

/// Primes up to this bound score the multipliers.
const SCORE_PRIMES: u32 = 2000;

/// `a * b mod p`, for `a, b < p < 2^32`.
#[inline]
pub(crate) fn mul_mod(a: u32, b: u32, p: u32) -> u32 {
    (u64::from(a) * u64::from(b) % u64::from(p)) as u32
}

/// `b^e mod p`.
pub(crate) fn pow_mod(b: u32, mut e: u64, p: u32) -> u32 {
    let mut result = 1 % p;
    let mut b = b % p;
    while e > 0 {
        if e & 1 == 1 {
            result = mul_mod(result, b, p);
        }
        b = mul_mod(b, b, p);
        e >>= 1;
    }
    result
}

/// Inverse of `a` modulo `p` (`a` coprime to `p`, `p > 1`).
pub(crate) fn inverse(a: u32, p: u32) -> u32 {
    let (mut r0, mut r1) = (i64::from(p), i64::from(a % p));
    let (mut t0, mut t1) = (0i64, 1i64);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (t0, t1) = (t1, t0 - q * t1);
    }
    debug_assert_eq!(r0, 1, "{a} is not invertible modulo {p}");
    t0.rem_euclid(i64::from(p)) as u32
}

/// Legendre symbol `(a/p)` for an odd prime `p`: 1, 0 or -1.
pub(crate) fn legendre(a: u32, p: u32) -> i32 {
    match pow_mod(a, u64::from((p - 1) / 2), p) {
        0 => 0,
        1 => 1,
        _ => -1,
    }
}

/// A square root of `a` modulo the prime `p` (Tonelli-Shanks), if `a` is a square.
pub(crate) fn sqrt_mod(a: u32, p: u32) -> Option<u32> {
    let a = a % p;
    if p == 2 || a == 0 {
        return Some(a);
    }
    if legendre(a, p) != 1 {
        return None;
    }
    if p % 4 == 3 {
        return Some(pow_mod(a, u64::from(p / 4 + 1), p));
    }
    let (mut q, mut s) = (p - 1, 0u32);
    while q.is_multiple_of(2) {
        q /= 2;
        s += 1;
    }
    let z = (2..p).find(|&z| legendre(z, p) == -1)?;
    let mut m = s;
    let mut c = pow_mod(z, u64::from(q), p);
    let mut t = pow_mod(a, u64::from(q), p);
    let mut r = pow_mod(a, u64::from(q.div_ceil(2)), p);
    while t != 1 {
        let mut i = 0;
        let mut t2 = t;
        while t2 != 1 {
            t2 = mul_mod(t2, t2, p);
            i += 1;
        }
        let b = pow_mod(c, 1u64 << (m - i - 1), p);
        m = i;
        c = mul_mod(b, b, p);
        t = mul_mod(t, c, p);
        r = mul_mod(r, b, p);
    }
    Some(r)
}

/// Inverse of `p` modulo 2^32 (0 for even `p`), by Newton's iteration.
pub(crate) fn inverse_2_32(p: u32) -> u32 {
    if p.is_multiple_of(2) {
        return 0;
    }
    // Correct to 3 bits, then doubling the correct bits.
    let mut x = p;
    for _ in 0..4 {
        x = x.wrapping_mul(2u32.wrapping_sub(p.wrapping_mul(x)));
    }
    x
}

/// The Knuth-Schroeppel multiplier of `n`: the `k` making the most small primes (and 2) divide
/// the values `x^2 - kn` on average, for the size of `kn`.
pub(crate) fn multiplier(n: &Integer) -> u32 {
    let ln2 = std::f64::consts::LN_2;
    let mut best = (f64::NEG_INFINITY, 1);
    for k in MULTIPLIERS {
        let kn = Integer::from(n * k);
        let mut score = -0.5 * f64::from(k).ln();
        score += match kn.mod_u(8) {
            1 => 2.0 * ln2,
            5 => ln2,
            _ => 0.5 * ln2,
        };
        for p in primal::Primes::all()
            .skip(1)
            .take_while(|&p| p < SCORE_PRIMES as usize)
        {
            let p = p as u32;
            let lnp = f64::from(p).ln();
            if k.is_multiple_of(p) {
                score += lnp / f64::from(p);
            } else if legendre(kn.mod_u(p), p) == 1 {
                score += 2.0 * lnp / f64::from(p - 1);
            }
        }
        if score > best.0 {
            best = (score, k);
        }
    }
    best.1
}

/// The factor base: 2, the odd primes dividing the multiplier, and the odd primes `p` modulo
/// which `kn` is a nonzero square, with a square root of `kn` and the rounded logarithm.
pub(crate) struct FactorBase {
    pub primes: Vec<u32>,
    /// A square root of `kn` modulo each prime (0 for the primes dividing `kn`).
    pub sqrt: Vec<u32>,
    /// `log2(p)`, scaled (see [`FactorBase::new`]) and rounded.
    pub logs: Vec<u8>,
    /// Whether each prime divides the multiplier (a single root: sieved only by trial division).
    pub divides_k: Vec<bool>,
    /// Inverse of each odd prime modulo 2^32, and `(2^32 - 1) / p`: `p` divides `d < 2^32` if
    /// and only if `d * inv32 mod 2^32 <= lim32`.
    pub inv32: Vec<u32>,
    pub lim32: Vec<u32>,
}

/// Why no factor base: a prime considered divides `n`.
pub(crate) struct Divides(pub u32);

impl FactorBase {
    /// The factor base of `size` primes for `kn` (`n` odd), with logarithms scaled by `scale`,
    /// or a prime of the factor base range dividing `n`.
    pub(crate) fn new(n: &Integer, k: u32, size: usize, scale: f64) -> Result<Self, Divides> {
        let kn = Integer::from(n * k);
        let mut fb = Self {
            primes: Vec::with_capacity(size),
            sqrt: Vec::with_capacity(size),
            logs: Vec::with_capacity(size),
            divides_k: Vec::with_capacity(size),
            inv32: Vec::with_capacity(size),
            lim32: Vec::with_capacity(size),
        };
        let push = |fb: &mut Self, p: u32, root: u32, divides_k: bool| {
            fb.primes.push(p);
            fb.sqrt.push(root);
            fb.logs
                .push((f64::from(p).log2() * scale).round().clamp(1.0, 255.0) as u8);
            fb.divides_k.push(divides_k);
            fb.inv32.push(inverse_2_32(p));
            fb.lim32.push(u32::MAX / p);
        };
        push(&mut fb, 2, kn.mod_u(2), false);
        for p in primal::Primes::all().skip(1) {
            if fb.primes.len() >= size {
                break;
            }
            let p = u32::try_from(p).expect("a factor base below 2^32");
            let r = n.mod_u(p);
            if r == 0 {
                return Err(Divides(p));
            }
            if k.is_multiple_of(p) {
                push(&mut fb, p, 0, true);
                continue;
            }
            if let Some(root) = sqrt_mod(kn.mod_u(p), p) {
                push(&mut fb, p, root, false);
            }
        }
        Ok(fb)
    }

    pub(crate) fn len(&self) -> usize {
        self.primes.len()
    }

    /// The largest prime.
    pub(crate) fn max(&self) -> u32 {
        *self.primes.last().expect("a nonempty factor base")
    }
}

/// Whether `n` is a probable prime (the test of the driver).
pub(crate) fn is_prime(n: &Integer) -> bool {
    n.is_probably_prime(25) != IsPrime::No
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::str::FromStr;

    #[test]
    fn modular_helpers() {
        for p in [3u32, 5, 7, 13, 17, 97, 65537, 1_000_003, 4_294_967_291] {
            for a in [1u32, 2, 3, 10, 12345, p - 1] {
                let a = a % p;
                if a == 0 {
                    continue;
                }
                assert_eq!(mul_mod(a, inverse(a, p), p), 1, "{a} mod {p}");
                if let Some(r) = sqrt_mod(a, p) {
                    assert_eq!(mul_mod(r, r, p), a, "sqrt {a} mod {p}");
                } else {
                    assert_eq!(legendre(a, p), -1);
                }
            }
            if p % 2 == 1 {
                let inv = inverse_2_32(p);
                assert_eq!(p.wrapping_mul(inv), 1);
                for d in [
                    0u32,
                    1,
                    p,
                    p.wrapping_mul(2),
                    p.wrapping_mul(3).wrapping_add(1),
                    65535,
                    p.wrapping_add(65536),
                    u32::MAX,
                ] {
                    assert_eq!(d.wrapping_mul(inv) <= u32::MAX / p, d % p == 0, "{d} {p}");
                }
            }
        }
        // All the residues modulo a prime 1 mod 8 (the Tonelli-Shanks loop).
        let p = 257;
        for a in 1..p {
            if let Some(r) = sqrt_mod(a, p) {
                assert_eq!(mul_mod(r, r, p), a);
            }
        }
    }

    #[test]
    fn factor_base() {
        let p = Integer::from_str("1000000000000000000000000000000000000000000000000")
            .unwrap()
            .next_prime();
        let n = Integer::from(&p * &Integer::from(&p + 2u32).next_prime());
        let k = multiplier(&n);
        assert!(MULTIPLIERS.contains(&k));
        let fb = FactorBase::new(&n, k, 200, 1.0).ok().unwrap();
        let kn = Integer::from(&n * k);
        assert_eq!(fb.len(), 200);
        assert_eq!(fb.primes[0], 2);
        for (i, &p) in fb.primes.iter().enumerate().skip(1) {
            let r = fb.sqrt[i];
            assert_eq!(mul_mod(r, r, p), kn.mod_u(p), "{p}");
            assert_eq!(fb.divides_k[i], k.is_multiple_of(p));
        }
        assert!(fb.primes.windows(2).all(|w| w[0] < w[1]));
        // A prime of the range dividing n.
        let n =
            Integer::from(1_000_003u32) * Integer::from_str("1000000000000000000000007").unwrap();
        assert!(matches!(
            FactorBase::new(&n, 1, 50_000, 1.0),
            Err(Divides(1_000_003))
        ));
    }
}
