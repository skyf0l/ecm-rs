//! Compile-time choices of the arithmetic and their defaults.
//!
//! The `ecm-rs` tool includes this file (with `#[path]`) for `--printconfig`: what it prints is
//! what the library does.

use gmp_mpfr_sys::gmp;

/// Largest number of 64-bit limbs handled by the Montgomery arithmetic on fixed-size arrays
/// (`Mont`); larger odd moduli use the Montgomery arithmetic of runtime length (`MontLarge`,
/// with GMP's `mpn` functions), or GMP integers (`Plain`) if [`MPN_ENABLED`] is `false`.
pub const MAX_LIMBS: usize = 16;

/// Smallest number of limbs from which `Mont` multiplies and squares with GMP's `mpn`
/// functions (a product, then a Montgomery reduction) instead of its own code.
///
/// Measured in stage 1: GMP's squaring is faster from 11 limbs up, and the multiplications cost
/// the same, but GMP's compiled code does not depend on the build profile, while our unrolled
/// CIOS got up to 30% slower from 13 limbs up with `lto = "fat"` and `codegen-units = 1`.
pub const GMP_LIMBS: usize = 11;

/// Smallest number of limbs from which `MontLarge` reduces with GMP's `mpn_redc_2` (two limbs
/// at a time) instead of `mpn_redc_1`.
///
/// Measured on chains of modular multiplications and squarings (GMP's product, then each
/// reduction; the fastest of 40 interleaved batches, on an i7-8750H): `mpn_redc_1` is 1-4%
/// faster at 17-19 limbs, they are equal at 20-22, `mpn_redc_2` is 3-6% faster from 24.
pub const REDC_2_LIMBS: usize = 20;

/// Smallest number of limbs from which `MontLarge` reduces with GMP's subquadratic
/// `mpn_redc_n` (a low half product and a wrap-around product) instead of `mpn_redc_2`.
///
/// Measured as [`REDC_2_LIMBS`]: `mpn_redc_n` is slower up to 44 limbs, equal at 48-52,
/// faster from 56 (by 3%, 6% at 64 limbs, 21% at 128, 35% at 256).
pub const REDC_N_LIMBS: usize = 52;

/// Whether GMP's limbs are our 64-bit limbs (and `mpn_redc_1` returns its carry, from GMP
/// 5.1): if not, GMP's `mpn` functions are never called.
pub const MPN_ENABLED: bool = gmp::LIMB_BITS == 64
    && gmp::NAIL_BITS == 0
    && (gmp::VERSION > 5 || (gmp::VERSION == 5 && gmp::VERSION_MINOR >= 1));

/// Largest memory (in bytes) used by the polynomial continuation of stage 2, by default.
pub const MAX_POLY_MEMORY: usize = 256 * 1024 * 1024;

/// Whether the CPU has BMI2 and ADX, for the copies of the hottest loops (the stage 1 ladder,
/// the pairs of stage 2) compiled with these extensions.
#[cfg(target_arch = "x86_64")]
#[inline]
pub fn bmi2_adx_detected() -> bool {
    std::is_x86_feature_detected!("bmi2") && std::is_x86_feature_detected!("adx")
}

/// GMP-ECM's `BASE2_THRESHOLD`: largest ratio of `k` to the size of `n` for the special
/// reduction.
pub const BASE2_THRESHOLD: f64 = 1.4;

/// GMP-ECM's `MOD_MINBASE2`: smallest `k` for the special reduction.
pub const BASE2_MIN_EXPONENT: u32 = 16;

/// Smallest number of limbs `k/64` from which the arithmetic modulo `2^k - 1` (`k` a multiple
/// of 64) multiplies with GMP's wrap-around products (`mpn_mulmod_bnm1`, `mpn_sqrmod_bnm1`)
/// instead of a full product and our fold.
///
/// Measured in instructions (Cachegrind, chains of 2000-5000 multiplications or squarings of
/// cofactors of `2^k - 1`, GMP 6.3 tuned for Skylake), multiplications / squarings: 8% / 12%
/// fewer at 8 to 17 limbs (GMP's own full product and addition, below its threshold of 15
/// limbs or for an odd `k/64`), 21% / 14% at 16, 24% / 18% at 18, 35% / 34% at 32, 43% at 64,
/// 45% at 128, 40-55% up to 2048 limbs (the gain depends on the powers of 2 dividing `k/64`:
/// GMP splits `B^rn - 1` into `(B^(rn/2) - 1)(B^(rn/2) + 1)` while `rn` is even). Below 8
/// limbs, the special reduction is rarely faster than Montgomery's.
pub const BASE2_WRAP_LIMBS: usize = 8;

/// Smallest number of limbs `k/64` from which the arithmetic modulo `2^k + 1` (`k` a multiple
/// of 64) multiplies with GMP's `mpn_mul_fft` (modulo `2^k + 1` directly) instead of a full
/// product and a fold.
///
/// Measured in instructions (Callgrind, GMP's products alone with GMP's best FFT size,
/// `mpn_fft_best_k`, against `mpn_mul_n` and an addition): 1% slower at 256 limbs, 3% at 384,
/// then faster: 14% at 448, 18% at 512, 32% at 1024, 45% at 2048 (GMP-ECM uses it from 512
/// limbs, `2^32768 + 1`). In chains of 2000 operations modulo a cofactor of `2^32768 + 1`: 19%
/// fewer instructions per multiplication, 12% per squaring. GMP's `mpn_mulmod_bknp1` (modulo `2^k + 1`
/// for `k/64` a multiple of 3, 5, 7...) is slower up to 36 limbs and 2-25% faster above, but
/// only helps numbers whose `k` is at least 1.5 times their size (a multiple of 3) or rare
/// ones: unused.
pub const BASE2_FFT_LIMBS: usize = 448;
