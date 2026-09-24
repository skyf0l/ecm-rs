//! Compile-time choices of the arithmetic and their defaults.
//!
//! The `ecm-rs` tool includes this file (with `#[path]`) for `--printconfig`: what it prints is
//! what the library does.

use gmp_mpfr_sys::gmp;

/// Largest number of 64-bit limbs handled by the Montgomery arithmetic on fixed-size arrays
/// (`Mont`); larger moduli use GMP integers (`Plain`).
pub const MAX_LIMBS: usize = 16;

/// Smallest number of limbs from which `Mont` multiplies and squares with GMP's `mpn`
/// functions (a product, then a Montgomery reduction) instead of its own code.
///
/// Measured in stage 1: GMP's squaring is faster from 11 limbs up, and the multiplications cost
/// the same, but GMP's compiled code does not depend on the build profile, while our unrolled
/// CIOS got up to 30% slower from 13 limbs up with `lto = "fat"` and `codegen-units = 1`.
pub const GMP_LIMBS: usize = 11;

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
