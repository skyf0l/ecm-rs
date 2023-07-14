#![doc = include_str!("../README.md")]
#![deny(rust_2018_idioms)]
#![warn(missing_docs)]

use ecm_sys::__mpz_struct;
use rug::Integer;
use std::{ffi::CStr, mem::MaybeUninit};

mod params;
pub use params::*;

/// Returns the version of the ECM library.
pub fn ecm_version() -> &'static str {
    unsafe { CStr::from_ptr(ecm_sys::ecm_version()).to_str().unwrap() }
}

/// Returns one factor of N using the Elliptic Curve Method.
pub fn ecm_factor(mut n: Integer) -> Integer {
    unsafe {
        let mut params = {
            let mut params = MaybeUninit::uninit();
            ecm_sys::ecm_init(params.as_mut_ptr());
            params.assume_init()
        };

        let mut factor = Integer::ZERO;
        let b1 = 10000f64;
        let raw_params = &mut params as *mut ecm_sys::__ecm_param_struct;

        // Args: factor, n, b1, params
        ecm_sys::ecm_factor(
            factor.as_raw_mut() as *mut __mpz_struct,
            n.as_raw_mut() as *mut __mpz_struct,
            b1,
            raw_params,
        );
        factor
    }
}

/// Returns all factors of N using the Elliptic Curve Method.
pub fn ecm(mut n: Integer) -> Vec<Integer> {
    let mut factors = Vec::new();

    while n != 1 {
        let factor = ecm_factor(n.clone());
        if factor == 1 {
            break;
        }
        factors.push(factor.clone());
        n /= factor;
    }

    factors.sort();
    factors
}
