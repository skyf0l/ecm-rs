#![doc = include_str!("../README.md")]
#![deny(rust_2018_idioms)]
#![warn(missing_docs)]

mod ecm;
mod point;

pub use crate::ecm::{ecm, ecm_one_factor, ecm_with_params, Error};

/// Internals exposed for benchmarks only. Not part of the public API, no stability guarantees.
#[cfg(feature = "bench")]
#[doc(hidden)]
pub mod bench {
    pub use crate::ecm::{
        optimal_params, run_curve, stage1, stage1_multiplier, stage2, suyama_curve, trial_division,
        CurveOutcome,
    };
    pub use crate::point::Point;
}
