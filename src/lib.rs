#![doc = include_str!("../README.md")]
#![deny(rust_2018_idioms)]
#![warn(missing_docs)]

mod arith;
mod cost;
mod curve;
mod driver;
mod ecm;
mod pm1;
mod point;
mod poly;
mod rho;
mod stage2;
mod stage2_poly;

pub use crate::ecm::{ecm, ecm_one_factor, ecm_with_params, Error};

/// Internals exposed for benchmarks only. Not part of the public API, no stability guarantees.
#[cfg(feature = "bench")]
#[doc(hidden)]
pub mod bench {
    pub use crate::driver::factor;
    pub use crate::ecm::{
        batch2_curve, curve, random_sigma, run_curve, square_curve, stage1, stage1_multiplier,
        stage2, suyama_curve, trial_division, CurveOutcome, Param,
    };
    pub use crate::pm1::Pm1;
    pub use crate::point::Point;
    pub use crate::rho::{ecm_prob, pm1_prob};
    pub use crate::stage2::Stage2Plan;
}
