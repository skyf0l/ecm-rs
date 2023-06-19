#![doc = include_str!("../README.md")]
#![deny(rust_2018_idioms)]
#![warn(missing_docs)]

mod ecm;
mod point;

pub use crate::ecm::*;
