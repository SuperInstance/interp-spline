//! # interp-spline
//!
//! Interpolation and spline methods: linear, cubic Hermite, B-spline, Akima.
//!
//! ## Modules
//! - `linear` — Linear interpolation
//! - `cubic` — Cubic Hermite spline interpolation
//! - `bspline` — B-spline basis functions and interpolation
//! - `akima` — Akima spline interpolation
//! - `knots` — Knot vector utilities

pub mod linear;
pub mod cubic;
pub mod bspline;
pub mod akima;
pub mod knots;
