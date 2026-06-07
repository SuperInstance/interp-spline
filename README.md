# interp-spline

Research-grade interpolation and spline methods in pure Rust.

## Features

- **Linear**: Piecewise linear interpolation
- **Cubic Hermite**: With finite-difference and natural cubic spline
- **B-spline**: Basis functions, curve evaluation, interpolation via collocation
- **Akima**: Robust interpolation avoiding overshoot artifacts
- **Knots**: Uniform/clamped knot vectors, knot insertion

## Usage

```rust
use interp_spline::cubic::NaturalCubicSpline;

fn main() {
    let xs = vec![0.0, 1.0, 2.0, 3.0];
    let ys = vec![0.0, 1.0, 0.0, 1.0];
    let spline = NaturalCubicSpline::new(xs, ys);
    println!("y(1.5) = {:.6}", spline.eval(1.5).unwrap());
}
```

## No Dependencies

This crate uses only `std`.

## License

MIT OR Apache-2.0
