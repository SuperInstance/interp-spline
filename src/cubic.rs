//! Cubic Hermite spline interpolation.

/// Cubic Hermite interpolation between two points with given derivatives.
///
/// Given (x0, y0, d0) and (x1, y1, d1) where d = dy/dx,
/// interpolate at point x.
pub fn hermite(x0: f64, y0: f64, d0: f64, x1: f64, y1: f64, d1: f64, x: f64) -> f64 {
    let h = x1 - x0;
    let t = (x - x0) / h;
    let t2 = t * t;
    let t3 = t2 * t;

    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;

    h00 * y0 + h10 * h * d0 + h01 * y1 + h11 * h * d1
}

/// Cubic Hermite spline with finite-difference derivatives.
///
/// Automatically computes derivatives using central differences
/// (forward/backward at endpoints).
pub fn cubic_hermite_spline(xs: &[f64], ys: &[f64], x: f64) -> Option<f64> {
    if xs.len() < 2 || x < xs[0] || x > *xs.last().unwrap() {
        return None;
    }

    let n = xs.len();
    // Compute derivatives
    let mut dy = vec![0.0; n];
    dy[0] = (ys[1] - ys[0]) / (xs[1] - xs[0]);
    dy[n - 1] = (ys[n - 1] - ys[n - 2]) / (xs[n - 1] - xs[n - 2]);
    for i in 1..n - 1 {
        dy[i] = (ys[i + 1] - ys[i - 1]) / (xs[i + 1] - xs[i - 1]);
    }

    // Find interval
    let idx = match xs.binary_search_by(|v| v.partial_cmp(&x).unwrap()) {
        Ok(i) => {
            if i == n - 1 {
                return Some(ys[n - 1]);
            }
            i
        }
        Err(i) => {
            if i == 0 || i >= n {
                return None;
            }
            i - 1
        }
    };

    Some(hermite(
        xs[idx], ys[idx], dy[idx],
        xs[idx + 1], ys[idx + 1], dy[idx + 1],
        x,
    ))
}

/// Natural cubic spline interpolation.
///
/// Computes the full tridiagonal system for natural boundary conditions (S''=0 at endpoints).
pub struct NaturalCubicSpline {
    xs: Vec<f64>,
    ys: Vec<f64>,
    /// Second derivatives at each knot.
    m: Vec<f64>,
}

impl NaturalCubicSpline {
    /// Create a new natural cubic spline from sorted data.
    pub fn new(xs: Vec<f64>, ys: Vec<f64>) -> Self {
        let n = xs.len();
        assert!(n >= 2, "Need at least 2 points");

        let mut m = vec![0.0; n];
        let mut h = vec![0.0; n - 1];
        for i in 0..n - 1 {
            h[i] = xs[i + 1] - xs[i];
        }

        if n == 2 {
            return Self { xs, ys, m };
        }

        // Solve tridiagonal system
        let mut alpha = vec![0.0; n];
        for i in 1..n - 1 {
            alpha[i] = 3.0 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1]);
        }

        let mut l = vec![0.0; n];
        let mut mu = vec![0.0; n];
        let mut z = vec![0.0; n];

        l[0] = 1.0;
        for i in 1..n - 1 {
            l[i] = 2.0 * (xs[i + 1] - xs[i - 1]) - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l[i];
            z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
        }

        l[n - 1] = 1.0;

        for i in (0..n - 1).rev() {
            m[i] = z[i] - mu[i] * m[i + 1];
        }

        Self { xs, ys, m }
    }

    /// Evaluate the spline at point x.
    pub fn eval(&self, x: f64) -> Option<f64> {
        if x < self.xs[0] || x > *self.xs.last().unwrap() {
            return None;
        }

        let n = self.xs.len();
        let idx = match self.xs.binary_search_by(|v| v.partial_cmp(&x).unwrap()) {
            Ok(i) => {
                if i == n - 1 {
                    return Some(self.ys[n - 1]);
                }
                i
            }
            Err(i) => {
                if i == 0 {
                    0
                } else if i >= n {
                    n - 2
                } else {
                    i - 1
                }
            }
        };

        let h = self.xs[idx + 1] - self.xs[idx];
        let a = (self.xs[idx + 1] - x) / h;
        let b = (x - self.xs[idx]) / h;

        Some(
            a * self.ys[idx]
                + b * self.ys[idx + 1]
                + (a * a * a - a) * self.m[idx] * h * h / 6.0
                + (b * b * b - b) * self.m[idx + 1] * h * h / 6.0,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hermite_at_nodes() {
        let v = hermite(0.0, 1.0, 0.0, 1.0, 2.0, 0.0, 0.0);
        assert!((v - 1.0).abs() < 1e-10);
        let v = hermite(0.0, 1.0, 0.0, 1.0, 2.0, 0.0, 1.0);
        assert!((v - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_cubic_hermite_spline_sin() {
        let xs: Vec<f64> = (0..=10).map(|i| i as f64 * 0.1 * std::f64::consts::PI).collect();
        let ys: Vec<f64> = xs.iter().map(|x| x.sin()).collect();
        let x = 0.25 * std::f64::consts::PI;
        let v = cubic_hermite_spline(&xs, &ys, x).unwrap();
        assert!((v - x.sin()).abs() < 0.01);
    }

    #[test]
    fn test_natural_cubic_exact_at_nodes() {
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let ys = vec![0.0, 1.0, 0.0, 1.0];
        let spline = NaturalCubicSpline::new(xs, ys);
        assert!((spline.eval(0.0).unwrap() - 0.0).abs() < 1e-10);
        assert!((spline.eval(1.0).unwrap() - 1.0).abs() < 1e-10);
        assert!((spline.eval(2.0).unwrap() - 0.0).abs() < 1e-10);
        assert!((spline.eval(3.0).unwrap() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_natural_cubic_smooth() {
        // Interpolate sin(x) with cubic spline
        let xs: Vec<f64> = (0..=10).map(|i| i as f64 * 0.1 * std::f64::consts::PI).collect();
        let ys: Vec<f64> = xs.iter().map(|x| x.sin()).collect();
        let spline = NaturalCubicSpline::new(xs, ys);
        let x = 0.15 * std::f64::consts::PI;
        let v = spline.eval(x).unwrap();
        assert!((v - x.sin()).abs() < 0.005);
    }
}
