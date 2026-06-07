//! Akima spline interpolation.
//!
//! The Akima spline produces a smooth curve that avoids overshooting
//! artifacts common with cubic splines.

/// Akima spline interpolation.
///
/// Given sorted data points, interpolate at x using the Akima method.
/// Returns None if x is outside the range or there are fewer than 2 points.
pub fn akima_interp(xs: &[f64], ys: &[f64], x: f64) -> Option<f64> {
    let n = xs.len();
    if n < 2 {
        return None;
    }
    if x < xs[0] || x > *xs.last().unwrap() {
        return None;
    }

    // Compute slopes
    let mut slopes = Vec::with_capacity(n);
    for i in 0..n - 1 {
        slopes.push((ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i]));
    }

    // Compute Akima weights
    let mut w = vec![0.0; n];
    for i in 2..n - 2 {
        let d_left = slopes[i - 1] - slopes[i - 2];
        let d_right = slopes[i] - slopes[i - 1];
        let wl = d_left.abs();
        let wr = d_right.abs();
        w[i] = if wl + wr > 1e-30 {
            (wr * slopes[i - 1] + wl * slopes[i]) / (wl + wr)
        } else {
            (slopes[i - 1] + slopes[i]) / 2.0
        };
    }

    // Boundary: mirror slopes
    if n >= 3 {
        w[0] = if slopes.len() > 1 {
            2.0 * slopes[0] - w.get(1).copied().unwrap_or(slopes[0])
        } else {
            slopes[0]
        };
    } else {
        w[0] = slopes[0];
    }
    if n >= 3 {
        w[1] = if slopes.len() > 1 {
            2.0 * slopes[0] - w[0]
        } else {
            slopes[0]
        };
    } else {
        w[1] = slopes.get(1).copied().unwrap_or(slopes[0]);
    }
    if n >= 3 {
        w[n - 2] = 2.0 * slopes[n - 2] - w.get(n - 3).copied().unwrap_or(slopes[n - 2]);
        w[n - 1] = 2.0 * slopes[n - 2] - w[n - 2];
    } else {
        w[n - 1] = slopes[n - 2];
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
            if i == 0 {
                0
            } else if i >= n {
                n - 2
            } else {
                i - 1
            }
        }
    };

    let h = xs[idx + 1] - xs[idx];
    let t = (x - xs[idx]) / h;
    let t2 = t * t;
    let t3 = t2 * t;

    let a = ys[idx];
    let b = w[idx];
    let c = 3.0 * slopes[idx] - 2.0 * w[idx] - w[idx + 1];
    let d = -2.0 * slopes[idx] + w[idx] + w[idx + 1];

    Some(a + b * t * h + c * t2 * h + d * t3 * h)
}

/// Evaluate Akima spline at multiple query points.
pub fn akima_interpolate(xs: &[f64], ys: &[f64], query: &[f64]) -> Vec<Option<f64>> {
    query.iter().map(|&x| akima_interp(xs, ys, x)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_akima_at_nodes() {
        let xs = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = vec![0.0, 1.0, 0.0, 1.0, 0.0];
        for i in 0..xs.len() {
            let v = akima_interp(&xs, &ys, xs[i]).unwrap();
            assert!((v - ys[i]).abs() < 1e-10, "Failed at node {}: got {} expected {}", i, v, ys[i]);
        }
    }

    #[test]
    fn test_akima_smooth() {
        let xs: Vec<f64> = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
        let ys: Vec<f64> = xs.iter().map(|&x| x.sin()).collect();
        let v = akima_interp(&xs, &ys, 2.5).unwrap();
        assert!((v - 2.5_f64.sin()).abs() < 0.1);
    }

    #[test]
    fn test_akima_no_overshoot() {
        // Step function - Akima should not overshoot much
        let xs = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = vec![0.0, 0.0, 1.0, 1.0, 1.0];
        let v = akima_interp(&xs, &ys, 1.5).unwrap();
        assert!(v >= -0.5 && v <= 1.5, "Akima overshoot: {}", v);
    }

    #[test]
    fn test_akima_out_of_range() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![0.0, 1.0, 2.0];
        assert!(akima_interp(&xs, &ys, -1.0).is_none());
        assert!(akima_interp(&xs, &ys, 3.0).is_none());
    }
}
