//! Linear interpolation.

/// Linearly interpolate at point x given data points (x0, y0) and (x1, y1).
pub fn lerp(x0: f64, y0: f64, x1: f64, y1: f64, x: f64) -> f64 {
    y0 + (y1 - y0) * (x - x0) / (x1 - x0)
}

/// Piecewise linear interpolation on sorted data.
///
/// `xs` must be sorted in ascending order.
/// Returns None if x is outside the range.
pub fn piecewise_linear(xs: &[f64], ys: &[f64], x: f64) -> Option<f64> {
    if xs.is_empty() || x < xs[0] || x > *xs.last().unwrap() {
        return if !xs.is_empty() && (x - xs[0]).abs() < 1e-14 {
            Some(ys[0])
        } else if !xs.is_empty() && (x - *xs.last().unwrap()).abs() < 1e-14 {
            Some(*ys.last().unwrap())
        } else {
            None
        };
    }

    // Binary search for the interval
    let idx = match xs.binary_search_by(|v| v.partial_cmp(&x).unwrap()) {
        Ok(i) => return Some(ys[i]),
        Err(i) => {
            if i == 0 {
                0
            } else if i >= xs.len() {
                xs.len() - 1
            } else {
                i
            }
        }
    };

    let i = idx - 1;
    let t = (x - xs[i]) / (xs[i + 1] - xs[i]);
    Some(ys[i] + t * (ys[i + 1] - ys[i]))
}

/// Compute linear interpolation for multiple query points.
pub fn interpolate(xs: &[f64], ys: &[f64], query: &[f64]) -> Vec<Option<f64>> {
    query.iter().map(|&x| piecewise_linear(xs, ys, x)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lerp_midpoint() {
        let v = lerp(0.0, 0.0, 1.0, 1.0, 0.5);
        assert!((v - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_lerp_endpoints() {
        assert!((lerp(0.0, 1.0, 2.0, 3.0, 0.0) - 1.0).abs() < 1e-10);
        assert!((lerp(0.0, 1.0, 2.0, 3.0, 2.0) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_piecewise_linear_simple() {
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let ys = vec![0.0, 1.0, 4.0, 9.0];
        let v = piecewise_linear(&xs, &ys, 1.5).unwrap();
        assert!((v - 2.5).abs() < 1e-10);
    }

    #[test]
    fn test_piecewise_linear_at_nodes() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![0.0, 2.0, 4.0];
        assert!((piecewise_linear(&xs, &ys, 0.0).unwrap() - 0.0).abs() < 1e-10);
        assert!((piecewise_linear(&xs, &ys, 1.0).unwrap() - 2.0).abs() < 1e-10);
        assert!((piecewise_linear(&xs, &ys, 2.0).unwrap() - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_piecewise_linear_out_of_range() {
        let xs = vec![0.0, 1.0];
        let ys = vec![0.0, 1.0];
        assert!(piecewise_linear(&xs, &ys, -1.0).is_none());
        assert!(piecewise_linear(&xs, &ys, 2.0).is_none());
    }
}
