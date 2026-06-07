//! Knot vector utilities for spline computations.

/// Create a uniform knot vector.
pub fn uniform(n: usize, degree: usize) -> Vec<f64> {
    let m = n + degree + 1;
    (0..m).map(|i| i as f64 / (m - 1) as f64).collect()
}

/// Create a clamped (open) uniform knot vector.
///
/// The first and last `degree+1` knots are repeated.
pub fn clamped(n: usize, degree: usize) -> Vec<f64> {
    let m = n + degree + 1;
    let mut knots = Vec::with_capacity(m);
    for i in 0..m {
        if i <= degree {
            knots.push(0.0);
        } else if i >= m - degree - 1 {
            knots.push(1.0);
        } else {
            knots.push((i - degree) as f64 / (n - degree) as f64);
        }
    }
    knots
}

/// Compute the knot span index for parameter u.
///
/// Returns the index i such that knots[i] <= u < knots[i+1].
pub fn find_span(degree: usize, knots: &[f64], u: f64) -> usize {
    let n = knots.len() - degree - 2;
    if u >= knots[n + 1] {
        return n;
    }
    if u <= knots[degree] {
        return degree;
    }

    let mut lo = degree;
    let mut hi = n + 1;
    let mut mid = (lo + hi) / 2;
    while u < knots[mid] || u >= knots[mid + 1] {
        if u < knots[mid] {
            hi = mid;
        } else {
            lo = mid;
        }
        mid = (lo + hi) / 2;
    }
    mid
}

/// Insert a knot into the knot vector and update control points.
///
/// Returns the new knot vector and updated control points.
#[allow(clippy::needless_range_loop)]
pub fn insert_knot(knots: &[f64], control_points: &[(f64, f64)], degree: usize, u: f64) -> (Vec<f64>, Vec<(f64, f64)>) {
    let span = find_span(degree, knots, u);

    let mut new_knots = knots.to_vec();
    new_knots.insert(span + 1, u);

    let n = control_points.len();
    let mut new_cp = Vec::with_capacity(n + 1);
    for i in 0..=span - degree {
        new_cp.push(control_points[i]);
    }
    for i in (span - degree + 1)..=span {
        let alpha = (u - knots[i]) / (knots[i + degree] - knots[i]);
        let x = alpha * control_points[i].0 + (1.0 - alpha) * control_points[i - 1].0;
        let y = alpha * control_points[i].1 + (1.0 - alpha) * control_points[i - 1].1;
        new_cp.push((x, y));
    }
    for i in span..n {
        new_cp.push(control_points[i]);
    }

    (new_knots, new_cp)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uniform_knots() {
        let knots = uniform(4, 2);
        assert_eq!(knots.len(), 7);
        assert!((knots[0] - 0.0).abs() < 1e-10);
        assert!((knots[6] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_clamped_knots() {
        let knots = clamped(4, 2);
        assert_eq!(knots.len(), 7);
        assert_eq!(knots[0], 0.0);
        assert_eq!(knots[1], 0.0);
        assert_eq!(knots[2], 0.0);
        assert_eq!(knots[4], 1.0);
        assert_eq!(knots[5], 1.0);
        assert_eq!(knots[6], 1.0);
    }

    #[test]
    fn test_find_span() {
        let knots = vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0];
        let span = find_span(2, &knots, 0.3);
        assert_eq!(span, 2);
        let span2 = find_span(2, &knots, 0.7);
        assert_eq!(span2, 3);
    }

    #[test]
    fn test_knot_insertion() {
        let knots = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let cps = vec![(0.0, 0.0), (0.5, 1.0), (1.0, 0.5), (2.0, 0.0)];
        let (new_knots, new_cp) = insert_knot(&knots, &cps, 2, 0.5);
        assert_eq!(new_knots.len(), 7);
        assert_eq!(new_cp.len(), 5);
    }
}
