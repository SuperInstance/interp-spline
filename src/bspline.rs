//! B-spline basis functions and interpolation.

/// Evaluate the B-spline basis function of degree p at knot span i, parameter u.
///
/// Uses the Cox-de Boor recursion formula.
pub fn basis(i: usize, p: usize, u: f64, knots: &[f64]) -> f64 {
    if p == 0 {
        if u >= knots[i] && u < knots[i + 1] {
            return 1.0;
        }
        // Handle right endpoint
        if u == knots[i + 1] && (i + 1 == knots.len() - 1 || knots[i + 2..].iter().all(|&k| k == u)) {
            return 1.0;
        }
        return 0.0;
    }

    let d1 = knots[i + p] - knots[i];
    let d2 = knots[i + p + 1] - knots[i + 1];

    let c1 = if d1.abs() > 1e-14 {
        (u - knots[i]) / d1 * basis(i, p - 1, u, knots)
    } else {
        0.0
    };

    let c2 = if d2.abs() > 1e-14 {
        (knots[i + p + 1] - u) / d2 * basis(i + 1, p - 1, u, knots)
    } else {
        0.0
    };

    c1 + c2
}

/// Evaluate all non-zero B-spline basis functions at parameter u.
///
/// Returns a vector of basis function values for the given degree and knot vector.
pub fn basis_all(u: f64, p: usize, knots: &[f64]) -> Vec<f64> {
    let n = knots.len() - p - 1;
    (0..n).map(|i| basis(i, p, u, knots)).collect()
}

/// B-spline curve evaluation.
///
/// Evaluates a B-spline curve defined by control points and a knot vector
/// at parameter u.
#[allow(clippy::needless_range_loop)]
pub fn eval_curve(control_points: &[(f64, f64)], degree: usize, knots: &[f64], u: f64) -> (f64, f64) {
    let n = control_points.len();
    let mut x = 0.0;
    let mut y = 0.0;
    for i in 0..n {
        let b = basis(i, degree, u, knots);
        x += b * control_points[i].0;
        y += b * control_points[i].1;
    }
    (x, y)
}

/// B-spline interpolation: compute control points from data points.
///
/// Uses a uniform knot vector and solves the collocation system.
#[allow(clippy::type_complexity)]
pub fn interpolate(xs: &[f64], ys: &[f64], degree: usize) -> Option<(Vec<(f64, f64)>, Vec<f64>)> {
    let n = xs.len();
    if n < degree + 1 {
        return None;
    }

    // Create uniform knot vector
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

    // Parameter values (uniform)
    let params: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();

    // Build and solve collocation system using Gaussian elimination
    let mut mat = vec![vec![0.0; n]; n];
    #[allow(clippy::needless_range_loop)]
    for (row, &u) in params.iter().enumerate() {
        for col in 0..n {
            mat[row][col] = basis(col, degree, u, &knots);
        }
    }

    // Solve for x and y control points
    let cx = solve_system(&mat, xs)?;
    let cy = solve_system(&mat, ys)?;

    let control_points: Vec<(f64, f64)> = cx.iter().zip(cy.iter()).map(|(&x, &y)| (x, y)).collect();

    Some((control_points, knots))
}

#[allow(clippy::needless_range_loop)]
fn solve_system(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    let mut aug = vec![vec![0.0; n + 1]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a[i][j];
        }
        aug[i][n] = b[i];
    }

    for col in 0..n {
        // Pivot
        let mut max_row = col;
        let mut max_val = aug[col][col].abs();
        for row in col + 1..n {
            if aug[row][col].abs() > max_val {
                max_val = aug[row][col].abs();
                max_row = row;
            }
        }
        if max_val < 1e-14 {
            return None;
        }
        aug.swap(col, max_row);

        let pivot = aug[col][col];
        for j in col..=n {
            aug[col][j] /= pivot;
        }
        for row in 0..n {
            if row != col {
                let factor = aug[row][col];
                for j in col..=n {
                    aug[row][j] -= factor * aug[col][j];
                }
            }
        }
    }

    Some((0..n).map(|i| aug[i][n]).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basis_degree0() {
        let knots = vec![0.0, 1.0, 2.0];
        assert!((basis(0, 0, 0.5, &knots) - 1.0).abs() < 1e-10);
        assert!(basis(1, 0, 0.5, &knots).abs() < 1e-10);
    }

    #[test]
    fn test_basis_partition_of_unity() {
        let knots = vec![0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0];
        for u in [0.1, 0.3, 0.5, 0.7, 0.9] {
            let sum: f64 = basis_all(u, 2, &knots).iter().sum();
            assert!((sum - 1.0).abs() < 1e-10, "Partition of unity failed at u={}: sum={}", u, sum);
        }
    }

    #[test]
    fn test_bspline_curve_linear() {
        // Linear B-spline is just a line
        let cps = vec![(0.0, 0.0), (1.0, 1.0)];
        let knots = vec![0.0, 0.0, 1.0, 1.0];
        let (x, y) = eval_curve(&cps, 1, &knots, 0.5);
        assert!((x - 0.5).abs() < 1e-10);
        assert!((y - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_bspline_interpolation() {
        let xs = vec![0.0, 1.0, 2.0, 3.0];
        let ys = vec![0.0, 1.0, 0.5, 2.0];
        let result = interpolate(&xs, &ys, 2);
        assert!(result.is_some());
        let (cps, knots) = result.unwrap();
        // Should have 4 control points
        assert_eq!(cps.len(), 4);
        assert_eq!(knots.len(), 7); // n + degree + 1 = 4 + 2 + 1
    }
}
