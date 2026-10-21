//! Native covariance-standardized residual screen (Baarda w statistic).
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
use crate::filter::{
    FilterError, GaussianState, Matrix, cholesky, multiply, solve_cholesky, transpose,
};

#[derive(Debug, Clone, PartialEq)]
pub struct ResidualScreen {
    /// (S^-1 r)_i / sqrt((S^-1)_ii), signed. These are not raw r_i/sigma_i.
    pub standardized_residuals: Vec<f64>,
    /// First row attaining the largest absolute statistic strictly above gate.
    pub worst_row: Option<usize>,
}
fn select(
    weighted: &[f64],
    inverse_diagonal: &[f64],
    threshold_sigma: f64,
) -> Result<ResidualScreen, FilterError> {
    if !threshold_sigma.is_finite() || threshold_sigma <= 0.0 {
        return Err(FilterError::InvalidGate);
    }
    let mut values = Vec::with_capacity(weighted.len());
    let mut worst = None;
    let mut largest = threshold_sigma;
    for (i, (&w, &v)) in weighted.iter().zip(inverse_diagonal).enumerate() {
        if !w.is_finite() || !v.is_finite() {
            return Err(FilterError::NumericalFailure);
        }
        if v <= 0.0 {
            return Err(FilterError::InvalidCovariance);
        }
        let value = w / v.sqrt();
        if !value.is_finite() {
            return Err(FilterError::NumericalFailure);
        }
        if value.abs() > largest {
            largest = value.abs();
            worst = Some(i);
        }
        values.push(value);
    }
    Ok(ResidualScreen {
        standardized_residuals: values,
        worst_row: worst,
    })
}

/// Actual native formula from caller-supplied S^-1. This low-level kernel
/// checks dimensions/finiteness/positive diagonal, not full positive definiteness.
/// Preserve entries of the derived inverse: native numerical inverses can have
/// tiny asymmetry and are not independently supplied measurement covariance.
/// Use `screen_measurements` to validate actual prior/noise covariance with SPD.
/// Rust returns errors for invalid rows which the native helper would skip.
/// For diagonal R this is the standardized post-fit residual R S^-1 r with
/// covariance R S^-1 R. With correlated R it remains the named native statistic.
pub fn screen_with_innovation_inverse(
    innovations: &[f64],
    inverse: &Matrix,
    threshold_sigma: f64,
) -> Result<ResidualScreen, FilterError> {
    let m = innovations.len();
    if !(1..=256).contains(&m) || inverse.len() != m {
        return Err(FilterError::InvalidDimensions);
    }
    if !innovations.iter().all(|v| v.is_finite()) {
        return Err(FilterError::NonFinite);
    }
    if inverse.iter().any(|row| row.len() != m) {
        return Err(FilterError::InvalidDimensions);
    }
    if !inverse.iter().flatten().all(|v| v.is_finite()) {
        return Err(FilterError::NonFinite);
    }
    let weighted: Vec<_> = inverse
        .iter()
        .map(|row| row.iter().zip(innovations).map(|(a, b)| a * b).sum())
        .collect();
    select(
        &weighted,
        &(0..m).map(|i| inverse[i][i]).collect::<Vec<_>>(),
        threshold_sigma,
    )
}

// Compensated products retain a high-accuracy residual for iterative refinement.
// This improves the new screening solve without changing Kalman update arithmetic.
fn refined_solve(s: &Matrix, l: &Matrix, b: &[f64]) -> Vec<f64> {
    let mut x = solve_cholesky(l, b);
    for _ in 0..2 {
        let residual: Vec<_> = s
            .iter()
            .zip(b)
            .map(|(row, &target)| {
                let mut sum: f64 = 0.0;
                let mut correction = 0.0;
                for (&a, &v) in row.iter().zip(&x) {
                    let product = a * v;
                    let product_error = a.mul_add(v, -product);
                    let next = sum + product;
                    let addition_error = if sum.abs() >= product.abs() {
                        (sum - next) + product
                    } else {
                        (product - next) + sum
                    };
                    correction += product_error + addition_error;
                    sum = next;
                }
                (target - sum) - correction
            })
            .collect();
        let delta = solve_cholesky(l, &residual);
        for (v, d) in x.iter_mut().zip(delta) {
            *v += d;
        }
    }
    x
}

/// Build S=H P H'+R from a checked Gaussian prior. Solve for S^-1 r and
/// inverse diagonal with Cholesky; no posterior is formed or committed.
pub fn screen_measurements(
    prior: &GaussianState,
    design: &Matrix,
    innovations: &[f64],
    noise: &Matrix,
    threshold_sigma: f64,
) -> Result<ResidualScreen, FilterError> {
    let n = prior.mean.len();
    let m = innovations.len();
    if !(1..=128).contains(&n)
        || !(1..=256).contains(&m)
        || design.len() != m
        || design.iter().any(|r| r.len() != n)
        || prior.covariance.len() != n
        || noise.len() != m
    {
        return Err(FilterError::InvalidDimensions);
    }
    if !prior
        .mean
        .iter()
        .chain(innovations)
        .chain(design.iter().flatten())
        .all(|v| v.is_finite())
    {
        return Err(FilterError::NonFinite);
    }
    cholesky(&prior.covariance)?;
    cholesky(noise)?;
    let mut s = multiply(design, &multiply(&prior.covariance, &transpose(design)));
    for (i, row) in s.iter_mut().enumerate() {
        for (j, v) in row.iter_mut().enumerate() {
            *v += noise[i][j];
        }
    }
    // Validate inputs before symmetrizing analytically symmetric derived S.
    for i in 1..m {
        let (previous, remaining) = s.split_at_mut(i);
        let row = &mut remaining[0];
        for (j, other) in previous.iter_mut().enumerate() {
            let v = (row[j] + other[i]) * 0.5;
            row[j] = v;
            other[i] = v;
        }
    }
    let l = cholesky(&s)?;
    let weighted = refined_solve(&s, &l, innovations);
    let diag = (0..m)
        .map(|i| {
            let mut e = vec![0.0; m];
            e[i] = 1.0;
            refined_solve(&s, &l, &e)[i]
        })
        .collect::<Vec<_>>();
    select(&weighted, &diag, threshold_sigma)
}
