//! Dense Gaussian measurement update with correlated measurement noise.
//!
//! Equivalent to the active-state update in libgnss++ `kalman.hpp` when all
//! states are explicitly active. Zero-valued states are valid here. Covariance
//! uses the Joseph form for stability; inputs are never modified on failure.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.

use std::fmt;

pub type Matrix = Vec<Vec<f64>>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterError {
    InvalidDimensions,
    NonFinite,
    InvalidCovariance,
    InvalidGate,
    InnovationRejected,
    NumericalFailure,
}
impl fmt::Display for FilterError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidDimensions => "invalid filter dimensions (states <= 128, rows <= 256)",
            Self::NonFinite => "filter input is not finite",
            Self::InvalidCovariance => "filter covariance must be symmetric positive definite",
            Self::InvalidGate => "innovation gate must be finite and positive",
            Self::InnovationRejected => "normalized innovation exceeded the configured gate",
            Self::NumericalFailure => "filter arithmetic failed",
        })
    }
}
impl std::error::Error for FilterError {}

#[derive(Debug, Clone, PartialEq)]
pub struct GaussianState {
    pub mean: Vec<f64>,
    pub covariance: Matrix,
}
#[derive(Debug, Clone, PartialEq)]
pub struct FilterUpdate {
    pub posterior: GaussianState,
    pub normalized_innovation_squared: f64,
    pub innovation_variances: Vec<f64>,
}

pub(crate) fn zeros(rows: usize, cols: usize) -> Matrix {
    vec![vec![0.0; cols]; rows]
}
pub(crate) fn transpose(a: &Matrix) -> Matrix {
    (0..a[0].len())
        .map(|j| a.iter().map(|r| r[j]).collect())
        .collect()
}
pub(crate) fn multiply(a: &Matrix, b: &Matrix) -> Matrix {
    let mut result = zeros(a.len(), b[0].len());
    for (row, result_row) in a.iter().zip(&mut result) {
        for (coefficient, b_row) in row.iter().zip(b) {
            for (value, b_value) in result_row.iter_mut().zip(b_row) {
                *value += coefficient * b_value;
            }
        }
    }
    result
}
fn symmetrize(a: &mut Matrix) {
    for i in 1..a.len() {
        let (previous, remaining) = a.split_at_mut(i);
        let row = &mut remaining[0];
        for (j, other) in previous.iter_mut().enumerate() {
            let value = (row[j] + other[i]) * 0.5;
            row[j] = value;
            other[i] = value;
        }
    }
}
fn square(a: &Matrix, n: usize) -> Result<(), FilterError> {
    if a.len() != n || a.iter().any(|r| r.len() != n) {
        return Err(FilterError::InvalidDimensions);
    }
    if !a.iter().flatten().all(|x| x.is_finite()) {
        return Err(FilterError::NonFinite);
    }
    Ok(())
}
pub(crate) fn cholesky(a: &Matrix) -> Result<Matrix, FilterError> {
    let n = a.len();
    square(a, n)?;
    for (i, row) in a.iter().enumerate() {
        for (j, other) in a.iter().enumerate().take(i) {
            let scale = row[j].abs().max(other[i].abs());
            if (row[j] - other[i]).abs() > scale * 1e-12 {
                return Err(FilterError::InvalidCovariance);
            }
        }
    }
    let mut l = zeros(n, n);
    for i in 0..n {
        for j in 0..=i {
            let value = a[i][j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if !value.is_finite() {
                return Err(FilterError::NumericalFailure);
            }
            if i == j {
                if value <= 0.0 {
                    return Err(FilterError::InvalidCovariance);
                }
                l[i][j] = value.sqrt();
            } else {
                l[i][j] = value / l[j][j];
            }
        }
    }
    Ok(l)
}
pub(crate) fn solve_cholesky(l: &Matrix, b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut x = vec![0.0; n];
    for i in 0..n {
        x[i] = (b[i] - (0..i).map(|j| l[i][j] * x[j]).sum::<f64>()) / l[i][i];
    }
    for i in (0..n).rev() {
        x[i] = (x[i] - (i + 1..n).map(|j| l[j][i] * x[j]).sum::<f64>()) / l[i][i];
    }
    x
}

/// Update from measurement-minus-prediction innovations. `design` maps states
/// to measurements. Optional gate is NIS / row count; rejection returns no
/// posterior. All states (including exact zeros) participate explicitly.
pub fn measurement_update(
    prior: &GaussianState,
    design: &Matrix,
    innovations: &[f64],
    noise: &Matrix,
    max_nis_per_row: Option<f64>,
) -> Result<FilterUpdate, FilterError> {
    let n = prior.mean.len();
    let m = innovations.len();
    if !(1..=128).contains(&n)
        || !(1..=256).contains(&m)
        || design.len() != m
        || design.iter().any(|r| r.len() != n)
    {
        return Err(FilterError::InvalidDimensions);
    }
    if !prior
        .mean
        .iter()
        .chain(innovations)
        .chain(design.iter().flatten())
        .all(|x| x.is_finite())
    {
        return Err(FilterError::NonFinite);
    }
    if max_nis_per_row.is_some_and(|v| !v.is_finite() || v <= 0.0) {
        return Err(FilterError::InvalidGate);
    }
    square(&prior.covariance, n)?;
    square(noise, m)?;
    cholesky(&prior.covariance)?;
    cholesky(noise)?;
    let f = multiply(&prior.covariance, &transpose(design));
    let mut s = multiply(design, &f);
    for i in 0..m {
        for j in 0..m {
            s[i][j] += noise[i][j];
        }
    }
    // H P H' is analytically symmetric. Floating-point cancellation between
    // correlated SD ambiguities can introduce asymmetry in this computed
    // matrix. Symmetrize derived products after validating the input matrices.
    symmetrize(&mut s);
    let l = cholesky(&s)?;
    let weighted = solve_cholesky(&l, innovations);
    let nis: f64 = innovations.iter().zip(&weighted).map(|(a, b)| a * b).sum();
    if !nis.is_finite() || nis < 0.0 {
        return Err(FilterError::NumericalFailure);
    }
    if max_nis_per_row.is_some_and(|limit| nis / m as f64 > limit) {
        return Err(FilterError::InnovationRejected);
    }
    let gain: Matrix = f.iter().map(|row| solve_cholesky(&l, row)).collect();
    let mean: Vec<_> = prior
        .mean
        .iter()
        .zip(&gain)
        .map(|(value, row)| value + row.iter().zip(innovations).map(|(a, b)| a * b).sum::<f64>())
        .collect();
    let mut ikh = multiply(&gain, design);
    for (i, row) in ikh.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            *value = if i == j { 1.0 } else { 0.0 } - *value;
        }
    }
    let mut p = multiply(&multiply(&ikh, &prior.covariance), &transpose(&ikh));
    let krk = multiply(&multiply(&gain, noise), &transpose(&gain));
    for i in 0..n {
        for j in 0..n {
            p[i][j] += krk[i][j];
        }
    }
    symmetrize(&mut p);
    if !mean.iter().chain(p.iter().flatten()).all(|x| x.is_finite()) {
        return Err(FilterError::NumericalFailure);
    }
    cholesky(&p)?;
    Ok(FilterUpdate {
        posterior: GaussianState {
            mean,
            covariance: p,
        },
        normalized_innovation_squared: nis,
        innovation_variances: (0..m).map(|i| s[i][i]).collect(),
    })
}
