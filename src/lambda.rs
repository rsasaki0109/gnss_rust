//! LAMBDA decorrelation and modified LAMBDA integer least-squares search.
//!
//! Port of libgnss++ `src/algorithms/lambda.cpp`, derived from RTKLIB.
//! Copyright (c) 2007-2013, T. Takasu. Copyright (c) 2024 LibGNSS++ Contributors.
//! See LICENSE and LICENSE-RTKLIB for redistribution terms.
//! Ambiguities are in cycles and covariance in squared cycles. This component
//! returns candidates, not a FIX decision or a position solution.

use std::fmt;

const EXACT_INTEGER_LIMIT: i64 = 1_i64 << 52;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LambdaError {
    InvalidDimensions,
    InvalidConfiguration,
    NonFinite,
    InvalidCovariance,
    NumericalRange,
    IterationLimit,
}
impl fmt::Display for LambdaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidDimensions => "expected 1..=128 ambiguities and a square covariance",
            Self::InvalidConfiguration => {
                "expected 1..=64 candidates and positive iteration limits"
            }
            Self::NonFinite => "LAMBDA input is not finite",
            Self::InvalidCovariance => "covariance must be symmetric and positive definite",
            Self::NumericalRange => "LAMBDA arithmetic exceeded its numerical range",
            Self::IterationLimit => "LAMBDA reduction or search exhausted its iteration limit",
        })
    }
}
impl std::error::Error for LambdaError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LambdaConfig {
    pub candidate_count: usize,
    pub max_search_iterations: usize,
    pub max_reduction_iterations: usize,
}
impl Default for LambdaConfig {
    fn default() -> Self {
        Self {
            candidate_count: 2,
            max_search_iterations: 10_000,
            max_reduction_iterations: 10_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct IntegerCandidate {
    pub ambiguities: Vec<i64>,
    /// Mahalanobis squared residual, ordered from smallest to largest.
    pub squared_residual: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LambdaResult {
    pub candidates: Vec<IntegerCandidate>,
    /// Diagonal D after reduction: Qz = L transpose * diag(D) * L.
    pub conditional_variances: Vec<f64>,
    /// Row-major Z, with decorrelated_float = Z transpose * float ambiguities.
    pub decorrelation_transform: Vec<Vec<i64>>,
    pub decorrelated_float: Vec<f64>,
    /// Row-major Z transpose * Q * Z.
    pub decorrelated_covariance: Vec<Vec<f64>>,
}
impl LambdaResult {
    /// Upstream second/best convention, including zero when best residual is
    /// exactly zero. None if only one candidate was requested. No FIX threshold
    /// is applied; a large ratio alone does not establish a valid RTK solution.
    pub fn ratio(&self) -> Option<f64> {
        let best = self.candidates.first()?.squared_residual;
        let second = self.candidates.get(1)?.squared_residual;
        Some(if best > 0.0 { second / best } else { 0.0 })
    }
}

/// Search the integer lattice using LAMBDA reduction and MLAMBDA enumeration.
/// Covariance is row-major, symmetric to relative tolerance 1e-12; its lower
/// triangle is authoritative within that tolerance, matching upstream LD.
/// Exact integer arithmetic is checked; no partial result is returned when a
/// limit is reached. FFRT, success-rate diagnostics and partial AR are pending.
pub fn search(
    float_ambiguities: &[f64],
    covariance: &[Vec<f64>],
    config: LambdaConfig,
) -> Result<LambdaResult, LambdaError> {
    let n = float_ambiguities.len();
    if !(1..=128).contains(&n) || covariance.len() != n || covariance.iter().any(|r| r.len() != n) {
        return Err(LambdaError::InvalidDimensions);
    }
    if !(1..=64).contains(&config.candidate_count)
        || config.max_search_iterations == 0
        || config.max_reduction_iterations == 0
    {
        return Err(LambdaError::InvalidConfiguration);
    }
    if !float_ambiguities.iter().all(|x| x.is_finite())
        || !covariance.iter().flatten().all(|x| x.is_finite())
    {
        return Err(LambdaError::NonFinite);
    }
    let mut q = vec![0.0; n * n];
    for i in 0..n {
        if covariance[i][i] <= 0.0 {
            return Err(LambdaError::InvalidCovariance);
        }
        for j in 0..n {
            let scale = covariance[i][j].abs().max(covariance[j][i].abs());
            if (covariance[i][j] - covariance[j][i]).abs() > 1e-12 * scale {
                return Err(LambdaError::InvalidCovariance);
            }
            q[i + j * n] = covariance[i.max(j)][i.min(j)];
        }
    }
    let (mut l, mut d) = factor(n, &q)?;
    let mut transform = vec![0_i64; n * n];
    let mut inverse_transpose = vec![0_i64; n * n];
    for i in 0..n {
        transform[i + i * n] = 1;
        inverse_transpose[i + i * n] = 1;
    }
    reduce(
        n,
        &mut l,
        &mut d,
        &mut transform,
        &mut inverse_transpose,
        config,
    )?;
    let z: Vec<f64> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| transform[j + i * n] as f64 * float_ambiguities[j])
                .sum()
        })
        .collect();
    for &value in &z {
        integer(value)?;
    }
    let reduced_candidates = enumerate(n, &l, &d, &z, config)?;
    let mut candidates = Vec::with_capacity(config.candidate_count);
    for (e, squared_residual) in reduced_candidates {
        let mut ambiguities = Vec::with_capacity(n);
        for i in 0..n {
            let mut sum = 0_i128;
            for j in 0..n {
                sum = sum
                    .checked_add(i128::from(inverse_transpose[i + j * n]) * i128::from(e[j]))
                    .ok_or(LambdaError::NumericalRange)?;
            }
            let value = i64::try_from(sum).map_err(|_| LambdaError::NumericalRange)?;
            check_integer(value)?;
            ambiguities.push(value);
        }
        candidates.push(IntegerCandidate {
            ambiguities,
            squared_residual,
        });
    }
    let mut qz = vec![vec![0.0; n]; n];
    // Two matrix products, preserving column-major reduction conventions.
    let mut q_times_z = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            q_times_z[i + j * n] = (0..n)
                .map(|k| q[i + k * n] * transform[k + j * n] as f64)
                .sum();
        }
    }
    for i in 0..n {
        for j in 0..n {
            qz[i][j] = (0..n)
                .map(|k| transform[k + i * n] as f64 * q_times_z[k + j * n])
                .sum();
            finite(qz[i][j])?;
        }
    }
    Ok(LambdaResult {
        candidates,
        conditional_variances: d,
        decorrelation_transform: (0..n)
            .map(|i| (0..n).map(|j| transform[i + j * n]).collect())
            .collect(),
        decorrelated_float: z,
        decorrelated_covariance: qz,
    })
}

fn finite(value: f64) -> Result<f64, LambdaError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(LambdaError::NumericalRange)
    }
}
fn check_integer(value: i64) -> Result<i64, LambdaError> {
    if value.unsigned_abs() < EXACT_INTEGER_LIMIT as u64 {
        Ok(value)
    } else {
        Err(LambdaError::NumericalRange)
    }
}
fn integer(value: f64) -> Result<i64, LambdaError> {
    let rounded = (finite(value)? + 0.5).floor();
    if rounded.abs() >= EXACT_INTEGER_LIMIT as f64 {
        return Err(LambdaError::NumericalRange);
    }
    Ok(rounded as i64)
}
fn sign(value: f64) -> f64 {
    if value <= 0.0 { -1.0 } else { 1.0 }
}

fn factor(n: usize, q: &[f64]) -> Result<(Vec<f64>, Vec<f64>), LambdaError> {
    let mut a = q.to_vec();
    let mut l = vec![0.0; n * n];
    let mut d = vec![0.0; n];
    for i in (0..n).rev() {
        d[i] = finite(a[i + i * n])?;
        if d[i] <= 0.0 {
            return Err(LambdaError::InvalidCovariance);
        }
        let root = d[i].sqrt();
        for j in 0..=i {
            l[i + j * n] = finite(a[i + j * n] / root)?;
        }
        for j in 0..i {
            for k in 0..=j {
                a[j + k * n] = finite(a[j + k * n] - l[i + k * n] * l[i + j * n])?;
            }
        }
        let diagonal = l[i + i * n];
        for j in 0..=i {
            l[i + j * n] = finite(l[i + j * n] / diagonal)?;
        }
    }
    Ok((l, d))
}

fn add_multiple(target: i64, source: i64, multiple: i64) -> Result<i64, LambdaError> {
    let value = i128::from(target) + i128::from(source) * i128::from(multiple);
    check_integer(i64::try_from(value).map_err(|_| LambdaError::NumericalRange)?)
}

fn reduce(
    n: usize,
    l: &mut [f64],
    d: &mut [f64],
    z: &mut [i64],
    inverse: &mut [i64],
    config: LambdaConfig,
) -> Result<(), LambdaError> {
    if n == 1 {
        return Ok(());
    }
    let mut j = (n - 2) as isize;
    let mut k = j;
    let mut iterations = 0;
    while j >= 0 {
        if iterations >= config.max_reduction_iterations {
            return Err(LambdaError::IterationLimit);
        }
        iterations += 1;
        let col = j as usize;
        if j <= k {
            for i in col + 1..n {
                let mu = integer(l[i + col * n])?;
                if mu != 0 {
                    for row in i..n {
                        l[row + col * n] = finite(l[row + col * n] - mu as f64 * l[row + i * n])?;
                    }
                    for row in 0..n {
                        z[row + col * n] = add_multiple(z[row + col * n], z[row + i * n], -mu)?;
                        inverse[row + i * n] =
                            add_multiple(inverse[row + i * n], inverse[row + col * n], mu)?;
                    }
                }
            }
        }
        let off = l[col + 1 + col * n];
        let delta = finite(d[col] + off * off * d[col + 1])?;
        if delta <= 0.0 {
            return Err(LambdaError::NumericalRange);
        }
        if delta + 1e-6 < d[col + 1] {
            let eta = d[col] / delta;
            let lam = d[col + 1] * off / delta;
            d[col] = finite(eta * d[col + 1])?;
            d[col + 1] = delta;
            if d[col] <= 0.0 {
                return Err(LambdaError::NumericalRange);
            }
            for row in 0..col {
                let a0 = l[col + row * n];
                let a1 = l[col + 1 + row * n];
                l[col + row * n] = finite(-off * a0 + a1)?;
                l[col + 1 + row * n] = finite(eta * a0 + lam * a1)?;
            }
            l[col + 1 + col * n] = finite(lam)?;
            for row in col + 2..n {
                l.swap(row + col * n, row + (col + 1) * n);
            }
            for row in 0..n {
                z.swap(row + col * n, row + (col + 1) * n);
                inverse.swap(row + col * n, row + (col + 1) * n);
            }
            k = j;
            j = (n - 2) as isize;
        } else {
            j -= 1;
        }
    }
    Ok(())
}

fn enumerate(
    n: usize,
    l: &[f64],
    d: &[f64],
    zs: &[f64],
    config: LambdaConfig,
) -> Result<Vec<(Vec<i64>, f64)>, LambdaError> {
    let m = config.candidate_count;
    let mut candidates: Vec<(Vec<i64>, f64)> = Vec::with_capacity(m);
    let mut imax = 0;
    let mut maxdist = 1e99;
    let mut s = vec![0.0; n * n];
    let mut dist = vec![0.0; n];
    let mut zb = vec![0.0; n];
    let mut z = vec![0_i64; n];
    let mut step = vec![0.0; n];
    let mut k = n - 1;
    zb[k] = zs[k];
    z[k] = integer(zb[k])?;
    let mut y = zb[k] - z[k] as f64;
    step[k] = sign(y);
    for _ in 0..config.max_search_iterations {
        let newdist = finite(dist[k] + y * y / d[k])?;
        if newdist < maxdist {
            if k != 0 {
                k -= 1;
                dist[k] = newdist;
                for i in 0..=k {
                    s[k + i * n] = finite(
                        s[k + 1 + i * n] + (z[k + 1] as f64 - zb[k + 1]) * l[k + 1 + i * n],
                    )?;
                }
                zb[k] = finite(zs[k] + s[k + k * n])?;
                z[k] = integer(zb[k])?;
                y = zb[k] - z[k] as f64;
                step[k] = sign(y);
            } else {
                if candidates.len() < m {
                    if candidates.is_empty() || newdist > candidates[imax].1 {
                        imax = candidates.len();
                    }
                    candidates.push((z.clone(), newdist));
                } else {
                    if newdist < candidates[imax].1 {
                        candidates[imax] = (z.clone(), newdist);
                        imax = 0;
                        for i in 0..m {
                            if candidates[imax].1 < candidates[i].1 {
                                imax = i;
                            }
                        }
                    }
                    maxdist = candidates[imax].1;
                }
                z[0] = add_multiple(z[0], 1, integer(step[0])?)?;
                y = zb[0] - z[0] as f64;
                step[0] = -step[0] - sign(step[0]);
            }
        } else {
            if k == n - 1 {
                if candidates.len() != m {
                    return Err(LambdaError::NumericalRange);
                }
                // Preserve upstream's tie order (it swaps equal residuals).
                for i in 0..m - 1 {
                    for j in i + 1..m {
                        if candidates[i].1 >= candidates[j].1 {
                            candidates.swap(i, j);
                        }
                    }
                }
                return Ok(candidates);
            }
            k += 1;
            z[k] = add_multiple(z[k], 1, integer(step[k])?)?;
            y = zb[k] - z[k] as f64;
            step[k] = -step[k] - sign(step[k]);
        }
    }
    Err(LambdaError::IterationLimit)
}
