/// The same Gaussian update and NIS gate using a whitened least-squares QR
/// in unit-prior coordinates. The identity rows represent the supplied prior,
/// not an extra state constraint or covariance floor. Column equilibration is
/// undone in the mean/covariance. Useful with diffuse ambiguity variances and
/// small carrier noise, without forming/inverting the innovation covariance.
pub fn measurement_update_whitened_qr(
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
        .all(|v| v.is_finite())
    {
        return Err(FilterError::NonFinite);
    }
    if max_nis_per_row.is_some_and(|v| !v.is_finite() || v <= 0.0) {
        return Err(FilterError::InvalidGate);
    }
    square(&prior.covariance, n)?;
    square(noise, m)?;
    let lp = cholesky(&prior.covariance)?;
    let lr = cholesky(noise)?;
    let hl = multiply(design, &lp);
    let innovation_variances: Vec<_> = hl
        .iter()
        .enumerate()
        .map(|(i, row)| row.iter().map(|v| v * v).sum::<f64>() + noise[i][i])
        .collect();
    let forward = |rhs: &[f64]| {
        let mut y = vec![0.0; m];
        for i in 0..m {
            y[i] = (rhs[i] - (0..i).map(|j| lr[i][j] * y[j]).sum::<f64>()) / lr[i][i];
        }
        y
    };
    let mut a = zeros(n + m, n);
    let mut rhs = vec![0.0; n];
    rhs.extend(forward(innovations));
    for (i, row) in a.iter_mut().take(n).enumerate() {
        row[i] = 1.0;
    }
    for j in 0..n {
        let column: Vec<_> = hl.iter().map(|r| r[j]).collect();
        let w = forward(&column);
        for i in 0..m {
            a[n + i][j] = w[i];
        }
    }
    let scales: Vec<f64> = (0..n)
        .map(|j| a.iter().fold(0.0_f64, |x, row| x.hypot(row[j])))
        .collect();
    if !scales.iter().all(|v| v.is_finite() && *v > 0.0)
        || !rhs
            .iter()
            .chain(innovation_variances.iter())
            .all(|v| v.is_finite())
    {
        return Err(FilterError::NumericalFailure);
    }
    for row in &mut a {
        for j in 0..n {
            row[j] /= scales[j];
        }
    }
    for k in 0..n {
        let norm = a[k..].iter().fold(0.0_f64, |x, row| x.hypot(row[k]));
        if !norm.is_finite() || norm <= 0.0 {
            return Err(FilterError::NumericalFailure);
        }
        let mut v: Vec<_> = a[k..].iter().map(|r| r[k]).collect();
        v[0] += norm.copysign(v[0]);
        let vnorm = v.iter().fold(0.0_f64, |x, y| x.hypot(*y));
        if !vnorm.is_finite() || vnorm <= 0.0 {
            return Err(FilterError::NumericalFailure);
        }
        for x in &mut v {
            *x /= vnorm;
        }
        for j in k..n {
            let dot = v
                .iter()
                .enumerate()
                .map(|(i, x)| x * a[k + i][j])
                .sum::<f64>();
            for (i, x) in v.iter().enumerate() {
                a[k + i][j] -= 2.0 * x * dot;
            }
        }
        let dot = v
            .iter()
            .enumerate()
            .map(|(i, x)| x * rhs[k + i])
            .sum::<f64>();
        for (i, x) in v.iter().enumerate() {
            rhs[k + i] -= 2.0 * x * dot;
        }
    }
    let nis = rhs[n..].iter().map(|v| v * v).sum::<f64>();
    if !nis.is_finite() {
        return Err(FilterError::NumericalFailure);
    }
    if max_nis_per_row.is_some_and(|limit| nis / m as f64 > limit) {
        return Err(FilterError::InnovationRejected);
    }
    let backsolve = |rhs: &[f64]| {
        let mut x = vec![0.0; n];
        for i in (0..n).rev() {
            x[i] = (rhs[i] - (i + 1..n).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
        }
        x
    };
    let u = backsolve(&rhs[..n]);
    let z: Vec<_> = u.iter().zip(&scales).map(|(u, s)| u / s).collect();
    let mean: Vec<_> = prior
        .mean
        .iter()
        .zip(&lp)
        .map(|(x, row)| x + row.iter().zip(&z).map(|(a, b)| a * b).sum::<f64>())
        .collect();
    let mut inverse_r = zeros(n, n);
    for j in 0..n {
        let mut e = vec![0.0; n];
        e[j] = 1.0;
        let column = backsolve(&e);
        for i in 0..n {
            inverse_r[i][j] = column[i];
        }
    }
    let mut scaled_lp = lp;
    for row in &mut scaled_lp {
        for j in 0..n {
            row[j] /= scales[j];
        }
    }
    let root = multiply(&scaled_lp, &inverse_r);
    let mut covariance = multiply(&root, &transpose(&root));
    symmetrize(&mut covariance);
    if !mean
        .iter()
        .chain(covariance.iter().flatten())
        .all(|v| v.is_finite())
    {
        return Err(FilterError::NumericalFailure);
    }
    cholesky(&covariance)?;
    Ok(FilterUpdate {
        posterior: GaussianState { mean, covariance },
        normalized_innovation_squared: nis,
        innovation_variances,
    })
}
