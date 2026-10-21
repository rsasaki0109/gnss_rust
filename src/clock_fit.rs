//! Weighted QR with independent constellation clock columns.
use crate::Error;
use crate::spp::Row;

pub(crate) struct ClockLinearFit {
    pub mean: Vec<f64>,
    pub covariance: Vec<Vec<f64>>,
    pub design: Vec<Vec<f64>>,
}
pub(crate) fn least_squares_with_clocks(
    rows: &[Row],
    systems: &[crate::GnssSystem],
) -> Result<ClockLinearFit, Error> {
    let n = 3 + systems.len();
    let m = rows.len();
    if m < n {
        return Err(Error::InsufficientObservations);
    }
    let design: Vec<Vec<f64>> = rows
        .iter()
        .map(|r| {
            let mut h = vec![0.0; n];
            h[..3].copy_from_slice(&r.design[..3]);
            h[3 + systems
                .binary_search(&r.satellite.system())
                .expect("collected system")] = 1.0;
            h
        })
        .collect();
    let mut a: Vec<Vec<f64>> = design
        .iter()
        .zip(rows)
        .map(|(h, r)| h.iter().map(|v| v * r.weight).collect())
        .collect();
    let mut b: Vec<f64> = rows.iter().map(|r| r.residual * r.weight).collect();
    let scale = rows.iter().map(|r| r.weight * r.weight).sum::<f64>().sqrt();
    for k in 0..n {
        let norm = a[k..].iter().fold(0.0_f64, |v, row| v.hypot(row[k]));
        if !norm.is_finite() || norm <= scale * 1e-10 {
            return Err(Error::SingularGeometry);
        }
        let mut v: Vec<f64> = a[k..].iter().map(|r| r[k]).collect();
        v[0] += norm.copysign(v[0]);
        let vnorm = v.iter().fold(0.0_f64, |v, x| v.hypot(*x));
        for x in &mut v {
            *x /= vnorm;
        }
        let dots: Vec<f64> = (0..n)
            .map(|j| v.iter().enumerate().map(|(i, x)| x * a[k + i][j]).sum())
            .collect();
        for (row, vi) in a[k..].iter_mut().zip(&v) {
            for j in k..n {
                row[j] -= 2.0 * vi * dots[j];
            }
        }
        let dot = v.iter().enumerate().map(|(i, x)| x * b[k + i]).sum::<f64>();
        for (i, x) in v.iter().enumerate() {
            b[k + i] -= 2.0 * x * dot;
        }
    }
    let backsolve = |rhs: &[f64]| {
        let mut x = vec![0.0; n];
        for i in (0..n).rev() {
            x[i] = (rhs[i] - (i + 1..n).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
        }
        x
    };
    let mean = backsolve(&b[..n]);
    let inverse: Vec<Vec<f64>> = (0..n)
        .map(|j| {
            let mut e = vec![0.0; n];
            e[j] = 1.0;
            backsolve(&e)
        })
        .collect();
    let covariance: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| (0..n).map(|k| inverse[k][i] * inverse[k][j]).sum())
                .collect()
        })
        .collect();
    if !mean
        .iter()
        .chain(covariance.iter().flatten())
        .all(|v| v.is_finite())
    {
        return Err(Error::NonFinite);
    }
    Ok(ClockLinearFit {
        mean,
        covariance,
        design,
    })
}
