//! Unweighted common-clock geometry, independent of a measurement/filter prior.
use crate::spp::{Row, least_squares};
use crate::{EcefCoord, Error, GnssSystem, SatelliteId};
use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq)]
pub struct DilutionOfPrecision {
    /// Inverse HᵀH in ECEF x/y/z and clock-range coordinates for unit row
    /// variance. With several clocks this is the leading 4x4 principal block;
    /// full_geometry_covariance retains the complete matrix. Not an RTK prior.
    pub geometry_covariance: [[f64; 4]; 4],
    pub gdop: f64,
    pub pdop: f64,
    pub tdop: f64,
    /// Empty denotes one common clock. Otherwise the covariance has one clock
    /// per listed system, and geometry_covariance is its leading 4x4 block.
    pub clock_systems: Vec<GnssSystem>,
    /// Full ECEF/clock covariance when more than one clock is estimated.
    pub full_geometry_covariance: Option<Vec<Vec<f64>>>,
}

/// Common receiver-clock geometry H=[-LOS_x,-LOS_y,-LOS_z,1]. Use one
/// constellation or already aligned clock datums. No noise/elevation weighting,
/// normal-matrix inversion, diagonal floor or prior position constraint is used.
/// Householder QR avoids squaring the condition number in the factorization.
pub fn dilution_of_precision(
    receiver: EcefCoord,
    satellites: &[(SatelliteId, EcefCoord)],
) -> Result<DilutionOfPrecision, Error> {
    if !(4..=128).contains(&satellites.len())
        || satellites
            .iter()
            .map(|(s, _)| *s)
            .collect::<BTreeSet<_>>()
            .len()
            != satellites.len()
    {
        return Err(Error::InvalidConfiguration);
    }
    if !receiver.is_finite() || satellites.iter().any(|(_, p)| !p.is_finite()) {
        return Err(Error::NonFinite);
    }
    let rows: Vec<_> = satellites
        .iter()
        .map(|(satellite, position)| {
            let delta = [
                position.x - receiver.x,
                position.y - receiver.y,
                position.z - receiver.z,
            ];
            let range = delta[0].hypot(delta[1]).hypot(delta[2]);
            if !range.is_finite() || range <= 0.0 {
                return Err(Error::SingularGeometry);
            }
            Ok(Row {
                satellite: *satellite,
                design: [-delta[0] / range, -delta[1] / range, -delta[2] / range, 1.0],
                residual: 0.0,
                weight: 1.0,
            })
        })
        .collect::<Result<_, _>>()?;
    let (_, covariance) = least_squares(&rows)?;
    let position_trace = (0..3).map(|i| covariance[i][i]).sum::<f64>();
    let time_variance = covariance[3][3];
    if position_trace <= 0.0 || time_variance <= 0.0 {
        return Err(Error::SingularGeometry);
    }
    let result = DilutionOfPrecision {
        geometry_covariance: covariance,
        gdop: (position_trace + time_variance).sqrt(),
        pdop: position_trace.sqrt(),
        tdop: time_variance.sqrt(),
        clock_systems: Vec::new(),
        full_geometry_covariance: None,
    };
    if ![result.gdop, result.pdop, result.tdop]
        .iter()
        .all(|v| v.is_finite())
    {
        return Err(Error::NonFinite);
    }
    Ok(result)
}

/// Geometry with an independent receiver-clock range for each constellation.
/// This matches separate within-system DD blocks; no cross-system clock or
/// inter-system code-bias alignment is assumed. No prior or diagonal floor.
pub fn dilution_of_precision_by_system(
    receiver: EcefCoord,
    satellites: &[(SatelliteId, EcefCoord)],
) -> Result<DilutionOfPrecision, Error> {
    let systems: Vec<_> = satellites
        .iter()
        .map(|(s, _)| s.system())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    if systems.len() == 1 {
        let mut result = dilution_of_precision(receiver, satellites)?;
        result.clock_systems = systems;
        return Ok(result);
    }
    let n = 3 + systems.len();
    let m = satellites.len();
    if systems.is_empty()
        || systems.contains(&GnssSystem::Unknown)
        || m < n
        || m > 128
        || satellites
            .iter()
            .map(|(s, _)| *s)
            .collect::<BTreeSet<_>>()
            .len()
            != m
    {
        return Err(Error::InvalidConfiguration);
    }
    if !receiver.is_finite() || satellites.iter().any(|(_, p)| !p.is_finite()) {
        return Err(Error::NonFinite);
    }
    let mut a = Vec::with_capacity(m);
    for (satellite, p) in satellites {
        let d = [p.x - receiver.x, p.y - receiver.y, p.z - receiver.z];
        let range = d[0].hypot(d[1]).hypot(d[2]);
        if !range.is_finite() || range <= 0.0 {
            return Err(Error::SingularGeometry);
        }
        let mut row = vec![0.0; n];
        for i in 0..3 {
            row[i] = -d[i] / range;
        }
        row[3 + systems
            .binary_search(&satellite.system())
            .map_err(|_| Error::InvalidConfiguration)?] = 1.0;
        a.push(row);
    }
    // Householder QR, using the same rank criterion as the common-clock code
    // solver. Invert R, not H^T H, so the condition number is not squared.
    let scale = (m as f64).sqrt();
    for k in 0..n {
        let norm = a[k..].iter().fold(0.0_f64, |v, r| v.hypot(r[k]));
        if !norm.is_finite() || norm <= scale * 1e-10 {
            return Err(Error::SingularGeometry);
        }
        let mut v: Vec<_> = a[k..].iter().map(|r| r[k]).collect();
        v[0] += norm.copysign(v[0]);
        let vnorm = v.iter().fold(0.0_f64, |a, b| a.hypot(*b));
        for x in &mut v {
            *x /= vnorm;
        }
        let dots: Vec<f64> = (0..n)
            .map(|j| v.iter().enumerate().map(|(i, x)| x * a[k + i][j]).sum())
            .collect();
        for (r, vi) in a[k..].iter_mut().zip(v) {
            for j in k..n {
                r[j] -= 2.0 * vi * dots[j];
            }
        }
    }
    let inverse: Vec<Vec<f64>> = (0..n)
        .map(|column| {
            let mut x = vec![0.0; n];
            for i in (0..n).rev() {
                x[i] = ((if i == column { 1.0 } else { 0.0 })
                    - (i + 1..n).map(|j| a[i][j] * x[j]).sum::<f64>())
                    / a[i][i];
            }
            x
        })
        .collect();
    let covariance: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| (0..n).map(|k| inverse[k][i] * inverse[k][j]).sum())
                .collect()
        })
        .collect();
    let p: f64 = (0..3).map(|i| covariance[i][i]).sum();
    let t: f64 = (3..n).map(|i| covariance[i][i]).sum();
    if !covariance.iter().flatten().all(|v| v.is_finite()) || p <= 0.0 || t <= 0.0 {
        return Err(Error::SingularGeometry);
    }
    Ok(DilutionOfPrecision {
        geometry_covariance: std::array::from_fn(|i| std::array::from_fn(|j| covariance[i][j])),
        gdop: (p + t).sqrt(),
        pdop: p.sqrt(),
        tdop: t.sqrt(),
        clock_systems: systems,
        full_geometry_covariance: Some(covariance),
    })
}
