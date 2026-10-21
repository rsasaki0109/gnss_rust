//! Unweighted common-clock geometry, independent of a measurement/filter prior.
use crate::spp::{Row, least_squares};
use crate::{EcefCoord, Error, SatelliteId};
use std::collections::BTreeSet;

#[derive(Debug, Clone, PartialEq)]
pub struct DilutionOfPrecision {
    /// Inverse HᵀH in ECEF x/y/z and clock-range coordinates, for unit
    /// independent row variance. This is not the RTK posterior covariance.
    pub geometry_covariance: [[f64; 4]; 4],
    pub gdop: f64,
    pub pdop: f64,
    pub tdop: f64,
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
    };
    if ![result.gdop, result.pdop, result.tdop]
        .iter()
        .all(|v| v.is_finite())
    {
        return Err(Error::NonFinite);
    }
    Ok(result)
}
