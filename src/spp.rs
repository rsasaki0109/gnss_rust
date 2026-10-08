//! Initial GPS/QZSS L1 code solver. RTK/PPP and full upstream SPP quality control
//! are not emulated by this implementation.

use crate::coordinates::{ecef_to_enu, ecef_to_geodetic, geometric_distance};
use crate::models::{KlobucharParameters, klobuchar, saastamoinen};
use crate::navigation::NavigationData;
use crate::observation::{Observation, ObservationEpoch};
use crate::{
    EcefCoord, Error, GnssSystem, PositionSolution, SatelliteId, SolutionStatus, constants as c,
};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SppConfig {
    pub max_iterations: usize,
    pub convergence_m: f64,
    pub elevation_mask_rad: f64,
    pub pseudorange_sigma_m: f64,
    pub max_ephemeris_age_s: f64,
    pub use_troposphere: bool,
    pub use_ionosphere: bool,
}
impl Default for SppConfig {
    fn default() -> Self {
        Self {
            max_iterations: 10,
            convergence_m: 1e-4,
            elevation_mask_rad: 10_f64.to_radians(),
            pseudorange_sigma_m: 3.0,
            max_ephemeris_age_s: 14_400.0,
            use_troposphere: true,
            use_ionosphere: true,
        }
    }
}
impl SppConfig {
    pub(crate) fn validate(self) -> Result<(), Error> {
        if self.max_iterations == 0
            || !self.convergence_m.is_finite()
            || self.convergence_m <= 0.0
            || !self.elevation_mask_rad.is_finite()
            || !(0.0..std::f64::consts::FRAC_PI_2).contains(&self.elevation_mask_rad)
            || !self.pseudorange_sigma_m.is_finite()
            || self.pseudorange_sigma_m <= 0.0
            || !self.max_ephemeris_age_s.is_finite()
            || self.max_ephemeris_age_s < 0.0
        {
            return Err(Error::InvalidConfiguration);
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RejectionReason {
    UnsupportedSignal,
    MissingCode,
    MissingEphemeris,
    BelowElevationMask,
}
#[derive(Debug, Clone, PartialEq)]
pub struct SppSolution {
    pub position: PositionSolution,
    pub iterations: usize,
    /// Unweighted post-fit code residual RMS in metres.
    pub residual_rms_m: f64,
    pub residuals_m: Vec<(SatelliteId, f64)>,
    pub rejected: Vec<(SatelliteId, RejectionReason)>,
}

pub(crate) struct Row {
    pub satellite: SatelliteId,
    pub design: [f64; 4],
    pub residual: f64,
    pub weight: f64,
}

type RejectedObservations = Vec<(SatelliteId, RejectionReason)>;

/// Solve one measurement epoch, choosing one exact L1 tracking code per GPS/
/// QZSS satellite. No phase measurements are consumed; flag 6 is not a code epoch.
pub fn solve_epoch(
    epoch: &ObservationEpoch,
    navigation: &NavigationData,
    config: SppConfig,
) -> Result<SppSolution, Error> {
    config.validate()?;
    let (selected, rejected) = select_observations(epoch)?;
    fit_epoch(epoch, config, rejected, |position, clock_m| {
        build_rows(epoch, navigation, &selected, position, clock_m, config)
    })
}

pub(crate) fn select_observations(
    epoch: &ObservationEpoch,
) -> Result<(BTreeMap<SatelliteId, &Observation>, RejectedObservations), Error> {
    if !matches!(epoch.flag, 0 | 1) {
        return Err(Error::InsufficientObservations);
    }
    let mut selected: BTreeMap<SatelliteId, &Observation> = BTreeMap::new();
    let mut rejected = Vec::new();
    for obs in &epoch.observations {
        if !matches!(obs.satellite.system(), GnssSystem::Gps | GnssSystem::Qzss)
            || !obs.tracking_code.starts_with('1')
        {
            rejected.push((obs.satellite, RejectionReason::UnsupportedSignal));
            continue;
        }
        if !obs.pseudorange_m.is_some_and(|p| p.is_finite() && p > 0.0) {
            rejected.push((obs.satellite, RejectionReason::MissingCode));
            continue;
        }
        // Prefer C1C; otherwise keep a deterministic lexicographic code.
        let rank = |obs: &Observation| {
            (
                obs.tracking_code.as_str() != "1C",
                obs.tracking_code.clone(),
            )
        };
        if selected
            .get(&obs.satellite)
            .is_none_or(|previous| rank(obs) < rank(previous))
        {
            selected.insert(obs.satellite, obs);
        }
    }
    if selected.len() < 4 {
        return Err(Error::InsufficientObservations);
    }
    Ok((selected, rejected))
}

/// Shared four-state QR iteration; measurement providers retain their own
/// orbit/clock provenance and errors. Every final residual uses the final state.
pub(crate) fn fit_epoch<E: From<Error>>(
    epoch: &ObservationEpoch,
    config: SppConfig,
    mut rejected: RejectedObservations,
    mut build: impl FnMut(EcefCoord, f64) -> Result<(Vec<Row>, RejectedObservations), E>,
) -> Result<SppSolution, E> {
    let initial = epoch.approximate_position.unwrap_or_default();
    if !initial.is_finite() {
        return Err(Error::NonFinite.into());
    }
    let mut position = initial;
    let mut clock_m = 0.0;
    for iteration in 1..=config.max_iterations {
        let (rows, _) = build(position, clock_m)?;
        let (update, _) = least_squares(&rows)?;
        position = EcefCoord::new(
            position.x + update[0],
            position.y + update[1],
            position.z + update[2],
        );
        clock_m += update[3];
        if !position.is_finite() || !clock_m.is_finite() {
            return Err(Error::NonFinite.into());
        }
        if update[0].hypot(update[1]).hypot(update[2]) <= config.convergence_m
            && update[3].abs() <= config.convergence_m
        {
            let (rows, final_rejected) = build(position, clock_m)?;
            let (_, covariance) = least_squares(&rows)?;
            rejected.extend(final_rejected);
            let residuals_m = rows
                .iter()
                .map(|row| (row.satellite, row.residual))
                .collect();
            let residual_rms_m = (rows
                .iter()
                .map(|row| row.residual * row.residual)
                .sum::<f64>()
                / rows.len() as f64)
                .sqrt();
            return Ok(SppSolution {
                position: PositionSolution {
                    time: epoch.time,
                    status: SolutionStatus::Spp,
                    position_ecef: Some(position),
                    position_covariance: Some(std::array::from_fn(|i| {
                        std::array::from_fn(|j| covariance[i][j])
                    })),
                    num_satellites: rows.len(),
                    receiver_clock_bias: clock_m / c::SPEED_OF_LIGHT,
                    ..PositionSolution::default()
                },
                iterations: iteration,
                residual_rms_m,
                residuals_m,
                rejected,
            });
        }
    }
    Err(Error::NonConvergent.into())
}

fn build_rows(
    epoch: &ObservationEpoch,
    navigation: &NavigationData,
    selected: &BTreeMap<SatelliteId, &Observation>,
    position: EcefCoord,
    clock_m: f64,
    config: SppConfig,
) -> Result<(Vec<Row>, RejectedObservations), Error> {
    let geo = ecef_to_geodetic(position).ok();
    let ionosphere = navigation
        .gps_ionosphere
        .map(|(alpha, beta)| KlobucharParameters { alpha, beta })
        .unwrap_or_default();
    let mut rows = Vec::new();
    let mut rejected = Vec::new();
    for (&satellite, obs) in selected {
        let pseudorange = obs.pseudorange_m.unwrap();
        let approximate_tx = epoch
            .time
            .checked_add_seconds(-pseudorange / c::SPEED_OF_LIGHT)?;
        let eph = match navigation.ephemeris(satellite, approximate_tx, config.max_ephemeris_age_s)
        {
            Ok(eph) => eph,
            Err(Error::MissingEphemeris) => {
                rejected.push((satellite, RejectionReason::MissingEphemeris));
                continue;
            }
            Err(e) => return Err(e),
        };
        let (_, initial_satellite_clock) = eph.position_clock(approximate_tx)?;
        let tx = approximate_tx.checked_add_seconds(-initial_satellite_clock)?;
        let (sat_pos, satellite_clock) = eph.position_clock(tx)?;
        let difference = EcefCoord::new(
            sat_pos.x - position.x,
            sat_pos.y - position.y,
            sat_pos.z - position.z,
        );
        let distance = difference.x.hypot(difference.y).hypot(difference.z);
        if !distance.is_finite() || distance <= 0.0 {
            return Err(Error::SingularGeometry);
        }
        let (mut elevation, mut azimuth) = (std::f64::consts::FRAC_PI_2, 0.0);
        if let Some(geo) = geo {
            let enu = ecef_to_enu(difference, geo)?;
            elevation = enu.up.atan2(enu.east.hypot(enu.north));
            azimuth = enu.east.atan2(enu.north);
            if elevation < config.elevation_mask_rad {
                rejected.push((satellite, RejectionReason::BelowElevationMask));
                continue;
            }
        }
        let mut atmosphere = 0.0;
        if let Some(geo) = geo {
            if config.use_troposphere {
                atmosphere += saastamoinen(geo, elevation)?;
            }
            if config.use_ionosphere && elevation >= 0.0 {
                atmosphere += klobuchar(geo, azimuth, elevation, epoch.time.tow(), ionosphere)?;
            }
        }
        let predicted = geometric_distance(sat_pos, position)? + clock_m
            - satellite_clock * c::SPEED_OF_LIGHT
            + eph.clock.tgd * c::SPEED_OF_LIGHT
            + atmosphere;
        // Derivative includes the analytical first-order Sagnac term.
        let design = [
            -difference.x / distance - c::OMEGA_E * sat_pos.y / c::SPEED_OF_LIGHT,
            -difference.y / distance + c::OMEGA_E * sat_pos.x / c::SPEED_OF_LIGHT,
            -difference.z / distance,
            1.0,
        ];
        let weight = elevation.sin().max(0.1) / config.pseudorange_sigma_m;
        let residual = pseudorange - predicted;
        if !residual.is_finite() {
            return Err(Error::NonFinite);
        }
        rows.push(Row {
            satellite,
            design,
            residual,
            weight,
        });
    }
    if rows.len() < 4 {
        return Err(Error::InsufficientObservations);
    }
    Ok((rows, rejected))
}

/// Householder QR for the four-column code problem. Avoids squaring the
/// condition number via normal equations. Covariance uses the configured code
/// variance and elevation weights, not a zero residual variance estimate.
fn least_squares(rows: &[Row]) -> Result<([f64; 4], [[f64; 4]; 4]), Error> {
    let mut a: Vec<[f64; 4]> = rows
        .iter()
        .map(|r| r.design.map(|x| x * r.weight))
        .collect();
    let mut b: Vec<f64> = rows.iter().map(|r| r.residual * r.weight).collect();
    let scale = a.iter().map(|r| r[3] * r[3]).sum::<f64>().sqrt();
    for k in 0..4 {
        let norm = a[k..].iter().fold(0.0_f64, |norm, row| norm.hypot(row[k]));
        if !norm.is_finite() || norm <= scale * 1e-10 {
            return Err(Error::SingularGeometry);
        }
        let mut v: Vec<f64> = a[k..].iter().map(|row| row[k]).collect();
        v[0] += norm.copysign(v[0]);
        let vnorm = v.iter().fold(0.0_f64, |norm, &x| norm.hypot(x));
        for value in &mut v {
            *value /= vnorm;
        }
        let dots: [f64; 4] = std::array::from_fn(|column| {
            v.iter()
                .enumerate()
                .map(|(i, &v)| v * a[k + i][column])
                .sum()
        });
        for (row, &vi) in a[k..].iter_mut().zip(&v) {
            for (column, value) in row.iter_mut().enumerate().skip(k) {
                *value -= 2.0 * vi * dots[column];
            }
        }
        let dot = v
            .iter()
            .enumerate()
            .map(|(i, &v)| v * b[k + i])
            .sum::<f64>();
        for (i, &v) in v.iter().enumerate() {
            b[k + i] -= 2.0 * v * dot;
        }
    }
    let backsolve = |rhs: [f64; 4]| {
        let mut result = [0.0; 4];
        for i in (0..4).rev() {
            result[i] = (rhs[i] - (i + 1..4).map(|j| a[i][j] * result[j]).sum::<f64>()) / a[i][i];
        }
        result
    };
    let update = backsolve([b[0], b[1], b[2], b[3]]);
    let inverse: [[f64; 4]; 4] = std::array::from_fn(|column| {
        backsolve(std::array::from_fn(|i| if i == column { 1.0 } else { 0.0 }))
    });
    let covariance = std::array::from_fn(|i| {
        std::array::from_fn(|j| (0..4).map(|k| inverse[k][i] * inverse[k][j]).sum())
    });
    if !update.iter().all(|x| x.is_finite())
        || !covariance.iter().flatten().all(|x: &f64| x.is_finite())
    {
        return Err(Error::NonFinite);
    }
    Ok((update, covariance))
}
