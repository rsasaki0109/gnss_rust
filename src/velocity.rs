//! Doppler velocity/clock drift, with optional independent constellation clocks.
//! The single-system model follows libgnss++ spp_velocity.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
use crate::clock_fit::{ClockLinearFit as VelocityLinearFit, least_squares_with_clocks};
use crate::constants::SPEED_OF_LIGHT;
use crate::spp::{Row, least_squares};
use crate::{EcefCoord, Error, SatelliteId};
use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VelocityClockModel {
    SingleSystem,
    /// One independent receiver clock-range-rate per observed constellation.
    PerSystem,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DopplerObservation {
    pub satellite: SatelliteId,
    /// Position and velocity already rotated into the reception ECEF frame.
    pub satellite_position: EcefCoord,
    pub satellite_velocity_mps: [f64; 3],
    pub satellite_clock_drift_sps: f64,
    pub doppler_hz: f64,
    pub frequency_hz: f64,
    pub elevation_rad: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VelocityConfig {
    pub clock_model: VelocityClockModel,
    pub doppler_sigma_mps: f64,
    pub min_satellites: usize,
    pub max_chi_square_per_dof: f64,
    pub max_unscaled_velocity_variance_m2s2: f64,
}
impl Default for VelocityConfig {
    fn default() -> Self {
        Self {
            clock_model: VelocityClockModel::SingleSystem,
            doppler_sigma_mps: 0.5,
            min_satellites: 4,
            max_chi_square_per_dof: 25.0,
            max_unscaled_velocity_variance_m2s2: 400.0,
        }
    }
}
impl VelocityConfig {
    fn validate(self) -> Result<(), VelocityError> {
        if !(4..=128).contains(&self.min_satellites)
            || [
                self.doppler_sigma_mps,
                self.max_chi_square_per_dof,
                self.max_unscaled_velocity_variance_m2s2,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(VelocityError::InvalidConfiguration);
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct VelocitySolution {
    pub velocity_mps: [f64; 3],
    /// First system's clock drift; all clocks are named below.
    pub receiver_clock_drift_sps: f64,
    pub clock_systems: Vec<crate::GnssSystem>,
    pub receiver_clock_drifts_sps: Vec<f64>,
    /// Covariance in [vx, vy, vz, clock-range rate] coordinates, (m/s)².
    /// Leading 4x4 principal block with multiple clocks. Full matrices below
    /// retain all velocity/clock and inter-clock cross terms. No position prior.
    pub covariance_m2s2: [[f64; 4]; 4],
    pub unscaled_covariance_m2s2: [[f64; 4]; 4],
    pub full_covariance_m2s2: Option<Vec<Vec<f64>>>,
    pub full_unscaled_covariance_m2s2: Option<Vec<Vec<f64>>>,
    pub satellites: Vec<SatelliteId>,
    pub residuals_mps: Vec<f64>,
    pub residual_rms_mps: f64,
    pub chi_square_per_dof: f64,
    pub covariance_scale: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VelocityError {
    InvalidConfiguration,
    InvalidObservation,
    DuplicateObservation,
    MultipleClockSystems,
    InsufficientObservations,
    InconsistentFit { chi_square_per_dof: f64 },
    WeakGeometry { maximum_unscaled_variance_m2s2: f64 },
    Basic(Error),
}
impl fmt::Display for VelocityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfiguration => f.write_str("invalid Doppler velocity configuration"),
            Self::InvalidObservation => f.write_str("invalid prepared Doppler observation"),
            Self::DuplicateObservation => f.write_str("duplicate Doppler satellite"),
            Self::MultipleClockSystems => {
                f.write_str("single-system Doppler velocity cannot mix receiver clocks")
            }
            Self::InsufficientObservations => f.write_str("insufficient Doppler observations"),
            Self::InconsistentFit { .. } => f.write_str("Doppler fit exceeds its chi-square bound"),
            Self::WeakGeometry { .. } => {
                f.write_str("Doppler geometry exceeds its velocity variance bound")
            }
            Self::Basic(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for VelocityError {}
impl From<Error> for VelocityError {
    fn from(error: Error) -> Self {
        Self::Basic(error)
    }
}

/// Weighted LS of `[-LOS, 1]` and
/// `-Doppler*c/f - satellite_velocity·LOS + c*satellite_clock_drift`.
/// Inputs are already admitted/prepared; invalid
/// rows fail explicitly. A measured zero Doppler is valid, unlike the native
/// numeric zero sentinel. SingleSystem is the default; PerSystem estimates
/// independent clock-range rates without aligning inter-system clock datums.
pub fn solve_velocity(
    observations: &[DopplerObservation],
    receiver: EcefCoord,
    config: VelocityConfig,
) -> Result<VelocitySolution, VelocityError> {
    config.validate()?;
    if !receiver.is_finite() {
        return Err(Error::NonFinite.into());
    }
    if observations.len() < config.min_satellites {
        return Err(VelocityError::InsufficientObservations);
    }
    if observations.len() > 128 {
        return Err(VelocityError::InvalidObservation);
    }
    let systems: BTreeSet<_> = observations.iter().map(|o| o.satellite.system()).collect();
    if systems.len() != 1 && config.clock_model == VelocityClockModel::SingleSystem {
        return Err(VelocityError::MultipleClockSystems);
    }
    if systems.contains(&crate::GnssSystem::Unknown) {
        return Err(VelocityError::InvalidObservation);
    }
    if observations.len() < 3 + systems.len() {
        return Err(VelocityError::InsufficientObservations);
    }
    let mut seen = BTreeSet::new();
    let mut rows = Vec::with_capacity(observations.len());
    for o in observations {
        if !seen.insert(o.satellite) {
            return Err(VelocityError::DuplicateObservation);
        }
        if !o.satellite_position.is_finite()
            || !o
                .satellite_velocity_mps
                .iter()
                .chain([
                    &o.satellite_clock_drift_sps,
                    &o.doppler_hz,
                    &o.frequency_hz,
                    &o.elevation_rad,
                ])
                .all(|v| v.is_finite())
            || o.frequency_hz <= 0.0
            || o.elevation_rad.abs() > std::f64::consts::FRAC_PI_2
        {
            return Err(VelocityError::InvalidObservation);
        }
        let d = [
            o.satellite_position.x - receiver.x,
            o.satellite_position.y - receiver.y,
            o.satellite_position.z - receiver.z,
        ];
        let range = d[0].hypot(d[1]).hypot(d[2]);
        if !range.is_finite() || range < 1.0 {
            return Err(VelocityError::InvalidObservation);
        }
        let los = d.map(|v| v / range);
        let design = [-los[0], -los[1], -los[2], 1.0];
        let residual = -o.doppler_hz * (SPEED_OF_LIGHT / o.frequency_hz)
            - o.satellite_velocity_mps
                .iter()
                .zip(los)
                .map(|(v, l)| v * l)
                .sum::<f64>()
            + SPEED_OF_LIGHT * o.satellite_clock_drift_sps;
        let weight = o.elevation_rad.sin().max(0.1) / config.doppler_sigma_mps;
        if !residual.is_finite() || !weight.is_finite() || weight <= 0.0 {
            return Err(Error::NonFinite.into());
        }
        rows.push(Row {
            satellite: o.satellite,
            design,
            residual,
            weight,
        });
    }
    let systems: Vec<_> = systems.into_iter().collect();
    let fit = if systems.len() == 1 {
        let (mean, covariance) = least_squares(&rows)?;
        VelocityLinearFit {
            mean: mean.to_vec(),
            covariance: covariance.iter().map(|r| r.to_vec()).collect(),
            design: rows.iter().map(|r| r.design.to_vec()).collect(),
        }
    } else {
        least_squares_with_clocks(&rows, &systems)?
    };
    let mean = fit.mean;
    let unscaled = fit.covariance;
    let residuals: Vec<_> = rows
        .iter()
        .zip(&fit.design)
        .map(|(r, h)| r.residual - h.iter().zip(&mean).map(|(h, x)| h * x).sum::<f64>())
        .collect();
    let weighted_sum: f64 = residuals
        .iter()
        .zip(&rows)
        .map(|(v, r)| (v * r.weight).powi(2))
        .sum();
    let chi = weighted_sum / (observations.len() - mean.len()).max(1) as f64;
    if !chi.is_finite() || chi > config.max_chi_square_per_dof {
        return Err(VelocityError::InconsistentFit {
            chi_square_per_dof: chi,
        });
    }
    let maximum = (0..3).map(|i| unscaled[i][i]).fold(0.0, f64::max);
    if maximum > config.max_unscaled_velocity_variance_m2s2 {
        return Err(VelocityError::WeakGeometry {
            maximum_unscaled_variance_m2s2: maximum,
        });
    }
    let scale = chi.max(1.0);
    let covariance: Vec<Vec<f64>> = unscaled
        .iter()
        .map(|row| row.iter().map(|v| v * scale).collect())
        .collect();
    let rms = (residuals.iter().map(|v| v * v).sum::<f64>() / residuals.len() as f64).sqrt();
    if !covariance.iter().flatten().all(|v| v.is_finite()) || !rms.is_finite() {
        return Err(Error::NonFinite.into());
    }
    Ok(VelocitySolution {
        velocity_mps: [mean[0], mean[1], mean[2]],
        receiver_clock_drift_sps: mean[3] / SPEED_OF_LIGHT,
        clock_systems: systems,
        receiver_clock_drifts_sps: mean[3..].iter().map(|v| v / SPEED_OF_LIGHT).collect(),
        covariance_m2s2: std::array::from_fn(|i| std::array::from_fn(|j| covariance[i][j])),
        unscaled_covariance_m2s2: std::array::from_fn(|i| std::array::from_fn(|j| unscaled[i][j])),
        full_covariance_m2s2: (mean.len() > 4).then_some(covariance),
        full_unscaled_covariance_m2s2: (mean.len() > 4).then_some(unscaled),
        satellites: observations.iter().map(|o| o.satellite).collect(),
        residuals_mps: residuals,
        residual_rms_mps: rms,
        chi_square_per_dof: chi,
        covariance_scale: scale,
    })
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VelocityEpochConfig {
    pub velocity: VelocityConfig,
    /// Optional admission rule: require a finite carrier phase and clear LLI
    /// loss-of-lock/half-cycle bits. Phase values are not fused into velocity.
    /// No lock duration or temporal continuity is inferred from this rule.
    pub require_carrier_lock: bool,
    pub max_ephemeris_age_s: f64,
    pub elevation_mask_rad: f64,
}
impl Default for VelocityEpochConfig {
    fn default() -> Self {
        Self {
            velocity: VelocityConfig::default(),
            require_carrier_lock: false,
            max_ephemeris_age_s: 14400.0,
            elevation_mask_rad: 10_f64.to_radians(),
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DopplerRejection {
    UnsupportedSignal,
    AlternativeTrackingCode,
    MissingDoppler,
    InvalidDoppler,
    MissingCarrierPhase,
    InvalidCarrierPhase,
    CarrierLossOfLock,
    InvalidCode,
    Navigation(Error),
    BelowElevationMask,
}
#[derive(Debug, Clone, PartialEq)]
pub struct PreparedDoppler {
    pub raw: crate::observation::Observation,
    pub measurement: DopplerObservation,
    pub transmit_time: crate::GnssTime,
    pub ephemeris_toe: crate::GnssTime,
    pub reception_frame_rotation_rad: f64,
    /// True only when code is absent, using the native 0.075 s initial guess.
    pub used_missing_code_travel_time: bool,
}
#[derive(Debug, Clone, PartialEq)]
pub struct VelocityEpochReport {
    pub prepared: Vec<PreparedDoppler>,
    pub rejected: Vec<(crate::observation::Observation, DopplerRejection)>,
    pub solution: Result<VelocitySolution, VelocityError>,
}

/// GPS 1C-only preparation from raw observations at a caller-supplied receiver
/// position. No phase, integer or truth value is used. Preserve raw fields and
/// every rejection. Rotate broadcast position AND velocity into reception ECEF
/// using the pinned native convention; two clock-refined transmit-time solves.
pub fn solve_gps_l1_epoch(
    epoch: &crate::observation::ObservationEpoch,
    navigation: &crate::navigation::NavigationData,
    receiver: EcefCoord,
    config: VelocityEpochConfig,
) -> Result<VelocityEpochReport, Error> {
    solve_epoch(epoch, navigation, receiver, config, false)
}

/// GPS 1C plus Galileo E1 (valid Doppler and code preferred, then 1C before
/// 1X). Requires explicit per-system clocks. All selected and alternate raw
/// tracking groups and missing values remain in the report.
pub fn solve_gps_galileo_e1_epoch(
    epoch: &crate::observation::ObservationEpoch,
    navigation: &crate::navigation::NavigationData,
    receiver: EcefCoord,
    config: VelocityEpochConfig,
) -> Result<VelocityEpochReport, Error> {
    if config.velocity.clock_model != VelocityClockModel::PerSystem {
        return Err(Error::InvalidConfiguration);
    }
    solve_epoch(epoch, navigation, receiver, config, true)
}

fn solve_epoch(
    epoch: &crate::observation::ObservationEpoch,
    navigation: &crate::navigation::NavigationData,
    receiver: EcefCoord,
    config: VelocityEpochConfig,
    include_galileo: bool,
) -> Result<VelocityEpochReport, Error> {
    use crate::constants::{GPS_L1_FREQ, OMEGA_E};
    use crate::coordinates::{ecef_to_enu, ecef_to_geodetic};
    if config.velocity.validate().is_err()
        || !config.max_ephemeris_age_s.is_finite()
        || config.max_ephemeris_age_s < 0.0
        || !config.elevation_mask_rad.is_finite()
        || !(0.0..std::f64::consts::FRAC_PI_2).contains(&config.elevation_mask_rad)
        || !matches!(epoch.flag, 0 | 1)
        || epoch.receiver_clock_offset_s.is_some_and(|v| v != 0.0)
    {
        return Err(Error::InvalidConfiguration);
    }
    let geo = ecef_to_geodetic(receiver)?;
    let eligible = |raw: &crate::observation::Observation| {
        (raw.satellite.system() == crate::GnssSystem::Gps && raw.tracking_code == "1C")
            || (include_galileo
                && raw.satellite.system() == crate::GnssSystem::Galileo
                && matches!(raw.tracking_code.as_str(), "1C" | "1X"))
    };
    let mut seen = BTreeSet::new();
    let mut selected = BTreeMap::<SatelliteId, &crate::observation::Observation>::new();
    for raw in epoch.observations.iter().filter(|r| eligible(r)) {
        if !seen.insert((raw.satellite, raw.tracking_code.as_str())) {
            return Err(Error::InvalidConfiguration);
        }
        let rank = |r: &crate::observation::Observation| {
            (
                !r.doppler_hz.is_some_and(f64::is_finite),
                !r.pseudorange_m.is_some_and(|v| v.is_finite() && v > 0.0),
                r.tracking_code != "1C",
            )
        };
        if selected
            .get(&raw.satellite)
            .is_none_or(|old| rank(raw) < rank(old))
        {
            selected.insert(raw.satellite, raw);
        }
    }
    let mut prepared = Vec::new();
    let mut rejected = Vec::new();
    for raw in &epoch.observations {
        let mut reject = |reason| rejected.push((raw.clone(), reason));
        if !eligible(raw) {
            reject(DopplerRejection::UnsupportedSignal);
            continue;
        }
        if selected[&raw.satellite].tracking_code != raw.tracking_code {
            reject(DopplerRejection::AlternativeTrackingCode);
            continue;
        }
        let Some(doppler) = raw.doppler_hz else {
            reject(DopplerRejection::MissingDoppler);
            continue;
        };
        if !doppler.is_finite() {
            reject(DopplerRejection::InvalidDoppler);
            continue;
        }
        if config.require_carrier_lock {
            match raw.carrier_phase_cycles {
                None => {
                    reject(DopplerRejection::MissingCarrierPhase);
                    continue;
                }
                Some(value) if !value.is_finite() => {
                    reject(DopplerRejection::InvalidCarrierPhase);
                    continue;
                }
                Some(_) => {}
            }
            if raw.loss_of_lock() || raw.half_cycle_ambiguity() {
                reject(DopplerRejection::CarrierLossOfLock);
                continue;
            }
        }
        let travel = match raw.pseudorange_m {
            None => 0.075,
            Some(value) if value.is_finite() && value > 0.0 => value / SPEED_OF_LIGHT,
            Some(_) => {
                reject(DopplerRejection::InvalidCode);
                continue;
            }
        };
        let approximate_tx = epoch.time.checked_add_seconds(-travel)?;
        let first = navigation
            .ephemeris(raw.satellite, approximate_tx, config.max_ephemeris_age_s)
            .and_then(|e| e.satellite_state(approximate_tx));
        let first = match first {
            Ok(value) => value,
            Err(error) => {
                reject(DopplerRejection::Navigation(error));
                continue;
            }
        };
        let tx = approximate_tx.checked_add_seconds(-first.clock_bias_s)?;
        let eph = match navigation.ephemeris(raw.satellite, tx, config.max_ephemeris_age_s) {
            Ok(value) => value,
            Err(error) => {
                reject(DopplerRejection::Navigation(error));
                continue;
            }
        };
        let state = match eph.satellite_state(tx) {
            Ok(value) => value,
            Err(error) => {
                reject(DopplerRejection::Navigation(error));
                continue;
            }
        };
        let range = (state.position.x - receiver.x)
            .hypot(state.position.y - receiver.y)
            .hypot(state.position.z - receiver.z);
        let angle = OMEGA_E * range / SPEED_OF_LIGHT;
        let (sin, cos) = angle.sin_cos();
        let rotate = |v: [f64; 3]| [cos * v[0] + sin * v[1], -sin * v[0] + cos * v[1], v[2]];
        let p = rotate([state.position.x, state.position.y, state.position.z]);
        let enu = ecef_to_enu(
            EcefCoord::new(p[0] - receiver.x, p[1] - receiver.y, p[2] - receiver.z),
            geo,
        )?;
        let elevation = enu.up.atan2(enu.east.hypot(enu.north));
        if elevation < config.elevation_mask_rad {
            reject(DopplerRejection::BelowElevationMask);
            continue;
        }
        prepared.push(PreparedDoppler {
            raw: raw.clone(),
            measurement: DopplerObservation {
                satellite: raw.satellite,
                satellite_position: EcefCoord::new(p[0], p[1], p[2]),
                satellite_velocity_mps: rotate(state.velocity_mps),
                satellite_clock_drift_sps: state.clock_drift,
                doppler_hz: doppler,
                frequency_hz: GPS_L1_FREQ,
                elevation_rad: elevation,
            },
            transmit_time: tx,
            ephemeris_toe: eph.toe,
            reception_frame_rotation_rad: angle,
            used_missing_code_travel_time: raw.pseudorange_m.is_none(),
        });
    }
    let measurements: Vec<_> = prepared.iter().map(|p| p.measurement).collect();
    let solution = solve_velocity(&measurements, receiver, config.velocity);
    Ok(VelocityEpochReport {
        prepared,
        rejected,
        solution,
    })
}
