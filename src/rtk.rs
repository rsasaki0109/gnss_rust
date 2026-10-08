//! Initial static-baseline GPS L1 C/A FLOAT filter.
//!
//! Ports DD signs, correlated noise, SD ambiguity states and ambiguity
//! transformation from libgnss++ RTK measurement/filter components.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
//! This is a restricted FLOAT path, not the complete upstream RTK processor.

use crate::coordinates::{ecef_to_enu, ecef_to_geodetic, geometric_distance};
use crate::filter::{self, FilterError, GaussianState, Matrix};
use crate::lambda::{self, LambdaConfig, LambdaError, LambdaResult};
use crate::models::saastamoinen;
use crate::navigation::NavigationData;
use crate::observation::{Observation, ObservationEpoch};
use crate::{
    EcefCoord, Error, GnssSystem, GnssTime, PositionSolution, SatelliteId, SolutionStatus,
    constants as c,
};
use std::collections::BTreeMap;
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RtkFloatConfig {
    pub elevation_mask_rad: f64,
    pub pseudorange_sigma_m: f64,
    pub carrier_phase_sigma_m: f64,
    pub initial_position_variance_m2: f64,
    pub initial_ambiguity_variance_cycles2: f64,
    pub position_process_noise_m2_per_s: f64,
    pub ambiguity_process_noise_cycles2_per_s: f64,
    pub max_ephemeris_age_s: f64,
    pub max_phase_gap_s: f64,
    pub max_baseline_m: f64,
    pub max_iterations: usize,
    pub convergence_m: f64,
    pub max_nis_per_row: f64,
    pub use_troposphere: bool,
}
impl Default for RtkFloatConfig {
    fn default() -> Self {
        Self {
            elevation_mask_rad: 10_f64.to_radians(),
            pseudorange_sigma_m: 0.3,
            carrier_phase_sigma_m: 0.002,
            initial_position_variance_m2: 900.0,
            initial_ambiguity_variance_cycles2: 900.0,
            position_process_noise_m2_per_s: 1e-4,
            ambiguity_process_noise_cycles2_per_s: 1e-8,
            max_ephemeris_age_s: 14_400.0,
            max_phase_gap_s: 120.0,
            max_baseline_m: 10_000.0,
            max_iterations: 8,
            convergence_m: 1e-4,
            max_nis_per_row: 100.0,
            use_troposphere: true,
        }
    }
}
impl RtkFloatConfig {
    fn validate(self) -> Result<(), RtkError> {
        let positive = [
            self.pseudorange_sigma_m,
            self.carrier_phase_sigma_m,
            self.initial_position_variance_m2,
            self.initial_ambiguity_variance_cycles2,
            self.max_phase_gap_s,
            self.max_baseline_m,
            self.convergence_m,
            self.max_nis_per_row,
        ];
        let nonnegative = [
            self.position_process_noise_m2_per_s,
            self.ambiguity_process_noise_cycles2_per_s,
            self.max_ephemeris_age_s,
        ];
        if positive.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || nonnegative.iter().any(|v| !v.is_finite() || *v < 0.0)
            || !self.elevation_mask_rad.is_finite()
            || !(0.0..std::f64::consts::FRAC_PI_2).contains(&self.elevation_mask_rad)
            || !(1..=32).contains(&self.max_iterations)
        {
            return Err(RtkError::InvalidConfiguration);
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RtkError {
    InvalidConfiguration,
    InvalidEpoch,
    UnsynchronizedEpochs,
    OutOfOrder,
    DuplicateObservation,
    InsufficientCommonSatellites,
    BaselineTooLong,
    NonConvergent,
    Gnss(Error),
    Filter(FilterError),
}
impl fmt::Display for RtkError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfiguration => f.write_str("invalid RTK FLOAT configuration"),
            Self::InvalidEpoch => f.write_str("RTK FLOAT requires normal/power-failure epochs and unapplied clock offsets of zero"),
            Self::UnsynchronizedEpochs => f.write_str("base and rover epoch times differ by more than 1e-7 seconds"),
            Self::OutOfOrder => f.write_str("RTK FLOAT epochs must increase strictly in GPS time"),
            Self::DuplicateObservation => f.write_str("duplicate satellite/tracking observation"),
            Self::InsufficientCommonSatellites => f.write_str("at least four common GPS L1 C/A code/phase satellites are required"),
            Self::BaselineTooLong => f.write_str("baseline exceeds the configured short-baseline limit"),
            Self::NonConvergent => f.write_str("RTK FLOAT iteration did not converge"),
            Self::Gnss(e) => e.fmt(f), Self::Filter(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for RtkError {}
impl From<Error> for RtkError {
    fn from(e: Error) -> Self {
        Self::Gnss(e)
    }
}
impl From<FilterError> for RtkError {
    fn from(e: FilterError) -> Self {
        Self::Filter(e)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RtkRejection {
    UnsupportedSignal,
    MissingPair,
    MissingCodeOrPhase,
    HalfCycle,
    MissingEphemeris,
    BelowElevationMask,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AmbiguityTransform {
    pub head_state: Vec<f64>,
    pub dd_float_cycles: Vec<f64>,
    pub ambiguity_covariance_cycles2: Matrix,
    pub head_ambiguity_covariance: Matrix,
}

/// Map SD states to reference-minus-satellite DD cycles and covariances. Zero
/// ambiguities remain valid; invalid indices never create a zero placeholder.
pub fn ambiguity_transform(
    state: &GaussianState,
    head_count: usize,
    differences: &[(usize, usize)],
) -> Result<AmbiguityTransform, FilterError> {
    let n = state.mean.len();
    if !(1..=128).contains(&n)
        || head_count > n
        || differences.is_empty()
        || differences.len() > 128
        || state.covariance.len() != n
        || state.covariance.iter().any(|r| r.len() != n)
        || differences.iter().any(|&(r, s)| r >= n || s >= n || r == s)
    {
        return Err(FilterError::InvalidDimensions);
    }
    if !state.mean.iter().all(|v| v.is_finite()) {
        return Err(FilterError::NonFinite);
    }
    filter::cholesky(&state.covariance)?;
    let q = &state.covariance;
    let dd_float_cycles: Vec<_> = differences
        .iter()
        .map(|&(r, s)| state.mean[r] - state.mean[s])
        .collect();
    let count = differences.len();
    let mut ambiguity_covariance_cycles2 = filter::zeros(count, count);
    for (i, &(r, s)) in differences.iter().enumerate() {
        for (j, &(t, u)) in differences.iter().enumerate().skip(i) {
            let value = q[r][t] - q[r][u] - q[s][t] + q[s][u];
            // Match upstream's symmetric upper-triangle construction.
            ambiguity_covariance_cycles2[i][j] = value;
            ambiguity_covariance_cycles2[j][i] = value;
        }
    }
    let head_ambiguity_covariance: Matrix = (0..head_count)
        .map(|i| {
            differences
                .iter()
                .map(|&(r, s)| q[i][r] - q[i][s])
                .collect()
        })
        .collect();
    if !dd_float_cycles
        .iter()
        .chain(ambiguity_covariance_cycles2.iter().flatten())
        .chain(head_ambiguity_covariance.iter().flatten())
        .all(|v| v.is_finite())
    {
        return Err(FilterError::NumericalFailure);
    }
    Ok(AmbiguityTransform {
        head_state: state.mean[..head_count].to_vec(),
        dd_float_cycles,
        ambiguity_covariance_cycles2,
        head_ambiguity_covariance,
    })
}

/// DD covariance within a shared-reference block: ref_variance * 11' + diag(sat).
/// The variances here are already rover-minus-base SD variances (metres squared).
pub fn double_difference_covariance(
    reference_variance: f64,
    satellite_variances: &[f64],
) -> Result<Matrix, FilterError> {
    if satellite_variances.is_empty() || satellite_variances.len() > 128 {
        return Err(FilterError::InvalidDimensions);
    }
    if !reference_variance.is_finite() || !satellite_variances.iter().all(|v| v.is_finite()) {
        return Err(FilterError::NonFinite);
    }
    if reference_variance < 0.0 || satellite_variances.iter().any(|v| *v <= 0.0) {
        return Err(FilterError::InvalidCovariance);
    }
    let covariance: Matrix = satellite_variances
        .iter()
        .enumerate()
        .map(|(i, &v)| {
            (0..satellite_variances.len())
                .map(|j| reference_variance + if i == j { v } else { 0.0 })
                .collect()
        })
        .collect();
    if !covariance.iter().flatten().all(|v| v.is_finite()) {
        return Err(FilterError::NumericalFailure);
    }
    Ok(covariance)
}

#[derive(Debug, Clone, PartialEq)]
pub struct RtkFloatSolution {
    pub position: PositionSolution,
    pub reference_satellite: SatelliteId,
    pub ambiguity_satellites: Vec<SatelliteId>,
    pub ambiguities: AmbiguityTransform,
    pub normalized_innovation_squared: f64,
    pub code_residual_rms_m: f64,
    pub phase_residual_rms_m: f64,
    pub iterations: usize,
    pub reset_ambiguities: Vec<SatelliteId>,
    pub rejected: Vec<(SatelliteId, RtkRejection)>,
    pub(crate) fix_context: FixContext,
}
impl RtkFloatSolution {
    /// Diagnostic candidates only. The filter never changes FLOAT to FIX.
    pub fn integer_candidates(&self, config: LambdaConfig) -> Result<LambdaResult, LambdaError> {
        lambda::search(
            &self.ambiguities.dd_float_cycles,
            &self.ambiguities.ambiguity_covariance_cycles2,
            config,
        )
    }
}

#[derive(Debug, Clone, PartialEq)]
struct SatelliteData {
    satellite: SatelliteId,
    rover_position: EcefCoord,
    base_position: EcefCoord,
    code_sd_m: f64,
    phase_sd_cycles: f64,
    elevation: f64,
    loss_of_lock: bool,
}
/// Immutable measurement snapshot for validating candidates against the same
/// satellites, transmit times and observations as the FLOAT update.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct FixContext {
    pub state: GaussianState,
    pub base_position: EcefCoord,
    pub config: RtkFloatConfig,
    pub time: GnssTime,
    pub reset_satellites: Vec<SatelliteId>,
    data: Vec<SatelliteData>,
    reference: usize,
}
pub(crate) struct FixedResiduals {
    pub phase_m: Vec<f64>,
    pub code_m: Vec<f64>,
    pub noise: Matrix,
    pub code_phase_innovations_m: Vec<f64>,
}
impl FixContext {
    pub fn reference_satellite(&self) -> SatelliteId {
        self.data[self.reference].satellite
    }
    pub fn satellites(&self) -> Vec<SatelliteId> {
        self.data.iter().map(|d| d.satellite).collect()
    }
    pub fn head_prior(&self) -> GaussianState {
        GaussianState {
            mean: self.state.mean[..3].to_vec(),
            covariance: self.state.covariance[..3]
                .iter()
                .map(|r| r[..3].to_vec())
                .collect(),
        }
    }
    pub fn ambiguities(&self) -> Result<AmbiguityTransform, FilterError> {
        let pairs: Vec<_> = (0..self.data.len())
            .filter(|&i| i != self.reference)
            .map(|i| (3 + self.reference, 3 + i))
            .collect();
        ambiguity_transform(&self.state, 3, &pairs)
    }
    pub fn residuals(
        &self,
        baseline: &[f64],
        integers: &[i64],
    ) -> Result<FixedResiduals, RtkError> {
        if baseline.len() != 3 || integers.len() + 1 != self.data.len() {
            return Err(FilterError::InvalidDimensions.into());
        }
        let mut fixed_state = self.state.mean.clone();
        fixed_state[..3].copy_from_slice(baseline);
        let mut pairs = Vec::new();
        for (sat, &integer) in (0..self.data.len())
            .filter(|&i| i != self.reference)
            .zip(integers)
        {
            fixed_state[3 + sat] = fixed_state[3 + self.reference] - integer as f64;
            pairs.push(
                (self.data[self.reference].code_sd_m - self.data[sat].code_sd_m)
                    - ((self.data[self.reference].phase_sd_cycles
                        - self.data[sat].phase_sd_cycles)
                        * c::GPS_L1_WAVELENGTH
                        - integer as f64 * c::GPS_L1_WAVELENGTH),
            );
        }
        let system = build_system(
            &self.data,
            self.reference,
            &fixed_state,
            self.base_position,
            self.config,
        )?;
        let rows = self.data.len() - 1;
        Ok(FixedResiduals {
            phase_m: system.residuals[..rows].to_vec(),
            code_m: system.residuals[rows..].to_vec(),
            noise: system.noise,
            code_phase_innovations_m: pairs,
        })
    }
}
type RejectedSatellites = Vec<(SatelliteId, RtkRejection)>;
struct MeasurementSystem {
    design: Matrix,
    residuals: Vec<f64>,
    noise: Matrix,
}

/// Static baseline and one SD ambiguity per tracked GPS 1C satellite. Satellite
/// loss, LLI loss of lock, power failure, a long gap or explicit invalidation
/// starts new ambiguity arcs. Unflagged slips and kinematic motion are pending.
pub struct RtkFloatFilter {
    base_position: EcefCoord,
    initial_rover_position: EcefCoord,
    config: RtkFloatConfig,
    posterior: Option<GaussianState>,
    satellites: Vec<SatelliteId>,
    last_time: Option<GnssTime>,
    continuity_lost: bool,
}
impl RtkFloatFilter {
    pub fn new(
        base_position: EcefCoord,
        initial_rover_position: EcefCoord,
        config: RtkFloatConfig,
    ) -> Result<Self, RtkError> {
        config.validate()?;
        ecef_to_geodetic(base_position)?;
        ecef_to_geodetic(initial_rover_position)?;
        if distance(base_position, initial_rover_position) > config.max_baseline_m {
            return Err(RtkError::BaselineTooLong);
        }
        Ok(Self {
            base_position,
            initial_rover_position,
            config,
            posterior: None,
            satellites: Vec::new(),
            last_time: None,
            continuity_lost: false,
        })
    }
    pub fn state(&self) -> Option<&GaussianState> {
        self.posterior.as_ref()
    }
    pub fn tracked_satellites(&self) -> &[SatelliteId] {
        &self.satellites
    }
    /// Call for skipped event records or externally detected slips. Numerical
    /// state is retained, but all ambiguity arcs restart at the next valid pair.
    pub fn invalidate_ambiguities(&mut self) {
        self.continuity_lost = true;
    }
    pub fn process_pair(
        &mut self,
        rover: &ObservationEpoch,
        base: &ObservationEpoch,
        nav: &NavigationData,
    ) -> Result<RtkFloatSolution, RtkError> {
        let result = self.process_pair_inner(rover, base, nav);
        if result.is_err() {
            self.continuity_lost = true;
        }
        result
    }
    fn process_pair_inner(
        &mut self,
        rover: &ObservationEpoch,
        base: &ObservationEpoch,
        nav: &NavigationData,
    ) -> Result<RtkFloatSolution, RtkError> {
        if !matches!(rover.flag, 0 | 1)
            || !matches!(base.flag, 0 | 1)
            || rover.receiver_clock_offset_s.is_some_and(|v| v != 0.0)
            || base.receiver_clock_offset_s.is_some_and(|v| v != 0.0)
        {
            return Err(RtkError::InvalidEpoch);
        }
        if rover.time.difference_seconds(base.time).abs() > 1e-7 {
            return Err(RtkError::UnsynchronizedEpochs);
        }
        let dt = self
            .last_time
            .map(|t| rover.time.difference_seconds(t))
            .unwrap_or(0.0);
        if self.last_time.is_some() && dt <= 0.0 {
            return Err(RtkError::OutOfOrder);
        }
        let position = self
            .posterior
            .as_ref()
            .map(|s| add(self.base_position, &s.mean[..3]))
            .unwrap_or(self.initial_rover_position);
        let (data, rejected) =
            prepare(rover, base, nav, position, self.base_position, self.config)?;
        if data.len() < 4 {
            return Err(RtkError::InsufficientCommonSatellites);
        }
        let reference = data.iter().enumerate().fold(0, |best, (i, d)| {
            if d.elevation > data[best].elevation {
                i
            } else {
                best
            }
        });
        let reset_all = self.continuity_lost
            || rover.flag == 1
            || base.flag == 1
            || dt > self.config.max_phase_gap_s;
        let (prior, resets) = self.prior(&data, dt, reset_all);
        let mut estimate = prior.mean.clone();
        let mut accepted = None;
        for iteration in 1..=self.config.max_iterations {
            let system =
                build_system(&data, reference, &estimate, self.base_position, self.config)?;
            validate_geometry(&system.design[..data.len() - 1])?;
            // Iterated EKF: relinearize with the SAME prior. Do not count one
            // epoch's observations again as independent measurements.
            let innovation: Vec<_> = system
                .residuals
                .iter()
                .zip(&system.design)
                .map(|(&v, row)| {
                    v + row
                        .iter()
                        .zip(estimate.iter().zip(&prior.mean))
                        .map(|(h, (x, x0))| h * (x - x0))
                        .sum::<f64>()
                })
                .collect();
            let update = filter::measurement_update(
                &prior,
                &system.design,
                &innovation,
                &system.noise,
                Some(self.config.max_nis_per_row),
            )?;
            let change = distance(
                EcefCoord::new(estimate[0], estimate[1], estimate[2]),
                EcefCoord::new(
                    update.posterior.mean[0],
                    update.posterior.mean[1],
                    update.posterior.mean[2],
                ),
            );
            estimate = update.posterior.mean.clone();
            if distance(self.base_position, add(self.base_position, &estimate[..3]))
                > self.config.max_baseline_m
            {
                return Err(RtkError::BaselineTooLong);
            }
            if change <= self.config.convergence_m {
                accepted = Some((update, iteration));
                break;
            }
        }
        let (update, iterations) = accepted.ok_or(RtkError::NonConvergent)?;
        let posterior = update.posterior;
        let final_system = build_system(
            &data,
            reference,
            &posterior.mean,
            self.base_position,
            self.config,
        )?;
        let rows = data.len() - 1;
        let ambiguity_satellites: Vec<_> = data
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != reference)
            .map(|(_, d)| d.satellite)
            .collect();
        let differences: Vec<_> = (0..data.len())
            .filter(|&i| i != reference)
            .map(|i| (3 + reference, 3 + i))
            .collect();
        let ambiguities = ambiguity_transform(&posterior, 3, &differences)?;
        let result = RtkFloatSolution {
            position: PositionSolution {
                time: rover.time,
                status: SolutionStatus::Float,
                position_ecef: Some(add(self.base_position, &posterior.mean[..3])),
                position_covariance: Some(std::array::from_fn(|i| {
                    std::array::from_fn(|j| posterior.covariance[i][j])
                })),
                num_satellites: data.len(),
                ..PositionSolution::default()
            },
            reference_satellite: data[reference].satellite,
            ambiguity_satellites,
            ambiguities,
            normalized_innovation_squared: update.normalized_innovation_squared,
            code_residual_rms_m: rms(&final_system.residuals[rows..]),
            phase_residual_rms_m: rms(&final_system.residuals[..rows]),
            iterations,
            reset_ambiguities: resets.clone(),
            rejected,
            fix_context: FixContext {
                state: posterior.clone(),
                base_position: self.base_position,
                config: self.config,
                time: rover.time,
                reset_satellites: resets.clone(),
                data,
                reference,
            },
        };
        self.posterior = Some(posterior);
        self.satellites = result.fix_context.satellites();
        self.last_time = Some(rover.time);
        self.continuity_lost = false;
        Ok(result)
    }
    fn prior(
        &self,
        data: &[SatelliteData],
        dt: f64,
        reset_all: bool,
    ) -> (GaussianState, Vec<SatelliteId>) {
        let n = 3 + data.len();
        let mut mean = vec![0.0; n];
        let mut covariance = filter::zeros(n, n);
        let mut old_indices: Vec<Option<usize>> = vec![None; n];
        let initial = [
            self.initial_rover_position.x - self.base_position.x,
            self.initial_rover_position.y - self.base_position.y,
            self.initial_rover_position.z - self.base_position.z,
        ];
        for i in 0..3 {
            mean[i] = initial[i];
            covariance[i][i] = self.config.initial_position_variance_m2;
            old_indices[i] = Some(i);
        }
        let mut resets = Vec::new();
        for (i, d) in data.iter().enumerate() {
            let index = 3 + i;
            mean[index] = d.phase_sd_cycles - d.code_sd_m / c::GPS_L1_WAVELENGTH;
            covariance[index][index] = self.config.initial_ambiguity_variance_cycles2;
            if !reset_all && !d.loss_of_lock {
                old_indices[index] = self
                    .satellites
                    .iter()
                    .position(|s| *s == d.satellite)
                    .map(|i| 3 + i);
            }
            if old_indices[index].is_none() {
                resets.push(d.satellite);
            }
        }
        if let Some(old) = &self.posterior {
            for i in 0..n {
                if let Some(oi) = old_indices[i] {
                    mean[i] = old.mean[oi];
                    for j in 0..n {
                        if let Some(oj) = old_indices[j] {
                            covariance[i][j] = old.covariance[oi][oj];
                        }
                    }
                }
            }
        }
        for (i, row) in covariance.iter_mut().enumerate() {
            row[i] += dt
                * if i < 3 {
                    self.config.position_process_noise_m2_per_s
                } else {
                    self.config.ambiguity_process_noise_cycles2_per_s
                };
        }
        (GaussianState { mean, covariance }, resets)
    }
}

fn distance(a: EcefCoord, b: EcefCoord) -> f64 {
    (a.x - b.x).hypot(a.y - b.y).hypot(a.z - b.z)
}
fn add(a: EcefCoord, b: &[f64]) -> EcefCoord {
    EcefCoord::new(a.x + b[0], a.y + b[1], a.z + b[2])
}
fn rms(a: &[f64]) -> f64 {
    (a.iter().map(|v| v * v).sum::<f64>() / a.len() as f64).sqrt()
}
fn elevation(sat: EcefCoord, receiver: EcefCoord) -> Result<f64, Error> {
    let enu = ecef_to_enu(
        EcefCoord::new(sat.x - receiver.x, sat.y - receiver.y, sat.z - receiver.z),
        ecef_to_geodetic(receiver)?,
    )?;
    Ok(enu.up.atan2(enu.east.hypot(enu.north)))
}
fn selected(epoch: &ObservationEpoch) -> Result<BTreeMap<SatelliteId, &Observation>, RtkError> {
    let mut seen = std::collections::BTreeSet::new();
    let mut out = BTreeMap::new();
    for obs in &epoch.observations {
        if !seen.insert((obs.satellite, obs.tracking_code.as_str())) {
            return Err(RtkError::DuplicateObservation);
        }
        if obs.satellite.system() == GnssSystem::Gps && obs.tracking_code == "1C" {
            out.insert(obs.satellite, obs);
        }
    }
    Ok(out)
}
fn satellite_position(
    obs: &Observation,
    time: GnssTime,
    nav: &NavigationData,
    age: f64,
) -> Result<EcefCoord, Error> {
    let tx = time.checked_add_seconds(-obs.pseudorange_m.unwrap() / c::SPEED_OF_LIGHT)?;
    let eph = nav.ephemeris(obs.satellite, tx, age)?;
    let (_, clock) = eph.position_clock(tx)?;
    let (p, _) = eph.position_clock(tx.checked_add_seconds(-clock)?)?;
    Ok(p)
}
fn prepare(
    rover: &ObservationEpoch,
    base: &ObservationEpoch,
    nav: &NavigationData,
    position: EcefCoord,
    base_position: EcefCoord,
    config: RtkFloatConfig,
) -> Result<(Vec<SatelliteData>, RejectedSatellites), RtkError> {
    let r = selected(rover)?;
    let b = selected(base)?;
    let mut data = Vec::new();
    let mut rejected = Vec::new();
    let ids: std::collections::BTreeSet<_> = rover
        .satellites()
        .union(&base.satellites())
        .copied()
        .collect();
    for satellite in ids {
        if satellite.system() != GnssSystem::Gps {
            rejected.push((satellite, RtkRejection::UnsupportedSignal));
            continue;
        }
        let (Some(ro), Some(ba)) = (r.get(&satellite), b.get(&satellite)) else {
            rejected.push((satellite, RtkRejection::MissingPair));
            continue;
        };
        if [*ro, *ba].iter().any(|o| {
            !o.pseudorange_m.is_some_and(|v| v.is_finite() && v > 0.0)
                || !o.carrier_phase_cycles.is_some_and(|v| v.is_finite())
        }) {
            rejected.push((satellite, RtkRejection::MissingCodeOrPhase));
            continue;
        }
        if ro.half_cycle_ambiguity() || ba.half_cycle_ambiguity() {
            rejected.push((satellite, RtkRejection::HalfCycle));
            continue;
        }
        let (rp, bp) = match (
            satellite_position(ro, rover.time, nav, config.max_ephemeris_age_s),
            satellite_position(ba, base.time, nav, config.max_ephemeris_age_s),
        ) {
            (Ok(rp), Ok(bp)) => (rp, bp),
            (Err(Error::MissingEphemeris), _) | (_, Err(Error::MissingEphemeris)) => {
                rejected.push((satellite, RtkRejection::MissingEphemeris));
                continue;
            }
            (Err(e), _) | (_, Err(e)) => return Err(e.into()),
        };
        let el = elevation(rp, position)?;
        if el < config.elevation_mask_rad
            || elevation(bp, base_position)? < config.elevation_mask_rad
        {
            rejected.push((satellite, RtkRejection::BelowElevationMask));
            continue;
        }
        data.push(SatelliteData {
            satellite,
            rover_position: rp,
            base_position: bp,
            code_sd_m: ro.pseudorange_m.unwrap() - ba.pseudorange_m.unwrap(),
            phase_sd_cycles: ro.carrier_phase_cycles.unwrap() - ba.carrier_phase_cycles.unwrap(),
            elevation: el,
            loss_of_lock: ro.loss_of_lock() || ba.loss_of_lock(),
        });
    }
    Ok((data, rejected))
}
fn build_system(
    data: &[SatelliteData],
    reference: usize,
    state: &[f64],
    base: EcefCoord,
    config: RtkFloatConfig,
) -> Result<MeasurementSystem, RtkError> {
    let rover = add(base, &state[..3]);
    let geo_r = ecef_to_geodetic(rover)?;
    let geo_b = ecef_to_geodetic(base)?;
    let mut geometry = Vec::new();
    let mut los = Vec::new();
    let mut elevations = Vec::new();
    for d in data {
        let el = elevation(d.rover_position, rover)?;
        let bel = elevation(d.base_position, base)?;
        let trop = if config.use_troposphere {
            saastamoinen(geo_r, el)? - saastamoinen(geo_b, bel)?
        } else {
            0.0
        };
        geometry.push(
            geometric_distance(d.rover_position, rover)?
                - geometric_distance(d.base_position, base)?
                + trop,
        );
        let diff = [
            d.rover_position.x - rover.x,
            d.rover_position.y - rover.y,
            d.rover_position.z - rover.z,
        ];
        let norm = diff[0].hypot(diff[1]).hypot(diff[2]);
        if !norm.is_finite() || norm <= 0.0 {
            return Err(Error::SingularGeometry.into());
        }
        los.push(diff.map(|v| v / norm));
        elevations.push(el);
    }
    let indices: Vec<_> = (0..data.len()).filter(|&i| i != reference).collect();
    let count = indices.len();
    let mut design = filter::zeros(2 * count, state.len());
    let mut residuals = vec![0.0; 2 * count];
    let mut noise = filter::zeros(2 * count, 2 * count);
    for (row, &sat) in indices.iter().enumerate() {
        for j in 0..3 {
            design[row][j] = -los[reference][j] + los[sat][j];
            design[count + row][j] = design[row][j];
        }
        design[row][3 + reference] = c::GPS_L1_WAVELENGTH;
        design[row][3 + sat] = -c::GPS_L1_WAVELENGTH;
        let geom = geometry[reference] - geometry[sat];
        residuals[row] = (data[reference].phase_sd_cycles - data[sat].phase_sd_cycles)
            * c::GPS_L1_WAVELENGTH
            - geom
            - c::GPS_L1_WAVELENGTH * (state[3 + reference] - state[3 + sat]);
        residuals[count + row] = data[reference].code_sd_m - data[sat].code_sd_m - geom;
    }
    for (block, sigma) in [config.carrier_phase_sigma_m, config.pseudorange_sigma_m]
        .into_iter()
        .enumerate()
    {
        let variance = |el: f64| {
            let sin = el.sin().max(0.1);
            2.0 * sigma * sigma * (1.0 + 1.0 / (sin * sin))
        };
        let cov = double_difference_covariance(
            variance(elevations[reference]),
            &indices
                .iter()
                .map(|&i| variance(elevations[i]))
                .collect::<Vec<_>>(),
        )?;
        for i in 0..count {
            for j in 0..count {
                noise[block * count + i][block * count + j] = cov[i][j];
            }
        }
    }
    if !design
        .iter()
        .flatten()
        .chain(&residuals)
        .chain(noise.iter().flatten())
        .all(|v| v.is_finite())
    {
        return Err(Error::NonFinite.into());
    }
    Ok(MeasurementSystem {
        design,
        residuals,
        noise,
    })
}
fn validate_geometry(design: &[Vec<f64>]) -> Result<(), Error> {
    let mut columns: Vec<Vec<f64>> = (0..3)
        .map(|j| design.iter().map(|r| r[j]).collect())
        .collect();
    for i in 0..3 {
        let (previous, remaining) = columns.split_at_mut(i);
        let column = &mut remaining[0];
        for previous_column in previous.iter() {
            let dot: f64 = column.iter().zip(previous_column).map(|(a, b)| a * b).sum();
            for (value, previous_value) in column.iter_mut().zip(previous_column) {
                *value -= dot * previous_value;
            }
        }
        let norm = column.iter().fold(0.0_f64, |norm, &v| norm.hypot(v));
        if !norm.is_finite() || norm < 1e-8 {
            return Err(Error::SingularGeometry);
        }
        for v in column {
            *v /= norm;
        }
    }
    Ok(())
}
