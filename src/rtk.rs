//! Initial static/kinematic GPS L1 C/A and optional Galileo E1 FLOAT filter.
//!
//! Ports DD signs, correlated noise, SD ambiguity states and ambiguity
//! transformation from libgnss++ RTK measurement/filter components.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
//! This is a restricted FLOAT path, not the complete upstream RTK processor.

use crate::coordinates::{ecef_to_enu, ecef_to_geodetic, geometric_distance};
use crate::filter::{self, FilterError, GaussianState, Matrix};
#[cfg(test)]
#[path = "rtk_locked_subset_tests.rs"]
mod locked_subset_tests;
use crate::geometry::{
    DilutionOfPrecision, dilution_of_precision, dilution_of_precision_by_system,
};
use crate::lambda::{self, LambdaConfig, LambdaError, LambdaResult};
use crate::models::saastamoinen;
use crate::navigation::NavigationData;
use crate::observation::{Observation, ObservationEpoch};
use crate::rtk_ambiguity::{
    DatumAlignmentDiagnostic, DatumAlignmentInput, align_single_difference_datums,
};
use crate::rtk_quality::{PrefitResidualScreen, screen_prefit_residuals};
use crate::rtk_slip::{self, RtkSlipConfig, SlipTest};
use crate::spp::{self, SppConfig};
use crate::{
    EcefCoord, Error, GnssSystem, GnssTime, PositionSolution, SatelliteId, SolutionStatus,
    constants as c,
};
use std::collections::BTreeMap;
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RtkPositionMode {
    Static,
    /// Reinitialize position uncertainty/cross-covariance each epoch, while
    /// retaining the SD ambiguity block and its continuity. The previous
    /// position is a linearization seed, not a static position constraint.
    Kinematic,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RtkConstellations {
    Gps,
    /// GPS 1C plus Galileo E1 (1C preferred, otherwise 1X per receiver).
    /// Each system has its own DD reference/receiver-clock datum. Kinematic
    /// mixed-system FIX is opt-in; the default FIX policy remains GPS static.
    GpsGalileo,
}

impl RtkConstellations {
    fn includes(self, system: GnssSystem) -> bool {
        system == GnssSystem::Gps || self == Self::GpsGalileo && system == GnssSystem::Galileo
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RtkFloatConfig {
    pub constellations: RtkConstellations,
    pub position_mode: RtkPositionMode,
    /// Translate retained SD ambiguity means by their system's mean raw
    /// code/phase offset before initializing new/reset arcs. No covariance is
    /// changed. False preserves the legacy ambiguity datum policy.
    pub align_ambiguity_datums: bool,
    /// Optional code-only position seed for each kinematic epoch.
    /// Its covariance is not fused. On SPP failure retain the previous RTK
    /// position (or initial position), and report the error in the solution.
    /// Static mode requires None; None preserves the existing seed policy.
    pub kinematic_spp_seed: Option<SppConfig>,
    /// Explicit SPP seed scope, independent of RTK phase selection. Default
    /// GPS 1C preserves the previous path; mixed uses separate GPS/Galileo clocks.
    pub kinematic_spp_constellations: RtkConstellations,
    /// Optional SD code/phase and Doppler continuity tests. Only accepted
    /// epochs enter their history; detected slips reset that satellite's SD
    /// ambiguity without discarding its current measurements.
    pub slip_detection: Option<RtkSlipConfig>,
    /// Optional absolute DD prefit residual gate, checked once at the initial
    /// position/ambiguity seed. Rejects the entire update; no rows or entries
    /// of the correlated measurement covariance are modified.
    pub max_prefit_residual_m: Option<f64>,
    /// Optional unweighted geometry gate, evaluated before committing a
    /// posterior. None preserves the existing admission policy.
    pub max_geometry_gdop: Option<f64>,
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
            constellations: RtkConstellations::Gps,
            position_mode: RtkPositionMode::Static,
            align_ambiguity_datums: false,
            kinematic_spp_seed: None,
            kinematic_spp_constellations: RtkConstellations::Gps,
            slip_detection: None,
            max_prefit_residual_m: None,
            max_geometry_gdop: None,
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
        if self
            .max_prefit_residual_m
            .is_some_and(|v| !v.is_finite() || v <= 0.0)
            || self
                .slip_detection
                .is_some_and(|config| config.validate().is_err())
            || self.kinematic_spp_seed.is_some_and(|config| {
                self.position_mode != RtkPositionMode::Kinematic || config.validate().is_err()
            })
            || self
                .max_geometry_gdop
                .is_some_and(|v| !v.is_finite() || v <= 0.0)
            || positive.iter().any(|v| !v.is_finite() || *v <= 0.0)
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
    GeometryTooWeak,
    PrefitResidualRejected,
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
            Self::InsufficientCommonSatellites => f.write_str("at least three independent within-system code/phase double differences are required"),
            Self::BaselineTooLong => f.write_str("baseline exceeds the configured short-baseline limit"),
            Self::NonConvergent => f.write_str("RTK FLOAT iteration did not converge"),
            Self::GeometryTooWeak => f.write_str("RTK geometry exceeds the configured GDOP limit"),
            Self::PrefitResidualRejected => f.write_str("RTK initial DD residual exceeds the configured prefit limit"),
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
    InsufficientSystemSatellites,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RtkSeedSource {
    InitialPosition,
    PreviousRtk,
    SinglePoint,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RtkSeedDiagnostic {
    pub source: RtkSeedSource,
    pub position: EcefCoord,
    /// Present only when an attempted SPP seed failed. The seed source then
    /// explains which fallback was used; no SPP output is labeled RTK.
    pub spp_failure: Option<Error>,
    pub spp_requested_constellations: Option<RtkConstellations>,
    pub spp_satellites: usize,
    pub spp_gps_clock_bias_s: Option<f64>,
    pub spp_galileo_clock_bias_s: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RtkPrefitDiagnostic {
    pub time: GnssTime,
    pub position_seed: RtkSeedDiagnostic,
    pub reference_satellite: SatelliteId,
    /// Maps each phase/code screen row to its non-reference DD satellite.
    pub row_satellites: Vec<SatelliteId>,
    pub row_reference_satellites: Vec<SatelliteId>,
    pub screen: PrefitResidualScreen,
    /// Evidence from this attempt, including a rejected update. These raw
    /// observations enter the continuity history only if the update succeeds.
    pub slip_diagnostics: Vec<RtkSlipDiagnostic>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RtkSlipDiagnostic {
    pub satellite: SatelliteId,
    /// Tests are absent across a known continuity break or without history.
    pub code: Option<SlipTest>,
    pub doppler: Option<SlipTest>,
    /// Missing Doppler remains unavailable; a measured zero remains valid.
    pub doppler_available: bool,
    pub doppler_error: Option<Error>,
}

#[derive(Debug, Clone, Copy)]
struct SlipHistory {
    code_minus_phase_m: f64,
    doppler_phase_m: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RtkFloatSolution {
    pub ambiguity_datum_alignment: Vec<DatumAlignmentDiagnostic>,
    pub position: PositionSolution,
    pub position_seed: RtkSeedDiagnostic,
    pub slip_diagnostics: Vec<RtkSlipDiagnostic>,
    pub reference_satellite: SatelliteId,
    /// One reference per participating system; the scalar above is the first.
    pub reference_satellites: Vec<SatelliteId>,
    pub ambiguity_satellites: Vec<SatelliteId>,
    pub ambiguities: AmbiguityTransform,
    pub normalized_innovation_squared: f64,
    pub code_residual_rms_m: f64,
    pub phase_residual_rms_m: f64,
    pub iterations: usize,
    pub reset_ambiguities: Vec<SatelliteId>,
    pub rejected: Vec<(SatelliteId, RtkRejection)>,
    /// Independent common-clock geometry at the final estimated position.
    /// An unavailable diagnostic is retained explicitly when its gate is off.
    pub geometry: Result<DilutionOfPrecision, Error>,
    pub(crate) fix_context: FixContext,
}
impl RtkFloatSolution {
    pub fn tracking_pairs(&self) -> Vec<RtkTrackingPair> {
        self.fix_context
            .data
            .iter()
            .map(|d| RtkTrackingPair {
                satellite: d.satellite,
                rover_tracking_code: d.rover_tracking_code,
                base_tracking_code: d.base_tracking_code,
            })
            .collect()
    }
    pub fn used_satellites(&self) -> Vec<SatelliteId> {
        self.fix_context.satellites()
    }
    /// Diagnostic candidates only. The filter never changes FLOAT to FIX.
    pub fn integer_candidates(&self, config: LambdaConfig) -> Result<LambdaResult, LambdaError> {
        lambda::search(
            &self.ambiguities.dd_float_cycles,
            &self.ambiguities.ambiguity_covariance_cycles2,
            config,
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RtkTrackingPair {
    pub satellite: SatelliteId,
    pub rover_tracking_code: &'static str,
    pub base_tracking_code: &'static str,
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
    rover_tracking_code: &'static str,
    base_tracking_code: &'static str,
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
    references: Vec<usize>,
}
pub(crate) struct FixedResiduals {
    pub phase_m: Vec<f64>,
    pub code_m: Vec<f64>,
    pub noise: Matrix,
    pub code_phase_innovations_m: Vec<f64>,
}
impl FixContext {
    /// Marginal snapshot for lock-qualified satellites. Omitted SD states are
    /// integrated out by taking the full principal covariance block; they are
    /// never assigned zero covariance or treated as fixed. The FLOAT posterior
    /// already used all observations and is not reprocessed here.
    pub(crate) fn locked_subset(
        &self,
        locks: &BTreeMap<SatelliteId, usize>,
        minimum_lock: usize,
    ) -> Result<Option<Self>, RtkError> {
        let qualified = self
            .data
            .iter()
            .filter_map(|d| {
                (locks.get(&d.satellite).copied().unwrap_or(0) >= minimum_lock)
                    .then_some(d.satellite)
            })
            .collect();
        self.satellite_subset(&qualified)
    }
    pub(crate) fn ambiguity_pair_subset(
        &self,
        indices: &[usize],
    ) -> Result<Option<Self>, RtkError> {
        let pairs = dd_pairs(&self.data, &self.references);
        let mut unique = std::collections::BTreeSet::new();
        let mut qualified = std::collections::BTreeSet::new();
        if indices.is_empty() {
            return Err(FilterError::InvalidDimensions.into());
        }
        for &i in indices {
            if i >= pairs.len() || !unique.insert(i) {
                return Err(FilterError::InvalidDimensions.into());
            }
            let (r, s) = pairs[i];
            qualified.insert(self.data[r].satellite);
            qualified.insert(self.data[s].satellite);
        }
        self.satellite_subset(&qualified)
    }
    fn satellite_subset(
        &self,
        qualified: &std::collections::BTreeSet<SatelliteId>,
    ) -> Result<Option<Self>, RtkError> {
        let mut groups = BTreeMap::<_, Vec<usize>>::new();
        for (i, d) in self.data.iter().enumerate() {
            if qualified.contains(&d.satellite) {
                groups.entry(d.satellite.system()).or_default().push(i);
            }
        }
        let mut kept = Vec::new();
        let mut references = Vec::new();
        for members in groups.values().filter(|members| members.len() >= 2) {
            let mut best = members[0];
            for &i in members {
                if self.data[i].elevation > self.data[best].elevation {
                    best = i;
                }
            }
            references.push(kept.len() + members.iter().position(|&i| i == best).unwrap());
            kept.extend_from_slice(members);
        }
        if kept.is_empty() {
            return Ok(None);
        }
        let indices: Vec<_> = (0..3).chain(kept.iter().map(|i| 3 + i)).collect();
        let state = GaussianState {
            mean: indices.iter().map(|&i| self.state.mean[i]).collect(),
            covariance: indices
                .iter()
                .map(|&i| {
                    indices
                        .iter()
                        .map(|&j| self.state.covariance[i][j])
                        .collect()
                })
                .collect(),
        };
        filter::cholesky(&state.covariance)?;
        Ok(Some(Self {
            state,
            data: kept.iter().map(|&i| self.data[i].clone()).collect(),
            reference: references[0],
            references,
            reset_satellites: self
                .reset_satellites
                .iter()
                .copied()
                .filter(|s| kept.iter().any(|&i| self.data[i].satellite == *s))
                .collect(),
            base_position: self.base_position,
            config: self.config,
            time: self.time,
        }))
    }
    pub(crate) fn with_reference(&self, reference: SatelliteId) -> Result<Self, RtkError> {
        if self.config.constellations != RtkConstellations::Gps {
            return Err(Error::UnsupportedConstellation.into());
        }
        self.with_references(&[reference])
    }
    pub(crate) fn with_references(&self, references: &[SatelliteId]) -> Result<Self, RtkError> {
        let systems: std::collections::BTreeSet<_> =
            self.data.iter().map(|d| d.satellite.system()).collect();
        let mut chosen = BTreeMap::new();
        for &satellite in references {
            let index = self
                .data
                .iter()
                .position(|d| d.satellite == satellite)
                .ok_or(RtkError::InvalidConfiguration)?;
            if chosen.insert(satellite.system(), index).is_some() {
                return Err(RtkError::InvalidConfiguration);
            }
        }
        if chosen
            .keys()
            .copied()
            .collect::<std::collections::BTreeSet<_>>()
            != systems
        {
            return Err(RtkError::InvalidConfiguration);
        }
        let mut view = self.clone();
        view.references = chosen.into_values().collect();
        view.reference = *view
            .references
            .first()
            .ok_or(RtkError::InvalidConfiguration)?;
        Ok(view)
    }
    pub fn reference_satellite(&self) -> SatelliteId {
        self.data[self.reference].satellite
    }
    pub fn reference_satellites(&self) -> Vec<SatelliteId> {
        self.references
            .iter()
            .map(|&i| self.data[i].satellite)
            .collect()
    }
    pub fn ambiguity_satellites(&self) -> Vec<SatelliteId> {
        dd_pairs(&self.data, &self.references)
            .iter()
            .map(|(_, s)| self.data[*s].satellite)
            .collect()
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
        let pairs: Vec<_> = dd_pairs(&self.data, &self.references)
            .iter()
            .map(|(r, s)| (3 + r, 3 + s))
            .collect();
        ambiguity_transform(&self.state, 3, &pairs)
    }
    pub fn residuals(
        &self,
        baseline: &[f64],
        integers: &[i64],
    ) -> Result<FixedResiduals, RtkError> {
        let dd = dd_pairs(&self.data, &self.references);
        if baseline.len() != 3 || integers.len() != dd.len() {
            return Err(FilterError::InvalidDimensions.into());
        }
        let mut fixed_state = self.state.mean.clone();
        fixed_state[..3].copy_from_slice(baseline);
        let mut pairs = Vec::new();
        for ((reference, sat), &integer) in dd.iter().copied().zip(integers) {
            fixed_state[3 + sat] = fixed_state[3 + reference] - integer as f64;
            pairs.push(
                (self.data[reference].code_sd_m - self.data[sat].code_sd_m)
                    - ((self.data[reference].phase_sd_cycles - self.data[sat].phase_sd_cycles)
                        * c::GPS_L1_WAVELENGTH
                        - integer as f64 * c::GPS_L1_WAVELENGTH),
            );
        }
        let system = build_system_with_references(
            &self.data,
            &self.references,
            &fixed_state,
            self.base_position,
            self.config,
        )?;
        let rows = dd.len();
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

/// Static or epoch-wise kinematic baseline and one SD ambiguity per selected L1/E1 satellite. Satellite
/// loss, LLI loss of lock, power failure, a long gap or explicit invalidation
/// starts new ambiguity arcs. Optional code/phase and Doppler tests can also
/// reset individual GPS 1C arcs; they do not identify every possible slip.
pub struct RtkFloatFilter {
    base_position: EcefCoord,
    initial_rover_position: EcefCoord,
    config: RtkFloatConfig,
    posterior: Option<GaussianState>,
    satellites: Vec<SatelliteId>,
    last_time: Option<GnssTime>,
    continuity_lost: bool,
    slip_history: BTreeMap<SatelliteId, SlipHistory>,
    last_prefit_diagnostic: Option<RtkPrefitDiagnostic>,
    tracking_pairs: BTreeMap<SatelliteId, (&'static str, &'static str)>,
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
            slip_history: BTreeMap::new(),
            last_prefit_diagnostic: None,
            tracking_pairs: BTreeMap::new(),
        })
    }
    pub fn state(&self) -> Option<&GaussianState> {
        self.posterior.as_ref()
    }
    pub fn tracked_satellites(&self) -> &[SatelliteId] {
        &self.satellites
    }
    /// Diagnostic from the most recent process_pair attempt, including a
    /// prefit rejection. None if the gate is disabled or was not reached.
    /// This diagnostic never commits a numerical state or continuity history.
    pub fn last_prefit_diagnostic(&self) -> Option<&RtkPrefitDiagnostic> {
        self.last_prefit_diagnostic.as_ref()
    }
    /// Call for skipped event records or externally detected slips. Numerical
    /// state is retained, but all ambiguity arcs restart at the next valid pair.
    pub fn invalidate_ambiguities(&mut self) {
        self.continuity_lost = true;
    }
    /// Preserve continuity evidence when an input is not used for a paired
    /// update (for example 5 Hz rover data between 1 Hz common epochs).
    /// A power/slip event or GPS 1C LLI loss/half-cycle conservatively resets
    /// all ambiguity arcs on the next update. No numeric state/time advances.
    pub fn observe_skipped_epoch(&mut self, epoch: &ObservationEpoch) -> bool {
        let broken = matches!(epoch.flag, 1 | 6)
            || epoch.observations.iter().any(|o| {
                self.config.constellations.includes(o.satellite.system())
                    && (o.tracking_code == "1C"
                        || o.satellite.system() == GnssSystem::Galileo && o.tracking_code == "1X")
                    && (o.loss_of_lock() || o.half_cycle_ambiguity())
            });
        if broken {
            self.invalidate_ambiguities();
        }
        broken
    }
    pub fn process_pair(
        &mut self,
        rover: &ObservationEpoch,
        base: &ObservationEpoch,
        nav: &NavigationData,
    ) -> Result<RtkFloatSolution, RtkError> {
        self.last_prefit_diagnostic = None;
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
        let previous_position = self
            .posterior
            .as_ref()
            .map(|s| add(self.base_position, &s.mean[..3]))
            .unwrap_or(self.initial_rover_position);
        let mut position_seed = RtkSeedDiagnostic {
            source: if self.posterior.is_some() {
                RtkSeedSource::PreviousRtk
            } else {
                RtkSeedSource::InitialPosition
            },
            position: previous_position,
            spp_failure: None,
            spp_requested_constellations: self
                .config
                .kinematic_spp_seed
                .map(|_| self.config.kinematic_spp_constellations),
            spp_satellites: 0,
            spp_gps_clock_bias_s: None,
            spp_galileo_clock_bias_s: None,
        };
        let mut spp_position = None;
        if let Some(config) = self.config.kinematic_spp_seed {
            let mut seed_epoch = rover.clone();
            seed_epoch.approximate_position = Some(previous_position);
            let result = match self.config.kinematic_spp_constellations {
                RtkConstellations::Gps => {
                    selected(rover, RtkConstellations::Gps)?;
                    seed_epoch.observations.retain(|o| {
                        o.satellite.system() == GnssSystem::Gps && o.tracking_code == "1C"
                    });
                    spp::solve_epoch(&seed_epoch, nav, config).map(|s| {
                        let clocks = vec![(GnssSystem::Gps, s.position.receiver_clock_bias)];
                        (s, clocks)
                    })
                }
                RtkConstellations::GpsGalileo => {
                    crate::spp_multiclock::solve_gps_galileo_epoch(&seed_epoch, nav, config)
                        .map(|s| (s.spp, s.clock_biases_s))
                }
            };
            match result {
                Ok((solution, clocks)) => {
                    // A successful SPP solve always returns a finite position.
                    let position = solution
                        .position
                        .position_ecef
                        .ok_or(RtkError::InvalidEpoch)?;
                    spp_position = Some(position);
                    position_seed.source = RtkSeedSource::SinglePoint;
                    position_seed.position = position;
                    position_seed.spp_satellites = solution.position.num_satellites;
                    position_seed.spp_gps_clock_bias_s = clocks
                        .iter()
                        .find(|(s, _)| *s == GnssSystem::Gps)
                        .map(|(_, v)| *v);
                    position_seed.spp_galileo_clock_bias_s = clocks
                        .iter()
                        .find(|(s, _)| *s == GnssSystem::Galileo)
                        .map(|(_, v)| *v);
                }
                Err(error) => position_seed.spp_failure = Some(error),
            }
        }
        let (mut data, rejected) = prepare(
            rover,
            base,
            nav,
            position_seed.position,
            self.base_position,
            self.config,
        )?;
        if data.len() < 4 {
            return Err(RtkError::InsufficientCommonSatellites);
        }
        for d in &mut data {
            if self
                .tracking_pairs
                .get(&d.satellite)
                .is_some_and(|old| *old != (d.rover_tracking_code, d.base_tracking_code))
            {
                d.loss_of_lock = true;
            }
        }
        let mut by_system = BTreeMap::<GnssSystem, usize>::new();
        for (i, d) in data.iter().enumerate() {
            let best = by_system.entry(d.satellite.system()).or_insert(i);
            if d.elevation > data[*best].elevation {
                *best = i;
            }
        }
        let references: Vec<_> = by_system.into_values().collect();
        let reference = references[0];
        let pairs = dd_pairs(&data, &references);
        if pairs.len() < 3 {
            return Err(RtkError::InsufficientCommonSatellites);
        }
        let reset_all = self.continuity_lost
            || rover.flag == 1
            || base.flag == 1
            || dt > self.config.max_phase_gap_s;
        let (slip_diagnostics, next_slip_history) =
            self.detect_slips(&mut data, rover, base, dt, reset_all)?;
        let (prior, resets, ambiguity_datum_alignment) =
            self.prior(&data, dt, reset_all, spp_position)?;
        let mut estimate = prior.mean.clone();
        let mut accepted = None;
        for iteration in 1..=self.config.max_iterations {
            let system = build_system_with_references(
                &data,
                &references,
                &estimate,
                self.base_position,
                self.config,
            )?;
            validate_geometry(&system.design[..pairs.len()])?;
            if iteration == 1
                && let Some(threshold) = self.config.max_prefit_residual_m
            {
                let count = pairs.len();
                let screen = screen_prefit_residuals(
                    &system.residuals[..count],
                    &system.residuals[count..],
                    threshold,
                )?;
                let accepted = screen.accepted();
                self.last_prefit_diagnostic = Some(RtkPrefitDiagnostic {
                    time: rover.time,
                    position_seed,
                    reference_satellite: data[reference].satellite,
                    row_satellites: pairs.iter().map(|(_, sat)| data[*sat].satellite).collect(),
                    row_reference_satellites: pairs
                        .iter()
                        .map(|(r, _)| data[*r].satellite)
                        .collect(),
                    screen,
                    slip_diagnostics: slip_diagnostics.clone(),
                });
                if !accepted {
                    return Err(RtkError::PrefitResidualRejected);
                }
            }
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
        let geometry_function = if references.len() == 1 {
            dilution_of_precision
        } else {
            dilution_of_precision_by_system
        };
        let geometry = geometry_function(
            add(self.base_position, &posterior.mean[..3]),
            &data
                .iter()
                .map(|d| (d.satellite, d.rover_position))
                .collect::<Vec<_>>(),
        );
        if let Some(maximum) = self.config.max_geometry_gdop {
            match &geometry {
                Ok(dop) if dop.gdop <= maximum => {}
                Ok(_) => return Err(RtkError::GeometryTooWeak),
                Err(error) => return Err((*error).into()),
            }
        }
        let final_system = build_system_with_references(
            &data,
            &references,
            &posterior.mean,
            self.base_position,
            self.config,
        )?;
        let rows = pairs.len();
        let ambiguity_satellites: Vec<_> =
            pairs.iter().map(|(_, sat)| data[*sat].satellite).collect();
        let differences: Vec<_> = pairs.iter().map(|(r, sat)| (3 + r, 3 + sat)).collect();
        let ambiguities = ambiguity_transform(&posterior, 3, &differences)?;
        let result = RtkFloatSolution {
            ambiguity_datum_alignment,
            position_seed,
            slip_diagnostics,
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
            reference_satellites: references.iter().map(|r| data[*r].satellite).collect(),
            ambiguity_satellites,
            ambiguities,
            normalized_innovation_squared: update.normalized_innovation_squared,
            code_residual_rms_m: rms(&final_system.residuals[rows..]),
            phase_residual_rms_m: rms(&final_system.residuals[..rows]),
            iterations,
            reset_ambiguities: resets.clone(),
            rejected,
            geometry,
            fix_context: FixContext {
                state: posterior.clone(),
                base_position: self.base_position,
                config: self.config,
                time: rover.time,
                reset_satellites: resets.clone(),
                data,
                reference,
                references,
            },
        };
        self.tracking_pairs = result
            .tracking_pairs()
            .into_iter()
            .map(|p| (p.satellite, (p.rover_tracking_code, p.base_tracking_code)))
            .collect();
        self.posterior = Some(posterior);
        self.satellites = result.fix_context.satellites();
        self.last_time = Some(rover.time);
        self.continuity_lost = false;
        self.slip_history = next_slip_history;
        Ok(result)
    }

    fn detect_slips(
        &self,
        data: &mut [SatelliteData],
        rover: &ObservationEpoch,
        base: &ObservationEpoch,
        dt: f64,
        reset_all: bool,
    ) -> Result<(Vec<RtkSlipDiagnostic>, BTreeMap<SatelliteId, SlipHistory>), RtkError> {
        let Some(config) = self.config.slip_detection else {
            return Ok((Vec::new(), BTreeMap::new()));
        };
        let r = selected(rover, self.config.constellations)?;
        let b = selected(base, self.config.constellations)?;
        let mut diagnostics = Vec::with_capacity(data.len());
        let mut next = BTreeMap::new();
        for d in data {
            let code_minus_phase_m = d.code_sd_m - d.phase_sd_cycles * c::GPS_L1_WAVELENGTH;
            let phase_m = d.phase_sd_cycles * c::GPS_L1_WAVELENGTH;
            if !code_minus_phase_m.is_finite() || !phase_m.is_finite() {
                return Err(Error::NonFinite.into());
            }
            let mut diagnostic = RtkSlipDiagnostic {
                satellite: d.satellite,
                code: None,
                doppler: None,
                doppler_available: false,
                doppler_error: None,
            };
            let mut rate = None;
            if config.doppler_threshold_m.is_some()
                && let (Some(rd), Some(bd)) =
                    (r[&d.satellite].doppler_hz, b[&d.satellite].doppler_hz)
            {
                match rtk_slip::single_difference_range_rate_mps(rd, bd, c::GPS_L1_WAVELENGTH) {
                    Ok(value) => {
                        rate = Some(value);
                        diagnostic.doppler_available = true;
                    }
                    Err(error) => diagnostic.doppler_error = Some(error),
                }
            }
            if !reset_all
                && !d.loss_of_lock
                && let Some(previous) = self.slip_history.get(&d.satellite)
            {
                if let Some(threshold) = config.code_threshold_m {
                    diagnostic.code = Some(rtk_slip::test_code_minus_phase(
                        previous.code_minus_phase_m,
                        code_minus_phase_m,
                        threshold,
                    )?);
                }
                if let (Some(threshold), Some(previous_phase), Some(rate)) =
                    (config.doppler_threshold_m, previous.doppler_phase_m, rate)
                    && dt <= config.max_doppler_interval_s
                {
                    diagnostic.doppler = Some(rtk_slip::test_doppler_phase(
                        previous_phase,
                        phase_m,
                        rate,
                        dt,
                        threshold,
                    )?);
                }
            }
            d.loss_of_lock |= diagnostic.code.is_some_and(|v| v.detected)
                || diagnostic.doppler.is_some_and(|v| v.detected);
            next.insert(
                d.satellite,
                SlipHistory {
                    code_minus_phase_m,
                    doppler_phase_m: rate.map(|_| phase_m),
                },
            );
            diagnostics.push(diagnostic);
        }
        Ok((diagnostics, next))
    }
    fn prior(
        &self,
        data: &[SatelliteData],
        dt: f64,
        reset_all: bool,
        spp_position: Option<EcefCoord>,
    ) -> Result<
        (
            GaussianState,
            Vec<SatelliteId>,
            Vec<DatumAlignmentDiagnostic>,
        ),
        RtkError,
    > {
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
                            if self.config.position_mode == RtkPositionMode::Static
                                || (i >= 3 && j >= 3)
                            {
                                covariance[i][j] = old.covariance[oi][oj];
                            }
                        }
                    }
                }
            }
        }
        let mut alignment = Vec::new();
        if self.config.align_ambiguity_datums {
            let inputs: Vec<_> = data
                .iter()
                .enumerate()
                .map(|(i, d)| DatumAlignmentInput {
                    satellite: d.satellite,
                    raw_bias_cycles: d.phase_sd_cycles - d.code_sd_m / c::GPS_L1_WAVELENGTH,
                    retained_cycles: old_indices[3 + i]
                        .and_then(|oi| self.posterior.as_ref().map(|old| old.mean[oi])),
                })
                .collect();
            let aligned = align_single_difference_datums(&inputs)?;
            mean[3..].copy_from_slice(&aligned.single_difference_cycles);
            alignment = aligned.systems;
        }
        if let Some(position) = spp_position {
            // Kinematic position covariance/cross-covariance has already been
            // reset above. Preserve the ambiguity block and its continuity.
            mean[..3].copy_from_slice(&[
                position.x - self.base_position.x,
                position.y - self.base_position.y,
                position.z - self.base_position.z,
            ]);
        }
        for (i, row) in covariance.iter_mut().enumerate() {
            row[i] += dt
                * if i < 3 {
                    self.config.position_process_noise_m2_per_s
                } else {
                    self.config.ambiguity_process_noise_cycles2_per_s
                };
        }
        Ok((GaussianState { mean, covariance }, resets, alignment))
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
fn selected(
    epoch: &ObservationEpoch,
    constellations: RtkConstellations,
) -> Result<BTreeMap<SatelliteId, &Observation>, RtkError> {
    let mut seen = std::collections::BTreeSet::new();
    let mut out = BTreeMap::<SatelliteId, &Observation>::new();
    for obs in &epoch.observations {
        if !seen.insert((obs.satellite, obs.tracking_code.as_str())) {
            return Err(RtkError::DuplicateObservation);
        }
        if constellations.includes(obs.satellite.system())
            && (obs.tracking_code == "1C"
                || obs.satellite.system() == GnssSystem::Galileo && obs.tracking_code == "1X")
        {
            let rank = |o: &Observation| {
                (
                    !o.pseudorange_m.is_some_and(|v| v.is_finite() && v > 0.0)
                        || !o.carrier_phase_cycles.is_some_and(f64::is_finite),
                    o.tracking_code.as_str() != "1C",
                )
            };
            if out
                .get(&obs.satellite)
                .is_none_or(|old| rank(obs) < rank(old))
            {
                out.insert(obs.satellite, obs);
            }
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
    let r = selected(rover, config.constellations)?;
    let b = selected(base, config.constellations)?;
    let mut data = Vec::new();
    let mut rejected = Vec::new();
    let ids: std::collections::BTreeSet<_> = rover
        .satellites()
        .union(&base.satellites())
        .copied()
        .collect();
    for satellite in ids {
        if !config.constellations.includes(satellite.system()) {
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
            rover_tracking_code: if ro.tracking_code == "1C" { "1C" } else { "1X" },
            base_tracking_code: if ba.tracking_code == "1C" { "1C" } else { "1X" },
        });
    }
    let mut counts = BTreeMap::<GnssSystem, usize>::new();
    for d in &data {
        *counts.entry(d.satellite.system()).or_default() += 1;
    }
    data.retain(|d| {
        if counts[&d.satellite.system()] >= 2 {
            true
        } else {
            rejected.push((d.satellite, RtkRejection::InsufficientSystemSatellites));
            false
        }
    });
    Ok((data, rejected))
}

fn dd_pairs(data: &[SatelliteData], references: &[usize]) -> Vec<(usize, usize)> {
    references
        .iter()
        .flat_map(|&reference| {
            data.iter()
                .enumerate()
                .filter(move |(i, d)| {
                    *i != reference && d.satellite.system() == data[reference].satellite.system()
                })
                .map(move |(i, _)| (reference, i))
        })
        .collect()
}

fn build_system_with_references(
    data: &[SatelliteData],
    references: &[usize],
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
    let pairs = dd_pairs(data, references);
    let count = pairs.len();
    let mut design = filter::zeros(2 * count, state.len());
    let mut residuals = vec![0.0; 2 * count];
    let mut noise = filter::zeros(2 * count, 2 * count);
    for (row, &(reference, sat)) in pairs.iter().enumerate() {
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
        let mut offset = 0;
        for &reference in references {
            let indices: Vec<_> = pairs
                .iter()
                .filter(|(r, _)| *r == reference)
                .map(|(_, sat)| *sat)
                .collect();
            let cov = double_difference_covariance(
                variance(elevations[reference]),
                &indices
                    .iter()
                    .map(|&i| variance(elevations[i]))
                    .collect::<Vec<_>>(),
            )?;
            for i in 0..indices.len() {
                for j in 0..indices.len() {
                    noise[block * count + offset + i][block * count + offset + j] = cov[i][j];
                }
            }
            offset += indices.len();
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
