//! Initial GPS ionosphere-free PPP FLOAT with precise products.
//! State is [ECEF XYZ, clock metres, zenith total delay metres, optional
//! ECEF velocity m/s, IF ambiguities]. Static remains the default.
//! Optional nominal-yaw phase wind-up; Bias-SINEX OSB is bound externally.
//! Optional fixed receiver calibration and nominal-yaw satellite IF-PCO.
//! Optional explicit legacy solid Earth tide. No satellite PCV or AR.
//! Equations follow pinned libgnss++ PPP; independent priors and transactional
//! IEKF updates are intentional initial Rust choices. See README and GOALS.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.

use crate::antenna::{
    AntennaError, ReceiverAntennaCorrection, ReceiverAntennaModel, SatelliteAntennaCorrection,
    SatelliteAntennaModel,
};
use crate::coordinates::{ecef_to_enu, ecef_to_geodetic, geometric_distance};
use crate::dual_frequency::{
    DualFrequencyConfig, IonosphereFreeObservation, PhaseBiasMode, SignalCorrections,
};
use crate::filter::{FilterError, GaussianState, Matrix, measurement_update, zeros};
use crate::models::{niell_hydrostatic_mapping, zenith_troposphere_climatology};
use crate::observation::ObservationEpoch;
use crate::phase_arc::{
    PhaseArcCommit, PhaseArcConfig, PhaseArcDecision, PhaseArcError, PhaseArcId, PhaseArcTracker,
};
use crate::precise::{OrbitReferencePoint, PreciseProducts};
use crate::precise_transmit::{
    PreciseTransmitConfig, PreciseTransmitError, PreciseTransmitState, evaluate_precise_transmit,
};
use crate::tides::{
    ReceiverPositionReference, SolidEarthTideCorrection, SolidEarthTideModel, TideError,
    solid_earth_tide,
};
use crate::time::day_of_year_gps_scale;
use crate::windup::{
    PppWindupModel, WindupCycles, WindupError, approximate_sun_position_ecef, nominal_yaw_windup,
};
use crate::{
    EcefCoord, Error, GnssSystem, GnssTime, PositionSolution, SatelliteId, SolutionStatus,
    constants as c,
};
use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PppKinematicConfig {
    pub initial_velocity_mps: [f64; 3],
    /// Independent initial variance for each ECEF velocity component.
    pub initial_velocity_variance_m2_s2: f64,
    /// Velocity random-walk spectral density, (m/s)^2 per second.
    pub velocity_process_noise_m2_s3: f64,
}
impl Default for PppKinematicConfig {
    fn default() -> Self {
        Self {
            initial_velocity_mps: [0.0; 3],
            initial_velocity_variance_m2_s2: 100.0,
            velocity_process_noise_m2_s3: 1e-4,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PppAccelerationConfig {
    pub initial_velocity_mps: [f64; 3],
    pub initial_velocity_variance_m2_s2: f64,
    pub initial_acceleration_mps2: [f64; 3],
    pub initial_acceleration_variance_m2_s4: f64,
    /// Local east/north acceleration random-walk sigma in m/s² per sqrt(s).
    pub horizontal_sigma_mps2_per_sqrt_s: f64,
    /// Local up acceleration random-walk sigma in m/s² per sqrt(s).
    pub vertical_sigma_mps2_per_sqrt_s: f64,
}
impl Default for PppAccelerationConfig {
    fn default() -> Self {
        Self {
            initial_velocity_mps: [0.0; 3],
            initial_velocity_variance_m2_s2: 100.0,
            initial_acceleration_mps2: [0.0; 3],
            initial_acceleration_variance_m2_s4: 100.0,
            horizontal_sigma_mps2_per_sqrt_s: 1.0,
            vertical_sigma_mps2_per_sqrt_s: 0.1,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum PppPositionModel {
    #[default]
    Static,
    /// Ordinary precise-PPP F(position,velocity)=dt and diagonal velocity Q.
    /// Does not fuse Doppler, SPP seeds, acceleration, inertial data or truth.
    ConstantVelocity(PppKinematicConfig),
    /// Native PVA F and local acceleration Q only; no SPP/Doppler fusion,
    /// reseeding or inertia input. XYZ/velocity process noise is unused.
    ConstantAcceleration(PppAccelerationConfig),
}
impl PppPositionModel {
    fn head_len(self) -> usize {
        match self {
            Self::Static => 5,
            Self::ConstantVelocity(_) => 8,
            Self::ConstantAcceleration(_) => 11,
        }
    }
    fn initial_velocity(self) -> Option<([f64; 3], f64)> {
        match self {
            Self::Static => None,
            Self::ConstantVelocity(c) => {
                Some((c.initial_velocity_mps, c.initial_velocity_variance_m2_s2))
            }
            Self::ConstantAcceleration(c) => {
                Some((c.initial_velocity_mps, c.initial_velocity_variance_m2_s2))
            }
        }
    }
    fn validate(self) -> bool {
        match self {
            Self::Static => true,
            Self::ConstantVelocity(c) => {
                c.initial_velocity_mps.iter().all(|v| v.is_finite())
                    && c.initial_velocity_variance_m2_s2.is_finite()
                    && c.initial_velocity_variance_m2_s2 > 0.0
                    && c.velocity_process_noise_m2_s3.is_finite()
                    && c.velocity_process_noise_m2_s3 >= 0.0
            }
            Self::ConstantAcceleration(c) => {
                c.initial_velocity_mps
                    .iter()
                    .chain(&c.initial_acceleration_mps2)
                    .all(|v| v.is_finite())
                    && [
                        c.initial_velocity_variance_m2_s2,
                        c.initial_acceleration_variance_m2_s4,
                    ]
                    .iter()
                    .all(|v| v.is_finite() && *v > 0.0)
                    && [
                        c.horizontal_sigma_mps2_per_sqrt_s,
                        c.vertical_sigma_mps2_per_sqrt_s,
                    ]
                    .iter()
                    .all(|v| v.is_finite() && *v >= 0.0 && (v * v).is_finite())
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PppFilterSettings {
    pub position_model: PppPositionModel,
    pub clock_model: PppClockModel,
    pub rejected_update_policy: PppRejectedUpdatePolicy,
    /// Optional caller-admitted SPP/Doppler motion initialization/recovery.
    pub motion_seed_policy: Option<PppMotionSeedPolicy>,
    /// Opt-in native covariance-standardized satellite exclusion (native sigma 4).
    pub standardized_residual_gate_sigma: Option<f64>,
    /// Strict calibrated phase remains the default. Explicit raw mode retains
    /// absent phase OSBs and absorbs constant biases only into FLOAT ambiguities.
    pub phase_bias_mode: PhaseBiasMode,
    /// Independent initial variance for XYZ, clock_m and zenith_m.
    pub initial_head_variance_m2: [f64; 5],
    /// Additive random-walk variance per second for the same five states.
    pub head_process_noise_m2_per_s: [f64; 5],
    /// Independent zero-mean prior for each new/reset IF ambiguity; no reuse
    /// of this epoch's measurements as an independent ambiguity prior.
    pub initial_ambiguity_variance_m2: f64,
    pub ambiguity_process_noise_m2_per_s: f64,
    pub elevation_mask_rad: f64,
    pub max_iterations: usize,
    pub convergence_m: f64,
    pub max_nis_per_row: f64,
    pub max_postfit_code_m: f64,
    pub max_postfit_phase_m: f64,
}
impl Default for PppFilterSettings {
    fn default() -> Self {
        Self {
            position_model: PppPositionModel::Static,
            clock_model: PppClockModel::RandomWalk,
            rejected_update_policy: PppRejectedUpdatePolicy::InterruptAllArcs,
            motion_seed_policy: None,
            standardized_residual_gate_sigma: None,
            phase_bias_mode: PhaseBiasMode::RequireSignalCorrections,
            initial_head_variance_m2: [1e6, 1e6, 1e6, 1e8, 100.0],
            head_process_noise_m2_per_s: [0.0, 0.0, 0.0, 100.0, 1e-8],
            initial_ambiguity_variance_m2: 1e14,
            ambiguity_process_noise_m2_per_s: 1e-8,
            elevation_mask_rad: 10_f64.to_radians(),
            max_iterations: 8,
            convergence_m: 1e-4,
            max_nis_per_row: 25.0,
            max_postfit_code_m: 2.0,
            max_postfit_phase_m: 0.05,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PppClockModel {
    #[default]
    RandomWalk,
    /// Fresh independent clock at each measurement update, with the configured
    /// initial clock variance. Centre this finite prior at the last committed
    /// clock mean; discard its correlations, preserving the other marginal.
    /// No SPP clock, Doppler or observed clock seed is fused. Clock Q is unused.
    IndependentEpoch,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PppRejectedUpdatePolicy {
    #[default]
    InterruptAllArcs,
    /// On measurement/geometry rejection, keep only previously accepted arcs
    /// with unchanged identity and no current slip/missing/gap/reset reason.
    /// Numeric state, accepted counts, GF/MW and wind-up anchors do not advance.
    /// Configuration, product, model and arithmetic errors still interrupt all.
    RetainContinuousArcs,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PppMotionSeedPolicy {
    pub recovery: crate::ppp_motion::PvaRecoveryConfig,
    pub initialize_from_seed: bool,
    pub require_seed: bool,
    pub update_velocity: bool,
    /// Total 3-D NIS; None matches the native ungated velocity update.
    pub velocity_nis_gate: Option<f64>,
}
impl Default for PppMotionSeedPolicy {
    fn default() -> Self {
        Self {
            recovery: Default::default(),
            initialize_from_seed: true,
            require_seed: true,
            update_velocity: true,
            velocity_nis_gate: None,
        }
    }
}
impl PppMotionSeedPolicy {
    fn valid(self) -> bool {
        [
            self.recovery.position_variance_m2,
            self.recovery.velocity_variance_m2_s2,
            self.recovery.acceleration_variance_m2_s4,
        ]
        .iter()
        .all(|x| x.is_finite() && *x > 0.0)
            && [
                self.recovery.position_disagreement_gate,
                self.velocity_nis_gate,
            ]
            .into_iter()
            .flatten()
            .all(|x| x.is_finite() && x > 0.0)
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct PppMotionSeed {
    pub solution: PositionSolution,
    /// Caller assertion of compatibility, not an automatic frame transform.
    pub coordinate_system: String,
    pub receiver_reference: ReceiverPositionReference,
    pub provenance: String,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PppMotionSeedError {
    Disabled,
    InvalidSolution,
    InvalidPosition,
    TimeMismatch,
    FrameMismatch,
    ReceiverReferenceMismatch,
    MissingProvenance,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PppSeedVelocityUpdate {
    Disabled,
    Unavailable,
    MissingCovariance,
    Applied { normalized_innovation_squared: f64 },
    Rejected(FilterError),
}
#[derive(Debug, Clone, PartialEq)]
pub struct PppMotionSeedDiagnostic {
    pub seed_time: GnssTime,
    pub provenance: String,
    pub initialized: bool,
    pub restarted: bool,
    pub variance_trigger: bool,
    pub normalized_position_disagreement: Option<f64>,
    pub seed_covariance_fallback: bool,
    pub velocity_update: PppSeedVelocityUpdate,
    /// All operations above act on a proposed prior until GNSS update success.
    pub committed: bool,
}
#[derive(Debug, Clone, PartialEq)]
pub struct PppResidualScreeningDiagnostic {
    pub satellite: SatelliteId,
    pub phase_row: bool,
    pub standardized_residual: f64,
    pub threshold_sigma: f64,
    pub remaining_code_satellites: usize,
    pub remaining_phase_satellites: usize,
    /// False means minimum code/phase counts prevent another update attempt.
    pub removed: bool,
}
#[derive(Debug, Clone, PartialEq)]
pub struct PppConfig {
    /// Independent initial position in the declared SP3 coordinate frame.
    pub initial_position: EcefCoord,
    pub initial_receiver_clock_s: f64,
    pub initial_zenith_delay_m: f64,
    pub expected_coordinate_system: String,
    pub expected_reference_point: OrbitReferencePoint,
    pub transmit: PreciseTransmitConfig,
    pub arcs: PhaseArcConfig,
    pub filter: PppFilterSettings,
    /// Explicit opt-in. Uncorrected synthetic data must keep this None.
    pub phase_windup: Option<PppWindupModel>,
    /// Fixed receiver calibration and marker-to-ARP geometry. No automatic
    /// RINEX antenna lookup; satellite APC/COM declaration is independent.
    pub receiver_antenna: Option<ReceiverAntennaModel>,
    /// Required for COM input, forbidden for APC input. Fixed GPS tracking
    /// pair/calibrations and explicit compatible clock-datum assertion.
    pub satellite_antenna: Option<SatelliteAntennaModel>,
    pub receiver_reference: ReceiverPositionReference,
    /// Explicit legacy alternative; not native default IERS Step-1+2.
    pub solid_earth_tide: Option<SolidEarthTideModel>,
}
#[derive(Debug, Clone, PartialEq)]
pub enum PppError {
    InvalidConfiguration,
    MotionSeed(PppMotionSeedError),
    MissingMotionSeed,
    ResidualScreeningExhausted {
        satellite: SatelliteId,
        code: usize,
        phase: usize,
    },
    ProductFrameMismatch,
    UnsupportedSelection(SatelliteId),
    InsufficientCode {
        available: usize,
    },
    InsufficientPhase {
        available: usize,
    },
    NonConvergent,
    PostfitRejected,
    NonPhysicalState,
    Gnss(Error),
    Arc(PhaseArcError),
    Filter(FilterError),
    Transmit(PreciseTransmitError),
    Windup(WindupError),
    Antenna(AntennaError),
    Tide(TideError),
}
impl fmt::Display for PppError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "initial PPP: {self:?}")
    }
}
impl std::error::Error for PppError {}
impl From<Error> for PppError {
    fn from(e: Error) -> Self {
        Self::Gnss(e)
    }
}
impl From<FilterError> for PppError {
    fn from(e: FilterError) -> Self {
        Self::Filter(e)
    }
}
impl From<PhaseArcError> for PppError {
    fn from(e: PhaseArcError) -> Self {
        Self::Arc(e)
    }
}
impl From<PreciseTransmitError> for PppError {
    fn from(e: PreciseTransmitError) -> Self {
        Self::Transmit(e)
    }
}
impl From<WindupError> for PppError {
    fn from(e: WindupError) -> Self {
        Self::Windup(e)
    }
}
impl From<AntennaError> for PppError {
    fn from(e: AntennaError) -> Self {
        Self::Antenna(e)
    }
}
impl From<TideError> for PppError {
    fn from(e: TideError) -> Self {
        Self::Tide(e)
    }
}
impl PppConfig {
    /// Opt in to the ordinary native precise-mode initial head prior.
    /// Position and clock seed means remain caller supplied; ZTD is seeded
    /// from native climatology. XYZ/clock/ZTD variances become 100/1e8/25 m².
    /// Constant-velocity mode starts at rest with 25 (m/s)², matching native
    /// non-PVA initialization. PVA starts V/A at zero with variances 100 each.
    /// Ambiguity priors, noise, gates, correction
    /// settings and subsequent dynamics remain independent and unchanged.
    /// This is not a whole-native PPP initialization or recovery policy.
    pub fn with_native_precise_head_prior(mut self, time: GnssTime) -> Result<Self, PppError> {
        self.initial_zenith_delay_m = zenith_troposphere_climatology(
            ecef_to_geodetic(self.initial_position)?,
            day_of_year_gps_scale(time)?,
        )?
        .total_delay_m();
        self.filter.initial_head_variance_m2 = [100.0, 100.0, 100.0, 1e8, 25.0];
        match &mut self.filter.position_model {
            PppPositionModel::ConstantVelocity(model) => {
                model.initial_velocity_mps = [0.0; 3];
                model.initial_velocity_variance_m2_s2 = 25.0;
            }
            PppPositionModel::ConstantAcceleration(model) => {
                model.initial_velocity_mps = [0.0; 3];
                model.initial_velocity_variance_m2_s2 = 100.0;
                model.initial_acceleration_mps2 = [0.0; 3];
                model.initial_acceleration_variance_m2_s4 = 100.0;
            }
            PppPositionModel::Static => (),
        }
        Ok(self)
    }

    fn validate(&self) -> Result<(), PppError> {
        self.check_references()?;
        if self.solid_earth_tide.is_some() {
            crate::tides::validate_station(self.initial_position)?;
        }
        let f = self.filter;
        ecef_to_geodetic(self.initial_position)?;
        self.transmit.validate()?;
        if f.standardized_residual_gate_sigma.is_some_and(|g| {
            !g.is_finite() || g <= 0.0 || matches!(f.position_model, PppPositionModel::Static)
        }) || !f.position_model.validate()
            || f.motion_seed_policy.is_some_and(|p| {
                !p.valid() || !matches!(f.position_model, PppPositionModel::ConstantAcceleration(_))
            })
            || !self.initial_receiver_clock_s.is_finite()
            || !(self.initial_receiver_clock_s * c::SPEED_OF_LIGHT).is_finite()
            || !self.initial_zenith_delay_m.is_finite()
            || self.initial_zenith_delay_m < 0.0
            || self.expected_coordinate_system.trim().is_empty()
            || f.initial_head_variance_m2
                .iter()
                .chain([&f.initial_ambiguity_variance_m2])
                .any(|v| !v.is_finite() || *v <= 0.0)
            || f.head_process_noise_m2_per_s
                .iter()
                .chain([&f.ambiguity_process_noise_m2_per_s])
                .any(|v| !v.is_finite() || *v < 0.0)
            || !f.elevation_mask_rad.is_finite()
            || !(0.0..std::f64::consts::FRAC_PI_2).contains(&f.elevation_mask_rad)
            || !(1..=32).contains(&f.max_iterations)
            || [
                f.convergence_m,
                f.max_nis_per_row,
                f.max_postfit_code_m,
                f.max_postfit_phase_m,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(PppError::InvalidConfiguration);
        }
        Ok(())
    }
    fn check_products(&self, products: &PreciseProducts) -> Result<(), PppError> {
        self.check_references()?;
        let header = products.orbit().header();
        if header.coordinate_system != self.expected_coordinate_system
            || header.reference_point != self.expected_reference_point
        {
            return Err(PppError::ProductFrameMismatch);
        }
        Ok(())
    }
    fn check_references(&self) -> Result<(), PppError> {
        if (self.receiver_reference == ReceiverPositionReference::NominalMarker)
            != self.solid_earth_tide.is_some()
        {
            return Err(PppError::InvalidConfiguration);
        }
        if (self.expected_reference_point == OrbitReferencePoint::CentreOfMass)
            != self.satellite_antenna.is_some()
        {
            return Err(PppError::InvalidConfiguration);
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct PppModel {
    /// Instantaneous marker used consistently for transmit, geometry and all
    /// receiver corrections. Equals the XYZ state when tide correction None.
    pub receiver_position_m: EcefCoord,
    pub solid_earth_tide: Option<SolidEarthTideCorrection>,
    /// Original precise transmit state, COM when satellite PCO is enabled.
    pub transmit: PreciseTransmitState,
    /// Geometry, antenna projection and wind-up use this APC position.
    pub satellite_position_m: EcefCoord,
    pub satellite_antenna: Option<SatelliteAntennaCorrection>,
    pub geometric_range_m: f64,
    pub elevation_rad: f64,
    pub troposphere_mapping: f64,
    pub code_prediction_m: f64,
    /// XYZ, clock metres, zenith total delay. Adds first-order Sagnac position
    /// partials to native -LOS; mapping/orbit/tide receiver partials are neglected.
    pub head_design: [f64; 5],
}
/// Native IF code model rho + receiver clock - c*satellite clock + NMF*ZTD.
/// Periodic relativity is already applied once in transmit.corrected_clock.
/// Carrier prediction is code_prediction + IF ambiguity metres.
pub fn evaluate_ppp_model(
    products: &PreciseProducts,
    satellite: SatelliteId,
    time: GnssTime,
    receiver: EcefCoord,
    receiver_clock_m: f64,
    zenith_delay_m: f64,
    config: &PppConfig,
) -> Result<PppModel, PppError> {
    config.check_products(products)?;
    if !receiver_clock_m.is_finite() || !zenith_delay_m.is_finite() {
        return Err(Error::NonFinite.into());
    }
    let tide = config
        .solid_earth_tide
        .map(|model| solid_earth_tide(time, receiver, model))
        .transpose()?;
    let receiver = tide.as_ref().map_or(receiver, |t| t.instantaneous_marker_m);
    let geo = ecef_to_geodetic(receiver)?;
    let tx = evaluate_precise_transmit(products, satellite, time, receiver, config.transmit)?;
    let satellite_antenna = if let Some(antenna) = &config.satellite_antenna {
        if tx.evaluated_state.clock_source != antenna.binding().clock_source {
            return Err(AntennaError::DatumMismatch.into());
        }
        Some(antenna.correction(
            time,
            satellite,
            tx.position_m,
            approximate_sun_position_ecef(time)?,
        )?)
    } else {
        None
    };
    let satellite_position = satellite_antenna
        .as_ref()
        .map_or(tx.position_m, |a| a.position_apc_m);
    let difference = EcefCoord::new(
        satellite_position.x - receiver.x,
        satellite_position.y - receiver.y,
        satellite_position.z - receiver.z,
    );
    let norm = difference.x.hypot(difference.y).hypot(difference.z);
    if !norm.is_finite() || norm <= 1.0 {
        return Err(Error::SingularGeometry.into());
    }
    let enu = ecef_to_enu(difference, geo)?;
    let elevation = enu.up.atan2(enu.east.hypot(enu.north));
    if elevation < 0.0 {
        return Err(PppError::NonPhysicalState);
    }
    let mapping = niell_hydrostatic_mapping(geo, elevation, day_of_year_gps_scale(time)?)?;
    let range = geometric_distance(satellite_position, receiver)?;
    let prediction = range + receiver_clock_m - c::SPEED_OF_LIGHT * tx.corrected_clock_bias_s
        + mapping * zenith_delay_m;
    if !prediction.is_finite() {
        return Err(Error::NonFinite.into());
    }
    let head_design = [
        -difference.x / norm - c::OMEGA_E * satellite_position.y / c::SPEED_OF_LIGHT,
        -difference.y / norm + c::OMEGA_E * satellite_position.x / c::SPEED_OF_LIGHT,
        -difference.z / norm,
        1.0,
        mapping,
    ];
    Ok(PppModel {
        receiver_position_m: receiver,
        solid_earth_tide: tide,
        transmit: tx,
        satellite_position_m: satellite_position,
        satellite_antenna,
        geometric_range_m: range,
        elevation_rad: elevation,
        troposphere_mapping: mapping,
        code_prediction_m: prediction,
        head_design,
    })
}

#[derive(Debug, Clone, PartialEq)]
pub enum PppRejection {
    StandardizedResidual { phase_row: bool, statistic: f64 },
    Model(PppError),
    ClockSourceMismatch,
    BelowElevationMask,
    NoMeasurement,
}
#[derive(Debug, Clone, PartialEq)]
pub struct PppMeasurement {
    pub satellite: SatelliteId,
    pub observation: IonosphereFreeObservation,
    pub model: PppModel,
    pub code_residual_m: f64,
    pub phase_residual_m: Option<f64>,
    pub windup: Option<PppWindupCorrection>,
    pub receiver_antenna: Option<ReceiverAntennaCorrection>,
}
#[derive(Debug, Clone, PartialEq)]
pub struct PppWindupCorrection {
    pub model: PppWindupModel,
    pub sun_position_m: EcefCoord,
    pub previous_cycles: f64,
    /// Fixed for this epoch's relinearizations. For a fresh arc, the first
    /// predicted-geometry value is unwrapped against zero, then held here.
    pub unwrap_reference_cycles: f64,
    pub angles: WindupCycles,
    /// Additive corrections: subtract cycles times each wavelength.
    pub phase_add_m: [f64; 2],
    pub phase_if_add_m: f64,
}
#[derive(Debug, Clone, PartialEq)]
pub struct PppWindupState {
    pub arc_id: PhaseArcId,
    pub time: GnssTime,
    pub cycles: f64,
}
#[derive(Debug, Clone, PartialEq)]
pub struct PppAmbiguity {
    pub satellite: SatelliteId,
    pub arc_id: PhaseArcId,
    pub metres: f64,
    pub variance_m2: f64,
}
#[derive(Debug, Clone, PartialEq)]
pub struct PppSolution {
    pub position_model: PppPositionModel,
    pub clock_model: PppClockModel,
    pub rejected_update_policy: PppRejectedUpdatePolicy,
    pub phase_bias_mode: PhaseBiasMode,
    /// XYZ state/report reference; model receiver_position_m is instantaneous.
    pub receiver_reference: ReceiverPositionReference,
    pub position: PositionSolution,
    pub zenith_delay_m: f64,
    pub zenith_variance_m2: f64,
    pub acceleration_ecef_mps2: Option<[f64; 3]>,
    pub acceleration_covariance_m2_s4: Option<[[f64; 3]; 3]>,
    pub motion_seed: Option<PppMotionSeedDiagnostic>,
    pub residual_screening: Vec<PppResidualScreeningDiagnostic>,
    pub ambiguities: Vec<PppAmbiguity>,
    pub measurements: Vec<PppMeasurement>,
    pub arc_decisions: BTreeMap<SatelliteId, PhaseArcDecision>,
    pub arc_commit: PhaseArcCommit,
    pub rejected: BTreeMap<SatelliteId, PppRejection>,
    pub iterations: usize,
    pub normalized_innovation_squared: f64,
    pub code_rms_m: f64,
    pub phase_rms_m: f64,
    pub coordinate_system: String,
    /// Input SP3 reference. With satellite PCO enabled the model geometry is
    /// APC, while the original transmit state and this declaration remain COM.
    pub reference_point: OrbitReferencePoint,
}
#[derive(Debug, Clone, PartialEq)]
pub enum PppEpoch {
    Solution(Box<PppSolution>),
    Event(PhaseArcCommit),
}

#[derive(Debug)]
pub struct PppFilter {
    config: PppConfig,
    state: GaussianState,
    state_time: Option<GnssTime>,
    /// Last successfully committed numeric state's ambiguity layout.
    layout: BTreeMap<SatelliteId, (usize, PhaseArcId)>,
    arcs: PhaseArcTracker,
    last_rejections: BTreeMap<SatelliteId, PppRejection>,
    last_retained_phase_arcs: BTreeSet<SatelliteId>,
    last_motion_seed: Option<PppMotionSeedDiagnostic>,
    last_residual_screening: Vec<PppResidualScreeningDiagnostic>,
    windup: BTreeMap<SatelliteId, PppWindupState>,
}
impl PppFilter {
    pub fn new(config: PppConfig) -> Result<Self, PppError> {
        config.validate()?;
        let arcs = PhaseArcTracker::new(config.arcs)?;
        let head_len = config.filter.position_model.head_len();
        let mut covariance = zeros(head_len, head_len);
        for (i, row) in covariance.iter_mut().take(5).enumerate() {
            row[i] = config.filter.initial_head_variance_m2[i];
        }
        let mut mean = vec![
            config.initial_position.x,
            config.initial_position.y,
            config.initial_position.z,
            config.initial_receiver_clock_s * c::SPEED_OF_LIGHT,
            config.initial_zenith_delay_m,
        ];
        if let Some((velocity, variance)) = config.filter.position_model.initial_velocity() {
            mean.extend(velocity);
            for (i, row) in covariance.iter_mut().enumerate().skip(5).take(3) {
                row[i] = variance;
            }
        }
        if let PppPositionModel::ConstantAcceleration(c) = config.filter.position_model {
            mean.extend(c.initial_acceleration_mps2);
            for (i, row) in covariance.iter_mut().enumerate().skip(8).take(3) {
                row[i] = c.initial_acceleration_variance_m2_s4;
            }
        }
        Ok(Self {
            config,
            state: GaussianState { mean, covariance },
            state_time: None,
            layout: BTreeMap::new(),
            arcs,
            last_rejections: BTreeMap::new(),
            last_retained_phase_arcs: BTreeSet::new(),
            last_motion_seed: None,
            last_residual_screening: Vec::new(),
            windup: BTreeMap::new(),
        })
    }
    /// Last committed numeric state is preserved on rejected updates/events.
    pub fn state(&self) -> &GaussianState {
        &self.state
    }
    pub fn state_time(&self) -> Option<GnssTime> {
        self.state_time
    }
    pub fn phase_arc(&self, sat: SatelliteId) -> Option<&crate::phase_arc::PhaseArc> {
        self.arcs.arc(sat)
    }
    pub fn last_rejections(&self) -> &BTreeMap<SatelliteId, PppRejection> {
        &self.last_rejections
    }
    /// Arcs retained after the last processed rejected measurement update.
    /// Empty after a successful update, event or global invalidation.
    pub fn last_retained_phase_arcs(&self) -> &BTreeSet<SatelliteId> {
        &self.last_retained_phase_arcs
    }
    pub fn last_residual_screening(&self) -> &[PppResidualScreeningDiagnostic] {
        &self.last_residual_screening
    }
    pub fn last_motion_seed(&self) -> Option<&PppMotionSeedDiagnostic> {
        self.last_motion_seed.as_ref()
    }
    pub fn windup_state(&self, sat: SatelliteId) -> Option<&PppWindupState> {
        self.windup.get(&sat)
    }
    pub fn invalidate_phase(&mut self) -> Result<(), PppError> {
        self.arcs.invalidate_all()?;
        self.windup.clear();
        self.last_retained_phase_arcs.clear();
        self.last_motion_seed = None;
        self.last_residual_screening.clear();
        Ok(())
    }

    /// GPS-only IF PPP in the configured position model. Every code/phase requires declared
    /// corrections, matching actual selected clock source, frame and reference
    /// point. At least five independent code rows and four phases are required.
    /// Failure preserves numeric mean/covariance. The default interrupts all
    /// phase continuity; optional retention never accepts rejected measurements.
    pub fn process_epoch(
        &mut self,
        epoch: &ObservationEpoch,
        products: &PreciseProducts,
        corrections: &SignalCorrections,
        selections: &BTreeMap<SatelliteId, DualFrequencyConfig>,
    ) -> Result<PppEpoch, PppError> {
        self.process_epoch_with_seed(epoch, products, corrections, selections, None)
    }

    /// Explicit SPP/Doppler motion seed. A same-epoch restart is a heuristic
    /// data-dependent prior, not an independent code-measurement factor.
    /// Seed receiver clock and clock drift are never imported. The legitimate
    /// full-state correlations of a velocity measurement remain active.
    pub fn process_epoch_with_seed(
        &mut self,
        epoch: &ObservationEpoch,
        products: &PreciseProducts,
        corrections: &SignalCorrections,
        selections: &BTreeMap<SatelliteId, DualFrequencyConfig>,
        seed: Option<&PppMotionSeed>,
    ) -> Result<PppEpoch, PppError> {
        if self.arcs.last_epoch().is_some_and(|t| epoch.time <= t) {
            return Err(PhaseArcError::OutOfOrder.into());
        }
        if let Err(error) = self.config.check_products(products) {
            self.invalidate_phase()?;
            return Err(error);
        }
        if epoch.flag != 6
            && let Some(seed) = seed
            && let Err(reason) = self.check_motion_seed(epoch.time, seed)
        {
            self.invalidate_phase()?;
            return Err(PppError::MotionSeed(reason));
        }
        for sat in selections.keys() {
            if sat.system() != GnssSystem::Gps || sat.prn() > 32 {
                self.invalidate_phase()?;
                return Err(PppError::UnsupportedSelection(*sat));
            }
        }
        if let Some(antenna) = &self.config.satellite_antenna {
            for selection in selections.values() {
                if let Err(error) = antenna.check_datum(
                    &selection.tracking_codes,
                    selection.clock_source,
                    &selection.clock_reference,
                ) {
                    self.invalidate_phase()?;
                    return Err(error.into());
                }
            }
        }
        let prepared = self.arcs.prepare_epoch_with_phase_bias_mode(
            epoch,
            corrections,
            selections,
            self.config.filter.phase_bias_mode,
        );
        let decisions = match prepared {
            Ok(prepared) => prepared.decisions().clone(),
            Err(error) => {
                self.windup.clear();
                self.last_retained_phase_arcs.clear();
                self.last_motion_seed = None;
                self.last_residual_screening.clear();
                return Err(error.into());
            }
        };
        self.last_rejections.clear();
        self.last_retained_phase_arcs.clear();
        self.last_motion_seed = None;
        self.last_residual_screening.clear();
        if epoch.flag == 6 {
            self.windup.clear();
            return Ok(PppEpoch::Event(self.arcs.finish_epoch(&BTreeSet::new())?));
        }
        let result = if seed.is_none()
            && self
                .config
                .filter
                .motion_seed_policy
                .is_some_and(|p| p.require_seed)
        {
            Err(PppError::MissingMotionSeed)
        } else {
            self.update(epoch, products, &decisions, seed)
        };
        match result {
            Ok((posterior, layout, measurements, iterations, nis)) => {
                let eligible = layout.keys().copied().collect();
                let commit = self.arcs.finish_epoch(&eligible)?;
                if let Some(d) = &mut self.last_motion_seed {
                    d.committed = true;
                }
                self.windup = measurements
                    .iter()
                    .filter_map(|m| {
                        m.windup.as_ref().map(|w| {
                            (
                                m.satellite,
                                PppWindupState {
                                    arc_id: layout[&m.satellite].1,
                                    time: epoch.time,
                                    cycles: w.angles.cycles,
                                },
                            )
                        })
                    })
                    .collect();
                let position = PositionSolution {
                    time: epoch.time,
                    status: SolutionStatus::PppFloat,
                    position_ecef: Some(receiver(&posterior.mean)),
                    position_covariance: Some(std::array::from_fn(|i| {
                        std::array::from_fn(|j| posterior.covariance[i][j])
                    })),
                    num_satellites: measurements.len(),
                    receiver_clock_bias: posterior.mean[3] / c::SPEED_OF_LIGHT,
                    velocity_ecef: self
                        .config
                        .filter
                        .position_model
                        .initial_velocity()
                        .is_some()
                        .then(|| std::array::from_fn(|i| posterior.mean[5 + i])),
                    velocity_covariance: self
                        .config
                        .filter
                        .position_model
                        .initial_velocity()
                        .is_some()
                        .then(|| {
                            std::array::from_fn(|i| {
                                std::array::from_fn(|j| posterior.covariance[5 + i][5 + j])
                            })
                        }),
                };
                let ambiguities = layout
                    .iter()
                    .map(|(&sat, &(i, id))| PppAmbiguity {
                        satellite: sat,
                        arc_id: id,
                        metres: posterior.mean[i],
                        variance_m2: posterior.covariance[i][i],
                    })
                    .collect();
                let code_rms_m = (measurements
                    .iter()
                    .map(|m| m.code_residual_m * m.code_residual_m)
                    .sum::<f64>()
                    / measurements.len() as f64)
                    .sqrt();
                let phase_rms_m = (measurements
                    .iter()
                    .filter_map(|m| m.phase_residual_m)
                    .map(|v| v * v)
                    .sum::<f64>()
                    / layout.len() as f64)
                    .sqrt();
                let solution = PppSolution {
                    position_model: self.config.filter.position_model,
                    clock_model: self.config.filter.clock_model,
                    rejected_update_policy: self.config.filter.rejected_update_policy,
                    phase_bias_mode: self.config.filter.phase_bias_mode,
                    receiver_reference: self.config.receiver_reference,
                    position,
                    zenith_delay_m: posterior.mean[4],
                    zenith_variance_m2: posterior.covariance[4][4],
                    acceleration_ecef_mps2: matches!(
                        self.config.filter.position_model,
                        PppPositionModel::ConstantAcceleration(_)
                    )
                    .then(|| std::array::from_fn(|i| posterior.mean[8 + i])),
                    acceleration_covariance_m2_s4: matches!(
                        self.config.filter.position_model,
                        PppPositionModel::ConstantAcceleration(_)
                    )
                    .then(|| {
                        std::array::from_fn(|i| {
                            std::array::from_fn(|j| posterior.covariance[8 + i][8 + j])
                        })
                    }),
                    motion_seed: self.last_motion_seed.clone(),
                    residual_screening: self.last_residual_screening.clone(),
                    ambiguities,
                    measurements,
                    arc_decisions: decisions,
                    arc_commit: commit,
                    rejected: self.last_rejections.clone(),
                    iterations,
                    normalized_innovation_squared: nis,
                    code_rms_m,
                    phase_rms_m,
                    coordinate_system: self.config.expected_coordinate_system.clone(),
                    reference_point: self.config.expected_reference_point,
                };
                self.state = posterior;
                self.layout = layout;
                self.state_time = Some(epoch.time);
                Ok(PppEpoch::Solution(Box::new(solution)))
            }
            Err(error) => {
                let retain = self.config.filter.rejected_update_policy
                    == PppRejectedUpdatePolicy::RetainContinuousArcs
                    && matches!(
                        error,
                        PppError::InsufficientCode { .. }
                            | PppError::InsufficientPhase { .. }
                            | PppError::ResidualScreeningExhausted { .. }
                            | PppError::NonConvergent
                            | PppError::PostfitRejected
                            | PppError::NonPhysicalState
                            | PppError::Filter(FilterError::InnovationRejected)
                            | PppError::Gnss(Error::SingularGeometry)
                    );
                if retain {
                    let screened = self
                        .last_residual_screening
                        .iter()
                        .map(|d| d.satellite)
                        .collect();
                    let rejection = self
                        .arcs
                        .reject_epoch_retaining_continuous_excluding(&screened)?;
                    self.last_retained_phase_arcs = rejection.retained.keys().copied().collect();
                    self.windup
                        .retain(|sat, _| self.last_retained_phase_arcs.contains(sat));
                } else {
                    self.arcs.finish_epoch(&BTreeSet::new())?;
                    self.windup.clear();
                }
                Err(error)
            }
        }
    }

    fn check_motion_seed(
        &self,
        time: GnssTime,
        seed: &PppMotionSeed,
    ) -> Result<(), PppMotionSeedError> {
        if self.config.filter.motion_seed_policy.is_none() {
            return Err(PppMotionSeedError::Disabled);
        }
        if seed.solution.status != SolutionStatus::Spp || !seed.solution.is_valid() {
            return Err(PppMotionSeedError::InvalidSolution);
        }
        if ecef_to_geodetic(seed.solution.position_ecef.unwrap()).is_err() {
            return Err(PppMotionSeedError::InvalidPosition);
        }
        if seed.solution.time.difference_seconds(time).abs() > 1e-7 {
            return Err(PppMotionSeedError::TimeMismatch);
        }
        if seed.coordinate_system != self.config.expected_coordinate_system {
            return Err(PppMotionSeedError::FrameMismatch);
        }
        if seed.receiver_reference != self.config.receiver_reference {
            return Err(PppMotionSeedError::ReceiverReferenceMismatch);
        }
        if seed.provenance.trim().is_empty() {
            return Err(PppMotionSeedError::MissingProvenance);
        }
        Ok(())
    }

    fn apply_motion_seed(
        &mut self,
        prior: &mut GaussianState,
        seed: Option<&PppMotionSeed>,
    ) -> Result<(), PppError> {
        let Some(policy) = self.config.filter.motion_seed_policy else {
            return Ok(());
        };
        let Some(seed) = seed else {
            return Ok(());
        };
        let position = seed.solution.position_ecef.unwrap();
        let mut diagnostic = PppMotionSeedDiagnostic {
            seed_time: seed.solution.time,
            provenance: seed.provenance.clone(),
            initialized: false,
            restarted: false,
            variance_trigger: false,
            normalized_position_disagreement: None,
            seed_covariance_fallback: false,
            velocity_update: if policy.update_velocity {
                PppSeedVelocityUpdate::Unavailable
            } else {
                PppSeedVelocityUpdate::Disabled
            },
            committed: false,
        };
        if self.state_time.is_none() && policy.initialize_from_seed {
            let PppPositionModel::ConstantAcceleration(model) = self.config.filter.position_model
            else {
                return Err(PppError::InvalidConfiguration);
            };
            let velocity = seed
                .solution
                .velocity_ecef
                .filter(|v| v.iter().all(|x| x.is_finite()))
                .unwrap_or([0.0; 3]);
            for (start, mean, variance) in [
                (
                    0,
                    [position.x, position.y, position.z],
                    self.config.filter.initial_head_variance_m2[..3]
                        .try_into()
                        .unwrap(),
                ),
                (5, velocity, [model.initial_velocity_variance_m2_s2; 3]),
                (8, [0.0; 3], [model.initial_acceleration_variance_m2_s4; 3]),
            ] {
                for j in 0..3 {
                    let i = start + j;
                    prior.mean[i] = mean[j];
                    for row in &mut prior.covariance {
                        row[i] = 0.0;
                    }
                    prior.covariance[i].fill(0.0);
                    prior.covariance[i][i] = variance[j];
                }
            }
            diagnostic.initialized = true;
        } else {
            let r = crate::ppp_motion::recover_pva_motion(
                prior,
                Some(crate::ppp_motion::MotionSeed {
                    position_ecef: position,
                    position_covariance_m2: seed.solution.position_covariance,
                    velocity_mps: seed.solution.velocity_ecef,
                }),
                policy.recovery,
            )?;
            diagnostic.restarted = r.reset;
            diagnostic.variance_trigger = r.variance_trigger;
            diagnostic.normalized_position_disagreement = r.normalized_position_disagreement;
            diagnostic.seed_covariance_fallback = r.seed_covariance_fallback;
            *prior = r.prior;
        }
        if policy.update_velocity
            && let Some(velocity) = seed.solution.velocity_ecef
        {
            if let Some(covariance) = seed.solution.velocity_covariance {
                diagnostic.velocity_update = match crate::ppp_motion::update_velocity(
                    prior,
                    5,
                    velocity,
                    covariance,
                    policy.velocity_nis_gate,
                ) {
                    Ok(update) => {
                        *prior = update.posterior;
                        PppSeedVelocityUpdate::Applied {
                            normalized_innovation_squared: update.normalized_innovation_squared,
                        }
                    }
                    Err(error) => PppSeedVelocityUpdate::Rejected(error),
                };
            } else {
                diagnostic.velocity_update = PppSeedVelocityUpdate::MissingCovariance;
            }
        }
        self.last_motion_seed = Some(diagnostic);
        Ok(())
    }

    fn predicted_prior(
        &self,
        time: GnssTime,
        layout: &BTreeMap<SatelliteId, (usize, PhaseArcId)>,
    ) -> Result<GaussianState, PppError> {
        let head_len = self.config.filter.position_model.head_len();
        let n = head_len + layout.len();
        let mut mean = vec![0.0; n];
        let mut covariance = zeros(n, n);
        let mut source: Vec<Option<usize>> = (0..head_len).map(Some).collect();
        for (&sat, &(_, id)) in layout {
            source.push(
                self.layout
                    .get(&sat)
                    .filter(|(_, old_id)| *old_id == id)
                    .map(|(i, _)| *i),
            );
        }
        let dt = self.state_time.map_or(0.0, |t| time.difference_seconds(t));
        for i in 0..n {
            if let Some(old) = source[i] {
                mean[i] = self.state.mean[old];
            }
            for j in 0..n {
                if let (Some(a), Some(b)) = (source[i], source[j]) {
                    covariance[i][j] = self.state.covariance[a][b];
                }
            }
            if i >= head_len && source[i].is_none() {
                covariance[i][i] = self.config.filter.initial_ambiguity_variance_m2;
            }
            if matches!(self.config.filter.position_model, PppPositionModel::Static)
                && !(i == 3 && self.config.filter.clock_model == PppClockModel::IndependentEpoch)
            {
                covariance[i][i] += dt
                    * if i < 5 {
                        self.config.filter.head_process_noise_m2_per_s[i]
                    } else {
                        self.config.filter.ambiguity_process_noise_m2_per_s
                    };
            }
        }
        if let PppPositionModel::ConstantVelocity(c) = self.config.filter.position_model {
            let mut head_noise = self.config.filter.head_process_noise_m2_per_s;
            if self.config.filter.clock_model == PppClockModel::IndependentEpoch {
                head_noise[3] = 0.0;
            }
            predict_kinematic_state(
                &mut mean,
                &mut covariance,
                dt,
                c,
                head_noise,
                self.config.filter.ambiguity_process_noise_m2_per_s,
            );
        }
        if let PppPositionModel::ConstantAcceleration(c) = self.config.filter.position_model {
            let mut head_noise = self.config.filter.head_process_noise_m2_per_s;
            if self.config.filter.clock_model == PppClockModel::IndependentEpoch {
                head_noise[3] = 0.0;
            }
            predict_acceleration_state(
                &mut mean,
                &mut covariance,
                dt,
                c,
                head_noise,
                self.config.filter.ambiguity_process_noise_m2_per_s,
            )?;
        }
        if self.config.filter.clock_model == PppClockModel::IndependentEpoch {
            reset_epoch_clock_covariance(
                &mut covariance,
                self.config.filter.initial_head_variance_m2[3],
            );
        }
        if !mean
            .iter()
            .chain(covariance.iter().flatten())
            .all(|v| v.is_finite())
        {
            return Err(Error::NonFinite.into());
        }
        Ok(GaussianState { mean, covariance })
    }
    fn update(
        &mut self,
        epoch: &ObservationEpoch,
        products: &PreciseProducts,
        decisions: &BTreeMap<SatelliteId, PhaseArcDecision>,
        seed: Option<&PppMotionSeed>,
    ) -> Result<UpdateResult, PppError> {
        let mut admission_receiver = receiver(&self.state.mean);
        if self
            .config
            .filter
            .position_model
            .initial_velocity()
            .is_some()
        {
            let dt = self
                .state_time
                .map_or(0.0, |t| epoch.time.difference_seconds(t));
            admission_receiver = EcefCoord::new(
                admission_receiver.x + dt * self.state.mean[5],
                admission_receiver.y + dt * self.state.mean[6],
                admission_receiver.z + dt * self.state.mean[7],
            );
            if matches!(
                self.config.filter.position_model,
                PppPositionModel::ConstantAcceleration(_)
            ) {
                admission_receiver.x += 0.5 * dt * dt * self.state.mean[8];
                admission_receiver.y += 0.5 * dt * dt * self.state.mean[9];
                admission_receiver.z += 0.5 * dt * dt * self.state.mean[10];
            }
            if !admission_receiver.is_finite() {
                return Err(Error::NonFinite.into());
            }
        }
        if self.config.filter.motion_seed_policy.is_some() && seed.is_some() {
            let mut preview = self.predicted_prior(epoch.time, &BTreeMap::new())?;
            self.apply_motion_seed(&mut preview, seed)?;
            admission_receiver = receiver(&preview.mean);
        }
        let mut observations = BTreeMap::new();
        let mut layout = BTreeMap::new();
        for (&sat, decision) in decisions {
            let Some(observation) = &decision.measurement else {
                self.last_rejections
                    .insert(sat, PppRejection::NoMeasurement);
                continue;
            };
            let model = match evaluate_ppp_model(
                products,
                sat,
                epoch.time,
                admission_receiver,
                self.state.mean[3],
                self.state.mean[4],
                &self.config,
            ) {
                Ok(m) => m,
                Err(error @ (PppError::Antenna(_) | PppError::Tide(_))) => return Err(error),
                Err(error) => {
                    self.last_rejections.insert(sat, PppRejection::Model(error));
                    continue;
                }
            };
            if observation.corrections[0].clock_source
                != model.transmit.evaluated_state.clock_source
            {
                self.last_rejections
                    .insert(sat, PppRejection::ClockSourceMismatch);
                continue;
            }
            if model.elevation_rad < self.config.filter.elevation_mask_rad {
                self.last_rejections
                    .insert(sat, PppRejection::BelowElevationMask);
                continue;
            }
            if let Some(id) = decision.proposed_arc {
                layout.insert(
                    sat,
                    (
                        self.config.filter.position_model.head_len() + layout.len(),
                        id,
                    ),
                );
            }
            observations.insert(sat, observation.clone());
        }
        'screening: loop {
            if observations.len() < 5 {
                return Err(PppError::InsufficientCode {
                    available: observations.len(),
                });
            }
            if layout.len() < 4 {
                return Err(PppError::InsufficientPhase {
                    available: layout.len(),
                });
            }
            let mut prior = self.predicted_prior(epoch.time, &layout)?;
            self.apply_motion_seed(&mut prior, seed)?;
            let windup_context = self.prepare_windup(epoch.time, products, &layout, &prior.mean)?;
            let mut iterate = prior.mean.clone();
            for iteration in 1..=self.config.filter.max_iterations {
                let system = self.measurement_system(
                    epoch.time,
                    products,
                    &observations,
                    &layout,
                    &iterate,
                    windup_context.as_ref(),
                )?;
                ensure_geometry(&system.code_design)?;
                // IEKF always relinearizes the SAME predicted prior. This correction
                // prevents repeated iterations treating one epoch as independent data.
                let innovation: Vec<_> = system
                    .innovations
                    .iter()
                    .zip(&system.design)
                    .map(|(r, h)| {
                        r + h
                            .iter()
                            .zip(iterate.iter().zip(&prior.mean))
                            .map(|(h, (x, p))| h * (x - p))
                            .sum::<f64>()
                    })
                    .collect();
                if iteration == 1
                    && let Some(gate) = self.config.filter.standardized_residual_gate_sigma
                {
                    let screen = crate::ppp_screening::screen_measurements(
                        &prior,
                        &system.design,
                        &innovation,
                        &system.noise,
                        gate,
                    )?;
                    if let Some(row) = screen.worst_row {
                        let (satellite, phase_row) = system.row_satellites[row];
                        let code = observations.len() - 1;
                        let phase = layout.len() - usize::from(layout.contains_key(&satellite));
                        let removed = code >= 5 && phase >= 4;
                        let statistic = screen.standardized_residuals[row];
                        self.last_residual_screening
                            .push(PppResidualScreeningDiagnostic {
                                satellite,
                                phase_row,
                                standardized_residual: statistic,
                                threshold_sigma: gate,
                                remaining_code_satellites: code,
                                remaining_phase_satellites: phase,
                                removed,
                            });
                        self.last_rejections.insert(
                            satellite,
                            PppRejection::StandardizedResidual {
                                phase_row,
                                statistic,
                            },
                        );
                        if !removed {
                            return Err(PppError::ResidualScreeningExhausted {
                                satellite,
                                code,
                                phase,
                            });
                        }
                        observations.remove(&satellite);
                        layout = layout
                            .iter()
                            .filter(|(sat, _)| **sat != satellite)
                            .enumerate()
                            .map(|(j, (&sat, &(_, id)))| {
                                (sat, (self.config.filter.position_model.head_len() + j, id))
                            })
                            .collect();
                        // Every pass starts from the same committed state. Never feed
                        // a rejected GNSS update or velocity-fused proposal back in.
                        continue 'screening;
                    }
                }
                let update = measurement_update(
                    &prior,
                    &system.design,
                    &innovation,
                    &system.noise,
                    Some(self.config.filter.max_nis_per_row),
                )?;
                let difference = update
                    .posterior
                    .mean
                    .iter()
                    .zip(&iterate)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max);
                iterate = update.posterior.mean.clone();
                if difference <= self.config.filter.convergence_m {
                    if iterate[4] < 0.0 || !iterate[4].is_finite() {
                        return Err(PppError::NonPhysicalState);
                    }
                    let final_system = self.measurement_system(
                        epoch.time,
                        products,
                        &observations,
                        &layout,
                        &iterate,
                        windup_context.as_ref(),
                    )?;
                    if final_system.measurements.iter().any(|m| {
                        m.code_residual_m.abs() > self.config.filter.max_postfit_code_m
                            || m.phase_residual_m
                                .is_some_and(|r| r.abs() > self.config.filter.max_postfit_phase_m)
                    }) {
                        return Err(PppError::PostfitRejected);
                    }
                    return Ok((
                        update.posterior,
                        layout,
                        final_system.measurements,
                        iteration,
                        update.normalized_innovation_squared,
                    ));
                }
            }
            return Err(PppError::NonConvergent);
        }
    }
    fn prepare_windup(
        &self,
        time: GnssTime,
        products: &PreciseProducts,
        layout: &Layout,
        mean: &[f64],
    ) -> Result<Option<WindupEpochContext>, PppError> {
        let Some(model) = self.config.phase_windup else {
            return Ok(None);
        };
        let sun = approximate_sun_position_ecef(time)?;
        let mut anchors = BTreeMap::new();
        for (&sat, &(_, id)) in layout {
            let (previous, reference) = if let Some(old) =
                self.windup.get(&sat).filter(|s| s.arc_id == id)
            {
                (old.cycles, old.cycles)
            } else {
                let m = evaluate_ppp_model(
                    products,
                    sat,
                    time,
                    receiver(mean),
                    mean[3],
                    mean[4],
                    &self.config,
                )?;
                (
                    0.0,
                    nominal_yaw_windup(m.receiver_position_m, m.satellite_position_m, sun, 0.0)?
                        .cycles,
                )
            };
            anchors.insert(sat, (previous, reference));
        }
        Ok(Some(WindupEpochContext {
            model,
            sun,
            anchors,
        }))
    }
    fn measurement_system(
        &self,
        time: GnssTime,
        products: &PreciseProducts,
        observations: &BTreeMap<SatelliteId, IonosphereFreeObservation>,
        layout: &BTreeMap<SatelliteId, (usize, PhaseArcId)>,
        mean: &[f64],
        windup_context: Option<&WindupEpochContext>,
    ) -> Result<MeasurementSystem, PppError> {
        let mut design = Vec::new();
        let mut innovations = Vec::new();
        let mut variances = Vec::new();
        let mut measurements = Vec::new();
        let mut code_design = Vec::new();
        let mut row_satellites = Vec::new();
        for (&sat, observation) in observations {
            let model = evaluate_ppp_model(
                products,
                sat,
                time,
                receiver(mean),
                mean[3],
                mean[4],
                &self.config,
            )?;
            if model.elevation_rad < self.config.filter.elevation_mask_rad {
                return Err(PppError::NonPhysicalState);
            }
            if model.transmit.evaluated_state.clock_source
                != observation.corrections[0].clock_source
            {
                return Err(PppError::InvalidConfiguration);
            }
            let mut row = vec![0.0; mean.len()];
            row[..5].copy_from_slice(&model.head_design);
            code_design.push(model.head_design);
            let antenna = self
                .config
                .receiver_antenna
                .as_ref()
                .map(|a| {
                    a.correction(
                        time,
                        sat,
                        &observation.raw.clone().map(|r| r.tracking_code),
                        observation.coefficients,
                        [model.receiver_position_m, model.satellite_position_m],
                    )
                })
                .transpose()?;
            let antenna_add = antenna.as_ref().map_or(0.0, |a| a.if_add_m);
            let code_residual = observation.code_if_m + antenna_add - model.code_prediction_m;
            design.push(row.clone());
            innovations.push(code_residual);
            row_satellites.push((sat, false));
            variances.push(observation.code_variance_m2);
            let mut windup = None;
            let phase_residual = if let Some(&(index, _)) = layout.get(&sat) {
                let phase = observation.phase.as_ref().unwrap();
                if let Some(context) = windup_context {
                    let (previous, reference) = context.anchors[&sat];
                    let sun = context.sun;
                    let angles = nominal_yaw_windup(
                        model.receiver_position_m,
                        model.satellite_position_m,
                        sun,
                        reference,
                    )?;
                    let phase_add_m = observation.wavelength_m.map(|w| -angles.cycles * w);
                    // a*lambda1 + b*lambda2 == c/(f1+f2); native IF contract.
                    let phase_if_add_m = -angles.cycles * c::SPEED_OF_LIGHT
                        / (observation.frequency_hz[0] + observation.frequency_hz[1]);
                    windup = Some(PppWindupCorrection {
                        model: context.model,
                        sun_position_m: sun,
                        previous_cycles: previous,
                        unwrap_reference_cycles: reference,
                        angles,
                        phase_add_m,
                        phase_if_add_m,
                    });
                }
                let residual = phase.phase_if_m
                    + antenna_add
                    + windup.as_ref().map_or(0.0, |w| w.phase_if_add_m)
                    - (model.code_prediction_m + mean[index]);
                row[index] = 1.0;
                design.push(row);
                innovations.push(residual);
                row_satellites.push((sat, true));
                variances.push(phase.variance_m2);
                Some(residual)
            } else {
                None
            };
            measurements.push(PppMeasurement {
                satellite: sat,
                observation: observation.clone(),
                model,
                code_residual_m: code_residual,
                phase_residual_m: phase_residual,
                windup,
                receiver_antenna: antenna,
            });
        }
        let mut noise = zeros(variances.len(), variances.len());
        for (i, v) in variances.iter().enumerate() {
            noise[i][i] = *v;
        }
        Ok(MeasurementSystem {
            design,
            innovations,
            noise,
            measurements,
            code_design,
            row_satellites,
        })
    }
}

fn predict_kinematic_state(
    mean: &mut [f64],
    covariance: &mut Matrix,
    dt: f64,
    config: PppKinematicConfig,
    head_noise: [f64; 5],
    ambiguity_noise: f64,
) {
    let n = mean.len();
    let mut transition = zeros(n, n);
    for (i, row) in transition.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for i in 0..3 {
        mean[i] += dt * mean[5 + i];
        transition[i][5 + i] = dt;
    }
    *covariance = crate::filter::multiply(
        &crate::filter::multiply(&transition, covariance),
        &crate::filter::transpose(&transition),
    );
    for (i, row) in covariance.iter_mut().enumerate() {
        row[i] += dt
            * if i < 5 {
                head_noise[i]
            } else if i < 8 {
                config.velocity_process_noise_m2_s3
            } else {
                ambiguity_noise
            };
    }
}

fn predict_acceleration_state(
    mean: &mut [f64],
    covariance: &mut Matrix,
    dt: f64,
    config: PppAccelerationConfig,
    head_noise: [f64; 5],
    ambiguity_noise: f64,
) -> Result<(), PppError> {
    let n = mean.len();
    let mut transition = zeros(n, n);
    for (i, row) in transition.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for i in 0..3 {
        mean[i] += dt * mean[5 + i] + 0.5 * dt * dt * mean[8 + i];
        mean[5 + i] += dt * mean[8 + i];
        transition[i][5 + i] = dt;
        transition[i][8 + i] = 0.5 * dt * dt;
        transition[5 + i][8 + i] = dt;
    }
    *covariance = crate::filter::multiply(
        &crate::filter::multiply(&transition, covariance),
        &crate::filter::transpose(&transition),
    );
    let geo = ecef_to_geodetic(receiver(mean))?;
    let rotation = [
        crate::EnuCoord::new(1.0, 0.0, 0.0),
        crate::EnuCoord::new(0.0, 1.0, 0.0),
        crate::EnuCoord::new(0.0, 0.0, 1.0),
    ]
    .map(|e| crate::coordinates::enu_to_ecef(e, geo))
    .into_iter()
    .collect::<Result<Vec<_>, _>>()?;
    let sigmas = [
        config.horizontal_sigma_mps2_per_sqrt_s,
        config.horizontal_sigma_mps2_per_sqrt_s,
        config.vertical_sigma_mps2_per_sqrt_s,
    ];
    for (axis, sigma) in rotation.iter().zip(sigmas) {
        let v = [axis.x, axis.y, axis.z];
        for i in 0..3 {
            for j in 0..3 {
                covariance[8 + i][8 + j] += v[i] * v[j] * sigma * sigma * dt;
            }
        }
    }
    covariance[3][3] += head_noise[3] * dt;
    covariance[4][4] += head_noise[4] * dt;
    for (i, row) in covariance.iter_mut().enumerate().skip(11) {
        row[i] += ambiguity_noise * dt;
    }
    Ok(())
}

/// Pinned native reinitializeScalarState covariance operation, applied to
/// the clock slot while retaining its last committed seed mean.
fn reset_epoch_clock_covariance(covariance: &mut Matrix, variance: f64) {
    for row in covariance.iter_mut() {
        row[3] = 0.0;
    }
    covariance[3].fill(0.0);
    covariance[3][3] = variance;
}

type Layout = BTreeMap<SatelliteId, (usize, PhaseArcId)>;
type UpdateResult = (GaussianState, Layout, Vec<PppMeasurement>, usize, f64);
struct WindupEpochContext {
    model: PppWindupModel,
    sun: EcefCoord,
    anchors: BTreeMap<SatelliteId, (f64, f64)>,
}
struct MeasurementSystem {
    design: Matrix,
    innovations: Vec<f64>,
    noise: Matrix,
    measurements: Vec<PppMeasurement>,
    code_design: Vec<[f64; 5]>,
    row_satellites: Vec<(SatelliteId, bool)>,
}
fn receiver(mean: &[f64]) -> EcefCoord {
    EcefCoord::new(mean[0], mean[1], mean[2])
}
/// Modified Gram-Schmidt checks that the five code-state columns are independent;
/// a finite prior alone must not disguise rank-deficient satellite geometry.
fn ensure_geometry(rows: &[[f64; 5]]) -> Result<(), PppError> {
    let mut basis: Vec<Vec<f64>> = Vec::new();
    for column in 0..5 {
        let mut v: Vec<_> = rows.iter().map(|r| r[column]).collect();
        let scale = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        for _ in 0..2 {
            for q in &basis {
                let dot = v.iter().zip(q).map(|(a, b)| a * b).sum::<f64>();
                for (a, b) in v.iter_mut().zip(q) {
                    *a -= dot * b;
                }
            }
        }
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if !norm.is_finite() || norm <= scale.max(1.0) * 1e-10 {
            return Err(Error::SingularGeometry.into());
        }
        basis.push(v.into_iter().map(|x| x / norm).collect());
    }
    Ok(())
}

#[cfg(test)]
mod kinematic_prediction_tests {
    use super::*;
    #[test]
    fn acceleration_prediction_matches_48_actual_native_pva_calls() {
        let mut count = 0;
        for line in include_str!("../tests/fixtures/upstream_ppp_acceleration_predict.csv")
            .lines()
            .skip(1)
        {
            let d: Vec<f64> = line.split(',').map(|s| s.parse().unwrap()).collect();
            let n = d[1] as usize;
            assert_eq!(d.len(), 8 + 2 * n + 2 * n * n);
            let mut mean = d[8..8 + n].to_vec();
            let mut covariance: Matrix = d[8 + n..8 + n + n * n]
                .chunks(n)
                .map(|r| r.to_vec())
                .collect();
            predict_acceleration_state(
                &mut mean,
                &mut covariance,
                d[2],
                PppAccelerationConfig {
                    horizontal_sigma_mps2_per_sqrt_s: d[6],
                    vertical_sigma_mps2_per_sqrt_s: d[7],
                    ..Default::default()
                },
                [d[3], d[3], d[3], 100.0, d[4]],
                d[5],
            )
            .unwrap();
            let expected_mean = &d[8 + n + n * n..8 + 2 * n + n * n];
            let expected_covariance = &d[8 + 2 * n + n * n..];
            for (a, b) in mean
                .iter()
                .zip(expected_mean)
                .chain(covariance.iter().flatten().zip(expected_covariance))
            {
                assert!(
                    (a - b).abs() <= 1e-10 * b.abs().max(1.0),
                    "native case {}: {a} versus {b}",
                    d[0]
                );
            }
            crate::filter::cholesky(&covariance).unwrap();
            count += 1;
        }
        assert_eq!(count, 48);
    }

    #[test]
    fn pva_local_noise_rotates_variances_and_does_not_add_position_or_velocity_noise() {
        let mut mean = vec![0.0; 12];
        mean[0] = c::WGS84_A;
        mean[4] = 2.3;
        let mut covariance = zeros(12, 12);
        for (i, row) in covariance.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        predict_acceleration_state(
            &mut mean,
            &mut covariance,
            2.0,
            PppAccelerationConfig {
                horizontal_sigma_mps2_per_sqrt_s: 2.0,
                vertical_sigma_mps2_per_sqrt_s: 0.5,
                ..Default::default()
            },
            [1000.0, 1000.0, 1000.0, 2.0, 0.01],
            0.03,
        )
        .unwrap();
        for i in 0..3 {
            assert_eq!(covariance[i][i], 9.0);
            assert_eq!(covariance[5 + i][5 + i], 5.0);
            assert_eq!(covariance[i][5 + i], 6.0);
            assert_eq!(covariance[i][8 + i], 2.0);
            assert_eq!(covariance[5 + i][8 + i], 2.0);
        }
        assert_eq!(
            [covariance[8][8], covariance[9][9], covariance[10][10]],
            [1.5, 9.0, 9.0]
        );
        assert_eq!(covariance[3][3], 5.0);
        assert_eq!(covariance[4][4], 1.02);
        assert_eq!(covariance[11][11], 1.06);
        crate::filter::cholesky(&covariance).unwrap();
    }
    #[test]
    fn clock_covariance_reset_matches_all_16_actual_native_scalar_reinitializations() {
        let mut count = 0;
        for line in include_str!("../tests/fixtures/upstream_ppp_epoch_clock.csv")
            .lines()
            .skip(1)
        {
            let data: Vec<f64> = line.split(',').map(|s| s.parse().unwrap()).collect();
            let n = data[1] as usize;
            assert_eq!(data.len(), 3 + 2 * n * n);
            let mut covariance: Matrix = data[3..3 + n * n].chunks(n).map(|r| r.to_vec()).collect();
            reset_epoch_clock_covariance(&mut covariance, data[2]);
            for (a, b) in covariance.iter().flatten().zip(&data[3 + n * n..]) {
                assert!((a - b).abs() < 1e-10);
            }
            crate::filter::cholesky(&covariance).unwrap();
            count += 1;
        }
        assert_eq!(count, 16);
    }
    #[test]
    fn full_covariance_prediction_matches_all_24_actual_native_ppp_calls() {
        let mut count = 0;
        for line in include_str!("../tests/fixtures/upstream_ppp_kinematic_predict.csv")
            .lines()
            .skip(1)
        {
            let data: Vec<f64> = line.split(',').map(|v| v.parse().unwrap()).collect();
            let n = data[1] as usize;
            assert_eq!(data.len(), 7 + 2 * n + 2 * n * n);
            let (dt, qpos, qtrop, qamb, qvel) = (data[2], data[3], data[4], data[5], data[6]);
            let mut mean = data[7..7 + n].to_vec();
            let mut covariance: Matrix = data[7 + n..7 + n + n * n]
                .chunks(n)
                .map(|r| r.to_vec())
                .collect();
            let expected_mean = &data[7 + n + n * n..7 + 2 * n + n * n];
            let expected_covariance = &data[7 + 2 * n + n * n..];
            predict_kinematic_state(
                &mut mean,
                &mut covariance,
                dt,
                PppKinematicConfig {
                    velocity_process_noise_m2_s3: qvel,
                    ..Default::default()
                },
                [qpos, qpos, qpos, 100.0, qtrop],
                qamb,
            );
            for (actual, expected) in mean
                .iter()
                .zip(expected_mean)
                .chain(covariance.iter().flatten().zip(expected_covariance))
            {
                assert!(
                    (actual - expected).abs() <= 1e-10 * expected.abs().max(1.0),
                    "case {}: {actual} versus {expected}",
                    data[0]
                );
            }
            crate::filter::cholesky(&covariance).unwrap();
            count += 1;
        }
        assert_eq!(count, 24);
    }
}
