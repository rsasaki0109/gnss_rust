//! Initial static GPS ionosphere-free PPP FLOAT with precise products.
//! State is [ECEF XYZ, clock metres, zenith total delay metres, IF ambiguities].
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
use crate::dual_frequency::{DualFrequencyConfig, IonosphereFreeObservation, SignalCorrections};
use crate::filter::{FilterError, GaussianState, Matrix, measurement_update, zeros};
use crate::models::niell_hydrostatic_mapping;
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
pub struct PppFilterSettings {
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
    ProductFrameMismatch,
    UnsupportedSelection(SatelliteId),
    InsufficientCode { available: usize },
    InsufficientPhase { available: usize },
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
    fn validate(&self) -> Result<(), PppError> {
        self.check_references()?;
        if self.solid_earth_tide.is_some() {
            crate::tides::validate_station(self.initial_position)?;
        }
        let f = self.filter;
        ecef_to_geodetic(self.initial_position)?;
        self.transmit.validate()?;
        if !self.initial_receiver_clock_s.is_finite()
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
    /// XYZ state/report reference; model receiver_position_m is instantaneous.
    pub receiver_reference: ReceiverPositionReference,
    pub position: PositionSolution,
    pub zenith_delay_m: f64,
    pub zenith_variance_m2: f64,
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
    windup: BTreeMap<SatelliteId, PppWindupState>,
}
impl PppFilter {
    pub fn new(config: PppConfig) -> Result<Self, PppError> {
        config.validate()?;
        let arcs = PhaseArcTracker::new(config.arcs)?;
        let mut covariance = zeros(5, 5);
        for (i, row) in covariance.iter_mut().enumerate() {
            row[i] = config.filter.initial_head_variance_m2[i];
        }
        let mean = vec![
            config.initial_position.x,
            config.initial_position.y,
            config.initial_position.z,
            config.initial_receiver_clock_s * c::SPEED_OF_LIGHT,
            config.initial_zenith_delay_m,
        ];
        Ok(Self {
            config,
            state: GaussianState { mean, covariance },
            state_time: None,
            layout: BTreeMap::new(),
            arcs,
            last_rejections: BTreeMap::new(),
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
    pub fn windup_state(&self, sat: SatelliteId) -> Option<&PppWindupState> {
        self.windup.get(&sat)
    }
    pub fn invalidate_phase(&mut self) -> Result<(), PppError> {
        self.arcs.invalidate_all()?;
        self.windup.clear();
        Ok(())
    }

    /// Static GPS-only IF PPP. Every selected code/phase requires declared
    /// corrections, matching actual selected clock source, frame and reference
    /// point. At least five independent code rows and four phases are required.
    /// Failure preserves numeric mean/covariance and interrupts phase continuity.
    pub fn process_epoch(
        &mut self,
        epoch: &ObservationEpoch,
        products: &PreciseProducts,
        corrections: &SignalCorrections,
        selections: &BTreeMap<SatelliteId, DualFrequencyConfig>,
    ) -> Result<PppEpoch, PppError> {
        if self.arcs.last_epoch().is_some_and(|t| epoch.time <= t) {
            return Err(PhaseArcError::OutOfOrder.into());
        }
        if let Err(error) = self.config.check_products(products) {
            self.invalidate_phase()?;
            return Err(error);
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
        let prepared = self.arcs.prepare_epoch(epoch, corrections, selections);
        let decisions = match prepared {
            Ok(prepared) => prepared.decisions().clone(),
            Err(error) => {
                self.windup.clear();
                return Err(error.into());
            }
        };
        self.last_rejections.clear();
        if epoch.flag == 6 {
            self.windup.clear();
            return Ok(PppEpoch::Event(self.arcs.finish_epoch(&BTreeSet::new())?));
        }
        let result = self.update(epoch, products, &decisions);
        match result {
            Ok((posterior, layout, measurements, iterations, nis)) => {
                let eligible = layout.keys().copied().collect();
                let commit = self.arcs.finish_epoch(&eligible)?;
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
                    ..Default::default()
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
                    receiver_reference: self.config.receiver_reference,
                    position,
                    zenith_delay_m: posterior.mean[4],
                    zenith_variance_m2: posterior.covariance[4][4],
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
                self.arcs.finish_epoch(&BTreeSet::new())?;
                self.windup.clear();
                Err(error)
            }
        }
    }

    fn predicted_prior(
        &self,
        time: GnssTime,
        layout: &BTreeMap<SatelliteId, (usize, PhaseArcId)>,
    ) -> Result<GaussianState, PppError> {
        let n = 5 + layout.len();
        let mut mean = vec![0.0; n];
        let mut covariance = zeros(n, n);
        let mut source: Vec<Option<usize>> = (0..5).map(Some).collect();
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
            if i >= 5 && source[i].is_none() {
                covariance[i][i] = self.config.filter.initial_ambiguity_variance_m2;
            }
            covariance[i][i] += dt
                * if i < 5 {
                    self.config.filter.head_process_noise_m2_per_s[i]
                } else {
                    self.config.filter.ambiguity_process_noise_m2_per_s
                };
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
    ) -> Result<UpdateResult, PppError> {
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
                receiver(&self.state.mean),
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
                layout.insert(sat, (5 + layout.len(), id));
            }
            observations.insert(sat, observation.clone());
        }
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
        let prior = self.predicted_prior(epoch.time, &layout)?;
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
        Err(PppError::NonConvergent)
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
        })
    }
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
