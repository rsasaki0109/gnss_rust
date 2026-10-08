//! Initial GPS/QZSS L1 code positioning with precise transmit states.
//! Requires caller-declared signal/clock-datum corrections, never broadcast TGD.
//! This is code SPP, not phase PPP, and does not implement antenna corrections.

use crate::coordinates::{ecef_to_enu, ecef_to_geodetic, geometric_distance};
use crate::models::{KlobucharParameters, klobuchar, saastamoinen};
use crate::observation::ObservationEpoch;
use crate::precise::{ClockSource, OrbitReferencePoint, PreciseProducts};
use crate::precise_transmit::{
    PreciseTransmitConfig, PreciseTransmitError, PreciseTransmitState, evaluate_precise_transmit,
};
use crate::spp::{self, RejectionReason, Row, SppConfig};
use crate::{EcefCoord, Error, GnssTime, PositionSolution, SatelliteId, constants as c};
use std::{collections::BTreeMap, fmt};

#[derive(Debug, Clone, PartialEq)]
pub struct PreciseCodeBias {
    pub satellite: SatelliteId,
    pub tracking_code: String,
    /// P_model_datum = P_raw + correction_m. No implicit sign inversion.
    pub correction_m: f64,
    pub clock_source: ClockSource,
    /// Caller-declared identity of the clock datum (not verified from a file).
    pub clock_reference: String,
    pub valid_from: GnssTime,
    pub valid_until: GnssTime,
}
#[derive(Debug, Clone)]
pub struct PreciseCodeBiases {
    entries: BTreeMap<(SatelliteId, String), Vec<PreciseCodeBias>>,
}
impl PreciseCodeBiases {
    pub fn new(entries: Vec<PreciseCodeBias>) -> Result<Self, PreciseSppError> {
        let mut table: BTreeMap<(SatelliteId, String), Vec<PreciseCodeBias>> = BTreeMap::new();
        for bias in entries {
            if !bias.correction_m.is_finite()
                || bias.clock_reference.trim().is_empty()
                || bias.valid_until < bias.valid_from
                || bias.tracking_code.len() != 2
                || !bias.tracking_code.as_bytes()[0].is_ascii_digit()
                || !bias.tracking_code.as_bytes()[1].is_ascii_uppercase()
            {
                return Err(PreciseSppError::InvalidBiasTable);
            }
            let rows = table
                .entry((bias.satellite, bias.tracking_code.clone()))
                .or_default();
            if rows.iter().any(|b| {
                b.clock_source == bias.clock_source
                    && b.valid_from <= bias.valid_until
                    && bias.valid_from <= b.valid_until
            }) {
                return Err(PreciseSppError::InvalidBiasTable);
            }
            rows.push(bias);
        }
        Ok(Self { entries: table })
    }
    pub fn entries(&self) -> &BTreeMap<(SatelliteId, String), Vec<PreciseCodeBias>> {
        &self.entries
    }
    fn select(
        &self,
        sat: SatelliteId,
        tracking: &str,
        source: ClockSource,
        time: GnssTime,
        reference: &str,
    ) -> Result<&PreciseCodeBias, PreciseSppRejection> {
        let rows = self
            .entries
            .get(&(sat, tracking.to_owned()))
            .ok_or(PreciseSppRejection::MissingCodeBias)?;
        if !rows.iter().any(|b| b.clock_source == source) {
            return Err(PreciseSppRejection::BiasClockSourceMismatch);
        }
        let bias = rows
            .iter()
            .find(|b| b.clock_source == source && b.valid_from <= time && time <= b.valid_until)
            .ok_or(PreciseSppRejection::ExpiredCodeBias)?;
        if bias.clock_reference != reference {
            return Err(PreciseSppRejection::BiasClockReferenceMismatch);
        }
        Ok(bias)
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct PreciseSppConfig {
    pub spp: SppConfig,
    pub transmit: PreciseTransmitConfig,
    pub expected_coordinate_system: String,
    pub expected_reference_point: OrbitReferencePoint,
    pub clock_reference: String,
    /// Required when spp.use_ionosphere is true. No broadcast/default lookup.
    pub ionosphere: Option<KlobucharParameters>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PreciseSppError {
    InvalidConfiguration,
    InvalidBiasTable,
    ProductFrameMismatch,
    InvalidEpoch,
    Gnss(Error),
    Transmit(PreciseTransmitError),
}
impl fmt::Display for PreciseSppError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
        Self::InvalidConfiguration => f.write_str("invalid precise code SPP configuration"),
        Self::InvalidBiasTable => f.write_str("invalid or overlapping precise code bias table"),
        Self::ProductFrameMismatch => f.write_str("precise coordinate system/reference point differs from the declared model"),
        Self::InvalidEpoch => f.write_str("precise code SPP requires normal/power-failure epoch and zero unapplied receiver clock offset"),
        Self::Gnss(e) => e.fmt(f), Self::Transmit(e) => e.fmt(f),
    }
    }
}
impl std::error::Error for PreciseSppError {}
impl From<Error> for PreciseSppError {
    fn from(e: Error) -> Self {
        Self::Gnss(e)
    }
}
impl From<PreciseTransmitError> for PreciseSppError {
    fn from(e: PreciseTransmitError) -> Self {
        Self::Transmit(e)
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PreciseSppRejection {
    Observation(RejectionReason),
    Transmit(PreciseTransmitError),
    MissingCodeBias,
    ExpiredCodeBias,
    BiasClockSourceMismatch,
    BiasClockReferenceMismatch,
}
#[derive(Debug, Clone, PartialEq)]
pub struct PreciseCodeMeasurement {
    pub satellite: SatelliteId,
    pub tracking_code: String,
    pub raw_pseudorange_m: f64,
    pub corrected_pseudorange_m: f64,
    pub applied_bias: PreciseCodeBias,
    pub transmit: PreciseTransmitState,
}
#[derive(Debug, Clone, PartialEq)]
pub struct PreciseSppSolution {
    pub position: PositionSolution,
    pub iterations: usize,
    pub residual_rms_m: f64,
    pub residuals_m: Vec<(SatelliteId, f64)>,
    pub rejected: Vec<(SatelliteId, PreciseSppRejection)>,
    pub measurements: Vec<PreciseCodeMeasurement>,
    pub coordinate_system: String,
    pub reference_point: OrbitReferencePoint,
}

/// Precise code-only SPP. Every admitted tracking code requires an explicit
/// correction (including zero), clock source/datum and validity interval.
/// Coordinates are in the declared SP3 frame; no frame/antenna transform runs.
pub fn solve_epoch_precise(
    epoch: &ObservationEpoch,
    products: &PreciseProducts,
    biases: &PreciseCodeBiases,
    config: &PreciseSppConfig,
) -> Result<PreciseSppSolution, PreciseSppError> {
    config.spp.validate()?;
    config.transmit.validate()?;
    if config.clock_reference.trim().is_empty()
        || config.expected_coordinate_system.trim().is_empty()
        || (config.spp.use_ionosphere && config.ionosphere.is_none())
        || config
            .ionosphere
            .is_some_and(|p| !p.alpha.iter().chain(&p.beta).all(|x| x.is_finite()))
    {
        return Err(PreciseSppError::InvalidConfiguration);
    }
    if !matches!(epoch.flag, 0 | 1)
        || epoch
            .receiver_clock_offset_s
            .is_some_and(|v| !v.is_finite() || v != 0.0)
    {
        return Err(PreciseSppError::InvalidEpoch);
    }
    let header = products.orbit().header();
    if header.coordinate_system != config.expected_coordinate_system
        || header.reference_point != config.expected_reference_point
    {
        return Err(PreciseSppError::ProductFrameMismatch);
    }
    let (selected, initial_rejections) = spp::select_observations(epoch)?;
    let mut measurements = Vec::new();
    let mut rejected = Vec::new();
    let solution = spp::fit_epoch::<PreciseSppError>(
        epoch,
        config.spp,
        initial_rejections,
        |receiver, clock_m| {
            measurements.clear();
            rejected.clear();
            let geo = ecef_to_geodetic(receiver).ok();
            let mut rows = Vec::new();
            for (&satellite, obs) in &selected {
                let transmit = match evaluate_precise_transmit(
                    products,
                    satellite,
                    epoch.time,
                    receiver,
                    config.transmit,
                ) {
                    Ok(s) => s,
                    Err(e) => {
                        rejected.push((satellite, PreciseSppRejection::Transmit(e)));
                        continue;
                    }
                };
                let bias = match biases.select(
                    satellite,
                    &obs.tracking_code,
                    transmit.evaluated_state.clock_source,
                    epoch.time,
                    &config.clock_reference,
                ) {
                    Ok(b) => b,
                    Err(e) => {
                        rejected.push((satellite, e));
                        continue;
                    }
                };
                let raw = obs.pseudorange_m.unwrap();
                let corrected = raw + bias.correction_m;
                if !corrected.is_finite() || corrected <= 0.0 {
                    return Err(Error::NonFinite.into());
                }
                let p = transmit.position_m;
                let difference =
                    EcefCoord::new(p.x - receiver.x, p.y - receiver.y, p.z - receiver.z);
                let norm = difference.x.hypot(difference.y).hypot(difference.z);
                if !norm.is_finite() || norm <= 0.0 {
                    return Err(Error::SingularGeometry.into());
                }
                let mut elevation = std::f64::consts::FRAC_PI_2;
                let mut atmosphere = 0.0;
                if let Some(geo) = geo {
                    let enu = ecef_to_enu(difference, geo)?;
                    elevation = enu.up.atan2(enu.east.hypot(enu.north));
                    if elevation < config.spp.elevation_mask_rad {
                        rejected.push((
                            satellite,
                            PreciseSppRejection::Observation(RejectionReason::BelowElevationMask),
                        ));
                        continue;
                    }
                    if config.spp.use_troposphere {
                        atmosphere += saastamoinen(geo, elevation)?;
                    }
                    if config.spp.use_ionosphere {
                        atmosphere += klobuchar(
                            geo,
                            enu.east.atan2(enu.north),
                            elevation,
                            epoch.time.tow(),
                            config.ionosphere.unwrap(),
                        )?;
                    }
                }
                let predicted = geometric_distance(p, receiver)? + clock_m
                    - c::SPEED_OF_LIGHT * transmit.corrected_clock_bias_s
                    + atmosphere;
                let residual = corrected - predicted;
                if !residual.is_finite() {
                    return Err(Error::NonFinite.into());
                }
                rows.push(Row {
                    satellite,
                    design: [
                        -difference.x / norm - c::OMEGA_E * p.y / c::SPEED_OF_LIGHT,
                        -difference.y / norm + c::OMEGA_E * p.x / c::SPEED_OF_LIGHT,
                        -difference.z / norm,
                        1.0,
                    ],
                    residual,
                    weight: elevation.sin().max(0.1) / config.spp.pseudorange_sigma_m,
                });
                measurements.push(PreciseCodeMeasurement {
                    satellite,
                    tracking_code: obs.tracking_code.clone(),
                    raw_pseudorange_m: raw,
                    corrected_pseudorange_m: corrected,
                    applied_bias: bias.clone(),
                    transmit,
                });
            }
            if rows.len() < 4 {
                return Err(Error::InsufficientObservations.into());
            }
            Ok((rows, Vec::new()))
        },
    )?;
    let mut all_rejected: Vec<_> = solution
        .rejected
        .into_iter()
        .map(|(s, r)| (s, PreciseSppRejection::Observation(r)))
        .collect();
    all_rejected.extend(rejected);
    Ok(PreciseSppSolution {
        position: solution.position,
        iterations: solution.iterations,
        residual_rms_m: solution.residual_rms_m,
        residuals_m: solution.residuals_m,
        rejected: all_rejected,
        measurements,
        coordinate_system: header.coordinate_system.clone(),
        reference_point: header.reference_point,
    })
}
