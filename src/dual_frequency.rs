//! Initial ionosphere-free code/phase observations, not a PPP estimator.
//!
//! Exact RINEX tracking identities are selected by the caller. Code and phase
//! corrections are additive metres in declared clock/phase datums; there is no
//! implicit TGD, antenna, wind-up or product-datum conversion. IF ambiguity is
//! a real-valued metre quantity, not an integer number of cycles.
//! Formula source: libgnss++ PPP utilities/observations at the pinned commit.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.

use crate::filter::cholesky;
use crate::observation::{Observation, ObservationEpoch};
use crate::precise::ClockSource;
use crate::{Error, GnssSystem, GnssTime, SatelliteId, constants as c};
use std::{collections::BTreeMap, fmt};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DualFrequencyError {
    InvalidConfiguration,
    UnsupportedTrackingCode,
    InvalidFrequencyPair,
    InvalidEpoch,
    MissingObservation { tracking_code: String },
    DuplicateObservation { tracking_code: String },
    InvalidObservation { tracking_code: String },
    MissingCode { tracking_code: String },
    InvalidCorrections,
    MissingCorrection { tracking_code: String },
    ExpiredCorrection { tracking_code: String },
    ClockDatumMismatch { tracking_code: String },
    PhaseDatumMismatch { tracking_code: String },
    InvalidNoise,
    Gnss(Error),
}
impl fmt::Display for DualFrequencyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "dual-frequency observation: {self:?}")
    }
}
impl std::error::Error for DualFrequencyError {}
impl From<Error> for DualFrequencyError {
    fn from(value: Error) -> Self {
        Self::Gnss(value)
    }
}

/// Frequency lookup for an explicit subset of RINEX 3 tracking codes.
/// GLONASS L1/L2 is FDMA and requires the observation's channel (-7..=6).
/// No GLONASS CDMA, Galileo E5ab, BeiDou B2ab, SBAS or NavIC mapping is implied.
pub fn tracking_frequency_hz(
    satellite: SatelliteId,
    tracking_code: &str,
    glonass_channel: Option<i8>,
) -> Result<f64, DualFrequencyError> {
    let bytes = tracking_code.as_bytes();
    if bytes.len() != 2 || !bytes[1].is_ascii_uppercase() {
        return Err(DualFrequencyError::UnsupportedTrackingCode);
    }
    let (frequency, tracking): (f64, &str) = match (satellite.system(), bytes[0]) {
        (GnssSystem::Gps, b'1') => (c::GPS_L1_FREQ, "CSLXPWYMN"),
        (GnssSystem::Gps, b'2') => (c::GPS_L2_FREQ, "CDSLXPWYMN"),
        (GnssSystem::Gps, b'5') => (c::GPS_L5_FREQ, "IQX"),
        (GnssSystem::Qzss, b'1') => (c::GPS_L1_FREQ, "CSLXZ"),
        (GnssSystem::Qzss, b'2') => (c::GPS_L2_FREQ, "SLX"),
        (GnssSystem::Qzss, b'5') => (c::GPS_L5_FREQ, "IQX"),
        (GnssSystem::Galileo, b'1') => (c::GAL_E1_FREQ, "ABCXZ"),
        (GnssSystem::Galileo, b'5') => (c::GAL_E5A_FREQ, "IQX"),
        (GnssSystem::Galileo, b'7') => (c::GAL_E5B_FREQ, "IQX"),
        (GnssSystem::Galileo, b'6') => (c::GAL_E6_FREQ, "ABCXZ"),
        (GnssSystem::BeiDou, b'2') => (c::BDS_B1I_FREQ, "IQX"),
        (GnssSystem::BeiDou, b'7') => (c::BDS_B2I_FREQ, "IQX"),
        (GnssSystem::BeiDou, b'6') => (c::BDS_B3I_FREQ, "IQX"),
        (GnssSystem::BeiDou, b'1') => (c::BDS_B1C_FREQ, "DPX"),
        (GnssSystem::BeiDou, b'5') => (c::BDS_B2A_FREQ, "DPX"),
        (GnssSystem::Glonass, band @ (b'1' | b'2')) => {
            let channel = glonass_channel.ok_or(Error::MissingGlonassChannel)?;
            if !(-7..=6).contains(&channel) {
                return Err(Error::InvalidGlonassChannel.into());
            }
            let frequency = if band == b'1' {
                c::GLO_L1_BASE_FREQ + f64::from(channel) * c::GLO_L1_STEP_FREQ
            } else {
                c::GLO_L2_BASE_FREQ + f64::from(channel) * c::GLO_L2_STEP_FREQ
            };
            (frequency, "CP")
        }
        _ => return Err(DualFrequencyError::UnsupportedTrackingCode),
    };
    if !tracking.as_bytes().contains(&bytes[1]) {
        return Err(DualFrequencyError::UnsupportedTrackingCode);
    }
    Ok(frequency)
}

/// Coefficients [f1²/(f1²-f2²), -f2²/(f1²-f2²)].
/// A degenerate pair is an error, never a single-frequency fallback.
pub fn ionosphere_free_coefficients(f1: f64, f2: f64) -> Result<[f64; 2], DualFrequencyError> {
    let denominator = f1 * f1 - f2 * f2;
    if !f1.is_finite()
        || !f2.is_finite()
        || f1 <= 0.0
        || f2 <= 0.0
        || (f1 - f2).abs() < 1.0
        || !denominator.is_finite()
        || denominator.abs() < 1.0
    {
        return Err(DualFrequencyError::InvalidFrequencyPair);
    }
    Ok([f1 * f1 / denominator, -f2 * f2 / denominator])
}

/// Melbourne-Wubbena in metres, with phases in cycles and codes in metres.
/// Order determines the sign of the wide-lane wavelength c/(f1-f2).
pub fn melbourne_wubbena_m(
    phase_cycles: [f64; 2],
    code_m: [f64; 2],
    frequency_hz: [f64; 2],
) -> Result<f64, DualFrequencyError> {
    let [f1, f2] = frequency_hz;
    ionosphere_free_coefficients(f1, f2)?;
    if !phase_cycles.iter().chain(&code_m).all(|v| v.is_finite()) {
        return Err(Error::NonFinite.into());
    }
    finite(
        (phase_cycles[0] - phase_cycles[1]) * c::SPEED_OF_LIGHT / (f1 - f2)
            - (f1 * code_m[0] + f2 * code_m[1]) / (f1 + f2),
    )
}
fn finite(value: f64) -> Result<f64, DualFrequencyError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::NonFinite.into())
    }
}

/// Caller-declared additive per-signal corrections in metres, inclusive GPS
/// reception-time validity. Zero corrections must also be declared explicitly.
/// Names are caller declarations, not identities verified from product files.
#[derive(Debug, Clone, PartialEq)]
pub struct SignalCorrection {
    pub satellite: SatelliteId,
    pub tracking_code: String,
    pub code_add_m: f64,
    pub phase_add_m: Option<f64>,
    pub clock_source: ClockSource,
    pub clock_reference: String,
    /// Required whenever phase_add_m is present. Includes any OSB gauge choice.
    pub phase_reference: Option<String>,
    pub valid_from: GnssTime,
    pub valid_until: GnssTime,
}

#[derive(Debug, Clone)]
pub struct SignalCorrections {
    entries: BTreeMap<(SatelliteId, String), Vec<SignalCorrection>>,
}
impl SignalCorrections {
    pub fn new(entries: Vec<SignalCorrection>) -> Result<Self, DualFrequencyError> {
        let mut table: BTreeMap<_, Vec<SignalCorrection>> = BTreeMap::new();
        for entry in entries {
            if !entry.code_add_m.is_finite()
                || entry.phase_add_m.is_some_and(|v| !v.is_finite())
                || entry.clock_reference.trim().is_empty()
                || entry.valid_until < entry.valid_from
                || entry.tracking_code.len() != 2
                || !matches!(entry.tracking_code.as_bytes()[0], b'1'..=b'9')
                || !entry.tracking_code.as_bytes()[1].is_ascii_uppercase()
                || match (&entry.phase_reference, entry.phase_add_m) {
                    (Some(name), Some(_)) => name.trim().is_empty(),
                    (None, None) => false,
                    _ => true,
                }
            {
                return Err(DualFrequencyError::InvalidCorrections);
            }
            let rows = table
                .entry((entry.satellite, entry.tracking_code.clone()))
                .or_default();
            if rows.iter().any(|previous| {
                previous.clock_source == entry.clock_source
                    && previous.valid_from <= entry.valid_until
                    && entry.valid_from <= previous.valid_until
            }) {
                return Err(DualFrequencyError::InvalidCorrections);
            }
            rows.push(entry);
        }
        Ok(Self { entries: table })
    }
    pub fn entries(&self) -> &BTreeMap<(SatelliteId, String), Vec<SignalCorrection>> {
        &self.entries
    }
    fn select(
        &self,
        satellite: SatelliteId,
        code: &str,
        time: GnssTime,
        config: &DualFrequencyConfig,
    ) -> Result<&SignalCorrection, DualFrequencyError> {
        let rows = self
            .entries
            .get(&(satellite, code.to_owned()))
            .ok_or_else(|| DualFrequencyError::MissingCorrection {
                tracking_code: code.to_owned(),
            })?;
        let current: Vec<_> = rows
            .iter()
            .filter(|b| b.valid_from <= time && time <= b.valid_until)
            .collect();
        if current.is_empty() {
            return Err(DualFrequencyError::ExpiredCorrection {
                tracking_code: code.to_owned(),
            });
        }
        current
            .into_iter()
            .find(|b| {
                b.clock_source == config.clock_source && b.clock_reference == config.clock_reference
            })
            .ok_or_else(|| DualFrequencyError::ClockDatumMismatch {
                tracking_code: code.to_owned(),
            })
    }
}

/// Metre² covariance of each frequency's code/phase. Correlations BETWEEN
/// frequencies are retained. Code/phase cross covariance is not represented.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DualFrequencyNoise {
    pub code_m2: [[f64; 2]; 2],
    pub phase_m2: [[f64; 2]; 2],
}
fn combined_variance(
    covariance: [[f64; 2]; 2],
    weights: [f64; 2],
) -> Result<f64, DualFrequencyError> {
    if covariance[0][1] != covariance[1][0] {
        return Err(DualFrequencyError::InvalidNoise);
    }
    let matrix: Vec<Vec<f64>> = covariance.iter().map(|r| r.to_vec()).collect();
    let lower = cholesky(&matrix).map_err(|_| DualFrequencyError::InvalidNoise)?;
    // ||L^T w||² avoids cancellation for correlated input noise.
    let first = lower[0][0] * weights[0] + lower[1][0] * weights[1];
    let second = lower[1][1] * weights[1];
    let variance = first * first + second * second;
    if !variance.is_finite() || variance <= 0.0 {
        return Err(DualFrequencyError::InvalidNoise);
    }
    Ok(variance)
}

#[derive(Debug, Clone, PartialEq)]
pub struct DualFrequencyConfig {
    /// Exact tracking identities; no preference/fallback search is performed.
    pub tracking_codes: [String; 2],
    pub clock_source: ClockSource,
    pub clock_reference: String,
    pub phase_reference: String,
    pub noise: DualFrequencyNoise,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PhaseUnavailable {
    MissingPhase { tracking_code: String },
    HalfCycle { tracking_code: String },
    MissingCorrection { tracking_code: String },
}
#[derive(Debug, Clone, PartialEq)]
pub struct IonosphereFreePhase {
    pub corrected_phase_m: [f64; 2],
    pub phase_if_m: f64,
    pub variance_m2: f64,
    /// Raw, uncorrected L1*lambda1 - L2*lambda2 diagnostic.
    pub raw_geometry_free_m: f64,
    /// Raw, uncorrected MW diagnostic. GLONASS is excluded, like native PPP.
    pub raw_melbourne_wubbena_m: Option<f64>,
    /// A valid numeric phase can start a new arc. Never continue old ambiguity
    /// across either signal's LLI loss-of-lock or a power-failure epoch.
    pub reset_required: bool,
}
#[derive(Debug, Clone, PartialEq)]
pub struct IonosphereFreeObservation {
    pub time: GnssTime,
    pub satellite: SatelliteId,
    pub raw: [Observation; 2],
    pub frequency_hz: [f64; 2],
    pub wavelength_m: [f64; 2],
    pub coefficients: [f64; 2],
    pub corrections: [SignalCorrection; 2],
    pub corrected_code_m: [f64; 2],
    pub code_if_m: f64,
    pub code_variance_m2: f64,
    pub phase: Option<IonosphereFreePhase>,
    pub phase_unavailability: Option<PhaseUnavailable>,
}

/// Construct one satellite's IF measurement at a reception epoch. This is a
/// stateless measurement builder. [`crate::phase_arc::PhaseArcTracker`] manages
/// missing/rejected phase, time gaps, signal/datum changes and failed updates.
/// It must also bind config.clock_source/reference to the actual satellite clock.
pub fn form_ionosphere_free(
    epoch: &ObservationEpoch,
    satellite: SatelliteId,
    corrections: &SignalCorrections,
    config: &DualFrequencyConfig,
) -> Result<IonosphereFreeObservation, DualFrequencyError> {
    if config.clock_reference.trim().is_empty() || config.phase_reference.trim().is_empty() {
        return Err(DualFrequencyError::InvalidConfiguration);
    }
    if !matches!(epoch.flag, 0 | 1)
        || epoch
            .receiver_clock_offset_s
            .is_some_and(|v| v != 0.0 || !v.is_finite())
    {
        return Err(DualFrequencyError::InvalidEpoch);
    }
    let mut observations = Vec::new();
    for code in &config.tracking_codes {
        let matches: Vec<_> = epoch
            .observations
            .iter()
            .filter(|o| o.satellite == satellite && o.tracking_code == *code)
            .collect();
        match matches.as_slice() {
            [] => {
                return Err(DualFrequencyError::MissingObservation {
                    tracking_code: code.clone(),
                });
            }
            [observation] => observations.push((*observation).clone()),
            _ => {
                return Err(DualFrequencyError::DuplicateObservation {
                    tracking_code: code.clone(),
                });
            }
        }
    }
    let raw = [observations.remove(0), observations.remove(0)];
    let mut frequency_hz = [0.0; 2];
    let mut code_m = [0.0; 2];
    for (i, observation) in raw.iter().enumerate() {
        frequency_hz[i] = tracking_frequency_hz(
            satellite,
            &observation.tracking_code,
            observation.glonass_channel,
        )?;
        let code = observation
            .pseudorange_m
            .ok_or_else(|| DualFrequencyError::MissingCode {
                tracking_code: observation.tracking_code.clone(),
            })?;
        if !code.is_finite()
            || code <= 0.0
            || observation
                .carrier_phase_cycles
                .is_some_and(|v| !v.is_finite())
            || observation.lli.is_some_and(|v| v > 7)
        {
            return Err(DualFrequencyError::InvalidObservation {
                tracking_code: observation.tracking_code.clone(),
            });
        }
        code_m[i] = code;
    }
    if satellite.system() == GnssSystem::Glonass && raw[0].glonass_channel != raw[1].glonass_channel
    {
        return Err(DualFrequencyError::InvalidFrequencyPair);
    }
    let coefficients = ionosphere_free_coefficients(frequency_hz[0], frequency_hz[1])?;
    let wavelength_m = frequency_hz.map(|f| c::SPEED_OF_LIGHT / f);
    let code_variance_m2 = combined_variance(config.noise.code_m2, coefficients)?;
    let phase_variance_m2 = combined_variance(config.noise.phase_m2, coefficients)?;
    let corrections = [
        corrections
            .select(satellite, &config.tracking_codes[0], epoch.time, config)?
            .clone(),
        corrections
            .select(satellite, &config.tracking_codes[1], epoch.time, config)?
            .clone(),
    ];
    for correction in &corrections {
        if correction
            .phase_reference
            .as_ref()
            .is_some_and(|name| *name != config.phase_reference)
        {
            return Err(DualFrequencyError::PhaseDatumMismatch {
                tracking_code: correction.tracking_code.clone(),
            });
        }
    }
    let corrected_code_m = [
        finite(code_m[0] + corrections[0].code_add_m)?,
        finite(code_m[1] + corrections[1].code_add_m)?,
    ];
    if corrected_code_m.iter().any(|p| *p <= 0.0) {
        return Err(Error::OutOfRange.into());
    }
    let combine =
        |values: [f64; 2]| finite(coefficients[0] * values[0] + coefficients[1] * values[1]);
    let mut phase_unavailability = None;
    for i in 0..2 {
        let code = raw[i].tracking_code.clone();
        let reason = if raw[i].carrier_phase_cycles.is_none() {
            Some(PhaseUnavailable::MissingPhase {
                tracking_code: code,
            })
        } else if raw[i].half_cycle_ambiguity() {
            Some(PhaseUnavailable::HalfCycle {
                tracking_code: code,
            })
        } else if corrections[i].phase_add_m.is_none() {
            Some(PhaseUnavailable::MissingCorrection {
                tracking_code: code,
            })
        } else {
            None
        };
        if reason.is_some() {
            phase_unavailability = reason;
            break;
        }
    }
    let phase = if phase_unavailability.is_none() {
        let cycles = [
            raw[0].carrier_phase_cycles.unwrap(),
            raw[1].carrier_phase_cycles.unwrap(),
        ];
        let phase_m = [
            finite(cycles[0] * wavelength_m[0])?,
            finite(cycles[1] * wavelength_m[1])?,
        ];
        let corrected_phase_m = [
            finite(phase_m[0] + corrections[0].phase_add_m.unwrap())?,
            finite(phase_m[1] + corrections[1].phase_add_m.unwrap())?,
        ];
        Some(IonosphereFreePhase {
            corrected_phase_m,
            phase_if_m: combine(corrected_phase_m)?,
            variance_m2: phase_variance_m2,
            raw_geometry_free_m: finite(phase_m[0] - phase_m[1])?,
            raw_melbourne_wubbena_m: if satellite.system() == GnssSystem::Glonass {
                None
            } else {
                Some(melbourne_wubbena_m(cycles, code_m, frequency_hz)?)
            },
            reset_required: epoch.flag == 1 || raw.iter().any(Observation::loss_of_lock),
        })
    } else {
        None
    };
    Ok(IonosphereFreeObservation {
        time: epoch.time,
        satellite,
        raw,
        frequency_hz,
        wavelength_m,
        coefficients,
        corrections,
        corrected_code_m,
        code_if_m: combine(corrected_code_m)?,
        code_variance_m2,
        phase,
        phase_unavailability,
    })
}
