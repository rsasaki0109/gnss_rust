//! Precise satellite states for code measurements at transmission time.
//!
//! Reception-time Taylor projection follows libgnss++ spp.cpp and
//! ppp_corrections.cpp (p -= v*tau, clock -= drift*tau, relativity once).
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
//! No broadcast orbit/clock, antenna or signal-bias fallback is performed.

use crate::precise::{
    ClockSource, InterpolationConfig, InterpolationError, PreciseProducts, PreciseState,
    precise_clock_relativistic_correction,
};
use crate::{EcefCoord, Error, GnssTime, SatelliteId, constants::SPEED_OF_LIGHT};
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PreciseTransmitMethod {
    /// Native upstream's receive-time state projected by one geometric tau.
    ReceptionTaylor,
    /// Reinterpolate at emission time and iterate Euclidean geometric tau.
    IteratedEmission,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PreciseTransmitConfig {
    pub interpolation: InterpolationConfig,
    pub method: PreciseTransmitMethod,
    pub max_iterations: usize,
    pub convergence_s: f64,
    pub max_light_time_s: f64,
}
impl Default for PreciseTransmitConfig {
    fn default() -> Self {
        Self {
            interpolation: Default::default(),
            method: PreciseTransmitMethod::ReceptionTaylor,
            max_iterations: 8,
            convergence_s: 1e-11,
            max_light_time_s: 1.0,
        }
    }
}
impl PreciseTransmitConfig {
    pub(crate) fn validate(self) -> Result<(), PreciseTransmitError> {
        self.interpolation.validate()?;
        if !(1..=64).contains(&self.max_iterations)
            || !self.convergence_s.is_finite()
            || self.convergence_s <= 0.0
            || !self.max_light_time_s.is_finite()
            || self.max_light_time_s <= 0.0
            || self.max_light_time_s > 1.0
        {
            return Err(PreciseTransmitError::InvalidConfiguration);
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PreciseTransmitError {
    InvalidConfiguration,
    Interpolation(InterpolationError),
    ClockUnavailable {
        source: ClockSource,
        reason: InterpolationError,
    },
    QualityBoundary,
    InvalidGeometry,
    NonConvergent,
    Gnss(Error),
}
impl fmt::Display for PreciseTransmitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfiguration => f.write_str("invalid precise transmit configuration"),
            Self::Interpolation(e) => e.fmt(f),
            Self::ClockUnavailable { source, reason } => {
                write!(f, "selected {source:?} clock unavailable: {reason}")
            }
            Self::QualityBoundary => {
                f.write_str("light-time interval crosses a precise-product quality boundary")
            }
            Self::InvalidGeometry => f.write_str("invalid precise light-time geometry"),
            Self::NonConvergent => f.write_str("precise emission-time iteration did not converge"),
            Self::Gnss(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for PreciseTransmitError {}
impl From<InterpolationError> for PreciseTransmitError {
    fn from(e: InterpolationError) -> Self {
        Self::Interpolation(e)
    }
}
impl From<Error> for PreciseTransmitError {
    fn from(e: Error) -> Self {
        Self::Gnss(e)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PreciseTransmitState {
    pub reception_time: GnssTime,
    pub emission_time: GnssTime,
    pub method: PreciseTransmitMethod,
    pub light_time_s: f64,
    pub iterations: usize,
    pub position_m: EcefCoord,
    pub velocity_m_per_s: [f64; 3],
    pub published_clock_bias_s: f64,
    /// Published-bias interpolant derivative, excluding periodic relativity.
    pub published_clock_drift_s_per_s: f64,
    pub periodic_relativity_s: f64,
    /// Published clock plus periodic relativity exactly once. No TGD/DCB/OSB.
    pub corrected_clock_bias_s: f64,
    /// State actually evaluated (reception for Taylor, emission otherwise).
    pub evaluated_state: PreciseState,
    /// Endpoint admission evidence, also for Taylor. Both endpoints must have
    /// usable orbit and clock; no quality boundary may lie between them.
    pub emission_state: PreciseState,
}

fn require_clock(s: &PreciseState) -> Result<(), PreciseTransmitError> {
    if s.clock.is_none() {
        return Err(PreciseTransmitError::ClockUnavailable {
            source: s.clock_source,
            reason: s
                .clock_unavailability
                .unwrap_or(InterpolationError::InsufficientSamples),
        });
    }
    Ok(())
}
fn light_time(p: EcefCoord, r: EcefCoord, max: f64) -> Result<f64, PreciseTransmitError> {
    let tau = (p.x - r.x).hypot(p.y - r.y).hypot(p.z - r.z) / SPEED_OF_LIGHT;
    if !tau.is_finite() || tau <= 0.0 || tau > max {
        return Err(PreciseTransmitError::InvalidGeometry);
    }
    Ok(tau)
}
fn check_interval(
    products: &PreciseProducts,
    satellite: SatelliteId,
    tx: GnssTime,
    rx: GnssTime,
    source: ClockSource,
) -> Result<(), PreciseTransmitError> {
    if products.orbit().samples()[&satellite].iter().any(|s| {
        s.time > tx
            && s.time <= rx
            && (s.flags.orbit_maneuver
                || s.position_m.is_none()
                || (source == ClockSource::Sp3
                    && (s.flags.clock_event || s.clock_bias_s.is_none())))
    }) {
        return Err(PreciseTransmitError::QualityBoundary);
    }
    Ok(())
}
/// Evaluate a precise state using a receiver position estimate. Tau uses
/// Euclidean range, matching upstream; apply Sagnac only in the subsequent
/// geometric-distance measurement model. The receiver may be zero for an
/// initial acquisition iteration. Both endpoint products must be admitted.
pub fn evaluate_precise_transmit(
    products: &PreciseProducts,
    satellite: SatelliteId,
    reception: GnssTime,
    receiver: EcefCoord,
    config: PreciseTransmitConfig,
) -> Result<PreciseTransmitState, PreciseTransmitError> {
    config.validate()?;
    if !receiver.is_finite() {
        return Err(PreciseTransmitError::InvalidGeometry);
    }
    let received = products.interpolate(satellite, reception, config.interpolation)?;
    require_clock(&received)?;
    let mut tau = light_time(received.position_m, receiver, config.max_light_time_s)?;
    let mut iterations = 1;
    let mut emission_time = reception.checked_add_seconds(-tau)?;
    check_interval(
        products,
        satellite,
        emission_time,
        reception,
        received.clock_source,
    )?;
    let mut emitted = products.interpolate(satellite, emission_time, config.interpolation)?;
    require_clock(&emitted)?;
    if config.method == PreciseTransmitMethod::IteratedEmission {
        let mut converged = false;
        for iteration in 1..=config.max_iterations {
            let next_tau = light_time(emitted.position_m, receiver, config.max_light_time_s)?;
            iterations = iteration;
            if (next_tau - tau).abs() <= config.convergence_s {
                converged = true;
                break;
            }
            tau = next_tau;
            emission_time = reception.checked_add_seconds(-tau)?;
            check_interval(
                products,
                satellite,
                emission_time,
                reception,
                received.clock_source,
            )?;
            emitted = products.interpolate(satellite, emission_time, config.interpolation)?;
            require_clock(&emitted)?;
        }
        if !converged {
            return Err(PreciseTransmitError::NonConvergent);
        }
    }
    let evaluated = if config.method == PreciseTransmitMethod::ReceptionTaylor {
        received
    } else {
        emitted.clone()
    };
    let clock = evaluated.clock.as_ref().unwrap();
    let (position, published_clock) = if config.method == PreciseTransmitMethod::ReceptionTaylor {
        let [vx, vy, vz] = evaluated.velocity_m_per_s;
        (
            EcefCoord::new(
                evaluated.position_m.x - vx * tau,
                evaluated.position_m.y - vy * tau,
                evaluated.position_m.z - vz * tau,
            ),
            clock.bias_s - clock.drift_s_per_s * tau,
        )
    } else {
        (evaluated.position_m, clock.bias_s)
    };
    let periodic = precise_clock_relativistic_correction(position, evaluated.velocity_m_per_s)?;
    let corrected = published_clock + periodic;
    if !position.is_finite() || !published_clock.is_finite() || !corrected.is_finite() {
        return Err(PreciseTransmitError::InvalidGeometry);
    }
    Ok(PreciseTransmitState {
        reception_time: reception,
        emission_time,
        method: config.method,
        light_time_s: tau,
        iterations,
        position_m: position,
        velocity_m_per_s: evaluated.velocity_m_per_s,
        published_clock_bias_s: published_clock,
        published_clock_drift_s_per_s: clock.drift_s_per_s,
        periodic_relativity_s: periodic,
        corrected_clock_bias_s: corrected,
        evaluated_state: evaluated,
        emission_state: emitted,
    })
}
