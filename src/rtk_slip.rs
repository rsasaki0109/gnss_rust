//! Single-frequency SD code/phase and Doppler slip evidence.
//!
//! Valid-input equations and strict threshold comparisons follow pinned
//! libgnss++ rtk_slip_detection.hpp. Invalid/missing input must not become a
//! numeric zero or a negative detection: this API returns explicit errors.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.

use crate::Error;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlipTest {
    /// Signed change, or observed phase change minus predicted change, in m.
    pub residual_m: f64,
    pub detected: bool,
}

fn finite(values: &[f64]) -> Result<(), Error> {
    if values.iter().all(|v| v.is_finite()) {
        Ok(())
    } else {
        Err(Error::NonFinite)
    }
}

fn positive(value: f64) -> Result<(), Error> {
    finite(&[value])?;
    if value > 0.0 {
        Ok(())
    } else {
        Err(Error::OutOfRange)
    }
}

pub fn single_difference_range_rate_mps(
    rover_doppler_hz: f64,
    base_doppler_hz: f64,
    wavelength_m: f64,
) -> Result<f64, Error> {
    finite(&[rover_doppler_hz, base_doppler_hz])?;
    positive(wavelength_m)?;
    let value = -(rover_doppler_hz - base_doppler_hz) * wavelength_m;
    finite(&[value])?;
    Ok(value)
}

pub fn single_difference_code_minus_phase_m(
    rover_code_m: f64,
    base_code_m: f64,
    rover_phase_cycles: f64,
    base_phase_cycles: f64,
    wavelength_m: f64,
) -> Result<f64, Error> {
    finite(&[
        rover_code_m,
        base_code_m,
        rover_phase_cycles,
        base_phase_cycles,
    ])?;
    positive(wavelength_m)?;
    let value =
        (rover_code_m - base_code_m) - (rover_phase_cycles - base_phase_cycles) * wavelength_m;
    finite(&[value])?;
    Ok(value)
}

pub fn test_code_minus_phase(
    previous_m: f64,
    current_m: f64,
    threshold_m: f64,
) -> Result<SlipTest, Error> {
    finite(&[previous_m, current_m])?;
    positive(threshold_m)?;
    let residual_m = current_m - previous_m;
    finite(&[residual_m])?;
    Ok(SlipTest {
        residual_m,
        detected: residual_m.abs() > threshold_m,
    })
}

pub fn test_doppler_phase(
    previous_sd_phase_m: f64,
    current_sd_phase_m: f64,
    sd_range_rate_mps: f64,
    dt_s: f64,
    threshold_m: f64,
) -> Result<SlipTest, Error> {
    finite(&[previous_sd_phase_m, current_sd_phase_m, sd_range_rate_mps])?;
    positive(dt_s)?;
    positive(threshold_m)?;
    let residual_m = (current_sd_phase_m - previous_sd_phase_m) - sd_range_rate_mps * dt_s;
    finite(&[residual_m])?;
    Ok(SlipTest {
        residual_m,
        detected: residual_m.abs() > threshold_m,
    })
}

/// Defaults match the pinned native GPS L1 thresholds, including its dynamic
/// floors. The filter enables this policy only when explicitly requested.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RtkSlipConfig {
    pub code_threshold_m: Option<f64>,
    pub doppler_threshold_m: Option<f64>,
    /// Native Doppler detection is restricted to intervals <= 5 seconds.
    pub max_doppler_interval_s: f64,
}

impl Default for RtkSlipConfig {
    fn default() -> Self {
        Self {
            code_threshold_m: Some(5.0),
            doppler_threshold_m: Some(0.20),
            max_doppler_interval_s: 5.0,
        }
    }
}

impl RtkSlipConfig {
    pub(crate) fn validate(self) -> Result<(), Error> {
        if self.code_threshold_m.is_none() && self.doppler_threshold_m.is_none() {
            return Err(Error::InvalidConfiguration);
        }
        for value in [self.code_threshold_m, self.doppler_threshold_m]
            .into_iter()
            .flatten()
        {
            positive(value)?;
        }
        positive(self.max_doppler_interval_s)
    }
}
