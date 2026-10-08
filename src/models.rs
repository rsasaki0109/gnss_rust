//! Initial atmosphere models ported from the upstream headers.

use crate::{Error, GeodeticCoord, constants::SPEED_OF_LIGHT};
use std::f64::consts::PI;

/// Niell (1996) hydrostatic mapping with seasonal/latitude interpolation and
/// height correction, ported from libgnss++ models/troposphere.hpp.
/// Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
/// Uses a 0.05 floor on sin(elevation), native coefficient clamps, and a
/// caller-supplied calendar day (1..=366). No wet/dry split is performed here.
pub fn niell_hydrostatic_mapping(
    position: GeodeticCoord,
    elevation: f64,
    day_of_year: u16,
) -> Result<f64, Error> {
    position.validate()?;
    if !elevation.is_finite()
        || !(0.0..=PI / 2.0).contains(&elevation)
        || !(1..=366).contains(&day_of_year)
    {
        return Err(Error::InvalidConfiguration);
    }
    let latitude = position.latitude.abs().to_degrees().clamp(0.0, 90.0);
    let interpolate = |values: [f64; 5]| {
        if latitude <= 15.0 {
            return values[0];
        }
        if latitude >= 75.0 {
            return values[4];
        }
        let i = ((latitude - 15.0) / 15.0).floor() as usize;
        let weight = (latitude - (15.0 + 15.0 * i as f64)) / 15.0;
        values[i] * (1.0 - weight) + values[i + 1] * weight
    };
    let phase = (2.0 * PI * (f64::from(day_of_year) - 28.0) / 365.25
        + if position.latitude < 0.0 { PI } else { 0.0 })
    .cos();
    let a = interpolate([
        1.2769934e-3,
        1.2683230e-3,
        1.2465397e-3,
        1.2196049e-3,
        1.2045996e-3,
    ]) - interpolate([0.0, 1.2709626e-5, 2.6523662e-5, 3.4000452e-5, 4.1202191e-5]) * phase;
    let b = interpolate([
        2.9153695e-3,
        2.9152299e-3,
        2.9288445e-3,
        2.9022565e-3,
        2.9024912e-3,
    ]) - interpolate([0.0, 2.1414979e-5, 3.0160779e-5, 7.2562722e-5, 11.723375e-5]) * phase;
    let c = interpolate([
        62.610505e-3,
        62.837393e-3,
        63.721774e-3,
        63.824265e-3,
        64.258455e-3,
    ]) - interpolate([0.0, 9.0128400e-5, 4.3497037e-5, 84.795348e-5, 170.37206e-5]) * phase;
    let sin = elevation.sin().max(0.05);
    let continued = |a: f64, b: f64, c: f64| {
        (1.0 + a / (1.0 + b / (1.0 + c))) / (sin + a / (sin + b / (sin + c)))
    };
    let mapping = continued(a, b, c)
        + (1.0 / sin - continued(2.53e-5, 5.49e-3, 1.14e-3)) * position.height.max(0.0) * 1e-3;
    if mapping.is_finite() {
        Ok(mapping)
    } else {
        Err(Error::NonFinite)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KlobucharParameters {
    pub alpha: [f64; 4],
    pub beta: [f64; 4],
}
impl Default for KlobucharParameters {
    fn default() -> Self {
        Self {
            alpha: [0.1118e-7, -0.7451e-8, -0.5961e-7, 0.1192e-6],
            beta: [0.1167e6, -0.2294e6, -0.1311e6, 0.1049e7],
        }
    }
}

/// L1 code ionospheric delay in metres. Elevation/azimuth are radians. The caller
/// chooses actual broadcast coefficients or the documented upstream defaults.
pub fn klobuchar(
    position: GeodeticCoord,
    azimuth: f64,
    elevation: f64,
    tow: f64,
    parameters: KlobucharParameters,
) -> Result<f64, Error> {
    position.validate()?;
    if ![azimuth, elevation, tow]
        .iter()
        .chain(parameters.alpha.iter())
        .chain(parameters.beta.iter())
        .all(|x| x.is_finite())
        || !(0.0..=PI / 2.0).contains(&elevation)
        || !(0.0..604_800.0).contains(&tow)
    {
        return Err(Error::InvalidConfiguration);
    }
    let phi_u = position.latitude / PI;
    let lam_u = position.longitude / PI;
    let el_sc = elevation / PI;
    let psi = 0.0137 / (el_sc + 0.11) - 0.022;
    let phi_i = (phi_u + psi * azimuth.cos()).clamp(-0.416, 0.416);
    let lam_i = lam_u + psi * azimuth.sin() / (phi_i * PI).cos();
    let phi_m = phi_i + 0.064 * ((lam_i - 1.617) * PI).cos();
    let local_time = (43_200.0 * lam_i + tow).rem_euclid(86_400.0);
    let obliquity = 1.0 + 16.0 * (0.53 - el_sc).powi(3);
    let a = parameters.alpha;
    let b = parameters.beta;
    let amplitude =
        (a[0] + a[1] * phi_m + a[2] * phi_m * phi_m + a[3] * phi_m * phi_m * phi_m).max(0.0);
    let period =
        (b[0] + b[1] * phi_m + b[2] * phi_m * phi_m + b[3] * phi_m * phi_m * phi_m).max(72_000.0);
    let x = 2.0 * PI * (local_time - 50_400.0) / period;
    let delay = SPEED_OF_LIGHT
        * obliquity
        * if x.abs() < 1.57 {
            5e-9 + amplitude * (1.0 - x * x / 2.0 + x.powi(4) / 24.0)
        } else {
            5e-9
        };
    if !delay.is_finite() {
        return Err(Error::NonFinite);
    }
    Ok(delay)
}

/// Standard-atmosphere Saastamoinen hydrostatic/wet delay with simple mapping.
/// Matches upstream: zero below elevation 0.05 rad or outside -100..10000 m.
pub fn saastamoinen(position: GeodeticCoord, elevation: f64) -> Result<f64, Error> {
    position.validate()?;
    if !elevation.is_finite() || !(-PI / 2.0..=PI / 2.0).contains(&elevation) {
        return Err(Error::InvalidConfiguration);
    }
    if elevation <= 0.05 || !(-100.0..=10_000.0).contains(&position.height) {
        return Ok(0.0);
    }
    let height = position.height.max(0.0);
    let temperature = 15.0 - 6.5e-3 * height + 273.16;
    let pressure = 1013.25 * (1.0 - 2.2557e-5 * height).powf(5.2568);
    let vapor = 6.108 * 0.7 * ((17.15 * temperature - 4684.0) / (temperature - 38.45)).exp();
    let mapping = (PI / 2.0 - elevation).cos().max(0.067);
    Ok(0.002_276_8 * pressure
        / (1.0 - 0.00266 * (2.0 * position.latitude).cos() - 0.00028 * height * 1e-3)
        / mapping
        + 0.002277 * (1255.0 / temperature + 0.05) * vapor / mapping)
}
