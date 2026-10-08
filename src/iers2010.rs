//! IERS 2010 Dehant solid Earth tide, Step-1 + Step-2, explicit-body kernel.
//! All positions must share an Earth-fixed frame, in metres. This module does
//! not generate Sun/Moon ephemerides or convert GPS/UTC/TT/UT1 or ICRS/ITRS.
//! Permanent-tide Step-3 is omitted, matching the pinned native kernel.
//!
//! SPDX-License-Identifier: Apache-2.0
//! Copyright 2020 Geoscience Australia. Copyright 2024 LibGNSS++ Contributors.
//! Rust translation of ginan-iers2010/dehanttideinel/dehanttide_all.cpp;
//! modifications: arrays instead of Eigen, one station per call, explicit
//! validated arguments, decomposed results, errors for singular/invalid input.
//! See LICENSE-IERS2010 and NOTICE-IERS2010 for provenance and attribution.

use crate::{EcefCoord, constants};
use std::fmt;

type Vector = [f64; 3];
const RADIUS: f64 = 6_378_136.6; // IERS radius; not WGS84_A.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IersTideError {
    NonFinite,
    InvalidArguments,
    InvalidStation,
    AxialStation,
    InvalidBody,
}
impl fmt::Display for IersTideError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "IERS 2010 solid Earth tide: {self:?}")
    }
}
impl std::error::Error for IersTideError {}

/// TT centuries from J2000 and fractional UT hour [0,24). Native wrapper
/// approximates UT by UTC for this kernel; neither field is a GPS clock label.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IersTideArguments {
    julian_centuries_tt: f64,
    fractional_hour_ut: f64,
}
impl IersTideArguments {
    pub fn new(julian_centuries_tt: f64, fractional_hour_ut: f64) -> Result<Self, IersTideError> {
        if !julian_centuries_tt.is_finite() || !fractional_hour_ut.is_finite() {
            return Err(IersTideError::NonFinite);
        }
        if !(0.0..24.0).contains(&fractional_hour_ut) {
            return Err(IersTideError::InvalidArguments);
        }
        Ok(Self {
            julian_centuries_tt,
            fractional_hour_ut,
        })
    }
    pub fn julian_centuries_tt(self) -> f64 {
        self.julian_centuries_tt
    }
    pub fn fractional_hour_ut(self) -> f64 {
        self.fractional_hour_ut
    }
}

/// Geocentric Sun and Moon in the station's Earth-fixed frame (m), evaluated
/// at the epoch of the supplied time arguments. Caller guarantees the frame,
/// units and epoch; magnitudes alone cannot verify an ICRS→ITRS rotation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EarthFixedSunMoon {
    sun_m: EcefCoord,
    moon_m: EcefCoord,
}
impl EarthFixedSunMoon {
    pub fn new(sun_m: EcefCoord, moon_m: EcefCoord) -> Result<Self, IersTideError> {
        for body in [sun_m, moon_m] {
            let norm = checked_norm(vector(body))?;
            if norm <= 1.0 {
                return Err(IersTideError::InvalidBody);
            }
        }
        Ok(Self { sun_m, moon_m })
    }
    pub fn sun_m(self) -> EcefCoord {
        self.sun_m
    }
    pub fn moon_m(self) -> EcefCoord {
        self.moon_m
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct IersSolidEarthTideCorrection {
    pub arguments: IersTideArguments,
    pub bodies: EarthFixedSunMoon,
    pub nominal_marker_m: EcefCoord,
    pub degree2_degree3_m: EcefCoord,
    pub step1_diurnal_m: EcefCoord,
    pub step1_semidiurnal_m: EcefCoord,
    pub step1_latitude_m: EcefCoord,
    pub step2_diurnal_m: EcefCoord,
    pub step2_long_period_m: EcefCoord,
    /// Add exactly once to nominal_marker_m. No permanent-tide conversion.
    pub displacement_m: EcefCoord,
    pub instantaneous_marker_m: EcefCoord,
}

fn vector(v: EcefCoord) -> Vector {
    [v.x, v.y, v.z]
}
fn coord(v: Vector) -> EcefCoord {
    EcefCoord::new(v[0], v[1], v[2])
}
fn checked_norm(v: Vector) -> Result<f64, IersTideError> {
    if v.iter().any(|x| !x.is_finite()) {
        return Err(IersTideError::NonFinite);
    }
    let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if !norm.is_finite() {
        return Err(IersTideError::NonFinite);
    }
    Ok(norm)
}
fn add(a: Vector, b: Vector) -> Vector {
    std::array::from_fn(|i| a[i] + b[i])
}

struct Site {
    position: Vector,
    radius: f64,
    sin_phi: f64,
    cos_phi: f64,
    sin_lambda: f64,
    cos_lambda: f64,
}
impl Site {
    fn new(position: Vector) -> Result<Self, IersTideError> {
        let radius = checked_norm(position)?;
        if radius < 0.5 * constants::WGS84_A {
            return Err(IersTideError::InvalidStation);
        }
        let cos_phi = (position[0] * position[0] + position[1] * position[1]).sqrt() / radius;
        if cos_phi == 0.0 {
            return Err(IersTideError::AxialStation);
        }
        Ok(Self {
            position,
            radius,
            sin_phi: position[2] / radius,
            cos_phi,
            sin_lambda: position[1] / cos_phi / radius,
            cos_lambda: position[0] / cos_phi / radius,
        })
    }
    /// Geocentric local radial, north, east → station ECEF.
    fn rotate(&self, radial: f64, north: f64, east: f64) -> Vector {
        [
            radial * self.cos_lambda * self.cos_phi
                - east * self.sin_lambda
                - north * self.sin_phi * self.cos_lambda,
            radial * self.sin_lambda * self.cos_phi + east * self.cos_lambda
                - north * self.sin_phi * self.sin_lambda,
            radial * self.sin_phi + north * self.cos_phi,
        ]
    }
}

/// Explicit epoch/body kernel. Does not select a PPP tide model or silently
/// replace missing IERS ephemerides with the legacy approximate Sun/Moon.
pub fn solid_earth_tide(
    nominal_marker_m: EcefCoord,
    bodies: EarthFixedSunMoon,
    arguments: IersTideArguments,
) -> Result<IersSolidEarthTideCorrection, IersTideError> {
    let site = Site::new(vector(nominal_marker_m))?;
    let mut degree23 = [0.0; 3];
    let mut diurnal = [0.0; 3];
    let mut semidiurnal = [0.0; 3];
    let mut latitude = [0.0; 3];
    let s = site.sin_phi;
    let c = site.cos_phi;
    let sl = site.sin_lambda;
    let cl = site.cos_lambda;
    let cos2phi = c * c - s * s;
    let cos2lambda = cl * cl - sl * sl;
    let sin2lambda = 2.0 * cl * sl;
    let h2 = 0.6078 - 0.0006 * (1.0 - 1.5 * c * c);
    let l2 = 0.0847 + 0.0002 * (1.0 - 1.5 * c * c);
    // Each body has the same Love/Shida expressions and distinct mass/scale.
    for (body, mass_ratio) in [(bodies.sun_m, 332946.0482), (bodies.moon_m, 0.0123000371)] {
        let b = vector(body);
        let radius = checked_norm(b)?;
        let radius2 = radius * radius;
        let factor2 = mass_ratio * RADIUS * (RADIUS / radius).powi(3);
        let factor3 = factor2 * (RADIUS / radius);
        let projection =
            site.position.iter().zip(b).map(|(a, b)| a * b).sum::<f64>() / site.radius / radius;
        let p2 = 3.0 * (h2 / 2.0 - l2) * projection * projection - h2 / 2.0;
        let p3 =
            2.5 * (0.292 - 3.0 * 0.015) * projection.powi(3) + 1.5 * (0.015 - 0.292) * projection;
        let x2 = 3.0 * l2 * projection;
        let x3 = 3.0 * 0.015 / 2.0 * (5.0 * projection * projection - 1.0);
        degree23 = add(
            degree23,
            std::array::from_fn(|i| {
                factor2 * (x2 * b[i] / radius + p2 * site.position[i] / site.radius)
                    + factor3 * (x3 * b[i] / radius + p3 * site.position[i] / site.radius)
            }),
        );

        let cross = b[0] * sl - b[1] * cl;
        let along = b[0] * cl + b[1] * sl;
        diurnal = add(
            diurnal,
            site.rotate(
                -3.0 * (-0.0025) * s * c * factor2 * b[2] * cross / radius2,
                -3.0 * (-0.0007) * cos2phi * factor2 * b[2] * cross / radius2,
                -3.0 * (-0.0007) * s * factor2 * b[2] * along / radius2,
            ),
        );
        let xy_difference = b[0] * b[0] - b[1] * b[1];
        let cross2 = xy_difference * sin2lambda - 2.0 * b[0] * b[1] * cos2lambda;
        let along2 = xy_difference * cos2lambda + 2.0 * b[0] * b[1] * sin2lambda;
        semidiurnal = add(
            semidiurnal,
            site.rotate(
                -3.0 / 4.0 * (-0.0022) * c * c * factor2 * cross2 / radius2,
                3.0 / 2.0 * (-0.0007) * s * c * factor2 * cross2 / radius2,
                -3.0 / 2.0 * (-0.0007) * c * factor2 * along2 / radius2,
            ),
        );
        let north_d = -0.0012 * s * s * factor2 * b[2] * along / radius2;
        let east_d = 0.0012 * s * cos2phi * factor2 * b[2] * cross / radius2;
        let north_sd = -0.0024 / 2.0 * s * c * factor2 * along2 / radius2;
        let east_sd = -0.0024 / 2.0 * s * s * c * factor2 * cross2 / radius2;
        latitude = add(latitude, site.rotate(0.0, 3.0 * north_d, 3.0 * east_d));
        latitude = add(latitude, site.rotate(0.0, 3.0 * north_sd, 3.0 * east_sd));
    }
    let angles = step2_angles(arguments);
    let mut step2_diurnal = [0.0; 3];
    let longitude = site.position[1].atan2(site.position[0]);
    for row in DIURNAL {
        let theta = (angles[5] + (0..5).map(|j| row[j] * angles[j]).sum::<f64>()).to_radians();
        let (st, ct) = (theta + longitude).sin_cos();
        let radial = row[5] * (2.0 * s * c) * st + row[6] * (2.0 * s * c) * ct;
        let north = row[7] * cos2phi * st + row[8] * cos2phi * ct;
        let east = row[7] * s * ct - row[8] * s * st;
        step2_diurnal = add(step2_diurnal, site.rotate(radial, north, east));
    }
    step2_diurnal = step2_diurnal.map(|x| x * 1e-3);
    let mut step2_long_period = [0.0; 3];
    for row in LONG_PERIOD {
        let theta = (0..5).map(|j| row[j] * angles[j]).sum::<f64>().to_radians();
        let (st, ct) = theta.sin_cos();
        let f1 = (3.0 * s * s - 1.0) / 2.0;
        let f2 = 2.0 * c * s;
        let radial = row[5] * f1 * ct + row[7] * f1 * st;
        let north = row[6] * f2 * ct + row[8] * f2 * st;
        step2_long_period = add(step2_long_period, site.rotate(radial, north, 0.0));
    }
    step2_long_period = step2_long_period.map(|x| x * 1e-3);
    let displacement = [
        diurnal,
        semidiurnal,
        latitude,
        step2_diurnal,
        step2_long_period,
    ]
    .into_iter()
    .fold(degree23, add);
    let instantaneous = add(site.position, displacement);
    if angles
        .iter()
        .chain(displacement.iter())
        .chain(instantaneous.iter())
        .any(|x| !x.is_finite())
    {
        return Err(IersTideError::NonFinite);
    }
    Ok(IersSolidEarthTideCorrection {
        arguments,
        bodies,
        nominal_marker_m,
        degree2_degree3_m: coord(degree23),
        step1_diurnal_m: coord(diurnal),
        step1_semidiurnal_m: coord(semidiurnal),
        step1_latitude_m: coord(latitude),
        step2_diurnal_m: coord(step2_diurnal),
        step2_long_period_m: coord(step2_long_period),
        displacement_m: coord(displacement),
        instantaneous_marker_m: coord(instantaneous),
    })
}

fn step2_angles(arguments: IersTideArguments) -> [f64; 6] {
    let t = arguments.julian_centuries_tt;
    let mut s = 218.31664563 + (481267.88194 + (-0.0014663889 + 0.00000185139 * t) * t) * t;
    let tau = arguments.fractional_hour_ut * 15.0
        + 280.4606184
        + (36000.7700536 + (0.00038793 - 0.0000000258 * t) * t) * t
        - s;
    s += (1.396971278 + (0.000308889 + (0.000000021 + 0.000000007 * t) * t) * t) * t;
    let h = 280.46645
        + (36000.7697489 + (0.00030322222 + (0.000000020 - 0.00000000654 * t) * t) * t) * t;
    let p = 83.35324312
        + (4069.01363525 + (-0.01032172222 + (-0.0000124991 + 0.00000005263 * t) * t) * t) * t;
    let zns = 234.95544499
        + (1934.13626197 + (-0.00207561111 + (-0.00000213944 + 0.00000001650 * t) * t) * t) * t;
    let ps = 282.93734098
        + (1.71945766667 + (0.00045688889 + (-0.00000001778 - 0.00000000334 * t) * t) * t) * t;
    // Rust remainder has the same signed convention as native fmod.
    [s, h, p, zns, ps, tau].map(|x| x % 360.0)
}

const DIURNAL: [[f64; 9]; 31] = [
    [-3e0, 0e0, 2e0, 0e0, 0e0, -0.01e0, 0e0, 0e0, 0e0],
    [-3e0, 2e0, 0e0, 0e0, 0e0, -0.01e0, 0e0, 0e0, 0e0],
    [-2e0, 0e0, 1e0, -1e0, 0e0, -0.02e0, 0e0, 0e0, 0e0],
    [-2e0, 0e0, 1e0, 0e0, 0e0, -0.08e0, 0e0, -0.01e0, 0.01e0],
    [-2e0, 2e0, -1e0, 0e0, 0e0, -0.02e0, 0e0, 0e0, 0e0],
    [-1e0, 0e0, 0e0, -1e0, 0e0, -0.10e0, 0e0, 0e0, 0e0],
    [-1e0, 0e0, 0e0, 0e0, 0e0, -0.51e0, 0e0, -0.02e0, 0.03e0],
    [-1e0, 2e0, 0e0, 0e0, 0e0, 0.01e0, 0e0, 0e0, 0e0],
    [0e0, -2e0, 1e0, 0e0, 0e0, 0.01e0, 0e0, 0e0, 0e0],
    [0e0, 0e0, -1e0, 0e0, 0e0, 0.02e0, 0e0, 0e0, 0e0],
    [0e0, 0e0, 1e0, 0e0, 0e0, 0.06e0, 0e0, 0e0, 0e0],
    [0e0, 0e0, 1e0, 1e0, 0e0, 0.01e0, 0e0, 0e0, 0e0],
    [0e0, 2e0, -1e0, 0e0, 0e0, 0.01e0, 0e0, 0e0, 0e0],
    [1e0, -3e0, 0e0, 0e0, 1e0, -0.06e0, 0e0, 0e0, 0e0],
    [1e0, -2e0, 0e0, -1e0, 0e0, 0.01e0, 0e0, 0e0, 0e0],
    [1e0, -2e0, 0e0, 0e0, 0e0, -1.23e0, -0.07e0, 0.06e0, 0.01e0],
    [1e0, -1e0, 0e0, 0e0, -1e0, 0.02e0, 0e0, 0e0, 0e0],
    [1e0, -1e0, 0e0, 0e0, 1e0, 0.04e0, 0e0, 0e0, 0e0],
    [1e0, 0e0, 0e0, -1e0, 0e0, -0.22e0, 0.01e0, 0.01e0, 0e0],
    [1e0, 0e0, 0e0, 0e0, 0e0, 12.00e0, -0.80e0, -0.67e0, -0.03e0],
    [1e0, 0e0, 0e0, 1e0, 0e0, 1.73e0, -0.12e0, -0.10e0, 0e0],
    [1e0, 0e0, 0e0, 2e0, 0e0, -0.04e0, 0e0, 0e0, 0e0],
    [1e0, 1e0, 0e0, 0e0, -1e0, -0.50e0, -0.01e0, 0.03e0, 0e0],
    [1e0, 1e0, 0e0, 0e0, 1e0, 0.01e0, 0e0, 0e0, 0e0],
    [0e0, 1e0, 0e0, 1e0, -1e0, -0.01e0, 0e0, 0e0, 0e0],
    [1e0, 2e0, -2e0, 0e0, 0e0, -0.01e0, 0e0, 0e0, 0e0],
    [1e0, 2e0, 0e0, 0e0, 0e0, -0.11e0, 0.01e0, 0.01e0, 0e0],
    [2e0, -2e0, 1e0, 0e0, 0e0, -0.01e0, 0e0, 0e0, 0e0],
    [2e0, 0e0, -1e0, 0e0, 0e0, -0.02e0, 0e0, 0e0, 0e0],
    [3e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0, 0e0],
    [3e0, 0e0, 0e0, 1e0, 0e0, 0e0, 0e0, 0e0, 0e0],
];

const LONG_PERIOD: [[f64; 9]; 5] = [
    [0e0, 0e0, 0e0, 1e0, 0e0, 0.47e0, 0.23e0, 0.16e0, 0.07e0],
    [0e0, 2e0, 0e0, 0e0, 0e0, -0.20e0, -0.12e0, -0.11e0, -0.05e0],
    [1e0, 0e0, -1e0, 0e0, 0e0, -0.11e0, -0.08e0, -0.09e0, -0.04e0],
    [2e0, 0e0, 0e0, 0e0, 0e0, -0.13e0, -0.11e0, -0.15e0, -0.07e0],
    [2e0, 0e0, 0e0, 1e0, 0e0, -0.05e0, -0.05e0, -0.06e0, -0.03e0],
];
