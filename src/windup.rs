//! Wu et al. (1993) phase wind-up, nominal yaw and approximate Sun in ECEF.
//! Port of pinned libgnss++ ppp.cpp and ppp_corrections.cpp helpers.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
//! No eclipse/noon/midnight yaw model, EOP/UT1 or receiver attitude estimation.

use crate::coordinates::ecef_to_geodetic;
use crate::time::day_of_year_gps_scale;
use crate::{EcefCoord, Error, GnssTime};
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindupError {
    NonFinite,
    DegenerateLineOfSight,
    InvalidFrame,
    MissingSun,
    DegenerateAttitude,
    DegenerateDipoles,
    Gnss(Error),
}
impl fmt::Display for WindupError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "phase wind-up: {self:?}")
    }
}
impl std::error::Error for WindupError {}
impl From<Error> for WindupError {
    fn from(error: Error) -> Self {
        Self::Gnss(error)
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DipoleFrame {
    /// Unit, orthogonal ECEF dipole axes; z is x cross y.
    pub x: [f64; 3],
    pub y: [f64; 3],
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindupCycles {
    pub principal_cycles: f64,
    /// Principal angle plus integer nearest the previous accepted cycle count.
    pub cycles: f64,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PppWindupModel {
    /// Native approximate Sun in GPS-labelled time, nominal yaw steering,
    /// fixed north-west-up receiver. Explicit opt-in, no fallback attitude.
    NominalYawApproximateSun,
}
type V = [f64; 3];
fn finite(v: V) -> bool {
    v.iter().all(|x| x.is_finite())
}
fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: V) -> f64 {
    a[0].hypot(a[1]).hypot(a[2])
}
fn scale(a: V, s: f64) -> V {
    a.map(|x| x * s)
}
fn sub(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] - b[i])
}
fn add(a: V, b: V) -> V {
    std::array::from_fn(|i| a[i] + b[i])
}
fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn vector(p: EcefCoord) -> V {
    [p.x, p.y, p.z]
}
fn validate_frame(frame: DipoleFrame) -> Result<(), WindupError> {
    if !finite(frame.x) || !finite(frame.y) {
        return Err(WindupError::NonFinite);
    }
    if (norm(frame.x) - 1.0).abs() > 1e-10
        || (norm(frame.y) - 1.0).abs() > 1e-10
        || dot(frame.x, frame.y).abs() > 1e-10
    {
        return Err(WindupError::InvalidFrame);
    }
    Ok(())
}
/// Explicit ECEF dipole frames, with line of sight pointing satellite to
/// receiver. Uses native acos/orientation and C++ round (ties away from zero).
pub fn phase_windup_from_frames(
    los_satellite_to_receiver: V,
    satellite: DipoleFrame,
    receiver: DipoleFrame,
    previous_cycles: f64,
) -> Result<WindupCycles, WindupError> {
    if !finite(los_satellite_to_receiver) || !previous_cycles.is_finite() {
        return Err(WindupError::NonFinite);
    }
    validate_frame(satellite)?;
    validate_frame(receiver)?;
    let length = norm(los_satellite_to_receiver);
    if !length.is_finite() {
        return Err(WindupError::NonFinite);
    }
    if length < 1.0 {
        return Err(WindupError::DegenerateLineOfSight);
    }
    let e = scale(los_satellite_to_receiver, 1.0 / length);
    let ds = sub(
        sub(satellite.x, scale(e, dot(e, satellite.x))),
        cross(e, satellite.y),
    );
    let dr = add(
        sub(receiver.x, scale(e, dot(e, receiver.x))),
        cross(e, receiver.y),
    );
    let ns = norm(ds);
    let nr = norm(dr);
    if ns < 1e-10 || nr < 1e-10 {
        return Err(WindupError::DegenerateDipoles);
    }
    let cosine = (dot(ds, dr) / (ns * nr)).clamp(-1.0, 1.0);
    let mut principal = cosine.acos() / std::f64::consts::TAU;
    if dot(e, cross(ds, dr)) < 0.0 {
        principal = -principal;
    }
    let cycles = principal + (previous_cycles - principal).round();
    if !cycles.is_finite() {
        return Err(WindupError::NonFinite);
    }
    Ok(WindupCycles {
        principal_cycles: principal,
        cycles,
    })
}
/// Nominal yaw-steering satellite axes and fixed north-west-up receiver axes.
/// Invalid/missing Sun or degenerate geometry is an error, unlike native's
/// fallback frame or returning the old correction without a quality signal.
pub fn nominal_yaw_windup(
    receiver: EcefCoord,
    satellite: EcefCoord,
    sun: EcefCoord,
    previous_cycles: f64,
) -> Result<WindupCycles, WindupError> {
    let frame = nominal_yaw_satellite_frame(satellite, sun)?;
    let geo = ecef_to_geodetic(receiver)?;
    let (sl, cl) = geo.latitude.sin_cos();
    let (so, co) = geo.longitude.sin_cos();
    phase_windup_from_frames(
        sub(vector(receiver), vector(satellite)),
        frame,
        DipoleFrame {
            x: [-sl * co, -sl * so, cl],
            y: [so, -co, 0.0],
        },
        previous_cycles,
    )
}
/// Nominal yaw axes shared by wind-up and satellite antenna PCO. ANTEX body
/// z = x cross y points toward Earth; no eclipse/yaw fallback is provided.
pub fn nominal_yaw_satellite_frame(
    satellite: EcefCoord,
    sun: EcefCoord,
) -> Result<DipoleFrame, WindupError> {
    let rs = vector(satellite);
    let sun = vector(sun);
    if !finite(rs) || !finite(sun) {
        return Err(WindupError::NonFinite);
    }
    if norm(sun) <= 1.0 {
        return Err(WindupError::MissingSun);
    }
    if norm(rs) < 1.0 {
        return Err(WindupError::DegenerateAttitude);
    }
    let z = scale(rs, -1.0 / norm(rs));
    let direction = sub(sun, rs);
    let length = norm(direction);
    if length < 1e-3 {
        return Err(WindupError::DegenerateAttitude);
    }
    let y = cross(z, scale(direction, 1.0 / length));
    let length = norm(y);
    if length < 1e-10 {
        return Err(WindupError::DegenerateAttitude);
    }
    let y = scale(y, 1.0 / length);
    let x = cross(y, z);
    Ok(DipoleFrame { x, y })
}
/// Native low-order solar ephemeris and GMST, using fractional GPS labels
/// directly (no GPS→UTC, UT1, precession/nutation or polar motion correction).
pub fn approximate_sun_position_ecef(time: GnssTime) -> Result<EcefCoord, WindupError> {
    day_of_year_gps_scale(time)?; // reject unsupported calendar range
    let jd = 2444244.5 + f64::from(time.week()) * 7.0 + time.tow() / 86400.0;
    let d = jd - 2451545.0;
    let t = d / 36525.0;
    let rad = std::f64::consts::PI / 180.0;
    let wrap = |x: f64| x.rem_euclid(std::f64::consts::TAU);
    let longitude = wrap((280.460 + 0.9856474 * d) * rad);
    let anomaly = wrap((357.528 + 0.9856003 * d) * rad);
    let ecliptic = wrap(longitude + (1.915 * anomaly.sin() + 0.020 * (2.0 * anomaly).sin()) * rad);
    let obliquity = (23.439 - 0.0000004 * d) * rad;
    let radius =
        (1.00014 - 0.01671 * anomaly.cos() - 0.00014 * (2.0 * anomaly).cos()) * 149597870700.0;
    let eci = [
        radius * ecliptic.cos(),
        radius * obliquity.cos() * ecliptic.sin(),
        radius * obliquity.sin() * ecliptic.sin(),
    ];
    let gmst = wrap(
        (280.46061837 + 360.98564736629 * d + 0.000387933 * t * t - t * t * t / 38710000.0) * rad,
    );
    let (s, c) = gmst.sin_cos();
    let result = EcefCoord::new(c * eci[0] + s * eci[1], -s * eci[0] + c * eci[1], eci[2]);
    if !finite(vector(result)) {
        return Err(WindupError::NonFinite);
    }
    Ok(result)
}
