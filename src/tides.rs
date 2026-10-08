//! Pinned libgnss++ legacy degree-2 Love-number solid Earth tide.
//! Explicit approximation: not the native default IERS 2010 Step-1+2 model.
//! GPS-labelled low-order Sun/Moon and GMST; no UTC/UT1/EOP or permanent-tide
//! convention conversion, ocean loading, pole tide or atmospheric loading.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.

use crate::coordinates::ecef_to_geodetic;
use crate::time::day_of_year_gps_scale;
use crate::windup::{WindupError, approximate_sun_position_ecef};
use crate::{EcefCoord, Error, GnssTime, constants as c};
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TideError {
    NonFinite,
    InvalidStation,
    InvalidBody,
    InvalidBodyGravity,
    Gnss(Error),
    Ephemeris(WindupError),
}
impl fmt::Display for TideError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "solid Earth tide: {self:?}")
    }
}
impl std::error::Error for TideError {}
impl From<Error> for TideError {
    fn from(e: Error) -> Self {
        Self::Gnss(e)
    }
}

/// Defines the XYZ states/report, independently from the SP3 satellite COM/APC.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReceiverPositionReference {
    /// Geometry uses the reported marker position directly; tide model None.
    InstantaneousMarker,
    /// Reference marker before adding the explicitly selected tide model.
    /// This does not assert a standard tide-free/zero-tide terrestrial datum.
    NominalMarker,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolidEarthTideModel {
    /// Native alternative when use_iers_solid_tide=false. Degree-2 fixed
    /// h2=.6078/l2=.0847; geodetic up and native approximate Sun/Moon in ECEF.
    LegacyLoveApproximateSunMoon,
}
#[derive(Debug, Clone, PartialEq)]
pub struct SolidEarthTideCorrection {
    pub model: SolidEarthTideModel,
    pub time: GnssTime,
    pub nominal_marker_m: EcefCoord,
    pub sun_position_m: EcefCoord,
    pub moon_position_m: EcefCoord,
    pub sun_displacement_m: EcefCoord,
    pub moon_displacement_m: EcefCoord,
    /// ADD to nominal marker, once, to obtain instantaneous geometry marker.
    pub displacement_m: EcefCoord,
    pub instantaneous_marker_m: EcefCoord,
}
fn vector(v: EcefCoord) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn finite(v: EcefCoord) -> bool {
    vector(v).iter().all(|x| x.is_finite())
}
fn length(v: EcefCoord) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
pub(crate) fn validate_station(position: EcefCoord) -> Result<(), TideError> {
    if !finite(position) || !length(position).is_finite() {
        return Err(TideError::NonFinite);
    }
    if length(position) < c::WGS84_A * 0.5 {
        return Err(TideError::InvalidStation);
    }
    Ok(())
}

/// Native low-order Moon in ECEF at fractional GPS-labelled epoch. The frame
/// and time convention match the existing approximate Sun, not SOFA/UT1.
pub fn approximate_moon_position_ecef(time: GnssTime) -> Result<EcefCoord, TideError> {
    day_of_year_gps_scale(time)?;
    let jd = 2444244.5 + f64::from(time.week()) * 7.0 + time.tow() / 86400.0;
    let d = jd - 2451545.0;
    let t = d / 36525.0;
    let rad = std::f64::consts::PI / 180.0;
    let wrap = |x: f64| x.rem_euclid(std::f64::consts::TAU);
    let longitude = wrap((218.316 + 13.176396 * d) * rad);
    let moon_anomaly = wrap((134.963 + 13.064993 * d) * rad);
    let sun_anomaly = wrap((357.529 + 0.98560028 * d) * rad);
    let elongation = wrap((297.850 + 12.190749 * d) * rad);
    let latitude_argument = wrap((93.272 + 13.229350 * d) * rad);
    let lon = wrap(
        longitude
            + (6.289 * moon_anomaly.sin()
                + 1.274 * (2.0 * elongation - moon_anomaly).sin()
                + 0.658 * (2.0 * elongation).sin()
                + 0.214 * (2.0 * moon_anomaly).sin()
                - 0.186 * sun_anomaly.sin())
                * rad,
    );
    let lat = (5.128 * latitude_argument.sin()
        + 0.280 * (moon_anomaly + latitude_argument).sin()
        + 0.277 * (moon_anomaly - latitude_argument).sin()
        + 0.173 * (2.0 * elongation - latitude_argument).sin())
        * rad;
    let radius = (385001.0
        - 20905.0 * moon_anomaly.cos()
        - 3699.0 * (2.0 * elongation - moon_anomaly).cos()
        - 2956.0 * (2.0 * elongation).cos()
        - 570.0 * (2.0 * moon_anomaly).cos())
        * 1000.0;
    let obliquity = (23.439 - 0.0000004 * d) * rad;
    let eci = [
        radius * lat.cos() * lon.cos(),
        radius * (obliquity.cos() * lat.cos() * lon.sin() - obliquity.sin() * lat.sin()),
        radius * (obliquity.sin() * lat.cos() * lon.sin() + obliquity.cos() * lat.sin()),
    ];
    let gmst = wrap(
        (280.46061837 + 360.98564736629 * d + 0.000387933 * t * t - t * t * t / 38710000.0) * rad,
    );
    let result = EcefCoord::new(
        gmst.cos() * eci[0] + gmst.sin() * eci[1],
        -gmst.sin() * eci[0] + gmst.cos() * eci[1],
        eci[2],
    );
    if !finite(result) {
        return Err(TideError::NonFinite);
    }
    Ok(result)
}

/// Explicit-body component. Earth-centred ECEF body vector and GM in m³/s²;
/// returns ECEF displacement in m. Native invalid-body zero fallback is an
/// error here. Uses geodetic up and fixed WGS84 radius, as the legacy helper.
pub fn legacy_body_tide_displacement(
    receiver: EcefCoord,
    body: EcefCoord,
    body_gm_m3_s2: f64,
) -> Result<EcefCoord, TideError> {
    validate_station(receiver)?;
    if !finite(body) || !body_gm_m3_s2.is_finite() {
        return Err(TideError::NonFinite);
    }
    if body_gm_m3_s2 <= 0.0 {
        return Err(TideError::InvalidBodyGravity);
    }
    let distance = length(body);
    if !distance.is_finite() {
        return Err(TideError::NonFinite);
    }
    if distance <= 1.0 {
        return Err(TideError::InvalidBody);
    }
    let geo = ecef_to_geodetic(receiver)?;
    let mut up = [
        geo.latitude.cos() * geo.longitude.cos(),
        geo.latitude.cos() * geo.longitude.sin(),
        geo.latitude.sin(),
    ];
    let norm = up[0].hypot(up[1]).hypot(up[2]);
    up = up.map(|x| x / norm);
    let direction = vector(body).map(|x| x / distance);
    let projection = up
        .iter()
        .zip(direction)
        .map(|(u, b)| u * b)
        .sum::<f64>()
        .clamp(-1.0, 1.0);
    let scale = body_gm_m3_s2 / 3.986004418e14 * (c::WGS84_A / distance).powi(3) * c::WGS84_A;
    let displacement: [f64; 3] = std::array::from_fn(|i| {
        scale
            * (0.6078 * (1.5 * projection * projection - 0.5) * up[i]
                + 3.0 * 0.0847 * projection * (direction[i] - projection * up[i]))
    });
    let result = EcefCoord::new(displacement[0], displacement[1], displacement[2]);
    if !finite(result) {
        return Err(TideError::NonFinite);
    }
    Ok(result)
}
pub fn solid_earth_tide(
    time: GnssTime,
    nominal_marker_m: EcefCoord,
    model: SolidEarthTideModel,
) -> Result<SolidEarthTideCorrection, TideError> {
    validate_station(nominal_marker_m)?;
    let sun_position_m = approximate_sun_position_ecef(time).map_err(TideError::Ephemeris)?;
    let moon_position_m = approximate_moon_position_ecef(time)?;
    let sun_displacement_m =
        legacy_body_tide_displacement(nominal_marker_m, sun_position_m, 1.32712440018e20)?;
    let moon_displacement_m =
        legacy_body_tide_displacement(nominal_marker_m, moon_position_m, 4.902801e12)?;
    let displacement_m = EcefCoord::new(
        sun_displacement_m.x + moon_displacement_m.x,
        sun_displacement_m.y + moon_displacement_m.y,
        sun_displacement_m.z + moon_displacement_m.z,
    );
    let instantaneous_marker_m = EcefCoord::new(
        nominal_marker_m.x + displacement_m.x,
        nominal_marker_m.y + displacement_m.y,
        nominal_marker_m.z + displacement_m.z,
    );
    if !finite(instantaneous_marker_m) {
        return Err(TideError::NonFinite);
    }
    Ok(SolidEarthTideCorrection {
        model,
        time,
        nominal_marker_m,
        sun_position_m,
        moon_position_m,
        sun_displacement_m,
        moon_displacement_m,
        displacement_m,
        instantaneous_marker_m,
    })
}
