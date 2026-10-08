//! WGS-84 coordinate transformations and Sagnac-corrected geometric range.
//! ENU rotations operate on difference vectors, not absolute positions.

use crate::{EcefCoord, EnuCoord, Error, GeodeticCoord, constants as c};

pub fn geodetic_to_ecef(position: GeodeticCoord) -> Result<EcefCoord, Error> {
    position.validate()?;
    let (sin_lat, cos_lat) = position.latitude.sin_cos();
    let (sin_lon, cos_lon) = position.longitude.sin_cos();
    let n = c::WGS84_A / (1.0 - c::WGS84_E2 * sin_lat * sin_lat).sqrt();
    let result = EcefCoord::new(
        (n + position.height) * cos_lat * cos_lon,
        (n + position.height) * cos_lat * sin_lon,
        (n * (1.0 - c::WGS84_E2) + position.height) * sin_lat,
    );
    if !result.is_finite() {
        return Err(Error::NonFinite);
    }
    Ok(result)
}

/// Inverse WGS-84 transformation. The Earth centre has no unique geodetic
/// coordinate. Nonconvergence for other interior points is reported explicitly.
pub fn ecef_to_geodetic(position: EcefCoord) -> Result<GeodeticCoord, Error> {
    if !position.is_finite() {
        return Err(Error::NonFinite);
    }
    let p = position.x.hypot(position.y);
    if !p.is_finite() || !p.hypot(position.z).is_finite() {
        return Err(Error::NonFinite);
    }
    if p == 0.0 {
        if position.z == 0.0 {
            return Err(Error::UndefinedGeodeticOrigin);
        }
        return GeodeticCoord::new(
            std::f64::consts::FRAC_PI_2.copysign(position.z),
            0.0,
            position.z.abs() - c::WGS84_B,
        );
    }
    let longitude = position.y.atan2(position.x);
    let mut z = position.z;
    for _ in 0..32 {
        let sin_lat = z / p.hypot(z);
        let n = c::WGS84_A / (1.0 - c::WGS84_E2 * sin_lat * sin_lat).sqrt();
        let next_z = position.z + n * c::WGS84_E2 * sin_lat;
        if (next_z - z).abs() <= 1e-8 {
            let latitude = next_z.atan2(p);
            let (sin_lat, cos_lat) = latitude.sin_cos();
            let height = p * cos_lat + position.z * sin_lat
                - c::WGS84_A * (1.0 - c::WGS84_E2 * sin_lat * sin_lat).sqrt();
            return GeodeticCoord::new(latitude, longitude, height);
        }
        z = next_z;
    }
    Err(Error::NonConvergent)
}

/// Rotate the ECEF difference `target - origin` into the origin's ENU frame.
pub fn ecef_to_enu(difference: EcefCoord, origin: GeodeticCoord) -> Result<EnuCoord, Error> {
    origin.validate()?;
    if !difference.is_finite() {
        return Err(Error::NonFinite);
    }
    let (sin_lat, cos_lat) = origin.latitude.sin_cos();
    let (sin_lon, cos_lon) = origin.longitude.sin_cos();
    let result = EnuCoord::new(
        -sin_lon * difference.x + cos_lon * difference.y,
        -sin_lat * cos_lon * difference.x - sin_lat * sin_lon * difference.y
            + cos_lat * difference.z,
        cos_lat * cos_lon * difference.x
            + cos_lat * sin_lon * difference.y
            + sin_lat * difference.z,
    );
    if !result.is_finite() {
        return Err(Error::NonFinite);
    }
    Ok(result)
}

/// Rotate local ENU components into an ECEF difference vector.
pub fn enu_to_ecef(enu: EnuCoord, origin: GeodeticCoord) -> Result<EcefCoord, Error> {
    origin.validate()?;
    if !enu.is_finite() {
        return Err(Error::NonFinite);
    }
    let (sin_lat, cos_lat) = origin.latitude.sin_cos();
    let (sin_lon, cos_lon) = origin.longitude.sin_cos();
    let result = EcefCoord::new(
        -sin_lon * enu.east - sin_lat * cos_lon * enu.north + cos_lat * cos_lon * enu.up,
        cos_lon * enu.east - sin_lat * sin_lon * enu.north + cos_lat * sin_lon * enu.up,
        cos_lat * enu.north + sin_lat * enu.up,
    );
    if !result.is_finite() {
        return Err(Error::NonFinite);
    }
    Ok(result)
}

/// Euclidean range plus the upstream first-order Earth-rotation (Sagnac)
/// correction. Inputs are satellite/receiver ECEF positions in metres.
pub fn geometric_distance(satellite: EcefCoord, receiver: EcefCoord) -> Result<f64, Error> {
    if !satellite.is_finite() || !receiver.is_finite() {
        return Err(Error::NonFinite);
    }
    let distance = (satellite.x - receiver.x)
        .hypot(satellite.y - receiver.y)
        .hypot(satellite.z - receiver.z);
    let correction =
        c::OMEGA_E / c::SPEED_OF_LIGHT * (satellite.x * receiver.y - satellite.y * receiver.x);
    let result = distance + correction;
    if !result.is_finite() {
        return Err(Error::NonFinite);
    }
    Ok(result)
}
