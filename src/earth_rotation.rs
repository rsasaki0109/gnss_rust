//! GPS-era time scales and IAU Earth rotation with explicit leap coverage/EOP.
//! IERS Conventions: ERA (2003, chapter 5, eq. 14), TIO secular locator s',
//! and polar-motion rotations. Mathematical implementation in Rust; no SOFA
//! source/FFI. Pinned libgnss++/unchanged SOFA provide independent references.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.

use crate::iers2010::{IersTideArguments, IersTideError};
use crate::time::{CalendarDateTime, TimeScale};
use crate::{Error, GnssTime};
use std::fmt;

const GPS_MJD: i32 = 44244;
const DAY: f64 = 86400.0;
const ARCSEC: f64 = std::f64::consts::PI / (180.0 * 3600.0);
pub type Matrix3 = [[f64; 3]; 3];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RotationError {
    NonFinite,
    OutOfRange,
    InvalidLeapTable,
    LeapCoverage,
    InvalidRotation,
    Gnss(Error),
    Tide(IersTideError),
}
impl fmt::Display for RotationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Earth rotation: {self:?}")
    }
}
impl std::error::Error for RotationError {}
impl From<Error> for RotationError {
    fn from(e: Error) -> Self {
        Self::Gnss(e)
    }
}

/// Uniform MJD, retaining integer day and SI seconds rather than a large JD.
/// UTC is a separate type because its days may contain a leap second.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniformMjd {
    day: i32,
    seconds: f64,
}
impl UniformMjd {
    pub fn new(day: i32, seconds: f64) -> Result<Self, RotationError> {
        if !seconds.is_finite() {
            return Err(RotationError::NonFinite);
        }
        let shift = (seconds / DAY).floor();
        if shift < f64::from(i32::MIN) || shift > f64::from(i32::MAX) {
            return Err(RotationError::OutOfRange);
        }
        let mut day = day
            .checked_add(shift as i32)
            .ok_or(RotationError::OutOfRange)?;
        let mut seconds = seconds.rem_euclid(DAY);
        if seconds >= DAY {
            day = day.checked_add(1).ok_or(RotationError::OutOfRange)?;
            seconds = 0.0;
        }
        Ok(Self { day, seconds })
    }
    pub fn day_mjd(self) -> i32 {
        self.day
    }
    pub fn seconds_of_day(self) -> f64 {
        self.seconds
    }
    pub fn mjd(self) -> f64 {
        f64::from(self.day) + self.seconds / DAY
    }
    pub fn julian_date_parts(self) -> [f64; 2] {
        [2400000.5 + f64::from(self.day), self.seconds / DAY]
    }
    pub fn julian_centuries_tt(self) -> f64 {
        (f64::from(self.day) - 51544.5 + self.seconds / DAY) / 36525.0
    }
    pub fn difference_seconds(self, other: Self) -> f64 {
        (i64::from(self.day) - i64::from(other.day)) as f64 * DAY + self.seconds - other.seconds
    }
    fn add(self, seconds: f64) -> Result<Self, RotationError> {
        Self::new(self.day, self.seconds + seconds)
    }
}

/// UTC quasi-MJD: fraction uses the actual day length (86400 or 86401).
/// During a positive leap second seconds_of_day is in [86400,86401).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UtcEpoch {
    day: i32,
    seconds: f64,
    day_length: u32,
    tai_minus_utc_s: i32,
}
impl UtcEpoch {
    pub fn day_mjd(self) -> i32 {
        self.day
    }
    pub fn seconds_of_day(self) -> f64 {
        self.seconds
    }
    pub fn day_length_s(self) -> u32 {
        self.day_length
    }
    pub fn tai_minus_utc_s(self) -> i32 {
        self.tai_minus_utc_s
    }
    pub fn is_leap_second(self) -> bool {
        self.seconds >= DAY
    }
    pub fn quasi_mjd(self) -> f64 {
        f64::from(self.day) + self.seconds / f64::from(self.day_length)
    }
    pub fn julian_date_parts(self) -> [f64; 2] {
        [
            2400000.5 + f64::from(self.day),
            self.seconds / f64::from(self.day_length),
        ]
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LeapSecondTransition {
    pub utc_day_mjd: i32,
    pub tai_minus_utc_s: i32,
}
/// GPS-era positive-leap table. The half-open validity end is mandatory;
/// caller declares that no additional transitions are missing before it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeapSecondTable {
    product_id: String,
    transitions: Vec<LeapSecondTransition>,
    valid_until_utc_mjd: i32,
}
impl LeapSecondTable {
    pub fn new(
        product_id: impl Into<String>,
        transitions: Vec<LeapSecondTransition>,
        valid_until_utc_mjd: i32,
    ) -> Result<Self, RotationError> {
        let product_id = product_id.into();
        if product_id.trim().is_empty()
            || transitions.first()
                != Some(&LeapSecondTransition {
                    utc_day_mjd: GPS_MJD,
                    tai_minus_utc_s: 19,
                })
            || transitions
                .last()
                .is_none_or(|t| valid_until_utc_mjd <= t.utc_day_mjd)
            || transitions.windows(2).any(|w| {
                w[1].utc_day_mjd <= w[0].utc_day_mjd
                    || w[0].tai_minus_utc_s.checked_add(1) != Some(w[1].tai_minus_utc_s)
            })
        {
            return Err(RotationError::InvalidLeapTable);
        }
        Ok(Self {
            product_id,
            transitions,
            valid_until_utc_mjd,
        })
    }
    /// Known transitions through 2017-01-01. The caller supplies the coverage
    /// end explicitly; this method does not verify current announcements.
    pub fn historical(
        product_id: impl Into<String>,
        valid_until_utc_mjd: i32,
    ) -> Result<Self, RotationError> {
        let mut transitions = vec![LeapSecondTransition {
            utc_day_mjd: GPS_MJD,
            tai_minus_utc_s: 19,
        }];
        for (i, &(year, month, day)) in crate::time::UTC_LEAP_TRANSITIONS.iter().enumerate() {
            let label = CalendarDateTime {
                year,
                month,
                day,
                hour: 0,
                minute: 0,
                second: 0.0,
            }
            .to_gnss_time(TimeScale::Gps)?;
            let mjd = GPS_MJD + label.week() * 7 + (label.tow() / DAY) as i32;
            transitions.push(LeapSecondTransition {
                utc_day_mjd: mjd,
                tai_minus_utc_s: 20 + i as i32,
            });
        }
        Self::new(product_id, transitions, valid_until_utc_mjd)
    }
    pub fn product_id(&self) -> &str {
        &self.product_id
    }
    pub fn transitions(&self) -> &[LeapSecondTransition] {
        &self.transitions
    }
    pub fn valid_until_utc_mjd(&self) -> i32 {
        self.valid_until_utc_mjd
    }
    fn utc_from_tai(&self, tai: UniformMjd) -> Result<UtcEpoch, RotationError> {
        let n = self.transitions.partition_point(|t| {
            tai.difference_seconds(UniformMjd {
                day: t.utc_day_mjd,
                seconds: f64::from(t.tai_minus_utc_s),
            }) >= 0.0
        });
        if n == 0 {
            return Err(RotationError::LeapCoverage);
        }
        let offset = self.transitions[n - 1].tai_minus_utc_s;
        let mut uniform = tai.add(-f64::from(offset))?;
        let mut seconds = uniform.seconds;
        if let Some(next) = self.transitions.get(n) {
            let remaining = UniformMjd {
                day: next.utc_day_mjd,
                seconds: f64::from(next.tai_minus_utc_s),
            }
            .difference_seconds(tai);
            if remaining <= 1.0 {
                uniform.day = next.utc_day_mjd - 1;
                seconds = DAY + 1.0 - remaining;
            }
        }
        if uniform.day < GPS_MJD || uniform.day >= self.valid_until_utc_mjd {
            return Err(RotationError::LeapCoverage);
        }
        let day_length = if self
            .transitions
            .iter()
            .any(|t| t.utc_day_mjd == uniform.day + 1)
        {
            86401
        } else {
            86400
        };
        Ok(UtcEpoch {
            day: uniform.day,
            seconds,
            day_length,
            tai_minus_utc_s: offset,
        })
    }
}

/// Explicit values at the evaluated epoch. No implicit zero-EOP policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EarthOrientationParams {
    pub ut1_minus_utc_s: f64,
    pub xp_arcsec: f64,
    pub yp_arcsec: f64,
}
impl EarthOrientationParams {
    pub fn new(
        ut1_minus_utc_s: f64,
        xp_arcsec: f64,
        yp_arcsec: f64,
    ) -> Result<Self, RotationError> {
        let value = Self {
            ut1_minus_utc_s,
            xp_arcsec,
            yp_arcsec,
        };
        value.validate()?;
        Ok(value)
    }
    fn validate(self) -> Result<(), RotationError> {
        if [self.ut1_minus_utc_s, self.xp_arcsec, self.yp_arcsec]
            .iter()
            .any(|x| !x.is_finite())
        {
            Err(RotationError::NonFinite)
        } else {
            Ok(())
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RotationTimeScales {
    pub gps: GnssTime,
    pub utc: UtcEpoch,
    pub tai: UniformMjd,
    pub tt: UniformMjd,
    pub ut1: UniformMjd,
    pub eop: EarthOrientationParams,
    pub leap_product_id: String,
    pub leap_valid_until_utc_mjd: i32,
}
impl RotationTimeScales {
    pub fn from_gps(
        gps: GnssTime,
        leaps: &LeapSecondTable,
        eop: EarthOrientationParams,
    ) -> Result<Self, RotationError> {
        eop.validate()?;
        let day = i64::from(GPS_MJD) + i64::from(gps.week()) * 7 + (gps.tow() / DAY).floor() as i64;
        let day = day.try_into().map_err(|_| RotationError::OutOfRange)?;
        let tai = UniformMjd::new(day, gps.tow() % DAY + 19.0)?;
        let utc = leaps.utc_from_tai(tai)?;
        let tt = tai.add(32.184)?;
        let ut1 = tai.add(eop.ut1_minus_utc_s - f64::from(utc.tai_minus_utc_s))?;
        Ok(Self {
            gps,
            utc,
            tai,
            tt,
            ut1,
            eop,
            leap_product_id: leaps.product_id.clone(),
            leap_valid_until_utc_mjd: leaps.valid_until_utc_mjd,
        })
    }
    /// Match native tide wrapper's UT≈UTC convention, using quasi-UTC hour.
    pub fn tide_arguments_utc_approx(&self) -> Result<IersTideArguments, RotationError> {
        IersTideArguments::new(
            self.tt.julian_centuries_tt(),
            self.utc.seconds / f64::from(self.utc.day_length) * 24.0,
        )
        .map_err(RotationError::Tide)
    }
}

/// ERA from uniform UT1, IERS 2003 ch.5 eq.14 / Capitaine et al. (2000).
pub fn earth_rotation_angle(ut1: UniformMjd) -> f64 {
    let fraction = ut1.seconds / DAY;
    let d = f64::from(ut1.day) - 51544.5 + fraction;
    (std::f64::consts::TAU * (0.7790572732640 + fraction + 0.5 + 0.00273781191135448 * d))
        .rem_euclid(std::f64::consts::TAU)
}
fn multiply(a: Matrix3, b: Matrix3) -> Matrix3 {
    std::array::from_fn(|i| std::array::from_fn(|j| (0..3).map(|k| a[i][k] * b[k][j]).sum()))
}
pub fn transpose(m: Matrix3) -> Matrix3 {
    std::array::from_fn(|i| std::array::from_fn(|j| m[j][i]))
}
fn spin_z(angle: f64) -> Matrix3 {
    let (s, c) = angle.sin_cos();
    [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]]
}
/// Polar motion including the IAU secular TIO locator s'=−47 µas/century.
pub fn polar_motion_matrix(
    tt: UniformMjd,
    eop: EarthOrientationParams,
) -> Result<Matrix3, RotationError> {
    eop.validate()?;
    let xp = eop.xp_arcsec * ARCSEC;
    let yp = eop.yp_arcsec * ARCSEC;
    let sp = -47e-6 * ARCSEC * tt.julian_centuries_tt();
    if [xp, yp, sp].iter().any(|x| !x.is_finite()) {
        return Err(RotationError::NonFinite);
    }
    let (sx, cx) = xp.sin_cos();
    let (sy, cy) = yp.sin_cos();
    let ry = [[cx, 0.0, sx], [0.0, 1.0, 0.0], [-sx, 0.0, cx]];
    let rx = [[1.0, 0.0, 0.0], [0.0, cy, -sy], [0.0, sy, cy]];
    Ok(multiply(rx, multiply(ry, spin_z(sp))))
}
/// Compose with a caller-supplied celestial→intermediate (CIO) matrix.
/// Orthogonality/det are checked; this cannot prove its frame/epoch/model.
/// IAU 2006/2000A precession-nutation is NOT generated by this function.
pub fn compose_celestial_to_terrestrial(
    times: &RotationTimeScales,
    celestial_to_intermediate: Matrix3,
) -> Result<Matrix3, RotationError> {
    validate_rotation(celestial_to_intermediate)?;
    let polar = polar_motion_matrix(times.tt, times.eop)?;
    Ok(multiply(
        polar,
        multiply(
            spin_z(earth_rotation_angle(times.ut1)),
            celestial_to_intermediate,
        ),
    ))
}
fn validate_rotation(m: Matrix3) -> Result<(), RotationError> {
    if m.iter().flatten().any(|x| !x.is_finite()) {
        return Err(RotationError::NonFinite);
    }
    let gram = multiply(m, transpose(m));
    if (0..3).any(|i| (0..3).any(|j| (gram[i][j] - f64::from(i == j)).abs() > 1e-12)) {
        return Err(RotationError::InvalidRotation);
    }
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if (det - 1.0).abs() > 1e-12 {
        return Err(RotationError::InvalidRotation);
    }
    Ok(())
}
