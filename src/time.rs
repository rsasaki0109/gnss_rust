//! GPS-scale time without UTC leap-second conversion.

use crate::{
    Error,
    constants::{GPS_EPOCH_UNIX_SECONDS, SECONDS_PER_WEEK},
};
use std::cmp::Ordering;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

// Shared historical GPST-UTC transitions; existing calendar policy unchanged.
pub(crate) const UTC_LEAP_TRANSITIONS: [(i32, u8, u8); 18] = [
    (1981, 7, 1),
    (1982, 7, 1),
    (1983, 7, 1),
    (1985, 7, 1),
    (1988, 1, 1),
    (1990, 1, 1),
    (1991, 1, 1),
    (1992, 7, 1),
    (1993, 7, 1),
    (1994, 7, 1),
    (1996, 1, 1),
    (1997, 7, 1),
    (1999, 1, 1),
    (2006, 1, 1),
    (2009, 1, 1),
    (2012, 7, 1),
    (2015, 7, 1),
    (2017, 1, 1),
];

/// Full signed GPS week and normalized seconds of week (`0 <= tow < 604800`).
///
/// Equality is exact. Use [`Self::difference_seconds`] for a tolerance check.
/// Week rollover is not guessed from a truncated receiver week number.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct GnssTime {
    week: i32,
    tow: f64,
}

impl GnssTime {
    /// Construct and normalize a GPS time, including negative seconds of week.
    pub fn new(week: i32, tow: f64) -> Result<Self, Error> {
        if !tow.is_finite() {
            return Err(Error::NonFinite);
        }
        let shift = (tow / SECONDS_PER_WEEK).floor();
        if shift < i64::MIN as f64 || shift >= i64::MAX as f64 {
            return Err(Error::OutOfRange);
        }
        let mut week = i64::from(week)
            .checked_add(shift as i64)
            .ok_or(Error::OutOfRange)?;
        let mut tow = tow.rem_euclid(SECONDS_PER_WEEK);
        // rem_euclid can round a tiny negative remainder to the divisor.
        if tow >= SECONDS_PER_WEEK {
            week = week.checked_add(1).ok_or(Error::OutOfRange)?;
            tow = 0.0;
        }
        Ok(Self {
            week: week.try_into().map_err(|_| Error::OutOfRange)?,
            tow: if tow == 0.0 { 0.0 } else { tow },
        })
    }

    pub fn week(self) -> i32 {
        self.week
    }
    pub fn tow(self) -> f64 {
        self.tow
    }

    /// Add signed seconds, rejecting invalid inputs and week overflow.
    pub fn checked_add_seconds(self, seconds: f64) -> Result<Self, Error> {
        if !seconds.is_finite() {
            return Err(Error::NonFinite);
        }
        let tow = self.tow + seconds;
        if !tow.is_finite() {
            return Err(Error::OutOfRange);
        }
        Self::new(self.week, tow)
    }

    /// `self - other` in seconds. Small fractions are preserved within a week.
    pub fn difference_seconds(self, other: Self) -> f64 {
        (i64::from(self.week) - i64::from(other.week)) as f64 * SECONDS_PER_WEEK
            + (self.tow - other.tow)
    }

    /// Encode the GPS clock label as a SystemTime, retaining fractional seconds.
    ///
    /// This matches the upstream epoch-offset convention, NOT a conversion to
    /// UTC. A real UTC timestamp requires a separate GPS-UTC leap-second model.
    pub fn to_system_time_gps_scale(self) -> Result<SystemTime, Error> {
        let seconds = i64::from(self.week) * 604_800 + GPS_EPOCH_UNIX_SECONDS as i64;
        let base = if seconds >= 0 {
            UNIX_EPOCH.checked_add(Duration::from_secs(seconds as u64))
        } else {
            UNIX_EPOCH.checked_sub(Duration::from_secs(seconds.unsigned_abs()))
        }
        .ok_or(Error::OutOfRange)?;
        base.checked_add(Duration::try_from_secs_f64(self.tow).map_err(|_| Error::OutOfRange)?)
            .ok_or(Error::OutOfRange)
    }

    /// Decode a GPS-scale SystemTime label; no UTC leap-second adjustment occurs.
    pub fn from_system_time_gps_scale(time: SystemTime) -> Result<Self, Error> {
        let epoch = UNIX_EPOCH
            .checked_add(Duration::from_secs(GPS_EPOCH_UNIX_SECONDS))
            .ok_or(Error::OutOfRange)?;
        let (duration, sign) = match time.duration_since(epoch) {
            Ok(d) => (d, 1_i64),
            Err(e) => (e.duration(), -1_i64),
        };
        let week =
            i64::try_from(duration.as_secs() / 604_800).map_err(|_| Error::OutOfRange)? * sign;
        let tow = ((duration.as_secs() % 604_800) as f64
            + f64::from(duration.subsec_nanos()) * 1e-9)
            * sign as f64;
        Self::new(week.try_into().map_err(|_| Error::OutOfRange)?, tow)
    }
}

impl PartialOrd for GnssTime {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(
            self.week
                .cmp(&other.week)
                .then_with(|| self.tow.total_cmp(&other.tow)),
        )
    }
}

/// Calendar-label time scales admitted by the RINEX 3 reader.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TimeScale {
    #[default]
    Gps,
    /// Galileo/QZSS calendar labels are aligned with GPS time.
    Galileo,
    Qzss,
    BeiDou,
    Utc,
}

/// Ordinal day (1..=366) of the GPS-scale calendar label, without a UTC leap
/// adjustment. Matches the upstream GPS-label SystemTime mapping used by NMF.
pub fn day_of_year_gps_scale(time: GnssTime) -> Result<u16, Error> {
    fn year_start(year: i64) -> i64 {
        let y = year - 1;
        365 * y + y / 4 - y / 100 + y / 400
    }
    let day =
        year_start(1980) + 5 + i64::from(time.week()) * 7 + (time.tow() / 86_400.0).floor() as i64;
    if day < year_start(1) || day >= year_start(10_000) {
        return Err(Error::OutOfRange);
    }
    let (mut low, mut high) = (1, 10_000);
    while low + 1 < high {
        let middle = (low + high) / 2;
        if year_start(middle) <= day {
            low = middle;
        } else {
            high = middle;
        }
    }
    Ok((day - year_start(low) + 1) as u16)
}

/// Gregorian calendar label. Leap-second instants (`second == 60`) are
/// explicitly unsupported; ordinary UTC labels use the historical leap table.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CalendarDateTime {
    pub year: i32,
    pub month: u8,
    pub day: u8,
    pub hour: u8,
    pub minute: u8,
    pub second: f64,
}

impl CalendarDateTime {
    pub fn to_gnss_time(self, scale: TimeScale) -> Result<GnssTime, Error> {
        let leap = self.year % 4 == 0 && (self.year % 100 != 0 || self.year % 400 == 0);
        let months = [
            31,
            if leap { 29 } else { 28 },
            31,
            30,
            31,
            30,
            31,
            31,
            30,
            31,
            30,
            31,
        ];
        if !(1..=9999).contains(&self.year)
            || !(1..=12).contains(&self.month)
            || self.day == 0
            || self.day > months[usize::from(self.month - 1)]
            || self.hour > 23
            || self.minute > 59
            || !self.second.is_finite()
            || !(0.0..60.0).contains(&self.second)
        {
            return Err(Error::InvalidCalendar);
        }
        fn days_before_year(year: i32) -> i32 {
            let y = year - 1;
            365 * y + y / 4 - y / 100 + y / 400
        }
        let ordinal = days_before_year(self.year)
            + months[..usize::from(self.month - 1)]
                .iter()
                .map(|&x| i32::from(x))
                .sum::<i32>()
            + i32::from(self.day)
            - 1;
        // 1980-01-06: five days after the start of the year.
        let days = ordinal - (days_before_year(1980) + 5);
        let adjustment = match scale {
            TimeScale::BeiDou => 14.0,
            TimeScale::Utc => {
                // GPST-UTC transitions through 2017-01-01. Future UTC labels
                // use the last announced offset, not a speculative leap second.

                if (self.year, self.month, self.day) < (1980, 1, 6) {
                    return Err(Error::InvalidCalendar);
                }
                UTC_LEAP_TRANSITIONS
                    .iter()
                    .filter(|&&date| (self.year, self.month, self.day) >= date)
                    .count() as f64
            }
            _ => 0.0,
        };
        GnssTime::new(
            days.div_euclid(7),
            f64::from(days.rem_euclid(7)) * 86_400.0
                + f64::from(self.hour) * 3600.0
                + f64::from(self.minute) * 60.0
                + self.second
                + adjustment,
        )
    }
}
