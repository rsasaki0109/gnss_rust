//! Satellite identities, signal frequencies and basic solution data.

use crate::{Error, GnssTime, constants as c};
use std::{fmt, str::FromStr};

/// Values retain the upstream constellation identifiers and ordering.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum GnssSystem {
    Unknown = 0x00,
    Gps = 0x01,
    Glonass = 0x02,
    Galileo = 0x04,
    BeiDou = 0x08,
    Qzss = 0x10,
    Sbas = 0x20,
    NavIC = 0x40,
}

impl GnssSystem {
    pub fn rinex_prefix(self) -> char {
        match self {
            Self::Unknown => 'U',
            Self::Gps => 'G',
            Self::Glonass => 'R',
            Self::Galileo => 'E',
            Self::BeiDou => 'C',
            Self::Qzss => 'J',
            Self::Sbas => 'S',
            Self::NavIC => 'I',
        }
    }
}

/// Known constellation and a nonzero PRN, ordered by constellation then PRN.
///
/// PRNs retain the upstream u8 range. Constellation-specific slot allocation
/// and QZSS/SBAS RINEX-number remapping are the responsibility of file readers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SatelliteId {
    system: GnssSystem,
    prn: u8,
}

impl SatelliteId {
    pub fn new(system: GnssSystem, prn: u8) -> Result<Self, Error> {
        if system == GnssSystem::Unknown || prn == 0 {
            return Err(Error::InvalidSatelliteId);
        }
        Ok(Self { system, prn })
    }
    pub fn system(self) -> GnssSystem {
        self.system
    }
    pub fn prn(self) -> u8 {
        self.prn
    }
}

impl fmt::Display for SatelliteId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{:02}", self.system.rinex_prefix(), self.prn)
    }
}

impl FromStr for SatelliteId {
    type Err = Error;
    fn from_str(value: &str) -> Result<Self, Error> {
        let bytes = value.as_bytes();
        if !(3..=4).contains(&bytes.len()) || !bytes[1..].iter().all(u8::is_ascii_digit) {
            return Err(Error::InvalidSatelliteId);
        }
        let system = match bytes[0] {
            b'G' => GnssSystem::Gps,
            b'R' => GnssSystem::Glonass,
            b'E' => GnssSystem::Galileo,
            b'C' => GnssSystem::BeiDou,
            b'J' => GnssSystem::Qzss,
            b'S' => GnssSystem::Sbas,
            b'I' => GnssSystem::NavIC,
            _ => return Err(Error::InvalidSatelliteId),
        };
        Self::new(
            system,
            value[1..].parse().map_err(|_| Error::InvalidSatelliteId)?,
        )
    }
}

/// Signals currently declared by the upstream core types (not all RINEX codes).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum SignalType {
    GpsL1Ca,
    GpsL1P,
    GpsL2P,
    GpsL2C,
    GpsL5,
    GloL1Ca,
    GloL1P,
    GloL2Ca,
    GloL2P,
    GalE1,
    GalE5A,
    GalE5B,
    GalE6,
    BdsB1I,
    BdsB2I,
    BdsB3I,
    BdsB1C,
    BdsB2A,
    QzsL1Ca,
    QzsL2C,
    QzsL5,
}

impl SignalType {
    pub fn system(self) -> GnssSystem {
        match self {
            Self::GpsL1Ca | Self::GpsL1P | Self::GpsL2P | Self::GpsL2C | Self::GpsL5 => {
                GnssSystem::Gps
            }
            Self::GloL1Ca | Self::GloL1P | Self::GloL2Ca | Self::GloL2P => GnssSystem::Glonass,
            Self::GalE1 | Self::GalE5A | Self::GalE5B | Self::GalE6 => GnssSystem::Galileo,
            Self::BdsB1I | Self::BdsB2I | Self::BdsB3I | Self::BdsB1C | Self::BdsB2A => {
                GnssSystem::BeiDou
            }
            Self::QzsL1Ca | Self::QzsL2C | Self::QzsL5 => GnssSystem::Qzss,
        }
    }

    /// Carrier frequency in hertz. FDMA GLONASS L1/L2 requires a channel -7..=6.
    /// The channel argument is ignored for fixed-frequency signals.
    pub fn frequency_hz(self, glonass_channel: Option<i8>) -> Result<f64, Error> {
        let frequency = match self {
            Self::GpsL1Ca | Self::GpsL1P | Self::QzsL1Ca => c::GPS_L1_FREQ,
            Self::GpsL2P | Self::GpsL2C | Self::QzsL2C => c::GPS_L2_FREQ,
            Self::GpsL5 | Self::QzsL5 => c::GPS_L5_FREQ,
            Self::GalE1 => c::GAL_E1_FREQ,
            Self::GalE5A => c::GAL_E5A_FREQ,
            Self::GalE5B => c::GAL_E5B_FREQ,
            Self::GalE6 => c::GAL_E6_FREQ,
            Self::BdsB1I => c::BDS_B1I_FREQ,
            Self::BdsB2I => c::BDS_B2I_FREQ,
            Self::BdsB3I => c::BDS_B3I_FREQ,
            Self::BdsB1C => c::BDS_B1C_FREQ,
            Self::BdsB2A => c::BDS_B2A_FREQ,
            Self::GloL1Ca | Self::GloL1P | Self::GloL2Ca | Self::GloL2P => {
                let channel = glonass_channel.ok_or(Error::MissingGlonassChannel)?;
                if !(-7..=6).contains(&channel) {
                    return Err(Error::InvalidGlonassChannel);
                }
                match self {
                    Self::GloL1Ca | Self::GloL1P => {
                        c::GLO_L1_BASE_FREQ + f64::from(channel) * c::GLO_L1_STEP_FREQ
                    }
                    _ => c::GLO_L2_BASE_FREQ + f64::from(channel) * c::GLO_L2_STEP_FREQ,
                }
            }
        };
        Ok(frequency)
    }

    pub fn wavelength_m(self, glonass_channel: Option<i8>) -> Result<f64, Error> {
        Ok(c::SPEED_OF_LIGHT / self.frequency_hz(glonass_channel)?)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PositioningMode {
    Spp,
    Dgps,
    RtkFloat,
    RtkFixed,
    Ppp,
    PppAr,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SolutionStatus {
    #[default]
    None,
    Spp,
    Dgps,
    Float,
    Fixed,
    PppFloat,
    PppFixed,
    Propagated,
}

/// Latitude/longitude in radians, ellipsoidal height in metres.
/// Public fields permit construction from data; conversion functions validate them.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct GeodeticCoord {
    pub latitude: f64,
    pub longitude: f64,
    pub height: f64,
}

impl GeodeticCoord {
    pub fn new(latitude: f64, longitude: f64, height: f64) -> Result<Self, Error> {
        let result = Self {
            latitude,
            longitude,
            height,
        };
        result.validate()?;
        Ok(result)
    }
    pub(crate) fn validate(self) -> Result<(), Error> {
        if ![self.latitude, self.longitude, self.height]
            .iter()
            .all(|x| x.is_finite())
        {
            return Err(Error::NonFinite);
        }
        if self.latitude.abs() > std::f64::consts::FRAC_PI_2 {
            return Err(Error::InvalidLatitude);
        }
        Ok(())
    }
}

/// Earth-centred Earth-fixed coordinates in metres.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct EcefCoord {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
impl EcefCoord {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }
    pub fn to_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}
impl From<[f64; 3]> for EcefCoord {
    fn from([x, y, z]: [f64; 3]) -> Self {
        Self { x, y, z }
    }
}

/// Local east, north and up components in metres.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct EnuCoord {
    pub east: f64,
    pub north: f64,
    pub up: f64,
}
impl EnuCoord {
    pub const fn new(east: f64, north: f64, up: f64) -> Self {
        Self { east, north, up }
    }
    pub fn is_finite(self) -> bool {
        self.east.is_finite() && self.north.is_finite() && self.up.is_finite()
    }
}

/// Initial solution contract. Solver-specific diagnostics are not ported yet.
/// Absent position, velocity and covariance are explicit `None` values.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PositionSolution {
    pub time: GnssTime,
    pub status: SolutionStatus,
    pub position_ecef: Option<EcefCoord>,
    pub position_covariance: Option<[[f64; 3]; 3]>,
    pub velocity_ecef: Option<[f64; 3]>,
    pub velocity_covariance: Option<[[f64; 3]; 3]>,
    pub num_satellites: usize,
    /// Receiver clock bias in seconds.
    pub receiver_clock_bias: f64,
}
impl PositionSolution {
    /// A propagated solution may have no current satellites, but must have a
    /// finite position. GNSS solutions additionally need at least four satellites.
    pub fn is_valid(&self) -> bool {
        self.position_ecef.is_some_and(EcefCoord::is_finite)
            && self.status != SolutionStatus::None
            && (self.status == SolutionStatus::Propagated || self.num_satellites >= 4)
    }
    /// Status classification only; use `is_valid` to check usability as well.
    pub fn is_fixed(&self) -> bool {
        matches!(
            self.status,
            SolutionStatus::Fixed | SolutionStatus::PppFixed
        )
    }
}
