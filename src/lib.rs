//! Core GNSS primitives ported from libgnss++.
//!
//! Angles are radians, distances are metres, frequencies are hertz and GPS
//! time is represented as a full week number and seconds of week. GPS time
//! does not include UTC leap seconds. This crate currently provides the core
//! foundation, streaming RINEX 3 input, broadcast navigation and initial
//! GPS/QZSS L1 SPP, LAMBDA integer ambiguity search, and an initial static
//! GPS L1 C/A FLOAT filter with initial full-set integer FIX validation.
//! SP3-c/d and RINEX satellite CLK products provide interpolated precise
//! positions, velocities and published clocks with explicit coverage/quality.
//! Precise transmit states feed an initial code-only SPP with explicit signal
//! corrections. No phase PPP is implied by a precise code position.
//! A separate dual-frequency code/phase observation builder retains signal
//! corrections, covariance and raw GF/MW diagnostics, without a PPP estimator.
//! Transactional phase-arc tracking commits only accepted measurements and
//! reports continuity resets; it does not estimate an ambiguity value.
//! An initial static GPS IF PPP FLOAT filter connects precise states, phase arcs
//! and Niell-mapped zenith delay. Full upstream RTK/PPP policies and real-data
//! accuracy/runtime parity are not implemented or established yet.
//! IERS 2010 Step-1+2 is an explicit-body component; automatic astronomical
//! inputs and its PPP integration are still separate pending dependencies.
//!
//! ```
//! use gnss_rust::{GeodeticCoord, GnssTime};
//! use gnss_rust::coordinates::{geodetic_to_ecef, ecef_to_geodetic};
//!
//! let origin = GeodeticCoord::new(35_f64.to_radians(), 139_f64.to_radians(), 45.0)?;
//! let ecef = geodetic_to_ecef(origin)?;
//! let round_trip = ecef_to_geodetic(ecef)?;
//! assert!((round_trip.height - 45.0).abs() < 1e-6);
//! let time = GnssTime::new(2300, 604_799.5)?.checked_add_seconds(1.0)?;
//! assert_eq!(time.week(), 2301);
//! assert_eq!(time.tow(), 0.5);
//! # Ok::<(), gnss_rust::Error>(())
//! ```

pub mod antenna;
pub mod bias;
pub mod constants;
pub mod coordinates;
pub mod dual_frequency;
pub mod earth_rotation;
pub mod filter;
pub mod iers2010;
pub mod lambda;
pub mod models;
pub mod navigation;
pub mod observation;
pub mod phase_arc;
pub mod ppp;
pub mod precise;
pub mod precise_spp;
pub mod precise_transmit;
pub mod rinex;
pub mod rtk;
pub mod rtk_fix;
pub mod spp;
pub mod tides;
pub mod time;
pub mod types;
pub mod windup;

pub use time::GnssTime;
pub use types::{
    EcefCoord, EnuCoord, GeodeticCoord, GnssSystem, PositionSolution, PositioningMode, SatelliteId,
    SignalType, SolutionStatus,
};

/// Invalid input or a result outside the representable domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    NonFinite,
    InvalidLatitude,
    UndefinedGeodeticOrigin,
    NonConvergent,
    OutOfRange,
    InvalidSatelliteId,
    MissingGlonassChannel,
    InvalidGlonassChannel,
    InvalidCalendar,
    UnsupportedConstellation,
    InvalidEphemeris,
    MissingEphemeris,
    InsufficientObservations,
    SingularGeometry,
    InvalidConfiguration,
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::NonFinite => "input or result is not finite",
            Self::InvalidLatitude => "latitude must be between -pi/2 and pi/2 radians",
            Self::UndefinedGeodeticOrigin => {
                "geodetic coordinates are undefined at the Earth centre"
            }
            Self::NonConvergent => "numerical iteration did not converge",
            Self::OutOfRange => "result is outside the supported range",
            Self::InvalidSatelliteId => "expected a known constellation and nonzero u8 PRN",
            Self::MissingGlonassChannel => "GLONASS FDMA frequency requires a channel",
            Self::InvalidGlonassChannel => "GLONASS channel must be in -7..=6",
            Self::InvalidCalendar => "invalid calendar date or time",
            Self::UnsupportedConstellation => "constellation is not supported by this operation",
            Self::InvalidEphemeris => "invalid broadcast ephemeris",
            Self::MissingEphemeris => "no healthy ephemeris within the requested age limit",
            Self::InsufficientObservations => "at least four usable satellites are required",
            Self::SingularGeometry => "satellite geometry is rank deficient",
            Self::InvalidConfiguration => "invalid solver configuration",
        })
    }
}

impl std::error::Error for Error {}
