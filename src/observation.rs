//! Exact tracking-code observations; missing measurements are explicit.

use crate::{EcefCoord, GnssTime, SatelliteId};
use std::collections::BTreeSet;

/// One code/phase/doppler/SNR group (e.g. RINEX C1C/L1C/D1C/S1C).
#[derive(Debug, Clone, PartialEq)]
pub struct Observation {
    pub satellite: SatelliteId,
    /// Frequency digit and tracking code, for example `1C` or `2W`.
    pub tracking_code: String,
    pub pseudorange_m: Option<f64>,
    pub carrier_phase_cycles: Option<f64>,
    pub doppler_hz: Option<f64>,
    pub snr_db_hz: Option<f64>,
    /// Raw RINEX loss-of-lock bits; preserved even when phase is missing.
    pub lli: Option<u8>,
    pub signal_strength: Option<u8>,
    pub glonass_channel: Option<i8>,
}

impl Observation {
    pub fn new(satellite: SatelliteId, tracking_code: String) -> Self {
        Self {
            satellite,
            tracking_code,
            pseudorange_m: None,
            carrier_phase_cycles: None,
            doppler_hz: None,
            snr_db_hz: None,
            lli: None,
            signal_strength: None,
            glonass_channel: None,
        }
    }
    pub fn loss_of_lock(&self) -> bool {
        self.lli.is_some_and(|lli| lli & 1 != 0)
    }
    pub fn half_cycle_ambiguity(&self) -> bool {
        self.lli.is_some_and(|lli| lli & 2 != 0)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ObservationEpoch {
    pub time: GnssTime,
    /// 0 normal, 1 power failure, 6 cycle-slip record. Header/event records are
    /// handled by the reader and never returned as measurement epochs.
    pub flag: u8,
    pub approximate_position: Option<EcefCoord>,
    /// Raw optional epoch clock offset in seconds, not silently applied to the data.
    pub receiver_clock_offset_s: Option<f64>,
    pub observations: Vec<Observation>,
}

impl ObservationEpoch {
    pub fn satellites(&self) -> BTreeSet<SatelliteId> {
        self.observations.iter().map(|obs| obs.satellite).collect()
    }
}
