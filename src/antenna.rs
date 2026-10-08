//! Absolute ANTEX 1.4 calibrations and receiver PCO/NOAZI-PCV.
//! ANTEX offsets are North/East/Up in millimetres; this API uses metres.
//! Nominal-yaw satellite IF-PCO; no satellite PCV, azimuth grid or RMS model.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.

use crate::coordinates::{ecef_to_enu, ecef_to_geodetic};
use crate::dual_frequency::{ionosphere_free_coefficients, tracking_frequency_hz};
use crate::precise::ClockSource;
use crate::time::{CalendarDateTime, TimeScale};
use crate::windup::{DipoleFrame, WindupError, nominal_yaw_satellite_frame};
use crate::{EcefCoord, EnuCoord, Error, GnssSystem, GnssTime, SatelliteId};
use std::{collections::BTreeMap, fmt, io::BufRead};

#[derive(Debug, Clone, PartialEq)]
pub enum AntennaError {
    Parse { line: usize, reason: String },
    Unsupported(String),
    MissingCalibration,
    AmbiguousCalibration,
    MissingFrequency(String),
    OutsideValidity,
    InvalidGeometry,
    DatumMismatch,
    Attitude(WindupError),
    Gnss(Error),
}
impl fmt::Display for AntennaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "antenna: {self:?}")
    }
}
impl std::error::Error for AntennaError {}
impl From<Error> for AntennaError {
    fn from(e: Error) -> Self {
        Self::Gnss(e)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AntennaFrequency {
    /// ANTEX North/East/Up, including satellite body-axis entries.
    pub offset_neu_m: [f64; 3],
    pub noazi_m: Vec<f64>,
}
#[derive(Debug, Clone, PartialEq)]
pub struct AntennaCalibration {
    pub antenna_type: String,
    pub serial: String,
    pub svn: String,
    pub source_line: usize,
    /// GPS calendar labels, inclusive endpoints; None means unbounded.
    pub valid_from: Option<GnssTime>,
    pub valid_until: Option<GnssTime>,
    pub zenith_start_deg: f64,
    pub zenith_end_deg: f64,
    pub zenith_step_deg: f64,
    pub frequencies: BTreeMap<String, AntennaFrequency>,
}
impl AntennaCalibration {
    pub fn contains(&self, time: GnssTime) -> bool {
        self.valid_from.is_none_or(|t| time >= t) && self.valid_until.is_none_or(|t| time <= t)
    }
    pub fn satellite(&self) -> Option<SatelliteId> {
        if self.serial.len() == 3 {
            self.serial.parse().ok()
        } else {
            None
        }
    }
    /// Native receiver interpvar convention: linear NOAZI, clamped endpoints.
    pub fn pcv_m(&self, frequency: &str, zenith_deg: f64) -> Result<f64, AntennaError> {
        if !zenith_deg.is_finite() {
            return Err(Error::NonFinite.into());
        }
        let f = self
            .frequencies
            .get(frequency)
            .ok_or_else(|| AntennaError::MissingFrequency(frequency.into()))?;
        if f.noazi_m.is_empty()
            || !self.zenith_step_deg.is_finite()
            || self.zenith_step_deg <= 0.0
            || !self.zenith_start_deg.is_finite()
            || f.noazi_m.iter().any(|v| !v.is_finite())
        {
            return Err(AntennaError::InvalidGeometry);
        }
        let last = f.noazi_m.len() - 1;
        let node =
            ((zenith_deg - self.zenith_start_deg) / self.zenith_step_deg).clamp(0.0, last as f64);
        let i = node.floor() as usize;
        if i == last {
            return Ok(f.noazi_m[last]);
        }
        let fraction = node - i as f64;
        Ok(f.noazi_m[i] * (1.0 - fraction) + f.noazi_m[i + 1] * fraction)
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct Antex {
    entries: Vec<AntennaCalibration>,
}
impl Antex {
    pub fn entries(&self) -> &[AntennaCalibration] {
        &self.entries
    }
    /// PRN, exact SVN and GPS validity must identify a single entry. Selection
    /// is frozen by SatelliteAntennaModel; no automatic PRN/SVN handover.
    pub fn satellite_calibration(
        &self,
        satellite: SatelliteId,
        svn: &str,
        time: GnssTime,
    ) -> Result<&AntennaCalibration, AntennaError> {
        let mut found = self.entries.iter().filter(|e| {
            e.satellite() == Some(satellite) && e.svn == svn.trim() && e.contains(time)
        });
        let result = found.next().ok_or(AntennaError::MissingCalibration)?;
        if found.next().is_some() {
            return Err(AntennaError::AmbiguousCalibration);
        }
        Ok(result)
    }
    /// Match the complete type/radome and an exact serial. A generic blank
    /// serial must be explicitly requested; there is no fallback selection.
    pub fn receiver(
        &self,
        antenna_type: &str,
        serial: &str,
        time: GnssTime,
    ) -> Result<&AntennaCalibration, AntennaError> {
        let name = normalize_type(antenna_type);
        let mut found = self.entries.iter().filter(|e| {
            e.satellite().is_none()
                && e.antenna_type == name
                && e.serial == serial.trim()
                && e.contains(time)
        });
        let result = found.next().ok_or(AntennaError::MissingCalibration)?;
        if found.next().is_some() {
            return Err(AntennaError::AmbiguousCalibration);
        }
        Ok(result)
    }
}

/// Caller assertion that the selected ANTEX model matches the product clock
/// datum. This initial GPS path uses one explicit tracking pair for all PRNs.
#[derive(Debug, Clone, PartialEq)]
pub struct SatelliteAntennaBinding {
    pub product_id: String,
    pub clock_source: ClockSource,
    pub clock_reference: String,
    pub tracking_codes: [String; 2],
}
#[derive(Debug, Clone, PartialEq)]
pub struct SatelliteAntennaModel {
    binding: SatelliteAntennaBinding,
    frequency_codes: [String; 2],
    coefficients: [f64; 2],
    calibrations: BTreeMap<SatelliteId, AntennaCalibration>,
}
#[derive(Debug, Clone, PartialEq)]
pub struct SatelliteAntennaCorrection {
    pub product_id: String,
    pub satellite: SatelliteId,
    pub svn: String,
    pub antenna_type: String,
    pub source_line: usize,
    pub frequency_codes: [String; 2],
    pub coefficients: [f64; 2],
    /// Native ANTEX x/y/z body offsets; never exchange North and East here.
    pub body_offset_m: [[f64; 3]; 2],
    pub body_if_m: [f64; 3],
    pub frame: DipoleFrame,
    pub sun_position_m: EcefCoord,
    pub offset_ecef_m: EcefCoord,
    pub position_apc_m: EcefCoord,
}
impl SatelliteAntennaModel {
    pub fn new(
        binding: SatelliteAntennaBinding,
        calibrations: Vec<AntennaCalibration>,
    ) -> Result<Self, AntennaError> {
        if binding.product_id.trim().is_empty()
            || binding.clock_reference.trim().is_empty()
            || calibrations.is_empty()
        {
            return Err(AntennaError::DatumMismatch);
        }
        let gps: SatelliteId = "G01".parse().unwrap();
        let mut hz = [0.0; 2];
        for (i, code) in binding.tracking_codes.iter().enumerate() {
            hz[i] = tracking_frequency_hz(gps, code, None)
                .map_err(|_| AntennaError::Unsupported("satellite tracking code".into()))?;
        }
        let coefficients = ionosphere_free_coefficients(hz[0], hz[1])
            .map_err(|_| AntennaError::Unsupported("satellite frequency pair".into()))?;
        let frequency_codes = binding
            .tracking_codes
            .clone()
            .map(|code| format!("G0{}", code.as_bytes()[0] as char));
        let mut selected = BTreeMap::new();
        for c in calibrations {
            let satellite = c
                .satellite()
                .filter(|s| s.system() == GnssSystem::Gps && s.prn() <= 32)
                .ok_or_else(|| AntennaError::Unsupported("satellite calibration PRN".into()))?;
            if c.antenna_type.trim().is_empty()
                || c.svn.trim().is_empty()
                || matches!((c.valid_from,c.valid_until),(Some(a),Some(b)) if a>b)
            {
                return Err(AntennaError::InvalidGeometry);
            }
            for code in &frequency_codes {
                let f = c
                    .frequencies
                    .get(code)
                    .ok_or_else(|| AntennaError::MissingFrequency(code.clone()))?;
                if f.offset_neu_m.iter().any(|v| !v.is_finite()) {
                    return Err(Error::NonFinite.into());
                }
            }
            if selected.insert(satellite, c).is_some() {
                return Err(AntennaError::AmbiguousCalibration);
            }
        }
        Ok(Self {
            binding,
            frequency_codes,
            coefficients,
            calibrations: selected,
        })
    }
    pub fn binding(&self) -> &SatelliteAntennaBinding {
        &self.binding
    }
    pub fn calibrations(&self) -> &BTreeMap<SatelliteId, AntennaCalibration> {
        &self.calibrations
    }
    pub fn check_datum(
        &self,
        codes: &[String; 2],
        source: ClockSource,
        reference: &str,
    ) -> Result<(), AntennaError> {
        if *codes != self.binding.tracking_codes
            || source != self.binding.clock_source
            || reference != self.binding.clock_reference
        {
            return Err(AntennaError::DatumMismatch);
        }
        Ok(())
    }
    /// Time is GPS reception time, as in native ANTEX selection. Rotate the
    /// IF body PCO at transmit COM; shift position only, leaving COM velocity,
    /// published clock and periodic relativity untouched. Satellite PCV unused.
    pub fn correction(
        &self,
        time: GnssTime,
        satellite: SatelliteId,
        position_com_m: EcefCoord,
        sun_position_m: EcefCoord,
    ) -> Result<SatelliteAntennaCorrection, AntennaError> {
        let c = self
            .calibrations
            .get(&satellite)
            .ok_or(AntennaError::MissingCalibration)?;
        if !c.contains(time) {
            return Err(AntennaError::OutsideValidity);
        }
        let body_offset_m = self
            .frequency_codes
            .clone()
            .map(|code| c.frequencies[&code].offset_neu_m);
        let body_if_m: [f64; 3] = std::array::from_fn(|i| {
            self.coefficients[0] * body_offset_m[0][i] + self.coefficients[1] * body_offset_m[1][i]
        });
        let frame = nominal_yaw_satellite_frame(position_com_m, sun_position_m)
            .map_err(AntennaError::Attitude)?;
        let z = [
            frame.x[1] * frame.y[2] - frame.x[2] * frame.y[1],
            frame.x[2] * frame.y[0] - frame.x[0] * frame.y[2],
            frame.x[0] * frame.y[1] - frame.x[1] * frame.y[0],
        ];
        let offset: [f64; 3] = std::array::from_fn(|i| {
            body_if_m[0] * frame.x[i] + body_if_m[1] * frame.y[i] + body_if_m[2] * z[i]
        });
        let offset_ecef_m = EcefCoord::new(offset[0], offset[1], offset[2]);
        let position_apc_m = EcefCoord::new(
            position_com_m.x + offset[0],
            position_com_m.y + offset[1],
            position_com_m.z + offset[2],
        );
        if body_if_m
            .iter()
            .chain(offset.iter())
            .chain([&position_apc_m.x, &position_apc_m.y, &position_apc_m.z])
            .any(|v| !v.is_finite())
        {
            return Err(Error::NonFinite.into());
        }
        Ok(SatelliteAntennaCorrection {
            product_id: self.binding.product_id.clone(),
            satellite,
            svn: c.svn.clone(),
            antenna_type: c.antenna_type.clone(),
            source_line: c.source_line,
            frequency_codes: self.frequency_codes.clone(),
            coefficients: self.coefficients,
            body_offset_m,
            body_if_m,
            frame,
            sun_position_m,
            offset_ecef_m,
            position_apc_m,
        })
    }
}
fn normalize_type(value: &str) -> String {
    value
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
        .to_ascii_uppercase()
}
fn parse_error(line: usize, reason: &str) -> AntennaError {
    AntennaError::Parse {
        line,
        reason: reason.into(),
    }
}
fn numbers(value: &str, line: usize) -> Result<Vec<f64>, AntennaError> {
    value
        .split_whitespace()
        .map(|v| {
            v.parse::<f64>()
                .ok()
                .filter(|x| x.is_finite())
                .ok_or_else(|| parse_error(line, "invalid finite number"))
        })
        .collect()
}
fn epoch(value: &str, line: usize) -> Result<GnssTime, AntennaError> {
    let n = numbers(value, line)?;
    if n.len() != 6
        || n[..5].iter().any(|v| v.fract() != 0.0)
        || !(1.0..=9999.0).contains(&n[0])
        || n[1..5].iter().any(|v| !(0.0..=255.0).contains(v))
    {
        return Err(parse_error(line, "invalid calendar fields"));
    }
    CalendarDateTime {
        year: n[0] as i32,
        month: n[1] as u8,
        day: n[2] as u8,
        hour: n[3] as u8,
        minute: n[4] as u8,
        second: n[5],
    }
    .to_gnss_time(TimeScale::Gps)
    .map_err(|_| parse_error(line, "invalid GPS calendar"))
}

type FrequencyFields = (String, Option<[f64; 3]>, Option<Vec<f64>>);

/// Strict initial subset. Mixed satellite/receiver files are retained, but
/// azimuth-dependent grids, relative PCV and frequency RMS blocks are rejected.
/// Missing PCO/NOAZI, malformed blocks and incomplete files never become zero.
pub fn read_antex(reader: impl BufRead) -> Result<Antex, AntennaError> {
    let mut entries = Vec::new();
    let mut header = true;
    let mut version = false;
    let mut absolute = false;
    let mut current: Option<AntennaCalibration> = None;
    let mut frequency: Option<FrequencyFields> = None;
    let mut expected_frequencies = None;
    let mut grid_seen = false;
    let mut dazi_seen = false;
    let mut type_seen = false;
    let mut line_number = 0;
    for (index, input) in reader.lines().enumerate() {
        line_number = index + 1;
        let line = input.map_err(|e| parse_error(line_number, &e.to_string()))?;
        if !line.is_ascii() {
            return Err(parse_error(line_number, "non-ASCII input"));
        }
        let label = if line.trim_start().starts_with("NOAZI") {
            ""
        } else {
            line.get(60..).unwrap_or("").trim()
        };
        let data = &line[..line.len().min(60)];
        if header {
            match label {
                "ANTEX VERSION / SYST" => {
                    if version || data.get(..8).unwrap_or("").trim() != "1.4" {
                        return Err(AntennaError::Unsupported("ANTEX version".into()));
                    }
                    version = true;
                }
                "PCV TYPE / REFANT" => {
                    if absolute || !data.starts_with('A') {
                        return Err(AntennaError::Unsupported(
                            "relative/duplicate PCV declaration".into(),
                        ));
                    }
                    absolute = true;
                }
                "END OF HEADER" if version && absolute => header = false,
                "COMMENT" => (),
                _ => return Err(parse_error(line_number, "invalid header")),
            }
            continue;
        }
        if label == "COMMENT" || line.trim().is_empty() {
            continue;
        }
        if label == "START OF ANTENNA" {
            if current.is_some() {
                return Err(parse_error(line_number, "nested antenna"));
            }
            current = Some(AntennaCalibration {
                antenna_type: String::new(),
                serial: String::new(),
                svn: String::new(),
                source_line: line_number,
                valid_from: None,
                valid_until: None,
                zenith_start_deg: 0.0,
                zenith_end_deg: 0.0,
                zenith_step_deg: 0.0,
                frequencies: BTreeMap::new(),
            });
            expected_frequencies = None;
            grid_seen = false;
            dazi_seen = false;
            type_seen = false;
            continue;
        }
        let entry = current
            .as_mut()
            .ok_or_else(|| parse_error(line_number, "record outside antenna"))?;
        match label {
            "TYPE / SERIAL NO" if !type_seen && frequency.is_none() => {
                entry.antenna_type = normalize_type(
                    line.get(..20)
                        .ok_or_else(|| parse_error(line_number, "short type"))?,
                );
                entry.serial = line[20..40].trim().into();
                entry.svn = line[40..50].trim().into();
                type_seen = !entry.antenna_type.is_empty();
            }
            "DAZI" if !dazi_seen && frequency.is_none() => {
                if numbers(data, line_number)? != [0.0] {
                    return Err(AntennaError::Unsupported("azimuth-dependent PCV".into()));
                }
                dazi_seen = true;
            }
            "ZEN1 / ZEN2 / DZEN" if !grid_seen && frequency.is_none() => {
                let n = numbers(data, line_number)?;
                if n.len() != 3
                    || n[0] < 0.0
                    || n[1] > 180.0
                    || n[1] < n[0]
                    || n[2] <= 0.0
                    || ((n[1] - n[0]) / n[2] - ((n[1] - n[0]) / n[2]).round()).abs() > 1e-9
                    || (n[1] - n[0]) / n[2] > 10000.0
                {
                    return Err(parse_error(line_number, "invalid zenith grid"));
                }
                entry.zenith_start_deg = n[0];
                entry.zenith_end_deg = n[1];
                entry.zenith_step_deg = n[2];
                grid_seen = true;
            }
            "# OF FREQUENCIES" if expected_frequencies.is_none() && frequency.is_none() => {
                expected_frequencies = Some(
                    data.trim()
                        .parse::<usize>()
                        .map_err(|_| parse_error(line_number, "invalid frequency count"))?,
                );
            }
            "VALID FROM" if entry.valid_from.is_none() && frequency.is_none() => {
                entry.valid_from = Some(epoch(data, line_number)?)
            }
            "VALID UNTIL" if entry.valid_until.is_none() && frequency.is_none() => {
                entry.valid_until = Some(epoch(data, line_number)?)
            }
            "START OF FREQUENCY" if frequency.is_none() && grid_seen => {
                let key = data.trim();
                if key.len() != 3
                    || !b"GRECJSI".contains(&key.as_bytes()[0])
                    || !key.as_bytes()[1..].iter().all(u8::is_ascii_digit)
                {
                    return Err(parse_error(line_number, "invalid frequency code"));
                }
                frequency = Some((key.into(), None, None));
            }
            "NORTH / EAST / UP" => {
                let (_, pco, _) = frequency
                    .as_mut()
                    .ok_or_else(|| parse_error(line_number, "PCO outside frequency"))?;
                let n = numbers(data, line_number)?;
                if pco.is_some() || n.len() != 3 {
                    return Err(parse_error(line_number, "duplicate/invalid PCO"));
                }
                *pco = Some([n[0] * 1e-3, n[1] * 1e-3, n[2] * 1e-3]);
            }
            "END OF FREQUENCY" => {
                let (key, pco, pcv) = frequency
                    .take()
                    .ok_or_else(|| parse_error(line_number, "unmatched frequency end"))?;
                if data.trim() != key {
                    return Err(parse_error(line_number, "frequency end mismatch"));
                }
                let f = AntennaFrequency {
                    offset_neu_m: pco.ok_or_else(|| parse_error(line_number, "missing PCO"))?,
                    noazi_m: pcv.ok_or_else(|| parse_error(line_number, "missing NOAZI"))?,
                };
                if entry.frequencies.insert(key, f).is_some() {
                    return Err(parse_error(line_number, "duplicate frequency"));
                }
            }
            "END OF ANTENNA" if frequency.is_none() => {
                if !type_seen
                    || !grid_seen
                    || !dazi_seen
                    || entry.frequencies.is_empty()
                    || expected_frequencies != Some(entry.frequencies.len())
                    || matches!((entry.valid_from, entry.valid_until), (Some(a),Some(b)) if a > b)
                {
                    return Err(parse_error(line_number, "incomplete/invalid antenna"));
                }
                entries.push(current.take().unwrap());
            }
            "METH / BY / # / DATE" | "SINEX CODE" if frequency.is_none() => (),
            "START OF FREQ RMS" => {
                return Err(AntennaError::Unsupported("frequency RMS block".into()));
            }
            "" if line.trim_start().starts_with("NOAZI") => {
                let (_, _, pcv) = frequency
                    .as_mut()
                    .ok_or_else(|| parse_error(line_number, "NOAZI outside frequency"))?;
                let n = numbers(
                    line.trim_start().strip_prefix("NOAZI").unwrap(),
                    line_number,
                )?;
                let expected = ((entry.zenith_end_deg - entry.zenith_start_deg)
                    / entry.zenith_step_deg)
                    .round() as usize
                    + 1;
                if pcv.is_some() || n.len() != expected {
                    return Err(parse_error(line_number, "duplicate/incorrect NOAZI count"));
                }
                *pcv = Some(n.into_iter().map(|v| v * 1e-3).collect());
            }
            _ => return Err(parse_error(line_number, "unsupported/unexpected record")),
        }
    }
    if header || current.is_some() || frequency.is_some() || entries.is_empty() {
        return Err(parse_error(line_number, "incomplete ANTEX"));
    }
    Ok(Antex { entries })
}

/// Immutable selection for one filter lifetime. Changing calibration or
/// antenna geometry requires a new filter; expiry interrupts the phase arcs.
#[derive(Debug, Clone, PartialEq)]
pub struct ReceiverAntennaModel {
    product_id: String,
    calibration: AntennaCalibration,
    marker_delta_enu_m: EnuCoord,
}
#[derive(Debug, Clone, PartialEq)]
pub struct ReceiverAntennaCorrection {
    pub product_id: String,
    pub antenna_type: String,
    pub serial: String,
    pub source_line: usize,
    pub frequency_codes: [String; 2],
    pub pcv_m: [f64; 2],
    pub pco_projection_m: [f64; 2],
    /// Add to code and phase, once, after OSB; marker range convention.
    pub observation_add_m: [f64; 2],
    pub if_add_m: f64,
}
impl ReceiverAntennaModel {
    pub fn new(
        product_id: String,
        calibration: &AntennaCalibration,
        marker_delta_enu_m: EnuCoord,
    ) -> Result<Self, AntennaError> {
        if product_id.trim().is_empty()
            || calibration.satellite().is_some()
            || ![
                marker_delta_enu_m.east,
                marker_delta_enu_m.north,
                marker_delta_enu_m.up,
            ]
            .iter()
            .all(|v| v.is_finite())
        {
            return Err(AntennaError::InvalidGeometry);
        }
        // Calibration fields are public for provenance, so validate mutations.
        let size = (calibration.zenith_end_deg - calibration.zenith_start_deg)
            / calibration.zenith_step_deg;
        if !size.is_finite()
            || !(0.0..=10000.0).contains(&size)
            || (size - size.round()).abs() > 1e-9
            || calibration.zenith_step_deg <= 0.0
            || calibration.antenna_type.trim().is_empty()
            || calibration.zenith_start_deg < 0.0
            || calibration.zenith_end_deg > 180.0
            || matches!((calibration.valid_from, calibration.valid_until),(Some(a),Some(b)) if a>b)
            || calibration.frequencies.is_empty()
            || calibration.frequencies.values().any(|f| {
                f.noazi_m.len() != size.round() as usize + 1
                    || f.noazi_m
                        .iter()
                        .chain(f.offset_neu_m.iter())
                        .any(|v| !v.is_finite())
            })
        {
            return Err(AntennaError::InvalidGeometry);
        }
        Ok(Self {
            product_id,
            calibration: calibration.clone(),
            marker_delta_enu_m,
        })
    }
    pub fn calibration(&self) -> &AntennaCalibration {
        &self.calibration
    }
    /// Receiver fixed to local North/East/Up. Initial GPS 1/2/5 mapping only.
    pub fn correction(
        &self,
        time: GnssTime,
        satellite: SatelliteId,
        tracking_codes: &[String; 2],
        coefficients: [f64; 2],
        positions: [EcefCoord; 2],
    ) -> Result<ReceiverAntennaCorrection, AntennaError> {
        let [marker, satellite_position] = positions;
        if !self.calibration.contains(time) {
            return Err(AntennaError::OutsideValidity);
        }
        if satellite.system() != GnssSystem::Gps {
            return Err(AntennaError::Unsupported("receiver signal system".into()));
        }
        if coefficients.iter().any(|v| !v.is_finite()) {
            return Err(Error::NonFinite.into());
        }
        let geo = ecef_to_geodetic(marker)?;
        let los = ecef_to_enu(
            EcefCoord::new(
                satellite_position.x - marker.x,
                satellite_position.y - marker.y,
                satellite_position.z - marker.z,
            ),
            geo,
        )?;
        let norm = los.east.hypot(los.north).hypot(los.up);
        if !norm.is_finite() || norm < 1.0 {
            return Err(AntennaError::InvalidGeometry);
        }
        let zenith = 90.0 - los.up.atan2(los.east.hypot(los.north)).to_degrees();
        let mut frequency_codes = [String::new(), String::new()];
        let mut pcv_m = [0.0; 2];
        let mut pco_projection_m = [0.0; 2];
        for i in 0..2 {
            let code = &tracking_codes[i];
            if code.len() != 2
                || !b"125".contains(&code.as_bytes()[0])
                || !code.as_bytes()[1].is_ascii_uppercase()
            {
                return Err(AntennaError::Unsupported("receiver tracking code".into()));
            }
            frequency_codes[i] = format!("G0{}", code.as_bytes()[0] as char);
            let f = self
                .calibration
                .frequencies
                .get(&frequency_codes[i])
                .ok_or_else(|| AntennaError::MissingFrequency(frequency_codes[i].clone()))?;
            let [north, east, up] = f.offset_neu_m;
            let d = self.marker_delta_enu_m;
            pco_projection_m[i] =
                ((east + d.east) * los.east + (north + d.north) * los.north + (up + d.up) * los.up)
                    / norm;
            pcv_m[i] = self.calibration.pcv_m(&frequency_codes[i], zenith)?;
        }
        let observation_add_m = [
            pco_projection_m[0] - pcv_m[0],
            pco_projection_m[1] - pcv_m[1],
        ];
        let if_add_m =
            coefficients[0] * observation_add_m[0] + coefficients[1] * observation_add_m[1];
        if !if_add_m.is_finite() {
            return Err(Error::NonFinite.into());
        }
        Ok(ReceiverAntennaCorrection {
            product_id: self.product_id.clone(),
            antenna_type: self.calibration.antenna_type.clone(),
            serial: self.calibration.serial.clone(),
            source_line: self.calibration.source_line,
            frequency_codes,
            pcv_m,
            pco_projection_m,
            observation_add_m,
            if_add_m,
        })
    }
}
