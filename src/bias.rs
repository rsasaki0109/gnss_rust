//! Satellite Bias-SINEX 1.00 OSB/DSB with explicit time/datum admission.
//! Native reference: core/ionex_dcb_products.cpp and PPP code-bias conversion.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
//! DSB is retained, never interpreted as an absolute OSB without a gauge.

use crate::dual_frequency::{
    DualFrequencyError, SignalCorrection, SignalCorrections, tracking_frequency_hz,
};
use crate::precise::ClockSource;
use crate::time::{CalendarDateTime, TimeScale};
use crate::{GnssTime, SatelliteId, constants::SPEED_OF_LIGHT};
use std::{collections::BTreeMap, fmt, io::BufRead};

#[derive(Debug)]
pub enum BiasError {
    Io(std::io::Error),
    Parse {
        line: usize,
        message: String,
    },
    Unsupported {
        line: usize,
        feature: String,
    },
    MissingTimeScale,
    ConflictingTimeScale,
    OverlappingEntries {
        first_line: usize,
        second_line: usize,
    },
    InvalidBinding,
    MissingOsb {
        satellite: SatelliteId,
        observation: String,
    },
    ExpiredOsb {
        satellite: SatelliteId,
        observation: String,
    },
    AmbiguousOsb {
        satellite: SatelliteId,
        observation: String,
    },
    RelativeProduct,
    Correction(DualFrequencyError),
}
impl fmt::Display for BiasError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Bias-SINEX: {self:?}")
    }
}
impl std::error::Error for BiasError {}
impl From<std::io::Error> for BiasError {
    fn from(value: std::io::Error) -> Self {
        Self::Io(value)
    }
}
impl From<DualFrequencyError> for BiasError {
    fn from(value: DualFrequencyError) -> Self {
        Self::Correction(value)
    }
}
fn bad(line: usize, message: impl Into<String>) -> BiasError {
    BiasError::Parse {
        line,
        message: message.into(),
    }
}
fn unsupported(line: usize, feature: impl Into<String>) -> BiasError {
    BiasError::Unsupported {
        line,
        feature: feature.into(),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BiasKind {
    Osb,
    Dsb,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BiasUnit {
    Nanoseconds,
    Metres,
    Cycles,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BiasEntry {
    pub line: usize,
    pub kind: BiasKind,
    pub svn: Option<String>,
    pub satellite: SatelliteId,
    /// Receiver-specific rows are retained but never used as satellite OSBs.
    pub station: Option<String>,
    pub observation_1: String,
    pub observation_2: Option<String>,
    /// Half-open [start, end). None is the explicit 0000:000:00000 sentinel.
    pub valid_from: Option<GnssTime>,
    pub valid_until: Option<GnssTime>,
    pub start_label: String,
    pub end_label: String,
    pub unit: BiasUnit,
    pub unit_label: String,
    pub value: f64,
    pub sigma: f64,
}
impl BiasEntry {
    pub fn contains(&self, time: GnssTime) -> bool {
        self.valid_from.is_none_or(|t| t <= time) && self.valid_until.is_none_or(|t| time < t)
    }
    /// Only cycle-valued phase biases require a frequency/channel declaration.
    pub fn metres(&self, glonass_channel: Option<i8>) -> Result<(f64, f64), BiasError> {
        if !self.value.is_finite() || !self.sigma.is_finite() || self.sigma < 0.0 {
            return Err(bad(self.line, "invalid bias value or sigma"));
        }
        if self.unit == BiasUnit::Cycles
            && (self.observation_1.len() != 3
                || !self.observation_1.is_ascii()
                || !self.observation_1.starts_with('L'))
        {
            return Err(bad(self.line, "invalid cycle-valued phase observable"));
        }
        let scale = match self.unit {
            BiasUnit::Nanoseconds => SPEED_OF_LIGHT * 1e-9,
            BiasUnit::Metres => 1.0,
            BiasUnit::Cycles => {
                SPEED_OF_LIGHT
                    / tracking_frequency_hz(
                        self.satellite,
                        &self.observation_1[1..],
                        glonass_channel,
                    )?
            }
        };
        let value = self.value * scale;
        let sigma = self.sigma * scale;
        if !value.is_finite() || !sigma.is_finite() {
            return Err(bad(self.line, "metre conversion overflow"));
        }
        Ok((value, sigma))
    }
}
#[derive(Debug, Clone)]
pub struct BiasProducts {
    header: String,
    description: BTreeMap<String, Vec<Vec<String>>>,
    time_scale: TimeScale,
    entries: Vec<BiasEntry>,
    osb_index: BTreeMap<(SatelliteId, String), Vec<usize>>,
}
impl BiasProducts {
    pub fn header(&self) -> &str {
        &self.header
    }
    pub fn description(&self) -> &BTreeMap<String, Vec<Vec<String>>> {
        &self.description
    }
    pub fn time_scale(&self) -> TimeScale {
        self.time_scale
    }
    pub fn entries(&self) -> &[BiasEntry] {
        &self.entries
    }

    pub fn select_osb(
        &self,
        satellite: SatelliteId,
        observation: &str,
        time: GnssTime,
    ) -> Result<&BiasEntry, BiasError> {
        let rows: Vec<_> = self
            .osb_index
            .get(&(satellite, observation.to_owned()))
            .into_iter()
            .flatten()
            .map(|&index| &self.entries[index])
            .collect();
        let current: Vec<_> = rows.iter().copied().filter(|r| r.contains(time)).collect();
        if current.len() > 1
            || current.iter().any(|r| {
                r.observation_2
                    .as_ref()
                    .is_some_and(|other| other != observation)
            })
        {
            return Err(BiasError::AmbiguousOsb {
                satellite,
                observation: observation.into(),
            });
        }
        current.first().copied().ok_or_else(|| {
            if rows.is_empty() {
                BiasError::MissingOsb {
                    satellite,
                    observation: observation.into(),
                }
            } else {
                BiasError::ExpiredOsb {
                    satellite,
                    observation: observation.into(),
                }
            }
        })
    }

    /// Materialize corrections only at this reception epoch. This avoids
    /// turning a half-open product interval into the inclusive ledger interval.
    /// Clock/phase identities are explicit caller bindings, not inferred proof
    /// of product compatibility. Missing code OSBs fail; missing phase OSBs
    /// retain code without inventing zero phase corrections.
    pub fn corrections_at(
        &self,
        time: GnssTime,
        binding: &OsbBinding,
        requests: &[BiasSignalRequest],
    ) -> Result<BiasCorrectionSnapshot, BiasError> {
        if binding.product_id.trim().is_empty()
            || binding.clock_reference.trim().is_empty()
            || binding
                .phase_reference
                .as_ref()
                .is_some_and(|s| s.trim().is_empty())
        {
            return Err(BiasError::InvalidBinding);
        }
        if self
            .description
            .get("BIAS_MODE")
            .is_some_and(|v| v.len() != 1 || v[0].as_slice() != ["ABSOLUTE"])
            || self.header.split_whitespace().nth(7) == Some("R")
        {
            return Err(BiasError::RelativeProduct);
        }
        let mut entries = Vec::new();
        let mut provenance = BTreeMap::new();
        for request in requests {
            let tracking = &request.tracking_code;
            if tracking.len() != 2
                || !tracking.is_ascii()
                || !matches!(tracking.as_bytes()[0], b'1'..=b'9')
                || !tracking.as_bytes()[1].is_ascii_uppercase()
            {
                return Err(BiasError::InvalidBinding);
            }
            let code = self
                .select_osb(request.satellite, &format!("C{tracking}"), time)?
                .clone();
            let (phase, phase_unavailability) = if binding.phase_reference.is_some() {
                match self.select_osb(request.satellite, &format!("L{tracking}"), time) {
                    Ok(entry) => (Some(entry.clone()), None),
                    Err(BiasError::MissingOsb { .. }) => (None, Some(PhaseOsbUnavailable::Missing)),
                    Err(BiasError::ExpiredOsb { .. }) => (None, Some(PhaseOsbUnavailable::Expired)),
                    Err(error) => return Err(error),
                }
            } else {
                (None, Some(PhaseOsbUnavailable::NotBound))
            };
            let (code_m, code_sigma_m) = code.metres(request.glonass_channel)?;
            let phase_m = phase
                .as_ref()
                .map(|p| p.metres(request.glonass_channel))
                .transpose()?;
            entries.push(SignalCorrection {
                satellite: request.satellite,
                tracking_code: tracking.clone(),
                code_add_m: -code_m,
                phase_add_m: phase_m.map(|p| -p.0),
                clock_source: binding.clock_source,
                clock_reference: binding.clock_reference.clone(),
                phase_reference: phase.as_ref().and(binding.phase_reference.clone()),
                valid_from: time,
                valid_until: time,
            });
            provenance.insert(
                (request.satellite, tracking.clone()),
                AppliedOsb {
                    product_id: binding.product_id.clone(),
                    code,
                    phase,
                    phase_unavailability,
                    code_sigma_m,
                    phase_sigma_m: phase_m.map(|p| p.1),
                },
            );
        }
        Ok(BiasCorrectionSnapshot {
            time,
            corrections: SignalCorrections::new(entries)?,
            provenance,
        })
    }
}
#[derive(Debug, Clone)]
pub struct OsbBinding {
    pub product_id: String,
    pub clock_source: ClockSource,
    pub clock_reference: String,
    pub phase_reference: Option<String>,
}
#[derive(Debug, Clone)]
pub struct BiasSignalRequest {
    pub satellite: SatelliteId,
    pub tracking_code: String,
    pub glonass_channel: Option<i8>,
}
#[derive(Debug, Clone)]
pub struct AppliedOsb {
    pub product_id: String,
    pub code: BiasEntry,
    pub phase: Option<BiasEntry>,
    pub phase_unavailability: Option<PhaseOsbUnavailable>,
    /// Uncertainty is retained as provenance, not independently added to each
    /// epoch's measurement noise: a product bias error can persist over time.
    pub code_sigma_m: f64,
    pub phase_sigma_m: Option<f64>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseOsbUnavailable {
    Missing,
    Expired,
    NotBound,
}
#[derive(Debug, Clone)]
pub struct BiasCorrectionSnapshot {
    pub time: GnssTime,
    pub corrections: SignalCorrections,
    pub provenance: BTreeMap<(SatelliteId, String), AppliedOsb>,
}

/// If TIME_SYSTEM is missing, a caller declaration is mandatory. If present,
/// the declaration must agree. No GPS time scale is silently assumed.
#[derive(Debug, Clone, Copy, Default)]
pub struct BiasReadOptions {
    pub time_scale: Option<TimeScale>,
}

fn epoch(text: &str, scale: TimeScale, line: usize) -> Result<Option<GnssTime>, BiasError> {
    if text == "0000:000:00000" {
        return Ok(None);
    }
    let b = text.as_bytes();
    if b.len() != 14
        || b[4] != b':'
        || b[8] != b':'
        || !b
            .iter()
            .enumerate()
            .all(|(i, c)| i == 4 || i == 8 || c.is_ascii_digit())
    {
        return Err(bad(line, "expected YYYY:DDD:SSSSS"));
    }
    let year: i32 = text[..4].parse().map_err(|_| bad(line, "invalid year"))?;
    let day: u16 = text[5..8].parse().map_err(|_| bad(line, "invalid day"))?;
    let seconds: u32 = text[9..]
        .parse()
        .map_err(|_| bad(line, "invalid seconds"))?;
    let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    if year == 0 || day == 0 || day > if leap { 366 } else { 365 } || seconds >= 86400 {
        return Err(bad(line, "out-of-range bias epoch"));
    }
    // Convert the ordinal to a calendar date before applying UTC/BDT offsets;
    // adding GPS days across a UTC leap-second boundary would be incorrect.
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
    let mut remaining = day;
    let mut month = 1;
    for length in months {
        if remaining <= length {
            break;
        }
        remaining -= length;
        month += 1;
    }
    CalendarDateTime {
        year,
        month,
        day: remaining as u8,
        hour: (seconds / 3600) as u8,
        minute: (seconds / 60 % 60) as u8,
        second: f64::from(seconds % 60),
    }
    .to_gnss_time(scale)
    .map(Some)
    .map_err(|_| bad(line, "invalid bias epoch"))
}
fn number(text: &str, line: usize) -> Result<f64, BiasError> {
    let n: f64 = text
        .replace(['D', 'd'], "E")
        .parse()
        .map_err(|_| bad(line, "invalid bias number"))?;
    if n.is_finite() {
        Ok(n)
    } else {
        Err(bad(line, "non-finite bias number"))
    }
}
fn observation(text: &str, line: usize) -> Result<String, BiasError> {
    let b = text.as_bytes();
    if b.len() != 3
        || !matches!(b[0], b'C' | b'L')
        || !matches!(b[1], b'1'..=b'9')
        || !b[2].is_ascii_uppercase()
    {
        return Err(unsupported(line, "expected exact RINEX 3 C/L observable"));
    }
    Ok(text.to_owned())
}
fn row(text: &str, line: usize, scale: TimeScale) -> Result<BiasEntry, BiasError> {
    let f: Vec<_> = text.split_whitespace().collect();
    let kind = match f.first().copied() {
        Some("OSB") => BiasKind::Osb,
        Some("DSB") => BiasKind::Dsb,
        _ => return Err(unsupported(line, "bias type (only OSB/DSB)")),
    };
    let start = f
        .iter()
        .position(|s| s.contains(':'))
        .ok_or_else(|| bad(line, "missing bias epochs"))?;
    if f.len() != start + 5 {
        return Err(unsupported(
            line,
            "truncated row or bias slope/extra columns",
        ));
    }
    let mut index = 1;
    let svn = if f
        .get(index)
        .is_some_and(|s| s.len() == 4 && s.as_bytes()[1..].iter().all(u8::is_ascii_digit))
    {
        let s = Some(f[index].to_owned());
        index += 1;
        s
    } else {
        None
    };
    let sat_token = f.get(index).ok_or_else(|| bad(line, "missing PRN"))?;
    if sat_token.len() != 3 || !sat_token.as_bytes()[1..].iter().all(u8::is_ascii_digit) {
        return Err(unsupported(
            line,
            "satellite PRN (receiver-only rows unsupported)",
        ));
    }
    let satellite: SatelliteId = sat_token
        .parse()
        .map_err(|_| bad(line, "invalid satellite"))?;
    if svn
        .as_ref()
        .is_some_and(|s| s.as_bytes()[0] != sat_token.as_bytes()[0])
    {
        return Err(bad(line, "SVN/PRN constellation mismatch"));
    }
    index += 1;
    let station = if f
        .get(index)
        .is_some_and(|s| s.len() != 3 || !matches!(s.as_bytes()[0], b'C' | b'L'))
    {
        let s = f[index];
        if s.len() > 9 || !s.bytes().all(|c| c.is_ascii_alphanumeric()) {
            return Err(bad(line, "invalid station"));
        }
        index += 1;
        Some(s.to_owned())
    } else {
        None
    };
    if start < index + 1 || start > index + 2 {
        return Err(bad(line, "invalid observable columns"));
    }
    let observation_1 = observation(f[index], line)?;
    let observation_2 = if start == index + 2 && f[index + 1] != "-" {
        Some(observation(f[index + 1], line)?)
    } else {
        None
    };
    if kind == BiasKind::Dsb
        && (observation_2.is_none() || observation_2.as_ref() == Some(&observation_1))
    {
        return Err(bad(line, "DSB needs two distinct observables"));
    }
    let valid_from = epoch(f[start], scale, line)?;
    let valid_until = epoch(f[start + 1], scale, line)?;
    if valid_from.zip(valid_until).is_some_and(|(a, b)| a >= b) {
        return Err(bad(line, "empty/reversed bias interval"));
    }
    let unit_label = f[start + 2].to_owned();
    let unit = match unit_label.to_ascii_uppercase().as_str() {
        "NS" => BiasUnit::Nanoseconds,
        "M" | "METER" | "METERS" => BiasUnit::Metres,
        "CYC" => BiasUnit::Cycles,
        _ => return Err(unsupported(line, "bias unit")),
    };
    if unit == BiasUnit::Cycles && (kind != BiasKind::Osb || !observation_1.starts_with('L')) {
        return Err(unsupported(line, "cycle unit requires a phase OSB"));
    }
    let value = number(f[start + 3], line)?;
    let sigma = number(f[start + 4], line)?;
    if sigma < 0.0 {
        return Err(bad(line, "negative bias sigma"));
    }
    Ok(BiasEntry {
        line,
        kind,
        svn,
        satellite,
        station,
        observation_1,
        observation_2,
        valid_from,
        valid_until,
        start_label: f[start].into(),
        end_label: f[start + 1].into(),
        unit,
        unit_label,
        value,
        sigma,
    })
}

pub fn read_bias_sinex(
    reader: impl BufRead,
    options: BiasReadOptions,
) -> Result<BiasProducts, BiasError> {
    let mut header = None;
    let mut block: Option<String> = None;
    let mut description: BTreeMap<String, Vec<Vec<String>>> = BTreeMap::new();
    let mut raw_rows = Vec::new();
    let mut ended = false;
    let mut seen_solution = false;
    for (i, text) in reader.lines().enumerate() {
        let line = i + 1;
        let text = text?;
        if !text.is_ascii() {
            return Err(bad(line, "records must be ASCII"));
        }
        let text = text.trim();
        if text.is_empty() || text.starts_with('*') {
            continue;
        }
        if ended {
            return Err(bad(line, "data after %=ENDBIA"));
        }
        if header.is_none() {
            let f: Vec<_> = text.split_whitespace().collect();
            if f.first() != Some(&"%=BIA") || f.get(1) != Some(&"1.00") {
                return Err(unsupported(line, "Bias-SINEX version (only 1.00)"));
            }
            header = Some(text.to_owned());
            continue;
        }
        if text == "%=ENDBIA" {
            if block.is_some() {
                return Err(bad(line, "unclosed block"));
            }
            ended = true;
            continue;
        }
        if let Some(name) = text.strip_prefix('+') {
            if block.is_some() {
                return Err(bad(line, "nested block"));
            }
            if name == "BIAS/SOLUTION" {
                if seen_solution {
                    return Err(bad(line, "duplicate solution block"));
                }
                seen_solution = true;
            }
            block = Some(name.to_owned());
            continue;
        }
        if let Some(name) = text.strip_prefix('-') {
            if block.as_deref() != Some(name) {
                return Err(bad(line, "mismatched block end"));
            }
            block = None;
            continue;
        }
        match block.as_deref() {
            Some("BIAS/SOLUTION") => raw_rows.push((line, text.to_owned())),
            Some("BIAS/DESCRIPTION") => {
                let mut fields = text.split_whitespace();
                let key = fields
                    .next()
                    .ok_or_else(|| bad(line, "empty description"))?;
                if matches!(key, "TIME_SYSTEM" | "BIAS_MODE") && description.contains_key(key) {
                    return Err(bad(line, "duplicate description field"));
                }
                description
                    .entry(key.into())
                    .or_default()
                    .push(fields.map(str::to_owned).collect());
            }
            Some(_) => {} // Non-solution metadata blocks do not contain bias estimates.
            None => return Err(bad(line, "record outside block")),
        }
    }
    if !ended || !seen_solution || raw_rows.is_empty() {
        return Err(bad(0, "missing end marker or solution data"));
    }
    let declared = match description.get("TIME_SYSTEM").map(Vec::as_slice) {
        None => None,
        Some([values]) if values.len() == 1 => Some(match values[0].as_str() {
            "G" => TimeScale::Gps,
            "U" => TimeScale::Utc,
            "E" => TimeScale::Galileo,
            "C" => TimeScale::BeiDou,
            "J" => TimeScale::Qzss,
            _ => return Err(unsupported(0, "TIME_SYSTEM")),
        }),
        _ => return Err(bad(0, "invalid TIME_SYSTEM description")),
    };
    if declared
        .zip(options.time_scale)
        .is_some_and(|(a, b)| a != b)
    {
        return Err(BiasError::ConflictingTimeScale);
    }
    let time_scale = declared
        .or(options.time_scale)
        .ok_or(BiasError::MissingTimeScale)?;
    let mut entries: Vec<BiasEntry> = Vec::new();
    let mut groups = BTreeMap::new();
    let mut osb_index: BTreeMap<_, Vec<usize>> = BTreeMap::new();
    for (line, text) in raw_rows {
        let entry = row(&text, line, time_scale)?;
        let other = entry
            .observation_2
            .as_ref()
            .filter(|s| *s != &entry.observation_1)
            .cloned();
        let key = (
            entry.satellite,
            entry.station.clone(),
            entry.kind as u8,
            entry.observation_1.clone(),
            other,
        );
        groups
            .entry(key)
            .or_insert_with(Vec::new)
            .push(entries.len());
        if entry.kind == BiasKind::Osb && entry.station.is_none() {
            osb_index
                .entry((entry.satellite, entry.observation_1.clone()))
                .or_default()
                .push(entries.len());
        }
        entries.push(entry);
    }
    for indices in groups.values_mut() {
        indices.sort_by(|&a, &b| {
            entries[a]
                .valid_from
                .partial_cmp(&entries[b].valid_from)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for pair in indices.windows(2) {
            let a = &entries[pair[0]];
            let b = &entries[pair[1]];
            if !a
                .valid_until
                .zip(b.valid_from)
                .is_some_and(|(end, start)| end <= start)
            {
                return Err(BiasError::OverlappingEntries {
                    first_line: a.line,
                    second_line: b.line,
                });
            }
        }
    }
    Ok(BiasProducts {
        header: header.unwrap_or_default(),
        description,
        time_scale,
        entries,
        osb_index,
    })
}
