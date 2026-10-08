//! Streaming RINEX 3 observations and Kepler navigation records.
//! Unsupported record families are errors, not silently discarded data.

use crate::navigation::{BroadcastClock, BroadcastEphemeris, KeplerOrbit, NavigationData};
use crate::observation::{Observation, ObservationEpoch};
use crate::time::{CalendarDateTime, TimeScale};
use crate::{EcefCoord, GnssSystem, SatelliteId};
use std::{collections::BTreeMap, fmt, io::BufRead};

#[derive(Debug)]
pub enum RinexError {
    Io(std::io::Error),
    Parse { line: usize, message: String },
    Unsupported { line: usize, feature: String },
}
impl fmt::Display for RinexError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "RINEX I/O: {e}"),
            Self::Parse { line, message } => write!(f, "RINEX line {line}: {message}"),
            Self::Unsupported { line, feature } => {
                write!(f, "RINEX line {line}: unsupported {feature}")
            }
        }
    }
}
impl std::error::Error for RinexError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            _ => None,
        }
    }
}
impl From<std::io::Error> for RinexError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}
fn bad(line: usize, message: impl Into<String>) -> RinexError {
    RinexError::Parse {
        line,
        message: message.into(),
    }
}
fn unsupported(line: usize, feature: impl Into<String>) -> RinexError {
    RinexError::Unsupported {
        line,
        feature: feature.into(),
    }
}

struct Input<R> {
    reader: R,
    line: usize,
    pending: Option<(usize, String)>,
}
impl<R: BufRead> Input<R> {
    fn next_line(&mut self) -> Result<Option<(usize, String)>, RinexError> {
        if let Some(line) = self.pending.take() {
            return Ok(Some(line));
        }
        let mut text = String::new();
        if self.reader.read_line(&mut text)? == 0 {
            return Ok(None);
        }
        self.line += 1;
        if !text.is_ascii() {
            return Err(bad(self.line, "RINEX records must be ASCII"));
        }
        Ok(Some((
            self.line,
            text.trim_end_matches(['\n', '\r']).to_owned(),
        )))
    }
    fn required(&mut self) -> Result<(usize, String), RinexError> {
        self.next_line()?
            .ok_or_else(|| bad(self.line + 1, "unexpected end of file"))
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RinexHeader {
    pub version: f64,
    pub observation_types: BTreeMap<GnssSystem, Vec<String>>,
    pub approximate_position: Option<EcefCoord>,
    pub time_scale: TimeScale,
    pub clock_offsets_applied: bool,
    /// Ordered ledger, retaining duplicate GLONASS channel declarations.
    pub glonass_channels: Vec<(SatelliteId, i8)>,
    pub gps_ionosphere_alpha: Option<[f64; 4]>,
    pub gps_ionosphere_beta: Option<[f64; 4]>,
}

fn system(prefix: u8, line: usize) -> Result<GnssSystem, RinexError> {
    match prefix {
        b'G' => Ok(GnssSystem::Gps),
        b'R' => Ok(GnssSystem::Glonass),
        b'E' => Ok(GnssSystem::Galileo),
        b'C' => Ok(GnssSystem::BeiDou),
        b'J' => Ok(GnssSystem::Qzss),
        b'S' => Ok(GnssSystem::Sbas),
        b'I' => Ok(GnssSystem::NavIC),
        _ => Err(unsupported(line, "satellite system")),
    }
}
fn numeric(text: &str, line: usize) -> Result<f64, RinexError> {
    let value: f64 = text
        .trim()
        .replace(['D', 'd'], "E")
        .parse()
        .map_err(|_| bad(line, format!("invalid number {text:?}")))?;
    if !value.is_finite() {
        return Err(bad(line, "non-finite measurement"));
    }
    Ok(value)
}
fn integer(text: &str, line: usize) -> Result<i32, RinexError> {
    text.parse()
        .map_err(|_| bad(line, format!("invalid integer {text:?}")))
}
fn field(text: &str, start: usize, width: usize) -> &str {
    text.get(start..text.len().min(start + width)).unwrap_or("")
}
fn optional_numeric(text: &str, line: usize) -> Result<Option<f64>, RinexError> {
    if text.trim().is_empty() {
        Ok(None)
    } else {
        numeric(text, line).map(Some)
    }
}
fn unsigned_field(value: f64, line: usize) -> Result<u16, RinexError> {
    if value.fract() != 0.0 || !(0.0..=65535.0).contains(&value) {
        return Err(bad(line, "out-of-range integer field"));
    }
    Ok(value as u16)
}

fn read_header<R: BufRead>(input: &mut Input<R>, kind: u8) -> Result<RinexHeader, RinexError> {
    let (line, text) = input.required()?;
    if !text
        .get(60..)
        .is_some_and(|x| x.starts_with("RINEX VERSION / TYPE"))
    {
        return Err(bad(line, "missing RINEX VERSION / TYPE"));
    }
    let version = numeric(field(&text, 0, 9), line)?;
    if !(3.0..4.0).contains(&version) {
        return Err(unsupported(line, format!("RINEX version {version}")));
    }
    if text.as_bytes().get(20) != Some(&kind) {
        return Err(bad(line, "unexpected RINEX file type"));
    }
    let mut header = RinexHeader {
        version,
        observation_types: BTreeMap::new(),
        approximate_position: None,
        time_scale: TimeScale::Gps,
        clock_offsets_applied: false,
        glonass_channels: Vec::new(),
        gps_ionosphere_alpha: None,
        gps_ionosphere_beta: None,
    };
    loop {
        let (line, text) = input.required()?;
        let label = text
            .get(60..)
            .ok_or_else(|| bad(line, "short header record"))?
            .trim();
        if label == "END OF HEADER" {
            break;
        }
        apply_header_line(input, &mut header, line, &text)?;
    }
    if kind == b'O' && header.observation_types.is_empty() {
        return Err(bad(input.line, "missing observation types"));
    }
    Ok(header)
}

fn apply_header_line<R: BufRead>(
    input: &mut Input<R>,
    header: &mut RinexHeader,
    line: usize,
    text: &str,
) -> Result<(), RinexError> {
    let label = text
        .get(60..)
        .ok_or_else(|| bad(line, "short header record"))?
        .trim();
    match label {
        "APPROX POSITION XYZ" => {
            let values: Vec<_> = text[..60]
                .split_whitespace()
                .map(|s| numeric(s, line))
                .collect::<Result<_, _>>()?;
            if values.len() != 3 {
                return Err(bad(line, "expected three approximate coordinates"));
            }
            header.approximate_position = Some(EcefCoord::new(values[0], values[1], values[2]));
        }
        "SYS / # / OBS TYPES" => {
            let sys = system(text.as_bytes()[0], line)?;
            let count: usize = field(text, 3, 3)
                .trim()
                .parse()
                .map_err(|_| bad(line, "invalid observation type count"))?;
            if !(1..=999).contains(&count) {
                return Err(bad(line, "invalid observation type count"));
            }
            let mut types: Vec<String> = field(text, 7, 53)
                .split_whitespace()
                .map(str::to_owned)
                .collect();
            while types.len() < count {
                let (next_line, continuation) = input.required()?;
                if !continuation.get(60..).is_some_and(|x| x.trim() == label)
                    || !continuation
                        .as_bytes()
                        .first()
                        .is_some_and(|&x| x == b' ' || x == text.as_bytes()[0])
                {
                    return Err(bad(next_line, "invalid observation type continuation"));
                }
                let extra: Vec<_> = field(&continuation, 7, 53)
                    .split_whitespace()
                    .map(str::to_owned)
                    .collect();
                if extra.is_empty() {
                    return Err(bad(next_line, "empty observation type continuation"));
                }
                types.extend(extra);
            }
            if types.len() != count {
                return Err(bad(line, "observation type count mismatch"));
            }
            let mut seen = std::collections::BTreeSet::new();
            for code in &types {
                let b = code.as_bytes();
                if b.len() != 3
                    || !matches!(b[0], b'C' | b'L' | b'D' | b'S')
                    || !b[1].is_ascii_digit()
                    || !b[2].is_ascii_uppercase()
                    || !seen.insert(code)
                {
                    return Err(bad(
                        line,
                        format!("invalid or duplicate observation type {code}"),
                    ));
                }
            }
            header.observation_types.insert(sys, types);
        }
        "TIME OF FIRST OBS" => {
            header.time_scale = match field(text, 48, 3).trim() {
                "" | "GPS" => TimeScale::Gps,
                "GAL" => TimeScale::Galileo,
                "QZS" => TimeScale::Qzss,
                "BDT" => TimeScale::BeiDou,
                "UTC" | "GLO" => TimeScale::Utc,
                scale => return Err(unsupported(line, format!("time scale {scale}"))),
            };
        }
        "RCV CLOCK OFFS APPL" => {
            header.clock_offsets_applied = match field(text, 0, 6).trim() {
                "0" => false,
                "1" => true,
                _ => return Err(bad(line, "invalid receiver clock flag")),
            };
        }
        "GLONASS SLOT / FRQ #" => {
            let tokens: Vec<_> = text[..60].split_whitespace().collect();
            let count: usize = tokens
                .first()
                .ok_or_else(|| bad(line, "missing GLONASS count"))?
                .parse()
                .map_err(|_| bad(line, "invalid GLONASS count"))?;
            if count > 999 {
                return Err(bad(line, "invalid GLONASS count"));
            }
            let mut entries: Vec<String> = tokens[1..].iter().map(|x| (*x).to_owned()).collect();
            while entries.len() < count * 2 {
                let (next_line, continuation) = input.required()?;
                if !continuation.get(60..).is_some_and(|x| x.trim() == label) {
                    return Err(bad(next_line, "invalid GLONASS continuation"));
                }
                let extra: Vec<_> = continuation[..60]
                    .split_whitespace()
                    .map(str::to_owned)
                    .collect();
                if extra.is_empty() {
                    return Err(bad(next_line, "empty GLONASS continuation"));
                }
                entries.extend(extra);
            }
            if entries.len() != count * 2 {
                return Err(bad(line, "GLONASS count mismatch"));
            }
            for entry in entries.chunks_exact(2) {
                let sat: SatelliteId = entry[0]
                    .parse()
                    .map_err(|_| bad(line, "invalid GLONASS satellite"))?;
                let channel: i8 = entry[1]
                    .parse()
                    .map_err(|_| bad(line, "invalid GLONASS channel"))?;
                if sat.system() != GnssSystem::Glonass || !(-7..=6).contains(&channel) {
                    return Err(bad(line, "invalid GLONASS channel entry"));
                }
                header.glonass_channels.push((sat, channel));
            }
        }
        "IONOSPHERIC CORR" if matches!(field(text, 0, 4), "GPSA" | "GPSB") => {
            let mut coefficients = [0.0; 4];
            for (i, value) in coefficients.iter_mut().enumerate() {
                *value = numeric(field(text, 5 + i * 12, 12), line)?;
            }
            if text.starts_with("GPSA") {
                header.gps_ionosphere_alpha = Some(coefficients);
            } else {
                header.gps_ionosphere_beta = Some(coefficients);
            }
        }
        "SYS / SCALE FACTOR" => {
            if integer(field(text, 2, 4).trim(), line)? != 1 {
                return Err(unsupported(line, "non-unit observation scale factor"));
            }
        }
        "SYS / PHASE SHIFT" if numeric(field(text, 6, 8), line)? != 0.0 => {
            return Err(unsupported(line, "nonzero carrier phase shift"));
        }
        _ => {}
    }
    Ok(())
}

fn calendar(tokens: &[&str], line: usize) -> Result<CalendarDateTime, RinexError> {
    if tokens.len() != 6 {
        return Err(bad(line, "expected six calendar fields"));
    }
    let small = |i| {
        integer(tokens[i], line)?
            .try_into()
            .map_err(|_| bad(line, "out-of-range calendar field"))
    };
    Ok(CalendarDateTime {
        year: integer(tokens[0], line)?,
        month: small(1)?,
        day: small(2)?,
        hour: small(3)?,
        minute: small(4)?,
        second: numeric(tokens[5], line)?,
    })
}

pub struct RinexObservationReader<R> {
    input: Input<R>,
    header: RinexHeader,
    done: bool,
}
impl<R: BufRead> RinexObservationReader<R> {
    pub fn new(reader: R) -> Result<Self, RinexError> {
        let mut input = Input {
            reader,
            line: 0,
            pending: None,
        };
        let header = read_header(&mut input, b'O')?;
        Ok(Self {
            input,
            header,
            done: false,
        })
    }
    pub fn header(&self) -> &RinexHeader {
        &self.header
    }
    fn next_epoch(&mut self) -> Result<Option<ObservationEpoch>, RinexError> {
        loop {
            let Some((line, text)) = self.input.next_line()? else {
                return Ok(None);
            };
            if text.trim().is_empty() {
                continue;
            }
            if !text.starts_with('>') {
                return Err(bad(line, "expected epoch marker"));
            }
            let tokens: Vec<_> = text[1..].split_whitespace().collect();
            if !(8..=9).contains(&tokens.len()) {
                return Err(bad(line, "invalid epoch header"));
            }
            let flag: u8 = tokens[6]
                .parse()
                .map_err(|_| bad(line, "invalid epoch flag"))?;
            let count: usize = tokens[7]
                .parse()
                .map_err(|_| bad(line, "invalid epoch record count"))?;
            if flag > 6 || count > 999 {
                return Err(bad(line, "invalid epoch flag or record count"));
            }
            if (2..=5).contains(&flag) {
                let start = self.input.line;
                while self.input.line - start < count {
                    let (event_line, event) = self.input.required()?;
                    if flag == 4 {
                        apply_header_line(&mut self.input, &mut self.header, event_line, &event)?;
                    }
                }
                if self.input.line - start != count {
                    return Err(bad(line, "header event count mismatch"));
                }
                continue;
            }
            let time = calendar(&tokens[..6], line)?
                .to_gnss_time(self.header.time_scale)
                .map_err(|e| bad(line, e.to_string()))?;
            let receiver_clock_offset_s = tokens.get(8).map(|s| numeric(s, line)).transpose()?;
            let mut observations = Vec::new();
            let mut seen = std::collections::BTreeSet::new();
            for _ in 0..count {
                let (sat_line, sat_text) = self.input.required()?;
                let sat: SatelliteId = field(&sat_text, 0, 3)
                    .parse()
                    .map_err(|_| bad(sat_line, "invalid satellite record"))?;
                if !seen.insert(sat) {
                    return Err(bad(sat_line, "duplicate satellite in epoch"));
                }
                let types = self
                    .header
                    .observation_types
                    .get(&sat.system())
                    .ok_or_else(|| bad(sat_line, "undeclared observation types"))?;
                let mut fields = observation_fields(field(&sat_text, 3, sat_text.len()), sat_line)?;
                while fields.len() < types.len() {
                    match self.input.next_line()? {
                        Some((continuation_line, continuation))
                            if continuation.starts_with("   ") =>
                        {
                            let extra = observation_fields(&continuation[3..], continuation_line)?;
                            if extra.is_empty() {
                                return Err(bad(
                                    continuation_line,
                                    "empty observation continuation",
                                ));
                            }
                            fields.extend(extra);
                        }
                        Some(next) => {
                            self.input.pending = Some(next);
                            break;
                        }
                        None => break,
                    }
                }
                if fields.len() > types.len() {
                    if fields[types.len()..].iter().any(|x| !x.trim().is_empty()) {
                        return Err(bad(sat_line, "too many observation fields"));
                    }
                    fields.truncate(types.len());
                }
                fields.resize(types.len(), String::new());
                let mut groups: BTreeMap<String, Observation> = BTreeMap::new();
                for (code, value) in types.iter().zip(&fields) {
                    let tracking = &code[1..];
                    let obs = groups
                        .entry(tracking.to_owned())
                        .or_insert_with(|| Observation::new(sat, tracking.to_owned()));
                    obs.glonass_channel = self
                        .header
                        .glonass_channels
                        .iter()
                        .rev()
                        .find(|&&(s, _)| s == sat)
                        .map(|&(_, c)| c);
                    let measurement = optional_numeric(field(value, 0, 14), sat_line)?;
                    let indicator = |index| match value.as_bytes().get(index).copied() {
                        None | Some(b' ') => Ok(None),
                        Some(b @ b'0'..=b'9') => Ok(Some(b - b'0')),
                        _ => Err(bad(sat_line, "invalid observation indicator")),
                    };
                    match code.as_bytes()[0] {
                        b'C' => obs.pseudorange_m = measurement,
                        b'L' => {
                            obs.carrier_phase_cycles = measurement;
                            obs.lli = indicator(14)?;
                            if obs.lli.is_some_and(|lli| lli > 7) {
                                return Err(bad(sat_line, "invalid LLI bits"));
                            }
                            obs.signal_strength = indicator(15)?;
                        }
                        b'D' => obs.doppler_hz = measurement,
                        b'S' => obs.snr_db_hz = measurement,
                        _ => unreachable!(),
                    }
                }
                observations.extend(groups.into_values());
            }
            return Ok(Some(ObservationEpoch {
                time,
                flag,
                approximate_position: self.header.approximate_position,
                receiver_clock_offset_s,
                observations,
            }));
        }
    }
}
impl<R: BufRead> Iterator for RinexObservationReader<R> {
    type Item = Result<ObservationEpoch, RinexError>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        match self.next_epoch() {
            Ok(Some(epoch)) => Some(Ok(epoch)),
            Ok(None) => {
                self.done = true;
                None
            }
            Err(e) => {
                self.done = true;
                Some(Err(e))
            }
        }
    }
}
impl<R: BufRead> std::iter::FusedIterator for RinexObservationReader<R> {}

fn observation_fields(payload: &str, line: usize) -> Result<Vec<String>, RinexError> {
    let mut fields = Vec::new();
    let mut start = 0;
    while start < payload.len() {
        let chunk = field(payload, start, 16);
        if chunk.len() < 14 && !chunk.trim().is_empty() {
            return Err(bad(line, "truncated observation field"));
        }
        if chunk.len() < 14 && chunk.trim().is_empty() {
            break;
        }
        fields.push(chunk.to_owned());
        start += 16;
    }
    Ok(fields)
}

/// Read RINEX 3 Kepler navigation. GLONASS/SBAS state-vector records and RINEX
/// 2/4 are deliberately rejected until their models and parsers are implemented.
pub fn read_navigation<R: BufRead>(reader: R) -> Result<NavigationData, RinexError> {
    let mut input = Input {
        reader,
        line: 0,
        pending: None,
    };
    let header = read_header(&mut input, b'N')?;
    let mut navigation = NavigationData::default();
    navigation.gps_ionosphere = header.gps_ionosphere_alpha.zip(header.gps_ionosphere_beta);
    while let Some((line, text)) = input.next_line()? {
        if text.trim().is_empty() {
            continue;
        }
        let satellite: SatelliteId = field(&text, 0, 3)
            .parse()
            .map_err(|_| bad(line, "invalid navigation satellite"))?;
        if !matches!(
            satellite.system(),
            GnssSystem::Gps | GnssSystem::Galileo | GnssSystem::Qzss | GnssSystem::BeiDou
        ) {
            return Err(unsupported(
                line,
                format!("navigation for {:?}", satellite.system()),
            ));
        }
        let scale = if satellite.system() == GnssSystem::BeiDou {
            TimeScale::BeiDou
        } else {
            TimeScale::Gps
        };
        let toc = calendar(
            &field(&text, 3, 20).split_whitespace().collect::<Vec<_>>(),
            line,
        )?
        .to_gnss_time(scale)
        .map_err(|e| bad(line, e.to_string()))?;
        let af0 = numeric(field(&text, 23, 19), line)?;
        let af1 = numeric(field(&text, 42, 19), line)?;
        let af2 = numeric(field(&text, 61, 19), line)?;
        let mut rows = [[None; 4]; 7];
        for row in &mut rows {
            let (row_line, record) = input.required()?;
            if !record.starts_with("    ") {
                return Err(bad(row_line, "invalid navigation continuation"));
            }
            for (i, value) in row.iter_mut().enumerate() {
                *value = optional_numeric(field(&record, 4 + 19 * i, 19), row_line)?;
            }
        }
        let required = |r: usize, c: usize| {
            rows[r][c].ok_or_else(|| bad(line + 1 + r, "missing navigation field"))
        };
        let week = required(4, 2)?;
        if week.fract() != 0.0 || week < 0.0 || week > f64::from(i32::MAX - 1356) {
            return Err(bad(line + 5, "invalid broadcast week"));
        }
        let toes = required(2, 0)?;
        let toe = if satellite.system() == GnssSystem::BeiDou {
            crate::GnssTime::new(week as i32 + 1356, toes + 14.0)
        } else {
            crate::GnssTime::new(week as i32, toes)
        }
        .map_err(|e| bad(line, e.to_string()))?;
        let orbit = KeplerOrbit {
            crs: required(0, 1)?,
            delta_n: required(0, 2)?,
            mean_anomaly: required(0, 3)?,
            cuc: required(1, 0)?,
            eccentricity: required(1, 1)?,
            cus: required(1, 2)?,
            sqrt_a: required(1, 3)?,
            cic: required(2, 1)?,
            omega0: required(2, 2)?,
            cis: required(2, 3)?,
            inclination: required(3, 0)?,
            crc: required(3, 1)?,
            argument_of_perigee: required(3, 2)?,
            omega_dot: required(3, 3)?,
            idot: required(4, 0)?,
        };
        let clock = BroadcastClock {
            af0,
            af1,
            af2,
            tgd: required(5, 2)?,
            secondary_tgd: if matches!(satellite.system(), GnssSystem::Galileo | GnssSystem::BeiDou)
            {
                rows[5][3]
            } else {
                None
            },
        };
        let eph = BroadcastEphemeris {
            satellite,
            toe,
            toc,
            toes,
            orbit,
            clock,
            health: unsigned_field(required(5, 1)?, line + 6)?,
            iode: unsigned_field(required(0, 0)?, line + 1)?,
            iodc: if matches!(satellite.system(), GnssSystem::Gps | GnssSystem::Qzss) {
                unsigned_field(required(5, 3)?, line + 6)?
            } else if satellite.system() == GnssSystem::Galileo {
                unsigned_field(required(0, 0)?, line + 1)?
            } else {
                unsigned_field(rows[6][1].unwrap_or(0.0), line + 7)?
            },
            data_source: unsigned_field(rows[4][1].unwrap_or(0.0), line + 5)?,
            accuracy_m: rows[5][0],
        };
        navigation
            .add_ephemeris(eph)
            .map_err(|e| bad(line, e.to_string()))?;
    }
    Ok(navigation)
}
