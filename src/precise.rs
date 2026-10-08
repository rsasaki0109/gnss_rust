//! SP3-c/d position products and RINEX 3 satellite clocks, in SI units.
//!
//! Neville interpolation and central-difference rates follow libgnss++
//! `precise_products.cpp`. Copyright (c) 2024 LibGNSS++ Contributors.
//! See LICENSE. Strict coverage/quality admission is intentional; this module
//! does not implement PPP or apply the periodic relativistic clock correction.

use crate::time::{CalendarDateTime, TimeScale};
use crate::{EcefCoord, GnssTime, SatelliteId, constants::SPEED_OF_LIGHT};
use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
    io::BufRead,
};

#[derive(Debug)]
pub enum ProductError {
    Io(std::io::Error),
    Parse { line: usize, message: String },
    Unsupported { line: usize, feature: String },
}
impl fmt::Display for ProductError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "precise product I/O: {e}"),
            Self::Parse { line, message } => write!(f, "precise product line {line}: {message}"),
            Self::Unsupported { line, feature } => {
                write!(f, "precise product line {line}: unsupported {feature}")
            }
        }
    }
}
impl std::error::Error for ProductError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            _ => None,
        }
    }
}
impl From<std::io::Error> for ProductError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}
fn bad(line: usize, message: impl Into<String>) -> ProductError {
    ProductError::Parse {
        line,
        message: message.into(),
    }
}
fn unsupported(line: usize, feature: impl Into<String>) -> ProductError {
    ProductError::Unsupported {
        line,
        feature: feature.into(),
    }
}
struct Input<R> {
    reader: R,
    line: usize,
}
impl<R: BufRead> Input<R> {
    fn next(&mut self) -> Result<Option<String>, ProductError> {
        let mut text = String::new();
        if self.reader.read_line(&mut text)? == 0 {
            return Ok(None);
        }
        self.line += 1;
        if !text.is_ascii() {
            return Err(bad(self.line, "records must be ASCII"));
        }
        Ok(Some(text.trim_end_matches(['\r', '\n']).to_owned()))
    }
    fn required(&mut self) -> Result<String, ProductError> {
        self.next()?
            .ok_or_else(|| bad(self.line + 1, "unexpected end of file"))
    }
}
fn field(text: &str, start: usize, width: usize) -> &str {
    text.get(start..text.len().min(start + width)).unwrap_or("")
}
fn number(text: &str, line: usize) -> Result<f64, ProductError> {
    let value: f64 = text
        .trim()
        .replace(['D', 'd'], "E")
        .parse()
        .map_err(|_| bad(line, "invalid number"))?;
    if !value.is_finite() {
        return Err(bad(line, "nonfinite number"));
    }
    Ok(value)
}
fn integer<T: std::str::FromStr>(text: &str, line: usize) -> Result<T, ProductError> {
    text.trim()
        .parse()
        .map_err(|_| bad(line, "invalid integer"))
}
fn satellite(text: &str, line: usize) -> Result<SatelliteId, ProductError> {
    text.parse()
        .map_err(|_| bad(line, format!("invalid satellite {text:?}")))
}
fn time_scale(text: &str, line: usize) -> Result<TimeScale, ProductError> {
    match text.trim() {
        "GPS" => Ok(TimeScale::Gps),
        "GAL" => Ok(TimeScale::Galileo),
        "QZS" => Ok(TimeScale::Qzss),
        "BDT" => Ok(TimeScale::BeiDou),
        "UTC" => Ok(TimeScale::Utc),
        value => Err(unsupported(line, format!("time system {value:?}"))),
    }
}
fn calendar(v: &[&str], scale: TimeScale, line: usize) -> Result<GnssTime, ProductError> {
    if v.len() != 6 {
        return Err(bad(line, "expected six calendar fields"));
    }
    CalendarDateTime {
        year: integer(v[0], line)?,
        month: integer(v[1], line)?,
        day: integer(v[2], line)?,
        hour: integer(v[3], line)?,
        minute: integer(v[4], line)?,
        second: number(v[5], line)?,
    }
    .to_gnss_time(scale)
    .map_err(|e| bad(line, format!("invalid calendar: {e}")))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrbitReferencePoint {
    CentreOfMass,
    AntennaPhaseCenter,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Sp3Header {
    pub version: char,
    pub start_time: GnssTime,
    pub declared_epochs: usize,
    pub interval_s: f64,
    pub coordinate_system: String,
    pub data_used: String,
    pub orbit_type: String,
    pub agency: String,
    pub time_scale: TimeScale,
    pub satellites: Vec<SatelliteId>,
    /// Per-satellite header accuracy exponents, in satellite-list order.
    pub accuracy_exponents: Vec<Option<u16>>,
    pub position_sigma_base: f64,
    pub clock_sigma_base: f64,
    pub reference_point: OrbitReferencePoint,
    /// Retains uninterpreted metadata, comments, and declared accuracy bases.
    pub raw_lines: Vec<String>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Sp3Flags {
    pub clock_event: bool,
    pub clock_predicted: bool,
    pub orbit_maneuver: bool,
    pub orbit_predicted: bool,
}
#[derive(Debug, Clone, PartialEq)]
pub struct Sp3Sample {
    pub satellite: SatelliteId,
    pub time: GnssTime,
    /// SP3 zero coordinate components denote missing orbit data.
    pub position_m: Option<EcefCoord>,
    /// Missing SP3 clocks (999999.999999 microseconds) remain None; zero is valid.
    pub clock_bias_s: Option<f64>,
    /// Raw record accuracy exponents (x, y, z, clock). Blank means unknown.
    pub accuracy_exponents: [Option<u16>; 4],
    pub flags: Sp3Flags,
}
#[derive(Debug, Clone)]
pub struct Sp3Data {
    header: Sp3Header,
    samples: BTreeMap<SatelliteId, Vec<Sp3Sample>>,
}
impl Sp3Data {
    pub fn header(&self) -> &Sp3Header {
        &self.header
    }
    pub fn samples(&self) -> &BTreeMap<SatelliteId, Vec<Sp3Sample>> {
        &self.samples
    }
}

/// Read complete SP3-c/d P products. V, EP and EV records are explicitly
/// unsupported. Epochs and satellite declarations must be complete/consistent.
pub fn read_sp3<R: BufRead>(reader: R) -> Result<Sp3Data, ProductError> {
    let mut input = Input { reader, line: 0 };
    let first = input.required()?;
    if !first.starts_with("#c") && !first.starts_with("#d") {
        return Err(unsupported(1, "SP3 version (requires c or d header)"));
    }
    if first.as_bytes().get(2) != Some(&b'P') {
        return Err(unsupported(1, "SP3 velocity product"));
    }
    let epoch_fields: Vec<_> = field(&first, 3, 28).split_whitespace().collect();
    let epochs = integer::<usize>(field(&first, 32, 7), 1)?;
    if epochs == 0 {
        return Err(bad(1, "zero declared epochs"));
    }
    let frame = field(&first, 46, 5).trim().to_owned();
    if frame.is_empty() {
        return Err(bad(1, "missing coordinate system"));
    }
    let mut raw_lines = vec![first.clone()];
    let second = input.required()?;
    if !second.starts_with("##") {
        return Err(bad(2, "missing SP3 interval record"));
    }
    let interval = number(field(&second, 24, 14), 2)?;
    if interval <= 0.0 {
        return Err(bad(2, "SP3 interval must be positive"));
    }
    raw_lines.push(second);
    let (mut declared_satellites, mut satellites, mut accuracy) = (None, Vec::new(), Vec::new());
    let (mut scale, mut sigma_bases) = (None, None);
    let mut reference_point = OrbitReferencePoint::CentreOfMass;
    let mut pending;
    loop {
        let text = input.required()?;
        if text.starts_with('*') {
            pending = Some(text);
            break;
        }
        let line = input.line;
        if text.starts_with("++") {
            for slot in (9..60).step_by(3) {
                let value = field(&text, slot, 3).trim();
                accuracy.push(if value.is_empty() {
                    None
                } else {
                    Some(integer(value, line)?)
                });
            }
        } else if text.starts_with('+') {
            let count = field(&text, 3, 3).trim();
            if !count.is_empty() {
                let count: usize = integer(count, line)?;
                if declared_satellites.replace(count).is_some() {
                    return Err(bad(line, "duplicate satellite count"));
                }
            }
            for slot in (9..60).step_by(3) {
                let value = field(&text, slot, 3).trim();
                if value.is_empty() || value == "0" {
                    continue;
                }
                satellites.push(satellite(value, line)?);
            }
        } else if text.starts_with("%c") {
            if scale.is_none() {
                scale = Some(time_scale(field(&text, 9, 3), line)?);
            }
        } else if let Some(content) = text.strip_prefix("%f") {
            if sigma_bases.is_none() {
                let values: Vec<_> = content.split_whitespace().collect();
                if values.len() < 2 {
                    return Err(bad(line, "missing accuracy bases"));
                }
                let bases = (number(values[0], line)?, number(values[1], line)?);
                if bases.0 <= 0.0 || bases.1 <= 0.0 {
                    return Err(bad(line, "invalid accuracy bases"));
                }
                sigma_bases = Some(bases);
            }
        } else if text.starts_with("/*") {
            if text.contains("ANTENNA PHASE CENTER") {
                reference_point = OrbitReferencePoint::AntennaPhaseCenter;
            }
        } else if !text.starts_with("%i") {
            return Err(unsupported(line, "SP3 header record"));
        }
        raw_lines.push(text);
    }
    if declared_satellites != Some(satellites.len()) || satellites.is_empty() {
        return Err(bad(input.line, "satellite count mismatch"));
    }
    let declared: BTreeSet<_> = satellites.iter().copied().collect();
    if declared.len() != satellites.len() {
        return Err(bad(input.line, "duplicate satellite declaration"));
    }
    if accuracy.len() < satellites.len() {
        return Err(bad(input.line, "missing satellite accuracy declarations"));
    }
    accuracy.truncate(satellites.len());
    let scale = scale.ok_or_else(|| bad(input.line, "missing SP3 time system"))?;
    let bases = sigma_bases.ok_or_else(|| bad(input.line, "missing SP3 accuracy bases"))?;
    let header = Sp3Header {
        version: first.as_bytes()[1] as char,
        start_time: calendar(&epoch_fields, scale, 1)?,
        declared_epochs: epochs,
        interval_s: interval,
        coordinate_system: frame,
        data_used: field(&first, 40, 5).trim().to_owned(),
        orbit_type: field(&first, 52, 3).trim().to_owned(),
        agency: field(&first, 56, 4).trim().to_owned(),
        time_scale: scale,
        satellites,
        accuracy_exponents: accuracy,
        position_sigma_base: bases.0,
        clock_sigma_base: bases.1,
        reference_point,
        raw_lines,
    };
    let mut samples: BTreeMap<SatelliteId, Vec<Sp3Sample>> = BTreeMap::new();
    let (mut current, mut seen, mut epoch_count) = (None, BTreeSet::new(), 0);
    let mut eof = false;
    while let Some(text) = match pending.take() {
        Some(t) => Some(t),
        None => input.next()?,
    } {
        let line = input.line;
        if let Some(content) = text.strip_prefix('*') {
            if current.is_some() && seen != declared {
                return Err(bad(line, "incomplete satellite epoch"));
            }
            let fields: Vec<_> = content.split_whitespace().collect();
            let time = calendar(&fields, scale, line)?;
            if let Some(previous) = current {
                let elapsed = time.difference_seconds(previous);
                if elapsed <= 0.0 || (elapsed - interval).abs() > 1e-7 {
                    return Err(bad(
                        line,
                        "epochs must increase by the declared SP3 interval",
                    ));
                }
            } else if time != header.start_time {
                return Err(bad(line, "first epoch differs from header"));
            }
            current = Some(time);
            seen.clear();
            epoch_count += 1;
        } else if text.starts_with('P') {
            let time = current.ok_or_else(|| bad(line, "position before epoch"))?;
            if text.len() < 60 {
                return Err(bad(line, "short SP3 position record"));
            }
            let sat = satellite(field(&text, 1, 3), line)?;
            if !declared.contains(&sat) || !seen.insert(sat) {
                return Err(bad(line, "undeclared or duplicate satellite"));
            }
            let x = number(field(&text, 4, 14), line)?;
            let y = number(field(&text, 18, 14), line)?;
            let z = number(field(&text, 32, 14), line)?;
            let clock = number(field(&text, 46, 14), line)?;
            let mut exponents = [None; 4];
            for (i, (start, width)) in [(61, 2), (64, 2), (67, 2), (70, 3)].into_iter().enumerate()
            {
                let value = field(&text, start, width).trim();
                if !value.is_empty() {
                    exponents[i] = Some(integer(value, line)?);
                }
            }
            let flag = |index, expected| -> Result<bool, ProductError> {
                match text.as_bytes().get(index).copied().unwrap_or(b' ') {
                    b' ' => Ok(false),
                    actual if actual == expected => Ok(true),
                    _ => Err(bad(line, "invalid SP3 quality flag")),
                }
            };
            let flags = Sp3Flags {
                clock_event: flag(74, b'E')?,
                clock_predicted: flag(75, b'P')?,
                orbit_maneuver: flag(78, b'M')?,
                orbit_predicted: flag(79, b'P')?,
            };
            let metres = [x * 1000.0, y * 1000.0, z * 1000.0];
            if !metres.iter().all(|v| v.is_finite()) {
                return Err(bad(line, "position unit conversion overflow"));
            }
            samples.entry(sat).or_default().push(Sp3Sample {
                satellite: sat,
                time,
                position_m: if x == 0.0 || y == 0.0 || z == 0.0 {
                    None
                } else {
                    Some(EcefCoord::new(metres[0], metres[1], metres[2]))
                },
                clock_bias_s: if clock.abs() >= 999999.0 {
                    None
                } else {
                    Some(clock * 1e-6)
                },
                accuracy_exponents: exponents,
                flags,
            });
        } else if text.trim() == "EOF" {
            eof = true;
            break;
        } else {
            return Err(unsupported(line, "SP3 data record (requires P)"));
        }
    }
    if !eof || epoch_count != epochs || seen != declared {
        return Err(bad(input.line, "missing EOF or incomplete SP3 epochs"));
    }
    while let Some(text) = input.next()? {
        if !text.trim().is_empty() {
            return Err(bad(input.line, "data after EOF"));
        }
    }
    Ok(Sp3Data { header, samples })
}

#[derive(Debug, Clone, PartialEq)]
pub struct ClockHeader {
    pub version: f64,
    pub time_scale: TimeScale,
    /// No TIME SYSTEM ID means the RINEX clock specification's GPS default.
    pub time_system_declared: bool,
    pub raw_lines: Vec<String>,
}
#[derive(Debug, Clone, PartialEq)]
pub struct ClockSample {
    pub satellite: SatelliteId,
    pub time: GnssTime,
    pub bias_s: f64,
    pub sigma_s: Option<f64>,
    pub drift_s_per_s: Option<f64>,
    pub drift_sigma_s_per_s: Option<f64>,
    pub rate_s_per_s2: Option<f64>,
    pub rate_sigma_s_per_s2: Option<f64>,
}
#[derive(Debug, Clone)]
pub struct ClockData {
    header: ClockHeader,
    samples: BTreeMap<SatelliteId, Vec<ClockSample>>,
}
impl ClockData {
    pub fn header(&self) -> &ClockHeader {
        &self.header
    }
    pub fn samples(&self) -> &BTreeMap<SatelliteId, Vec<ClockSample>> {
        &self.samples
    }
}
/// Read RINEX 3.00–3.04 satellite AS clocks, including 1–6 values and
/// continuation lines. Other clock record types fail explicitly.
pub fn read_clk<R: BufRead>(reader: R) -> Result<ClockData, ProductError> {
    let mut input = Input { reader, line: 0 };
    let first = input.required()?;
    if first.get(60..).unwrap_or("").trim() != "RINEX VERSION / TYPE" {
        return Err(bad(1, "missing RINEX VERSION / TYPE"));
    }
    let version = number(field(&first, 0, 9), 1)?;
    if !(3.0..=3.04).contains(&version) {
        return Err(unsupported(1, "RINEX clock version"));
    }
    if first.as_bytes().get(20) != Some(&b'C') {
        return Err(bad(1, "expected clock file type C"));
    }
    let mut header = ClockHeader {
        version,
        time_scale: TimeScale::Gps,
        time_system_declared: false,
        raw_lines: vec![first],
    };
    loop {
        let text = input.required()?;
        // RINEX CLK 3.04 uses labels starting at column 66, older versions 61.
        let start = if version >= 3.04 { 65 } else { 60 };
        let label = text
            .get(start..)
            .ok_or_else(|| bad(input.line, "short clock header"))?
            .trim();
        if label == "TIME SYSTEM ID" {
            if header.time_system_declared {
                return Err(bad(input.line, "duplicate clock time system"));
            }
            header.time_scale = time_scale(text[..start].trim(), input.line)?;
            header.time_system_declared = true;
        }
        let end = label == "END OF HEADER";
        header.raw_lines.push(text);
        if end {
            break;
        }
    }
    let mut samples: BTreeMap<SatelliteId, Vec<ClockSample>> = BTreeMap::new();
    let mut last_epoch = None;
    while let Some(text) = input.next()? {
        if text.trim().is_empty() {
            continue;
        }
        let line = input.line;
        let tokens: Vec<_> = text.split_whitespace().collect();
        if tokens.first().copied() != Some("AS") {
            return Err(unsupported(line, "non-AS clock record"));
        }
        if tokens.len() < 10 {
            return Err(bad(line, "short clock record"));
        }
        let sat = satellite(tokens[1], line)?;
        let time = calendar(&tokens[2..8], header.time_scale, line)?;
        if last_epoch.is_some_and(|t| time < t) {
            return Err(bad(line, "clock epochs must not decrease"));
        }
        last_epoch = Some(time);
        let count: usize = integer(tokens[8], line)?;
        if !(1..=6).contains(&count) {
            return Err(unsupported(line, "clock value count (requires 1–6)"));
        }
        let mut values: Vec<_> = tokens[9..]
            .iter()
            .map(|v| number(v, line))
            .collect::<Result<_, _>>()?;
        while values.len() < count {
            let continuation = input.required()?;
            if !continuation.starts_with("   ") {
                return Err(bad(input.line, "expected clock value continuation"));
            }
            for value in continuation.split_whitespace() {
                values.push(number(value, input.line)?);
            }
        }
        if values.len() != count {
            return Err(bad(line, "clock value count mismatch"));
        }
        for index in [1, 3, 5] {
            if values.get(index).is_some_and(|v| *v < 0.0) {
                return Err(bad(line, "negative clock sigma"));
            }
        }
        let rows = samples.entry(sat).or_default();
        if rows.last().is_some_and(|s| s.time == time) {
            return Err(bad(line, "duplicate satellite clock epoch"));
        }
        rows.push(ClockSample {
            satellite: sat,
            time,
            bias_s: values[0],
            sigma_s: values.get(1).copied(),
            drift_s_per_s: values.get(2).copied(),
            drift_sigma_s_per_s: values.get(3).copied(),
            rate_s_per_s2: values.get(4).copied(),
            rate_sigma_s_per_s2: values.get(5).copied(),
        });
    }
    if samples.is_empty() {
        return Err(bad(input.line, "no AS clocks"));
    }
    Ok(ClockData { header, samples })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InterpolationError {
    InvalidConfiguration,
    MissingSatellite,
    InsufficientSamples,
    OutOfRange,
    GapOrDiscontinuity,
    NumericalFailure,
}
impl fmt::Display for InterpolationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidConfiguration => "invalid precise interpolation configuration",
            Self::MissingSatellite => "satellite absent from precise product",
            Self::InsufficientSamples => "at least two valid precise samples are required",
            Self::OutOfRange => "query outside precise product coverage",
            Self::GapOrDiscontinuity => {
                "precise samples unavailable across a gap or quality boundary"
            }
            Self::NumericalFailure => "precise interpolation arithmetic failed",
        })
    }
}
impl std::error::Error for InterpolationError {}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterpolationConfig {
    pub max_samples: usize,
    pub max_orbit_gap_s: f64,
    pub max_clock_gap_s: f64,
    /// Bounded extrapolation at the product's outer ends only; default zero.
    pub max_extrapolation_s: f64,
    pub allow_predicted: bool,
}
impl Default for InterpolationConfig {
    fn default() -> Self {
        Self {
            max_samples: 10,
            max_orbit_gap_s: 900.0,
            max_clock_gap_s: 900.0,
            max_extrapolation_s: 0.0,
            allow_predicted: false,
        }
    }
}
impl InterpolationConfig {
    pub(crate) fn validate(self) -> Result<(), InterpolationError> {
        if !(2..=10).contains(&self.max_samples)
            || ![self.max_orbit_gap_s, self.max_clock_gap_s]
                .iter()
                .all(|v| v.is_finite() && *v > 0.0)
            || !self.max_extrapolation_s.is_finite()
            || !(0.0..=900.0).contains(&self.max_extrapolation_s)
        {
            return Err(InterpolationError::InvalidConfiguration);
        }
        Ok(())
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct InterpolationSupport {
    pub sample_times: Vec<GnssTime>,
    pub extrapolated: bool,
    pub includes_predicted: bool,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ClockSource {
    Sp3,
    RinexClk,
}
#[derive(Debug, Clone, PartialEq)]
pub struct InterpolatedClock {
    /// Published clock, without periodic relativity or signal-bias corrections.
    pub bias_s: f64,
    pub drift_s_per_s: f64,
    pub source: ClockSource,
    pub support: InterpolationSupport,
}
#[derive(Debug, Clone, PartialEq)]
pub struct PreciseState {
    pub satellite: SatelliteId,
    pub time: GnssTime,
    pub position_m: EcefCoord,
    pub velocity_m_per_s: [f64; 3],
    /// Selected for this satellite even when the clock is unavailable.
    pub clock_source: ClockSource,
    pub clock: Option<InterpolatedClock>,
    pub clock_unavailability: Option<InterpolationError>,
    pub orbit_support: InterpolationSupport,
    pub coordinate_system: String,
    pub reference_point: OrbitReferencePoint,
}
#[derive(Debug, Clone)]
pub struct PreciseProducts {
    orbit: Sp3Data,
    clocks: Option<ClockData>,
}
impl PreciseProducts {
    pub fn new(orbit: Sp3Data, clocks: Option<ClockData>) -> Self {
        Self { orbit, clocks }
    }
    pub fn orbit(&self) -> &Sp3Data {
        &self.orbit
    }
    pub fn clocks(&self) -> Option<&ClockData> {
        self.clocks.as_ref()
    }
    pub fn interpolate(
        &self,
        satellite: SatelliteId,
        time: GnssTime,
        config: InterpolationConfig,
    ) -> Result<PreciseState, InterpolationError> {
        config.validate()?;
        let rows = self
            .orbit
            .samples
            .get(&satellite)
            .ok_or(InterpolationError::MissingSatellite)?;
        let orbit: Vec<_> = rows
            .iter()
            .map(|s| Point {
                time: s.time,
                values: s.position_m.map(|p| [p.x, p.y, p.z]),
                predicted: s.flags.orbit_predicted,
                boundary: s.flags.orbit_maneuver,
            })
            .collect();
        let (position, velocity, support) =
            interpolate_points(&orbit, time, config, config.max_orbit_gap_s)?;
        // CLK is selected for the entire satellite arc, never blended with SP3
        // or replaced by SP3 in a CLK gap (product clock datums may differ).
        let (source, clock_points) =
            if let Some(rows) = self.clocks.as_ref().and_then(|c| c.samples.get(&satellite)) {
                (
                    ClockSource::RinexClk,
                    rows.iter()
                        .map(|s| Point {
                            time: s.time,
                            values: Some([s.bias_s]),
                            predicted: false,
                            boundary: false,
                        })
                        .collect::<Vec<_>>(),
                )
            } else {
                (
                    ClockSource::Sp3,
                    rows.iter()
                        .map(|s| Point {
                            time: s.time,
                            values: s.clock_bias_s.map(|v| [v]),
                            predicted: s.flags.clock_predicted,
                            boundary: s.flags.clock_event,
                        })
                        .collect::<Vec<_>>(),
                )
            };
        let (clock, clock_unavailability) =
            match interpolate_points(&clock_points, time, config, config.max_clock_gap_s) {
                Ok((bias, drift, support)) => (
                    Some(InterpolatedClock {
                        bias_s: bias[0],
                        drift_s_per_s: drift[0],
                        source,
                        support,
                    }),
                    None,
                ),
                Err(e) => (None, Some(e)),
            };
        Ok(PreciseState {
            satellite,
            time,
            position_m: EcefCoord::new(position[0], position[1], position[2]),
            velocity_m_per_s: velocity,
            clock_source: source,
            clock,
            clock_unavailability,
            orbit_support: support,
            coordinate_system: self.orbit.header.coordinate_system.clone(),
            reference_point: self.orbit.header.reference_point,
        })
    }
}
#[derive(Clone, Copy)]
struct Point<const N: usize> {
    time: GnssTime,
    values: Option<[f64; N]>,
    predicted: bool,
    boundary: bool,
}
fn interpolate_points<const N: usize>(
    points: &[Point<N>],
    query: GnssTime,
    config: InterpolationConfig,
    max_gap: f64,
) -> Result<([f64; N], [f64; N], InterpolationSupport), InterpolationError> {
    let first = points
        .first()
        .ok_or(InterpolationError::InsufficientSamples)?
        .time;
    let last = points.last().unwrap().time;
    let outside = query < first || query > last;
    if first.difference_seconds(query) > config.max_extrapolation_s
        || query.difference_seconds(last) > config.max_extrapolation_s
    {
        return Err(InterpolationError::OutOfRange);
    }
    let mut arcs: Vec<Vec<Point<N>>> = Vec::new();
    let mut arc: Vec<Point<N>> = Vec::new();
    for &p in points {
        let admitted = p.values.is_some() && (config.allow_predicted || !p.predicted);
        let split = !admitted
            || p.boundary
            || arc
                .last()
                .is_some_and(|a| p.time.difference_seconds(a.time) > max_gap);
        if split && !arc.is_empty() {
            arcs.push(std::mem::take(&mut arc));
        }
        if admitted {
            arc.push(p);
        }
    }
    if !arc.is_empty() {
        arcs.push(arc);
    }
    let arc = arcs
        .iter()
        .find(|a| {
            let start = a[0].time;
            let end = a.last().unwrap().time;
            (query >= start && query <= end)
                || (query < first && start == first)
                || (query > last && end == last)
        })
        .ok_or(InterpolationError::GapOrDiscontinuity)?;
    if arc.len() < 2 {
        return Err(InterpolationError::InsufficientSamples);
    }
    let n = arc.len().min(config.max_samples);
    let before = arc.partition_point(|p| p.time <= query).saturating_sub(1);
    let start = before.saturating_sub((n - 1) / 2).min(arc.len() - n);
    let window = &arc[start..start + n];
    let xs: Vec<_> = window
        .iter()
        .map(|p| p.time.difference_seconds(query))
        .collect();
    let mut value = [0.0; N];
    let mut rate = [0.0; N];
    for axis in 0..N {
        let ys: Vec<_> = window.iter().map(|p| p.values.unwrap()[axis]).collect();
        value[axis] = neville(&xs, ys.clone())?;
        let plus: Vec<_> = xs.iter().map(|x| x - 1.0).collect();
        let minus: Vec<_> = xs.iter().map(|x| x + 1.0).collect();
        rate[axis] = (neville(&plus, ys.clone())? - neville(&minus, ys)?) / 2.0;
        if !rate[axis].is_finite() {
            return Err(InterpolationError::NumericalFailure);
        }
    }
    Ok((
        value,
        rate,
        InterpolationSupport {
            sample_times: window.iter().map(|p| p.time).collect(),
            extrapolated: outside,
            includes_predicted: window.iter().any(|p| p.predicted),
        },
    ))
}
fn neville(xs: &[f64], mut ys: Vec<f64>) -> Result<f64, InterpolationError> {
    for j in 1..xs.len() {
        for i in 0..xs.len() - j {
            let denominator = xs[i + j] - xs[i];
            if !denominator.is_finite() || denominator.abs() < 1e-12 {
                return Err(InterpolationError::NumericalFailure);
            }
            ys[i] = (xs[i + j] * ys[i] - xs[i] * ys[i + 1]) / denominator;
            if !ys[i].is_finite() {
                return Err(InterpolationError::NumericalFailure);
            }
        }
    }
    Ok(ys[0])
}

/// Periodic relativistic correction in seconds, to be added once to published
/// precise clocks at transmission time. Does not apply antenna/signal biases.
pub fn precise_clock_relativistic_correction(
    position_m: EcefCoord,
    velocity_m_per_s: [f64; 3],
) -> Result<f64, InterpolationError> {
    let position = [position_m.x, position_m.y, position_m.z];
    if !position
        .iter()
        .chain(&velocity_m_per_s)
        .all(|v| v.is_finite())
    {
        return Err(InterpolationError::NumericalFailure);
    }
    let dot: f64 = position
        .iter()
        .zip(velocity_m_per_s)
        .map(|(p, v)| p * v)
        .sum();
    let result = -2.0 * dot / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    if !result.is_finite() {
        return Err(InterpolationError::NumericalFailure);
    }
    Ok(result)
}
