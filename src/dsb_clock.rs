//! Relative GPS code DSBs under an explicit IF(C1W,C2W) clock constraint.
//! This is a caller-declared gauge, not an absolute OSB product or proof of
//! compatibility with a particular precise clock. No phase bias is invented.
use crate::bias::{BiasEntry, BiasError, BiasKind, BiasProducts, BiasSignalRequest};
use crate::dual_frequency::{
    DualFrequencyError, SignalCorrection, SignalCorrections, ionosphere_free_coefficients,
};
use crate::precise::ClockSource;
use crate::{GnssSystem, GnssTime, SatelliteId, constants as c};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug)]
pub enum DsbClockError {
    InvalidBinding,
    NotRelativeProduct,
    UnsupportedRequest,
    MissingPair {
        satellite: SatelliteId,
        first: String,
        second: String,
    },
    ExpiredPair {
        satellite: SatelliteId,
        first: String,
        second: String,
    },
    AmbiguousPair {
        satellite: SatelliteId,
        first: String,
        second: String,
    },
    NonFinite,
    Bias(BiasError),
    Correction(DualFrequencyError),
}
impl std::fmt::Display for DsbClockError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "GPS IF DSB clock binding: {self:?}")
    }
}
impl std::error::Error for DsbClockError {}
impl From<BiasError> for DsbClockError {
    fn from(e: BiasError) -> Self {
        Self::Bias(e)
    }
}
impl From<DualFrequencyError> for DsbClockError {
    fn from(e: DualFrequencyError) -> Self {
        Self::Correction(e)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GpsIfClockBinding {
    pub product_id: String,
    pub clock_source: ClockSource,
    /// Caller asserts the clock uses a*b_C1W + b*b_C2W = 0, where [a,b]
    /// are the GPS L1/L2 ionosphere-free coefficients. Not inferred from BIA.
    pub clock_reference: String,
}
#[derive(Debug, Clone, PartialEq)]
pub struct DsbCodeTerm {
    /// Original row: direction, units, sigma, SVN and half-open validity.
    pub entry: BiasEntry,
    /// Additive correction in metres = sum(coefficient * entry.metres()).
    pub coefficient: f64,
}
#[derive(Debug, Clone, PartialEq)]
pub struct DerivedDsbCode {
    pub code_add_m: f64,
    pub source_terms: Vec<DsbCodeTerm>,
    /// sum(abs(coefficient)*source_sigma). A conservative standard-deviation
    /// bound permitting unknown correlations, not an independent-error RMS.
    /// It is evidence only and is not silently fused as measurement noise.
    pub sigma_upper_bound_m: f64,
}
#[derive(Debug, Clone)]
pub struct DsbClockSnapshot {
    pub time: GnssTime,
    pub binding: GpsIfClockBinding,
    pub corrections: SignalCorrections,
    pub provenance: BTreeMap<(SatelliteId, String), DerivedDsbCode>,
}

pub struct GpsIfDsbClock<'a> {
    products: &'a BiasProducts,
    pairs: BTreeMap<(SatelliteId, String, String), Vec<usize>>,
}
impl<'a> GpsIfDsbClock<'a> {
    pub fn new(products: &'a BiasProducts) -> Result<Self, DsbClockError> {
        if products.header().split_whitespace().nth(7) != Some("R")
            || !products
                .description()
                .get("BIAS_MODE")
                .is_some_and(|m| m.len() == 1 && m[0].as_slice() == ["RELATIVE"])
        {
            return Err(DsbClockError::NotRelativeProduct);
        }
        let mut pairs: BTreeMap<_, Vec<usize>> = BTreeMap::new();
        for (i, e) in products.entries().iter().enumerate() {
            if e.kind == BiasKind::Dsb
                && e.station.is_none()
                && e.satellite.system() == GnssSystem::Gps
                && e.observation_1.starts_with('C')
                && e.observation_2.as_ref().is_some_and(|s| s.starts_with('C'))
            {
                pairs
                    .entry((
                        e.satellite,
                        e.observation_1.clone(),
                        e.observation_2.clone().unwrap(),
                    ))
                    .or_default()
                    .push(i);
            }
        }
        Ok(Self { products, pairs })
    }

    fn pair(
        &self,
        sat: SatelliteId,
        first: &str,
        second: &str,
        time: GnssTime,
    ) -> Result<(&'a BiasEntry, f64), DsbClockError> {
        let mut found = Vec::new();
        let mut exists = false;
        for (a, b, sign) in [(first, second, 1.0), (second, first, -1.0)] {
            if let Some(indices) = self.pairs.get(&(sat, a.to_owned(), b.to_owned())) {
                exists = true;
                found.extend(
                    indices
                        .iter()
                        .map(|&i| &self.products.entries()[i])
                        .filter(|e| e.contains(time))
                        .map(|e| (e, sign)),
                );
            }
        }
        match found.len() {
            1 => Ok(found[0]),
            0 if exists => Err(DsbClockError::ExpiredPair {
                satellite: sat,
                first: first.into(),
                second: second.into(),
            }),
            0 => Err(DsbClockError::MissingPair {
                satellite: sat,
                first: first.into(),
                second: second.into(),
            }),
            _ => Err(DsbClockError::AmbiguousPair {
                satellite: sat,
                first: first.into(),
                second: second.into(),
            }),
        }
    }

    /// Materialize only the requested epoch, preserving BIA's half-open end.
    /// Uses exactly C1W-C2W plus C1C-C1W or C2W-C2L where required. It never
    /// silently substitutes an alternate graph path or fits inconsistent edges.
    pub fn corrections_at(
        &self,
        time: GnssTime,
        binding: &GpsIfClockBinding,
        requests: &[BiasSignalRequest],
    ) -> Result<DsbClockSnapshot, DsbClockError> {
        if binding.product_id.trim().is_empty() || binding.clock_reference.trim().is_empty() {
            return Err(DsbClockError::InvalidBinding);
        }
        let [a, b] = ionosphere_free_coefficients(c::GPS_L1_FREQ, c::GPS_L2_FREQ)?;
        let mut keys = BTreeSet::new();
        let mut entries = Vec::new();
        let mut provenance = BTreeMap::new();
        for request in requests {
            let sat = request.satellite;
            let tracking = request.tracking_code.as_str();
            if sat.system() != GnssSystem::Gps
                || !matches!(tracking, "1C" | "1W" | "2W" | "2L")
                || request.glonass_channel.is_some()
            {
                return Err(DsbClockError::UnsupportedRequest);
            }
            if !keys.insert((sat, tracking)) {
                return Err(DsbClockError::InvalidBinding);
            }
            let (datum, sign) = self.pair(sat, "C1W", "C2W", time)?;
            let mut terms = vec![DsbCodeTerm {
                entry: datum.clone(),
                coefficient: sign * if tracking.starts_with('1') { -b } else { a },
            }];
            let extra = match tracking {
                "1C" => Some(("C1C", "C1W", -1.0)),
                "2L" => Some(("C2W", "C2L", 1.0)),
                _ => None,
            };
            if let Some((first, second, coefficient)) = extra {
                let (entry, sign) = self.pair(sat, first, second, time)?;
                terms.push(DsbCodeTerm {
                    entry: entry.clone(),
                    coefficient: coefficient * sign,
                });
            }
            let (mut correction, mut bound) = (0.0, 0.0);
            for t in &terms {
                let (value, sigma) = t.entry.metres(None)?;
                correction += t.coefficient * value;
                bound += t.coefficient.abs() * sigma;
            }
            if !correction.is_finite() || !bound.is_finite() {
                return Err(DsbClockError::NonFinite);
            }
            entries.push(SignalCorrection {
                satellite: sat,
                tracking_code: tracking.into(),
                code_add_m: correction,
                phase_add_m: None,
                phase_reference: None,
                clock_source: binding.clock_source,
                clock_reference: binding.clock_reference.clone(),
                valid_from: time,
                valid_until: time,
            });
            provenance.insert(
                (sat, tracking.to_owned()),
                DerivedDsbCode {
                    code_add_m: correction,
                    source_terms: terms,
                    sigma_upper_bound_m: bound,
                },
            );
        }
        Ok(DsbClockSnapshot {
            time,
            binding: binding.clone(),
            corrections: SignalCorrections::new(entries)?,
            provenance,
        })
    }
}
