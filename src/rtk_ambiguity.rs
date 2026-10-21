//! Per-system SD ambiguity datum alignment, ported from libgnss++ updateBias.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
use crate::{Error, GnssSystem, SatelliteId};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DatumAlignmentInput {
    pub satellite: SatelliteId,
    /// SD carrier cycles minus SD code metres / L1/E1 wavelength.
    pub raw_bias_cycles: f64,
    /// None for a new, lost-lock, changed-tracking or globally reset arc.
    /// A retained value of zero is valid and contributes to the mean.
    pub retained_cycles: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DatumAlignmentDiagnostic {
    pub system: GnssSystem,
    pub retained_satellites: Vec<SatelliteId>,
    pub initialized_satellites: Vec<SatelliteId>,
    pub offset_cycles: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DatumAlignmentResult {
    /// Aligned values in the input order. New/reset arcs keep their raw bias;
    /// retained arcs gain their system's mean (raw bias - retained value).
    pub single_difference_cycles: Vec<f64>,
    pub systems: Vec<DatumAlignmentDiagnostic>,
}

/// Align persistent GPS L1 / Galileo E1 states to the current raw code/phase
/// datum before joining new arcs. Within-system differences between retained
/// states are unchanged. This translates means only, as the upstream does;
/// it neither estimates a clock nor treats the offset as a new measurement.
/// No covariance, outlier threshold or arc-continuity decision is changed.
pub fn align_single_difference_datums(
    inputs: &[DatumAlignmentInput],
) -> Result<DatumAlignmentResult, Error> {
    let mut seen = BTreeSet::new();
    let mut groups = BTreeMap::<GnssSystem, (f64, usize)>::new();
    for input in inputs {
        if !matches!(
            input.satellite.system(),
            GnssSystem::Gps | GnssSystem::Galileo
        ) {
            return Err(Error::UnsupportedConstellation);
        }
        if !seen.insert(input.satellite) {
            return Err(Error::InvalidConfiguration);
        }
        if !input.raw_bias_cycles.is_finite()
            || input.retained_cycles.is_some_and(|v| !v.is_finite())
        {
            return Err(Error::NonFinite);
        }
        let group = groups.entry(input.satellite.system()).or_default();
        if let Some(old) = input.retained_cycles {
            group.0 += input.raw_bias_cycles - old;
            group.1 += 1;
        }
    }
    let mut systems = Vec::with_capacity(groups.len());
    for (system, (sum, count)) in groups {
        let offset = if count == 0 { 0.0 } else { sum / count as f64 };
        if !offset.is_finite() {
            return Err(Error::NonFinite);
        }
        let members: Vec<_> = inputs
            .iter()
            .filter(|i| i.satellite.system() == system)
            .collect();
        systems.push(DatumAlignmentDiagnostic {
            system,
            retained_satellites: members
                .iter()
                .filter(|i| i.retained_cycles.is_some())
                .map(|i| i.satellite)
                .collect(),
            initialized_satellites: members
                .iter()
                .filter(|i| i.retained_cycles.is_none())
                .map(|i| i.satellite)
                .collect(),
            offset_cycles: offset,
        });
    }
    let values: Vec<_> = inputs
        .iter()
        .map(|input| {
            input.retained_cycles.map_or(input.raw_bias_cycles, |old| {
                old + systems
                    .iter()
                    .find(|s| s.system == input.satellite.system())
                    .expect("validated system")
                    .offset_cycles
            })
        })
        .collect();
    if values.iter().any(|v| !v.is_finite()) {
        return Err(Error::NonFinite);
    }
    Ok(DatumAlignmentResult {
        single_difference_cycles: values,
        systems,
    })
}
