//! GPS L1 C/A + Galileo E1 SPP with independent receiver clocks.
//! Broadcast group delays use the pinned native legacy-primary convention.
//! Full upstream SPP FDE/QC and performance equivalence are not implied.

use crate::clock_fit::{ClockLinearFit, least_squares_with_clocks};
use crate::navigation::NavigationData;
use crate::observation::{Observation, ObservationEpoch};
use crate::spp::{self, RejectionReason, SppConfig, SppSolution};
use crate::{
    EcefCoord, Error, GnssSystem, GnssTime, PositionSolution, SatelliteId, SolutionStatus,
    constants as c,
};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObservationDisposition {
    Used,
    UnsupportedTracking,
    MissingCode,
    InvalidCode,
    AlternateTracking,
    MissingEphemeris,
    BelowElevationMask,
    OutlierRejected,
    RaimFdeRejected,
}
#[derive(Debug, Clone, PartialEq)]
pub struct CodeObservationEvidence {
    /// The original receiver row, including phase, Doppler, LLI and flags.
    pub raw_observation: Observation,
    pub disposition: ObservationDisposition,
}
#[derive(Debug, Clone, PartialEq)]
pub struct BroadcastCodeEvidence {
    pub satellite: SatelliteId,
    pub tracking_code: String,
    pub toe: GnssTime,
    pub toc: GnssTime,
    pub data_source: u16,
    /// Legacy-primary TGD/BGD is added to the predicted code, matching the
    /// pinned default SPP path. No I/NAV/F/NAV source-specific switch is made.
    pub primary_group_delay_s: f64,
    pub secondary_group_delay_s: Option<f64>,
    pub variance_m2: f64,
}
#[derive(Debug, Clone, PartialEq)]
pub struct MultiClockSppSolution {
    /// Scalar receiver_clock_bias labels the first clock in clock_biases_s.
    pub spp: SppSolution,
    pub clock_biases_s: Vec<(GnssSystem, f64)>,
    /// [ECEF metres, then clock range biases in metres], including all
    /// position/clock and cross-clock covariances; configured noise scaling.
    pub full_covariance_m2: Vec<Vec<f64>>,
    pub observations: Vec<CodeObservationEvidence>,
    pub broadcast_models: Vec<BroadcastCodeEvidence>,
    pub input_epoch_flag: u8,
    pub input_receiver_clock_offset_s: Option<f64>,
    pub variance_model: crate::spp_quality::SppVarianceModel,
}

pub(crate) fn fit(rows: &[spp::Row], systems: &[GnssSystem]) -> Result<ClockLinearFit, Error> {
    if systems.len() == 1 {
        let (mean, covariance) = spp::least_squares(rows)?;
        Ok(ClockLinearFit {
            mean: mean.to_vec(),
            covariance: covariance.iter().map(|r| r.to_vec()).collect(),
            design: rows.iter().map(|r| r.design.to_vec()).collect(),
        })
    } else {
        least_squares_with_clocks(rows, systems)
    }
}

/// Code-only GPS 1C / Galileo 1C or 1X; valid code is required, phase is not.
/// Preserve every input row and deterministic 1C-before-1X choice as evidence.
/// Every represented system has its own clock; singleton clocks provide no
/// position information and no prior/floor is used to repair deficient rank.
pub fn solve_gps_galileo_epoch(
    epoch: &ObservationEpoch,
    navigation: &NavigationData,
    config: SppConfig,
) -> Result<MultiClockSppSolution, Error> {
    solve_gps_galileo_epoch_with_variance_model(
        epoch,
        navigation,
        config,
        crate::spp_quality::SppVarianceModel::ElevationOnly,
    )
}

/// Explicit alternative measurement variance; orbit, clock and correction
/// models are unchanged. The default API above remains elevation-only.
pub fn solve_gps_galileo_epoch_with_variance_model(
    epoch: &ObservationEpoch,
    navigation: &NavigationData,
    config: SppConfig,
    variance_model: crate::spp_quality::SppVarianceModel,
) -> Result<MultiClockSppSolution, Error> {
    solve_with_initial_clocks(epoch, navigation, config, variance_model, &[])
}

pub(crate) fn solve_with_initial_clocks(
    epoch: &ObservationEpoch,
    navigation: &NavigationData,
    config: SppConfig,
    variance_model: crate::spp_quality::SppVarianceModel,
    initial_clock_biases_s: &[(GnssSystem, f64)],
) -> Result<MultiClockSppSolution, Error> {
    config.validate()?;
    if let Some(offset) = epoch.receiver_clock_offset_s {
        if !offset.is_finite() {
            return Err(Error::NonFinite);
        }
        if offset != 0.0 {
            return Err(Error::UnsupportedReceiverClockOffset);
        }
    }
    if !matches!(epoch.flag, 0 | 1) {
        return Err(Error::InsufficientObservations);
    }
    let mut selected = BTreeMap::<SatelliteId, usize>::new();
    let mut observations = Vec::new();
    for (i, obs) in epoch.observations.iter().enumerate() {
        let supported = match obs.satellite.system() {
            GnssSystem::Gps => obs.tracking_code == "1C",
            GnssSystem::Galileo => matches!(obs.tracking_code.as_str(), "1C" | "1X"),
            _ => false,
        };
        let disposition = if !supported {
            ObservationDisposition::UnsupportedTracking
        } else {
            match obs.pseudorange_m {
                None => ObservationDisposition::MissingCode,
                Some(v) if !v.is_finite() || v <= 0.0 => ObservationDisposition::InvalidCode,
                Some(_) => ObservationDisposition::Used,
            }
        };
        observations.push(CodeObservationEvidence {
            raw_observation: obs.clone(),
            disposition,
        });
        if disposition != ObservationDisposition::Used {
            continue;
        }
        if let Some(&old) = selected.get(&obs.satellite) {
            let rank = |j: usize| epoch.observations[j].tracking_code != "1C";
            if !rank(i) && rank(old) {
                observations[old].disposition = ObservationDisposition::AlternateTracking;
                selected.insert(obs.satellite, i);
            } else {
                observations[i].disposition = ObservationDisposition::AlternateTracking;
            }
        } else {
            selected.insert(obs.satellite, i);
        }
    }
    let selected: BTreeMap<_, _> = selected
        .iter()
        .map(|(&s, &i)| (s, &epoch.observations[i]))
        .collect();
    let mut position = epoch.approximate_position.unwrap_or_default();
    if !position.is_finite() {
        return Err(Error::NonFinite);
    }
    let mut clocks: BTreeMap<_, _> = initial_clock_biases_s
        .iter()
        .map(|&(s, v)| (s, v * c::SPEED_OF_LIGHT))
        .collect();
    for iteration in 1..=config.max_iterations {
        let (rows, _) = spp::build_rows_with_variance(
            epoch,
            navigation,
            &selected,
            position,
            |s| clocks.get(&s).copied().unwrap_or(0.0),
            config,
            variance_model,
        )?;
        let systems: Vec<_> = rows
            .iter()
            .map(|r| r.satellite.system())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let linear = fit(&rows, &systems)?;
        position = EcefCoord::new(
            position.x + linear.mean[0],
            position.y + linear.mean[1],
            position.z + linear.mean[2],
        );
        for (i, &system) in systems.iter().enumerate() {
            *clocks.entry(system).or_insert(0.0) += linear.mean[3 + i];
        }
        if !position.is_finite() || !clocks.values().all(|v| v.is_finite()) {
            return Err(Error::NonFinite);
        }
        let converged = linear.mean[0].hypot(linear.mean[1]).hypot(linear.mean[2])
            <= config.convergence_m
            && linear.mean[3..]
                .iter()
                .all(|v| v.abs() <= config.convergence_m);
        if !converged {
            continue;
        }
        let (rows, rejected) = spp::build_rows_with_variance(
            epoch,
            navigation,
            &selected,
            position,
            |s| clocks.get(&s).copied().unwrap_or(0.0),
            config,
            variance_model,
        )?;
        let final_systems: Vec<_> = rows
            .iter()
            .map(|r| r.satellite.system())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        let final_fit = fit(&rows, &final_systems)?;
        if final_systems != systems
            || final_fit.mean[0]
                .hypot(final_fit.mean[1])
                .hypot(final_fit.mean[2])
                > config.convergence_m
            || final_fit.mean[3..]
                .iter()
                .any(|v| v.abs() > config.convergence_m)
        {
            continue;
        }
        for evidence in &mut observations {
            if evidence.disposition != ObservationDisposition::Used {
                continue;
            }
            if let Some((_, reason)) = rejected
                .iter()
                .find(|(s, _)| *s == evidence.raw_observation.satellite)
            {
                evidence.disposition = match reason {
                    RejectionReason::MissingEphemeris => ObservationDisposition::MissingEphemeris,
                    RejectionReason::BelowElevationMask => {
                        ObservationDisposition::BelowElevationMask
                    }
                    RejectionReason::MissingCode => ObservationDisposition::MissingCode,
                    RejectionReason::UnsupportedSignal => {
                        ObservationDisposition::UnsupportedTracking
                    }
                };
            }
        }
        let models = rows
            .iter()
            .map(|r| {
                let raw = selected[&r.satellite];
                let approximate_tx = epoch
                    .time
                    .checked_add_seconds(-raw.pseudorange_m.unwrap() / c::SPEED_OF_LIGHT)?;
                let eph = navigation.ephemeris(
                    r.satellite,
                    approximate_tx,
                    config.max_ephemeris_age_s,
                )?;
                Ok(BroadcastCodeEvidence {
                    satellite: r.satellite,
                    tracking_code: raw.tracking_code.clone(),
                    toe: eph.toe,
                    toc: eph.toc,
                    data_source: eph.data_source,
                    primary_group_delay_s: eph.clock.tgd,
                    secondary_group_delay_s: eph.clock.secondary_tgd,
                    variance_m2: 1.0 / (r.weight * r.weight),
                })
            })
            .collect::<Result<Vec<_>, Error>>()?;
        let clock_biases_s: Vec<_> = systems
            .iter()
            .map(|s| (*s, clocks[s] / c::SPEED_OF_LIGHT))
            .collect();
        return Ok(MultiClockSppSolution {
            spp: SppSolution {
                position: PositionSolution {
                    time: epoch.time,
                    status: SolutionStatus::Spp,
                    position_ecef: Some(position),
                    position_covariance: Some(std::array::from_fn(|i| {
                        std::array::from_fn(|j| final_fit.covariance[i][j])
                    })),
                    num_satellites: rows.len(),
                    receiver_clock_bias: clock_biases_s[0].1,
                    ..Default::default()
                },
                iterations: iteration,
                residual_rms_m: (rows.iter().map(|r| r.residual * r.residual).sum::<f64>()
                    / rows.len() as f64)
                    .sqrt(),
                residuals_m: rows.iter().map(|r| (r.satellite, r.residual)).collect(),
                rejected,
            },
            clock_biases_s,
            full_covariance_m2: final_fit.covariance,
            observations,
            broadcast_models: models,
            input_epoch_flag: epoch.flag,
            input_receiver_clock_offset_s: epoch.receiver_clock_offset_s,
            variance_model,
        });
    }
    Err(Error::NonConvergent)
}
