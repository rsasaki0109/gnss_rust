//! One-pass MAD refit and leave-one-out RAIM FDE for mixed-clock code SPP.
//! Control flow and thresholds follow pinned libgnss++ spp.cpp.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
//!
//! Geometry deliberately uses each actual system clock, unlike the pinned
//! native final common-clock DOP. This is an explicit, separately reported
//! model difference, not whole-native SPP equivalence.
use crate::navigation::NavigationData;
use crate::observation::ObservationEpoch;
use crate::spp::{self, SppConfig};
use crate::spp_multiclock::{
    MultiClockSppSolution, ObservationDisposition, fit, solve_with_initial_clocks,
};
use crate::spp_quality::{SppVarianceModel, residual_inlier_mask, select_raim_fde_candidate};
use crate::{Error, GnssSystem, SatelliteId};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SppQualityConfig {
    pub outlier_threshold_sigma: Option<f64>,
    pub enable_raim_fde: bool,
    pub raim_minimum_improvement_ratio: f64,
    pub raim_minimum_improvement_m: f64,
    /// Unit-noise geometry with one clock per represented system. Native's
    /// numerical default is 50 but its mixed-system DOP uses a common clock.
    pub max_gdop: Option<f64>,
}
impl Default for SppQualityConfig {
    fn default() -> Self {
        Self {
            outlier_threshold_sigma: Some(3.0),
            enable_raim_fde: true,
            raim_minimum_improvement_ratio: 0.25,
            raim_minimum_improvement_m: 1.0,
            max_gdop: Some(50.0),
        }
    }
}
impl SppQualityConfig {
    pub(crate) fn validate(self) -> Result<(), Error> {
        if [self.outlier_threshold_sigma, self.max_gdop]
            .into_iter()
            .flatten()
            .any(|v| !v.is_finite() || v <= 0.0)
            || [
                self.raim_minimum_improvement_ratio,
                self.raim_minimum_improvement_m,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(Error::InvalidConfiguration);
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SppRefitEvidence {
    pub excluded_satellites: Vec<SatelliteId>,
    pub satellites: Option<usize>,
    pub residual_rms_m: Option<f64>,
    pub independent_clock_gdop: Option<f64>,
    pub error: Option<Error>,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SppQualitySummary {
    pub pre_qc_satellites: usize,
    pub pre_qc_residual_rms_m: f64,
    pub final_residual_rms_m: f64,
    pub outlier_rejections: usize,
    pub raim_attempts: usize,
    pub raim_rejections: usize,
    pub independent_clock_gdop: f64,
}
#[derive(Debug, Clone, PartialEq)]
pub struct SppQualityEvidence {
    pub config: SppQualityConfig,
    pub summary: SppQualitySummary,
    pub pre_qc_max_abs_residual_m: f64,
    pub pre_qc_clock_systems: Vec<GnssSystem>,
    pub pre_qc_degrees_of_freedom: usize,
    /// A proposed mask can fail its refit. Rejected proposals do not remove
    /// observations and do not alter the baseline used by RAIM candidates.
    pub mad_refit: Option<SppRefitEvidence>,
    pub mad_refit_accepted: bool,
    pub raim_candidates: Vec<SppRefitEvidence>,
    pub raim_excluded_satellite: Option<SatelliteId>,
}
#[derive(Debug, Clone, PartialEq)]
pub struct QualityControlledSppSolution {
    pub solution: MultiClockSppSolution,
    pub evidence: SppQualityEvidence,
}

fn gdop(
    epoch: &ObservationEpoch,
    navigation: &NavigationData,
    config: SppConfig,
    solution: &MultiClockSppSolution,
) -> Result<f64, Error> {
    let selected: BTreeMap<_, _> = solution
        .observations
        .iter()
        .filter(|o| o.disposition == ObservationDisposition::Used)
        .map(|o| (o.raw_observation.satellite, &o.raw_observation))
        .collect();
    let clocks: BTreeMap<_, _> = solution.clock_biases_s.iter().copied().collect();
    let (rows, _) = spp::build_rows_with_variance(
        epoch,
        navigation,
        &selected,
        solution
            .spp
            .position
            .position_ecef
            .ok_or(Error::NonFinite)?,
        |s| clocks[&s] * crate::constants::SPEED_OF_LIGHT,
        config,
        solution.variance_model,
    )?;
    let unit_rows: Vec<_> = rows
        .iter()
        .map(|r| spp::Row {
            satellite: r.satellite,
            design: r.design,
            residual: 0.0,
            weight: 1.0,
        })
        .collect();
    let systems: Vec<_> = rows
        .iter()
        .map(|r| r.satellite.system())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let covariance = fit(&unit_rows, &systems)?.covariance;
    let trace = (0..covariance.len()).map(|i| covariance[i][i]).sum::<f64>();
    if trace.is_finite() && trace > 0.0 {
        Ok(trace.sqrt())
    } else {
        Err(Error::NonFinite)
    }
}

fn gate(gdop: f64, config: SppQualityConfig) -> Result<(), Error> {
    if config.max_gdop.is_some_and(|limit| gdop > limit) {
        Err(Error::SppGeometryTooWeak)
    } else {
        Ok(())
    }
}

fn refit(
    epoch: &ObservationEpoch,
    nav: &NavigationData,
    config: SppConfig,
    quality: SppQualityConfig,
    baseline: &MultiClockSppSolution,
    excluded: &[SatelliteId],
) -> (Option<MultiClockSppSolution>, SppRefitEvidence) {
    let mut evidence = SppRefitEvidence {
        excluded_satellites: excluded.to_vec(),
        satellites: None,
        residual_rms_m: None,
        independent_clock_gdop: None,
        error: None,
    };
    let used: BTreeSet<_> = baseline
        .observations
        .iter()
        .filter(|o| {
            o.disposition == ObservationDisposition::Used
                && !excluded.contains(&o.raw_observation.satellite)
        })
        .map(|o| {
            (
                o.raw_observation.satellite,
                o.raw_observation.tracking_code.as_str(),
            )
        })
        .collect();
    let mut child = epoch.clone();
    child.approximate_position = baseline.spp.position.position_ecef;
    child
        .observations
        .retain(|o| used.contains(&(o.satellite, o.tracking_code.as_str())));
    let result = solve_with_initial_clocks(
        &child,
        nav,
        config,
        baseline.variance_model,
        &baseline.clock_biases_s,
    )
    .and_then(|s| {
        evidence.satellites = Some(s.spp.position.num_satellites);
        evidence.residual_rms_m = Some(s.spp.residual_rms_m);
        let geometry = gdop(&child, nav, config, &s)?;
        evidence.independent_clock_gdop = Some(geometry);
        gate(geometry, quality)?;
        Ok(s)
    });
    match result {
        Ok(s) => (Some(s), evidence),
        Err(e) => {
            evidence.error = Some(e);
            (None, evidence)
        }
    }
}

fn restore_raw_evidence(
    child: &mut MultiClockSppSolution,
    baseline: &MultiClockSppSolution,
    excluded: &[SatelliteId],
    reason: ObservationDisposition,
) {
    let mut observations = baseline.observations.clone();
    for o in &mut observations {
        if o.disposition != ObservationDisposition::Used {
            continue;
        }
        o.disposition = if excluded.contains(&o.raw_observation.satellite) {
            reason
        } else {
            child
                .observations
                .iter()
                .find(|c| {
                    c.raw_observation.satellite == o.raw_observation.satellite
                        && c.raw_observation.tracking_code == o.raw_observation.tracking_code
                })
                .expect("child retains every eligible raw row")
                .disposition
        };
    }
    child.observations = observations;
    for r in &baseline.spp.rejected {
        if !child.spp.rejected.contains(r) {
            child.spp.rejected.push(*r);
        }
    }
}

/// Run a baseline, at most one MAD refit, then (only if MAD was not adopted)
/// all leave-one-out RAIM refits. Children never recurse into QC. Candidate
/// failures leave the baseline intact and clocks stay paired with each fit.
pub fn solve_gps_galileo_epoch_with_quality_control(
    epoch: &ObservationEpoch,
    navigation: &NavigationData,
    config: SppConfig,
    variance_model: SppVarianceModel,
    quality: SppQualityConfig,
) -> Result<QualityControlledSppSolution, Error> {
    quality.validate()?;
    let baseline = solve_with_initial_clocks(epoch, navigation, config, variance_model, &[])?;
    let m = baseline.spp.residuals_m.len();
    let n = 3 + baseline.clock_biases_s.len();
    let mut evidence = SppQualityEvidence {
        config: quality,
        summary: SppQualitySummary {
            pre_qc_satellites: m,
            pre_qc_residual_rms_m: baseline.spp.residual_rms_m,
            final_residual_rms_m: baseline.spp.residual_rms_m,
            outlier_rejections: 0,
            raim_attempts: 0,
            raim_rejections: 0,
            independent_clock_gdop: 0.0,
        },
        pre_qc_max_abs_residual_m: baseline
            .spp
            .residuals_m
            .iter()
            .map(|(_, v)| v.abs())
            .fold(0.0, f64::max),
        pre_qc_clock_systems: baseline.clock_biases_s.iter().map(|(s, _)| *s).collect(),
        pre_qc_degrees_of_freedom: m.saturating_sub(n),
        mad_refit: None,
        mad_refit_accepted: false,
        raim_candidates: Vec::new(),
        raim_excluded_satellite: None,
    };
    if let Some(threshold) = quality.outlier_threshold_sigma.filter(|_| m > n) {
        let residuals: Vec<_> = baseline.spp.residuals_m.iter().map(|(_, v)| *v).collect();
        let mask = residual_inlier_mask(&residuals, threshold, config.pseudorange_sigma_m.max(1.0));
        let excluded: Vec<_> = baseline
            .spp
            .residuals_m
            .iter()
            .zip(mask)
            .filter(|(_, keep)| !*keep)
            .map(|((s, _), _)| *s)
            .collect();
        if !excluded.is_empty() {
            let (solution, attempt) =
                refit(epoch, navigation, config, quality, &baseline, &excluded);
            evidence.mad_refit = Some(attempt);
            if let Some(mut solution) = solution {
                evidence.mad_refit_accepted = true;
                evidence.summary.outlier_rejections = excluded.len();
                evidence.summary.final_residual_rms_m = solution.spp.residual_rms_m;
                evidence.summary.independent_clock_gdop = evidence
                    .mad_refit
                    .as_ref()
                    .unwrap()
                    .independent_clock_gdop
                    .unwrap();
                restore_raw_evidence(
                    &mut solution,
                    &baseline,
                    &excluded,
                    ObservationDisposition::OutlierRejected,
                );
                return Ok(QualityControlledSppSolution { solution, evidence });
            }
        }
    }
    if quality.enable_raim_fde && m.saturating_sub(n) >= 2 && baseline.spp.residual_rms_m > 0.0 {
        let mut solutions = Vec::new();
        let mut rms = Vec::new();
        for &(s, _) in &baseline.spp.residuals_m {
            let (solution, attempt) = refit(epoch, navigation, config, quality, &baseline, &[s]);
            rms.push(
                solution
                    .as_ref()
                    .map_or(f64::INFINITY, |s| s.spp.residual_rms_m),
            );
            solutions.push(solution);
            evidence.raim_candidates.push(attempt);
        }
        evidence.summary.raim_attempts = solutions.len();
        if let Some(index) = select_raim_fde_candidate(
            baseline.spp.residual_rms_m,
            &rms,
            quality.raim_minimum_improvement_ratio,
            quality.raim_minimum_improvement_m,
        ) {
            let mut solution = solutions[index]
                .take()
                .expect("selected candidate is valid");
            let excluded = baseline.spp.residuals_m[index].0;
            evidence.raim_excluded_satellite = Some(excluded);
            evidence.summary.raim_rejections = 1;
            evidence.summary.final_residual_rms_m = solution.spp.residual_rms_m;
            evidence.summary.independent_clock_gdop = evidence.raim_candidates[index]
                .independent_clock_gdop
                .unwrap();
            restore_raw_evidence(
                &mut solution,
                &baseline,
                &[excluded],
                ObservationDisposition::RaimFdeRejected,
            );
            return Ok(QualityControlledSppSolution { solution, evidence });
        }
    }
    let geometry = gdop(epoch, navigation, config, &baseline)?;
    gate(geometry, quality)?;
    evidence.summary.independent_clock_gdop = geometry;
    Ok(QualityControlledSppSolution {
        solution: baseline,
        evidence,
    })
}
