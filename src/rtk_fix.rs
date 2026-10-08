//! Initial full-set GPS L1 integer validation on immutable FLOAT snapshots.
//!
//! Head-state conditioning and position-independent code/phase validation are
//! derived from libgnss++ rtk_ar_evaluation / rtk_cp_pr_gate components.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
//! This policy does not hold integers or feed fixed values into the FLOAT filter.

use crate::filter::{self, FilterError, GaussianState, Matrix};
use crate::lambda::LambdaConfig;
use crate::rtk::{RtkError, RtkFloatSolution};
use crate::{EcefCoord, GnssTime, PositionSolution, SatelliteId, SolutionStatus};
use std::{collections::BTreeMap, fmt};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FixError {
    InvalidConfiguration,
    InvalidInteger,
    OutOfOrder,
    Filter(FilterError),
    Rtk(RtkError),
}
impl fmt::Display for FixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfiguration => f.write_str("invalid integer validation configuration"),
            Self::InvalidInteger => {
                f.write_str("integer candidate exceeds the exact arithmetic range")
            }
            Self::OutOfOrder => {
                f.write_str("FIX policy evaluations must increase strictly in GPS time")
            }
            Self::Filter(e) => e.fmt(f),
            Self::Rtk(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for FixError {}
impl From<FilterError> for FixError {
    fn from(e: FilterError) -> Self {
        Self::Filter(e)
    }
}
impl From<RtkError> for FixError {
    fn from(e: RtkError) -> Self {
        Self::Rtk(e)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ConditionedHead {
    /// Same units and coordinates as the input head (RTK: ECEF baseline metres).
    pub state: GaussianState,
    pub squared_ambiguity_residual: f64,
}

/// x_fixed = x_float - Qab Qbb^-1 (b_float - b_integer).
/// Q_fixed = Qaa - Qab Qbb^-1 Qba. Computes the covariance from the joint
/// Cholesky factor's conditional block, avoiding subtractive cancellation.
/// The joint covariance must be positive definite; no diagonal floor is added.
pub fn condition_head(
    head: &GaussianState,
    dd_float: &[f64],
    qbb: &Matrix,
    qab: &Matrix,
    integers: &[i64],
) -> Result<ConditionedHead, FixError> {
    let h = head.mean.len();
    let k = dd_float.len();
    if !(1..=128).contains(&h)
        || !(1..=128).contains(&k)
        || integers.len() != k
        || head.covariance.len() != h
        || head.covariance.iter().any(|r| r.len() != h)
        || qbb.len() != k
        || qbb.iter().any(|r| r.len() != k)
        || qab.len() != h
        || qab.iter().any(|r| r.len() != k)
    {
        return Err(FilterError::InvalidDimensions.into());
    }
    if !head
        .mean
        .iter()
        .chain(dd_float)
        .chain(qab.iter().flatten())
        .all(|v| v.is_finite())
    {
        return Err(FilterError::NonFinite.into());
    }
    if integers.iter().any(|i| i.unsigned_abs() >= (1_u64 << 52)) {
        return Err(FixError::InvalidInteger);
    }
    filter::cholesky(&head.covariance)?;
    filter::cholesky(qbb)?;
    let mut joint = filter::zeros(k + h, k + h);
    for i in 0..k {
        joint[i][..k].copy_from_slice(&qbb[i]);
    }
    for (i, row) in qab.iter().enumerate() {
        for (j, &value) in row.iter().enumerate() {
            joint[k + i][j] = value;
            joint[j][k + i] = value;
        }
        joint[k + i][k..].copy_from_slice(&head.covariance[i]);
    }
    let l = filter::cholesky(&joint)?;
    let lb: Matrix = l[..k].iter().map(|r| r[..k].to_vec()).collect();
    let delta: Vec<_> = dd_float
        .iter()
        .zip(integers)
        .map(|(f, i)| f - *i as f64)
        .collect();
    let weighted = filter::solve_cholesky(&lb, &delta);
    let mean: Vec<_> = head
        .mean
        .iter()
        .zip(qab)
        .map(|(x, row)| x - row.iter().zip(&weighted).map(|(a, b)| a * b).sum::<f64>())
        .collect();
    let lc: Matrix = l[k..].iter().map(|r| r[k..].to_vec()).collect();
    let covariance = filter::multiply(&lc, &filter::transpose(&lc));
    let cost: f64 = delta.iter().zip(&weighted).map(|(a, b)| a * b).sum();
    if !mean
        .iter()
        .chain(covariance.iter().flatten())
        .all(|v| v.is_finite())
        || !cost.is_finite()
        || cost < 0.0
    {
        return Err(FilterError::NumericalFailure.into());
    }
    filter::cholesky(&covariance)?;
    Ok(ConditionedHead {
        state: GaussianState { mean, covariance },
        squared_ambiguity_residual: cost,
    })
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CodePhaseObservation {
    pub dd_code_m: f64,
    pub dd_carrier_m: f64,
    pub fixed_ambiguity_m: f64,
}
#[derive(Debug, Clone, PartialEq)]
pub struct CodePhaseCheck {
    pub consistent: bool,
    pub innovations_m: Vec<f64>,
    pub rms_innovation_m: f64,
    pub max_abs_innovation_m: f64,
    pub bad_pairs: usize,
}

/// Position-independent DD_PR - (DD_CP - fixed ambiguity in metres).
/// Checks every pair; nonfinite inputs fail rather than being dropped.
pub fn code_phase_check(
    observations: &[CodePhaseObservation],
    threshold_m: f64,
) -> Result<CodePhaseCheck, FixError> {
    if observations.is_empty() || observations.len() > 128 {
        return Err(FilterError::InvalidDimensions.into());
    }
    if !threshold_m.is_finite() || threshold_m <= 0.0 {
        return Err(FixError::InvalidConfiguration);
    }
    let mut innovations = Vec::with_capacity(observations.len());
    for o in observations {
        if ![o.dd_code_m, o.dd_carrier_m, o.fixed_ambiguity_m]
            .iter()
            .all(|v| v.is_finite())
        {
            return Err(FilterError::NonFinite.into());
        }
        let innovation = o.dd_code_m - (o.dd_carrier_m - o.fixed_ambiguity_m);
        if !innovation.is_finite() {
            return Err(FilterError::NumericalFailure.into());
        }
        innovations.push(innovation);
    }
    summarize_code_phase(innovations, threshold_m)
}
fn summarize_code_phase(
    innovations: Vec<f64>,
    threshold_m: f64,
) -> Result<CodePhaseCheck, FixError> {
    let rms = rms(&innovations);
    let max = maximum(&innovations);
    let bad_pairs = innovations.iter().filter(|v| v.abs() > threshold_m).count();
    if !rms.is_finite() || !max.is_finite() {
        return Err(FilterError::NumericalFailure.into());
    }
    Ok(CodePhaseCheck {
        consistent: bad_pairs == 0,
        innovations_m: innovations,
        rms_innovation_m: rms,
        max_abs_innovation_m: max,
        bad_pairs,
    })
}
fn rms(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |sum, v| sum.hypot(*v)) / (values.len() as f64).sqrt()
}
fn maximum(values: &[f64]) -> f64 {
    values.iter().map(|v| v.abs()).fold(0.0, f64::max)
}
fn position(base: EcefCoord, baseline: &[f64]) -> EcefCoord {
    EcefCoord::new(
        base.x + baseline[0],
        base.y + baseline[1],
        base.z + baseline[2],
    )
}
fn distance(a: EcefCoord, b: EcefCoord) -> f64 {
    (a.x - b.x).hypot(a.y - b.y).hypot(a.z - b.z)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RtkFixConfig {
    pub ratio_threshold: f64,
    pub min_lock_epochs: usize,
    pub min_ambiguity_pairs: usize,
    pub min_confirmation_epochs: usize,
    pub max_float_correction_m: f64,
    pub max_position_std_m: f64,
    pub max_phase_residual_rms_m: f64,
    pub max_phase_residual_abs_m: f64,
    pub max_code_residual_rms_m: f64,
    pub max_code_phase_innovation_m: f64,
    pub max_postfit_nis_per_row: f64,
    pub max_history_jump_m: f64,
    pub max_history_age_s: f64,
}
impl Default for RtkFixConfig {
    fn default() -> Self {
        Self {
            ratio_threshold: 3.0,
            min_lock_epochs: 5,
            min_ambiguity_pairs: 4,
            min_confirmation_epochs: 2,
            max_float_correction_m: 1.0,
            max_position_std_m: 0.05,
            max_phase_residual_rms_m: 0.01,
            max_phase_residual_abs_m: 0.03,
            max_code_residual_rms_m: 1.0,
            max_code_phase_innovation_m: 0.5,
            max_postfit_nis_per_row: 9.0,
            max_history_jump_m: 0.1,
            max_history_age_s: 120.0,
        }
    }
}
impl RtkFixConfig {
    fn validate(self) -> Result<(), FixError> {
        let positive = [
            self.max_float_correction_m,
            self.max_position_std_m,
            self.max_phase_residual_rms_m,
            self.max_phase_residual_abs_m,
            self.max_code_residual_rms_m,
            self.max_code_phase_innovation_m,
            self.max_postfit_nis_per_row,
            self.max_history_jump_m,
            self.max_history_age_s,
        ];
        if positive.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || !self.ratio_threshold.is_finite()
            || self.ratio_threshold <= 1.0
            || !(1..=10_000).contains(&self.min_lock_epochs)
            || !(4..=128).contains(&self.min_ambiguity_pairs)
            || !(1..=10_000).contains(&self.min_confirmation_epochs)
        {
            return Err(FixError::InvalidConfiguration);
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FixDecision {
    Accepted,
    InsufficientPairs,
    InsufficientLock,
    SearchFailed,
    RatioRejected,
    ConditioningFailed,
    PositionCorrection,
    CovarianceTooLarge,
    PhaseResiduals,
    CodeResiduals,
    CodePhaseInconsistency,
    PostfitInnovation,
    BaselineTooLong,
    HistoryJump,
    Confirming,
}
impl fmt::Display for FixDecision {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Accepted => "accepted",
            Self::InsufficientPairs => "insufficient_pairs",
            Self::InsufficientLock => "insufficient_lock",
            Self::SearchFailed => "search_failed",
            Self::RatioRejected => "ratio_rejected",
            Self::ConditioningFailed => "conditioning_failed",
            Self::PositionCorrection => "position_correction",
            Self::CovarianceTooLarge => "covariance_too_large",
            Self::PhaseResiduals => "phase_residuals",
            Self::CodeResiduals => "code_residuals",
            Self::CodePhaseInconsistency => "code_phase_inconsistency",
            Self::PostfitInnovation => "postfit_innovation",
            Self::BaselineTooLong => "baseline_too_long",
            Self::HistoryJump => "history_jump",
            Self::Confirming => "confirming",
        })
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct CandidateValidation {
    pub conditioned: ConditionedHead,
    pub phase_residual_rms_m: f64,
    pub phase_residual_max_abs_m: f64,
    pub code_residual_rms_m: f64,
    pub code_phase: CodePhaseCheck,
    pub postfit_nis_per_row: f64,
    pub float_correction_m: f64,
    pub max_position_std_m: f64,
    /// None means local numerical/residual tests passed. Eligibility, ratio,
    /// repeat confirmation and history are separate policy checks.
    pub rejection: Option<FixDecision>,
}

/// Validate supplied integers without declaring FIX. Useful for independent
/// candidate diagnostics and wrong-integer tests. Only the immutable FLOAT
/// observation snapshot is used; no state is modified.
pub fn validate_candidate(
    float: &RtkFloatSolution,
    integers: &[i64],
    config: RtkFixConfig,
) -> Result<CandidateValidation, FixError> {
    config.validate()?;
    let context = &float.fix_context;
    let ambiguity = context.ambiguities()?;
    let conditioned = condition_head(
        &context.head_prior(),
        &ambiguity.dd_float_cycles,
        &ambiguity.ambiguity_covariance_cycles2,
        &ambiguity.head_ambiguity_covariance,
        integers,
    )?;
    let residuals = context.residuals(&conditioned.state.mean, integers)?;
    let code_phase = summarize_code_phase(
        residuals.code_phase_innovations_m,
        config.max_code_phase_innovation_m,
    )?;
    let mut v = residuals.phase_m.clone();
    v.extend_from_slice(&residuals.code_m);
    let l = filter::cholesky(&residuals.noise)?;
    let weighted = filter::solve_cholesky(&l, &v);
    let nis: f64 = v.iter().zip(&weighted).map(|(a, b)| a * b).sum::<f64>() / v.len() as f64;
    let phase_rms = rms(&residuals.phase_m);
    let phase_max = maximum(&residuals.phase_m);
    let code_rms = rms(&residuals.code_m);
    let fixed_position = position(context.base_position, &conditioned.state.mean);
    let correction = distance(
        fixed_position,
        position(context.base_position, &context.state.mean),
    );
    let std = (0..3)
        .map(|i| conditioned.state.covariance[i][i].sqrt())
        .fold(0.0, f64::max);
    if ![nis, phase_rms, phase_max, code_rms, correction, std]
        .iter()
        .all(|v| v.is_finite())
        || nis < 0.0
    {
        return Err(FilterError::NumericalFailure.into());
    }
    let rejection =
        if distance(context.base_position, fixed_position) > context.config.max_baseline_m {
            Some(FixDecision::BaselineTooLong)
        } else if correction > config.max_float_correction_m {
            Some(FixDecision::PositionCorrection)
        } else if std > config.max_position_std_m {
            Some(FixDecision::CovarianceTooLarge)
        } else if phase_rms > config.max_phase_residual_rms_m
            || phase_max > config.max_phase_residual_abs_m
        {
            Some(FixDecision::PhaseResiduals)
        } else if code_rms > config.max_code_residual_rms_m {
            Some(FixDecision::CodeResiduals)
        } else if !code_phase.consistent {
            Some(FixDecision::CodePhaseInconsistency)
        } else if nis > config.max_postfit_nis_per_row {
            Some(FixDecision::PostfitInnovation)
        } else {
            None
        };
    Ok(CandidateValidation {
        conditioned,
        phase_residual_rms_m: phase_rms,
        phase_residual_max_abs_m: phase_max,
        code_residual_rms_m: code_rms,
        code_phase,
        postfit_nis_per_row: nis,
        float_correction_m: correction,
        max_position_std_m: std,
        rejection,
    })
}

#[derive(Debug, Clone, PartialEq)]
pub struct RtkFixAttempt {
    pub decision: FixDecision,
    pub ratio: Option<f64>,
    pub integers: Vec<i64>,
    pub validation: Option<CandidateValidation>,
    pub confirmation_count: usize,
    /// Present only after every policy and local validation check passes.
    pub fixed_solution: Option<PositionSolution>,
}
struct Pending {
    key: Vec<(SatelliteId, i64)>,
    position: EcefCoord,
    time: GnssTime,
    count: usize,
}
/// Conservative initial full-set AR policy. Accepted integers are checked anew
/// every epoch; no fix-and-hold or conditional state feeds back to FLOAT.
pub struct RtkFixPolicy {
    config: RtkFixConfig,
    locks: BTreeMap<SatelliteId, usize>,
    last_time: Option<GnssTime>,
    pending: Option<Pending>,
    anchor: Option<(EcefCoord, GnssTime)>,
}
impl RtkFixPolicy {
    pub fn new(config: RtkFixConfig) -> Result<Self, FixError> {
        config.validate()?;
        Ok(Self {
            config,
            locks: BTreeMap::new(),
            last_time: None,
            pending: None,
            anchor: None,
        })
    }
    pub fn reset(&mut self) {
        self.locks.clear();
        self.last_time = None;
        self.pending = None;
        self.anchor = None;
    }
    fn reject(&mut self, mut attempt: RtkFixAttempt, decision: FixDecision) -> RtkFixAttempt {
        self.pending = None;
        attempt.decision = decision;
        attempt
    }
    pub fn evaluate(&mut self, float: &RtkFloatSolution) -> Result<RtkFixAttempt, FixError> {
        let context = &float.fix_context;
        let time = context.time;
        if self
            .last_time
            .is_some_and(|t| time.difference_seconds(t) <= 0.0)
        {
            return Err(FixError::OutOfOrder);
        }
        let satellites = context.satellites();
        let history_gap = self
            .last_time
            .is_some_and(|t| time.difference_seconds(t) > self.config.max_history_age_s);
        if history_gap {
            self.locks.clear();
        }
        if !context.reset_satellites.is_empty() || history_gap {
            self.pending = None;
            self.anchor = None;
        }
        self.locks = satellites
            .iter()
            .map(|&s| {
                let count = if context.reset_satellites.contains(&s) {
                    1
                } else {
                    self.locks.get(&s).copied().unwrap_or(0).saturating_add(1)
                };
                (s, count)
            })
            .collect();
        self.last_time = Some(time);
        let mut attempt = RtkFixAttempt {
            decision: FixDecision::InsufficientPairs,
            ratio: None,
            integers: Vec::new(),
            validation: None,
            confirmation_count: 0,
            fixed_solution: None,
        };
        if satellites.len() - 1 < self.config.min_ambiguity_pairs {
            return Ok(self.reject(attempt, FixDecision::InsufficientPairs));
        }
        if self
            .locks
            .values()
            .any(|count| *count < self.config.min_lock_epochs)
        {
            return Ok(self.reject(attempt, FixDecision::InsufficientLock));
        }
        let ambiguity = match context.ambiguities() {
            Ok(a) => a,
            Err(_) => return Ok(self.reject(attempt, FixDecision::ConditioningFailed)),
        };
        let candidates = match crate::lambda::search(
            &ambiguity.dd_float_cycles,
            &ambiguity.ambiguity_covariance_cycles2,
            LambdaConfig::default(),
        ) {
            Ok(c) => c,
            Err(_) => return Ok(self.reject(attempt, FixDecision::SearchFailed)),
        };
        attempt.ratio = candidates.ratio();
        attempt.integers = candidates.candidates[0].ambiguities.clone();
        if !attempt
            .ratio
            .is_some_and(|r| r.is_finite() && r >= self.config.ratio_threshold)
        {
            return Ok(self.reject(attempt, FixDecision::RatioRejected));
        }
        let validation = match validate_candidate(float, &attempt.integers, self.config) {
            Ok(v) => v,
            Err(_) => return Ok(self.reject(attempt, FixDecision::ConditioningFailed)),
        };
        let rejection = validation.rejection;
        let fixed_position = position(context.base_position, &validation.conditioned.state.mean);
        attempt.validation = Some(validation);
        if let Some(reason) = rejection {
            return Ok(self.reject(attempt, reason));
        }
        if let Some((anchor, at)) = self.anchor {
            if time.difference_seconds(at) > self.config.max_history_age_s {
                self.anchor = None;
            } else if distance(fixed_position, anchor) > self.config.max_history_jump_m {
                return Ok(self.reject(attempt, FixDecision::HistoryJump));
            }
        }
        // Canonical reference-independent DD labels for repeated confirmation.
        let mut dd = BTreeMap::new();
        let reference = context.reference_satellite();
        dd.insert(reference, 0_i64);
        for (&sat, &integer) in satellites
            .iter()
            .filter(|s| **s != reference)
            .zip(&attempt.integers)
        {
            dd.insert(sat, integer);
        }
        let root = *dd.values().next().ok_or(FilterError::InvalidDimensions)?;
        let key: Vec<_> = dd.into_iter().map(|(s, n)| (s, n - root)).collect();
        let count = match &self.pending {
            Some(p)
                if p.key == key
                    && time.difference_seconds(p.time) <= self.config.max_history_age_s
                    && distance(fixed_position, p.position) <= self.config.max_history_jump_m =>
            {
                p.count.saturating_add(1)
            }
            _ => 1,
        };
        self.pending = Some(Pending {
            key,
            position: fixed_position,
            time,
            count,
        });
        attempt.confirmation_count = count;
        if count < self.config.min_confirmation_epochs {
            attempt.decision = FixDecision::Confirming;
            return Ok(attempt);
        }
        let conditioned = &attempt.validation.as_ref().unwrap().conditioned;
        attempt.fixed_solution = Some(PositionSolution {
            time,
            status: SolutionStatus::Fixed,
            position_ecef: Some(fixed_position),
            position_covariance: Some(std::array::from_fn(|i| {
                std::array::from_fn(|j| conditioned.state.covariance[i][j])
            })),
            num_satellites: satellites.len(),
            ..PositionSolution::default()
        });
        attempt.decision = FixDecision::Accepted;
        self.anchor = Some((fixed_position, time));
        Ok(attempt)
    }
}
