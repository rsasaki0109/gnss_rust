//! GPS L1 / Galileo E1 integer validation on immutable FLOAT snapshots.
//! Default FIX uses GPS static snapshots. Kinematic GPS/Galileo is opt-in and
//! requires timestamped velocity evidence for displacement-history checks.
//!
//! Head-state conditioning and position-independent code/phase validation are
//! derived from libgnss++ rtk_ar_evaluation / rtk_cp_pr_gate components.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
//! This policy does not hold integers or feed fixed values into the FLOAT filter.

use crate::filter::{self, FilterError, GaussianState, Matrix};
use crate::lambda::LambdaConfig;
use crate::rtk::{FixContext, RtkError, RtkFloatSolution};
use crate::{EcefCoord, GnssTime, PositionSolution, SatelliteId, SolutionStatus};
use std::{collections::BTreeMap, fmt};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FixError {
    InvalidConfiguration,
    UnsupportedPositionMode,
    UnsupportedConstellations,
    MotionEvidenceRequired,
    InvalidMotionEvidence,
    ReferenceSetLimit,
    InvalidInteger,
    OutOfOrder,
    Filter(FilterError),
    Rtk(RtkError),
}
impl fmt::Display for FixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfiguration => f.write_str("invalid integer validation configuration"),
            Self::UnsupportedPositionMode => {
                f.write_str("FIX policy position mode does not match the FLOAT snapshot")
            }
            Self::UnsupportedConstellations => {
                f.write_str("initial FIX validation requires GPS-only observations")
            }
            Self::MotionEvidenceRequired => {
                f.write_str("kinematic FIX requires explicit motion evaluation")
            }
            Self::InvalidMotionEvidence => {
                f.write_str("invalid or mismatched kinematic motion evidence")
            }
            Self::ReferenceSetLimit => {
                f.write_str("candidate audit exceeds the reference-set limit")
            }
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
    MotionUnavailable,
    MotionUncertain,
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
            Self::MotionUnavailable => "motion_unavailable",
            Self::MotionUncertain => "motion_uncertain",
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
    validate_context(&float.fix_context, integers, config)
}

fn validate_context(
    context: &FixContext,
    integers: &[i64],
    config: RtkFixConfig,
) -> Result<CandidateValidation, FixError> {
    config.validate()?;
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

/// Local validation of the same physical integer candidate in one DD basis.
/// This is not a fresh LAMBDA search or a policy-level FIX decision.
#[derive(Debug, Clone, PartialEq)]
pub struct ReferenceCandidateValidation {
    pub reference_satellite: SatelliteId,
    pub ambiguity_satellites: Vec<SatelliteId>,
    pub integers: Vec<i64>,
    pub validation: CandidateValidation,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CandidateReferenceAudit {
    pub time: GnssTime,
    pub source_reference: SatelliteId,
    pub references: Vec<ReferenceCandidateValidation>,
}
impl CandidateReferenceAudit {
    pub fn local_rejection_depends_on_reference(&self) -> bool {
        self.references.first().is_some_and(|first| {
            self.references
                .iter()
                .any(|entry| entry.validation.rejection != first.validation.rejection)
        })
    }

    /// Require every reference to pass the existing local gates. If several
    /// gates reject, report the first in the same gate order as validate_candidate.
    /// This intersection selects no favourable reference and is independent
    /// of the order of the audited references. It is not a policy FIX decision.
    pub fn conservative_rejection(&self) -> Option<FixDecision> {
        self.references
            .iter()
            .filter_map(|entry| entry.validation.rejection)
            .min_by_key(|reason| match reason {
                FixDecision::BaselineTooLong => 0,
                FixDecision::PositionCorrection => 1,
                FixDecision::CovarianceTooLarge => 2,
                FixDecision::PhaseResiduals => 3,
                FixDecision::CodeResiduals => 4,
                FixDecision::CodePhaseInconsistency => 5,
                FixDecision::PostfitInnovation => 6,
                _ => 7,
            })
    }
}

/// Re-express one supplied candidate under every reference in the immutable
/// FLOAT observation snapshot and run the unchanged local validation gates.
///
/// Integers use the snapshot's original reference-minus-target ordering.
/// Rebased values use checked integer differences and the same strict 2^52
/// bound as condition_head. The complete audit fails if any basis cannot
/// represent the candidate in that range.
///
/// Raw observations, SD state/covariance, filter and FIX-policy history are
/// unchanged. Public display fields do not override the immutable snapshot.
/// Passing local validation is distinct from eligibility, ratio, confirmation
/// and history checks. This function selects no reference or winning result.
pub fn audit_candidate_references(
    float: &RtkFloatSolution,
    integers: &[i64],
    config: RtkFixConfig,
) -> Result<CandidateReferenceAudit, FixError> {
    if float.fix_context.config.constellations != crate::rtk::RtkConstellations::Gps {
        return Err(FixError::UnsupportedConstellations);
    }
    config.validate()?;
    let source = &float.fix_context;
    let satellites = source.satellites();
    if integers.len() + 1 != satellites.len() {
        return Err(FilterError::InvalidDimensions.into());
    }
    if integers.iter().any(|i| i.unsigned_abs() >= (1_u64 << 52)) {
        return Err(FixError::InvalidInteger);
    }
    let original_reference = source.reference_satellite();
    let mut original_dd = BTreeMap::from([(original_reference, 0_i64)]);
    for (satellite, integer) in satellites
        .iter()
        .filter(|s| **s != original_reference)
        .zip(integers)
    {
        original_dd.insert(*satellite, *integer);
    }
    // Check all bases before numerical evaluation; no partially evaluated
    // result is returned for an out-of-range equivalent representation.
    let minimum = original_dd
        .values()
        .min()
        .ok_or(FilterError::InvalidDimensions)?;
    let maximum = original_dd
        .values()
        .max()
        .ok_or(FilterError::InvalidDimensions)?;
    let span = maximum
        .checked_sub(*minimum)
        .ok_or(FixError::InvalidInteger)?;
    if span.unsigned_abs() >= (1_u64 << 52) {
        return Err(FixError::InvalidInteger);
    }
    let mut references = Vec::with_capacity(satellites.len());
    for &reference in &satellites {
        let targets: Vec<_> = satellites
            .iter()
            .copied()
            .filter(|s| *s != reference)
            .collect();
        let equivalent = targets
            .iter()
            .map(|satellite| {
                original_dd[satellite]
                    .checked_sub(original_dd[&reference])
                    .ok_or(FixError::InvalidInteger)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let context = source.with_reference(reference)?;
        let validation = validate_context(&context, &equivalent, config)?;
        references.push(ReferenceCandidateValidation {
            reference_satellite: reference,
            ambiguity_satellites: targets,
            integers: equivalent,
            validation,
        });
    }
    Ok(CandidateReferenceAudit {
        time: source.time,
        source_reference: original_reference,
        references,
    })
}

/// The same candidate under one combination of within-system references.
#[derive(Debug, Clone, PartialEq)]
pub struct SystemReferenceCandidateValidation {
    pub reference_satellites: Vec<SatelliteId>,
    pub ambiguity_satellites: Vec<SatelliteId>,
    pub integers: Vec<i64>,
    pub validation: CandidateValidation,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SystemCandidateReferenceAudit {
    pub time: GnssTime,
    pub source_references: Vec<SatelliteId>,
    pub reference_sets: Vec<SystemReferenceCandidateValidation>,
}
impl SystemCandidateReferenceAudit {
    pub fn local_rejection_depends_on_reference(&self) -> bool {
        self.reference_sets.first().is_some_and(|first| {
            self.reference_sets
                .iter()
                .any(|v| v.validation.rejection != first.validation.rejection)
        })
    }
    /// Intersection of the unchanged local gates across every reference set.
    /// Does not perform lock, ratio, motion-history or FIX confirmation checks.
    pub fn conservative_rejection(&self) -> Option<FixDecision> {
        self.reference_sets
            .iter()
            .filter_map(|v| v.validation.rejection)
            .min_by_key(|reason| match reason {
                FixDecision::BaselineTooLong => 0,
                FixDecision::PositionCorrection => 1,
                FixDecision::CovarianceTooLarge => 2,
                FixDecision::PhaseResiduals => 3,
                FixDecision::CodeResiduals => 4,
                FixDecision::CodePhaseInconsistency => 5,
                FixDecision::PostfitInnovation => 6,
                _ => 7,
            })
    }
}

/// Re-express one supplied GPS/Galileo candidate under the Cartesian product
/// of within-system reference choices, preserving the full joint covariance.
/// Never forms a difference between constellations or runs a new LAMBDA search.
///
/// The caller bounds work with max_reference_sets (1..=4096). If the complete
/// product exceeds this limit, fail before evaluating any basis; no partial
/// audit or favourable reference subset is returned. Numerical conditioning
/// and all local thresholds are identical to validate_candidate. No filter or
/// policy history changes, and a passing audit does not declare FIX.
pub fn audit_candidate_system_references(
    float: &RtkFloatSolution,
    integers: &[i64],
    config: RtkFixConfig,
    max_reference_sets: usize,
) -> Result<SystemCandidateReferenceAudit, FixError> {
    config.validate()?;
    if !(1..=4096).contains(&max_reference_sets) {
        return Err(FixError::InvalidConfiguration);
    }
    let source = &float.fix_context;
    let references = source.reference_satellites();
    let targets = source.ambiguity_satellites();
    if integers.len() != targets.len() {
        return Err(FilterError::InvalidDimensions.into());
    }
    if integers.iter().any(|n| n.unsigned_abs() >= (1_u64 << 52)) {
        return Err(FixError::InvalidInteger);
    }
    let mut dd = BTreeMap::new();
    for reference in &references {
        dd.insert(*reference, 0_i64);
    }
    for (&satellite, &integer) in targets.iter().zip(integers) {
        dd.insert(satellite, integer);
    }
    let mut groups = BTreeMap::<_, Vec<_>>::new();
    for satellite in source.satellites() {
        groups
            .entry(satellite.system())
            .or_default()
            .push(satellite);
    }
    let mut product = 1_usize;
    for members in groups.values() {
        product = product
            .checked_mul(members.len())
            .ok_or(FixError::ReferenceSetLimit)?;
        if product > max_reference_sets {
            return Err(FixError::ReferenceSetLimit);
        }
        let min = members
            .iter()
            .map(|s| dd[s])
            .min()
            .ok_or(FilterError::InvalidDimensions)?;
        let max = members
            .iter()
            .map(|s| dd[s])
            .max()
            .ok_or(FilterError::InvalidDimensions)?;
        let span = max.checked_sub(min).ok_or(FixError::InvalidInteger)?;
        if span.unsigned_abs() >= (1_u64 << 52) {
            return Err(FixError::InvalidInteger);
        }
    }
    let mut combinations = vec![Vec::new()];
    for members in groups.values() {
        let mut next = Vec::with_capacity(combinations.len() * members.len());
        for partial in &combinations {
            for &satellite in members {
                let mut set = partial.clone();
                set.push(satellite);
                next.push(set);
            }
        }
        combinations = next;
    }
    let mut reference_sets = Vec::with_capacity(product);
    for set in combinations {
        let context = source.with_references(&set)?;
        let targets = context.ambiguity_satellites();
        let roots: BTreeMap<_, _> = set.iter().map(|s| (s.system(), dd[s])).collect();
        let equivalent: Vec<_> = targets
            .iter()
            .map(|s| {
                dd[s]
                    .checked_sub(roots[&s.system()])
                    .ok_or(FixError::InvalidInteger)
            })
            .collect::<Result<_, _>>()?;
        let validation = validate_context(&context, &equivalent, config)?;
        reference_sets.push(SystemReferenceCandidateValidation {
            reference_satellites: set,
            ambiguity_satellites: targets,
            integers: equivalent,
            validation,
        });
    }
    Ok(SystemCandidateReferenceAudit {
        time: source.time,
        source_references: references,
        reference_sets,
    })
}

#[derive(Debug, Clone, PartialEq)]
pub struct RtkFixAttempt {
    pub decision: FixDecision,
    pub ratio: Option<f64>,
    pub integers: Vec<i64>,
    pub validation: Option<CandidateValidation>,
    /// Present in AllReferences mode once local validation succeeds numerically.
    /// validation above remains the original-reference view; decision applies
    /// the intersection of all references' local gates.
    pub reference_audit: Option<CandidateReferenceAudit>,
    pub system_reference_audit: Option<SystemCandidateReferenceAudit>,
    pub motion_max_velocity_std_mps: Option<f64>,
    pub motion_checks: Vec<MotionHistoryCheck>,
    pub confirmation_count: usize,
    /// Present only after every policy and local validation check passes.
    pub fixed_solution: Option<PositionSolution>,
}
/// Scope of local candidate gates. Search, lock, ratio, confirmation and history
/// remain separate. AllReferences can reject candidates that pass Selected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceGateMode {
    Selected,
    AllReferences,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KinematicMotionConfig {
    pub max_velocity_std_mps: f64,
    pub max_interval_s: f64,
}
impl Default for KinematicMotionConfig {
    fn default() -> Self {
        Self {
            max_velocity_std_mps: 1.0,
            max_interval_s: 2.0,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RtkMotionEvidence {
    pub time: GnssTime,
    pub velocity_mps: [f64; 3],
    pub covariance_m2s2: [[f64; 3]; 3],
}
impl RtkMotionEvidence {
    pub fn from_velocity_solution(
        time: GnssTime,
        velocity: &crate::velocity::VelocitySolution,
    ) -> Self {
        Self {
            time,
            velocity_mps: velocity.velocity_mps,
            covariance_m2s2: std::array::from_fn(|i| {
                std::array::from_fn(|j| velocity.covariance_m2s2[i][j])
            }),
        }
    }
    fn validate(self, time: GnssTime) -> Result<f64, FixError> {
        if self.time.difference_seconds(time).abs() > 1e-7
            || !self.velocity_mps.iter().all(|v| v.is_finite())
        {
            return Err(FixError::InvalidMotionEvidence);
        }
        let covariance: Matrix = self.covariance_m2s2.iter().map(|r| r.to_vec()).collect();
        filter::cholesky(&covariance).map_err(|_| FixError::InvalidMotionEvidence)?;
        Ok((0..3)
            .map(|i| self.covariance_m2s2[i][i].sqrt())
            .fold(0.0, f64::max))
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MotionHistoryKind {
    LastFixed,
    PreviousCandidate,
}
#[derive(Debug, Clone, PartialEq)]
pub struct MotionHistoryCheck {
    pub kind: MotionHistoryKind,
    pub interval_s: f64,
    pub predicted_displacement_m: [f64; 3],
    pub displacement_residual_m: [f64; 3],
    pub residual_norm_m: f64,
    pub maximum_m: f64,
}
#[derive(Clone, Copy)]
struct HistoryPoint {
    position: EcefCoord,
    time: GnssTime,
    velocity: Option<[f64; 3]>,
}
struct Pending {
    key: Vec<(SatelliteId, i64)>,
    point: HistoryPoint,
    count: usize,
}
/// Conservative initial full-set AR policy. Accepted integers are checked anew
/// every epoch; no fix-and-hold or conditional state feeds back to FLOAT.
pub struct RtkFixPolicy {
    config: RtkFixConfig,
    motion_config: Option<KinematicMotionConfig>,
    reference_gate_mode: ReferenceGateMode,
    locks: BTreeMap<SatelliteId, usize>,
    last_time: Option<GnssTime>,
    pending: Option<Pending>,
    anchor: Option<HistoryPoint>,
}
impl RtkFixPolicy {
    pub fn new(config: RtkFixConfig) -> Result<Self, FixError> {
        Self::with_reference_gates(config, ReferenceGateMode::Selected)
    }
    /// Opt in to conservative validation under every reference. AllReferences
    /// performs additional conditioning and residual evaluations, with no
    /// threshold changes, new searches or feedback into FLOAT.
    pub fn with_reference_gates(
        config: RtkFixConfig,
        reference_gate_mode: ReferenceGateMode,
    ) -> Result<Self, FixError> {
        config.validate()?;
        Ok(Self {
            config,
            motion_config: None,
            reference_gate_mode,
            locks: BTreeMap::new(),
            last_time: None,
            pending: None,
            anchor: None,
        })
    }
    /// Opt-in full-set kinematic GPS/Galileo policy. Requires at least two
    /// confirmations and explicit velocity evidence; no integer holding or
    /// velocity/position fusion. Existing local thresholds stay unchanged.
    pub fn for_kinematic(
        config: RtkFixConfig,
        motion: KinematicMotionConfig,
        reference_gate_mode: ReferenceGateMode,
    ) -> Result<Self, FixError> {
        if config.min_confirmation_epochs < 2
            || [motion.max_velocity_std_mps, motion.max_interval_s]
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(FixError::InvalidConfiguration);
        }
        let mut policy = Self::with_reference_gates(config, reference_gate_mode)?;
        policy.motion_config = Some(motion);
        Ok(policy)
    }
    pub fn reset(&mut self) {
        self.locks.clear();
        self.last_time = None;
        self.pending = None;
        self.anchor = None;
    }
    fn reject(&mut self, mut attempt: RtkFixAttempt, decision: FixDecision) -> RtkFixAttempt {
        self.pending = None;
        if self.motion_config.is_some()
            && matches!(
                decision,
                FixDecision::MotionUnavailable
                    | FixDecision::MotionUncertain
                    | FixDecision::HistoryJump
            )
        {
            self.anchor = None;
        }
        attempt.decision = decision;
        attempt
    }
    fn history_age_limit(&self) -> f64 {
        self.motion_config
            .map_or(self.config.max_history_age_s, |m| {
                m.max_interval_s.min(self.config.max_history_age_s)
            })
    }
    fn history_difference(
        &self,
        current: HistoryPoint,
        previous: HistoryPoint,
        kind: MotionHistoryKind,
    ) -> (f64, Option<MotionHistoryCheck>) {
        if self.motion_config.is_none() {
            return (distance(current.position, previous.position), None);
        }
        let dt = current.time.difference_seconds(previous.time);
        let a = previous.velocity.expect("kinematic history velocity");
        let b = current.velocity.expect("admitted kinematic velocity");
        let predicted = std::array::from_fn(|i| (0.5 * a[i] + 0.5 * b[i]) * dt);
        let observed = [
            current.position.x - previous.position.x,
            current.position.y - previous.position.y,
            current.position.z - previous.position.z,
        ];
        let residual: [f64; 3] = std::array::from_fn(|i| observed[i] - predicted[i]);
        let norm = residual[0].hypot(residual[1]).hypot(residual[2]);
        (
            norm,
            Some(MotionHistoryCheck {
                kind,
                interval_s: dt,
                predicted_displacement_m: predicted,
                displacement_residual_m: residual,
                residual_norm_m: norm,
                maximum_m: self.config.max_history_jump_m,
            }),
        )
    }
    pub fn evaluate(&mut self, float: &RtkFloatSolution) -> Result<RtkFixAttempt, FixError> {
        if self.motion_config.is_some() {
            return Err(FixError::MotionEvidenceRequired);
        }
        self.evaluate_inner(float, None)
    }
    pub fn evaluate_kinematic(
        &mut self,
        float: &RtkFloatSolution,
        motion: Option<&RtkMotionEvidence>,
    ) -> Result<RtkFixAttempt, FixError> {
        if self.motion_config.is_none() {
            return Err(FixError::UnsupportedPositionMode);
        }
        self.evaluate_inner(float, motion)
    }
    fn evaluate_inner(
        &mut self,
        float: &RtkFloatSolution,
        motion: Option<&RtkMotionEvidence>,
    ) -> Result<RtkFixAttempt, FixError> {
        let context = &float.fix_context;
        if self.motion_config.is_none()
            && context.config.constellations != crate::rtk::RtkConstellations::Gps
        {
            return Err(FixError::UnsupportedConstellations);
        }
        let mode = if self.motion_config.is_some() {
            crate::rtk::RtkPositionMode::Kinematic
        } else {
            crate::rtk::RtkPositionMode::Static
        };
        if context.config.position_mode != mode {
            return Err(FixError::UnsupportedPositionMode);
        }
        let time = context.time;
        let velocity_std = motion.map(|m| m.validate(time)).transpose()?;
        if let Some(motion) = motion {
            let v = motion.velocity_mps;
            if !(v[0].hypot(v[1]).hypot(v[2]) * self.history_age_limit()).is_finite() {
                return Err(FixError::InvalidMotionEvidence);
            }
        }
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
            reference_audit: None,
            system_reference_audit: None,
            motion_max_velocity_std_mps: velocity_std,
            motion_checks: Vec::new(),
            confirmation_count: 0,
            fixed_solution: None,
        };
        if let Some(config) = self.motion_config {
            let Some(std) = velocity_std else {
                return Ok(self.reject(attempt, FixDecision::MotionUnavailable));
            };
            if std > config.max_velocity_std_mps {
                return Ok(self.reject(attempt, FixDecision::MotionUncertain));
            }
        }
        if satellites.len() - context.reference_satellites().len() < self.config.min_ambiguity_pairs
        {
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
        let local_result = match self.reference_gate_mode {
            ReferenceGateMode::Selected => {
                validate_candidate(float, &attempt.integers, self.config)
                    .map(|validation| (validation.rejection, validation))
            }
            ReferenceGateMode::AllReferences if self.motion_config.is_some() => {
                audit_candidate_system_references(float, &attempt.integers, self.config, 4096)
                    .and_then(|audit| {
                        let rejection = audit.conservative_rejection();
                        let original = audit
                            .reference_sets
                            .iter()
                            .find(|entry| entry.reference_satellites == audit.source_references)
                            .ok_or(FilterError::InvalidDimensions)?
                            .validation
                            .clone();
                        attempt.system_reference_audit = Some(audit);
                        Ok((rejection, original))
                    })
            }
            ReferenceGateMode::AllReferences => {
                audit_candidate_references(float, &attempt.integers, self.config).and_then(
                    |audit| {
                        let rejection = audit.conservative_rejection();
                        let original = audit
                            .references
                            .iter()
                            .find(|entry| entry.reference_satellite == audit.source_reference)
                            .ok_or(FilterError::InvalidDimensions)?
                            .validation
                            .clone();
                        attempt.reference_audit = Some(audit);
                        Ok((rejection, original))
                    },
                )
            }
        };
        let (rejection, validation) = match local_result {
            Ok(v) => v,
            Err(_) => return Ok(self.reject(attempt, FixDecision::ConditioningFailed)),
        };
        let fixed_position = position(context.base_position, &validation.conditioned.state.mean);
        attempt.validation = Some(validation);
        if let Some(reason) = rejection {
            return Ok(self.reject(attempt, reason));
        }
        let point = HistoryPoint {
            position: fixed_position,
            time,
            velocity: motion.map(|m| m.velocity_mps),
        };
        if let Some(anchor) = self.anchor {
            if time.difference_seconds(anchor.time) > self.history_age_limit() {
                self.anchor = None;
            } else {
                let (residual, check) =
                    self.history_difference(point, anchor, MotionHistoryKind::LastFixed);
                if let Some(check) = check {
                    attempt.motion_checks.push(check);
                }
                if residual > self.config.max_history_jump_m {
                    return Ok(self.reject(attempt, FixDecision::HistoryJump));
                }
            }
        }
        // Canonical reference-independent DD labels for repeated confirmation.
        let mut dd = BTreeMap::new();
        for reference in context.reference_satellites() {
            dd.insert(reference, 0_i64);
        }
        for (&sat, &integer) in context.ambiguity_satellites().iter().zip(&attempt.integers) {
            dd.insert(sat, integer);
        }
        let mut roots = BTreeMap::new();
        for (&sat, &n) in &dd {
            roots.entry(sat.system()).or_insert(n);
        }
        let key: Vec<_> = dd
            .into_iter()
            .map(|(s, n)| (s, n - roots[&s.system()]))
            .collect();
        let mut count = 1;
        if let Some(p) = &self.pending {
            if p.key == key && time.difference_seconds(p.point.time) <= self.history_age_limit() {
                let (residual, check) =
                    self.history_difference(point, p.point, MotionHistoryKind::PreviousCandidate);
                if let Some(check) = check {
                    attempt.motion_checks.push(check);
                }
                if residual <= self.config.max_history_jump_m {
                    count = p.count.saturating_add(1);
                } else if self.motion_config.is_some() {
                    return Ok(self.reject(attempt, FixDecision::HistoryJump));
                }
            }
        }
        self.pending = Some(Pending { key, point, count });
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
        self.anchor = Some(point);
        Ok(attempt)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rinex::{RinexObservationReader, read_navigation};
    use crate::rtk::{RtkFloatConfig, RtkFloatFilter};
    use std::io::Cursor;

    #[test]
    fn kinematic_confirmations_survive_independent_reference_changes_in_both_systems() {
        use crate::rtk::{RtkConstellations, RtkPositionMode};
        let mut r = RinexObservationReader::new(Cursor::new(include_str!(
            "../tests/fixtures/synthetic_rtk_mixed_rover.obs"
        )))
        .unwrap();
        let mut b = RinexObservationReader::new(Cursor::new(include_str!(
            "../tests/fixtures/synthetic_rtk_mixed_base.obs"
        )))
        .unwrap();
        let seed = r.header().approximate_position.unwrap();
        let base = b.header().approximate_position.unwrap();
        let nav = read_navigation(Cursor::new(include_str!(
            "../tests/fixtures/synthetic_rtk_mixed.nav"
        )))
        .unwrap();
        let truth: Vec<Vec<f64>> = include_str!("../tests/fixtures/synthetic_rtk_mixed_truth.csv")
            .lines()
            .skip(1)
            .map(|l| l.split(',').map(|v| v.parse().unwrap()).collect())
            .collect();
        let velocity = std::array::from_fn(|j| truth[1][2 + j] - truth[0][2 + j]);
        let mut engine = RtkFloatFilter::new(
            base,
            seed,
            RtkFloatConfig {
                constellations: RtkConstellations::GpsGalileo,
                position_mode: RtkPositionMode::Kinematic,
                ..Default::default()
            },
        )
        .unwrap();
        let mut original = RtkFixPolicy::for_kinematic(
            RtkFixConfig::default(),
            KinematicMotionConfig::default(),
            ReferenceGateMode::AllReferences,
        )
        .unwrap();
        let mut changed = RtkFixPolicy::for_kinematic(
            RtkFixConfig::default(),
            KinematicMotionConfig::default(),
            ReferenceGateMode::AllReferences,
        )
        .unwrap();
        for (i, (ro, ba)) in r.by_ref().zip(b.by_ref()).enumerate() {
            let ro = ro.unwrap();
            let float = engine.process_pair(&ro, &ba.unwrap(), &nav).unwrap();
            let sample = RtkMotionEvidence {
                time: ro.time,
                velocity_mps: velocity,
                covariance_m2s2: [[1e-4, 0.0, 0.0], [0.0, 1e-4, 0.0], [0.0, 0.0, 1e-4]],
            };
            let mut rebased = float.clone();
            let satellites = float.fix_context.satellites();
            let mut systems = BTreeMap::<_, Vec<_>>::new();
            for s in satellites {
                systems.entry(s.system()).or_default().push(s);
            }
            let references: Vec<_> = systems
                .values()
                .enumerate()
                .map(|(j, s)| s[(i + j) % s.len()])
                .collect();
            rebased.fix_context = float.fix_context.with_references(&references).unwrap();
            let a = original.evaluate_kinematic(&float, Some(&sample)).unwrap();
            let z = changed.evaluate_kinematic(&rebased, Some(&sample)).unwrap();
            assert_eq!(a.decision, z.decision, "epoch {i}");
            assert_eq!(a.confirmation_count, z.confirmation_count);
            if let (Some(a), Some(z)) = (a.fixed_solution, z.fixed_solution) {
                assert!(distance(a.position_ecef.unwrap(), z.position_ecef.unwrap()) < 1e-6);
            }
        }
    }

    #[test]
    fn all_reference_gates_survive_actual_snapshot_reference_changes() {
        let coordinates: Vec<f64> = include_str!("../tests/fixtures/synthetic_rtk_truth.csv")
            .lines()
            .find(|line| !line.starts_with('#'))
            .unwrap()
            .split(',')
            .map(|v| v.parse().unwrap())
            .collect();
        let coord = |i| EcefCoord::new(coordinates[i], coordinates[i + 1], coordinates[i + 2]);
        let nav = read_navigation(Cursor::new(include_str!(
            "../tests/fixtures/synthetic_rtk.nav"
        )))
        .unwrap();
        for biased in [false, true] {
            let rover = RinexObservationReader::new(Cursor::new(include_str!(
                "../tests/fixtures/synthetic_rtk_rover.obs"
            )))
            .unwrap();
            let base = RinexObservationReader::new(Cursor::new(include_str!(
                "../tests/fixtures/synthetic_rtk_base.obs"
            )))
            .unwrap();
            let mut engine =
                RtkFloatFilter::new(coord(0), coord(6), RtkFloatConfig::default()).unwrap();
            let mut last = None;
            for (i, (rover, base)) in rover.zip(base).take(11).enumerate() {
                let mut rover = rover.unwrap();
                if biased && i == 10 {
                    *rover.observations[0].carrier_phase_cycles.as_mut().unwrap() += 0.15;
                }
                last = Some(engine.process_pair(&rover, &base.unwrap(), &nav).unwrap());
            }
            let original = last.unwrap();
            let candidates = original
                .integer_candidates(LambdaConfig::default())
                .unwrap();
            for wrong in [false, true] {
                let mut integers = candidates.candidates[0].ambiguities.clone();
                if wrong {
                    integers[0] += 1;
                }
                let audit =
                    audit_candidate_references(&original, &integers, RtkFixConfig::default())
                        .unwrap();
                assert_eq!(audit.conservative_rejection().is_none(), !biased && !wrong);
                for source in &audit.references {
                    let mut rebased = original.clone();
                    rebased.fix_context = original
                        .fix_context
                        .with_reference(source.reference_satellite)
                        .unwrap();
                    let next = audit_candidate_references(
                        &rebased,
                        &source.integers,
                        RtkFixConfig::default(),
                    )
                    .unwrap();
                    assert_eq!(next.source_reference, source.reference_satellite);
                    assert_eq!(
                        next.conservative_rejection(),
                        audit.conservative_rejection()
                    );
                    for (actual, expected) in next.references.iter().zip(&audit.references) {
                        assert_eq!(actual.reference_satellite, expected.reference_satellite);
                        assert_eq!(actual.integers, expected.integers);
                        assert_eq!(actual.validation, expected.validation);
                    }
                }
            }
        }
    }
}
