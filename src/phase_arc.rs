//! Transactional IF phase-arc continuity for a future PPP estimator.
//! No ambiguity value, covariance, receiver position or FIX is estimated here.
//! GF/MW thresholds follow the ordinary pinned libgnss++ PPP state path.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.

use crate::dual_frequency::{
    DualFrequencyConfig, DualFrequencyError, IonosphereFreeObservation, PhaseUnavailable,
    SignalCorrections, form_ionosphere_free,
};
use crate::observation::ObservationEpoch;
use crate::precise::ClockSource;
use crate::{GnssTime, SatelliteId};
use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhaseArcConfig {
    pub cycle_slip_threshold_m: f64,
    /// Strictly greater gaps reset continuity; equality is admitted.
    pub max_gap_s: f64,
}
impl Default for PhaseArcConfig {
    fn default() -> Self {
        Self {
            cycle_slip_threshold_m: 0.05,
            max_gap_s: 120.0,
        }
    }
}
impl PhaseArcConfig {
    pub fn geometry_free_threshold_m(self) -> f64 {
        self.cycle_slip_threshold_m.max(0.5)
    }
    pub fn melbourne_wubbena_threshold_m(self) -> f64 {
        (self.cycle_slip_threshold_m * 100.0).max(10.0)
    }
    fn validate(self) -> Result<(), PhaseArcError> {
        if !self.cycle_slip_threshold_m.is_finite()
            || self.cycle_slip_threshold_m <= 0.0
            || !(self.cycle_slip_threshold_m * 100.0).is_finite()
            || !self.max_gap_s.is_finite()
            || self.max_gap_s <= 0.0
        {
            return Err(PhaseArcError::InvalidConfiguration);
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CombinationSlip {
    pub geometry_free_delta_m: f64,
    pub melbourne_wubbena_delta_m: Option<f64>,
    pub geometry_free: bool,
    pub melbourne_wubbena: bool,
}
/// Raw GF/MW differences with strict native thresholds. Missing MW at either
/// endpoint cannot vote for a slip. Nonfinite diagnostics are errors.
pub fn detect_combination_slip(
    previous_gf_m: f64,
    current_gf_m: f64,
    previous_mw_m: Option<f64>,
    current_mw_m: Option<f64>,
    config: PhaseArcConfig,
) -> Result<CombinationSlip, PhaseArcError> {
    config.validate()?;
    if ![previous_gf_m, current_gf_m].iter().all(|v| v.is_finite())
        || previous_mw_m
            .into_iter()
            .chain(current_mw_m)
            .any(|v| !v.is_finite())
    {
        return Err(PhaseArcError::InvalidDiagnostic);
    }
    let gf = current_gf_m - previous_gf_m;
    let mw = previous_mw_m.zip(current_mw_m).map(|(p, c)| c - p);
    if !gf.is_finite() || mw.is_some_and(|v| !v.is_finite()) {
        return Err(PhaseArcError::InvalidDiagnostic);
    }
    Ok(CombinationSlip {
        geometry_free_delta_m: gf,
        melbourne_wubbena_delta_m: mw,
        geometry_free: gf.abs() > config.geometry_free_threshold_m(),
        melbourne_wubbena: mw.is_some_and(|v| v.abs() > config.melbourne_wubbena_threshold_m()),
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PhaseArcError {
    InvalidConfiguration,
    InvalidEpoch,
    OutOfOrder,
    PendingUpdate,
    NoPendingUpdate,
    IneligibleAcceptance { satellite: SatelliteId },
    InvalidDiagnostic,
    CounterOverflow,
}
impl fmt::Display for PhaseArcError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "phase arc: {self:?}")
    }
}
impl std::error::Error for PhaseArcError {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ArcResetReason {
    FirstObservation,
    MissingSatellite,
    NotSelected,
    PhaseUnavailable(PhaseUnavailable),
    MeasurementRejected(DualFrequencyError),
    LossOfLock,
    PowerFailure,
    CycleSlipEvent,
    TimeGap,
    SignalChanged,
    ClockDatumChanged,
    PhaseDatumChanged,
    CodeCorrectionChanged,
    PhaseCorrectionChanged,
    GeometryFreeJump,
    MelbourneWubbenaJump,
    UpdateRejected,
    InvalidEpoch,
    ExternalInvalidation,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PhaseArcId(u64);
impl PhaseArcId {
    pub fn value(self) -> u64 {
        self.0
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PhaseArc {
    pub id: PhaseArcId,
    pub started_at: GnssTime,
    pub last_accepted_at: GnssTime,
    pub accepted_epochs: u64,
}
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseArcDecision {
    pub satellite: SatelliteId,
    pub measurement: Option<IonosphereFreeObservation>,
    /// Proposed ID; not committed until finish_epoch accepts this satellite.
    pub proposed_arc: Option<PhaseArcId>,
    /// Committed epochs in a continuing arc; zero for a fresh proposal.
    pub previous_accepted_epochs: u64,
    pub reset_reasons: Vec<ArcResetReason>,
    /// Only compared within an otherwise uninterrupted, unchanged signal arc.
    pub combination_slip: Option<CombinationSlip>,
}
#[derive(Debug)]
pub struct PreparedPhaseEpoch {
    time: GnssTime,
    event_only: bool,
    decisions: BTreeMap<SatelliteId, PhaseArcDecision>,
}
impl PreparedPhaseEpoch {
    pub fn time(&self) -> GnssTime {
        self.time
    }
    pub fn event_only(&self) -> bool {
        self.event_only
    }
    pub fn decisions(&self) -> &BTreeMap<SatelliteId, PhaseArcDecision> {
        &self.decisions
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct PhaseArcCommit {
    pub time: GnssTime,
    pub accepted: BTreeMap<SatelliteId, PhaseArc>,
    pub interrupted: BTreeMap<SatelliteId, Vec<ArcResetReason>>,
}

#[derive(Debug, Clone, PartialEq)]
struct Identity {
    tracking: [String; 2],
    frequency: [f64; 2],
    channel: [Option<i8>; 2],
    clock_source: ClockSource,
    clock_reference: String,
    phase_reference: String,
    code_correction: [f64; 2],
    phase_correction: [f64; 2],
}
impl Identity {
    fn from_measurement(m: &IonosphereFreeObservation) -> Self {
        Self {
            tracking: [
                m.raw[0].tracking_code.clone(),
                m.raw[1].tracking_code.clone(),
            ],
            frequency: m.frequency_hz,
            channel: [m.raw[0].glonass_channel, m.raw[1].glonass_channel],
            clock_source: m.corrections[0].clock_source,
            clock_reference: m.corrections[0].clock_reference.clone(),
            phase_reference: m.corrections[0]
                .phase_reference
                .clone()
                .expect("builder admits calibrated phase"),
            code_correction: [m.corrections[0].code_add_m, m.corrections[1].code_add_m],
            phase_correction: [
                m.corrections[0].phase_add_m.unwrap(),
                m.corrections[1].phase_add_m.unwrap(),
            ],
        }
    }
    fn changes(&self, next: &Self) -> Vec<ArcResetReason> {
        let mut r = Vec::new();
        if self.tracking != next.tracking
            || self.frequency != next.frequency
            || self.channel != next.channel
        {
            r.push(ArcResetReason::SignalChanged);
        }
        if self.clock_source != next.clock_source || self.clock_reference != next.clock_reference {
            r.push(ArcResetReason::ClockDatumChanged);
        }
        if self.phase_reference != next.phase_reference {
            r.push(ArcResetReason::PhaseDatumChanged);
        }
        if self.code_correction != next.code_correction {
            r.push(ArcResetReason::CodeCorrectionChanged);
        }
        if self.phase_correction != next.phase_correction {
            r.push(ArcResetReason::PhaseCorrectionChanged);
        }
        r
    }
}
#[derive(Debug, Clone)]
struct ActiveArc {
    arc: PhaseArc,
    identity: Identity,
    gf_m: f64,
    mw_m: Option<f64>,
}
#[derive(Debug, Clone, Default)]
struct History {
    active: Option<ActiveArc>,
    ever_accepted: bool,
    break_reasons: Vec<ArcResetReason>,
}
#[derive(Debug)]
struct Pending {
    view: PreparedPhaseEpoch,
    candidates: BTreeMap<SatelliteId, ActiveArc>,
    next_id: u64,
}

/// One selected IF pair per satellite. Prepare computes a read-only proposal;
/// finish records only phases actually accepted by the estimator. Calling
/// finish with an empty set rejects all phase updates and breaks continuity.
/// Pending proposals must be finished before another epoch can be prepared.
#[derive(Debug)]
pub struct PhaseArcTracker {
    config: PhaseArcConfig,
    last_epoch: Option<GnssTime>,
    next_id: u64,
    histories: BTreeMap<SatelliteId, History>,
    pending: Option<Pending>,
}
impl PhaseArcTracker {
    pub fn new(config: PhaseArcConfig) -> Result<Self, PhaseArcError> {
        config.validate()?;
        Ok(Self {
            config,
            last_epoch: None,
            next_id: 1,
            histories: BTreeMap::new(),
            pending: None,
        })
    }
    pub fn arc(&self, satellite: SatelliteId) -> Option<&PhaseArc> {
        self.histories
            .get(&satellite)?
            .active
            .as_ref()
            .map(|a| &a.arc)
    }
    pub fn last_epoch(&self) -> Option<GnssTime> {
        self.last_epoch
    }
    pub fn pending(&self) -> Option<&PreparedPhaseEpoch> {
        self.pending.as_ref().map(|p| &p.view)
    }
    pub fn invalidate_all(&mut self) -> Result<(), PhaseArcError> {
        if self.pending.is_some() {
            return Err(PhaseArcError::PendingUpdate);
        }
        self.break_all(ArcResetReason::ExternalInvalidation);
        Ok(())
    }
    fn break_all(&mut self, reason: ArcResetReason) {
        for h in self.histories.values_mut() {
            h.active = None;
            h.break_reasons = vec![reason.clone()];
        }
    }

    /// Normal/power-failure epochs and flag 6 event epochs are supported.
    /// Forward invalid epochs fail closed by breaking existing arcs. Out-of-order
    /// epochs and duplicate calls leave history unchanged. Global errors do not
    /// advance time; a repaired forward epoch can be retried with new arcs.
    pub fn prepare_epoch(
        &mut self,
        epoch: &ObservationEpoch,
        corrections: &SignalCorrections,
        selections: &BTreeMap<SatelliteId, DualFrequencyConfig>,
    ) -> Result<&PreparedPhaseEpoch, PhaseArcError> {
        if self.pending.is_some() {
            return Err(PhaseArcError::PendingUpdate);
        }
        if self.last_epoch.is_some_and(|t| epoch.time <= t) {
            return Err(PhaseArcError::OutOfOrder);
        }
        if !matches!(epoch.flag, 0 | 1 | 6)
            || epoch
                .receiver_clock_offset_s
                .is_some_and(|v| !v.is_finite() || v != 0.0)
        {
            self.break_all(ArcResetReason::InvalidEpoch);
            return Err(PhaseArcError::InvalidEpoch);
        }
        let satellites: BTreeSet<_> = self
            .histories
            .keys()
            .chain(selections.keys())
            .copied()
            .collect();
        let present = epoch.satellites();
        let mut next_id = self.next_id;
        let mut candidates = BTreeMap::new();
        let mut decisions = BTreeMap::new();
        for sat in satellites {
            let history = self.histories.get(&sat);
            let active = history.and_then(|h| h.active.as_ref());
            let mut reasons = Vec::new();
            if epoch.flag == 1 {
                reasons.push(ArcResetReason::PowerFailure);
            }
            let mut measurement = None;
            if epoch.flag == 6 {
                reasons.push(ArcResetReason::CycleSlipEvent);
            } else if let Some(selection) = selections.get(&sat) {
                if !present.contains(&sat) {
                    reasons.push(ArcResetReason::MissingSatellite);
                } else {
                    match form_ionosphere_free(epoch, sat, corrections, selection) {
                        Ok(m) => {
                            if let Some(reason) = &m.phase_unavailability {
                                reasons.push(ArcResetReason::PhaseUnavailable(reason.clone()));
                            }
                            measurement = Some(m);
                        }
                        Err(error) => reasons.push(ArcResetReason::MeasurementRejected(error)),
                    }
                }
            } else {
                reasons.push(ArcResetReason::NotSelected);
            }
            let mut proposed_arc = None;
            let mut previous_accepted_epochs = 0;
            let mut combination_slip = None;
            if let Some(m) = measurement.as_ref().filter(|m| m.phase.is_some()) {
                let phase = m.phase.as_ref().unwrap();
                let identity = Identity::from_measurement(m);
                if m.raw.iter().any(|o| o.loss_of_lock()) {
                    reasons.push(ArcResetReason::LossOfLock);
                }
                if let Some(old) = active {
                    reasons.extend(old.identity.changes(&identity));
                    if epoch.time.difference_seconds(old.arc.last_accepted_at)
                        > self.config.max_gap_s
                    {
                        reasons.push(ArcResetReason::TimeGap);
                    }
                    if reasons.is_empty() {
                        let slip = match detect_combination_slip(
                            old.gf_m,
                            phase.raw_geometry_free_m,
                            old.mw_m,
                            phase.raw_melbourne_wubbena_m,
                            self.config,
                        ) {
                            Ok(slip) => slip,
                            Err(error) => {
                                self.break_all(ArcResetReason::InvalidEpoch);
                                return Err(error);
                            }
                        };
                        if slip.geometry_free {
                            reasons.push(ArcResetReason::GeometryFreeJump);
                        }
                        if slip.melbourne_wubbena {
                            reasons.push(ArcResetReason::MelbourneWubbenaJump);
                        }
                        combination_slip = Some(slip);
                    }
                } else if history.is_some_and(|h| h.ever_accepted) {
                    for reason in &history.unwrap().break_reasons {
                        if !reasons.contains(reason) {
                            reasons.push(reason.clone());
                        }
                    }
                } else {
                    reasons.push(ArcResetReason::FirstObservation);
                }
                let arc = if reasons.is_empty() {
                    let old = &active.unwrap().arc;
                    previous_accepted_epochs = old.accepted_epochs;
                    PhaseArc {
                        id: old.id,
                        started_at: old.started_at,
                        last_accepted_at: epoch.time,
                        accepted_epochs: old
                            .accepted_epochs
                            .checked_add(1)
                            .ok_or(PhaseArcError::CounterOverflow)?,
                    }
                } else {
                    let id = PhaseArcId(next_id);
                    next_id = next_id
                        .checked_add(1)
                        .ok_or(PhaseArcError::CounterOverflow)?;
                    PhaseArc {
                        id,
                        started_at: epoch.time,
                        last_accepted_at: epoch.time,
                        accepted_epochs: 1,
                    }
                };
                proposed_arc = Some(arc.id);
                candidates.insert(
                    sat,
                    ActiveArc {
                        arc,
                        identity,
                        gf_m: phase.raw_geometry_free_m,
                        mw_m: phase.raw_melbourne_wubbena_m,
                    },
                );
            }
            decisions.insert(
                sat,
                PhaseArcDecision {
                    satellite: sat,
                    measurement,
                    proposed_arc,
                    previous_accepted_epochs,
                    reset_reasons: reasons,
                    combination_slip,
                },
            );
        }
        self.pending = Some(Pending {
            view: PreparedPhaseEpoch {
                time: epoch.time,
                event_only: epoch.flag == 6,
                decisions,
            },
            candidates,
            next_id,
        });
        Ok(&self.pending.as_ref().unwrap().view)
    }

    /// Commit only eligible satellites whose phase update actually succeeded.
    /// Code-only acceptance must not be included. Invalid acceptance leaves the
    /// pending proposal intact for correction. IDs are never reused after finish,
    /// even if a proposed phase was rejected; numeric estimator state is external.
    pub fn finish_epoch(
        &mut self,
        accepted: &BTreeSet<SatelliteId>,
    ) -> Result<PhaseArcCommit, PhaseArcError> {
        let pending = self
            .pending
            .as_ref()
            .ok_or(PhaseArcError::NoPendingUpdate)?;
        for sat in accepted {
            if !pending.candidates.contains_key(sat) {
                return Err(PhaseArcError::IneligibleAcceptance { satellite: *sat });
            }
        }
        let mut pending = self.pending.take().unwrap();
        let mut committed = BTreeMap::new();
        let mut interrupted = BTreeMap::new();
        for (sat, decision) in &pending.view.decisions {
            let h = self.histories.entry(*sat).or_default();
            if accepted.contains(sat) {
                let candidate = pending.candidates.remove(sat).unwrap();
                committed.insert(*sat, candidate.arc.clone());
                h.active = Some(candidate);
                h.ever_accepted = true;
                h.break_reasons.clear();
            } else {
                let mut reasons = decision.reset_reasons.clone();
                if decision.proposed_arc.is_some() {
                    reasons.push(ArcResetReason::UpdateRejected);
                }
                if reasons.is_empty() {
                    reasons.push(ArcResetReason::UpdateRejected);
                }
                h.active = None;
                h.break_reasons = reasons.clone();
                interrupted.insert(*sat, reasons);
            }
        }
        self.last_epoch = Some(pending.view.time);
        self.next_id = pending.next_id;
        Ok(PhaseArcCommit {
            time: pending.view.time,
            accepted: committed,
            interrupted,
        })
    }
}
