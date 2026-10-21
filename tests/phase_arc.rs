use gnss_rust::dual_frequency::*;
use gnss_rust::observation::{Observation, ObservationEpoch};
use gnss_rust::phase_arc::*;
use gnss_rust::precise::ClockSource;
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::{GnssTime, SatelliteId, constants as c};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Cursor,
};

fn sat(prn: u8) -> SatelliteId {
    format!("G{prn:02}").parse().unwrap()
}
fn selections() -> BTreeMap<SatelliteId, DualFrequencyConfig> {
    (1..=2)
        .map(|prn| {
            (
                sat(prn),
                DualFrequencyConfig {
                    tracking_codes: ["1C".to_owned(), "2W".to_owned()],
                    clock_source: ClockSource::RinexClk,
                    clock_reference: "clock-v1".to_owned(),
                    phase_reference: "phase-v1".to_owned(),
                    noise: DualFrequencyNoise {
                        code_m2: [[0.25, 0.0], [0.0, 0.25]],
                        phase_m2: [[0.0001, 0.0], [0.0, 0.0001]],
                    },
                },
            )
        })
        .collect()
}
fn entries(selections: &BTreeMap<SatelliteId, DualFrequencyConfig>) -> Vec<SignalCorrection> {
    selections
        .iter()
        .flat_map(|(&sat, config)| {
            config
                .tracking_codes
                .iter()
                .map(move |code| SignalCorrection {
                    satellite: sat,
                    tracking_code: code.clone(),
                    clock_source: config.clock_source,
                    clock_reference: config.clock_reference.clone(),
                    phase_reference: Some(config.phase_reference.clone()),
                    code_add_m: 0.0,
                    phase_add_m: Some(0.0),
                    valid_from: GnssTime::new(2300, 0.0).unwrap(),
                    valid_until: GnssTime::new(2301, 1000.0).unwrap(),
                })
        })
        .collect()
}
fn table(s: &BTreeMap<SatelliteId, DualFrequencyConfig>) -> SignalCorrections {
    SignalCorrections::new(entries(s)).unwrap()
}
fn epoch(seconds: f64) -> ObservationEpoch {
    let observations = (1..=2)
        .flat_map(|prn| {
            (0..2).map(move |i| {
                let mut o = Observation::new(sat(prn), if i == 0 { "1C" } else { "2W" }.to_owned());
                let f = if i == 0 {
                    c::GPS_L1_FREQ
                } else {
                    c::GPS_L2_FREQ
                };
                let range = 23000000.0 + f64::from(prn) * 1000.0 + seconds * 5.0;
                o.pseudorange_m = Some(range);
                o.carrier_phase_cycles =
                    Some(range * f / c::SPEED_OF_LIGHT + if i == 0 { 12.0 } else { -7.0 });
                o.lli = Some(0);
                o
            })
        })
        .collect();
    ObservationEpoch {
        time: GnssTime::new(2300, 345600.0 + seconds).unwrap(),
        flag: 0,
        approximate_position: None,
        receiver_clock_offset_s: None,
        observations,
    }
}
fn tracker() -> PhaseArcTracker {
    PhaseArcTracker::new(PhaseArcConfig::default()).unwrap()
}
fn all_eligible(t: &PhaseArcTracker) -> BTreeSet<SatelliteId> {
    t.pending()
        .unwrap()
        .decisions()
        .iter()
        .filter_map(|(&sat, d)| d.proposed_arc.map(|_| sat))
        .collect()
}
fn accept_all(t: &mut PhaseArcTracker) -> PhaseArcCommit {
    t.finish_epoch(&all_eligible(t)).unwrap()
}
fn init(t: &mut PhaseArcTracker, s: &BTreeMap<SatelliteId, DualFrequencyConfig>) {
    t.prepare_epoch(&epoch(0.0), &table(s), s).unwrap();
    accept_all(t);
}
fn decision(t: &PhaseArcTracker, prn: u8) -> &PhaseArcDecision {
    &t.pending().unwrap().decisions()[&sat(prn)]
}

#[test]
fn rejected_retention_preserves_accepted_counts_times_and_gf_mw_anchors() {
    let s = selections();
    let mut t = tracker();
    init(&mut t, &s);
    let old = t.arc(sat(1)).unwrap().clone();
    let mut rejected = epoch(1.0);
    *rejected.observations[0]
        .carrier_phase_cycles
        .as_mut()
        .unwrap() += 0.3 * c::GPS_L1_FREQ / c::SPEED_OF_LIGHT;
    t.prepare_epoch(&rejected, &table(&s), &s).unwrap();
    let r = t.reject_epoch_retaining_continuous().unwrap();
    assert_eq!(r.retained.len(), 2);
    assert!(r.interrupted.is_empty());
    assert_eq!(t.arc(sat(1)), Some(&old));
    assert_eq!(t.last_epoch(), Some(rejected.time));
    let mut next = epoch(2.0);
    *next.observations[0].carrier_phase_cycles.as_mut().unwrap() +=
        0.6 * c::GPS_L1_FREQ / c::SPEED_OF_LIGHT;
    t.prepare_epoch(&next, &table(&s), &s).unwrap();
    assert!(
        decision(&t, 1)
            .reset_reasons
            .contains(&ArcResetReason::GeometryFreeJump)
    );
    assert!(
        (decision(&t, 1)
            .combination_slip
            .unwrap()
            .geometry_free_delta_m
            - 0.6)
            .abs()
            < 1e-7
    );
    assert_eq!(decision(&t, 2).previous_accepted_epochs, 1);
    let accepted = accept_all(&mut t);
    assert_eq!(accepted.accepted[&sat(2)].accepted_epochs, 2);
}

#[test]
fn repeated_rejected_retention_cannot_extend_the_last_accepted_gap_limit() {
    let s = selections();
    let mut t = tracker();
    init(&mut t, &s);
    let old = t.arc(sat(1)).unwrap().clone();
    for time in [60.0, 119.0, 120.0] {
        t.prepare_epoch(&epoch(time), &table(&s), &s).unwrap();
        assert_eq!(
            t.reject_epoch_retaining_continuous()
                .unwrap()
                .retained
                .len(),
            2
        );
        assert_eq!(t.arc(sat(1)), Some(&old));
    }
    t.prepare_epoch(&epoch(120.001), &table(&s), &s).unwrap();
    let r = t.reject_epoch_retaining_continuous().unwrap();
    assert!(r.retained.is_empty());
    assert!(r.interrupted[&sat(1)].contains(&ArcResetReason::TimeGap));
    assert!(t.arc(sat(1)).is_none());
}

#[test]
fn intrinsic_discontinuities_interrupt_even_when_the_estimator_rejects_the_epoch() {
    for case in 0..12 {
        let mut s = selections();
        let mut t = tracker();
        init(&mut t, &s);
        let mut e = epoch(1.0);
        let reason = match case {
            0 => {
                e.observations[0].lli = Some(1);
                ArcResetReason::LossOfLock
            }
            1 => {
                e.observations[0].lli = Some(2);
                ArcResetReason::PhaseUnavailable(PhaseUnavailable::HalfCycle {
                    tracking_code: "1C".into(),
                })
            }
            2 => {
                e.observations[0].carrier_phase_cycles = None;
                ArcResetReason::PhaseUnavailable(PhaseUnavailable::MissingPhase {
                    tracking_code: "1C".into(),
                })
            }
            3 => {
                e.observations.retain(|o| o.satellite != sat(1));
                ArcResetReason::MissingSatellite
            }
            4 => {
                s.remove(&sat(1));
                ArcResetReason::NotSelected
            }
            5 => {
                e.flag = 1;
                ArcResetReason::PowerFailure
            }
            6 => {
                e.flag = 6;
                ArcResetReason::CycleSlipEvent
            }
            7 => {
                s.get_mut(&sat(1)).unwrap().clock_reference = "changed-clock".into();
                ArcResetReason::ClockDatumChanged
            }
            8 => {
                s.get_mut(&sat(1)).unwrap().phase_reference = "changed-phase".into();
                ArcResetReason::PhaseDatumChanged
            }
            9 => {
                *e.observations[0].carrier_phase_cycles.as_mut().unwrap() +=
                    c::GPS_L1_FREQ / c::SPEED_OF_LIGHT;
                ArcResetReason::GeometryFreeJump
            }
            10 => {
                *e.observations[0].pseudorange_m.as_mut().unwrap() += 1000.0;
                ArcResetReason::MelbourneWubbenaJump
            }
            11 => {
                e.observations[1].tracking_code = "2L".into();
                s.get_mut(&sat(1)).unwrap().tracking_codes[1] = "2L".into();
                ArcResetReason::SignalChanged
            }
            _ => unreachable!(),
        };
        t.prepare_epoch(&e, &table(&s), &s).unwrap();
        let r = t.reject_epoch_retaining_continuous().unwrap();
        assert!(
            r.interrupted[&sat(1)].contains(&reason),
            "case {case}: {r:?}"
        );
        assert!(!r.retained.contains_key(&sat(1)));
        assert!(t.arc(sat(1)).is_none());
        assert_eq!(r.retained.contains_key(&sat(2)), !matches!(case, 5 | 6));
    }
}

#[test]
fn rejected_new_arcs_are_never_accepted_and_default_rejection_still_interrupts() {
    let s = selections();
    let mut t = tracker();
    assert_eq!(
        t.reject_epoch_retaining_continuous(),
        Err(PhaseArcError::NoPendingUpdate)
    );
    t.prepare_epoch(&epoch(0.0), &table(&s), &s).unwrap();
    let proposed = decision(&t, 1).proposed_arc.unwrap();
    assert!(
        t.reject_epoch_retaining_continuous()
            .unwrap()
            .retained
            .is_empty()
    );
    assert!(t.arc(sat(1)).is_none());
    t.prepare_epoch(&epoch(1.0), &table(&s), &s).unwrap();
    assert!(decision(&t, 1).proposed_arc.unwrap() > proposed);
    accept_all(&mut t);
    t.prepare_epoch(&epoch(2.0), &table(&s), &s).unwrap();
    assert!(
        t.finish_epoch(&BTreeSet::new())
            .unwrap()
            .accepted
            .is_empty()
    );
    assert!(t.arc(sat(1)).is_none());
}

#[test]
fn native_48_predicates_match_floor_and_strict_positive_negative_boundaries() {
    let mut count = 0;
    for line in include_str!("fixtures/upstream_phase_arc.csv")
        .lines()
        .skip(1)
    {
        let v: Vec<_> = line.split(',').collect();
        let n = |i: usize| v[i].parse::<f64>().unwrap();
        let config = PhaseArcConfig {
            cycle_slip_threshold_m: n(0),
            ..Default::default()
        };
        assert_eq!(config.geometry_free_threshold_m(), n(1));
        assert_eq!(config.melbourne_wubbena_threshold_m(), n(2));
        let result = detect_combination_slip(
            n(3),
            n(4),
            if v[7] == "1" { Some(n(5)) } else { None },
            if v[8] == "1" { Some(n(6)) } else { None },
            config,
        )
        .unwrap();
        assert_eq!(result.geometry_free, v[9] == "1");
        assert_eq!(result.melbourne_wubbena, v[10] == "1");
        assert_eq!(result.geometry_free_delta_m, n(4) - n(3));
        count += 1;
    }
    assert_eq!(count, 48);
}

#[test]
fn native_generated_continuous_rinex_trace_handles_all_lifecycle_events() {
    let epochs = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_phase_arc.obs"
    )))
    .unwrap()
    .collect::<Result<Vec<_>, _>>()
    .unwrap();
    assert_eq!(epochs.len(), 16);
    let s = selections();
    let table = table(&s);
    let mut t = tracker();
    let mut accepted = 0;
    let mut g2_previous = None;
    for (i, e) in epochs.iter().enumerate() {
        t.prepare_epoch(e, &table, &s).unwrap();
        let reason = match i {
            0 => Some(ArcResetReason::FirstObservation),
            2 => Some(ArcResetReason::LossOfLock),
            5 => Some(ArcResetReason::PhaseUnavailable(
                PhaseUnavailable::MissingPhase {
                    tracking_code: "2W".to_owned(),
                },
            )),
            6 => Some(ArcResetReason::GeometryFreeJump),
            9 => Some(ArcResetReason::MissingSatellite),
            10 => Some(ArcResetReason::PowerFailure),
            12 => Some(ArcResetReason::CycleSlipEvent),
            13 => Some(ArcResetReason::CycleSlipEvent),
            14 => Some(ArcResetReason::TimeGap),
            15 => Some(ArcResetReason::MelbourneWubbenaJump),
            _ => None,
        };
        if let Some(reason) = reason {
            assert!(
                decision(&t, 1).reset_reasons.contains(&reason),
                "epoch {i}: {:?}",
                decision(&t, 1).reset_reasons
            );
        }
        if [4, 8, 11, 12].contains(&i) {
            assert!(decision(&t, 1).proposed_arc.is_none());
        }
        if i == 4 || i == 11 {
            assert!(
                decision(&t, 1)
                    .measurement
                    .as_ref()
                    .unwrap()
                    .phase
                    .is_none()
            );
            assert!(
                decision(&t, 1)
                    .measurement
                    .as_ref()
                    .unwrap()
                    .code_if_m
                    .is_finite()
            );
        }
        if i == 12 {
            assert!(t.pending().unwrap().event_only());
            assert!(
                decision(&t, 2)
                    .reset_reasons
                    .contains(&ArcResetReason::CycleSlipEvent)
            );
        }
        if i == 6 {
            let slip = decision(&t, 1).combination_slip.unwrap();
            assert!(slip.geometry_free);
            assert!(!slip.melbourne_wubbena);
        }
        if i == 15 {
            let slip = decision(&t, 1).combination_slip.unwrap();
            assert!(!slip.geometry_free);
            assert!(slip.melbourne_wubbena);
        }
        let mut eligible = all_eligible(&t);
        if i == 7 {
            eligible.remove(&sat(1));
        }
        let committed = t.finish_epoch(&eligible).unwrap();
        accepted += committed.accepted.len();
        if i == 7 {
            assert!(committed.interrupted[&sat(1)].contains(&ArcResetReason::UpdateRejected));
            assert!(t.arc(sat(1)).is_none());
        }
        if let Some(g2) = committed.accepted.get(&sat(2)) {
            if i > 0 && ![10, 13, 14].contains(&i) {
                assert_eq!(Some(g2.id), g2_previous, "G02 continuity at {i}");
            }
            g2_previous = Some(g2.id);
        }
    }
    assert_eq!(accepted, 26);
    assert_eq!(t.arc(sat(1)).unwrap().id.value(), 13);
    assert_eq!(t.arc(sat(2)).unwrap().accepted_epochs, 2);
}

#[test]
fn proposals_do_not_commit_or_inflate_lock_and_invalid_acceptance_is_retryable() {
    let mut t = tracker();
    let s = selections();
    let table = table(&s);
    t.prepare_epoch(&epoch(0.0), &table, &s).unwrap();
    assert!(t.arc(sat(1)).is_none());
    assert_eq!(decision(&t, 1).previous_accepted_epochs, 0);
    assert_eq!(
        t.prepare_epoch(&epoch(30.0), &table, &s).unwrap_err(),
        PhaseArcError::PendingUpdate
    );
    assert_eq!(t.invalidate_all(), Err(PhaseArcError::PendingUpdate));
    assert_eq!(
        t.finish_epoch(&BTreeSet::from([sat(3)])),
        Err(PhaseArcError::IneligibleAcceptance { satellite: sat(3) })
    );
    assert!(t.pending().is_some());
    accept_all(&mut t);
    assert_eq!(t.arc(sat(1)).unwrap().accepted_epochs, 1);
    assert_eq!(
        t.finish_epoch(&BTreeSet::new()),
        Err(PhaseArcError::NoPendingUpdate)
    );
    t.prepare_epoch(&epoch(30.0), &table, &s).unwrap();
    assert_eq!(decision(&t, 1).previous_accepted_epochs, 1);
    assert!(decision(&t, 1).reset_reasons.is_empty());
    assert_eq!(t.arc(sat(1)).unwrap().accepted_epochs, 1);
    accept_all(&mut t);
    assert_eq!(t.arc(sat(1)).unwrap().accepted_epochs, 2);
}

#[test]
fn failed_update_breaks_continuity_and_never_reuses_a_rejected_new_arc_id() {
    let mut t = tracker();
    let s = selections();
    init(&mut t, &s);
    let original = t.arc(sat(1)).unwrap().id;
    t.prepare_epoch(&epoch(30.0), &table(&s), &s).unwrap();
    t.finish_epoch(&BTreeSet::new()).unwrap();
    assert!(t.arc(sat(1)).is_none());
    assert_eq!(t.last_epoch(), Some(epoch(30.0).time));
    t.prepare_epoch(&epoch(60.0), &table(&s), &s).unwrap();
    let rejected = decision(&t, 1).proposed_arc.unwrap();
    assert_ne!(rejected, original);
    assert!(
        decision(&t, 1)
            .reset_reasons
            .contains(&ArcResetReason::UpdateRejected)
    );
    t.finish_epoch(&BTreeSet::new()).unwrap();
    t.prepare_epoch(&epoch(90.0), &table(&s), &s).unwrap();
    let final_id = decision(&t, 1).proposed_arc.unwrap();
    assert!(final_id > rejected);
    accept_all(&mut t);
    assert_eq!(t.arc(sat(1)).unwrap().accepted_epochs, 1);
}

#[test]
fn subset_acceptance_breaks_only_rejected_satellite() {
    let mut t = tracker();
    let s = selections();
    init(&mut t, &s);
    let old = t.arc(sat(2)).unwrap().id;
    t.prepare_epoch(&epoch(30.0), &table(&s), &s).unwrap();
    t.finish_epoch(&BTreeSet::from([sat(2)])).unwrap();
    t.prepare_epoch(&epoch(60.0), &table(&s), &s).unwrap();
    assert!(
        decision(&t, 1)
            .reset_reasons
            .contains(&ArcResetReason::UpdateRejected)
    );
    assert_eq!(decision(&t, 2).proposed_arc, Some(old));
    assert!(decision(&t, 2).reset_reasons.is_empty());
    accept_all(&mut t);
    assert_eq!(t.arc(sat(2)).unwrap().accepted_epochs, 3);
}

#[test]
fn either_signals_lli_and_receiver_power_failure_reset_without_rejecting_phase() {
    for index in 0..2 {
        let mut t = tracker();
        let s = selections();
        init(&mut t, &s);
        let mut e = epoch(30.0);
        e.observations[index].lli = Some(1);
        t.prepare_epoch(&e, &table(&s), &s).unwrap();
        assert!(
            decision(&t, 1)
                .reset_reasons
                .contains(&ArcResetReason::LossOfLock)
        );
        assert!(decision(&t, 2).reset_reasons.is_empty());
        assert!(
            decision(&t, 1)
                .measurement
                .as_ref()
                .unwrap()
                .phase
                .is_some()
        );
        accept_all(&mut t);
        assert_eq!(t.arc(sat(1)).unwrap().accepted_epochs, 1);
    }
    let mut t = tracker();
    let s = selections();
    init(&mut t, &s);
    let mut e = epoch(30.0);
    e.flag = 1;
    t.prepare_epoch(&e, &table(&s), &s).unwrap();
    for prn in 1..=2 {
        assert!(
            decision(&t, prn)
                .reset_reasons
                .contains(&ArcResetReason::PowerFailure)
        );
    }
    accept_all(&mut t);
    for prn in 1..=2 {
        assert_eq!(t.arc(sat(prn)).unwrap().accepted_epochs, 1);
    }
}

#[test]
fn missing_phase_half_cycle_missing_satellite_and_removed_selection_break_arcs() {
    let s = selections();
    for kind in 0..4 {
        let mut t = tracker();
        init(&mut t, &s);
        let old = t.arc(sat(1)).unwrap().id;
        let mut e = epoch(30.0);
        let mut chosen = s.clone();
        match kind {
            0 => e.observations[1].carrier_phase_cycles = None,
            1 => e.observations[1].lli = Some(2),
            2 => e.observations.retain(|o| o.satellite != sat(1)),
            _ => {
                chosen.remove(&sat(1));
            }
        }
        t.prepare_epoch(&e, &table(&s), &chosen).unwrap();
        assert!(decision(&t, 1).proposed_arc.is_none());
        assert!(decision(&t, 1).reset_reasons.len() == 1);
        if kind < 2 {
            assert!(decision(&t, 1).measurement.is_some());
            assert!(
                decision(&t, 1)
                    .measurement
                    .as_ref()
                    .unwrap()
                    .phase
                    .is_none()
            );
        }
        accept_all(&mut t);
        assert!(t.arc(sat(1)).is_none());
        t.prepare_epoch(&epoch(60.0), &table(&s), &s).unwrap();
        assert_ne!(decision(&t, 1).proposed_arc, Some(old));
        assert!(
            !decision(&t, 1)
                .reset_reasons
                .contains(&ArcResetReason::FirstObservation)
        );
        accept_all(&mut t);
        assert_eq!(t.arc(sat(1)).unwrap().accepted_epochs, 1);
    }
}

#[test]
fn rejected_measurement_is_isolated_and_does_not_accept_code_as_phase() {
    let mut t = tracker();
    let s = selections();
    init(&mut t, &s);
    let mut e = epoch(30.0);
    e.observations[0].pseudorange_m = Some(f64::NAN);
    t.prepare_epoch(&e, &table(&s), &s).unwrap();
    assert!(matches!(
        decision(&t, 1).reset_reasons.as_slice(),
        [ArcResetReason::MeasurementRejected(
            DualFrequencyError::InvalidObservation { .. }
        )]
    ));
    assert_eq!(
        t.finish_epoch(&BTreeSet::from([sat(1)])),
        Err(PhaseArcError::IneligibleAcceptance { satellite: sat(1) })
    );
    accept_all(&mut t);
    assert!(t.arc(sat(1)).is_none());
    assert_eq!(t.arc(sat(2)).unwrap().accepted_epochs, 2);
}

#[test]
fn signal_clock_phase_datum_and_correction_changes_reset_before_diagnostic_comparison() {
    for kind in 0..6 {
        let mut t = tracker();
        let s = selections();
        init(&mut t, &s);
        let mut chosen = s.clone();
        let mut e = epoch(30.0);
        match kind {
            0 => {
                chosen.get_mut(&sat(1)).unwrap().tracking_codes[0] = "1W".to_owned();
                e.observations[0].tracking_code = "1W".to_owned();
            }
            1 => chosen.get_mut(&sat(1)).unwrap().clock_reference = "clock-v2".to_owned(),
            2 => chosen.get_mut(&sat(1)).unwrap().clock_source = ClockSource::Sp3,
            3 => chosen.get_mut(&sat(1)).unwrap().phase_reference = "phase-v2".to_owned(),
            _ => {}
        }
        let mut correction = entries(&chosen);
        if kind == 4 {
            correction[0].code_add_m = 1.0;
        }
        if kind == 5 {
            correction[1].phase_add_m = Some(0.01);
        }
        let expected = match kind {
            0 => ArcResetReason::SignalChanged,
            1 | 2 => ArcResetReason::ClockDatumChanged,
            3 => ArcResetReason::PhaseDatumChanged,
            4 => ArcResetReason::CodeCorrectionChanged,
            _ => ArcResetReason::PhaseCorrectionChanged,
        };
        t.prepare_epoch(&e, &SignalCorrections::new(correction).unwrap(), &chosen)
            .unwrap();
        assert!(decision(&t, 1).reset_reasons.contains(&expected));
        assert!(decision(&t, 1).combination_slip.is_none());
        assert!(decision(&t, 2).reset_reasons.is_empty());
        accept_all(&mut t);
        assert_eq!(t.arc(sat(1)).unwrap().accepted_epochs, 1);
    }
}

#[test]
fn validity_extension_and_noise_change_preserve_an_unchanged_phase_datum() {
    let mut t = tracker();
    let s = selections();
    init(&mut t, &s);
    let old = t.arc(sat(1)).unwrap().id;
    let mut changed = s.clone();
    changed.get_mut(&sat(1)).unwrap().noise.code_m2 = [[1.0, 0.2], [0.2, 1.0]];
    let mut correction = entries(&changed);
    for row in &mut correction {
        row.valid_from = epoch(30.0).time;
        row.valid_until = epoch(90.0).time;
    }
    t.prepare_epoch(
        &epoch(30.0),
        &SignalCorrections::new(correction).unwrap(),
        &changed,
    )
    .unwrap();
    assert_eq!(decision(&t, 1).proposed_arc, Some(old));
    assert!(decision(&t, 1).reset_reasons.is_empty());
    accept_all(&mut t);
}

#[test]
fn strict_time_gap_and_full_week_time_are_used() {
    for gap in [120.0, 120.000001] {
        let mut t = tracker();
        let s = selections();
        let mut first = epoch(0.0);
        first.time = GnssTime::new(2300, 604740.0).unwrap();
        t.prepare_epoch(&first, &table(&s), &s).unwrap();
        accept_all(&mut t);
        let old = t.arc(sat(1)).unwrap().id;
        let mut next = epoch(gap);
        next.time = first.time.checked_add_seconds(gap).unwrap();
        t.prepare_epoch(&next, &table(&s), &s).unwrap();
        if gap == 120.0 {
            assert_eq!(decision(&t, 1).proposed_arc, Some(old));
            assert!(decision(&t, 1).reset_reasons.is_empty());
        } else {
            assert!(
                decision(&t, 1)
                    .reset_reasons
                    .contains(&ArcResetReason::TimeGap)
            );
            assert!(decision(&t, 1).combination_slip.is_none());
        }
        accept_all(&mut t);
    }
}

#[test]
fn duplicate_or_backward_epochs_do_not_modify_committed_history() {
    let mut t = tracker();
    let s = selections();
    init(&mut t, &s);
    let before = t.arc(sat(1)).unwrap().clone();
    for seconds in [0.0, -1.0] {
        assert_eq!(
            t.prepare_epoch(&epoch(seconds), &table(&s), &s)
                .unwrap_err(),
            PhaseArcError::OutOfOrder
        );
        assert_eq!(t.arc(sat(1)), Some(&before));
        assert!(t.pending().is_none());
    }
}

#[test]
fn forward_invalid_epoch_fails_closed_and_can_be_repaired_at_same_time() {
    for kind in 0..3 {
        let mut t = tracker();
        let s = selections();
        init(&mut t, &s);
        let mut bad = epoch(30.0);
        match kind {
            0 => bad.flag = 2,
            1 => bad.receiver_clock_offset_s = Some(1e-6),
            _ => bad.receiver_clock_offset_s = Some(f64::NAN),
        }
        assert_eq!(
            t.prepare_epoch(&bad, &table(&s), &s).unwrap_err(),
            PhaseArcError::InvalidEpoch
        );
        assert!(t.arc(sat(1)).is_none());
        assert_eq!(t.last_epoch(), Some(epoch(0.0).time));
        t.prepare_epoch(&epoch(30.0), &table(&s), &s).unwrap();
        assert!(
            decision(&t, 1)
                .reset_reasons
                .contains(&ArcResetReason::InvalidEpoch)
        );
        accept_all(&mut t);
    }
}

#[test]
fn sparse_flag6_event_interrupts_even_satellites_absent_from_event_records() {
    let mut t = tracker();
    let s = selections();
    init(&mut t, &s);
    let mut event = epoch(30.0);
    event.flag = 6;
    event.observations.retain(|o| o.satellite == sat(1));
    t.prepare_epoch(&event, &table(&s), &s).unwrap();
    assert!(all_eligible(&t).is_empty());
    let committed = accept_all(&mut t);
    for prn in 1..=2 {
        assert!(committed.interrupted[&sat(prn)].contains(&ArcResetReason::CycleSlipEvent));
        assert!(t.arc(sat(prn)).is_none());
    }
    t.prepare_epoch(&epoch(60.0), &table(&s), &s).unwrap();
    for prn in 1..=2 {
        assert!(
            decision(&t, prn)
                .reset_reasons
                .contains(&ArcResetReason::CycleSlipEvent)
        );
    }
    accept_all(&mut t);
}

#[test]
fn glonass_channel_change_resets_and_mw_does_not_vote() {
    let r: SatelliteId = "R01".parse().unwrap();
    let mut s = selections();
    let mut config = s.remove(&sat(1)).unwrap();
    config.tracking_codes[1] = "2C".to_owned();
    s = BTreeMap::from([(r, config)]);
    let mut first = epoch(0.0);
    first.observations.retain(|o| o.satellite == sat(1));
    for (i, o) in first.observations.iter_mut().enumerate() {
        o.satellite = r;
        o.glonass_channel = Some(0);
        o.tracking_code = if i == 0 { "1C" } else { "2C" }.to_owned();
    }
    let mut t = tracker();
    t.prepare_epoch(&first, &table(&s), &s).unwrap();
    accept_all(&mut t);
    let mut next = first.clone();
    next.time = next.time.checked_add_seconds(30.0).unwrap();
    next.observations[0].pseudorange_m = Some(30000000.0);
    t.prepare_epoch(&next, &table(&s), &s).unwrap();
    let d = &t.pending().unwrap().decisions()[&r];
    assert!(!d.combination_slip.unwrap().melbourne_wubbena);
    assert!(d.reset_reasons.is_empty());
    accept_all(&mut t);
    next.time = next.time.checked_add_seconds(30.0).unwrap();
    for o in &mut next.observations {
        o.glonass_channel = Some(1);
    }
    t.prepare_epoch(&next, &table(&s), &s).unwrap();
    let d = &t.pending().unwrap().decisions()[&r];
    assert!(d.reset_reasons.contains(&ArcResetReason::SignalChanged));
    assert!(d.combination_slip.is_none());
    accept_all(&mut t);
}

#[test]
fn external_invalidation_requires_fresh_arcs() {
    let mut t = tracker();
    let s = selections();
    init(&mut t, &s);
    t.invalidate_all().unwrap();
    assert!(t.arc(sat(1)).is_none());
    t.prepare_epoch(&epoch(30.0), &table(&s), &s).unwrap();
    assert!(
        decision(&t, 1)
            .reset_reasons
            .contains(&ArcResetReason::ExternalInvalidation)
    );
    accept_all(&mut t);
    assert_eq!(t.arc(sat(1)).unwrap().accepted_epochs, 1);
}

#[test]
fn nonfinite_diagnostics_and_invalid_configuration_are_errors() {
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(matches!(
            PhaseArcTracker::new(PhaseArcConfig {
                cycle_slip_threshold_m: value,
                ..Default::default()
            }),
            Err(PhaseArcError::InvalidConfiguration)
        ));
    }
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            PhaseArcTracker::new(PhaseArcConfig {
                max_gap_s: value,
                ..Default::default()
            }),
            Err(PhaseArcError::InvalidConfiguration)
        ));
    }
    assert_eq!(
        detect_combination_slip(f64::NAN, 0.0, None, None, PhaseArcConfig::default()),
        Err(PhaseArcError::InvalidDiagnostic)
    );
    assert_eq!(
        detect_combination_slip(
            0.0,
            0.0,
            Some(f64::INFINITY),
            None,
            PhaseArcConfig::default()
        ),
        Err(PhaseArcError::InvalidDiagnostic)
    );
    assert_eq!(
        detect_combination_slip(-f64::MAX, f64::MAX, None, None, PhaseArcConfig::default()),
        Err(PhaseArcError::InvalidDiagnostic)
    );
}

#[test]
fn quality_exclusion_breaks_only_named_arcs_without_advancing_other_rejected_history() {
    let s = selections();
    let b = table(&s);
    let mut t = PhaseArcTracker::new(Default::default()).unwrap();
    t.prepare_epoch(&epoch(0.0), &b, &s).unwrap();
    t.finish_epoch(&BTreeSet::from([sat(1), sat(2)])).unwrap();
    let one = t.arc(sat(1)).unwrap().clone();
    let two = t.arc(sat(2)).unwrap().clone();
    t.prepare_epoch(&epoch(1.0), &b, &s).unwrap();
    let r = t
        .reject_epoch_retaining_continuous_excluding(&BTreeSet::from([sat(1)]))
        .unwrap();
    assert_eq!(r.retained[&sat(2)], two);
    assert!(!r.retained.contains_key(&sat(1)));
    assert!(t.arc(sat(1)).is_none());
    assert_eq!(t.arc(sat(2)), Some(&two));
    let p = t.prepare_epoch(&epoch(2.0), &b, &s).unwrap();
    assert_ne!(p.decisions()[&sat(1)].proposed_arc, Some(one.id));
    assert_eq!(p.decisions()[&sat(2)].proposed_arc, Some(two.id));
}
