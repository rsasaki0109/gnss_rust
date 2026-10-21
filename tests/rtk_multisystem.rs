use gnss_rust::observation::ObservationEpoch;
use gnss_rust::rinex::{RinexObservationReader, read_navigation};
use gnss_rust::rtk::{
    RtkConstellations, RtkError, RtkFloatConfig, RtkFloatFilter, RtkPositionMode, RtkRejection,
};
use gnss_rust::{EcefCoord, GnssSystem, SatelliteId};
use std::io::Cursor;

fn data() -> (
    Vec<ObservationEpoch>,
    Vec<ObservationEpoch>,
    gnss_rust::navigation::NavigationData,
    EcefCoord,
    EcefCoord,
) {
    let mut r = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_mixed_rover.obs"
    )))
    .unwrap();
    let mut b = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_mixed_base.obs"
    )))
    .unwrap();
    let seed = r.header().approximate_position.unwrap();
    let base = b.header().approximate_position.unwrap();
    let rover = r.by_ref().collect::<Result<Vec<_>, _>>().unwrap();
    let bases = b.by_ref().collect::<Result<Vec<_>, _>>().unwrap();
    let nav = read_navigation(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_mixed.nav"
    )))
    .unwrap();
    (rover, bases, nav, base, seed)
}
fn config() -> RtkFloatConfig {
    RtkFloatConfig {
        constellations: RtkConstellations::GpsGalileo,
        position_mode: RtkPositionMode::Kinematic,
        ..Default::default()
    }
}
fn truth(index: usize) -> EcefCoord {
    let v: Vec<f64> = include_str!("fixtures/synthetic_rtk_mixed_truth.csv")
        .lines()
        .nth(index + 1)
        .unwrap()
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect();
    EcefCoord::new(v[2], v[3], v[4])
}
fn distance(a: EcefCoord, b: EcefCoord) -> f64 {
    (a.x - b.x).hypot(a.y - b.y).hypot(a.z - b.z)
}
fn gal_target(epoch: &ObservationEpoch) -> SatelliteId {
    epoch
        .observations
        .iter()
        .find(|o| o.satellite.system() == GnssSystem::Galileo)
        .unwrap()
        .satellite
}

fn assert_native_component_startup(result: &gnss_rust::rtk::RtkFloatSolution, name: &str) {
    let line = include_str!("fixtures/upstream_rtk_mixed_startup.csv")
        .lines()
        .find(|l| l.starts_with(&format!("{name},")))
        .unwrap();
    let v: Vec<f64> = line
        .split(',')
        .skip(2)
        .map(|s| s.parse().unwrap())
        .collect();
    assert!(
        distance(
            result.position.position_ecef.unwrap(),
            EcefCoord::new(v[0], v[1], v[2])
        ) < 1e-6,
        "{name}: native startup position"
    );
    let covariance = result.position.position_covariance.unwrap();
    for i in 0..3 {
        for j in 0..3 {
            let expected = v[3 + 3 * i + j];
            assert!(
                (covariance[i][j] - expected).abs() < 1e-10 * expected.abs().max(1.0),
                "{name}: startup Q[{i},{j}]"
            );
        }
    }
}

#[test]
fn mixed_moving_truth_uses_separate_references_and_clock_geometry() {
    let (r, b, nav, base, seed) = data();
    let mut filter = RtkFloatFilter::new(base, seed, config()).unwrap();
    for i in 0..r.len() {
        let result = filter.process_pair(&r[i], &b[i], &nav).unwrap();
        if i == 0 {
            assert_native_component_startup(&result, "full");
        }
        assert!(
            distance(result.position.position_ecef.unwrap(), truth(i)) < 0.003,
            "epoch {i}"
        );
        assert_eq!(result.position.num_satellites, 8);
        assert_eq!(
            result
                .reference_satellites
                .iter()
                .map(|s| s.system())
                .collect::<Vec<_>>(),
            vec![GnssSystem::Gps, GnssSystem::Galileo]
        );
        assert_eq!(result.reference_satellite, result.reference_satellites[0]);
        assert_eq!(result.ambiguity_satellites.len(), 6);
        assert_eq!(result.reset_ambiguities.len(), if i == 0 { 8 } else { 0 });
        let geometry = result.geometry.as_ref().unwrap();
        assert_eq!(
            geometry.clock_systems,
            vec![GnssSystem::Gps, GnssSystem::Galileo]
        );
        assert_eq!(geometry.full_geometry_covariance.as_ref().unwrap().len(), 5);
        if i == 31 {
            let candidate = result.integer_candidates(Default::default()).unwrap();
            let expected: Vec<i64> = result
                .ambiguity_satellites
                .iter()
                .map(|s| {
                    let reference = result
                        .reference_satellites
                        .iter()
                        .find(|r| r.system() == s.system())
                        .unwrap();
                    i64::from(reference.prn()) - i64::from(s.prn())
                })
                .collect();
            assert_eq!(candidate.candidates[0].ambiguities, expected);
        }
    }
}

#[test]
fn three_gps_and_three_galileo_recover_position_with_independent_clocks() {
    let (mut r, mut b, nav, base, seed) = data();
    for e in r.iter_mut().chain(&mut b) {
        let mut counts = std::collections::BTreeMap::new();
        e.observations.retain(|o| {
            let n = counts.entry(o.satellite.system()).or_insert(0);
            *n += 1;
            *n <= 3
        });
    }
    let mut gps = RtkFloatFilter::new(
        base,
        seed,
        RtkFloatConfig {
            position_mode: RtkPositionMode::Kinematic,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        gps.process_pair(&r[0], &b[0], &nav),
        Err(RtkError::InsufficientCommonSatellites)
    );
    let mut mixed = RtkFloatFilter::new(base, seed, config()).unwrap();
    for i in 0..r.len() {
        let result = mixed.process_pair(&r[i], &b[i], &nav).unwrap();
        assert_eq!(result.position.num_satellites, 6);
        assert_eq!(result.ambiguity_satellites.len(), 4);
        if i == 0 {
            // The finite initial position prior produces a 3.865 mm startup
            // bias also present in the independent native component result.
            assert_native_component_startup(&result, "three_each");
        } else {
            assert!(distance(result.position.position_ecef.unwrap(), truth(i)) < 0.003);
        }
    }
}

#[test]
fn a_common_galileo_phase_offset_cancels_inside_its_dd_block() {
    let (mut r, b, nav, base, seed) = data();
    let mut control = RtkFloatFilter::new(base, seed, config()).unwrap();
    let mut shifted = RtkFloatFilter::new(base, seed, config()).unwrap();
    for i in 0..3 {
        let expected = control.process_pair(&r[i], &b[i], &nav).unwrap();
        for o in &mut r[i].observations {
            if o.satellite.system() == GnssSystem::Galileo {
                *o.carrier_phase_cycles.as_mut().unwrap() += 0.125;
            }
        }
        let result = shifted.process_pair(&r[i], &b[i], &nav).unwrap();
        assert!(
            distance(
                expected.position.position_ecef.unwrap(),
                result.position.position_ecef.unwrap()
            ) < 1e-6
        );
    }
}

#[test]
fn galileo_reference_loss_and_reacquisition_preserve_other_sd_arcs() {
    let (r, b, nav, base, seed) = data();
    let mut filter = RtkFloatFilter::new(base, seed, config()).unwrap();
    let first = filter.process_pair(&r[0], &b[0], &nav).unwrap();
    let target = first.reference_satellites[1];
    let mut missing = r[1].clone();
    missing.observations.retain(|o| o.satellite != target);
    let second = filter.process_pair(&missing, &b[1], &nav).unwrap();
    assert_ne!(second.reference_satellites[1], target);
    assert!(second.reset_ambiguities.is_empty());
    assert_eq!(second.reference_satellite, first.reference_satellite);
    let third = filter.process_pair(&r[2], &b[2], &nav).unwrap();
    assert_eq!(third.reset_ambiguities, vec![target]);
    assert!(distance(third.position.position_ecef.unwrap(), truth(2)) < 0.003);
}

#[test]
fn tracking_code_changes_start_new_arcs_without_lli() {
    let (mut r, b, nav, base, seed) = data();
    let target = gal_target(&r[0]);
    let mut filter = RtkFloatFilter::new(base, seed, config()).unwrap();
    filter.process_pair(&r[0], &b[0], &nav).unwrap();
    for i in 1..3 {
        let o = r[i]
            .observations
            .iter_mut()
            .find(|o| o.satellite == target)
            .unwrap();
        o.tracking_code = "1X".into();
        *o.carrier_phase_cycles.as_mut().unwrap() += 0.4;
        let result = filter.process_pair(&r[i], &b[i], &nav).unwrap();
        assert_eq!(
            result.reset_ambiguities,
            if i == 1 { vec![target] } else { vec![] }
        );
        let pair = result
            .tracking_pairs()
            .into_iter()
            .find(|p| p.satellite == target)
            .unwrap();
        assert_eq!(
            (pair.rover_tracking_code, pair.base_tracking_code),
            ("1X", "1X")
        );
        assert!(distance(result.position.position_ecef.unwrap(), truth(i)) < 0.003);
    }
}

#[test]
fn singleton_constellations_are_reported_and_not_counted_as_measurements() {
    let (r, b, nav, base, seed) = data();
    let target = gal_target(&r[0]);
    let mut initial = r[0].clone();
    initial
        .observations
        .retain(|o| o.satellite.system() == GnssSystem::Gps || o.satellite == target);
    let mut filter = RtkFloatFilter::new(base, seed, config()).unwrap();
    let first = filter.process_pair(&initial, &b[0], &nav).unwrap();
    assert_eq!(first.position.num_satellites, 4);
    assert_eq!(first.reference_satellites.len(), 1);
    assert!(
        first
            .rejected
            .contains(&(target, RtkRejection::InsufficientSystemSatellites))
    );
    let second = filter.process_pair(&r[1], &b[1], &nav).unwrap();
    assert_eq!(second.reset_ambiguities.len(), 4);
    assert!(
        second
            .reset_ambiguities
            .iter()
            .all(|s| s.system() == GnssSystem::Galileo)
    );
}

#[test]
fn satellite_count_alone_does_not_replace_the_number_of_clock_differences() {
    let (mut r, b, nav, base, seed) = data();
    let mut counts = std::collections::BTreeMap::new();
    r[0].observations.retain(|o| {
        let n = counts.entry(o.satellite.system()).or_insert(0);
        *n += 1;
        *n <= 2
    });
    let mut filter = RtkFloatFilter::new(base, seed, config()).unwrap();
    assert_eq!(
        filter.process_pair(&r[0], &b[0], &nav),
        Err(RtkError::InsufficientCommonSatellites)
    );
    assert!(filter.state().is_none());
    let mut gal = r[1].clone();
    gal.observations
        .retain(|o| o.satellite.system() == GnssSystem::Galileo);
    let result = filter.process_pair(&gal, &b[1], &nav).unwrap();
    assert_eq!(result.reference_satellites.len(), 1);
    assert_eq!(result.reference_satellite.system(), GnssSystem::Galileo);
    // Four Galileo satellites retain a measurable finite-prior startup bias;
    // compare it to native components rather than assuming zero prior weight.
    assert_native_component_startup(&result, "galileo_only");
}

#[test]
fn initial_fix_policy_rejects_mixed_constellations_even_with_local_candidate_validation() {
    use gnss_rust::rtk_fix::{
        FixError, RtkFixConfig, RtkFixPolicy, audit_candidate_references, validate_candidate,
    };
    let (r, b, nav, base, seed) = data();
    let mut mixed = RtkFloatFilter::new(
        base,
        seed,
        RtkFloatConfig {
            position_mode: RtkPositionMode::Static,
            ..config()
        },
    )
    .unwrap();
    let result = mixed.process_pair(&r[0], &b[0], &nav).unwrap();
    let mut policy = RtkFixPolicy::new(RtkFixConfig::default()).unwrap();
    assert_eq!(
        policy.evaluate(&result),
        Err(FixError::UnsupportedConstellations)
    );
    let integers = vec![0; result.ambiguity_satellites.len()];
    assert!(validate_candidate(&result, &integers, RtkFixConfig::default()).is_ok());
    assert_eq!(
        audit_candidate_references(&result, &integers, RtkFixConfig::default()),
        Err(FixError::UnsupportedConstellations)
    );
    let mut gps = RtkFloatFilter::new(base, seed, Default::default()).unwrap();
    let gps = gps.process_pair(&r[0], &b[0], &nav).unwrap();
    assert!(policy.evaluate(&gps).is_ok());
}

#[test]
fn complete_galileo_tracking_is_preferred_and_duplicate_groups_remain_invalid() {
    let (r, b, nav, base, seed) = data();
    let target = gal_target(&r[0]);
    let mut modified = r[0].clone();
    let mut alternate = modified
        .observations
        .iter()
        .find(|o| o.satellite == target)
        .unwrap()
        .clone();
    alternate.tracking_code = "1X".into();
    modified.observations.push(alternate);
    let mut control = RtkFloatFilter::new(base, seed, config()).unwrap();
    control.process_pair(&r[0], &b[0], &nav).unwrap();
    let mut preferred = RtkFloatFilter::new(base, seed, config()).unwrap();
    let first = preferred.process_pair(&modified, &b[0], &nav).unwrap();
    assert_eq!(preferred.state(), control.state());
    assert_eq!(
        first
            .tracking_pairs()
            .into_iter()
            .find(|p| p.satellite == target)
            .unwrap()
            .rover_tracking_code,
        "1C"
    );
    modified
        .observations
        .iter_mut()
        .find(|o| o.satellite == target && o.tracking_code == "1C")
        .unwrap()
        .pseudorange_m = None;
    let mut fallback = RtkFloatFilter::new(base, seed, config()).unwrap();
    let result = fallback.process_pair(&modified, &b[0], &nav).unwrap();
    assert_eq!(
        result
            .tracking_pairs()
            .into_iter()
            .find(|p| p.satellite == target)
            .unwrap()
            .rover_tracking_code,
        "1X"
    );
    modified.observations.push(
        modified
            .observations
            .iter()
            .find(|o| o.satellite == target && o.tracking_code == "1X")
            .unwrap()
            .clone(),
    );
    let before = fallback.state().unwrap().clone();
    assert_eq!(
        fallback.process_pair(&modified, &b[0], &nav),
        Err(RtkError::OutOfOrder)
    );
    modified.time = modified.time.checked_add_seconds(1.0).unwrap();
    assert_eq!(
        fallback.process_pair(&modified, &b[1], &nav),
        Err(RtkError::DuplicateObservation)
    );
    assert_eq!(fallback.state(), Some(&before));
}

#[test]
fn skipped_galileo_lli_is_continuity_evidence_only_for_enabled_systems() {
    let (r, b, nav, base, seed) = data();
    let mut filter = RtkFloatFilter::new(base, seed, config()).unwrap();
    filter.process_pair(&r[0], &b[0], &nav).unwrap();
    let before = filter.state().unwrap().clone();
    let mut skipped = r[0].clone();
    skipped.time = skipped.time.checked_add_seconds(0.2).unwrap();
    skipped
        .observations
        .retain(|o| o.satellite.system() == GnssSystem::Galileo);
    skipped.observations[0].lli = Some(1);
    assert!(filter.observe_skipped_epoch(&skipped));
    assert_eq!(filter.state(), Some(&before));
    let result = filter.process_pair(&r[1], &b[1], &nav).unwrap();
    assert_eq!(result.reset_ambiguities.len(), 8);
    let mut gps = RtkFloatFilter::new(base, seed, Default::default()).unwrap();
    assert!(!gps.observe_skipped_epoch(&skipped));
}

#[test]
fn galileo_unflagged_slip_resets_the_galileo_satellite_only() {
    use gnss_rust::rtk_slip::RtkSlipConfig;
    let (mut r, b, nav, base, seed) = data();
    let target = gal_target(&r[0]);
    let mut filter = RtkFloatFilter::new(
        base,
        seed,
        RtkFloatConfig {
            slip_detection: Some(RtkSlipConfig::default()),
            ..config()
        },
    )
    .unwrap();
    filter.process_pair(&r[0], &b[0], &nav).unwrap();
    *r[1]
        .observations
        .iter_mut()
        .find(|o| o.satellite == target)
        .unwrap()
        .carrier_phase_cycles
        .as_mut()
        .unwrap() += 50.0;
    let result = filter.process_pair(&r[1], &b[1], &nav).unwrap();
    assert_eq!(result.reset_ambiguities, vec![target]);
    assert!(
        result
            .slip_diagnostics
            .iter()
            .find(|d| d.satellite == target)
            .unwrap()
            .code
            .unwrap()
            .detected
    );
    assert!(distance(result.position.position_ecef.unwrap(), truth(1)) < 0.003);
}

#[test]
fn mixed_prefit_rows_map_to_references_in_the_same_constellation() {
    let (mut r, b, nav, base, seed) = data();
    let target = gal_target(&r[0]);
    let mut filter = RtkFloatFilter::new(
        base,
        seed,
        RtkFloatConfig {
            max_prefit_residual_m: Some(3.0),
            ..config()
        },
    )
    .unwrap();
    filter.process_pair(&r[0], &b[0], &nav).unwrap();
    let before = filter.state().unwrap().clone();
    *r[1]
        .observations
        .iter_mut()
        .find(|o| o.satellite == target)
        .unwrap()
        .pseudorange_m
        .as_mut()
        .unwrap() += 200.0;
    assert_eq!(
        filter.process_pair(&r[1], &b[1], &nav),
        Err(RtkError::PrefitResidualRejected)
    );
    assert_eq!(filter.state(), Some(&before));
    let diagnostic = filter.last_prefit_diagnostic().unwrap();
    assert_eq!(diagnostic.row_satellites.len(), 6);
    assert!(
        diagnostic
            .row_satellites
            .iter()
            .zip(&diagnostic.row_reference_satellites)
            .all(|(s, r)| s.system() == r.system())
    );
    assert!(
        diagnostic
            .screen
            .rejected_code_rows
            .iter()
            .all(|i| diagnostic.row_satellites[*i].system() == GnssSystem::Galileo)
    );
}

#[test]
fn datum_alignment_survives_reacquisition_with_independent_system_phase_offsets() {
    let (r, b, nav, base, seed) = data();
    let cfg = RtkFloatConfig {
        align_ambiguity_datums: true,
        ..config()
    };
    let mut control = RtkFloatFilter::new(base, seed, cfg).unwrap();
    let mut shifted = RtkFloatFilter::new(base, seed, cfg).unwrap();
    let gps = r[0]
        .observations
        .iter()
        .find(|o| o.satellite.system() == GnssSystem::Gps)
        .unwrap()
        .satellite;
    let gal = gal_target(&r[0]);
    for i in 0..r.len() {
        let mut plain = r[i].clone();
        if i == 1 {
            plain
                .observations
                .retain(|o| o.satellite != gps && o.satellite != gal);
        }
        let mut moved = plain.clone();
        for o in &mut moved.observations {
            let cycles = if o.satellite.system() == GnssSystem::Gps {
                1.25
            } else {
                -2.5
            };
            *o.carrier_phase_cycles.as_mut().unwrap() += cycles * i as f64;
        }
        let expected = control.process_pair(&plain, &b[i], &nav).unwrap();
        let actual = shifted.process_pair(&moved, &b[i], &nav).unwrap();
        assert!(
            distance(
                actual.position.position_ecef.unwrap(),
                expected.position.position_ecef.unwrap()
            ) < 1e-6,
            "epoch {i}"
        );
        assert_eq!(actual.reset_ambiguities, expected.reset_ambiguities);
        let a = actual.position.position_covariance.unwrap();
        let e = expected.position.position_covariance.unwrap();
        for j in 0..3 {
            for k in 0..3 {
                assert!((a[j][k] - e[j][k]).abs() < 1e-10 * e[j][k].abs().max(1.0));
            }
        }
        if i == 2 {
            for d in &actual.ambiguity_datum_alignment {
                assert_eq!(d.retained_satellites.len(), 3);
                assert_eq!(d.initialized_satellites.len(), 1);
            }
        }
    }
}

#[test]
fn datum_alignment_never_commits_a_failed_update_and_global_recovery_has_no_old_arcs() {
    let (mut r, b, nav, base, seed) = data();
    let mut filter = RtkFloatFilter::new(
        base,
        seed,
        RtkFloatConfig {
            align_ambiguity_datums: true,
            max_prefit_residual_m: Some(3.0),
            ..config()
        },
    )
    .unwrap();
    filter.process_pair(&r[0], &b[0], &nav).unwrap();
    let before = filter.state().unwrap().clone();
    let target = gal_target(&r[0]);
    let recovery = r[1].clone();
    for o in &mut r[1].observations {
        if o.satellite.system() == GnssSystem::Galileo {
            *o.carrier_phase_cycles.as_mut().unwrap() += 12.5;
        }
        if o.satellite == target {
            *o.pseudorange_m.as_mut().unwrap() += 200.0;
        }
    }
    assert_eq!(
        filter.process_pair(&r[1], &b[1], &nav),
        Err(RtkError::PrefitResidualRejected)
    );
    assert_eq!(filter.state(), Some(&before));
    // A failed attempt has not advanced time. Retry its clean measurements;
    // two moving epochs later can legitimately exceed this same 3 m gate.
    let result = filter.process_pair(&recovery, &b[1], &nav).unwrap();
    assert_eq!(result.reset_ambiguities.len(), 8);
    for d in result.ambiguity_datum_alignment {
        assert!(d.retained_satellites.is_empty());
        assert_eq!(d.initialized_satellites.len(), 4);
        assert_eq!(d.offset_cycles, 0.0);
    }
}

#[test]
fn a_tracking_change_is_excluded_from_the_datum_mean() {
    let (mut r, b, nav, base, seed) = data();
    let target = gal_target(&r[0]);
    let mut filter = RtkFloatFilter::new(
        base,
        seed,
        RtkFloatConfig {
            align_ambiguity_datums: true,
            ..config()
        },
    )
    .unwrap();
    filter.process_pair(&r[0], &b[0], &nav).unwrap();
    let observation = r[1]
        .observations
        .iter_mut()
        .find(|o| o.satellite == target)
        .unwrap();
    observation.tracking_code = "1X".into();
    *observation.carrier_phase_cycles.as_mut().unwrap() += 0.25;
    let result = filter.process_pair(&r[1], &b[1], &nav).unwrap();
    assert_eq!(result.reset_ambiguities, vec![target]);
    let d = result
        .ambiguity_datum_alignment
        .iter()
        .find(|d| d.system == GnssSystem::Galileo)
        .unwrap();
    assert!(!d.retained_satellites.contains(&target));
    assert_eq!(d.initialized_satellites, vec![target]);
    assert_eq!(d.retained_satellites.len(), 3);
}

#[test]
fn aligned_zero_baseline_retains_zero_ambiguities_as_active_states() {
    let (_, b, nav, base, _) = data();
    let mut filter = RtkFloatFilter::new(
        base,
        base,
        RtkFloatConfig {
            align_ambiguity_datums: true,
            ..config()
        },
    )
    .unwrap();
    for (i, epoch) in b.iter().enumerate() {
        let result = filter.process_pair(epoch, epoch, &nav).unwrap();
        assert_eq!(result.position.position_ecef, Some(base));
        assert!(filter.state().unwrap().mean.iter().all(|v| *v == 0.0));
        assert_eq!(result.reset_ambiguities.len(), if i == 0 { 8 } else { 0 });
        for d in &result.ambiguity_datum_alignment {
            assert_eq!(d.retained_satellites.len(), if i == 0 { 0 } else { 4 });
            assert_eq!(d.offset_cycles, 0.0);
        }
    }
}

fn known_integers(float: &gnss_rust::rtk::RtkFloatSolution) -> Vec<i64> {
    float
        .ambiguity_satellites
        .iter()
        .map(|s| {
            let reference = float
                .reference_satellites
                .iter()
                .find(|r| r.system() == s.system())
                .unwrap();
            i64::from(reference.prn()) - i64::from(s.prn())
        })
        .collect()
}

#[test]
fn mixed_known_integers_validate_all_32_moving_snapshots_without_declaring_fix() {
    use gnss_rust::rtk_fix::{RtkFixConfig, validate_candidate};
    let (r, b, nav, base, seed) = data();
    let mut filter = RtkFloatFilter::new(base, seed, config()).unwrap();
    for i in 0..r.len() {
        let float = filter.process_pair(&r[i], &b[i], &nav).unwrap();
        let before = filter.state().unwrap().clone();
        let result =
            validate_candidate(&float, &known_integers(&float), RtkFixConfig::default()).unwrap();
        assert_eq!(result.rejection, None, "epoch {i}");
        let fixed = EcefCoord::new(
            base.x + result.conditioned.state.mean[0],
            base.y + result.conditioned.state.mean[1],
            base.z + result.conditioned.state.mean[2],
        );
        assert!(distance(fixed, truth(i)) < 0.003, "epoch {i}");
        assert_eq!(float.position.status, gnss_rust::SolutionStatus::Float);
        assert_eq!(filter.state(), Some(&before));
    }
}

#[test]
fn mixed_audit_checks_all_16_reference_sets_with_the_same_physical_integers() {
    use gnss_rust::rtk_fix::{RtkFixConfig, audit_candidate_system_references};
    let (r, b, nav, base, seed) = data();
    let mut filter = RtkFloatFilter::new(base, seed, config()).unwrap();
    let float = filter.process_pair(&r[0], &b[0], &nav).unwrap();
    let before = filter.state().unwrap().clone();
    let audit = audit_candidate_system_references(
        &float,
        &known_integers(&float),
        RtkFixConfig::default(),
        16,
    )
    .unwrap();
    assert_eq!(audit.reference_sets.len(), 16);
    assert_eq!(audit.conservative_rejection(), None);
    assert!(!audit.local_rejection_depends_on_reference());
    let first = &audit.reference_sets[0].validation.conditioned;
    for entry in &audit.reference_sets {
        assert_eq!(entry.reference_satellites.len(), 2);
        for (&target, &integer) in entry.ambiguity_satellites.iter().zip(&entry.integers) {
            let reference = entry
                .reference_satellites
                .iter()
                .find(|r| r.system() == target.system())
                .unwrap();
            assert_eq!(
                integer,
                i64::from(reference.prn()) - i64::from(target.prn())
            );
        }
        for i in 0..3 {
            assert!(
                (entry.validation.conditioned.state.mean[i] - first.state.mean[i]).abs() < 1e-6
            );
            for j in 0..3 {
                assert!(
                    (entry.validation.conditioned.state.covariance[i][j]
                        - first.state.covariance[i][j])
                        .abs()
                        < 1e-10
                );
            }
        }
    }
    assert_eq!(filter.state(), Some(&before));
}

#[test]
fn wrong_gps_or_galileo_integers_are_rejected_under_all_reference_combinations() {
    use gnss_rust::rtk_fix::{RtkFixConfig, audit_candidate_system_references, validate_candidate};
    let (r, b, nav, base, seed) = data();
    let mut filter = RtkFloatFilter::new(base, seed, config()).unwrap();
    let float = filter.process_pair(&r[0], &b[0], &nav).unwrap();
    let before = filter.state().unwrap().clone();
    for system in [GnssSystem::Gps, GnssSystem::Galileo] {
        let mut integers = known_integers(&float);
        let i = float
            .ambiguity_satellites
            .iter()
            .position(|s| s.system() == system)
            .unwrap();
        integers[i] += 1;
        assert!(
            validate_candidate(&float, &integers, RtkFixConfig::default())
                .unwrap()
                .rejection
                .is_some()
        );
        let audit =
            audit_candidate_system_references(&float, &integers, RtkFixConfig::default(), 16)
                .unwrap();
        assert!(
            audit
                .reference_sets
                .iter()
                .all(|v| v.validation.rejection.is_some())
        );
    }
    assert_eq!(filter.state(), Some(&before));
}

#[test]
fn mixed_audit_checks_work_limits_and_integer_ranges_before_returning_any_basis() {
    use gnss_rust::filter::FilterError;
    use gnss_rust::rtk_fix::{FixError, RtkFixConfig, audit_candidate_system_references as audit};
    let (r, b, nav, base, seed) = data();
    let mut filter = RtkFloatFilter::new(base, seed, config()).unwrap();
    let float = filter.process_pair(&r[0], &b[0], &nav).unwrap();
    let integers = known_integers(&float);
    let cfg = RtkFixConfig::default();
    assert_eq!(
        audit(&float, &integers, cfg, 15),
        Err(FixError::ReferenceSetLimit)
    );
    assert_eq!(
        audit(&float, &integers, cfg, 0),
        Err(FixError::InvalidConfiguration)
    );
    assert_eq!(
        audit(&float, &integers, cfg, 4097),
        Err(FixError::InvalidConfiguration)
    );
    assert_eq!(
        audit(&float, &integers[..5], cfg, 16),
        Err(FixError::Filter(FilterError::InvalidDimensions))
    );
    let mut invalid = integers.clone();
    invalid[0] = 1_i64 << 52;
    assert_eq!(
        audit(&float, &invalid, cfg, 16),
        Err(FixError::InvalidInteger)
    );
    invalid[0] = -(1_i64 << 51);
    invalid[1] = 1_i64 << 51;
    assert_eq!(
        audit(&float, &invalid, cfg, 16),
        Err(FixError::InvalidInteger)
    );
}

#[test]
fn single_system_audit_is_identical_to_the_existing_gps_audit() {
    use gnss_rust::rtk_fix::{
        RtkFixConfig, audit_candidate_references, audit_candidate_system_references,
    };
    let (r, b, nav, base, seed) = data();
    let mut filter = RtkFloatFilter::new(base, seed, RtkFloatConfig::default()).unwrap();
    let float = filter.process_pair(&r[0], &b[0], &nav).unwrap();
    let integers = known_integers(&float);
    let cfg = RtkFixConfig::default();
    let old = audit_candidate_references(&float, &integers, cfg).unwrap();
    let new = audit_candidate_system_references(&float, &integers, cfg, 4).unwrap();
    assert_eq!(old.conservative_rejection(), new.conservative_rejection());
    assert_eq!(old.references.len(), new.reference_sets.len());
    for (a, b) in old.references.iter().zip(new.reference_sets) {
        assert_eq!(vec![a.reference_satellite], b.reference_satellites);
        assert_eq!(a.ambiguity_satellites, b.ambiguity_satellites);
        assert_eq!(a.integers, b.integers);
        assert_eq!(a.validation, b.validation);
    }
}

#[test]
fn mixed_validation_and_reference_audit_use_immutable_snapshot_identity() {
    use gnss_rust::rtk_fix::{
        RtkFixConfig, audit_candidate_system_references as audit, validate_candidate,
    };
    let (r, b, nav, base, seed) = data();
    let mut filter = RtkFloatFilter::new(base, seed, config()).unwrap();
    let float = filter.process_pair(&r[0], &b[0], &nav).unwrap();
    let integers = known_integers(&float);
    let cfg = RtkFixConfig::default();
    let expected = validate_candidate(&float, &integers, cfg).unwrap();
    let expected_audit = audit(&float, &integers, cfg, 16).unwrap();
    let mut edited = float.clone();
    edited.reference_satellites.reverse();
    edited.reference_satellite = "E36".parse().unwrap();
    edited.ambiguity_satellites.reverse();
    edited.ambiguities.dd_float_cycles.fill(42.0);
    edited.position.position_ecef = Some(EcefCoord::new(1.0, 2.0, 3.0));
    assert_eq!(
        validate_candidate(&edited, &integers, cfg).unwrap(),
        expected
    );
    assert_eq!(audit(&edited, &integers, cfg, 16).unwrap(), expected_audit);
}
