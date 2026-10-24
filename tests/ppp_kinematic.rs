use gnss_rust::bias::{BiasSignalRequest, OsbBinding, read_bias_sinex};
use gnss_rust::dual_frequency::{DualFrequencyConfig, DualFrequencyNoise};
use gnss_rust::ppp::*;
use gnss_rust::precise::{ClockSource, OrbitReferencePoint, PreciseProducts, read_clk, read_sp3};
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::tides::ReceiverPositionReference;
use gnss_rust::{EcefCoord, SatelliteId, SolutionStatus};
use std::{collections::BTreeMap, io::Cursor};

fn products() -> PreciseProducts {
    PreciseProducts::new(
        read_sp3(Cursor::new(include_str!("fixtures/synthetic_ppp.sp3"))).unwrap(),
        Some(read_clk(Cursor::new(include_str!("fixtures/synthetic_ppp.clk"))).unwrap()),
    )
}
fn epochs(noisy: bool) -> Vec<gnss_rust::observation::ObservationEpoch> {
    RinexObservationReader::new(Cursor::new(if noisy {
        include_str!("fixtures/synthetic_ppp_kinematic_noisy.obs")
    } else {
        include_str!("fixtures/synthetic_ppp_kinematic.obs")
    }))
    .unwrap()
    .collect::<Result<_, _>>()
    .unwrap()
}
fn truth() -> Vec<Vec<f64>> {
    include_str!("fixtures/synthetic_ppp_kinematic_truth.csv")
        .lines()
        .skip(1)
        .map(|l| l.split(',').map(|v| v.parse().unwrap()).collect())
        .collect()
}
fn selections() -> BTreeMap<SatelliteId, DualFrequencyConfig> {
    [1, 2, 7, 8, 9, 18, 19, 25, 26]
        .into_iter()
        .map(|n| {
            (
                format!("G{n:02}").parse().unwrap(),
                DualFrequencyConfig {
                    tracking_codes: ["1C".into(), "2W".into()],
                    clock_source: ClockSource::RinexClk,
                    clock_reference: "synthetic-precise-clock-v1".into(),
                    phase_reference: "synthetic-phase-v1".into(),
                    noise: DualFrequencyNoise {
                        code_m2: [[0.0025, 0.0], [0.0, 0.0025]],
                        phase_m2: [[0.000001, 0.0], [0.0, 0.000001]],
                    },
                },
            )
        })
        .collect()
}
fn biases(t: gnss_rust::GnssTime) -> gnss_rust::dual_frequency::SignalCorrections {
    let p = read_bias_sinex(
        Cursor::new(include_str!("fixtures/synthetic_ppp_osb.bsx")),
        Default::default(),
    )
    .unwrap();
    let requests: Vec<_> = selections()
        .into_iter()
        .flat_map(|(sat, c)| {
            c.tracking_codes
                .into_iter()
                .map(move |code| BiasSignalRequest {
                    satellite: sat,
                    tracking_code: code,
                    glonass_channel: None,
                })
        })
        .collect();
    p.corrections_at(
        t,
        &OsbBinding {
            product_id: "synthetic_ppp_osb.bsx".into(),
            clock_source: ClockSource::RinexClk,
            clock_reference: "synthetic-precise-clock-v1".into(),
            phase_reference: Some("synthetic-phase-v1".into()),
        },
        &requests,
    )
    .unwrap()
    .corrections
}
fn config() -> PppConfig {
    PppConfig {
        initial_position: epochs(false)[0].approximate_position.unwrap(),
        initial_receiver_clock_s: 0.0,
        initial_zenith_delay_m: 2.3,
        expected_coordinate_system: "IGS20".into(),
        expected_reference_point: OrbitReferencePoint::AntennaPhaseCenter,
        transmit: Default::default(),
        arcs: Default::default(),
        filter: PppFilterSettings {
            position_model: PppPositionModel::ConstantVelocity(Default::default()),
            ..Default::default()
        },
        phase_windup: None,
        receiver_antenna: None,
        satellite_antenna: None,
        receiver_reference: ReceiverPositionReference::InstantaneousMarker,
        solid_earth_tide: None,
    }
}

#[test]
fn optional_native_head_prior_matches_six_actual_initializers_and_keeps_other_policies() {
    for line in include_str!("fixtures/upstream_ppp_initial.csv")
        .lines()
        .skip(1)
    {
        let fields: Vec<_> = line.split(',').collect();
        let v: Vec<f64> = fields[1..].iter().map(|s| s.parse().unwrap()).collect();
        let time = gnss_rust::GnssTime::new(v[0] as i32, v[1]).unwrap();
        let mut original = config();
        original.initial_position = EcefCoord::new(v[2], v[3], v[4]);
        let changed = original
            .clone()
            .with_native_precise_head_prior(time)
            .unwrap();
        assert_eq!(changed.initial_position, original.initial_position);
        assert_eq!(
            changed.initial_receiver_clock_s,
            original.initial_receiver_clock_s
        );
        assert!((changed.initial_zenith_delay_m - v[6]).abs() < 1e-10);
        let filter = PppFilter::new(changed.clone()).unwrap();
        let state = filter.state();
        let expected = [v[7], v[8], v[9], v[11], v[12], v[10], v[10], v[10]];
        for (i, row) in state.covariance.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                assert_eq!(value, if i == j { expected[i] } else { 0.0 });
            }
        }
        let mut restored = changed;
        restored.initial_zenith_delay_m = original.initial_zenith_delay_m;
        restored.filter.initial_head_variance_m2 = original.filter.initial_head_variance_m2;
        restored.filter.position_model = original.filter.position_model;
        assert_eq!(
            restored, original,
            "only declared initial head means/variances may change"
        );
        let mut static_config = original;
        static_config.filter.position_model = PppPositionModel::Static;
        let static_prior = static_config.with_native_precise_head_prior(time).unwrap();
        assert_eq!(static_prior.filter.position_model, PppPositionModel::Static);
        assert_eq!(PppFilter::new(static_prior).unwrap().state().mean.len(), 5);
    }
    let mut invalid = config();
    invalid.initial_position = EcefCoord::new(f64::NAN, 0.0, 0.0);
    assert!(
        invalid
            .with_native_precise_head_prior(epochs(false)[0].time)
            .is_err()
    );
}
fn solution(result: PppEpoch) -> Box<PppSolution> {
    match result {
        PppEpoch::Solution(s) => s,
        _ => panic!("expected solution"),
    }
}

fn inject_receiver_clock_metres(epoch: &mut gnss_rust::observation::ObservationEpoch, offset: f64) {
    for observation in &mut epoch.observations {
        if let Some(code) = &mut observation.pseudorange_m {
            *code += offset;
        }
        if let Some(phase) = &mut observation.carrier_phase_cycles {
            let hz = gnss_rust::dual_frequency::tracking_frequency_hz(
                observation.satellite,
                &observation.tracking_code,
                observation.glonass_channel,
            )
            .unwrap();
            *phase += offset * hz / gnss_rust::constants::SPEED_OF_LIGHT;
        }
    }
}

#[test]
fn independent_epoch_clock_preserves_clean_noisy_accuracy_and_tracks_injected_common_jumps() {
    let products = products();
    let selections = selections();
    let truth = truth();
    assert_eq!(
        PppFilterSettings::default().clock_model,
        PppClockModel::RandomWalk
    );
    for noisy in [false, true] {
        let mut cfg = config();
        cfg.filter.clock_model = PppClockModel::IndependentEpoch;
        let mut control = PppFilter::new(cfg.clone()).unwrap();
        let mut jumped = PppFilter::new(cfg).unwrap();
        let mut legacy = PppFilter::new(config()).unwrap();
        let mut legacy_rejected = 0;
        for (i, e) in epochs(noisy).iter().enumerate() {
            let offset = if i < 8 {
                0.0
            } else {
                [30000.0, -20000.0, 50000.0, 0.0][(i - 8) % 4]
            };
            let mut shifted = e.clone();
            inject_receiver_clock_metres(&mut shifted, offset);
            let a = solution(
                control
                    .process_epoch(e, &products, &biases(e.time), &selections)
                    .unwrap(),
            );
            let b = solution(
                jumped
                    .process_epoch(&shifted, &products, &biases(e.time), &selections)
                    .unwrap(),
            );
            assert_eq!(b.clock_model, PppClockModel::IndependentEpoch);
            assert_eq!(b.position.status, SolutionStatus::PppFloat);
            let aa = a.position.position_ecef.unwrap();
            let bb = b.position.position_ecef.unwrap();
            assert!(
                (aa.x - bb.x).hypot(aa.y - bb.y).hypot(aa.z - bb.z) < 0.001,
                "clock-injection equivalence, epoch {i}"
            );
            assert!(
                ((b.position.receiver_clock_bias - a.position.receiver_clock_bias)
                    * gnss_rust::constants::SPEED_OF_LIGHT
                    - offset)
                    .abs()
                    < 0.001
            );
            assert!(distance(bb, &truth[i][3..6]) < if noisy { 1.0 } else { 0.01 });
            if legacy
                .process_epoch(&shifted, &products, &biases(e.time), &selections)
                .is_err()
            {
                legacy_rejected += 1;
            }
        }
        assert!(
            legacy_rejected > 10,
            "original random-walk gate must retain clock-jump failures"
        );
    }
}

#[test]
fn independent_clock_rejection_and_event_preserve_committed_state_before_recovery() {
    let mut cfg = config();
    cfg.filter.clock_model = PppClockModel::IndependentEpoch;
    let mut filter = PppFilter::new(cfg).unwrap();
    let p = products();
    let selection = selections();
    let mut epochs = epochs(false);
    for e in epochs.iter().take(6) {
        filter
            .process_epoch(e, &p, &biases(e.time), &selection)
            .unwrap();
    }
    let state = filter.state().clone();
    let time = filter.state_time();
    *epochs[6].observations[0].pseudorange_m.as_mut().unwrap() += 1000.0;
    assert!(
        filter
            .process_epoch(&epochs[6], &p, &biases(epochs[6].time), &selection)
            .is_err()
    );
    assert_eq!(filter.state(), &state);
    assert_eq!(filter.state_time(), time);
    epochs[7].flag = 6;
    assert!(matches!(
        filter
            .process_epoch(&epochs[7], &p, &biases(epochs[7].time), &selection)
            .unwrap(),
        PppEpoch::Event(_)
    ));
    assert_eq!(filter.state(), &state);
    assert_eq!(filter.state_time(), time);
    inject_receiver_clock_metres(&mut epochs[8], 30000.0);
    let result = solution(
        filter
            .process_epoch(&epochs[8], &p, &biases(epochs[8].time), &selection)
            .unwrap(),
    );
    assert_eq!(result.clock_model, PppClockModel::IndependentEpoch);
    assert!(distance(result.position.position_ecef.unwrap(), &truth()[8][3..6]) < 0.01);
    assert!(
        result
            .arc_commit
            .accepted
            .values()
            .all(|a| a.accepted_epochs == 1)
    );
}
fn distance(a: EcefCoord, t: &[f64]) -> f64 {
    (a.x - t[0]).hypot(a.y - t[1]).hypot(a.z - t[2])
}

fn acceleration_epochs(noisy: bool) -> Vec<gnss_rust::observation::ObservationEpoch> {
    RinexObservationReader::new(Cursor::new(if noisy {
        include_str!("fixtures/synthetic_ppp_acceleration_noisy.obs")
    } else {
        include_str!("fixtures/synthetic_ppp_acceleration.obs")
    }))
    .unwrap()
    .collect::<Result<_, _>>()
    .unwrap()
}
fn acceleration_truth() -> Vec<Vec<f64>> {
    include_str!("fixtures/synthetic_ppp_acceleration_truth.csv")
        .lines()
        .skip(1)
        .map(|l| l.split(',').map(|v| v.parse().unwrap()).collect())
        .collect()
}
fn acceleration_config() -> PppConfig {
    let mut c = config();
    c.filter.position_model = PppPositionModel::ConstantAcceleration(Default::default());
    c
}

// Independent code-only seed: no truth position/velocity is fed to either solver.
fn code_motion_seed(e: &gnss_rust::observation::ObservationEpoch) -> PppMotionSeed {
    use gnss_rust::precise_spp::*;
    let corrections = biases(e.time);
    let table = PreciseCodeBiases::new(
        corrections
            .entries()
            .values()
            .flatten()
            .filter(|b| b.tracking_code == "1C")
            .map(|b| PreciseCodeBias {
                satellite: b.satellite,
                tracking_code: b.tracking_code.clone(),
                correction_m: b.code_add_m,
                clock_source: b.clock_source,
                clock_reference: b.clock_reference.clone(),
                valid_from: b.valid_from,
                valid_until: b.valid_until,
            })
            .collect(),
    )
    .unwrap();
    let result = solve_epoch_precise(
        e,
        &products(),
        &table,
        &PreciseSppConfig {
            spp: gnss_rust::spp::SppConfig {
                use_ionosphere: false,
                ..Default::default()
            },
            expected_coordinate_system: "IGS20".into(),
            expected_reference_point: OrbitReferencePoint::AntennaPhaseCenter,
            clock_reference: "synthetic-precise-clock-v1".into(),
            transmit: Default::default(),
            ionosphere: None,
        },
    )
    .unwrap();
    PppMotionSeed {
        solution: result.position,
        coordinate_system: result.coordinate_system,
        receiver_reference: ReceiverPositionReference::InstantaneousMarker,
        provenance: "synthetic-precise-1C-code-SPP-no-truth-no-Doppler".into(),
    }
}

#[test]
fn pva_code_seed_tracks_the_independent_car_without_velocity_or_truth_input() {
    let p = products();
    let s = selections();
    for noisy in [false, true] {
        let mut c = acceleration_config();
        c.filter.motion_seed_policy = Some(PppMotionSeedPolicy {
            update_velocity: false,
            ..Default::default()
        });
        let mut f = PppFilter::new(c).unwrap();
        for (i, e) in acceleration_epochs(noisy).iter().enumerate() {
            let seed = code_motion_seed(e);
            assert_eq!(seed.solution.velocity_ecef, None);
            let r = solution(
                f.process_epoch_with_seed(e, &p, &biases(e.time), &s, Some(&seed))
                    .unwrap(),
            );
            assert!(
                distance(
                    r.position.position_ecef.unwrap(),
                    &acceleration_truth()[i][3..6]
                ) < if noisy { 1.0 } else { 0.01 },
                "car seed epoch {i}"
            );
            let d = r.motion_seed.unwrap();
            assert_eq!(d.initialized, i == 0);
            assert!(d.committed);
            assert_eq!(d.velocity_update, PppSeedVelocityUpdate::Disabled);
        }
    }
}

#[test]
fn seed_contract_rejects_time_frame_reference_status_and_missing_source_transactionally() {
    let e = &acceleration_epochs(false)[0];
    let p = products();
    let s = selections();
    let seed = code_motion_seed(e);
    let mut c = acceleration_config();
    c.filter.motion_seed_policy = Some(Default::default());
    for defect in 0..7 {
        let mut f = PppFilter::new(c.clone()).unwrap();
        let old = f.state().clone();
        let mut bad = seed.clone();
        let reason = match defect {
            0 => {
                bad.solution.time = e.time.checked_add_seconds(0.2).unwrap();
                PppMotionSeedError::TimeMismatch
            }
            1 => {
                bad.coordinate_system = "ITRF?".into();
                PppMotionSeedError::FrameMismatch
            }
            2 => {
                bad.receiver_reference = ReceiverPositionReference::NominalMarker;
                PppMotionSeedError::ReceiverReferenceMismatch
            }
            3 => {
                bad.solution.status = SolutionStatus::PppFloat;
                PppMotionSeedError::InvalidSolution
            }
            4 => {
                bad.provenance = " ".into();
                PppMotionSeedError::MissingProvenance
            }
            5 => {
                bad.solution.position_ecef = Some(EcefCoord::new(f64::NAN, 0.0, 0.0));
                PppMotionSeedError::InvalidSolution
            }
            _ => {
                bad.solution.position_ecef = Some(EcefCoord::new(0.0, 0.0, 0.0));
                PppMotionSeedError::InvalidPosition
            }
        };
        assert_eq!(
            f.process_epoch_with_seed(e, &p, &biases(e.time), &s, Some(&bad)),
            Err(PppError::MotionSeed(reason))
        );
        assert_eq!(f.state(), &old);
        assert_eq!(f.state_time(), None);
        assert_eq!(f.last_motion_seed(), None);
    }
    let mut f = PppFilter::new(c).unwrap();
    assert_eq!(
        f.process_epoch(e, &p, &biases(e.time), &s),
        Err(PppError::MissingMotionSeed)
    );
    let mut event = e.clone();
    event.time = e.time.checked_add_seconds(0.1).unwrap();
    event.flag = 6;
    assert!(matches!(
        f.process_epoch(&event, &p, &biases(e.time), &s),
        Ok(PppEpoch::Event(_))
    ));
    let mut f = PppFilter::new(acceleration_config()).unwrap();
    assert_eq!(
        f.process_epoch_with_seed(e, &p, &biases(e.time), &s, Some(&seed)),
        Err(PppError::MotionSeed(PppMotionSeedError::Disabled))
    );
}

#[test]
fn first_velocity_seed_is_fused_once_and_never_imports_seed_clock() {
    let e = &acceleration_epochs(false)[0];
    let mut seed = code_motion_seed(e);
    // Injected unit input for the covariance calculation, not a synthetic Doppler claim.
    seed.provenance.push_str(";injected-unit-velocity");
    seed.solution.velocity_ecef = Some([1.0, -2.0, 3.0]);
    seed.solution.velocity_covariance = Some([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]]);
    let mut c = acceleration_config();
    c.filter.motion_seed_policy = Some(Default::default());
    let mut left = PppFilter::new(c.clone()).unwrap();
    let mut right = PppFilter::new(c).unwrap();
    seed.solution.receiver_clock_bias = 1e100;
    let r = solution(
        left.process_epoch_with_seed(e, &products(), &biases(e.time), &selections(), Some(&seed))
            .unwrap(),
    );
    seed.solution.receiver_clock_bias = -1e100;
    right
        .process_epoch_with_seed(e, &products(), &biases(e.time), &selections(), Some(&seed))
        .unwrap();
    assert_eq!(left.state(), right.state());
    assert_eq!(r.position.velocity_ecef, seed.solution.velocity_ecef);
    for i in 5..8 {
        assert!((left.state().covariance[i][i] - 100.0 * 4.0 / 104.0).abs() < 1e-12);
    }
    assert!(matches!(
        r.motion_seed.unwrap().velocity_update,
        PppSeedVelocityUpdate::Applied { .. }
    ));
}

#[test]
fn rejected_seeded_gnss_update_does_not_commit_any_motion_update() {
    let ee = acceleration_epochs(false);
    let p = products();
    let s = selections();
    let mut c = acceleration_config();
    c.filter.motion_seed_policy = Some(Default::default());
    let mut f = PppFilter::new(c).unwrap();
    let seed = code_motion_seed(&ee[0]);
    f.process_epoch_with_seed(&ee[0], &p, &biases(ee[0].time), &s, Some(&seed))
        .unwrap();
    let old = f.state().clone();
    let old_time = f.state_time();
    let seed = code_motion_seed(&ee[1]);
    let mut bad = ee[1].clone();
    inject_receiver_clock_metres(&mut bad, 1e7);
    assert!(
        f.process_epoch_with_seed(&bad, &p, &biases(bad.time), &s, Some(&seed))
            .is_err()
    );
    assert_eq!(f.state(), &old);
    assert_eq!(f.state_time(), old_time);
    assert!(!f.last_motion_seed().unwrap().committed);
    let seed = code_motion_seed(&ee[2]);
    assert!(
        f.process_epoch_with_seed(&ee[2], &p, &biases(ee[2].time), &s, Some(&seed))
            .is_ok()
    );
}

#[test]
fn velocity_gate_and_missing_or_bad_covariance_are_explicit_without_zero_fallback() {
    let e = &acceleration_epochs(false)[0];
    let mut c = acceleration_config();
    c.filter.motion_seed_policy = Some(PppMotionSeedPolicy {
        initialize_from_seed: false,
        velocity_nis_gate: Some(16.27),
        recovery: gnss_rust::ppp_motion::PvaRecoveryConfig {
            position_variance_m2: 1e7,
            position_disagreement_gate: None,
            ..Default::default()
        },
        ..Default::default()
    });
    let mut baseline = PppFilter::new(acceleration_config()).unwrap();
    baseline
        .process_epoch(e, &products(), &biases(e.time), &selections())
        .unwrap();
    for defect in 0..4 {
        let mut seed = code_motion_seed(e);
        let expected = match defect {
            0 => PppSeedVelocityUpdate::Unavailable,
            1 => {
                seed.solution.velocity_ecef = Some([1.0; 3]);
                PppSeedVelocityUpdate::MissingCovariance
            }
            2 => {
                seed.solution.velocity_ecef = Some([1.0; 3]);
                seed.solution.velocity_covariance = Some([[0.0; 3]; 3]);
                PppSeedVelocityUpdate::Rejected(gnss_rust::filter::FilterError::InvalidCovariance)
            }
            _ => {
                seed.solution.velocity_ecef = Some([10000.0; 3]);
                seed.solution.velocity_covariance =
                    Some([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]]);
                PppSeedVelocityUpdate::Rejected(gnss_rust::filter::FilterError::InnovationRejected)
            }
        };
        let mut f = PppFilter::new(c.clone()).unwrap();
        let r = solution(
            f.process_epoch_with_seed(e, &products(), &biases(e.time), &selections(), Some(&seed))
                .unwrap(),
        );
        assert_eq!(r.motion_seed.unwrap().velocity_update, expected);
        assert_eq!(f.state(), baseline.state());
    }
}

#[test]
fn optional_pva_head_prior_matches_six_native_initializations_and_preserves_noise() {
    for line in include_str!("fixtures/upstream_ppp_acceleration_initial.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = line.split(',').collect();
        let d: Vec<f64> = f[1..].iter().map(|s| s.parse().unwrap()).collect();
        assert_eq!(d.len(), 18);
        let time = gnss_rust::GnssTime::new(d[0] as i32, d[1]).unwrap();
        let mut c = acceleration_config();
        c.initial_position = EcefCoord::new(d[2], d[3], d[4]);
        let original = c.clone();
        let changed = c.with_native_precise_head_prior(time).unwrap();
        let f = PppFilter::new(changed.clone()).unwrap();
        assert!((f.state().mean[4] - d[6]).abs() < 1e-10);
        assert_eq!(f.state().mean[3], d[5]);
        assert_eq!(&f.state().mean[5..], &d[12..]);
        let expected = [
            d[7], d[7], d[7], d[10], d[11], d[8], d[8], d[8], d[9], d[9], d[9],
        ];
        for (i, row) in f.state().covariance.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                assert_eq!(*v, if i == j { expected[i] } else { 0.0 });
            }
        }
        let mut restored = changed;
        restored.initial_zenith_delay_m = original.initial_zenith_delay_m;
        restored.filter.initial_head_variance_m2 = original.filter.initial_head_variance_m2;
        assert_eq!(restored, original);
    }
}

#[test]
fn independent_five_hz_car_trace_tracks_acceleration_under_predeclared_accuracy_bounds() {
    let p = products();
    let s = selections();
    let truth = acceleration_truth();
    for noisy in [false, true] {
        let mut f = PppFilter::new(acceleration_config()).unwrap();
        for (i, e) in acceleration_epochs(noisy).iter().enumerate() {
            let r = solution(f.process_epoch(e, &p, &biases(e.time), &s).unwrap());
            assert!(
                distance(r.position.position_ecef.unwrap(), &truth[i][3..6])
                    < if noisy { 1.0 } else { 0.01 },
                "car position epoch {i}"
            );
            assert_eq!(
                r.position_model,
                acceleration_config().filter.position_model
            );
            assert!(
                r.acceleration_ecef_mps2
                    .unwrap()
                    .iter()
                    .all(|v| v.is_finite())
            );
            assert_eq!(f.state().mean.len(), 11 + r.ambiguities.len());
            if i >= 8 {
                let v = r.position.velocity_ecef.unwrap();
                let delta = (v[0] - truth[i][6])
                    .hypot(v[1] - truth[i][7])
                    .hypot(v[2] - truth[i][8]);
                assert!(delta < 1.0, "car velocity epoch {i}: {delta}");
            }
            if i == 0 {
                assert_eq!(r.position.velocity_ecef, Some([0.0; 3]));
                assert_eq!(r.acceleration_ecef_mps2, Some([0.0; 3]));
                for j in 0..3 {
                    assert_eq!(r.acceleration_covariance_m2_s4.unwrap()[j][j], 100.0);
                }
            }
        }
    }
}

#[test]
fn pva_rejection_keeps_full_numeric_state_and_recovery_predicts_over_the_full_gap() {
    let p = products();
    let s = selections();
    let mut ee = acceleration_epochs(false);
    let mut c = acceleration_config();
    c.filter.clock_model = PppClockModel::IndependentEpoch;
    c.filter.rejected_update_policy = PppRejectedUpdatePolicy::RetainContinuousArcs;
    let mut f = PppFilter::new(c.clone()).unwrap();
    let mut skipped = PppFilter::new(c).unwrap();
    for e in &ee[..12] {
        f.process_epoch(e, &p, &biases(e.time), &s).unwrap();
        skipped.process_epoch(e, &p, &biases(e.time), &s).unwrap();
    }
    let old = f.state().clone();
    let time = f.state_time();
    inject_receiver_clock_metres(&mut ee[12], 1e7);
    // Even an independent finite clock prior rejects an enormous common jump.
    assert_eq!(
        f.process_epoch(&ee[12], &p, &biases(ee[12].time), &s),
        Err(PppError::Filter(
            gnss_rust::filter::FilterError::InnovationRejected
        ))
    );
    assert_eq!(f.state(), &old);
    assert_eq!(f.state_time(), time);
    let r = solution(
        f.process_epoch(&ee[13], &p, &biases(ee[13].time), &s)
            .unwrap(),
    );
    skipped
        .process_epoch(&ee[13], &p, &biases(ee[13].time), &s)
        .unwrap();
    assert_eq!(f.state(), skipped.state());
    assert!(
        distance(
            r.position.position_ecef.unwrap(),
            &acceleration_truth()[13][3..6]
        ) < 0.01
    );
    assert_eq!(
        r.acceleration_covariance_m2_s4.unwrap(),
        std::array::from_fn(|i| std::array::from_fn(|j| f.state().covariance[8 + i][8 + j]))
    );
}

#[test]
fn pva_phase_outage_remaps_ambiguities_after_eleven_head_states_and_events_reset_arcs() {
    let p = products();
    let s = selections();
    let mut ee = acceleration_epochs(false);
    let mut f = PppFilter::new(acceleration_config()).unwrap();
    for e in &ee[..8] {
        f.process_epoch(e, &p, &biases(e.time), &s).unwrap();
    }
    let sat = ee[8].observations[0].satellite;
    let id = f.phase_arc(sat).unwrap().id;
    for o in &mut ee[8].observations {
        if o.satellite == sat {
            o.carrier_phase_cycles = None;
        }
    }
    let r = solution(
        f.process_epoch(&ee[8], &p, &biases(ee[8].time), &s)
            .unwrap(),
    );
    assert!(r.ambiguities.iter().all(|a| a.satellite != sat));
    assert_eq!(f.state().mean.len(), 11 + r.ambiguities.len());
    let r = solution(
        f.process_epoch(&ee[9], &p, &biases(ee[9].time), &s)
            .unwrap(),
    );
    assert_ne!(r.arc_commit.accepted[&sat].id, id);
    assert_eq!(r.arc_commit.accepted[&sat].accepted_epochs, 1);
    ee[10].flag = 6;
    let state = f.state().clone();
    assert!(matches!(
        f.process_epoch(&ee[10], &p, &biases(ee[10].time), &s)
            .unwrap(),
        PppEpoch::Event(_)
    ));
    assert_eq!(f.state(), &state);
    let r = solution(
        f.process_epoch(&ee[11], &p, &biases(ee[11].time), &s)
            .unwrap(),
    );
    assert!(
        r.arc_commit
            .accepted
            .values()
            .all(|a| a.accepted_epochs == 1)
    );
}

#[test]
fn pva_requires_finite_initial_states_positive_variances_and_finite_squared_noise() {
    assert_eq!(
        PppFilterSettings::default().position_model,
        PppPositionModel::Static
    );
    for case in 0..8 {
        let mut pva = PppAccelerationConfig::default();
        match case {
            0 => pva.initial_velocity_mps[0] = f64::NAN,
            1 => pva.initial_acceleration_mps2[1] = f64::INFINITY,
            2 => pva.initial_velocity_variance_m2_s2 = 0.0,
            3 => pva.initial_acceleration_variance_m2_s4 = -1.0,
            4 => pva.horizontal_sigma_mps2_per_sqrt_s = -0.1,
            5 => pva.vertical_sigma_mps2_per_sqrt_s = f64::NAN,
            6 => pva.horizontal_sigma_mps2_per_sqrt_s = 1e308,
            7 => pva.initial_acceleration_variance_m2_s4 = f64::INFINITY,
            _ => unreachable!(),
        }
        let mut c = acceleration_config();
        c.filter.position_model = PppPositionModel::ConstantAcceleration(pva);
        assert!(matches!(
            PppFilter::new(c),
            Err(PppError::InvalidConfiguration)
        ));
    }
}

#[test]
fn retained_moving_arcs_keep_exact_state_and_only_successful_lock_counts() {
    let p = products();
    let selection = selections();
    let mut cfg = config();
    cfg.filter.rejected_update_policy = PppRejectedUpdatePolicy::RetainContinuousArcs;
    let mut f = PppFilter::new(cfg).unwrap();
    let mut ee = epochs(false);
    for e in &ee[..6] {
        f.process_epoch(e, &p, &biases(e.time), &selection).unwrap();
    }
    let old: BTreeMap<_, _> = selection
        .keys()
        .filter_map(|s| f.phase_arc(*s).map(|a| (*s, a.clone())))
        .collect();
    let state = f.state().clone();
    let time = f.state_time();
    let broken = ee[6].observations[0].satellite;
    *ee[6].observations[0].pseudorange_m.as_mut().unwrap() += 1000.0;
    assert!(
        f.process_epoch(&ee[6], &p, &biases(ee[6].time), &selection)
            .is_err()
    );
    assert_eq!(f.state(), &state);
    assert_eq!(f.state_time(), time);
    assert!(!f.last_retained_phase_arcs().is_empty());
    assert!(f.phase_arc(broken).is_none());
    for s in f.last_retained_phase_arcs() {
        assert_eq!(f.phase_arc(*s), old.get(s));
    }
    let retained = f.last_retained_phase_arcs().clone();
    let r = solution(
        f.process_epoch(&ee[7], &p, &biases(ee[7].time), &selection)
            .unwrap(),
    );
    assert_eq!(
        r.rejected_update_policy,
        PppRejectedUpdatePolicy::RetainContinuousArcs
    );
    assert!(f.last_retained_phase_arcs().is_empty());
    for s in retained {
        let a = &r.arc_commit.accepted[&s];
        assert_eq!(a.id, old[&s].id);
        assert_eq!(a.accepted_epochs, old[&s].accepted_epochs + 1);
    }
    assert!(distance(r.position.position_ecef.unwrap(), &truth()[7][3..6]) < 0.002);
}

#[test]
fn retention_cannot_accept_failed_initialization_or_survive_events_and_invalid_epochs() {
    let p = products();
    let s = selections();
    let mut cfg = config();
    cfg.filter.rejected_update_policy = PppRejectedUpdatePolicy::RetainContinuousArcs;
    let mut f = PppFilter::new(cfg).unwrap();
    let mut ee = epochs(false);
    *ee[0].observations[0].pseudorange_m.as_mut().unwrap() += 1000.0;
    assert!(
        f.process_epoch(&ee[0], &p, &biases(ee[0].time), &s)
            .is_err()
    );
    assert!(f.state_time().is_none());
    assert!(f.last_retained_phase_arcs().is_empty());
    for e in &ee[1..5] {
        f.process_epoch(e, &p, &biases(e.time), &s).unwrap();
    }
    let old = f.state().clone();
    ee[5].flag = 6;
    assert!(matches!(
        f.process_epoch(&ee[5], &p, &biases(ee[5].time), &s)
            .unwrap(),
        PppEpoch::Event(_)
    ));
    assert_eq!(f.state(), &old);
    assert!(s.keys().all(|s| f.phase_arc(*s).is_none()));
    let recovered = solution(
        f.process_epoch(&ee[6], &p, &biases(ee[6].time), &s)
            .unwrap(),
    );
    assert!(
        recovered
            .arc_commit
            .accepted
            .values()
            .all(|a| a.accepted_epochs == 1)
    );
    let old = f.state().clone();
    ee[7].flag = 9;
    assert!(
        f.process_epoch(&ee[7], &p, &biases(ee[7].time), &s)
            .is_err()
    );
    assert_eq!(f.state(), &old);
    assert!(f.last_retained_phase_arcs().is_empty());
    assert!(s.keys().all(|s| f.phase_arc(*s).is_none()));
}

#[test]
fn rounded_clean_trace_retains_full_trace_two_mm_failures_and_meets_the_bound_afterward() {
    let products = products();
    let mut f = PppFilter::new(config()).unwrap();
    let truth = truth();
    let s = selections();
    let mut previous = BTreeMap::new();
    let mut over_two_mm = Vec::new();
    for (i, e) in epochs(false).iter().enumerate() {
        let result = solution(f.process_epoch(e, &products, &biases(e.time), &s).unwrap());
        assert_eq!(result.position.status, SolutionStatus::PppFloat);
        assert_eq!(result.position_model, config().filter.position_model);
        let error = distance(result.position.position_ecef.unwrap(), &truth[i][3..6]);
        if error >= 0.002 {
            over_two_mm.push(i);
        }
        // The original full-trace 2 mm criterion is NOT met on rounded RINEX.
        // Keep all failures visible, and preserve the bound after convergence.
        if i >= 4 {
            assert!(error < 0.002, "epoch {i}: {error}");
            assert!(result.code_rms_m < 0.0001);
        }
        assert!(result.phase_rms_m < 0.0001);
        let velocity = result.position.velocity_ecef.unwrap();
        if i == 0 {
            assert_eq!(velocity, [0.0; 3]);
            for j in 0..3 {
                assert_eq!(result.position.velocity_covariance.unwrap()[j][j], 100.0);
            }
        }
        if i >= 2 {
            for j in 0..3 {
                assert!((velocity[j] - truth[i][6 + j]).abs() < 0.02);
            }
        }
        assert_eq!(f.state().mean.len(), 8 + result.ambiguities.len());
        let current: BTreeMap<_, _> = result
            .arc_commit
            .accepted
            .iter()
            .map(|(&sat, a)| {
                if let Some((id, count)) = previous.get(&sat) {
                    if *id == a.id {
                        assert_eq!(a.accepted_epochs, count + 1);
                    }
                }
                (sat, (a.id, a.accepted_epochs))
            })
            .collect();
        previous = current;
    }
    assert_eq!(over_two_mm, vec![1, 2, 3]);
}

#[test]
fn independently_generated_unrounded_observables_isolate_carrier_quantization() {
    let p = products();
    let t = truth();
    let selections = selections();
    for (phase, code) in [(true, false), (false, true), (true, true)] {
        let mut epochs = epochs(false);
        for line in include_str!("fixtures/upstream_ppp_kinematic_unrounded.csv")
            .lines()
            .skip(1)
        {
            let v: Vec<_> = line.split(',').collect();
            let i: usize = v[0].parse().unwrap();
            let sat: SatelliteId = v[1].parse().unwrap();
            for (j, tracking) in ["1C", "2W"].iter().enumerate() {
                let o = epochs[i]
                    .observations
                    .iter_mut()
                    .find(|o| o.satellite == sat && o.tracking_code == *tracking)
                    .unwrap();
                if code {
                    o.pseudorange_m = Some(v[2 + 2 * j].parse().unwrap());
                }
                if phase {
                    o.carrier_phase_cycles = Some(v[3 + 2 * j].parse().unwrap());
                }
            }
        }
        let mut f = PppFilter::new(config()).unwrap();
        let mut maximum = 0.0_f64;
        for (i, e) in epochs.iter().enumerate() {
            let result = solution(
                f.process_epoch(e, &p, &biases(e.time), &selections)
                    .unwrap(),
            );
            maximum = maximum.max(distance(
                result.position.position_ecef.unwrap(),
                &t[i][3..6],
            ));
        }
        if phase {
            assert!(maximum < 0.002, "phase={phase} code={code}: {maximum}");
        } else {
            assert!(
                maximum >= 0.002,
                "code precision alone does not remove the observed carrier-quantization effect"
            );
        }
    }
}

#[test]
fn independent_noisy_moving_trace_meets_predeclared_position_and_velocity_bounds() {
    let p = products();
    let truth = truth();
    let mut f = PppFilter::new(config()).unwrap();
    let s = selections();
    let mut last = None;
    for (i, e) in epochs(true).iter().enumerate() {
        let result = solution(f.process_epoch(e, &p, &biases(e.time), &s).unwrap());
        assert!(
            distance(result.position.position_ecef.unwrap(), &truth[i][3..6]) < 1.0,
            "epoch {i}"
        );
        if i >= 2 {
            for j in 0..3 {
                assert!((result.position.velocity_ecef.unwrap()[j] - truth[i][6 + j]).abs() < 0.02);
            }
        }
        assert!(result.code_rms_m < 0.2);
        assert!(result.phase_rms_m < 0.005);
        last = Some(result);
    }
    assert!(
        distance(
            last.unwrap().position.position_ecef.unwrap(),
            &truth[63][3..6]
        ) < 0.1
    );
}

#[test]
fn default_static_path_exposes_motion_failure_without_committing_predicted_state() {
    let mut c = config();
    c.filter.position_model = PppPositionModel::Static;
    let mut f = PppFilter::new(c).unwrap();
    let p = products();
    let e = epochs(false);
    let s = selections();
    let first = solution(f.process_epoch(&e[0], &p, &biases(e[0].time), &s).unwrap());
    assert!(first.position.velocity_ecef.is_none() && first.position.velocity_covariance.is_none());
    assert_eq!(first.position_model, PppPositionModel::Static);
    let before = f.state().clone();
    let time = f.state_time();
    assert!(f.process_epoch(&e[1], &p, &biases(e[1].time), &s).is_err());
    assert_eq!(f.state(), &before);
    assert_eq!(f.state_time(), time);
}

#[test]
fn rejected_moving_update_preserves_velocity_covariance_and_recovers_over_full_gap() {
    let mut f = PppFilter::new(config()).unwrap();
    let p = products();
    let mut e = epochs(false);
    let s = selections();
    f.process_epoch(&e[0], &p, &biases(e[0].time), &s).unwrap();
    let before = f.state().clone();
    let time = f.state_time();
    *e[1].observations[0].pseudorange_m.as_mut().unwrap() += 1000.0;
    assert!(f.process_epoch(&e[1], &p, &biases(e[1].time), &s).is_err());
    assert_eq!(f.state(), &before);
    assert_eq!(f.state_time(), time);
    let recovered = solution(f.process_epoch(&e[2], &p, &biases(e[2].time), &s).unwrap());
    assert!(
        distance(recovered.position.position_ecef.unwrap(), &truth()[2][3..6]) < 0.002,
        "error {}",
        distance(recovered.position.position_ecef.unwrap(), &truth()[2][3..6])
    );
    assert!(
        recovered
            .arc_commit
            .accepted
            .values()
            .all(|a| a.accepted_epochs == 1)
    );
}

#[test]
fn phase_outage_and_reacquisition_remap_ambiguities_after_the_velocity_head() {
    let mut f = PppFilter::new(config()).unwrap();
    let p = products();
    let mut e = epochs(false);
    let s = selections();
    let mut first = None;
    for epoch in e.iter().take(10) {
        first = Some(solution(
            f.process_epoch(epoch, &p, &biases(epoch.time), &s).unwrap(),
        ));
    }
    let first = first.unwrap();
    let sat: SatelliteId = "G01".parse().unwrap();
    let other: SatelliteId = "G02".parse().unwrap();
    for o in &mut e[10].observations {
        if o.satellite == sat {
            o.carrier_phase_cycles = None;
        }
    }
    let outage = solution(
        f.process_epoch(&e[10], &p, &biases(e[10].time), &s)
            .unwrap(),
    );
    assert_eq!(f.state().mean.len(), 8 + outage.ambiguities.len());
    assert!(!outage.arc_commit.accepted.contains_key(&sat));
    let recovered = solution(
        f.process_epoch(&e[11], &p, &biases(e[11].time), &s)
            .unwrap(),
    );
    assert_ne!(
        first.arc_commit.accepted[&sat].id,
        recovered.arc_commit.accepted[&sat].id
    );
    assert_eq!(
        first.arc_commit.accepted[&other].id,
        recovered.arc_commit.accepted[&other].id
    );
    assert!(
        distance(
            recovered.position.position_ecef.unwrap(),
            &truth()[11][3..6]
        ) < 0.002
    );
    assert_eq!(f.state().mean.len(), 8 + recovered.ambiguities.len());
}

#[test]
fn kinematic_configuration_requires_finite_velocity_and_valid_noise_without_changing_default() {
    assert_eq!(
        PppFilterSettings::default().position_model,
        PppPositionModel::Static
    );
    let original = PppKinematicConfig::default();
    for k in [
        PppKinematicConfig {
            initial_velocity_mps: [f64::NAN, 0.0, 0.0],
            ..original
        },
        PppKinematicConfig {
            initial_velocity_variance_m2_s2: 0.0,
            ..original
        },
        PppKinematicConfig {
            initial_velocity_variance_m2_s2: f64::INFINITY,
            ..original
        },
        PppKinematicConfig {
            velocity_process_noise_m2_s3: -1.0,
            ..original
        },
        PppKinematicConfig {
            velocity_process_noise_m2_s3: f64::NAN,
            ..original
        },
    ] {
        let mut c = config();
        c.filter.position_model = PppPositionModel::ConstantVelocity(k);
        assert!(matches!(
            PppFilter::new(c),
            Err(PppError::InvalidConfiguration)
        ));
    }
}

fn code_only_biases(t: gnss_rust::GnssTime) -> gnss_rust::dual_frequency::SignalCorrections {
    let b = biases(t);
    let entries = b
        .entries()
        .values()
        .flatten()
        .cloned()
        .map(|mut c| {
            c.phase_add_m = None;
            c.phase_reference = None;
            c
        })
        .collect();
    gnss_rust::dual_frequency::SignalCorrections::new(entries).unwrap()
}
fn raw_selections() -> BTreeMap<SatelliteId, DualFrequencyConfig> {
    selections()
        .into_iter()
        .map(|(s, mut c)| {
            c.phase_reference = "raw-carrier-FLOAT-gauge-v1".into();
            (s, c)
        })
        .collect()
}
fn dsb_biases(t: gnss_rust::GnssTime) -> gnss_rust::dual_frequency::SignalCorrections {
    let p = read_bias_sinex(
        Cursor::new(include_str!("fixtures/synthetic_ppp_raw_phase_dsb.bsx")),
        Default::default(),
    )
    .unwrap();
    let m = gnss_rust::dsb_clock::GpsIfDsbClock::new(&p).unwrap();
    let requests: Vec<_> = selections()
        .into_iter()
        .flat_map(|(sat, c)| {
            c.tracking_codes
                .into_iter()
                .map(move |code| BiasSignalRequest {
                    satellite: sat,
                    tracking_code: code,
                    glonass_channel: None,
                })
        })
        .collect();
    m.corrections_at(
        t,
        &gnss_rust::dsb_clock::GpsIfClockBinding {
            product_id: "synthetic_ppp_raw_phase_dsb.bsx".into(),
            clock_source: ClockSource::RinexClk,
            clock_reference: "synthetic-precise-clock-v1".into(),
        },
        &requests,
    )
    .unwrap()
    .corrections
}

#[test]
fn explicit_raw_float_keeps_absent_phase_osbs_and_default_still_requires_them() {
    use gnss_rust::dual_frequency::*;
    let e = &epochs(false)[0];
    let before = e.clone();
    let sat = "G01".parse().unwrap();
    let selection = raw_selections();
    let code = code_only_biases(e.time);
    let strict = form_ionosphere_free(e, sat, &code, &selection[&sat]).unwrap();
    assert!(strict.phase.is_none());
    assert!(matches!(
        strict.phase_unavailability,
        Some(PhaseUnavailable::MissingCorrection { .. })
    ));
    let raw = form_ionosphere_free_with_phase_bias_mode(
        e,
        sat,
        &code,
        &selection[&sat],
        PhaseBiasMode::UncalibratedFloat,
    )
    .unwrap();
    assert_eq!(raw.phase_bias_mode, PhaseBiasMode::UncalibratedFloat);
    assert_eq!(raw.phase_reference, selection[&sat].phase_reference);
    assert!(
        raw.corrections
            .iter()
            .all(|c| c.phase_add_m.is_none() && c.phase_reference.is_none())
    );
    for i in 0..2 {
        assert_eq!(
            raw.phase.as_ref().unwrap().corrected_phase_m[i],
            raw.raw[i].carrier_phase_cycles.unwrap() * raw.wavelength_m[i]
        );
    }
    assert_eq!(*e, before);
    let mut cfg = config();
    let mut filter = PppFilter::new(cfg.clone()).unwrap();
    let state = filter.state().clone();
    assert!(matches!(
        filter.process_epoch(e, &products(), &code, &selection),
        Err(PppError::InsufficientPhase { available: 0 })
    ));
    assert_eq!(filter.state(), &state);
    cfg.filter.phase_bias_mode = PhaseBiasMode::UncalibratedFloat;
    let result = solution(
        PppFilter::new(cfg)
            .unwrap()
            .process_epoch(e, &products(), &code, &selection)
            .unwrap(),
    );
    assert_eq!(result.phase_bias_mode, PhaseBiasMode::UncalibratedFloat);
    assert!(!result.position.is_fixed());
}

#[test]
fn raw_float_rejects_supplied_phase_corrections_and_preserves_missing_or_half_cycle_reasons() {
    use gnss_rust::dual_frequency::*;
    let mut e = epochs(false)[0].clone();
    let sat = "G01".parse().unwrap();
    let selections = raw_selections();
    assert!(matches!(
        form_ionosphere_free_with_phase_bias_mode(
            &e,
            sat,
            &biases(e.time),
            &selections[&sat],
            PhaseBiasMode::UncalibratedFloat
        ),
        Err(DualFrequencyError::PhaseBiasModeMismatch { .. })
    ));
    let code = code_only_biases(e.time);
    let supplied_zero = SignalCorrections::new(
        code.entries()
            .values()
            .flatten()
            .cloned()
            .map(|mut c| {
                c.phase_add_m = Some(0.0);
                c.phase_reference = Some("declared-zero-OSB".into());
                c
            })
            .collect(),
    )
    .unwrap();
    assert!(matches!(
        form_ionosphere_free_with_phase_bias_mode(
            &e,
            sat,
            &supplied_zero,
            &selections[&sat],
            PhaseBiasMode::UncalibratedFloat
        ),
        Err(DualFrequencyError::PhaseBiasModeMismatch { .. })
    ));
    let mut wrong_clock = selections[&sat].clone();
    wrong_clock.clock_reference = "different-clock".into();
    assert!(matches!(
        form_ionosphere_free_with_phase_bias_mode(
            &e,
            sat,
            &code,
            &wrong_clock,
            PhaseBiasMode::UncalibratedFloat
        ),
        Err(DualFrequencyError::ClockDatumMismatch { .. })
    ));
    let mut offset_epoch = e.clone();
    offset_epoch.receiver_clock_offset_s = Some(1e-6);
    assert!(matches!(
        form_ionosphere_free_with_phase_bias_mode(
            &offset_epoch,
            sat,
            &code,
            &selections[&sat],
            PhaseBiasMode::UncalibratedFloat
        ),
        Err(DualFrequencyError::InvalidEpoch)
    ));
    e.observations[0].carrier_phase_cycles = None;
    let missing = form_ionosphere_free_with_phase_bias_mode(
        &e,
        sat,
        &code,
        &selections[&sat],
        PhaseBiasMode::UncalibratedFloat,
    )
    .unwrap();
    assert!(matches!(
        missing.phase_unavailability,
        Some(PhaseUnavailable::MissingPhase { .. })
    ));
    e = epochs(false)[0].clone();
    e.observations[0].lli = Some(2);
    let half = form_ionosphere_free_with_phase_bias_mode(
        &e,
        sat,
        &code,
        &selections[&sat],
        PhaseBiasMode::UncalibratedFloat,
    )
    .unwrap();
    assert!(matches!(
        half.phase_unavailability,
        Some(PhaseUnavailable::HalfCycle { .. })
    ));
}

#[test]
fn changing_raw_phase_model_or_reference_starts_a_new_arc_without_fake_zero_osbs() {
    use gnss_rust::dual_frequency::PhaseBiasMode;
    use gnss_rust::phase_arc::{ArcResetReason, PhaseArcTracker};
    use std::collections::BTreeSet;
    let e = epochs(false);
    let mut tracker = PhaseArcTracker::new(Default::default()).unwrap();
    let mut s = selections();
    for c in s.values_mut() {
        c.phase_reference = "same-name".into();
    }
    let mut calibrated = biases(e[0].time)
        .entries()
        .values()
        .flatten()
        .cloned()
        .collect::<Vec<_>>();
    for c in &mut calibrated {
        c.phase_reference = Some("same-name".into());
    }
    let calibrated = gnss_rust::dual_frequency::SignalCorrections::new(calibrated).unwrap();
    tracker.prepare_epoch(&e[0], &calibrated, &s).unwrap();
    let used: BTreeSet<_> = s.keys().copied().collect();
    let first = tracker.finish_epoch(&used).unwrap();
    let prepared = tracker
        .prepare_epoch_with_phase_bias_mode(
            &e[1],
            &code_only_biases(e[1].time),
            &s,
            PhaseBiasMode::UncalibratedFloat,
        )
        .unwrap();
    assert!(prepared.decisions().values().all(|d| {
        d.reset_reasons
            .contains(&ArcResetReason::PhaseBiasModelChanged)
    }));
    let second = tracker.finish_epoch(&used).unwrap();
    for sat in &used {
        assert_ne!(first.accepted[sat].id, second.accepted[sat].id);
    }
    for c in s.values_mut() {
        c.phase_reference = "new-raw-reference".into();
    }
    let prepared = tracker
        .prepare_epoch_with_phase_bias_mode(
            &e[2],
            &code_only_biases(e[2].time),
            &s,
            PhaseBiasMode::UncalibratedFloat,
        )
        .unwrap();
    assert!(
        prepared
            .decisions()
            .values()
            .all(|d| d.reset_reasons.contains(&ArcResetReason::PhaseDatumChanged))
    );
    tracker.finish_epoch(&used).unwrap();
}

#[test]
fn raw_float_and_relative_code_dsb_gauges_preserve_moving_positions_and_shift_real_ambiguities() {
    use gnss_rust::dual_frequency::PhaseBiasMode;
    let p = products();
    let mut expected = BTreeMap::new();
    for line in include_str!("fixtures/upstream_ppp_raw_phase.csv")
        .lines()
        .skip(1)
    {
        let v: Vec<_> = line.split(',').collect();
        expected.insert(
            v[0].parse::<SatelliteId>().unwrap(),
            v[2].parse::<f64>().unwrap(),
        );
    }
    for noisy in [false, true] {
        let mut calibrated = PppFilter::new(config()).unwrap();
        let mut c = config();
        c.filter.phase_bias_mode = PhaseBiasMode::UncalibratedFloat;
        let mut raw = PppFilter::new(c.clone()).unwrap();
        let mut dsb = PppFilter::new(c).unwrap();
        let s = selections();
        let rs = raw_selections();
        for e in epochs(noisy) {
            let base = solution(
                calibrated
                    .process_epoch(&e, &p, &biases(e.time), &s)
                    .unwrap(),
            );
            let raw = solution(
                raw.process_epoch(&e, &p, &code_only_biases(e.time), &rs)
                    .unwrap(),
            );
            let dsb = solution(dsb.process_epoch(&e, &p, &dsb_biases(e.time), &rs).unwrap());
            let reference = base.position.position_ecef.unwrap();
            let xyz = [reference.x, reference.y, reference.z];
            for r in [&raw, &dsb] {
                assert!(distance(r.position.position_ecef.unwrap(), &xyz) < 1e-6);
                assert!(!r.position.is_fixed());
                assert_eq!(r.phase_bias_mode, PhaseBiasMode::UncalibratedFloat);
                assert!(r.measurements.iter().all(|m| {
                    m.observation
                        .corrections
                        .iter()
                        .all(|c| c.phase_add_m.is_none() && c.phase_reference.is_none())
                }));
                for (b, a) in base.ambiguities.iter().zip(&r.ambiguities) {
                    assert_eq!(b.satellite, a.satellite);
                    assert!((a.metres - b.metres - expected[&a.satellite]).abs() < 1e-6);
                }
            }
        }
    }
}

fn add_code_fault(e: &mut gnss_rust::observation::ObservationEpoch, sat: SatelliteId, metres: f64) {
    for o in e.observations.iter_mut().filter(|o| o.satellite == sat) {
        if let Some(p) = &mut o.pseudorange_m {
            *p += metres;
        }
    }
}
#[test]
fn standardized_screening_is_opt_in_for_motion_with_finite_positive_sigma_only() {
    let original = acceleration_config();
    assert_eq!(original.filter.standardized_residual_gate_sigma, None);
    for gate in [0.0, -4.0, f64::NAN, f64::INFINITY] {
        let mut c = original.clone();
        c.filter.standardized_residual_gate_sigma = Some(gate);
        assert_eq!(
            PppFilter::new(c).unwrap_err(),
            PppError::InvalidConfiguration
        );
    }
    let mut c = original.clone();
    c.filter.position_model = PppPositionModel::Static;
    c.filter.standardized_residual_gate_sigma = Some(4.0);
    assert_eq!(
        PppFilter::new(c).unwrap_err(),
        PppError::InvalidConfiguration
    );
    for model in [
        original.filter.position_model,
        PppPositionModel::ConstantVelocity(Default::default()),
    ] {
        let mut c = original.clone();
        c.filter.position_model = model;
        c.filter.standardized_residual_gate_sigma = Some(4.0);
        PppFilter::new(c).unwrap();
    }
}
#[test]
fn native_sigma_screen_keeps_clean_and_noisy_car_states_identical_without_false_exclusions() {
    let p = products();
    let s = selections();
    for noisy in [false, true] {
        let base = acceleration_config();
        let mut c = base.clone();
        c.filter.standardized_residual_gate_sigma = Some(4.0);
        let mut actual = PppFilter::new(c).unwrap();
        let mut control = PppFilter::new(base).unwrap();
        for (i, e) in acceleration_epochs(noisy).iter().enumerate() {
            let r = solution(actual.process_epoch(e, &p, &biases(e.time), &s).unwrap());
            control.process_epoch(e, &p, &biases(e.time), &s).unwrap();
            assert!(
                r.residual_screening.is_empty(),
                "car sigma screen noisy={noisy} epoch {i}: {:?}",
                r.residual_screening
            );
            assert_eq!(actual.state(), control.state());
        }
    }
}
#[test]
fn single_and_multiple_code_faults_remove_whole_satellites_and_recompute_from_the_same_prior() {
    let ee = acceleration_epochs(false);
    let p = products();
    let s = selections();
    let bad1: SatelliteId = "G09".parse().unwrap();
    let bad2: SatelliteId = "G18".parse().unwrap();
    for multiple in [false, true] {
        let base = acceleration_config();
        let mut c = base.clone();
        c.filter.standardized_residual_gate_sigma = Some(4.0);
        let mut actual = PppFilter::new(c).unwrap();
        let mut control = PppFilter::new(base).unwrap();
        for e in &ee[..12] {
            actual.process_epoch(e, &p, &biases(e.time), &s).unwrap();
            control.process_epoch(e, &p, &biases(e.time), &s).unwrap();
        }
        let mut faulty = ee[12].clone();
        add_code_fault(&mut faulty, bad1, 40.0);
        if multiple {
            add_code_fault(&mut faulty, bad2, -80.0);
        }
        let mut omitted = ee[12].clone();
        omitted
            .observations
            .retain(|o| o.satellite != bad1 && (!multiple || o.satellite != bad2));
        let r = solution(
            actual
                .process_epoch(&faulty, &p, &biases(faulty.time), &s)
                .unwrap(),
        );
        control
            .process_epoch(&omitted, &p, &biases(omitted.time), &s)
            .unwrap();
        let excluded: std::collections::BTreeSet<_> =
            r.residual_screening.iter().map(|d| d.satellite).collect();
        assert_eq!(
            excluded,
            if multiple {
                [bad1, bad2].into()
            } else {
                [bad1].into()
            }
        );
        assert!(
            r.residual_screening
                .iter()
                .all(|d| d.removed && d.standardized_residual.abs() > 4.0)
        );
        assert!(
            r.measurements
                .iter()
                .all(|m| !excluded.contains(&m.satellite))
        );
        assert_eq!(actual.state(), control.state());
        assert!(actual.phase_arc(bad1).is_none());
        assert!(
            distance(
                r.position.position_ecef.unwrap(),
                &acceleration_truth()[12][3..6]
            ) < 0.01
        );
        let next = solution(
            actual
                .process_epoch(&ee[13], &p, &biases(ee[13].time), &s)
                .unwrap(),
        );
        assert_eq!(next.arc_commit.accepted[&bad1].accepted_epochs, 1);
    }
}
#[test]
fn screening_exhaustion_is_transactional_and_retention_never_keeps_identified_bad_arcs() {
    let ee = acceleration_epochs(false);
    let p = products();
    let s = selections();
    let bad: SatelliteId = "G09".parse().unwrap();
    let mut c = acceleration_config();
    c.filter.standardized_residual_gate_sigma = Some(4.0);
    c.filter.rejected_update_policy = PppRejectedUpdatePolicy::RetainContinuousArcs;
    let mut f = PppFilter::new(c).unwrap();
    for e in &ee[..12] {
        f.process_epoch(e, &p, &biases(e.time), &s).unwrap();
    }
    let old = f.state().clone();
    let time = f.state_time();
    let mut faulty = ee[12].clone();
    faulty.observations.retain(|o| o.satellite.prn() <= 9);
    add_code_fault(&mut faulty, bad, 400.0);
    assert!(matches!(
        f.process_epoch(&faulty, &p, &biases(faulty.time), &s),
        Err(PppError::ResidualScreeningExhausted { .. })
    ));
    assert_eq!(f.state(), &old);
    assert_eq!(f.state_time(), time);
    assert!(!f.last_residual_screening().last().unwrap().removed);
    assert!(!f.last_retained_phase_arcs().contains(&bad));
    assert!(f.phase_arc(bad).is_none());
    assert!(
        f.last_retained_phase_arcs()
            .iter()
            .all(|&sat| f.phase_arc(sat).unwrap().accepted_epochs == 12)
    );
    assert!(
        f.process_epoch(&ee[13], &p, &biases(ee[13].time), &s)
            .is_ok()
    );
}
#[test]
fn screened_retry_does_not_fuse_the_velocity_seed_more_than_once() {
    let e = &acceleration_epochs(false)[0];
    let p = products();
    let s = selections();
    let mut seed = code_motion_seed(e);
    // Deliberate covariance unit input, not measured Doppler or a truth seed.
    seed.solution.velocity_ecef = Some([1.0, -2.0, 3.0]);
    seed.solution.velocity_covariance = Some([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]]);
    let mut faulty = e.clone();
    add_code_fault(&mut faulty, "G09".parse().unwrap(), 400.0);
    let mut c = acceleration_config();
    c.filter.motion_seed_policy = Some(Default::default());
    c.filter.standardized_residual_gate_sigma = Some(4.0);
    let mut f = PppFilter::new(c).unwrap();
    let r = solution(
        f.process_epoch_with_seed(&faulty, &p, &biases(e.time), &s, Some(&seed))
            .unwrap(),
    );
    assert!(!r.residual_screening.is_empty());
    assert!(r.residual_screening.iter().all(|d| d.removed));
    assert_eq!(r.position.velocity_ecef, seed.solution.velocity_ecef);
    for i in 5..8 {
        assert!((f.state().covariance[i][i] - 100.0 * 4.0 / 104.0).abs() < 1e-12);
    }
    assert!(r.motion_seed.unwrap().committed);
}

#[test]
fn reference_if_initialization_is_used_once_per_new_or_reset_arc() {
    let p = products();
    let s = selections();
    let mut ee = acceleration_epochs(false);
    let mut c = acceleration_config();
    assert_eq!(
        c.filter.ambiguity_initialization,
        PppAmbiguityInitialization::IndependentZeroMean
    );
    c.filter.ambiguity_initialization = PppAmbiguityInitialization::ReferencePrecise;
    c.filter.motion_seed_policy = Some(PppMotionSeedPolicy {
        initialize_clock_from_seed: true,
        ..Default::default()
    });
    let mut f = PppFilter::new(c).unwrap();
    let r = solution(
        f.process_epoch_with_seed(
            &ee[0],
            &p,
            &biases(ee[0].time),
            &s,
            Some(&code_motion_seed(&ee[0])),
        )
        .unwrap(),
    );
    assert_eq!(r.ambiguity_seeds.len(), 9);
    for seed in &r.ambiguity_seeds {
        assert_eq!(
            seed.mean_m,
            seed.corrected_phase_if_m - seed.predicted_range_m
        );
        assert_eq!(seed.variance_m2, 25.0);
        assert!(seed.committed);
    }
    assert_eq!(f.last_ambiguity_seeds(), r.ambiguity_seeds);
    assert_eq!(
        r.motion_seed.as_ref().unwrap().initial_clock_seed_s,
        Some(code_motion_seed(&ee[0]).solution.receiver_clock_bias)
    );
    let r = solution(
        f.process_epoch_with_seed(
            &ee[1],
            &p,
            &biases(ee[1].time),
            &s,
            Some(&code_motion_seed(&ee[1])),
        )
        .unwrap(),
    );
    assert!(r.ambiguity_seeds.is_empty());
    assert_eq!(r.motion_seed.as_ref().unwrap().initial_clock_seed_s, None);
    let reset: SatelliteId = "G09".parse().unwrap();
    for o in ee[2]
        .observations
        .iter_mut()
        .filter(|o| o.satellite == reset)
    {
        o.lli = Some(1);
    }
    let r = solution(
        f.process_epoch_with_seed(
            &ee[2],
            &p,
            &biases(ee[2].time),
            &s,
            Some(&code_motion_seed(&ee[2])),
        )
        .unwrap(),
    );
    assert_eq!(r.ambiguity_seeds.len(), 1);
    assert_eq!(r.ambiguity_seeds[0].satellite, reset);
    assert_eq!(r.arc_commit.accepted[&reset].accepted_epochs, 1);
}
#[test]
fn failed_reference_ambiguity_initialization_epoch_never_commits_the_proposed_state() {
    let p = products();
    let s = selections();
    let mut ee = acceleration_epochs(false);
    let mut c = acceleration_config();
    c.filter.ambiguity_initialization = PppAmbiguityInitialization::ReferencePrecise;
    c.filter.motion_seed_policy = Some(PppMotionSeedPolicy {
        initialize_clock_from_seed: true,
        ..Default::default()
    });
    let mut f = PppFilter::new(c).unwrap();
    f.process_epoch_with_seed(
        &ee[0],
        &p,
        &biases(ee[0].time),
        &s,
        Some(&code_motion_seed(&ee[0])),
    )
    .unwrap();
    let old = f.state().clone();
    let time = f.state_time();
    let reset: SatelliteId = "G09".parse().unwrap();
    for o in ee[1]
        .observations
        .iter_mut()
        .filter(|o| o.satellite == reset)
    {
        o.lli = Some(1);
    }
    let seed = code_motion_seed(&ee[1]);
    inject_receiver_clock_metres(&mut ee[1], 1e7);
    assert!(
        f.process_epoch_with_seed(&ee[1], &p, &biases(ee[1].time), &s, Some(&seed))
            .is_err()
    );
    assert_eq!(f.state(), &old);
    assert_eq!(f.state_time(), time);
    assert!(!f.last_ambiguity_seeds().is_empty());
    assert!(f.last_ambiguity_seeds().iter().all(|seed| !seed.committed));
}

#[test]
fn explicit_reference_initial_clock_rejects_nonfinite_seed_before_mutating_state() {
    let e = &acceleration_epochs(false)[0];
    let p = products();
    let s = selections();
    for clock in [f64::NAN, f64::INFINITY, f64::MAX] {
        let mut c = acceleration_config();
        c.filter.motion_seed_policy = Some(PppMotionSeedPolicy {
            initialize_clock_from_seed: true,
            ..Default::default()
        });
        let mut f = PppFilter::new(c).unwrap();
        let old = f.state().clone();
        let mut seed = code_motion_seed(e);
        seed.solution.receiver_clock_bias = clock;
        assert_eq!(
            f.process_epoch_with_seed(e, &p, &biases(e.time), &s, Some(&seed)),
            Err(PppError::MotionSeed(PppMotionSeedError::InvalidClock))
        );
        assert_eq!(f.state(), &old);
        assert_eq!(f.state_time(), None);
    }
}
