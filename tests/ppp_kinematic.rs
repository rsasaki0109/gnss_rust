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
