use gnss_rust::observation::ObservationEpoch;
use gnss_rust::rinex::{RinexObservationReader, read_navigation};
use gnss_rust::rtk::{RtkConstellations, RtkFloatConfig, RtkFloatFilter, RtkPositionMode};
use gnss_rust::rtk_fix::{
    FixDecision, FixError, KinematicMotionConfig, ReferenceGateMode, RtkFixConfig, RtkFixPolicy,
    RtkMotionEvidence,
};
use gnss_rust::{EcefCoord, GnssTime, SolutionStatus};
use std::io::Cursor;

struct Data {
    rover: Vec<ObservationEpoch>,
    base: Vec<ObservationEpoch>,
    nav: gnss_rust::navigation::NavigationData,
    base_position: EcefCoord,
    seed: EcefCoord,
}
fn data() -> Data {
    let mut r = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_mixed_rover.obs"
    )))
    .unwrap();
    let mut b = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_mixed_base.obs"
    )))
    .unwrap();
    let seed = r.header().approximate_position.unwrap();
    let base_position = b.header().approximate_position.unwrap();
    let rover = r.by_ref().collect::<Result<_, _>>().unwrap();
    let base = b.by_ref().collect::<Result<_, _>>().unwrap();
    let nav = read_navigation(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_mixed.nav"
    )))
    .unwrap();
    Data {
        rover,
        base,
        nav,
        base_position,
        seed,
    }
}
fn truth(i: usize) -> EcefCoord {
    let v: Vec<f64> = include_str!("fixtures/synthetic_rtk_mixed_truth.csv")
        .lines()
        .nth(i + 1)
        .unwrap()
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect();
    EcefCoord::new(v[2], v[3], v[4])
}
fn motion(time: GnssTime) -> RtkMotionEvidence {
    let a = truth(0);
    let b = truth(1);
    RtkMotionEvidence {
        time,
        velocity_mps: [b.x - a.x, b.y - a.y, b.z - a.z],
        covariance_m2s2: [[1e-4, 0.0, 0.0], [0.0, 1e-4, 0.0], [0.0, 0.0, 1e-4]],
    }
}
fn filter(d: &Data) -> RtkFloatFilter {
    RtkFloatFilter::new(
        d.base_position,
        d.seed,
        RtkFloatConfig {
            constellations: RtkConstellations::GpsGalileo,
            position_mode: RtkPositionMode::Kinematic,
            ..Default::default()
        },
    )
    .unwrap()
}
fn policy(mode: ReferenceGateMode) -> RtkFixPolicy {
    RtkFixPolicy::for_kinematic(
        RtkFixConfig::default(),
        KinematicMotionConfig::default(),
        mode,
    )
    .unwrap()
}
fn distance(a: EcefCoord, b: EcefCoord) -> f64 {
    (a.x - b.x).hypot(a.y - b.y).hypot(a.z - b.z)
}

#[test]
fn moving_mixed_policy_confirms_correct_integers_in_both_reference_modes_without_feedback() {
    for mode in [
        ReferenceGateMode::Selected,
        ReferenceGateMode::AllReferences,
    ] {
        let d = data();
        let mut engine = filter(&d);
        let mut policy = policy(mode);
        let mut fixed = 0;
        for i in 0..d.rover.len() {
            let float = engine
                .process_pair(&d.rover[i], &d.base[i], &d.nav)
                .unwrap();
            let before = engine.state().unwrap().clone();
            let attempt = policy
                .evaluate_kinematic(&float, Some(&motion(d.rover[i].time)))
                .unwrap();
            assert_eq!(
                attempt.fixed_solution.is_some(),
                i >= 5,
                "epoch {i}: {:?}",
                attempt.decision
            );
            if let Some(solution) = attempt.fixed_solution {
                fixed += 1;
                assert_eq!(solution.status, SolutionStatus::Fixed);
                assert!(
                    distance(solution.position_ecef.unwrap(), truth(i)) < 0.003,
                    "epoch {i}"
                );
                assert!(
                    attempt
                        .motion_checks
                        .iter()
                        .all(|c| c.residual_norm_m <= c.maximum_m)
                );
                assert_eq!(
                    attempt.system_reference_audit.is_some(),
                    mode == ReferenceGateMode::AllReferences
                );
            }
            assert_eq!(float.position.status, SolutionStatus::Float);
            assert_eq!(engine.state(), Some(&before));
        }
        assert_eq!(fixed, 27);
    }
}

#[test]
fn missing_uncertain_and_reset_motion_require_fresh_confirmation() {
    let d = data();
    let mut engine = filter(&d);
    let mut policy = policy(ReferenceGateMode::AllReferences);
    for i in 0..12 {
        let float = engine
            .process_pair(&d.rover[i], &d.base[i], &d.nav)
            .unwrap();
        let mut v = motion(d.rover[i].time);
        if i == 9 {
            v.covariance_m2s2 = [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]];
        }
        let attempt = policy
            .evaluate_kinematic(&float, if i == 7 { None } else { Some(&v) })
            .unwrap();
        match i {
            5 | 6 | 11 => assert_eq!(attempt.decision, FixDecision::Accepted),
            7 => assert_eq!(attempt.decision, FixDecision::MotionUnavailable),
            9 => assert_eq!(attempt.decision, FixDecision::MotionUncertain),
            8 | 10 => {
                assert_eq!(attempt.decision, FixDecision::Confirming);
                assert_eq!(attempt.confirmation_count, 1);
            }
            _ => assert!(attempt.fixed_solution.is_none()),
        }
    }
    policy.reset();
    let float = engine
        .process_pair(&d.rover[12], &d.base[12], &d.nav)
        .unwrap();
    assert_eq!(
        policy
            .evaluate_kinematic(&float, Some(&motion(d.rover[12].time)))
            .unwrap()
            .decision,
        FixDecision::InsufficientLock
    );
}

#[test]
fn inconsistent_velocity_rejects_a_position_jump_and_recovers_after_reconfirmation() {
    let d = data();
    let mut engine = filter(&d);
    let mut policy = policy(ReferenceGateMode::Selected);
    for i in 0..9 {
        let float = engine
            .process_pair(&d.rover[i], &d.base[i], &d.nav)
            .unwrap();
        let before = engine.state().unwrap().clone();
        let mut v = motion(d.rover[i].time);
        if i == 6 {
            v.velocity_mps[0] += 1.0;
        }
        let attempt = policy.evaluate_kinematic(&float, Some(&v)).unwrap();
        if i == 6 {
            assert_eq!(attempt.decision, FixDecision::HistoryJump);
            assert!(attempt.fixed_solution.is_none());
            assert!(attempt.motion_checks[0].residual_norm_m > 0.4);
        }
        if i == 7 {
            assert_eq!(attempt.decision, FixDecision::Confirming);
            assert_eq!(attempt.confirmation_count, 1);
        }
        if i == 8 {
            assert_eq!(attempt.decision, FixDecision::Accepted);
        }
        assert_eq!(engine.state(), Some(&before));
    }
}

#[test]
fn invalid_motion_and_wrong_entry_points_never_advance_a_policy_epoch() {
    let d = data();
    let mut engine = filter(&d);
    let float = engine
        .process_pair(&d.rover[0], &d.base[0], &d.nav)
        .unwrap();
    let mut kinematic = policy(ReferenceGateMode::Selected);
    let mut static_policy = RtkFixPolicy::new(RtkFixConfig::default()).unwrap();
    assert_eq!(
        kinematic.evaluate(&float),
        Err(FixError::MotionEvidenceRequired)
    );
    assert_eq!(
        static_policy.evaluate_kinematic(&float, Some(&motion(d.rover[0].time))),
        Err(FixError::UnsupportedPositionMode)
    );
    let mut invalid = motion(d.rover[0].time);
    invalid.time = invalid.time.checked_add_seconds(1.0).unwrap();
    assert_eq!(
        kinematic.evaluate_kinematic(&float, Some(&invalid)),
        Err(FixError::InvalidMotionEvidence)
    );
    for value in [f64::NAN, f64::INFINITY, f64::MAX] {
        invalid = motion(d.rover[0].time);
        invalid.velocity_mps[0] = value;
        assert_eq!(
            kinematic.evaluate_kinematic(&float, Some(&invalid)),
            Err(FixError::InvalidMotionEvidence)
        );
    }
    for value in [-1.0, f64::NAN] {
        invalid = motion(d.rover[0].time);
        invalid.covariance_m2s2[0][0] = value;
        assert_eq!(
            kinematic.evaluate_kinematic(&float, Some(&invalid)),
            Err(FixError::InvalidMotionEvidence)
        );
    }
    assert_eq!(
        kinematic
            .evaluate_kinematic(&float, Some(&motion(d.rover[0].time)))
            .unwrap()
            .decision,
        FixDecision::InsufficientLock
    );
    assert_eq!(
        kinematic.evaluate_kinematic(&float, Some(&motion(d.rover[0].time))),
        Err(FixError::OutOfOrder)
    );
}

#[test]
fn motion_gaps_restart_confirmation_and_lli_requires_fresh_lock() {
    let mut d = data();
    let mut engine = filter(&d);
    let mut policy = policy(ReferenceGateMode::AllReferences);
    for i in 0..6 {
        let float = engine
            .process_pair(&d.rover[i], &d.base[i], &d.nav)
            .unwrap();
        policy
            .evaluate_kinematic(&float, Some(&motion(d.rover[i].time)))
            .unwrap();
    }
    for i in 8..16 {
        if i == 10 {
            d.rover[i].observations[0].lli = Some(1);
        }
        let float = engine
            .process_pair(&d.rover[i], &d.base[i], &d.nav)
            .unwrap();
        let attempt = policy
            .evaluate_kinematic(&float, Some(&motion(d.rover[i].time)))
            .unwrap();
        match i {
            8 | 14 => {
                assert_eq!(attempt.decision, FixDecision::Confirming);
                assert_eq!(attempt.confirmation_count, 1);
            }
            9 | 15 => assert_eq!(attempt.decision, FixDecision::Accepted),
            10..=13 => assert_eq!(attempt.decision, FixDecision::InsufficientLock),
            _ => unreachable!(),
        }
    }
}

#[test]
fn kinematic_configuration_and_uncertainty_boundary_are_explicit() {
    let cfg = RtkFixConfig::default();
    let mode = ReferenceGateMode::Selected;
    assert!(matches!(
        RtkFixPolicy::for_kinematic(
            RtkFixConfig {
                min_confirmation_epochs: 1,
                ..cfg
            },
            KinematicMotionConfig::default(),
            mode
        ),
        Err(FixError::InvalidConfiguration)
    ));
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            RtkFixPolicy::for_kinematic(
                cfg,
                KinematicMotionConfig {
                    max_velocity_std_mps: bad,
                    ..Default::default()
                },
                mode
            ),
            Err(FixError::InvalidConfiguration)
        ));
        assert!(matches!(
            RtkFixPolicy::for_kinematic(
                cfg,
                KinematicMotionConfig {
                    max_interval_s: bad,
                    ..Default::default()
                },
                mode
            ),
            Err(FixError::InvalidConfiguration)
        ));
    }
    let d = data();
    let mut engine = filter(&d);
    let float = engine
        .process_pair(&d.rover[0], &d.base[0], &d.nav)
        .unwrap();
    let mut v = motion(d.rover[0].time);
    v.covariance_m2s2 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut exact = policy(mode);
    assert_eq!(
        exact.evaluate_kinematic(&float, Some(&v)).unwrap().decision,
        FixDecision::InsufficientLock
    );
    let mut lower = RtkFixPolicy::for_kinematic(
        cfg,
        KinematicMotionConfig {
            max_velocity_std_mps: f64::from_bits(1.0_f64.to_bits() - 1),
            ..Default::default()
        },
        mode,
    )
    .unwrap();
    assert_eq!(
        lower.evaluate_kinematic(&float, Some(&v)).unwrap().decision,
        FixDecision::MotionUncertain
    );
}

#[test]
fn gps_only_kinematic_path_preserves_the_original_moving_truth_tolerance() {
    let mut r = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_moving_rover.obs"
    )))
    .unwrap();
    let mut b = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_moving_base.obs"
    )))
    .unwrap();
    let seed = r.header().approximate_position.unwrap();
    let base = b.header().approximate_position.unwrap();
    let nav = read_navigation(Cursor::new(include_str!("fixtures/synthetic_rtk.nav"))).unwrap();
    let points: Vec<Vec<f64>> = include_str!("fixtures/synthetic_rtk_moving_truth.csv")
        .lines()
        .skip(1)
        .map(|l| l.split(',').map(|v| v.parse().unwrap()).collect())
        .collect();
    let velocity = std::array::from_fn(|j| points[1][2 + j] - points[0][2 + j]);
    let mut engine = RtkFloatFilter::new(
        base,
        seed,
        RtkFloatConfig {
            position_mode: RtkPositionMode::Kinematic,
            ..Default::default()
        },
    )
    .unwrap();
    let mut policy = policy(ReferenceGateMode::Selected);
    let mut count = 0;
    for (i, (ro, ba)) in r.by_ref().zip(b.by_ref()).enumerate() {
        let ro = ro.unwrap();
        let float = engine.process_pair(&ro, &ba.unwrap(), &nav).unwrap();
        let motion = RtkMotionEvidence {
            time: ro.time,
            velocity_mps: velocity,
            covariance_m2s2: [[1e-4, 0.0, 0.0], [0.0, 1e-4, 0.0], [0.0, 0.0, 1e-4]],
        };
        let attempt = policy.evaluate_kinematic(&float, Some(&motion)).unwrap();
        if let Some(fixed) = attempt.fixed_solution {
            count += 1;
            assert!(
                distance(
                    fixed.position_ecef.unwrap(),
                    EcefCoord::new(points[i][2], points[i][3], points[i][4])
                ) < 0.003
            );
        }
    }
    assert_eq!(count, 59);
}
