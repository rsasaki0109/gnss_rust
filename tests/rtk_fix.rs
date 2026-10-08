use gnss_rust::filter::{FilterError, GaussianState, Matrix};
use gnss_rust::lambda::LambdaConfig;
use gnss_rust::navigation::NavigationData;
use gnss_rust::observation::ObservationEpoch;
use gnss_rust::rinex::{RinexObservationReader, read_navigation};
use gnss_rust::rtk::{RtkFloatConfig, RtkFloatFilter};
use gnss_rust::rtk_fix::{
    CodePhaseObservation, FixDecision, FixError, RtkFixConfig, RtkFixPolicy, code_phase_check,
    condition_head, validate_candidate,
};
use gnss_rust::{EcefCoord, SatelliteId, SolutionStatus};
use std::{collections::BTreeMap, io::Cursor};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1.0),
        "{a:.17} != {b:.17}"
    );
}
fn take(v: &[f64], offset: &mut usize, n: usize) -> Vec<f64> {
    let values = v[*offset..*offset + n].to_vec();
    *offset += n;
    values
}
fn matrix(v: &[f64], offset: &mut usize, rows: usize, cols: usize) -> Matrix {
    take(v, offset, rows * cols)
        .chunks(cols)
        .map(|r| r.to_vec())
        .collect()
}

#[test]
fn head_conditioning_matches_native_cpp_and_independent_eigen_covariance() {
    let mut cases = 0;
    for line in include_str!("fixtures/upstream_rtk_fix.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let mut fields = line.split(',');
        let name = fields.next().unwrap();
        let h = fields.next().unwrap().parse().unwrap();
        let k = fields.next().unwrap().parse().unwrap();
        let values: Vec<f64> = fields.map(|v| v.parse().unwrap()).collect();
        let mut offset = 0;
        let head = GaussianState {
            mean: take(&values, &mut offset, h),
            covariance: matrix(&values, &mut offset, h, h),
        };
        let qab = matrix(&values, &mut offset, h, k);
        let qbb = matrix(&values, &mut offset, k, k);
        let float = take(&values, &mut offset, k);
        let integers: Vec<_> = take(&values, &mut offset, k)
            .iter()
            .map(|v| *v as i64)
            .collect();
        let expected_mean = take(&values, &mut offset, h);
        let expected_cov = matrix(&values, &mut offset, h, h);
        let expected_cost = take(&values, &mut offset, 1)[0];
        let result = condition_head(&head, &float, &qbb, &qab, &integers).unwrap();
        for (a, b) in result.state.mean.iter().zip(expected_mean) {
            close(*a, b, 1e-11);
        }
        for (a, b) in result
            .state
            .covariance
            .iter()
            .flatten()
            .zip(expected_cov.iter().flatten())
        {
            close(*a, *b, 1e-11);
        }
        close(result.squared_ambiguity_residual, expected_cost, 1e-11);
        assert_eq!(offset, values.len(), "{name}");
        cases += 1;
    }
    assert_eq!(cases, 5);
}

#[test]
fn code_phase_gate_matches_native_cpp_including_exact_boundary() {
    let mut cases = 0;
    for line in include_str!("fixtures/upstream_rtk_code_phase.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let mut fields = line.split(',');
        fields.next().unwrap();
        let k: usize = fields.next().unwrap().parse().unwrap();
        let threshold = fields.next().unwrap().parse().unwrap();
        let values: Vec<f64> = fields.map(|v| v.parse().unwrap()).collect();
        let observations: Vec<_> = values[..3 * k]
            .chunks(3)
            .map(|v| CodePhaseObservation {
                dd_code_m: v[0],
                dd_carrier_m: v[1],
                fixed_ambiguity_m: v[2],
            })
            .collect();
        let check = code_phase_check(&observations, threshold).unwrap();
        assert_eq!(check.consistent, values[3 * k] != 0.0);
        assert_eq!(check.bad_pairs, values[3 * k + 1] as usize);
        close(check.rms_innovation_m, values[3 * k + 2], 1e-14);
        close(check.max_abs_innovation_m, values[3 * k + 3], 1e-14);
        cases += 1;
    }
    assert_eq!(cases, 3);
}

#[test]
fn scalar_conditioning_and_invalid_joint_inputs() {
    let head = GaussianState {
        mean: vec![0.0],
        covariance: vec![vec![2.0]],
    };
    let out = condition_head(&head, &[0.25], &vec![vec![1.0]], &vec![vec![0.5]], &[0]).unwrap();
    close(out.state.mean[0], -0.125, 1e-15);
    close(out.state.covariance[0][0], 1.75, 1e-15);
    close(out.squared_ambiguity_residual, 0.0625, 1e-15);
    assert!(condition_head(&head, &[0.0], &vec![vec![1.0]], &vec![vec![2.0]], &[0]).is_err());
    assert_eq!(
        condition_head(&head, &[f64::NAN], &vec![vec![1.0]], &vec![vec![0.5]], &[0]),
        Err(FixError::Filter(FilterError::NonFinite))
    );
    assert_eq!(
        condition_head(
            &head,
            &[0.0],
            &vec![vec![1.0]],
            &vec![vec![0.5]],
            &[1_i64 << 52]
        ),
        Err(FixError::InvalidInteger)
    );
    assert!(condition_head(&head, &[0.0], &vec![vec![1.0]], &vec![vec![0.5]], &[]).is_err());
    let observation = CodePhaseObservation {
        dd_code_m: f64::NAN,
        dd_carrier_m: 0.0,
        fixed_ambiguity_m: 0.0,
    };
    assert_eq!(
        code_phase_check(&[observation], 0.5),
        Err(FixError::Filter(FilterError::NonFinite))
    );
    assert!(code_phase_check(&[], 0.5).is_err());
    assert!(
        RtkFixPolicy::new(RtkFixConfig {
            min_confirmation_epochs: 0,
            ..Default::default()
        })
        .is_err()
    );
}

struct Fixture {
    nav: NavigationData,
    rover: Vec<ObservationEpoch>,
    base: Vec<ObservationEpoch>,
    base_position: EcefCoord,
    seed: EcefCoord,
    truth: EcefCoord,
    sd: BTreeMap<SatelliteId, i64>,
}
fn fixture(noisy: bool) -> Fixture {
    let mut lines = include_str!("fixtures/synthetic_rtk_truth.csv")
        .lines()
        .filter(|l| !l.starts_with('#'));
    let p: Vec<f64> = lines
        .next()
        .unwrap()
        .split(',')
        .map(|v| v.parse().unwrap())
        .collect();
    let coord = |i| EcefCoord::new(p[i], p[i + 1], p[i + 2]);
    let read = |s| {
        RinexObservationReader::new(Cursor::new(s))
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap()
    };
    Fixture {
        nav: read_navigation(Cursor::new(include_str!("fixtures/synthetic_rtk.nav"))).unwrap(),
        rover: read(if noisy {
            include_str!("fixtures/synthetic_rtk_rover_noisy.obs")
        } else {
            include_str!("fixtures/synthetic_rtk_rover.obs")
        }),
        base: read(if noisy {
            include_str!("fixtures/synthetic_rtk_base_noisy.obs")
        } else {
            include_str!("fixtures/synthetic_rtk_base.obs")
        }),
        base_position: coord(0),
        truth: coord(3),
        seed: coord(6),
        sd: lines
            .map(|l| {
                let (s, n) = l.split_once(',').unwrap();
                (s.parse().unwrap(), n.parse().unwrap())
            })
            .collect(),
    }
}
fn filter(f: &Fixture) -> RtkFloatFilter {
    RtkFloatFilter::new(f.base_position, f.seed, RtkFloatConfig::default()).unwrap()
}
fn distance(a: EcefCoord, b: EcefCoord) -> f64 {
    (a.x - b.x).hypot(a.y - b.y).hypot(a.z - b.z)
}

#[test]
fn known_integers_produce_validated_fix_without_modifying_float_posterior() {
    for noisy in [false, true] {
        let f = fixture(noisy);
        let mut engine = filter(&f);
        let mut policy = RtkFixPolicy::new(RtkFixConfig::default()).unwrap();
        let mut maximum_error = 0.0_f64;
        let mut fixed = 0;
        for (i, (r, b)) in f.rover.iter().zip(&f.base).enumerate() {
            let float = engine.process_pair(r, b, &f.nav).unwrap();
            let prior = engine.state().unwrap().clone();
            let attempt = policy.evaluate(&float).unwrap();
            assert_eq!(*engine.state().unwrap(), prior);
            assert_eq!(float.position.status, SolutionStatus::Float);
            if i < 4 {
                assert_eq!(attempt.decision, FixDecision::InsufficientLock);
            } else if i == 4 {
                assert_eq!(
                    attempt.decision,
                    FixDecision::Confirming,
                    "noisy={noisy}: {attempt:?}"
                );
            } else {
                assert_eq!(
                    attempt.decision,
                    FixDecision::Accepted,
                    "epoch={i} noisy={noisy}: {attempt:?}"
                );
                let solution = attempt.fixed_solution.unwrap();
                assert_eq!(solution.status, SolutionStatus::Fixed);
                let expected: Vec<_> = float
                    .ambiguity_satellites
                    .iter()
                    .map(|s| f.sd[&float.reference_satellite] - f.sd[s])
                    .collect();
                assert_eq!(attempt.integers, expected);
                let error = distance(solution.position_ecef.unwrap(), f.truth);
                maximum_error = maximum_error.max(error);
                assert!(error < if noisy { 0.005 } else { 0.003 }, "{error}");
                for j in 0..3 {
                    let covariance = solution.position_covariance.unwrap()[j][j];
                    assert!(
                        covariance > 0.0
                            && covariance < float.position.position_covariance.unwrap()[j][j]
                    );
                }
                fixed += 1;
            }
        }
        assert_eq!(fixed, 27);
        eprintln!("noisy={noisy}: 27 FIX, maximum 3-D error {maximum_error:.9} m");
    }
}

#[test]
fn wrong_integer_and_strict_gates_reject_without_destroying_float() {
    let f = fixture(true);
    let mut engine = filter(&f);
    let mut last = None;
    for i in 0..10 {
        last = Some(
            engine
                .process_pair(&f.rover[i], &f.base[i], &f.nav)
                .unwrap(),
        );
    }
    let float = last.unwrap();
    let state = engine.state().unwrap().clone();
    let mut integers = float
        .integer_candidates(LambdaConfig::default())
        .unwrap()
        .candidates[0]
        .ambiguities
        .clone();
    assert_eq!(
        validate_candidate(&float, &integers, RtkFixConfig::default())
            .unwrap()
            .rejection,
        None
    );
    integers[0] += 1;
    let bad = validate_candidate(&float, &integers, RtkFixConfig::default()).unwrap();
    assert!(bad.rejection.is_some(), "wrong candidate {bad:?}");
    assert!(bad.phase_residual_max_abs_m > 0.03 || bad.float_correction_m > 1.0);
    assert_eq!(*engine.state().unwrap(), state);
    let configurations = [
        (
            RtkFixConfig {
                ratio_threshold: 1e30,
                ..Default::default()
            },
            FixDecision::RatioRejected,
        ),
        (
            RtkFixConfig {
                max_position_std_m: 1e-9,
                ..Default::default()
            },
            FixDecision::CovarianceTooLarge,
        ),
        (
            RtkFixConfig {
                max_phase_residual_rms_m: 1e-10,
                ..Default::default()
            },
            FixDecision::PhaseResiduals,
        ),
        (
            RtkFixConfig {
                max_code_phase_innovation_m: 1e-10,
                ..Default::default()
            },
            FixDecision::CodePhaseInconsistency,
        ),
    ];
    for (config, expected) in configurations {
        let mut policy = RtkFixPolicy::new(RtkFixConfig {
            min_lock_epochs: 1,
            ..config
        })
        .unwrap();
        let attempt = policy.evaluate(&float).unwrap();
        assert_eq!(attempt.decision, expected);
        assert!(attempt.fixed_solution.is_none());
        assert_eq!(*engine.state().unwrap(), state);
    }
}

#[test]
fn duplicate_evaluation_does_not_confirm_and_lli_slip_returns_to_float() {
    let f = fixture(false);
    let mut engine = filter(&f);
    let mut policy = RtkFixPolicy::new(RtkFixConfig::default()).unwrap();
    for i in 0..6 {
        let float = engine
            .process_pair(&f.rover[i], &f.base[i], &f.nav)
            .unwrap();
        let attempt = policy.evaluate(&float).unwrap();
        assert_eq!(policy.evaluate(&float), Err(FixError::OutOfOrder));
        if i == 4 {
            assert_eq!(attempt.confirmation_count, 1);
            assert!(attempt.fixed_solution.is_none());
        }
        if i == 5 {
            assert!(attempt.fixed_solution.is_some());
        }
    }
    let satellite = f.rover[6].observations[0].satellite;
    for i in 6..18 {
        let mut rover = f.rover[i].clone();
        *rover.observations[0].carrier_phase_cycles.as_mut().unwrap() += 7.0;
        if i == 6 {
            rover.observations[0].lli = Some(1);
        }
        let float = engine.process_pair(&rover, &f.base[i], &f.nav).unwrap();
        let attempt = policy.evaluate(&float).unwrap();
        if i < 10 {
            assert_eq!(attempt.decision, FixDecision::InsufficientLock);
        } else if i == 10 {
            assert_eq!(attempt.decision, FixDecision::Confirming);
        } else {
            assert_eq!(
                attempt.decision,
                FixDecision::Accepted,
                "epoch {i}: {attempt:?}"
            );
            let sd = |s| f.sd[&s] + if s == satellite { 7 } else { 0 };
            let expected: Vec<_> = float
                .ambiguity_satellites
                .iter()
                .map(|s| sd(float.reference_satellite) - sd(*s))
                .collect();
            assert_eq!(attempt.integers, expected);
            assert!(
                distance(
                    attempt.fixed_solution.unwrap().position_ecef.unwrap(),
                    f.truth
                ) < 0.003
            );
        }
    }
}

#[test]
fn high_ratio_with_inconsistent_phase_still_does_not_fix() {
    let f = fixture(false);
    let mut engine = filter(&f);
    let mut policy = RtkFixPolicy::new(RtkFixConfig::default()).unwrap();
    let mut high_ratio_rejections = 0;
    for i in 0..32 {
        let mut rover = f.rover[i].clone();
        // A transient outlier after convergence: integer history may remain
        // strong while the current carrier-phase observation is inconsistent.
        if i == 6 {
            *rover.observations[0].carrier_phase_cycles.as_mut().unwrap() += 0.25;
        }
        let float = engine.process_pair(&rover, &f.base[i], &f.nav).unwrap();
        let attempt = policy.evaluate(&float).unwrap();
        if i == 6 {
            assert!(attempt.fixed_solution.is_none(), "epoch {i}: {attempt:?}");
            assert!(attempt.ratio.is_some_and(|r| r >= 3.0), "{attempt:?}");
            assert_eq!(attempt.decision, FixDecision::PhaseResiduals);
            high_ratio_rejections += 1;
        } else if i == 7 {
            assert_eq!(attempt.decision, FixDecision::Confirming);
        } else if i >= 8 {
            assert_eq!(attempt.decision, FixDecision::Accepted);
        }
    }
    assert!(high_ratio_rejections > 0);
}

#[test]
fn persistent_fractional_phase_bias_is_not_promoted_to_fix() {
    let f = fixture(false);
    let mut engine = filter(&f);
    let mut policy = RtkFixPolicy::new(RtkFixConfig::default()).unwrap();
    for i in 0..32 {
        let mut rover = f.rover[i].clone();
        *rover.observations[0].carrier_phase_cycles.as_mut().unwrap() += 0.25;
        let float = engine.process_pair(&rover, &f.base[i], &f.nav).unwrap();
        let attempt = policy.evaluate(&float).unwrap();
        assert!(attempt.fixed_solution.is_none());
        if i >= 4 {
            assert_eq!(attempt.decision, FixDecision::RatioRejected);
        }
    }
}

#[test]
fn static_history_discontinuity_rejects_an_otherwise_valid_candidate() {
    let f = fixture(true);
    let mut engine = filter(&f);
    // Enforce a tighter static-position history limit than the noisy epochs
    // satisfy, while keeping all local residual/covariance tests unchanged.
    let mut policy = RtkFixPolicy::new(RtkFixConfig {
        max_history_jump_m: 1e-5,
        min_confirmation_epochs: 1,
        ..Default::default()
    })
    .unwrap();
    let mut rejected = 0;
    for i in 0..32 {
        let float = engine
            .process_pair(&f.rover[i], &f.base[i], &f.nav)
            .unwrap();
        let attempt = policy.evaluate(&float).unwrap();
        if attempt.decision == FixDecision::HistoryJump {
            assert!(attempt.fixed_solution.is_none());
            assert_eq!(attempt.validation.unwrap().rejection, None);
            rejected += 1;
        }
    }
    assert!(rejected > 0);
}

#[test]
fn reference_loss_reconfirms_and_reacquisition_waits_for_a_new_lock() {
    let f = fixture(false);
    let mut engine = filter(&f);
    let mut policy = RtkFixPolicy::new(RtkFixConfig::default()).unwrap();
    let mut reference = None;
    for i in 0..16 {
        let mut rover = f.rover[i].clone();
        if (6..8).contains(&i) {
            rover
                .observations
                .retain(|o| Some(o.satellite) != reference);
        }
        let float = engine.process_pair(&rover, &f.base[i], &f.nav).unwrap();
        if i == 5 {
            reference = Some(float.reference_satellite);
        }
        let attempt = policy.evaluate(&float).unwrap();
        if i == 6 {
            assert_ne!(Some(float.reference_satellite), reference);
            assert_eq!(attempt.decision, FixDecision::Confirming);
        } else if i == 7 {
            assert_eq!(attempt.decision, FixDecision::Accepted);
        } else if (8..12).contains(&i) {
            assert_eq!(attempt.decision, FixDecision::InsufficientLock);
        } else if i == 12 {
            assert_eq!(attempt.decision, FixDecision::Confirming);
        } else if i >= 13 {
            assert_eq!(attempt.decision, FixDecision::Accepted);
        }
        if let Some(solution) = attempt.fixed_solution {
            let expected: Vec<_> = float
                .ambiguity_satellites
                .iter()
                .map(|s| f.sd[&float.reference_satellite] - f.sd[s])
                .collect();
            assert_eq!(attempt.integers, expected);
            assert!(distance(solution.position_ecef.unwrap(), f.truth) < 0.003);
        }
    }
}

#[test]
fn power_gap_event_and_failed_update_require_fresh_confirmation() {
    for event in 0..4 {
        let f = fixture(false);
        let mut engine = filter(&f);
        let mut policy = RtkFixPolicy::new(RtkFixConfig::default()).unwrap();
        for i in 0..6 {
            let float = engine
                .process_pair(&f.rover[i], &f.base[i], &f.nav)
                .unwrap();
            let attempt = policy.evaluate(&float).unwrap();
            if i == 5 {
                assert!(attempt.fixed_solution.is_some());
            }
        }
        let start = match event {
            1 => 10, // A 150-second gap also resets the FLOAT ambiguity arcs.
            2 => {
                engine.invalidate_ambiguities();
                policy.reset();
                6
            }
            3 => {
                let previous = engine.state().unwrap().clone();
                let mut rover = f.rover[6].clone();
                *rover.observations[0].carrier_phase_cycles.as_mut().unwrap() += 50.0;
                assert!(engine.process_pair(&rover, &f.base[6], &f.nav).is_err());
                assert_eq!(*engine.state().unwrap(), previous);
                7
            }
            _ => 6,
        };
        for i in start..start + 6 {
            let mut rover = f.rover[i].clone();
            if event == 0 && i == start {
                rover.flag = 1;
            }
            let float = engine.process_pair(&rover, &f.base[i], &f.nav).unwrap();
            let attempt = policy.evaluate(&float).unwrap();
            assert_eq!(
                attempt.decision,
                if i < start + 4 {
                    FixDecision::InsufficientLock
                } else if i == start + 4 {
                    FixDecision::Confirming
                } else {
                    FixDecision::Accepted
                },
                "event {event}, epoch {i}"
            );
        }
    }
}

#[test]
fn a_shorter_policy_history_gap_also_restarts_lock_counts() {
    let f = fixture(false);
    let mut engine = filter(&f);
    let mut policy = RtkFixPolicy::new(RtkFixConfig {
        max_history_age_s: 45.0,
        ..Default::default()
    })
    .unwrap();
    for i in 0..6 {
        let float = engine
            .process_pair(&f.rover[i], &f.base[i], &f.nav)
            .unwrap();
        policy.evaluate(&float).unwrap();
    }
    let float = engine
        .process_pair(&f.rover[7], &f.base[7], &f.nav)
        .unwrap();
    assert!(float.reset_ambiguities.is_empty());
    let attempt = policy.evaluate(&float).unwrap();
    assert_eq!(attempt.decision, FixDecision::InsufficientLock);
    assert!(attempt.fixed_solution.is_none());
}

#[test]
fn exact_zero_residual_uses_upstream_zero_ratio_and_remains_float() {
    let f = fixture(false);
    let mut engine =
        RtkFloatFilter::new(f.base_position, f.base_position, RtkFloatConfig::default()).unwrap();
    let mut policy = RtkFixPolicy::new(RtkFixConfig::default()).unwrap();
    for i in 0..8 {
        let float = engine.process_pair(&f.base[i], &f.base[i], &f.nav).unwrap();
        let attempt = policy.evaluate(&float).unwrap();
        assert!(attempt.fixed_solution.is_none());
        if i >= 4 {
            assert_eq!(attempt.ratio, Some(0.0));
            assert_eq!(attempt.decision, FixDecision::RatioRejected);
            assert_eq!(attempt.integers, vec![0; 8]);
            let validation =
                validate_candidate(&float, &attempt.integers, RtkFixConfig::default()).unwrap();
            assert_eq!(validation.rejection, None);
            assert_eq!(validation.phase_residual_rms_m, 0.0);
        }
    }
}
