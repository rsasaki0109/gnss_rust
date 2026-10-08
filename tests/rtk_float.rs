use gnss_rust::filter::{FilterError, GaussianState, Matrix, measurement_update};
use gnss_rust::lambda::LambdaConfig;
use gnss_rust::navigation::NavigationData;
use gnss_rust::observation::ObservationEpoch;
use gnss_rust::rinex::{RinexObservationReader, read_navigation};
use gnss_rust::rtk::{
    RtkError, RtkFloatConfig, RtkFloatFilter, RtkRejection, ambiguity_transform,
    double_difference_covariance,
};
use gnss_rust::{EcefCoord, Error, SatelliteId, SolutionStatus};
use std::io::Cursor;

fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * b.abs().max(1.0),
        "{a:.17} != {b:.17}"
    );
}
fn vector(v: &[f64], offset: &mut usize, n: usize) -> Vec<f64> {
    let out = v[*offset..*offset + n].to_vec();
    *offset += n;
    out
}
fn matrix(v: &[f64], offset: &mut usize, n: usize, m: usize) -> Matrix {
    vector(v, offset, n * m)
        .chunks(m)
        .map(|r| r.to_vec())
        .collect()
}

#[test]
fn filter_and_ambiguity_transform_match_native_cpp_components() {
    let mut cases = 0;
    for line in include_str!("fixtures/upstream_rtk_update.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let mut fields = line.split(',');
        let name = fields.next().unwrap();
        let n: usize = fields.next().unwrap().parse().unwrap();
        let m: usize = fields.next().unwrap().parse().unwrap();
        let v: Vec<f64> = fields.map(|s| s.parse().unwrap()).collect();
        let mut offset = 0;
        let prior = GaussianState {
            mean: vector(&v, &mut offset, n),
            covariance: matrix(&v, &mut offset, n, n),
        };
        let h = matrix(&v, &mut offset, m, n);
        let residual = vector(&v, &mut offset, m);
        let noise = matrix(&v, &mut offset, m, m);
        let expected_mean = vector(&v, &mut offset, n);
        let expected_cov = matrix(&v, &mut offset, n, n);
        let expected_nis = vector(&v, &mut offset, 1)[0];
        let expected_variances = vector(&v, &mut offset, m);
        let out = measurement_update(&prior, &h, &residual, &noise, None).unwrap();
        for i in 0..n {
            close(out.posterior.mean[i], expected_mean[i], 1e-10);
            for (actual, expected) in out.posterior.covariance[i].iter().zip(&expected_cov[i]) {
                close(*actual, *expected, 1e-10);
            }
        }
        close(out.normalized_innovation_squared, expected_nis, 1e-10);
        for (actual, expected) in out.innovation_variances.iter().zip(&expected_variances) {
            close(*actual, *expected, 1e-12);
        }
        let k = vector(&v, &mut offset, 1)[0] as usize;
        let pairs: Vec<_> = vector(&v, &mut offset, 2 * k)
            .chunks(2)
            .map(|p| (p[0] as usize, p[1] as usize))
            .collect();
        let expected_dd = vector(&v, &mut offset, k);
        let expected_q = matrix(&v, &mut offset, k, k);
        let expected_cross = matrix(&v, &mut offset, 3, k);
        let transform = ambiguity_transform(&out.posterior, 3, &pairs).unwrap();
        for i in 0..k {
            close(transform.dd_float_cycles[i], expected_dd[i], 1e-10);
            for (actual, expected) in transform.ambiguity_covariance_cycles2[i]
                .iter()
                .zip(&expected_q[i])
            {
                close(*actual, *expected, 1e-10);
            }
            for (actual_row, expected_row) in transform
                .head_ambiguity_covariance
                .iter()
                .zip(&expected_cross)
            {
                close(actual_row[i], expected_row[i], 1e-10);
            }
        }
        assert_eq!(offset, v.len(), "{name}");
        cases += 1;
    }
    assert_eq!(cases, 3);
}

#[test]
fn zero_state_scalar_update_matches_closed_form_and_remains_positive() {
    let prior = GaussianState {
        mean: vec![0.],
        covariance: vec![vec![2.]],
    };
    let out = measurement_update(&prior, &vec![vec![1.]], &[1.], &vec![vec![3.]], None).unwrap();
    close(out.posterior.mean[0], 0.4, 1e-14);
    close(out.posterior.covariance[0][0], 1.2, 1e-14);
    close(out.normalized_innovation_squared, 0.2, 1e-14);
    let mut state = GaussianState {
        mean: vec![0.],
        covariance: vec![vec![1000.]],
    };
    for _ in 0..500 {
        state = measurement_update(&state, &vec![vec![1.]], &[0.], &vec![vec![1e-8]], None)
            .unwrap()
            .posterior;
        assert!(state.covariance[0][0] > 0.0);
        assert_eq!(state.mean[0], 0.0);
    }
    close(
        state.covariance[0][0],
        1.0 / (1.0 / 1000.0 + 500.0 / 1e-8),
        1e-20,
    );
}

#[test]
fn covariance_keeps_shared_reference_correlation_and_checks_inputs() {
    assert_eq!(
        double_difference_covariance(0.1, &[1., 2., 3.]).unwrap(),
        vec![
            vec![1.1, 0.1, 0.1],
            vec![0.1, 2.1, 0.1],
            vec![0.1, 0.1, 3.1]
        ]
    );
    assert_eq!(
        double_difference_covariance(0.1, &[]),
        Err(FilterError::InvalidDimensions)
    );
    assert_eq!(
        double_difference_covariance(-0.1, &[1.]),
        Err(FilterError::InvalidCovariance)
    );
    assert_eq!(
        double_difference_covariance(0.1, &[0.]),
        Err(FilterError::InvalidCovariance)
    );
    assert_eq!(
        double_difference_covariance(0.1, &[f64::NAN]),
        Err(FilterError::NonFinite)
    );
}

#[test]
fn malformed_covariances_dimensions_gates_and_transforms_fail() {
    let prior = GaussianState {
        mean: vec![0.],
        covariance: vec![vec![2.]],
    };
    let h = vec![vec![1.]];
    let r = vec![vec![1.]];
    assert_eq!(
        measurement_update(&prior, &h, &[], &r, None),
        Err(FilterError::InvalidDimensions)
    );
    assert_eq!(
        measurement_update(&prior, &vec![vec![]], &[1.], &r, None),
        Err(FilterError::InvalidDimensions)
    );
    assert_eq!(
        measurement_update(&prior, &h, &[f64::NAN], &r, None),
        Err(FilterError::NonFinite)
    );
    assert_eq!(
        measurement_update(&prior, &h, &[1.], &vec![vec![-1.]], None),
        Err(FilterError::InvalidCovariance)
    );
    assert_eq!(
        measurement_update(&prior, &h, &[1.], &r, Some(0.)),
        Err(FilterError::InvalidGate)
    );
    assert_eq!(
        measurement_update(&prior, &h, &[100.], &r, Some(1.)),
        Err(FilterError::InnovationRejected)
    );
    assert_eq!(
        ambiguity_transform(&prior, 1, &[(0, 0)]),
        Err(FilterError::InvalidDimensions)
    );
    assert_eq!(
        ambiguity_transform(&prior, 1, &[(0, 2)]),
        Err(FilterError::InvalidDimensions)
    );
    assert_eq!(
        ambiguity_transform(&prior, 2, &[(0, 1)]),
        Err(FilterError::InvalidDimensions)
    );
    let invalid = GaussianState {
        mean: vec![0., 0.],
        covariance: vec![vec![1., 0.2], vec![0.1, 1.]],
    };
    assert_eq!(
        measurement_update(&invalid, &vec![vec![1., 1.]], &[0.], &r, None),
        Err(FilterError::InvalidCovariance)
    );
}

struct Fixture {
    nav: NavigationData,
    rover: Vec<ObservationEpoch>,
    base: Vec<ObservationEpoch>,
    truth: EcefCoord,
    base_position: EcefCoord,
    seed: EcefCoord,
    sd: std::collections::BTreeMap<SatelliteId, f64>,
}
fn fixture(noisy: bool) -> Fixture {
    let file = include_str!("fixtures/synthetic_rtk_truth.csv");
    let mut lines = file.lines().filter(|l| !l.starts_with('#'));
    let p: Vec<f64> = lines
        .next()
        .unwrap()
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect();
    let coord = |i| EcefCoord::new(p[i], p[i + 1], p[i + 2]);
    let sd = lines
        .map(|s| {
            let (sat, n) = s.split_once(',').unwrap();
            (sat.parse().unwrap(), n.parse().unwrap())
        })
        .collect();
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
        truth: coord(3),
        base_position: coord(0),
        seed: coord(6),
        sd,
    }
}
fn distance(a: EcefCoord, b: EcefCoord) -> f64 {
    (a.x - b.x).hypot(a.y - b.y).hypot(a.z - b.z)
}
fn processor(f: &Fixture) -> RtkFloatFilter {
    RtkFloatFilter::new(f.base_position, f.seed, RtkFloatConfig::default()).unwrap()
}

#[test]
fn rinex_static_float_pipeline_recovers_baseline_and_passes_dd_states_to_lambda() {
    let f = fixture(false);
    let mut filter = processor(&f);
    let mut max_error = 0.0_f64;
    assert_eq!(f.rover.len(), 32);
    assert_eq!(f.base.len(), 32);
    for (i, (rover, base)) in f.rover.iter().zip(&f.base).enumerate() {
        let result = filter.process_pair(rover, base, &f.nav).unwrap();
        assert_eq!(result.position.status, SolutionStatus::Float);
        assert_eq!(result.position.num_satellites, 9);
        assert!(!result.position.is_fixed());
        let error = distance(result.position.position_ecef.unwrap(), f.truth);
        max_error = max_error.max(error);
        assert!(error < 0.003, "epoch {i} error {error}");
        assert!(result.phase_residual_rms_m < 0.0001);
        assert_eq!(result.reset_ambiguities.len(), if i == 0 { 9 } else { 0 });
        let candidates = result.integer_candidates(LambdaConfig::default()).unwrap();
        let expected: Vec<_> = result
            .ambiguity_satellites
            .iter()
            .map(|s| (f.sd[&result.reference_satellite] - f.sd[s]) as i64)
            .collect();
        assert_eq!(candidates.candidates[0].ambiguities, expected);
    }
    eprintln!("noiseless maximum 3-D FLOAT error: {max_error} m");
}

#[test]
fn noisy_static_epochs_keep_float_status_and_finite_positive_covariance() {
    let f = fixture(true);
    let mut filter = processor(&f);
    let mut max_error = 0.0_f64;
    for (rover, base) in f.rover.iter().zip(&f.base) {
        let result = filter.process_pair(rover, base, &f.nav).unwrap();
        assert_eq!(result.position.status, SolutionStatus::Float);
        let error = distance(result.position.position_ecef.unwrap(), f.truth);
        max_error = max_error.max(error);
        assert!(error < 0.15, "{error}");
        for i in 0..3 {
            assert!(result.position.position_covariance.unwrap()[i][i] > 0.0);
        }
        assert!(
            result
                .ambiguities
                .dd_float_cycles
                .iter()
                .all(|v| v.is_finite())
        );
        let candidates = result.integer_candidates(LambdaConfig::default()).unwrap();
        assert_eq!(candidates.candidates.len(), 2);
        assert!(candidates.ratio().unwrap().is_finite());
    }
    eprintln!("noisy maximum 3-D FLOAT error: {max_error} m");
}

#[test]
fn lli_base_slip_half_cycles_missing_satellites_and_reappearance_restart_arcs() {
    let f = fixture(false);
    let mut filter = processor(&f);
    filter
        .process_pair(&f.rover[0], &f.base[0], &f.nav)
        .unwrap();
    let slipped = f.rover[1].observations[0].satellite;
    let mut rover = f.rover[1].clone();
    rover.observations[0].lli = Some(1);
    *rover.observations[0].carrier_phase_cycles.as_mut().unwrap() += 7.0;
    let result = filter.process_pair(&rover, &f.base[1], &f.nav).unwrap();
    assert_eq!(result.reset_ambiguities, vec![slipped]);
    let mut base = f.base[2].clone();
    base.observations[1].lli = Some(1);
    // Returning from the injected shifted phase also starts a new flagged arc.
    let mut rover = f.rover[2].clone();
    rover.observations[0].lli = Some(1);
    let result = filter.process_pair(&rover, &base, &f.nav).unwrap();
    assert_eq!(
        result.reset_ambiguities,
        vec![slipped, base.observations[1].satellite]
    );
    let excluded = f.rover[3].observations[2].satellite;
    let mut rover = f.rover[3].clone();
    rover.observations[2].lli = Some(2);
    let result = filter.process_pair(&rover, &f.base[3], &f.nav).unwrap();
    assert!(
        result
            .rejected
            .contains(&(excluded, RtkRejection::HalfCycle))
    );
    assert_eq!(result.position.num_satellites, 8);
    let result = filter
        .process_pair(&f.rover[4], &f.base[4], &f.nav)
        .unwrap();
    assert!(result.reset_ambiguities.contains(&excluded));
    let absent = f.rover[5].observations[3].satellite;
    let mut rover = f.rover[5].clone();
    rover.observations.remove(3);
    let result = filter.process_pair(&rover, &f.base[5], &f.nav).unwrap();
    assert!(
        result
            .rejected
            .contains(&(absent, RtkRejection::MissingPair))
    );
    let result = filter
        .process_pair(&f.rover[6], &f.base[6], &f.nav)
        .unwrap();
    assert_eq!(result.reset_ambiguities, vec![absent]);
}

#[test]
fn power_failure_gap_and_external_events_reset_all_ambiguities() {
    let f = fixture(false);
    let mut filter = processor(&f);
    filter
        .process_pair(&f.rover[0], &f.base[0], &f.nav)
        .unwrap();
    let mut rover = f.rover[1].clone();
    rover.flag = 1;
    assert_eq!(
        filter
            .process_pair(&rover, &f.base[1], &f.nav)
            .unwrap()
            .reset_ambiguities
            .len(),
        9
    );
    assert_eq!(
        filter
            .process_pair(&f.rover[8], &f.base[8], &f.nav)
            .unwrap()
            .reset_ambiguities
            .len(),
        9
    );
    filter.invalidate_ambiguities();
    assert_eq!(
        filter
            .process_pair(&f.rover[9], &f.base[9], &f.nav)
            .unwrap()
            .reset_ambiguities
            .len(),
        9
    );
}

#[test]
fn failed_or_out_of_order_epochs_do_not_commit_state_and_break_phase_continuity() {
    let f = fixture(false);
    let mut filter = processor(&f);
    filter
        .process_pair(&f.rover[0], &f.base[0], &f.nav)
        .unwrap();
    let previous = filter.state().unwrap().clone();
    let mut base = f.base[1].clone();
    base.time = base.time.checked_add_seconds(0.01).unwrap();
    assert_eq!(
        filter.process_pair(&f.rover[1], &base, &f.nav),
        Err(RtkError::UnsynchronizedEpochs)
    );
    assert_eq!(filter.state(), Some(&previous));
    assert_eq!(
        filter
            .process_pair(&f.rover[1], &f.base[1], &f.nav)
            .unwrap()
            .reset_ambiguities
            .len(),
        9
    );
    assert_eq!(
        filter.process_pair(&f.rover[1], &f.base[1], &f.nav),
        Err(RtkError::OutOfOrder)
    );
}

#[test]
fn reference_changes_preserve_sd_states_and_dd_candidate_signs() {
    let f = fixture(false);
    let mut filter = processor(&f);
    let first = filter
        .process_pair(&f.rover[0], &f.base[0], &f.nav)
        .unwrap();
    let mut rover = f.rover[1].clone();
    rover
        .observations
        .retain(|o| o.satellite != first.reference_satellite);
    let second = filter.process_pair(&rover, &f.base[1], &f.nav).unwrap();
    assert_ne!(first.reference_satellite, second.reference_satellite);
    assert!(second.reset_ambiguities.is_empty());
    assert_eq!(second.position.num_satellites, 8);
    let integers = second.integer_candidates(LambdaConfig::default()).unwrap();
    let expected: Vec<_> = second
        .ambiguity_satellites
        .iter()
        .map(|s| (f.sd[&second.reference_satellite] - f.sd[s]) as i64)
        .collect();
    assert_eq!(integers.candidates[0].ambiguities, expected);
    let third = filter
        .process_pair(&f.rover[2], &f.base[2], &f.nav)
        .unwrap();
    assert_eq!(third.reset_ambiguities, vec![first.reference_satellite]);
    assert_eq!(third.reference_satellite, first.reference_satellite);
}

#[test]
fn large_phase_jump_rejected_by_innovation_gate_does_not_commit_and_restarts_arcs() {
    let f = fixture(false);
    let mut filter = processor(&f);
    filter
        .process_pair(&f.rover[0], &f.base[0], &f.nav)
        .unwrap();
    let previous = filter.state().unwrap().clone();
    let mut rover = f.rover[1].clone();
    *rover.observations[0].carrier_phase_cycles.as_mut().unwrap() += 50.0;
    assert_eq!(
        filter.process_pair(&rover, &f.base[1], &f.nav),
        Err(RtkError::Filter(FilterError::InnovationRejected))
    );
    assert_eq!(filter.state(), Some(&previous));
    assert_eq!(
        filter
            .process_pair(&f.rover[2], &f.base[2], &f.nav)
            .unwrap()
            .reset_ambiguities
            .len(),
        9
    );
}

#[test]
fn insufficient_bad_tracking_clock_metadata_duplicate_and_stale_navigation_fail() {
    let f = fixture(false);
    let mut rover = f.rover[0].clone();
    rover.observations.truncate(3);
    assert_eq!(
        processor(&f).process_pair(&rover, &f.base[0], &f.nav),
        Err(RtkError::InsufficientCommonSatellites)
    );
    let mut base = f.base[0].clone();
    for obs in &mut base.observations {
        obs.tracking_code = "1W".into();
    }
    assert_eq!(
        processor(&f).process_pair(&f.rover[0], &base, &f.nav),
        Err(RtkError::InsufficientCommonSatellites)
    );
    let mut rover = f.rover[0].clone();
    rover.receiver_clock_offset_s = Some(0.001);
    assert_eq!(
        processor(&f).process_pair(&rover, &f.base[0], &f.nav),
        Err(RtkError::InvalidEpoch)
    );
    rover.receiver_clock_offset_s = None;
    rover.flag = 6;
    assert_eq!(
        processor(&f).process_pair(&rover, &f.base[0], &f.nav),
        Err(RtkError::InvalidEpoch)
    );
    rover.flag = 0;
    rover.observations.push(rover.observations[0].clone());
    assert_eq!(
        processor(&f).process_pair(&rover, &f.base[0], &f.nav),
        Err(RtkError::DuplicateObservation)
    );
    assert_eq!(
        processor(&f).process_pair(&f.rover[0], &f.base[0], &NavigationData::default()),
        Err(RtkError::InsufficientCommonSatellites)
    );
    let mut rover = f.rover[0].clone();
    let mut base = f.base[0].clone();
    rover.time = rover.time.checked_add_seconds(604800.).unwrap();
    base.time = rover.time;
    assert_eq!(
        processor(&f).process_pair(&rover, &base, &f.nav),
        Err(RtkError::InsufficientCommonSatellites)
    );
}

#[test]
fn zero_baseline_and_zero_ambiguities_remain_active() {
    let f = fixture(false);
    let mut filter =
        RtkFloatFilter::new(f.base_position, f.base_position, RtkFloatConfig::default()).unwrap();
    for base in &f.base[..3] {
        let result = filter.process_pair(base, base, &f.nav).unwrap();
        assert_eq!(result.position.position_ecef, Some(f.base_position));
        assert!(result.ambiguities.dd_float_cycles.iter().all(|v| *v == 0.0));
        assert_eq!(result.position.status, SolutionStatus::Float);
        assert_eq!(
            result
                .integer_candidates(LambdaConfig::default())
                .unwrap()
                .ratio(),
            Some(0.0)
        );
    }
}

#[test]
fn invalid_initialization_configuration_nonconvergence_and_degenerate_geometry_fail() {
    let f = fixture(false);
    assert!(matches!(
        RtkFloatFilter::new(EcefCoord::default(), f.seed, RtkFloatConfig::default()),
        Err(RtkError::Gnss(Error::UndefinedGeodeticOrigin))
    ));
    assert!(matches!(
        RtkFloatFilter::new(
            f.base_position,
            EcefCoord::new(
                f.base_position.x + 10001.,
                f.base_position.y,
                f.base_position.z
            ),
            RtkFloatConfig::default()
        ),
        Err(RtkError::BaselineTooLong)
    ));
    assert!(matches!(
        RtkFloatFilter::new(
            f.base_position,
            f.seed,
            RtkFloatConfig {
                carrier_phase_sigma_m: 0.,
                ..RtkFloatConfig::default()
            }
        ),
        Err(RtkError::InvalidConfiguration)
    ));
    let mut filter = RtkFloatFilter::new(
        f.base_position,
        f.seed,
        RtkFloatConfig {
            max_iterations: 1,
            ..RtkFloatConfig::default()
        },
    )
    .unwrap();
    assert_eq!(
        filter.process_pair(&f.rover[0], &f.base[0], &f.nav),
        Err(RtkError::NonConvergent)
    );
    assert_eq!(filter.state(), None);
    let ephemeris = f.nav.records().next().unwrap().clone();
    let mut nav = NavigationData::default();
    for obs in &f.rover[0].observations {
        let mut eph = ephemeris.clone();
        eph.satellite = obs.satellite;
        nav.add_ephemeris(eph).unwrap();
    }
    let config = RtkFloatConfig {
        elevation_mask_rad: 0.0,
        ..RtkFloatConfig::default()
    };
    let mut filter = RtkFloatFilter::new(f.base_position, f.seed, config).unwrap();
    assert!(matches!(
        filter.process_pair(&f.rover[0], &f.base[0], &nav),
        Err(RtkError::Gnss(Error::SingularGeometry)) | Err(RtkError::InsufficientCommonSatellites)
    ));
}
