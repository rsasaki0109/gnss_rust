use gnss_rust::EcefCoord;
use gnss_rust::filter::{FilterError, GaussianState};
use gnss_rust::ppp_motion::*;

fn data(line: &str) -> Vec<f64> {
    line.split(',').map(|s| s.parse().unwrap()).collect()
}
fn state(d: &[f64], start: usize, n: usize) -> GaussianState {
    GaussianState {
        mean: d[start..start + n].to_vec(),
        covariance: d[start + n..start + n + n * n]
            .chunks(n)
            .map(|r| r.to_vec())
            .collect(),
    }
}
fn close(actual: &GaussianState, expected: &GaussianState) {
    for (a, b) in actual.mean.iter().zip(&expected.mean).chain(
        actual
            .covariance
            .iter()
            .flatten()
            .zip(expected.covariance.iter().flatten()),
    ) {
        assert!((a - b).abs() <= 1e-10 * b.abs().max(1.0), "{a} versus {b}");
    }
}
fn identity(n: usize) -> GaussianState {
    GaussianState {
        mean: vec![0.0; n],
        covariance: (0..n)
            .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect(),
    }
}

#[test]
fn native_pva_seed_restart_matches_all_32_full_cross_covariance_cases() {
    let mut count = 0;
    let mut resets = 0;
    for line in include_str!("fixtures/upstream_ppp_motion_reset.csv")
        .lines()
        .skip(1)
    {
        let d = data(line);
        let n = d[1] as usize;
        assert_eq!(d.len(), 22 + 2 * n + 2 * n * n);
        let prior = state(&d, 22, n);
        let expected = state(&d, 22 + n + n * n, n);
        let seed = (d[4] == 1.0).then(|| MotionSeed {
            position_ecef: EcefCoord::new(d[7], d[8], d[9]),
            position_covariance_m2: (d[6] == 1.0)
                .then(|| std::array::from_fn(|i| std::array::from_fn(|j| d[13 + 3 * i + j]))),
            velocity_mps: (d[5] == 1.0).then(|| [d[10], d[11], d[12]]),
        });
        let r = recover_pva_motion(
            &prior,
            seed,
            PvaRecoveryConfig {
                position_variance_m2: d[2],
                position_disagreement_gate: (d[3] > 0.0).then_some(d[3]),
                ..Default::default()
            },
        )
        .unwrap();
        close(&r.prior, &expected);
        assert_eq!(r.reset, expected.mean[8] == 0.0);
        resets += usize::from(r.reset);
        count += 1;
    }
    assert_eq!(count, 32);
    assert!(resets > 0 && resets < 32);
}

#[test]
fn native_velocity_update_matches_32_dense_cases_and_rejects_the_same_outliers() {
    let mut count = 0;
    let mut rejected = 0;
    for line in include_str!("fixtures/upstream_ppp_motion_velocity.csv")
        .lines()
        .skip(1)
    {
        let d = data(line);
        let n = d[1] as usize;
        assert_eq!(d.len(), 16 + 2 * n + 2 * n * n);
        let prior = state(&d, 15, n);
        let expected = state(&d, 16 + n + n * n, n);
        let result = update_velocity(
            &prior,
            5,
            [d[3], d[4], d[5]],
            std::array::from_fn(|i| std::array::from_fn(|j| d[6 + 3 * i + j])),
            (d[2] > 0.0).then_some(d[2]),
        );
        if d[15 + n + n * n] == 1.0 {
            close(&result.unwrap().posterior, &expected);
        } else {
            assert_eq!(result, Err(FilterError::InnovationRejected));
            close(&prior, &expected);
            rejected += 1;
        }
        count += 1;
    }
    assert_eq!(count, 32);
    assert!(rejected > 0 && rejected < 32);
}

#[test]
fn native_seed_disagreement_matches_32_spd_and_failure_cases() {
    let mut count = 0;
    for line in include_str!("fixtures/upstream_ppp_motion_disagreement.csv")
        .lines()
        .skip(1)
    {
        let d = data(line);
        assert_eq!(d.len(), 26);
        let result = position_disagreement(
            EcefCoord::new(d[1], d[2], d[3]),
            std::array::from_fn(|i| std::array::from_fn(|j| d[7 + 3 * i + j])),
            EcefCoord::new(d[4], d[5], d[6]),
            std::array::from_fn(|i| std::array::from_fn(|j| d[16 + 3 * i + j])),
        );
        if d[25].is_infinite() {
            assert!(result.is_err());
        } else {
            let v = result.unwrap();
            assert!((v - d[25]).abs() < 1e-10 * d[25].abs().max(1.0));
        }
        count += 1;
    }
    assert_eq!(count, 32);
}

#[test]
fn strict_restart_boundary_preserves_unreset_marginal_and_missing_velocity_is_not_fused() {
    let mut p = identity(13);
    for i in 0..3 {
        p.covariance[i][i] = 3600.0;
    }
    p.mean[3] = 10.0;
    p.mean[4] = 2.4;
    p.covariance[3][4] = 0.2;
    p.covariance[4][3] = 0.2;
    let seed = MotionSeed {
        position_ecef: EcefCoord::new(1.0, 2.0, 3.0),
        position_covariance_m2: None,
        velocity_mps: None,
    };
    let c = PvaRecoveryConfig {
        position_disagreement_gate: None,
        ..Default::default()
    };
    assert!(!recover_pva_motion(&p, Some(seed), c).unwrap().reset);
    for i in 0..3 {
        p.covariance[i][i] = f64::from_bits(3600.0_f64.to_bits() + 1);
    }
    let r = recover_pva_motion(&p, Some(seed), c).unwrap();
    assert!(r.reset && r.variance_trigger);
    assert_eq!(&r.prior.mean[..3], &[1.0, 2.0, 3.0]);
    assert_eq!(&r.prior.mean[5..11], &[0.0; 6]);
    for i in [3, 4, 11, 12] {
        assert_eq!(r.prior.mean[i], p.mean[i]);
        for j in [3, 4, 11, 12] {
            assert_eq!(r.prior.covariance[i][j], p.covariance[i][j]);
        }
    }
    assert_eq!(recover_pva_motion(&p, None, c).unwrap().prior, p);
}

#[test]
fn invalid_motion_inputs_and_unphysical_velocity_covariance_fail_explicitly() {
    let p = identity(13);
    let zero = [[0.0; 3]; 3];
    let one = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    assert_eq!(
        update_velocity(&p, 5, [0.0; 3], zero, None),
        Err(FilterError::InvalidCovariance)
    );
    assert_eq!(
        update_velocity(&p, usize::MAX, [0.0; 3], one, None),
        Err(FilterError::InvalidDimensions)
    );
    assert_eq!(
        update_velocity(&p, 5, [f64::NAN, 0.0, 0.0], one, None),
        Err(FilterError::NonFinite)
    );
    assert_eq!(
        update_velocity(&p, 5, [0.0; 3], one, Some(0.0)),
        Err(FilterError::InvalidGate)
    );
    assert!(recover_pva_motion(&identity(8), None, Default::default()).is_err());
    let mut asymmetric = one;
    asymmetric[0][1] = 0.1;
    assert!(
        position_disagreement(
            EcefCoord::new(0.0, 0.0, 0.0),
            one,
            EcefCoord::new(1.0, 2.0, 3.0),
            asymmetric
        )
        .is_err()
    );
}
