use gnss_rust::filter::{FilterError, GaussianState};
use gnss_rust::ppp_screening::*;

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 1e-10 * a.abs().max(b.abs()).max(1.0),
        "{a:.17e} != {b:.17e}"
    );
}
#[test]
fn standardized_residuals_match_128_actual_native_calls_and_128_gaussian_reconstructions() {
    for (data, production) in [
        (include_str!("fixtures/upstream_ppp_screening.csv"), false),
        (
            include_str!("fixtures/upstream_ppp_screening_gaussian.csv"),
            true,
        ),
    ] {
        for line in data.lines().skip(1) {
            let f: Vec<_> = line.split(',').collect();
            let m: usize = f[1].parse().unwrap();
            let gate = f[2].parse().unwrap();
            let worst: i32 = f[3].parse().unwrap();
            let x: Vec<f64> = f[4..].iter().map(|v| v.parse().unwrap()).collect();
            assert_eq!(x.len(), 2 * m * m + 2 * m);
            let r = &x[..m];
            let s: Vec<_> = x[m..m + m * m].chunks(m).map(|r| r.to_vec()).collect();
            let inverse: Vec<_> = x[m + m * m..m + 2 * m * m]
                .chunks(m)
                .map(|r| r.to_vec())
                .collect();
            let expected = &x[m + 2 * m * m..];
            let actual = screen_with_innovation_inverse(r, &inverse, gate)
                .unwrap_or_else(|e| panic!("native case {}: {e:?}", f[0]));
            assert_eq!(
                actual.worst_row,
                worst.try_into().ok(),
                "native case {}",
                f[0]
            );
            for (&a, &b) in actual.standardized_residuals.iter().zip(expected) {
                close(a, b);
            }
            // Independent covariance construction P=.01 I, H=I, R=S-P is SPD.
            let prior = GaussianState {
                mean: vec![0.0; m],
                covariance: (0..m)
                    .map(|i| (0..m).map(|j| if i == j { 0.01 } else { 0.0 }).collect())
                    .collect(),
            };
            let h = (0..m)
                .map(|i| (0..m).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
                .collect();
            let mut noise = s.clone();
            for (i, row) in noise.iter_mut().enumerate() {
                row[i] -= 0.01;
            }
            let gaussian = screen_measurements(&prior, &h, r, &noise, gate).unwrap();
            assert_eq!(gaussian.worst_row, actual.worst_row);
            for (&a, &b) in gaussian.standardized_residuals.iter().zip(expected) {
                assert!(
                    (a - b).abs() <= 1e-10 * a.abs().max(b.abs()).max(1.0),
                    "Gaussian production={production} case {}: {a:.17e} != {b:.17e}",
                    f[0]
                );
            }
        }
    }
}
#[test]
fn strict_threshold_first_tie_and_shared_prior_are_not_raw_sigma_screening() {
    let identity = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    assert_eq!(
        screen_with_innovation_inverse(&[4.0, -4.0], &identity, 4.0)
            .unwrap()
            .worst_row,
        None
    );
    assert_eq!(
        screen_with_innovation_inverse(&[6.0, -6.0], &identity, 4.0)
            .unwrap()
            .worst_row,
        Some(0)
    );
    let p = GaussianState {
        mean: vec![0.0],
        covariance: vec![vec![1e8]],
    };
    let result =
        screen_measurements(&p, &vec![vec![1.0]; 2], &[1000.0; 2], &identity, 4.0).unwrap();
    assert_eq!(result.worst_row, None);
    assert!(result.standardized_residuals.iter().all(|x| x.abs() < 0.01));
}
#[test]
fn invalid_screen_inputs_are_errors_not_clean_data_and_do_not_mutate_the_prior() {
    let p = GaussianState {
        mean: vec![0.0],
        covariance: vec![vec![1.0]],
    };
    let old = p.clone();
    for gate in [0.0, -4.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            screen_with_innovation_inverse(&[1.0], &vec![vec![1.0]], gate),
            Err(FilterError::InvalidGate)
        );
    }
    assert_eq!(
        screen_with_innovation_inverse(&[], &vec![], 4.0),
        Err(FilterError::InvalidDimensions)
    );
    assert_eq!(
        screen_with_innovation_inverse(&[1.0], &vec![vec![0.0]], 4.0),
        Err(FilterError::InvalidCovariance)
    );
    assert_eq!(
        screen_with_innovation_inverse(&[f64::NAN], &vec![vec![1.0]], 4.0),
        Err(FilterError::NonFinite)
    );
    assert_eq!(
        screen_with_innovation_inverse(&[1.0, 2.0], &vec![vec![1.0]], 4.0),
        Err(FilterError::InvalidDimensions)
    );
    assert!(screen_measurements(&p, &vec![vec![1.0]], &[1.0], &vec![vec![0.0]], 4.0).is_err());
    assert!(screen_measurements(&p, &vec![vec![f64::NAN]], &[1.0], &vec![vec![1.0]], 4.0).is_err());
    assert!(screen_measurements(&p, &vec![vec![]], &[1.0], &vec![vec![1.0]], 4.0).is_err());
    assert_eq!(p, old);
}
