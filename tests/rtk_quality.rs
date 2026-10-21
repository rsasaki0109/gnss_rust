use gnss_rust::Error;
use gnss_rust::filter::{GaussianState, measurement_update};
use gnss_rust::rtk_quality::screen_prefit_residuals;

#[test]
fn classifications_match_ten_pinned_native_row_suppression_cases() {
    let mut count = 0;
    for line in include_str!("fixtures/upstream_rtk_prefit.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split(',').collect();
        let n: usize = fields[1].parse().unwrap();
        let threshold: f64 = fields[2].parse().unwrap();
        let v: Vec<f64> = fields[3..].iter().map(|s| s.parse().unwrap()).collect();
        assert_eq!(v.len(), 4 * n + 1);
        let screen = screen_prefit_residuals(&v[..n], &v[n..2 * n], threshold).unwrap();
        for i in 0..n {
            assert_eq!(
                screen.rejected_phase_rows.contains(&i),
                v[2 * n + i] != 0.0,
                "{} phase {i}",
                fields[0]
            );
            assert_eq!(
                screen.rejected_code_rows.contains(&i),
                v[3 * n + i] != 0.0,
                "{} code {i}",
                fields[0]
            );
        }
        assert_eq!(
            screen.rejected_phase_rows.len() + screen.rejected_code_rows.len(),
            v[4 * n] as usize
        );
        assert_eq!(screen.accepted(), v[4 * n] == 0.0);
        count += 1;
    }
    assert_eq!(count, 10);
}

#[test]
fn invalid_shapes_values_and_limits_remain_explicit() {
    assert_eq!(
        screen_prefit_residuals(&[], &[], 3.0),
        Err(Error::InsufficientObservations)
    );
    assert_eq!(
        screen_prefit_residuals(&[0.0], &[0.0, 1.0], 3.0),
        Err(Error::InsufficientObservations)
    );
    for v in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            screen_prefit_residuals(&[v], &[0.0], 3.0),
            Err(Error::NonFinite)
        );
        assert_eq!(
            screen_prefit_residuals(&[0.0], &[v], 3.0),
            Err(Error::NonFinite)
        );
        assert_eq!(
            screen_prefit_residuals(&[0.0], &[0.0], v),
            Err(Error::NonFinite)
        );
    }
    for v in [0.0, -1.0] {
        assert_eq!(
            screen_prefit_residuals(&[0.0], &[0.0], v),
            Err(Error::OutOfRange)
        );
    }
}

#[test]
fn threshold_equality_passes_and_one_ulp_exceedance_is_retained() {
    let above = f64::from_bits(3.0_f64.to_bits() + 1);
    let result = screen_prefit_residuals(&[-3.0, 3.0], &[0.0, above], 3.0).unwrap();
    assert!(!result.accepted());
    assert!(result.rejected_phase_rows.is_empty());
    assert_eq!(result.rejected_code_rows, vec![1]);
    assert_eq!(result.maximum_phase_abs_m, 3.0);
    assert_eq!(result.maximum_code_abs_m, above);
}

#[test]
fn zeroing_a_correlated_row_and_marginalizing_it_are_different_updates() {
    let prior = GaussianState {
        mean: vec![0.0],
        covariance: vec![vec![1.0]],
    };
    let zeroed = measurement_update(
        &prior,
        &vec![vec![0.0], vec![1.0]],
        &[0.0, 1.0],
        &vec![vec![1.0, 0.5], vec![0.5, 1.0]],
        None,
    )
    .unwrap();
    let marginal =
        measurement_update(&prior, &vec![vec![1.0]], &[1.0], &vec![vec![1.0]], None).unwrap();
    assert!((zeroed.posterior.mean[0] - 4.0 / 7.0).abs() < 1e-14);
    assert!((zeroed.posterior.covariance[0][0] - 3.0 / 7.0).abs() < 1e-14);
    assert!((marginal.posterior.mean[0] - 0.5).abs() < 1e-14);
    assert!((marginal.posterior.covariance[0][0] - 0.5).abs() < 1e-14);
    assert!(
        !screen_prefit_residuals(&[4.0], &[1.0], 3.0)
            .unwrap()
            .accepted()
    );
    assert_eq!(prior.mean, vec![0.0]);
    assert_eq!(prior.covariance, vec![vec![1.0]]);
}
