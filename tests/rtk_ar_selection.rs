use gnss_rust::rtk_ar_selection::{
    RelativeVarianceConfig, SelectionError, filter_pairs_by_relative_variance,
    progressive_variance_drop_subsets,
};

#[test]
fn relative_filter_and_progressive_paths_match_all_40_native_cases() {
    let mut cases = 0;
    for line in include_str!("fixtures/upstream_rtk_ar_selection.csv")
        .lines()
        .skip(1)
    {
        let mut fields = line.split(',').map(|s| s.parse::<f64>().unwrap());
        let id = fields.next().unwrap();
        let n = fields.next().unwrap() as usize;
        let minimum = fields.next().unwrap() as usize;
        let steps = fields.next().unwrap() as usize;
        let multiplier = fields.next().unwrap();
        let floor = fields.next().unwrap();
        let variances: Vec<_> = fields.by_ref().take(n).collect();
        let selected_count = fields.next().unwrap() as usize;
        let expected: Vec<_> = fields
            .by_ref()
            .take(selected_count)
            .map(|v| v as usize)
            .collect();
        let count = fields.next().unwrap() as usize;
        let paths: Vec<Vec<usize>> = (0..count)
            .map(|_| {
                let n = fields.next().unwrap() as usize;
                fields.by_ref().take(n).map(|v| v as usize).collect()
            })
            .collect();
        let r = filter_pairs_by_relative_variance(
            &variances,
            RelativeVarianceConfig {
                multiplier,
                minimum_variance_cycles2: floor,
                minimum_pairs: minimum,
            },
        )
        .unwrap();
        assert_eq!(r.selected_indices, expected, "case {id}");
        assert_eq!(
            progressive_variance_drop_subsets(&variances, minimum, steps).unwrap(),
            paths,
            "case {id}"
        );
        assert!(fields.next().is_none());
        cases += 1;
    }
    assert_eq!(cases, 40);
}

#[test]
fn invalid_config_variances_and_overflow_fail_without_selecting_any_subset() {
    let defaults = RelativeVarianceConfig::default();
    for bad in [f64::NAN, f64::INFINITY, -1.0] {
        assert_eq!(
            filter_pairs_by_relative_variance(&[bad], defaults),
            Err(SelectionError::InvalidVariance)
        );
        assert_eq!(
            progressive_variance_drop_subsets(&[bad], 1, 6),
            Err(SelectionError::InvalidVariance)
        );
        assert_eq!(
            filter_pairs_by_relative_variance(
                &[1.0],
                RelativeVarianceConfig {
                    minimum_variance_cycles2: bad,
                    ..defaults
                }
            ),
            Err(SelectionError::InvalidConfiguration)
        );
    }
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            filter_pairs_by_relative_variance(
                &[1.0],
                RelativeVarianceConfig {
                    multiplier: bad,
                    ..defaults
                }
            ),
            Err(SelectionError::InvalidConfiguration)
        );
    }
    for n in [0, 129] {
        assert_eq!(
            filter_pairs_by_relative_variance(
                &[1.0],
                RelativeVarianceConfig {
                    minimum_pairs: n,
                    ..defaults
                }
            ),
            Err(SelectionError::InvalidConfiguration)
        );
        assert_eq!(
            progressive_variance_drop_subsets(&[1.0], n, 6),
            Err(SelectionError::InvalidConfiguration)
        );
    }
    assert_eq!(
        filter_pairs_by_relative_variance(&[f64::MAX; 4], defaults),
        Err(SelectionError::NumericalRange)
    );
    assert_eq!(
        progressive_variance_drop_subsets(&[1.0], 1, 129),
        Err(SelectionError::InvalidConfiguration)
    );
    assert_eq!(
        filter_pairs_by_relative_variance(&[1.0; 129], defaults),
        Err(SelectionError::InvalidVariance)
    );
}

#[test]
fn upper_median_exact_threshold_minimum_pair_fallback_and_tie_order_are_explicit() {
    let cfg = RelativeVarianceConfig::default();
    let r = filter_pairs_by_relative_variance(&[1.0, 2.0, 3.0, 29.0], cfg).unwrap();
    assert_eq!(r.median_variance_cycles2, Some(3.0));
    assert_eq!(r.threshold_cycles2, Some(30.0));
    assert_eq!(r.selected_indices, vec![0, 1, 2, 3]);
    let r = filter_pairs_by_relative_variance(&[0.0, 0.0, 0.0, 1.0], cfg).unwrap();
    assert_eq!(r.threshold_cycles2, Some(1e-4));
    assert!(r.selected_indices.is_empty());
    let r = filter_pairs_by_relative_variance(
        &[0.0, 0.0, 0.0, 1e-4, f64::from_bits(1e-4_f64.to_bits() + 1)],
        cfg,
    )
    .unwrap();
    assert_eq!(r.selected_indices, vec![0, 1, 2, 3]);
    assert_eq!(
        progressive_variance_drop_subsets(&[4.0, 4.0, 3.0, 2.0, 1.0, 1.0], 4, 6).unwrap(),
        vec![vec![1, 2, 3, 4, 5], vec![2, 3, 4, 5]]
    );
}
