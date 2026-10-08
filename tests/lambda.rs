use gnss_rust::lambda::{LambdaConfig, LambdaError, search};

fn config(candidate_count: usize) -> LambdaConfig {
    LambdaConfig {
        candidate_count,
        ..LambdaConfig::default()
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:.17} != {expected:.17}"
    );
}

#[test]
fn matches_actual_cpp_candidates_reduction_residuals_and_ratio() {
    let mut rows = 0;
    let mut names = std::collections::BTreeSet::new();
    for line in include_str!("fixtures/upstream_lambda.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split(',').collect();
        names.insert(fields[0]);
        let n: usize = fields[1].parse().unwrap();
        let m: usize = fields[2].parse().unwrap();
        let index: usize = fields[3].parse().unwrap();
        let values: Vec<f64> = fields[4..].iter().map(|s| s.parse().unwrap()).collect();
        let a = &values[..n];
        let q: Vec<_> = values[n..n + n * n].chunks(n).map(|r| r.to_vec()).collect();
        let out = search(a, &q, config(m)).unwrap();
        let mut offset = n + n * n;
        let expected = &values[offset..offset + n];
        offset += n;
        for (actual, expected) in out.candidates[index].ambiguities.iter().zip(expected) {
            // Upstream uses floating-point LU to back-transform integers.
            assert!((*expected - expected.round()).abs() < 1e-8);
            assert_eq!(*actual, expected.round() as i64);
        }
        close(
            out.candidates[index].squared_residual,
            values[offset],
            1e-12,
        );
        offset += 1;
        for i in 0..n {
            close(out.conditional_variances[i], values[offset + i], 1e-12);
        }
        offset += n;
        for i in 0..n {
            for j in 0..n {
                assert_eq!(
                    out.decorrelation_transform[i][j] as f64,
                    values[offset + i * n + j]
                );
            }
        }
        offset += n * n;
        for i in 0..n {
            close(out.decorrelated_float[i], values[offset + i], 1e-12);
        }
        offset += n;
        for i in 0..n {
            for j in 0..n {
                close(
                    out.decorrelated_covariance[i][j],
                    values[offset + i * n + j],
                    1e-12,
                );
            }
        }
        offset += n * n;
        if m >= 2 {
            close(out.ratio().unwrap(), values[offset], 1e-12);
        } else {
            assert_eq!(out.ratio(), None);
        }
        assert_eq!(offset + 1, values.len());
        rows += 1;
    }
    assert_eq!(names.len(), 11);
    assert_eq!(rows, 34);
}

// Independent exhaustive 2-D Mahalanobis scoring, with analytic matrix inverse.
fn exhaustive(a: &[f64; 2], q: &[Vec<f64>]) -> Vec<([i64; 2], f64)> {
    let det = q[0][0] * q[1][1] - q[0][1] * q[1][0];
    let mut candidates = Vec::new();
    for i in -12..=12 {
        for j in -12..=12 {
            let x = a[0] - i as f64;
            let y = a[1] - j as f64;
            let residual = (q[1][1] * x * x - 2.0 * q[0][1] * x * y + q[0][0] * y * y) / det;
            candidates.push(([i, j], residual));
        }
    }
    candidates.sort_by(|a, b| a.1.total_cmp(&b.1));
    candidates
}

#[test]
fn correlated_integer_lattice_matches_exhaustive_search() {
    let mut seed = 0x7839_u64;
    let mut random = || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        (seed >> 11) as f64 / (1_u64 << 53) as f64
    };
    for _ in 0..100 {
        let a = [6.0 * random() - 3.0, 6.0 * random() - 3.0];
        let u = 0.05 + random();
        let v = 0.05 + random();
        let rho = 1.9 * random() - 0.95;
        let q = vec![vec![u * u, rho * u * v], vec![rho * u * v, v * v]];
        let expected = exhaustive(&a, &q);
        let actual = search(&a, &q, config(4)).unwrap();
        for (actual, expected) in actual.candidates.iter().zip(expected.iter().take(4)) {
            assert_eq!(actual.ambiguities, expected.0);
            close(actual.squared_residual, expected.1, 1e-10);
        }
    }
    // Correlation means independent rounding can miss the best lattice point.
    let a = [0.49, -0.49];
    let q = vec![vec![1.0, 0.99], vec![0.99, 1.0]];
    let out = search(&a, &q, config(2)).unwrap();
    assert_ne!(out.candidates[0].ambiguities, vec![0, 0]);
    let expected = exhaustive(&a, &q);
    close(out.candidates[0].squared_residual, expected[0].1, 1e-10);
}

#[test]
fn one_dimensional_ties_and_zero_residual_follow_upstream_conventions() {
    let out = search(&[-1.5], &[vec![0.04]], config(4)).unwrap();
    assert_eq!(
        out.candidates
            .iter()
            .map(|c| c.ambiguities[0])
            .collect::<Vec<_>>(),
        vec![-2, -1, -3, 0]
    );
    assert_eq!(out.ratio(), Some(1.0));
    let exact = search(&[7.0], &[vec![0.04]], config(2)).unwrap();
    assert_eq!(exact.candidates[0].ambiguities, vec![7]);
    assert_eq!(exact.candidates[0].squared_residual, 0.0);
    assert_eq!(exact.ratio(), Some(0.0));
}

#[test]
fn translation_covariance_scale_and_permutation_preserve_integer_solutions() {
    let a = [3.16, -2.73];
    let q = vec![vec![0.13, -0.057], vec![-0.057, 0.04]];
    let base = search(&a, &q, config(3)).unwrap();
    let shifted = search(&[a[0] + 100., a[1] - 200.], &q, config(3)).unwrap();
    let scaled = search(
        &a,
        &q.iter()
            .map(|r| r.iter().map(|v| v * 5.0).collect())
            .collect::<Vec<_>>(),
        config(3),
    )
    .unwrap();
    let permuted = search(
        &[a[1], a[0]],
        &[vec![q[1][1], q[1][0]], vec![q[0][1], q[0][0]]],
        config(3),
    )
    .unwrap();
    for i in 0..3 {
        assert_eq!(
            shifted.candidates[i].ambiguities,
            vec![
                base.candidates[i].ambiguities[0] + 100,
                base.candidates[i].ambiguities[1] - 200
            ]
        );
        assert_eq!(
            scaled.candidates[i].ambiguities,
            base.candidates[i].ambiguities
        );
        assert_eq!(
            permuted.candidates[i].ambiguities,
            base.candidates[i]
                .ambiguities
                .iter()
                .copied()
                .rev()
                .collect::<Vec<_>>()
        );
        close(
            scaled.candidates[i].squared_residual * 5.0,
            base.candidates[i].squared_residual,
            1e-12,
        );
    }
    close(scaled.ratio().unwrap(), base.ratio().unwrap(), 1e-12);
}

#[test]
fn malformed_nonfinite_asymmetric_and_nonpositive_covariances_are_errors() {
    use LambdaError::*;
    assert_eq!(search(&[], &[], config(2)), Err(InvalidDimensions));
    assert_eq!(search(&[0.], &[], config(2)), Err(InvalidDimensions));
    assert_eq!(
        search(&[0., 0.], &[vec![1.], vec![1.]], config(2)),
        Err(InvalidDimensions)
    );
    assert_eq!(search(&[f64::NAN], &[vec![1.]], config(2)), Err(NonFinite));
    assert_eq!(
        search(&[0.], &[vec![f64::INFINITY]], config(2)),
        Err(NonFinite)
    );
    for q in [
        vec![vec![0.]],
        vec![vec![-1.]],
        vec![vec![1., 1.], vec![1., 1.]],
        vec![vec![1., 2.], vec![2., 1.]],
        vec![vec![1., 0.1], vec![0.2, 1.]],
    ] {
        assert_eq!(
            search(&vec![0.; q.len()], &q, config(2)),
            Err(InvalidCovariance)
        );
    }
    assert_eq!(
        search(&[0.], &[vec![1.]], config(0)),
        Err(InvalidConfiguration)
    );
    assert_eq!(
        search(&[0.], &[vec![1.]], config(65)),
        Err(InvalidConfiguration)
    );
    assert_eq!(
        search(
            &[0.],
            &[vec![1.]],
            LambdaConfig {
                max_search_iterations: 0,
                ..config(2)
            }
        ),
        Err(InvalidConfiguration)
    );
}

#[test]
fn exhausted_search_and_reduction_do_not_return_partial_candidates() {
    let a = [0.2, -0.3, 0.4];
    let q = vec![vec![1., 0., 0.], vec![0., 1., 0.], vec![0., 0., 1.]];
    assert_eq!(
        search(
            &a,
            &q,
            LambdaConfig {
                max_search_iterations: 1,
                ..config(2)
            }
        ),
        Err(LambdaError::IterationLimit)
    );
    assert_eq!(
        search(
            &a,
            &q,
            LambdaConfig {
                max_reduction_iterations: 1,
                ..config(2)
            }
        ),
        Err(LambdaError::IterationLimit)
    );
    assert!(search(&a, &q, config(2)).is_ok());
}

#[test]
fn unrepresentable_integer_or_arithmetic_range_is_rejected() {
    assert_eq!(
        search(&[2_f64.powi(52)], &[vec![1.]], config(2)),
        Err(LambdaError::NumericalRange)
    );
    assert_eq!(
        search(&[-2_f64.powi(52)], &[vec![1.]], config(2)),
        Err(LambdaError::NumericalRange)
    );
    assert_eq!(
        search(&[0.1], &[vec![f64::MIN_POSITIVE]], config(2)),
        Err(LambdaError::NumericalRange)
    );
    assert_eq!(
        search(&[0.1], &[vec![1e-320]], config(2)),
        Err(LambdaError::NumericalRange)
    );
}
