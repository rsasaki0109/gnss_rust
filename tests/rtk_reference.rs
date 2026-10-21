use gnss_rust::SatelliteId;
use gnss_rust::filter::{FilterError, GaussianState, Matrix, measurement_update};
use gnss_rust::rtk::double_difference_covariance;
use gnss_rust::rtk_fix::RtkFixConfig;
use gnss_rust::rtk_reference::{DoubleDifferenceBasis, ReferenceError};
fn gps(prn: usize) -> SatelliteId {
    format!("G{prn:02}").parse().unwrap()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * b.abs().max(1.0),
        "{a:.17} != {b:.17}"
    );
}
fn vector(v: &[f64], offset: &mut usize, n: usize) -> Vec<f64> {
    let result = v[*offset..*offset + n].to_vec();
    *offset += n;
    result
}
fn matrix(v: &[f64], offset: &mut usize, rows: usize, columns: usize) -> Matrix {
    vector(v, offset, rows * columns)
        .chunks(columns)
        .map(<[f64]>::to_vec)
        .collect()
}
struct Case {
    id: usize,
    satellites: usize,
    reference: usize,
    prior: GaussianState,
    h: Matrix,
    v: Vec<f64>,
    r: Matrix,
    posterior: GaussianState,
    nis: f64,
    energy: f64,
    phase_rms: f64,
    code_rms: f64,
}
fn cases() -> Vec<Case> {
    include_str!("fixtures/upstream_rtk_reference_basis.csv")
        .lines()
        .filter(|s| !s.starts_with('#'))
        .map(|line| {
            let f: Vec<_> = line.split(',').collect();
            let id = f[0].parse().unwrap();
            let satellites = f[1].parse().unwrap();
            let n = f[2].parse().unwrap();
            let rows = f[3].parse().unwrap();
            let reference = f[4].parse().unwrap();
            let v: Vec<f64> = f[5..].iter().map(|s| s.parse().unwrap()).collect();
            let mut offset = 0;
            let prior = GaussianState {
                mean: vector(&v, &mut offset, n),
                covariance: matrix(&v, &mut offset, n, n),
            };
            let h = matrix(&v, &mut offset, rows, n);
            let residual = vector(&v, &mut offset, rows);
            let r = matrix(&v, &mut offset, rows, rows);
            let posterior = GaussianState {
                mean: vector(&v, &mut offset, n),
                covariance: matrix(&v, &mut offset, n, n),
            };
            let tail = vector(&v, &mut offset, 4);
            assert_eq!(offset, v.len());
            Case {
                id,
                satellites,
                reference,
                prior,
                h,
                v: residual,
                r,
                posterior,
                nis: tail[0],
                energy: tail[1],
                phase_rms: tail[2],
                code_rms: tail[3],
            }
        })
        .collect()
}

#[test]
fn all_22_native_reference_systems_match_design_noise_posterior_and_energy() {
    let cases = cases();
    assert_eq!(cases.len(), 22);
    for c in &cases {
        let first = cases
            .iter()
            .find(|x| x.id == c.id && x.reference == 1)
            .unwrap();
        let basis =
            DoubleDifferenceBasis::new((1..=c.satellites).map(gps).collect(), gps(1)).unwrap();
        let out = basis
            .rebase(gps(c.reference), 2, &first.h, &first.v, &first.r)
            .unwrap();
        assert_eq!(out.basis.reference(), gps(c.reference));
        for (a, b) in out.design.iter().flatten().zip(c.h.iter().flatten()) {
            close(*a, *b, 1e-12);
        }
        for (a, b) in out.residual.iter().zip(&c.v) {
            close(*a, *b, 1e-12);
        }
        for (a, b) in out.covariance.iter().flatten().zip(c.r.iter().flatten()) {
            close(*a, *b, 1e-12);
        }
        close(
            out.nominal_weighted_sum_of_squares().unwrap(),
            c.energy,
            1e-10,
        );
        let updated =
            measurement_update(&c.prior, &out.design, &out.residual, &out.covariance, None)
                .unwrap();
        close(updated.normalized_innovation_squared, c.nis, 1e-10);
        for i in 0..c.prior.mean.len() {
            close(updated.posterior.mean[i], c.posterior.mean[i], 1e-10);
            close(updated.posterior.mean[i], first.posterior.mean[i], 1e-10);
            for j in 0..c.prior.mean.len() {
                close(
                    updated.posterior.covariance[i][j],
                    c.posterior.covariance[i][j],
                    1e-10,
                );
                close(
                    updated.posterior.covariance[i][j],
                    first.posterior.covariance[i][j],
                    1e-10,
                );
            }
        }
    }
}

#[test]
fn same_information_has_reference_dependent_unweighted_phase_and_code_rms() {
    let cases = cases();
    for id in 0..4 {
        let group: Vec<_> = cases.iter().filter(|x| x.id == id).collect();
        let phase: Vec<_> = group.iter().map(|x| x.phase_rms).collect();
        let code: Vec<_> = group.iter().map(|x| x.code_rms).collect();
        let range = |v: &[f64]| {
            v.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                - v.iter().copied().fold(f64::INFINITY, f64::min)
        };
        assert!(range(&phase) > 0.001);
        assert!(range(&code) > 0.1);
        for c in group {
            close(c.nis, cases.iter().find(|x| x.id == id).unwrap().nis, 1e-10);
            close(
                c.energy,
                cases.iter().find(|x| x.id == id).unwrap().energy,
                1e-10,
            );
        }
    }
}

#[test]
fn phase_code_cross_covariance_is_material_and_preserved() {
    let c = cases().remove(0);
    let basis = DoubleDifferenceBasis::new((1..=c.satellites).map(gps).collect(), gps(1)).unwrap();
    let complete = basis.rebase(gps(2), 2, &c.h, &c.v, &c.r).unwrap();
    let k = c.satellites - 1;
    assert!(
        complete
            .covariance
            .iter()
            .take(k)
            .any(|r| r[k..].iter().any(|x| *x != 0.0))
    );
    let mut wrong = complete.clone();
    for i in 0..k {
        for j in k..2 * k {
            wrong.covariance[i][j] = 0.0;
            wrong.covariance[j][i] = 0.0;
        }
    }
    assert!(
        (complete.nominal_weighted_sum_of_squares().unwrap()
            - wrong.nominal_weighted_sum_of_squares().unwrap())
        .abs()
            > 1e-5
    );
}

#[test]
fn independent_dd_process_noise_is_basis_dependent_but_uniform_sd_noise_is_not() {
    let basis = DoubleDifferenceBasis::new(vec![gps(1), gps(2), gps(3)], gps(1)).unwrap();
    let h = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let v = vec![0.0, 0.0];
    let independent = basis.rebase(gps(2), 1, &h, &v, &h).unwrap();
    assert_eq!(independent.covariance, vec![vec![1.0, 1.0], vec![1.0, 2.0]]);
    assert_ne!(independent.covariance, h);
    let sd = double_difference_covariance(1.0, &[1.0, 1.0]).unwrap();
    assert_eq!(basis.rebase(gps(2), 1, &h, &v, &sd).unwrap().covariance, sd);
}

#[test]
fn existing_phase_rms_gate_can_cross_its_bound_under_equivalent_references() {
    let basis = DoubleDifferenceBasis::new(vec![gps(1), gps(2), gps(3)], gps(1)).unwrap();
    let h = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    // SD residuals [0, -.010, -.015] m, expressed as reference-minus-target.
    let v = vec![0.010, 0.015];
    let r = double_difference_covariance(4e-6, &[4e-6, 4e-6]).unwrap();
    let old = basis.rebase(gps(1), 1, &h, &v, &r).unwrap();
    let new = basis.rebase(gps(2), 1, &h, &v, &r).unwrap();
    let rms = |v: &[f64]| (v.iter().map(|x| x * x).sum::<f64>() / v.len() as f64).sqrt();
    let limit = RtkFixConfig::default().max_phase_residual_rms_m;
    assert!(rms(&old.residual) > limit);
    assert!(rms(&new.residual) < limit);
    close(
        old.nominal_weighted_sum_of_squares().unwrap(),
        new.nominal_weighted_sum_of_squares().unwrap(),
        1e-12,
    );
    // This proves one gate is basis dependent, not that all FIX gates pass.
}

#[test]
fn arbitrary_satellite_order_and_round_trip_keep_all_blocks() {
    let basis = DoubleDifferenceBasis::new(vec![gps(4), gps(1), gps(3), gps(2)], gps(4)).unwrap();
    let h: Matrix = (0..6).map(|i| vec![i as f64, -(i as f64)]).collect();
    let v: Vec<_> = (0..6).map(|i| i as f64 * 0.1).collect();
    let r: Matrix = (0..6)
        .map(|i| (0..6).map(|j| if i == j { 2.0 } else { 0.2 }).collect())
        .collect();
    let a = basis.rebase(gps(3), 2, &h, &v, &r).unwrap();
    assert_eq!(
        a.basis.targets().collect::<Vec<_>>(),
        vec![gps(4), gps(1), gps(2)]
    );
    let b = a
        .basis
        .rebase(gps(4), 2, &a.design, &a.residual, &a.covariance)
        .unwrap();
    assert_eq!(b.basis, basis);
    for (a, b) in b.design.iter().flatten().zip(h.iter().flatten()) {
        close(*a, *b, 1e-12);
    }
    for (a, b) in b.residual.iter().zip(v) {
        close(*a, b, 1e-12);
    }
    for (a, b) in b.covariance.iter().flatten().zip(r.iter().flatten()) {
        close(*a, *b, 1e-12);
    }
}

#[test]
fn invalid_identity_shape_and_covariance_do_not_mutate_inputs() {
    assert_eq!(
        DoubleDifferenceBasis::new(vec![gps(1)], gps(1)),
        Err(ReferenceError::InvalidBasis)
    );
    assert_eq!(
        DoubleDifferenceBasis::new(vec![gps(1), gps(1)], gps(1)),
        Err(ReferenceError::InvalidBasis)
    );
    assert_eq!(
        DoubleDifferenceBasis::new(vec![gps(1), "E01".parse().unwrap()], gps(1)),
        Err(ReferenceError::InvalidBasis)
    );
    assert_eq!(
        DoubleDifferenceBasis::new(vec![gps(1), gps(2)], gps(3)),
        Err(ReferenceError::ReferenceUnavailable)
    );
    let basis = DoubleDifferenceBasis::new(vec![gps(1), gps(2), gps(3)], gps(1)).unwrap();
    assert_eq!(
        basis.transformation(gps(2), 0),
        Err(ReferenceError::Filter(FilterError::InvalidDimensions))
    );
    assert_eq!(
        basis.transformation(gps(2), usize::MAX),
        Err(ReferenceError::Filter(FilterError::InvalidDimensions))
    );
    assert_eq!(
        basis.transformation(gps(4), 1),
        Err(ReferenceError::ReferenceUnavailable)
    );
    let h = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let v = vec![0.1, -0.1];
    let before = (h.clone(), v.clone());
    for r in [
        vec![vec![1.0]],
        vec![vec![0.0, 0.0], vec![0.0, 0.0]],
        vec![vec![1.0, 0.1], vec![0.0, 1.0]],
        vec![vec![f64::NAN, 0.0], vec![0.0, 1.0]],
    ] {
        assert!(basis.rebase(gps(2), 1, &h, &v, &r).is_err());
        assert_eq!((h.clone(), v.clone()), before);
    }
    assert!(basis.rebase(gps(2), 1, &Vec::new(), &v, &h).is_err());
    assert!(
        basis
            .rebase(gps(2), 1, &h, &[f64::INFINITY, 0.0], &h)
            .is_err()
    );
}
