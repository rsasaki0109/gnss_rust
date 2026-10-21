use gnss_rust::rinex::{RinexObservationReader, read_navigation};
use gnss_rust::spp::SppConfig;
use gnss_rust::spp_multiclock::{
    solve_gps_galileo_epoch, solve_gps_galileo_epoch_with_variance_model,
};
use gnss_rust::spp_quality::{
    MeasurementVarianceInputs, SppVarianceModel, measurement_variance, residual_inlier_mask,
    select_raim_fde_candidate,
};
use gnss_rust::{EcefCoord, Error};
use std::io::Cursor;

#[test]
fn noise_mask_and_fde_selectors_match_all_80_actual_native_cases() {
    let mut counts = [0; 3];
    for line in include_str!("fixtures/upstream_spp_quality.csv")
        .lines()
        .skip(1)
    {
        let mut fields = line.split(',');
        let tag = fields.next().unwrap();
        let mut v = fields.map(|s| s.parse::<f64>().unwrap());
        let id = v.next().unwrap();
        match tag {
            "V" => {
                let elevation = v.next().unwrap();
                let snr = v.next().unwrap();
                let sigma = v.next().unwrap();
                let iono = v.next().unwrap();
                let trop = v.next().unwrap();
                let ic = v.next().unwrap() != 0.0;
                let tc = v.next().unwrap() != 0.0;
                let reference = v.next().unwrap();
                let iscale = v.next().unwrap();
                let tscale = v.next().unwrap();
                let expected = v.next().unwrap();
                let actual = measurement_variance(MeasurementVarianceInputs {
                    elevation_rad: elevation,
                    snr_db_hz: Some(snr),
                    pseudorange_sigma_m: sigma,
                    ionosphere_delay_m: iono,
                    troposphere_delay_m: trop,
                    ionosphere_corrected: ic,
                    troposphere_corrected: tc,
                    snr_reference_db_hz: reference,
                    ionosphere_sigma_scale: iscale,
                    troposphere_sigma_scale: tscale,
                })
                .unwrap();
                assert!(
                    (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
                    "case {id}: {actual} vs {expected}"
                );
                counts[0] += 1;
            }
            "M" => {
                let n = v.next().unwrap() as usize;
                let threshold = v.next().unwrap();
                let floor = v.next().unwrap();
                let residuals: Vec<_> = v.by_ref().take(n).collect();
                let mask: Vec<_> = v.by_ref().take(n).map(|x| x != 0.0).collect();
                assert_eq!(
                    residual_inlier_mask(&residuals, threshold, floor),
                    mask,
                    "case {id}"
                );
                counts[1] += 1;
            }
            "F" => {
                let baseline = v.next().unwrap();
                let ratio = v.next().unwrap();
                let minimum = v.next().unwrap();
                let n = v.next().unwrap() as usize;
                let rms: Vec<_> = v.by_ref().take(n).collect();
                let expected = v.next().unwrap() as i32;
                assert_eq!(
                    select_raim_fde_candidate(baseline, &rms, ratio, minimum)
                        .map(|i| i as i32)
                        .unwrap_or(-1),
                    expected,
                    "case {id}"
                );
                counts[2] += 1;
            }
            _ => panic!("unexpected tag"),
        }
        assert!(v.next().is_none());
    }
    assert_eq!(counts, [32, 24, 24]);
}

#[test]
fn noise_clamps_missing_cn0_and_invalid_active_corrections_are_explicit() {
    let p = MeasurementVarianceInputs {
        elevation_rad: std::f64::consts::FRAC_PI_2,
        ionosphere_corrected: true,
        troposphere_corrected: true,
        ..Default::default()
    };
    let normal = measurement_variance(p).unwrap();
    assert!((normal - 9.065).abs() < 1e-12);
    for snr in [Some(0.0), Some(f64::NAN), Some(45.0), Some(50.0), None] {
        assert_eq!(
            measurement_variance(MeasurementVarianceInputs {
                snr_db_hz: snr,
                ..p
            })
            .unwrap(),
            normal
        );
    }
    let weak = measurement_variance(MeasurementVarianceInputs {
        snr_db_hz: Some(10.0),
        ..p
    })
    .unwrap();
    assert!((weak - 144.065).abs() < 1e-12);
    assert_eq!(
        measurement_variance(MeasurementVarianceInputs {
            ionosphere_delay_m: f64::NAN,
            ..p
        }),
        Err(Error::NonFinite)
    );
    assert_eq!(
        measurement_variance(MeasurementVarianceInputs {
            pseudorange_sigma_m: f64::MAX,
            ..p
        }),
        Err(Error::NonFinite)
    );
    for sigma in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            measurement_variance(MeasurementVarianceInputs {
                pseudorange_sigma_m: sigma,
                ..p
            }),
            Err(Error::InvalidConfiguration)
        );
    }
}

#[test]
fn robust_masks_and_fde_gates_keep_exact_boundaries_and_first_ties() {
    assert_eq!(
        residual_inlier_mask(&[0.0, 0.0, 0.0, 0.0, 9.0], 3.0, 3.0),
        vec![true; 5]
    );
    assert_eq!(
        residual_inlier_mask(
            &[0.0, 0.0, 0.0, 0.0, f64::from_bits(9.0_f64.to_bits() + 1)],
            3.0,
            3.0
        ),
        vec![true, true, true, true, false]
    );
    assert_eq!(
        residual_inlier_mask(&[0.0, 0.0, 0.0, 2.0], 3.0, 3.0),
        vec![true; 4]
    );
    assert_eq!(
        select_raim_fde_candidate(10.0, &[7.5, 7.5], 0.25, 1.0),
        Some(0)
    );
    assert_eq!(
        select_raim_fde_candidate(10.0, &[f64::from_bits(7.5_f64.to_bits() + 1)], 0.25, 1.0),
        None
    );
    assert_eq!(select_raim_fde_candidate(2.0, &[1.0], 0.25, 1.0), Some(0));
    assert_eq!(
        select_raim_fde_candidate(10.0, &[-1.0, f64::NAN, f64::INFINITY], 0.25, 1.0),
        None
    );
}

fn data() -> (
    Vec<gnss_rust::observation::ObservationEpoch>,
    gnss_rust::navigation::NavigationData,
) {
    (
        RinexObservationReader::new(Cursor::new(include_str!(
            "fixtures/synthetic_rtk_mixed_rover.obs"
        )))
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap(),
        read_navigation(Cursor::new(include_str!(
            "fixtures/synthetic_rtk_mixed.nav"
        )))
        .unwrap(),
    )
}
fn truth(i: usize) -> EcefCoord {
    let v: Vec<f64> = include_str!("fixtures/synthetic_rtk_mixed_truth.csv")
        .lines()
        .nth(i + 1)
        .unwrap()
        .split(',')
        .map(|v| v.parse().unwrap())
        .collect();
    EcefCoord::new(v[2], v[3], v[4])
}
fn distance(a: EcefCoord, b: EcefCoord) -> f64 {
    (a.x - b.x).hypot(a.y - b.y).hypot(a.z - b.z)
}
fn config() -> SppConfig {
    SppConfig {
        use_ionosphere: false,
        ..Default::default()
    }
}

#[test]
fn native_variance_preserves_all_32_clean_positions_clocks_and_raw_rows() {
    let (epochs, nav) = data();
    for (i, e) in epochs.iter().enumerate() {
        let legacy = solve_gps_galileo_epoch(e, &nav, config()).unwrap();
        let before = e.clone();
        let s = solve_gps_galileo_epoch_with_variance_model(
            e,
            &nav,
            config(),
            SppVarianceModel::NativeSnrAtmosphere,
        )
        .unwrap();
        assert!(distance(s.spp.position.position_ecef.unwrap(), truth(i)) < 0.003);
        assert_eq!(s.variance_model, SppVarianceModel::NativeSnrAtmosphere);
        assert!(
            (s.clock_biases_s[0].1 - 2e-6).abs() < 1e-11
                && (s.clock_biases_s[1].1 + 3e-6).abs() < 1e-11
        );
        for j in 0..5 {
            assert!(s.full_covariance_m2[j][j] >= legacy.full_covariance_m2[j][j]);
        }
        assert!(
            s.broadcast_models
                .iter()
                .all(|m| m.variance_m2.is_finite() && m.variance_m2 > 0.0)
        );
        assert_eq!(*e, before);
    }
}

#[test]
fn a_weak_cn0_code_outlier_is_downweighted_without_dropping_or_modifying_it() {
    let (epochs, nav) = data();
    let mut epoch = epochs[0].clone();
    for o in &mut epoch.observations {
        o.snr_db_hz = Some(45.0);
    }
    epoch.observations[0].pseudorange_m = epoch.observations[0].pseudorange_m.map(|v| v + 100.0);
    epoch.observations[0].snr_db_hz = Some(10.0);
    let before = epoch.clone();
    let a = solve_gps_galileo_epoch(&epoch, &nav, config()).unwrap();
    let b = solve_gps_galileo_epoch_with_variance_model(
        &epoch,
        &nav,
        config(),
        SppVarianceModel::NativeSnrAtmosphere,
    )
    .unwrap();
    assert!(
        distance(b.spp.position.position_ecef.unwrap(), truth(0))
            < distance(a.spp.position.position_ecef.unwrap(), truth(0))
    );
    assert_eq!(a.spp.position.num_satellites, b.spp.position.num_satellites);
    assert_eq!(epoch, before);
    assert_eq!(b.observations[0].raw_observation, epoch.observations[0]);
}
