use gnss_rust::observation::ObservationEpoch;
use gnss_rust::rinex::{RinexObservationReader, read_navigation};
use gnss_rust::spp::SppConfig;
use gnss_rust::spp_multiclock::{
    ObservationDisposition, solve_gps_galileo_epoch_with_variance_model,
};
use gnss_rust::spp_qc::{SppQualityConfig, solve_gps_galileo_epoch_with_quality_control as solve};
use gnss_rust::spp_quality::SppVarianceModel;
use gnss_rust::{EcefCoord, Error, GnssSystem};
use std::io::Cursor;

fn data() -> (Vec<ObservationEpoch>, gnss_rust::navigation::NavigationData) {
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
fn config() -> SppConfig {
    SppConfig {
        use_ionosphere: false,
        ..Default::default()
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
fn distance(a: EcefCoord, b: EcefCoord) -> f64 {
    (a.x - b.x).hypot(a.y - b.y).hypot(a.z - b.z)
}
const VARIANCE: SppVarianceModel = SppVarianceModel::NativeSnrAtmosphere;

#[test]
fn clean_physical_epochs_keep_exact_baseline_and_test_every_raim_candidate() {
    let (epochs, nav) = data();
    for (i, e) in epochs.iter().enumerate() {
        let baseline =
            solve_gps_galileo_epoch_with_variance_model(e, &nav, config(), VARIANCE).unwrap();
        let result = solve(e, &nav, config(), VARIANCE, SppQualityConfig::default()).unwrap();
        assert_eq!(result.solution, baseline);
        let q = result.evidence;
        assert_eq!(q.pre_qc_degrees_of_freedom, 3);
        assert_eq!(q.summary.raim_attempts, 8);
        assert_eq!(q.raim_candidates.len(), 8);
        assert_eq!(q.summary.outlier_rejections + q.summary.raim_rejections, 0);
        assert!(q.mad_refit.is_none());
        assert!(q.summary.independent_clock_gdop > 0.0 && q.summary.independent_clock_gdop < 50.0);
        assert!(
            distance(
                result.solution.spp.position.position_ecef.unwrap(),
                truth(i)
            ) < 0.003
        );
    }
}

#[test]
fn injected_outlier_is_removed_by_one_raim_refit_with_raw_and_clocks_preserved() {
    let (mut epochs, nav) = data();
    let e = &mut epochs[0];
    let bad = e.observations[0].satellite;
    *e.observations[0].pseudorange_m.as_mut().unwrap() += 100.0;
    e.observations[0].lli = Some(3);
    let before = e.clone();
    let result = solve(
        e,
        &nav,
        config(),
        VARIANCE,
        SppQualityConfig {
            outlier_threshold_sigma: None,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(*e, before);
    assert_eq!(result.evidence.raim_excluded_satellite, Some(bad));
    assert_eq!(result.evidence.summary.raim_attempts, 8);
    assert_eq!(result.evidence.summary.raim_rejections, 1);
    assert_eq!(result.solution.spp.position.num_satellites, 7);
    assert!(
        distance(
            result.solution.spp.position.position_ecef.unwrap(),
            truth(0)
        ) < 0.003
    );
    assert!((result.solution.clock_biases_s[0].1 - 2e-6).abs() < 1e-11);
    assert!((result.solution.clock_biases_s[1].1 + 3e-6).abs() < 1e-11);
    assert_eq!(
        result.solution.observations[0].disposition,
        ObservationDisposition::RaimFdeRejected
    );
    assert_eq!(result.solution.observations.len(), e.observations.len());
    for (a, b) in result
        .solution
        .observations
        .iter()
        .zip(&before.observations)
    {
        assert_eq!(a.raw_observation, *b);
    }
}

#[test]
fn native_mad_short_circuit_preserves_raw_but_can_swamp_a_good_satellite() {
    let (mut epochs, nav) = data();
    let e = &mut epochs[0];
    let bad = e.observations[0].satellite;
    *e.observations[0].pseudorange_m.as_mut().unwrap() += 500.0;
    let mut alternate = e.observations[4].clone();
    alternate.tracking_code = "1X".to_string();
    e.observations.push(alternate);
    let mut unsupported = e.observations[0].clone();
    unsupported.tracking_code = "2W".to_string();
    e.observations.push(unsupported);
    let before = e.clone();
    let result = solve(e, &nav, config(), VARIANCE, SppQualityConfig::default()).unwrap();
    assert!(result.evidence.mad_refit_accepted, "{:?}", result.evidence);
    assert!(result.evidence.summary.outlier_rejections > 0);
    assert_eq!(result.evidence.summary.raim_attempts, 0);
    assert!(result.evidence.raim_candidates.is_empty());
    // Native MAD control flow is not an integrity guarantee: leverage masks
    // G01's injected error and swamps the otherwise good G07 observation.
    assert_eq!(
        result
            .evidence
            .mad_refit
            .as_ref()
            .unwrap()
            .excluded_satellites,
        vec!["G07".parse().unwrap()]
    );
    assert_ne!(
        result
            .evidence
            .mad_refit
            .as_ref()
            .unwrap()
            .excluded_satellites[0],
        bad
    );
    assert_eq!(
        result.solution.observations[0].disposition,
        ObservationDisposition::Used
    );
    assert_eq!(
        result.solution.observations[2].disposition,
        ObservationDisposition::OutlierRejected
    );
    assert!(
        distance(
            result.solution.spp.position.position_ecef.unwrap(),
            truth(0)
        ) > 700.0
    );
    assert_eq!(
        result.solution.observations[8].disposition,
        ObservationDisposition::AlternateTracking
    );
    assert_eq!(
        result.solution.observations[9].disposition,
        ObservationDisposition::UnsupportedTracking
    );
    assert_eq!(
        result.solution.observations.len(),
        before.observations.len()
    );
    for (a, b) in result
        .solution
        .observations
        .iter()
        .zip(&before.observations)
    {
        assert_eq!(a.raw_observation, *b);
    }
    assert_eq!(*e, before);
}

#[test]
fn geometry_rejected_raim_candidates_leave_baseline_and_each_clock_intact() {
    let (mut epochs, nav) = data();
    let e = &mut epochs[0];
    *e.observations[0].pseudorange_m.as_mut().unwrap() += 100.0;
    let no_refit = SppQualityConfig {
        outlier_threshold_sigma: None,
        enable_raim_fde: false,
        max_gdop: None,
        ..Default::default()
    };
    let baseline = solve(e, &nav, config(), VARIANCE, no_refit).unwrap();
    let result = solve(
        e,
        &nav,
        config(),
        VARIANCE,
        SppQualityConfig {
            enable_raim_fde: true,
            max_gdop: Some(baseline.evidence.summary.independent_clock_gdop),
            ..no_refit
        },
    )
    .unwrap();
    assert_eq!(result.solution, baseline.solution);
    assert_eq!(result.evidence.summary.raim_attempts, 8);
    assert_eq!(result.evidence.summary.raim_rejections, 0);
    assert!(
        result
            .evidence
            .raim_candidates
            .iter()
            .all(|c| c.error == Some(Error::SppGeometryTooWeak)
                && c.residual_rms_m.is_some()
                && c.independent_clock_gdop.is_some())
    );
}

#[test]
fn mixed_clock_degrees_of_freedom_prevent_fde_at_only_one_redundant_row() {
    let (mut epochs, nav) = data();
    let e = &mut epochs[0];
    let mut counts = std::collections::BTreeMap::new();
    e.observations.retain(|o| {
        let n = counts.entry(o.satellite.system()).or_insert(0);
        *n += 1;
        *n <= 3
    });
    let q = solve(e, &nav, config(), VARIANCE, SppQualityConfig::default()).unwrap();
    assert_eq!(
        q.evidence.pre_qc_clock_systems,
        vec![GnssSystem::Gps, GnssSystem::Galileo]
    );
    assert_eq!(q.evidence.pre_qc_degrees_of_freedom, 1);
    assert_eq!(q.evidence.summary.raim_attempts, 0);
    assert!(q.evidence.raim_candidates.is_empty());
}

#[test]
fn gdop_and_invalid_qc_are_explicit_failures_and_calls_are_stateless() {
    let (epochs, nav) = data();
    let e = &epochs[0];
    let mut q = SppQualityConfig {
        outlier_threshold_sigma: None,
        enable_raim_fde: false,
        ..Default::default()
    };
    let good = solve(e, &nav, config(), VARIANCE, q).unwrap();
    let limit = good.evidence.summary.independent_clock_gdop;
    q.max_gdop = Some(limit);
    assert!(solve(e, &nav, config(), VARIANCE, q).is_ok());
    q.max_gdop = Some(f64::from_bits(limit.to_bits() - 1));
    assert_eq!(
        solve(e, &nav, config(), VARIANCE, q),
        Err(Error::SppGeometryTooWeak)
    );
    for value in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        q.max_gdop = Some(value);
        assert_eq!(
            solve(e, &nav, config(), VARIANCE, q),
            Err(Error::InvalidConfiguration)
        );
    }
    assert_eq!(
        solve(e, &nav, config(), VARIANCE, SppQualityConfig::default())
            .unwrap()
            .solution,
        good.solution
    );
}

#[test]
fn rtk_qc_seed_keeps_clean_float_covariance_and_reports_geometry_fallback() {
    use gnss_rust::rtk::{
        RtkConstellations, RtkFloatConfig, RtkFloatFilter, RtkPositionMode, RtkSeedSource,
    };
    let (epochs, nav) = data();
    let base = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_mixed_base.obs"
    )))
    .unwrap();
    let base_position = base.header().approximate_position.unwrap();
    let bases = base.collect::<Result<Vec<_>, _>>().unwrap();
    let seed = epochs[0].approximate_position.unwrap();
    let cfg = RtkFloatConfig {
        constellations: RtkConstellations::GpsGalileo,
        position_mode: RtkPositionMode::Kinematic,
        kinematic_spp_seed: Some(config()),
        kinematic_spp_constellations: RtkConstellations::GpsGalileo,
        kinematic_spp_variance: VARIANCE,
        ..Default::default()
    };
    let mut baseline = RtkFloatFilter::new(base_position, seed, cfg).unwrap();
    let mut quality = RtkFloatFilter::new(
        base_position,
        seed,
        RtkFloatConfig {
            kinematic_spp_quality: Some(SppQualityConfig::default()),
            ..cfg
        },
    )
    .unwrap();
    let mut fallback = RtkFloatFilter::new(
        base_position,
        seed,
        RtkFloatConfig {
            kinematic_spp_quality: Some(SppQualityConfig {
                max_gdop: Some(1.0),
                ..Default::default()
            }),
            ..cfg
        },
    )
    .unwrap();
    for (i, (e, b)) in epochs.iter().zip(&bases).enumerate() {
        let old = baseline.process_pair(e, b, &nav).unwrap();
        let new = quality.process_pair(e, b, &nav).unwrap();
        assert_eq!(new.position, old.position);
        assert_eq!(new.position_seed.source, RtkSeedSource::SinglePoint);
        assert!(new.position_seed.spp_quality.is_some());
        let failed = fallback.process_pair(e, b, &nav).unwrap();
        assert_eq!(
            failed.position_seed.spp_failure,
            Some(Error::SppGeometryTooWeak)
        );
        assert_eq!(
            failed.position_seed.source,
            if i == 0 {
                RtkSeedSource::InitialPosition
            } else {
                RtkSeedSource::PreviousRtk
            }
        );
        assert!(failed.position_seed.spp_quality.is_none());
        assert!(failed.position_seed.spp_gps_clock_bias_s.is_none());
        assert!(failed.position_seed.spp_galileo_clock_bias_s.is_none());
    }
    assert!(
        RtkFloatFilter::new(
            base_position,
            seed,
            RtkFloatConfig {
                kinematic_spp_quality: Some(SppQualityConfig::default()),
                ..Default::default()
            }
        )
        .is_err()
    );
}
