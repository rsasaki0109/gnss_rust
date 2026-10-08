use gnss_rust::coordinates::geodetic_to_ecef;
use gnss_rust::models::{KlobucharParameters, klobuchar, saastamoinen};
use gnss_rust::rinex::{RinexObservationReader, read_navigation};
use gnss_rust::spp::{RejectionReason, SppConfig, solve_epoch};
use gnss_rust::{EcefCoord, Error, GeodeticCoord, GnssTime, SatelliteId};
use std::io::Cursor;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} (tolerance {tolerance})"
    );
}
fn vacuum() -> SppConfig {
    SppConfig {
        use_troposphere: false,
        use_ionosphere: false,
        ..SppConfig::default()
    }
}
fn epoch() -> gnss_rust::observation::ObservationEpoch {
    RinexObservationReader::new(Cursor::new(include_str!("fixtures/synthetic_spp.obs")))
        .unwrap()
        .next()
        .unwrap()
        .unwrap()
}
fn nav() -> gnss_rust::navigation::NavigationData {
    read_navigation(Cursor::new(include_str!("fixtures/synthetic_spp.nav"))).unwrap()
}

#[test]
fn atmosphere_models_match_four_unmodified_cpp_reference_cases() {
    let rows: Vec<_> = include_str!("fixtures/upstream_atmosphere.csv")
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(rows.len(), 4);
    for row in rows {
        let v: Vec<f64> = row.split(',').map(|s| s.parse().unwrap()).collect();
        assert_eq!(v.len(), 8);
        let geo = GeodeticCoord::new(v[0], v[1], v[2]).unwrap();
        near(
            klobuchar(geo, v[3], v[4], v[5], KlobucharParameters::default()).unwrap(),
            v[6],
            1e-9,
        );
        near(saastamoinen(geo, v[4]).unwrap(), v[7], 1e-8);
    }
    near(
        klobuchar(
            GeodeticCoord::new(35_f64.to_radians(), 139_f64.to_radians(), 45.0).unwrap(),
            120_f64.to_radians(),
            30_f64.to_radians(),
            16000.0,
            KlobucharParameters::default(),
        )
        .unwrap(),
        7.691021650363509,
        1e-9,
    );
}

#[test]
fn atmosphere_bounds_and_invalid_inputs() {
    let geo = GeodeticCoord::default();
    assert_eq!(saastamoinen(geo, 0.01).unwrap(), 0.0);
    assert_eq!(
        saastamoinen(GeodeticCoord::new(0.0, 0.0, 11000.0).unwrap(), 1.0).unwrap(),
        0.0
    );
    assert!(saastamoinen(geo, f64::NAN).is_err());
    assert!(klobuchar(geo, 0.0, -0.1, 0.0, KlobucharParameters::default()).is_err());
    let p = KlobucharParameters {
        alpha: [f64::NAN; 4],
        ..KlobucharParameters::default()
    };
    assert!(klobuchar(geo, 0.0, 1.0, 0.0, p).is_err());
}

#[test]
fn rinex_to_spp_recovers_all_eight_cpp_generated_vacuum_epochs() {
    let expected = geodetic_to_ecef(
        GeodeticCoord::new(35_f64.to_radians(), 139_f64.to_radians(), 45.0).unwrap(),
    )
    .unwrap();
    let nav = nav();
    let epochs =
        RinexObservationReader::new(Cursor::new(include_str!("fixtures/synthetic_spp.obs")))
            .unwrap();
    let mut count = 0;
    for epoch in epochs {
        let epoch = epoch.unwrap();
        let solution = solve_epoch(&epoch, &nav, vacuum()).unwrap();
        let position = solution.position.position_ecef.unwrap();
        assert!(solution.position.is_valid());
        assert_eq!(solution.position.num_satellites, 9);
        assert!(solution.iterations <= 5);
        assert!(solution.residual_rms_m < 0.001);
        for (actual, expected) in position.to_array().into_iter().zip(expected.to_array()) {
            near(actual, expected, 0.003);
        }
        near(solution.position.receiver_clock_bias, 2e-6, 1e-11);
        let covariance = solution.position.position_covariance.unwrap();
        for (i, row) in covariance.iter().enumerate() {
            assert!(row[i] > 0.0);
            for (j, &value) in row.iter().enumerate() {
                near(value, covariance[j][i], 1e-12);
            }
        }
        count += 1;
    }
    assert_eq!(count, 8);
}

#[test]
fn default_atmosphere_spp_recovers_cpp_generated_ionosphere_and_troposphere() {
    let expected = geodetic_to_ecef(
        GeodeticCoord::new(35_f64.to_radians(), 139_f64.to_radians(), 45.0).unwrap(),
    )
    .unwrap();
    let nav = nav();
    let epochs = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_spp_atmosphere.obs"
    )))
    .unwrap();
    let mut count = 0;
    for epoch in epochs {
        let epoch = epoch.unwrap();
        let solution = solve_epoch(&epoch, &nav, SppConfig::default()).unwrap();
        for (actual, expected) in solution
            .position
            .position_ecef
            .unwrap()
            .to_array()
            .into_iter()
            .zip(expected.to_array())
        {
            near(actual, expected, 0.003);
        }
        near(solution.position.receiver_clock_bias, 2e-6, 1e-11);
        assert!(solution.residual_rms_m < 0.001);
        count += 1;
    }
    assert_eq!(count, 8);
}

#[test]
fn spp_can_start_without_approximate_position() {
    let mut epoch = epoch();
    epoch.approximate_position = None;
    let solution = solve_epoch(&epoch, &nav(), vacuum()).unwrap();
    let expected = geodetic_to_ecef(
        GeodeticCoord::new(35_f64.to_radians(), 139_f64.to_radians(), 45.0).unwrap(),
    )
    .unwrap();
    for (actual, expected) in solution
        .position
        .position_ecef
        .unwrap()
        .to_array()
        .into_iter()
        .zip(expected.to_array())
    {
        near(actual, expected, 0.003);
    }
}

#[test]
fn duplicate_tracking_signals_do_not_inflate_satellite_count() {
    let mut epoch = epoch();
    let mut alternate = epoch.observations[0].clone();
    alternate.tracking_code = "1W".to_owned();
    alternate.pseudorange_m = Some(alternate.pseudorange_m.unwrap() + 1000.0);
    epoch.observations.insert(0, alternate);
    let result = solve_epoch(&epoch, &nav(), vacuum()).unwrap();
    assert_eq!(result.position.num_satellites, 9);
    assert!(result.residual_rms_m < 0.001);
}

#[test]
fn missing_code_or_navigation_is_reported_in_successful_epoch() {
    let mut epoch = epoch();
    let sat = epoch.observations[0].satellite;
    epoch
        .observations
        .iter_mut()
        .filter(|o| o.satellite == sat)
        .for_each(|o| o.pseudorange_m = None);
    let solution = solve_epoch(&epoch, &nav(), vacuum()).unwrap();
    assert_eq!(solution.position.num_satellites, 8);
    assert!(
        solution
            .rejected
            .contains(&(sat, RejectionReason::MissingCode))
    );
    let mut missing = epoch
        .observations
        .iter()
        .find(|o| o.tracking_code == "1C" && o.pseudorange_m.is_some())
        .unwrap()
        .clone();
    missing.satellite = "G32".parse().unwrap();
    epoch.observations.push(missing);
    let solution = solve_epoch(&epoch, &nav(), vacuum()).unwrap();
    assert!(
        solution
            .rejected
            .contains(&("G32".parse().unwrap(), RejectionReason::MissingEphemeris))
    );
}

#[test]
fn insufficient_stale_phase_only_or_invalid_geometry_cannot_succeed() {
    let nav = nav();
    let mut fewer = epoch();
    let satellites: Vec<SatelliteId> = fewer.satellites().into_iter().take(3).collect();
    fewer
        .observations
        .retain(|o| satellites.contains(&o.satellite));
    assert_eq!(
        solve_epoch(&fewer, &nav, vacuum()).unwrap_err(),
        Error::InsufficientObservations
    );
    let mut phase_only = epoch();
    phase_only
        .observations
        .iter_mut()
        .for_each(|o| o.pseudorange_m = None);
    assert_eq!(
        solve_epoch(&phase_only, &nav, vacuum()).unwrap_err(),
        Error::InsufficientObservations
    );
    let mut stale = epoch();
    stale.time = GnssTime::new(2301, 345600.0).unwrap();
    assert_eq!(
        solve_epoch(&stale, &nav, vacuum()).unwrap_err(),
        Error::InsufficientObservations
    );
    let mut invalid = epoch();
    invalid.approximate_position = Some(EcefCoord::new(f64::NAN, 0.0, 0.0));
    assert_eq!(
        solve_epoch(&invalid, &nav, vacuum()).unwrap_err(),
        Error::NonFinite
    );
    let mut degenerate = gnss_rust::navigation::NavigationData::default();
    let eph = nav.records().next().unwrap();
    for original in nav.records() {
        let mut clone = eph.clone();
        clone.satellite = original.satellite;
        degenerate.add_ephemeris(clone).unwrap();
    }
    // Remove clock-dependent transmit-time differences so identical satellite
    // orbits produce genuinely rank-deficient geometry.
    let mut epoch = epoch();
    epoch
        .observations
        .iter_mut()
        .for_each(|o| o.pseudorange_m = Some(26_500_000.0));
    assert_eq!(
        solve_epoch(&epoch, &degenerate, vacuum()).unwrap_err(),
        Error::SingularGeometry
    );
}

#[test]
fn bad_solver_configuration_or_nonconvergence_are_errors() {
    let epoch = epoch();
    let nav = nav();
    for config in [
        SppConfig {
            max_iterations: 0,
            ..vacuum()
        },
        SppConfig {
            convergence_m: f64::NAN,
            ..vacuum()
        },
        SppConfig {
            pseudorange_sigma_m: 0.0,
            ..vacuum()
        },
        SppConfig {
            max_ephemeris_age_s: -1.0,
            ..vacuum()
        },
        SppConfig {
            elevation_mask_rad: std::f64::consts::FRAC_PI_2,
            ..vacuum()
        },
    ] {
        assert_eq!(
            solve_epoch(&epoch, &nav, config).unwrap_err(),
            Error::InvalidConfiguration
        );
    }
    assert_eq!(
        solve_epoch(
            &epoch,
            &nav,
            SppConfig {
                max_iterations: 1,
                ..vacuum()
            }
        )
        .unwrap_err(),
        Error::NonConvergent
    );
}
