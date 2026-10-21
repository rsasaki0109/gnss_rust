use gnss_rust::velocity::{DopplerObservation, VelocityConfig, VelocityError, solve_velocity};
use gnss_rust::{EcefCoord, Error};

type NativeCase = (
    String,
    EcefCoord,
    VelocityConfig,
    Vec<DopplerObservation>,
    Vec<f64>,
);
fn cases() -> Vec<NativeCase> {
    include_str!("fixtures/upstream_velocity.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|line| {
            let mut fields = line.split(',');
            let name = fields.next().unwrap().to_string();
            let n: usize = fields.next().unwrap().parse().unwrap();
            let v: Vec<f64> = fields.map(|s| s.parse().unwrap()).collect();
            let receiver = EcefCoord::new(v[0], v[1], v[2]);
            let config = VelocityConfig {
                doppler_sigma_mps: v[3],
                min_satellites: v[4] as usize,
                ..Default::default()
            };
            let rows = (0..n)
                .map(|i| {
                    let a = &v[5 + 10 * i..5 + 10 * (i + 1)];
                    DopplerObservation {
                        satellite: format!("G{:02}", i + 1).parse().unwrap(),
                        satellite_position: EcefCoord::new(a[0], a[1], a[2]),
                        satellite_velocity_mps: [a[3], a[4], a[5]],
                        satellite_clock_drift_sps: a[6],
                        doppler_hz: a[7],
                        frequency_hz: a[8],
                        elevation_rad: a[9],
                    }
                })
                .collect();
            (name, receiver, config, rows, v[5 + 10 * n..].to_vec())
        })
        .collect()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * b.abs().max(1.0),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn doppler_velocity_matches_eleven_native_cases_with_qc_enabled() {
    let mut count = 0;
    for (name, receiver, config, rows, expected) in cases() {
        if name == "native_zero_sentinel" {
            continue;
        }
        let result = solve_velocity(&rows, receiver, config);
        assert_eq!(result.is_ok(), expected[0] == 1.0, "{name}: {result:?}");
        if let Ok(solution) = result {
            assert_eq!(solution.satellites.len(), expected[1] as usize);
            for i in 0..3 {
                close(solution.velocity_mps[i], expected[2 + i], 1e-9);
            }
            close(solution.receiver_clock_drift_sps, expected[5], 1e-16);
            for i in 0..3 {
                for j in 0..3 {
                    close(
                        solution.covariance_m2s2[i][j],
                        expected[6 + 3 * i + j],
                        1e-9,
                    );
                }
            }
            close(
                solution.covariance_m2s2[3][3] / gnss_rust::constants::SPEED_OF_LIGHT.powi(2),
                expected[15],
                1e-23,
            );
            close(solution.residual_rms_mps, expected[16], 1e-9);
        }
        count += 1;
    }
    assert_eq!(count, 11);
}

#[test]
fn measured_zero_is_active_and_a_priori_covariance_never_shrinks() {
    let (_, receiver, cfg, rows, expected) = cases()
        .into_iter()
        .find(|c| c.0 == "native_zero_sentinel")
        .unwrap();
    assert_eq!(expected[0], 0.0); // native drops its numeric-zero sentinel
    let result = solve_velocity(&rows, receiver, cfg).unwrap();
    assert_eq!(result.velocity_mps, [0.0; 3]);
    assert_eq!(result.receiver_clock_drift_sps, 0.0);
    assert_eq!(result.satellites.len(), 8);
    assert_eq!(result.covariance_scale, 1.0);
    assert_eq!(result.covariance_m2s2, result.unscaled_covariance_m2s2);
    assert!(result.covariance_m2s2[0][0] > 0.0);
}

#[test]
fn exactly_four_observations_can_fail_geometry_despite_near_zero_residual() {
    let (_, receiver, cfg, rows, _) = cases()
        .into_iter()
        .find(|c| c.0 == "poor_geometry")
        .unwrap();
    let error = solve_velocity(&rows, receiver, cfg).unwrap_err();
    assert!(matches!(error, VelocityError::WeakGeometry { .. }));
    let (_, receiver, cfg, rows, _) = cases().into_iter().find(|c| c.0 == "bad_fit").unwrap();
    assert!(matches!(
        solve_velocity(&rows, receiver, cfg),
        Err(VelocityError::InconsistentFit { .. })
    ));
}

#[test]
fn quality_boundaries_are_inclusive_and_one_ulp_below_is_rejected() {
    let (_, receiver, cfg, rows, _) = cases()
        .into_iter()
        .find(|c| c.0 == "noise_1.000000")
        .unwrap();
    let result = solve_velocity(&rows, receiver, cfg).unwrap();
    let chi = result.chi_square_per_dof;
    let bound = VelocityConfig {
        max_chi_square_per_dof: chi,
        ..cfg
    };
    assert!(solve_velocity(&rows, receiver, bound).is_ok());
    assert!(matches!(
        solve_velocity(
            &rows,
            receiver,
            VelocityConfig {
                max_chi_square_per_dof: f64::from_bits(chi.to_bits() - 1),
                ..bound
            }
        ),
        Err(VelocityError::InconsistentFit { .. })
    ));
    let maximum = (0..3)
        .map(|i| result.unscaled_covariance_m2s2[i][i])
        .fold(0.0, f64::max);
    let bound = VelocityConfig {
        max_unscaled_velocity_variance_m2s2: maximum,
        ..cfg
    };
    assert!(solve_velocity(&rows, receiver, bound).is_ok());
    assert!(matches!(
        solve_velocity(
            &rows,
            receiver,
            VelocityConfig {
                max_unscaled_velocity_variance_m2s2: f64::from_bits(maximum.to_bits() - 1),
                ..bound
            }
        ),
        Err(VelocityError::WeakGeometry { .. })
    ));
}

#[test]
fn invalid_config_rows_duplicates_and_independent_clocks_fail_explicitly() {
    let (_, receiver, cfg, rows, _) = cases().remove(0);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            solve_velocity(
                &rows,
                receiver,
                VelocityConfig {
                    doppler_sigma_mps: bad,
                    ..cfg
                }
            ),
            Err(VelocityError::InvalidConfiguration)
        );
        let mut invalid = rows.clone();
        invalid[0].frequency_hz = bad;
        assert_eq!(
            solve_velocity(&invalid, receiver, cfg),
            Err(VelocityError::InvalidObservation)
        );
    }
    assert_eq!(
        solve_velocity(&rows, EcefCoord::new(f64::NAN, 0.0, 0.0), cfg),
        Err(VelocityError::Basic(Error::NonFinite))
    );
    let mut duplicate = rows.clone();
    duplicate[1].satellite = duplicate[0].satellite;
    assert_eq!(
        solve_velocity(&duplicate, receiver, cfg),
        Err(VelocityError::DuplicateObservation)
    );
    let mut mixed = rows.clone();
    mixed[0].satellite = "E01".parse().unwrap();
    assert_eq!(
        solve_velocity(&mixed, receiver, cfg),
        Err(VelocityError::MultipleClockSystems)
    );
    let mut invalid = rows.clone();
    invalid[0].doppler_hz = f64::NAN;
    assert_eq!(
        solve_velocity(&invalid, receiver, cfg),
        Err(VelocityError::InvalidObservation)
    );
    invalid[0] = rows[0];
    invalid[0].satellite_position = receiver;
    assert_eq!(
        solve_velocity(&invalid, receiver, cfg),
        Err(VelocityError::InvalidObservation)
    );
}

fn epoch_data() -> (
    gnss_rust::observation::ObservationEpoch,
    gnss_rust::navigation::NavigationData,
    EcefCoord,
) {
    use gnss_rust::rinex::{RinexObservationReader, read_navigation};
    use std::io::Cursor;
    let mut reader = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_rover.obs"
    )))
    .unwrap();
    let position = reader.header().approximate_position.unwrap();
    let mut epoch = reader.next().unwrap().unwrap();
    for o in &mut epoch.observations {
        o.doppler_hz = Some(0.0);
    }
    let nav = read_navigation(Cursor::new(include_str!("fixtures/synthetic_rtk.nav"))).unwrap();
    (epoch, nav, position)
}

#[test]
fn raw_preparation_preserves_zero_missing_code_phase_tracking_and_rejection_evidence() {
    use gnss_rust::velocity::{DopplerRejection, VelocityEpochConfig, solve_gps_l1_epoch};
    let (mut epoch, nav, position) = epoch_data();
    let cfg = VelocityEpochConfig::default();
    epoch.observations[0].doppler_hz = None;
    epoch.observations[1].pseudorange_m = None;
    epoch.observations[2].carrier_phase_cycles = None;
    epoch.observations[2].lli = Some(3);
    let mut unsupported = epoch.observations[3].clone();
    unsupported.tracking_code = "1W".into();
    epoch.observations.push(unsupported);
    let report = solve_gps_l1_epoch(&epoch, &nav, position, cfg).unwrap();
    assert_eq!(report.prepared.len(), 8);
    assert_eq!(report.rejected.len(), 2);
    assert_eq!(report.rejected[0].1, DopplerRejection::MissingDoppler);
    assert_eq!(report.rejected[1].1, DopplerRejection::UnsupportedSignal);
    assert!(
        report
            .prepared
            .iter()
            .all(|p| p.raw.doppler_hz == Some(0.0) && p.measurement.doppler_hz == 0.0)
    );
    let absent = &report.prepared[0];
    assert!(absent.used_missing_code_travel_time);
    assert_eq!(absent.raw.pseudorange_m, None);
    let no_phase = &report.prepared[1];
    assert_eq!(no_phase.raw.carrier_phase_cycles, None);
    assert_eq!(no_phase.raw.lli, Some(3));
    for p in &report.prepared {
        assert!(p.transmit_time.difference_seconds(epoch.time) < 0.0);
        assert!(p.reception_frame_rotation_rad > 0.0);
        let state = nav
            .ephemeris(p.raw.satellite, p.transmit_time, 14400.0)
            .unwrap()
            .satellite_state(p.transmit_time)
            .unwrap();
        let norm = |v: [f64; 3]| v[0].hypot(v[1]).hypot(v[2]);
        close(
            norm(p.measurement.satellite_velocity_mps),
            norm(state.velocity_mps),
            1e-12,
        );
    }
}

#[test]
fn raw_preparation_rejects_stale_navigation_bad_code_duplicates_and_epoch_metadata() {
    use gnss_rust::velocity::{DopplerRejection, VelocityEpochConfig, solve_gps_l1_epoch as solve};
    let (mut epoch, nav, position) = epoch_data();
    let cfg = VelocityEpochConfig::default();
    epoch.observations[0].pseudorange_m = Some(0.0);
    epoch.observations[1].doppler_hz = Some(f64::NAN);
    let report = solve(&epoch, &nav, position, cfg).unwrap();
    assert_eq!(report.rejected[0].1, DopplerRejection::InvalidCode);
    assert_eq!(report.rejected[1].1, DopplerRejection::InvalidDoppler);
    let stale = solve(
        &epoch,
        &nav,
        position,
        VelocityEpochConfig {
            max_ephemeris_age_s: 0.0,
            ..cfg
        },
    )
    .unwrap();
    assert!(stale.prepared.is_empty());
    assert!(
        stale
            .rejected
            .iter()
            .skip(2)
            .all(|r| r.1 == DopplerRejection::Navigation(Error::MissingEphemeris))
    );
    assert_eq!(stale.solution, Err(VelocityError::InsufficientObservations));
    epoch.receiver_clock_offset_s = Some(0.001);
    assert_eq!(
        solve(&epoch, &nav, position, cfg),
        Err(Error::InvalidConfiguration)
    );
    epoch.receiver_clock_offset_s = None;
    epoch.flag = 6;
    assert_eq!(
        solve(&epoch, &nav, position, cfg),
        Err(Error::InvalidConfiguration)
    );
    epoch.flag = 0;
    epoch.observations.push(epoch.observations[2].clone());
    assert_eq!(
        solve(&epoch, &nav, position, cfg),
        Err(Error::InvalidConfiguration)
    );
}

#[test]
fn optional_carrier_lock_admission_retains_zero_phase_and_reports_each_failure() {
    use gnss_rust::velocity::{DopplerRejection, VelocityEpochConfig, solve_gps_l1_epoch as solve};
    let (mut epoch, nav, position) = epoch_data();
    epoch.observations[0].carrier_phase_cycles = None;
    epoch.observations[1].lli = Some(1);
    epoch.observations[2].lli = Some(2);
    epoch.observations[3].carrier_phase_cycles = Some(f64::NAN);
    epoch.observations[4].carrier_phase_cycles = Some(0.0);
    let legacy = solve(&epoch, &nav, position, VelocityEpochConfig::default()).unwrap();
    assert_eq!(legacy.prepared.len(), 9);
    let locked = solve(
        &epoch,
        &nav,
        position,
        VelocityEpochConfig {
            require_carrier_lock: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(locked.prepared.len(), 5);
    assert_eq!(
        locked
            .rejected
            .iter()
            .map(|(_, why)| *why)
            .collect::<Vec<_>>(),
        vec![
            DopplerRejection::MissingCarrierPhase,
            DopplerRejection::CarrierLossOfLock,
            DopplerRejection::CarrierLossOfLock,
            DopplerRejection::InvalidCarrierPhase
        ]
    );
    assert!(
        locked
            .prepared
            .iter()
            .any(|p| p.raw.satellite == epoch.observations[4].satellite
                && p.raw.carrier_phase_cycles == Some(0.0))
    );
    assert!(legacy.rejected.is_empty());
}
