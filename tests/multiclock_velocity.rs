use gnss_rust::velocity::{
    DopplerObservation, VelocityClockModel, VelocityConfig, VelocityError, solve_velocity,
};
use gnss_rust::{EcefCoord, GnssSystem};

struct ReferenceCase {
    name: String,
    clocks: usize,
    receiver: EcefCoord,
    rows: Vec<DopplerObservation>,
    expected: Vec<f64>,
}
fn cases() -> Vec<ReferenceCase> {
    include_str!("fixtures/upstream_multiclock_velocity.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|line| {
            let f: Vec<_> = line.split(',').collect();
            let n: usize = f[1].parse().unwrap();
            let clocks = f[2].parse().unwrap();
            let mut i = 3;
            let take = |i: &mut usize| {
                let v = f[*i].parse::<f64>().unwrap();
                *i += 1;
                v
            };
            let receiver = EcefCoord::new(take(&mut i), take(&mut i), take(&mut i));
            let mut rows = Vec::new();
            for _ in 0..n {
                let satellite = f[i].parse().unwrap();
                i += 1;
                rows.push(DopplerObservation {
                    satellite,
                    satellite_position: EcefCoord::new(take(&mut i), take(&mut i), take(&mut i)),
                    satellite_velocity_mps: [take(&mut i), take(&mut i), take(&mut i)],
                    satellite_clock_drift_sps: take(&mut i),
                    doppler_hz: take(&mut i),
                    frequency_hz: take(&mut i),
                    elevation_rad: take(&mut i),
                });
            }
            ReferenceCase {
                name: f[0].into(),
                clocks,
                receiver,
                rows,
                expected: f[i..].iter().map(|s| s.parse().unwrap()).collect(),
            }
        })
        .collect()
}
fn config() -> VelocityConfig {
    VelocityConfig {
        clock_model: VelocityClockModel::PerSystem,
        ..Default::default()
    }
}
fn close(a: f64, b: f64, t: f64) {
    assert!(
        (a - b).abs() <= t * b.abs().max(1.0),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn two_and_three_clock_velocity_matches_sixteen_independent_eigen_cases() {
    let all = cases();
    assert_eq!(all.len(), 16);
    for case in all {
        let actual = solve_velocity(&case.rows, case.receiver, config());
        assert_eq!(
            actual.is_ok(),
            case.expected[0] == 1.0,
            "{}: {actual:?}",
            case.name
        );
        if let Ok(result) = actual {
            let n = 3 + case.clocks;
            assert_eq!(result.clock_systems.len(), case.clocks);
            for i in 0..3 {
                close(result.velocity_mps[i], case.expected[1 + i], 1e-9);
            }
            for i in 0..case.clocks {
                close(
                    result.receiver_clock_drifts_sps[i],
                    case.expected[4 + i] / gnss_rust::constants::SPEED_OF_LIGHT,
                    1e-16,
                );
            }
            let covariance = result.full_covariance_m2s2.as_ref().unwrap();
            assert_eq!(covariance.len(), n);
            for (i, row) in covariance.iter().enumerate() {
                for (j, value) in row.iter().enumerate() {
                    close(*value, case.expected[1 + n + i * n + j], 1e-9);
                }
            }
            close(
                result.chi_square_per_dof,
                case.expected[1 + n + n * n],
                1e-9,
            );
            close(result.covariance_scale, case.expected[2 + n + n * n], 1e-9);
            close(result.residual_rms_mps, case.expected[3 + n + n * n], 1e-9);
        }
    }
}

#[test]
fn changing_one_clock_rate_leaves_velocity_and_other_clocks_unchanged() {
    let mut case = cases().remove(0);
    let before = solve_velocity(&case.rows, case.receiver, config()).unwrap();
    for o in &mut case.rows {
        if o.satellite.system() == GnssSystem::Galileo {
            o.doppler_hz += 5.25;
        }
    }
    let after = solve_velocity(&case.rows, case.receiver, config()).unwrap();
    for i in 0..3 {
        close(after.velocity_mps[i], before.velocity_mps[i], 1e-10);
    }
    let gal = after
        .clock_systems
        .iter()
        .position(|s| *s == GnssSystem::Galileo)
        .unwrap();
    close(
        (after.receiver_clock_drifts_sps[gal] - before.receiver_clock_drifts_sps[gal])
            * gnss_rust::constants::SPEED_OF_LIGHT,
        -5.25 * gnss_rust::constants::GPS_L1_WAVELENGTH,
        1e-10,
    );
    close(
        after.receiver_clock_drifts_sps[0],
        before.receiver_clock_drifts_sps[0],
        1e-18,
    );
    assert_eq!(after.full_covariance_m2s2, before.full_covariance_m2s2);
}

#[test]
fn a_single_system_keeps_the_old_path_and_a_singleton_clock_adds_no_velocity_information() {
    let case = cases().remove(0);
    // Use eight broadly distributed GPS directions: the four alternating
    // GPS rows alone legitimately fail the unchanged 400 (m/s)^2 bound.
    let gps: Vec<_> = case
        .rows
        .iter()
        .copied()
        .enumerate()
        .map(|(i, mut row)| {
            if row.satellite.system() == GnssSystem::Galileo {
                row.doppler_hz -= 2.4 / gnss_rust::constants::GPS_L1_WAVELENGTH;
            }
            row.satellite = format!("G{:02}", i + 1).parse().unwrap();
            row
        })
        .collect();
    let before = solve_velocity(&gps, case.receiver, VelocityConfig::default()).unwrap();
    assert_eq!(
        solve_velocity(&gps, case.receiver, config()).unwrap(),
        before
    );
    let mut mixed = gps.clone();
    mixed.push(
        *case
            .rows
            .iter()
            .find(|r| r.satellite.system() == GnssSystem::Galileo)
            .unwrap(),
    );
    let after = solve_velocity(&mixed, case.receiver, config()).unwrap();
    for i in 0..3 {
        close(after.velocity_mps[i], before.velocity_mps[i], 1e-9);
        for j in 0..3 {
            close(
                after.covariance_m2s2[i][j],
                before.covariance_m2s2[i][j],
                1e-9,
            );
        }
    }
    assert_eq!(
        after.clock_systems,
        vec![GnssSystem::Gps, GnssSystem::Galileo]
    );
    assert_eq!(after.full_covariance_m2s2.unwrap().len(), 5);
}

#[test]
fn row_order_preserves_full_joint_solution_and_clock_labels() {
    let mut case = cases().into_iter().find(|c| c.clocks == 3).unwrap();
    let before = solve_velocity(&case.rows, case.receiver, config()).unwrap();
    case.rows.reverse();
    let after = solve_velocity(&case.rows, case.receiver, config()).unwrap();
    assert_eq!(before.clock_systems, after.clock_systems);
    for i in 0..3 {
        close(after.velocity_mps[i], before.velocity_mps[i], 1e-10);
    }
    for (a, b) in after
        .full_covariance_m2s2
        .unwrap()
        .iter()
        .flatten()
        .zip(before.full_covariance_m2s2.unwrap().iter().flatten())
    {
        close(*a, *b, 1e-10);
    }
}

#[test]
fn raw_galileo_selection_prefers_usable_doppler_and_retains_alternate_groups() {
    use gnss_rust::rinex::{RinexObservationReader, read_navigation};
    use gnss_rust::velocity::{
        DopplerRejection, VelocityEpochConfig, solve_gps_galileo_e1_epoch as solve,
    };
    use std::io::Cursor;
    let mut reader = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_mixed_rover.obs"
    )))
    .unwrap();
    let position = reader.header().approximate_position.unwrap();
    let mut epoch = reader.next().unwrap().unwrap();
    let nav = read_navigation(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_mixed.nav"
    )))
    .unwrap();
    for o in &mut epoch.observations {
        o.doppler_hz = Some(0.0);
    }
    let index = epoch
        .observations
        .iter()
        .position(|o| o.satellite.system() == GnssSystem::Galileo)
        .unwrap();
    let satellite = epoch.observations[index].satellite;
    let mut alternate = epoch.observations[index].clone();
    alternate.tracking_code = "1X".into();
    alternate.doppler_hz = Some(5.0);
    epoch.observations.push(alternate);
    let cfg = VelocityEpochConfig {
        velocity: config(),
        ..Default::default()
    };
    let report = solve(&epoch, &nav, position, cfg).unwrap();
    assert_eq!(report.prepared.len(), 8);
    assert_eq!(
        report
            .prepared
            .iter()
            .find(|p| p.raw.satellite == satellite)
            .unwrap()
            .raw
            .tracking_code,
        "1C"
    );
    assert!(
        report
            .rejected
            .iter()
            .any(|(r, why)| r.satellite == satellite
                && r.tracking_code == "1X"
                && *why == DopplerRejection::AlternativeTrackingCode)
    );
    epoch.observations[index].doppler_hz = None;
    let report = solve(&epoch, &nav, position, cfg).unwrap();
    assert_eq!(
        report
            .prepared
            .iter()
            .find(|p| p.raw.satellite == satellite)
            .unwrap()
            .raw
            .tracking_code,
        "1X"
    );
    assert_eq!(
        solve(&epoch, &nav, position, VelocityEpochConfig::default()),
        Err(gnss_rust::Error::InvalidConfiguration)
    );
    assert_eq!(
        solve_velocity(
            &case_rows_for_rejection(),
            position,
            VelocityConfig::default()
        ),
        Err(VelocityError::MultipleClockSystems)
    );
}
fn case_rows_for_rejection() -> Vec<DopplerObservation> {
    cases().remove(0).rows
}

#[test]
fn three_satellites_per_system_recover_velocity_without_assuming_a_shared_clock() {
    let case = cases()
        .into_iter()
        .find(|c| c.name == "three_each")
        .unwrap();
    for system in [GnssSystem::Gps, GnssSystem::Galileo] {
        let rows: Vec<_> = case
            .rows
            .iter()
            .filter(|r| r.satellite.system() == system)
            .copied()
            .collect();
        assert_eq!(rows.len(), 3);
        assert_eq!(
            solve_velocity(&rows, case.receiver, VelocityConfig::default()),
            Err(VelocityError::InsufficientObservations)
        );
    }
    let solution = solve_velocity(&case.rows, case.receiver, config()).unwrap();
    for (actual, expected) in solution.velocity_mps.into_iter().zip([-3.2, 4.5, 1.2]) {
        close(actual, expected, 1e-9);
    }
    assert_ne!(
        solution.receiver_clock_drifts_sps[0],
        solution.receiver_clock_drifts_sps[1]
    );
}
