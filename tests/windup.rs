use gnss_rust::windup::*;
use gnss_rust::{EcefCoord, Error, GnssTime};
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        a.is_finite() && (a - b).abs() <= tol,
        "{a:.17} vs {b:.17}, delta={}",
        (a - b).abs()
    );
}

#[test]
fn native_12_solar_positions_include_fractional_time_and_week_boundaries() {
    let mut count = 0;
    for line in include_str!("fixtures/upstream_windup_sun.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = line.split(',').collect();
        let time = GnssTime::new(f[0].parse().unwrap(), f[1].parse().unwrap()).unwrap();
        let actual = approximate_sun_position_ecef(time).unwrap();
        for (i, x) in actual.to_array().iter().enumerate() {
            close(*x, f[2 + i].parse().unwrap(), 0.001);
        }
        count += 1;
    }
    assert_eq!(count, 12);
    assert!(matches!(
        approximate_sun_position_ecef(GnssTime::new(i32::MAX, 0.0).unwrap()),
        Err(WindupError::Gnss(Error::OutOfRange))
    ));
}
#[test]
fn native_360_nominal_yaw_windup_cycles_match_both_hemispheres_and_gauges() {
    let mut count = 0;
    for line in include_str!("fixtures/upstream_windup.csv").lines().skip(1) {
        let n: Vec<f64> = line.split(',').map(|s| s.parse().unwrap()).collect();
        let p = |i| EcefCoord::new(n[i], n[i + 1], n[i + 2]);
        let actual = nominal_yaw_windup(p(0), p(3), p(6), n[9]).unwrap();
        close(actual.cycles, n[10], 2e-12);
        assert!((-0.5..=0.5).contains(&actual.principal_cycles));
        assert!((actual.cycles - n[9]).abs() <= 0.5);
        count += 1;
    }
    assert_eq!(count, 360);
}
fn frame(angle: f64) -> DipoleFrame {
    let (s, c) = angle.sin_cos();
    DipoleFrame {
        x: [c, s, 0.0],
        y: [s, -c, 0.0],
    }
}
const RECEIVER: DipoleFrame = DipoleFrame {
    x: [1.0, 0.0, 0.0],
    y: [0.0, 1.0, 0.0],
};
#[test]
fn complete_rotations_unwrap_branch_cuts_without_cycle_jumps() {
    for direction in [-1.0, 1.0] {
        let mut prior = 0.0;
        let mut branch_crossings = 0;
        let mut last_principal = 0.0;
        for i in 0..=72 {
            let true_cycles = direction * f64::from(i) / 24.0;
            let actual = phase_windup_from_frames(
                [0.0, 0.0, -2e7],
                frame(true_cycles * std::f64::consts::TAU),
                RECEIVER,
                prior,
            )
            .unwrap();
            close(actual.cycles, true_cycles, 2e-12);
            assert!((actual.cycles - prior).abs() < 0.05);
            if (actual.principal_cycles - last_principal).abs() > 0.5 {
                branch_crossings += 1;
            }
            prior = actual.cycles;
            last_principal = actual.principal_cycles;
        }
        assert_eq!(branch_crossings, 3);
    }
}
#[test]
fn fresh_arc_reference_keeps_integer_branch_fixed_across_relinearizations() {
    let first = phase_windup_from_frames(
        [0.0, 0.0, -2e7],
        frame(0.4999 * std::f64::consts::TAU),
        RECEIVER,
        0.0,
    )
    .unwrap();
    let relinearized = phase_windup_from_frames(
        [0.0, 0.0, -2e7],
        frame(0.5001 * std::f64::consts::TAU),
        RECEIVER,
        first.cycles,
    )
    .unwrap();
    close(first.cycles, 0.4999, 2e-12);
    close(relinearized.cycles, 0.5001, 2e-12);
    assert!(relinearized.principal_cycles < 0.0);
}

#[test]
fn exact_half_cycle_ties_follow_native_round_away_from_zero() {
    for (previous, expected) in [(-1.5, -2.0), (-0.5, -1.0), (0.5, 1.0), (1.5, 2.0)] {
        let actual =
            phase_windup_from_frames([0.0, 0.0, -2e7], frame(0.0), RECEIVER, previous).unwrap();
        assert_eq!(actual.principal_cycles, 0.0);
        assert_eq!(actual.cycles, expected);
    }
}
#[test]
fn invalid_geometry_frames_and_missing_sun_never_reuse_old_corrections() {
    let los = [0.0, 0.0, -2e7];
    assert_eq!(
        phase_windup_from_frames([0.0; 3], frame(0.0), RECEIVER, 5.0),
        Err(WindupError::DegenerateLineOfSight)
    );
    assert_eq!(
        phase_windup_from_frames(
            los,
            DipoleFrame {
                x: [2.0, 0.0, 0.0],
                ..frame(0.0)
            },
            RECEIVER,
            5.0
        ),
        Err(WindupError::InvalidFrame)
    );
    assert_eq!(
        phase_windup_from_frames(
            los,
            DipoleFrame {
                x: [1.0, 0.0, 0.0],
                y: [1.0, 0.0, 0.0]
            },
            RECEIVER,
            5.0
        ),
        Err(WindupError::InvalidFrame)
    );
    assert_eq!(
        phase_windup_from_frames(los, RECEIVER, RECEIVER, 5.0),
        Err(WindupError::DegenerateDipoles)
    );
    assert_eq!(
        phase_windup_from_frames(los, frame(0.0), RECEIVER, f64::NAN),
        Err(WindupError::NonFinite)
    );
    assert_eq!(
        phase_windup_from_frames([f64::MAX; 3], frame(0.0), RECEIVER, 0.0),
        Err(WindupError::NonFinite)
    );
    let receiver = EcefCoord::new(6378137.0, 0.0, 0.0);
    let satellite = EcefCoord::new(2.6e7, 0.0, 0.0);
    assert_eq!(
        nominal_yaw_windup(receiver, satellite, EcefCoord::default(), 5.0),
        Err(WindupError::MissingSun)
    );
    assert_eq!(
        nominal_yaw_windup(receiver, satellite, EcefCoord::new(1.5e11, 0.0, 0.0), 5.0),
        Err(WindupError::DegenerateAttitude)
    );
    assert!(
        nominal_yaw_windup(
            EcefCoord::default(),
            satellite,
            EcefCoord::new(0.0, 1.5e11, 0.0),
            5.0
        )
        .is_err()
    );
}
