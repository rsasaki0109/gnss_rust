use gnss_rust::GnssTime;
use gnss_rust::earth_rotation::*;
use gnss_rust::time::{CalendarDateTime, TimeScale};

fn table() -> LeapSecondTable {
    LeapSecondTable::historical("frozen-reference-through-2017", 61000).unwrap()
}
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() < t, "{a:.17} != {b:.17}, tolerance={t}");
}
fn identity() -> Matrix3 {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
fn matrix(v: &[f64]) -> Matrix3 {
    std::array::from_fn(|i| std::array::from_fn(|j| v[3 * i + j]))
}
fn apply(m: Matrix3, v: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| (0..3).map(|j| m[i][j] * v[j]).sum())
}
fn midnight_2017() -> GnssTime {
    CalendarDateTime {
        year: 2017,
        month: 1,
        day: 1,
        hour: 0,
        minute: 0,
        second: 0.0,
    }
    .to_gnss_time(TimeScale::Utc)
    .unwrap()
}

#[test]
fn actual_sofa_times_era_polar_and_composition_match_all_reference_rows() {
    let leaps = table();
    let mut count = 0;
    let mut native_mismatches = 0;
    let mut largest_native: f64 = 0.0;
    for line in include_str!("fixtures/upstream_earth_rotation.csv")
        .lines()
        .skip(1)
    {
        let v: Vec<f64> = line.split(',').map(|x| x.parse().unwrap()).collect();
        let gps = GnssTime::new(v[1] as i32, v[2]).unwrap();
        let eop = EarthOrientationParams::new(v[13], v[14], v[15]).unwrap();
        let times = RotationTimeScales::from_gps(gps, &leaps, eop).unwrap();
        assert_eq!(times.utc.day_mjd(), v[3] as i32);
        close(times.utc.seconds_of_day(), v[4], 3e-9);
        assert_eq!(times.utc.day_length_s(), v[5] as u32);
        assert_eq!(times.utc.tai_minus_utc_s(), v[6] as i32);
        close(times.utc.quasi_mjd(), v[7], 1e-11);
        assert_eq!(times.tt.day_mjd(), v[9] as i32);
        close(times.tt.seconds_of_day(), v[10], 3e-9);
        assert_eq!(times.ut1.day_mjd(), v[11] as i32);
        close(times.ut1.seconds_of_day(), v[12], 3e-9);
        close(earth_rotation_angle(times.ut1), v[16], 2e-12);
        let pom = polar_motion_matrix(times.tt, eop).unwrap();
        let composed = compose_celestial_to_terrestrial(&times, matrix(&v[18..27])).unwrap();
        for (actual, reference) in [(pom, &v[27..36]), (composed, &v[36..45])] {
            for (a, b) in actual.iter().flatten().zip(reference) {
                close(*a, *b, 2e-12);
            }
        }
        let args = times.tide_arguments_utc_approx().unwrap();
        close(args.julian_centuries_tt(), v[45], 1e-16);
        close(args.fractional_hour_ut(), v[46], 1e-12);
        assert_eq!(times.leap_product_id, leaps.product_id());
        let difference = (v[8] - v[7]).abs() * 86400.0;
        largest_native = largest_native.max(difference);
        if v[5] == 86400.0 {
            close(times.utc.quasi_mjd(), v[8], 1e-11);
        } else if difference > 1e-6 {
            native_mismatches += 1;
        }
        count += 1;
    }
    eprintln!(
        "native one-pass UTC vs SOFA quasi-UTC: {native_mismatches} differing rows; largest equivalent-seconds difference={largest_native:.9}"
    );
    assert_eq!(count, 558);
    assert!(native_mismatches > 100);
    assert!(largest_native > 0.49);
}

#[test]
fn positive_leap_second_is_explicit_and_tt_ut1_remain_continuous() {
    let gps = midnight_2017();
    let leaps = table();
    let mut previous: Option<RotationTimeScales> = None;
    for delta in [-1.5, -1.0, -0.5, 0.0, 0.5] {
        let eop =
            EarthOrientationParams::new(if delta < 0.0 { -0.4 } else { 0.6 }, 0.1, -0.2).unwrap();
        let times =
            RotationTimeScales::from_gps(gps.checked_add_seconds(delta).unwrap(), &leaps, eop)
                .unwrap();
        assert_eq!(times.utc.is_leap_second(), (-1.0..0.0).contains(&delta));
        if delta < 0.0 {
            assert_eq!(times.utc.day_mjd(), 57753);
            assert_eq!(times.utc.day_length_s(), 86401);
            assert_eq!(times.utc.tai_minus_utc_s(), 36);
        } else {
            assert_eq!(times.utc.day_mjd(), 57754);
            assert_eq!(times.utc.day_length_s(), 86400);
            assert_eq!(times.utc.tai_minus_utc_s(), 37);
        }
        if let Some(prev) = previous {
            close(times.tt.difference_seconds(prev.tt), 0.5, 1e-10);
            close(times.ut1.difference_seconds(prev.ut1), 0.5, 1e-10);
            close(
                (earth_rotation_angle(times.ut1) - earth_rotation_angle(prev.ut1))
                    .rem_euclid(std::f64::consts::TAU),
                0.5 * 7.29211514670698e-5,
                2e-12,
            );
        }
        previous = Some(times);
    }
}

#[test]
fn coverage_is_half_open_and_future_leap_updates_are_explicit() {
    let baseline = table();
    let eop = EarthOrientationParams::new(0.0, 0.0, 0.0).unwrap();
    assert_eq!(
        RotationTimeScales::from_gps(GnssTime::new(0, -0.1).unwrap(), &baseline, eop),
        Err(RotationError::LeapCoverage)
    );
    let endgps = GnssTime::new(
        (61000 - 44244) / 7,
        f64::from((61000 - 44244) % 7) * 86400.0 + 18.0,
    )
    .unwrap();
    assert!(
        RotationTimeScales::from_gps(endgps.checked_add_seconds(-0.1).unwrap(), &baseline, eop)
            .is_ok()
    );
    assert_eq!(
        RotationTimeScales::from_gps(endgps, &baseline, eop),
        Err(RotationError::LeapCoverage)
    );
    let mut entries = baseline.transitions().to_vec();
    entries.push(LeapSecondTransition {
        utc_day_mjd: 60000,
        tai_minus_utc_s: 38,
    });
    let updated = LeapSecondTable::new("synthetic-future-leap", entries, 61000).unwrap();
    let next = GnssTime::new(
        (60000 - 44244) / 7,
        f64::from((60000 - 44244) % 7) * 86400.0 + 19.0,
    )
    .unwrap();
    let before =
        RotationTimeScales::from_gps(next.checked_add_seconds(-0.5).unwrap(), &updated, eop)
            .unwrap();
    let after = RotationTimeScales::from_gps(next, &updated, eop).unwrap();
    assert!(before.utc.is_leap_second());
    assert_eq!(before.utc.day_mjd(), 59999);
    assert_eq!(before.utc.tai_minus_utc_s(), 37);
    assert_eq!(after.utc.day_mjd(), 60000);
    assert_eq!(after.utc.tai_minus_utc_s(), 38);
    assert_eq!(after.leap_product_id, "synthetic-future-leap");
}

#[test]
fn invalid_tables_offsets_rotations_and_numeric_ranges_are_rejected() {
    let good = table();
    let first = good.transitions()[0];
    for entries in [
        vec![],
        vec![LeapSecondTransition {
            utc_day_mjd: 44244,
            tai_minus_utc_s: 18,
        }],
        vec![first, first],
        vec![
            first,
            LeapSecondTransition {
                utc_day_mjd: 45000,
                tai_minus_utc_s: 18,
            },
        ],
        vec![
            first,
            LeapSecondTransition {
                utc_day_mjd: 45000,
                tai_minus_utc_s: 21,
            },
        ],
    ] {
        assert_eq!(
            LeapSecondTable::new("invalid", entries, 61000),
            Err(RotationError::InvalidLeapTable)
        );
    }
    assert_eq!(
        LeapSecondTable::new("", vec![first], 61000),
        Err(RotationError::InvalidLeapTable)
    );
    assert_eq!(
        LeapSecondTable::new("bad-end", vec![first], 44244),
        Err(RotationError::InvalidLeapTable)
    );
    assert_eq!(UniformMjd::new(1, f64::NAN), Err(RotationError::NonFinite));
    assert_eq!(UniformMjd::new(1, f64::MAX), Err(RotationError::OutOfRange));
    assert_eq!(
        UniformMjd::new(i32::MAX, 86400.0),
        Err(RotationError::OutOfRange)
    );
    assert_eq!(
        EarthOrientationParams::new(f64::NAN, 0.0, 0.0),
        Err(RotationError::NonFinite)
    );
    let eop = EarthOrientationParams::new(0.0, 0.0, 0.0).unwrap();
    assert_eq!(
        RotationTimeScales::from_gps(GnssTime::new(i32::MAX, 0.0).unwrap(), &good, eop),
        Err(RotationError::OutOfRange)
    );
    let times = RotationTimeScales::from_gps(midnight_2017(), &good, eop).unwrap();
    let mut invalid = identity();
    invalid[0][0] = f64::NAN;
    assert_eq!(
        compose_celestial_to_terrestrial(&times, invalid),
        Err(RotationError::NonFinite)
    );
    invalid = identity();
    invalid[0][0] = -1.0;
    assert_eq!(
        compose_celestial_to_terrestrial(&times, invalid),
        Err(RotationError::InvalidRotation)
    );
    invalid = identity();
    invalid[0][1] = 0.001;
    assert_eq!(
        compose_celestial_to_terrestrial(&times, invalid),
        Err(RotationError::InvalidRotation)
    );
}

#[test]
fn matrix_direction_inverse_and_units_preserve_earth_fixed_geometry() {
    let leaps = table();
    let times = RotationTimeScales::from_gps(
        midnight_2017(),
        &leaps,
        EarthOrientationParams::new(0.3345, 0.1, -0.2).unwrap(),
    )
    .unwrap();
    let rotation = compose_celestial_to_terrestrial(&times, identity()).unwrap();
    let input = [137859926952.015, 54228127881.435, 23509422341.696];
    let terrestrial = apply(rotation, input);
    let returned = apply(transpose(rotation), terrestrial);
    for (a, b) in returned.iter().zip(input) {
        close(*a, b, 5e-5);
    }
    let zero = compose_celestial_to_terrestrial(&times, identity()).unwrap();
    let angle = earth_rotation_angle(times.ut1);
    let got = apply(zero, [1.0, 0.0, 0.0]);
    close(got[0], angle.cos(), 2e-6);
    close(got[1], -angle.sin(), 2e-6);
    let polar = polar_motion_matrix(times.tt, times.eop).unwrap();
    assert!(polar[0][2] > 0.0);
    assert!(polar[2][0] < 0.0);
    close(polar[0][2], 0.1_f64.to_radians() / 3600.0, 1e-15);
}

#[test]
fn fractional_gps_week_rollover_keeps_tt_offset_and_ut1_small_differences() {
    let leaps = table();
    let eop = EarthOrientationParams::new(-0.1234, 0.0, 0.0).unwrap();
    let a = RotationTimeScales::from_gps(GnssTime::new(2300, 604799.999).unwrap(), &leaps, eop)
        .unwrap();
    let b = RotationTimeScales::from_gps(GnssTime::new(2301, 0.001).unwrap(), &leaps, eop).unwrap();
    close(b.tt.difference_seconds(a.tt), 0.002, 1e-10);
    close(b.ut1.difference_seconds(a.ut1), 0.002, 1e-10);
    close(a.tt.difference_seconds(a.tai), 32.184, 1e-10);
    close(a.tai.difference_seconds(a.ut1), 37.1234, 1e-10);
}
