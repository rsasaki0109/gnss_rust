use gnss_rust::iers2010::*;
use gnss_rust::{EcefCoord, constants};

fn vector(values: &[f64]) -> EcefCoord {
    EcefCoord::new(values[0], values[1], values[2])
}
fn close(a: EcefCoord, b: EcefCoord, tolerance: f64) {
    for (x, y) in [a.x, a.y, a.z].iter().zip([b.x, b.y, b.z]) {
        assert!(
            (x - y).abs() < tolerance,
            "{x:.17} != {y:.17} (tolerance {tolerance})"
        );
    }
}
fn length(v: EcefCoord) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn difference(a: EcefCoord, b: EcefCoord) -> EcefCoord {
    EcefCoord::new(a.x - b.x, a.y - b.y, a.z - b.z)
}
fn reference_inputs() -> (EcefCoord, EarthFixedSunMoon, IersTideArguments) {
    (
        EcefCoord::new(4075578.385, 931852.890, 4801570.154),
        EarthFixedSunMoon::new(
            EcefCoord::new(137859926952.015, 54228127881.4350, 23509422341.6960),
            EcefCoord::new(-179996231.920342, -312468450.131567, -169288918.592160),
        )
        .unwrap(),
        IersTideArguments::new((54934.0 + 66.184 / 86400.0 - 51544.5) / 36525.0, 0.0).unwrap(),
    )
}

#[test]
fn published_dehant_reference_case_matches_without_permanent_tide_step3() {
    let (station, bodies, args) = reference_inputs();
    let result = solid_earth_tide(station, bodies, args).unwrap();
    close(
        result.displacement_m,
        EcefCoord::new(
            0.07700420357108126,
            0.06304056321824968,
            0.05516568152597247,
        ),
        1e-6,
    );
    assert_eq!(result.nominal_marker_m, station);
    assert_eq!(result.bodies, bodies);
    assert_eq!(result.arguments, args);
    assert_eq!(
        result.instantaneous_marker_m,
        EcefCoord::new(
            station.x + result.displacement_m.x,
            station.y + result.displacement_m.y,
            station.z + result.displacement_m.z
        )
    );
    assert!(length(result.step1_diurnal_m) > 1e-4);
    assert!(length(result.step2_diurnal_m) > 1e-3);
}

#[test]
fn actual_native_kernel_and_each_step_match_505_independent_cases() {
    let mut count = 0;
    for line in include_str!("fixtures/upstream_iers2010_kernel.csv")
        .lines()
        .skip(1)
    {
        let v: Vec<f64> = line.split(',').map(|x| x.parse().unwrap()).collect();
        let station = vector(&v[3..6]);
        let bodies = EarthFixedSunMoon::new(vector(&v[6..9]), vector(&v[9..12])).unwrap();
        let result =
            solid_earth_tide(station, bodies, IersTideArguments::new(v[1], v[2]).unwrap()).unwrap();
        for (actual, reference) in [
            result.degree2_degree3_m,
            result.step1_diurnal_m,
            result.step1_semidiurnal_m,
            result.step1_latitude_m,
            result.step2_diurnal_m,
            result.step2_long_period_m,
            result.displacement_m,
        ]
        .into_iter()
        .zip(v[12..33].chunks_exact(3))
        {
            close(actual, vector(reference), 2e-12);
        }
        close(result.instantaneous_marker_m, vector(&v[33..36]), 1e-8);
        count += 1;
    }
    assert_eq!(count, 505);
}

#[test]
fn native_sofa_body_inputs_reproduce_iers_wrapper_for_every_hour_and_eop_case() {
    let mut count = 0;
    let mut first = None;
    let mut previous_day = None;
    for line in include_str!("fixtures/upstream_iers2010_epochs.csv")
        .lines()
        .skip(1)
    {
        let v: Vec<f64> = line.split(',').map(|x| x.parse().unwrap()).collect();
        let result = solid_earth_tide(
            vector(&v[7..10]),
            EarthFixedSunMoon::new(vector(&v[10..13]), vector(&v[13..16])).unwrap(),
            IersTideArguments::new(v[2], v[3]).unwrap(),
        )
        .unwrap();
        close(result.displacement_m, vector(&v[16..19]), 2e-12);
        close(result.instantaneous_marker_m, vector(&v[19..22]), 1e-8);
        assert!(length(result.displacement_m) < 0.5);
        if count % 50 == 0 {
            first = Some(result.displacement_m);
        }
        if count % 50 == 24 {
            assert!(length(difference(result.displacement_m, first.unwrap())) > 0.03);
        }
        if count % 50 == 48 {
            previous_day = Some(v[1].floor());
        }
        if count % 50 == 49 {
            assert_eq!(Some(v[1].floor()), previous_day);
            assert_eq!(v[3], 0.0);
        }
        count += 1;
    }
    assert_eq!(count, 150);
}

#[test]
fn tt_and_ut_arguments_affect_frequency_terms_without_mutating_body_terms() {
    let (station, bodies, args) = reference_inputs();
    let zero = solid_earth_tide(station, bodies, args).unwrap();
    let hour = solid_earth_tide(
        station,
        bodies,
        IersTideArguments::new(args.julian_centuries_tt(), 6.0).unwrap(),
    )
    .unwrap();
    assert_eq!(zero.degree2_degree3_m, hour.degree2_degree3_m);
    assert_eq!(zero.step1_diurnal_m, hour.step1_diurnal_m);
    assert_eq!(zero.step1_semidiurnal_m, hour.step1_semidiurnal_m);
    assert_eq!(zero.step1_latitude_m, hour.step1_latitude_m);
    assert_eq!(zero.step2_long_period_m, hour.step2_long_period_m);
    assert!(length(difference(zero.step2_diurnal_m, hour.step2_diurnal_m)) > 0.001);
    let year = solid_earth_tide(
        station,
        bodies,
        IersTideArguments::new(args.julian_centuries_tt() + 0.01, 0.0).unwrap(),
    )
    .unwrap();
    assert!(
        length(difference(
            zero.step2_long_period_m,
            year.step2_long_period_m
        )) > 1e-5
    );
}

#[test]
fn step2_is_continuous_at_ut_midnight_and_includes_more_than_legacy_love_terms() {
    let (station, bodies, args) = reference_inputs();
    let a = solid_earth_tide(
        station,
        bodies,
        IersTideArguments::new(
            args.julian_centuries_tt() - 0.001 / 86400.0 / 36525.0,
            24.0 - 0.001 / 3600.0,
        )
        .unwrap(),
    )
    .unwrap();
    let b = solid_earth_tide(
        station,
        bodies,
        IersTideArguments::new(
            args.julian_centuries_tt() + 0.001 / 86400.0 / 36525.0,
            0.001 / 3600.0,
        )
        .unwrap(),
    )
    .unwrap();
    close(a.displacement_m, b.displacement_m, 1e-7);
    let legacy_sun =
        gnss_rust::tides::legacy_body_tide_displacement(station, bodies.sun_m(), 1.32712440018e20)
            .unwrap();
    let legacy_moon =
        gnss_rust::tides::legacy_body_tide_displacement(station, bodies.moon_m(), 4.902801e12)
            .unwrap();
    let legacy = EcefCoord::new(
        legacy_sun.x + legacy_moon.x,
        legacy_sun.y + legacy_moon.y,
        legacy_sun.z + legacy_moon.z,
    );
    assert!(length(difference(b.displacement_m, legacy)) > 0.001);
}

#[test]
fn invalid_body_station_and_axis_fail_explicitly_instead_of_returning_zero() {
    let (station, bodies, args) = reference_inputs();
    for invalid in [EcefCoord::new(0.0, 0.0, 0.0), EcefCoord::new(1.0, 0.0, 0.0)] {
        assert_eq!(
            EarthFixedSunMoon::new(invalid, bodies.moon_m()),
            Err(IersTideError::InvalidBody)
        );
        assert_eq!(
            EarthFixedSunMoon::new(bodies.sun_m(), invalid),
            Err(IersTideError::InvalidBody)
        );
    }
    for invalid in [EcefCoord::new(0.0, 0.0, 0.0), EcefCoord::new(1e6, 0.0, 0.0)] {
        assert_eq!(
            solid_earth_tide(invalid, bodies, args),
            Err(IersTideError::InvalidStation)
        );
    }
    for z in [-constants::WGS84_A, constants::WGS84_A] {
        assert_eq!(
            solid_earth_tide(EcefCoord::new(0.0, 0.0, z), bodies, args),
            Err(IersTideError::AxialStation)
        );
    }
    for value in [f64::NAN, f64::INFINITY, f64::MAX] {
        let bad = EcefCoord::new(value, station.y, station.z);
        assert_eq!(
            solid_earth_tide(bad, bodies, args),
            Err(IersTideError::NonFinite)
        );
        assert_eq!(
            EarthFixedSunMoon::new(bad, bodies.moon_m()),
            Err(IersTideError::NonFinite)
        );
    }
}

#[test]
fn invalid_time_arguments_are_not_wrapped_or_treated_as_gps_labels() {
    for ut in [-0.01, 24.0, 25.0] {
        assert_eq!(
            IersTideArguments::new(0.0, ut),
            Err(IersTideError::InvalidArguments)
        );
    }
    for (tt, ut) in [(f64::NAN, 0.0), (f64::INFINITY, 0.0), (0.0, f64::NAN)] {
        assert_eq!(
            IersTideArguments::new(tt, ut),
            Err(IersTideError::NonFinite)
        );
    }
    let (station, bodies, _) = reference_inputs();
    assert_eq!(
        solid_earth_tide(
            station,
            bodies,
            IersTideArguments::new(f64::MAX, 0.0).unwrap()
        ),
        Err(IersTideError::NonFinite)
    );
}
