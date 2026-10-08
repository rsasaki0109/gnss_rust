use gnss_rust::coordinates::{enu_to_ecef, geodetic_to_ecef};
use gnss_rust::tides::*;
use gnss_rust::{EcefCoord, EnuCoord, GeodeticCoord, GnssTime, constants as c};

fn v(n: &[f64]) -> EcefCoord {
    EcefCoord::new(n[0], n[1], n[2])
}
fn close(a: EcefCoord, b: EcefCoord, tolerance: f64) {
    for (x, y) in [a.x, a.y, a.z].iter().zip([b.x, b.y, b.z]) {
        assert!((x - y).abs() < tolerance, "{x} != {y}");
    }
}
#[test]
fn actual_native_moon_helper_matches_fractional_epochs_and_week_boundaries() {
    let mut count = 0;
    for line in include_str!("fixtures/upstream_legacy_moon.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = line.split(',').collect();
        let n: Vec<f64> = f[2..].iter().map(|x| x.parse().unwrap()).collect();
        close(
            approximate_moon_position_ecef(
                GnssTime::new(f[0].parse().unwrap(), f[1].parse().unwrap()).unwrap(),
            )
            .unwrap(),
            v(&n),
            1e-6,
        );
        count += 1;
    }
    assert_eq!(count, 12);
}
#[test]
fn actual_native_legacy_body_tides_match_sun_moon_and_total_displacement_at_252_stations() {
    let mut count = 0;
    for line in include_str!("fixtures/upstream_legacy_solid_tide.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = line.split(',').collect();
        let n: Vec<f64> = f[2..].iter().map(|x| x.parse().unwrap()).collect();
        let time = GnssTime::new(f[0].parse().unwrap(), f[1].parse().unwrap()).unwrap();
        let rx = v(&n[..3]);
        let sun = v(&n[3..6]);
        let moon = v(&n[6..9]);
        close(
            legacy_body_tide_displacement(rx, sun, 1.32712440018e20).unwrap(),
            v(&n[9..12]),
            2e-12,
        );
        close(
            legacy_body_tide_displacement(rx, moon, 4.902801e12).unwrap(),
            v(&n[12..15]),
            2e-12,
        );
        let tide =
            solid_earth_tide(time, rx, SolidEarthTideModel::LegacyLoveApproximateSunMoon).unwrap();
        assert_eq!(tide.nominal_marker_m, rx);
        assert_eq!(tide.time, time);
        close(tide.displacement_m, v(&n[15..18]), 2e-12);
        close(tide.instantaneous_marker_m, v(&n[18..21]), 1e-8);
        assert_eq!(
            tide.instantaneous_marker_m,
            EcefCoord::new(
                rx.x + tide.displacement_m.x,
                rx.y + tide.displacement_m.y,
                rx.z + tide.displacement_m.z
            )
        );
        count += 1;
    }
    assert_eq!(count, 252);
}
#[test]
fn explicit_body_love_terms_have_correct_radial_sign_and_inverse_cube_scaling() {
    let rx = EcefCoord::new(c::WGS84_A, 0.0, 0.0);
    let d = 1e11;
    let gm = 1.32712440018e20;
    let scale = gm / 3.986004418e14 * (c::WGS84_A / d).powi(3) * c::WGS84_A;
    let toward = legacy_body_tide_displacement(rx, EcefCoord::new(d, 0.0, 0.0), gm).unwrap();
    close(toward, EcefCoord::new(scale * 0.6078, 0.0, 0.0), 1e-15);
    let away = legacy_body_tide_displacement(rx, EcefCoord::new(-d, 0.0, 0.0), gm).unwrap();
    close(toward, away, 1e-15);
    let perpendicular = legacy_body_tide_displacement(rx, EcefCoord::new(0.0, d, 0.0), gm).unwrap();
    close(
        perpendicular,
        EcefCoord::new(-scale * 0.6078 * 0.5, 0.0, 0.0),
        1e-15,
    );
    let far = legacy_body_tide_displacement(rx, EcefCoord::new(2.0 * d, 0.0, 0.0), gm).unwrap();
    close(far, EcefCoord::new(toward.x / 8.0, 0.0, 0.0), 1e-15);
    let diagonal = legacy_body_tide_displacement(
        rx,
        EcefCoord::new(d / 2_f64.sqrt(), d / 2_f64.sqrt(), 0.0),
        gm,
    )
    .unwrap();
    assert!(diagonal.x > 0.0 && diagonal.y > 0.0);
}
#[test]
fn legacy_body_uses_geodetic_up_and_is_continuous_across_gps_week_rollover() {
    let geo = GeodeticCoord::new(45_f64.to_radians(), 139_f64.to_radians(), 1500.0).unwrap();
    let rx = geodetic_to_ecef(geo).unwrap();
    let up = enu_to_ecef(
        EnuCoord {
            east: 0.0,
            north: 0.0,
            up: 1e11,
        },
        geo,
    )
    .unwrap();
    let tide = legacy_body_tide_displacement(rx, up, 1.32712440018e20).unwrap();
    let normal = enu_to_ecef(
        EnuCoord {
            east: 0.0,
            north: 0.0,
            up: tide.x.hypot(tide.y).hypot(tide.z),
        },
        geo,
    )
    .unwrap();
    close(tide, normal, 1e-15);
    let a = solid_earth_tide(
        GnssTime::new(2300, 604799.999).unwrap(),
        rx,
        SolidEarthTideModel::LegacyLoveApproximateSunMoon,
    )
    .unwrap();
    let b = solid_earth_tide(
        GnssTime::new(2301, 0.001).unwrap(),
        rx,
        SolidEarthTideModel::LegacyLoveApproximateSunMoon,
    )
    .unwrap();
    close(a.displacement_m, b.displacement_m, 1e-6);
}
#[test]
fn invalid_body_station_gravity_and_calendar_never_become_zero_tides() {
    let rx = EcefCoord::new(c::WGS84_A, 0.0, 0.0);
    let body = EcefCoord::new(1e11, 0.0, 0.0);
    assert_eq!(
        legacy_body_tide_displacement(rx, EcefCoord::new(0.0, 0.0, 0.0), 1e20),
        Err(TideError::InvalidBody)
    );
    for gm in [0.0, -1.0] {
        assert_eq!(
            legacy_body_tide_displacement(rx, body, gm),
            Err(TideError::InvalidBodyGravity)
        );
    }
    assert_eq!(
        legacy_body_tide_displacement(rx, body, f64::NAN),
        Err(TideError::NonFinite)
    );
    for rx in [EcefCoord::new(0.0, 0.0, 0.0), EcefCoord::new(1e6, 0.0, 0.0)] {
        assert_eq!(
            solid_earth_tide(
                GnssTime::default(),
                rx,
                SolidEarthTideModel::LegacyLoveApproximateSunMoon
            ),
            Err(TideError::InvalidStation)
        );
    }
    assert_eq!(
        legacy_body_tide_displacement(rx, EcefCoord::new(f64::INFINITY, 0.0, 0.0), 1e20),
        Err(TideError::NonFinite)
    );
    assert!(approximate_moon_position_ecef(GnssTime::new(i32::MAX, 0.0).unwrap()).is_err());
}
