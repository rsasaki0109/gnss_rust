use gnss_rust::coordinates::*;
use gnss_rust::{EcefCoord, EnuCoord, Error, GeodeticCoord, constants as c};
use std::f64::consts::{FRAC_PI_2, PI};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.15} != {expected:.15}, tolerance {tolerance}"
    );
}

#[test]
fn wgs84_equatorial_axes_and_exact_poles() {
    for (input, expected) in [
        (
            GeodeticCoord::new(0.0, 0.0, 0.0).unwrap(),
            [c::WGS84_A, 0.0, 0.0],
        ),
        (
            GeodeticCoord::new(0.0, FRAC_PI_2, 100.0).unwrap(),
            [0.0, c::WGS84_A + 100.0, 0.0],
        ),
        (
            GeodeticCoord::new(FRAC_PI_2, 0.0, 0.0).unwrap(),
            [0.0, 0.0, c::WGS84_B],
        ),
        (
            GeodeticCoord::new(-FRAC_PI_2, PI, 20.0).unwrap(),
            [0.0, 0.0, -c::WGS84_B - 20.0],
        ),
    ] {
        for (actual, expected) in geodetic_to_ecef(input)
            .unwrap()
            .to_array()
            .into_iter()
            .zip(expected)
        {
            near(actual, expected, 1e-8);
        }
    }
    for sign in [-1.0, 1.0] {
        let geo = ecef_to_geodetic(EcefCoord::new(0.0, 0.0, sign * (c::WGS84_B + 100.0))).unwrap();
        assert_eq!(geo.latitude, sign * FRAC_PI_2);
        assert_eq!(geo.longitude, 0.0);
        near(geo.height, 100.0, 1e-9);
    }
}

#[test]
fn synthetic_ppp_receiver_matches_independent_fixture_coordinates() {
    let position = GeodeticCoord::new(35_f64.to_radians(), 139_f64.to_radians(), 45.0).unwrap();
    let ecef = geodetic_to_ecef(position).unwrap();
    // The fixture builder offsets APPROX POSITION XYZ by (+8, -5, +3) m
    // from truth. Undo that offset; the RINEX header is rounded to 0.1 mm.
    for (actual, expected) in
        ecef.to_array()
            .into_iter()
            .zip([-3_947_481.064_9, 3_431_492.937_5, 3_637_892.720_3])
    {
        near(actual, expected, 5e-5);
    }
}

#[test]
fn coordinate_round_trips_across_hemispheres_poles_and_orbital_heights() {
    for latitude in [-90_f64, -89.999_999, -45.0, 0.0, 35.0, 89.999_999, 90.0] {
        for longitude in [-180_f64, -139.0, -1.0, 0.0, 139.0, 180.0] {
            for height in [-500.0, 0.0, 45.0, 10_000.0, 20_200_000.0] {
                let input =
                    GeodeticCoord::new(latitude.to_radians(), longitude.to_radians(), height)
                        .unwrap();
                let ecef = geodetic_to_ecef(input).unwrap();
                let geo = ecef_to_geodetic(ecef).unwrap();
                near(geo.latitude, input.latitude, 1e-12);
                near(geo.height, height, 1e-7);
                for (actual, expected) in geodetic_to_ecef(geo)
                    .unwrap()
                    .to_array()
                    .into_iter()
                    .zip(ecef.to_array())
                {
                    near(actual, expected, 1e-7);
                }
            }
        }
    }
}

#[test]
fn invalid_geodetic_inputs_are_checked_even_with_public_fields() {
    assert_eq!(
        GeodeticCoord::new(FRAC_PI_2 + 0.01, 0.0, 0.0),
        Err(Error::InvalidLatitude)
    );
    for input in [
        GeodeticCoord {
            latitude: -PI,
            longitude: 0.0,
            height: 0.0,
        },
        GeodeticCoord {
            latitude: 0.0,
            longitude: f64::NAN,
            height: 0.0,
        },
        GeodeticCoord {
            latitude: f64::INFINITY,
            longitude: 0.0,
            height: 0.0,
        },
        GeodeticCoord {
            latitude: 0.0,
            longitude: 0.0,
            height: f64::NEG_INFINITY,
        },
    ] {
        assert!(geodetic_to_ecef(input).is_err());
        assert!(ecef_to_enu(EcefCoord::default(), input).is_err());
        assert!(enu_to_ecef(EnuCoord::default(), input).is_err());
    }
}

#[test]
fn earth_centre_and_nonfinite_ecef_are_not_valid_geodetic_positions() {
    assert_eq!(
        ecef_to_geodetic(EcefCoord::default()),
        Err(Error::UndefinedGeodeticOrigin)
    );
    for input in [
        EcefCoord::new(f64::NAN, 1.0, 1.0),
        EcefCoord::new(0.0, f64::INFINITY, 0.0),
        EcefCoord::new(0.0, 0.0, f64::NEG_INFINITY),
        EcefCoord::new(f64::MAX, f64::MAX, f64::MAX),
    ] {
        assert_eq!(ecef_to_geodetic(input), Err(Error::NonFinite));
    }
}

#[test]
fn enu_at_equator_has_east_y_north_z_up_x() {
    let origin = GeodeticCoord::default();
    let enu = ecef_to_enu(EcefCoord::new(1.0, 2.0, 3.0), origin).unwrap();
    assert_eq!(enu, EnuCoord::new(2.0, 3.0, 1.0));
    assert_eq!(
        enu_to_ecef(enu, origin).unwrap(),
        EcefCoord::new(1.0, 2.0, 3.0)
    );
}

#[test]
fn enu_inverse_preserves_vector_and_length_across_origin_frames() {
    let difference = EcefCoord::new(1200.0, -300.0, 45.0);
    for latitude in [-90_f64, -45.0, 0.0, 35.0, 90.0] {
        for longitude in [-180_f64, -90.0, 0.0, 139.0, 180.0] {
            let origin =
                GeodeticCoord::new(latitude.to_radians(), longitude.to_radians(), 100.0).unwrap();
            let enu = ecef_to_enu(difference, origin).unwrap();
            near(
                enu.east.hypot(enu.north).hypot(enu.up),
                difference.x.hypot(difference.y).hypot(difference.z),
                1e-12,
            );
            for (actual, expected) in enu_to_ecef(enu, origin)
                .unwrap()
                .to_array()
                .into_iter()
                .zip(difference.to_array())
            {
                near(actual, expected, 1e-12);
            }
        }
    }
}

#[test]
fn enu_rejects_nonfinite_vectors() {
    let origin = GeodeticCoord::default();
    assert_eq!(
        ecef_to_enu(EcefCoord::new(0.0, f64::NAN, 1.0), origin),
        Err(Error::NonFinite)
    );
    assert_eq!(
        enu_to_ecef(EnuCoord::new(f64::INFINITY, 0.0, 1.0), origin),
        Err(Error::NonFinite)
    );
}

#[test]
fn geometric_range_matches_analytic_sagnac_and_sign_convention() {
    let satellite = EcefCoord::new(20_000_000.0, 0.0, 0.0);
    let receiver = EcefCoord::new(0.0, 6_000_000.0, 0.0);
    let euclidean = 20_000_000_f64.hypot(6_000_000.0);
    let forward = geometric_distance(satellite, receiver).unwrap();
    let reverse = geometric_distance(receiver, satellite).unwrap();
    // Independently evaluated first-order Sagnac correction for these axes.
    near(forward - euclidean, 29.188_653_491_876_704, 1e-8);
    near(reverse - euclidean, -29.188_653_491_876_704, 1e-8);
    near((forward + reverse) / 2.0, euclidean, 1e-8);
    assert_eq!(geometric_distance(receiver, receiver).unwrap(), 0.0);
}

#[test]
fn geometric_range_rejects_nonfinite_and_overflow() {
    assert_eq!(
        geometric_distance(EcefCoord::new(f64::NAN, 0.0, 0.0), EcefCoord::default()),
        Err(Error::NonFinite)
    );
    assert_eq!(
        geometric_distance(
            EcefCoord::new(f64::MAX, f64::MAX, 0.0),
            EcefCoord::new(-f64::MAX, 1.0, 0.0)
        ),
        Err(Error::NonFinite)
    );
}
