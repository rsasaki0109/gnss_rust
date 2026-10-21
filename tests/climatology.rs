use gnss_rust::GeodeticCoord;
use gnss_rust::models::zenith_troposphere_climatology;

#[test]
fn climatology_matches_all_312_actual_native_latitude_season_height_cases() {
    let mut count = 0;
    for line in include_str!("fixtures/upstream_zenith_climatology.csv")
        .lines()
        .skip(1)
    {
        let v: Vec<f64> = line.split(',').map(|s| s.parse().unwrap()).collect();
        let p = GeodeticCoord::new(v[0].to_radians(), 0.0, v[1]).unwrap();
        let d = zenith_troposphere_climatology(p, v[2] as u16).unwrap();
        for (actual, expected) in [
            d.pressure_hpa,
            d.temperature_k,
            d.water_vapor_pressure_hpa,
            d.hydrostatic_delay_m,
            d.wet_delay_m,
            d.total_delay_m(),
        ]
        .into_iter()
        .zip(&v[3..])
        {
            assert!(
                (actual - expected).abs() < 1e-10,
                "{line}: {actual} != {expected}"
            );
        }
        count += 1;
    }
    assert_eq!(count, 312);
}

#[test]
fn negative_height_clamp_and_seasonal_hemisphere_and_latitude_boundaries_are_explicit() {
    let p = |latitude: f64, height| GeodeticCoord::new(latitude.to_radians(), 0.0, height).unwrap();
    for latitude in [-90., -75., -15., 0., 15., 75., 90.] {
        assert_eq!(
            zenith_troposphere_climatology(p(latitude, -500.), 216).unwrap(),
            zenith_troposphere_climatology(p(latitude, 0.), 216).unwrap()
        );
    }
    let north = zenith_troposphere_climatology(p(45., 0.), 28).unwrap();
    let south = zenith_troposphere_climatology(p(-45., 0.), 28).unwrap();
    assert!(south.temperature_k > north.temperature_k);
    assert!(south.wet_delay_m > north.wet_delay_m);
    assert!(
        zenith_troposphere_climatology(p(45., 10000.), 28)
            .unwrap()
            .total_delay_m()
            < north.total_delay_m()
    );
    for boundary in [15., 30., 45., 60., 75.] {
        let a = zenith_troposphere_climatology(p(boundary - 1e-8, 45.), 216)
            .unwrap()
            .total_delay_m();
        let b = zenith_troposphere_climatology(p(boundary + 1e-8, 45.), 216)
            .unwrap()
            .total_delay_m();
        assert!((a - b).abs() < 1e-8);
    }
}

#[test]
fn invalid_climatology_input_does_not_silently_return_a_constant_or_zero_delay() {
    let p = GeodeticCoord::new(0.0, 0.0, 0.0).unwrap();
    for day in [0, 367, u16::MAX] {
        assert!(zenith_troposphere_climatology(p, day).is_err());
    }
    for bad in [
        GeodeticCoord {
            latitude: f64::NAN,
            ..p
        },
        GeodeticCoord {
            latitude: std::f64::consts::PI,
            ..p
        },
        GeodeticCoord {
            height: f64::INFINITY,
            ..p
        },
        GeodeticCoord { height: 1e9, ..p },
    ] {
        assert!(zenith_troposphere_climatology(bad, 216).is_err());
    }
}
