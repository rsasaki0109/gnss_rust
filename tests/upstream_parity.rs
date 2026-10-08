use gnss_rust::coordinates::*;
use gnss_rust::{EcefCoord, EnuCoord, GeodeticCoord};

fn near(actual: f64, expected: f64, tolerance: f64, row: usize, column: &str) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "reference row {row}, {column}: {actual:.17} != {expected:.17} (tolerance {tolerance})"
    );
}

#[test]
fn agrees_with_unmodified_upstream_cpp_coordinate_header() {
    let rows: Vec<_> = include_str!("fixtures/upstream_coordinates.csv")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(rows.len(), 8, "the reference test must execute all cases");
    for (row, line) in rows.iter().enumerate() {
        let values: Vec<f64> = line.split(',').map(|s| s.parse().unwrap()).collect();
        assert_eq!(values.len(), 16);
        let geo = GeodeticCoord::new(values[0], values[1], values[2]).unwrap();
        let ecef = geodetic_to_ecef(geo).unwrap();
        for i in 0..3 {
            near(ecef.to_array()[i], values[3 + i], 1e-7, row, "ECEF");
        }
        // Supply the exact C++ ECEF vector to isolate the inverse conversion.
        let reference_ecef = EcefCoord::new(values[3], values[4], values[5]);
        let inverse = ecef_to_geodetic(reference_ecef).unwrap();
        near(inverse.latitude, values[6], 1e-12, row, "latitude");
        near(inverse.longitude, values[7], 1e-12, row, "longitude");
        near(inverse.height, values[8], 1e-5, row, "height");
        let enu = ecef_to_enu(EcefCoord::new(100.0, -200.0, 300.0), geo).unwrap();
        for (i, actual) in [enu.east, enu.north, enu.up].into_iter().enumerate() {
            near(actual, values[9 + i], 1e-12, row, "ENU");
        }
        let delta = enu_to_ecef(EnuCoord::new(100.0, -200.0, 300.0), geo).unwrap();
        for i in 0..3 {
            near(
                delta.to_array()[i],
                values[12 + i],
                1e-12,
                row,
                "ECEF difference",
            );
        }
        near(
            geometric_distance(
                EcefCoord::new(20_200_000.0, 14_000_000.0, 21_700_000.0),
                reference_ecef,
            )
            .unwrap(),
            values[15],
            1e-7,
            row,
            "geometric range",
        );
    }
}
