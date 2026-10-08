use gnss_rust::models::niell_hydrostatic_mapping;
use gnss_rust::time::{CalendarDateTime, TimeScale, day_of_year_gps_scale};
use gnss_rust::{Error, GeodeticCoord};

#[test]
fn native_96_mapping_cases_cover_season_hemisphere_height_and_low_elevation() {
    let mut count = 0;
    for line in include_str!("fixtures/upstream_niell.csv").lines().skip(1) {
        let v: Vec<_> = line.split(',').collect();
        let n = |i: usize| v[i].parse::<f64>().unwrap();
        let actual = niell_hydrostatic_mapping(
            GeodeticCoord::new(n(0), 0.0, n(1)).unwrap(),
            n(2),
            v[3].parse().unwrap(),
        )
        .unwrap();
        assert!((actual - n(4)).abs() < 1e-12, "{actual} vs {}", n(4));
        count += 1;
    }
    assert_eq!(count, 96);
}
#[test]
fn calendar_day_is_gps_label_and_handles_leap_year_and_boundaries() {
    for (year, month, day, expected) in [
        (2024, 2, 8, 39),
        (2024, 2, 29, 60),
        (2024, 12, 31, 366),
        (2025, 1, 1, 1),
        (2025, 3, 1, 60),
        (2000, 3, 1, 61),
        (1900, 3, 1, 60),
    ] {
        let time = CalendarDateTime {
            year,
            month,
            day,
            hour: 23,
            minute: 59,
            second: 59.9,
        }
        .to_gnss_time(TimeScale::Gps)
        .unwrap();
        assert_eq!(day_of_year_gps_scale(time).unwrap(), expected);
    }
    let last = CalendarDateTime {
        year: 2024,
        month: 12,
        day: 31,
        hour: 23,
        minute: 59,
        second: 59.9,
    }
    .to_gnss_time(TimeScale::Gps)
    .unwrap();
    assert_eq!(
        day_of_year_gps_scale(last.checked_add_seconds(0.2).unwrap()).unwrap(),
        1
    );
    let last = CalendarDateTime {
        year: 9999,
        month: 12,
        day: 31,
        hour: 23,
        minute: 59,
        second: 59.9,
    }
    .to_gnss_time(TimeScale::Gps)
    .unwrap();
    assert_eq!(
        day_of_year_gps_scale(last.checked_add_seconds(0.2).unwrap()),
        Err(Error::OutOfRange)
    );
}
#[test]
fn mapping_validates_geometry_angles_and_day() {
    let geo = GeodeticCoord::new(0.0, 0.0, 0.0).unwrap();
    for day in [0, 367] {
        assert_eq!(
            niell_hydrostatic_mapping(geo, 0.5, day),
            Err(Error::InvalidConfiguration)
        );
    }
    for elevation in [-0.1, 2.0, f64::NAN] {
        assert_eq!(
            niell_hydrostatic_mapping(geo, elevation, 39),
            Err(Error::InvalidConfiguration)
        );
    }
    assert!(
        niell_hydrostatic_mapping(
            GeodeticCoord {
                height: f64::INFINITY,
                ..geo
            },
            0.5,
            39
        )
        .is_err()
    );
}
