use gnss_rust::rinex::{RinexObservationReader, read_navigation};
use gnss_rust::time::{CalendarDateTime, TimeScale};
use gnss_rust::{Error, GnssSystem, GnssTime, SatelliteId};
use std::io::Cursor;

fn header(content: &str, label: &str) -> String {
    format!("{content:<60}{label}\n")
}
fn observation_file(types: &str, body: &str) -> String {
    header(
        "     3.04           O                   M",
        "RINEX VERSION / TYPE",
    ) + &header(types, "SYS / # / OBS TYPES")
        + &header("", "END OF HEADER")
        + body
}
fn epoch(count: usize) -> String {
    format!("> 2024 02 08 00 00  0.0000000  0 {count:2}\n")
}
fn calendar(year: i32, month: u8, day: u8, hour: u8, minute: u8, second: f64) -> CalendarDateTime {
    CalendarDateTime {
        year,
        month,
        day,
        hour,
        minute,
        second,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn calendar_labels_have_correct_epochs_and_time_scales() {
    assert_eq!(
        calendar(1980, 1, 6, 0, 0, 0.0)
            .to_gnss_time(TimeScale::Gps)
            .unwrap(),
        GnssTime::default()
    );
    assert_eq!(
        calendar(2006, 1, 1, 0, 0, 0.0)
            .to_gnss_time(TimeScale::BeiDou)
            .unwrap(),
        GnssTime::new(1356, 14.0).unwrap()
    );
    let modern = calendar(2017, 1, 1, 0, 0, 0.125);
    near(
        modern
            .to_gnss_time(TimeScale::Utc)
            .unwrap()
            .difference_seconds(modern.to_gnss_time(TimeScale::Gps).unwrap()),
        18.0,
        1e-10,
    );
    let historic = calendar(1981, 7, 1, 0, 0, 0.0);
    assert_eq!(
        historic
            .to_gnss_time(TimeScale::Utc)
            .unwrap()
            .difference_seconds(historic.to_gnss_time(TimeScale::Gps).unwrap()),
        1.0
    );
    assert_eq!(
        modern.to_gnss_time(TimeScale::Galileo).unwrap(),
        modern.to_gnss_time(TimeScale::Gps).unwrap()
    );
    assert_eq!(
        modern.to_gnss_time(TimeScale::Qzss).unwrap(),
        modern.to_gnss_time(TimeScale::Gps).unwrap()
    );
    let before = calendar(1980, 1, 5, 23, 59, 59.5)
        .to_gnss_time(TimeScale::Gps)
        .unwrap();
    assert_eq!((before.week(), before.tow()), (-1, 604799.5));
}

#[test]
fn calendar_checks_dates_leap_years_and_fractional_seconds() {
    for date in [
        calendar(0, 1, 1, 0, 0, 0.0),
        calendar(2023, 2, 29, 0, 0, 0.0),
        calendar(1900, 2, 29, 0, 0, 0.0),
        calendar(2024, 0, 1, 0, 0, 0.0),
        calendar(2024, 13, 1, 0, 0, 0.0),
        calendar(2024, 4, 31, 0, 0, 0.0),
        calendar(2024, 1, 0, 0, 0, 0.0),
        calendar(2024, 1, 1, 24, 0, 0.0),
        calendar(2024, 1, 1, 0, 60, 0.0),
        calendar(2024, 1, 1, 0, 0, 60.0),
        calendar(2024, 1, 1, 0, 0, f64::NAN),
    ] {
        assert_eq!(
            date.to_gnss_time(TimeScale::Gps),
            Err(Error::InvalidCalendar)
        );
    }
    for year in [2000, 2024] {
        assert!(
            calendar(year, 2, 29, 0, 0, 0.125)
                .to_gnss_time(TimeScale::Gps)
                .is_ok()
        );
    }
}

#[test]
fn upstream_synthetic_fixture_preserves_both_carriers_and_all_epochs() {
    let reader =
        RinexObservationReader::new(Cursor::new(include_str!("fixtures/synthetic_spp.obs")))
            .unwrap();
    assert!(reader.header().approximate_position.unwrap().is_finite());
    let epochs = reader.collect::<Result<Vec<_>, _>>().unwrap();
    assert_eq!(epochs.len(), 8);
    let satellite_count = epochs[0].satellites().len();
    assert!(satellite_count >= 6);
    for (i, epoch) in epochs.iter().enumerate() {
        assert_eq!(
            epoch.time,
            GnssTime::new(2300, 345600.0 + i as f64 * 30.0).unwrap()
        );
        assert_eq!(epoch.satellites().len(), satellite_count);
        assert_eq!(epoch.observations.len(), 2 * satellite_count);
        for obs in &epoch.observations {
            assert!(obs.pseudorange_m.unwrap() > 1e7);
            assert!(obs.carrier_phase_cycles.unwrap() > 1e7);
            assert!(obs.doppler_hz.is_none());
            assert!(obs.snr_db_hz.is_none());
        }
    }
}

#[test]
fn missing_values_zero_values_lli_and_tracking_codes_are_distinct() {
    let record = format!(
        "G01{:14.3}  {:14.3}36{:16}{:14.3}  \n",
        0.0, 123.5, "", 45.0
    );
    let text = observation_file("G    4 C1C L1C D1C S1C", &(epoch(1) + &record));
    let mut reader = RinexObservationReader::new(Cursor::new(text)).unwrap();
    let e = reader.next().unwrap().unwrap();
    let obs = &e.observations[0];
    assert_eq!(obs.pseudorange_m, Some(0.0));
    assert_eq!(obs.carrier_phase_cycles, Some(123.5));
    assert_eq!(obs.doppler_hz, None);
    assert_eq!(obs.snr_db_hz, Some(45.0));
    assert!(obs.loss_of_lock());
    assert!(obs.half_cycle_ambiguity());
    assert_eq!(obs.signal_strength, Some(6));
    assert!(reader.next().is_none());
}

#[test]
fn observation_continuations_and_trailing_missing_fields_do_not_consume_next_satellite() {
    let types = "G    6 C1C L1C D1C S1C C2W L2W";
    let first = format!(
        "G01{:14.3}  {:14.3}  {:14.3}  {:14.3}  \n   {:14.3}  {:14.3}  \n",
        20e6, 1e8, -2.0, 45.0, 21e6, 8e7
    );
    let second = format!("G02{:14.3}  \n", 22e6);
    let text = observation_file(types, &(epoch(2) + &first + &second));
    let mut reader = RinexObservationReader::new(Cursor::new(text)).unwrap();
    let e = reader.next().unwrap().unwrap();
    assert_eq!(e.observations.len(), 4);
    let secondary = e
        .observations
        .iter()
        .find(|o| o.satellite.to_string() == "G01" && o.tracking_code == "2W")
        .unwrap();
    assert_eq!(secondary.carrier_phase_cycles, Some(8e7));
    let missing = e
        .observations
        .iter()
        .find(|o| o.satellite.to_string() == "G02" && o.tracking_code == "2W")
        .unwrap();
    assert!(missing.pseudorange_m.is_none());
    assert!(reader.next().is_none());
}

#[test]
fn observation_header_continuation_and_mid_file_header_update() {
    let codes = "C1C L1C D1C S1C C1W L1W D1W S1W C2W L2W D2W S2W C5Q";
    let text = header(
        "     3.04           O                   M",
        "RINEX VERSION / TYPE",
    ) + &header(&format!("G   14 {codes}"), "SYS / # / OBS TYPES")
        + &header("G      L5Q", "SYS / # / OBS TYPES")
        + &header("", "END OF HEADER")
        + "> 2024 02 08 00 00  0.0000000  4  1\n"
        + &header("G    1 C1C", "SYS / # / OBS TYPES")
        + &epoch(1)
        + &format!("G01{:14.3}  \n", 2e7);
    let mut reader = RinexObservationReader::new(Cursor::new(text)).unwrap();
    assert_eq!(
        reader.header().observation_types[&GnssSystem::Gps].len(),
        14
    );
    let e = reader.next().unwrap().unwrap();
    assert_eq!(e.observations.len(), 1);
    assert_eq!(
        reader.header().observation_types[&GnssSystem::Gps],
        vec!["C1C"]
    );
}

#[test]
fn glonass_channel_zero_is_present_and_duplicates_are_retained() {
    let text = header(
        "     3.04           O                   M",
        "RINEX VERSION / TYPE",
    ) + &header("R    1 C1C", "SYS / # / OBS TYPES")
        + &header("  3 R01 -7 R02 0 R02 6", "GLONASS SLOT / FRQ #")
        + &header("", "END OF HEADER")
        + &epoch(2)
        + &format!("R01{:14.3}  \nR02{:14.3}  \n", 2e7, 2.1e7);
    let mut reader = RinexObservationReader::new(Cursor::new(text)).unwrap();
    assert_eq!(reader.header().glonass_channels.len(), 3);
    assert_eq!(reader.header().glonass_channels[1].1, 0);
    let e = reader.next().unwrap().unwrap();
    assert_eq!(e.observations[0].glonass_channel, Some(-7));
    assert_eq!(e.observations[1].glonass_channel, Some(6));
}

#[test]
fn malformed_observations_fail_and_reader_is_fused_after_error() {
    let mut cases = vec![
        observation_file("G    1 C1C", &(epoch(2) + &format!("G01{:14.3}  \n", 2e7))),
        observation_file("G    1 C1C", &(epoch(1) + "G01           NaN  \n")),
        observation_file("G    1 C1C", &(epoch(1) + "G01x\n")),
        observation_file(
            "G    1 C1C",
            &(epoch(2) + &format!("G01{:14.3}  \nG01{:14.3}  \n", 2e7, 2e7)),
        ),
        observation_file("G    1 L1C", &(epoch(1) + &format!("G01{:14.3}86\n", 1e8))),
        observation_file("G    1 C1C", "> 2024 02 30 00 00  0.0000000  0  0\n"),
    ];
    cases.push(observation_file("G    1 C1C", &(epoch(1) + "é01\n")));
    for text in cases {
        let mut reader = RinexObservationReader::new(Cursor::new(text)).unwrap();
        assert!(reader.next().unwrap().is_err());
        assert!(reader.next().is_none());
        assert!(reader.next().is_none());
    }
}

#[test]
fn unsupported_or_malformed_headers_and_navigation_are_explicit_errors() {
    for text in [
        observation_file("G    2 C1C C1C", ""),
        observation_file("G    1 X1C", ""),
        observation_file("G    2 C1C", ""),
        observation_file("G    1 C1C", "").replacen("3.04", "2.11", 1),
    ] {
        assert!(RinexObservationReader::new(Cursor::new(text)).is_err());
    }
    let unsupported = header(
        "     3.04           N                   M",
        "RINEX VERSION / TYPE",
    ) + &header("", "END OF HEADER")
        + "R01 2024 02 08 00 00 00\n";
    assert!(
        read_navigation(Cursor::new(unsupported))
            .unwrap_err()
            .to_string()
            .contains("unsupported")
    );
    let valid = include_str!("fixtures/upstream_broadcast.nav");
    assert!(
        read_navigation(Cursor::new(
            valid.lines().take(6).collect::<Vec<_>>().join("\n")
        ))
        .is_err()
    );
}

#[test]
fn scaled_and_phase_shifted_data_are_not_silently_misinterpreted() {
    for (content, label) in [
        ("G   10", "SYS / SCALE FACTOR"),
        ("G L1C  0.25000", "SYS / PHASE SHIFT"),
        ("9999999999999999999", "GLONASS SLOT / FRQ #"),
    ] {
        let text = header(
            "     3.04           O                   M",
            "RINEX VERSION / TYPE",
        ) + &header("G    1 C1C", "SYS / # / OBS TYPES")
            + &header(content, label)
            + &header("", "END OF HEADER");
        assert!(RinexObservationReader::new(Cursor::new(text)).is_err());
    }
    let unit_scale = header(
        "     3.04           O                   M",
        "RINEX VERSION / TYPE",
    ) + &header("G    1 C1C", "SYS / # / OBS TYPES")
        + &header("G    1", "SYS / SCALE FACTOR")
        + &header("", "END OF HEADER");
    assert!(RinexObservationReader::new(Cursor::new(unit_scale)).is_ok());
}

#[test]
fn navigation_fortran_exponents_preserve_fields() {
    let original = include_str!("fixtures/upstream_broadcast.nav");
    let fortran = original.replace('e', "D");
    assert_eq!(
        read_navigation(Cursor::new(original)).unwrap(),
        read_navigation(Cursor::new(fortran)).unwrap()
    );
}

#[test]
fn receiver_clock_metadata_and_cycle_slip_epochs_are_retained() {
    let text = header(
        "     3.04           O                   M",
        "RINEX VERSION / TYPE",
    ) + &header("G    1 L1C", "SYS / # / OBS TYPES")
        + &header("     1", "RCV CLOCK OFFS APPL")
        + &header("", "END OF HEADER")
        + "> 2024 02 08 00 00  0.1250000  6  1 0.000002\n"
        + &format!("G01{:14.3}16\n", 1e8);
    let mut reader = RinexObservationReader::new(Cursor::new(text)).unwrap();
    assert!(reader.header().clock_offsets_applied);
    let epoch = reader.next().unwrap().unwrap();
    assert_eq!(epoch.flag, 6);
    assert_eq!(epoch.time, GnssTime::new(2300, 345600.125).unwrap());
    assert_eq!(epoch.receiver_clock_offset_s, Some(2e-6));
    assert!(epoch.observations[0].loss_of_lock());
    assert!(reader.next().is_none());
}

#[test]
fn rinex_kepler_navigation_agrees_with_twenty_actual_cpp_state_cases() {
    let nav =
        read_navigation(Cursor::new(include_str!("fixtures/upstream_broadcast.nav"))).unwrap();
    assert_eq!(nav.len(), 5);
    let rows: Vec<_> = include_str!("fixtures/upstream_navigation.csv")
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(rows.len(), 20);
    for line in rows {
        let values: Vec<_> = line.split(',').collect();
        assert_eq!(values.len(), 10);
        let satellite: SatelliteId = values[0].parse().unwrap();
        let v: Vec<f64> = values[1..].iter().map(|s| s.parse().unwrap()).collect();
        let time = GnssTime::new(2300, 345600.0 + v[0]).unwrap();
        let eph = nav.ephemeris(satellite, time, 14400.0).unwrap();
        assert_eq!(eph.toe, GnssTime::new(2300, 345600.0).unwrap());
        assert_eq!(eph.toc, eph.toe);
        let state = eph.satellite_state(time).unwrap();
        for i in 0..3 {
            near(state.position.to_array()[i], v[1 + i], 1e-4);
            near(state.velocity_mps[i], v[4 + i], 1e-6);
        }
        near(state.clock_bias_s, v[7], 1e-15);
        near(state.clock_drift, v[8], 1e-15);
    }
}

#[test]
fn unhealthy_stale_and_missing_navigation_cannot_be_selected() {
    let mut nav =
        read_navigation(Cursor::new(include_str!("fixtures/upstream_broadcast.nav"))).unwrap();
    let satellite: SatelliteId = "G01".parse().unwrap();
    let time = GnssTime::new(2300, 345600.0).unwrap();
    assert_eq!(
        nav.ephemeris(
            satellite,
            time.checked_add_seconds(14401.0).unwrap(),
            14400.0
        )
        .unwrap_err(),
        Error::MissingEphemeris
    );
    assert_eq!(
        nav.ephemeris("G31".parse().unwrap(), time, 14400.0)
            .unwrap_err(),
        Error::MissingEphemeris
    );
    assert_eq!(
        nav.ephemeris(satellite, time, f64::NAN).unwrap_err(),
        Error::InvalidConfiguration
    );
    let mut unhealthy = nav.ephemeris(satellite, time, 14400.0).unwrap().clone();
    unhealthy.health = 1;
    unhealthy.satellite = "G02".parse().unwrap();
    nav.add_ephemeris(unhealthy.clone()).unwrap();
    assert_eq!(
        nav.ephemeris(unhealthy.satellite, time, 14400.0)
            .unwrap_err(),
        Error::MissingEphemeris
    );
    unhealthy.orbit.eccentricity = 1.0;
    assert_eq!(
        unhealthy.satellite_state(time),
        Err(Error::InvalidEphemeris)
    );
    unhealthy.orbit.eccentricity = f64::NAN;
    assert_eq!(nav.add_ephemeris(unhealthy), Err(Error::InvalidEphemeris));
}

#[test]
fn navigation_uses_full_week_for_age_admission_not_modulo_week() {
    let nav =
        read_navigation(Cursor::new(include_str!("fixtures/upstream_broadcast.nav"))).unwrap();
    let satellite: SatelliteId = "G01".parse().unwrap();
    assert_eq!(
        nav.ephemeris(satellite, GnssTime::new(2301, 345600.0).unwrap(), 14400.0)
            .unwrap_err(),
        Error::MissingEphemeris
    );
    let eph = nav.records().find(|e| e.satellite == satellite).unwrap();
    assert!(
        eph.satellite_state(GnssTime::new(2301, 345600.0).unwrap())
            .is_ok()
    );
    // The raw broadcast evaluator intentionally wraps orbital time, but the
    // navigation selector must never admit a week-old record for a real solve.
}
