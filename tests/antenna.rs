use gnss_rust::antenna::*;
use gnss_rust::coordinates::{enu_to_ecef, geodetic_to_ecef};
use gnss_rust::time::{CalendarDateTime, TimeScale};
use gnss_rust::{EcefCoord, EnuCoord, GeodeticCoord, GnssTime};
use std::io::Cursor;

const INPUT: &str = include_str!("fixtures/synthetic_receiver.atx");
fn time() -> GnssTime {
    GnssTime::new(2300, 346200.0).unwrap()
}
fn delta() -> EnuCoord {
    EnuCoord {
        east: 0.025,
        north: -0.035,
        up: 0.8,
    }
}
fn model() -> ReceiverAntennaModel {
    let atx = read_antex(Cursor::new(INPUT)).unwrap();
    ReceiverAntennaModel::new(
        "synthetic_receiver.atx".into(),
        atx.receiver("TESTANT NONE", "RX001", time()).unwrap(),
        delta(),
    )
    .unwrap()
}
#[test]
fn native_loader_axes_pcv_and_range_correction_match() {
    let m = model();
    let geo = GeodeticCoord::new(35_f64.to_radians(), 139_f64.to_radians(), 45.0).unwrap();
    let rx = geodetic_to_ecef(geo).unwrap();
    let mut count = 0;
    for line in include_str!("fixtures/upstream_receiver_antenna.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = line.split(',').collect();
        let v: Vec<f64> = f[1..].iter().map(|x| x.parse().unwrap()).collect();
        let calibration = m.calibration();
        assert_eq!(
            calibration.frequencies[f[0]].offset_neu_m,
            [v[1], v[0], v[2]]
        );
        let pcv = calibration.pcv_m(f[0], 90.0 - v[4].to_degrees()).unwrap();
        assert!((pcv - v[5]).abs() < 1e-15);
        if v[4] <= std::f64::consts::FRAC_PI_2 {
            let offset = enu_to_ecef(
                EnuCoord {
                    east: 1e7 * v[3].sin() * v[4].cos(),
                    north: 1e7 * v[3].cos() * v[4].cos(),
                    up: 1e7 * v[4].sin(),
                },
                geo,
            )
            .unwrap();
            let sat = EcefCoord::new(rx.x + offset.x, rx.y + offset.y, rx.z + offset.z);
            let c = m
                .correction(
                    time(),
                    "G01".parse().unwrap(),
                    &["1C".into(), "2W".into()],
                    [2.5, -1.5],
                    [rx, sat],
                )
                .unwrap();
            let i = usize::from(f[0] == "G02");
            assert!((c.observation_add_m[i] - v[7]).abs() < 3e-12);
            assert!(
                (c.if_add_m - (2.5 * c.observation_add_m[0] - 1.5 * c.observation_add_m[1])).abs()
                    < 1e-15
            );
        }
        count += 1;
    }
    assert_eq!(count, 80);
}
#[test]
fn complete_type_radome_serial_and_inclusive_validity_are_required() {
    let a = read_antex(Cursor::new(INPUT)).unwrap();
    assert_eq!(a.entries().len(), 2);
    assert_eq!(a.entries()[1].satellite(), Some("G01".parse().unwrap()));
    let c = a.receiver(" testant   none ", "RX001", time()).unwrap();
    for t in [c.valid_from.unwrap(), c.valid_until.unwrap()] {
        assert!(a.receiver("TESTANT NONE", "RX001", t).is_ok());
    }
    for t in [
        c.valid_from.unwrap().checked_add_seconds(-0.001).unwrap(),
        c.valid_until.unwrap().checked_add_seconds(0.001).unwrap(),
    ] {
        assert_eq!(
            a.receiver("TESTANT NONE", "RX001", t),
            Err(AntennaError::MissingCalibration)
        );
    }
    for (name, serial) in [
        ("TESTANT SCIS", "RX001"),
        ("TESTANT NONE", ""),
        ("BLOCK IIF", "G01"),
    ] {
        assert_eq!(
            a.receiver(name, serial, time()),
            Err(AntennaError::MissingCalibration)
        );
    }
    let expected = CalendarDateTime {
        year: 2020,
        month: 1,
        day: 1,
        hour: 0,
        minute: 0,
        second: 0.0,
    }
    .to_gnss_time(TimeScale::Gps)
    .unwrap();
    assert_eq!(c.valid_from, Some(expected));
}
#[test]
fn overlapping_calibrations_are_ambiguous_without_first_match_fallback() {
    let start = INPUT.find("START OF ANTENNA").unwrap() - 60;
    let end = INPUT.find("END OF ANTENNA").unwrap() + "END OF ANTENNA\n".len();
    let duplicate = format!("{}{}", INPUT, &INPUT[start..end]);
    let a = read_antex(Cursor::new(duplicate)).unwrap();
    assert_eq!(
        a.receiver("TESTANT NONE", "RX001", time()),
        Err(AntennaError::AmbiguousCalibration)
    );
}
#[test]
fn malformed_or_unsupported_calibrations_never_become_zero() {
    let mutations = [
        INPUT.replace("     1.4", "     1.3"),
        INPUT.replacen(
            "A                                                           PCV TYPE",
            "R                                                           PCV TYPE",
            1,
        ),
        INPUT.replacen(
            "     0.0                                                    DAZI",
            "     5.0                                                    DAZI",
            1,
        ),
        INPUT.replacen(
            "   G02                                                      END OF FREQUENCY",
            "   G01                                                      END OF FREQUENCY",
            1,
        ),
        INPUT.replacen(
            "     2                                                      # OF FREQUENCIES",
            "     3                                                      # OF FREQUENCIES",
            1,
        ),
        INPUT
            .lines()
            .filter(|l| !l.contains("NORTH / EAST / UP"))
            .collect::<Vec<_>>()
            .join("\n"),
        INPUT
            .lines()
            .filter(|l| !l.contains("NOAZI"))
            .collect::<Vec<_>>()
            .join("\n"),
        INPUT.replacen("    0.0000000", "   60.0000000", 1),
        INPUT.replacen("    90.0    10.0", "    90.0     7.0", 1),
        INPUT.replacen("   NOAZI", "   0.000", 1),
        INPUT.replacen("   NOAZI", "   NOAZI NaN", 1),
        format!(
            "{}\n",
            INPUT
                .lines()
                .take(INPUT.lines().count() - 1)
                .collect::<Vec<_>>()
                .join("\n")
        ),
    ];
    for (i, bad) in mutations.iter().enumerate() {
        assert_ne!(bad, INPUT, "ineffective mutation {i}");
        assert!(read_antex(Cursor::new(bad)).is_err(), "mutation {i}");
    }
}
#[test]
fn missing_band_expiry_and_invalid_geometry_are_explicit() {
    let m = model();
    let rx = EcefCoord::new(4e6, 3e6, 4e6);
    let sat = EcefCoord::new(2e7, 1e7, 2e7);
    let c = |t, codes, pos| m.correction(t, "G01".parse().unwrap(), &codes, [2.5, -1.5], pos);
    assert_eq!(
        c(time(), ["1C".into(), "5Q".into()], [rx, sat]),
        Err(AntennaError::MissingFrequency("G05".into()))
    );
    assert_eq!(
        c(
            m.calibration()
                .valid_until
                .unwrap()
                .checked_add_seconds(1.0)
                .unwrap(),
            ["1C".into(), "2W".into()],
            [rx, sat]
        ),
        Err(AntennaError::OutsideValidity)
    );
    assert_eq!(
        c(time(), ["1C".into(), "2W".into()], [rx, rx]),
        Err(AntennaError::InvalidGeometry)
    );
    assert!(
        m.correction(
            time(),
            "E01".parse().unwrap(),
            &["1C".into(), "5Q".into()],
            [2.5, -1.5],
            [rx, sat]
        )
        .is_err()
    );
    assert!(m.calibration().pcv_m("G01", f64::NAN).is_err());
}
#[test]
fn modified_public_calibration_is_validated_before_filter_binding() {
    let m = model();
    let c = m.calibration();
    for i in 0..5 {
        let mut bad = c.clone();
        match i {
            0 => bad.frequencies.get_mut("G01").unwrap().noazi_m.clear(),
            1 => bad.zenith_step_deg = 0.0,
            2 => bad.frequencies.get_mut("G01").unwrap().offset_neu_m[0] = f64::INFINITY,
            3 => bad.valid_until = Some(bad.valid_from.unwrap().checked_add_seconds(-1.0).unwrap()),
            _ => bad.serial = "G01".into(),
        }
        assert!(ReceiverAntennaModel::new("calibration".into(), &bad, delta()).is_err());
    }
}

fn satellite_binding(pair: [String; 2]) -> SatelliteAntennaBinding {
    SatelliteAntennaBinding {
        product_id: "synthetic_satellite.atx".into(),
        clock_source: gnss_rust::precise::ClockSource::RinexClk,
        clock_reference: "synthetic-precise-clock-v1".into(),
        tracking_codes: pair,
    }
}
fn satellite_model(pair: [String; 2]) -> SatelliteAntennaModel {
    let atx = read_antex(Cursor::new(include_str!(
        "fixtures/synthetic_satellite.atx"
    )))
    .unwrap();
    SatelliteAntennaModel::new(satellite_binding(pair), atx.entries().to_vec()).unwrap()
}
#[test]
fn native_satellite_loader_and_nominal_yaw_if_pco_match_all_three_axes_and_pairs() {
    let mut count = 0;
    for line in include_str!("fixtures/upstream_satellite_antenna.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = line.split(',').collect();
        let n: Vec<f64> = f[5..].iter().map(|v| v.parse().unwrap()).collect();
        let m = satellite_model([f[3].into(), f[4].into()]);
        let sat = f[0].parse().unwrap();
        let time = GnssTime::new(f[1].parse().unwrap(), f[2].parse().unwrap()).unwrap();
        let c = m
            .correction(
                time,
                sat,
                EcefCoord::new(n[0], n[1], n[2]),
                EcefCoord::new(n[3], n[4], n[5]),
            )
            .unwrap();
        assert_eq!(c.satellite, sat);
        assert_eq!(c.frequency_codes[0], format!("G0{}", &f[3][..1]));
        for (a, b) in c.body_if_m.iter().zip(&n[6..9]) {
            assert!((a - b).abs() < 2e-14);
        }
        for (a, b) in [c.offset_ecef_m.x, c.offset_ecef_m.y, c.offset_ecef_m.z]
            .iter()
            .zip(&n[9..12])
        {
            assert!((a - b).abs() < 2e-12);
        }
        for (a, b) in [c.position_apc_m.x, c.position_apc_m.y, c.position_apc_m.z]
            .iter()
            .zip(&n[12..15])
        {
            assert!((a - b).abs() < 1e-6);
        }
        // No satellite PCV is applied even though fixture grids are nonzero.
        assert!(m.calibrations()[&sat].frequencies["G01"].noazi_m[0] > 0.0);
        count += 1;
    }
    assert_eq!(count, 81);
}
#[test]
fn satellite_selection_requires_prn_svn_validity_and_unique_fixed_calibration() {
    let atx = read_antex(Cursor::new(include_str!(
        "fixtures/synthetic_satellite.atx"
    )))
    .unwrap();
    let sat = "G01".parse().unwrap();
    let c = atx.satellite_calibration(sat, "G101", time()).unwrap();
    assert_eq!(c.serial, "G01");
    assert_eq!(c.svn, "G101");
    assert_eq!(
        atx.satellite_calibration(sat, "G102", time()),
        Err(AntennaError::MissingCalibration)
    );
    let pair = ["1C".into(), "2W".into()];
    let m = SatelliteAntennaModel::new(satellite_binding(pair.clone()), vec![c.clone()]).unwrap();
    let com = EcefCoord::new(2e7, 1e7, 1e7);
    let sun = EcefCoord::new(1e11, -3e10, 2e10);
    for t in [c.valid_from.unwrap(), c.valid_until.unwrap()] {
        assert!(m.correction(t, sat, com, sun).is_ok());
    }
    assert_eq!(
        m.correction(
            c.valid_until.unwrap().checked_add_seconds(0.001).unwrap(),
            sat,
            com,
            sun
        ),
        Err(AntennaError::OutsideValidity)
    );
    assert_eq!(
        m.correction(time(), "G02".parse().unwrap(), com, sun),
        Err(AntennaError::MissingCalibration)
    );
    assert_eq!(
        SatelliteAntennaModel::new(satellite_binding(pair.clone()), vec![c.clone(), c.clone()]),
        Err(AntennaError::AmbiguousCalibration)
    );
    let start = include_str!("fixtures/synthetic_satellite.atx")
        .find("START OF ANTENNA")
        .unwrap()
        - 60;
    let duplicate = format!(
        "{}{}",
        include_str!("fixtures/synthetic_satellite.atx"),
        &include_str!("fixtures/synthetic_satellite.atx")[start..]
    );
    let overlap = read_antex(Cursor::new(duplicate)).unwrap();
    assert_eq!(
        overlap.satellite_calibration(sat, "G101", time()),
        Err(AntennaError::AmbiguousCalibration)
    );
    let receiver = read_antex(Cursor::new(INPUT)).unwrap();
    assert!(
        SatelliteAntennaModel::new(satellite_binding(pair), vec![receiver.entries()[0].clone()])
            .is_err()
    );
}
#[test]
fn satellite_binding_rejects_missing_frequency_unsupported_pair_and_clock_datum() {
    let model = satellite_model(["1C".into(), "2W".into()]);
    let src = gnss_rust::precise::ClockSource::RinexClk;
    assert!(
        model
            .check_datum(
                &["1C".into(), "2W".into()],
                src,
                "synthetic-precise-clock-v1"
            )
            .is_ok()
    );
    for (pair, source, name) in [
        (["1C".into(), "2W".into()], src, "wrong"),
        (
            ["1C".into(), "2W".into()],
            gnss_rust::precise::ClockSource::Sp3,
            "synthetic-precise-clock-v1",
        ),
        (
            ["1C".into(), "2C".into()],
            src,
            "synthetic-precise-clock-v1",
        ),
    ] {
        assert_eq!(
            model.check_datum(&pair, source, name),
            Err(AntennaError::DatumMismatch)
        );
    }
    let entries = model.calibrations().values().cloned().collect::<Vec<_>>();
    for pair in [
        ["1C".into(), "1W".into()],
        ["1Z".into(), "2W".into()],
        ["1C".into(), "6C".into()],
    ] {
        assert!(SatelliteAntennaModel::new(satellite_binding(pair), entries.clone()).is_err());
    }
    let mut bad = entries[0].clone();
    bad.frequencies.remove("G02");
    assert_eq!(
        SatelliteAntennaModel::new(satellite_binding(["1C".into(), "2W".into()]), vec![bad]),
        Err(AntennaError::MissingFrequency("G02".into()))
    );
    let mut bad = entries[0].clone();
    bad.frequencies.get_mut("G01").unwrap().offset_neu_m[0] = f64::NAN;
    assert!(
        SatelliteAntennaModel::new(satellite_binding(["1C".into(), "2W".into()]), vec![bad])
            .is_err()
    );
    let mut binding = model.binding().clone();
    binding.clock_reference.clear();
    assert_eq!(
        SatelliteAntennaModel::new(binding, entries),
        Err(AntennaError::DatumMismatch)
    );
}
#[test]
fn satellite_pco_has_no_fallback_for_missing_or_degenerate_sun_attitude() {
    use gnss_rust::windup::WindupError;
    let m = satellite_model(["1C".into(), "2W".into()]);
    let sat = "G01".parse().unwrap();
    let com = EcefCoord::new(2e7, 0.0, 0.0);
    assert_eq!(
        m.correction(time(), sat, com, EcefCoord::new(0.0, 0.0, 0.0)),
        Err(AntennaError::Attitude(WindupError::MissingSun))
    );
    assert_eq!(
        m.correction(time(), sat, com, EcefCoord::new(1e11, 0.0, 0.0)),
        Err(AntennaError::Attitude(WindupError::DegenerateAttitude))
    );
    assert!(
        m.correction(
            time(),
            sat,
            EcefCoord::new(f64::NAN, 0.0, 0.0),
            EcefCoord::new(1e11, 1e11, 0.0)
        )
        .is_err()
    );
}
