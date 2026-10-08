use gnss_rust::bias::*;
use gnss_rust::dual_frequency::DualFrequencyError;
use gnss_rust::precise::ClockSource;
use gnss_rust::time::TimeScale;
use gnss_rust::{Error, GnssTime, constants as c};
use std::io::{self, Cursor};

fn parse(text: &str) -> Result<BiasProducts, BiasError> {
    read_bias_sinex(Cursor::new(text), Default::default())
}
fn file(rows: &str) -> String {
    format!(
        "%=BIA 1.00 TEST 2024:039:00000 TEST\n+BIAS/DESCRIPTION\n TIME_SYSTEM G\n BIAS_MODE ABSOLUTE\n-BIAS/DESCRIPTION\n+BIAS/SOLUTION\n{rows}\n-BIAS/SOLUTION\n%=ENDBIA\n"
    )
}
const ROW: &str = " OSB G080 G01           C1C      2024:039:00600 2024:039:02520 ns 1.234 0.1";
fn time() -> GnssTime {
    GnssTime::new(2300, 346200.0).unwrap()
}
fn binding() -> OsbBinding {
    OsbBinding {
        product_id: "test-product".into(),
        clock_source: ClockSource::RinexClk,
        clock_reference: "clock-v1".into(),
        phase_reference: Some("phase-v1".into()),
    }
}
fn request(tracking: &str) -> BiasSignalRequest {
    BiasSignalRequest {
        satellite: "G01".parse().unwrap(),
        tracking_code: tracking.into(),
        glonass_channel: None,
    }
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(a.is_finite() && (a - b).abs() <= tol, "{a} vs {b}");
}

#[test]
fn native_loader_39_rows_and_literal_conversion_sign_match() {
    let p = parse(include_str!("fixtures/upstream_bias.bsx")).unwrap();
    assert_eq!(p.entries().len(), 39);
    let reference: Vec<_> = include_str!("fixtures/upstream_bias.csv")
        .lines()
        .skip(1)
        .collect();
    for (entry, row) in p.entries().iter().zip(&reference) {
        let f: Vec<_> = row.split(',').collect();
        assert_eq!(
            entry.kind,
            if f[0] == "OSB" {
                BiasKind::Osb
            } else {
                BiasKind::Dsb
            }
        );
        assert_eq!(entry.satellite.to_string(), f[1]);
        assert_eq!(entry.observation_1, f[2]);
        assert_eq!(entry.observation_2.as_deref(), Some(f[3]));
        assert_eq!(entry.unit_label, f[4]);
        assert_eq!(entry.value, f[5].parse::<f64>().unwrap());
        assert_eq!(entry.sigma, f[6].parse::<f64>().unwrap());
        let (m, sigma) = entry.metres(None).unwrap();
        close(m, f[7].parse().unwrap(), 2e-15);
        close(-m, f[8].parse().unwrap(), 2e-15);
        assert!(sigma >= 0.0);
        if entry.kind == BiasKind::Osb
            && entry.observation_1.starts_with('C')
            && entry.contains(time())
        {
            let s = p
                .corrections_at(
                    time(),
                    &binding(),
                    &[BiasSignalRequest {
                        satellite: entry.satellite,
                        tracking_code: entry.observation_1[1..].into(),
                        glonass_channel: None,
                    }],
                )
                .unwrap();
            close(
                s.corrections.entries().values().next().unwrap()[0].code_add_m,
                f[8].parse().unwrap(),
                2e-15,
            );
        }
    }
    assert_eq!(reference.len(), 39);
}
#[test]
fn standard_osb_blank_second_observable_retains_svn_time_sigma_and_binding() {
    let p = parse(include_str!("fixtures/synthetic_ppp_osb.bsx")).unwrap();
    assert_eq!(p.entries().len(), 36);
    assert_eq!(p.time_scale(), TimeScale::Gps);
    assert_eq!(p.description()["BIAS_MODE"], [vec!["ABSOLUTE"]]);
    assert!(p.header().starts_with("%=BIA 1.00 SYN"));
    assert!(
        p.entries()
            .iter()
            .all(|e| e.observation_2.is_none() && e.station.is_none())
    );
    let snapshot = p
        .corrections_at(time(), &binding(), &[request("1C"), request("2W")])
        .unwrap();
    let row = &snapshot.corrections.entries()[&("G01".parse().unwrap(), "1C".into())][0];
    close(row.code_add_m, -1.2, 1e-14);
    close(row.phase_add_m.unwrap(), -0.008, 1e-16);
    assert_eq!(row.valid_from, time());
    assert_eq!(row.valid_until, time());
    assert_eq!(row.clock_source, ClockSource::RinexClk);
    assert_eq!(row.clock_reference, "clock-v1");
    assert_eq!(row.phase_reference.as_deref(), Some("phase-v1"));
    let ledger = &snapshot.provenance[&(row.satellite, row.tracking_code.clone())];
    assert_eq!(ledger.product_id, "test-product");
    assert_eq!(ledger.code.svn.as_deref(), Some("G081"));
    assert_eq!(ledger.code.start_label, "2024:039:00600");
    assert_eq!(ledger.phase.as_ref().unwrap().observation_1, "L1C");
    close(ledger.code_sigma_m, 0.001 * c::SPEED_OF_LIGHT * 1e-9, 1e-18);
}
#[test]
fn exact_half_open_boundary_selects_next_bias_and_never_reuses_expired_osb() {
    let next = ROW
        .replace("00600 2024:039:02520", "02520 2024:040:00000")
        .replace("1.234", "2.345");
    let p = parse(&file(&format!("{next}\n{ROW}"))).unwrap(); // input need not be sorted
    let sat = "G01".parse().unwrap();
    assert_eq!(p.select_osb(sat, "C1C", time()).unwrap().value, 1.234);
    let boundary = time().checked_add_seconds(1920.0).unwrap();
    assert_eq!(
        p.select_osb(sat, "C1C", boundary.checked_add_seconds(-1e-6).unwrap())
            .unwrap()
            .value,
        1.234
    );
    assert_eq!(p.select_osb(sat, "C1C", boundary).unwrap().value, 2.345);
    assert!(matches!(
        p.select_osb(sat, "C1C", time().checked_add_seconds(-1.0).unwrap()),
        Err(BiasError::ExpiredOsb { .. })
    ));
    let end = GnssTime::new(2300, 432000.0).unwrap();
    assert!(matches!(
        p.select_osb(sat, "C1C", end),
        Err(BiasError::ExpiredOsb { .. })
    ));
    assert!(matches!(
        p.select_osb(sat, "C1W", time()),
        Err(BiasError::MissingOsb { .. })
    ));
}
#[test]
fn overlap_duplicate_alias_and_open_ended_intervals_are_validated() {
    for second in [
        ROW.to_owned(),
        ROW.replace("C1C      ", "C1C C1C "),
        ROW.replace("00600", "00601"),
    ] {
        assert!(matches!(
            parse(&file(&format!("{ROW}\n{second}"))),
            Err(BiasError::OverlappingEntries { .. })
        ));
    }
    let open = ROW.replace(
        "2024:039:00600 2024:039:02520",
        "0000:000:00000 0000:000:00000",
    );
    let p = parse(&file(&open)).unwrap();
    assert_eq!(p.entries()[0].valid_from, None);
    assert_eq!(p.entries()[0].valid_until, None);
    assert!(p.entries()[0].contains(time().checked_add_seconds(-604800.0).unwrap()));
    assert!(matches!(
        parse(&file(&format!("{open}\n{ROW}"))),
        Err(BiasError::OverlappingEntries { .. })
    ));
}
#[test]
fn station_rows_and_dsb_do_not_supply_absolute_satellite_corrections() {
    let station = ROW.replace("G01           ", "G01 ABCD00USA ");
    let dsb = ROW.replace("OSB", "DSB").replace("C1C      ", "C1C C2W ");
    let p = parse(&file(&format!("{station}\n{dsb}"))).unwrap();
    assert_eq!(p.entries()[0].station.as_deref(), Some("ABCD00USA"));
    assert_eq!(p.entries()[1].kind, BiasKind::Dsb);
    assert!(matches!(
        p.corrections_at(time(), &binding(), &[request("1C")]),
        Err(BiasError::MissingOsb { .. })
    ));
    let paired = ROW.replace("C1C      ", "C1C C2W ");
    let p = parse(&file(&paired)).unwrap();
    assert!(matches!(
        p.corrections_at(time(), &binding(), &[request("1C")]),
        Err(BiasError::AmbiguousOsb { .. })
    ));
    let p = parse(&file(ROW).replace("ABSOLUTE", "RELATIVE")).unwrap();
    assert!(matches!(
        p.corrections_at(time(), &binding(), &[request("1C")]),
        Err(BiasError::RelativeProduct)
    ));
    let relative_header =
        "%=BIA 1.00 SYN 2024:039:00000 SYN 2024:039:00000 2024:040:00000 R 00000001";
    let text = file(ROW).replacen("%=BIA 1.00 TEST 2024:039:00000 TEST", relative_header, 1);
    let p = parse(&text).unwrap();
    assert!(matches!(
        p.corrections_at(time(), &binding(), &[request("1C")]),
        Err(BiasError::RelativeProduct)
    ));
}
#[test]
fn missing_phase_or_unbound_phase_preserves_code_without_zero_phase_osb() {
    let p = parse(&file(ROW)).unwrap();
    let s = p
        .corrections_at(time(), &binding(), &[request("1C")])
        .unwrap();
    let e = &s.corrections.entries().values().next().unwrap()[0];
    assert_eq!(e.phase_add_m, None);
    assert_eq!(e.phase_reference, None);
    assert!(s.provenance.values().next().unwrap().phase.is_none());
    assert_eq!(
        s.provenance.values().next().unwrap().phase_unavailability,
        Some(PhaseOsbUnavailable::Missing)
    );
    let phase = ROW.replace("C1C", "L1C");
    let p = parse(&file(&format!("{ROW}\n{phase}"))).unwrap();
    let mut b = binding();
    b.phase_reference = None;
    let s = p.corrections_at(time(), &b, &[request("1C")]).unwrap();
    assert_eq!(
        s.provenance.values().next().unwrap().phase_unavailability,
        Some(PhaseOsbUnavailable::NotBound)
    );
    assert_eq!(
        s.corrections.entries().values().next().unwrap()[0].phase_add_m,
        None
    );
    let expired = phase.replace(
        "2024:039:00600 2024:039:02520",
        "2024:038:00600 2024:039:00600",
    );
    let p = parse(&file(&format!("{ROW}\n{expired}"))).unwrap();
    let s = p
        .corrections_at(time(), &binding(), &[request("1C")])
        .unwrap();
    assert_eq!(
        s.provenance.values().next().unwrap().phase_unavailability,
        Some(PhaseOsbUnavailable::Expired)
    );
}
#[test]
fn phase_cycle_units_use_exact_frequency_and_require_glonass_channel() {
    let gps = ROW.replace("C1C", "L1C").replace("ns 1.234", "cyc 2.5");
    let glo = gps.replace("G080 G01", "R080 R01");
    let p = parse(&file(&format!("{gps}\n{glo}"))).unwrap();
    close(
        p.entries()[0].metres(None).unwrap().0,
        2.5 * c::GPS_L1_WAVELENGTH,
        1e-15,
    );
    close(
        p.entries()[0].metres(None).unwrap().1,
        0.1 * c::GPS_L1_WAVELENGTH,
        1e-16,
    );
    assert!(matches!(
        p.entries()[1].metres(None),
        Err(BiasError::Correction(DualFrequencyError::Gnss(
            Error::MissingGlonassChannel
        )))
    ));
    assert!(matches!(
        p.entries()[1].metres(Some(7)),
        Err(BiasError::Correction(DualFrequencyError::Gnss(
            Error::InvalidGlonassChannel
        )))
    ));
    close(
        p.entries()[1].metres(Some(-7)).unwrap().0,
        2.5 * c::SPEED_OF_LIGHT / (c::GLO_L1_BASE_FREQ - 7.0 * c::GLO_L1_STEP_FREQ),
        1e-15,
    );
    assert!(parse(&file(&gps.replace("L1C", "C1C"))).is_err());
}
#[test]
fn epoch_scale_requires_declaration_and_converts_utc_and_bdt() {
    let gps = parse(&file(ROW)).unwrap().entries()[0].valid_from.unwrap();
    for (code, offset) in [("U", 18.0), ("C", 14.0), ("E", 0.0), ("J", 0.0)] {
        let p = parse(&file(ROW).replace("TIME_SYSTEM G", &format!("TIME_SYSTEM {code}"))).unwrap();
        assert_eq!(
            p.entries()[0].valid_from.unwrap().difference_seconds(gps),
            offset
        );
    }
    let absent = file(ROW).replace(" TIME_SYSTEM G\n", "");
    assert!(matches!(parse(&absent), Err(BiasError::MissingTimeScale)));
    assert_eq!(
        read_bias_sinex(
            Cursor::new(&absent),
            BiasReadOptions {
                time_scale: Some(TimeScale::Gps)
            }
        )
        .unwrap()
        .entries()[0]
            .valid_from,
        Some(gps)
    );
    assert!(matches!(
        read_bias_sinex(
            Cursor::new(file(ROW)),
            BiasReadOptions {
                time_scale: Some(TimeScale::Utc)
            }
        ),
        Err(BiasError::ConflictingTimeScale)
    ));
    assert!(parse(&file(ROW).replace("TIME_SYSTEM G", "TIME_SYSTEM R")).is_err());
}
#[test]
fn utc_leap_boundary_and_leap_year_ordinals_are_converted_as_calendar_dates() {
    let text = file(" OSB G01 C1C 2016:366:86399 2017:001:00000 m 1 0")
        .replace("TIME_SYSTEM G", "TIME_SYSTEM U");
    let p = parse(&text).unwrap();
    let e = &p.entries()[0];
    assert_eq!(
        e.valid_until
            .unwrap()
            .difference_seconds(e.valid_from.unwrap()),
        2.0
    );
    let text = file(" OSB G01 C1C 2024:366:86399 2025:001:00000 m 1 0");
    let p = parse(&text).unwrap();
    let e = &p.entries()[0];
    assert_eq!(
        e.valid_until
            .unwrap()
            .difference_seconds(e.valid_from.unwrap()),
        1.0
    );
    assert!(parse(&text.replace("2024:366", "1900:366")).is_err());
}
#[test]
fn malformed_unsupported_nonfinite_and_truncated_inputs_fail_with_line_numbers() {
    let variants = [
        ROW.replace("ns", "s"),
        ROW.replace("1.234", "NaN"),
        ROW.replace("0.1", "-0.1"),
        ROW.replace("039:00600", "000:00600"),
        ROW.replace("039:00600", "039:89999"),
        ROW.replace("039:02520", "039:00600"),
        ROW.replace("G080", "E080"),
        ROW.replace("G01", "G00"),
        ROW.replace("C1C", "C1"),
        format!("{ROW} 0.01 0.001"),
        ROW.replace("OSB", "ISB"),
        ROW.replace("C1C", "D1C"),
        "OSB G01".into(),
        "OSB 2024:039:00600".into(),
    ];
    for row in variants {
        match parse(&file(&row)).unwrap_err() {
            BiasError::Parse { line, .. } | BiasError::Unsupported { line, .. } => {
                assert!(line > 0)
            }
            e => panic!("unexpected error {e} for {row}"),
        }
    }
    for text in [
        file(ROW).replace("%=ENDBIA\n", ""),
        file(ROW).replace("-BIAS/SOLUTION", "-BIAS/DESCRIPTION"),
        file(ROW).replace("1.00", "2.00"),
        file(ROW).replace("+BIAS/SOLUTION", "+BIAS/SOLUTION\n+BIAS/SOLUTION"),
        file(ROW).replace("1.234", "é"),
        file(ROW) + ROW,
    ] {
        assert!(parse(&text).is_err());
    }
    let p = parse(&file(&ROW.replace("1.234", "1.234D+00"))).unwrap();
    assert_eq!(p.entries()[0].value, 1.234);
}
#[test]
fn invalid_binding_duplicate_requests_and_invalid_external_entries_fail() {
    let p = parse(&file(ROW)).unwrap();
    for field in 0..3 {
        let mut b = binding();
        match field {
            0 => b.product_id.clear(),
            1 => b.clock_reference.clear(),
            _ => b.phase_reference = Some(" ".into()),
        }
        assert!(matches!(
            p.corrections_at(time(), &b, &[request("1C")]),
            Err(BiasError::InvalidBinding)
        ));
    }
    assert!(
        p.corrections_at(time(), &binding(), &[request("1C"), request("1C")])
            .is_err()
    );
    assert!(
        p.corrections_at(time(), &binding(), &[request("C1C")])
            .is_err()
    );
    let p = parse(&file(
        &ROW.replace("1.234", "1e308")
            .replace("C1C", "L1C")
            .replace("ns", "cyc"),
    ))
    .unwrap();
    // Frequency conversion is finite for this value; deliberately corrupt a
    // public entry to verify conversion validates external construction too.
    let mut e = p.entries()[0].clone();
    e.value = f64::INFINITY;
    assert!(e.metres(None).is_err());
    e.observation_1 = "".into();
    assert!(e.metres(None).is_err());
}
#[test]
fn repeated_clock_reference_metadata_is_preserved_without_datum_guessing() {
    let text=file(ROW).replace(" TIME_SYSTEM G", " SATELLITE_CLOCK_REFERENCE_OBSERVABLES G C1W C2W\n SATELLITE_CLOCK_REFERENCE_OBSERVABLES R C1P C2P\n TIME_SYSTEM G");
    let p = parse(&text).unwrap();
    assert_eq!(
        p.description()["SATELLITE_CLOCK_REFERENCE_OBSERVABLES"],
        [vec!["G", "C1W", "C2W"], vec!["R", "C1P", "C2P"]]
    );
    assert!(parse(&file(ROW).replace("TIME_SYSTEM G", "TIME_SYSTEM G\n TIME_SYSTEM G")).is_err());
    assert!(
        parse(&file(ROW).replace(
            "BIAS_MODE ABSOLUTE",
            "BIAS_MODE ABSOLUTE\n BIAS_MODE RELATIVE"
        ))
        .is_err()
    );
}

#[test]
fn io_errors_are_not_reported_as_empty_products() {
    struct Broken;
    impl io::Read for Broken {
        fn read(&mut self, _: &mut [u8]) -> io::Result<usize> {
            Err(io::Error::other("broken"))
        }
    }
    impl io::BufRead for Broken {
        fn fill_buf(&mut self) -> io::Result<&[u8]> {
            Err(io::Error::other("broken"))
        }
        fn consume(&mut self, _: usize) {}
    }
    assert!(matches!(
        read_bias_sinex(Broken, Default::default()),
        Err(BiasError::Io(_))
    ));
}
