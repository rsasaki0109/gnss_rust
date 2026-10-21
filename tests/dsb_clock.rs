use gnss_rust::bias::{BiasProducts, BiasSignalRequest, OsbBinding, read_bias_sinex};
use gnss_rust::dsb_clock::{DsbClockError, GpsIfClockBinding, GpsIfDsbClock};
use gnss_rust::precise::ClockSource;
use gnss_rust::{GnssTime, constants as c};
use std::io::Cursor;

fn binding() -> GpsIfClockBinding {
    GpsIfClockBinding {
        product_id: "declared-relative-product".into(),
        clock_source: ClockSource::RinexClk,
        clock_reference: "caller-asserted-GPS-IF-C1W-C2W".into(),
    }
}
fn request(code: &str) -> BiasSignalRequest {
    BiasSignalRequest {
        satellite: "G01".parse().unwrap(),
        tracking_code: code.into(),
        glonass_channel: None,
    }
}
fn time() -> GnssTime {
    GnssTime::new(2300, 346200.0).unwrap()
}
fn file(rows: &str) -> String {
    format!(
        "%=BIA 1.00 SYN 2024:039:00000 SYN 2024:039:00000 2024:040:00000 R 00000003\n+BIAS/DESCRIPTION\n TIME_SYSTEM G\n BIAS_MODE RELATIVE\n-BIAS/DESCRIPTION\n+BIAS/SOLUTION\n{rows}-BIAS/SOLUTION\n%=ENDBIA\n"
    )
}
fn row(a: &str, b: &str, value: f64, sigma: f64, ns: bool) -> String {
    let scale = if ns { c::SPEED_OF_LIGHT * 1e-9 } else { 1.0 };
    format!(
        " DSB G080 G01           {a} {b} 2024:039:00600 2024:039:02520 {} {:.17e} {:.17e}\n",
        if ns { "ns" } else { "m" },
        value / scale,
        sigma / scale
    )
}
fn product(d12: f64, d1c: f64, d2l: f64, reverse: bool) -> BiasProducts {
    let mut rows = String::new();
    for (a, b, d, s, ns) in [
        ("C1W", "C2W", d12, 0.01, true),
        ("C1C", "C1W", d1c, 0.02, false),
        ("C2W", "C2L", d2l, 0.03, true),
    ] {
        rows += &if reverse {
            row(b, a, -d, s, ns)
        } else {
            row(a, b, d, s, ns)
        };
    }
    read_bias_sinex(Cursor::new(file(&rows)), Default::default()).unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-11, "{a} versus {b}");
}

#[test]
fn all_24_independent_eigen_gauges_and_96_code_corrections_match() {
    for line in include_str!("fixtures/upstream_dsb_clock.csv")
        .lines()
        .skip(1)
    {
        let v: Vec<f64> = line.split(',').map(|s| s.parse().unwrap()).collect();
        assert_eq!(v.len(), 12);
        for reverse in [false, true] {
            let p = product(v[1], v[2], v[3], reverse);
            let model = GpsIfDsbClock::new(&p).unwrap();
            let snap = model
                .corrections_at(
                    time(),
                    &binding(),
                    &[request("1W"), request("2W"), request("1C"), request("2L")],
                )
                .unwrap();
            for (j, code) in ["1W", "2W", "1C", "2L"].iter().enumerate() {
                let d = &snap.provenance[&("G01".parse().unwrap(), code.to_string())];
                close(d.code_add_m, v[4 + j]);
                close(d.sigma_upper_bound_m, v[8 + j]);
                let value = d
                    .source_terms
                    .iter()
                    .map(|t| t.coefficient * t.entry.metres(None).unwrap().0)
                    .sum::<f64>();
                close(value, d.code_add_m);
            }
            let a = c::IFLC_F1;
            let b = c::IFLC_F2;
            close(
                a * snap.provenance[&("G01".parse().unwrap(), "1W".into())].code_add_m
                    + b * snap.provenance[&("G01".parse().unwrap(), "2W".into())].code_add_m,
                0.0,
            );
        }
    }
}

#[test]
fn explicit_gauge_corrects_injected_raw_codes_but_never_invents_phase_osbs() {
    let a = c::IFLC_F1;
    let b = c::IFLC_F2;
    let biases = [2.0, -a / b * 2.0, 3.5, -1.2];
    let p = product(
        biases[0] - biases[1],
        biases[2] - biases[0],
        biases[1] - biases[3],
        false,
    );
    let snap = GpsIfDsbClock::new(&p)
        .unwrap()
        .corrections_at(
            time(),
            &binding(),
            &[request("1W"), request("2W"), request("1C"), request("2L")],
        )
        .unwrap();
    for (j, code) in ["1W", "2W", "1C", "2L"].iter().enumerate() {
        let rows = &snap.corrections.entries()[&("G01".parse().unwrap(), code.to_string())];
        let correction = &rows[0];
        close(
            20_000_000.0 + biases[j] + correction.code_add_m,
            20_000_000.0,
        );
        assert_eq!(correction.phase_add_m, None);
        assert_eq!(correction.phase_reference, None);
        assert_eq!(correction.valid_from, time());
        assert_eq!(correction.valid_until, time());
        assert_eq!(correction.clock_reference, binding().clock_reference);
    }
    let osb = OsbBinding {
        product_id: "relative".into(),
        clock_source: ClockSource::RinexClk,
        clock_reference: "IF".into(),
        phase_reference: None,
    };
    assert!(matches!(
        p.corrections_at(time(), &osb, &[request("1C")]),
        Err(gnss_rust::bias::BiasError::RelativeProduct)
    ));
}

#[test]
fn product_half_open_end_is_not_admitted_as_an_inclusive_snapshot() {
    let p = product(1.0, 2.0, 3.0, false);
    let m = GpsIfDsbClock::new(&p).unwrap();
    assert!(
        m.corrections_at(time(), &binding(), &[request("1C")])
            .is_ok()
    );
    let end = time().checked_add_seconds(1920.0).unwrap();
    assert!(matches!(
        m.corrections_at(end, &binding(), &[request("1C")]),
        Err(DsbClockError::ExpiredPair { .. })
    ));
    assert!(
        m.corrections_at(
            end.checked_add_seconds(-1e-6).unwrap(),
            &binding(),
            &[request("1C")]
        )
        .is_ok()
    );
}

#[test]
fn alternate_graph_paths_are_not_silently_substituted_or_averaged() {
    let rows = row("C1W", "C2W", 1.0, 0.01, false) + &row("C1C", "C2W", 2.0, 0.02, false);
    let p = read_bias_sinex(Cursor::new(file(&rows)), Default::default()).unwrap();
    let m = GpsIfDsbClock::new(&p).unwrap();
    assert!(matches!(
        m.corrections_at(time(), &binding(), &[request("1C")]),
        Err(DsbClockError::MissingPair { .. })
    ));
    let p = product(1.0, 2.0, 3.0, false);
    let mut rows = String::new();
    for e in p.entries() {
        rows += &row(
            &e.observation_1,
            e.observation_2.as_ref().unwrap(),
            e.metres(None).unwrap().0,
            e.metres(None).unwrap().1,
            false,
        );
    }
    rows += &row("C1C", "C2W", 999.0, 0.01, false);
    let p = read_bias_sinex(Cursor::new(file(&rows)), Default::default()).unwrap();
    let s = GpsIfDsbClock::new(&p)
        .unwrap()
        .corrections_at(time(), &binding(), &[request("1C")])
        .unwrap();
    close(
        s.provenance[&("G01".parse().unwrap(), "1C".into())].code_add_m,
        -c::IFLC_F2 - 2.0,
    );
}

#[test]
fn conflicting_reverse_pair_is_ambiguous_even_if_its_value_agrees() {
    let rows = row("C1W", "C2W", 1.0, 0.01, false) + &row("C2W", "C1W", -1.0, 0.01, false);
    let p = read_bias_sinex(Cursor::new(file(&rows)), Default::default()).unwrap();
    assert!(matches!(
        GpsIfDsbClock::new(&p)
            .unwrap()
            .corrections_at(time(), &binding(), &[request("1W")]),
        Err(DsbClockError::AmbiguousPair { .. })
    ));
}

#[test]
fn wrong_product_modes_unsupported_requests_and_duplicate_bindings_are_explicit() {
    let p = product(1.0, 2.0, 3.0, false);
    let m = GpsIfDsbClock::new(&p).unwrap();
    let mut q = request("1C");
    q.satellite = "E01".parse().unwrap();
    assert!(matches!(
        m.corrections_at(time(), &binding(), &[q]),
        Err(DsbClockError::UnsupportedRequest)
    ));
    assert!(matches!(
        m.corrections_at(time(), &binding(), &[request("1X")]),
        Err(DsbClockError::UnsupportedRequest)
    ));
    assert!(matches!(
        m.corrections_at(time(), &binding(), &[request("1C"), request("1C")]),
        Err(DsbClockError::InvalidBinding)
    ));
    let mut b = binding();
    b.clock_reference.clear();
    assert!(matches!(
        m.corrections_at(time(), &b, &[request("1C")]),
        Err(DsbClockError::InvalidBinding)
    ));
    let text = file(&row("C1W", "C2W", 1.0, 0.01, false))
        .replace(" BIAS_MODE RELATIVE", " BIAS_MODE ABSOLUTE");
    let p = read_bias_sinex(Cursor::new(text), Default::default()).unwrap();
    assert!(matches!(
        GpsIfDsbClock::new(&p),
        Err(DsbClockError::NotRelativeProduct)
    ));
}

#[test]
fn unknown_source_correlations_use_a_sigma_bound_and_receiver_rows_do_not_set_gauge() {
    let rows = row("C1W", "C2W", 1.0, 0.01, false)
        + &row("C1C", "C1W", 2.0, 0.02, false)
        + " DSB G080 G01 ABCD C1W C2W 2024:039:00600 2024:039:02520 m 999.0 0.01\n";
    let p = read_bias_sinex(Cursor::new(file(&rows)), Default::default()).unwrap();
    let before = p.entries().to_vec();
    let snap = GpsIfDsbClock::new(&p)
        .unwrap()
        .corrections_at(time(), &binding(), &[request("1C")])
        .unwrap();
    let d = &snap.provenance[&("G01".parse().unwrap(), "1C".into())];
    close(d.sigma_upper_bound_m, c::IFLC_F2.abs() * 0.01 + 0.02);
    assert!(d.sigma_upper_bound_m > ((c::IFLC_F2 * 0.01).powi(2) + 0.02_f64.powi(2)).sqrt());
    assert!(d.source_terms.iter().all(|t| t.entry.station.is_none()));
    assert_eq!(p.entries(), before);
}
