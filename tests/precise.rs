use gnss_rust::precise::{
    ClockSource, InterpolationConfig, InterpolationError, OrbitReferencePoint, PreciseProducts,
    ProductError, precise_clock_relativistic_correction, read_clk, read_sp3,
};
use gnss_rust::time::TimeScale;
use gnss_rust::{EcefCoord, GnssTime};
use std::io::Cursor;

const SP3: &str = include_str!("fixtures/synthetic_precise.sp3");
const TWO: &str = include_str!("fixtures/synthetic_precise_two.sp3");
const CLK: &str = include_str!("fixtures/synthetic_precise.clk");
fn products(sp3: &str, clk: Option<&str>) -> PreciseProducts {
    PreciseProducts::new(
        read_sp3(Cursor::new(sp3)).unwrap(),
        clk.map(|s| read_clk(Cursor::new(s)).unwrap()),
    )
}
fn time(offset: f64) -> GnssTime {
    GnssTime::new(2300, 345600.0 + offset).unwrap()
}
fn query(
    p: &PreciseProducts,
    offset: f64,
) -> Result<gnss_rust::precise::PreciseState, InterpolationError> {
    p.interpolate(
        "G01".parse().unwrap(),
        time(offset),
        InterpolationConfig::default(),
    )
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, difference {:.4e}",
        (a - b).abs()
    );
}
fn patch_record(input: &str, index: usize, start: usize, width: usize, value: &str) -> String {
    assert_eq!(value.len(), width);
    let mut seen = 0;
    let mut out = String::new();
    for line in input.lines() {
        let mut line = line.to_owned();
        if line.starts_with("PG01") {
            if seen == index {
                if line.len() < start + width {
                    line.push_str(&" ".repeat(start + width - line.len()));
                }
                line.replace_range(start..start + width, value);
            }
            seen += 1;
        }
        out.push_str(&line);
        out.push('\n');
    }
    out
}

#[test]
fn sp3_and_clk_preserve_metadata_units_zero_clock_and_distinct_grids() {
    let p = products(SP3, Some(CLK));
    let h = p.orbit().header();
    assert_eq!(h.version, 'd');
    assert_eq!(h.declared_epochs, 12);
    assert_eq!(h.start_time, time(0.0));
    assert_eq!(h.interval_s, 300.0);
    assert_eq!(h.coordinate_system, "IGS20");
    assert_eq!(h.agency, "TEST");
    assert_eq!(h.orbit_type, "FIT");
    assert_eq!(h.time_scale, TimeScale::Gps);
    assert_eq!(h.reference_point, OrbitReferencePoint::AntennaPhaseCenter);
    assert_eq!(h.accuracy_exponents, vec![Some(5), Some(6)]);
    close(h.position_sigma_base, 1.25, 0.0);
    close(h.clock_sigma_base, 1.025, 0.0);
    let rows = &p.orbit().samples()[&"G01".parse().unwrap()];
    assert_eq!(rows.len(), 12);
    assert_eq!(
        rows[0].position_m,
        Some(EcefCoord::new(20200000.0, -14100000.0, 21700000.0))
    );
    close(rows[0].clock_bias_s.unwrap(), 2e-5, 1e-20);
    assert_eq!(
        rows[0].accuracy_exponents,
        [Some(5), Some(6), Some(7), Some(8)]
    );
    assert_eq!(
        p.orbit().samples()[&"E12".parse().unwrap()][0].clock_bias_s,
        Some(0.0)
    );
    let clocks = p.clocks().unwrap();
    assert!(clocks.header().time_system_declared);
    assert_eq!(clocks.samples()[&"G01".parse().unwrap()].len(), 40);
    assert_eq!(
        clocks.samples()[&"G01".parse().unwrap()][0].sigma_s,
        Some(1e-12)
    );
    assert_eq!(
        query(&p, 150.0).unwrap().clock.unwrap().source,
        ClockSource::RinexClk
    );
    let e12 = p
        .interpolate("E12".parse().unwrap(), time(150.0), Default::default())
        .unwrap();
    assert_eq!(e12.clock.as_ref().unwrap().source, ClockSource::Sp3);
    assert_eq!(e12.clock.unwrap().bias_s, 0.0);
}

#[test]
fn precise_states_match_native_neville_on_matching_grids_and_expose_mixed_clock_difference() {
    let sp3 = products(SP3, None);
    let clk = products(SP3, Some(CLK));
    let two = products(TWO, None);
    let mut matched = 0;
    let mut mixed_differences = 0;
    for line in include_str!("fixtures/upstream_precise.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split(',').collect();
        let mode: usize = fields[0].parse().unwrap();
        let satellite = fields[1].parse().unwrap();
        let week = fields[2].parse().unwrap();
        let v: Vec<f64> = fields[3..].iter().map(|v| v.parse().unwrap()).collect();
        let p = match mode {
            0 => &sp3,
            2 => &two,
            _ => &clk,
        };
        let state = p
            .interpolate(
                satellite,
                GnssTime::new(week, v[0]).unwrap(),
                Default::default(),
            )
            .unwrap();
        for (actual, expected) in [state.position_m.x, state.position_m.y, state.position_m.z]
            .into_iter()
            .zip(&v[1..4])
        {
            close(actual, *expected, 1e-6);
        }
        for (actual, expected) in state.velocity_m_per_s.into_iter().zip(&v[4..7]) {
            close(actual, *expected, 2e-7);
        }
        close(
            precise_clock_relativistic_correction(state.position_m, state.velocity_m_per_s)
                .unwrap(),
            v[10],
            1e-16,
        );
        assert_eq!(state.clock.is_some(), v[7] != 0.0);
        let clock = state.clock.unwrap();
        if mode == 3 && fields[1] == "G01" {
            if (clock.bias_s - v[8]).abs() > 1e-7 {
                mixed_differences += 1;
            }
            continue;
        }
        close(clock.bias_s, v[8], 1e-16);
        close(clock.drift_s_per_s, v[9], 1e-16);
        assert!(!state.orbit_support.extrapolated);
        assert_eq!(
            state.orbit_support.sample_times.len(),
            if mode == 2 { 2 } else { 10 }
        );
        matched += 1;
    }
    assert_eq!(matched, 55); // 46 matching-grid cases plus 9 unchanged E12 cases.
    assert!(mixed_differences > 0);
}

#[test]
fn two_point_interpolation_has_analytic_linear_position_and_velocity() {
    let p = products(TWO, None);
    let rows = &p.orbit().samples()[&"G01".parse().unwrap()];
    let a = rows[0].position_m.unwrap();
    let b = rows[1].position_m.unwrap();
    let state = query(&p, 150.0).unwrap();
    for (i, (a, b)) in [a.x, a.y, a.z].into_iter().zip([b.x, b.y, b.z]).enumerate() {
        close(
            [state.position_m.x, state.position_m.y, state.position_m.z][i],
            (a + b) / 2.0,
            1e-8,
        );
        close(state.velocity_m_per_s[i], (b - a) / 300.0, 1e-8);
    }
    let corrected =
        precise_clock_relativistic_correction(EcefCoord::new(2e7, 0.0, 0.0), [1000.0, 0.0, 0.0])
            .unwrap();
    close(
        corrected,
        -4e10 / gnss_rust::constants::SPEED_OF_LIGHT.powi(2),
        1e-22,
    );
    assert!(
        precise_clock_relativistic_correction(EcefCoord::new(f64::NAN, 0.0, 0.0), [0.0; 3])
            .is_err()
    );
}

#[test]
fn coverage_limits_week_boundaries_and_explicit_extrapolation() {
    let p = products(SP3, None);
    assert_eq!(query(&p, -0.125), Err(InterpolationError::OutOfRange));
    assert_eq!(query(&p, 3300.125), Err(InterpolationError::OutOfRange));
    assert_eq!(query(&p, 604800.0), Err(InterpolationError::OutOfRange));
    let cfg = InterpolationConfig {
        max_extrapolation_s: 30.0,
        ..Default::default()
    };
    for offset in [-30.0, 3330.0] {
        let s = p
            .interpolate("G01".parse().unwrap(), time(offset), cfg)
            .unwrap();
        assert!(s.orbit_support.extrapolated);
        assert!(s.clock.unwrap().support.extrapolated);
    }
    assert_eq!(
        p.interpolate("G01".parse().unwrap(), time(-30.001), cfg),
        Err(InterpolationError::OutOfRange)
    );
    assert_eq!(
        p.interpolate("G99".parse().unwrap(), time(0.0), Default::default()),
        Err(InterpolationError::MissingSatellite)
    );
    for cfg in [
        InterpolationConfig {
            max_samples: 1,
            ..Default::default()
        },
        InterpolationConfig {
            max_clock_gap_s: f64::NAN,
            ..Default::default()
        },
        InterpolationConfig {
            max_extrapolation_s: 901.0,
            ..Default::default()
        },
    ] {
        assert_eq!(
            p.interpolate("G01".parse().unwrap(), time(0.0), cfg),
            Err(InterpolationError::InvalidConfiguration)
        );
    }
}

#[test]
fn orbit_missing_values_quality_boundaries_and_gaps_do_not_get_bridged() {
    let missing = patch_record(
        SP3,
        5,
        4,
        42,
        &format!("{:14.6}{:14.6}{:14.6}", 0.0, 0.0, 0.0),
    );
    let p = products(&missing, None);
    assert!(
        p.orbit().samples()[&"G01".parse().unwrap()][5]
            .position_m
            .is_none()
    );
    for offset in [1350.0, 1500.0, 1650.0] {
        assert_eq!(
            query(&p, offset),
            Err(InterpolationError::GapOrDiscontinuity)
        );
    }
    assert_eq!(
        query(&p, 1200.0).unwrap().orbit_support.sample_times.len(),
        5
    );
    assert!(query(&p, 2100.0).is_ok());
    let maneuver = patch_record(SP3, 5, 78, 1, "M");
    let p = products(&maneuver, None);
    assert_eq!(
        query(&p, 1350.0),
        Err(InterpolationError::GapOrDiscontinuity)
    );
    let new_arc = query(&p, 1500.0).unwrap();
    assert!(
        new_arc
            .orbit_support
            .sample_times
            .iter()
            .all(|t| *t >= time(1500.0))
    );
    let predicted = patch_record(SP3, 5, 79, 1, "P");
    let p = products(&predicted, None);
    assert_eq!(
        query(&p, 1500.0),
        Err(InterpolationError::GapOrDiscontinuity)
    );
    let s = p
        .interpolate(
            "G01".parse().unwrap(),
            time(1500.0),
            InterpolationConfig {
                allow_predicted: true,
                ..Default::default()
            },
        )
        .unwrap();
    assert!(s.orbit_support.includes_predicted);
    let p = products(SP3, None);
    assert_eq!(
        p.interpolate(
            "G01".parse().unwrap(),
            time(150.0),
            InterpolationConfig {
                max_orbit_gap_s: 299.0,
                ..Default::default()
            }
        ),
        Err(InterpolationError::GapOrDiscontinuity)
    );
    // Even opt-in outer extrapolation cannot cross an internal missing epoch.
    let p = products(&missing, None);
    assert_eq!(
        p.interpolate(
            "G01".parse().unwrap(),
            time(1500.0),
            InterpolationConfig {
                max_extrapolation_s: 900.0,
                ..Default::default()
            }
        ),
        Err(InterpolationError::GapOrDiscontinuity)
    );
}

#[test]
fn missing_clock_and_clock_events_keep_orbit_but_never_return_zero_as_substitute() {
    let missing = patch_record(SP3, 5, 46, 14, " 999999.999999");
    let p = products(&missing, None);
    let s = query(&p, 1500.0).unwrap();
    assert!(s.clock.is_none());
    assert_eq!(
        s.clock_unavailability,
        Some(InterpolationError::GapOrDiscontinuity)
    );
    let event = patch_record(SP3, 5, 74, 1, "E");
    let p = products(&event, None);
    assert!(query(&p, 1350.0).unwrap().clock.is_none());
    assert!(
        query(&p, 1500.0)
            .unwrap()
            .clock
            .unwrap()
            .support
            .sample_times
            .iter()
            .all(|t| *t >= time(1500.0))
    );
    let predicted = patch_record(SP3, 5, 75, 1, "P");
    let p = products(&predicted, None);
    assert!(query(&p, 1500.0).unwrap().clock.is_none());
    let s = p
        .interpolate(
            "G01".parse().unwrap(),
            time(1500.0),
            InterpolationConfig {
                allow_predicted: true,
                ..Default::default()
            },
        )
        .unwrap();
    assert!(s.clock.unwrap().support.includes_predicted);
}

#[test]
fn clk_gaps_and_short_coverage_do_not_fall_back_to_a_different_sp3_datum() {
    // Use the first AS record only; it cannot provide a rate/interpolation arc.
    let mut one = CLK
        .lines()
        .take_while(|l| !l.starts_with("AS"))
        .collect::<Vec<_>>()
        .join("\n");
    one.push('\n');
    one.push_str(CLK.lines().find(|l| l.starts_with("AS")).unwrap());
    one.push('\n');
    let p = products(SP3, Some(&one));
    let s = query(&p, 0.0).unwrap();
    assert!(s.clock.is_none());
    assert_eq!(
        s.clock_unavailability,
        Some(InterpolationError::InsufficientSamples)
    );
    assert_eq!(s.clock_source, ClockSource::RinexClk);
    assert!(query(&p, 150.0).unwrap().clock.is_none());
    let gapped = CLK
        .lines()
        .filter(|l| {
            !l.starts_with("AS")
                || l.split_whitespace().nth(6).unwrap().parse::<u8>().unwrap() < 10
                || l.split_whitespace().nth(6).unwrap().parse::<u8>().unwrap() > 40
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let p = products(SP3, Some(&gapped));
    let s = query(&p, 1500.0).unwrap();
    assert!(s.clock.is_none());
    assert_eq!(
        s.clock_unavailability,
        Some(InterpolationError::GapOrDiscontinuity)
    );
}

#[test]
fn time_systems_are_converted_and_unknown_scales_are_rejected() {
    for (label, scale, offset) in [
        ("UTC", TimeScale::Utc, 18.0),
        ("BDT", TimeScale::BeiDou, 14.0),
        ("GAL", TimeScale::Galileo, 0.0),
        ("QZS", TimeScale::Qzss, 0.0),
    ] {
        let sp3 = SP3.replace("%c M  cc GPS", &format!("%c M  cc {label}"));
        let data = read_sp3(Cursor::new(sp3)).unwrap();
        assert_eq!(data.header().time_scale, scale);
        assert_eq!(data.header().start_time, time(offset));
        let clocks = read_clk(Cursor::new(CLK.replace("GPS ", &format!("{label} ")))).unwrap();
        assert_eq!(clocks.header().time_scale, scale);
        assert_eq!(
            clocks.samples()[&"G01".parse().unwrap()][0].time,
            time(offset)
        );
    }
    for label in ["TAI", "GLO", "xxx"] {
        assert!(matches!(
            read_sp3(Cursor::new(
                SP3.replace("%c M  cc GPS", &format!("%c M  cc {label}"))
            )),
            Err(ProductError::Unsupported { .. })
        ));
    }
    let no_id = CLK
        .lines()
        .filter(|l| !l.ends_with("TIME SYSTEM ID"))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let c = read_clk(Cursor::new(no_id)).unwrap();
    assert_eq!(c.header().time_scale, TimeScale::Gps);
    assert!(!c.header().time_system_declared);
}

#[test]
fn clock_304_header_continuations_fortran_exponents_and_all_six_values() {
    let header = |value: &str, label: &str| format!("{value:<65}{label}\n");
    let text = header("     3.04           C", "RINEX VERSION / TYPE")
        + &header("GPS", "TIME SYSTEM ID")
        + &header("", "END OF HEADER")
        + "AS G01 2024 02 08 00 00 0.00000000 6 0.0D+00 0.0D+00\n   -2.0D-11 3.0D-12 -4.0D-15 5.0D-16\n";
    let c = read_clk(Cursor::new(text)).unwrap();
    let row = &c.samples()[&"G01".parse().unwrap()][0];
    assert_eq!(row.bias_s, 0.0);
    assert_eq!(row.sigma_s, Some(0.0));
    assert_eq!(row.drift_s_per_s, Some(-2e-11));
    assert_eq!(row.drift_sigma_s_per_s, Some(3e-12));
    assert_eq!(row.rate_s_per_s2, Some(-4e-15));
    assert_eq!(row.rate_sigma_s_per_s2, Some(5e-16));
}

#[test]
fn malformed_incomplete_or_unsupported_products_fail_instead_of_losing_records() {
    for text in [
        SP3.replacen("#dP", "#aP", 1),
        SP3.replacen("#dP", "#dV", 1),
        SP3.replace("EOF", ""),
        SP3.replacen("      12 ORBIT", "      13 ORBIT", 1),
        SP3.replacen("G01E12", "G01G01", 1),
        SP3.replacen("PG01", "VG01", 1),
        patch_record(SP3, 0, 4, 14, "           NaN"),
        patch_record(SP3, 0, 4, 14, "         1e308"),
        patch_record(SP3, 0, 79, 1, "X"),
        SP3.to_owned() + "PG01\n",
    ] {
        assert!(read_sp3(Cursor::new(text)).is_err());
    }
    let first = SP3.lines().find(|l| l.starts_with("PG01")).unwrap();
    assert!(read_sp3(Cursor::new(SP3.replacen(first, "", 1))).is_err());
    assert!(
        read_sp3(Cursor::new(SP3.replacen(
            first,
            &format!("{first}\n{first}"),
            1
        )))
        .is_err()
    );
    for text in [
        CLK.replacen("AS G01", "AR G01", 1),
        CLK.replacen("  2 ", "  7 ", 1),
        CLK.replacen("1.000000000000E-12", "-1.000000000000E-12", 1),
        CLK.replacen("5.000000000000e-05", "NaN", 1),
    ] {
        assert!(read_clk(Cursor::new(text)).is_err());
    }
    let first = CLK.lines().find(|l| l.starts_with("AS")).unwrap();
    assert!(
        read_clk(Cursor::new(CLK.replacen(
            first,
            &format!("{first}\n{first}"),
            1
        )))
        .is_err()
    );
    assert!(read_clk(Cursor::new(CLK.replace("END OF HEADER", "COMMENT"))).is_err());
}

#[test]
fn interpolation_spans_a_gps_week_boundary_with_fractional_queries() {
    let mut epoch = 0;
    let text = SP3
        .lines()
        .map(|line| {
            if line.starts_with("#dP") {
                return line.replacen("2024 02 08 00 00", "2024 02 10 23 45", 1);
            }
            if line.starts_with("##") {
                return line.replacen("345600.00000000", "603900.00000000", 1);
            }
            if line.starts_with('*') {
                let minutes = 45 + 5 * epoch;
                epoch += 1;
                return if minutes < 60 {
                    format!("*  2024 02 10 23 {minutes:02}  0.00000000")
                } else {
                    format!("*  2024 02 11 00 {:02}  0.00000000", minutes - 60)
                };
            }
            line.to_owned()
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let p = products(&text, None);
    assert_eq!(
        p.orbit().header().start_time,
        GnssTime::new(2300, 603900.0).unwrap()
    );
    let original = products(SP3, None);
    for offset in [899.875, 900.0, 900.125, 1350.0] {
        let s = p
            .interpolate(
                "G01".parse().unwrap(),
                GnssTime::new(2300, 603900.0 + offset).unwrap(),
                Default::default(),
            )
            .unwrap();
        let expected = query(&original, offset).unwrap();
        assert_eq!(s.position_m, expected.position_m);
        assert_eq!(s.velocity_m_per_s, expected.velocity_m_per_s);
        assert_eq!(s.clock.unwrap().bias_s, expected.clock.unwrap().bias_s);
        assert!(
            s.orbit_support
                .sample_times
                .iter()
                .any(|t| t.week() == 2300)
        );
        assert!(
            s.orbit_support
                .sample_times
                .iter()
                .any(|t| t.week() == 2301)
        );
    }
}

#[test]
fn satellite_list_continuations_are_supported_in_sp3_c_and_d() {
    let first_list = (1..=17).map(|i| format!("G{i:02}")).collect::<String>();
    let mut records = String::new();
    for line in TWO.lines() {
        if line.starts_with("++") {
            records.push_str(&format!(
                "++       {}\n++       {:3}{}\n",
                "  0".repeat(17),
                0,
                "   ".repeat(16)
            ));
        } else if line.starts_with('+') {
            records.push_str(&format!(
                "+   18   {first_list}\n+        G18{}\n",
                "  0".repeat(16)
            ));
        } else if line.starts_with("PG01") {
            for i in 1..=18 {
                records.push_str(&line.replacen("PG01", &format!("PG{i:02}"), 1));
                records.push('\n');
            }
        } else if !line.starts_with("PE12") && !line.starts_with("/*") {
            records.push_str(line);
            records.push('\n');
        }
    }
    for version in ["#cP", "#dP"] {
        let p = products(&records.replacen("#dP", version, 1), None);
        assert_eq!(p.orbit().header().satellites.len(), 18);
        assert_eq!(
            p.orbit().header().reference_point,
            OrbitReferencePoint::CentreOfMass
        );
        let s = p
            .interpolate("G18".parse().unwrap(), time(150.0), Default::default())
            .unwrap();
        assert_eq!(s.orbit_support.sample_times.len(), 2);
    }
}
