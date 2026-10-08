use gnss_rust::coordinates::geometric_distance;
use gnss_rust::precise::{ClockSource, InterpolationError, PreciseProducts, read_clk, read_sp3};
use gnss_rust::precise_transmit::{
    PreciseTransmitConfig, PreciseTransmitError, PreciseTransmitMethod, evaluate_precise_transmit,
};
use gnss_rust::{EcefCoord, GnssTime, constants::SPEED_OF_LIGHT};
use std::io::Cursor;

const ORBIT: &str = include_str!("fixtures/synthetic_precise.sp3");
fn products(orbit: &str, clk: bool) -> PreciseProducts {
    PreciseProducts::new(
        read_sp3(Cursor::new(orbit)).unwrap(),
        if clk {
            Some(read_clk(Cursor::new(include_str!("fixtures/synthetic_precise.clk"))).unwrap())
        } else {
            None
        },
    )
}
fn receiver() -> EcefCoord {
    EcefCoord::new(-3947481.0649388283, 3431492.9375319984, 3637892.720317735)
}
fn time(offset: f64) -> GnssTime {
    GnssTime::new(2300, 345600.0 + offset).unwrap()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "{a:.17e} != {b:.17e} (difference {:.3e})",
        (a - b).abs()
    );
}
fn flag(index: usize, column: usize, value: char) -> String {
    let mut n = 0;
    ORBIT
        .lines()
        .map(|line| {
            let mut text = line.to_owned();
            if line.starts_with("PG01") {
                if n == index {
                    text.replace_range(column..column + 1, &value.to_string());
                }
                n += 1;
            }
            text
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

#[test]
fn precise_transmit_matches_native_components_and_literal_taylor_statements() {
    let sp3 = products(ORBIT, false);
    let clk = products(ORBIT, true);
    let mut count = 0;
    for line in include_str!("fixtures/upstream_precise_transmit.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let f: Vec<_> = line.split(',').collect();
        let p = if f[0] == "0" { &sp3 } else { &clk };
        let method = if f[1] == "0" {
            PreciseTransmitMethod::ReceptionTaylor
        } else {
            PreciseTransmitMethod::IteratedEmission
        };
        let sat = f[2].parse().unwrap();
        let week = f[3].parse().unwrap();
        let v: Vec<f64> = f[4..].iter().map(|v| v.parse().unwrap()).collect();
        let rx = EcefCoord::new(v[1], v[2], v[3]);
        let r = evaluate_precise_transmit(
            p,
            sat,
            GnssTime::new(week, v[0]).unwrap(),
            rx,
            PreciseTransmitConfig {
                method,
                ..Default::default()
            },
        )
        .unwrap();
        close(r.emission_time.tow(), v[4], 1e-10);
        close(r.light_time_s, v[5], 1e-13);
        assert_eq!(r.iterations, v[6] as usize);
        for (a, b) in r.position_m.to_array().iter().zip(&v[7..10]) {
            close(*a, *b, 1e-6);
        }
        for (a, b) in r.velocity_m_per_s.iter().zip(&v[10..13]) {
            close(*a, *b, 2e-7);
        }
        close(r.published_clock_bias_s, v[13], 1e-16);
        close(r.published_clock_drift_s_per_s, v[14], 1e-16);
        close(r.periodic_relativity_s, v[15], 1e-16);
        close(r.corrected_clock_bias_s, v[16], 1e-16);
        close(geometric_distance(r.position_m, rx).unwrap(), v[17], 1e-6);
        assert_eq!(
            r.corrected_clock_bias_s,
            r.published_clock_bias_s + r.periodic_relativity_s
        );
        assert_eq!(r.emission_state.time, r.emission_time);
        assert_eq!(
            r.evaluated_state.time,
            if method == PreciseTransmitMethod::ReceptionTaylor {
                r.reception_time
            } else {
                r.emission_time
            }
        );
        count += 1;
    }
    assert_eq!(count, 80);
}

#[test]
fn nonzero_motion_clock_drift_and_relativity_are_applied_once_with_the_right_sign() {
    let p = products(ORBIT, true);
    let r = evaluate_precise_transmit(
        &p,
        "G01".parse().unwrap(),
        time(150.0),
        receiver(),
        Default::default(),
    )
    .unwrap();
    let receive = &r.evaluated_state;
    let c = receive.clock.as_ref().unwrap();
    assert_eq!(c.source, ClockSource::RinexClk);
    for ((position, original), velocity) in r
        .position_m
        .to_array()
        .iter()
        .zip(receive.position_m.to_array())
        .zip(receive.velocity_m_per_s)
    {
        close(*position, original - velocity * r.light_time_s, 1e-9);
    }
    close(
        r.published_clock_bias_s,
        c.bias_s - c.drift_s_per_s * r.light_time_s,
        1e-20,
    );
    assert!(r.periodic_relativity_s.abs() * SPEED_OF_LIGHT > 50.0);
    assert!((r.corrected_clock_bias_s - r.published_clock_bias_s).abs() * SPEED_OF_LIGHT > 50.0);
    let e = evaluate_precise_transmit(
        &p,
        "E12".parse().unwrap(),
        time(150.0),
        receiver(),
        Default::default(),
    )
    .unwrap();
    assert_eq!(e.published_clock_bias_s, 0.0);
    assert_ne!(e.corrected_clock_bias_s, 0.0);
    assert_eq!(e.corrected_clock_bias_s, e.periodic_relativity_s);
}

#[test]
fn iterated_emission_is_consistent_and_the_taylor_approximation_difference_is_visible() {
    let p = products(ORBIT, true);
    let a = evaluate_precise_transmit(
        &p,
        "G01".parse().unwrap(),
        time(1350.125),
        receiver(),
        Default::default(),
    )
    .unwrap();
    let b = evaluate_precise_transmit(
        &p,
        "G01".parse().unwrap(),
        time(1350.125),
        receiver(),
        PreciseTransmitConfig {
            method: PreciseTransmitMethod::IteratedEmission,
            ..Default::default()
        },
    )
    .unwrap();
    let r = receiver();
    let x = b.position_m;
    let geometric = (x.x - r.x).hypot(x.y - r.y).hypot(x.z - r.z) / SPEED_OF_LIGHT;
    close(b.light_time_s, geometric, 1e-11);
    assert!(b.iterations >= 2);
    let difference = (a.position_m.x - b.position_m.x)
        .hypot(a.position_m.y - b.position_m.y)
        .hypot(a.position_m.z - b.position_m.z);
    assert!(difference > 1e-4 && difference < 0.01, "{difference}");
    let failure = evaluate_precise_transmit(
        &p,
        "G01".parse().unwrap(),
        time(1350.125),
        receiver(),
        PreciseTransmitConfig {
            method: PreciseTransmitMethod::IteratedEmission,
            max_iterations: 1,
            ..Default::default()
        },
    );
    assert_eq!(failure, Err(PreciseTransmitError::NonConvergent));
}

#[test]
fn emission_coverage_is_required_and_quality_boundaries_are_not_projected_across() {
    for method in [
        PreciseTransmitMethod::ReceptionTaylor,
        PreciseTransmitMethod::IteratedEmission,
    ] {
        let p = products(ORBIT, true);
        assert_eq!(
            evaluate_precise_transmit(
                &p,
                "G01".parse().unwrap(),
                time(0.0),
                receiver(),
                PreciseTransmitConfig {
                    method,
                    ..Default::default()
                }
            ),
            Err(PreciseTransmitError::Interpolation(
                InterpolationError::OutOfRange
            ))
        );
        for (column, value) in [(78, 'M'), (74, 'E')] {
            let p = products(&flag(3, column, value), false);
            let r = evaluate_precise_transmit(
                &p,
                "G01".parse().unwrap(),
                time(900.0),
                receiver(),
                PreciseTransmitConfig {
                    method,
                    ..Default::default()
                },
            );
            assert_eq!(r, Err(PreciseTransmitError::QualityBoundary));
        }
    }
}

#[test]
fn missing_selected_clock_bad_geometry_and_invalid_configuration_are_errors() {
    let one = include_str!("fixtures/synthetic_precise.clk")
        .lines()
        .take_while(|l| !l.starts_with("AS"))
        .chain(std::iter::once("AS G01 2024 02 08 00 00 0.00000000 1 0.0"))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let p = PreciseProducts::new(
        read_sp3(Cursor::new(ORBIT)).unwrap(),
        Some(read_clk(Cursor::new(one)).unwrap()),
    );
    assert!(matches!(
        evaluate_precise_transmit(
            &p,
            "G01".parse().unwrap(),
            time(150.0),
            receiver(),
            Default::default()
        ),
        Err(PreciseTransmitError::ClockUnavailable {
            source: ClockSource::RinexClk,
            ..
        })
    ));
    let p = products(ORBIT, false);
    assert_eq!(
        evaluate_precise_transmit(
            &p,
            "G01".parse().unwrap(),
            time(150.0),
            EcefCoord::new(f64::NAN, 0.0, 0.0),
            Default::default()
        ),
        Err(PreciseTransmitError::InvalidGeometry)
    );
    for config in [
        PreciseTransmitConfig {
            max_iterations: 0,
            ..Default::default()
        },
        PreciseTransmitConfig {
            convergence_s: f64::NAN,
            ..Default::default()
        },
        PreciseTransmitConfig {
            max_light_time_s: 1.01,
            ..Default::default()
        },
    ] {
        assert_eq!(
            evaluate_precise_transmit(&p, "G01".parse().unwrap(), time(150.0), receiver(), config),
            Err(PreciseTransmitError::InvalidConfiguration)
        );
    }
    assert_eq!(
        evaluate_precise_transmit(
            &p,
            "G01".parse().unwrap(),
            time(150.0),
            receiver(),
            PreciseTransmitConfig {
                max_light_time_s: 0.01,
                ..Default::default()
            }
        ),
        Err(PreciseTransmitError::InvalidGeometry)
    );
}

#[test]
fn receive_and_emission_times_can_cross_a_gps_week_without_stale_modulo_admission() {
    let mut n = 0;
    let text = ORBIT
        .lines()
        .map(|line| {
            if line.starts_with("#dP") {
                return line.replacen("2024 02 08 00 00", "2024 02 10 23 45", 1);
            }
            if line.starts_with("##") {
                return line.replacen("345600.00000000", "603900.00000000", 1);
            }
            if line.starts_with('*') {
                let minute = 45 + n * 5;
                n += 1;
                return if minute < 60 {
                    format!("*  2024 02 10 23 {minute:02}  0.00000000")
                } else {
                    format!("*  2024 02 11 00 {:02}  0.00000000", minute - 60)
                };
            }
            line.to_owned()
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let shifted = products(&text, false);
    let original = products(ORBIT, false);
    for method in [
        PreciseTransmitMethod::ReceptionTaylor,
        PreciseTransmitMethod::IteratedEmission,
    ] {
        let cfg = PreciseTransmitConfig {
            method,
            ..Default::default()
        };
        let s = evaluate_precise_transmit(
            &shifted,
            "G01".parse().unwrap(),
            GnssTime::new(2301, 0.005).unwrap(),
            receiver(),
            cfg,
        )
        .unwrap();
        let expected = evaluate_precise_transmit(
            &original,
            "G01".parse().unwrap(),
            time(900.005),
            receiver(),
            cfg,
        )
        .unwrap();
        assert_eq!(s.emission_time.week(), 2300);
        assert!(s.emission_time.tow() > 604799.0);
        close(
            s.reception_time.difference_seconds(s.emission_time),
            s.light_time_s,
            1e-10,
        );
        for (a, b) in s
            .position_m
            .to_array()
            .into_iter()
            .zip(expected.position_m.to_array())
        {
            close(a, b, 1e-6);
        }
        assert!(
            evaluate_precise_transmit(
                &shifted,
                "G01".parse().unwrap(),
                GnssTime::new(2302, 0.005).unwrap(),
                receiver(),
                cfg
            )
            .is_err()
        );
    }
}
