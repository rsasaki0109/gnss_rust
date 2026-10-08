use gnss_rust::precise::{ClockSource, OrbitReferencePoint, PreciseProducts, read_clk, read_sp3};
use gnss_rust::precise_spp::{
    PreciseCodeBias, PreciseCodeBiases, PreciseSppConfig, PreciseSppError, PreciseSppRejection,
    solve_epoch_precise,
};
use gnss_rust::precise_transmit::PreciseTransmitMethod;
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::spp::SppConfig;
use gnss_rust::{EcefCoord, Error, GnssTime, SolutionStatus};
use std::io::Cursor;

fn products() -> PreciseProducts {
    PreciseProducts::new(
        read_sp3(Cursor::new(include_str!(
            "fixtures/synthetic_precise_spp.sp3"
        )))
        .unwrap(),
        Some(
            read_clk(Cursor::new(include_str!(
                "fixtures/synthetic_precise_spp.clk"
            )))
            .unwrap(),
        ),
    )
}
fn biases() -> Vec<PreciseCodeBias> {
    include_str!("fixtures/synthetic_precise_spp_bias.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| {
            let v: Vec<_> = l.split(',').collect();
            PreciseCodeBias {
                satellite: v[0].parse().unwrap(),
                tracking_code: v[1].to_owned(),
                clock_source: ClockSource::RinexClk,
                valid_from: GnssTime::new(v[3].parse().unwrap(), v[4].parse().unwrap()).unwrap(),
                valid_until: GnssTime::new(v[5].parse().unwrap(), v[6].parse().unwrap()).unwrap(),
                correction_m: v[7].parse().unwrap(),
                clock_reference: v[8].to_owned(),
            }
        })
        .collect()
}
fn config(exact: bool) -> PreciseSppConfig {
    PreciseSppConfig {
        spp: SppConfig {
            use_ionosphere: false,
            ..Default::default()
        },
        transmit: gnss_rust::precise_transmit::PreciseTransmitConfig {
            method: if exact {
                PreciseTransmitMethod::IteratedEmission
            } else {
                PreciseTransmitMethod::ReceptionTaylor
            },
            ..Default::default()
        },
        expected_coordinate_system: "IGS20".to_owned(),
        expected_reference_point: OrbitReferencePoint::AntennaPhaseCenter,
        clock_reference: "synthetic-precise-clock-v1".to_owned(),
        ionosphere: None,
    }
}
fn epochs(exact: bool) -> Vec<gnss_rust::observation::ObservationEpoch> {
    RinexObservationReader::new(Cursor::new(if exact {
        include_str!("fixtures/synthetic_precise_spp_emission.obs")
    } else {
        include_str!("fixtures/synthetic_precise_spp_taylor.obs")
    }))
    .unwrap()
    .collect::<Result<_, _>>()
    .unwrap()
}
fn truth() -> EcefCoord {
    let v: Vec<f64> = include_str!("fixtures/synthetic_precise_spp_truth.csv")
        .lines()
        .find(|l| !l.starts_with('#'))
        .unwrap()
        .split(',')
        .map(|v| v.parse().unwrap())
        .collect();
    EcefCoord::new(v[0], v[1], v[2])
}
fn distance(p: EcefCoord, q: EcefCoord) -> f64 {
    (p.x - q.x).hypot(p.y - q.y).hypot(p.z - q.z)
}

#[test]
fn native_generated_precise_code_epochs_recover_truth_for_both_light_time_methods() {
    let p = products();
    let table = PreciseCodeBiases::new(biases()).unwrap();
    for exact in [false, true] {
        let mut max = 0.0_f64;
        for epoch in epochs(exact) {
            let result = solve_epoch_precise(&epoch, &p, &table, &config(exact)).unwrap();
            let error = distance(result.position.position_ecef.unwrap(), truth());
            max = max.max(error);
            assert!(error < 0.001, "{error}");
            assert!((result.position.receiver_clock_bias - 2e-6).abs() < 1e-12);
            assert_eq!(result.position.status, SolutionStatus::Spp);
            assert_eq!(result.position.num_satellites, 9);
            assert_eq!(result.measurements.len(), 9);
            assert!(result.residual_rms_m < 1e-4);
            assert!(result.rejected.is_empty());
            assert_eq!(result.coordinate_system, "IGS20");
            for m in result.measurements {
                assert_eq!(m.applied_bias.clock_reference, "synthetic-precise-clock-v1");
                assert_eq!(
                    m.transmit.evaluated_state.clock_source,
                    ClockSource::RinexClk
                );
                assert_eq!(
                    m.corrected_pseudorange_m,
                    m.raw_pseudorange_m + m.applied_bias.correction_m
                );
                assert_eq!(
                    m.transmit.corrected_clock_bias_s,
                    m.transmit.published_clock_bias_s + m.transmit.periodic_relativity_s
                );
                assert!(m.transmit.emission_time < epoch.time);
            }
        }
        eprintln!("precise code SPP exact={exact}: maximum 3-D error {max:.9} m");
    }
}

#[test]
fn missing_expired_or_mismatched_biases_are_reported_per_satellite() {
    let p = products();
    let epoch = epochs(false).remove(0);
    let cfg = config(false);
    for defect in 0..4 {
        let mut entries = biases();
        let sat = entries[0].satellite;
        let expected = match defect {
            0 => {
                entries.remove(0);
                PreciseSppRejection::MissingCodeBias
            }
            1 => {
                entries[0].valid_until = entries[0].valid_from.checked_add_seconds(-1.0).unwrap();
                entries[0].valid_from = entries[0].valid_until;
                PreciseSppRejection::ExpiredCodeBias
            }
            2 => {
                entries[0].clock_source = ClockSource::Sp3;
                PreciseSppRejection::BiasClockSourceMismatch
            }
            _ => {
                entries[0].clock_reference = "another-clock-datum".to_owned();
                PreciseSppRejection::BiasClockReferenceMismatch
            }
        };
        let result =
            solve_epoch_precise(&epoch, &p, &PreciseCodeBiases::new(entries).unwrap(), &cfg)
                .unwrap();
        assert_eq!(result.position.num_satellites, 8);
        assert!(result.rejected.contains(&(sat, expected)));
        assert!(!result.measurements.iter().any(|m| m.satellite == sat));
        assert!(distance(result.position.position_ecef.unwrap(), truth()) < 0.001);
    }
    assert_eq!(
        solve_epoch_precise(
            &epoch,
            &p,
            &PreciseCodeBiases::new(Vec::new()).unwrap(),
            &cfg
        ),
        Err(PreciseSppError::Gnss(Error::InsufficientObservations))
    );
}

#[test]
fn missing_precise_clock_is_rejected_even_when_sp3_has_a_clock_datum() {
    let orbit = include_str!("fixtures/synthetic_precise_spp.sp3")
        .replace(" 999999.999999", "     20.000000");
    let mut one = include_str!("fixtures/synthetic_precise_spp.clk")
        .lines()
        .take_while(|l| !l.starts_with("AS"))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    for line in include_str!("fixtures/synthetic_precise_spp.clk")
        .lines()
        .filter(|l| l.starts_with("AS"))
    {
        if !line.starts_with("AS G01") || !one.contains("AS G01") {
            one.push_str(line);
            one.push('\n');
        }
    }
    let p = PreciseProducts::new(
        read_sp3(Cursor::new(orbit)).unwrap(),
        Some(read_clk(Cursor::new(one)).unwrap()),
    );
    let result = solve_epoch_precise(
        &epochs(false)[0],
        &p,
        &PreciseCodeBiases::new(biases()).unwrap(),
        &config(false),
    )
    .unwrap();
    assert_eq!(result.position.num_satellites, 8);
    assert!(
        result
            .rejected
            .iter()
            .any(|(s, r)| *s == "G01".parse().unwrap()
                && matches!(
                    r,
                    PreciseSppRejection::Transmit(
                        gnss_rust::precise_transmit::PreciseTransmitError::ClockUnavailable {
                            source: ClockSource::RinexClk,
                            ..
                        }
                    )
                ))
    );
}

#[test]
fn declared_zero_bias_is_allowed_but_omitting_real_signal_corrections_changes_the_solution() {
    let mut entries = biases();
    for b in &mut entries {
        b.correction_m = 0.0;
    }
    let result = solve_epoch_precise(
        &epochs(false)[0],
        &products(),
        &PreciseCodeBiases::new(entries).unwrap(),
        &config(false),
    )
    .unwrap();
    assert!(distance(result.position.position_ecef.unwrap(), truth()) > 0.1);
    assert!(result.residual_rms_m > 0.1);
}

#[test]
fn clock_offsets_product_frames_and_undeclared_ionosphere_models_are_not_assumed() {
    let p = products();
    let table = PreciseCodeBiases::new(biases()).unwrap();
    let epoch = epochs(false).remove(0);
    let mut cfg = config(false);
    cfg.expected_coordinate_system = "ITRF-mismatch".to_owned();
    assert_eq!(
        solve_epoch_precise(&epoch, &p, &table, &cfg),
        Err(PreciseSppError::ProductFrameMismatch)
    );
    cfg = config(false);
    cfg.expected_reference_point = OrbitReferencePoint::CentreOfMass;
    assert_eq!(
        solve_epoch_precise(&epoch, &p, &table, &cfg),
        Err(PreciseSppError::ProductFrameMismatch)
    );
    cfg = config(false);
    cfg.spp.use_ionosphere = true;
    assert_eq!(
        solve_epoch_precise(&epoch, &p, &table, &cfg),
        Err(PreciseSppError::InvalidConfiguration)
    );
    let mut bad = epoch.clone();
    bad.receiver_clock_offset_s = Some(1e-6);
    assert_eq!(
        solve_epoch_precise(&bad, &p, &table, &config(false)),
        Err(PreciseSppError::InvalidEpoch)
    );
    bad = epoch;
    bad.flag = 6;
    assert_eq!(
        solve_epoch_precise(&bad, &p, &table, &config(false)),
        Err(PreciseSppError::InvalidEpoch)
    );
}

#[test]
fn invalid_or_overlapping_bias_intervals_and_tracking_names_fail() {
    for defect in 0..6 {
        let mut entries = biases();
        match defect {
            0 => entries.push(entries[0].clone()),
            1 => entries[0].correction_m = f64::NAN,
            2 => entries[0].clock_reference.clear(),
            3 => entries[0].valid_until = entries[0].valid_from.checked_add_seconds(-1.0).unwrap(),
            4 => entries[0].tracking_code = "1c".to_owned(),
            _ => entries[0].tracking_code = "".to_owned(),
        }
        assert!(matches!(
            PreciseCodeBiases::new(entries),
            Err(PreciseSppError::InvalidBiasTable)
        ));
    }
}

#[test]
fn precise_code_spp_can_acquire_without_receiver_seed_and_checks_iteration_failure() {
    let p = products();
    let table = PreciseCodeBiases::new(biases()).unwrap();
    let mut epoch = epochs(false).remove(0);
    epoch.approximate_position = None;
    let result = solve_epoch_precise(&epoch, &p, &table, &config(false)).unwrap();
    assert!(distance(result.position.position_ecef.unwrap(), truth()) < 0.001);
    let mut cfg = config(false);
    cfg.spp.max_iterations = 1;
    assert_eq!(
        solve_epoch_precise(&epoch, &p, &table, &cfg),
        Err(PreciseSppError::Gnss(Error::NonConvergent))
    );
}
