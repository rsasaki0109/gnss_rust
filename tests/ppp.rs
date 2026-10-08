use gnss_rust::bias::{BiasSignalRequest, OsbBinding, read_bias_sinex};
use gnss_rust::dual_frequency::*;
use gnss_rust::filter::{FilterError, measurement_update};
use gnss_rust::phase_arc::{ArcResetReason, PhaseArcError};
use gnss_rust::ppp::*;
use gnss_rust::precise::{ClockSource, OrbitReferencePoint, PreciseProducts, read_clk, read_sp3};
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::windup::{PppWindupModel, approximate_sun_position_ecef, nominal_yaw_windup};
use gnss_rust::{EcefCoord, Error, GnssTime, SatelliteId, SolutionStatus, constants as c};
use std::{collections::BTreeMap, io::Cursor};

fn tide_config() -> PppConfig {
    let mut cfg = satellite_config();
    cfg.receiver_reference = gnss_rust::tides::ReceiverPositionReference::NominalMarker;
    cfg.solid_earth_tide =
        Some(gnss_rust::tides::SolidEarthTideModel::LegacyLoveApproximateSunMoon);
    cfg
}
fn tide_epochs(noisy: bool) -> Vec<gnss_rust::observation::ObservationEpoch> {
    RinexObservationReader::new(Cursor::new(if noisy {
        include_str!("fixtures/synthetic_ppp_legacy_tide_noisy.obs")
    } else {
        include_str!("fixtures/synthetic_ppp_legacy_tide.obs")
    }))
    .unwrap()
    .collect::<Result<_, _>>()
    .unwrap()
}
#[test]
fn native_legacy_tide_geometry_uses_instantaneous_receiver_before_precise_transmit() {
    let p = satellite_products();
    let cfg = tide_config();
    let mut count = 0;
    let mut transmit_changed = false;
    for line in include_str!("fixtures/upstream_ppp_legacy_tide_model.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = line.split(',').collect();
        let n: Vec<f64> = f[3..].iter().map(|v| v.parse().unwrap()).collect();
        let time = GnssTime::new(f[1].parse().unwrap(), f[2].parse().unwrap()).unwrap();
        let i = (time.tow() - 346200.0) / 30.0;
        let nominal = EcefCoord::new(n[0], n[1], n[2]);
        let sat = f[0].parse().unwrap();
        let m = evaluate_ppp_model(
            &p,
            sat,
            time,
            nominal,
            (2e-6 + i * 5e-10) * c::SPEED_OF_LIGHT,
            2.4,
            &cfg,
        )
        .unwrap();
        let tide = m.solid_earth_tide.as_ref().unwrap();
        assert_eq!(tide.nominal_marker_m, nominal);
        assert_eq!(m.receiver_position_m, tide.instantaneous_marker_m);
        for (v, expected) in [
            (m.receiver_position_m, &n[3..6]),
            (m.transmit.position_m, &n[6..9]),
            (m.satellite_position_m, &n[9..12]),
        ] {
            for (a, b) in [v.x, v.y, v.z].iter().zip(expected) {
                close(*a, *b, 1e-6);
            }
        }
        close(m.geometric_range_m, n[12], 1e-6);
        close(m.elevation_rad, n[13], 1e-11);
        close(m.troposphere_mapping, n[14], 1e-10);
        close(m.code_prediction_m, n[15], 1e-6);
        for (a, b) in m.head_design[..3].iter().zip(&n[16..19]) {
            close(*a, *b, 1e-10);
        }
        close(m.transmit.corrected_clock_bias_s, n[19], 1e-16);
        close(m.transmit.light_time_s, n[20], 1e-13);
        let no_tide =
            evaluate_ppp_model(&p, sat, time, nominal, 0.0, 2.4, &satellite_config()).unwrap();
        transmit_changed |= (m.transmit.light_time_s - no_tide.transmit.light_time_s).abs() > 1e-12;
        count += 1;
    }
    assert_eq!(count, 54);
    assert!(transmit_changed);
}
#[test]
fn complete_native_legacy_tide_injection_recovers_nominal_marker_and_continuous_arcs() {
    use gnss_rust::tides::ReceiverPositionReference;
    let p = satellite_products();
    let refs: BTreeMap<(usize, SatelliteId), Vec<f64>> =
        include_str!("fixtures/upstream_ppp_legacy_tide.csv")
            .lines()
            .skip(1)
            .map(|l| {
                let v: Vec<_> = l.split(',').collect();
                (
                    (v[0].parse().unwrap(), v[1].parse().unwrap()),
                    v[2..].iter().map(|x| x.parse().unwrap()).collect(),
                )
            })
            .collect();
    let truth: Vec<Vec<f64>> = include_str!("fixtures/synthetic_ppp_legacy_tide_truth.csv")
        .lines()
        .skip(1)
        .map(|l| l.split(',').map(|x| x.parse().unwrap()).collect())
        .collect();
    for noisy in [false, true] {
        let cfg = tide_config();
        let mut filter = PppFilter::new(cfg.clone()).unwrap();
        let mut max_error: f64 = 0.0;
        let mut rows = 0;
        for (i, epoch) in tide_epochs(noisy).iter().enumerate() {
            let result = solution(
                filter
                    .process_epoch(epoch, &p, &biases(), &selections())
                    .unwrap(),
            );
            let nominal = result.position.position_ecef.unwrap();
            assert_eq!(result.position.status, SolutionStatus::PppFloat);
            assert_eq!(
                result.receiver_reference,
                ReceiverPositionReference::NominalMarker
            );
            max_error = max_error.max(distance(nominal));
            let expected = refs
                .iter()
                .filter_map(|(&(j, s), v)| (i == j && v[2] == 1.0).then_some(s))
                .collect::<std::collections::BTreeSet<_>>();
            assert_eq!(
                result
                    .arc_commit
                    .accepted
                    .keys()
                    .copied()
                    .collect::<std::collections::BTreeSet<_>>(),
                expected
            );
            for m in &result.measurements {
                let tide = m.model.solid_earth_tide.as_ref().unwrap();
                let v = &truth[i];
                assert_eq!(tide.nominal_marker_m, nominal);
                assert_eq!(tide.instantaneous_marker_m, m.model.receiver_position_m);
                assert_eq!(tide.time, epoch.time);
                for (a, b) in [
                    tide.displacement_m.x,
                    tide.displacement_m.y,
                    tide.displacement_m.z,
                ]
                .iter()
                .zip(&v[6..9])
                {
                    close(*a, *b, if noisy { 1e-7 } else { 2e-10 });
                }
                let antenna = m.receiver_antenna.as_ref().unwrap();
                let windup = m.windup.as_ref().unwrap();
                close(
                    windup.angles.cycles,
                    refs[&(i, m.satellite)][0],
                    if noisy { 1e-7 } else { 1e-9 },
                );
                close(
                    antenna.if_add_m,
                    refs[&(i, m.satellite)][1],
                    if noisy { 1e-6 } else { 1e-8 },
                );
                let recomputed = evaluate_ppp_model(
                    &p,
                    m.satellite,
                    epoch.time,
                    nominal,
                    result.position.receiver_clock_bias * c::SPEED_OF_LIGHT,
                    result.zenith_delay_m,
                    &cfg,
                )
                .unwrap();
                assert_eq!(m.model, recomputed);
                assert_eq!(
                    windup.angles,
                    nominal_yaw_windup(
                        m.model.receiver_position_m,
                        m.model.satellite_position_m,
                        windup.sun_position_m,
                        windup.unwrap_reference_cycles
                    )
                    .unwrap()
                );
                let n = result
                    .ambiguities
                    .iter()
                    .find(|a| a.satellite == m.satellite)
                    .unwrap()
                    .metres;
                close(
                    m.phase_residual_m.unwrap(),
                    m.observation.phase.as_ref().unwrap().phase_if_m
                        + antenna.if_add_m
                        + windup.phase_if_add_m
                        - (m.model.code_prediction_m + n),
                    1e-9,
                );
                assert_eq!(
                    result.arc_commit.accepted[&m.satellite].accepted_epochs,
                    i as u64 + 1
                );
                assert_eq!(
                    result.arc_decisions[&m.satellite].reset_reasons.is_empty(),
                    i != 0
                );
                rows += 1;
            }
            if i == 63 {
                let final_error = distance(nominal);
                let max_ambiguity = result
                    .ambiguities
                    .iter()
                    .map(|a| (a.metres - ambiguities()[&a.satellite]).abs())
                    .fold(0.0, f64::max);
                eprintln!(
                    "legacy tide noisy={noisy} max_position={max_error:.9} final_position={final_error:.9} max_ambiguity={max_ambiguity:.9} phase_rms={:.9}",
                    result.phase_rms_m
                );
                // Rounded carriers (1e-4 cycles) amplify the early height/ZTD
                // transient. This 5 mm synthetic regression bound includes
                // input quantization; it is not a native solver parity bound.
                // The separate unrounded control retains a 0.5 mm bound.
                assert!(max_error < if noisy { 1.0 } else { 0.005 });
                assert!(final_error < if noisy { 0.01 } else { 0.001 });
                assert!(max_ambiguity < if noisy { 0.01 } else { 0.001 });
            }
        }
        assert_eq!(rows, 568);
    }
}
#[test]
fn tide_coordinate_reference_contradictions_and_invalid_station_are_rejected() {
    use gnss_rust::tides::*;
    let mut cfg = tide_config();
    cfg.receiver_reference = ReceiverPositionReference::InstantaneousMarker;
    assert_eq!(
        PppFilter::new(cfg).unwrap_err(),
        PppError::InvalidConfiguration
    );
    let mut cfg = tide_config();
    cfg.solid_earth_tide = None;
    assert_eq!(
        PppFilter::new(cfg).unwrap_err(),
        PppError::InvalidConfiguration
    );
    let mut cfg = tide_config();
    cfg.initial_position = EcefCoord::new(1e6, 0.0, 0.0);
    assert_eq!(
        PppFilter::new(cfg).unwrap_err(),
        PppError::Tide(TideError::InvalidStation)
    );
    let cfg = tide_config();
    let p = satellite_products();
    assert_eq!(
        evaluate_ppp_model(
            &p,
            "G01".parse().unwrap(),
            GnssTime::new(2300, 346200.0).unwrap(),
            EcefCoord::new(0.0, 0.0, 0.0),
            0.0,
            2.4,
            &cfg
        ),
        Err(PppError::Tide(TideError::InvalidStation))
    );
}
#[test]
fn rejecting_tide_epoch_preserves_numeric_state_and_recovery_starts_fresh_arcs() {
    let epochs = tide_epochs(false);
    let p = satellite_products();
    let mut f = PppFilter::new(tide_config()).unwrap();
    let sat = "G01".parse().unwrap();
    for e in &epochs[..8] {
        f.process_epoch(e, &p, &biases(), &selections()).unwrap();
    }
    let state = f.state().clone();
    let time = f.state_time();
    let arc = f.phase_arc(sat).unwrap().clone();
    let windup = f.windup_state(sat).unwrap().clone();
    assert_eq!(
        f.process_epoch(&epochs[7], &p, &biases(), &selections()),
        Err(PppError::Arc(PhaseArcError::OutOfOrder))
    );
    assert_eq!(f.state(), &state);
    assert_eq!(f.phase_arc(sat), Some(&arc));
    assert_eq!(f.windup_state(sat), Some(&windup));
    let mut bad = epochs[8].clone();
    for observation in &mut bad.observations {
        observation.pseudorange_m = observation.pseudorange_m.map(|p| p + 5000.0);
    }
    assert_eq!(
        f.process_epoch(&bad, &p, &biases(), &selections()),
        Err(PppError::Filter(FilterError::InnovationRejected))
    );
    assert_eq!(f.state(), &state);
    assert_eq!(f.state_time(), time);
    assert!(f.phase_arc(sat).is_none());
    assert!(f.windup_state(sat).is_none());
    let r = solution(
        f.process_epoch(&epochs[9], &p, &biases(), &selections())
            .unwrap(),
    );
    for (sat, arc) in &r.arc_commit.accepted {
        assert_eq!(arc.accepted_epochs, 1);
        assert!(!r.arc_decisions[sat].reset_reasons.is_empty());
    }
    assert!(
        r.measurements
            .iter()
            .all(|m| m.windup.as_ref().unwrap().previous_cycles == 0.0)
    );
}
#[test]
fn tide_disabled_filter_reports_the_instantaneous_marker_at_the_first_epoch() {
    use gnss_rust::tides::ReceiverPositionReference;
    let epoch = tide_epochs(false).remove(0);
    let p = satellite_products();
    let mut filter = PppFilter::new(satellite_config()).unwrap();
    let result = solution(
        filter
            .process_epoch(&epoch, &p, &biases(), &selections())
            .unwrap(),
    );
    assert_eq!(
        result.receiver_reference,
        ReceiverPositionReference::InstantaneousMarker
    );
    let v: Vec<f64> = include_str!("fixtures/synthetic_ppp_legacy_tide_truth.csv")
        .lines()
        .nth(1)
        .unwrap()
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect();
    let position = result.position.position_ecef.unwrap();
    let error = (position.x - v[9])
        .hypot(position.y - v[10])
        .hypot(position.z - v[11]);
    assert!(error < 0.003);
    assert!(distance(position) > 0.08);
    for m in &result.measurements {
        assert!(m.model.solid_earth_tide.is_none());
        assert_eq!(m.model.receiver_position_m, position);
    }
}

#[test]
fn native_unrounded_carriers_isolate_rinex_quantization_from_tide_model_error() {
    let exact: BTreeMap<(usize, SatelliteId), Vec<f64>> =
        include_str!("fixtures/upstream_ppp_legacy_tide_unrounded.csv")
            .lines()
            .skip(1)
            .map(|l| {
                let v: Vec<_> = l.split(',').collect();
                (
                    (v[0].parse().unwrap(), v[1].parse().unwrap()),
                    v[2..].iter().map(|x| x.parse().unwrap()).collect(),
                )
            })
            .collect();
    let p = satellite_products();
    let mut f = PppFilter::new(tide_config()).unwrap();
    let mut max_error: f64 = 0.0;
    let mut max_phase_rounding: f64 = 0.0;
    for (i, mut epoch) in tide_epochs(false).into_iter().enumerate() {
        for observation in &mut epoch.observations {
            let row = &exact[&(i, observation.satellite)];
            let j = usize::from(observation.tracking_code == "2W");
            close(observation.pseudorange_m.unwrap(), row[2 * j], 5.01e-6);
            let delta = (observation.carrier_phase_cycles.unwrap() - row[2 * j + 1]).abs();
            assert!(delta < 5.01e-5);
            max_phase_rounding = max_phase_rounding.max(delta);
            observation.pseudorange_m = Some(row[2 * j]);
            observation.carrier_phase_cycles = Some(row[2 * j + 1]);
        }
        let r = solution(
            f.process_epoch(&epoch, &p, &biases(), &selections())
                .unwrap(),
        );
        max_error = max_error.max(distance(r.position.position_ecef.unwrap()));
        if i == 63 {
            assert!(distance(r.position.position_ecef.unwrap()) < 1e-6);
            eprintln!(
                "unrounded legacy tide max_position={max_error:.9} final_position={:.9} max_phase_rounding_cycles={max_phase_rounding:.9}",
                distance(r.position.position_ecef.unwrap())
            );
        }
    }
    assert!(max_error < 0.0005);
    assert!(max_phase_rounding > 4.9e-5);
}

fn satellite_products() -> PreciseProducts {
    PreciseProducts::new(
        read_sp3(Cursor::new(include_str!("fixtures/synthetic_ppp_com.sp3"))).unwrap(),
        Some(read_clk(Cursor::new(include_str!("fixtures/synthetic_ppp.clk"))).unwrap()),
    )
}
fn satellite_model(pair: [String; 2]) -> gnss_rust::antenna::SatelliteAntennaModel {
    use gnss_rust::antenna::*;
    let atx = read_antex(Cursor::new(include_str!(
        "fixtures/synthetic_satellite.atx"
    )))
    .unwrap();
    let entries = atx
        .entries()
        .iter()
        .map(|entry| {
            atx.satellite_calibration(
                entry.satellite().unwrap(),
                &entry.svn,
                GnssTime::new(2300, 346200.0).unwrap(),
            )
            .unwrap()
            .clone()
        })
        .collect();
    SatelliteAntennaModel::new(
        SatelliteAntennaBinding {
            product_id: "synthetic_satellite.atx".into(),
            clock_source: ClockSource::RinexClk,
            clock_reference: "synthetic-precise-clock-v1".into(),
            tracking_codes: pair,
        },
        entries,
    )
    .unwrap()
}
fn satellite_config() -> PppConfig {
    let mut cfg = config();
    cfg.expected_reference_point = OrbitReferencePoint::CentreOfMass;
    cfg.satellite_antenna = Some(satellite_model(["1C".into(), "2W".into()]));
    cfg.receiver_antenna = Some(receiver_antenna_model(None));
    cfg.phase_windup = Some(PppWindupModel::NominalYawApproximateSun);
    cfg
}
fn satellite_epochs(noisy: bool) -> Vec<gnss_rust::observation::ObservationEpoch> {
    RinexObservationReader::new(Cursor::new(if noisy {
        include_str!("fixtures/synthetic_ppp_satellite_antenna_noisy.obs")
    } else {
        include_str!("fixtures/synthetic_ppp_satellite_antenna.obs")
    }))
    .unwrap()
    .collect::<Result<_, _>>()
    .unwrap()
}
#[test]
fn native_com_to_apc_geometry_uses_shifted_position_and_unchanged_precise_clock() {
    let p = satellite_products();
    let truth = truth();
    let receiver = EcefCoord::new(truth[0], truth[1], truth[2]);
    let mut count = 0;
    for line in include_str!("fixtures/upstream_satellite_antenna.csv")
        .lines()
        .skip(1)
    {
        let f: Vec<_> = line.split(',').collect();
        let n: Vec<f64> = f[5..].iter().map(|v| v.parse().unwrap()).collect();
        let mut cfg = satellite_config();
        cfg.satellite_antenna = Some(satellite_model([f[3].into(), f[4].into()]));
        let tow = f[2].parse::<f64>().unwrap();
        let i = (tow - 346200.0) / 30.0;
        let m = evaluate_ppp_model(
            &p,
            f[0].parse().unwrap(),
            GnssTime::new(2300, tow).unwrap(),
            receiver,
            (2e-6 + i * 5e-10) * c::SPEED_OF_LIGHT,
            2.4,
            &cfg,
        )
        .unwrap();
        assert_eq!(
            m.transmit.evaluated_state.reference_point,
            OrbitReferencePoint::CentreOfMass
        );
        for (a, b) in [
            m.transmit.position_m.x,
            m.transmit.position_m.y,
            m.transmit.position_m.z,
        ]
        .iter()
        .zip(&n[..3])
        {
            close(*a, *b, 1e-6);
        }
        for (a, b) in [
            m.satellite_position_m.x,
            m.satellite_position_m.y,
            m.satellite_position_m.z,
        ]
        .iter()
        .zip(&n[12..15])
        {
            close(*a, *b, 1e-6);
        }
        assert_eq!(
            m.satellite_position_m,
            m.satellite_antenna.as_ref().unwrap().position_apc_m
        );
        close(m.geometric_range_m, n[15], 1e-6);
        close(m.elevation_rad, n[16], 1e-11);
        close(m.troposphere_mapping, n[17], 1e-10);
        close(m.code_prediction_m, n[18], 1e-6);
        for (a, b) in m.head_design[..3].iter().zip(&n[19..22]) {
            close(*a, *b, 1e-10);
        }
        close(m.transmit.corrected_clock_bias_s, n[22], 1e-16);
        count += 1;
    }
    assert_eq!(count, 81);
}
#[test]
fn native_satellite_receiver_and_windup_injection_recovers_marker_and_continuous_arcs() {
    let reference: BTreeMap<(usize, SatelliteId), Vec<f64>> =
        include_str!("fixtures/upstream_ppp_satellite_antenna.csv")
            .lines()
            .skip(1)
            .map(|l| {
                let v: Vec<_> = l.split(',').collect();
                (
                    (v[0].parse().unwrap(), v[1].parse().unwrap()),
                    v[2..].iter().map(|s| s.parse().unwrap()).collect(),
                )
            })
            .collect();
    let p = satellite_products();
    for noisy in [false, true] {
        let cfg = satellite_config();
        let mut f = PppFilter::new(cfg.clone()).unwrap();
        let mut max_position: f64 = 0.0;
        let mut count = 0;
        for (i, epoch) in satellite_epochs(noisy).iter().enumerate() {
            let r = solution(
                f.process_epoch(epoch, &p, &biases(), &selections())
                    .unwrap(),
            );
            let pos = r.position.position_ecef.unwrap();
            assert_eq!(r.position.status, SolutionStatus::PppFloat);
            assert_eq!(r.reference_point, OrbitReferencePoint::CentreOfMass);
            let expected = reference
                .iter()
                .filter_map(|(&(j, sat), v)| (i == j && v[9] == 1.0).then_some(sat))
                .collect::<std::collections::BTreeSet<_>>();
            assert_eq!(
                r.arc_commit
                    .accepted
                    .keys()
                    .copied()
                    .collect::<std::collections::BTreeSet<_>>(),
                expected
            );
            max_position = max_position.max(distance(pos));
            for m in &r.measurements {
                let v = &reference[&(i, m.satellite)];
                let a = m.model.satellite_antenna.as_ref().unwrap();
                assert_eq!(a.product_id, "synthetic_satellite.atx");
                assert_eq!(a.satellite, m.satellite);
                for (x, y) in [
                    m.model.satellite_position_m.x,
                    m.model.satellite_position_m.y,
                    m.model.satellite_position_m.z,
                ]
                .iter()
                .zip(&v[..3])
                {
                    close(*x, *y, if noisy { 1e-4 } else { 1e-6 });
                }
                close(
                    m.windup.as_ref().unwrap().angles.cycles,
                    v[7],
                    if noisy { 1e-7 } else { 1e-9 },
                );
                close(
                    m.receiver_antenna.as_ref().unwrap().if_add_m,
                    v[8],
                    if noisy { 1e-6 } else { 1e-8 },
                );
                let recomputed = evaluate_ppp_model(
                    &p,
                    m.satellite,
                    epoch.time,
                    pos,
                    r.position.receiver_clock_bias * c::SPEED_OF_LIGHT,
                    r.zenith_delay_m,
                    &cfg,
                )
                .unwrap();
                assert_eq!(m.model, recomputed);
                let windup = m.windup.as_ref().unwrap();
                assert_eq!(
                    windup.angles,
                    nominal_yaw_windup(
                        pos,
                        m.model.satellite_position_m,
                        windup.sun_position_m,
                        windup.unwrap_reference_cycles
                    )
                    .unwrap()
                );
                let n = r
                    .ambiguities
                    .iter()
                    .find(|a| a.satellite == m.satellite)
                    .unwrap()
                    .metres;
                close(
                    m.phase_residual_m.unwrap(),
                    m.observation.phase.as_ref().unwrap().phase_if_m
                        + m.receiver_antenna.as_ref().unwrap().if_add_m
                        + windup.phase_if_add_m
                        - (m.model.code_prediction_m + n),
                    1e-9,
                );
                assert_eq!(
                    r.arc_commit.accepted[&m.satellite].accepted_epochs,
                    i as u64 + 1
                );
                assert_eq!(
                    r.arc_decisions[&m.satellite].reset_reasons.is_empty(),
                    i != 0
                );
                count += 1;
            }
            if i == 63 {
                let final_position = distance(pos);
                let max_ambiguity = r
                    .ambiguities
                    .iter()
                    .map(|a| (a.metres - ambiguities()[&a.satellite]).abs())
                    .fold(0.0, f64::max);
                eprintln!(
                    "satellite antenna noisy={noisy} max_position={max_position:.9} final_position={final_position:.9} max_ambiguity={max_ambiguity:.9} phase_rms={:.9}",
                    r.phase_rms_m
                );
                assert!(max_position < if noisy { 1.0 } else { 0.003 });
                assert!(final_position < if noisy { 0.01 } else { 0.001 });
                assert!(max_ambiguity < if noisy { 0.01 } else { 0.001 });
            }
        }
        assert_eq!(count, 568);
    }
}
#[test]
fn ppp_requires_satellite_pco_for_com_and_forbids_it_for_apc() {
    let mut cfg = config();
    cfg.satellite_antenna = Some(satellite_model(["1C".into(), "2W".into()]));
    assert_eq!(
        PppFilter::new(cfg).unwrap_err(),
        PppError::InvalidConfiguration
    );
    let mut cfg = satellite_config();
    cfg.satellite_antenna = None;
    assert_eq!(
        PppFilter::new(cfg).unwrap_err(),
        PppError::InvalidConfiguration
    );
    let mut f = PppFilter::new(satellite_config()).unwrap();
    let e = satellite_epochs(false);
    f.process_epoch(&e[0], &satellite_products(), &biases(), &selections())
        .unwrap();
    let state = f.state().clone();
    let time = f.state_time();
    assert_eq!(
        f.process_epoch(&e[1], &products(), &biases(), &selections()),
        Err(PppError::ProductFrameMismatch)
    );
    assert_eq!(f.state(), &state);
    assert_eq!(f.state_time(), time);
    assert!(f.phase_arc("G01".parse().unwrap()).is_none());
}
#[test]
fn satellite_clock_and_signal_binding_failure_interrupts_phase_without_numeric_update() {
    let e = satellite_epochs(false);
    let p = satellite_products();
    let sat = "G01".parse().unwrap();
    for kind in 0..3 {
        let mut f = PppFilter::new(satellite_config()).unwrap();
        f.process_epoch(&e[0], &p, &biases(), &selections())
            .unwrap();
        let state = f.state().clone();
        let time = f.state_time();
        let mut selected = selections();
        let s = selected.get_mut(&sat).unwrap();
        match kind {
            0 => s.clock_reference = "wrong-clock".into(),
            1 => s.clock_source = ClockSource::Sp3,
            _ => s.tracking_codes[1] = "2C".into(),
        }
        assert_eq!(
            f.process_epoch(&e[1], &p, &biases(), &selected),
            Err(PppError::Antenna(
                gnss_rust::antenna::AntennaError::DatumMismatch
            ))
        );
        assert_eq!(f.state(), &state);
        assert_eq!(f.state_time(), time);
        assert!(f.windup_state(sat).is_none());
        assert!(f.phase_arc(sat).is_none());
    }
    let mut cfg = satellite_config();
    let a = cfg.satellite_antenna.take().unwrap();
    let mut binding = a.binding().clone();
    binding.clock_source = ClockSource::Sp3;
    cfg.satellite_antenna = Some(
        gnss_rust::antenna::SatelliteAntennaModel::new(
            binding,
            a.calibrations().values().cloned().collect(),
        )
        .unwrap(),
    );
    let mut selected = selections();
    for s in selected.values_mut() {
        s.clock_source = ClockSource::Sp3;
    }
    // Public geometry also checks the actual CLK-selected source before PCO.
    assert_eq!(
        evaluate_ppp_model(&p, sat, e[0].time, cfg.initial_position, 0.0, 2.3, &cfg),
        Err(PppError::Antenna(
            gnss_rust::antenna::AntennaError::DatumMismatch
        ))
    );
}
#[test]
fn expired_or_missing_satellite_calibration_is_a_hard_failure_with_no_zero_fallback() {
    use gnss_rust::antenna::*;
    let epochs = satellite_epochs(false);
    let p = satellite_products();
    let sat = "G01".parse().unwrap();
    let a = satellite_model(["1C".into(), "2W".into()]);
    let mut entries = a.calibrations().values().cloned().collect::<Vec<_>>();
    entries
        .iter_mut()
        .find(|c| c.satellite() == Some(sat))
        .unwrap()
        .valid_until = Some(epochs[1].time);
    let mut cfg = satellite_config();
    cfg.satellite_antenna = Some(SatelliteAntennaModel::new(a.binding().clone(), entries).unwrap());
    let mut f = PppFilter::new(cfg).unwrap();
    for e in &epochs[..2] {
        f.process_epoch(e, &p, &biases(), &selections()).unwrap();
    }
    let state = f.state().clone();
    let time = f.state_time();
    let windup = f.windup_state(sat).unwrap().clone();
    assert!(matches!(
        f.process_epoch(&epochs[1], &p, &biases(), &selections()),
        Err(PppError::Arc(PhaseArcError::OutOfOrder))
    ));
    assert_eq!(f.windup_state(sat), Some(&windup));
    assert_eq!(
        f.process_epoch(&epochs[2], &p, &biases(), &selections()),
        Err(PppError::Antenna(AntennaError::OutsideValidity))
    );
    assert_eq!(f.state(), &state);
    assert_eq!(f.state_time(), time);
    assert!(f.phase_arc(sat).is_none());
    assert!(f.windup_state(sat).is_none());
    let mut cfg = satellite_config();
    cfg.satellite_antenna = Some(
        SatelliteAntennaModel::new(
            a.binding().clone(),
            a.calibrations()
                .iter()
                .filter_map(|(s, c)| (*s != sat).then_some(c.clone()))
                .collect(),
        )
        .unwrap(),
    );
    let mut f = PppFilter::new(cfg).unwrap();
    let state = f.state().clone();
    assert_eq!(
        f.process_epoch(&epochs[0], &p, &biases(), &selections()),
        Err(PppError::Antenna(AntennaError::MissingCalibration))
    );
    assert_eq!(f.state(), &state);
    assert_eq!(f.state_time(), None);
}

#[test]
fn degenerate_satellite_attitude_rejects_epoch_even_when_windup_is_disabled() {
    let epochs = satellite_epochs(false);
    let p = satellite_products();
    let mut cfg = satellite_config();
    cfg.phase_windup = None;
    let mut f = PppFilter::new(cfg).unwrap();
    f.process_epoch(&epochs[0], &p, &biases(), &selections())
        .unwrap();
    let state = f.state().clone();
    let time = f.state_time();
    let sun = approximate_sun_position_ecef(epochs[1].time).unwrap();
    let scale = 2.6e7 / sun.x.hypot(sun.y).hypot(sun.z);
    let orbit = include_str!("fixtures/synthetic_ppp_com.sp3")
        .lines()
        .map(|line| {
            let mut line = line.to_owned();
            if line.starts_with("PG01") {
                line.replace_range(
                    4..46,
                    &format!(
                        "{:14.6}{:14.6}{:14.6}",
                        sun.x * scale / 1000.0,
                        sun.y * scale / 1000.0,
                        sun.z * scale / 1000.0
                    ),
                );
            }
            line
        })
        .collect::<Vec<_>>()
        .join("\n");
    let degenerate = PreciseProducts::new(
        read_sp3(Cursor::new(orbit)).unwrap(),
        Some(read_clk(Cursor::new(include_str!("fixtures/synthetic_ppp.clk"))).unwrap()),
    );
    assert_eq!(
        f.process_epoch(&epochs[1], &degenerate, &biases(), &selections()),
        Err(PppError::Antenna(
            gnss_rust::antenna::AntennaError::Attitude(
                gnss_rust::windup::WindupError::DegenerateAttitude
            )
        ))
    );
    assert_eq!(f.state(), &state);
    assert_eq!(f.state_time(), time);
    assert!(f.phase_arc("G01".parse().unwrap()).is_none());
}

fn receiver_antenna_model(
    valid_until: Option<GnssTime>,
) -> gnss_rust::antenna::ReceiverAntennaModel {
    use gnss_rust::antenna::*;
    let atx = read_antex(Cursor::new(include_str!("fixtures/synthetic_receiver.atx"))).unwrap();
    let mut calibration = atx
        .receiver(
            "TESTANT NONE",
            "RX001",
            GnssTime::new(2300, 346200.0).unwrap(),
        )
        .unwrap()
        .clone();
    if let Some(t) = valid_until {
        calibration.valid_until = Some(t);
    }
    ReceiverAntennaModel::new(
        "synthetic_receiver.atx".into(),
        &calibration,
        gnss_rust::EnuCoord {
            east: 0.025,
            north: -0.035,
            up: 0.8,
        },
    )
    .unwrap()
}
fn antenna_epochs(noisy: bool) -> Vec<gnss_rust::observation::ObservationEpoch> {
    RinexObservationReader::new(Cursor::new(if noisy {
        include_str!("fixtures/synthetic_ppp_antenna_noisy.obs")
    } else {
        include_str!("fixtures/synthetic_ppp_antenna.obs")
    }))
    .unwrap()
    .collect::<Result<_, _>>()
    .unwrap()
}
#[test]
fn native_receiver_injection_recovers_marker_without_resetting_smooth_phase_arcs() {
    let references: BTreeMap<(usize, SatelliteId), f64> =
        include_str!("fixtures/upstream_ppp_receiver_antenna.csv")
            .lines()
            .skip(1)
            .map(|l| {
                let v: Vec<_> = l.split(',').collect();
                (
                    (v[0].parse().unwrap(), v[1].parse().unwrap()),
                    v[6].parse().unwrap(),
                )
            })
            .collect();
    let products = products();
    let expected = admitted();
    for noisy in [false, true] {
        let mut cfg = config();
        cfg.phase_windup = Some(PppWindupModel::NominalYawApproximateSun);
        cfg.receiver_antenna = Some(receiver_antenna_model(None));
        let mut filter = PppFilter::new(cfg).unwrap();
        let mut maximum: f64 = 0.0;
        let mut code_rows = 0;
        let mut phase_rows = 0;
        for (i, epoch) in antenna_epochs(noisy).iter().enumerate() {
            let result = solution(
                filter
                    .process_epoch(epoch, &products, &biases(), &selections())
                    .unwrap(),
            );
            assert_eq!(result.position.status, SolutionStatus::PppFloat);
            assert_eq!(
                result
                    .arc_commit
                    .accepted
                    .keys()
                    .copied()
                    .collect::<std::collections::BTreeSet<_>>(),
                expected[i]
            );
            maximum = maximum.max(distance(result.position.position_ecef.unwrap()));
            for m in &result.measurements {
                let antenna = m.receiver_antenna.as_ref().unwrap();
                assert_eq!(antenna.product_id, "synthetic_receiver.atx");
                assert_eq!(antenna.frequency_codes, ["G01", "G02"]);
                close(
                    antenna.if_add_m,
                    references[&(i, m.satellite)],
                    if noisy { 1e-6 } else { 1e-8 },
                );
                close(
                    m.code_residual_m,
                    m.observation.code_if_m + antenna.if_add_m - m.model.code_prediction_m,
                    1e-9,
                );
                let n = result
                    .ambiguities
                    .iter()
                    .find(|a| a.satellite == m.satellite)
                    .unwrap()
                    .metres;
                close(
                    m.phase_residual_m.unwrap(),
                    m.observation.phase.as_ref().unwrap().phase_if_m
                        + antenna.if_add_m
                        + m.windup.as_ref().unwrap().phase_if_add_m
                        - (m.model.code_prediction_m + n),
                    1e-9,
                );
                assert_eq!(
                    result.arc_commit.accepted[&m.satellite].accepted_epochs,
                    i as u64 + 1
                );
                assert_eq!(
                    result.arc_decisions[&m.satellite].reset_reasons.is_empty(),
                    i != 0
                );
                code_rows += 1;
                phase_rows += usize::from(m.phase_residual_m.is_some());
            }
            if i == 63 {
                let final_error = distance(result.position.position_ecef.unwrap());
                let max_ambiguity = result
                    .ambiguities
                    .iter()
                    .map(|a| (a.metres - ambiguities()[&a.satellite]).abs())
                    .fold(0.0, f64::max);
                eprintln!(
                    "receiver antenna noisy={noisy} max_position={maximum:.9} final_position={final_error:.9} max_ambiguity={max_ambiguity:.9} phase_rms={:.9}",
                    result.phase_rms_m
                );
                assert!(maximum < if noisy { 1.0 } else { 0.003 });
                assert!(final_error < if noisy { 0.01 } else { 0.001 });
                assert!(max_ambiguity < if noisy { 0.01 } else { 0.001 });
            }
        }
        assert_eq!((code_rows, phase_rows), (568, 568));
    }
}
#[test]
fn antenna_expiry_rejects_update_preserves_numeric_state_and_interrupts_phase() {
    let epochs = antenna_epochs(false);
    let p = products();
    let mut cfg = config();
    cfg.phase_windup = Some(PppWindupModel::NominalYawApproximateSun);
    cfg.receiver_antenna = Some(receiver_antenna_model(Some(epochs[1].time)));
    let mut f = PppFilter::new(cfg).unwrap();
    for epoch in &epochs[..2] {
        f.process_epoch(epoch, &p, &biases(), &selections())
            .unwrap();
    }
    let state = f.state().clone();
    let time = f.state_time();
    let sat = "G01".parse().unwrap();
    let windup = f.windup_state(sat).unwrap().clone();
    assert!(matches!(
        f.process_epoch(&epochs[1], &p, &biases(), &selections()),
        Err(PppError::Arc(PhaseArcError::OutOfOrder))
    ));
    assert_eq!(f.windup_state(sat), Some(&windup));
    assert_eq!(
        f.process_epoch(&epochs[2], &p, &biases(), &selections()),
        Err(PppError::Antenna(
            gnss_rust::antenna::AntennaError::OutsideValidity
        ))
    );
    assert_eq!(f.state(), &state);
    assert_eq!(f.state_time(), time);
    assert!(f.windup_state(sat).is_none());
    assert!(f.phase_arc(sat).is_none());
}
#[test]
fn disabling_antenna_on_injected_data_measures_a_different_receiver_reference_point() {
    let mut f = windup_filter();
    let p = products();
    let mut final_position = None;
    for epoch in antenna_epochs(false) {
        let result = solution(
            f.process_epoch(&epoch, &p, &biases(), &selections())
                .unwrap(),
        );
        assert!(
            result
                .measurements
                .iter()
                .all(|m| m.receiver_antenna.is_none())
        );
        final_position = result.position.position_ecef;
    }
    assert!(distance(final_position.unwrap()) > 0.5);
}

#[test]
fn missing_receiver_frequency_is_not_silently_replaced_by_zero_in_ppp() {
    let m = receiver_antenna_model(None);
    let mut c = m.calibration().clone();
    c.frequencies.remove("G02");
    let mut cfg = config();
    cfg.receiver_antenna = Some(
        gnss_rust::antenna::ReceiverAntennaModel::new(
            "missing-frequency".into(),
            &c,
            gnss_rust::EnuCoord {
                east: 0.025,
                north: -0.035,
                up: 0.8,
            },
        )
        .unwrap(),
    );
    let mut f = PppFilter::new(cfg).unwrap();
    let state = f.state().clone();
    assert_eq!(
        f.process_epoch(
            &antenna_epochs(false)[0],
            &products(),
            &biases(),
            &selections()
        ),
        Err(PppError::Antenna(
            gnss_rust::antenna::AntennaError::MissingFrequency("G02".into())
        ))
    );
    assert_eq!(f.state(), &state);
    assert_eq!(f.state_time(), None);
    assert!(f.phase_arc("G01".parse().unwrap()).is_none());
}

fn products() -> PreciseProducts {
    PreciseProducts::new(
        read_sp3(Cursor::new(include_str!("fixtures/synthetic_ppp.sp3"))).unwrap(),
        Some(read_clk(Cursor::new(include_str!("fixtures/synthetic_ppp.clk"))).unwrap()),
    )
}

fn osb_binding() -> OsbBinding {
    OsbBinding {
        product_id: "synthetic_ppp_osb.bsx".into(),
        clock_source: ClockSource::RinexClk,
        clock_reference: "synthetic-precise-clock-v1".into(),
        phase_reference: Some("synthetic-phase-v1".into()),
    }
}
fn osb_requests() -> Vec<BiasSignalRequest> {
    selections()
        .into_iter()
        .flat_map(|(satellite, selection)| {
            selection
                .tracking_codes
                .into_iter()
                .map(move |tracking_code| BiasSignalRequest {
                    satellite,
                    tracking_code,
                    glonass_channel: None,
                })
        })
        .collect()
}

fn windup_epochs(noisy: bool) -> Vec<gnss_rust::observation::ObservationEpoch> {
    RinexObservationReader::new(Cursor::new(if noisy {
        include_str!("fixtures/synthetic_ppp_windup_noisy.obs")
    } else {
        include_str!("fixtures/synthetic_ppp_windup.obs")
    }))
    .unwrap()
    .collect::<Result<_, _>>()
    .unwrap()
}
fn windup_filter() -> PppFilter {
    let mut config = config();
    config.phase_windup = Some(PppWindupModel::NominalYawApproximateSun);
    PppFilter::new(config).unwrap()
}
fn windup_reference() -> BTreeMap<(usize, SatelliteId), (f64, f64)> {
    include_str!("fixtures/upstream_ppp_windup.csv")
        .lines()
        .skip(1)
        .map(|line| {
            let f: Vec<_> = line.split(',').collect();
            (
                (f[0].parse().unwrap(), f[1].parse().unwrap()),
                (f[8].parse().unwrap(), f[9].parse().unwrap()),
            )
        })
        .collect()
}
#[test]
fn windup_final_iteration_recomputes_geometry_and_degenerate_attitude_rejects_epoch() {
    let p = products();
    let mut f = windup_filter();
    let epochs = windup_epochs(false);
    for epoch in &epochs[..8] {
        let r = solution(
            f.process_epoch(epoch, &p, &biases(), &selections())
                .unwrap(),
        );
        let receiver = r.position.position_ecef.unwrap();
        for m in &r.measurements {
            let w = m.windup.as_ref().unwrap();
            let expected = nominal_yaw_windup(
                receiver,
                m.model.transmit.position_m,
                w.sun_position_m,
                w.unwrap_reference_cycles,
            )
            .unwrap();
            assert_eq!(w.angles, expected);
            if epoch.time == epochs[0].time {
                let seed = config();
                let initial = evaluate_ppp_model(
                    &p,
                    m.satellite,
                    epoch.time,
                    seed.initial_position,
                    seed.initial_receiver_clock_s * c::SPEED_OF_LIGHT,
                    seed.initial_zenith_delay_m,
                    &seed,
                )
                .unwrap();
                let anchor = nominal_yaw_windup(
                    seed.initial_position,
                    initial.transmit.position_m,
                    w.sun_position_m,
                    0.0,
                )
                .unwrap();
                assert_eq!(w.previous_cycles, 0.0);
                assert_eq!(w.unwrap_reference_cycles, anchor.cycles);
            } else {
                assert_eq!(w.unwrap_reference_cycles, w.previous_cycles);
            }
            for (i, &add) in w.phase_add_m.iter().enumerate() {
                close(add, -w.angles.cycles * m.observation.wavelength_m[i], 1e-15);
            }
            assert_eq!(
                m.observation.code_if_m,
                m.observation.coefficients[0] * m.observation.corrected_code_m[0]
                    + m.observation.coefficients[1] * m.observation.corrected_code_m[1]
            );
        }
    }
    let state = f.state().clone();
    let time = f.state_time();
    let sun = approximate_sun_position_ecef(epochs[8].time).unwrap();
    let scale = 2.6e7 / (sun.x * sun.x + sun.y * sun.y + sun.z * sun.z).sqrt();
    let mut orbit = String::new();
    for line in include_str!("fixtures/synthetic_ppp.sp3").lines() {
        let mut line = line.to_owned();
        if line.starts_with("PG01") {
            let xyz = format!(
                "{:14.6}{:14.6}{:14.6}",
                sun.x * scale / 1000.0,
                sun.y * scale / 1000.0,
                sun.z * scale / 1000.0
            );
            line.replace_range(4..46, &xyz);
        }
        orbit.push_str(&line);
        orbit.push('\n');
    }
    let degenerate = PreciseProducts::new(
        read_sp3(Cursor::new(orbit)).unwrap(),
        Some(read_clk(Cursor::new(include_str!("fixtures/synthetic_ppp.clk"))).unwrap()),
    );
    assert!(matches!(
        f.process_epoch(&epochs[8], &degenerate, &biases(), &selections()),
        Err(PppError::Windup(
            gnss_rust::windup::WindupError::DegenerateAttitude
        ))
    ));
    assert_eq!(f.state(), &state);
    assert_eq!(f.state_time(), time);
    assert!(f.windup_state("G01".parse().unwrap()).is_none());
}

#[test]
fn native_windup_injected_phase_traces_recover_truth_and_keep_smooth_arcs() {
    let p = products();
    let refs = windup_reference();
    assert_eq!(refs.len(), 576);
    for noisy in [false, true] {
        let mut f = windup_filter();
        let mut count = 0;
        let mut maximum: f64 = 0.0;
        let expected = admitted();
        let mut previous = BTreeMap::new();
        for (i, epoch) in windup_epochs(noisy).iter().enumerate() {
            let result = solution(
                f.process_epoch(epoch, &p, &biases(), &selections())
                    .unwrap(),
            );
            assert_eq!(result.position.status, SolutionStatus::PppFloat);
            assert_eq!(
                result
                    .arc_commit
                    .accepted
                    .keys()
                    .copied()
                    .collect::<std::collections::BTreeSet<_>>(),
                expected[i]
            );
            maximum = maximum.max(distance(result.position.position_ecef.unwrap()));
            for m in &result.measurements {
                let w = m.windup.as_ref().unwrap();
                let cache = f.windup_state(m.satellite).unwrap();
                assert_eq!(cache.cycles, w.angles.cycles);
                assert_eq!(cache.time, epoch.time);
                assert_eq!(cache.arc_id, result.arc_commit.accepted[&m.satellite].id);
                assert_eq!(
                    w.previous_cycles,
                    *previous.get(&m.satellite).unwrap_or(&0.0)
                );
                previous.insert(m.satellite, w.angles.cycles);
                assert_eq!(
                    result.arc_commit.accepted[&m.satellite].accepted_epochs,
                    i as u64 + 1
                );
                close(
                    w.angles.cycles,
                    refs[&(i, m.satellite)].0,
                    if noisy { 1e-7 } else { 1e-9 },
                );
                close(
                    w.phase_if_add_m,
                    refs[&(i, m.satellite)].1,
                    if noisy { 1e-8 } else { 1e-10 },
                );
                close(
                    w.phase_if_add_m,
                    m.observation.coefficients[0] * w.phase_add_m[0]
                        + m.observation.coefficients[1] * w.phase_add_m[1],
                    1e-15,
                );
                close(
                    m.phase_residual_m.unwrap(),
                    m.observation.phase.as_ref().unwrap().phase_if_m + w.phase_if_add_m
                        - (m.model.code_prediction_m
                            + result
                                .ambiguities
                                .iter()
                                .find(|a| a.satellite == m.satellite)
                                .unwrap()
                                .metres),
                    1e-9,
                );
            }
            count += 1;
            if i == 63 {
                let final_error = distance(result.position.position_ecef.unwrap());
                assert!(maximum < if noisy { 1.0 } else { 0.002 });
                assert!(final_error < if noisy { 0.01 } else { 0.001 });
                let max_ambiguity = result
                    .ambiguities
                    .iter()
                    .map(|a| (a.metres - ambiguities()[&a.satellite]).abs())
                    .fold(0.0, f64::max);
                assert!(max_ambiguity < if noisy { 0.01 } else { 0.001 });
                eprintln!(
                    "windup noisy={noisy} max_position={maximum:.9} final_position={final_error:.9} max_ambiguity={max_ambiguity:.9} phase_rms={:.9}",
                    result.phase_rms_m
                );
                assert!(f.windup_state("G09".parse().unwrap()).is_none());
            }
        }
        assert_eq!(count, 64);
    }
}
#[test]
fn rejected_and_duplicate_epochs_do_not_commit_windup_and_reacquisition_resets_gauge() {
    let p = products();
    let mut f = windup_filter();
    let epochs = windup_epochs(false);
    for epoch in &epochs[..8] {
        f.process_epoch(epoch, &p, &biases(), &selections())
            .unwrap();
    }
    let state = f.state().clone();
    let time = f.state_time();
    let sat = "G01".parse().unwrap();
    let cache = f.windup_state(sat).unwrap().clone();
    assert!(matches!(
        f.process_epoch(&epochs[7], &p, &biases(), &selections()),
        Err(PppError::Arc(PhaseArcError::OutOfOrder))
    ));
    assert_eq!(f.windup_state(sat), Some(&cache));
    let mut bad = epochs[8].clone();
    for obs in &mut bad.observations {
        if let Some(code) = obs.pseudorange_m.as_mut() {
            *code += 5000.0;
        }
    }
    assert!(matches!(
        f.process_epoch(&bad, &p, &biases(), &selections()),
        Err(PppError::Filter(FilterError::InnovationRejected))
    ));
    assert_eq!(f.state(), &state);
    assert_eq!(f.state_time(), time);
    assert!(f.windup_state(sat).is_none());
    let result = solution(
        f.process_epoch(&epochs[9], &p, &biases(), &selections())
            .unwrap(),
    );
    assert!(
        result
            .measurements
            .iter()
            .all(|m| m.windup.as_ref().unwrap().previous_cycles == 0.0)
    );
    assert!(
        result
            .arc_commit
            .accepted
            .values()
            .all(|arc| arc.accepted_epochs == 1)
    );
}
#[test]
fn slip_phase_outage_event_and_external_invalidation_clear_matching_windup_history() {
    let p = products();
    let mut f = windup_filter();
    let mut epochs = windup_epochs(false);
    let sat = "G01".parse().unwrap();
    for epoch in &epochs[..4] {
        f.process_epoch(epoch, &p, &biases(), &selections())
            .unwrap();
    }
    let other = f.windup_state("G02".parse().unwrap()).unwrap().cycles;
    for e in &mut epochs[4..] {
        for o in &mut e.observations {
            if o.satellite == sat && o.tracking_code == "1C" {
                o.carrier_phase_cycles = o.carrier_phase_cycles.map(|x| x + 4.0);
                if e.time == GnssTime::new(2300, 346320.0).unwrap() {
                    o.lli = Some(1);
                }
            }
        }
    }
    let result = solution(
        f.process_epoch(&epochs[4], &p, &biases(), &selections())
            .unwrap(),
    );
    assert_eq!(
        result
            .measurements
            .iter()
            .find(|m| m.satellite == sat)
            .unwrap()
            .windup
            .as_ref()
            .unwrap()
            .previous_cycles,
        0.0
    );
    assert_eq!(
        result
            .measurements
            .iter()
            .find(|m| m.satellite.to_string() == "G02")
            .unwrap()
            .windup
            .as_ref()
            .unwrap()
            .previous_cycles,
        other
    );
    let mut outage = epochs[5].clone();
    for o in &mut outage.observations {
        if o.satellite == sat && o.tracking_code == "2W" {
            o.carrier_phase_cycles = None;
        }
    }
    let result = solution(
        f.process_epoch(&outage, &p, &biases(), &selections())
            .unwrap(),
    );
    assert_eq!(result.ambiguities.len(), 8);
    assert!(f.windup_state(sat).is_none());
    assert!(
        result
            .measurements
            .iter()
            .find(|m| m.satellite == sat)
            .unwrap()
            .windup
            .is_none()
    );
    f.process_epoch(&epochs[6], &p, &biases(), &selections())
        .unwrap();
    let state = f.state().clone();
    let mut event = epochs[7].clone();
    event.flag = 6;
    event.observations.truncate(1);
    assert!(matches!(
        f.process_epoch(&event, &p, &biases(), &selections())
            .unwrap(),
        PppEpoch::Event(_)
    ));
    assert_eq!(f.state(), &state);
    assert!(f.windup_state(sat).is_none());
    f.process_epoch(&epochs[8], &p, &biases(), &selections())
        .unwrap();
    f.invalidate_phase().unwrap();
    assert!(f.windup_state(sat).is_none());
    f.process_epoch(&epochs[9], &p, &biases(), &selections())
        .unwrap();
    let mut power = epochs[10].clone();
    power.flag = 1;
    let r = solution(
        f.process_epoch(&power, &p, &biases(), &selections())
            .unwrap(),
    );
    assert!(
        r.measurements
            .iter()
            .all(|m| m.windup.as_ref().unwrap().previous_cycles == 0.0)
    );
    assert!(
        r.arc_commit
            .accepted
            .values()
            .all(|arc| arc.accepted_epochs == 1)
    );
}

#[test]
fn sinex_code_expiry_before_filter_requires_external_invalidation_and_recovers() {
    let text = include_str!("fixtures/synthetic_ppp_osb.bsx");
    let full = read_bias_sinex(Cursor::new(text), Default::default()).unwrap();
    let expired = read_bias_sinex(
        Cursor::new(text.replace("2024:039:02520", "2024:039:00720")),
        Default::default(),
    )
    .unwrap();
    let p = products();
    let mut f = PppFilter::new(config()).unwrap();
    let epochs = epochs(false);
    let requests = osb_requests();
    for epoch in &epochs[..4] {
        let s = full
            .corrections_at(epoch.time, &osb_binding(), &requests)
            .unwrap();
        f.process_epoch(epoch, &p, &s.corrections, &selections())
            .unwrap();
    }
    let state = f.state().clone();
    let state_time = f.state_time();
    assert!(matches!(
        expired.corrections_at(epochs[4].time, &osb_binding(), &requests),
        Err(gnss_rust::bias::BiasError::ExpiredOsb { .. })
    ));
    // The adapter is external to PppFilter: it cannot alter the filter on its
    // own. Explicitly interrupt continuity when skipping this failed epoch.
    f.invalidate_phase().unwrap();
    assert_eq!(f.state(), &state);
    assert_eq!(f.state_time(), state_time);
    assert!(f.phase_arc("G01".parse().unwrap()).is_none());
    let s = full
        .corrections_at(epochs[5].time, &osb_binding(), &requests)
        .unwrap();
    let r = solution(
        f.process_epoch(&epochs[5], &p, &s.corrections, &selections())
            .unwrap(),
    );
    assert_eq!(r.position.status, SolutionStatus::PppFloat);
    assert!(
        r.ambiguities
            .iter()
            .all(|a| f.phase_arc(a.satellite).unwrap().accepted_epochs == 1)
    );
    assert!(distance(r.position.position_ecef.unwrap()) < 0.001);
}

#[test]
fn sinex_phase_outage_keeps_code_and_datum_mismatch_rejects_numeric_update() {
    let text = include_str!("fixtures/synthetic_ppp_osb.bsx");
    let full = read_bias_sinex(Cursor::new(text), Default::default()).unwrap();
    let without_phase = text
        .lines()
        .filter(|line| !line.contains("G081 G01") || !line.contains("L2W"))
        .collect::<Vec<_>>()
        .join("\n");
    let outage = read_bias_sinex(Cursor::new(without_phase), Default::default()).unwrap();
    let p = products();
    let mut f = PppFilter::new(config()).unwrap();
    let epochs = epochs(false);
    let requests = osb_requests();
    for epoch in &epochs[..4] {
        let s = full
            .corrections_at(epoch.time, &osb_binding(), &requests)
            .unwrap();
        f.process_epoch(epoch, &p, &s.corrections, &selections())
            .unwrap();
    }
    let sat = "G01".parse().unwrap();
    let old = f.phase_arc(sat).unwrap().id;
    let s = outage
        .corrections_at(epochs[4].time, &osb_binding(), &requests)
        .unwrap();
    let result = solution(
        f.process_epoch(&epochs[4], &p, &s.corrections, &selections())
            .unwrap(),
    );
    assert_eq!(result.position.num_satellites, 9);
    assert_eq!(result.ambiguities.len(), 8);
    assert!(f.phase_arc(sat).is_none());
    let measurement = result
        .measurements
        .iter()
        .find(|m| m.satellite == sat)
        .unwrap();
    assert!(matches!(
        measurement.observation.phase_unavailability,
        Some(PhaseUnavailable::MissingCorrection { .. })
    ));
    let s = full
        .corrections_at(epochs[5].time, &osb_binding(), &requests)
        .unwrap();
    f.process_epoch(&epochs[5], &p, &s.corrections, &selections())
        .unwrap();
    assert_ne!(f.phase_arc(sat).unwrap().id, old);
    assert_eq!(f.phase_arc(sat).unwrap().accepted_epochs, 1);
    assert_eq!(
        f.phase_arc("G02".parse().unwrap()).unwrap().accepted_epochs,
        6
    );
    let state = f.state().clone();
    let state_time = f.state_time();
    let mut wrong = osb_binding();
    wrong.clock_reference = "incompatible-clock".into();
    let s = full
        .corrections_at(epochs[6].time, &wrong, &requests)
        .unwrap();
    assert!(matches!(
        f.process_epoch(&epochs[6], &p, &s.corrections, &selections()),
        Err(PppError::InsufficientCode { available: 0 })
    ));
    assert_eq!(f.state(), &state);
    assert_eq!(f.state_time(), state_time);
    assert!(f.phase_arc(sat).is_none());
}

#[test]
fn sinex_adjacent_phase_bias_change_reinitializes_only_affected_satellite() {
    let text = include_str!("fixtures/synthetic_ppp_osb.bsx");
    let mut rows = Vec::new();
    for line in text.lines() {
        if line.contains("G081 G01") && line.contains("L1C") {
            rows.push(line.replace("2024:039:02520", "2024:039:00840"));
            rows.push(" OSB G081 G01 L1C 2024:039:00840 2024:039:02520 m 0.028 0.00001".into());
        } else {
            rows.push(line.to_owned());
        }
    }
    let osbs = read_bias_sinex(Cursor::new(rows.join("\n")), Default::default()).unwrap();
    let p = products();
    let mut f = PppFilter::new(config()).unwrap();
    let epochs = epochs(false);
    let requests = osb_requests();
    for epoch in &epochs[..8] {
        let s = osbs
            .corrections_at(epoch.time, &osb_binding(), &requests)
            .unwrap();
        f.process_epoch(epoch, &p, &s.corrections, &selections())
            .unwrap();
    }
    let old: BTreeMap<_, _> = selections()
        .keys()
        .map(|&sat| (sat, f.phase_arc(sat).unwrap().id))
        .collect();
    let sat = "G01".parse().unwrap();
    let s = osbs
        .corrections_at(epochs[8].time, &osb_binding(), &requests)
        .unwrap();
    let result = solution(
        f.process_epoch(&epochs[8], &p, &s.corrections, &selections())
            .unwrap(),
    );
    assert_eq!(result.position.status, SolutionStatus::PppFloat);
    assert_ne!(f.phase_arc(sat).unwrap().id, old[&sat]);
    assert_eq!(f.phase_arc(sat).unwrap().accepted_epochs, 1);
    for (&other, &arc) in &old {
        if other != sat {
            assert_eq!(f.phase_arc(other).unwrap().id, arc);
            assert_eq!(f.phase_arc(other).unwrap().accepted_epochs, 9);
        }
    }
    let a = result
        .ambiguities
        .iter()
        .find(|a| a.satellite == sat)
        .unwrap();
    close(a.metres, ambiguities()[&sat] - 0.02 * c::IFLC_F1, 0.001);
    assert!(distance(result.position.position_ecef.unwrap()) < 0.001);
}

#[test]
fn sinex_osbs_recover_both_complete_ppp_traces_without_breaking_continuous_arcs() {
    let p = products();
    let osbs = read_bias_sinex(
        Cursor::new(include_str!("fixtures/synthetic_ppp_osb.bsx")),
        Default::default(),
    )
    .unwrap();
    let requests = osb_requests();
    for noisy in [false, true] {
        let mut loaded = PppFilter::new(config()).unwrap();
        let mut ledger = PppFilter::new(config()).unwrap();
        let mut count = 0;
        for (i, epoch) in epochs(noisy).iter().enumerate() {
            let snapshot = osbs
                .corrections_at(epoch.time, &osb_binding(), &requests)
                .unwrap();
            assert_eq!(snapshot.provenance.len(), 18);
            let a = solution(
                loaded
                    .process_epoch(epoch, &p, &snapshot.corrections, &selections())
                    .unwrap(),
            );
            let b = solution(
                ledger
                    .process_epoch(epoch, &p, &biases(), &selections())
                    .unwrap(),
            );
            assert_eq!(a.position.status, SolutionStatus::PppFloat);
            assert_eq!(a.arc_commit, b.arc_commit);
            // Snapshot validity is this single epoch; the old ledger spans
            // all epochs. Compare actual lifecycle decisions independently
            // of this intentionally different correction provenance.
            for (sat, decision) in &a.arc_decisions {
                let other = &b.arc_decisions[sat];
                assert_eq!(decision.proposed_arc, other.proposed_arc);
                assert_eq!(decision.reset_reasons, other.reset_reasons);
                assert_eq!(
                    decision.previous_accepted_epochs,
                    other.previous_accepted_epochs
                );
                assert_eq!(decision.combination_slip, other.combination_slip);
            }
            for (x, y) in loaded.state().mean.iter().zip(&ledger.state().mean) {
                close(*x, *y, 1e-7);
            }
            for (x, y) in loaded
                .state()
                .covariance
                .iter()
                .flatten()
                .zip(ledger.state().covariance.iter().flatten())
            {
                close(*x, *y, 1e-7);
            }
            for used in &a.measurements {
                assert_eq!(used.observation.corrections[0].valid_from, epoch.time);
                let sat = used.observation.satellite;
                assert_eq!(loaded.phase_arc(sat).unwrap().accepted_epochs, i as u64 + 1);
            }
            count += 1;
            if i == 63 {
                assert!(
                    distance(a.position.position_ecef.unwrap()) < if noisy { 0.01 } else { 0.001 }
                );
            }
        }
        assert_eq!(count, 64);
    }
}
fn truth() -> Vec<f64> {
    include_str!("fixtures/synthetic_ppp_truth.csv")
        .lines()
        .nth(1)
        .unwrap()
        .split(',')
        .map(|v| v.parse().unwrap())
        .collect()
}
fn ambiguities() -> BTreeMap<SatelliteId, f64> {
    include_str!("fixtures/synthetic_ppp_truth.csv")
        .lines()
        .skip(3)
        .map(|line| {
            let v: Vec<_> = line.split(',').collect();
            (v[0].parse().unwrap(), v[3].parse().unwrap())
        })
        .collect()
}
fn bias_entries() -> Vec<SignalCorrection> {
    include_str!("fixtures/synthetic_ppp_bias.csv")
        .lines()
        .skip(1)
        .map(|line| {
            let v: Vec<_> = line.split(',').collect();
            SignalCorrection {
                satellite: v[0].parse().unwrap(),
                tracking_code: v[1].to_owned(),
                clock_source: ClockSource::RinexClk,
                valid_from: GnssTime::new(v[3].parse().unwrap(), v[4].parse().unwrap()).unwrap(),
                valid_until: GnssTime::new(v[5].parse().unwrap(), v[6].parse().unwrap()).unwrap(),
                code_add_m: v[7].parse().unwrap(),
                phase_add_m: Some(v[8].parse().unwrap()),
                clock_reference: v[9].to_owned(),
                phase_reference: Some(v[10].to_owned()),
            }
        })
        .collect()
}
fn biases() -> SignalCorrections {
    SignalCorrections::new(bias_entries()).unwrap()
}
fn selections() -> BTreeMap<SatelliteId, DualFrequencyConfig> {
    ambiguities()
        .keys()
        .map(|&sat| {
            (
                sat,
                DualFrequencyConfig {
                    tracking_codes: ["1C".to_owned(), "2W".to_owned()],
                    clock_source: ClockSource::RinexClk,
                    clock_reference: "synthetic-precise-clock-v1".to_owned(),
                    phase_reference: "synthetic-phase-v1".to_owned(),
                    noise: DualFrequencyNoise {
                        code_m2: [[0.0025, 0.0], [0.0, 0.0025]],
                        phase_m2: [[0.000001, 0.0], [0.0, 0.000001]],
                    },
                },
            )
        })
        .collect()
}
fn config() -> PppConfig {
    let t = truth();
    PppConfig {
        initial_position: EcefCoord::new(t[3], t[4], t[5]),
        initial_receiver_clock_s: t[6],
        initial_zenith_delay_m: 2.3,
        expected_coordinate_system: "IGS20".to_owned(),
        expected_reference_point: OrbitReferencePoint::AntennaPhaseCenter,
        transmit: Default::default(),
        arcs: Default::default(),
        filter: Default::default(),
        phase_windup: None,
        receiver_antenna: None,
        satellite_antenna: None,
        receiver_reference: gnss_rust::tides::ReceiverPositionReference::InstantaneousMarker,
        solid_earth_tide: None,
    }
}
fn epochs(noisy: bool) -> Vec<gnss_rust::observation::ObservationEpoch> {
    RinexObservationReader::new(Cursor::new(if noisy {
        include_str!("fixtures/synthetic_ppp_noisy.obs")
    } else {
        include_str!("fixtures/synthetic_ppp.obs")
    }))
    .unwrap()
    .collect::<Result<_, _>>()
    .unwrap()
}
fn solution(result: PppEpoch) -> Box<PppSolution> {
    match result {
        PppEpoch::Solution(s) => s,
        PppEpoch::Event(_) => panic!("expected solution"),
    }
}
fn distance(p: EcefCoord) -> f64 {
    let t = truth();
    (p.x - t[0]).hypot(p.y - t[1]).hypot(p.z - t[2])
}
fn admitted() -> Vec<std::collections::BTreeSet<SatelliteId>> {
    let mut epochs = vec![std::collections::BTreeSet::new(); 64];
    for line in include_str!("fixtures/upstream_ppp_admission.csv")
        .lines()
        .skip(1)
    {
        let v: Vec<_> = line.split(',').collect();
        if v[3] == "1" {
            epochs[v[0].parse::<usize>().unwrap()].insert(v[1].parse().unwrap());
        }
    }
    epochs
}
fn close(actual: f64, expected: f64, tol: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tol,
        "{actual:.17e} vs {expected:.17e}, tol={tol}"
    );
}

#[test]
fn native_54_precise_models_match_geometry_clock_mapping_and_rows() {
    let p = products();
    let config = config();
    let mut count = 0;
    for line in include_str!("fixtures/upstream_ppp_model.csv")
        .lines()
        .skip(1)
    {
        let v: Vec<_> = line.split(',').collect();
        let n = |i: usize| v[i].parse::<f64>().unwrap();
        let model = evaluate_ppp_model(
            &p,
            v[0].parse().unwrap(),
            GnssTime::new(v[1].parse().unwrap(), n(2)).unwrap(),
            EcefCoord::new(n(3), n(4), n(5)),
            n(6),
            n(7),
            &config,
        )
        .unwrap();
        close(model.geometric_range_m, n(8), 1e-6);
        close(model.transmit.corrected_clock_bias_s, n(9), 1e-16);
        close(model.elevation_rad, n(10), 1e-11);
        close(model.troposphere_mapping, n(11), 1e-10);
        close(model.code_prediction_m, n(12), 1e-6);
        for i in 0..5 {
            close(model.head_design[i], n(13 + i), 1e-10);
        }
        count += 1;
    }
    assert_eq!(count, 54);
}

#[test]
fn native_generated_static_noiseless_and_noisy_phase_epochs_recover_truth() {
    let p = products();
    let b = biases();
    let s = selections();
    let truth_amb = ambiguities();
    let admission = admitted();
    for noisy in [false, true] {
        let mut f = PppFilter::new(config()).unwrap();
        let mut max_position = 0.0_f64;
        let mut last = None;
        for (i, e) in epochs(noisy).iter().enumerate() {
            let result = solution(f.process_epoch(e, &p, &b, &s).unwrap());
            let error = distance(result.position.position_ecef.unwrap());
            max_position = max_position.max(error);
            assert_eq!(result.position.status, SolutionStatus::PppFloat);
            assert!(result.position.is_valid());
            assert!(!result.position.is_fixed());
            assert_eq!(result.measurements.len(), admission[i].len());
            assert_eq!(result.ambiguities.len(), admission[i].len());
            assert_eq!(result.position.num_satellites, admission[i].len());
            assert_eq!(
                result
                    .measurements
                    .iter()
                    .map(|m| m.satellite)
                    .collect::<std::collections::BTreeSet<_>>(),
                admission[i]
            );
            assert!(
                result
                    .rejected
                    .values()
                    .all(|r| *r == PppRejection::BelowElevationMask)
            );
            assert!(
                result
                    .position
                    .position_covariance
                    .unwrap()
                    .iter()
                    .flatten()
                    .all(|v| v.is_finite())
            );
            for ambiguity in &result.ambiguities {
                assert!(ambiguity.variance_m2 > 0.0);
                let length = admission[..=i]
                    .iter()
                    .rev()
                    .take_while(|set| set.contains(&ambiguity.satellite))
                    .count();
                assert_eq!(
                    f.phase_arc(ambiguity.satellite).unwrap().accepted_epochs,
                    length as u64
                );
            }
            last = Some(result);
        }
        let result = last.unwrap();
        let final_position = distance(result.position.position_ecef.unwrap());
        let final_clock =
            (result.position.receiver_clock_bias - (2e-6 + 63.0 * 5e-10)) * c::SPEED_OF_LIGHT;
        let max_ambiguity = result
            .ambiguities
            .iter()
            .map(|a| (a.metres - truth_amb[&a.satellite]).abs())
            .fold(0.0, f64::max);
        eprintln!(
            "noisy={noisy} max_pos={max_position:.9} final_pos={final_position:.9} clock_m={final_clock:.9} zenith={:.9} amb={max_ambiguity:.9} code_rms={:.9} phase_rms={:.9}",
            result.zenith_delay_m, result.code_rms_m, result.phase_rms_m
        );
        let final_bound = if noisy { 0.01 } else { 0.001 };
        assert!(max_position < if noisy { 1.0 } else { 0.002 });
        assert!(final_position < final_bound);
        assert!(final_clock.abs() < final_bound);
        assert!((result.zenith_delay_m - 2.4).abs() < final_bound);
        assert!(max_ambiguity < final_bound);
        assert!(result.phase_rms_m < if noisy { 0.005 } else { 0.0001 });
        assert!(result.code_rms_m < if noisy { 0.2 } else { 0.0001 });
    }
}

#[test]
fn initial_iteration_matches_single_linear_gaussian_information_update() {
    let p = products();
    let b = biases();
    let s = selections();
    let mut cfg = config();
    let t = truth();
    cfg.initial_position = EcefCoord::new(t[0], t[1], t[2]);
    cfg.initial_receiver_clock_s = 2e-6;
    cfg.initial_zenith_delay_m = 2.4;
    let mut f = PppFilter::new(cfg.clone()).unwrap();
    let initial = f.state().clone();
    let e = &epochs(false)[0];
    let result = solution(f.process_epoch(e, &p, &b, &s).unwrap());
    let mut prior = initial;
    let count = 5 + s.len();
    prior.mean.resize(count, 0.0);
    for r in &mut prior.covariance {
        r.resize(count, 0.0);
    }
    prior.covariance.resize(count, vec![0.0; count]);
    for i in 5..count {
        prior.covariance[i][i] = cfg.filter.initial_ambiguity_variance_m2;
    }
    let mut h = Vec::new();
    let mut residuals = Vec::new();
    let mut variance = Vec::new();
    for (j, (&sat, config)) in s.iter().enumerate() {
        let obs = form_ionosphere_free(e, sat, &b, config).unwrap();
        let model = evaluate_ppp_model(
            &p,
            sat,
            e.time,
            cfg.initial_position,
            2e-6 * c::SPEED_OF_LIGHT,
            2.4,
            &cfg,
        )
        .unwrap();
        let mut row = vec![0.0; count];
        row[..5].copy_from_slice(&model.head_design);
        h.push(row.clone());
        residuals.push(obs.code_if_m - model.code_prediction_m);
        variance.push(obs.code_variance_m2);
        row[5 + j] = 1.0;
        h.push(row);
        residuals.push(obs.phase.unwrap().phase_if_m - model.code_prediction_m);
        variance.push(
            config.noise.phase_m2[0][0] * (c::IFLC_F1 * c::IFLC_F1 + c::IFLC_F2 * c::IFLC_F2),
        );
    }
    let mut noise = vec![vec![0.0; variance.len()]; variance.len()];
    for (i, v) in variance.iter().enumerate() {
        noise[i][i] = *v;
    }
    let once = measurement_update(
        &prior,
        &h,
        &residuals,
        &noise,
        Some(cfg.filter.max_nis_per_row),
    )
    .unwrap()
    .posterior;
    // Nonlinear relinearization can change the final mean slightly; it must not
    // shrink covariance as though the same epoch was measured twice.
    for i in 0..count {
        close(f.state().covariance[i][i], once.covariance[i][i], 1e-5);
    }
    assert!(result.iterations >= 2);
}

fn warm(config: PppConfig) -> PppFilter {
    let mut f = PppFilter::new(config).unwrap();
    let p = products();
    let b = biases();
    let s = selections();
    for e in epochs(false).iter().take(8) {
        solution(f.process_epoch(e, &p, &b, &s).unwrap());
    }
    f
}
fn first_satellite() -> SatelliteId {
    *selections().keys().next().unwrap()
}

#[test]
fn innovation_rejection_preserves_numeric_state_and_requires_fresh_ambiguities() {
    let mut f = warm(config());
    let before = f.state().clone();
    let time = f.state_time();
    let mut bad = epochs(false)[8].clone();
    for obs in &mut bad.observations {
        obs.pseudorange_m = obs.pseudorange_m.map(|p| p + 5000.0);
    }
    assert_eq!(
        f.process_epoch(&bad, &products(), &biases(), &selections()),
        Err(PppError::Filter(FilterError::InnovationRejected))
    );
    assert_eq!(f.state(), &before);
    assert_eq!(f.state_time(), time);
    for sat in selections().keys() {
        assert!(f.phase_arc(*sat).is_none());
    }
    assert_eq!(
        f.process_epoch(&epochs(false)[8], &products(), &biases(), &selections()),
        Err(PppError::Arc(PhaseArcError::OutOfOrder))
    );
    let recovered = solution(
        f.process_epoch(&epochs(false)[9], &products(), &biases(), &selections())
            .unwrap(),
    );
    assert!(distance(recovered.position.position_ecef.unwrap()) < 0.01);
    for a in &recovered.ambiguities {
        assert_eq!(f.phase_arc(a.satellite).unwrap().accepted_epochs, 1);
    }
}

#[test]
fn postfit_rejection_does_not_commit_a_numerically_successful_update() {
    let mut cfg = config();
    cfg.filter.max_postfit_code_m = 0.1;
    let mut f = warm(cfg);
    let before = f.state().clone();
    let mut bad = epochs(false)[8].clone();
    bad.observations[0].pseudorange_m = bad.observations[0].pseudorange_m.map(|p| p + 0.25);
    assert_eq!(
        f.process_epoch(&bad, &products(), &biases(), &selections()),
        Err(PppError::PostfitRejected)
    );
    assert_eq!(f.state(), &before);
    assert!(f.phase_arc(first_satellite()).is_none());
}

#[test]
fn announced_and_unannounced_slips_reset_only_the_affected_ambiguity() {
    for announced in [false, true] {
        let mut f = warm(config());
        let p = products();
        let b = biases();
        let s = selections();
        let sat = first_satellite();
        let old = f.phase_arc(sat).unwrap().id;
        let cycles = if announced { 4.0 } else { 5.0 };
        for i in 8..12 {
            let mut e = epochs(false)[i].clone();
            let obs = e
                .observations
                .iter_mut()
                .find(|o| o.satellite == sat && o.tracking_code == "1C")
                .unwrap();
            obs.carrier_phase_cycles = obs.carrier_phase_cycles.map(|v| v + cycles);
            if announced && i == 8 {
                obs.lli = Some(1);
            }
            let result = solution(f.process_epoch(&e, &p, &b, &s).unwrap());
            if i == 8 {
                assert!(
                    result.arc_decisions[&sat]
                        .reset_reasons
                        .contains(&if announced {
                            ArcResetReason::LossOfLock
                        } else {
                            ArcResetReason::GeometryFreeJump
                        })
                );
            }
            let a = result
                .ambiguities
                .iter()
                .find(|a| a.satellite == sat)
                .unwrap();
            assert_ne!(a.arc_id, old);
            close(
                a.metres,
                ambiguities()[&sat] + c::IFLC_F1 * c::GPS_L1_WAVELENGTH * cycles,
                0.01,
            );
            assert_eq!(f.phase_arc(sat).unwrap().accepted_epochs, (i - 7) as u64);
            for other in result.ambiguities.iter().filter(|a| a.satellite != sat) {
                assert_eq!(
                    f.phase_arc(other.satellite).unwrap().accepted_epochs,
                    i as u64 + 1
                );
            }
        }
    }
}

#[test]
fn phase_outage_preserves_code_and_marginalizes_old_ambiguity() {
    let mut f = warm(config());
    let sat = first_satellite();
    let old = f.phase_arc(sat).unwrap().id;
    let mut e = epochs(false)[8].clone();
    e.observations
        .iter_mut()
        .find(|o| o.satellite == sat && o.tracking_code == "2W")
        .unwrap()
        .carrier_phase_cycles = None;
    let result = solution(
        f.process_epoch(&e, &products(), &biases(), &selections())
            .unwrap(),
    );
    assert_eq!(result.measurements.len(), 9);
    assert_eq!(result.ambiguities.len(), 8);
    assert_eq!(f.state().mean.len(), 13);
    let snapshot = result
        .measurements
        .iter()
        .find(|m| m.satellite == sat)
        .unwrap();
    assert!(snapshot.phase_residual_m.is_none());
    assert!(snapshot.code_residual_m.is_finite());
    assert!(f.phase_arc(sat).is_none());
    let next = solution(
        f.process_epoch(&epochs(false)[9], &products(), &biases(), &selections())
            .unwrap(),
    );
    let a = next
        .ambiguities
        .iter()
        .find(|a| a.satellite == sat)
        .unwrap();
    assert_ne!(a.arc_id, old);
    assert_eq!(f.phase_arc(sat).unwrap().accepted_epochs, 1);
    assert_eq!(f.state().mean.len(), 14);
}

#[test]
fn code_only_epochs_do_not_get_a_ppp_status_or_change_numeric_state() {
    let mut f = warm(config());
    let before = f.state().clone();
    let mut e = epochs(false)[8].clone();
    for o in &mut e.observations {
        o.carrier_phase_cycles = None;
    }
    assert_eq!(
        f.process_epoch(&e, &products(), &biases(), &selections()),
        Err(PppError::InsufficientPhase { available: 0 })
    );
    assert_eq!(f.state(), &before);
    assert!(f.phase_arc(first_satellite()).is_none());
}

#[test]
fn sparse_events_and_receiver_power_failure_reinitialize_arcs() {
    let mut f = warm(config());
    let before = f.state().clone();
    let mut event = epochs(false)[8].clone();
    event.flag = 6;
    event
        .observations
        .retain(|o| o.satellite == first_satellite());
    let result = f
        .process_epoch(&event, &products(), &biases(), &selections())
        .unwrap();
    match result {
        PppEpoch::Event(commit) => assert_eq!(commit.interrupted.len(), 9),
        _ => panic!("expected event"),
    }
    assert_eq!(f.state(), &before);
    assert_eq!(f.state_time(), Some(epochs(false)[7].time));
    let next = solution(
        f.process_epoch(&epochs(false)[9], &products(), &biases(), &selections())
            .unwrap(),
    );
    for a in &next.ambiguities {
        assert_eq!(f.phase_arc(a.satellite).unwrap().accepted_epochs, 1);
    }
    let mut power = epochs(false)[10].clone();
    power.flag = 1;
    let result = solution(
        f.process_epoch(&power, &products(), &biases(), &selections())
            .unwrap(),
    );
    for a in &result.ambiguities {
        assert!(
            result.arc_decisions[&a.satellite]
                .reset_reasons
                .contains(&ArcResetReason::PowerFailure)
        );
        assert_eq!(f.phase_arc(a.satellite).unwrap().accepted_epochs, 1);
    }
}

#[test]
fn actual_clock_source_mismatch_rejects_the_satellite_without_broadcast_fallback() {
    let mut f = warm(config());
    let sat = first_satellite();
    let mut s = selections();
    s.get_mut(&sat).unwrap().clock_source = ClockSource::Sp3;
    let mut entries = bias_entries();
    for row in entries.iter_mut().filter(|r| r.satellite == sat) {
        row.clock_source = ClockSource::Sp3;
    }
    let result = solution(
        f.process_epoch(
            &epochs(false)[8],
            &products(),
            &SignalCorrections::new(entries).unwrap(),
            &s,
        )
        .unwrap(),
    );
    assert_eq!(result.rejected[&sat], PppRejection::ClockSourceMismatch);
    assert_eq!(result.measurements.len(), 8);
    assert!(f.phase_arc(sat).is_none());
}

#[test]
fn missing_selected_clk_is_isolated_even_when_sp3_clocks_are_valid() {
    let mut orbit = String::new();
    for line in include_str!("fixtures/synthetic_ppp.sp3").lines() {
        let mut text = line.to_owned();
        if text.starts_with("PG01") {
            text.replace_range(46..60, &format!("{:14.6}", 20.0));
        }
        orbit.push_str(&text);
        orbit.push('\n');
    }
    let clock = include_str!("fixtures/synthetic_ppp.clk")
        .lines()
        .filter(|line| {
            if !line.starts_with("AS G01 ") {
                return true;
            }
            let v: Vec<_> = line.split_whitespace().collect();
            let seconds = v[5].parse::<i32>().unwrap() * 3600
                + v[6].parse::<i32>().unwrap() * 60
                + v[7].parse::<f64>().unwrap() as i32;
            !(900..=2250).contains(&seconds)
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let p = PreciseProducts::new(
        read_sp3(Cursor::new(orbit)).unwrap(),
        Some(read_clk(Cursor::new(clock)).unwrap()),
    );
    let mut f = warm(config());
    let result = solution(
        f.process_epoch(&epochs(false)[8], &p, &biases(), &selections())
            .unwrap(),
    );
    let sat = first_satellite();
    assert!(matches!(
        result.rejected[&sat],
        PppRejection::Model(PppError::Transmit(_))
    ));
    assert_eq!(result.measurements.len(), 8);
    assert!(f.phase_arc(sat).is_none());
}

#[test]
fn global_frame_and_epoch_errors_preserve_numeric_state_and_fail_closed() {
    for kind in 0..2 {
        let mut f = warm(config());
        let before = f.state().clone();
        let mut e = epochs(false)[8].clone();
        let p = if kind == 0 {
            PreciseProducts::new(
                read_sp3(Cursor::new(
                    include_str!("fixtures/synthetic_ppp.sp3").replace("IGS20", "ITRF0"),
                ))
                .unwrap(),
                Some(read_clk(Cursor::new(include_str!("fixtures/synthetic_ppp.clk"))).unwrap()),
            )
        } else {
            e.receiver_clock_offset_s = Some(1e-6);
            products()
        };
        let error = f
            .process_epoch(&e, &p, &biases(), &selections())
            .unwrap_err();
        assert_eq!(
            error,
            if kind == 0 {
                PppError::ProductFrameMismatch
            } else {
                PppError::Arc(PhaseArcError::InvalidEpoch)
            }
        );
        assert_eq!(f.state(), &before);
        assert!(f.phase_arc(first_satellite()).is_none());
        solution(
            f.process_epoch(&epochs(false)[8], &products(), &biases(), &selections())
                .unwrap(),
        );
    }
}

#[test]
fn insufficient_or_singular_geometry_is_not_hidden_by_a_finite_prior() {
    let mut e = epochs(false)[0].clone();
    let s = selections();
    let keep: Vec<_> = s.keys().take(4).copied().collect();
    e.observations.retain(|o| keep.contains(&o.satellite));
    let mut f = PppFilter::new(config()).unwrap();
    let before = f.state().clone();
    assert_eq!(
        f.process_epoch(&e, &products(), &biases(), &s),
        Err(PppError::InsufficientCode { available: 4 })
    );
    assert_eq!(f.state(), &before);
    let mut orbit = String::new();
    let mut first = None;
    for line in include_str!("fixtures/synthetic_ppp.sp3").lines() {
        let mut text = line.to_owned();
        if text.starts_with("PG01") {
            first = Some(text[4..46].to_owned());
        }
        if text.starts_with('P') {
            text.replace_range(4..46, first.as_ref().unwrap());
        }
        orbit.push_str(&text);
        orbit.push('\n');
    }
    let p = PreciseProducts::new(
        read_sp3(Cursor::new(orbit)).unwrap(),
        Some(read_clk(Cursor::new(include_str!("fixtures/synthetic_ppp.clk"))).unwrap()),
    );
    let mut f = PppFilter::new(config()).unwrap();
    assert_eq!(
        f.process_epoch(&epochs(false)[0], &p, &biases(), &s),
        Err(PppError::Gnss(Error::SingularGeometry))
    );
}

#[test]
fn nonconvergence_invalid_settings_and_unsupported_constellations_are_explicit() {
    let mut cfg = config();
    cfg.filter.max_iterations = 1;
    let mut f = PppFilter::new(cfg).unwrap();
    let before = f.state().clone();
    assert_eq!(
        f.process_epoch(&epochs(false)[0], &products(), &biases(), &selections()),
        Err(PppError::NonConvergent)
    );
    assert_eq!(f.state(), &before);
    assert!(f.phase_arc(first_satellite()).is_none());
    for kind in 0..5 {
        let mut cfg = config();
        match kind {
            0 => cfg.filter.initial_head_variance_m2[0] = 0.0,
            1 => cfg.filter.head_process_noise_m2_per_s[0] = -1.0,
            2 => cfg.initial_zenith_delay_m = -1.0,
            3 => cfg.filter.max_postfit_code_m = f64::NAN,
            _ => cfg.filter.max_iterations = 0,
        }
        assert!(matches!(
            PppFilter::new(cfg),
            Err(PppError::InvalidConfiguration)
        ));
    }
    let mut f = warm(config());
    let before = f.state().clone();
    let mut s = selections();
    let foreign: SatelliteId = "J01".parse().unwrap();
    s.insert(foreign, s[&first_satellite()].clone());
    assert_eq!(
        f.process_epoch(&epochs(false)[8], &products(), &biases(), &s),
        Err(PppError::UnsupportedSelection(foreign))
    );
    assert_eq!(f.state(), &before);
}

#[test]
fn duplicate_or_backward_calls_cannot_double_count_a_successful_epoch() {
    let mut f = warm(config());
    let before = f.state().clone();
    let arc = f.phase_arc(first_satellite()).unwrap().clone();
    for i in [6, 7] {
        assert_eq!(
            f.process_epoch(&epochs(false)[i], &products(), &biases(), &selections()),
            Err(PppError::Arc(PhaseArcError::OutOfOrder))
        );
        assert_eq!(f.state(), &before);
        assert_eq!(f.phase_arc(first_satellite()), Some(&arc));
    }
}

#[test]
fn calibration_changes_reset_bias_state_and_expose_applied_corrections() {
    let mut f = warm(config());
    let sat = first_satellite();
    let old = f.phase_arc(sat).unwrap().id;
    let mut entries = bias_entries();
    entries
        .iter_mut()
        .find(|r| r.satellite == sat && r.tracking_code == "1C")
        .unwrap()
        .phase_add_m = Some(0.02);
    let table = SignalCorrections::new(entries).unwrap();
    let result = solution(
        f.process_epoch(&epochs(false)[8], &products(), &table, &selections())
            .unwrap(),
    );
    let a = result
        .ambiguities
        .iter()
        .find(|a| a.satellite == sat)
        .unwrap();
    assert_ne!(a.arc_id, old);
    assert_eq!(f.phase_arc(sat).unwrap().accepted_epochs, 1);
    assert!(
        result.arc_decisions[&sat]
            .reset_reasons
            .contains(&ArcResetReason::PhaseCorrectionChanged)
    );
    let m = result
        .measurements
        .iter()
        .find(|m| m.satellite == sat)
        .unwrap();
    assert_eq!(m.observation.corrections[0].phase_add_m, Some(0.02));
    close(
        m.code_residual_m,
        m.observation.code_if_m - m.model.code_prediction_m,
        0.0,
    );
    close(
        m.phase_residual_m.unwrap(),
        m.observation.phase.as_ref().unwrap().phase_if_m - m.model.code_prediction_m - a.metres,
        1e-8,
    );
}

#[test]
fn large_initial_carrier_gauge_is_absorbed_without_biasing_position() {
    let p = products();
    let b = biases();
    let s = selections();
    let sat = first_satellite();
    let mut f = PppFilter::new(config()).unwrap();
    // A receiver may start its accumulated phase near zero rather than the
    // geometric range. A new IF ambiguity must support that large metre gauge.
    let shift_m = -20_000_000.0;
    for e in epochs(false).iter().take(12) {
        let mut changed = e.clone();
        for obs in changed
            .observations
            .iter_mut()
            .filter(|o| o.satellite == sat)
        {
            let lambda = if obs.tracking_code == "1C" {
                c::GPS_L1_WAVELENGTH
            } else {
                c::GPS_L2_WAVELENGTH
            };
            obs.carrier_phase_cycles = obs.carrier_phase_cycles.map(|v| v + shift_m / lambda);
        }
        let result = solution(f.process_epoch(&changed, &p, &b, &s).unwrap());
        assert!(distance(result.position.position_ecef.unwrap()) < 0.002);
        let a = result
            .ambiguities
            .iter()
            .find(|a| a.satellite == sat)
            .unwrap();
        close(a.metres, ambiguities()[&sat] + shift_m, 0.002);
    }
}
