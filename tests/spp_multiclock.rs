use gnss_rust::observation::ObservationEpoch;
use gnss_rust::rinex::{RinexObservationReader, read_navigation};
use gnss_rust::spp::{self, SppConfig};
use gnss_rust::spp_multiclock::{ObservationDisposition, solve_gps_galileo_epoch};
use gnss_rust::{EcefCoord, Error, GnssSystem};
use std::io::Cursor;

fn config() -> SppConfig {
    SppConfig {
        use_ionosphere: false,
        ..Default::default()
    }
}
fn data() -> (Vec<ObservationEpoch>, gnss_rust::navigation::NavigationData) {
    let r = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_mixed_rover.obs"
    )))
    .unwrap();
    (
        r.collect::<Result<_, _>>().unwrap(),
        read_navigation(Cursor::new(include_str!(
            "fixtures/synthetic_rtk_mixed.nav"
        )))
        .unwrap(),
    )
}
fn truth(i: usize) -> EcefCoord {
    let v: Vec<f64> = include_str!("fixtures/synthetic_rtk_mixed_truth.csv")
        .lines()
        .nth(i + 1)
        .unwrap()
        .split(',')
        .map(|v| v.parse().unwrap())
        .collect();
    EcefCoord::new(v[2], v[3], v[4])
}
fn distance(a: EcefCoord, b: EcefCoord) -> f64 {
    (a.x - b.x).hypot(a.y - b.y).hypot(a.z - b.z)
}

#[test]
fn mixed_spp_recovers_all_32_native_physical_positions_and_distinct_clocks() {
    let (epochs, nav) = data();
    for (i, epoch) in epochs.iter().enumerate() {
        let before = epoch.clone();
        let s = solve_gps_galileo_epoch(epoch, &nav, config()).unwrap();
        assert!(
            distance(s.spp.position.position_ecef.unwrap(), truth(i)) < 0.003,
            "epoch {i}"
        );
        assert_eq!(s.clock_biases_s.len(), 2);
        assert_eq!(s.full_covariance_m2.len(), 5);
        assert!((s.clock_biases_s[0].1 - 2e-6).abs() < 1e-11);
        assert!((s.clock_biases_s[1].1 + 3e-6).abs() < 1e-11);
        assert_eq!(s.clock_biases_s[0].0, GnssSystem::Gps);
        assert_eq!(s.clock_biases_s[1].0, GnssSystem::Galileo);
        assert!(s.spp.residual_rms_m < 0.001);
        assert_eq!(s.broadcast_models.len(), 8);
        assert_eq!(s.observations.len(), epoch.observations.len());
        assert_eq!(*epoch, before);
        for (raw, evidence) in epoch.observations.iter().zip(&s.observations) {
            assert_eq!(*raw, evidence.raw_observation);
            assert_eq!(evidence.disposition, ObservationDisposition::Used);
        }
    }
}

#[test]
fn three_gps_plus_three_galileo_supply_five_state_position_information() {
    let (epochs, nav) = data();
    let mut epoch = epochs[0].clone();
    let mut counts = std::collections::BTreeMap::new();
    epoch.observations.retain(|o| {
        let n = counts.entry(o.satellite.system()).or_insert(0);
        *n += 1;
        *n <= 3
    });
    assert_eq!(
        spp::solve_epoch(&epoch, &nav, config()),
        Err(Error::InsufficientObservations)
    );
    let s = solve_gps_galileo_epoch(&epoch, &nav, config()).unwrap();
    assert_eq!(s.spp.position.num_satellites, 6);
    assert!(distance(s.spp.position.position_ecef.unwrap(), truth(0)) < 0.003);
}

#[test]
fn a_singleton_constellation_adds_only_a_clock_and_no_position_information() {
    let (epochs, nav) = data();
    let mut gps = epochs[0].clone();
    gps.observations
        .retain(|o| o.satellite.system() == GnssSystem::Gps);
    let baseline = solve_gps_galileo_epoch(&gps, &nav, config()).unwrap();
    let mut mixed = gps.clone();
    mixed.observations.push(
        epochs[0]
            .observations
            .iter()
            .find(|o| o.satellite.system() == GnssSystem::Galileo)
            .unwrap()
            .clone(),
    );
    let s = solve_gps_galileo_epoch(&mixed, &nav, config()).unwrap();
    assert_eq!(s.clock_biases_s.len(), 2);
    assert!(
        distance(
            s.spp.position.position_ecef.unwrap(),
            baseline.spp.position.position_ecef.unwrap()
        ) < 1e-7
    );
    for i in 0..3 {
        for j in 0..3 {
            let a = s.full_covariance_m2[i][j];
            let b = baseline.full_covariance_m2[i][j];
            assert!((a - b).abs() < 1e-8 * b.abs().max(1.0));
        }
    }
}

#[test]
fn single_system_gps_is_exactly_the_existing_four_state_solver() {
    let (epochs, nav) = data();
    let mut epoch = epochs[0].clone();
    epoch
        .observations
        .retain(|o| o.satellite.system() == GnssSystem::Gps);
    assert_eq!(
        solve_gps_galileo_epoch(&epoch, &nav, config()).unwrap().spp,
        spp::solve_epoch(&epoch, &nav, config()).unwrap()
    );
}

#[test]
fn tracking_selection_retains_alternates_invalid_codes_and_raw_phase_flags() {
    let (epochs, nav) = data();
    let mut epoch = epochs[0].clone();
    let index = epoch
        .observations
        .iter()
        .position(|o| o.satellite.system() == GnssSystem::Galileo)
        .unwrap();
    let mut alternate = epoch.observations[index].clone();
    alternate.tracking_code = "1X".into();
    alternate.lli = Some(3);
    alternate.carrier_phase_cycles = None;
    epoch.observations[index].pseudorange_m = None;
    epoch.observations.push(alternate.clone());
    let s = solve_gps_galileo_epoch(&epoch, &nav, config()).unwrap();
    assert_eq!(
        s.observations[index].disposition,
        ObservationDisposition::MissingCode
    );
    assert_eq!(
        s.observations.last().unwrap().disposition,
        ObservationDisposition::Used
    );
    assert_eq!(s.observations.last().unwrap().raw_observation, alternate);
    epoch.observations[index].pseudorange_m = alternate.pseudorange_m;
    let s = solve_gps_galileo_epoch(&epoch, &nav, config()).unwrap();
    assert_eq!(
        s.observations.last().unwrap().disposition,
        ObservationDisposition::AlternateTracking
    );
    epoch.observations[0].pseudorange_m = Some(f64::INFINITY);
    let mut unsupported = epoch.observations[1].clone();
    unsupported.tracking_code = "2W".into();
    epoch.observations.push(unsupported);
    let s = solve_gps_galileo_epoch(&epoch, &nav, config()).unwrap();
    assert_eq!(
        s.observations[0].disposition,
        ObservationDisposition::InvalidCode
    );
    assert_eq!(
        s.observations.last().unwrap().disposition,
        ObservationDisposition::UnsupportedTracking
    );
}

#[test]
fn invalid_config_events_and_stale_navigation_fail_and_repeat_evaluation_is_stateless() {
    let (epochs, nav) = data();
    let epoch = &epochs[0];
    let a = solve_gps_galileo_epoch(epoch, &nav, config()).unwrap();
    assert_eq!(a, solve_gps_galileo_epoch(epoch, &nav, config()).unwrap());
    assert_eq!(
        solve_gps_galileo_epoch(
            epoch,
            &nav,
            SppConfig {
                max_iterations: 0,
                ..config()
            }
        ),
        Err(Error::InvalidConfiguration)
    );
    let mut stale = epoch.clone();
    stale.time = stale.time.checked_add_seconds(100000.0).unwrap();
    assert_eq!(
        solve_gps_galileo_epoch(&stale, &nav, config()),
        Err(Error::InsufficientObservations)
    );
    let mut event = epoch.clone();
    event.flag = 6;
    assert_eq!(
        solve_gps_galileo_epoch(&event, &nav, config()),
        Err(Error::InsufficientObservations)
    );
    let mut offset = epoch.clone();
    offset.receiver_clock_offset_s = Some(0.1);
    assert_eq!(
        solve_gps_galileo_epoch(&offset, &nav, config()),
        Err(Error::UnsupportedReceiverClockOffset)
    );
    offset.receiver_clock_offset_s = Some(f64::NAN);
    assert_eq!(
        solve_gps_galileo_epoch(&offset, &nav, config()),
        Err(Error::NonFinite)
    );
    offset.receiver_clock_offset_s = Some(0.0);
    assert_eq!(
        solve_gps_galileo_epoch(&offset, &nav, config())
            .unwrap()
            .input_receiver_clock_offset_s,
        Some(0.0)
    );
    let mut bad = epoch.clone();
    bad.approximate_position = Some(EcefCoord::new(f64::NAN, 0.0, 0.0));
    assert_eq!(
        solve_gps_galileo_epoch(&bad, &nav, config()),
        Err(Error::NonFinite)
    );
}

#[test]
fn independent_clock_columns_do_not_repair_rank_deficient_geometry() {
    let (epochs, nav) = data();
    let mut epoch = epochs[0].clone();
    let mut broken = gnss_rust::navigation::NavigationData::default();
    for original in nav.records() {
        let mut template = nav
            .records()
            .find(|e| e.satellite.system() == original.satellite.system())
            .unwrap()
            .clone();
        template.satellite = original.satellite;
        broken.add_ephemeris(template).unwrap();
    }
    let mut codes = std::collections::BTreeMap::new();
    for o in &epoch.observations {
        codes.entry(o.satellite.system()).or_insert(o.pseudorange_m);
    }
    for o in &mut epoch.observations {
        o.pseudorange_m = codes[&o.satellite.system()];
    }
    assert_eq!(
        solve_gps_galileo_epoch(&epoch, &broken, config()),
        Err(Error::SingularGeometry)
    );
}

#[test]
fn mixed_spp_seed_works_with_three_gps_and_reports_both_clocks_without_fusing_covariance() {
    use gnss_rust::rtk::{
        RtkConstellations, RtkFloatConfig, RtkFloatFilter, RtkPositionMode, RtkSeedSource,
    };
    let (mut epochs, nav) = data();
    let mut base = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_mixed_base.obs"
    )))
    .unwrap();
    let base_position = base.header().approximate_position.unwrap();
    let gps = epochs[0]
        .observations
        .iter()
        .find(|o| o.satellite.system() == GnssSystem::Gps)
        .unwrap()
        .satellite;
    for e in &mut epochs {
        e.observations.retain(|o| o.satellite != gps);
    }
    let origin = epochs[0].approximate_position.unwrap();
    let seed = EcefCoord::new(origin.x + 500.0, origin.y - 300.0, origin.z + 100.0);
    let mut filter = RtkFloatFilter::new(
        base_position,
        seed,
        RtkFloatConfig {
            constellations: RtkConstellations::GpsGalileo,
            position_mode: RtkPositionMode::Kinematic,
            kinematic_spp_seed: Some(config()),
            kinematic_spp_constellations: RtkConstellations::GpsGalileo,
            ..Default::default()
        },
    )
    .unwrap();
    for (i, (rover, base)) in epochs.iter().zip(base.by_ref()).enumerate() {
        let s = filter.process_pair(rover, &base.unwrap(), &nav).unwrap();
        let d = s.position_seed;
        assert_eq!(d.source, RtkSeedSource::SinglePoint);
        assert_eq!(
            d.spp_requested_constellations,
            Some(RtkConstellations::GpsGalileo)
        );
        assert_eq!(d.spp_satellites, 7);
        assert!(d.spp_failure.is_none());
        assert!(distance(d.position, truth(i)) < 0.003);
        assert!((d.spp_gps_clock_bias_s.unwrap() - 2e-6).abs() < 1e-11);
        assert!((d.spp_galileo_clock_bias_s.unwrap() + 3e-6).abs() < 1e-11);
        assert_eq!(s.position.status, gnss_rust::SolutionStatus::Float);
    }
}

#[test]
fn failed_mixed_spp_seed_reports_the_existing_position_fallback_without_any_clock_output() {
    use gnss_rust::rtk::{
        RtkConstellations, RtkFloatConfig, RtkFloatFilter, RtkPositionMode, RtkSeedSource,
    };
    let (epochs, nav) = data();
    let mut base = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_rtk_mixed_base.obs"
    )))
    .unwrap();
    let p = base.header().approximate_position.unwrap();
    let mut filter = RtkFloatFilter::new(
        p,
        epochs[0].approximate_position.unwrap(),
        RtkFloatConfig {
            constellations: RtkConstellations::GpsGalileo,
            position_mode: RtkPositionMode::Kinematic,
            kinematic_spp_seed: Some(SppConfig {
                max_iterations: 1,
                ..config()
            }),
            kinematic_spp_constellations: RtkConstellations::GpsGalileo,
            ..Default::default()
        },
    )
    .unwrap();
    for (i, (r, b)) in epochs.iter().zip(base.by_ref()).take(2).enumerate() {
        let s = filter.process_pair(r, &b.unwrap(), &nav).unwrap();
        let d = s.position_seed;
        assert_eq!(d.spp_failure, Some(Error::NonConvergent));
        assert_eq!(
            d.source,
            if i == 0 {
                RtkSeedSource::InitialPosition
            } else {
                RtkSeedSource::PreviousRtk
            }
        );
        assert_eq!(d.spp_satellites, 0);
        assert!(d.spp_gps_clock_bias_s.is_none() && d.spp_galileo_clock_bias_s.is_none());
    }
}
