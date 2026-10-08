//! Initial static GPS IF PPP demonstration on calibrated synthetic products.
//! `noisy` selects deterministic noise; `windup` selects native phase injection
//! and enables the matching nominal-yaw correction. These are synthetic data.
//! `antenna` adds a fixed receiver ANTEX model and also enables phase wind-up.
//! `satellite` selects COM products plus satellite/receiver PCO and wind-up.
//! `legacy-tide` adds the explicit legacy Love-number solid Earth tide.
use gnss_rust::SatelliteId;
use gnss_rust::antenna::{
    ReceiverAntennaModel, SatelliteAntennaBinding, SatelliteAntennaModel, read_antex,
};
use gnss_rust::bias::{BiasKind, BiasSignalRequest, OsbBinding, read_bias_sinex};
use gnss_rust::dual_frequency::{DualFrequencyConfig, DualFrequencyNoise};
use gnss_rust::ppp::{PppConfig, PppEpoch, PppFilter};
use gnss_rust::precise::{ClockSource, OrbitReferencePoint, PreciseProducts, read_clk, read_sp3};
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::tides::{ReceiverPositionReference, SolidEarthTideModel};
use gnss_rust::windup::PppWindupModel;
use std::{collections::BTreeMap, io::Cursor};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let arguments: Vec<_> = std::env::args().skip(1).collect();
    let mut noisy = false;
    let mut windup = false;
    let mut antenna = false;
    let mut satellite = false;
    let mut tide = false;
    for arg in &arguments {
        match arg.as_str() {
            "noisy" if !noisy => noisy = true,
            "windup" if !windup => windup = true,
            "antenna" if !antenna => antenna = true,
            "satellite" if !satellite => satellite = true,
            "legacy-tide" if !tide => tide = true,
            _ => {
                return Err(
                    "usage: ppp [noisy] [windup] [antenna] [satellite] [legacy-tide]".into(),
                );
            }
        }
    }
    // Antenna fixtures also carry native phase wind-up injection.
    satellite |= tide;
    antenna |= satellite;
    windup |= antenna;
    let products = PreciseProducts::new(
        read_sp3(Cursor::new(if satellite {
            include_str!("../tests/fixtures/synthetic_ppp_com.sp3")
        } else {
            include_str!("../tests/fixtures/synthetic_ppp.sp3")
        }))?,
        Some(read_clk(Cursor::new(include_str!(
            "../tests/fixtures/synthetic_ppp.clk"
        )))?),
    );
    let osbs = read_bias_sinex(
        Cursor::new(include_str!("../tests/fixtures/synthetic_ppp_osb.bsx")),
        Default::default(),
    )?;
    // This binding is known from the fixture generator. Real products need
    // independently established compatible clock/phase datum declarations.
    let binding = OsbBinding {
        product_id: "synthetic_ppp_osb.bsx".into(),
        clock_source: ClockSource::RinexClk,
        clock_reference: "synthetic-precise-clock-v1".into(),
        phase_reference: Some("synthetic-phase-v1".into()),
    };
    let selections: BTreeMap<SatelliteId, DualFrequencyConfig> = osbs
        .entries()
        .iter()
        .filter(|r| {
            r.kind == BiasKind::Osb && r.station.is_none() && r.observation_1.starts_with('C')
        })
        .map(|row| {
            (
                row.satellite,
                DualFrequencyConfig {
                    tracking_codes: ["1C".to_owned(), "2W".to_owned()],
                    clock_source: ClockSource::RinexClk,
                    clock_reference: binding.clock_reference.clone(),
                    phase_reference: binding.phase_reference.clone().unwrap(),
                    noise: DualFrequencyNoise {
                        code_m2: [[0.0025, 0.0], [0.0, 0.0025]],
                        phase_m2: [[0.000001, 0.0], [0.0, 0.000001]],
                    },
                },
            )
        })
        .collect();
    let requests: Vec<_> = selections
        .iter()
        .flat_map(|(&satellite, selection)| {
            selection
                .tracking_codes
                .iter()
                .map(move |tracking_code| BiasSignalRequest {
                    satellite,
                    tracking_code: tracking_code.clone(),
                    glonass_channel: None,
                })
        })
        .collect();
    let reader = RinexObservationReader::new(Cursor::new(
        match (noisy, windup, antenna, satellite, tide) {
            (false, _, _, _, true) => {
                include_str!("../tests/fixtures/synthetic_ppp_legacy_tide.obs")
            }
            (true, _, _, _, true) => {
                include_str!("../tests/fixtures/synthetic_ppp_legacy_tide_noisy.obs")
            }
            (false, _, _, true, false) => {
                include_str!("../tests/fixtures/synthetic_ppp_satellite_antenna.obs")
            }
            (true, _, _, true, false) => {
                include_str!("../tests/fixtures/synthetic_ppp_satellite_antenna_noisy.obs")
            }
            (false, _, true, false, false) => {
                include_str!("../tests/fixtures/synthetic_ppp_antenna.obs")
            }
            (true, _, true, false, false) => {
                include_str!("../tests/fixtures/synthetic_ppp_antenna_noisy.obs")
            }
            (false, false, false, false, false) => {
                include_str!("../tests/fixtures/synthetic_ppp.obs")
            }
            (true, false, false, false, false) => {
                include_str!("../tests/fixtures/synthetic_ppp_noisy.obs")
            }
            (false, true, false, false, false) => {
                include_str!("../tests/fixtures/synthetic_ppp_windup.obs")
            }
            (true, true, false, false, false) => {
                include_str!("../tests/fixtures/synthetic_ppp_windup_noisy.obs")
            }
        },
    ))?;
    let receiver_antenna = if antenna {
        let atx = read_antex(Cursor::new(include_str!(
            "../tests/fixtures/synthetic_receiver.atx"
        )))?;
        Some(ReceiverAntennaModel::new(
            "synthetic_receiver.atx".into(),
            atx.receiver(
                "TESTANT NONE",
                "RX001",
                gnss_rust::GnssTime::new(2300, 346200.0)?,
            )?,
            gnss_rust::EnuCoord {
                east: 0.025,
                north: -0.035,
                up: 0.8,
            },
        )?)
    } else {
        None
    };
    let satellite_antenna = if satellite {
        let atx = read_antex(Cursor::new(include_str!(
            "../tests/fixtures/synthetic_satellite.atx"
        )))?;
        let entries = atx
            .entries()
            .iter()
            .map(|c| {
                atx.satellite_calibration(
                    c.satellite().ok_or("satellite PRN missing")?,
                    &c.svn,
                    gnss_rust::GnssTime::new(2300, 346200.0)?,
                )
                .cloned()
                .map_err(Box::<dyn std::error::Error>::from)
            })
            .collect::<Result<Vec<_>, Box<dyn std::error::Error>>>()?;
        Some(SatelliteAntennaModel::new(
            SatelliteAntennaBinding {
                product_id: "synthetic_satellite.atx".into(),
                clock_source: binding.clock_source,
                clock_reference: binding.clock_reference.clone(),
                tracking_codes: ["1C".into(), "2W".into()],
            },
            entries,
        )?)
    } else {
        None
    };
    let config = PppConfig {
        initial_position: reader
            .header()
            .approximate_position
            .ok_or("fixture seed missing")?,
        initial_receiver_clock_s: 0.0,
        initial_zenith_delay_m: 2.3,
        expected_coordinate_system: "IGS20".to_owned(),
        expected_reference_point: if satellite {
            OrbitReferencePoint::CentreOfMass
        } else {
            OrbitReferencePoint::AntennaPhaseCenter
        },
        transmit: Default::default(),
        arcs: Default::default(),
        filter: Default::default(),
        phase_windup: windup.then_some(PppWindupModel::NominalYawApproximateSun),
        receiver_antenna,
        satellite_antenna,
        receiver_reference: if tide {
            ReceiverPositionReference::NominalMarker
        } else {
            ReceiverPositionReference::InstantaneousMarker
        },
        solid_earth_tide: tide.then_some(SolidEarthTideModel::LegacyLoveApproximateSunMoon),
    };
    let mut filter = PppFilter::new(config)?;
    println!(
        "week,tow,x_m,y_m,z_m,receiver_clock_s,satellites,status,zenith_delay_m,code_rms_m,phase_rms_m,iterations,phase_satellites,nis,reset_arcs,windup_satellites,windup_max_abs_cycles,receiver_antenna_satellites,receiver_antenna_max_abs_if_m,satellite_antenna_satellites,satellite_antenna_max_offset_m,solid_tide_satellites,solid_tide_displacement_m"
    );
    let mut count = 0;
    for epoch in reader {
        let epoch = epoch?;
        let snapshot = match osbs.corrections_at(epoch.time, &binding, &requests) {
            Ok(snapshot) => snapshot,
            Err(error) => {
                filter.invalidate_phase()?;
                return Err(error.into());
            }
        };
        match filter.process_epoch(&epoch, &products, &snapshot.corrections, &selections)? {
            PppEpoch::Solution(result) => {
                let p = result
                    .position
                    .position_ecef
                    .ok_or("PPP position missing")?;
                let resets = result
                    .arc_commit
                    .accepted
                    .keys()
                    .filter(|sat| !result.arc_decisions[sat].reset_reasons.is_empty())
                    .count();
                let cycles: Vec<_> = result
                    .measurements
                    .iter()
                    .filter_map(|m| m.windup.as_ref().map(|w| w.angles.cycles))
                    .collect();
                let antenna: Vec<_> = result
                    .measurements
                    .iter()
                    .filter_map(|m| m.receiver_antenna.as_ref().map(|a| a.if_add_m))
                    .collect();
                let satellite: Vec<_> = result
                    .measurements
                    .iter()
                    .filter_map(|m| {
                        m.model.satellite_antenna.as_ref().map(|a| {
                            a.offset_ecef_m
                                .x
                                .hypot(a.offset_ecef_m.y)
                                .hypot(a.offset_ecef_m.z)
                        })
                    })
                    .collect();
                let tides: Vec<_> = result
                    .measurements
                    .iter()
                    .filter_map(|m| {
                        m.model.solid_earth_tide.as_ref().map(|t| {
                            t.displacement_m
                                .x
                                .hypot(t.displacement_m.y)
                                .hypot(t.displacement_m.z)
                        })
                    })
                    .collect();
                println!(
                    "{},{:.7},{:.9},{:.9},{:.9},{:.15e},{},PPP_FLOAT,{:.9},{:.9},{:.9},{},{},{:.9},{},{},{:.9},{},{:.9},{},{:.9},{},{:.9}",
                    epoch.time.week(),
                    epoch.time.tow(),
                    p.x,
                    p.y,
                    p.z,
                    result.position.receiver_clock_bias,
                    result.position.num_satellites,
                    result.zenith_delay_m,
                    result.code_rms_m,
                    result.phase_rms_m,
                    result.iterations,
                    result.ambiguities.len(),
                    result.normalized_innovation_squared,
                    resets,
                    cycles.len(),
                    cycles.iter().map(|x| x.abs()).fold(0.0, f64::max),
                    antenna.len(),
                    antenna.iter().map(|x| x.abs()).fold(0.0, f64::max),
                    satellite.len(),
                    satellite.iter().copied().fold(0.0, f64::max),
                    tides.len(),
                    tides.iter().copied().fold(0.0, f64::max)
                );
                count += 1;
            }
            PppEpoch::Event(_) => eprintln!(
                "phase event at GPS {} {}",
                epoch.time.week(),
                epoch.time.tow()
            ),
        }
    }
    eprintln!(
        "completed {count} initial PPP FLOAT epochs (synthetic, noisy={noisy}, windup={windup})"
    );
    Ok(())
}
