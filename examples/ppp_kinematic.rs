//! Calibrated synthetic moving IF PPP; no real-data performance claim.
use gnss_rust::SatelliteId;
use gnss_rust::bias::{BiasSignalRequest, OsbBinding, read_bias_sinex};
use gnss_rust::dual_frequency::{DualFrequencyConfig, DualFrequencyNoise};
use gnss_rust::ppp::*;
use gnss_rust::precise::{ClockSource, OrbitReferencePoint, PreciseProducts, read_clk, read_sp3};
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::tides::ReceiverPositionReference;
use std::{collections::BTreeMap, io::Cursor};

fn products() -> PreciseProducts {
    PreciseProducts::new(
        read_sp3(Cursor::new(include_str!(
            "../tests/fixtures/synthetic_ppp.sp3"
        )))
        .unwrap(),
        Some(
            read_clk(Cursor::new(include_str!(
                "../tests/fixtures/synthetic_ppp.clk"
            )))
            .unwrap(),
        ),
    )
}
fn epochs(noisy: bool) -> Vec<gnss_rust::observation::ObservationEpoch> {
    RinexObservationReader::new(Cursor::new(if noisy {
        include_str!("../tests/fixtures/synthetic_ppp_kinematic_noisy.obs")
    } else {
        include_str!("../tests/fixtures/synthetic_ppp_kinematic.obs")
    }))
    .unwrap()
    .collect::<Result<_, _>>()
    .unwrap()
}
fn input_epochs(kind: &str) -> Vec<gnss_rust::observation::ObservationEpoch> {
    let mut epochs = epochs(kind == "noisy");
    let phase = matches!(kind, "unrounded-phase" | "unrounded-both");
    let code = matches!(kind, "unrounded-code" | "unrounded-both");
    if phase || code {
        for line in include_str!("../tests/fixtures/upstream_ppp_kinematic_unrounded.csv")
            .lines()
            .skip(1)
        {
            let v: Vec<_> = line.split(',').collect();
            let i: usize = v[0].parse().unwrap();
            let sat: SatelliteId = v[1].parse().unwrap();
            for (j, tracking) in ["1C", "2W"].iter().enumerate() {
                let o = epochs[i]
                    .observations
                    .iter_mut()
                    .find(|o| o.satellite == sat && o.tracking_code == *tracking)
                    .unwrap();
                if code {
                    o.pseudorange_m = Some(v[2 + 2 * j].parse().unwrap());
                }
                if phase {
                    o.carrier_phase_cycles = Some(v[3 + 2 * j].parse().unwrap());
                }
            }
        }
    }
    epochs
}
fn selections() -> BTreeMap<SatelliteId, DualFrequencyConfig> {
    [1, 2, 7, 8, 9, 18, 19, 25, 26]
        .into_iter()
        .map(|n| {
            (
                format!("G{n:02}").parse().unwrap(),
                DualFrequencyConfig {
                    tracking_codes: ["1C".into(), "2W".into()],
                    clock_source: ClockSource::RinexClk,
                    clock_reference: "synthetic-precise-clock-v1".into(),
                    phase_reference: "synthetic-phase-v1".into(),
                    noise: DualFrequencyNoise {
                        code_m2: [[0.0025, 0.0], [0.0, 0.0025]],
                        phase_m2: [[0.000001, 0.0], [0.0, 0.000001]],
                    },
                },
            )
        })
        .collect()
}
fn biases(t: gnss_rust::GnssTime) -> gnss_rust::dual_frequency::SignalCorrections {
    let p = read_bias_sinex(
        Cursor::new(include_str!("../tests/fixtures/synthetic_ppp_osb.bsx")),
        Default::default(),
    )
    .unwrap();
    let requests: Vec<_> = selections()
        .into_iter()
        .flat_map(|(sat, c)| {
            c.tracking_codes
                .into_iter()
                .map(move |code| BiasSignalRequest {
                    satellite: sat,
                    tracking_code: code,
                    glonass_channel: None,
                })
        })
        .collect();
    p.corrections_at(
        t,
        &OsbBinding {
            product_id: "synthetic_ppp_osb.bsx".into(),
            clock_source: ClockSource::RinexClk,
            clock_reference: "synthetic-precise-clock-v1".into(),
            phase_reference: Some("synthetic-phase-v1".into()),
        },
        &requests,
    )
    .unwrap()
    .corrections
}
fn config() -> PppConfig {
    PppConfig {
        initial_position: epochs(false)[0].approximate_position.unwrap(),
        initial_receiver_clock_s: 0.0,
        initial_zenith_delay_m: 2.3,
        expected_coordinate_system: "IGS20".into(),
        expected_reference_point: OrbitReferencePoint::AntennaPhaseCenter,
        transmit: Default::default(),
        arcs: Default::default(),
        filter: PppFilterSettings {
            position_model: PppPositionModel::ConstantVelocity(Default::default()),
            ..Default::default()
        },
        phase_windup: None,
        receiver_antenna: None,
        satellite_antenna: None,
        receiver_reference: ReceiverPositionReference::InstantaneousMarker,
        solid_earth_tide: None,
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    use std::io::Write;
    let args: Vec<_> = std::env::args().skip(1).collect();
    if !matches!(args.len(), 2 | 3)
        || !matches!(args[0].as_str(), "static" | "kinematic")
        || !matches!(
            args[1].as_str(),
            "clean" | "noisy" | "unrounded-phase" | "unrounded-code" | "unrounded-both"
        )
    {
        return Err(
            "usage: ppp_kinematic static|kinematic clean|noisy [INITIAL_VELOCITY_VARIANCE_M2_S2]"
                .into(),
        );
    }
    let mut c = config();
    if args[0] == "static" {
        c.filter.position_model = PppPositionModel::Static;
    }
    if args.len() == 3 {
        if args[0] != "kinematic" {
            return Err("velocity variance requires the kinematic model".into());
        }
        c.filter.position_model = PppPositionModel::ConstantVelocity(PppKinematicConfig {
            initial_velocity_variance_m2_s2: args[2].parse()?,
            ..Default::default()
        });
    }
    let selected_model = c.filter.position_model;
    let mut filter = PppFilter::new(c)?;
    let products = products();
    let selections = selections();
    let mut out = std::io::BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "week,tow,status,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps,clock_s,zenith_m,position_covariance_trace_m2,velocity_covariance_trace_m2_s2,code_rms_m,phase_rms_m,satellites,ambiguities,arc_resets,nis,error"
    )?;
    for epoch in input_epochs(&args[1]) {
        let mut row = vec![String::new(); 20];
        row[0] = epoch.time.week().to_string();
        row[1] = format!("{:.17}", epoch.time.tow());
        match filter.process_epoch(&epoch, &products, &biases(epoch.time), &selections) {
            Ok(PppEpoch::Solution(s)) => {
                row[2] = "PPP_FLOAT".into();
                let p = s.position.position_ecef.unwrap();
                for (i, v) in [p.x, p.y, p.z].into_iter().enumerate() {
                    row[3 + i] = format!("{v:.17e}");
                }
                if let Some(v) = s.position.velocity_ecef {
                    for i in 0..3 {
                        row[6 + i] = format!("{:.17e}", v[i]);
                    }
                }
                row[9] = format!("{:.17e}", s.position.receiver_clock_bias);
                row[10] = format!("{:.17e}", s.zenith_delay_m);
                let p = s.position.position_covariance.unwrap();
                row[11] = format!("{:.17e}", p[0][0] + p[1][1] + p[2][2]);
                if let Some(p) = s.position.velocity_covariance {
                    row[12] = format!("{:.17e}", p[0][0] + p[1][1] + p[2][2]);
                }
                row[13] = format!("{:.17e}", s.code_rms_m);
                row[14] = format!("{:.17e}", s.phase_rms_m);
                row[15] = s.measurements.len().to_string();
                row[16] = s.ambiguities.len().to_string();
                row[17] = s
                    .arc_decisions
                    .values()
                    .filter(|d| !d.reset_reasons.is_empty())
                    .count()
                    .to_string();
                row[18] = format!("{:.17e}", s.normalized_innovation_squared);
            }
            Ok(PppEpoch::Event(_)) => row[2] = "EVENT".into(),
            Err(e) => {
                row[2] = "INVALID".into();
                row[19] = format!("\"{}\"", e.to_string().replace('"', "\"\""));
            }
        }
        writeln!(out, "{}", row.join(","))?;
    }
    out.flush()?;
    eprintln!(
        "synthetic_moving_IF_PPP=true default_position_model=Static selected_model={selected_model:?} injected_noise={} truth_used_by_solver=false whole_native_PPP_executed=false real_PPC_performance_accepted=false",
        args[1]
    );
    Ok(())
}
