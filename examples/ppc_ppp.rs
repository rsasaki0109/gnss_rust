//! Diagnostic GPS moving PPP on unmodified PPC observations and matching CODE products.
//! No reference.csv, base observations, SPP, Doppler or truth enters the filter.
use gnss_rust::antenna::{SatelliteAntennaBinding, SatelliteAntennaModel, read_antex};
use gnss_rust::bias::{BiasSignalRequest, OsbBinding, read_satellite_biases_with_byte_comments};
use gnss_rust::dual_frequency::{
    DualFrequencyConfig, DualFrequencyNoise, PhaseBiasMode, SignalCorrections,
};
use gnss_rust::ppp::{PppConfig, PppEpoch, PppFilter, PppFilterSettings, PppPositionModel};
use gnss_rust::precise::{
    ClockSource, OrbitReferencePoint, PreciseProducts, read_satellite_clocks, read_sp3,
};
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::tides::ReceiverPositionReference;
use gnss_rust::windup::PppWindupModel;
use gnss_rust::{GnssSystem, SatelliteId};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufReader, BufWriter, Write},
};

fn csv(value: &str) -> String {
    if value.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.into()
    }
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 8
        || !matches!(args[6].as_str(), "1C2W" | "1C2L")
        || !matches!(args[7].as_str(), "calibrated" | "raw-float")
    {
        return Err("usage: ppc_ppp ROVER_OBS CODE_SP3 CODE_CLK CODE_OSB_BIA GPS_IGS20_ATX CASE_LABEL 1C2W|1C2L calibrated|raw-float".into());
    }
    let open = |name: &str| File::open(name).map(BufReader::new);
    let orbit = read_sp3(open(&args[1])?)?;
    let selected_clocks = read_satellite_clocks(open(&args[2])?)?;
    let selected_bias =
        read_satellite_biases_with_byte_comments(open(&args[3])?, Default::default())?;
    let bias = &selected_bias.selection.products;
    if orbit.header().coordinate_system != "IGS20"
        || !bias
            .description()
            .get("SATELLITE_CLOCK_REFERENCE_OBSERVABLES")
            .is_some_and(|r| r.iter().any(|v| v.as_slice() == ["G", "C1W", "C2W"]))
        || !bias
            .description()
            .get("APC_MODEL")
            .is_some_and(|r| r.len() == 1 && r[0].as_slice() == ["IGS20"])
        || !selected_clocks
            .clocks
            .header()
            .raw_lines
            .iter()
            .any(|l| l.starts_with('G') && l.contains("IGS20") && l.contains("SYS / PCVS APPLIED"))
    {
        return Err("required CODE GPS IF(C1W,C2W)/IGS20 APC declarations unavailable".into());
    }
    let products = PreciseProducts::new(orbit, Some(selected_clocks.clocks));
    let mut reader = RinexObservationReader::new(open(&args[0])?)?;
    let initial_position = reader
        .header()
        .approximate_position
        .ok_or("RINEX approximate position missing")?;
    let phase_shifts = reader.header().phase_shifts.len();
    let first = reader.next().transpose()?.ok_or("no rover epochs")?;
    let atx = read_antex(open(&args[4])?)?;
    let codes: [String; 2] = ["1C".into(), args[6][2..].into()];
    let clock_reference = format!(
        "COD0OPSRAP:{}:GPS-IF(C1W,C2W)",
        std::path::Path::new(&args[3])
            .file_name()
            .ok_or("bias filename")?
            .to_string_lossy()
    );
    let raw = args[7] == "raw-float";
    let phase_reference = if raw {
        "raw-carrier-FLOAT-gauge-v1".to_owned()
    } else {
        format!("CODE-absolute-phase:{}", clock_reference)
    };
    let binding = OsbBinding {
        product_id: args[3].clone(),
        clock_source: ClockSource::RinexClk,
        clock_reference: clock_reference.clone(),
        phase_reference: (!raw).then(|| phase_reference.clone()),
    };
    let mut selections = BTreeMap::new();
    let mut calibrations = Vec::new();
    for prn in 1..=32 {
        let sat: SatelliteId = format!("G{prn:02}").parse()?;
        let found: Vec<_> = atx
            .entries()
            .iter()
            .filter(|c| c.satellite() == Some(sat) && c.contains(first.time))
            .collect();
        if found.len() != 1 {
            return Err(format!("nonunique ANTEX at start: {sat}").into());
        }
        let c = atx.satellite_calibration(sat, &found[0].svn, first.time)?;
        if bias
            .entries()
            .iter()
            .filter(|e| e.satellite == sat && e.station.is_none() && e.contains(first.time))
            .any(|e| e.svn.as_deref() != Some(c.svn.as_str()))
        {
            return Err(format!("CODE/ANTEX SVN mismatch: {sat}").into());
        }
        calibrations.push(c.clone());
        selections.insert(
            sat,
            DualFrequencyConfig {
                tracking_codes: codes.clone(),
                clock_source: ClockSource::RinexClk,
                clock_reference: clock_reference.clone(),
                phase_reference: phase_reference.clone(),
                noise: DualFrequencyNoise {
                    code_m2: [[1.0, 0.0], [0.0, 1.0]],
                    phase_m2: [[1e-4, 0.0], [0.0, 1e-4]],
                },
            },
        );
        eprintln!(
            "calibration satellite={sat} svn={} source_line={} first_week={} first_tow={:.7}",
            c.svn,
            c.source_line,
            first.time.week(),
            first.time.tow()
        );
    }
    let antenna = SatelliteAntennaModel::new(
        SatelliteAntennaBinding {
            product_id: args[4].clone(),
            clock_source: ClockSource::RinexClk,
            clock_reference,
            tracking_codes: codes.clone(),
        },
        calibrations,
    )?;
    let config = PppConfig {
        initial_position,
        initial_receiver_clock_s: 0.0,
        initial_zenith_delay_m: 2.3,
        expected_coordinate_system: "IGS20".into(),
        expected_reference_point: OrbitReferencePoint::CentreOfMass,
        transmit: Default::default(),
        arcs: Default::default(),
        filter: PppFilterSettings {
            position_model: PppPositionModel::ConstantVelocity(Default::default()),
            phase_bias_mode: if raw {
                PhaseBiasMode::UncalibratedFloat
            } else {
                PhaseBiasMode::RequireSignalCorrections
            },
            ..Default::default()
        },
        phase_windup: Some(PppWindupModel::NominalYawApproximateSun),
        receiver_antenna: None,
        satellite_antenna: Some(antenna),
        receiver_reference: ReceiverPositionReference::InstantaneousMarker,
        solid_earth_tide: None,
    };
    eprintln!(
        "config={config:?} case={} pair={} phase={} RINEX_phase_shift_declarations={phase_shifts} excluded_receiver_clocks={} excluded_receiver_bias_rows={} byte_comment_lines={} receiver_calibration_known=false satellite_PCV_applied=false solid_earth_tide_applied=false truth_used_by_solver=false whole_native_PPP_executed=false practical_parity_accepted=false",
        args[5],
        args[6],
        args[7],
        selected_clocks.receiver_records.len(),
        selected_bias.selection.excluded_receiver_records.len(),
        selected_bias.preserved_free_text.len()
    );
    let mut filter = PppFilter::new(config)?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "week,tow,status,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps,clock_s,zenith_m,position_covariance_trace_m2,velocity_covariance_trace_m2_s2,code_rms_m,phase_rms_m,satellites,ambiguities,arc_resets,nis,error,phase_bias_mode,raw_epoch_flag,code_corrections,phase_OSBs,correction_unavailable,last_rejections,committed_week,committed_tow,phase_unavailable"
    )?;
    for epoch in std::iter::once(Ok(first)).chain(reader) {
        let epoch = epoch?;
        let mut corrections = Vec::new();
        let mut unavailable = Vec::new();
        for observation in epoch.observations.iter().filter(|o| {
            o.satellite.system() == GnssSystem::Gps
                && o.satellite.prn() <= 32
                && codes.contains(&o.tracking_code)
        }) {
            let request = BiasSignalRequest {
                satellite: observation.satellite,
                tracking_code: observation.tracking_code.clone(),
                glonass_channel: None,
            };
            match bias.corrections_at(epoch.time, &binding, &[request]) {
                Ok(snapshot) => {
                    corrections.extend(snapshot.corrections.entries().values().flatten().cloned())
                }
                Err(e) => unavailable.push(format!(
                    "{}:{}:{e:?}",
                    observation.satellite, observation.tracking_code
                )),
            }
        }
        let code_count = corrections.len();
        let phase_count = corrections
            .iter()
            .filter(|c| c.phase_add_m.is_some())
            .count();
        let corrections = SignalCorrections::new(corrections)?;
        let mut row = vec![String::new(); 29];
        row[0] = epoch.time.week().to_string();
        row[1] = format!("{:.17}", epoch.time.tow());
        row[20] = format!(
            "{:?}",
            if raw {
                PhaseBiasMode::UncalibratedFloat
            } else {
                PhaseBiasMode::RequireSignalCorrections
            }
        );
        row[21] = epoch.flag.to_string();
        row[22] = code_count.to_string();
        row[23] = phase_count.to_string();
        row[24] = unavailable.join(";");
        match filter.process_epoch(&epoch, &products, &corrections, &selections) {
            Ok(PppEpoch::Solution(s)) => {
                row[2] = "PPP_FLOAT".into();
                let p = s.position.position_ecef.unwrap();
                for (i, v) in [p.x, p.y, p.z].into_iter().enumerate() {
                    row[3 + i] = format!("{v:.17e}");
                }
                for (i, v) in s.position.velocity_ecef.unwrap().into_iter().enumerate() {
                    row[6 + i] = format!("{v:.17e}");
                }
                row[9] = format!("{:.17e}", s.position.receiver_clock_bias);
                row[10] = format!("{:.17e}", s.zenith_delay_m);
                let p = s.position.position_covariance.unwrap();
                row[11] = format!("{:.17e}", p[0][0] + p[1][1] + p[2][2]);
                let p = s.position.velocity_covariance.unwrap();
                row[12] = format!("{:.17e}", p[0][0] + p[1][1] + p[2][2]);
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
                row[28] = s
                    .arc_decisions
                    .iter()
                    .filter_map(|(sat, d)| {
                        d.measurement
                            .as_ref()
                            .and_then(|m| m.phase_unavailability.as_ref())
                            .map(|e| format!("{sat}:{e:?}"))
                    })
                    .collect::<Vec<_>>()
                    .join(";");
            }
            Ok(PppEpoch::Event(_)) => row[2] = "EVENT".into(),
            Err(e) => {
                row[2] = "INVALID".into();
                row[19] = e.to_string();
            }
        }
        row[25] = filter
            .last_rejections()
            .iter()
            .map(|(s, r)| format!("{s}:{r:?}"))
            .collect::<Vec<_>>()
            .join(";");
        if let Some(t) = filter.state_time() {
            row[26] = t.week().to_string();
            row[27] = format!("{:.17}", t.tow());
        }
        writeln!(
            out,
            "{}",
            row.iter().map(|s| csv(s)).collect::<Vec<_>>().join(",")
        )?;
    }
    out.flush()?;
    Ok(())
}
