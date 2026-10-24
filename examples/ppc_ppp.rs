//! Diagnostic GPS moving PPP on unmodified PPC observations and matching CODE products.
//! Optional explicit broadcast SPP/Doppler seeds; no reference.csv/IMU/truth input.
#[path = "support/ppp_seed.rs"]
mod ppp_seed;
#[path = "support/ppp_ztd_prior.rs"]
mod ppp_ztd_prior;
use gnss_rust::antenna::{SatelliteAntennaBinding, SatelliteAntennaModel, read_antex};
use gnss_rust::bias::{BiasSignalRequest, OsbBinding, read_satellite_biases_with_byte_comments};
use gnss_rust::dual_frequency::{
    DualFrequencyConfig, DualFrequencyNoise, PhaseBiasMode, SignalCorrections,
};
use gnss_rust::ppp::{
    PppAmbiguityInitialization, PppClockModel, PppConfig, PppEpoch, PppFilter, PppFilterSettings,
    PppMotionSeedPolicy, PppPositionModel, PppRejectedUpdatePolicy,
};
use gnss_rust::precise::{
    ClockSource, OrbitReferencePoint, PreciseProducts, read_satellite_clocks, read_sp3,
};
use gnss_rust::rinex::{RinexObservationReader, read_navigation_for_systems};
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
    let mut args: Vec<_> = std::env::args().skip(1).collect();
    let timing = args.last().is_some_and(|s| s == "-timing");
    if timing {
        args.pop();
    }
    let retain_arcs = args.last().is_some_and(|s| s == "-retain-arcs");
    if retain_arcs {
        args.pop();
    }
    let native_if_prior = args.last().is_some_and(|s| s == "-native-if-prior");
    if native_if_prior {
        args.pop();
    }
    let native_screen = args.last().is_some_and(|s| s == "-native-screen");
    if native_screen {
        args.pop();
    }
    let ztd_sigma = ppp_ztd_prior::take_sigma_flag(&mut args)?;
    let seed_profile = if args.len() >= 2
        && matches!(
            args[args.len() - 2].as_str(),
            "-seed-shadow" | "-seed-position" | "-seed-doppler"
        ) {
        let nav = args.pop().unwrap();
        Some((args.pop().unwrap(), nav))
    } else {
        None
    };
    let pva = args.last().is_some_and(|s| s == "-pva");
    if pva {
        args.pop();
    }
    let epoch_clock = args.last().is_some_and(|s| s == "-epoch-clock");
    if epoch_clock {
        args.pop();
    }
    if native_if_prior
        && !seed_profile
            .as_ref()
            .is_some_and(|(mode, _)| mode != "-seed-shadow")
    {
        return Err(
            "-native-if-prior requires an explicit -seed-position or -seed-doppler SPP source"
                .into(),
        );
    }
    if (seed_profile.is_some() && !pva)
        || !matches!(args.len(), 8 | 9)
        || !matches!(args[6].as_str(), "1C2W" | "1C2L")
        || !matches!(args[7].as_str(), "calibrated" | "raw-float")
        || (args.len() == 9
            && !matches!(
                args[8].as_str(),
                "native-variances" | "climatology" | "native-head"
            ))
    {
        return Err("usage: ppc_ppp ROVER_OBS CODE_SP3 CODE_CLK CODE_OSB_BIA GPS_IGS20_ATX CASE_LABEL 1C2W|1C2L calibrated|raw-float [native-variances|climatology|native-head] [-epoch-clock] [-pva] [-seed-shadow|-seed-position|-seed-doppler BASE_NAV] [-ztd-sigma SIGMA_M] [-native-screen] [-native-if-prior] [-retain-arcs] [-timing]".into());
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
    let mut config = PppConfig {
        initial_position,
        initial_receiver_clock_s: 0.0,
        initial_zenith_delay_m: 2.3,
        expected_coordinate_system: "IGS20".into(),
        expected_reference_point: OrbitReferencePoint::CentreOfMass,
        transmit: Default::default(),
        arcs: Default::default(),
        filter: PppFilterSettings {
            position_model: if pva {
                PppPositionModel::ConstantAcceleration(Default::default())
            } else {
                PppPositionModel::ConstantVelocity(Default::default())
            },
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
    if args.len() == 9 {
        let original = config.clone();
        config = config.with_native_precise_head_prior(first.time)?;
        match args[8].as_str() {
            "native-variances" => config.initial_zenith_delay_m = original.initial_zenith_delay_m,
            "climatology" => config.filter = original.filter,
            "native-head" => (),
            _ => unreachable!(),
        }
        eprintln!(
            "initial_head_profile={} independent_RINEX_header_seed=true thresholds_changed=false ambiguity_prior_changed=false",
            args[8]
        );
    }
    if epoch_clock {
        config.filter.clock_model = PppClockModel::IndependentEpoch;
        eprintln!(
            "clock_model=IndependentEpoch per_epoch_clock_measurement_fused=false thresholds_changed=false"
        );
    }
    if retain_arcs {
        config.filter.rejected_update_policy = PppRejectedUpdatePolicy::RetainContinuousArcs;
        eprintln!(
            "rejected_update_policy=RetainContinuousArcs rejected_measurements_accepted=false thresholds_changed=false"
        );
    }
    ppp_ztd_prior::apply(&mut config, ztd_sigma);
    if native_if_prior {
        config.filter.ambiguity_initialization = PppAmbiguityInitialization::ReferencePrecise;
    }
    eprintln!(
        "config={config:?} case={} pair={} phase={} RINEX_phase_shift_declarations={phase_shifts} excluded_receiver_clocks={} excluded_receiver_bias_rows={} byte_comment_lines={} receiver_calibration_known=false satellite_PCV_applied=false solid_earth_tide_applied=false truth_used_by_solver=false whole_native_PPP_executed=false practical_parity_accepted=false",
        args[5],
        args[6],
        args[7],
        selected_clocks.receiver_records.len(),
        selected_bias.selection.excluded_receiver_records.len(),
        selected_bias.preserved_free_text.len()
    );
    let navigation = seed_profile
        .as_ref()
        .map(|(_, path)| {
            read_navigation_for_systems(open(path)?, &[GnssSystem::Gps, GnssSystem::Galileo])
                .map_err(Box::<dyn std::error::Error>::from)
        })
        .transpose()?;
    if let Some((mode, path)) = &seed_profile {
        if mode != "-seed-shadow" {
            config.filter.motion_seed_policy = Some(PppMotionSeedPolicy {
                update_velocity: mode == "-seed-doppler",
                initialize_clock_from_seed: native_if_prior,
                ..Default::default()
            });
        }
        eprintln!(
            "seed_mode={mode} broadcast_navigation={path} same_epoch_data_dependent_prior=true position_factor_independent=false broadcast_WGS84_IGS20_compatibility_is_caller_assertion=true initial_clock_from_seed={native_if_prior} per_epoch_clock_measurement_fused=false Doppler_clock_imported=false SPP_fallback_reported_separately=true seed_policy={:?}",
            config.filter.motion_seed_policy
        );
    }
    if native_if_prior {
        eprintln!(
            "ambiguity_initialization=ReferencePrecise initial_mean=corrected_L_IF_minus_prediction variance_m2=25 initial_clock_from_SPP=true source=ppp_state.cpp:initializeAmbiguityState"
        );
    }
    if native_screen {
        config.filter.standardized_residual_gate_sigma = Some(4.0);
        eprintln!(
            "standardized_residual_screen_sigma=4 raw_sigma_screen=false whole_satellite_removed=true fresh_prior_every_retry=true existing_GNSS_gates_unchanged=true"
        );
    }
    let mut filter = PppFilter::new(config)?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "week,tow,status,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps,clock_s,zenith_m,position_covariance_trace_m2,velocity_covariance_trace_m2_s2,code_rms_m,phase_rms_m,satellites,ambiguities,arc_resets,nis,error,phase_bias_mode,raw_epoch_flag,code_corrections,phase_OSBs,correction_unavailable,last_rejections,committed_week,committed_tow,phase_unavailable{}{}{}{}{}{}",
        if timing {
            ",filter_ms,epoch_stage_ms"
        } else {
            ""
        },
        if retain_arcs {
            ",retained_arcs,retained_satellites"
        } else {
            ""
        },
        if pva {
            ",ax_mps2,ay_mps2,az_mps2,acceleration_qxx_m2_s4,acceleration_qyy_m2_s4,acceleration_qzz_m2_s4"
        } else {
            ""
        },
        if seed_profile.is_some() {
            ppp_seed::HEADER
        } else {
            ""
        },
        if native_screen {
            ",screen_attempts,screened_satellites,screen_diagnostics"
        } else {
            ""
        },
        if native_if_prior {
            ",ambiguity_seeds,ambiguity_seed_diagnostics"
        } else {
            ""
        }
    )?;
    for epoch in std::iter::once(Ok(first)).chain(reader) {
        let epoch = epoch?;
        let epoch_started = timing.then(std::time::Instant::now);
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
        let mut row =
            vec![
                String::new();
                29 + 2 * usize::from(timing) + 2 * usize::from(retain_arcs) + 6 * usize::from(pva)
            ];
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
        let seed_record = navigation
            .as_ref()
            .map(|nav| ppp_seed::Record::solve(&epoch, &nav.navigation));
        let mut supplied = seed_record.as_ref().and_then(|r| r.seed.clone());
        if seed_profile
            .as_ref()
            .is_some_and(|(m, _)| m == "-seed-position")
            && let Some(s) = &mut supplied
        {
            s.solution.velocity_ecef = None;
            s.solution.velocity_covariance = None;
        }
        let fuse = seed_profile
            .as_ref()
            .is_some_and(|(m, _)| m != "-seed-shadow");
        let filter_started = timing.then(std::time::Instant::now);
        let result = filter.process_epoch_with_seed(
            &epoch,
            &products,
            &corrections,
            &selections,
            if fuse { supplied.as_ref() } else { None },
        );
        let service_ppp = match &result {
            Ok(PppEpoch::Solution(s)) => Some(s.position.clone()),
            _ => None,
        };
        if let Some(started) = filter_started {
            row[29] = format!("{:.17e}", started.elapsed().as_secs_f64() * 1000.0);
        }
        if let Some(started) = epoch_started {
            row[30] = format!("{:.17e}", started.elapsed().as_secs_f64() * 1000.0);
        }
        if retain_arcs {
            let index = 29 + 2 * usize::from(timing);
            row[index] = filter.last_retained_phase_arcs().len().to_string();
            row[index + 1] = filter
                .last_retained_phase_arcs()
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(";");
        }
        match result {
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
                if pva {
                    let index = 29 + 2 * usize::from(timing) + 2 * usize::from(retain_arcs);
                    let acceleration = s.acceleration_ecef_mps2.unwrap();
                    let covariance = s.acceleration_covariance_m2_s4.unwrap();
                    for i in 0..3 {
                        row[index + i] = format!("{:.17e}", acceleration[i]);
                        row[index + 3 + i] = format!("{:.17e}", covariance[i][i]);
                    }
                }
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
        if let Some(r) = seed_record {
            row.extend(r.append(filter.last_motion_seed(), service_ppp.as_ref(), epoch.time));
        }
        if native_screen {
            let d = filter.last_residual_screening();
            row.extend([
                d.len().to_string(),
                d.iter()
                    .filter(|d| d.removed)
                    .map(|d| d.satellite.to_string())
                    .collect::<Vec<_>>()
                    .join(";"),
                format!("{d:?}"),
            ]);
        }
        if native_if_prior {
            let seeds = filter.last_ambiguity_seeds();
            row.extend([seeds.len().to_string(), format!("{seeds:?}")]);
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
