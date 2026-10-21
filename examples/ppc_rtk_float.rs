//! L1/E1 FLOAT on common epochs; every unmatched/invalid epoch is reported.
use gnss_rust::epoch_match::{EpochMatch, EpochMatcher, EpochSide};
use gnss_rust::rinex::{RinexObservationReader, read_navigation_for_systems};
use gnss_rust::rtk::{
    RtkConstellations, RtkFloatConfig, RtkFloatFilter, RtkPositionMode, RtkPrefitDiagnostic,
    RtkSeedDiagnostic, RtkSlipDiagnostic,
};
use gnss_rust::rtk_slip::RtkSlipConfig;
use gnss_rust::spp::SppConfig;
use gnss_rust::{EcefCoord, GnssSystem};
use std::{
    error::Error,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

fn write_seed(row: &mut [String], offset: usize, seed: RtkSeedDiagnostic) {
    row[offset] = format!("{:?}", seed.source);
    for (i, v) in [seed.position.x, seed.position.y, seed.position.z]
        .into_iter()
        .enumerate()
    {
        row[offset + 1 + i] = format!("{v:.17e}");
    }
    if let Some(error) = seed.spp_failure {
        row[offset + 4] = format!("{error:?}");
    }
}
fn write_spp_scope(row: &mut [String], seed: RtkSeedDiagnostic) {
    let i = row.len() - 4;
    row[i] = seed
        .spp_requested_constellations
        .map_or_else(String::new, |c| format!("{c:?}"));
    row[i + 1] = seed.spp_satellites.to_string();
    row[i + 2] = seed
        .spp_gps_clock_bias_s
        .map_or_else(String::new, |v| format!("{v:.17e}"));
    row[i + 3] = seed
        .spp_galileo_clock_bias_s
        .map_or_else(String::new, |v| format!("{v:.17e}"));
}

fn write_slips(row: &mut [String], offset: usize, evidence: &[RtkSlipDiagnostic]) {
    for (i, value) in [
        evidence
            .iter()
            .filter(|d| d.code.is_some_and(|t| t.detected))
            .count(),
        evidence
            .iter()
            .filter(|d| d.doppler.is_some_and(|t| t.detected))
            .count(),
        evidence.iter().filter(|d| d.code.is_some()).count(),
        evidence.iter().filter(|d| d.doppler.is_some()).count(),
    ]
    .into_iter()
    .enumerate()
    {
        row[offset + i] = value.to_string();
    }
    row[offset + 4] = evidence
        .iter()
        .filter(|d| d.code.is_some_and(|t| t.detected) || d.doppler.is_some_and(|t| t.detected))
        .map(|d| d.satellite.to_string())
        .collect::<Vec<_>>()
        .join(";");
    row[offset + 5] = evidence
        .iter()
        .filter_map(|d| d.doppler_error.map(|e| format!("{}:{e:?}", d.satellite)))
        .collect::<Vec<_>>()
        .join(";");
}

fn write_prefit(row: &mut [String], diagnostic: &RtkPrefitDiagnostic) {
    let screen = &diagnostic.screen;
    row[26] = format!("{:.17e}", screen.maximum_phase_abs_m);
    row[27] = format!("{:.17e}", screen.maximum_code_abs_m);
    row[28] = format!("{:.17e}", screen.threshold_m);
    row[29] = screen.rejected_phase_rows.len().to_string();
    row[30] = screen.rejected_code_rows.len().to_string();
    let mut rejected = Vec::new();
    for (kind, indices) in [
        ("phase", &screen.rejected_phase_rows),
        ("code", &screen.rejected_code_rows),
    ] {
        rejected.extend(indices.iter().map(|i| {
            format!(
                "{kind}:{}-{}",
                diagnostic.row_reference_satellites[*i], diagnostic.row_satellites[*i]
            )
        }));
    }
    row[31] = rejected.join(";");
    let seed = diagnostic.position_seed;
    row[32] = format!("{:?}", seed.source);
    for (i, v) in [seed.position.x, seed.position.y, seed.position.z]
        .into_iter()
        .enumerate()
    {
        row[33 + i] = format!("{v:.17e}");
    }
    if let Some(error) = seed.spp_failure {
        row[36] = format!("{error:?}");
    }
    let evidence = &diagnostic.slip_diagnostics;
    for (i, value) in [
        evidence
            .iter()
            .filter(|d| d.code.is_some_and(|t| t.detected))
            .count(),
        evidence
            .iter()
            .filter(|d| d.doppler.is_some_and(|t| t.detected))
            .count(),
        evidence.iter().filter(|d| d.code.is_some()).count(),
        evidence.iter().filter(|d| d.doppler.is_some()).count(),
    ]
    .into_iter()
    .enumerate()
    {
        row[37 + i] = value.to_string();
    }
    row[41] = evidence
        .iter()
        .filter(|d| d.code.is_some_and(|t| t.detected) || d.doppler.is_some_and(|t| t.detected))
        .map(|d| d.satellite.to_string())
        .collect::<Vec<_>>()
        .join(";");
    row[42] = evidence
        .iter()
        .filter_map(|d| d.doppler_error.map(|e| format!("{}:{e:?}", d.satellite)))
        .collect::<Vec<_>>()
        .join(";");
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() != 3 {
        return Err(
            "usage: ppc_rtk_float CASE_DIRECTORY MODE BASE_X,BASE_Y,BASE_Z (MODE: static, kinematic, kinematic-gdop50, kinematic-spp, kinematic-cmc, kinematic-doppler, kinematic-slips; Galileo: kinematic-galileo with optional -spp / -slips; append -prefit3 to supported gate modes; append -align to enable per-system SD ambiguity datum alignment)".into(),
        );
    }
    let requested_mode = args[1].to_str().ok_or("non-UTF8 positioning mode")?;
    let mixed_spp = requested_mode.ends_with("-ggspp");
    let requested_mode = requested_mode
        .strip_suffix("-ggspp")
        .unwrap_or(requested_mode);
    let align_datums = requested_mode.ends_with("-align");
    let mode = requested_mode
        .strip_suffix("-align")
        .unwrap_or(requested_mode);
    let position_mode = match Some(mode) {
        Some("static") => RtkPositionMode::Static,
        Some(
            "kinematic"
            | "kinematic-gdop50"
            | "kinematic-spp"
            | "kinematic-cmc"
            | "kinematic-doppler"
            | "kinematic-slips"
            | "kinematic-prefit3"
            | "kinematic-spp-prefit3"
            | "kinematic-slips-prefit3"
            | "kinematic-spp-slips-prefit3",
        ) => RtkPositionMode::Kinematic,
        Some(
            "kinematic-galileo"
            | "kinematic-galileo-slips"
            | "kinematic-galileo-spp"
            | "kinematic-galileo-spp-slips"
            | "kinematic-galileo-spp-slips-prefit3",
        ) => RtkPositionMode::Kinematic,
        _ => return Err("unsupported positioning mode; see ppc_rtk_float usage".into()),
    };
    let mixed = Some(mode).is_some_and(|s| s.contains("galileo"));
    let constellations = if mixed {
        RtkConstellations::GpsGalileo
    } else {
        RtkConstellations::Gps
    };
    let use_spp_seed = matches!(
        Some(mode),
        Some(
            "kinematic-spp"
                | "kinematic-spp-prefit3"
                | "kinematic-spp-slips-prefit3"
                | "kinematic-galileo-spp"
                | "kinematic-galileo-spp-slips"
                | "kinematic-galileo-spp-slips-prefit3"
        )
    );
    let use_prefit = Some(mode).is_some_and(|s| s.ends_with("-prefit3"));
    if mixed_spp && !use_spp_seed {
        return Err("-ggspp requires an explicit -spp mode".into());
    }
    let slip_detection = match Some(mode) {
        Some("kinematic-cmc") => Some(RtkSlipConfig {
            doppler_threshold_m: None,
            ..Default::default()
        }),
        Some("kinematic-doppler") => Some(RtkSlipConfig {
            code_threshold_m: None,
            ..Default::default()
        }),
        Some(
            "kinematic-slips"
            | "kinematic-slips-prefit3"
            | "kinematic-spp-slips-prefit3"
            | "kinematic-galileo-slips"
            | "kinematic-galileo-spp-slips"
            | "kinematic-galileo-spp-slips-prefit3",
        ) => Some(RtkSlipConfig::default()),
        _ => None,
    };
    let xyz: Vec<f64> = args[2]
        .to_str()
        .ok_or("non-UTF8 base position")?
        .split(',')
        .map(str::parse)
        .collect::<Result<_, _>>()?;
    if xyz.len() != 3 {
        return Err("expected three base ECEF components".into());
    }
    let base_position = EcefCoord::new(xyz[0], xyz[1], xyz[2]);
    let root = Path::new(&args[0]);
    let open = |name| File::open(root.join(name)).map(BufReader::new);
    let rover = RinexObservationReader::new(open("rover.obs")?)?;
    let seed = rover
        .header()
        .approximate_position
        .ok_or("missing rover header position")?;
    let base = RinexObservationReader::new(open("base.obs")?)?;
    let systems = if mixed || mixed_spp {
        vec![GnssSystem::Gps, GnssSystem::Galileo]
    } else {
        vec![GnssSystem::Gps]
    };
    let navigation = read_navigation_for_systems(open("base.nav")?, &systems)?;
    let excluded = navigation.excluded_records.len();
    let mut filter = RtkFloatFilter::new(
        base_position,
        seed,
        RtkFloatConfig {
            constellations,
            position_mode,
            align_ambiguity_datums: align_datums,
            kinematic_spp_seed: use_spp_seed.then(SppConfig::default),
            kinematic_spp_constellations: if mixed_spp {
                RtkConstellations::GpsGalileo
            } else {
                RtkConstellations::Gps
            },
            slip_detection,
            max_prefit_residual_m: use_prefit.then_some(3.0),
            max_geometry_gdop: if mode == "kinematic-gdop50" {
                Some(50.0)
            } else {
                None
            },
            ..Default::default()
        },
    )?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    let header = "kind,week,tow,side,rover_flag,base_flag,base_dt_s,status,x_m,y_m,z_m,satellites,reference,nis,code_rms_m,phase_rms_m,ambiguity_resets,reason,geometry_gdop,geometry_pdop,position_covariance_trace_m2,position_qxx_m2,position_qyy_m2,position_qzz_m2,used_satellites,geometry_error";
    let suffix = if use_prefit && mixed {
        ",prefit_phase_max_m,prefit_code_max_m,prefit_threshold_m,prefit_rejected_phase_rows,prefit_rejected_code_rows,prefit_rejected_pairs,seed_source,seed_x_m,seed_y_m,seed_z_m,spp_failure,code_slips,doppler_slips,code_tests,doppler_tests,slip_satellites,doppler_errors,reference_satellites,tracking_pairs,geometry_clock_systems,dd_phase_rows"
    } else if mixed {
        ",seed_source,seed_x_m,seed_y_m,seed_z_m,spp_failure,code_slips,doppler_slips,code_tests,doppler_tests,slip_satellites,doppler_errors,reference_satellites,tracking_pairs,geometry_clock_systems,dd_phase_rows"
    } else if use_prefit {
        ",prefit_phase_max_m,prefit_code_max_m,prefit_threshold_m,prefit_rejected_phase_rows,prefit_rejected_code_rows,prefit_rejected_pairs,seed_source,seed_x_m,seed_y_m,seed_z_m,spp_failure,code_slips,doppler_slips,code_tests,doppler_tests,slip_satellites,doppler_errors"
    } else if use_spp_seed {
        ",seed_source,seed_x_m,seed_y_m,seed_z_m,spp_failure"
    } else if slip_detection.is_some() {
        ",code_slips,doppler_slips,code_tests,doppler_tests,slip_satellites,doppler_errors"
    } else {
        ""
    };
    let datum_suffix = if align_datums {
        ",datum_offsets_cycles,datum_retained_satellites,datum_initialized_satellites"
    } else {
        ""
    };
    let spp_suffix = if mixed_spp {
        ",spp_scope,spp_satellites,spp_gps_clock_bias_s,spp_galileo_clock_bias_s"
    } else {
        ""
    };
    writeln!(out, "{header}{suffix}{datum_suffix}{spp_suffix}")?;
    let (
        mut pairs,
        mut floated,
        mut failed,
        mut events,
        mut rover_only,
        mut base_only,
        mut continuity,
    ) = (0, 0, 0, 0, 0, 0, 0);
    for event in EpochMatcher::new(rover, base) {
        let columns = if mixed && use_prefit {
            47
        } else if mixed {
            41
        } else if use_prefit {
            43
        } else if use_spp_seed {
            31
        } else if slip_detection.is_some() {
            32
        } else {
            26
        };
        let columns = columns + if align_datums { 3 } else { 0 } + if mixed_spp { 4 } else { 0 };
        let mut row = vec![String::new(); columns];
        match event? {
            EpochMatch::Unmatched { side, epoch } => {
                row[0] = "unmatched".into();
                row[1] = epoch.time.week().to_string();
                row[2] = format!("{:.17}", epoch.time.tow());
                row[7] = "UNMATCHED".into();
                match side {
                    EpochSide::Rover => {
                        row[3] = "rover".into();
                        row[4] = epoch.flag.to_string();
                        rover_only += 1;
                    }
                    EpochSide::Base => {
                        row[3] = "base".into();
                        row[5] = epoch.flag.to_string();
                        base_only += 1;
                    }
                }
                let broken = filter.observe_skipped_epoch(&epoch);
                continuity += usize::from(broken);
                if broken {
                    row[17] = "continuity_reset_pending".into();
                }
            }
            EpochMatch::Pair { rover, base } => {
                pairs += 1;
                row[0] = "pair".into();
                row[1] = rover.time.week().to_string();
                row[2] = format!("{:.17}", rover.time.tow());
                row[3] = "both".into();
                row[4] = rover.flag.to_string();
                row[5] = base.flag.to_string();
                row[6] = format!("{:.17e}", base.time.difference_seconds(rover.time));
                if rover.flag == 6 || base.flag == 6 {
                    filter.invalidate_ambiguities();
                    events += 1;
                    row[7] = "EVENT".into();
                    row[17] = "cycle_slip_record".into();
                } else {
                    match filter.process_pair(&rover, &base, &navigation.navigation) {
                        Ok(result) => {
                            if mixed_spp {
                                write_spp_scope(&mut row, result.position_seed);
                            }
                            if align_datums {
                                let last = columns - if mixed_spp { 4 } else { 0 } - 3;
                                row[last] = result
                                    .ambiguity_datum_alignment
                                    .iter()
                                    .map(|d| {
                                        format!(
                                            "{}:{:.17e}",
                                            d.system.rinex_prefix(),
                                            d.offset_cycles
                                        )
                                    })
                                    .collect::<Vec<_>>()
                                    .join(";");
                                row[last + 1] = result
                                    .ambiguity_datum_alignment
                                    .iter()
                                    .flat_map(|d| d.retained_satellites.iter())
                                    .map(ToString::to_string)
                                    .collect::<Vec<_>>()
                                    .join(";");
                                row[last + 2] = result
                                    .ambiguity_datum_alignment
                                    .iter()
                                    .flat_map(|d| d.initialized_satellites.iter())
                                    .map(ToString::to_string)
                                    .collect::<Vec<_>>()
                                    .join(";");
                            }
                            if mixed {
                                if !use_prefit {
                                    write_seed(&mut row, 26, result.position_seed);
                                    write_slips(&mut row, 31, &result.slip_diagnostics);
                                }
                                let extra = if use_prefit { 43 } else { 37 };
                                row[extra] = result
                                    .reference_satellites
                                    .iter()
                                    .map(ToString::to_string)
                                    .collect::<Vec<_>>()
                                    .join(";");
                                row[extra + 1] = result
                                    .tracking_pairs()
                                    .iter()
                                    .map(|p| {
                                        format!(
                                            "{}:{}:{}",
                                            p.satellite,
                                            p.rover_tracking_code,
                                            p.base_tracking_code
                                        )
                                    })
                                    .collect::<Vec<_>>()
                                    .join(";");
                                row[extra + 2] = result
                                    .reference_satellites
                                    .iter()
                                    .map(|s| s.system().rinex_prefix().to_string())
                                    .collect::<Vec<_>>()
                                    .join(";");
                                row[extra + 3] = result.ambiguity_satellites.len().to_string();
                            }
                            if slip_detection.is_some() && !use_prefit && !mixed {
                                let evidence = &result.slip_diagnostics;
                                row[26] = evidence
                                    .iter()
                                    .filter(|d| d.code.is_some_and(|t| t.detected))
                                    .count()
                                    .to_string();
                                row[27] = evidence
                                    .iter()
                                    .filter(|d| d.doppler.is_some_and(|t| t.detected))
                                    .count()
                                    .to_string();
                                row[28] = evidence
                                    .iter()
                                    .filter(|d| d.code.is_some())
                                    .count()
                                    .to_string();
                                row[29] = evidence
                                    .iter()
                                    .filter(|d| d.doppler.is_some())
                                    .count()
                                    .to_string();
                                row[30] = evidence
                                    .iter()
                                    .filter(|d| {
                                        d.code.is_some_and(|t| t.detected)
                                            || d.doppler.is_some_and(|t| t.detected)
                                    })
                                    .map(|d| d.satellite.to_string())
                                    .collect::<Vec<_>>()
                                    .join(";");
                                row[31] = evidence
                                    .iter()
                                    .filter_map(|d| {
                                        d.doppler_error.map(|e| format!("{}:{e:?}", d.satellite))
                                    })
                                    .collect::<Vec<_>>()
                                    .join(";");
                            }
                            if use_spp_seed && !use_prefit && !mixed {
                                let seed = result.position_seed;
                                row[26] = format!("{:?}", seed.source);
                                for (i, v) in [seed.position.x, seed.position.y, seed.position.z]
                                    .into_iter()
                                    .enumerate()
                                {
                                    row[27 + i] = format!("{v:.17e}");
                                }
                                if let Some(error) = seed.spp_failure {
                                    row[30] = format!("{error:?}");
                                }
                            }
                            floated += 1;
                            row[7] = "FLOAT".into();
                            let p = result
                                .position
                                .position_ecef
                                .ok_or("missing FLOAT position")?;
                            for (index, value) in [p.x, p.y, p.z].into_iter().enumerate() {
                                row[8 + index] = format!("{value:.9}");
                            }
                            row[11] = result.position.num_satellites.to_string();
                            row[12] = result.reference_satellite.to_string();
                            row[13] = format!("{:.17e}", result.normalized_innovation_squared);
                            row[14] = format!("{:.17e}", result.code_residual_rms_m);
                            row[15] = format!("{:.17e}", result.phase_residual_rms_m);
                            row[16] = result.reset_ambiguities.len().to_string();
                            match &result.geometry {
                                Ok(dop) => {
                                    row[18] = format!("{:.17e}", dop.gdop);
                                    row[19] = format!("{:.17e}", dop.pdop);
                                }
                                Err(error) => row[25] = format!("{error:?}"),
                            }
                            let covariance = result
                                .position
                                .position_covariance
                                .ok_or("missing FLOAT covariance")?;
                            row[20] =
                                format!("{:.17e}", (0..3).map(|i| covariance[i][i]).sum::<f64>());
                            for i in 0..3 {
                                row[21 + i] = format!("{:.17e}", covariance[i][i]);
                            }
                            row[24] = result
                                .used_satellites()
                                .iter()
                                .map(ToString::to_string)
                                .collect::<Vec<_>>()
                                .join(";");
                        }
                        Err(error) => {
                            failed += 1;
                            row[7] = "INVALID".into();
                            row[17] = format!("{error:?}");
                        }
                    }
                }
                if use_prefit
                    && let Some(diagnostic) = filter
                        .last_prefit_diagnostic()
                        .filter(|d| d.time == rover.time)
                {
                    write_prefit(&mut row, diagnostic);
                    if mixed_spp {
                        write_spp_scope(&mut row, diagnostic.position_seed);
                    }
                }
            }
        }
        // Quote the diagnostic field in case a future error contains commas.
        row[17] = format!("\"{}\"", row[17].replace('"', "\"\""));
        writeln!(out, "{}", row.join(","))?;
    }
    out.flush()?;
    eprintln!(
        "position_mode={position_mode:?} pairs={pairs} FLOAT={floated} INVALID={failed} event_pairs={events} rover_only={rover_only} base_only={base_only} unmatched_continuity_notices={continuity} excluded_navigation_records={excluded} FIX_policy_evaluated=false whole_native_comparison=false performance_accepted=false"
    );
    if align_datums {
        eprintln!(
            "ambiguity_datum_alignment=per_system_retained_mean default=false covariance_changed=false"
        );
    }
    if use_spp_seed {
        if mixed_spp {
            eprintln!(
                "kinematic_spp_seed=GPS_1C_Galileo_E1 independent_clocks=true legacy_primary_BGD=true default_SppConfig no_SPP_covariance_fusion=true"
            );
        } else {
            eprintln!(
                "kinematic_spp_seed=GPS_1C_only default_SppConfig no_SPP_covariance_fusion=true"
            );
        }
    }
    if let Some(config) = slip_detection {
        eprintln!(
            "slip_detection={config:?} accepted_epoch_history=true thresholds_from_pinned_native_defaults=true"
        );
    }
    if use_prefit {
        eprintln!("max_prefit_residual_m=3 whole_epoch_gate=true no_matrix_mutation=true");
    }
    if mixed {
        eprintln!(
            "constellations=GPS,Galileo separate_references=true SPP_seed_GPS_only={} kinematic_FIX=false",
            !mixed_spp
        );
    }
    Ok(())
}
