//! Read-only mixed-system integer candidate diagnostics, never a FIX policy.
use gnss_rust::epoch_match::{EpochMatch, EpochMatcher, EpochSide};
use gnss_rust::lambda::LambdaConfig;
use gnss_rust::rinex::{RinexObservationReader, read_navigation_for_systems};
use gnss_rust::rtk::{RtkConstellations, RtkFloatConfig, RtkFloatFilter, RtkPositionMode};
use gnss_rust::rtk_fix::{RtkFixConfig, audit_candidate_system_references, validate_candidate};
use gnss_rust::rtk_slip::RtkSlipConfig;
use gnss_rust::velocity::{
    VelocityClockModel, VelocityConfig, VelocityEpochConfig, VelocityEpochReport,
    solve_gps_galileo_e1_epoch, solve_gps_l1_epoch,
};
use gnss_rust::{EcefCoord, GnssSystem};
use std::{
    error::Error,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

fn write_velocity(row: &mut [String], report: Result<VelocityEpochReport, gnss_rust::Error>) {
    let report = match report {
        Ok(value) => value,
        Err(error) => {
            row[29] = "INVALID_EPOCH".into();
            row[47] = format!("{error:?}");
            return;
        }
    };
    row[30] = report.prepared.len().to_string();
    row[31] = report.rejected.len().to_string();
    row[49] = report
        .prepared
        .iter()
        .filter(|p| p.measurement.doppler_hz == 0.0)
        .count()
        .to_string();
    row[50] = report
        .prepared
        .iter()
        .filter(|p| p.used_missing_code_travel_time)
        .count()
        .to_string();
    if row.len() > 51 {
        row[53] = report
            .prepared
            .iter()
            .map(|p| format!("{}:{}", p.raw.satellite, p.raw.tracking_code))
            .collect::<Vec<_>>()
            .join(";");
        row[54] = report
            .rejected
            .iter()
            .map(|(raw, why)| format!("{}:{}:{why:?}", raw.satellite, raw.tracking_code))
            .collect::<Vec<_>>()
            .join(";");
    }
    match report.solution {
        Err(error) => {
            row[29] = "INVALID_VELOCITY".into();
            row[47] = format!("{error:?}");
        }
        Ok(solution) => {
            row[29] = "VELOCITY".into();
            for (i, v) in solution.velocity_mps.into_iter().enumerate() {
                row[32 + i] = format!("{v:.17e}");
            }
            row[35] = format!("{:.17e}", solution.receiver_clock_drift_sps);
            row[36] = format!("{:.17e}", solution.chi_square_per_dof);
            row[37] = format!("{:.17e}", solution.residual_rms_mps);
            for i in 0..3 {
                for j in 0..3 {
                    row[38 + 3 * i + j] = format!("{:.17e}", solution.covariance_m2s2[i][j]);
                }
            }
            row[48] = solution
                .satellites
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(";");
            if row.len() > 51 {
                row[51] = solution
                    .clock_systems
                    .iter()
                    .map(|s| s.rinex_prefix().to_string())
                    .collect::<Vec<_>>()
                    .join(";");
                row[52] = solution
                    .receiver_clock_drifts_sps
                    .iter()
                    .map(|v| format!("{v:.17e}"))
                    .collect::<Vec<_>>()
                    .join(";");
            }
        }
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 3 {
        return Err(
            "usage: ppc_rtk_candidates CASE_DIRECTORY galileo|galileo-slips BASE_X,BASE_Y,BASE_Z (append -velocity for GPS 1C Doppler diagnostics, or -velocity-mixed for GPS/Galileo with separate clocks)"
                .into(),
        );
    }
    let mixed_velocity = args[1].ends_with("-velocity-mixed");
    let requested = if mixed_velocity {
        args[1].strip_suffix("-mixed").unwrap()
    } else {
        &args[1]
    };
    let use_velocity = requested.ends_with("-velocity");
    let mode = requested.strip_suffix("-velocity").unwrap_or(requested);
    let slip_detection = match mode {
        "galileo" => None,
        "galileo-slips" => Some(RtkSlipConfig::default()),
        _ => return Err("unknown diagnostic mode".into()),
    };
    let xyz: Vec<f64> = args[2]
        .split(',')
        .map(str::parse)
        .collect::<Result<_, _>>()?;
    if xyz.len() != 3 {
        return Err("expected three base ECEF components".into());
    }
    let position = EcefCoord::new(xyz[0], xyz[1], xyz[2]);
    let root = Path::new(&args[0]);
    let open = |name| File::open(root.join(name)).map(BufReader::new);
    let rover = RinexObservationReader::new(open("rover.obs")?)?;
    let seed = rover
        .header()
        .approximate_position
        .ok_or("missing rover position")?;
    let base = RinexObservationReader::new(open("base.obs")?)?;
    let navigation =
        read_navigation_for_systems(open("base.nav")?, &[GnssSystem::Gps, GnssSystem::Galileo])?;
    let mut filter = RtkFloatFilter::new(
        position,
        seed,
        RtkFloatConfig {
            constellations: RtkConstellations::GpsGalileo,
            position_mode: RtkPositionMode::Kinematic,
            slip_detection,
            ..Default::default()
        },
    )?;
    let config = RtkFixConfig::default();
    let mut out = BufWriter::new(std::io::stdout().lock());
    let suffix = if use_velocity {
        ",velocity_status,velocity_prepared_count,velocity_rejected_count,vx_mps,vy_mps,vz_mps,velocity_clock_drift_sps,velocity_chi_square_per_dof,velocity_residual_rms_mps,vqxx,vqxy,vqxz,vqyx,vqyy,vqyz,vqzx,vqzy,vqzz,velocity_reason,velocity_used_satellites,velocity_zero_rows,velocity_missing_code_rows"
    } else {
        ""
    };
    let clock_suffix = if mixed_velocity {
        ",velocity_clock_systems,velocity_clock_drifts_sps,velocity_tracking_codes,velocity_rejections"
    } else {
        ""
    };
    writeln!(
        out,
        "kind,week,tow,side,status,x_m,y_m,z_m,used_satellites,reference_satellites,ambiguity_pairs,ambiguity_resets,search_status,ratio,integers,local_status,local_rejection,conditioned_x_m,conditioned_y_m,conditioned_z_m,correction_m,phase_rms_m,code_rms_m,postfit_nis_per_row,max_position_std_m,audit_status,reference_sets,conservative_rejection,reason{suffix}{clock_suffix}"
    )?;
    for event in EpochMatcher::new(rover, base) {
        let mut row = vec![
            String::new();
            if mixed_velocity {
                55
            } else if use_velocity {
                51
            } else {
                29
            }
        ];
        match event? {
            EpochMatch::Unmatched { side, epoch } => {
                row[0] = "unmatched".into();
                row[1] = epoch.time.week().to_string();
                row[2] = format!("{:.17}", epoch.time.tow());
                row[3] = match side {
                    EpochSide::Rover => "rover",
                    EpochSide::Base => "base",
                }
                .into();
                row[4] = "UNMATCHED".into();
                if filter.observe_skipped_epoch(&epoch) {
                    row[28] = "continuity_reset_pending".into();
                }
            }
            EpochMatch::Pair { rover, base } => {
                row[0] = "pair".into();
                row[1] = rover.time.week().to_string();
                row[2] = format!("{:.17}", rover.time.tow());
                row[3] = "both".into();
                if rover.flag == 6 || base.flag == 6 {
                    filter.invalidate_ambiguities();
                    row[4] = "EVENT".into();
                    row[28] = "cycle_slip_record".into();
                } else {
                    match filter.process_pair(&rover, &base, &navigation.navigation) {
                        Err(error) => {
                            row[4] = "INVALID".into();
                            row[28] = format!("{error:?}");
                        }
                        Ok(float) => {
                            row[4] = "FLOAT".into();
                            let p = float
                                .position
                                .position_ecef
                                .ok_or("missing FLOAT position")?;
                            for (i, v) in [p.x, p.y, p.z].into_iter().enumerate() {
                                row[5 + i] = format!("{v:.9}");
                            }
                            row[8] = float
                                .used_satellites()
                                .iter()
                                .map(ToString::to_string)
                                .collect::<Vec<_>>()
                                .join(";");
                            row[9] = float
                                .reference_satellites
                                .iter()
                                .map(ToString::to_string)
                                .collect::<Vec<_>>()
                                .join(";");
                            row[10] = float.ambiguity_satellites.len().to_string();
                            row[11] = float.reset_ambiguities.len().to_string();
                            let before = filter.state().cloned();
                            if use_velocity {
                                let velocity = if mixed_velocity {
                                    solve_gps_galileo_e1_epoch(
                                        &rover,
                                        &navigation.navigation,
                                        p,
                                        VelocityEpochConfig {
                                            velocity: VelocityConfig {
                                                clock_model: VelocityClockModel::PerSystem,
                                                ..Default::default()
                                            },
                                            ..Default::default()
                                        },
                                    )
                                } else {
                                    solve_gps_l1_epoch(
                                        &rover,
                                        &navigation.navigation,
                                        p,
                                        VelocityEpochConfig::default(),
                                    )
                                };
                                write_velocity(&mut row, velocity);
                            }
                            match float.integer_candidates(LambdaConfig::default()) {
                                Err(error) => {
                                    row[12] = format!("{error:?}");
                                }
                                Ok(search) => {
                                    row[12] = "OK".into();
                                    let ratio = search.ratio();
                                    if let Some(ratio) = ratio {
                                        row[13] = format!("{ratio:.17e}");
                                    }
                                    let integers = &search.candidates[0].ambiguities;
                                    row[14] = integers
                                        .iter()
                                        .map(ToString::to_string)
                                        .collect::<Vec<_>>()
                                        .join(";");
                                    match validate_candidate(&float, integers, config) {
                                        Err(error) => {
                                            row[15] = "ERROR".into();
                                            row[16] = format!("{error:?}");
                                        }
                                        Ok(validation) => {
                                            row[15] = if validation.rejection.is_none() {
                                                "PASS_LOCAL_GATES"
                                            } else {
                                                "REJECTED"
                                            }
                                            .into();
                                            if let Some(reason) = validation.rejection {
                                                row[16] = reason.to_string();
                                            }
                                            for (i, v) in [position.x, position.y, position.z]
                                                .into_iter()
                                                .enumerate()
                                            {
                                                row[17 + i] = format!(
                                                    "{:.9}",
                                                    v + validation.conditioned.state.mean[i]
                                                );
                                            }
                                            for (i, v) in [
                                                validation.float_correction_m,
                                                validation.phase_residual_rms_m,
                                                validation.code_residual_rms_m,
                                                validation.postfit_nis_per_row,
                                                validation.max_position_std_m,
                                            ]
                                            .into_iter()
                                            .enumerate()
                                            {
                                                row[20 + i] = format!("{v:.17e}");
                                            }
                                            if validation.rejection.is_none()
                                                && ratio.is_some_and(|r| {
                                                    r.is_finite() && r >= config.ratio_threshold
                                                })
                                            {
                                                match audit_candidate_system_references(
                                                    &float, integers, config, 4096,
                                                ) {
                                                    Err(error) => {
                                                        row[25] = "ERROR".into();
                                                        row[27] = format!("{error:?}");
                                                    }
                                                    Ok(audit) => {
                                                        row[25] = "COMPLETE".into();
                                                        row[26] =
                                                            audit.reference_sets.len().to_string();
                                                        row[27] = audit
                                                            .conservative_rejection()
                                                            .map_or_else(
                                                                || "PASS_LOCAL_GATES".into(),
                                                                |r| r.to_string(),
                                                            );
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            if filter.state() != before.as_ref() {
                                return Err("candidate diagnostic modified the FLOAT state".into());
                            }
                        }
                    }
                }
            }
        }
        for field in &mut row {
            if field.contains([',', '"', '\n']) {
                *field = format!("\"{}\"", field.replace('"', "\"\""));
            }
        }
        writeln!(out, "{}", row.join(","))?;
    }
    out.flush()?;
    eprintln!(
        "candidate_diagnostics=true FIX_policy_evaluated=false no_lock_confirmation_or_motion_history=true no_feedback_to_FLOAT=true audit_only_if_ratio_and_local_gates_pass=true reference_used=false"
    );
    if use_velocity {
        let model = if mixed_velocity {
            "GPS_Galileo_separate_clocks_current_FLOAT_position"
        } else {
            "GPS_1C_current_FLOAT_position"
        };
        eprintln!(
            "velocity_diagnostics={model} no_velocity_feedback=true no_motion_gate_evaluated=true"
        );
    }
    Ok(())
}
