//! Read-only mixed-system integer candidate diagnostics, never a FIX policy.
use gnss_rust::epoch_match::{EpochMatch, EpochMatcher, EpochSide};
use gnss_rust::lambda::LambdaConfig;
use gnss_rust::rinex::{RinexObservationReader, read_navigation_for_systems};
use gnss_rust::rtk::{RtkConstellations, RtkFloatConfig, RtkFloatFilter, RtkPositionMode};
use gnss_rust::rtk_fix::{RtkFixConfig, audit_candidate_system_references, validate_candidate};
use gnss_rust::rtk_slip::RtkSlipConfig;
use gnss_rust::{EcefCoord, GnssSystem};
use std::{
    error::Error,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 3 {
        return Err(
            "usage: ppc_rtk_candidates CASE_DIRECTORY galileo|galileo-slips BASE_X,BASE_Y,BASE_Z"
                .into(),
        );
    }
    let slip_detection = match args[1].as_str() {
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
    writeln!(
        out,
        "kind,week,tow,side,status,x_m,y_m,z_m,used_satellites,reference_satellites,ambiguity_pairs,ambiguity_resets,search_status,ratio,integers,local_status,local_rejection,conditioned_x_m,conditioned_y_m,conditioned_z_m,correction_m,phase_rms_m,code_rms_m,postfit_nis_per_row,max_position_std_m,audit_status,reference_sets,conservative_rejection,reason"
    )?;
    for event in EpochMatcher::new(rover, base) {
        let mut row = vec![String::new(); 29];
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
    Ok(())
}
