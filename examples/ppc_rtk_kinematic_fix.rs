//! Opt-in full-set kinematic FIX; FLOAT state is never modified by policy output.
use gnss_rust::epoch_match::{EpochMatch, EpochMatcher, EpochSide};
use gnss_rust::rinex::{RinexObservationReader, read_navigation_for_systems};
use gnss_rust::rtk::{RtkConstellations, RtkFloatConfig, RtkFloatFilter, RtkPositionMode};
use gnss_rust::rtk_fix::{
    KinematicMotionConfig, ReferenceGateMode, RtkFixConfig, RtkFixPolicy, RtkMotionEvidence,
};
use gnss_rust::rtk_slip::RtkSlipConfig;
use gnss_rust::velocity::{
    VelocityClockModel, VelocityConfig, VelocityEpochConfig, solve_gps_galileo_e1_epoch,
};
use gnss_rust::{EcefCoord, GnssSystem};
use std::{
    error::Error,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 3 && args.len() != 4 {
        return Err("usage: ppc_rtk_kinematic_fix CASE_DIRECTORY selected|all|selected-slips|all-slips[with -locked suffix] BASE_X,BASE_Y,BASE_Z [MAX_MOTION_STD_MPS]".into());
    }
    let locked_subset = args[1].ends_with("-locked");
    let mode = args[1].strip_suffix("-locked").unwrap_or(&args[1]);
    let (gates, slips) = match mode {
        "selected" => (ReferenceGateMode::Selected, false),
        "all" => (ReferenceGateMode::AllReferences, false),
        "selected-slips" => (ReferenceGateMode::Selected, true),
        "all-slips" => (ReferenceGateMode::AllReferences, true),
        _ => return Err("unknown FIX mode".into()),
    };
    let xyz: Vec<f64> = args[2]
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
        .ok_or("missing rover position")?;
    let base = RinexObservationReader::new(open("base.obs")?)?;
    let nav =
        read_navigation_for_systems(open("base.nav")?, &[GnssSystem::Gps, GnssSystem::Galileo])?;
    let mut filter = RtkFloatFilter::new(
        base_position,
        seed,
        RtkFloatConfig {
            constellations: RtkConstellations::GpsGalileo,
            position_mode: RtkPositionMode::Kinematic,
            slip_detection: slips.then(RtkSlipConfig::default),
            ..Default::default()
        },
    )?;
    let motion_std = if args.len() == 4 {
        args[3].parse()?
    } else {
        1.0
    };
    let constructor = if locked_subset {
        RtkFixPolicy::for_kinematic_locked_subset
    } else {
        RtkFixPolicy::for_kinematic
    };
    let mut policy = constructor(
        RtkFixConfig::default(),
        KinematicMotionConfig {
            max_velocity_std_mps: motion_std,
            ..Default::default()
        },
        gates,
    )?;
    let velocity_config = VelocityEpochConfig {
        velocity: VelocityConfig {
            clock_model: VelocityClockModel::PerSystem,
            ..Default::default()
        },
        require_carrier_lock: true,
        ..Default::default()
    };
    let mut out = BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "kind,week,tow,side,status,x_m,y_m,z_m,float_x_m,float_y_m,float_z_m,used_satellites,reference_satellites,ambiguity_pairs,ambiguity_resets,decision,ratio,confirmation_count,velocity_status,vx_mps,vy_mps,vz_mps,velocity_max_std_mps,velocity_reason,fixed_correction_m,fixed_phase_rms_m,fixed_code_rms_m,motion_max_residual_m,motion_checks,integers,float_code_rms_m,float_phase_rms_m,float_qxx_m2,float_qyy_m2,float_qzz_m2,solution_qxx_m2,solution_qyy_m2,solution_qzz_m2,system_reference_sets,reason,velocity_qxx,velocity_qxy,velocity_qxz,velocity_qyx,velocity_qyy,velocity_qyz,velocity_qzx,velocity_qzy,velocity_qzz,ar_satellites,ar_reference_satellites,ar_ambiguity_pairs,satellite_lock_epochs"
    )?;
    for event in EpochMatcher::new(rover, base) {
        let mut row = vec![String::new(); 53];
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
                    row[39] = "continuity_reset_pending".into();
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
                    row[39] = "cycle_slip_record".into();
                } else {
                    match filter.process_pair(&rover, &base, &nav.navigation) {
                        Err(error) => {
                            row[4] = "INVALID".into();
                            row[39] = format!("{error:?}");
                        }
                        Ok(float) => {
                            let p = float
                                .position
                                .position_ecef
                                .ok_or("missing FLOAT position")?;
                            row[4] = "FLOAT".into();
                            for (i, v) in [p.x, p.y, p.z].into_iter().enumerate() {
                                row[5 + i] = format!("{v:.9}");
                                row[8 + i] = format!("{v:.9}");
                            }
                            row[11] = float
                                .used_satellites()
                                .iter()
                                .map(ToString::to_string)
                                .collect::<Vec<_>>()
                                .join(";");
                            row[12] = float
                                .reference_satellites
                                .iter()
                                .map(ToString::to_string)
                                .collect::<Vec<_>>()
                                .join(";");
                            row[13] = float.ambiguity_satellites.len().to_string();
                            row[14] = float.reset_ambiguities.len().to_string();
                            row[30] = format!("{:.17e}", float.code_residual_rms_m);
                            row[31] = format!("{:.17e}", float.phase_residual_rms_m);
                            let covariance = float
                                .position
                                .position_covariance
                                .ok_or("missing FLOAT covariance")?;
                            for i in 0..3 {
                                row[32 + i] = format!("{:.17e}", covariance[i][i]);
                                row[35 + i] = row[32 + i].clone();
                            }
                            let before = filter.state().cloned();
                            let mut evidence = None;
                            match solve_gps_galileo_e1_epoch(
                                &rover,
                                &nav.navigation,
                                p,
                                velocity_config,
                            ) {
                                Err(error) => {
                                    row[18] = "INVALID_EPOCH".into();
                                    row[23] = format!("{error:?}");
                                }
                                Ok(report) => match report.solution {
                                    Err(error) => {
                                        row[18] = "INVALID_VELOCITY".into();
                                        row[23] = format!("{error:?}");
                                    }
                                    Ok(velocity) => {
                                        row[18] = "VELOCITY".into();
                                        for (i, v) in velocity.velocity_mps.into_iter().enumerate()
                                        {
                                            row[19 + i] = format!("{v:.17e}");
                                        }
                                        for i in 0..3 {
                                            for j in 0..3 {
                                                row[40 + 3 * i + j] = format!(
                                                    "{:.17e}",
                                                    velocity.covariance_m2s2[i][j]
                                                );
                                            }
                                        }
                                        evidence = Some(RtkMotionEvidence::from_velocity_solution(
                                            rover.time, &velocity,
                                        ));
                                    }
                                },
                            }
                            let attempt = policy.evaluate_kinematic(&float, evidence.as_ref())?;
                            row[49] = attempt
                                .ar_satellites
                                .iter()
                                .map(ToString::to_string)
                                .collect::<Vec<_>>()
                                .join(";");
                            row[50] = attempt
                                .ar_reference_satellites
                                .iter()
                                .map(ToString::to_string)
                                .collect::<Vec<_>>()
                                .join(";");
                            row[51] = (attempt.ar_satellites.len()
                                - attempt.ar_reference_satellites.len())
                            .to_string();
                            row[52] = attempt
                                .satellite_lock_epochs
                                .iter()
                                .map(|(s, n)| format!("{s}:{n}"))
                                .collect::<Vec<_>>()
                                .join(";");
                            row[15] = attempt.decision.to_string();
                            row[16] = attempt
                                .ratio
                                .map_or_else(String::new, |v| format!("{v:.17e}"));
                            row[17] = attempt.confirmation_count.to_string();
                            row[22] = attempt
                                .motion_max_velocity_std_mps
                                .map_or_else(String::new, |v| format!("{v:.17e}"));
                            if let Some(validation) = &attempt.validation {
                                row[24] = format!("{:.17e}", validation.float_correction_m);
                                row[25] = format!("{:.17e}", validation.phase_residual_rms_m);
                                row[26] = format!("{:.17e}", validation.code_residual_rms_m);
                            }
                            if !attempt.motion_checks.is_empty() {
                                row[27] = format!(
                                    "{:.17e}",
                                    attempt
                                        .motion_checks
                                        .iter()
                                        .map(|c| c.residual_norm_m)
                                        .fold(0.0, f64::max)
                                );
                                row[28]=attempt.motion_checks.iter().map(|c|format!("{:?}@{:.17e}@{:.17e}@{:.17e}@{:.17e}@{:.17e}@{:.17e}@{:.17e}@{:.17e}@{:.17e}",c.kind,c.interval_s,c.predicted_displacement_m[0],c.predicted_displacement_m[1],c.predicted_displacement_m[2],c.displacement_residual_m[0],c.displacement_residual_m[1],c.displacement_residual_m[2],c.residual_norm_m,c.maximum_m)).collect::<Vec<_>>().join(";");
                            }
                            row[29] = attempt
                                .integers
                                .iter()
                                .map(ToString::to_string)
                                .collect::<Vec<_>>()
                                .join(";");
                            if let Some(audit) = &attempt.system_reference_audit {
                                row[38] = audit.reference_sets.len().to_string();
                            }
                            if let Some(fixed) = attempt.fixed_solution {
                                row[4] = "FIX".into();
                                let p = fixed.position_ecef.ok_or("missing FIX position")?;
                                for (i, v) in [p.x, p.y, p.z].into_iter().enumerate() {
                                    row[5 + i] = format!("{v:.9}");
                                }
                                let q =
                                    fixed.position_covariance.ok_or("missing FIX covariance")?;
                                for i in 0..3 {
                                    row[35 + i] = format!("{:.17e}", q[i][i]);
                                }
                            }
                            if filter.state() != before.as_ref() {
                                return Err("FIX policy modified the FLOAT posterior".into());
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
        "kinematic_FIX_policy=true locked_subset={locked_subset} full_set_only={} integer_holding=false fixed_or_velocity_feedback=false GPS_Galileo=true carrier_lock_required=true motion_max_std_mps={motion_std} motion_max_gap_s=2 displacement_residual_limit_m=0.1 confirmation_epochs=2 reference_gates={gates:?} performance_accepted=false",
        !locked_subset
    );
    Ok(())
}
