//! GPS L1 FLOAT on common epochs; every unmatched/invalid epoch is reported.
use gnss_rust::epoch_match::{EpochMatch, EpochMatcher, EpochSide};
use gnss_rust::rinex::{RinexObservationReader, read_navigation_for_systems};
use gnss_rust::rtk::{RtkFloatConfig, RtkFloatFilter, RtkPositionMode};
use gnss_rust::rtk_slip::RtkSlipConfig;
use gnss_rust::spp::SppConfig;
use gnss_rust::{EcefCoord, GnssSystem};
use std::{
    error::Error,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() != 3 {
        return Err(
            "usage: ppc_rtk_float CASE_DIRECTORY static|kinematic|kinematic-gdop50|kinematic-spp|kinematic-cmc|kinematic-doppler|kinematic-slips BASE_X,BASE_Y,BASE_Z".into(),
        );
    }
    let position_mode = match args[1].to_str() {
        Some("static") => RtkPositionMode::Static,
        Some(
            "kinematic" | "kinematic-gdop50" | "kinematic-spp" | "kinematic-cmc"
            | "kinematic-doppler" | "kinematic-slips",
        ) => RtkPositionMode::Kinematic,
        _ => return Err("unsupported positioning mode; see ppc_rtk_float usage".into()),
    };
    let use_spp_seed = args[1] == "kinematic-spp";
    let slip_detection = match args[1].to_str() {
        Some("kinematic-cmc") => Some(RtkSlipConfig {
            doppler_threshold_m: None,
            ..Default::default()
        }),
        Some("kinematic-doppler") => Some(RtkSlipConfig {
            code_threshold_m: None,
            ..Default::default()
        }),
        Some("kinematic-slips") => Some(RtkSlipConfig::default()),
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
    let navigation = read_navigation_for_systems(open("base.nav")?, &[GnssSystem::Gps])?;
    let excluded = navigation.excluded_records.len();
    let mut filter = RtkFloatFilter::new(
        base_position,
        seed,
        RtkFloatConfig {
            position_mode,
            kinematic_spp_seed: use_spp_seed.then(SppConfig::default),
            slip_detection,
            max_geometry_gdop: if args[1] == "kinematic-gdop50" {
                Some(50.0)
            } else {
                None
            },
            ..Default::default()
        },
    )?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    let header = "kind,week,tow,side,rover_flag,base_flag,base_dt_s,status,x_m,y_m,z_m,satellites,reference,nis,code_rms_m,phase_rms_m,ambiguity_resets,reason,geometry_gdop,geometry_pdop,position_covariance_trace_m2,position_qxx_m2,position_qyy_m2,position_qzz_m2,used_satellites,geometry_error";
    let suffix = if use_spp_seed {
        ",seed_source,seed_x_m,seed_y_m,seed_z_m,spp_failure"
    } else if slip_detection.is_some() {
        ",code_slips,doppler_slips,code_tests,doppler_tests,slip_satellites,doppler_errors"
    } else {
        ""
    };
    writeln!(out, "{header}{suffix}")?;
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
        let columns = if use_spp_seed {
            31
        } else if slip_detection.is_some() {
            32
        } else {
            26
        };
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
                            if slip_detection.is_some() {
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
                            if use_spp_seed {
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
    if use_spp_seed {
        eprintln!("kinematic_spp_seed=GPS_1C_only default_SppConfig no_SPP_covariance_fusion=true");
    }
    if let Some(config) = slip_detection {
        eprintln!(
            "slip_detection={config:?} accepted_epoch_history=true thresholds_from_pinned_native_defaults=true"
        );
    }
    Ok(())
}
