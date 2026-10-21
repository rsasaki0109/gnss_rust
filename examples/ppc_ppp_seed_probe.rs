//! Read-only SPP/Doppler probe at a specified raw PPC epoch; no reference input.
use gnss_rust::rinex::{RinexObservationReader, read_navigation_for_systems};
use gnss_rust::spp::SppConfig;
use gnss_rust::spp_qc::{SppQualityConfig, solve_gps_galileo_epoch_with_quality_control};
use gnss_rust::spp_quality::SppVarianceModel;
use gnss_rust::velocity::{VelocityClockModel, VelocityEpochConfig, solve_gps_galileo_e1_epoch};
use gnss_rust::{GnssSystem, GnssTime};
use std::{fs::File, io::BufReader, path::Path};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 3 {
        return Err("usage: ppc_ppp_seed_probe CASE_DIRECTORY GPS_WEEK GPS_TOW".into());
    }
    let root = Path::new(&args[0]);
    let target = GnssTime::new(args[1].parse()?, args[2].parse()?)?;
    let nav = read_navigation_for_systems(
        BufReader::new(File::open(root.join("base.nav"))?),
        &[GnssSystem::Gps, GnssSystem::Galileo],
    )?;
    let reader = RinexObservationReader::new(BufReader::new(File::open(root.join("rover.obs"))?))?;
    let mut found = None;
    for epoch in reader {
        let e = epoch?;
        if e.time.difference_seconds(target).abs() <= 1e-6 {
            found = Some(e);
            break;
        }
    }
    let epoch = found.ok_or("target raw epoch absent")?;
    eprintln!(
        "raw_epoch={epoch:?} truth_used_by_solver=false PPP_filter_updated=false Doppler_not_fused=true"
    );
    println!(
        "profile,week,tow,status,x_m,y_m,z_m,code_rms_m,position_covariance_trace_m2,satellites,gdop,MAD_rejected,RAIM_rejected,vx_mps,vy_mps,vz_mps,velocity_covariance_trace_m2_s2,velocity_satellites,velocity_residual_rms_mps,error"
    );
    for gps_only in [true, false] {
        for qc in [false, true] {
            let mut input = epoch.clone();
            if gps_only {
                input
                    .observations
                    .retain(|o| o.satellite.system() == GnssSystem::Gps);
            }
            let profile = format!(
                "{}-{}",
                if gps_only { "GPS" } else { "GPS-Galileo" },
                if qc {
                    "native-variance-QC"
                } else {
                    "elevation-no-QC"
                }
            );
            let mut row = vec![String::new(); 20];
            row[0] = profile.clone();
            row[1] = epoch.time.week().to_string();
            row[2] = format!("{:.17}", epoch.time.tow());
            let quality = if qc {
                SppQualityConfig::default()
            } else {
                SppQualityConfig {
                    outlier_threshold_sigma: None,
                    enable_raim_fde: false,
                    max_gdop: None,
                    ..Default::default()
                }
            };
            match solve_gps_galileo_epoch_with_quality_control(
                &input,
                &nav.navigation,
                SppConfig::default(),
                if qc {
                    SppVarianceModel::NativeSnrAtmosphere
                } else {
                    SppVarianceModel::ElevationOnly
                },
                quality,
            ) {
                Ok(report) => {
                    eprintln!(
                        "profile={profile} quality={:?} SPP={:?}",
                        report.evidence, report.solution
                    );
                    let solution = &report.solution.spp;
                    let position = solution.position.position_ecef.unwrap();
                    let q = solution.position.position_covariance.unwrap();
                    row[3] = "SPP".into();
                    for (i, x) in [position.x, position.y, position.z].iter().enumerate() {
                        row[4 + i] = format!("{x:.17e}");
                    }
                    row[7] = format!("{:.17e}", solution.residual_rms_m);
                    row[8] = format!("{:.17e}", q[0][0] + q[1][1] + q[2][2]);
                    row[9] = solution.position.num_satellites.to_string();
                    row[10] = format!("{:.17e}", report.evidence.summary.independent_clock_gdop);
                    row[11] = report.evidence.summary.outlier_rejections.to_string();
                    row[12] = report.evidence.summary.raim_rejections.to_string();
                    let velocity_config = VelocityEpochConfig {
                        require_carrier_lock: true,
                        velocity: gnss_rust::velocity::VelocityConfig {
                            clock_model: VelocityClockModel::PerSystem,
                            ..Default::default()
                        },
                        ..Default::default()
                    };
                    match solve_gps_galileo_e1_epoch(
                        &input,
                        &nav.navigation,
                        position,
                        velocity_config,
                    ) {
                        Ok(v) => {
                            eprintln!("profile={profile} Doppler={v:?}");
                            match v.solution {
                                Ok(v) => {
                                    for (i, x) in v.velocity_mps.iter().enumerate() {
                                        row[13 + i] = format!("{x:.17e}");
                                    }
                                    row[16] = format!(
                                        "{:.17e}",
                                        v.covariance_m2s2[0][0]
                                            + v.covariance_m2s2[1][1]
                                            + v.covariance_m2s2[2][2]
                                    );
                                    row[17] = v.satellites.len().to_string();
                                    row[18] = format!("{:.17e}", v.residual_rms_mps);
                                }
                                Err(e) => row[19] = format!("Doppler:{e:?}"),
                            }
                        }
                        Err(e) => row[19] = format!("Doppler preparation:{e:?}"),
                    }
                }
                Err(e) => {
                    row[3] = "INVALID".into();
                    row[19] = format!("SPP:{e:?}");
                }
            }
            println!(
                "{}",
                row.iter()
                    .map(|v| format!("\"{}\"", v.replace('"', "\"\"")))
                    .collect::<Vec<_>>()
                    .join(",")
            );
        }
    }
    Ok(())
}
