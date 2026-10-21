//! GPS 1C velocity diagnostic at the fixed, surveyed PPC base station.
use gnss_rust::rinex::{RinexObservationReader, read_navigation_for_systems};
use gnss_rust::velocity::{VelocityEpochConfig, solve_gps_l1_epoch};
use gnss_rust::{EcefCoord, GnssSystem};
use std::{
    error::Error,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 2 {
        return Err("usage: ppc_base_velocity CASE_DIRECTORY BASE_X,BASE_Y,BASE_Z".into());
    }
    let xyz: Vec<f64> = args[1]
        .split(',')
        .map(str::parse)
        .collect::<Result<_, _>>()?;
    if xyz.len() != 3 {
        return Err("expected three base ECEF components".into());
    }
    let position = EcefCoord::new(xyz[0], xyz[1], xyz[2]);
    let root = Path::new(&args[0]);
    let open = |name| File::open(root.join(name)).map(BufReader::new);
    let reader = RinexObservationReader::new(open("base.obs")?)?;
    let nav = read_navigation_for_systems(open("base.nav")?, &[GnssSystem::Gps])?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "week,tow,flag,status,prepared_count,rejected_count,zero_doppler_count,missing_code_count,vx_mps,vy_mps,vz_mps,speed_mps,clock_drift_sps,residual_rms_mps,chi_square_per_dof,covariance_scale,velocity_covariance_trace_m2s2,max_unscaled_velocity_variance_m2s2,used_satellites,rejections,reason"
    )?;
    for epoch in reader {
        let epoch = epoch?;
        let mut row = vec![String::new(); 21];
        row[0] = epoch.time.week().to_string();
        row[1] = format!("{:.17}", epoch.time.tow());
        row[2] = epoch.flag.to_string();
        match solve_gps_l1_epoch(
            &epoch,
            &nav.navigation,
            position,
            VelocityEpochConfig::default(),
        ) {
            Err(error) => {
                row[3] = "INVALID_EPOCH".into();
                row[20] = format!("{error:?}");
            }
            Ok(report) => {
                row[4] = report.prepared.len().to_string();
                row[5] = report.rejected.len().to_string();
                row[6] = report
                    .prepared
                    .iter()
                    .filter(|p| p.measurement.doppler_hz == 0.0)
                    .count()
                    .to_string();
                row[7] = report
                    .prepared
                    .iter()
                    .filter(|p| p.used_missing_code_travel_time)
                    .count()
                    .to_string();
                row[19] = report
                    .rejected
                    .iter()
                    .map(|(raw, reason)| {
                        format!("{}:{}:{reason:?}", raw.satellite, raw.tracking_code)
                    })
                    .collect::<Vec<_>>()
                    .join(";");
                match report.solution {
                    Err(error) => {
                        row[3] = "INVALID_VELOCITY".into();
                        row[20] = format!("{error:?}");
                    }
                    Ok(solution) => {
                        row[3] = "VELOCITY".into();
                        for (i, v) in solution.velocity_mps.into_iter().enumerate() {
                            row[8 + i] = format!("{v:.17e}");
                        }
                        let v = solution.velocity_mps;
                        row[11] = format!("{:.17e}", v[0].hypot(v[1]).hypot(v[2]));
                        for (i, value) in [
                            solution.receiver_clock_drift_sps,
                            solution.residual_rms_mps,
                            solution.chi_square_per_dof,
                            solution.covariance_scale,
                            (0..3).map(|i| solution.covariance_m2s2[i][i]).sum::<f64>(),
                            (0..3)
                                .map(|i| solution.unscaled_covariance_m2s2[i][i])
                                .fold(0.0, f64::max),
                        ]
                        .into_iter()
                        .enumerate()
                        {
                            row[12 + i] = format!("{value:.17e}");
                        }
                        row[18] = solution
                            .satellites
                            .iter()
                            .map(ToString::to_string)
                            .collect::<Vec<_>>()
                            .join(";");
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
        "GPS_1C_only=true base_stationary_diagnostic=true velocity_not_fused=true no_FIX_policy=true surveyed_base_position=true receiver_clock_offsets_not_applied=true excluded_navigation_records={}",
        nav.excluded_records.len()
    );
    Ok(())
}
