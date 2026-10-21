//! Read-only single-epoch GPS versus GPS/QZSS code SPP diagnostics; no truth input.
use gnss_rust::geometry::dilution_of_precision;
use gnss_rust::rinex::{RinexObservationReader, read_navigation_for_systems};
use gnss_rust::spp::{SppConfig, solve_epoch};
use gnss_rust::{GnssSystem, GnssTime, constants};
use std::{collections::BTreeMap, error::Error, fs::File, io::BufReader, path::Path};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 3 {
        return Err("usage: ppc_spp_diagnostics CASE_DIR GPS_WEEK GPS_TOW".into());
    }
    let time = GnssTime::new(args[1].parse()?, args[2].parse()?)?;
    let root = Path::new(&args[0]);
    let open = |name| File::open(root.join(name)).map(BufReader::new);
    let mut rover = RinexObservationReader::new(open("rover.obs")?)?;
    let seed = rover
        .header()
        .approximate_position
        .ok_or("no header seed")?;
    let epoch = rover
        .find_map(|r| match r {
            Ok(r) if r.time == time => Some(Ok(r)),
            Ok(_) => None,
            Err(e) => Some(Err(e)),
        })
        .ok_or("target rover epoch absent")??;
    let mut base = RinexObservationReader::new(open("base.obs")?)?;
    let base = base
        .find_map(|r| match r {
            Ok(r) if r.time == time => Some(Ok(r)),
            Ok(_) => None,
            Err(e) => Some(Err(e)),
        })
        .ok_or("target base epoch absent")??;
    let nav = read_navigation_for_systems(open("base.nav")?, &[GnssSystem::Gps, GnssSystem::Qzss])?;
    println!(
        "kind,profile,week,tow,satellite,x_m,y_m,z_m,satellites,code_rms_m,gdop,position_std_norm_m,code_m,phase_cycles,doppler_hz,snr_db_hz,lli,used,code_residual_m,common_base_code_phase,error"
    );
    for qzss in [false, true] {
        let profile = if qzss { "gps_qzss" } else { "gps" };
        let mut input = epoch.clone();
        input.approximate_position = Some(seed);
        input.observations.retain(|o| {
            o.tracking_code == "1C"
                && (o.satellite.system() == GnssSystem::Gps
                    || qzss && o.satellite.system() == GnssSystem::Qzss)
        });
        let mut residuals = BTreeMap::new();
        let mut rejections = BTreeMap::new();
        let mut summary = vec![String::new(); 21];
        summary[0] = "summary".into();
        summary[1] = profile.into();
        summary[2] = time.week().to_string();
        summary[3] = format!("{:.17e}", time.tow());
        match solve_epoch(&input, &nav.navigation, SppConfig::default()) {
            Ok(solution) => {
                let position = solution
                    .position
                    .position_ecef
                    .ok_or("SPP position absent")?;
                for (i, v) in [position.x, position.y, position.z].into_iter().enumerate() {
                    summary[5 + i] = format!("{v:.17e}");
                }
                summary[8] = solution.position.num_satellites.to_string();
                summary[9] = format!("{:.17e}", solution.residual_rms_m);
                let q = solution
                    .position
                    .position_covariance
                    .ok_or("SPP covariance absent")?;
                summary[11] = format!("{:.17e}", (0..3).map(|i| q[i][i]).sum::<f64>().sqrt());
                residuals.extend(solution.residuals_m);
                rejections.extend(solution.rejected);
                let mut positions = Vec::new();
                for o in &input.observations {
                    if !residuals.contains_key(&o.satellite) {
                        continue;
                    }
                    let tx = time.checked_add_seconds(
                        -o.pseudorange_m.ok_or("used code absent")? / constants::SPEED_OF_LIGHT,
                    )?;
                    let eph = nav.navigation.ephemeris(
                        o.satellite,
                        tx,
                        SppConfig::default().max_ephemeris_age_s,
                    )?;
                    let (_, clock) = eph.position_clock(tx)?;
                    let (p, _) = eph.position_clock(tx.checked_add_seconds(-clock)?)?;
                    positions.push((o.satellite, p));
                }
                match dilution_of_precision(position, &positions) {
                    Ok(d) => summary[10] = format!("{:.17e}", d.gdop),
                    Err(e) => summary[20] = format!("DOP:{e:?}"),
                }
            }
            Err(error) => summary[20] = format!("{error:?}"),
        }
        println!("{}", summary.join(","));
        for o in &input.observations {
            let mut row = vec![String::new(); 21];
            row[0] = "observation".into();
            row[1] = profile.into();
            row[2] = time.week().to_string();
            row[3] = format!("{:.17e}", time.tow());
            row[4] = o.satellite.to_string();
            for (i, value) in [
                o.pseudorange_m,
                o.carrier_phase_cycles,
                o.doppler_hz,
                o.snr_db_hz,
            ]
            .into_iter()
            .enumerate()
            {
                if let Some(v) = value {
                    row[12 + i] = format!("{v:.17e}");
                }
            }
            if let Some(lli) = o.lli {
                row[16] = lli.to_string();
            }
            row[17] = u8::from(residuals.contains_key(&o.satellite)).to_string();
            if let Some(residual) = residuals.get(&o.satellite) {
                row[18] = format!("{residual:.17e}");
            }
            row[19] = u8::from(
                base.observations.iter().any(|b| {
                    b.satellite == o.satellite
                        && b.tracking_code == "1C"
                        && b.pseudorange_m.is_some()
                        && b.carrier_phase_cycles.is_some()
                }) && o.pseudorange_m.is_some()
                    && o.carrier_phase_cycles.is_some(),
            )
            .to_string();
            if let Some(reason) = rejections.get(&o.satellite) {
                row[20] = format!("{reason:?}");
            }
            println!("{}", row.join(","));
        }
    }
    Ok(())
}
