//! GPS L1 FLOAT on common epochs; every unmatched/invalid epoch is reported.
use gnss_rust::epoch_match::{EpochMatch, EpochMatcher, EpochSide};
use gnss_rust::rinex::{RinexObservationReader, read_navigation_for_systems};
use gnss_rust::rtk::{RtkFloatConfig, RtkFloatFilter, RtkPositionMode};
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
            "usage: ppc_rtk_float CASE_DIRECTORY static|kinematic BASE_X,BASE_Y,BASE_Z".into(),
        );
    }
    let position_mode = match args[1].to_str() {
        Some("static") => RtkPositionMode::Static,
        Some("kinematic") => RtkPositionMode::Kinematic,
        _ => return Err("expected static or kinematic".into()),
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
            ..Default::default()
        },
    )?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "kind,week,tow,side,rover_flag,base_flag,base_dt_s,status,x_m,y_m,z_m,satellites,reference,nis,code_rms_m,phase_rms_m,ambiguity_resets,reason"
    )?;
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
        let mut row = vec![String::new(); 18];
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
    Ok(())
}
