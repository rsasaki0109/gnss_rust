//! Read-only diagnostic of raw SD ambiguity datums before successful updates.
//! A reset arc never contributes to the retained-arc offset. No truth is read.
use gnss_rust::constants::GPS_L1_WAVELENGTH;
use gnss_rust::epoch_match::{EpochMatch, EpochMatcher, EpochSide};
use gnss_rust::rinex::{RinexObservationReader, read_navigation_for_systems};
use gnss_rust::rtk::{RtkConstellations, RtkFloatConfig, RtkFloatFilter, RtkPositionMode};
use gnss_rust::{EcefCoord, GnssSystem};
use std::{
    collections::BTreeMap,
    error::Error,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 3 {
        return Err("usage: ppc_rtk_gauge CASE_DIRECTORY gps|galileo BASE_X,BASE_Y,BASE_Z".into());
    }
    let constellations = match args[1].as_str() {
        "gps" => RtkConstellations::Gps,
        "galileo" => RtkConstellations::GpsGalileo,
        _ => return Err("expected gps or galileo".into()),
    };
    let systems = if constellations == RtkConstellations::Gps {
        vec![GnssSystem::Gps]
    } else {
        vec![GnssSystem::Gps, GnssSystem::Galileo]
    };
    let xyz: Vec<f64> = args[2]
        .split(',')
        .map(str::parse)
        .collect::<Result<_, _>>()?;
    if xyz.len() != 3 {
        return Err("expected three base ECEF components".into());
    }
    let root = Path::new(&args[0]);
    let open = |name| File::open(root.join(name)).map(BufReader::new);
    let rover = RinexObservationReader::new(open("rover.obs")?)?;
    let seed = rover
        .header()
        .approximate_position
        .ok_or("missing rover position")?;
    let base = RinexObservationReader::new(open("base.obs")?)?;
    let nav = read_navigation_for_systems(open("base.nav")?, &systems)?;
    let mut filter = RtkFloatFilter::new(
        EcefCoord::new(xyz[0], xyz[1], xyz[2]),
        seed,
        RtkFloatConfig {
            constellations,
            position_mode: RtkPositionMode::Kinematic,
            ..Default::default()
        },
    )?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "kind,week,tow,side,status,system,retained_arcs,reset_arcs,mean_offset_cycles,mean_offset_m,retained_raw_minus_state_cycles,reset_satellites,x_m,y_m,z_m,reason"
    )?;
    for event in EpochMatcher::new(rover, base) {
        let mut row = vec![String::new(); 16];
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
                    row[15] = "continuity_reset_pending".into();
                }
                writeln!(out, "{}", row.join(","))?;
            }
            EpochMatch::Pair { rover, base } => {
                row[0] = "pair".into();
                row[1] = rover.time.week().to_string();
                row[2] = format!("{:.17}", rover.time.tow());
                row[3] = "both".into();
                let old: BTreeMap<_, _> = filter
                    .state()
                    .map(|state| {
                        filter
                            .tracked_satellites()
                            .iter()
                            .enumerate()
                            .map(|(i, sat)| (*sat, state.mean[3 + i]))
                            .collect()
                    })
                    .unwrap_or_default();
                if rover.flag == 6 || base.flag == 6 {
                    filter.invalidate_ambiguities();
                    row[4] = "EVENT".into();
                    row[15] = "cycle_slip_record".into();
                    writeln!(out, "{}", row.join(","))?;
                    continue;
                }
                match filter.process_pair(&rover, &base, &nav.navigation) {
                    Err(error) => {
                        row[4] = "INVALID".into();
                        row[15] = format!("{error:?}").replace(',', ";");
                        writeln!(out, "{}", row.join(","))?;
                    }
                    Ok(result) => {
                        row[4] = "FLOAT".into();
                        let position = result
                            .position
                            .position_ecef
                            .ok_or("missing FLOAT position")?;
                        for (i, value) in
                            [position.x, position.y, position.z].into_iter().enumerate()
                        {
                            row[12 + i] = format!("{value:.9}");
                        }
                        let pairs = result.tracking_pairs();
                        for system in &systems {
                            let mut offsets = Vec::new();
                            let mut resets = Vec::new();
                            for pair in pairs.iter().filter(|p| p.satellite.system() == *system) {
                                if result.reset_ambiguities.contains(&pair.satellite) {
                                    resets.push(pair.satellite.to_string());
                                    continue;
                                }
                                let ro = rover
                                    .observations
                                    .iter()
                                    .find(|o| {
                                        o.satellite == pair.satellite
                                            && o.tracking_code == pair.rover_tracking_code
                                    })
                                    .ok_or("missing selected rover observation")?;
                                let ba = base
                                    .observations
                                    .iter()
                                    .find(|o| {
                                        o.satellite == pair.satellite
                                            && o.tracking_code == pair.base_tracking_code
                                    })
                                    .ok_or("missing selected base observation")?;
                                let raw = (ro.carrier_phase_cycles.ok_or("missing rover phase")?
                                    - ba.carrier_phase_cycles.ok_or("missing base phase")?)
                                    - (ro.pseudorange_m.ok_or("missing rover code")?
                                        - ba.pseudorange_m.ok_or("missing base code")?)
                                        / GPS_L1_WAVELENGTH;
                                offsets.push((
                                    pair.satellite,
                                    raw - old
                                        .get(&pair.satellite)
                                        .ok_or("retained arc without old state")?,
                                ));
                            }
                            if offsets.is_empty() && resets.is_empty() {
                                continue;
                            }
                            row[5] = system.rinex_prefix().to_string();
                            row[6] = offsets.len().to_string();
                            row[7] = resets.len().to_string();
                            row[8].clear();
                            row[9].clear();
                            if !offsets.is_empty() {
                                let offset = offsets.iter().map(|(_, v)| v).sum::<f64>()
                                    / offsets.len() as f64;
                                row[8] = format!("{offset:.17e}");
                                row[9] = format!("{:.17e}", offset * GPS_L1_WAVELENGTH);
                            }
                            row[10] = offsets
                                .iter()
                                .map(|(s, v)| format!("{s}:{v:.17e}"))
                                .collect::<Vec<_>>()
                                .join(";");
                            row[11] = resets.join(";");
                            writeln!(out, "{}", row.join(","))?;
                        }
                    }
                }
            }
        }
    }
    out.flush()?;
    Ok(())
}
