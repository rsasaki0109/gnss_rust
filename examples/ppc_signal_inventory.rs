//! Read-only GPS/Galileo code/phase and healthy-navigation availability at common epochs.
use gnss_rust::epoch_match::{EpochMatch, EpochMatcher};
use gnss_rust::observation::{Observation, ObservationEpoch};
use gnss_rust::rinex::{RinexObservationReader, read_navigation_for_systems};
use gnss_rust::{GnssSystem, SatelliteId, constants};
use std::{
    collections::BTreeMap,
    error::Error,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

fn selected(epoch: &ObservationEpoch, system: GnssSystem) -> BTreeMap<SatelliteId, &Observation> {
    let mut result = BTreeMap::<SatelliteId, &Observation>::new();
    for o in &epoch.observations {
        let allowed =
            o.tracking_code == "1C" || system == GnssSystem::Galileo && o.tracking_code == "1X";
        if o.satellite.system() != system
            || !allowed
            || o.half_cycle_ambiguity()
            || !o.pseudorange_m.is_some_and(|v| v.is_finite() && v > 0.0)
            || !o.carrier_phase_cycles.is_some_and(f64::is_finite)
        {
            continue;
        }
        if result
            .get(&o.satellite)
            .is_none_or(|old| o.tracking_code < old.tracking_code)
        {
            result.insert(o.satellite, o);
        }
    }
    result
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() != 1 {
        return Err("usage: ppc_signal_inventory CASE_DIRECTORY".into());
    }
    let root = Path::new(&args[0]);
    let open = |n| File::open(root.join(n)).map(BufReader::new);
    let rover = RinexObservationReader::new(open("rover.obs")?)?;
    let base = RinexObservationReader::new(open("base.obs")?)?;
    let nav =
        read_navigation_for_systems(open("base.nav")?, &[GnssSystem::Gps, GnssSystem::Galileo])?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "week,tow,gps_complete,gps_healthy_nav,galileo_complete,galileo_healthy_nav,pair_tracking_codes"
    )?;
    for event in EpochMatcher::new(rover, base) {
        let EpochMatch::Pair { rover, base } = event? else {
            continue;
        };
        let mut row = vec![
            rover.time.week().to_string(),
            format!("{:.17e}", rover.time.tow()),
        ];
        let mut tracking = Vec::new();
        for system in [GnssSystem::Gps, GnssSystem::Galileo] {
            let r = selected(&rover, system);
            let b = selected(&base, system);
            let (mut complete, mut healthy) = (0, 0);
            for (sat, ro) in r {
                let Some(ba) = b.get(&sat) else {
                    continue;
                };
                complete += 1;
                let rt = rover
                    .time
                    .checked_add_seconds(-ro.pseudorange_m.unwrap() / constants::SPEED_OF_LIGHT)?;
                let bt = base
                    .time
                    .checked_add_seconds(-ba.pseudorange_m.unwrap() / constants::SPEED_OF_LIGHT)?;
                if nav.navigation.ephemeris(sat, rt, 14400.0).is_ok()
                    && nav.navigation.ephemeris(sat, bt, 14400.0).is_ok()
                {
                    healthy += 1;
                    tracking.push(format!("{sat}:{}:{}", ro.tracking_code, ba.tracking_code));
                }
            }
            row.extend([complete.to_string(), healthy.to_string()]);
        }
        row.push(tracking.join(";"));
        writeln!(out, "{}", row.join(","))?;
    }
    eprintln!(
        "navigation_systems=GPS,Galileo excluded_navigation_records={} elevation_mask_not_applied=true no_position_solver=true",
        nav.excluded_records.len()
    );
    Ok(())
}
