//! Independent decoded GPS 1C/Galileo E1 selected-input ledger, without a solver.
use gnss_rust::epoch_match::{EpochMatch, EpochMatcher};
use gnss_rust::observation::{Observation, ObservationEpoch};
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::{GnssSystem, SatelliteId};
use std::{
    collections::BTreeMap,
    error::Error,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

fn eligible(o: &Observation) -> bool {
    o.satellite.system() == GnssSystem::Gps && o.tracking_code == "1C"
        || o.satellite.system() == GnssSystem::Galileo
            && matches!(o.tracking_code.as_str(), "1C" | "1X")
}
fn selected(epoch: &ObservationEpoch) -> BTreeMap<SatelliteId, &Observation> {
    let mut out = BTreeMap::<SatelliteId, &Observation>::new();
    let rank = |o: &Observation| {
        (
            !o.pseudorange_m.is_some_and(|v| v.is_finite() && v > 0.0)
                || !o.carrier_phase_cycles.is_some_and(f64::is_finite),
            o.tracking_code != "1C",
        )
    };
    for o in &epoch.observations {
        if eligible(o) && out.get(&o.satellite).is_none_or(|old| rank(o) < rank(old)) {
            out.insert(o.satellite, o);
        }
    }
    out
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() != 1 {
        return Err("usage: ppc_mixed_input_ledger CASE_DIRECTORY".into());
    }
    let root = Path::new(&args[0]);
    let open = |n| File::open(root.join(n)).map(BufReader::new);
    let rover = RinexObservationReader::new(open("rover.obs")?)?;
    let base = RinexObservationReader::new(open("base.obs")?)?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "side,week,tow,satellite,tracking_code,has_code,code_m,has_phase,phase_cycles,has_doppler,doppler_hz,snr_db_hz,lli,derived_reset"
    )?;
    let mut reset = false;
    for event in EpochMatcher::new(rover, base) {
        match event? {
            EpochMatch::Unmatched { epoch, .. } => {
                reset |= matches!(epoch.flag, 1 | 6)
                    || epoch
                        .observations
                        .iter()
                        .any(|o| eligible(o) && (o.loss_of_lock() || o.half_cycle_ambiguity()))
            }
            EpochMatch::Pair { rover, base } => {
                for (side, epoch) in [("rover", &rover), ("base", &base)] {
                    if epoch.flag != 0 {
                        return Err("ledger requires normal PPC epochs".into());
                    }
                    for o in selected(epoch).into_values() {
                        if o.pseudorange_m.is_none()
                            && o.carrier_phase_cycles.is_none()
                            && o.doppler_hz.is_none()
                            && o.snr_db_hz.unwrap_or(0.0) == 0.0
                            && o.lli.unwrap_or(0) == 0
                        {
                            continue;
                        }
                        writeln!(
                            out,
                            "{side},{},{:.17e},{},{},{},{:.17e},{},{:.17e},{},{:.17e},{:.17e},{},{}",
                            epoch.time.week(),
                            epoch.time.tow(),
                            o.satellite,
                            o.tracking_code,
                            u8::from(o.pseudorange_m.is_some()),
                            o.pseudorange_m.unwrap_or(0.0),
                            u8::from(o.carrier_phase_cycles.is_some()),
                            o.carrier_phase_cycles.unwrap_or(0.0),
                            u8::from(o.doppler_hz.is_some()),
                            o.doppler_hz.unwrap_or(0.0),
                            o.snr_db_hz.unwrap_or(0.0),
                            o.lli.unwrap_or(0),
                            u8::from(reset)
                        )?;
                    }
                }
                reset = false;
            }
        }
    }
    out.flush()?;
    Ok(())
}
