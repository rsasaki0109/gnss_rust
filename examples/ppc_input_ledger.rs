//! Canonical matched GPS 1C input ledger for independent decoder comparison.
use gnss_rust::epoch_match::{EpochMatch, EpochMatcher};
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::{GnssSystem, observation::ObservationEpoch};
use std::{
    error::Error,
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

fn broken(epoch: &ObservationEpoch) -> bool {
    matches!(epoch.flag, 1 | 6)
        || epoch.observations.iter().any(|o| {
            o.satellite.system() == GnssSystem::Gps
                && o.tracking_code == "1C"
                && (o.loss_of_lock() || o.half_cycle_ambiguity())
        })
}
fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() != 1 {
        return Err("usage: ppc_input_ledger CASE_DIRECTORY".into());
    }
    let root = Path::new(&args[0]);
    let open = |name| File::open(root.join(name)).map(BufReader::new);
    let rover = RinexObservationReader::new(open("rover.obs")?)?;
    let base = RinexObservationReader::new(open("base.obs")?)?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "side,week,tow,satellite,has_code,code_m,has_phase,phase_cycles,has_doppler,doppler_hz,snr_db_hz,lli,derived_reset"
    )?;
    let mut reset = false;
    for event in EpochMatcher::new(rover, base) {
        match event? {
            EpochMatch::Unmatched { epoch, .. } => reset |= broken(&epoch),
            EpochMatch::Pair { rover, base } => {
                for (side, epoch) in [("rover", &rover), ("base", &base)] {
                    if epoch.flag != 0 {
                        return Err("comparison ledger currently requires normal PPC epochs".into());
                    }
                    let mut observations: Vec<_> = epoch
                        .observations
                        .iter()
                        .filter(|o| {
                            o.satellite.system() == GnssSystem::Gps && o.tracking_code == "1C"
                        })
                        .filter(|o| {
                            o.pseudorange_m.is_some()
                                || o.carrier_phase_cycles.is_some()
                                || o.doppler_hz.is_some()
                                || o.snr_db_hz.unwrap_or(0.0) != 0.0
                                || o.lli.unwrap_or(0) != 0
                        })
                        .collect();
                    observations.sort_by_key(|o| o.satellite);
                    for o in observations {
                        writeln!(
                            out,
                            "{side},{},{:.17e},{},{},{:.17e},{},{:.17e},{},{:.17e},{:.17e},{},{}",
                            epoch.time.week(),
                            epoch.time.tow(),
                            o.satellite,
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
