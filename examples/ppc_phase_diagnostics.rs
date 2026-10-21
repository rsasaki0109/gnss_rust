//! Read-only raw SD code/phase/doppler diagnostics, not an applied slip policy.
use gnss_rust::constants::GPS_L1_WAVELENGTH as LAMBDA;
use gnss_rust::epoch_match::{EpochMatch, EpochMatcher};
use gnss_rust::observation::Observation;
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::{GnssSystem, GnssTime, SatelliteId};
use std::{collections::BTreeMap, error::Error, fs::File, io::BufReader, path::Path};

fn selected(observations: &[Observation]) -> BTreeMap<SatelliteId, &Observation> {
    observations
        .iter()
        .filter(|o| o.satellite.system() == GnssSystem::Gps && o.tracking_code == "1C")
        .map(|o| (o.satellite, o))
        .collect()
}
fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 3 {
        return Err("usage: ppc_phase_diagnostics CASE_DIRECTORY FROM_TOW TO_TOW".into());
    }
    let root = Path::new(&args[0]);
    let from: f64 = args[1].parse()?;
    let until: f64 = args[2].parse()?;
    let open = |name| File::open(root.join(name)).map(BufReader::new);
    let rover = RinexObservationReader::new(open("rover.obs")?)?;
    let base = RinexObservationReader::new(open("base.obs")?)?;
    let mut previous: BTreeMap<SatelliteId, (GnssTime, f64, f64)> = BTreeMap::new();
    println!(
        "week,tow,satellite,sd_code_m,sd_phase_m,cmc_m,dt_s,cmc_jump_m,doppler_error_m,rover_lli,base_lli,rover_snr,base_snr"
    );
    for event in EpochMatcher::new(rover, base) {
        if let EpochMatch::Pair { rover, base } = event? {
            let (ro, ba) = (selected(&rover.observations), selected(&base.observations));
            for (satellite, r) in ro {
                let Some(b) = ba.get(&satellite) else {
                    continue;
                };
                let (Some(rp), Some(bp), Some(rc), Some(bc)) = (
                    r.pseudorange_m,
                    b.pseudorange_m,
                    r.carrier_phase_cycles,
                    b.carrier_phase_cycles,
                ) else {
                    previous.remove(&satellite);
                    continue;
                };
                let code = rp - bp;
                let phase = (rc - bc) * LAMBDA;
                let cmc = code - phase;
                let mut fields = vec![String::new(); 13];
                fields[0] = rover.time.week().to_string();
                fields[1] = format!("{:.17}", rover.time.tow());
                fields[2] = satellite.to_string();
                for (i, v) in [code, phase, cmc].into_iter().enumerate() {
                    fields[3 + i] = format!("{v:.17e}");
                }
                if let Some((time, previous_cmc, previous_phase)) = previous.get(&satellite) {
                    let dt = rover.time.difference_seconds(*time);
                    fields[6] = format!("{dt:.17e}");
                    fields[7] = format!("{:.17e}", cmc - previous_cmc);
                    if let (Some(rd), Some(bd)) = (r.doppler_hz, b.doppler_hz) {
                        fields[8] =
                            format!("{:.17e}", phase - previous_phase + (rd - bd) * LAMBDA * dt);
                    }
                }
                fields[9] = r.lli.map(|v| v.to_string()).unwrap_or_default();
                fields[10] = b.lli.map(|v| v.to_string()).unwrap_or_default();
                fields[11] = r.snr_db_hz.map(|v| v.to_string()).unwrap_or_default();
                fields[12] = b.snr_db_hz.map(|v| v.to_string()).unwrap_or_default();
                if (from..=until).contains(&rover.time.tow()) {
                    println!("{}", fields.join(","));
                }
                previous.insert(satellite, (rover.time, cmc, phase));
            }
            if rover.time.tow() > until {
                break;
            }
        }
    }
    eprintln!("raw_pair_diagnostics_only=true slip_policy_applied=false");
    Ok(())
}
