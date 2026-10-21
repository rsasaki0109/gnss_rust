//! Raw, read-only L1/E1 phase increment versus current Doppler at one epoch.
use gnss_rust::constants::GPS_L1_WAVELENGTH;
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::{GnssSystem, SatelliteId};
use std::{
    collections::BTreeMap,
    error::Error,
    fs::File,
    io::{BufReader, BufWriter, Write},
};

fn number(value: Option<f64>) -> String {
    value.map_or_else(String::new, |v| format!("{v:.17e}"))
}
fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 3 {
        return Err("usage: ppc_phase_doppler_probe ROVER_OBS GPS_WEEK GPS_TOW".into());
    }
    let week: i32 = args[1].parse()?;
    let tow: f64 = args[2].parse()?;
    let reader = RinexObservationReader::new(BufReader::new(File::open(&args[0])?))?;
    let mut previous = BTreeMap::<(SatelliteId, String), (f64, Option<f64>, Option<u8>, u8)>::new();
    let mut previous_time = None;
    let mut out = BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "week,tow,flag,satellite,tracking,phase_cycles,doppler_hz,snr_db_hz,lli,previous_flag,previous_phase_cycles,previous_doppler_hz,previous_lli,dt_s,raw_phase_minus_Doppler_m"
    )?;
    for epoch in reader {
        let epoch = epoch?;
        let target = epoch.time.week() == week && (epoch.time.tow() - tow).abs() < 1e-7;
        let dt = previous_time.map(|t| epoch.time.difference_seconds(t));
        let mut current = BTreeMap::new();
        for raw in &epoch.observations {
            if raw.satellite.system() != GnssSystem::Gps
                && raw.satellite.system() != GnssSystem::Galileo
            {
                continue;
            }
            if raw.tracking_code != "1C"
                && (raw.satellite.system() != GnssSystem::Galileo || raw.tracking_code != "1X")
            {
                continue;
            }
            let key = (raw.satellite, raw.tracking_code.clone());
            let old = previous.get(&key);
            if target {
                let mut row = vec![String::new(); 15];
                row[0] = week.to_string();
                row[1] = format!("{:.17}", epoch.time.tow());
                row[2] = epoch.flag.to_string();
                row[3] = raw.satellite.to_string();
                row[4] = raw.tracking_code.clone();
                row[5] = number(raw.carrier_phase_cycles);
                row[6] = number(raw.doppler_hz);
                row[7] = number(raw.snr_db_hz);
                row[8] = raw.lli.map_or_else(String::new, |v| v.to_string());
                if let Some((phase, doppler, lli, flag)) = old {
                    row[9] = flag.to_string();
                    row[10] = number(Some(*phase));
                    row[11] = number(*doppler);
                    row[12] = lli.map_or_else(String::new, |v| v.to_string());
                    if let (Some(current_phase), Some(doppler), Some(dt)) =
                        (raw.carrier_phase_cycles, raw.doppler_hz, dt)
                    {
                        row[14] = number(Some(
                            ((current_phase - phase) + doppler * dt) * GPS_L1_WAVELENGTH,
                        ));
                    }
                }
                row[13] = number(dt);
                writeln!(out, "{}", row.join(","))?;
            }
            if let Some(phase) = raw.carrier_phase_cycles {
                current.insert(key, (phase, raw.doppler_hz, raw.lli, epoch.flag));
            }
        }
        if target {
            out.flush()?;
            return Ok(());
        }
        previous = current;
        previous_time = Some(epoch.time);
    }
    Err("requested epoch not found".into())
}
