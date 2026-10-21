//! Exact decoded GPS 1C/2W input ledger and explicit native OSB loader view.
use gnss_rust::bias::read_satellite_biases_with_byte_comments;
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::{GnssSystem, GnssTime, SatelliteId};
use std::{
    fs::File,
    io::{BufReader, BufWriter, Write},
};

fn optional(value: Option<f64>) -> String {
    value.map_or_else(String::new, |v| format!("{v:.17e}"))
}
fn observations(path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let reader = RinexObservationReader::new(BufReader::new(File::open(path)?))?;
    eprintln!(
        "RINEX_header={:?} selection=GPS_1C_2W_only raw_values_unchanged=true source_options_retained=true",
        reader.header()
    );
    let mut out = BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "kind,week,tow,satellite,tracking,code_m,phase_cycles,doppler_hz,snr_db_hz,lli,ssi,flag,receiver_offset_s,header_x_m,header_y_m,header_z_m"
    )?;
    let (mut count, mut groups) = (0, 0);
    for epoch in reader {
        let epoch = epoch?;
        let xyz = epoch
            .approximate_position
            .ok_or("missing RINEX approximate position")?;
        let mut row = vec![String::new(); 16];
        row[0] = "epoch".into();
        row[1] = epoch.time.week().to_string();
        row[2] = format!("{:.17}", epoch.time.tow());
        row[11] = epoch.flag.to_string();
        row[12] = optional(epoch.receiver_clock_offset_s);
        for (i, v) in [xyz.x, xyz.y, xyz.z].into_iter().enumerate() {
            row[13 + i] = format!("{v:.17e}");
        }
        writeln!(out, "{}", row.join(","))?;
        count += 1;
        for o in epoch.observations.iter().filter(|o| {
            o.satellite.system() == GnssSystem::Gps
                && o.satellite.prn() <= 32
                && matches!(o.tracking_code.as_str(), "1C" | "2W")
        }) {
            let mut row = vec![String::new(); 16];
            row[0] = "observation".into();
            row[1] = epoch.time.week().to_string();
            row[2] = format!("{:.17}", epoch.time.tow());
            row[3] = o.satellite.to_string();
            row[4] = o.tracking_code.clone();
            row[5] = optional(o.pseudorange_m);
            row[6] = optional(o.carrier_phase_cycles);
            row[7] = optional(o.doppler_hz);
            row[8] = optional(o.snr_db_hz);
            row[9] = o.lli.map_or_else(String::new, |v| v.to_string());
            row[10] = o
                .signal_strength
                .map_or_else(String::new, |v| v.to_string());
            writeln!(out, "{}", row.join(","))?;
            groups += 1;
        }
    }
    out.flush()?;
    eprintln!(
        "raw_epochs={count} selected_groups={groups} receiver_clock_unit=seconds phase_unit=cycles code_unit=metres solver_executed=false"
    );
    Ok(())
}
fn bias_view(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    let selected = read_satellite_biases_with_byte_comments(
        BufReader::new(File::open(&args[1])?),
        Default::default(),
    )?;
    let p = &selected.selection.products;
    let week = args[2].parse()?;
    let start = GnssTime::new(week, args[3].parse()?)?;
    let end = GnssTime::new(week, args[4].parse()?)?;
    if end < start
        || !p
            .description()
            .get("SATELLITE_CLOCK_REFERENCE_OBSERVABLES")
            .is_some_and(|v| v.iter().any(|x| x.as_slice() == ["G", "C1W", "C2W"]))
        || !p
            .description()
            .get("BIAS_MODE")
            .is_some_and(|v| v.len() == 1 && v[0].as_slice() == ["ABSOLUTE"])
    {
        return Err("requires absolute GPS IF(C1W,C2W) product and ordered case interval".into());
    }
    let mut rows = Vec::new();
    for prn in 1..=32 {
        let sat: SatelliteId = format!("G{prn:02}").parse()?;
        for code in ["C1C", "C2W"] {
            let a = p.select_osb(sat, code, start)?;
            let b = p.select_osb(sat, code, end)?;
            if a.line != b.line || a.observation_2.is_some() {
                return Err("OSB changes across the requested case interval".into());
            }
            rows.push(a.clone());
        }
    }
    let mut destination = BufWriter::new(
        std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&args[5])?,
    );
    let mut header: Vec<_> = p.header().split_whitespace().map(str::to_string).collect();
    if header.len() != 9 {
        return Err("unexpected Bias-SINEX header".into());
    }
    header[8] = format!("{:08}", rows.len());
    writeln!(destination, "{}", header.join(" "))?;
    writeln!(
        destination,
        "* Explicit GPS code-only selection; native loader requires repeated OBS2 for OSB."
    )?;
    writeln!(destination, "+BIAS/DESCRIPTION")?;
    for (key, values) in p.description() {
        for value in values {
            writeln!(destination, " {key} {}", value.join(" "))?;
        }
    }
    writeln!(destination, "-BIAS/DESCRIPTION\n+BIAS/SOLUTION")?;
    println!(
        "source_line,svn,satellite,observation,source_OBS2_present,unit,value,sigma,metres,code_add_m,start,end"
    );
    for entry in rows {
        let svn = entry.svn.as_deref().ok_or("missing SVN")?;
        let (metres, _) = entry.metres(None)?;
        writeln!(
            destination,
            " OSB {svn} {} {} {} {} {} {} {:.17e} {:.17e}",
            entry.satellite,
            entry.observation_1,
            entry.observation_1,
            entry.start_label,
            entry.end_label,
            entry.unit_label,
            entry.value,
            entry.sigma
        )?;
        println!(
            "{},{svn},{},{},false,{},{:.17e},{:.17e},{metres:.17e},{:.17e},{},{}",
            entry.line,
            entry.satellite,
            entry.observation_1,
            entry.unit_label,
            entry.value,
            entry.sigma,
            -metres,
            entry.start_label,
            entry.end_label
        );
    }
    writeln!(destination, "-BIAS/SOLUTION\n%=ENDBIA")?;
    destination.flush()?;
    eprintln!(
        "selected_GPS_code_OSBs=64 repeated_OBS2_only_in_compatibility_view=true phase_OSBs_selected=false synthetic_biases_created=false source_values_and_sigmas_unchanged=true full_case_half_open_validity_checked=true preserved_byte_comment_lines={} excluded_receiver_rows={}",
        selected.preserved_free_text.len(),
        selected.selection.excluded_receiver_records.len()
    );
    Ok(())
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() == 2 && args[0] == "observations" {
        observations(&args[1])
    } else if args.len() == 6 && args[0] == "native-bias" {
        bias_view(&args)
    } else {
        Err("usage: ppc_ppp_inputs observations ROVER_OBS | native-bias CODE_BIA GPS_WEEK FIRST_TOW LAST_TOW NEW_COMPAT_BIA".into())
    }
}
