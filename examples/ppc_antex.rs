//! Strict GPS satellite calibration admission; no antenna/clock compatibility inferred.
use gnss_rust::antenna::read_antex;
use gnss_rust::{GnssTime, SatelliteId};
use std::{
    fs::File,
    io::{BufReader, BufWriter, Write},
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 3 {
        return Err("usage: ppc_antex ANTEX_FILE GPS_WEEK GPS_TOW".into());
    }
    let time = GnssTime::new(args[1].parse()?, args[2].parse()?)?;
    let atx = read_antex(BufReader::new(File::open(&args[0])?))?;
    let mut out = BufWriter::new(std::io::stdout().lock());
    writeln!(
        out,
        "week,tow,satellite,status,svn,type,source_line,g01_n_m,g01_e_m,g01_u_m,g02_n_m,g02_e_m,g02_u_m,error"
    )?;
    for prn in 1..=32 {
        let sat: SatelliteId = format!("G{prn:02}").parse()?;
        let found: Vec<_> = atx
            .entries()
            .iter()
            .filter(|c| c.satellite() == Some(sat) && c.contains(time))
            .collect();
        let mut row = vec![String::new(); 14];
        row[0] = time.week().to_string();
        row[1] = format!("{:.17}", time.tow());
        row[2] = sat.to_string();
        if found.len() == 1 {
            let c = found[0];
            atx.satellite_calibration(sat, &c.svn, time)?;
            if let (Some(g1), Some(g2)) = (c.frequencies.get("G01"), c.frequencies.get("G02")) {
                row[3] = "AVAILABLE".into();
                row[4] = c.svn.clone();
                row[5] = c.antenna_type.clone();
                row[6] = c.source_line.to_string();
                for (i, v) in g1
                    .offset_neu_m
                    .into_iter()
                    .chain(g2.offset_neu_m)
                    .enumerate()
                {
                    row[7 + i] = format!("{v:.17e}");
                }
            } else {
                row[3] = "MISSING_FREQUENCY".into();
            }
        } else {
            row[3] = if found.is_empty() {
                "MISSING"
            } else {
                "AMBIGUOUS"
            }
            .into();
        }
        writeln!(out, "{}", row.join(","))?;
    }
    out.flush()?;
    eprintln!(
        "parsed_calibrations={} selection=unique_PRN_and_validity_explicitly_records_SVN receiver_calibration_known=false clock_compatibility_verified=false satellite_PCV_applied=false PPP_executed=false",
        atx.entries().len()
    );
    Ok(())
}
