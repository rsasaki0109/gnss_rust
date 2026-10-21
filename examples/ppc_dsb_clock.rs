//! Inspect explicit GPS IF-clock DSB gauges and raw PPC dual-frequency coverage.
//! No clock compatibility is inferred and no positioning/phase correction runs.
use gnss_rust::bias::{BiasSignalRequest, read_satellite_biases};
use gnss_rust::dsb_clock::{DsbClockError, GpsIfClockBinding, GpsIfDsbClock};
use gnss_rust::precise::ClockSource;
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::{GnssSystem, GnssTime, SatelliteId};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{BufReader, BufWriter, Write},
    path::Path,
};

fn request(satellite: SatelliteId, code: &str) -> BiasSignalRequest {
    BiasSignalRequest {
        satellite,
        tracking_code: code.into(),
        glonass_channel: None,
    }
}
fn error_kind(e: &DsbClockError) -> &'static str {
    match e {
        DsbClockError::MissingPair { .. } => "MISSING",
        DsbClockError::ExpiredPair { .. } => "EXPIRED",
        DsbClockError::AmbiguousPair { .. } => "AMBIGUOUS",
        _ => "ERROR",
    }
}
fn quoted(value: &str) -> String {
    format!("\"{}\"", value.replace('"', "\"\""))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() != 4 && args.len() != 3 {
        return Err(
            "usage: ppc_dsb_clock BIAS_FILE --at GPS_WEEK GPS_TOW | BIAS_FILE --rover ROVER_OBS"
                .into(),
        );
    }
    let selected =
        read_satellite_biases(BufReader::new(File::open(&args[0])?), Default::default())?;
    let products = selected.products;
    let model = GpsIfDsbClock::new(&products)?;
    let binding = GpsIfClockBinding {
        product_id: Path::new(&args[0])
            .file_name()
            .unwrap()
            .to_string_lossy()
            .into_owned(),
        clock_source: ClockSource::RinexClk,
        clock_reference: "UNVERIFIED-caller-IF-C1W-C2W-gauge-for-admission-audit".into(),
    };
    let mut out = BufWriter::new(std::io::stdout().lock());
    if args[1] == "--at" && args.len() == 4 {
        let time = GnssTime::new(args[2].parse()?, args[3].parse()?)?;
        writeln!(
            out,
            "week,tow,satellite,tracking,status,code_add_m,sigma_upper_bound_m,source_terms,error"
        )?;
        for prn in 1..=32 {
            let sat: SatelliteId = format!("G{prn:02}").parse()?;
            for code in ["1W", "2W", "1C", "2L"] {
                match model.corrections_at(time, &binding, &[request(sat, code)]) {
                    Ok(snapshot) => {
                        let d = &snapshot.provenance[&(sat, code.into())];
                        let terms = d
                            .source_terms
                            .iter()
                            .map(|t| {
                                format!(
                                    "{}-{}@line{}*{:.17e}*{:.17e}{}",
                                    t.entry.observation_1,
                                    t.entry.observation_2.as_ref().unwrap(),
                                    t.entry.line,
                                    t.coefficient,
                                    t.entry.value,
                                    t.entry.unit_label
                                )
                            })
                            .collect::<Vec<_>>()
                            .join(";");
                        writeln!(
                            out,
                            "{},{:.17},{sat},{code},AVAILABLE,{:.17e},{:.17e},{terms},",
                            time.week(),
                            time.tow(),
                            d.code_add_m,
                            d.sigma_upper_bound_m
                        )?;
                    }
                    Err(e) => writeln!(
                        out,
                        "{},{:.17},{sat},{code},{},,,,{}",
                        time.week(),
                        time.tow(),
                        error_kind(&e),
                        quoted(&e.to_string())
                    )?,
                }
            }
        }
    } else if args[1] == "--rover" && args.len() == 3 {
        let reader = RinexObservationReader::new(BufReader::new(File::open(&args[2])?))?;
        writeln!(
            out,
            "week,tow,flag,kind,gps_satellites,complete_1c2w,unflagged_1c2w,dsb_available_1c2w,dsb_missing_1c2w,dsb_expired_1c2w,dsb_ambiguous_1c2w,dsb_other_1c2w,complete_1c2l,unflagged_1c2l,dsb_available_1c2l,dsb_missing_1c2l,dsb_expired_1c2l,dsb_ambiguous_1c2l,dsb_other_1c2l,missing_lli_pair_rows"
        )?;
        for e in reader {
            let e = e?;
            let mut rows = BTreeMap::new();
            let mut satellites = BTreeSet::new();
            for o in &e.observations {
                if o.satellite.system() == GnssSystem::Gps {
                    satellites.insert(o.satellite);
                    if rows
                        .insert((o.satellite, o.tracking_code.as_str()), o)
                        .is_some()
                    {
                        return Err("duplicate GPS tracking row".into());
                    }
                }
            }
            let kind = if !matches!(e.flag, 0 | 1) {
                "EVENT"
            } else if e
                .receiver_clock_offset_s
                .is_some_and(|v| !v.is_finite() || v != 0.0)
            {
                "UNSUPPORTED_OFFSET"
            } else {
                "OBSERVATION"
            };
            let mut counts = [0usize; 14];
            let mut missing_lli = 0;
            if kind == "OBSERVATION" {
                for &sat in &satellites {
                    for (i, second) in ["2W", "2L"].iter().enumerate() {
                        let Some(first) = rows.get(&(sat, "1C")) else {
                            continue;
                        };
                        let Some(second_row) = rows.get(&(sat, *second)) else {
                            continue;
                        };
                        if ![first, second_row].iter().all(|o| {
                            o.pseudorange_m.is_some_and(|v| v.is_finite() && v > 0.0)
                                && o.carrier_phase_cycles.is_some_and(f64::is_finite)
                        }) {
                            continue;
                        }
                        let start = i * 7;
                        counts[start] += 1;
                        missing_lli += [first, second_row]
                            .iter()
                            .filter(|o| o.lli.is_none())
                            .count();
                        if [first, second_row]
                            .iter()
                            .all(|o| o.lli.is_none_or(|v| v & 3 == 0))
                        {
                            counts[start + 1] += 1;
                        }
                        let index = match model.corrections_at(
                            e.time,
                            &binding,
                            &[request(sat, "1C"), request(sat, second)],
                        ) {
                            Ok(_) => 2,
                            Err(DsbClockError::MissingPair { .. }) => 3,
                            Err(DsbClockError::ExpiredPair { .. }) => 4,
                            Err(DsbClockError::AmbiguousPair { .. }) => 5,
                            Err(_) => 6,
                        };
                        counts[start + index] += 1;
                    }
                }
            }
            write!(
                out,
                "{},{:.17},{},{kind},{}",
                e.time.week(),
                e.time.tow(),
                e.flag,
                satellites.len()
            )?;
            for count in counts {
                write!(out, ",{count}")?;
            }
            writeln!(out, ",{missing_lli}")?;
        }
    } else {
        return Err("invalid ppc_dsb_clock mode".into());
    }
    out.flush()?;
    eprintln!(
        "GPS_IF_C1W_C2W_gauge=explicit clock_compatibility_verified=false phase_OSB_generated=false DSB_sigma_correlations=unknown upper_bound_only=true PPP_executed=false positioning_executed=false performance_accepted=false excluded_receiver_rows={}",
        selected.excluded_receiver_records.len()
    );
    Ok(())
}
