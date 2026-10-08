//! Run the documented synthetic precise code-SPP model (IGS20/APC/no ionosphere).
//! The bias CSV is a demo ledger, not a RINEX/Bias-SINEX DCB/OSB reader.
use gnss_rust::GnssTime;
use gnss_rust::precise::{ClockSource, OrbitReferencePoint, PreciseProducts, read_clk, read_sp3};
use gnss_rust::precise_spp::{
    PreciseCodeBias, PreciseCodeBiases, PreciseSppConfig, solve_epoch_precise,
};
use gnss_rust::precise_transmit::{PreciseTransmitConfig, PreciseTransmitMethod};
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::spp::SppConfig;
use std::{
    env,
    error::Error,
    fs::{self, File},
    io::BufReader,
};

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = env::args_os().skip(1);
    let usage = "usage: precise_spp OBS SP3 CLK BIAS_CSV [emission] (synthetic model)";
    let observations =
        RinexObservationReader::new(BufReader::new(File::open(args.next().ok_or(usage)?)?))?;
    let orbit = read_sp3(BufReader::new(File::open(args.next().ok_or(usage)?)?))?;
    let clocks = read_clk(BufReader::new(File::open(args.next().ok_or(usage)?)?))?;
    let text = fs::read_to_string(args.next().ok_or(usage)?)?;
    let mut entries = Vec::new();
    for (line, row) in text.lines().enumerate() {
        if row.starts_with('#') || row.trim().is_empty() {
            continue;
        }
        let f: Vec<_> = row.split(',').collect();
        if f.len() != 9 {
            return Err(
                format!("bias ledger line {}: expected 9 unquoted fields", line + 1).into(),
            );
        }
        let source = match f[2] {
            "RinexClk" => ClockSource::RinexClk,
            "Sp3" => ClockSource::Sp3,
            _ => return Err("unknown bias clock source".into()),
        };
        entries.push(PreciseCodeBias {
            satellite: f[0].parse()?,
            tracking_code: f[1].to_owned(),
            clock_source: source,
            valid_from: GnssTime::new(f[3].parse()?, f[4].parse()?)?,
            valid_until: GnssTime::new(f[5].parse()?, f[6].parse()?)?,
            correction_m: f[7].parse()?,
            clock_reference: f[8].to_owned(),
        });
    }
    let method = match args.next() {
        None => PreciseTransmitMethod::ReceptionTaylor,
        Some(s) if s == "emission" => PreciseTransmitMethod::IteratedEmission,
        _ => return Err(usage.into()),
    };
    if args.next().is_some() {
        return Err(usage.into());
    }
    let products = PreciseProducts::new(orbit, Some(clocks));
    let biases = PreciseCodeBiases::new(entries)?;
    let config = PreciseSppConfig {
        spp: SppConfig {
            use_ionosphere: false,
            ..Default::default()
        },
        transmit: PreciseTransmitConfig {
            method,
            ..Default::default()
        },
        expected_coordinate_system: "IGS20".to_owned(),
        expected_reference_point: OrbitReferencePoint::AntennaPhaseCenter,
        clock_reference: "synthetic-precise-clock-v1".to_owned(),
        ionosphere: None,
    };
    println!("week,tow,x_m,y_m,z_m,receiver_clock_s,satellites,status,residual_rms_m,iterations");
    let mut count = 0;
    for epoch in observations {
        let epoch = epoch?;
        let result = solve_epoch_precise(&epoch, &products, &biases, &config)?;
        let p = result
            .position
            .position_ecef
            .ok_or("no precise code position")?;
        for (sat, reason) in &result.rejected {
            eprintln!("{sat} at {}: rejected {reason:?}", epoch.time.tow());
        }
        println!(
            "{},{:.7},{:.9},{:.9},{:.9},{:.15e},{},SPP,{:.9},{}",
            epoch.time.week(),
            epoch.time.tow(),
            p.x,
            p.y,
            p.z,
            result.position.receiver_clock_bias,
            result.position.num_satellites,
            result.residual_rms_m,
            result.iterations
        );
        count += 1;
    }
    if count == 0 {
        return Err("no precise code epochs".into());
    }
    eprintln!("processed {count} precise code SPP epochs ({method:?})");
    Ok(())
}
