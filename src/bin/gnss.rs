use gnss_rust::rinex::{RinexObservationReader, read_navigation};
use gnss_rust::rtk::{RtkFloatConfig, RtkFloatFilter};
use gnss_rust::rtk_fix::{RtkFixConfig, RtkFixPolicy};
use gnss_rust::spp::{SppConfig, solve_epoch};
use gnss_rust::{EcefCoord, lambda::LambdaConfig};
use std::{
    env,
    error::Error,
    fs::{File, OpenOptions},
    io::{self, BufReader, Write},
    path::PathBuf,
};

const HELP: &str = "gnss-rust spp --obs OBSERVATIONS --nav NAVIGATION [--output CSV] [--no-atmosphere]\ngnss-rust rtk-float --rover OBS --base OBS --nav NAV --base-position X,Y,Z [--rover-position X,Y,Z] [--output CSV]\ngnss-rust rtk --rover OBS --base OBS --nav NAV --base-position X,Y,Z [--rover-position X,Y,Z] [--output CSV]\n\nInitial SPP and experimental static-baseline GPS L1 C/A RTK from RINEX 3.\nrtk validates full-set integer candidates; rtk-float outputs only FLOAT.\nPPP is not implemented. Base coordinates must be surveyed ECEF metres.\nRTK pairs must have matching epoch sequences; no interpolation is performed.\n--no-atmosphere is for vacuum SPP inputs; normally leave models enabled.\n";

fn run() -> Result<(), Box<dyn Error>> {
    let mut args = env::args_os().skip(1);
    let Some(command) = args.next() else {
        return Err(HELP.into());
    };
    if command == "--help" || command == "-h" {
        print!("{HELP}");
        return Ok(());
    }
    if command == "rtk-float" || command == "rtk" {
        return run_rtk(args, command == "rtk");
    }
    if command != "spp" {
        return Err("expected 'spp', 'rtk-float' or 'rtk'; PPP is not implemented".into());
    }
    let (mut observations, mut navigation, mut output) = (None, None, None);
    let mut config = SppConfig::default();
    while let Some(arg) = args.next() {
        if arg == "--help" || arg == "-h" {
            print!("{HELP}");
            return Ok(());
        }
        if arg == "--no-atmosphere" {
            config.use_troposphere = false;
            config.use_ionosphere = false;
            continue;
        }
        let target = if arg == "--obs" {
            &mut observations
        } else if arg == "--nav" {
            &mut navigation
        } else if arg == "--output" {
            &mut output
        } else {
            return Err(format!("unknown option {}", arg.to_string_lossy()).into());
        };
        if target.is_some() {
            return Err(format!("duplicate option {}", arg.to_string_lossy()).into());
        }
        *target = Some(PathBuf::from(args.next().ok_or("missing option value")?));
    }
    let obs_path = observations.ok_or("--obs is required")?;
    let nav_path = navigation.ok_or("--nav is required")?;
    let obs_file = File::open(&obs_path)?;
    let nav = read_navigation(BufReader::new(File::open(&nav_path)?))?;
    let epochs = RinexObservationReader::new(BufReader::new(obs_file))?;
    if let Some(path) = output.as_ref() {
        // Prevent accidentally truncating either input, including symlink aliases.
        if path.exists() {
            let resolved = path.canonicalize()?;
            if resolved == obs_path.canonicalize()? || resolved == nav_path.canonicalize()? {
                return Err("output must differ from input files".into());
            }
        }
    }
    let mut out: Box<dyn Write> = match output {
        Some(path) => Box::new(OpenOptions::new().write(true).create_new(true).open(path)?),
        None => Box::new(io::stdout().lock()),
    };
    writeln!(
        out,
        "week,tow,x_m,y_m,z_m,receiver_clock_bias_s,satellites,status,residual_rms_m,iterations"
    )?;
    let (mut processed, mut skipped) = (0, 0);
    for epoch in epochs {
        let epoch = epoch?;
        if epoch.flag == 6 {
            skipped += 1;
            continue;
        }
        let solution = solve_epoch(&epoch, &nav, config).map_err(|e| {
            format!(
                "GPS week {} tow {}: {e}",
                epoch.time.week(),
                epoch.time.tow()
            )
        })?;
        let position = solution
            .position
            .position_ecef
            .ok_or("solver did not return a position")?;
        writeln!(
            out,
            "{},{:.7},{:.6},{:.6},{:.6},{:.12},{},SPP,{:.6},{}",
            epoch.time.week(),
            epoch.time.tow(),
            position.x,
            position.y,
            position.z,
            solution.position.receiver_clock_bias,
            solution.position.num_satellites,
            solution.residual_rms_m,
            solution.iterations
        )?;
        processed += 1;
    }
    if processed == 0 {
        return Err("no measurement epochs were solved".into());
    }
    out.flush()?;
    eprintln!(
        "processed {processed} measurement epochs; valid {processed}; skipped {skipped} cycle-slip records"
    );
    Ok(())
}

fn parse_position(text: &std::ffi::OsStr) -> Result<EcefCoord, Box<dyn Error>> {
    let values = text
        .to_str()
        .ok_or("coordinates must be UTF-8")?
        .split(',')
        .map(str::parse::<f64>)
        .collect::<Result<Vec<_>, _>>()?;
    if values.len() != 3 || !values.iter().all(|v| v.is_finite()) {
        return Err("expected finite ECEF X,Y,Z in metres".into());
    }
    Ok(EcefCoord::new(values[0], values[1], values[2]))
}

fn run_rtk(
    mut args: impl Iterator<Item = std::ffi::OsString>,
    enable_fix: bool,
) -> Result<(), Box<dyn Error>> {
    let (mut rover_path, mut base_path, mut nav_path, mut output) = (None, None, None, None);
    let (mut base_position, mut rover_position) = (None, None);
    let mut seen = std::collections::BTreeSet::new();
    while let Some(arg) = args.next() {
        if arg == "--help" || arg == "-h" {
            print!("{HELP}");
            return Ok(());
        }
        if !seen.insert(arg.clone()) {
            return Err(format!("duplicate option {}", arg.to_string_lossy()).into());
        }
        let value = args.next().ok_or("missing option value")?;
        if arg == "--base-position" {
            base_position = Some(parse_position(&value)?);
        } else if arg == "--rover-position" {
            rover_position = Some(parse_position(&value)?);
        } else {
            let target = if arg == "--rover" {
                &mut rover_path
            } else if arg == "--base" {
                &mut base_path
            } else if arg == "--nav" {
                &mut nav_path
            } else if arg == "--output" {
                &mut output
            } else {
                return Err(format!("unknown option {}", arg.to_string_lossy()).into());
            };
            *target = Some(PathBuf::from(value));
        }
    }
    let base_position = base_position.ok_or("--base-position surveyed ECEF X,Y,Z is required")?;
    let mut rover = RinexObservationReader::new(BufReader::new(File::open(
        rover_path.ok_or("--rover is required")?,
    )?))?;
    let mut base = RinexObservationReader::new(BufReader::new(File::open(
        base_path.ok_or("--base is required")?,
    )?))?;
    let nav = read_navigation(BufReader::new(File::open(
        nav_path.ok_or("--nav is required")?,
    )?))?;
    let seed = rover_position
        .or(rover.header().approximate_position)
        .ok_or("--rover-position or a RINEX approximate rover position is required")?;
    let mut filter = RtkFloatFilter::new(base_position, seed, RtkFloatConfig::default())?;
    let mut policy = RtkFixPolicy::new(RtkFixConfig::default())?;
    let mut out: Box<dyn Write> = match output {
        Some(path) => Box::new(OpenOptions::new().write(true).create_new(true).open(path)?),
        None => Box::new(io::stdout().lock()),
    };
    write!(
        out,
        "week,tow,x_m,y_m,z_m,satellites,status,reference,code_residual_rms_m,phase_residual_rms_m,{},iterations,ambiguity_resets,candidate_ratio",
        if enable_fix {
            "float_update_nis"
        } else {
            "nis"
        }
    )?;
    if enable_fix {
        write!(
            out,
            ",fix_decision,candidate_pairs,confirmations,candidate_postfit_nis_per_row"
        )?;
    }
    writeln!(out)?;
    let (mut processed, mut skipped, mut diagnostic_failures, mut fixed) = (0, 0, 0, 0);
    loop {
        let (ro,ba)=match (rover.next(),base.next()) {
            (None,None)=>break,
            (Some(r),Some(b))=>(r?,b?),
            _=>return Err("base and rover files have different epoch counts; matching epoch sequences are required".into()),
        };
        if ro.time.difference_seconds(ba.time).abs() > 1e-7 {
            return Err(format!(
                "unmatched rover/base times: GPS week {} tow {} / week {} tow {}",
                ro.time.week(),
                ro.time.tow(),
                ba.time.week(),
                ba.time.tow()
            )
            .into());
        }
        if ro.flag == 6 || ba.flag == 6 {
            filter.invalidate_ambiguities();
            policy.reset();
            skipped += 1;
            continue;
        }
        let result = filter
            .process_pair(&ro, &ba, &nav)
            .map_err(|e| format!("GPS week {} tow {}: {e}", ro.time.week(), ro.time.tow()))?;
        let attempt = if enable_fix {
            Some(policy.evaluate(&result)?)
        } else {
            None
        };
        let selected = attempt
            .as_ref()
            .and_then(|a| a.fixed_solution.as_ref())
            .unwrap_or(&result.position);
        let is_fixed = selected.is_fixed();
        if is_fixed {
            fixed += 1;
        }
        let position = selected.position_ecef.ok_or("RTK returned no position")?;
        let ratio = if let Some(a) = &attempt {
            if matches!(
                a.decision,
                gnss_rust::rtk_fix::FixDecision::SearchFailed
                    | gnss_rust::rtk_fix::FixDecision::ConditioningFailed
            ) {
                diagnostic_failures += 1;
                eprintln!(
                    "GPS week {} tow {}: FIX unavailable: {}",
                    ro.time.week(),
                    ro.time.tow(),
                    a.decision
                );
            }
            a.ratio.map(|r| format!("{r:.6}")).unwrap_or_default()
        } else {
            match result.integer_candidates(LambdaConfig::default()) {
                Ok(candidates) => candidates
                    .ratio()
                    .map(|r| format!("{r:.6}"))
                    .unwrap_or_default(),
                Err(e) => {
                    diagnostic_failures += 1;
                    eprintln!(
                        "GPS week {} tow {}: ambiguity candidates unavailable: {e}",
                        ro.time.week(),
                        ro.time.tow()
                    );
                    String::new()
                }
            }
        };
        let validated = attempt
            .as_ref()
            .filter(|_| is_fixed)
            .and_then(|a| a.validation.as_ref());
        let code_rms = validated.map_or(result.code_residual_rms_m, |v| v.code_residual_rms_m);
        let phase_rms = validated.map_or(result.phase_residual_rms_m, |v| v.phase_residual_rms_m);
        write!(
            out,
            "{},{:.7},{:.6},{:.6},{:.6},{},{},{},{:.6},{:.6},{:.6},{},{},{}",
            ro.time.week(),
            ro.time.tow(),
            position.x,
            position.y,
            position.z,
            result.position.num_satellites,
            if is_fixed { "FIX" } else { "FLOAT" },
            result.reference_satellite,
            code_rms,
            phase_rms,
            result.normalized_innovation_squared,
            result.iterations,
            result.reset_ambiguities.len(),
            ratio
        )?;
        if let Some(a) = &attempt {
            let postfit = a
                .validation
                .as_ref()
                .map(|v| format!("{:.6}", v.postfit_nis_per_row))
                .unwrap_or_default();
            write!(
                out,
                ",{},{},{},{}",
                a.decision,
                a.integers.len(),
                a.confirmation_count,
                postfit
            )?;
        }
        writeln!(out)?;
        processed += 1;
    }
    if processed == 0 {
        return Err("no measurement epoch pairs were solved".into());
    }
    out.flush()?;
    if enable_fix {
        eprintln!(
            "processed {processed} matched epoch pairs; valid FLOAT {}; valid FIX {fixed}; skipped {skipped} cycle-slip pairs; unavailable ambiguity diagnostics {diagnostic_failures}",
            processed - fixed
        );
    } else {
        eprintln!(
            "processed {processed} matched epoch pairs; valid FLOAT {processed}; skipped {skipped} cycle-slip pairs; unavailable ambiguity diagnostics {diagnostic_failures}"
        );
    }
    Ok(())
}

fn main() -> std::process::ExitCode {
    match run() {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("gnss-rust: {e}");
            std::process::ExitCode::FAILURE
        }
    }
}
