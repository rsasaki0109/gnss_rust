//! Inspect the synthetic GPS L1/L2 fixture with explicitly declared zero biases.
//! This example does not calibrate real observations or solve a PPP position.
use gnss_rust::dual_frequency::{
    DualFrequencyConfig, DualFrequencyNoise, SignalCorrection, SignalCorrections,
    form_ionosphere_free,
};
use gnss_rust::precise::ClockSource;
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::{GnssTime, SatelliteId};
use std::io::Cursor;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let satellite: SatelliteId = "G01".parse()?;
    let config = DualFrequencyConfig {
        tracking_codes: ["1C".to_owned(), "2W".to_owned()],
        clock_source: ClockSource::RinexClk,
        clock_reference: "synthetic-clock".to_owned(),
        phase_reference: "synthetic-phase".to_owned(),
        noise: DualFrequencyNoise {
            code_m2: [[0.25, 0.0], [0.0, 0.25]],
            phase_m2: [[0.0001, 0.0], [0.0, 0.0001]],
        },
    };
    let corrections = SignalCorrections::new(
        config
            .tracking_codes
            .iter()
            .map(|code| {
                Ok(SignalCorrection {
                    satellite,
                    tracking_code: code.clone(),
                    code_add_m: 0.0,
                    phase_add_m: Some(0.0),
                    clock_source: config.clock_source,
                    clock_reference: config.clock_reference.clone(),
                    phase_reference: Some(config.phase_reference.clone()),
                    valid_from: GnssTime::new(2300, 345600.0)?,
                    valid_until: GnssTime::new(2300, 345690.0)?,
                })
            })
            .collect::<Result<Vec<_>, gnss_rust::Error>>()?,
    )?;
    let reader = RinexObservationReader::new(Cursor::new(include_str!(
        "../tests/fixtures/synthetic_dual_frequency.obs"
    )))?;
    println!("sat,week,tow,code_if_m,phase_if_m,phase_minus_code_m,gf_m,mw_m,reset_required");
    for epoch in reader {
        let epoch = epoch?;
        let measurement = form_ionosphere_free(&epoch, satellite, &corrections, &config)?;
        let phase = measurement.phase.ok_or("synthetic phase missing")?;
        println!(
            "{satellite},{},{:.3},{:.9},{:.9},{:.9},{:.9},{:.9},{}",
            epoch.time.week(),
            epoch.time.tow(),
            measurement.code_if_m,
            phase.phase_if_m,
            phase.phase_if_m - measurement.code_if_m,
            phase.raw_geometry_free_m,
            phase.raw_melbourne_wubbena_m.ok_or("MW missing")?,
            phase.reset_required,
        );
    }
    Ok(())
}
