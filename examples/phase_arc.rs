//! Synthetic continuity demonstration; accepts phases without estimating PPP.
//! One satellite update is deliberately rejected to demonstrate recovery.
use gnss_rust::dual_frequency::{
    DualFrequencyConfig, DualFrequencyNoise, SignalCorrection, SignalCorrections,
};
use gnss_rust::phase_arc::{PhaseArcConfig, PhaseArcTracker};
use gnss_rust::precise::ClockSource;
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::{GnssTime, SatelliteId};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Cursor,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let selections: BTreeMap<SatelliteId, DualFrequencyConfig> = ["G01", "G02"]
        .into_iter()
        .map(|name| {
            Ok((
                name.parse()?,
                DualFrequencyConfig {
                    tracking_codes: ["1C".to_owned(), "2W".to_owned()],
                    clock_source: ClockSource::RinexClk,
                    clock_reference: "synthetic-clock".to_owned(),
                    phase_reference: "synthetic-phase".to_owned(),
                    noise: DualFrequencyNoise {
                        code_m2: [[0.25, 0.0], [0.0, 0.25]],
                        phase_m2: [[0.0001, 0.0], [0.0, 0.0001]],
                    },
                },
            ))
        })
        .collect::<Result<_, gnss_rust::Error>>()?;
    let mut entries = Vec::new();
    for (&sat, selection) in &selections {
        for code in &selection.tracking_codes {
            entries.push(SignalCorrection {
                satellite: sat,
                tracking_code: code.clone(),
                code_add_m: 0.0,
                phase_add_m: Some(0.0),
                clock_source: selection.clock_source,
                clock_reference: selection.clock_reference.clone(),
                phase_reference: Some(selection.phase_reference.clone()),
                valid_from: GnssTime::new(2300, 345600.0)?,
                valid_until: GnssTime::new(2300, 346200.0)?,
            });
        }
    }
    let corrections = SignalCorrections::new(entries)?;
    let reader = RinexObservationReader::new(Cursor::new(include_str!(
        "../tests/fixtures/synthetic_phase_arc.obs"
    )))?;
    let mut tracker = PhaseArcTracker::new(PhaseArcConfig::default())?;
    println!("sat,week,tow,proposed_arc,accepted_arc,accepted_epochs,reasons");
    for (index, epoch) in reader.enumerate() {
        let epoch = epoch?;
        let decisions = tracker
            .prepare_epoch(&epoch, &corrections, &selections)?
            .decisions()
            .clone();
        let mut accepted: BTreeSet<_> = decisions
            .iter()
            .filter_map(|(&sat, d)| d.proposed_arc.map(|_| sat))
            .collect();
        if index == 7 {
            accepted.remove(&"G01".parse()?);
        }
        let commit = tracker.finish_epoch(&accepted)?;
        for (&sat, d) in &decisions {
            let arc = commit.accepted.get(&sat);
            let reasons = commit.interrupted.get(&sat).unwrap_or(&d.reset_reasons);
            println!(
                "{sat},{},{:.3},{},{},{},\"{}\"",
                epoch.time.week(),
                epoch.time.tow(),
                d.proposed_arc
                    .map(|id| id.value().to_string())
                    .unwrap_or_default(),
                arc.map(|a| a.id.value().to_string()).unwrap_or_default(),
                arc.map(|a| a.accepted_epochs.to_string())
                    .unwrap_or_default(),
                format!("{reasons:?}").replace('"', "\"\"")
            );
        }
    }
    Ok(())
}
