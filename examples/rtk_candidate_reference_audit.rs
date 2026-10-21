//! Synthetic diagnostic only: audit one original-basis LAMBDA candidate.
use gnss_rust::lambda::LambdaConfig;
use gnss_rust::rinex::{RinexObservationReader, read_navigation};
use gnss_rust::rtk::{RtkFloatConfig, RtkFloatFilter};
use gnss_rust::rtk_fix::{RtkFixConfig, audit_candidate_references};
use gnss_rust::{EcefCoord, SatelliteId};
use std::collections::BTreeMap;
use std::io::{BufWriter, Cursor, Write};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let arguments: Vec<_> = std::env::args().collect();
    let mode = arguments.get(1).map_or("clean", String::as_str);
    if arguments.len() > 2 || !matches!(mode, "clean" | "noisy" | "phase-bias") {
        return Err("usage: rtk_candidate_reference_audit [clean|noisy|phase-bias]".into());
    }
    let mut truth = include_str!("../tests/fixtures/synthetic_rtk_truth.csv")
        .lines()
        .filter(|line| !line.starts_with('#'));
    let coordinates: Vec<f64> = truth
        .next()
        .ok_or("missing synthetic coordinates")?
        .split(',')
        .map(str::parse)
        .collect::<Result<_, _>>()?;
    let base_position = EcefCoord::new(coordinates[0], coordinates[1], coordinates[2]);
    let seed = EcefCoord::new(coordinates[6], coordinates[7], coordinates[8]);
    let sd: BTreeMap<SatelliteId, i64> = truth
        .map(|line| {
            let (satellite, integer) = line.split_once(',').ok_or("invalid synthetic SD")?;
            Ok((satellite.parse()?, integer.parse()?))
        })
        .collect::<Result<_, Box<dyn std::error::Error>>>()?;
    let target = *sd.keys().next().ok_or("missing synthetic satellite")?;
    let nav = read_navigation(Cursor::new(include_str!(
        "../tests/fixtures/synthetic_rtk.nav"
    )))?;
    let rover_text = if mode == "noisy" {
        include_str!("../tests/fixtures/synthetic_rtk_rover_noisy.obs")
    } else {
        include_str!("../tests/fixtures/synthetic_rtk_rover.obs")
    };
    let base_text = if mode == "noisy" {
        include_str!("../tests/fixtures/synthetic_rtk_base_noisy.obs")
    } else {
        include_str!("../tests/fixtures/synthetic_rtk_base.obs")
    };
    let rover =
        RinexObservationReader::new(Cursor::new(rover_text))?.collect::<Result<Vec<_>, _>>()?;
    let base =
        RinexObservationReader::new(Cursor::new(base_text))?.collect::<Result<Vec<_>, _>>()?;
    if rover.len() != base.len() {
        return Err("synthetic input length mismatch".into());
    }
    let mut engine = RtkFloatFilter::new(base_position, seed, RtkFloatConfig::default())?;
    let mut output = BufWriter::new(std::io::stdout().lock());
    writeln!(
        output,
        "case,week,tow,source_reference,reference,candidate_correct_for_synthetic_truth,targets,integers,baseline_x_m,baseline_y_m,baseline_z_m,qxx,qxy,qxz,qyy,qyz,qzz,ambiguity_cost,phase_rms_m,phase_max_m,code_rms_m,code_phase_rms_m,code_phase_max_m,code_phase_bad_pairs,nominal_postfit_score_per_row,local_rejection,reference_sensitive_local_rejection"
    )?;
    let mut sensitive = 0;
    let mut epochs = 0;
    let mut rows = 0;
    for (index, (mut rover, base)) in rover.into_iter().zip(base).enumerate() {
        if mode == "phase-bias" && index >= 10 {
            *rover
                .observations
                .iter_mut()
                .find(|o| o.satellite == target)
                .and_then(|o| o.carrier_phase_cycles.as_mut())
                .ok_or("missing synthetic fault target")? += 0.15;
        }
        let float = engine.process_pair(&rover, &base, &nav)?;
        let original = float.integer_candidates(LambdaConfig::default())?;
        let integers = &original
            .candidates
            .first()
            .ok_or("no original-basis candidate")?
            .ambiguities;
        let correct = float
            .ambiguity_satellites
            .iter()
            .zip(integers)
            .all(|(satellite, integer)| sd[&float.reference_satellite] - sd[satellite] == *integer);
        let before = engine.state().ok_or("no FLOAT state")?.clone();
        let audit = audit_candidate_references(&float, integers, RtkFixConfig::default())?;
        if engine.state() != Some(&before) {
            return Err("read-only audit changed FLOAT posterior".into());
        }
        let changed = audit.local_rejection_depends_on_reference();
        sensitive += usize::from(changed);
        epochs += 1;
        for entry in &audit.references {
            let validation = &entry.validation;
            let mean = &validation.conditioned.state.mean;
            let q = &validation.conditioned.state.covariance;
            let values = [
                mean[0],
                mean[1],
                mean[2],
                q[0][0],
                q[0][1],
                q[0][2],
                q[1][1],
                q[1][2],
                q[2][2],
                validation.conditioned.squared_ambiguity_residual,
                validation.phase_residual_rms_m,
                validation.phase_residual_max_abs_m,
                validation.code_residual_rms_m,
                validation.code_phase.rms_innovation_m,
                validation.code_phase.max_abs_innovation_m,
            ];
            let numeric = values
                .iter()
                .map(|v| format!("{v:.17e}"))
                .collect::<Vec<_>>()
                .join(",");
            let targets = entry
                .ambiguity_satellites
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(";");
            let integers = entry
                .integers
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(";");
            let rejection = entry
                .validation
                .rejection
                .map_or_else(|| "pass".to_owned(), |v| v.to_string());
            writeln!(
                output,
                "{mode},{},{:.17},{},{},{correct},{targets},{integers},{numeric},{},{:.17e},{rejection},{changed}",
                audit.time.week(),
                audit.time.tow(),
                audit.source_reference,
                entry.reference_satellite,
                validation.code_phase.bad_pairs,
                validation.postfit_nis_per_row
            )?;
            rows += 1;
        }
    }
    output.flush()?;
    eprintln!(
        "synthetic_only=true mode={mode} epochs={epochs} reference_rows={rows} reference_sensitive_epochs={sensitive} candidate_source=original_basis_LAMBDA best_reference_selected=false FIX_policy_evaluated=false phase_bias_applied={} phase_bias_cycles=0.15 phase_bias_first_epoch=10 phase_bias_target={target} performance_accepted=false",
        mode == "phase-bias"
    );
    Ok(())
}
