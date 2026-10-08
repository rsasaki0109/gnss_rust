use gnss_rust::dual_frequency::*;
use gnss_rust::observation::{Observation, ObservationEpoch};
use gnss_rust::precise::ClockSource;
use gnss_rust::rinex::RinexObservationReader;
use gnss_rust::{Error, GnssSystem, GnssTime, SatelliteId, constants as c};
use std::io::Cursor;

fn rows() -> Vec<Vec<String>> {
    include_str!("fixtures/upstream_dual_frequency.csv")
        .lines()
        .skip(1)
        .map(|l| l.split(',').map(str::to_owned).collect())
        .collect()
}
fn number(row: &[String], i: usize) -> f64 {
    row[i].parse().unwrap()
}
fn config(row: &[String]) -> DualFrequencyConfig {
    DualFrequencyConfig {
        tracking_codes: [row[1].clone(), row[2].clone()],
        clock_source: ClockSource::RinexClk,
        clock_reference: "synthetic-clock".to_owned(),
        phase_reference: "synthetic-phase".to_owned(),
        noise: DualFrequencyNoise {
            code_m2: [[0.25, 0.0], [0.0, 0.25]],
            phase_m2: [[0.0001, 0.0], [0.0, 0.0001]],
        },
    }
}
fn epoch(row: &[String]) -> ObservationEpoch {
    let sat = row[0].parse().unwrap();
    let observations = (0..2)
        .map(|i| {
            let mut o = Observation::new(sat, row[i + 1].clone());
            o.pseudorange_m = Some(number(row, 7 + i));
            o.carrier_phase_cycles = Some(number(row, 9 + i));
            o.lli = Some(0);
            o.glonass_channel = if row[3] == "99" {
                None
            } else {
                Some(row[3].parse().unwrap())
            };
            o
        })
        .collect();
    ObservationEpoch {
        time: GnssTime::new(2300, 345600.0 + number(row, 4) * 30.0).unwrap(),
        flag: 0,
        approximate_position: None,
        receiver_clock_offset_s: None,
        observations,
    }
}
fn correction_rows(epoch: &ObservationEpoch) -> Vec<SignalCorrection> {
    epoch
        .observations
        .iter()
        .map(|o| SignalCorrection {
            satellite: o.satellite,
            tracking_code: o.tracking_code.clone(),
            code_add_m: 0.0,
            phase_add_m: Some(0.0),
            clock_source: ClockSource::RinexClk,
            clock_reference: "synthetic-clock".to_owned(),
            phase_reference: Some("synthetic-phase".to_owned()),
            valid_from: GnssTime::new(2300, 345600.0).unwrap(),
            valid_until: GnssTime::new(2300, 345690.0).unwrap(),
        })
        .collect()
}
fn form(
    epoch: &ObservationEpoch,
    entries: Vec<SignalCorrection>,
    config: &DualFrequencyConfig,
) -> Result<IonosphereFreeObservation, DualFrequencyError> {
    form_ionosphere_free(
        epoch,
        epoch.observations[0].satellite,
        &SignalCorrections::new(entries)?,
        config,
    )
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:.17e} vs {expected:.17e}, tolerance {tolerance}"
    );
}

#[test]
fn native_48_cases_cover_five_systems_frequencies_units_and_noise() {
    let rows = rows();
    assert_eq!(rows.len(), 48);
    for row in rows {
        let e = epoch(&row);
        let r = form(&e, correction_rows(&e), &config(&row)).unwrap();
        for i in 0..2 {
            close(r.frequency_hz[i], number(&row, 5 + i), 0.0);
            close(r.coefficients[i], number(&row, 11 + i), 1e-14);
            close(
                r.wavelength_m[i],
                c::SPEED_OF_LIGHT / number(&row, 5 + i),
                0.0,
            );
        }
        close(r.code_if_m, number(&row, 13), 1e-7);
        let p = r.phase.unwrap();
        close(p.phase_if_m, number(&row, 14), 1e-7);
        close(p.raw_geometry_free_m, number(&row, 15), 1e-7);
        let mw = melbourne_wubbena_m(
            [number(&row, 9), number(&row, 10)],
            [number(&row, 7), number(&row, 8)],
            r.frequency_hz,
        )
        .unwrap();
        close(mw, number(&row, 16), 1e-7);
        if r.satellite.system() == GnssSystem::Glonass {
            assert!(p.raw_melbourne_wubbena_m.is_none());
        } else {
            close(p.raw_melbourne_wubbena_m.unwrap(), number(&row, 16), 1e-7);
        }
        close(r.code_variance_m2, number(&row, 17), 1e-12);
        close(p.variance_m2, number(&row, 18), 1e-14);
        assert_eq!(r.raw, e.observations.as_slice());
        assert!(r.phase_unavailability.is_none());
        assert!(!p.reset_required);
    }
}

#[test]
fn dispersive_ionosphere_cancels_but_if_ambiguity_remains_in_metres() {
    for row in rows() {
        let e = epoch(&row);
        let r = form(&e, correction_rows(&e), &config(&row)).unwrap();
        let range = 23000000.125 + number(&row, 4) * 17.5;
        close(r.code_if_m, range, 1e-6);
        close(r.phase.unwrap().phase_if_m - range, number(&row, 19), 1e-6);
        // There is no integer rounding or common wavelength applied to IF phase.
        assert!(number(&row, 19).abs() > 0.1);
    }
}

#[test]
fn native_rinex_file_flows_through_reader_to_48_measurements() {
    let mut reader = RinexObservationReader::new(Cursor::new(include_str!(
        "fixtures/synthetic_dual_frequency.obs"
    )))
    .unwrap();
    assert_eq!(reader.header().glonass_channels.len(), 2);
    let epochs = reader.by_ref().collect::<Result<Vec<_>, _>>().unwrap();
    assert_eq!(epochs.len(), 4);
    let reference = rows();
    let mut count = 0;
    for row in reference {
        let e = &epochs[number(&row, 4) as usize];
        let sat: SatelliteId = row[0].parse().unwrap();
        let selected = ObservationEpoch {
            observations: e
                .observations
                .iter()
                .filter(|o| {
                    o.satellite == sat && (o.tracking_code == row[1] || o.tracking_code == row[2])
                })
                .cloned()
                .collect(),
            ..e.clone()
        };
        let table = SignalCorrections::new(correction_rows(&selected)).unwrap();
        let r = form_ionosphere_free(e, sat, &table, &config(&row)).unwrap();
        close(r.code_if_m, number(&row, 13), 0.001);
        let phase = r.phase.unwrap();
        close(phase.phase_if_m, number(&row, 14), 0.001);
        close(phase.raw_geometry_free_m, number(&row, 15), 0.0001);
        if let Some(mw) = phase.raw_melbourne_wubbena_m {
            close(mw, number(&row, 16), 0.001);
        }
        assert_eq!(r.raw[0].lli, Some(0));
        assert_eq!(r.raw[0].signal_strength, Some(7));
        count += 1;
    }
    assert_eq!(count, 48);
}

#[test]
fn additive_corrections_apply_before_if_and_preserve_raw_slip_diagnostics() {
    let row = &rows()[1];
    let e = epoch(row);
    let conf = config(row);
    let base = form(&e, correction_rows(&e), &conf).unwrap();
    let mut entries = correction_rows(&e);
    entries[0].code_add_m = -3.2;
    entries[1].code_add_m = 1.7;
    entries[0].phase_add_m = Some(0.023);
    entries[1].phase_add_m = Some(-0.047);
    let r = form(&e, entries.clone(), &conf).unwrap();
    let [a, b] = r.coefficients;
    close(r.code_if_m - base.code_if_m, a * (-3.2) + b * 1.7, 1e-7);
    let p = r.phase.unwrap();
    let bp = base.phase.unwrap();
    close(p.phase_if_m - bp.phase_if_m, a * 0.023 - b * 0.047, 1e-7);
    assert_eq!(p.raw_geometry_free_m, bp.raw_geometry_free_m);
    assert_eq!(p.raw_melbourne_wubbena_m, bp.raw_melbourne_wubbena_m);
    assert_eq!(r.raw, e.observations.as_slice());
    assert_eq!(r.corrections, entries.as_slice());
}

#[test]
fn correlated_frequency_noise_is_transformed_without_diagonalizing() {
    let row = &rows()[0];
    let e = epoch(row);
    let mut conf = config(row);
    conf.noise.code_m2 = [[0.25, 0.1], [0.1, 0.36]];
    conf.noise.phase_m2 = [[0.0001, -0.00002], [-0.00002, 0.0004]];
    let r = form(&e, correction_rows(&e), &conf).unwrap();
    let [a, b] = r.coefficients;
    close(
        r.code_variance_m2,
        a * a * 0.25 + 2.0 * a * b * 0.1 + b * b * 0.36,
        1e-14,
    );
    close(
        r.phase.unwrap().variance_m2,
        a * a * 0.0001 - 2.0 * a * b * 0.00002 + b * b * 0.0004,
        1e-17,
    );
    for covariance in [
        [[1.0, 0.2], [0.3, 1.0]],
        [[1.0, 2.0], [2.0, 1.0]],
        [[0.0, 0.0], [0.0, 1.0]],
        [[f64::NAN, 0.0], [0.0, 1.0]],
    ] {
        conf.noise.code_m2 = covariance;
        assert_eq!(
            form(&e, correction_rows(&e), &conf),
            Err(DualFrequencyError::InvalidNoise)
        );
    }
    conf.noise.code_m2 = [[0.25, 0.0], [0.0, 0.25]];
    conf.noise.phase_m2 = [[1.0, 1.0], [1.0, 1.0]];
    assert_eq!(
        form(&e, correction_rows(&e), &conf),
        Err(DualFrequencyError::InvalidNoise)
    );
}

#[test]
fn missing_half_cycle_or_uncorrected_phase_keeps_if_code_and_reports_reason() {
    let row = &rows()[0];
    let e = epoch(row);
    let conf = config(row);
    for i in 0..2 {
        let mut changed = e.clone();
        changed.observations[i].carrier_phase_cycles = None;
        let r = form(&changed, correction_rows(&e), &conf).unwrap();
        assert!(r.phase.is_none());
        close(r.code_if_m, number(row, 13), 1e-7);
        assert_eq!(
            r.phase_unavailability,
            Some(PhaseUnavailable::MissingPhase {
                tracking_code: conf.tracking_codes[i].clone()
            })
        );
        let mut changed = e.clone();
        changed.observations[i].lli = Some(2);
        let r = form(&changed, correction_rows(&e), &conf).unwrap();
        assert!(r.phase.is_none());
        assert_eq!(
            r.phase_unavailability,
            Some(PhaseUnavailable::HalfCycle {
                tracking_code: conf.tracking_codes[i].clone()
            })
        );
        let mut entries = correction_rows(&e);
        entries[i].phase_add_m = None;
        entries[i].phase_reference = None;
        let r = form(&e, entries, &conf).unwrap();
        assert!(r.phase.is_none());
        assert_eq!(
            r.phase_unavailability,
            Some(PhaseUnavailable::MissingCorrection {
                tracking_code: conf.tracking_codes[i].clone()
            })
        );
    }
}

#[test]
fn loss_of_lock_and_power_failure_require_new_phase_arc() {
    let row = &rows()[0];
    let e = epoch(row);
    let conf = config(row);
    for i in 0..2 {
        let mut changed = e.clone();
        changed.observations[i].lli = Some(1);
        let r = form(&changed, correction_rows(&e), &conf).unwrap();
        assert!(r.phase.unwrap().reset_required);
    }
    let mut changed = e.clone();
    changed.flag = 1;
    assert!(
        form(&changed, correction_rows(&e), &conf)
            .unwrap()
            .phase
            .unwrap()
            .reset_required
    );
    changed = e.clone();
    changed.observations[0].lli = Some(4);
    let r = form(&changed, correction_rows(&e), &conf).unwrap();
    assert!(!r.phase.unwrap().reset_required);
    assert_eq!(r.raw[0].lli, Some(4));
}

#[test]
fn exact_tracking_selection_never_uses_alternative_or_duplicate() {
    let row = &rows()[0];
    let e = epoch(row);
    let conf = config(row);
    let entries = correction_rows(&e);
    let mut changed = e.clone();
    changed.observations[1].tracking_code = "2L".to_owned();
    assert!(matches!(
        form(&changed, entries.clone(), &conf),
        Err(DualFrequencyError::MissingObservation { .. })
    ));
    changed = e.clone();
    changed.observations.push(changed.observations[1].clone());
    assert!(matches!(
        form(&changed, entries.clone(), &conf),
        Err(DualFrequencyError::DuplicateObservation { .. })
    ));
    changed = e.clone();
    changed.observations[1].pseudorange_m = None;
    assert!(matches!(
        form(&changed, entries, &conf),
        Err(DualFrequencyError::MissingCode { .. })
    ));
}

#[test]
fn corrections_require_current_matching_clock_and_phase_datums() {
    let row = &rows()[0];
    let e = epoch(row);
    let conf = config(row);
    for i in 0..2 {
        let mut entries = correction_rows(&e);
        entries.remove(i);
        assert!(matches!(
            form(&e, entries, &conf),
            Err(DualFrequencyError::MissingCorrection { .. })
        ));
        let mut entries = correction_rows(&e);
        entries[i].valid_from = e.time.checked_add_seconds(1.0).unwrap();
        assert!(matches!(
            form(&e, entries, &conf),
            Err(DualFrequencyError::ExpiredCorrection { .. })
        ));
        let mut entries = correction_rows(&e);
        entries[i].clock_source = ClockSource::Sp3;
        assert!(matches!(
            form(&e, entries, &conf),
            Err(DualFrequencyError::ClockDatumMismatch { .. })
        ));
        let mut entries = correction_rows(&e);
        entries[i].clock_reference = "other-clock".to_owned();
        assert!(matches!(
            form(&e, entries, &conf),
            Err(DualFrequencyError::ClockDatumMismatch { .. })
        ));
        let mut entries = correction_rows(&e);
        entries[i].phase_reference = Some("other-phase".to_owned());
        assert!(matches!(
            form(&e, entries, &conf),
            Err(DualFrequencyError::PhaseDatumMismatch { .. })
        ));
    }
}

#[test]
fn invalid_or_ambiguous_correction_ledgers_are_rejected() {
    let e = epoch(&rows()[0]);
    let entries = correction_rows(&e);
    for kind in 0..8 {
        let mut changed = entries.clone();
        match kind {
            0 => changed.push(changed[0].clone()),
            1 => changed[0].code_add_m = f64::NAN,
            2 => changed[0].phase_add_m = Some(f64::INFINITY),
            3 => changed[0].valid_until = changed[0].valid_from.checked_add_seconds(-1.0).unwrap(),
            4 => changed[0].clock_reference = " ".to_owned(),
            5 => changed[0].phase_reference = None,
            6 => changed[0].phase_reference = Some(" ".to_owned()),
            _ => changed[0].tracking_code = "0C".to_owned(),
        }
        assert!(matches!(
            SignalCorrections::new(changed),
            Err(DualFrequencyError::InvalidCorrections)
        ));
    }
    let mut changed = entries;
    let mut next = changed[0].clone();
    next.valid_from = next.valid_until;
    next.valid_until = next.valid_until.checked_add_seconds(30.0).unwrap();
    changed.push(next);
    assert!(matches!(
        SignalCorrections::new(changed.clone()),
        Err(DualFrequencyError::InvalidCorrections)
    ));
    changed[2].valid_from = changed[2].valid_from.checked_add_seconds(0.001).unwrap();
    assert!(SignalCorrections::new(changed).is_ok());
}

#[test]
fn glonass_fdma_requires_consistent_valid_channels_and_retains_provenance() {
    let row = &rows()[40];
    let e = epoch(row);
    let conf = config(row);
    assert_eq!(e.observations[0].satellite.system(), GnssSystem::Glonass);
    let mut changed = e.clone();
    changed.observations[1].glonass_channel = None;
    assert_eq!(
        form(&changed, correction_rows(&e), &conf),
        Err(DualFrequencyError::Gnss(Error::MissingGlonassChannel))
    );
    changed = e.clone();
    changed.observations[1].glonass_channel = Some(7);
    assert_eq!(
        form(&changed, correction_rows(&e), &conf),
        Err(DualFrequencyError::Gnss(Error::InvalidGlonassChannel))
    );
    changed = e.clone();
    changed.observations[1].glonass_channel = Some(0);
    assert_eq!(
        form(&changed, correction_rows(&e), &conf),
        Err(DualFrequencyError::InvalidFrequencyPair)
    );
}

#[test]
fn invalid_frequencies_tracking_and_nonfinite_measurements_do_not_fallback() {
    for [f1, f2] in [
        [c::GPS_L1_FREQ, c::GPS_L1_FREQ],
        [c::GPS_L1_FREQ, c::GPS_L1_FREQ + 0.5],
        [0.0, 1e9],
        [-1e9, 1e9],
        [f64::INFINITY, 1e9],
        [f64::MAX, 1e9],
    ] {
        assert_eq!(
            ionosphere_free_coefficients(f1, f2),
            Err(DualFrequencyError::InvalidFrequencyPair)
        );
    }
    let sat: SatelliteId = "G01".parse().unwrap();
    for code in ["1Q", "9C", "1", "β", "C1C"] {
        assert_eq!(
            tracking_frequency_hz(sat, code, None),
            Err(DualFrequencyError::UnsupportedTrackingCode)
        );
    }
    for sat in ["E01", "C01", "S01", "I01"] {
        assert!(tracking_frequency_hz(sat.parse().unwrap(), "8X", None).is_err());
    }
    let row = &rows()[0];
    let e = epoch(row);
    let conf = config(row);
    for value in [f64::NAN, f64::INFINITY, -1.0, 0.0] {
        let mut changed = e.clone();
        changed.observations[1].pseudorange_m = Some(value);
        assert!(matches!(
            form(&changed, correction_rows(&e), &conf),
            Err(DualFrequencyError::InvalidObservation { .. })
        ));
    }
    let mut changed = e.clone();
    changed.observations[1].carrier_phase_cycles = Some(f64::NAN);
    assert!(matches!(
        form(&changed, correction_rows(&e), &conf),
        Err(DualFrequencyError::InvalidObservation { .. })
    ));
    changed = e.clone();
    changed.observations[1].lli = Some(8);
    assert!(matches!(
        form(&changed, correction_rows(&e), &conf),
        Err(DualFrequencyError::InvalidObservation { .. })
    ));
    changed = e.clone();
    changed.observations[1].carrier_phase_cycles = Some(0.0);
    assert!(
        form(&changed, correction_rows(&e), &conf)
            .unwrap()
            .phase
            .is_some()
    );
}

#[test]
fn epochs_offsets_and_empty_datums_are_checked() {
    let row = &rows()[0];
    let e = epoch(row);
    let mut conf = config(row);
    for flag in [2, 6, 255] {
        let mut changed = e.clone();
        changed.flag = flag;
        assert_eq!(
            form(&changed, correction_rows(&e), &conf),
            Err(DualFrequencyError::InvalidEpoch)
        );
    }
    for offset in [1e-7, f64::NAN, f64::INFINITY] {
        let mut changed = e.clone();
        changed.receiver_clock_offset_s = Some(offset);
        assert_eq!(
            form(&changed, correction_rows(&e), &conf),
            Err(DualFrequencyError::InvalidEpoch)
        );
    }
    conf.clock_reference = " ".to_owned();
    assert_eq!(
        form(&e, correction_rows(&e), &conf),
        Err(DualFrequencyError::InvalidConfiguration)
    );
    conf = config(row);
    conf.phase_reference.clear();
    assert_eq!(
        form(&e, correction_rows(&e), &conf),
        Err(DualFrequencyError::InvalidConfiguration)
    );
}

#[test]
fn reversed_pair_preserves_if_and_mw_but_reverses_geometry_free() {
    let row = &rows()[0];
    let e = epoch(row);
    let conf = config(row);
    let base = form(&e, correction_rows(&e), &conf).unwrap();
    let mut reverse = conf;
    reverse.tracking_codes.swap(0, 1);
    let r = form(&e, correction_rows(&e), &reverse).unwrap();
    close(r.code_if_m, base.code_if_m, 1e-7);
    let p = r.phase.unwrap();
    let b = base.phase.unwrap();
    close(p.phase_if_m, b.phase_if_m, 1e-7);
    close(p.raw_geometry_free_m, -b.raw_geometry_free_m, 1e-7);
    close(
        p.raw_melbourne_wubbena_m.unwrap(),
        b.raw_melbourne_wubbena_m.unwrap(),
        1e-7,
    );
}

#[test]
fn correction_validity_uses_full_gps_time_across_week_boundary() {
    let row = &rows()[0];
    let mut e = epoch(row);
    let conf = config(row);
    let mut entries = correction_rows(&e);
    for entry in &mut entries {
        entry.valid_from = GnssTime::new(2300, 604790.0).unwrap();
        entry.valid_until = GnssTime::new(2301, 10.0).unwrap();
    }
    for time in [
        GnssTime::new(2300, 604790.0).unwrap(),
        GnssTime::new(2301, 0.0).unwrap(),
        GnssTime::new(2301, 10.0).unwrap(),
    ] {
        e.time = time;
        assert!(form(&e, entries.clone(), &conf).is_ok());
    }
    e.time = GnssTime::new(2302, 0.0).unwrap();
    assert!(matches!(
        form(&e, entries, &conf),
        Err(DualFrequencyError::ExpiredCorrection { .. })
    ));
}
