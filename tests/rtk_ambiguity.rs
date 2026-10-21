use gnss_rust::rtk_ambiguity::{DatumAlignmentInput, align_single_difference_datums as align};
use gnss_rust::{Error, GnssSystem, SatelliteId};
use std::collections::BTreeMap;

fn input(satellite: &str, raw: f64, old: Option<f64>) -> DatumAlignmentInput {
    DatumAlignmentInput {
        satellite: satellite.parse().unwrap(),
        raw_bias_cycles: raw,
        retained_cycles: old,
    }
}

#[test]
fn alignment_matches_the_actual_native_update_bias_in_nine_cases() {
    let mut cases = BTreeMap::<&str, Vec<Vec<&str>>>::new();
    for line in include_str!("fixtures/upstream_rtk_datum.csv")
        .lines()
        .skip(1)
    {
        let fields: Vec<_> = line.split(',').collect();
        if fields[0] != "native_zero_sentinel" {
            cases.entry(fields[0]).or_default().push(fields);
        }
    }
    assert_eq!(cases.len(), 9);
    for (name, rows) in cases {
        let inputs: Vec<_> = rows
            .iter()
            .map(|r| {
                input(
                    r[1],
                    r[2].parse().unwrap(),
                    if r[3].is_empty() || r[4] == "1" {
                        None
                    } else {
                        Some(r[3].parse().unwrap())
                    },
                )
            })
            .collect();
        let actual = align(&inputs).unwrap();
        for (i, row) in rows.iter().enumerate() {
            let expected: f64 = row[5].parse().unwrap();
            assert!(
                (actual.single_difference_cycles[i] - expected).abs()
                    < 1e-12 * expected.abs().max(1.0),
                "{name} {}",
                row[1]
            );
            // The native means are translated, while retained covariance and
            // new/slipped-arc independence stay unchanged.
            let retained = inputs[i].retained_cycles.is_some();
            assert_eq!(
                row[6].parse::<f64>().unwrap(),
                if retained { 4.0 } else { 900.0 }
            );
            let peers = inputs
                .iter()
                .enumerate()
                .filter(|(j, v)| *j != i && v.retained_cycles.is_some())
                .count();
            assert_eq!(
                row[7].parse::<f64>().unwrap(),
                if retained { 0.125 * peers as f64 } else { 0.0 }
            );
        }
    }
}

#[test]
fn zero_is_a_retained_value_instead_of_the_native_uninitialized_sentinel() {
    let result = align(&[
        input("G01", 10.0, Some(0.0)),
        input("G02", 5.0, Some(2.0)),
        input("G03", 8.0, None),
    ])
    .unwrap();
    assert_eq!(result.single_difference_cycles, vec![6.5, 8.5, 8.0]);
    assert_eq!(result.systems[0].retained_satellites.len(), 2);
    assert_eq!(
        result.systems[0].initialized_satellites,
        vec!["G03".parse::<SatelliteId>().unwrap()]
    );
    let native: Vec<_> = include_str!("fixtures/upstream_rtk_datum.csv")
        .lines()
        .filter(|r| r.starts_with("native_zero_sentinel,"))
        .map(|r| r.split(',').nth(5).unwrap().parse::<f64>().unwrap())
        .collect();
    assert_eq!(native, vec![10.0, 5.0, 8.0]);
}

#[test]
fn per_system_shifts_preserve_every_retained_double_difference() {
    let inputs = vec![
        input("E01", -50.0, Some(1.0)),
        input("G01", 15.0, Some(2.0)),
        input("G02", 16.0, Some(-4.0)),
        input("E02", -48.0, Some(5.0)),
        input("E03", 22.0, None),
    ];
    let result = align(&inputs).unwrap();
    assert_eq!(
        result.systems.iter().map(|d| d.system).collect::<Vec<_>>(),
        vec![GnssSystem::Gps, GnssSystem::Galileo]
    );
    for (i, a) in inputs.iter().enumerate() {
        for (j, b) in inputs.iter().enumerate() {
            if a.satellite.system() == b.satellite.system() {
                if let (Some(x), Some(y)) = (a.retained_cycles, b.retained_cycles) {
                    assert_eq!(
                        result.single_difference_cycles[i] - result.single_difference_cycles[j],
                        x - y
                    );
                }
            }
        }
    }
    let mut reversed = inputs.clone();
    reversed.reverse();
    let mut values = align(&reversed).unwrap().single_difference_cycles;
    values.reverse();
    assert_eq!(result.single_difference_cycles, values);
    assert!(align(&[]).unwrap().systems.is_empty());
}

#[test]
fn malformed_inputs_and_overflow_are_rejected_without_changing_inputs() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(align(&[input("G01", value, None)]), Err(Error::NonFinite));
        assert_eq!(
            align(&[input("G01", 0.0, Some(value))]),
            Err(Error::NonFinite)
        );
    }
    let duplicate = input("G01", 1.0, Some(0.0));
    assert_eq!(
        align(&[duplicate, duplicate]),
        Err(Error::InvalidConfiguration)
    );
    assert_eq!(
        align(&[input("R01", 1.0, None)]),
        Err(Error::UnsupportedConstellation)
    );
    let overflowing = [input("G01", f64::MAX, Some(-f64::MAX))];
    assert_eq!(align(&overflowing), Err(Error::NonFinite));
    assert_eq!(overflowing[0].retained_cycles, Some(-f64::MAX));
}
