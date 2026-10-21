use super::*;

#[test]
fn locked_marginal_matches_native_selection_subset_and_independent_schur_covariance() {
    let mut cases = 0;
    for line in include_str!("../tests/fixtures/upstream_rtk_locked_subset.csv")
        .lines()
        .skip(1)
    {
        let mut values = line.split(',').map(|s| s.parse::<f64>().unwrap());
        let id = values.next().unwrap() as usize;
        let count = values.next().unwrap() as usize;
        let n = values.next().unwrap() as usize;
        let k = values.next().unwrap() as usize;
        let locks: Vec<_> = values.by_ref().take(count).map(|v| v as usize).collect();
        let elevations: Vec<_> = values.by_ref().take(count).collect();
        let mean: Vec<_> = values.by_ref().take(n).collect();
        let covariance: Matrix = (0..n).map(|_| values.by_ref().take(n).collect()).collect();
        let expected_pairs: Vec<_> = (0..k)
            .map(|_| {
                (
                    values.next().unwrap() as usize,
                    values.next().unwrap() as usize,
                )
            })
            .collect();
        let expected_dd: Vec<_> = values.by_ref().take(k).collect();
        let expected_qbb: Matrix = (0..k).map(|_| values.by_ref().take(k).collect()).collect();
        let expected_qab: Matrix = (0..3).map(|_| values.by_ref().take(k).collect()).collect();
        let data: Vec<_> = (0..count)
            .map(|i| SatelliteData {
                satellite: SatelliteId::new(
                    if i < 5 {
                        GnssSystem::Gps
                    } else {
                        GnssSystem::Galileo
                    },
                    if i < 5 { i + 1 } else { i - 4 } as u8,
                )
                .unwrap(),
                rover_position: EcefCoord::new(20e6, 10e6, 15e6),
                base_position: EcefCoord::new(20e6, 10e6, 15e6),
                code_sd_m: 0.0,
                phase_sd_cycles: 0.0,
                elevation: elevations[i],
                loss_of_lock: false,
                rover_tracking_code: "1C",
                base_tracking_code: "1C",
            })
            .collect();
        let lock_map = data
            .iter()
            .zip(locks)
            .map(|(d, n)| (d.satellite, n))
            .collect();
        let context = FixContext {
            state: GaussianState { mean, covariance },
            data,
            base_position: EcefCoord::new(4e6, 3e6, 4e6),
            config: RtkFloatConfig {
                constellations: RtkConstellations::GpsGalileo,
                ..Default::default()
            },
            time: GnssTime::new(2325, 100.0).unwrap(),
            reset_satellites: Vec::new(),
            reference: 4,
            references: vec![4, 8],
        };
        let before = context.clone();
        let subset = context.locked_subset(&lock_map, 5).unwrap();
        assert_eq!(context, before);
        cases += 1;
        if k == 0 {
            assert!(subset.is_none(), "case {id}");
            assert!(values.next().is_none());
            continue;
        }
        let subset = subset.unwrap();
        let selected_pairs: Vec<_> = dd_pairs(&subset.data, &subset.references)
            .iter()
            .map(|&(r, s)| (subset.data[r].satellite, subset.data[s].satellite))
            .collect();
        assert_eq!(
            selected_pairs,
            expected_pairs
                .iter()
                .map(|&(r, s)| (context.data[r].satellite, context.data[s].satellite))
                .collect::<Vec<_>>(),
            "case {id}"
        );
        let a = subset.ambiguities().unwrap();
        let close = |actual: f64, expected: f64| {
            assert!(
                (actual - expected).abs() <= 1e-11 * expected.abs().max(1.0),
                "case {id}: {actual} vs {expected}"
            )
        };
        for (&a, &b) in a.dd_float_cycles.iter().zip(&expected_dd) {
            close(a, b);
        }
        for (&a, &b) in a
            .ambiguity_covariance_cycles2
            .iter()
            .flatten()
            .zip(expected_qbb.iter().flatten())
        {
            close(a, b);
        }
        for (&a, &b) in a
            .head_ambiguity_covariance
            .iter()
            .flatten()
            .zip(expected_qab.iter().flatten())
        {
            close(a, b);
        }
        assert!(subset.ambiguity_pair_subset(&[]).is_err());
        assert!(subset.ambiguity_pair_subset(&[k]).is_err());
        assert!(subset.ambiguity_pair_subset(&[0, 0]).is_err());
        if k > 1 {
            let view = subset
                .ambiguity_pair_subset(&(0..k - 1).collect::<Vec<_>>())
                .unwrap()
                .unwrap();
            let reduced = view.ambiguities().unwrap();
            for i in 0..k - 1 {
                close(reduced.dd_float_cycles[i], expected_dd[i]);
                for (j, &value) in expected_qbb[i].iter().take(k - 1).enumerate() {
                    close(reduced.ambiguity_covariance_cycles2[i][j], value);
                }
                for (j, row) in expected_qab.iter().enumerate() {
                    close(reduced.head_ambiguity_covariance[j][i], row[i]);
                }
            }
            assert_eq!(subset.head_prior(), view.head_prior());
        }
        let integers: Vec<_> = expected_dd.iter().map(|v| v.round() as i64).collect();
        let fixed = crate::rtk_fix::condition_head(
            &subset.head_prior(),
            &a.dd_float_cycles,
            &a.ambiguity_covariance_cycles2,
            &a.head_ambiguity_covariance,
            &integers,
        )
        .unwrap();
        for &v in &fixed.state.mean {
            close(v, values.next().unwrap());
        }
        for &v in fixed.state.covariance.iter().flatten() {
            close(v, values.next().unwrap());
        }
        assert!(values.next().is_none());
    }
    assert_eq!(cases, 24);
}
