use gnss_rust::GnssTime;
use gnss_rust::epoch_match::{EpochMatch, EpochMatcher, EpochSide, MatchError};
use gnss_rust::observation::{Observation, ObservationEpoch};

fn epoch(tow: f64) -> ObservationEpoch {
    let mut observation = Observation::new("G01".parse().unwrap(), "1C".into());
    observation.carrier_phase_cycles = Some(123.5);
    ObservationEpoch {
        time: GnssTime::new(2300, tow).unwrap(),
        flag: 0,
        approximate_position: None,
        receiver_clock_offset_s: None,
        observations: vec![observation],
    }
}
fn stream(
    epochs: Vec<ObservationEpoch>,
) -> impl Iterator<Item = Result<ObservationEpoch, &'static str>> {
    epochs.into_iter().map(Ok)
}

#[test]
fn five_to_one_rates_preserve_every_epoch_tail_and_unmatched_continuity_evidence() {
    let mut rover: Vec<_> = (0..=10).map(|i| epoch(100.0 + i as f64 * 0.2)).collect();
    rover[2].observations[0].lli = Some(7);
    rover[3].flag = 6;
    rover[4].flag = 1;
    let base = vec![
        epoch(99.0),
        epoch(100.0),
        epoch(101.0),
        epoch(102.0),
        epoch(103.0),
    ];
    let events = EpochMatcher::new(stream(rover.clone()), stream(base.clone()))
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    let (mut recovered_rover, mut recovered_base, mut matched) = (Vec::new(), Vec::new(), 0);
    for event in events {
        match event {
            EpochMatch::Pair { rover, base } => {
                assert_eq!(rover.time, base.time);
                recovered_rover.push(rover);
                recovered_base.push(base);
                matched += 1;
            }
            EpochMatch::Unmatched {
                side: EpochSide::Rover,
                epoch,
            } => recovered_rover.push(epoch),
            EpochMatch::Unmatched {
                side: EpochSide::Base,
                epoch,
            } => recovered_base.push(epoch),
        }
    }
    assert_eq!(matched, 3);
    assert_eq!(recovered_rover, rover);
    assert_eq!(recovered_base, base);
}

#[test]
fn matching_uses_full_gps_week_and_a_fixed_small_time_tolerance() {
    let rover = vec![epoch(604799.8), epoch(604800.0), epoch(604800.2)];
    let base = vec![epoch(604800.0)];
    let events = EpochMatcher::new(stream(rover), stream(base))
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    assert!(
        matches!(&events[1], EpochMatch::Pair { rover, base } if rover.time.week()==2301 && rover.time == base.time)
    );
    for (difference, expected_pairs) in [(0.5e-7, 1), (2e-7, 0)] {
        let events = EpochMatcher::new(
            stream(vec![epoch(100.0)]),
            stream(vec![epoch(100.0 + difference)]),
        )
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
        assert_eq!(
            events
                .iter()
                .filter(|e| matches!(e, EpochMatch::Pair { .. }))
                .count(),
            expected_pairs
        );
    }
}

#[test]
fn duplicates_backwards_times_and_input_errors_are_side_specific_and_fused() {
    for side in [EpochSide::Rover, EpochSide::Base] {
        for final_time in [100.0, 99.0] {
            let bad = vec![epoch(100.0), epoch(final_time)];
            let good = vec![epoch(100.0), epoch(101.0)];
            let (r, b) = if side == EpochSide::Rover {
                (bad, good)
            } else {
                (good, bad)
            };
            let mut matcher = EpochMatcher::new(stream(r), stream(b));
            assert!(matcher.next().unwrap().is_ok());
            assert!(
                matches!(matcher.next().unwrap(), Err(MatchError::OutOfOrder { side: actual, .. }) if actual == side)
            );
            assert!(matcher.next().is_none() && matcher.next().is_none());
        }
        let good = vec![Ok(epoch(100.0))];
        let bad = vec![Err("parse failure")];
        let (r, b) = if side == EpochSide::Rover {
            (bad, good)
        } else {
            (good, bad)
        };
        let mut matcher = EpochMatcher::new(r.into_iter(), b.into_iter());
        assert!(
            matches!(matcher.next().unwrap(), Err(MatchError::Input { side: actual, source: "parse failure" }) if actual == side)
        );
        assert!(matcher.next().is_none());
    }
    assert_eq!(EpochMatcher::new(stream(vec![]), stream(vec![])).count(), 0);
    assert_eq!(
        EpochMatcher::new(stream(vec![epoch(100.0)]), stream(vec![])).count(),
        1
    );
}
