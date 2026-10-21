use gnss_rust::Error;
use gnss_rust::rtk_slip::{
    single_difference_code_minus_phase_m as cmc, single_difference_range_rate_mps as rate,
    test_code_minus_phase as code, test_doppler_phase as doppler,
};

#[test]
fn valid_single_frequency_slip_components_match_28_native_cases() {
    let mut count = 0;
    for line in include_str!("fixtures/upstream_rtk_slip.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split(',').collect();
        let v: Vec<f64> = fields[2..].iter().map(|s| s.parse().unwrap()).collect();
        let actual = match fields[1] {
            "rate" => rate(v[0], v[1], v[2]).unwrap(),
            "cmc" => cmc(v[0], v[1], v[2], v[3], v[4]).unwrap(),
            "code" => f64::from(code(v[0], v[1], v[2]).unwrap().detected),
            "doppler" => f64::from(doppler(v[0], v[1], v[2], v[3], v[4]).unwrap().detected),
            _ => panic!("unknown case"),
        };
        let expected = v[v.len() - 1];
        assert!(
            (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
            "{}: {actual} != {expected}",
            fields[0]
        );
        count += 1;
    }
    assert_eq!(count, 28);
}

#[test]
fn invalid_values_and_arithmetic_overflow_never_become_zero_or_negative_evidence() {
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for index in 0..3 {
            let mut v = [1.0; 3];
            v[index] = invalid;
            assert_eq!(rate(v[0], v[1], v[2]), Err(Error::NonFinite));
            assert_eq!(code(v[0], v[1], v[2]), Err(Error::NonFinite));
        }
        for index in 0..5 {
            let mut v = [1.0; 5];
            v[index] = invalid;
            assert_eq!(cmc(v[0], v[1], v[2], v[3], v[4]), Err(Error::NonFinite));
            assert_eq!(doppler(v[0], v[1], v[2], v[3], v[4]), Err(Error::NonFinite));
        }
    }
    for invalid in [0.0, -1.0] {
        assert_eq!(rate(0.0, 0.0, invalid), Err(Error::OutOfRange));
        assert_eq!(cmc(0.0, 0.0, 0.0, 0.0, invalid), Err(Error::OutOfRange));
        assert_eq!(code(0.0, 0.0, invalid), Err(Error::OutOfRange));
        assert_eq!(doppler(0.0, 0.0, 0.0, invalid, 1.0), Err(Error::OutOfRange));
        assert_eq!(doppler(0.0, 0.0, 0.0, 1.0, invalid), Err(Error::OutOfRange));
    }
    assert_eq!(rate(f64::MAX, -f64::MAX, 1.0), Err(Error::NonFinite));
    assert_eq!(
        cmc(f64::MAX, -f64::MAX, 0.0, 0.0, 1.0),
        Err(Error::NonFinite)
    );
    assert_eq!(code(f64::MAX, -f64::MAX, 1.0), Err(Error::NonFinite));
    assert_eq!(doppler(0.0, 0.0, f64::MAX, 2.0, 1.0), Err(Error::NonFinite));
}

#[test]
fn signed_residuals_zero_measurements_and_strict_boundary_remain_explicit() {
    assert_eq!(rate(0.0, 0.0, 0.19), Ok(0.0));
    assert_eq!(cmc(0.0, 0.0, 0.0, 0.0, 0.19), Ok(0.0));
    assert_eq!(code(0.0, -5.0, 5.0).unwrap().residual_m, -5.0);
    assert!(!code(0.0, -5.0, 5.0).unwrap().detected);
    assert!(
        code(0.0, f64::from_bits(5.0_f64.to_bits() + 1), 5.0)
            .unwrap()
            .detected
    );
    let test = doppler(10.0, 9.0, 0.0, 1.0, 0.2).unwrap();
    assert_eq!(test.residual_m, -1.0);
    assert!(test.detected);
    assert!(!doppler(10.0, 9.0, -1.0, 1.0, 0.2).unwrap().detected);
}
