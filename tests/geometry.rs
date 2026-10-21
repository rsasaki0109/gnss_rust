use gnss_rust::geometry::dilution_of_precision;
use gnss_rust::{EcefCoord, Error};

#[test]
fn common_clock_geometry_matches_native_dop_and_eigen_covariance() {
    let mut count = 0;
    for line in include_str!("fixtures/upstream_geometry.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let mut fields = line.split(',');
        let name = fields.next().unwrap();
        let satellites: usize = fields.next().unwrap().parse().unwrap();
        let values: Vec<f64> = fields.map(|v| v.parse().unwrap()).collect();
        let receiver = EcefCoord::new(values[0], values[1], values[2]);
        let mut data = Vec::new();
        for index in 0..satellites {
            let k = 3 + 3 * index;
            data.push((
                format!("G{:02}", index + 1).parse().unwrap(),
                EcefCoord::new(values[k], values[k + 1], values[k + 2]),
            ));
        }
        let expected = &values[3 + 3 * satellites..];
        let result = dilution_of_precision(receiver, &data).unwrap();
        let tolerance = if name.starts_with("weak") {
            1e-6
        } else {
            1e-10
        };
        let close = |a: f64, b: f64| {
            assert!(
                (a - b).abs() <= tolerance * b.abs().max(1.0),
                "{name}: {a} vs {b}"
            )
        };
        close(result.gdop, expected[0]);
        close(result.pdop, expected[1]);
        close(result.tdop, expected[2]);
        for i in 0..4 {
            for j in 0..4 {
                close(result.geometry_covariance[i][j], expected[3 + 4 * i + j]);
            }
        }
        data.reverse();
        let reversed = dilution_of_precision(receiver, &data).unwrap();
        assert!((result.gdop - reversed.gdop).abs() <= 1e-9 * result.gdop);
        assert!(result.gdop > result.pdop && result.gdop > result.tdop);
        count += 1;
    }
    assert_eq!(count, 16);
}

#[test]
fn invalid_and_rank_deficient_geometry_fails_without_regularization() {
    let receiver = EcefCoord::new(1e6, 2e6, 3e6);
    let points: Vec<_> = (1..=4)
        .map(|id| {
            (
                format!("G{id:02}").parse().unwrap(),
                EcefCoord::new(2e7, 1e7, 3e7),
            )
        })
        .collect();
    assert_eq!(
        dilution_of_precision(receiver, &points),
        Err(Error::SingularGeometry)
    );
    assert!(dilution_of_precision(receiver, &points[..3]).is_err());
    let mut duplicate = points.clone();
    duplicate[1].0 = duplicate[0].0;
    assert_eq!(
        dilution_of_precision(receiver, &duplicate),
        Err(Error::InvalidConfiguration)
    );
    let mut coincident = points.clone();
    coincident[0].1 = receiver;
    assert_eq!(
        dilution_of_precision(receiver, &coincident),
        Err(Error::SingularGeometry)
    );
    assert_eq!(
        dilution_of_precision(EcefCoord::new(f64::NAN, 0.0, 0.0), &points),
        Err(Error::NonFinite)
    );
}
