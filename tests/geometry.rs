use gnss_rust::geometry::{dilution_of_precision, dilution_of_precision_by_system};
use gnss_rust::{EcefCoord, Error};

#[test]
fn per_system_clock_geometry_matches_24_independent_eigen_cases() {
    let mut count = 0;
    for line in include_str!("fixtures/upstream_multiclock_geometry.csv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split(',').collect();
        let groups: usize = fields[1].parse().unwrap();
        let per: usize = fields[2].parse().unwrap();
        let m = groups * per;
        let n = 3 + groups;
        let v: Vec<f64> = fields[3..].iter().map(|s| s.parse().unwrap()).collect();
        assert_eq!(v.len(), 3 + 3 * m + n * n);
        let receiver = EcefCoord::new(v[0], v[1], v[2]);
        let data: Vec<_> = (0..m)
            .map(|i| {
                let prefix = ["G", "E", "J"][i / per];
                let k = 3 + 3 * i;
                (
                    format!("{prefix}{:02}", i % per + 1).parse().unwrap(),
                    EcefCoord::new(v[k], v[k + 1], v[k + 2]),
                )
            })
            .collect();
        let result = dilution_of_precision_by_system(receiver, &data).unwrap();
        assert_eq!(result.clock_systems.len(), groups);
        let actual = result.full_geometry_covariance.clone().unwrap_or_else(|| {
            result
                .geometry_covariance
                .iter()
                .map(|r| r.to_vec())
                .collect()
        });
        let expected = &v[3 + 3 * m..];
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (actual[i][j] - expected[i * n + j]).abs()
                        < 1e-10 * expected[i * n + j].abs().max(1.0),
                    "{} Q[{i},{j}]",
                    fields[0]
                );
            }
        }
        let p: f64 = (0..3).map(|i| expected[i * n + i]).sum();
        let t: f64 = (3..n).map(|i| expected[i * n + i]).sum();
        assert!((result.pdop - p.sqrt()).abs() < 1e-10 * p.sqrt().max(1.0));
        assert!((result.tdop - t.sqrt()).abs() < 1e-10 * t.sqrt().max(1.0));
        let mut reversed = data.clone();
        reversed.reverse();
        let reverse = dilution_of_precision_by_system(receiver, &reversed).unwrap();
        assert!((result.gdop - reverse.gdop).abs() < 1e-9 * result.gdop);
        if groups > 1 {
            assert!(result.pdop >= dilution_of_precision(receiver, &data).unwrap().pdop);
        }
        count += 1;
    }
    assert_eq!(count, 24);
}

#[test]
fn per_system_geometry_requires_enough_independent_rows_without_a_prior() {
    let receiver = EcefCoord::new(1e6, 2e6, 3e6);
    let data: Vec<_> = [
        ("G01", [2e7, 1e7, 3e7]),
        ("G02", [-1e7, 2e7, 3e7]),
        ("G03", [1e7, -2e7, 3e7]),
        ("E01", [2e7, 3e7, -1e7]),
    ]
    .into_iter()
    .map(|(s, p)| (s.parse().unwrap(), EcefCoord::new(p[0], p[1], p[2])))
    .collect();
    assert_eq!(
        dilution_of_precision_by_system(receiver, &data),
        Err(Error::InvalidConfiguration)
    );
    assert_eq!(
        dilution_of_precision_by_system(receiver, &[]),
        Err(Error::InvalidConfiguration)
    );
    let mut rank_deficient = data.clone();
    rank_deficient.extend([
        ("G04".parse().unwrap(), data[0].1),
        ("E02".parse().unwrap(), data[0].1),
    ]);
    for (_, p) in &mut rank_deficient {
        *p = data[0].1;
    }
    assert_eq!(
        dilution_of_precision_by_system(receiver, &rank_deficient),
        Err(Error::SingularGeometry)
    );
    let mut duplicate = rank_deficient.clone();
    duplicate[1].0 = duplicate[0].0;
    assert_eq!(
        dilution_of_precision_by_system(receiver, &duplicate),
        Err(Error::InvalidConfiguration)
    );
}

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
