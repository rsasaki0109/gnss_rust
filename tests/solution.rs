use gnss_rust::{EcefCoord, PositionSolution, SolutionStatus};

#[test]
fn absent_solution_state_cannot_be_mistaken_for_a_position() {
    let solution = PositionSolution::default();
    assert_eq!(solution.status, SolutionStatus::None);
    assert!(solution.position_ecef.is_none());
    assert!(solution.position_covariance.is_none());
    assert!(solution.velocity_ecef.is_none());
    assert!(solution.velocity_covariance.is_none());
    assert!(!solution.is_valid());
    assert!(!solution.is_fixed());
}

#[test]
fn propagated_position_needs_no_current_satellites_and_is_never_fixed() {
    let mut solution = PositionSolution {
        status: SolutionStatus::Propagated,
        position_ecef: Some(EcefCoord::new(1.0, 2.0, 3.0)),
        ..PositionSolution::default()
    };
    assert!(solution.is_valid());
    assert!(!solution.is_fixed());
    solution.position_ecef = Some(EcefCoord::new(f64::NAN, 2.0, 3.0));
    assert!(!solution.is_valid());
}

#[test]
fn gnss_solution_requires_finite_position_and_at_least_four_satellites() {
    for status in [
        SolutionStatus::Spp,
        SolutionStatus::Dgps,
        SolutionStatus::Float,
        SolutionStatus::Fixed,
        SolutionStatus::PppFloat,
        SolutionStatus::PppFixed,
    ] {
        let mut solution = PositionSolution {
            status,
            num_satellites: 4,
            ..PositionSolution::default()
        };
        assert!(!solution.is_valid());
        solution.position_ecef = Some(EcefCoord::new(1.0, 2.0, 3.0));
        assert!(solution.is_valid());
        solution.num_satellites = 3;
        assert!(!solution.is_valid());
        solution.num_satellites = 4;
        solution.position_ecef = Some(EcefCoord::new(1.0, f64::INFINITY, 3.0));
        assert!(!solution.is_valid());
    }
}

#[test]
fn fixed_status_is_distinct_from_solution_validity() {
    for status in [SolutionStatus::Fixed, SolutionStatus::PppFixed] {
        let solution = PositionSolution {
            status,
            ..PositionSolution::default()
        };
        assert!(solution.is_fixed());
        assert!(!solution.is_valid());
    }
    let none = PositionSolution {
        position_ecef: Some(EcefCoord::new(1.0, 2.0, 3.0)),
        num_satellites: 4,
        ..PositionSolution::default()
    };
    assert!(!none.is_valid());
}
