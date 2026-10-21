//! Native PVA motion-recovery primitives; caller admits seed provenance/time/frame.
//! These functions do not obtain SPP/Doppler observations or update a PPP filter.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
use crate::EcefCoord;
use crate::filter::{
    FilterError, FilterUpdate, GaussianState, cholesky, measurement_update, solve_cholesky, zeros,
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PvaRecoveryConfig {
    /// Both the mean XYZ variance trigger and independent reset XYZ variance.
    pub position_variance_m2: f64,
    pub velocity_variance_m2_s2: f64,
    pub acceleration_variance_m2_s4: f64,
    /// Total three-dimensional NIS, not NIS per row. None disables this trigger.
    pub position_disagreement_gate: Option<f64>,
}
impl Default for PvaRecoveryConfig {
    fn default() -> Self {
        Self {
            position_variance_m2: 3600.0,
            velocity_variance_m2_s2: 100.0,
            acceleration_variance_m2_s4: 100.0,
            position_disagreement_gate: Some(16.27),
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MotionSeed {
    pub position_ecef: EcefCoord,
    /// None/nonfinite uses native 100 m² identity only for the disagreement test.
    /// It is not fused as an independent position measurement.
    pub position_covariance_m2: Option<[[f64; 3]; 3]>,
    /// Missing/nonfinite velocity starts a reset from rest, without a measurement.
    pub velocity_mps: Option<[f64; 3]>,
}
#[derive(Debug, Clone, PartialEq)]
pub struct MotionRecovery {
    pub prior: GaussianState,
    pub reset: bool,
    pub variance_trigger: bool,
    pub normalized_position_disagreement: Option<f64>,
    pub seed_covariance_fallback: bool,
}

/// d' (P + C)^-1 d. Requires a finite, symmetric positive-definite sum.
/// Native failure returns infinity; Rust exposes its checked error explicitly.
pub fn position_disagreement(
    predicted: EcefCoord,
    predicted_covariance: [[f64; 3]; 3],
    seed: EcefCoord,
    seed_covariance: [[f64; 3]; 3],
) -> Result<f64, FilterError> {
    if !predicted.is_finite() || !seed.is_finite() {
        return Err(FilterError::NonFinite);
    }
    let sum = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| predicted_covariance[i][j] + seed_covariance[i][j])
                .collect()
        })
        .collect();
    let d = [
        predicted.x - seed.x,
        predicted.y - seed.y,
        predicted.z - seed.z,
    ];
    let solved = solve_cholesky(&cholesky(&sum)?, &d);
    let nis = d.iter().zip(solved).map(|(a, b)| a * b).sum::<f64>();
    if !nis.is_finite() || nis < 0.0 {
        return Err(FilterError::NumericalFailure);
    }
    Ok(nis)
}

/// Apply the native PVA variance/disagreement restart to a copied predicted
/// Rust [XYZ,clock,ZTD,VXYZ,AXYZ,IFN...] state. Remaining marginal covariance is
/// retained, while all motion rows/columns are independently reinitialized.
/// This is a data-dependent restart, not an independent seed-measurement fusion.
pub fn recover_pva_motion(
    prior: &GaussianState,
    seed: Option<MotionSeed>,
    config: PvaRecoveryConfig,
) -> Result<MotionRecovery, FilterError> {
    if prior.mean.len() < 11 || prior.covariance.len() != prior.mean.len() {
        return Err(FilterError::InvalidDimensions);
    }
    if !prior.mean.iter().all(|x| x.is_finite()) {
        return Err(FilterError::NonFinite);
    }
    cholesky(&prior.covariance)?;
    if [
        config.position_variance_m2,
        config.velocity_variance_m2_s2,
        config.acceleration_variance_m2_s4,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x <= 0.0)
        || config
            .position_disagreement_gate
            .is_some_and(|x| !x.is_finite() || x <= 0.0)
    {
        return Err(FilterError::InvalidGate);
    }
    let mut result = MotionRecovery {
        prior: prior.clone(),
        reset: false,
        variance_trigger: false,
        normalized_position_disagreement: None,
        seed_covariance_fallback: false,
    };
    let Some(seed) = seed else {
        return Ok(result);
    };
    if !seed.position_ecef.is_finite() {
        return Err(FilterError::NonFinite);
    }
    result.variance_trigger =
        (prior.covariance[0][0] + prior.covariance[1][1] + prior.covariance[2][2]) / 3.0
            > config.position_variance_m2;
    if !result.variance_trigger
        && let Some(gate) = config.position_disagreement_gate
    {
        let selected = seed
            .position_covariance_m2
            .filter(|c| c.iter().flatten().all(|x| x.is_finite()));
        result.seed_covariance_fallback = selected.is_none();
        let seed_covariance =
            selected.unwrap_or([[100.0, 0.0, 0.0], [0.0, 100.0, 0.0], [0.0, 0.0, 100.0]]);
        let nis = position_disagreement(
            EcefCoord::new(prior.mean[0], prior.mean[1], prior.mean[2]),
            std::array::from_fn(|i| std::array::from_fn(|j| prior.covariance[i][j])),
            seed.position_ecef,
            seed_covariance,
        )
        .unwrap_or(f64::INFINITY);
        result.normalized_position_disagreement = Some(nis);
        result.reset = nis > gate;
    }
    result.reset |= result.variance_trigger;
    if result.reset {
        let position = [
            seed.position_ecef.x,
            seed.position_ecef.y,
            seed.position_ecef.z,
        ];
        let velocity = seed
            .velocity_mps
            .filter(|v| v.iter().all(|x| x.is_finite()))
            .unwrap_or([0.0; 3]);
        for (start, mean, variance) in [
            (0, position, config.position_variance_m2),
            (5, velocity, config.velocity_variance_m2_s2),
            (8, [0.0; 3], config.acceleration_variance_m2_s4),
        ] {
            for (j, value) in mean.into_iter().enumerate() {
                let i = start + j;
                result.prior.mean[i] = value;
                for row in &mut result.prior.covariance {
                    row[i] = 0.0;
                }
                result.prior.covariance[i].fill(0.0);
                result.prior.covariance[i][i] = variance;
            }
        }
    }
    Ok(result)
}

/// Direct velocity measurement on the full Gaussian state. Uses checked SPD
/// noise and Joseph covariance update; native uses a subtractive update.
/// Gate is total 3-D NIS; None reproduces the native ungated normal path.
pub fn update_velocity(
    prior: &GaussianState,
    velocity_index: usize,
    velocity: [f64; 3],
    covariance_m2_s2: [[f64; 3]; 3],
    gate: Option<f64>,
) -> Result<FilterUpdate, FilterError> {
    let n = prior.mean.len();
    if velocity_index.checked_add(3).is_none_or(|end| end > n) {
        return Err(FilterError::InvalidDimensions);
    }
    let mut design = zeros(3, n);
    for (i, row) in design.iter_mut().enumerate() {
        row[velocity_index + i] = 1.0;
    }
    let residuals: Vec<_> = (0..3)
        .map(|i| velocity[i] - prior.mean[velocity_index + i])
        .collect();
    let noise = covariance_m2_s2.iter().map(|r| r.to_vec()).collect();
    measurement_update(prior, &design, &residuals, &noise, gate.map(|x| x / 3.0))
}
