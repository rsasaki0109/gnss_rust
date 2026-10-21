//! Absolute prefit residual screening for optional RTK epoch admission.
//!
//! The strict row classification matches pinned native suppressOutlierRows.
//! An epoch gate rejects the update instead of zeroing H/v while retaining
//! correlated R. No residual, design matrix or covariance is mutated here.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.

use crate::Error;

#[derive(Debug, Clone, PartialEq)]
pub struct PrefitResidualScreen {
    pub threshold_m: f64,
    pub maximum_phase_abs_m: f64,
    pub maximum_code_abs_m: f64,
    pub rejected_phase_rows: Vec<usize>,
    pub rejected_code_rows: Vec<usize>,
}

impl PrefitResidualScreen {
    pub fn accepted(&self) -> bool {
        self.rejected_phase_rows.is_empty() && self.rejected_code_rows.is_empty()
    }
}

/// Compare raw measurement-minus-prediction residuals at the initial seed.
/// Equality passes; any strict exceedance rejects the entire epoch update.
pub fn screen_prefit_residuals(
    phase_m: &[f64],
    code_m: &[f64],
    threshold_m: f64,
) -> Result<PrefitResidualScreen, Error> {
    if phase_m.is_empty() || phase_m.len() != code_m.len() {
        return Err(Error::InsufficientObservations);
    }
    if !threshold_m.is_finite() || !phase_m.iter().chain(code_m).all(|v| v.is_finite()) {
        return Err(Error::NonFinite);
    }
    if threshold_m <= 0.0 {
        return Err(Error::OutOfRange);
    }
    let rejected = |rows: &[f64]| {
        rows.iter()
            .enumerate()
            .filter_map(|(i, v)| (v.abs() > threshold_m).then_some(i))
            .collect()
    };
    Ok(PrefitResidualScreen {
        threshold_m,
        maximum_phase_abs_m: phase_m.iter().fold(0.0_f64, |a, b| a.max(b.abs())),
        maximum_code_abs_m: code_m.iter().fold(0.0_f64, |a, b| a.max(b.abs())),
        rejected_phase_rows: rejected(phase_m),
        rejected_code_rows: rejected(code_m),
    })
}
