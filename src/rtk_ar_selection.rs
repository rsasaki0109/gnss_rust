//! Deterministic ambiguity subset preselection from libgnss++ rtk_ar_selection.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RelativeVarianceConfig {
    pub multiplier: f64,
    pub minimum_variance_cycles2: f64,
    pub minimum_pairs: usize,
}
impl Default for RelativeVarianceConfig {
    fn default() -> Self {
        Self {
            multiplier: 10.0,
            minimum_variance_cycles2: 1e-4,
            minimum_pairs: 4,
        }
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct RelativeVarianceSelection {
    pub variances_cycles2: Vec<f64>,
    /// The upper middle value for even length, matching the native helper.
    pub median_variance_cycles2: Option<f64>,
    pub threshold_cycles2: Option<f64>,
    /// Empty also denotes fewer than minimum_pairs surviving the threshold.
    /// Native callers retain the full input when this vector is empty.
    pub selected_indices: Vec<usize>,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelectionError {
    InvalidConfiguration,
    InvalidVariance,
    NumericalRange,
}
impl fmt::Display for SelectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidConfiguration => "invalid ambiguity selection configuration",
            Self::InvalidVariance => "invalid ambiguity variance",
            Self::NumericalRange => "ambiguity variance threshold overflows",
        })
    }
}
impl std::error::Error for SelectionError {}

/// Native largest-variance-first nested drop path. Ties drop the first
/// remaining input index. Generates subsets without declaring FIX.
pub fn progressive_variance_drop_subsets(
    variances: &[f64],
    minimum_pairs: usize,
    max_drop_steps: usize,
) -> Result<Vec<Vec<usize>>, SelectionError> {
    if !(1..=128).contains(&minimum_pairs) || max_drop_steps > 128 {
        return Err(SelectionError::InvalidConfiguration);
    }
    if variances.len() > 128 || variances.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(SelectionError::InvalidVariance);
    }
    let mut current: Vec<_> = (0..variances.len()).collect();
    let mut subsets = Vec::new();
    for _ in 0..variances
        .len()
        .saturating_sub(minimum_pairs)
        .min(max_drop_steps)
    {
        let mut worst = current[0];
        for &i in &current {
            if variances[i] > variances[worst] {
                worst = i;
            }
        }
        current.retain(|&i| i != worst);
        subsets.push(current.clone());
    }
    Ok(subsets)
}

/// Native filterPairsByRelativeVariance for finite nonnegative variances.
/// No diagonal floor or covariance update is performed. Invalid inputs fail
/// explicitly instead of relying on native NaN sorting/comparison behavior.
pub fn filter_pairs_by_relative_variance(
    variances: &[f64],
    config: RelativeVarianceConfig,
) -> Result<RelativeVarianceSelection, SelectionError> {
    if !config.multiplier.is_finite()
        || config.multiplier <= 0.0
        || !config.minimum_variance_cycles2.is_finite()
        || config.minimum_variance_cycles2 < 0.0
        || !(1..=128).contains(&config.minimum_pairs)
    {
        return Err(SelectionError::InvalidConfiguration);
    }
    if variances.len() > 128 || variances.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(SelectionError::InvalidVariance);
    }
    let mut result = RelativeVarianceSelection {
        variances_cycles2: variances.to_vec(),
        median_variance_cycles2: None,
        threshold_cycles2: None,
        selected_indices: Vec::new(),
    };
    if variances.is_empty() {
        return Ok(result);
    }
    let mut sorted = variances.to_vec();
    sorted.sort_by(f64::total_cmp);
    let median = sorted[sorted.len() / 2];
    let product = median * config.multiplier;
    if !product.is_finite() {
        return Err(SelectionError::NumericalRange);
    }
    let threshold = product.max(config.minimum_variance_cycles2);
    let indices: Vec<_> = variances
        .iter()
        .enumerate()
        .filter_map(|(i, &v)| (v <= threshold).then_some(i))
        .collect();
    result.median_variance_cycles2 = Some(median);
    result.threshold_cycles2 = Some(threshold);
    if indices.len() >= config.minimum_pairs {
        result.selected_indices = indices;
    }
    Ok(result)
}
