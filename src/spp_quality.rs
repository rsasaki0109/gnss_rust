//! SPP noise and QC utilities ported from pinned libgnss++ spp_utils.
//! Copyright (c) 2024 LibGNSS++ Contributors. See LICENSE.
use crate::Error;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SppVarianceModel {
    #[default]
    ElevationOnly,
    NativeSnrAtmosphere,
}
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeasurementVarianceInputs {
    pub elevation_rad: f64,
    pub snr_db_hz: Option<f64>,
    pub pseudorange_sigma_m: f64,
    pub ionosphere_delay_m: f64,
    pub troposphere_delay_m: f64,
    pub ionosphere_corrected: bool,
    pub troposphere_corrected: bool,
    pub snr_reference_db_hz: f64,
    pub ionosphere_sigma_scale: f64,
    pub troposphere_sigma_scale: f64,
}
impl Default for MeasurementVarianceInputs {
    fn default() -> Self {
        Self {
            elevation_rad: 0.0,
            snr_db_hz: None,
            pseudorange_sigma_m: 3.0,
            ionosphere_delay_m: 0.0,
            troposphere_delay_m: 0.0,
            ionosphere_corrected: false,
            troposphere_corrected: false,
            snr_reference_db_hz: 45.0,
            ionosphere_sigma_scale: 0.5,
            troposphere_sigma_scale: 0.1,
        }
    }
}

/// Native code/elevation/CN0/atmosphere variance in square metres. Missing,
/// nonpositive or nonfinite CN0 causes no CN0 penalty, as in native code;
/// the caller retains its raw distinction. Nonfinite active corrections or
/// overflow fail explicitly instead of falling back to an optimistic variance.
pub fn measurement_variance(v: MeasurementVarianceInputs) -> Result<f64, Error> {
    if !v.pseudorange_sigma_m.is_finite() || v.pseudorange_sigma_m < 0.0 {
        return Err(Error::InvalidConfiguration);
    }
    if v.ionosphere_corrected && !v.ionosphere_delay_m.is_finite()
        || v.troposphere_corrected && !v.troposphere_delay_m.is_finite()
    {
        return Err(Error::NonFinite);
    }
    let base = v.pseudorange_sigma_m.max(0.1);
    let sin = v.elevation_rad.sin();
    let sin = if !sin.is_finite() || sin < 0.1 {
        0.1
    } else {
        sin
    };
    let reference = if v.snr_reference_db_hz.is_finite() && v.snr_reference_db_hz > 0.0 {
        v.snr_reference_db_hz
    } else {
        45.0
    };
    let mut code = base / sin;
    if let Some(snr) = v
        .snr_db_hz
        .filter(|s| s.is_finite() && *s > 0.0 && *s < reference)
    {
        code *= 10.0_f64.powf((reference - snr) / 20.0).clamp(1.0, 4.0);
    }
    let mut variance = code * code;
    let iono = if v.ionosphere_corrected {
        let scale = if v.ionosphere_sigma_scale.is_finite() && v.ionosphere_sigma_scale >= 0.0 {
            v.ionosphere_sigma_scale
        } else {
            0.5
        };
        (v.ionosphere_delay_m.abs() * scale).clamp(0.25, 8.0)
    } else {
        2.5 / sin
    };
    variance += iono * iono;
    let trop = if v.troposphere_corrected {
        let scale = if v.troposphere_sigma_scale.is_finite() && v.troposphere_sigma_scale >= 0.0 {
            v.troposphere_sigma_scale
        } else {
            0.1
        };
        (v.troposphere_delay_m.abs() * scale).clamp(0.05, 3.0)
    } else {
        0.5 / sin
    };
    variance += trop * trop;
    if !variance.is_finite() || variance <= 0.0 {
        return Err(Error::NonFinite);
    }
    Ok(variance)
}

fn median(values: &mut [f64]) -> f64 {
    values.sort_by(f64::total_cmp);
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        0.5 * (values[n / 2 - 1] + values[n / 2])
    }
}

/// Native robust mask. An all-true return for insufficient redundancy or fewer
/// than four surviving inliers means "do not filter", not a valid SPP solution.
pub fn residual_inlier_mask(
    residuals: &[f64],
    threshold_sigma: f64,
    sigma_floor: f64,
) -> Vec<bool> {
    let unchanged = vec![true; residuals.len()];
    if residuals.len() <= 4 || !threshold_sigma.is_finite() || threshold_sigma <= 0.0 {
        return unchanged;
    }
    let floor = if sigma_floor.is_finite() && sigma_floor > 0.0 {
        sigma_floor
    } else {
        1.0
    };
    let mut finite: Vec<_> = residuals
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .collect();
    if finite.len() <= 4 {
        return unchanged;
    }
    let center = median(&mut finite);
    let mut deviations: Vec<_> = finite.iter().map(|v| (v - center).abs()).collect();
    let mut sigma = (1.4826 * median(&mut deviations)).max(floor);
    if !sigma.is_finite() || sigma <= 0.0 {
        sigma =
            (finite.iter().map(|v| (v - center).powi(2)).sum::<f64>() / finite.len() as f64).sqrt();
    }
    if !sigma.is_finite() || sigma <= 0.0 {
        return unchanged;
    }
    let limit = threshold_sigma * sigma;
    let mask: Vec<_> = residuals
        .iter()
        .map(|v| v.is_finite() && (v - center).abs() <= limit)
        .collect();
    if mask.iter().filter(|&&v| v).count() >= 4 {
        mask
    } else {
        unchanged
    }
}

/// Native first-best leave-one-out RMS choice with both fractional and
/// absolute improvement requirements. Invalid candidate slots are skipped.
/// This helper neither performs a refit nor declares measurement validity.
pub fn select_raim_fde_candidate(
    baseline_rms: f64,
    candidates: &[f64],
    minimum_ratio: f64,
    minimum_m: f64,
) -> Option<usize> {
    if !baseline_rms.is_finite() || baseline_rms <= 0.0 {
        return None;
    }
    let mut best = baseline_rms;
    let mut index = None;
    for (i, &v) in candidates.iter().enumerate() {
        if v.is_finite() && v >= 0.0 && v < best {
            best = v;
            index = Some(i);
        }
    }
    let ratio = if minimum_ratio.is_finite() && minimum_ratio >= 0.0 {
        minimum_ratio
    } else {
        0.0
    };
    let absolute = if minimum_m.is_finite() && minimum_m >= 0.0 {
        minimum_m
    } else {
        0.0
    };
    index.filter(|_| baseline_rms - best >= absolute.max(baseline_rms * ratio))
}
