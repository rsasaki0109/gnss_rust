//! Read-only changes of double-difference measurement basis.
//!
//! A reference change must transform the design, residual and full covariance
//! together. This preserves Gaussian measurement information, including
//! correlations between phase and code blocks. Unweighted RMS and maximum
//! residual gates generally do not have this invariance.
use crate::SatelliteId;
use crate::filter::{self, FilterError, Matrix};
use std::collections::BTreeSet;
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReferenceError {
    InvalidBasis,
    ReferenceUnavailable,
    Filter(FilterError),
}
impl fmt::Display for ReferenceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBasis => {
                f.write_str("DD basis requires distinct satellites of one system")
            }
            Self::ReferenceUnavailable => {
                f.write_str("reference satellite is absent from the basis")
            }
            Self::Filter(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for ReferenceError {}
impl From<FilterError> for ReferenceError {
    fn from(value: FilterError) -> Self {
        Self::Filter(value)
    }
}

/// Row order is the supplied satellite order excluding the reference.
/// Every row means reference minus target. Each measurement block uses this
/// same ordering; satellites, raw observations and time tags are not changed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DoubleDifferenceBasis {
    satellites: Vec<SatelliteId>,
    reference: SatelliteId,
}
impl DoubleDifferenceBasis {
    pub fn new(
        satellites: Vec<SatelliteId>,
        reference: SatelliteId,
    ) -> Result<Self, ReferenceError> {
        let seen: BTreeSet<_> = satellites.iter().copied().collect();
        if !(2..=128).contains(&satellites.len())
            || seen.len() != satellites.len()
            || satellites.iter().any(|s| s.system() != reference.system())
        {
            return Err(ReferenceError::InvalidBasis);
        }
        if !seen.contains(&reference) {
            return Err(ReferenceError::ReferenceUnavailable);
        }
        Ok(Self {
            satellites,
            reference,
        })
    }

    pub fn satellites(&self) -> &[SatelliteId] {
        &self.satellites
    }
    pub fn reference(&self) -> SatelliteId {
        self.reference
    }
    pub fn targets(&self) -> impl Iterator<Item = SatelliteId> + '_ {
        self.satellites
            .iter()
            .copied()
            .filter(|s| *s != self.reference)
    }

    /// Block-diagonal M with new_DD = M old_DD. It acts on complete DD
    /// blocks, not a selected subset of rows. At most 256 output rows.
    pub fn transformation(
        &self,
        next_reference: SatelliteId,
        blocks: usize,
    ) -> Result<Matrix, ReferenceError> {
        if !self.satellites.contains(&next_reference) {
            return Err(ReferenceError::ReferenceUnavailable);
        }
        let k = self.satellites.len() - 1;
        let rows = k
            .checked_mul(blocks)
            .ok_or(FilterError::InvalidDimensions)?;
        if !(1..=256).contains(&rows) {
            return Err(FilterError::InvalidDimensions.into());
        }
        let old_targets: Vec<_> = self.targets().collect();
        let new_targets: Vec<_> = self
            .satellites
            .iter()
            .copied()
            .filter(|s| *s != next_reference)
            .collect();
        let mut result = vec![vec![0.0; rows]; rows];
        for block in 0..blocks {
            for (i, satellite) in new_targets.iter().enumerate() {
                if let Some(j) = old_targets.iter().position(|s| s == satellite) {
                    result[block * k + i][block * k + j] = 1.0;
                }
                if let Some(j) = old_targets.iter().position(|s| *s == next_reference) {
                    result[block * k + i][block * k + j] -= 1.0;
                }
            }
        }
        Ok(result)
    }

    /// Transform H, v and R together: H'=MH, v'=Mv, R'=MRMᵀ.
    /// R must be finite symmetric positive definite. All cross-block terms
    /// are retained. Inputs are borrowed and remain intact on any error.
    /// State coordinates remain fixed. A caller storing DD states must also
    /// transform its state, prior and process covariance separately.
    pub fn rebase(
        &self,
        next_reference: SatelliteId,
        blocks: usize,
        design: &Matrix,
        residual: &[f64],
        covariance: &Matrix,
    ) -> Result<RebasedMeasurements, ReferenceError> {
        let m = self.transformation(next_reference, blocks)?;
        let rows = m.len();
        let n = design.first().map_or(0, Vec::len);
        if !(1..=128).contains(&n)
            || design.len() != rows
            || design.iter().any(|row| row.len() != n)
            || residual.len() != rows
            || covariance.len() != rows
            || covariance.iter().any(|row| row.len() != rows)
        {
            return Err(FilterError::InvalidDimensions.into());
        }
        if !design
            .iter()
            .flatten()
            .chain(residual)
            .all(|v| v.is_finite())
        {
            return Err(FilterError::NonFinite.into());
        }
        filter::cholesky(covariance)?;
        let design = filter::multiply(&m, design);
        let residual: Vec<f64> = m
            .iter()
            .map(|row| row.iter().zip(residual).map(|(a, b)| a * b).sum())
            .collect();
        let mut covariance =
            filter::multiply(&filter::multiply(&m, covariance), &filter::transpose(&m));
        // Symmetrize only floating-point roundoff of the exact congruence.
        for i in 1..rows {
            let (previous, remaining) = covariance.split_at_mut(i);
            let row = &mut remaining[0];
            for (j, other) in previous.iter_mut().enumerate() {
                let value = (row[j] + other[i]) * 0.5;
                row[j] = value;
                other[i] = value;
            }
        }
        if !design
            .iter()
            .flatten()
            .chain(&residual)
            .chain(covariance.iter().flatten())
            .all(|v| v.is_finite())
        {
            return Err(FilterError::NumericalFailure.into());
        }
        filter::cholesky(&covariance)?;
        Ok(RebasedMeasurements {
            basis: Self::new(self.satellites.clone(), next_reference)?,
            design,
            residual,
            covariance,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RebasedMeasurements {
    pub basis: DoubleDifferenceBasis,
    pub design: Matrix,
    pub residual: Vec<f64>,
    pub covariance: Matrix,
}
impl RebasedMeasurements {
    /// Nominal vᵀR⁻¹v. For postfit residuals this is not innovation NIS,
    /// a calibrated chi-square statistic or an accuracy certificate.
    pub fn nominal_weighted_sum_of_squares(&self) -> Result<f64, FilterError> {
        if !(1..=256).contains(&self.residual.len()) || self.residual.len() != self.covariance.len()
        {
            return Err(FilterError::InvalidDimensions);
        }
        if !self.residual.iter().all(|v| v.is_finite()) {
            return Err(FilterError::NonFinite);
        }
        let l = filter::cholesky(&self.covariance)?;
        let weighted = filter::solve_cholesky(&l, &self.residual);
        let value = self
            .residual
            .iter()
            .zip(weighted)
            .map(|(a, b)| a * b)
            .sum::<f64>();
        if !value.is_finite() || value < 0.0 {
            return Err(FilterError::NumericalFailure);
        }
        Ok(value)
    }
}
