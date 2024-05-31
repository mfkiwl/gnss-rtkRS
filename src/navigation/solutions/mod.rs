//! PVT Solutions
use crate::prelude::{Candidate, Duration, Error, TimeScale, Vector3, SV};
use std::collections::HashMap;

use super::SVInput;
use nalgebra::base::{DVector, Matrix3, Matrix4};

pub(crate) mod validator;

#[cfg(feature = "serde")]
use serde::Deserialize;

#[derive(Debug, Copy, Clone, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(Deserialize))]
pub enum PVTSolutionType {
    /// Default, complete solution with Position,
    /// Velocity and Time components. Requires either
    /// 4 vehicles in sight, or 3 if you're working in fixed altitude
    /// (provided ahead of time).
    #[default]
    PositionVelocityTime,
    /// Resolve Time component only. Only requires 1 vehicle in sight.
    TimeOnly,
}

impl std::fmt::Display for PVTSolutionType {
    /*
     * Prints self
     */
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            Self::PositionVelocityTime => write!(f, "PVT"),
            Self::TimeOnly => write!(f, "TimeOnly"),
        }
    }
}

/// Dilution of Precision
#[derive(Debug, Clone)]
pub struct DOP {
    /// Geometric Dilution of Precision
    pub gdop: f64,
    /// Position Dilution of Precision
    pub pdop: f64,
    /// Time Dilution of Precision
    pub tdop: f64,
    /// 4x4 matrix
    q: Matrix4<f64>,
}

impl DOP {
    /// Builds new DOP
    pub(crate) fn new(pool: &[Candidate], solution: &DVector<f64>) -> Result<Self, Error> {
        let mut q = Matrix4::<f64>::zeros();
        let (x, y, z) = (solution[0], solution[1], solution[2]);
        for i in 0..4 {
            let position = pool[i].state.unwrap().position;
            let (x_i, y_i, z_i) = (position[0], position[1], position[2]);
            let r = ((x_i - x).powi(2) + (y_i - y).powi(2) + (z_i - z).powi(2)).sqrt();
            q[(i, 0)] = (x_i - x) / r;
            q[(i, 1)] = (y_i - y) / r;
            q[(i, 2)] = (z_i - z) / r;
            q[(i, 3)] = 1.0;
        }
        let q = (q.transpose() * q)
            .try_inverse()
            .ok_or(Error::MatrixInversionError)?;

        let pdop = (q[(0, 0)] + q[(1, 1)] + q[(2, 2)]).sqrt();
        let tdop = q[(3, 3)].sqrt();
        let gdop = (pdop.powi(2) + tdop.powi(2)).sqrt();
        Ok(Self {
            q,
            tdop,
            gdop,
            pdop,
        })
    }
    fn q_3x3(&self) -> Matrix3<f64> {
        Matrix3::<f64>::new(
            self.q[(0, 0)],
            self.q[(0, 1)],
            self.q[(0, 2)],
            self.q[(1, 0)],
            self.q[(1, 1)],
            self.q[(1, 2)],
            self.q[(2, 0)],
            self.q[(2, 1)],
            self.q[(2, 2)],
        )
    }
    fn q_enu(&self, lat: f64, lon: f64) -> Matrix3<f64> {
        let r = Matrix3::<f64>::new(
            -lon.sin(),
            -lon.cos() * lat.sin(),
            lat.cos() * lon.cos(),
            lon.cos(),
            -lat.sin() * lon.sin(),
            lat.cos() * lon.sin(),
            0.0_f64,
            lat.cos(),
            lon.sin(),
        );
        r.clone().transpose() * self.q_3x3() * r
    }
    pub fn hdop(&self, lat: f64, lon: f64) -> f64 {
        let q = self.q_enu(lat, lon);
        (q[(0, 0)] + q[(1, 1)]).sqrt()
    }
    pub fn vdop(&self, lat: f64, lon: f64) -> f64 {
        self.q_enu(lat, lon)[(2, 2)].sqrt()
    }
}

/// PVT Solution, always expressed as the correction to apply
/// to an Apriori / static position.
#[derive(Debug, Clone)]
// #[cfg_attr(feature = "serde", derive(Serialize))]
pub struct PVTSolution {
    /// Position error (in [m] ECEF)
    pub position: Vector3<f64>,
    /// Absolute Velocity (in [m/s] ECEF).
    pub velocity: Vector3<f64>,
    /// Timescale
    pub timescale: TimeScale,
    /// Offset to timescale
    pub dt: Duration,
    /// Space Vehicles that helped form this solution
    /// and data associated to each individual SV
    pub sv: HashMap<SV, SVInput>,
    /// Dilution of Precision
    pub dop: DOP,
}

impl PVTSolution {
    /// Returns list of Space Vehicles (SV) that help form this solution.
    pub fn sv(&self) -> Vec<SV> {
        self.sv.keys().copied().collect()
    }
}
