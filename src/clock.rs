use gnss::prelude::{Constellation, SV};
use hifitime::{Duration, Epoch};

/// Implement this trait to provide SV Clock States
pub trait ClockIter {
    /// Provide Clock states in chronological order.
    /// In order to simplify internal logic, we will panic if Clock states
    /// are not provided in chronological order.
    /// Unevenly spaced states (in time) will delay PVT solution production,
    /// as we maintain steady interpolation errors.
    fn next(&mut self) -> Option<Clock>;
}

/// Satellite Vehicle Clock state
#[derive(Debug, Clone)]
pub struct Clock {
    /// Satellite Vehicle
    pub(crate) sv: SV,
    /// Instant of this State snapshot
    pub(crate) epoch: Epoch,
    /// Clock offset [s]
    pub(crate) offset: f64,
    /// Clock drift [s/s]
    pub(crate) drift: f64,
    /// Clock drift rate [s/s^2]
    pub(crate) drift_rate: f64,
}

impl Clock {
    /// Builds new Clock State from (Clock Offset [s],
    /// Clock drift [s/s] and Clock drift rate [s/s^2]) estimates.
    pub fn new(
        sv: SV,
        epoch: Epoch,
        offset: f64,
        drift: Option<f64>,
        drift_rate: Option<f64>,
    ) -> Self {
        Self {
            sv,
            epoch,
            offset,
            drift: drift.unwrap_or(0.0),
            drift_rate: drift_rate.unwrap_or(0.0),
        }
    }
}
