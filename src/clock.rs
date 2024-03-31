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
    /// Clock correction (to be determined)
    pub(crate) correction: Option<Duration>,
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
            correction: None,
            drift: drift.unwrap_or(0.0),
            drift_rate: drift_rate.unwrap_or(0.0),
        }
    }
    pub(crate) fn clock_correction(&mut self, t: Epoch, toe: Epoch) {
        match self.sv.constellation {
            Constellation::Glonass => {
                panic!("glo_sv_clock_corr");
            },
            _ => {
                let mut dt = (t - toe).to_seconds();
                // TODO: does this apply to GST/BDT correctly as well ?
                const WEEK_SECONDS: f64 = 604800.0;
                if dt > WEEK_SECONDS / 2.0 {
                    dt -= WEEK_SECONDS;
                } else if dt < -WEEK_SECONDS / 2.0 {
                    dt += WEEK_SECONDS;
                }
                self.correction = Some(Duration::from_seconds(
                    self.offset + self.drift * dt + self.drift_rate * dt.powi(2),
                ));
            },
        }
    }
}
