use gnss::prelude::SV;
use hifitime::Epoch;

/// Implement this trait to provide SV signal observations.
pub trait ObservationIter {
    /// Provide SV signal observation in chronological order.
    /// On this returns None,
    fn next(&mut self) -> Option<Observation>;
}

/// Signal observation
#[derive(Debug, Default, Clone)]
pub struct Observation {
    /// SV (signal emitter)
    pub(crate) sv: SV,
    /// Actual observation
    pub(crate) value: f64,
    /// Signal sampling Epoch
    pub(crate) epoch: Epoch,
    /// carrier frequency [Hz]
    pub(crate) frequency_hz: f64,
    /// Optional (but recommended) SNR in [dB]
    pub(crate) snr_db: Option<f64>,
}

impl Observation {
    /// Builds new [Observation] of one value observed from one carrier signal
    /// of given frequency. We recommend providing the SNR estimate (in dB) if that is feasible.
    pub fn new(sv: SV, epoch: Epoch, value: f64, frequency_hz: f64, snr_db: Option<f64>) -> Self {
        Self {
            sv,
            epoch,
            value,
            frequency_hz,
            snr_db,
        }
    }
}
