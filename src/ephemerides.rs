use crate::prelude::{
    SV, Epoch, TimeScale,
};

/// Ephemerides
#[derive(Debug)]
pub struct Ephemerides {
    pub(crate) t: Epoch,
    pub(crate) sv: SV,
    pub(crate) a0: f64,
    pub(crate) a1: f64,
    pub(crate) a2: f64,
    pub(crate) toe: Epoch,
}

impl Ephemerides {
    pub fn new(t: Epoch, sv: SV, a: (f64, f64, f64), toe: Epoch) -> Self {
        // TODO: GPST
        assert!(toe.time_scale == TimeScale::GPST,
            "only GPST is handled at the moment..");
        Self {
            t,
            sv,
            a0: a.0,
            a1: a.1,
            a2: a.2,
            toe,
        }
    }
}

/// Ephemerides data source.
pub trait EphemeridesIter {
    /// Provide Time of Issue of Ephemerides, in chronological order,
    /// with the temporary limitation to the GPST timescale.
    fn next(&mut self) -> Option<Ephemerides>;
}
