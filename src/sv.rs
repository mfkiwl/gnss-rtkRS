use crate::prelude::{Duration, SV};

/// In-depth static SV information used to obtain precise results
pub struct SVInfo {
    /// SV
    pub sv: SV,
    /// Onboard total (group) delay
    pub tgd: Duration,
}

/// Implement this trait to provide more info on SV
pub trait SVInfoIter {
    /// Provide in depth information on SV for precise results
    fn next(&mut self) -> Option<SVInfo>;
}
