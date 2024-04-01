pub(crate) mod tropo;

pub(crate) mod iono;
pub use iono::{BdModel, IonosphereBiasModel, IonosphereBiasModelIter, KbModel, NgModel};

use crate::prelude::Epoch;

pub(crate) trait BiasModel {
    /// Confirms Bias model validity at desired Epoch
    fn valid(&self, t: Epoch) -> bool;
}
