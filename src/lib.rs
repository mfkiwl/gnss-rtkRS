#![doc = include_str!("../README.md")]
#![cfg_attr(docrs, feature(doc_cfg))]

extern crate gnss_rs as gnss;
extern crate nyx_space as nyx;

// private modules
mod apriori;
mod bias;
mod cfg;
mod clock;
mod ephemerides;
mod interp;
mod observation;
mod orbit;
mod solver;

pub(crate) mod navigation;

// prelude
pub mod prelude {
    pub use crate::apriori::AprioriPosition;
    pub use crate::bias::{
        BdModel, IonosphereBiasModel, IonosphereBiasModelIter, KbModel, NgModel,
    };
    pub use crate::cfg::{Config, Filter, Method};
    pub use crate::clock::{Clock, ClockIter};
    pub use crate::ephemerides::EphemeridesIter;
    pub use crate::navigation::solutions::{PVTSolution, PVTSolutionType};
    pub use crate::observation::{Observation, ObservationIter};
    pub use crate::orbit::{Orbit, OrbitIter};
    pub use crate::solver::Solver;
    // re-export
    pub use gnss::prelude::{Constellation, SV};
    pub use hifitime::{Duration, Epoch, TimeScale};
    pub use nalgebra::Vector3;
}
