#![doc = include_str!("../README.md")]
#![cfg_attr(docrs, feature(doc_cfg))]

extern crate gnss_rs as gnss;
extern crate nyx_space as nyx;

// private modules
mod apriori;
mod bias;
mod cfg;
mod clock;
pub(crate) mod interp;
mod observation;
mod orbit;
mod solutions;
mod solver;

// pub export
pub use solver::Error;

// prelude
pub mod prelude {
    pub use crate::apriori::AprioriPosition;
    pub use crate::bias::{BdModel, IonosphericBias, KbModel, NgModel, TroposphericBias};
    pub use crate::cfg::{Config, Filter, Method};
    pub use crate::solutions::{PVTSolution, PVTSolutionType};
    pub use crate::solver::Solver;
    // re-export
    pub use crate::clock::{Clock, ClockIter};
    pub use crate::observation::{Observation, ObservationIter};
    pub use crate::orbit::{Keplerian, Orbit, OrbitIter, Perturbations};
    pub use gnss::prelude::{Constellation, SV};
    pub use hifitime::{Duration, Epoch, TimeScale};
    pub use nalgebra::Vector3;
}
