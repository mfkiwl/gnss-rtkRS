use log::debug;
use nalgebra::{DMatrix, DVector, Vector3};
use nyx::cosmic::SPEED_OF_LIGHT;
use thiserror::Error;

use crate::{
    cfg::SolverOpts,
    navigation::filter::FilterState,
    navigation::{Input, PVTSolution},
    prelude::Candidate,
};

#[derive(Clone, Debug, Error)]
pub enum Error {
    #[error("gdop {0}: limit exceeded")]
    GDOPOutlier(f64),
    #[error("tdop limit exceeded {0}")]
    TDOPOutlier(f64),
    #[error("innovation outlier |{0}|")]
    InnovationOutlier(f64),
    #[error("coderes limit exceeded {0}")]
    CodeResidual(f64),
}

/// PVT Solution validator
pub struct Validator {}

impl Validator {
    pub fn validate(
        apriori_ecef: Vector3<f64>,
        pool: &[Candidate],
        input: &Input,
        solution: &PVTSolution,
        w: &DMatrix<f64>,
        state: &FilterState,
        opts: &SolverOpts,
    ) -> Result<(), Error> {
        let x = state.estimate();

        let (x, y, z, dt) = (
            apriori_ecef[0] + x[0],
            apriori_ecef[1] + x[1],
            apriori_ecef[2] + x[2],
            x[3] / SPEED_OF_LIGHT,
        );

        let mut residuals = DVector::<f64>::zeros(pool.len());

        for (idx, cd) in pool.iter().enumerate() {
            let sv = input
                .sv
                .iter()
                .filter_map(|(sv, data)| if *sv == cd.sv { Some(data) } else { None })
                .reduce(|k, _| k)
                .unwrap();

            let pr = cd.prefered_pseudorange().unwrap().value;

            let sv_pos = cd.state.unwrap().position;
            let (sv_x, sv_y, sv_z) = (sv_pos[0], sv_pos[1], sv_pos[2]);

            let rho = ((sv_x - x).powi(2) + (sv_y - y).powi(2) + (sv_z - z).powi(2)).sqrt();

            let dt = cd.clock_corr.to_seconds() - dt;

            residuals[idx] = pr;
            residuals[idx] -= rho;
            residuals[idx] += dt * SPEED_OF_LIGHT;
            residuals[idx] -= sv.tropo_bias.value().unwrap_or(0.0);
            residuals[idx] -= sv.iono_bias.value().unwrap_or(0.0);
            residuals[idx] /= w[(idx, idx)];
            debug!(
                "{} ({}): coderes={}/w={}",
                cd.t,
                cd.sv,
                residuals[idx],
                w[(idx, idx)]
            );
        }
        if let Some(max_gdop) = opts.gdop_threshold {
            if solution.dop.gdop > max_gdop {
                return Err(Error::GDOPOutlier(solution.dop.gdop));
            }
        }
        if let Some(max_tdop) = opts.tdop_threshold {
            if solution.dop.tdop > max_tdop {
                return Err(Error::TDOPOutlier(solution.dop.tdop));
            }
        }
        Ok(())
    }
}
