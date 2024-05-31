pub mod solutions;
pub use solutions::{PVTSolution, PVTSolutionType};

mod filter;

pub use filter::{Filter, FilterState};

use log::{debug, error};
use std::collections::HashMap;

use crate::{
    ambiguity::Ambiguities,
    bias::{Bias, IonosphereBias, RuntimeParam as BiasRuntimeParams, TropoModel, TroposphereBias},
    candidate::Candidate,
    cfg::Config,
    prelude::{Error, Method, SV},
};

use nalgebra::{base::dimension::U4, DMatrix, DVector, OMatrix};

use nyx::cosmic::SPEED_OF_LIGHT;

/// SV Navigation information
#[derive(Debug, Clone, Default)]
pub struct SVInput {
    /// SV azimuth angle in degrees
    pub azimuth: f64,
    /// SV elevation angle in degrees
    pub elevation: f64,
    /// Ionospheric bias in meters of delay
    pub iono_bias: Bias,
    /// Tropospheric bias in meters of delay
    pub tropo_bias: Bias,
}

/// Navigation Input
#[derive(Debug, Clone)]
pub struct Input {
    /// Measurement vector
    pub y: DVector<f64>,
    /// NAV Matrix
    pub g: DMatrix<f64>,
    /// SV dependent data
    pub sv: HashMap<SV, SVInput>,
}

/// Navigation Output
#[derive(Debug, Clone)]
pub struct Output {
    /// Time Dilution of Precision
    pub tdop: f64,
    /// Geometric Dilution of Precision
    pub gdop: f64,
    /// Position Dilution of Precision
    pub pdop: f64,
    /// Q covariance matrix
    pub q: DMatrix<f64>,
}

impl Output {
    pub(crate) fn q_covar4x4(&self) -> OMatrix<f64, U4, U4> {
        OMatrix::<f64, U4, U4>::new(
            self.q[(0, 0)],
            self.q[(0, 1)],
            self.q[(0, 2)],
            self.q[(0, 3)],
            self.q[(1, 0)],
            self.q[(1, 1)],
            self.q[(1, 2)],
            self.q[(1, 3)],
            self.q[(2, 0)],
            self.q[(2, 1)],
            self.q[(2, 2)],
            self.q[(2, 3)],
            self.q[(3, 0)],
            self.q[(3, 1)],
            self.q[(3, 2)],
            self.q[(3, 3)],
        )
    }
}

impl Input {
    /// Forms new Navigation Input
    pub fn new(
        apriori: (f64, f64, f64),
        apriori_geo: (f64, f64, f64),
        cfg: &Config,
        cd: &[Candidate],
        ambiguities: &Ambiguities,
        iono_bias: &IonosphereBias,
        tropo_bias: &TroposphereBias,
    ) -> Result<Self, Error> {
        let mut y = DVector::<f64>::zeros(cd.len() * 2);
        let mut g = DMatrix::<f64>::zeros(cd.len() * 2, 4 + cd.len());
        let mut sv = HashMap::<SV, SVInput>::with_capacity(cd.len());
        /*
         * Compensate for ARP (if possible)
         */
        let apriori = match cfg.arp_enu {
            Some(offset) => (
                apriori.0 + offset.0,
                apriori.1 + offset.1,
                apriori.2 + offset.2,
            ),
            None => apriori,
        };

        let (x0, y0, z0) = apriori;

        for row in 0..cd.len() {
            let mut sv_input = SVInput::default();

            let state = cd[row].state.ok_or(Error::UnresolvedState)?;
            let clock_corr = cd[row].clock_corr.to_seconds();

            let (azimuth, elevation) = (state.azimuth, state.elevation);
            sv_input.azimuth = azimuth;
            sv_input.elevation = elevation;

            let (sv_x, sv_y, sv_z) = (state.position[0], state.position[1], state.position[2]);
            let rho = ((sv_x - x0).powi(2) + (sv_y - y0).powi(2) + (sv_z - z0).powi(2)).sqrt();
            let (x_i, y_i, z_i) = ((x0 - sv_x) / rho, (y0 - sv_y) / rho, (z0 - sv_z) / rho);

            g[(2 * row, 0)] = x_i;
            g[(2 * row, 1)] = y_i;
            g[(2 * row, 2)] = z_i;
            g[(2 * row, 3)] = 1.0_f64;

            g[(2 * row + 1, 0)] = x_i;
            g[(2 * row + 1, 1)] = y_i;
            g[(2 * row + 1, 2)] = z_i;
            g[(2 * row + 1, 3)] = 1.0_f64;
            g[(2 * row + 1, 4 + row)] = 1.0_f64;

            let mut models = 0.0_f64;

            if cfg.modeling.sv_clock_bias {
                models -= clock_corr * SPEED_OF_LIGHT;
            }
            if let Some(delay) = cfg.externalref_delay {
                models -= delay * SPEED_OF_LIGHT;
            }

            let (pr, frequency) = match cfg.method {
                Method::SPP => {
                    let pr = cd[row]
                        .prefered_pseudorange()
                        .ok_or(Error::MissingPseudoRange)?;
                    (pr.value, pr.carrier.frequency())
                },
                Method::CPP | Method::PPP => {
                    let pr = cd[row]
                        .code_if_combination()
                        .ok_or(Error::PseudoRangeCombination)?;
                    (pr.value, pr.reference.frequency())
                },
            };

            // frequency dependent delay
            for delay in &cfg.int_delay {
                if delay.frequency == frequency {
                    models += delay.delay * SPEED_OF_LIGHT;
                }
            }

            /*
             * IONO + TROPO biases
             */
            let rtm = BiasRuntimeParams {
                t: cd[row].t,
                elevation,
                azimuth,
                frequency,
                apriori_geo,
            };

            /*
             * TROPO
             */
            if cfg.modeling.tropo_delay {
                if tropo_bias.needs_modeling() {
                    let bias = TroposphereBias::model(TropoModel::Niel, &rtm);
                    debug!(
                        "{}({}): modeled tropo delay {:.3E}[m]",
                        cd[row].t, cd[row].sv, bias
                    );
                    models += bias;
                    sv_input.tropo_bias = Bias::modeled(bias);
                } else if let Some(bias) = tropo_bias.bias(&rtm) {
                    debug!(
                        "{}({}): measured tropo delay {:.3E}[m]",
                        cd[row].t, cd[row].sv, bias
                    );
                    models += bias;
                    sv_input.tropo_bias = Bias::measured(bias);
                }
            }

            /*
             * IONO
             */
            if cfg.method == Method::SPP && cfg.modeling.iono_delay {
                if let Some(bias) = iono_bias.bias(&rtm) {
                    debug!(
                        "{} : modeled iono delay (f={:.3E}Hz) {:.3E}[m]",
                        cd[row].t, rtm.frequency, bias
                    );
                    models += bias;
                    sv_input.iono_bias = Bias::modeled(bias);
                }
            }

            y[2 * row] = pr - rho - models;

            if cfg.method == Method::PPP {
                let cmb = cd[row]
                    .phase_if_combination()
                    .ok_or(Error::PseudoRangeCombination)?;

                let f_1 = cmb.reference.frequency();
                let lambda_j = cmb.lhs.wavelength();
                let f_j = cmb.lhs.frequency();

                let (lambda_n, lambda_w) =
                    (SPEED_OF_LIGHT / (f_1 + f_j), SPEED_OF_LIGHT / (f_1 - f_j));

                let bias = if let Some(ambiguity) = ambiguities.get(&(cd[row].sv, cmb.reference)) {
                    let (n_1, n_w) = (ambiguity.n_1, ambiguity.n_w);
                    let b_c = lambda_n * (n_1 + (lambda_w / lambda_j) * n_w);
                    debug!(
                        "{} ({}/{}) b_c: {}",
                        cd[row].t, cd[row].sv, cmb.reference, b_c
                    );
                    b_c
                } else {
                    error!(
                        "{} ({}/{}): unresolved ambiguity",
                        cd[row].t, cd[row].sv, cmb.reference
                    );
                    return Err(Error::UnresolvedAmbiguity);
                };

                // TODO: conclude windup
                let windup = 0.0_f64;

                y[2 * row + 1] = cmb.value - rho - models - windup + bias;
            } else {
                y[2 * row + 1] = y[2 * row];
            }
            sv.insert(cd[row].sv, sv_input);
        }

        assert!(g.nrows() == y.nrows(), "invalid g matrix formulation");
        assert!(g.ncols() == cd.len() + 4, "invalid g matrix formulation");
        debug!(
            "y({}): {} g({},{}): {}",
            y.nrows(),
            y,
            g.nrows(),
            g.ncols(),
            g
        );
        Ok(Self { y, g, sv })
    }
}
