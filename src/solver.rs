//! PVT solver

use hifitime::{Epoch, Unit};
use log::{debug, error, warn};
use map_3d::deg2rad;
use thiserror::Error;

use nyx::cosmic::eclipse::{eclipse_state, EclipseState};
use nyx::cosmic::SPEED_OF_LIGHT;
use nyx::md::prelude::{Arc, Cosm};
use nyx::md::prelude::{Bodies, Frame, LightTimeCalc};

use gnss::prelude::SV;

use nalgebra::{DVector, Matrix3, Matrix4, Matrix4x1, MatrixXx4, Vector3};

use crate::{
    apriori::AprioriPosition,
    bias::{IonosphericBias, TroposphericBias},
    cfg::{Config, Filter, Method},
    clock::{Clock, ClockIter},
    observation::{Observation, ObservationIter},
    orbit::{Orbit, OrbitIter},
    solutions::{PVTSVData, PVTSolution, PVTSolutionType},
};

use itertools::Itertools;

/// Solver Error
#[derive(Debug, Clone, Error)]
pub enum Error {
    #[error("observations do not match solving strategy/method")]
    UnfitObservations,
    #[error("need more candidates to resolve a {0} a solution")]
    NotEnoughInputCandidates(PVTSolutionType),
    #[error("not enough candidates fit criteria")]
    NotEnoughFittingCandidates,
    #[error("failed to invert navigation matrix")]
    MatrixInversionError,
    #[error("failed to invert covar matrix")]
    CovarMatrixInversionError,
    #[error("reolved NaN: invalid input matrix")]
    TimeIsNan,
    #[error("undefined apriori position")]
    UndefinedAprioriPosition,
    #[error("missing pseudo range observation")]
    MissingPseudoRange,
    #[error("at least one pseudo range observation is mandatory")]
    NeedsAtLeastOnePseudoRange,
    #[error("failed to model or measure ionospheric delay")]
    MissingIonosphericDelayValue,
    #[error("unresolved state: interpolation should have passed")]
    UnresolvedState,
    #[error("unable to form signal combination")]
    SignalRecombination,
    #[error("physical non sense: rx prior tx")]
    PhysicalNonSenseRxPriorTx,
    #[error("physical non sense: t_rx is too late")]
    PhysicalNonSenseRxTooLate,
    // #[error("invalidated solution: {0}")]
    // InvalidatedSolution(SolutionInvalidation),
    // // Kalman filter bad op: should never happen
    // #[error("uninitialized kalman filter!")]
    // UninitializedKalmanFilter,
}

#[derive(Debug, Clone)]
pub(crate) struct LSQState {
    /* p matrix */
    pub(crate) p: Matrix4<f64>,
    /* x estimate */
    pub(crate) x: Matrix4x1<f64>,
}

#[derive(Debug, Clone)]
pub(crate) struct KfState {
    /* p matrix */
    pub(crate) p: Matrix4<f64>,
    /* x estimate */
    pub(crate) x: Matrix4x1<f64>,
}

// Filter state
#[derive(Debug, Clone)]
pub(crate) enum FilterState {
    /// LSQ state
    LSQState(LSQState),
    /// KF state
    KfState(KfState),
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
enum State {
    /// Need carrier signal measurements
    #[default]
    SignalAcquisition,
    /// Signal quality filter
    SignalFilter,
    /// SV clock interpolation
    ClockInterpolation,
    /// SV orbital state interpolation
    OrbitInterpolation,
    /// Processing
    Processing,
}

/// PVT Solver
#[derive(Debug, Clone)]
pub struct Solver {
    /// Solver configuration
    cfg: Config,
    /// Solver state
    state: State,
    /// apriori position
    apriori: AprioriPosition,
    /// Type of solutions to resolve
    solutions_type: PVTSolutionType,
    /// Clock states buffer
    clocks: Vec<Clock>,
    /// Orbital states buffer
    orbits: Vec<Orbit>,
    /// Observations buffer
    observations: Vec<Observation>,
    /// Cosmic model
    cosmic: Arc<Cosm>,
    /// (Reference) Earth frame.
    earth_frame: Frame,
    /// Sun Body frame
    sun_frame: Frame,
}

impl Solver {
    /// Evaluates Sun/Earth vector in meter ECEF at given Epoch
    fn sun_unit_vector(reference: &Frame, cosmic: &Arc<Cosm>, t: Epoch) -> Vector3<f64> {
        let orbit =
            cosmic.celestial_state(Bodies::Sun.ephem_path(), t, *reference, LightTimeCalc::None);
        Vector3::new(
            orbit.x_km * 1000.0,
            orbit.y_km * 1000.0,
            orbit.z_km * 1000.0,
        )
    }
    /// Minimum SV number to gather
    fn min_sv_required(solution: PVTSolutionType, fixed_alt: bool) -> usize {
        match solution {
            PVTSolutionType::TimeOnly => 1,
            _ => {
                if fixed_alt {
                    3
                } else {
                    4
                }
            },
        }
    }
    /// Builds new Position solver using given Configuration settings.
    /// [AprioriPosition] describes location this should converge to,
    /// usually resulting from a geodetic survey. If not known, simply use null.  
    /// solution_type: [PVTSolutionType] we want this solver to generate.
    pub fn new(cfg: Config, apriori: AprioriPosition, solutions_type: PVTSolutionType) -> Self {
        let cosmic = Cosm::de438();
        let sun_frame = cosmic.frame("Sun J2000");
        let earth_frame = cosmic.frame("EME2000");
        /*
         * print some infos on latched config
         */
        if cfg.method == Method::SPP && cfg.min_sv_sunlight_rate.is_some() {
            warn!("eclipse filter is not meaningful when using spp strategy");
        }
        if cfg.modeling.relativistic_path_range {
            warn!("relativistic path range modeling is not supported at the moment");
        }
        let interp_order = cfg.interp_order;
        Self {
            cfg,
            cosmic,
            apriori,
            sun_frame,
            earth_frame,
            solutions_type,
            state: State::default(),
            observations: Vec::<Observation>::with_capacity(64),
            clocks: Vec::<Clock>::with_capacity(interp_order),
            orbits: Vec::<Orbit>::with_capacity(interp_order),
        }
    }
    /// Runs signal quality filter
    fn signal_quality_filter(&mut self) {
        if let Some(min_snr) = self.cfg.min_snr {
            self.observations
                .retain(|obs| obs.snr_db.unwrap_or(-200.0) > min_snr);
        }
    }
    /// Runs attitude filter
    fn sv_attitude_filter(&mut self) {
        if let Some(min_elev) = self.cfg.min_sv_elev {
            let before = self.orbits.len();
            self.orbits.retain(|orb| orb.elevation > min_elev);
            let dropped = self.orbits.len() - before;
            debug!("dropped {} states for low elevation angles", dropped);
        }
    }
    /// Performs orbital interpolation
    fn orbit_interpolation<O: OrbitIter>(&mut self, orbit: O) -> State {
        State::default()
    }

    /// Performs clock interpolation
    fn clock_interpolation<C: ClockIter>(&mut self, clock: C) -> State {
        State::default()
    }

    /// Try to resolve one PVTSolution by sourcing [Observation]s, [Clock] states and
    /// [Orbit]al states. Solver's behavior highly depends on [Config] preset and the desired [PVTSolutionType].
    /// Returns Ok(None) if Solver is not in position at the moment to resolve.
    pub fn resolve<O: ObservationIter, CK: ClockIter>(
        &mut self,
        mut observation: O,
        mut clock: CK,
    ) -> Result<Vec<PVTSolution>, Error> {
        let min_sv_required =
            Self::min_sv_required(self.solutions_type, self.cfg.fixed_altitude.is_some());

        let mut current_samp_t = Epoch::default();

        loop {
            debug!("{} - {:?}", self.cfg.method, self.state);
            match self.state {
                State::SignalAcquisition => {
                    /*
                     * Gather signal observations
                     */
                    loop {
                        if let Some(obs) = observation.next() {
                            debug!("{} ({}) - new observation", obs.epoch, obs.sv);
                            self.observations.push(obs);
                        }
                        // can always process one signal sampling at a time..
                        if self
                            .observations
                            .iter()
                            .map(|obs| obs.epoch)
                            .unique()
                            .count()
                            > 1
                        {
                            self.state = State::SignalFilter;
                            break;
                        }
                    }
                },
                State::SignalFilter => {
                    current_samp_t = self
                        .observations
                        .iter()
                        .min_by(|a, b| a.epoch.cmp(&b.epoch))
                        .unwrap()
                        .epoch;

                    self.signal_quality_filter();

                    // TODO:
                    // needs to be more complex for PPP
                    let count = self
                        .observations
                        .iter()
                        .filter(|obs| obs.epoch == current_samp_t)
                        .count();

                    if count >= min_sv_required {
                        debug!(
                            "{} - {} observations passed quality checks",
                            current_samp_t, count
                        );
                        self.state = State::ClockInterpolation;
                    } else {
                        warn!(
                            "{} - too many samples below quality criteria",
                            current_samp_t
                        );
                        self.state = State::SignalAcquisition;
                    }
                },
                State::ClockInterpolation => {
                    while let Some(ck) = clock.next() {
                        self.clocks.push(ck);
                    }
                },
                _ => return Err(Error::UnfitObservations),
            }
        }
    }
    //    let (x0, y0, z0) = (
    //        self.apriori.ecef.x,
    //        self.apriori.ecef.y,
    //        self.apriori.ecef.z,
    //    );

    //    let (lat_ddeg, lon_ddeg, altitude_above_sea_m) = (
    //        self.apriori.geodetic.x,
    //        self.apriori.geodetic.y,
    //        self.apriori.geodetic.z,
    //    );

    //    let method = self.cfg.method;
    //    let modeling = self.cfg.modeling;
    //    let solver_opts = &self.cfg.solver;
    //    let filter = solver_opts.filter;
    //    let interp_order = self.cfg.interp_order;

    //    let _cosmic = &self.cosmic;
    //    let _earth_frame = self.earth_frame;

    //    /* interpolate positions */
    //    let mut pool: Vec<Candidate> = pool
    //        .iter()
    //        .filter_map(|c| {
    //            match c.transmission_time(&self.cfg) {
    //                Ok((t_tx, dt_ttx)) => {
    //                    debug!("{:?} ({}) : signal travel time: {}", t_tx, c.sv, dt_ttx);
    //                    if let Some(mut interpolated) =
    //                        (self.interpolator)(t_tx, c.sv, interp_order)
    //                    {
    //                        let mut c = c.clone();

    //                        let rot = match modeling.earth_rotation {
    //                            true => {
    //                                const EARTH_OMEGA_E_WGS84: f64 = 7.2921151467E-5;
    //                                let we = EARTH_OMEGA_E_WGS84 * dt_ttx;
    //                                let (we_cos, we_sin) = (we.cos(), we.sin());
    //                                Matrix3::<f64>::new(
    //                                    we_cos, we_sin, 0.0_f64, -we_sin, we_cos, 0.0_f64, 0.0_f64,
    //                                    0.0_f64, 1.0_f64,
    //                                )
    //                            },
    //                            false => Matrix3::<f64>::new(
    //                                1.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64, 0.0_f64, 0.0_f64,
    //                                0.0_f64, 1.0_f64,
    //                            ),
    //                        };

    //                        interpolated.position = InterpolatedPosition::AntennaPhaseCenter(
    //                            rot * interpolated.position(),
    //                        );

    //                        /* determine velocity */
    //                        if let Some((p_ttx, p_position)) = self.prev_sv_state.get(&c.sv) {
    //                            let dt = (t_tx - *p_ttx).to_seconds();
    //                            interpolated.velocity =
    //                                Some((rot * interpolated.position() - rot * p_position) / dt);
    //                        }

    //                        self.prev_sv_state
    //                            .insert(c.sv, (t_tx, interpolated.position()));

    //                        if modeling.relativistic_clock_bias {
    //                            /*
    //                             * following calculations need inst. velocity
    //                             */
    //                            if interpolated.velocity.is_some() {
    //                                let _orbit = interpolated.orbit(t_tx, self.earth_frame);
    //                                const EARTH_SEMI_MAJOR_AXIS_WGS84: f64 = 6378137.0_f64;
    //                                const EARTH_GRAVITATIONAL_CONST: f64 = 3986004.418 * 10.0E8;
    //                                let orbit = interpolated.orbit(t_tx, self.earth_frame);
    //                                let ea_rad = deg2rad(orbit.ea_deg());
    //                                let gm = (EARTH_SEMI_MAJOR_AXIS_WGS84
    //                                    * EARTH_GRAVITATIONAL_CONST)
    //                                    .sqrt();
    //                                let bias = -2.0_f64 * orbit.ecc() * ea_rad.sin() * gm
    //                                    / SPEED_OF_LIGHT
    //                                    / SPEED_OF_LIGHT
    //                                    * Unit::Second;
    //                                debug!(
    //                                    "{:?} ({}) : relativistic clock bias: {}",
    //                                    t_tx, c.sv, bias
    //                                );
    //                                c.clock_corr += bias;
    //                            }
    //                        }

    //                        debug!(
    //                            "{:?} ({}) : interpolated state: {:?}",
    //                            t_tx, c.sv, interpolated.position
    //                        );

    //                        c.state = Some(interpolated);
    //                        Some(c)
    //                    } else {
    //                        warn!("{:?} ({}) : interpolation failed", t_tx, c.sv);
    //                        None
    //                    }
    //                },
    //                Err(e) => {
    //                    error!("{} - transsmision time error: {:?}", c.sv, e);
    //                    None
    //                },
    //            }
    //        })
    //        .collect();

    //    /* apply elevation filter (if any) */
    //    if let Some(min_elev) = self.cfg.min_sv_elev {
    //        let mut idx: usize = 0;
    //        let mut nb_removed: usize = 0;
    //        while idx < pool.len() {
    //            if let Some(state) = pool[idx - nb_removed].state {
    //                if state.elevation < min_elev {
    //                    debug!(
    //                        "{:?} ({}) : below elevation mask",
    //                        pool[idx - nb_removed].t,
    //                        pool[idx - nb_removed].sv
    //                    );
    //                    let _ = pool.swap_remove(idx - nb_removed);
    //                    nb_removed += 1;
    //                }
    //            }
    //            idx += 1;
    //        }
    //    }

    //    /* remove observed signals above snr mask (if any) */
    //    if let Some(min_snr) = self.cfg.min_snr {
    //        let mut nb_removed: usize = 0;
    //        for idx in 0..pool.len() {
    //            let (init_code, init_phase) = (
    //                pool[idx - nb_removed].code.len(),
    //                pool[idx - nb_removed].phase.len(),
    //            );
    //            pool[idx - nb_removed].min_snr_mask(min_snr);
    //            let delta_code = init_code - pool[idx - nb_removed].code.len();
    //            let delta_phase = init_phase - pool[idx - nb_removed].phase.len();
    //            if delta_code > 0 || delta_phase > 0 {
    //                debug!(
    //                    "{:?} ({}) : {} code | {} phase below snr mask",
    //                    pool[idx - nb_removed].t,
    //                    pool[idx - nb_removed].sv,
    //                    delta_code,
    //                    delta_phase
    //                );
    //            }
    //            /* make sure we're still compliant */
    //            match method {
    //                Method::SPP => {
    //                    if pool[idx - nb_removed].code.is_empty() {
    //                        debug!(
    //                            "{:?} ({}) dropped on bad snr",
    //                            pool[idx - nb_removed].t,
    //                            pool[idx - nb_removed].sv
    //                        );
    //                        let _ = pool.swap_remove(idx - nb_removed);
    //                        nb_removed += 1;
    //                    }
    //                },
    //                Method::PPP => {
    //                    let mut drop = !pool[idx - nb_removed].dual_freq_pseudorange();
    //                    drop |= !pool[idx - nb_removed].dual_freq_phase();
    //                    if drop {
    //                        let _ = pool.swap_remove(idx - nb_removed);
    //                        nb_removed += 1;
    //                    }
    //                },
    //            }
    //        }
    //    }

    //    /* apply eclipse filter (if need be) */
    //    if let Some(min_rate) = self.cfg.min_sv_sunlight_rate {
    //        let mut nb_removed: usize = 0;
    //        for idx in 0..pool.len() {
    //            let state = pool[idx - nb_removed].state.unwrap(); // infaillible
    //            let orbit = state.orbit(pool[idx - nb_removed].t, self.earth_frame);
    //            let state = eclipse_state(&orbit, self.sun_frame, self.earth_frame, &self.cosmic);
    //            let eclipsed = match state {
    //                EclipseState::Umbra => true,
    //                EclipseState::Visibilis => false,
    //                EclipseState::Penumbra(r) => r < min_rate,
    //            };
    //            if eclipsed {
    //                debug!(
    //                    "{:?} ({}): dropped - eclipsed by earth",
    //                    pool[idx - nb_removed].t,
    //                    pool[idx - nb_removed].sv
    //                );
    //                let _ = pool.swap_remove(idx - nb_removed);
    //                nb_removed += 1;
    //            }
    //        }
    //    }

    //    /* make sure we still have enough SV */
    //    let nb_candidates = pool.len();
    //    if nb_candidates < min_required {
    //        return Err(Error::NotEnoughFittingCandidates);
    //    } else {
    //        debug!("{:?}: {} elected sv", t, nb_candidates);
    //    }

    //    /* form matrix */
    //    let mut y = DVector::<f64>::zeros(nb_candidates);
    //    let mut g = MatrixXx4::<f64>::zeros(nb_candidates);
    //    let mut pvt_sv_data = HashMap::<SV, PVTSVData>::with_capacity(nb_candidates);

    //    let r_sun = Self::sun_unit_vector(&self.earth_frame, &self.cosmic, t);

    //    for (row_index, cd) in pool.iter().enumerate() {
    //        if let Ok(sv_data) = cd.resolve(
    //            t,
    //            &self.cfg,
    //            (x0, y0, z0),
    //            (lat_ddeg, lon_ddeg, altitude_above_sea_m),
    //            iono_bias,
    //            tropo_bias,
    //            row_index,
    //            &mut y,
    //            &mut g,
    //            &r_sun,
    //        ) {
    //            pvt_sv_data.insert(cd.sv, sv_data);
    //        }
    //    }

    //    let w = self.cfg.solver.weight_matrix(
    //        nb_candidates,
    //        pvt_sv_data
    //            .values()
    //            .map(|sv_data| sv_data.elevation)
    //            .collect(),
    //    );

    //    let (mut pvt_solution, new_state) = PVTSolution::new(
    //        g.clone(),
    //        w.clone(),
    //        y.clone(),
    //        pvt_sv_data.clone(),
    //        filter,
    //        self.filter_state.clone(),
    //    )?;

    //    let validator = SolutionValidator::new(&self.apriori.ecef, &pool, &w, &pvt_solution);

    //    let valid = validator.valid(solver_opts);
    //    if valid.is_err() {
    //        return Err(Error::InvalidatedSolution(valid.err().unwrap()));
    //    }

    //    if let Some((prev_t, prev_pvt)) = &self.prev_pvt {
    //        pvt_solution.vel = (pvt_solution.pos - prev_pvt.pos) / (t - *prev_t).to_seconds();
    //    }

    //    self.prev_pvt = Some((t, pvt_solution.clone()));

    //    if filter != Filter::None {
    //        self.filter_state = new_state;
    //    }

    //    /*
    //     * slightly rework the solution so it ""looks"" like
    //     * what we expect based on the defined setup.
    //     */
    //    if let Some(alt) = self.cfg.fixed_altitude {
    //        pvt_solution.pos.z = self.apriori.ecef.z - alt;
    //        pvt_solution.vel.z = 0.0_f64;
    //    }

    //    match solution {
    //        PVTSolutionType::TimeOnly => {
    //            pvt_solution.pos = Vector3::<f64>::default();
    //            pvt_solution.vel = Vector3::<f64>::default();
    //        },
    //        _ => {},
    //    }

    //    Ok((t, pvt_solution))
    //}
}
