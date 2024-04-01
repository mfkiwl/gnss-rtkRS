//! PVT solver

use log::{debug, error, info, warn};
use map_3d::deg2rad;
use std::collections::HashMap;
use thiserror::Error;

use gnss::prelude::{Constellation, SV};
use hifitime::{Duration, Epoch, TimeScale, Unit};

use nalgebra::{DMatrix, DVector, Matrix3, Matrix4, Matrix4x1, MatrixXx4, Vector3};

use nyx::{
    cosmic::{
        eclipse::{eclipse_state, EclipseState},
        SPEED_OF_LIGHT,
    },
    md::prelude::{Arc, Bodies, Cosm, Frame, LightTimeCalc},
};

use crate::{
    apriori::AprioriPosition,
    bias::{IonosphericBias, TroposphericBias},
    cfg::{Config, Filter, Method},
    clock::{Clock, ClockIter},
    ephemerides::EphemeridesIter,
    interp::{ClockInterpolator, OrbitInterpolator},
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
    SignalAcquisition,
    /// Signal quality filter
    SignalFilter,
    #[default]
    EphemeridesGathering,
    /// Clock data gathering
    ClockGathering,
    /// Orbital state gathering
    OrbitGathering,
    /// SV clock interpolation
    ClockInterpolation,
    /// Orbital state interpolation
    OrbitInterpolation,
    /// Navigation
    Navigation,
    /// Solution Validation
    Validation,
}

/// PVT Solver
#[derive(Debug)]
pub struct Solver {
    /// Solver configuration
    cfg: Config,
    /// Solver state
    state: State,
    /// Apriori position
    apriori: AprioriPosition,
    /// Minimum required
    min_sv_required: usize,
    /// Type of solutions to resolve
    solutions_type: PVTSolutionType,
    /// Clock interpolator
    clock: ClockInterpolator,
    /// Orbit interpolator
    orbit: OrbitInterpolator,
    /// Time of issue of Ephemeris
    toe: Vec<Epoch>,
    /// Observations buffer
    signals: Vec<Observation>,
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
        if interp_order % 2 == 0 {
            panic!("only odd interpolation orders are currently supported");
        }

        let min_sv_required = Self::min_sv_required(solutions_type, cfg.fixed_altitude.is_some());

        Self {
            cfg,
            cosmic,
            apriori,
            sun_frame,
            earth_frame,
            solutions_type,
            min_sv_required,
            state: State::default(),
            toe: Vec::<Epoch>::with_capacity(16),
            signals: Vec::<Observation>::with_capacity(64),
            clock: ClockInterpolator::malloc(128),
            orbit: OrbitInterpolator::malloc(interp_order, 128),
        }
    }
    /// Runs signal quality filter
    fn signal_quality_filter(&mut self) {
        if let Some(min_snr) = self.cfg.min_snr {
            self.signals
                .retain(|obs| obs.snr_db.unwrap_or(-200.0) > min_snr);
        }
    }
    /// Returns TOE
    fn find_toe(&self, t: Epoch) -> Option<&Epoch> {
        self.toe.iter().filter(|toe| **toe > t).reduce(|k, _| k)
    }
    /// Code filter + total candidates filter
    fn signal_filter(&mut self) {
        match self.cfg.method {
            Method::SPP => {
                /*
                 * Retain 1 L1 PR for up to 4 SV
                 */
                let mut index = 0;
                let mut sv = Vec::<SV>::new();

                //TODO propose cfg::prefered frequency
                let freq_hz = 1575.42E6_f64;

                self.signals.retain(|sig| {
                    let retain = sig.frequency_hz == freq_hz
                        && index < self.min_sv_required
                        && !sv.contains(&sig.sv);
                    if retain {
                        index += 1;
                        sv.push(sig.sv);
                    }
                    retain
                });
            },
        }
    }
    // /// Runs attitude filter
    // fn sv_attitude_filter(&mut self) {
    //     if let Some(min_elev) = self.cfg.min_sv_elev {
    //         let before = self.orbits.len();
    //         self.orbits.retain(|orb| orb.elevation > min_elev);
    //         let dropped = self.orbits.len() - before;
    //         debug!("dropped {} states for low elevation angles", dropped);
    //     }
    // }
    /// Try to resolve PVTSolution by exploiting [Observation]s, [Clock] and
    /// [Orbit]al states sources. Solver's behavior highly depends on [Config] preset
    /// and the desired [PVTSolutionType].
    pub fn resolve<E: EphemeridesIter, O: ObservationIter, OR: OrbitIter, CK: ClockIter>(
        &mut self,
        mut ephemerides: E,
        mut orbit: OR,
        mut clock: CK,
        mut observation: O,
    ) -> Result<Vec<PVTSolution>, Error> {
        const WEEK_SECONDS: f64 = 604800.0;
        let mut solutions = Vec::<PVTSolution>::new();

        let apriori = self.apriori.clone();
        let apriori_ecef = apriori.ecef();
        let (x0, y0, z0) = (apriori_ecef[0], apriori_ecef[1], apriori_ecef[2]);

        let interp_ord = self.cfg.interp_order;
        let min_sv_required = self.min_sv_required;

        // interation logic
        let mut t_tx = Epoch::default();
        let mut t_rx = Epoch::default();

        // interpolated results
        let mut interpolated_ck = Vec::<Clock>::with_capacity(min_sv_required);
        let mut interpolated_orb = Vec::<Orbit>::with_capacity(min_sv_required);
        let mut clock_corr = HashMap::<SV, f64>::with_capacity(min_sv_required);

        // matrix calc.
        let mut y = DVector::<f64>::zeros(min_sv_required);
        let mut g = MatrixXx4::<f64>::zeros(min_sv_required);
        let mut w = DMatrix::<f64>::identity(min_sv_required, min_sv_required);

        loop {
            debug!("{} - {:?}", self.cfg.method, self.state);
            match self.state {
                State::EphemeridesGathering => {
                    while let Some(toe) = ephemerides.next() {
                        assert!(
                            toe.time_scale == TimeScale::GPST,
                            "We only tolerate TOE expressed in GPST at the moment.."
                        );
                        self.toe.push(toe);
                    }
                    self.state = State::ClockGathering;
                },
                State::ClockGathering => {
                    while let Some(ck) = clock.next() {
                        self.clock.new_clock(ck);
                    }
                    self.state = State::OrbitGathering;
                },
                State::OrbitGathering => {
                    while let Some(orb) = orbit.next() {
                        self.orbit.new_orbit(orb);
                    }
                    self.state = State::SignalAcquisition;
                },
                State::SignalAcquisition => loop {
                    if let Some(obs) = observation.next() {
                        debug!("{:?} ({}) - new observation", obs.epoch, obs.sv);
                        self.signals.push(obs);
                    } else {
                        info!("{:?} - consumed all signals", t_rx);
                        return Ok(solutions);
                    }
                    if self.signals.iter().map(|obs| obs.epoch).unique().count() > 1
                    // next Epoch is starting: we'll drop 1st SV..
                    {
                        // simple logic, but 1st SV of each observation Epoch is lost..
                        self.signals.pop(); //drop

                        let samp_t = self
                            .signals
                            .iter()
                            .min_by(|a, b| a.epoch.cmp(&b.epoch))
                            .unwrap()
                            .epoch;

                        if samp_t != t_rx {
                            // moving to new Epoch
                            t_rx = samp_t;
                            t_tx = samp_t;
                            interpolated_ck.clear();
                            interpolated_orb.clear();
                            info!("new Epoch {:?}", t_rx);
                        }

                        self.state = State::SignalFilter;
                        break;
                    }
                },
                State::SignalFilter => {
                    self.signal_filter();
                    self.signal_quality_filter();

                    let count = self.signals.len();
                    if count < min_sv_required {
                        warn!(
                            "{:?} - too many samples below quality criteria {}/{}",
                            t_rx, count, min_sv_required
                        );
                        self.signals.clear();
                        self.state = State::SignalAcquisition;
                    } else {
                        debug!("{:?} - {} observations passed quality checks", t_rx, count);

                        self.state = State::ClockInterpolation;
                    }
                },
                State::ClockInterpolation => {
                    for sv in self.signals.iter().map(|sig| sig.sv) {
                        if let Some(ck) = self.clock.interpolate(sv, t_rx) {
                            debug!("{:?} - Interpolated {:?}", t_rx, ck);
                            interpolated_ck.push(ck);
                        }
                    }

                    if interpolated_ck.len() < min_sv_required {
                        self.signals.clear();
                        interpolated_ck.clear();
                        interpolated_orb.clear();
                        self.state = State::SignalAcquisition;
                        error!("{:?} - too many interpolation failures", t_rx);
                    } else {
                        //TODO deviation criteria ?
                        self.state = State::OrbitInterpolation;
                    }
                },
                State::OrbitInterpolation => {
                    for (sv, pr) in self.signals.iter().map(|sig| (sig.sv, sig.value)) {
                        if self.cfg.modeling.signal_propagation {
                            let ts = t_rx.time_scale;
                            let seconds_ts = t_rx.to_duration().to_seconds();
                            let dt_tx = seconds_ts - pr / SPEED_OF_LIGHT;
                            t_tx = Epoch::from_duration(dt_tx * Unit::Second, ts);

                            if self.cfg.modeling.sv_clock_bias {
                                if let Some(toe) = self.find_toe(t_tx) {
                                    let clock =
                                        interpolated_ck.iter().find(|ck| ck.sv == sv).unwrap();

                                    assert!(
                                        clock.sv.constellation == Constellation::GPS,
                                        "Can only process GPS vehicles at the moment.."
                                    );

                                    let mut dt = (t_rx - *toe).to_seconds();
                                    if dt > WEEK_SECONDS / 2.0 {
                                        dt -= WEEK_SECONDS;
                                    } else if dt < -WEEK_SECONDS / 2.0 {
                                        dt += WEEK_SECONDS;
                                    }

                                    let correction = clock.offset
                                        + clock.drift * dt
                                        + clock.drift_rate * dt.powi(2);

                                    //TODO: il y a aussi une correction d'un effet relativiste là dessus

                                    clock_corr.insert(sv, correction);
                                    t_tx -= correction * Unit::Second;
                                    debug!(
                                        "{:?} ({}) - clock correction {}",
                                        t_rx,
                                        sv,
                                        Duration::from_seconds(correction)
                                    );
                                } else {
                                    error!(
                                        "{:?} ({}) - undetermined ephemeris: can't proceed",
                                        t_rx, sv
                                    );
                                }
                            }

                            if self.cfg.modeling.sv_total_group_delay {
                                //TODO
                            }

                            let dt = t_rx - t_tx;
                            let dt_secs = dt.to_seconds();
                            assert!(
                                dt_secs.is_sign_positive(),
                                "Physical non sense: RX {:?} prior TX {:?}",
                                t_rx,
                                t_tx,
                            );
                            assert!(
                                dt_secs < 0.1,
                                "Physical non sense: {:?} signal propagation looks suspicious",
                                dt
                            );
                            debug!("{:?} ({}) - signal propagation {}", t_rx, sv, dt);
                        }
                        if let Some(mut orb) = self.orbit.interpolate(sv, t_tx, &apriori) {
                            if self.cfg.modeling.earth_rotation {
                                let we = Orbit::EARTH_OMEGA_E_WGS84 * (t_rx - t_tx).to_seconds();
                                let (we_cos, we_sin) = (we.cos(), we.sin());
                                let rot = Matrix3::<f64>::new(
                                    we_cos, we_sin, 0.0_f64, -we_sin, we_cos, 0.0_f64, 0.0_f64,
                                    0.0_f64, 1.0_f64,
                                );

                                let rotated = rot
                                    * Vector3::<f64>::new(
                                        orb.position.0,
                                        orb.position.1,
                                        orb.position.2,
                                    );

                                orb.position = (rotated[0], rotated[1], rotated[2]);

                                //TODO: correct velocities as well
                            }

                            debug!("{:?} - Interpolated {:?}", t_rx, orb);
                            interpolated_orb.push(orb);
                        }
                    }

                    if interpolated_orb.len() < min_sv_required {
                        self.signals.clear();
                        interpolated_ck.clear();
                        interpolated_orb.clear();
                        self.state = State::SignalAcquisition;
                        error!("{:?} - too many interpolation failures", t_rx);
                    } else {
                        // attitude filter
                        if let Some(min) = self.cfg.min_elevation {
                            interpolated_orb.retain(|orb| orb.elevation > min);
                        }
                        if let Some(min) = self.cfg.min_azimuth {
                            interpolated_orb.retain(|orb| orb.azimuth > min);
                        }
                        if let Some(max) = self.cfg.max_azimuth {
                            interpolated_orb.retain(|orb| orb.azimuth < max);
                        }
                        if interpolated_orb.len() < min_sv_required {
                            self.signals.clear();
                            interpolated_ck.clear();
                            interpolated_orb.clear();
                            self.state = State::SignalAcquisition;
                            error!("{:?} - too many SV below attitude criteria", t_rx);
                        } else {
                            self.state = State::Navigation;
                        }
                    }
                },
                State::Navigation => {
                    for (index, signal) in self.signals.iter().enumerate() {
                        let orbit = interpolated_orb
                            .iter()
                            .find(|orb| orb.sv == signal.sv)
                            .unwrap();
                        let clock = interpolated_ck
                            .iter()
                            .find(|ck| ck.sv == signal.sv)
                            .unwrap();

                        let (sv_x, sv_y, sv_z) = orbit.position;
                        let rho = ((sv_x - x0).powi(2) + (sv_y - y0).powi(2) + (sv_z - z0).powi(2))
                            .sqrt();

                        g[(index, 0)] = (x0 - sv_x) / rho;
                        g[(index, 1)] = (y0 - sv_y) / rho;
                        g[(index, 2)] = (z0 - sv_z) / rho;
                        g[(index, 3)] = 1.0_f64;

                        let mut biases = 0.0_f64;
                        if self.cfg.modeling.sv_clock_bias {
                            if let Some((_, clock_corr)) =
                                clock_corr.iter().find(|(k, _)| **k == signal.sv)
                            {
                                biases -= clock_corr * SPEED_OF_LIGHT;
                            } else {
                                panic!(
                                    "{} ({}) - undetermined ephemeris: can't proceed",
                                    t_rx, orbit.sv
                                );
                            }
                        }

                        /*
                         * Possible delay compensation
                         */
                        if let Some(delay) = self.cfg.externalref_delay {
                            y[index] -= delay * SPEED_OF_LIGHT;
                        }

                        //let rtm = bias::RuntimeParam {
                        //    t,
                        //    elevation,
                        //    azimuth,
                        //    apriori_geo,
                        //    frequency,
                        //};

                        if self.cfg.modeling.tropo_delay {
                            //let bias = TroposphericBias::model(TropoModel::Niel, &rtm);
                            //debug!("{:?} : modeled tropo delay {:.3E}[m]", t, bias);
                            //biases += bias;
                            //sv_data.tropo_bias = PVTBias::modeled(bias);
                        }

                        if self.cfg.modeling.iono_delay {
                            //if let Some(bias) = iono_bias.bias(&rtm) {
                            //    debug!(
                            //        "{:?} : modeled iono delay (f={:.3E}Hz) {:.3E}[m]",
                            //        t, rtm.frequency, bias
                            //    );
                            //    biases += bias;
                            //    sv_data.iono_bias = PVTBias::modeled(bias);
                            //}
                        }

                        for delay in &self.cfg.int_delay {
                            if delay.frequency == signal.frequency_hz {
                                y[index] += delay.delay * SPEED_OF_LIGHT;
                            }
                        }

                        y[index] = signal.value - rho - biases;
                    }
                    debug!("{:?} - G: {} Y: {}", t_rx, g, y);

                    let mut pvt = match PVTSolution::new(t_rx, g.clone(), w.clone(), y.clone()) {
                        Ok(pvt) => {
                            debug!("{:?} - new {:?}", t_rx, pvt);
                            solutions.push(pvt);
                        },
                        Err(Error::TimeIsNan) => error!("time is nan"),
                        Err(e) => panic!("solver error {:?}", e),
                    };

                    self.signals.clear();
                    self.state = State::SignalAcquisition;
                },
                _ => return Err(Error::UnfitObservations),
            }
        }
    }

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
