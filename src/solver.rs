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
    bias::{tropo::TroposphereModel, BiasModel, IonosphereBiasModel, IonosphereBiasModelIter},
    cfg::{Config, Method},
    clock::{Clock, ClockIter},
    ephemerides::{Ephemerides, EphemeridesIter},
    interp::{ClockInterpolator, OrbitInterpolator},
    navigation::{
        lsq::LSQ,
        solutions::{PVTSolution, PVTSolutionType},
        Error, Filter, Navigation,
    },
    observation::{Observation, ObservationIter},
    orbit::{Orbit, OrbitIter},
    sv::{SVInfo, SVInfoIter},
};

use itertools::Itertools;

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
enum State {
    /// Need carrier signal measurements
    SignalAcquisition,
    /// Signal quality filter
    SignalFilter,
    #[default]
    EphemeridesGathering,
    /// SV Info gathering
    SVInfoGathering,
    /// Orbital state gathering
    OrbitGathering,
    /// [IonosphereBiasModel] gathering
    IonoModelsGathering,
    /// SV clock correction
    ClockCorrection,
    /// Orbital state interpolation
    OrbitInterpolation,
    /// Navigation
    Navigation,
    /// Solution Validation
    Validation,
    /// Aborting resolution
    Abort,
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
    /// Orbit interpolator
    orbit: OrbitInterpolator,
    /// Ionosphere Bias model(s)
    iono_models: Vec<IonosphereBiasModel>,
    /// Ephemerides
    eph: Vec<Ephemerides>,
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
        if !cfg.modeling.signal_propagation {
            if cfg.modeling.earth_rotation {
                warn!("cannot compensate Earth rotation without compensating signal propagation");
            }
        }

        let interp_order = cfg.interp_order;
        assert!(
            interp_order % 2 > 0,
            "only odd interpolation orders are supported"
        );

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
            eph: Vec::with_capacity(16),
            iono_models: Vec::with_capacity(8),
            signals: Vec::with_capacity(64),
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
    fn clock_correction(&self, sv: SV, t: Epoch) -> Option<Duration> {
        const WEEK_SECONDS: f64 = 604800.0;
        const HALF_WEEK_SECONDS :f64 = WEEK_SECONDS /2.0;
        let eph = self.eph.iter().filter(|eph| {
            eph.sv == sv && t >= eph.t && (t - eph.t).to_seconds() < 4000.0
        }).reduce(|k, _| k)?;
        debug!("IDENTIFIED: {:?}", eph);
        let mut dt = (t - eph.toe).to_seconds();
        // TODO: GPST
        if dt > HALF_WEEK_SECONDS {
            dt -= WEEK_SECONDS;
        } else if dt < -HALF_WEEK_SECONDS {
            dt += WEEK_SECONDS;
        }
        Some(Duration::from_seconds(eph.a0 + eph.a1 * dt + eph.a2 * dt.powi(2)))
    }
    /// Returns [IonosphereBiasModel] currently valid
    fn iono_bias_model(&self, t: Epoch) -> Option<&IonosphereBiasModel> {
        self.iono_models.iter().filter(|iono| iono.valid(t)).last()
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
                    let retain = sig.frequency == freq_hz
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
    /// Try to resolve PVTSolution by exploiting [Observation]s, [Clock] and
    /// [Orbit]al states sources. Solver's behavior highly depends on [Config] preset
    /// and the desired [PVTSolutionType].
    pub fn resolve<
        E: EphemeridesIter,
        O: ObservationIter,
        OR: OrbitIter,
        IO: IonosphereBiasModelIter,
        S: SVInfoIter,
    >(
        &mut self,
        mut ephemerides: E,
        mut orbit: OR,
        mut observation: O,
        mut sv_infos: S,
        mut ionosphere_models: IO,
    ) -> Vec<PVTSolution> {
        let mut solutions = Vec::<PVTSolution>::new();

        let apriori = self.apriori.clone();
        let apriori_ecef = apriori.ecef();
        let geo0_rad = apriori.geodetic_rad();

        let (x0, y0, z0) = (apriori_ecef[0], apriori_ecef[1], apriori_ecef[2]);
        let (lat0_rad, long0_rad, alt0) = (geo0_rad[0], geo0_rad[1], geo0_rad[2]);

        let filter = LSQ::new();
        let interp_ord = self.cfg.interp_order;
        let min_sv_required = self.min_sv_required;

        // interation logic
        let mut t_rx = Epoch::default();

        let mut pvt = PVTSolution::default();
        let mut past_pvt = Option::<PVTSolution>::None;
        let mut navigation = Navigation::new(filter);

        // interpolated results
        let mut sv_info = Vec::<SVInfo>::with_capacity(32);
        let mut interpolated_orb = Vec::<Orbit>::with_capacity(min_sv_required);
        let mut clock_corr = HashMap::<SV, Duration>::with_capacity(min_sv_required);

        loop {
            debug!("{} - {:?}", self.cfg.method, self.state);
            match self.state {
                State::EphemeridesGathering => {
                    while let Some(eph) = ephemerides.next() {
                        self.eph.push(eph);
                    }
                    debug!("DATASET: {:?}", self.eph);
                    // TODO: skip if PPP
                    self.state = State::IonoModelsGathering;
                },
                State::IonoModelsGathering => {
                    while let Some(model) = ionosphere_models.next() {
                        self.iono_models.push(model);
                    }
                    self.state = State::SVInfoGathering;
                },
                State::SVInfoGathering => {
                    while let Some(sv) = sv_infos.next() {
                        sv_info.push(sv);
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
                        return solutions;
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

                        self.state = State::ClockCorrection;
                    }
                },
                State::ClockCorrection => {
                    //TODO deviation criteria ?
                    self.state = State::OrbitInterpolation;
                },
                State::OrbitInterpolation => {
                    interpolated_orb.clear();

                    for (sv, pr) in self.signals.iter().map(|sig| (sig.sv, sig.value)) {
                        let mut t_tx = match self.cfg.modeling.signal_propagation {
                            false => t_rx,
                            true => {
                                let ts = t_rx.time_scale;
                                let seconds_ts = t_rx.to_duration().to_seconds();
                                let dt_tx_sec = seconds_ts - pr / SPEED_OF_LIGHT;
                                let t_tx =
                                    Epoch::from_duration(dt_tx_sec * Unit::Second, TimeScale::GPST); //TODO GPST
                                let dt = t_rx - t_tx;
                                let dt_sec = dt.to_seconds();
                                assert!(
                                    dt_sec.is_sign_positive(),
                                    "Physical non sense: RX {:?} prior TX {:?}",
                                    t_rx,
                                    t_tx,
                                );
                                assert!(
                                    dt_sec < 0.1,
                                    "Physical non sense: {} signal propagation looks suspicious",
                                    dt
                                );
                                debug!("{:?} ({}) - signal propagation {}", t_rx, sv, dt);
                                t_tx
                            },
                        };

                        if self.cfg.modeling.sv_total_group_delay {
                            if let Some(info) =
                                sv_info.iter().filter(|info| info.sv == sv).reduce(|k, _| k)
                            {
                                t_tx -= info.tgd;
                                debug!("{:?} ({}) - tgd: {}", t_rx, sv, info.tgd);
                            } else {
                                warn!(
                                    "{:?} ({}) - tgd not provided, therefore not compensated!",
                                    t_rx, sv
                                );
                            }
                        }

                        if self.cfg.modeling.sv_clock_bias {
                            if let Some(ck) = self.clock_correction(sv, t_tx) {
                                t_tx -= ck;
                                clock_corr.insert(sv, ck);
                                debug!(
                                    "{:?} ({}) - clock correction {}",
                                    t_rx,
                                    sv,
                                    ck,
                                );
                            } else {
                                error!(
                                    "{:?} ({}) - undetermined ephemeris: can't proceed",
                                    t_rx, sv
                                );
                            }
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
                                //TODO: velocity
                            }

                            debug!("{:?} ({}) - Interpolated {:?}", t_rx, sv, orb);
                            interpolated_orb.push(orb);
                        }
                    }

                    if interpolated_orb.len() < min_sv_required {
                        self.signals.clear();
                        self.state = State::SignalAcquisition;
                        error!("{:?} - too many (orbit) interpolation failures", t_rx);
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
                        // TODO: eclipse filter
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
                        if interpolated_orb.len() < min_sv_required {
                            self.signals.clear();
                            self.state = State::SignalAcquisition;
                            error!("{:?} - too many SV below attitude criteria", t_rx);
                        } else {
                            self.state = State::Navigation;
                        }
                    }
                },
                State::Navigation => {
                    for (index, signal) in self.signals.iter().enumerate() {
                        let mut biases = 0.0_f64;

                        let orbit = interpolated_orb
                            .iter()
                            .find(|orb| orb.sv == signal.sv)
                            .unwrap();

                        let (sv_x, sv_y, sv_z) = orbit.position;
                        let rho = ((sv_x - x0).powi(2) + (sv_y - y0).powi(2) + (sv_z - z0).powi(2))
                            .sqrt();

                        if self.cfg.modeling.sv_clock_bias {
                            if let Some((_, clock_corr)) =
                                clock_corr.iter().find(|(k, _)| **k == signal.sv)
                            {
                                biases -= clock_corr.to_seconds() * SPEED_OF_LIGHT;
                            } else {
                                error!(
                                    "{} ({}) - undetermined ephemeris: aborting",
                                    t_rx, orbit.sv
                                );
                                self.state = State::Abort;
                                break;
                            }
                        }

                        /*
                         * Possible delay compensation
                         */
                        if let Some(delay) = self.cfg.externalref_delay {
                            biases -= delay * SPEED_OF_LIGHT;
                        }

                        if self.cfg.modeling.tropo_delay {
                            let delay = TroposphereModel::Niel.bias(orbit.elevation, alt0);
                            debug!(
                                "{:?} ({}) - modeled tropo delay {:.3E}[m]",
                                t_rx, signal.sv, delay
                            );
                            biases += delay;
                        }

                        if self.cfg.modeling.iono_delay {
                            //TODO: only in SPP
                            if let Some(model) = self.iono_bias_model(t_rx) {
                                if let Some(bias) = model.bias(
                                    t_rx,
                                    signal.frequency,
                                    orbit.elevation,
                                    orbit.azimuth,
                                    lat0_rad,
                                    long0_rad,
                                ) {
                                    debug!(
                                        "{:?} ({}) - modeled iono delay (f={:.3E}Hz) {:.3E}[m]",
                                        t_rx, signal.sv, signal.frequency, bias
                                    );
                                    biases += bias;
                                    //sv_data.iono_bias = PVTBias::modeled(bias);
                                }
                            } else {
                                warn!("{:?} ({}) - no Ionosphere model!", t_rx, signal.sv);
                            }
                        }

                        for delay in &self.cfg.int_delay {
                            if delay.frequency == signal.frequency {
                                biases += delay.delay * SPEED_OF_LIGHT;
                            }
                        }

                        navigation.load(
                            index,
                            ((x0 - sv_x) / rho, (y0 - sv_y) / rho, (z0 - sv_z) / rho),
                            signal.value - rho - biases,
                        );
                    }
                    if self.state != State::Abort {
                        pvt = navigation.resolve(t_rx);
                        debug!("{:?} - new {:?}", t_rx, pvt);
                        self.state = State::Validation;
                    }
                },
                State::Validation => {
                    //if pvt.gdop() < 5.0 {
                    if let Some(ref past_pvt) = past_pvt {
                        let dt = (pvt.epoch - past_pvt.epoch).to_seconds();
                        pvt.vel = Vector3::<f64>::new(
                            (pvt.pos[0] - past_pvt.pos[0]) / dt,
                            (pvt.pos[1] - past_pvt.pos[1]) / dt,
                            (pvt.pos[2] - past_pvt.pos[2]) / dt,
                        );
                    }
                    solutions.push(pvt.clone());
                    past_pvt = Some(pvt.clone());
                    navigation.confirm();

                    debug!("{:?} - confirmed solution", t_rx);
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
                    //}
                    self.signals.clear();
                    self.state = State::SignalAcquisition;
                },
                State::Abort => {
                    self.signals.clear();
                    self.state = State::SignalAcquisition;
                },
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
}
