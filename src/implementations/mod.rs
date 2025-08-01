use anise::prelude::{Almanac, Frame};

use crate::{
    bias::Bias,
    candidate::Candidate,
    cfg::Config,
    ephemeris::EphemerisSource,
    navigation::PVTSolution,
    orbit::OrbitSource,
    prelude::{Epoch, Error, Rc},
    rtk::RTKBase,
    solver::Solver,
    time::AbsoluteTime,
    user::UserParameters,
};

// kinematic implementations
mod kinematic;
pub use kinematic::KinematicSolver;

/// [StaticSolver] to resolve [PVTSolution]s without dynamics.
/// This solver is particularly suited for static applications.
/// It can still be used to locate moving targets, but the dynamics
/// of the system are not modeled nor predicted.
///
/// ## Generics:
/// - EPH: [EphemerisSource], possible Ephemeris provider.
/// - ORB: [OrbitSource], possible Orbit state provider.
/// - B: [Bias] model.
/// - TIM: [AbsoluteTime] source for correct absolute time
pub struct StaticSolver<EPH: EphemerisSource, ORB: OrbitSource, B: Bias, TIM: AbsoluteTime> {
    solver: Solver<EPH, ORB, B, TIM>,
}

impl<EPH: EphemerisSource, ORB: OrbitSource, B: Bias, TIM: AbsoluteTime>
    StaticSolver<EPH, ORB, B, TIM>
{
    /// Creates a new [StaticSolver] with possible initial position.
    ///
    /// ## Input
    /// - almanac: provided valid [Almanac]
    /// - earth_cef: [Frame] that must be an ECEF
    /// - cfg: solver [Config]uration
    /// - eph_source: [EphemerisSource] implementation, serves as "raw" / indirect
    /// orbit provider.
    /// - orbit_source: [OrbitSource] implementation to provide [Orbit]al states directly
    /// - time_source: [AbsoluteTime] implementation
    /// - bias: external [Bias] model implementation, to improve overall accuracy.
    /// - initial_ecef_m: possible initial state (ECEF, meters)
    ///
    /// Usually you want to use this with a position preset. Otherwise,
    /// you can simply move to [Self::new_survey].
    pub fn new(
        almanac: Almanac,
        earth_cef: Frame,
        cfg: Config,
        eph_source: Rc<EPH>,
        orbit_source: Rc<ORB>,
        time_source: TIM,
        bias: B,
        initial_ecef_m: Option<(f64, f64, f64)>,
    ) -> Self {
        Self {
            solver: {
                Solver::<EPH, ORB, B, TIM>::new(
                    almanac,
                    earth_cef,
                    cfg,
                    eph_source,
                    orbit_source,
                    time_source,
                    bias,
                    initial_ecef_m,
                )
            },
        }
    }

    /// Creates a new [StaticSolver] with no a apriori knowledge.
    /// In this case, the solver will have to initialize itself.
    ///
    /// ## Input
    /// - almanac: provided valid [Almanac]
    /// - earth_cef: [Frame] that must be an ECEF
    /// - cfg: solver [Config]uration
    /// - eph_source: [EphemerisSource] implementation, serves as "raw" / indirect
    /// orbit provider.
    /// - orbit_source: [OrbitSource] implementation to provide [Orbit]al states directly
    /// - time_source: [AbsoluteTime] implementation
    /// - bias: external [Bias] model implementation, to improve overall accuracy.
    pub fn new_survey(
        almanac: Almanac,
        earth_cef: Frame,
        cfg: Config,
        eph_source: Rc<EPH>,
        orbit_source: Rc<ORB>,
        time_source: TIM,
        bias: B,
    ) -> Self {
        Self {
            solver: {
                Solver::<EPH, ORB, B, TIM>::new(
                    almanac,
                    earth_cef,
                    cfg,
                    eph_source,
                    orbit_source,
                    time_source,
                    bias,
                    None,
                )
            },
        }
    }

    /// [PVTSolution] solving attempt using PPP technique (no reference).
    /// Use this when no [RTKBase] may be accessed.
    /// Switch to RTK at any point in your session, when at least one [RTKBase] becomes
    /// accessible.
    ///
    /// ## Input
    /// - epoch: [Epoch] of measurement
    /// - params: [UserParameters]
    /// - candidates: proposed [Candidate]s (= measurements)
    /// - rtk_base: possible [RTKBase] we will connect to
    ///
    /// ## Output
    /// - [PVTSolution].
    pub fn ppp_solving(
        &mut self,
        epoch: Epoch,
        params: UserParameters,
        candidates: &[Candidate],
    ) -> Result<PVTSolution, Error> {
        self.solver.ppp_solving(epoch, params, candidates)
    }

    /// [PVTSolution] solving attempt using RTK technique and a single reference
    /// site. Switch to PPP at any point in your session, when access to remote
    /// site is lost.
    ///
    /// ## Input
    /// - epoch: [Epoch] of measurement
    /// - params: [UserParameters]
    /// - candidates: proposed [Candidate]s (= measurements)
    /// - base: [RTKBase] implementation, that must provide enough information
    /// for this to proceed. You may catch RTK related issues and
    /// retry using PPP technique.
    ///
    /// ## Output
    /// - [PVTSolution].
    pub fn rtk_solving<RTK: RTKBase>(
        &mut self,
        epoch: Epoch,
        params: UserParameters,
        candidates: &[Candidate],
        base: &RTK,
    ) -> Result<PVTSolution, Error> {
        self.solver.rtk_solving(epoch, params, candidates, base)
    }

    /// Reset this [StaticSolver].
    pub fn reset(&mut self) {
        self.solver.reset();
    }
}
