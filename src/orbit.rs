use map_3d::rad2deg;
use thiserror::Error;

use nyx::{cosmic::Orbit as NyxOrbit, md::prelude::Frame};

use crate::prelude::AprioriPosition;
use gnss::prelude::{Constellation, SV};
use hifitime::{Duration, Epoch, TimeScale, Unit, GPST_REF_EPOCH};
use std::f64::consts::PI;

/// Implement this trait to provide SV Orbital States
pub trait OrbitIter {
    /// Provide Orbital states in chronological order.
    fn next(&self) -> Option<Orbit>;
    /// If you're in position to provide Orbital states
    /// at Epoch we want them (without interpolation),
    /// use this entry that the internal core will prefer
    /// for performance issues. Otherwise, simply return
    /// None here.
    fn next_at(&self, t: Epoch) -> Option<Orbit>;
}

// /// Satellite Vehicle Position
// #[derive(Copy, Clone, Debug, PartialEq)]
// pub enum Position {
//     /// Mass Center position, in ECEF [m]
//     MassCenter(Vector3<f64>),
//     /// Antenna Phase Center position, in ECEF [m]
//     AntennaPhaseCenter(Vector3<f64>),
// }
//
// impl Default for Position {
//     fn default() -> Self {
//         Self::AntennaPhaseCenter(Vector3::<f64>::default())
//     }
// }
//
// impl Position {
//     /// Builds new [Position] from Mass Center (MC) position in ECEF [m]
//     pub fn mass_center(mc: (f64, f64, f64)) -> Self {
//         Self::MassCenter(Vector3::<f64>::new(mc.0, mc.1, mc.2))
//     }
//     /// Builds new [Position] from Antenna Phase Center (APC) position in ECEF [m]
//     pub fn antenna_phase_center(apc: (f64, f64, f64)) -> Self {
//         Self::AntennaPhaseCenter(Vector3::<f64>::new(apc.0, apc.1, apc.2))
//     }
// }

/// Satellite Vehicle Orbital state
#[derive(Debug, Clone)]
pub struct Orbit {
    /// Satellite Vehicle
    pub(crate) sv: SV,
    /// Instant of this State snapshot
    pub(crate) epoch: Epoch,
    /// 3D [Position] vector
    pub(crate) position: (f64, f64, f64),
    /// Azimuth angle
    pub(crate) azimuth: f64,
    /// Elevation angle
    pub(crate) elevation: f64,
}

/// Keplerian parameters
pub struct Keplerian {
    /// Semi major axis (m)
    pub a: f64,
    /// Eccentricity (n.a)
    pub e: f64,
    /// Inclination angle at reference time (semicircles)
    pub i_0: f64,
    /// Longitude of ascending node at reference time (semicircles)
    pub omega_0: f64,
    /// Mean anomaly at reference time (semicircles)
    pub m_0: f64,
    /// argument of perigee (semicircles)
    pub omega: f64,
}

/// Keplerian perturbations
pub struct Perturbations {
    /// Mean motion difference from computed value [semicircles.s-1]
    pub dn: f64,
    /// Inclination rate of change [semicircles.s-1]
    pub i_dot: f64,
    /// Right ascension rate of change [semicircles.s^-1]
    pub omega_dot: f64,
    /// Amplitude of sine harmonic correction term of the argument
    /// of latitude [rad]
    pub cus: f64,
    /// Amplitude of cosine harmonic correction term of the argument
    /// of latitude [rad]
    pub cuc: f64,
    /// Amplitude of sine harmonic correction term of the angle of inclination [rad]
    pub cis: f64,
    /// Amplitude of cosine harmonic correction term of the angle of inclination [rad]
    pub cic: f64,
    /// Amplitude of sine harmonic correction term of the orbit radius [m]
    pub crs: f64,
    /// Amplitude of cosine harmonic correction term of the orbit radius [m]
    pub crc: f64,
}

impl Orbit {
    /// Computes Elevation and Azimuth angles between given position
    /// in the Sky and apriori position on ground.
    fn elevation_azimuth(position: (f64, f64, f64), apriori: AprioriPosition) -> (f64, f64) {
        let (sv_x, sv_y, sv_z) = position;

        let ecef = apriori.ecef();
        let (ref_x, ref_y, ref_z) = (ecef[0], ecef[1], ecef[2]);

        let geodetic_rad = apriori.geodetic();
        let (ref_lat, ref_lon) = (geodetic_rad[0], geodetic_rad[1]);

        // || sv - ref_pos || pseudo range
        let a_i = (sv_x - ref_x, sv_y - ref_y, sv_z - ref_z);
        let norm = (a_i.0.powf(2.0) + a_i.1.powf(2.0) + a_i.2.powf(2.0)).sqrt();
        let a_i = (a_i.0 / norm, a_i.1 / norm, a_i.2 / norm);

        // ECEF to VEN 3X3 transform matrix
        let ecef_to_ven = (
            (
                ref_lat.cos() * ref_lon.cos(),
                ref_lat.cos() * ref_lon.sin(),
                ref_lat.sin(),
            ),
            (-ref_lon.sin(), ref_lon.cos(), 0.0_f64),
            (
                -ref_lat.sin() * ref_lon.cos(),
                -ref_lat.sin() * ref_lon.sin(),
                ref_lat.cos(),
            ),
        );
        // ECEF to VEN transform
        let ven = (
            (ecef_to_ven.0 .0 * a_i.0 + ecef_to_ven.0 .1 * a_i.1 + ecef_to_ven.0 .2 * a_i.2),
            (ecef_to_ven.1 .0 * a_i.0 + ecef_to_ven.1 .1 * a_i.1 + ecef_to_ven.1 .2 * a_i.2),
            (ecef_to_ven.2 .0 * a_i.0 + ecef_to_ven.2 .1 * a_i.1 + ecef_to_ven.2 .2 * a_i.2),
        );
        let el = rad2deg(PI / 2.0 - ven.0.acos());
        let mut az = rad2deg(ven.1.atan2(ven.2));
        if az < 0.0 {
            az += 360.0;
        }
        (el, az)
    }
    /// Generates Epoch of TOE
    fn toe_gpst(week: u16, secs: f64) -> Epoch {
        Epoch::from_duration(
            (week as f64) * 7.0 * Unit::Day + secs * Unit::Second,
            TimeScale::GPST,
        )
    }
    /// Builds new [Orbit] from Keplerian parameters estimated at given Epoch.
    /// TOE: Time of Issue of Ephemeris, should be past "epoch" otherwise this will panic.
    /// "toe" week is rolling counter in said Timescale.
    /// "toe_secs" is week elapsed within that week, of said Timescale.
    pub fn kepler(
        sv: SV,
        epoch: Epoch,
        toe_week: u16,
        toe_secs: f64,
        keplerian: Keplerian,
        perturbations: Perturbations,
        apriori: AprioriPosition,
    ) -> Self {
        /// Eearth mass * Gravitationnal field constant [m^3/s^2]
        const EARTH_GM_CONSTANT: f64 = 3.986004418E14_f64;
        /// Earth rotation rate in WGS84 frame [rad]
        const EARTH_OMEGA_E_WGS84: f64 = 7.2921151467E-5;

        let ts = sv
            .timescale()
            .unwrap_or_else(|| panic!("kepler(): failed to determine timescale"));

        let epoch_gpst = match epoch.time_scale {
            TimeScale::UTC => Epoch::from_duration(
                epoch - GPST_REF_EPOCH - epoch.leap_seconds(false).unwrap_or(0.0) * Unit::Second
                    + 19.0 * Unit::Second,
                TimeScale::GPST,
            ),
            _ => epoch,
        };

        let toe_week_gpst = match ts {
            TimeScale::GST => toe_week - 1024,
            _ => toe_week,
        };

        let toe = Self::toe_gpst(toe_week_gpst, toe_secs);

        println!("GPST EPOCH: {:?} | TOE (GPST): {:?} ", epoch_gpst, toe);

        let t_k = (epoch_gpst - toe).to_seconds();
        //if t_k < 0.0 {
        //    panic!("kepler(): invalid toe {:?} should be past epoch {:?}", toe, epoch);
        //}

        let n0 = (EARTH_GM_CONSTANT / keplerian.a.powf(3.0)).sqrt();
        let n = n0 + perturbations.dn;
        let m_k = keplerian.m_0 + n * t_k;
        let e_k = m_k + keplerian.e * m_k.sin();
        let nu_k =
            ((1.0 - keplerian.e.powf(2.0)).sqrt() * e_k.sin()).atan2(e_k.cos() - keplerian.e);
        let phi_k = nu_k + keplerian.omega;

        let du_k =
            perturbations.cuc * (2.0 * phi_k).cos() + perturbations.cus * (2.0 * phi_k).sin();
        let u_k = phi_k + du_k;

        let di_k =
            perturbations.cic * (2.0 * phi_k).cos() + perturbations.cis * (2.0 * phi_k).sin();
        let i_k = keplerian.i_0 + perturbations.i_dot * t_k + di_k;

        let dr_k =
            perturbations.crc * (2.0 * phi_k).cos() + perturbations.crs * (2.0 * phi_k).sin();
        let r_k = keplerian.a * (1.0 - keplerian.e * e_k.cos()) + dr_k;

        let omega_k = keplerian.omega_0 + (perturbations.omega_dot - EARTH_OMEGA_E_WGS84) * t_k
            - EARTH_OMEGA_E_WGS84 * toe_secs;

        let xp_k = r_k * u_k.cos();
        let yp_k = r_k * u_k.sin();

        let x_k = xp_k * omega_k.cos() - yp_k * omega_k.sin() * i_k.cos();
        let y_k = xp_k * omega_k.sin() + yp_k * omega_k.cos() * i_k.cos();
        let z_k = yp_k * i_k.sin();

        //let position = Position::mass_center((x_k, y_k, z_k));
        let position = (x_k, y_k, z_k);
        let (elevation, azimuth) = Self::elevation_azimuth(position, apriori);
        Self {
            sv,
            epoch,
            azimuth,
            position,
            elevation,
        }
    }
    /// Build Self from given position that must be expressed in ECEF [m]
    pub fn position(
        sv: SV,
        epoch: Epoch,
        position: (f64, f64, f64),
        apriori: AprioriPosition,
    ) -> Self {
        let (elevation, azimuth) = Self::elevation_azimuth(position, apriori);
        Self {
            sv,
            epoch,
            position,
            azimuth,
            elevation,
        }
    }
    /// Creates Nyx [Orbit] from Self
    pub(crate) fn to_nyx(
        &self,
        velocity: Option<(f64, f64, f64)>,
        dt: Epoch,
        frame: Frame,
    ) -> NyxOrbit {
        let (x, y, z) = self.position;
        let (v_x, v_y, v_z) = velocity.unwrap_or((0.0, 0.0, 0.0));
        NyxOrbit::cartesian(
            x / 1000.0,
            y / 1000.0,
            z / 1000.0,
            v_x / 1000.0,
            v_y / 1000.0,
            v_z / 1000.0,
            dt,
            frame,
        )
    }
}