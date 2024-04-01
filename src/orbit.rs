use hifitime::{Duration, Epoch, TimeScale, Unit, GPST_REF_EPOCH};
use map_3d::rad2deg;
use nyx::{cosmic::Orbit as NyxOrbit, md::prelude::Frame};
use std::f64::consts::PI;

use crate::prelude::AprioriPosition;
use gnss::prelude::{Constellation, SV};

/// Implement this trait to provide SV Orbital States (direct positions).
pub trait OrbitIter {
    /// Provide Orbital positions in chronological order.
    /// In order to simplify internal logic, we will panic if Orbital states
    /// are not provided in chronological order.
    /// Unevenly spaced states (in time) will delay PVT solution production,
    /// as we maintain steady interpolation errors.
    /// Tie this to None in case you're solely using [EphemeridesIter].
    fn next(&mut self) -> Option<Orbit>;
}

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

impl Orbit {
    /// Eearth mass * Gravitationnal field constant [m^3/s^2]
    pub(crate) const EARTH_GM_CONSTANT: f64 = 3.986004418E14_f64;
    /// Earth rotation rate in WGS84 frame [rad]
    pub(crate) const EARTH_OMEGA_E_WGS84: f64 = 7.2921151467E-5;

    /// Computes Elevation and Azimuth angles between given position
    /// in the Sky and apriori position on ground.
    fn elevation_azimuth(position: (f64, f64, f64), apriori: &AprioriPosition) -> (f64, f64) {
        let (sv_x, sv_y, sv_z) = position;

        let ecef = apriori.ecef();
        let (ref_x, ref_y, ref_z) = (ecef[0], ecef[1], ecef[2]);

        let geodetic_rad = apriori.geodetic_rad();
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
    /// Build Self from given position that must be expressed in ECEF [m]
    pub fn new(sv: SV, epoch: Epoch, position: (f64, f64, f64), apriori: &AprioriPosition) -> Self {
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
