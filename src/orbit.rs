use thiserror::Error;

use nyx::{cosmic::Orbit as NyxOrbit, md::prelude::Frame};

use crate::prelude::AprioriPosition;
use gnss::prelude::{Constellation, SV};
use hifitime::{Duration, Epoch, TimeScale, Unit, GPST_REF_EPOCH};

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
        apriori: &AprioriPosition,
    ) -> Self {
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
        // println!("GPST EPOCH: {:?} | TOE (GPST): {:?} ", epoch_gpst, toe);

        let t_k = (epoch_gpst - toe).to_seconds();
        //if t_k < 0.0 {
        //    panic!("kepler(): invalid toe {:?} should be past epoch {:?}", toe, epoch);
        //}

        let n0 = (Self::EARTH_GM_CONSTANT / keplerian.a.powf(3.0)).sqrt();
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

        let omega_k = keplerian.omega_0
            + (perturbations.omega_dot - Self::EARTH_OMEGA_E_WGS84) * t_k
            - Self::EARTH_OMEGA_E_WGS84 * toe_secs;

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
    pub fn new(
        sv: SV,
        epoch: Epoch,
        position: (f64, f64, f64),
        apriori: &AprioriPosition,
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

#[cfg(test)]
mod test {
    use crate::{
        orbit::{Keplerian, Orbit, Perturbations},
        prelude::AprioriPosition,
    };
    use gnss::prelude::{Constellation, SV};
    use hifitime::{Epoch, TimeScale};
    use nalgebra::Vector3;
    use std::str::FromStr;
    #[test]
    fn keplerian_orbit() {
        let sv = SV {
            constellation: Constellation::GPS,
            prn: 1,
        };
        for (epoch, week, secs, kepler, perturb, apriori, elev, azim, ecef) in [(
            Epoch::from_str("2022-01-01T00:00:00 UTC").unwrap(),
            2190,
            518400.0,
            Keplerian {
                a: 26561110.712759566,
                e: 0.00534839148168,
                i_0: 0.957537602313,
                omega_0: 1.03791041521,
                m_0: 2.30316624652,
                omega: -2.3834050415,
            },
            Perturbations {
                dn: 2.3949035344821167e-17,
                i_dot: 5.11807041192e-10,
                omega_dot: -8.0467641439e-09,
                cus: 6.09830021858e-06,
                cuc: 9.85339283943e-07,
                cis: -1.54599547386e-07,
                cic: -1.04308128357e-07,
                crs: 17.3125,
                crc: 258.34375,
            },
            AprioriPosition::from_ecef(Vector3::<f64>::new(3628427.9118, 562059.0936, 5197872.215)),
            8.386332281745226,
            133.44087594021298,
            (16685968.411769923, 20728763.631397538, -1574846.006229475),
        )] {
            let orbit = Orbit::kepler(sv, epoch, week, secs, kepler, perturb, &apriori);

            assert_eq!(orbit.sv, sv);

            let err = (
                (ecef.0 - orbit.position.0).abs(),
                (ecef.1 - orbit.position.1).abs(),
                (ecef.2 - orbit.position.2).abs(),
            );

            assert!(err.0 < 1.0E-6, "x(ecef) error too large {}", err.0);
            assert!(err.1 < 1.0E-6, "y(ecef) error too large {}", err.1);
            assert!(err.2 < 1.0E-6, "z(ecef) error too large {}", err.2);

            let elev_err = (orbit.elevation - elev).abs();
            let azim_err = (orbit.azimuth - azim).abs();

            assert!(elev_err < 1.0E-6, "sv_elev err too large {}", elev_err);
            assert!(azim_err < 1.0E-6, "sv_azim error too large {}", azim_err);

            //todo: test other API
            //let orbit = Orbit::position(position);
        }
    }
}
