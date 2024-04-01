use crate::prelude::Vector3;
use map_3d::{deg2rad, ecef2geodetic, geodetic2ecef, rad2deg, Ellipsoid};

#[derive(Default, Debug, Clone, PartialEq)]
pub struct AprioriPosition {
    /// ECEF coordinates in meters
    ecef: Vector3<f64>,
    /// Geodetic coordinates in radians
    geodetic: Vector3<f64>,
}

impl AprioriPosition {
    /// Returns Geodetic coordinates in decimal degrees
    pub fn geodetic_ddeg(&self) -> Vector3<f64> {
        Vector3::<f64>::new(
            rad2deg(self.geodetic[0]),
            rad2deg(self.geodetic[1]),
            rad2deg(self.geodetic[2]),
        )
    }
    /// Returns Geodetic coordinates in radians
    pub fn geodetic_rad(&self) -> Vector3<f64> {
        self.geodetic
    }
    /// Returns coordinates in ECEF [m]
    pub fn ecef(&self) -> Vector3<f64> {
        self.ecef
    }
    /// Builds Self from ECEF position [m]
    pub fn from_ecef(ecef: Vector3<f64>) -> Self {
        let (x, y, z) = (ecef[0], ecef[1], ecef[2]);
        let (lat, lon, h) = ecef2geodetic(x, y, z, Ellipsoid::WGS84);
        Self {
            ecef,
            geodetic: Vector3::new(lat, lon, h),
        }
    }
    /// Builds Self from Geodetic coordinates (latitude [ddeg], longitude [ddeg], altitude above sea [m])
    pub fn from_geo_ddeg(coords: Vector3<f64>) -> Self {
        let rad = Vector3::<f64>::new(deg2rad(coords[0]), deg2rad(coords[1]), coords[2]);
        Self::from_geo_rad(rad)
    }
    /// Builds Self from Geodetic coordinates (latitude [rad], longitude [rad], altitude above sea [m])
    pub fn from_geo_rad(coords: Vector3<f64>) -> Self {
        let (x, y, z) = geodetic2ecef(coords[0], coords[1], coords[2], Ellipsoid::WGS84);
        Self {
            geodetic: coords,
            ecef: Vector3::new(x, y, z),
        }
    }
}
