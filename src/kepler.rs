use hifitime::Epoch;
use map_3d::rad2deg;
use std::f64::consts::PI;

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
    /// Argument of perigee (semicircles)
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

/// Implement this trait to provide SV Ephemerides
pub trait EphemeridesIter {
    /// Provide Time of Issue of Ephemerides, in chronological order,
    /// expressed as Epoch in whatever timescale you want.
    fn next_toe(&mut self) -> Option<Epoch>;
    /// Provide SV Keplerian parameters that we would resolve
    /// to orbital state internally.
    /// Tie this to None if case you're solely using [OrbitIter].
    fn next_kepler(&mut self) -> Option<(Keplerian, Perturbations)>;
}
