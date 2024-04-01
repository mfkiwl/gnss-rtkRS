use crate::bias::RuntimeParam;
use crate::prelude::TimeScale;
use hifitime::Epoch;
use map_3d::deg2rad;
use std::f64::consts::PI;

/// Ionosphere Bias source
pub trait IonosphereBiasModelIter {
    /// Provide Bias model(s)
    fn next(&mut self) -> Option<IonosphereBiasModel>;
}

/// Ionosphere bias.
#[derive(Debug, Default)]
pub struct IonosphereBiasModel {
    /// (Start of) Validity of model
    pub validity: Epoch,
    /// Klobuchar Model to apply
    pub kb_model: Option<KbModel>,
    /// Nequick-G model to apply
    pub ng_model: Option<NgModel>,
    /// BDGIM model to apply
    pub bd_model: Option<BdModel>,
    /// Slan Total Electron Density estimate
    pub stec_meas: Option<f64>,
}

/// Klobuchar Model
#[derive(Clone, Copy, Default, Debug)]
pub struct KbModel {
    /// alpha coefficients
    pub alpha: (f64, f64, f64, f64),
    /// beta coefficients
    pub beta: (f64, f64, f64, f64),
    /// ionosphere layer in [km] for this model
    pub h_km: f64,
}

impl KbModel {
    pub(crate) fn bias(
        &self,
        t: Epoch,
        frequency: f64,
        elevation: f64,
        azimuth: f64,
        lat_rad: f64,
        lon_rad: f64,
    ) -> Option<f64> {
        const PHI_P: f64 = 78.3;
        const R_EARTH: f64 = 6378.0;
        const LAMBDA_P: f64 = 291.0;
        const L1_F: f64 = 1575.42E6;

        let fract = R_EARTH / (R_EARTH + self.h_km);
        let (elev, azim) = (deg2rad(elevation), deg2rad(azimuth));
        let (phi_u, lambda_u) = (lat_rad, lon_rad);

        let t_gps = t.to_duration_in_time_scale(TimeScale::GPST).to_seconds();

        let psi = PI / 2.0 - elev - (fract * elev.cos()).asin();
        let phi_i = (phi_u.sin() * psi.cos() + phi_u.cos() * psi.sin() * azim.cos()).asin();
        let lambda_i = lambda_u + azim.sin() * psi / phi_i.cos();
        let phi_m = (phi_i.sin() * PHI_P.sin()
            + phi_i.cos() * PHI_P.cos() * (lambda_i - LAMBDA_P).cos())
        .asin();

        let mut t_s = 43.2E3 * lambda_i / PI + t_gps;
        if t_s > 86.4E3 {
            t_s -= 86.4E3;
        } else if t_s < 0.0 {
            t_s += 86.4E3;
        }

        let mut a_i = self.alpha.0 * (phi_m / PI).powi(0)
            + self.alpha.1 * (phi_m / PI).powi(1)
            + self.alpha.2 * (phi_m / PI).powi(2)
            + self.alpha.3 * (phi_m / PI).powi(3);
        if a_i < 0.0 {
            a_i = 0.0_f64;
        }

        let mut p_i = self.beta.0 * (phi_m / PI).powi(0)
            + self.beta.1 * (phi_m / PI).powi(1)
            + self.beta.2 * (phi_m / PI).powi(2)
            + self.beta.3 * (phi_m / PI).powi(3);
        if p_i < 72.0E3 {
            p_i = 72.0E3;
        }

        let x_i = 2.0 * PI * (t_s - 50400.0) / p_i;
        let f = 1.0 / ((1.0 - fract * elev.cos()).powi(2)).sqrt();
        let i_1 = match x_i < PI / 2.0 {
            true => 5.0 * 10E-9 + a_i * x_i.cos(),
            false => f * 5.0 * 10E-9,
        };

        Some(i_1 * (L1_F / frequency).powi(2))
    }
}

/// Nequick-G Model
#[derive(Clone, Copy, Default, Debug)]
pub struct NgModel {
    /// alpha coefficients
    pub a: (f64, f64, f64),
}

impl NgModel {
    pub(crate) fn bias(&self) -> Option<f64> {
        //let phi = deg2rad(rtm.apriori_geo.0);
        //let mu = inclination / phi.cos().sqrt();
        panic!("We only support the Klobuchar model at the moment..");
    }
}

/// BDGIM Model
#[derive(Clone, Copy, Default, Debug)]
pub struct BdModel {
    /// Alpha coefficients in TECu
    pub alpha: (f64, f64, f64, f64, f64, f64, f64, f64, f64),
}

impl BdModel {
    pub(crate) fn bias(&self) -> Option<f64> {
        //let phi = deg2rad(rtm.apriori_geo.0);
        //let mu = inclination / phi.cos().sqrt();
        panic!("We only support the Klobuchar model at the moment..");
    }
}

impl IonosphereBiasModel {
    /// Creates a new Klobuchar model
    pub fn klobuchar(start: Epoch, model: KbModel) -> Self {
        let mut s = Self::default();
        s.validity = start;
        s.kb_model = Some(model);
        s
    }
    /// Creates a new Nequick-G model
    pub fn nequick_g(start: Epoch, model: NgModel) -> Self {
        let mut s = Self::default();
        s.validity = start;
        s.ng_model = Some(model);
        s
    }
    /// Creates a new BDGIM model
    pub fn bdgim(start: Epoch, model: BdModel) -> Self {
        let mut s = Self::default();
        s.validity = start;
        s.bd_model = Some(model);
        s
    }
    /// Returns bias as meters of delay
    pub(crate) fn bias(
        &self,
        t: Epoch,
        frequency: f64,
        elevation: f64,
        azimuth: f64,
        lat_rad: f64,
        lon_rad: f64,
    ) -> Option<f64> {
        if let Some(_stec) = self.stec_meas {
            // TODO
            // let alpha = 40.3 * 10E16 / frequency / frequency;
            None
        } else if let Some(kb) = self.kb_model {
            kb.bias(t, frequency, elevation, azimuth, lat_rad, lon_rad)
        } else if let Some(ng) = self.ng_model {
            ng.bias()
        } else if let Some(bd) = self.bd_model {
            bd.bias()
        } else {
            None
        }
    }
    pub(crate) fn valid(&self, t: Epoch) -> bool {
        self.validity <= t
    }
}
