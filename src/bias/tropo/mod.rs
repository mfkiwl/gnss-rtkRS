use crate::cfg::Error;
use map_3d::deg2rad;

// mod unb3;

#[derive(Default, Copy, Clone, Debug)]
pub(crate) enum TroposphereModel {
    #[default]
    Niel,
    //UNB3,
}

impl TroposphereModel {
    pub fn bias(&self, elevation: f64, h: f64) -> f64 {
        match self {
            Self::Niel => Self::niel_model(elevation, h),
            //TropoModel::UNB3 => {
            //    let (zwd, zdd) = unb3::unb3_model(rtm);
            //    (zwd + zdd) * 1.001_f64
            //        / (0.002001_f64 + deg2rad(rtm.elevation).sin().powi(2)).sqrt()
            //},
        }
    }
    fn niel_model(elevation: f64, h: f64) -> f64 {
        const NS: f64 = 324.8;

        let elev = deg2rad(elevation);
        let h_km = h / 1000.0;

        let f = match elev < 90.0 {
            true => 1.0_f64 / (elev.sin() + 0.00143 / (elev.tan() + 0.0455)),
            false => 1.0,
        };

        let delta_n = -7.32 * (0.005577 * NS).exp();

        let delta_r =
            (NS + 0.5 * delta_n - NS * h_km - 0.5 * delta_n * h_km.powi(2) + 1430.0 + 732.0)
                * 0.001;

        f * delta_r
    }
}

impl std::str::FromStr for TroposphereModel {
    type Err = Error;
    fn from_str(s: &str) -> Result<TroposphereModel, Error> {
        let c = s.trim().to_lowercase();
        match c.as_str() {
            "niel" => Ok(TroposphereModel::Niel),
            //"unb3" => Ok(TroposphereModel::UNB3),
            _ => Err(Error::UnknownTropoModel(c.to_string())),
        }
    }
}

/// Tropospheric Components to attach any
/// resolution attempt. Fill as much as you can.
/// An empty structure will not impeach the solver to compensate
/// for this effect.
#[derive(Default, Copy, Clone, Debug)]
pub struct TroposphericBias {
    /// Undifferentiated total Zenith delay (Dry + Wet components), in meters of delay
    pub total: Option<f64>,
    /// Zenith Dry and Zenith Wet delay components, in meters of delay
    pub zwd_zdd: Option<(f64, f64)>,
}

impl TroposphericBias {
    pub(crate) fn needs_modeling(&self) -> bool {
        self.total.is_none() && self.zwd_zdd.is_none()
    }
    pub(crate) fn bias(&self, elevation: f64, h: f64) -> Option<f64> {
        if let Some((zwd, zdd)) = self.zwd_zdd {
            Some((zdd + zwd) * 1.001_f64 / (0.002001_f64 + deg2rad(elevation).sin().powi(2)).sqrt())
        } else {
            self.total.map(|total| {
                total * 1.001_f64 / (0.002001_f64 + deg2rad(elevation).sin().powi(2)).sqrt()
            })
        }
    }
}
