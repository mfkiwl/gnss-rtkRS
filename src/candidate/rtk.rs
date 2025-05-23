use crate::prelude::{Bias, Candidate, Config, Error, Method, Vector3};

impl Candidate {
    /// Measurement vector contribution.
    pub(crate) fn rtk_vector_contribution<B: Bias>(
        &self,
        cfg: &Config,
        remote: &Self,
    ) -> Result<(f64, f64, f64), Error> {
        let mut rho_m = match cfg.method {
            Method::SPP => {
                let (_, range_m) = self.best_snr_range_m().ok_or(Error::MissingPseudoRange)?;

                range_m
            },
            Method::PPP | Method::CPP => {
                let combination = self
                    .code_if_combination()
                    .ok_or(Error::PseudoRangeCombination)?;

                combination.value
            },
        };

        match cfg.method {
            Method::SPP => {
                let (_, range_m) = remote.best_snr_range_m().ok_or(Error::MissingPseudoRange)?;
                rho_m -= range_m;
            },
            Method::PPP | Method::CPP => {
                let combination = remote
                    .code_if_combination()
                    .ok_or(Error::PseudoRangeCombination)?;

                rho_m -= combination.value;
            },
        }

        Ok((rho_m, 1.0, 0.0))
    }

    /// RTK Matrix contribution.
    pub(crate) fn rtk_matrix_contribution(
        &self,
        //x0_y0_z0_m: Vector3<f64>,
        base_r0: Vector3<f64>,
    ) -> Vector3<f64> {
        //let (x0_m, y0_m, z0_m) = (x0_y0_z0_m[0], x0_y0_z0_m[1], x0_y0_z0_m[2]);

        let orbit = self.orbit.unwrap_or_else(|| {
            panic!("internal error: matrix contribution prior vector contribution")
        });

        let pos_vel_m = orbit.to_cartesian_pos_vel() * 1.0E3;

        let sv_r = Vector3::new(pos_vel_m[0], pos_vel_m[1], pos_vel_m[2]);

        let dr = sv_r - base_r0;
        let rho = dr.norm();
        dr / rho
    }
}
