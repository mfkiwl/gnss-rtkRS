        let mut e_tx = Epoch::from_duration(dt_tx * Unit::Second, ts);

        if cfg.modeling.sv_clock_bias {
            debug!("{:?} ({}) clock_corr: {}", t, self.sv, self.clock_corr);
            e_tx -= self.clock_corr;
        }

        if cfg.modeling.sv_total_group_delay {
            if let Some(tgd) = self.tgd {
                debug!("{:?} ({}) tgd   : {}", t, self.sv, tgd);
                e_tx -= tgd;
            }
        }
    }
    /*
     * Resolves Self
     */
    pub(crate) fn resolve(
        &self,
        t: Epoch,
        cfg: &Config,
        apriori: (f64, f64, f64),
        apriori_geo: (f64, f64, f64),
        iono_bias: &IonosphericBias,
        tropo_bias: &TroposphericBias,
        row_index: usize,
        y: &mut DVector<f64>,
        g: &mut MatrixXx4<f64>,
        r_sun: &Vector3<f64>,
    ) -> Result<PVTSVData, Error> {
        // state
        let state = self.state.ok_or(Error::UnresolvedState)?;
        let clock_corr = self.clock_corr.to_seconds();
        let (azimuth, elevation) = (state.azimuth, state.elevation);

        /*
         * compensate for ARP (if possible)
         */
        let apriori = match cfg.arp_enu {
            Some(offset) => (
                apriori.0 + offset.0,
                apriori.1 + offset.1,
                apriori.2 + offset.2,
            ),
            None => apriori,
        };

        let (x0, y0, z0) = apriori;

        /*
         * Compensate for APC (if desired)
         */
        let sv_p = match cfg.modeling.sv_apc {
            true => {
                match state.position {
                    InterpolatedPosition::MassCenter(r_sat) => {
                        let delta_apc = 0.0_f64; //TODO
                        let k = -r_sat
                            / (r_sat[0].powi(2) + r_sat[1].powi(2) + r_sat[3].powi(2)).sqrt();
                        let norm = ((r_sun[0] - r_sat[0]).powi(2)
                            + (r_sun[1] - r_sat[1]).powi(2)
                            + (r_sun[2] - r_sat[2]).powi(2))
                        .sqrt();

                        let e = (r_sun - r_sat) / norm;
                        let j = Vector3::<f64>::new(k[0] * e[0], k[1] * e[1], k[2] * e[2]);
                        let i = Vector3::<f64>::new(j[0] * k[0], j[1] * k[1], j[2] * k[2]);
                        let r_dot = Vector3::<f64>::new(
                            (i[0] + j[0] + k[0]) * delta_apc,
                            (i[1] + j[1] + k[1]) * delta_apc,
                            (i[2] + j[2] + k[2]) * delta_apc,
                        );
                        r_sat + r_dot
                    },
                    InterpolatedPosition::AntennaPhaseCenter(r_sat) => r_sat,
                }
            },
            false => state.position(),
        };

        let pr = self
            .prefered_pseudorange()
            .ok_or(Error::MissingPseudoRange)?;

        let (pr, frequency) = (pr.value, pr.frequency);
    }
}
