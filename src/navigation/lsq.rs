use log::debug;

use crate::{
    navigation::{solutions::PVTSolution, Error, Filter},
    prelude::Epoch,
};
use nalgebra::{Matrix4, Matrix4x1, Vector3, Vector4};
use nyx::cosmic::SPEED_OF_LIGHT;

#[derive(Debug, Clone)]
pub struct LSQ {
    init: bool,
    p: Matrix4<f64>,
    x: Matrix4x1<f64>,
}

impl Filter for LSQ {
    fn new() -> Self {
        Self {
            p: Matrix4::zeros(),
            x: Matrix4x1::zeros(),
            init: false,
        }
    }
    fn reset(&mut self) {
        self.p = Matrix4::zeros();
        self.x = Matrix4x1::zeros();
        self.init = false;
    }
    fn initialized(&self) -> bool {
        self.init
    }
    fn confirm(&mut self) {
        // self.x = x.clone();
        // self.p = p.clone();
        // self.init = true;
    }
    fn resolve(
        &mut self,
        t: Epoch,
        y: Vector4<f64>,
        g: Matrix4<f64>,
        w: Matrix4<f64>,
    ) -> PVTSolution {
        if !self.init {
            debug!("{:?} - \nG: {} Y: {}: W: {}", t, g, y, w);

            let g_prime = g.clone().transpose();
            let q = (g_prime.clone() * g.clone())
                .try_inverse()
                .unwrap_or_else(|| panic!("failed to invert matrix"));

            self.p = (g_prime.clone() * w.clone() * g.clone())
                .try_inverse()
                .unwrap_or_else(|| panic!("failed to invert matrix"));

            self.x = self.p * (g_prime.clone() * w.clone() * y);

            let dt = self.x[3] / SPEED_OF_LIGHT;
            assert!(!dt.is_nan(), "dt: {}", dt);

            PVTSolution::new(
                t,
                Vector3::new(self.x[0], self.x[1], self.x[2]),
                Vector3::<f64>::default(), //TODO vel
                dt,
                q,
            )
        } else {
            debug!(
                "{:?} - G: {} Y: {}: W: {} | P: {} x: {}",
                t, g, y, w, self.p, self.x
            );

            let p_1 = self
                .p
                .try_inverse()
                .unwrap_or_else(|| panic!("failed to matrix"));

            let g_prime = g.clone().transpose();
            let q = (g_prime.clone() * g.clone())
                .try_inverse()
                .unwrap_or_else(|| panic!("failed to matrix"));

            let p = g_prime.clone() * w.clone() * g.clone();
            let p = (p_1 + p)
                .try_inverse()
                .unwrap_or_else(|| panic!("failed to matrix"));

            let x = p * (p_1 * self.x + (g_prime.clone() * w.clone() * y));
            let dt = x[3] / SPEED_OF_LIGHT;
            assert!(!dt.is_nan(), "dt: {}", dt);

            PVTSolution::new(
                t,
                Vector3::new(x[0], x[1], x[2]),
                Vector3::<f64>::default(), //TODO vel
                dt,
                q,
            )
        }
    }
}
