pub mod solutions;
//mod validator;
pub mod lsq;

use nalgebra::{Matrix4, Matrix4x1, Vector3, Vector4};

use crate::prelude::Epoch;
use log::debug;
use nyx::cosmic::SPEED_OF_LIGHT;
use solutions::PVTSolution;

#[derive(Debug, Clone)]
pub enum Error {
    TimeIsNan,
    MatrixInversion,
}

pub trait Filter {
    fn new() -> Self;
    fn reset(&mut self);
    fn initialized(&self) -> bool;
    fn confirm(&mut self, x: Matrix4x1<f64>, p: Matrix4<f64>);
    fn resolve(
        &mut self,
        t: Epoch,
        y: Vector4<f64>,
        g: Matrix4<f64>,
        w: Matrix4<f64>,
    ) -> PVTSolution;
}

#[derive(Debug, Clone)]
pub struct Navigation<F: Filter> {
    pub y: Vector4<f64>,
    pub g: Matrix4<f64>,
    pub w: Matrix4<f64>,
    p: Matrix4<f64>,
    x: Matrix4x1<f64>,
    pub filter: F,
}

impl<F: Filter> Navigation<F> {
    pub fn new(filter: F) -> Self {
        Self {
            filter,
            y: Vector4::<f64>::zeros(),
            g: Matrix4::<f64>::zeros(),
            w: Matrix4::<f64>::identity(),
            p: Matrix4::<f64>::zeros(),
            x: Matrix4x1::<f64>::zeros(),
        }
    }
    pub fn load(&mut self, row: usize, xyz: (f64, f64, f64), meas: f64) {
        self.y[row] = meas;
        self.g[(row, 0)] = xyz.0;
        self.g[(row, 1)] = xyz.1;
        self.g[(row, 2)] = xyz.2;
        self.g[(row, 3)] = 1.0_f64;
    }
    pub fn resolve(&mut self, t: Epoch) -> PVTSolution {
        //if self.filter.initialized() {
        //    self.filter.resolve(t, self.y, self.g, self.w)
        //} else {
        debug!("{:?} - G: {} Y: {}: W: {}", t, self.g, self.y, self.w);

        let g_prime = self.g.clone().transpose();
        let q = (g_prime.clone() * self.g.clone())
            .try_inverse()
            .unwrap_or_else(|| panic!("failed to invert matrix"));

        self.p = (g_prime.clone() * self.w.clone() * self.g.clone())
            .try_inverse()
            .unwrap_or_else(|| panic!("failed to invert matrix"));

        self.x = self.p * (g_prime.clone() * self.w.clone() * self.y);
        let dt = self.x[3] / SPEED_OF_LIGHT;
        assert!(!dt.is_nan(), "dt: {}", dt);

        PVTSolution::new(
            t,
            Vector3::new(self.x[0], self.x[1], self.x[2]),
            Vector3::<f64>::default(), //TODO vel
            dt,
            q,
        )
        //}
    }
    pub fn confirm(&mut self) {
        // self.filter.confirm(self.x, self.p);
    }
}
