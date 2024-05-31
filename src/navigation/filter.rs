use nalgebra::{DMatrix, DVector, Vector3};

#[cfg(feature = "serde")]
use serde::Deserialize;

use super::{Input, Output};
use crate::prelude::{Epoch, Error};

/// Navigation Filter.
#[derive(Default, Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(Deserialize))]
pub enum Filter {
    /// None: solver filter completely bypassed. Lighter calculations, no iterative behavior.
    None,
    #[default]
    /// LSQ Filter. Heavy computation.
    LSQ,
    /// Kalman Filter. Heavy+ computations. Compared to LSQ, the Kalman filter
    /// converges faster and has the ability to "improve" models
    Kalman,
}

#[derive(Debug, Clone)]
struct LSQState {
    /// number of iterations
    nth: u64,
    pub p: DMatrix<f64>,
    pub x: DVector<f64>,
}

impl LSQState {
    pub fn init() -> Self {
        Self {
            nth: 0,
            p: DMatrix::<f64>::zeros(8, 8),
            x: DVector::<f64>::zeros(8),
        }
    }
    pub fn update(&mut self, input: &Input, w: &DMatrix<f64>) -> Result<Output, Error> {
        let q = if self.nth == 0 {
            let g_prime = input.g.clone().transpose();
            self.p = g_prime.clone() * w * input.g.clone();
            self.x = self.p.clone() * (g_prime.clone() * w * input.y.clone());
            (g_prime.clone() * input.g.clone())
                .try_inverse()
                .ok_or(Error::MatrixInversionError)?
        } else {
            let p_1 = self
                .p
                .clone()
                .try_inverse()
                .ok_or(Error::MatrixInversionError)?;

            let g_prime = input.g.clone().transpose();
            let p = g_prime.clone() * w * input.g.clone();

            self.p = (p_1.clone() + p.clone())
                .try_inverse()
                .ok_or(Error::MatrixInversionError)?;

            self.x =
                self.p.clone() * (p_1 * self.x.clone() + (g_prime.clone() * w * input.y.clone()));

            (g_prime.clone() * input.g.clone())
                .try_inverse()
                .ok_or(Error::MatrixInversionError)?
        };
        self.nth += 1;
        Ok(Output {
            gdop: (q[(0, 0)] + q[(1, 1)] + q[(2, 2)] + q[(3, 3)]).sqrt(),
            pdop: (q[(0, 0)] + q[(1, 1)] + q[(2, 2)]).sqrt(),
            tdop: q[(4, 3)].sqrt(),
            q,
        })
    }
}

#[derive(Debug, Clone)]
struct KFState {
    pub nth: u64,
    pub q: DMatrix<f64>,
    pub p: DMatrix<f64>,
    pub x: DVector<f64>,
    pub phi: DMatrix<f64>,
}

impl KFState {
    pub fn init() -> Self {
        Self {
            nth: 0,
            x: DVector::<f64>::zeros(8),
            p: DMatrix::<f64>::zeros(8, 8),
            q: {
                let mut diag = DVector::<f64>::zeros(8);
                diag[7] = 1.0;
                DMatrix::<f64>::from_diagonal(&diag)
            },
            phi: {
                let mut diag = DVector::<f64>::zeros(8);
                for i in 0..8 {
                    diag[i] = 1.0;
                }
                DMatrix::<f64>::from_diagonal(&diag)
            },
        }
    }
    pub fn update(&mut self, input: &Input, w: &DMatrix<f64>) -> Result<Output, Error> {
        let x_bn = self.phi.clone() * self.x.clone();
        let p_bn = self.phi.clone() * self.p.clone() * self.phi.transpose() + self.q.clone();

        let p_bn_inv = p_bn.try_inverse().ok_or(Error::MatrixInversionError)?;

        let w_g = input.g.transpose() * w * input.y.clone();
        let w_gy_pbn = w_g + (p_bn_inv.clone() * x_bn);

        let q_n = input.g.transpose() * input.g.clone();

        // update
        self.p = (input.g.transpose() * w * input.g.clone() + p_bn_inv.clone())
            .try_inverse()
            .ok_or(Error::MatrixInversionError)?;

        self.x = self.p.clone() * w_gy_pbn;

        Ok(Output {
            gdop: (q_n[(0, 0)] + q_n[(1, 1)] + q_n[(2, 2)] + q_n[(3, 3)]).sqrt(),
            pdop: (q_n[(0, 0)] + q_n[(1, 1)] + q_n[(2, 2)]).sqrt(),
            tdop: q_n[(4, 3)].sqrt(),
            q: q_n,
        })
    }
}

#[derive(Debug, Clone)]
pub enum FilterState {
    Kf(KFState),
    Lsq(LSQState),
}

impl FilterState {
    pub fn init(filter: Filter) -> Self {
        match filter {
            Filter::Kalman => Self::Kf(KFState::init()),
            Filter::LSQ | Filter::None => Self::Lsq(LSQState::init()),
        }
    }
    pub(crate) fn estimate(&self) -> &DVector<f64> {
        match self {
            Self::Kf(state) => &state.x,
            Self::Lsq(state) => &state.x,
        }
    }
    pub fn update(&mut self, input: &Input, w: &DMatrix<f64>) -> Result<Output, Error> {
        assert!(input.g.nrows() == w.ncols(), "invalid w matrix formulation");
        assert!(input.g.nrows() == w.nrows(), "invalid w matrix formulation");
        match self {
            Self::Kf(kf) => kf.update(input, w),
            Self::Lsq(lsq) => lsq.update(input, w),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Copy, Default)]
pub(crate) struct State3D {
    pub t: Epoch,
    pub inner: Vector3<f64>,
}

impl std::fmt::LowerExp for State3D {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        let (x, y, z) = (self.inner[0], self.inner[1], self.inner[2]);
        write!(f, "({:.6E},{:.6E},{:.6E})", x, y, z)
    }
}

impl std::fmt::Display for State3D {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        let (x, y, z) = (self.inner[0], self.inner[1], self.inner[2]);
        write!(f, "({:.6E},{:.6E},{:.6E})", x, y, z)
    }
}

// impl NyxState for State3D {
//     type Size = U3;
//     type VecLength = U3;
//     fn as_vector(&self) -> OVector<f64, U3>, NyxError {
//         self.inner.into()
//     }
//     fn unset_stm(&mut self) {}
//     fn set(&mut self, t: Epoch, vector: &OVector<f64, U3>) -> () {
//         self.t = t;
//         self.inner = vector.clone();
//     }
//     fn epoch(&self) -> Epoch {
//         self.t
//     }
//     fn set_epoch(&mut self, t: Epoch) {
//         self.t = t;
//     }
// }
