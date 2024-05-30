use nalgebra::{base::dimension::U8, DMatrix, DVector, OMatrix, OVector, Vector3};

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
    pub p: DMatrix<f64>,
    pub x: DVector<f64>,
}

#[derive(Debug, Clone)]
struct KFState {
    pub q: DMatrix<f64>,
    pub p: DMatrix<f64>,
    pub x: DVector<f64>,
    pub phi: DMatrix<f64>,
}

#[derive(Debug, Clone)]
pub enum FilterState {
    Lsq(LSQState),
    Kf(KFState),
}

impl FilterState {
    fn lsq(state: LSQState) -> Self {
        Self::Lsq(state)
    }
    fn kf(state: KFState) -> Self {
        Self::Kf(state)
    }
    pub(crate) fn estimate(&self) -> &DVector<f64> {
        match self {
            Self::Lsq(state) => &state.x,
            Self::Kf(state) => &state.x,
        }
    }
}

impl Filter {
    fn lsq_resolve(
        input: &Input,
        w: &DMatrix<f64>,
        p_state: Option<FilterState>,
    ) -> Result<Output, Error> {
        match p_state {
            Some(FilterState::Lsq(p_state)) => {
                let p_1 = p_state.p.try_inverse().ok_or(Error::MatrixInversionError)?;

                let g_prime = input.g.clone().transpose();
                let q = (g_prime.clone() * input.g.clone())
                    .try_inverse()
                    .ok_or(Error::MatrixInversionError)?;

                let p = g_prime.clone() * w * input.g.clone();
                let p = (p_1.clone() + p)
                    .try_inverse()
                    .ok_or(Error::MatrixInversionError)?;

                let x = p.clone() * (p_1 * p_state.x + (g_prime * w * input.y.clone()));

                Ok(Output {
                    gdop: (q[(0, 0)] + q[(1, 1)] + q[(2, 2)] + q[(3, 3)]).sqrt(),
                    pdop: (q[(0, 0)] + q[(1, 1)] + q[(2, 2)]).sqrt(),
                    tdop: q[(4, 3)].sqrt(),
                    q,
                    state: FilterState::lsq(LSQState { p, x }),
                })
            },
            _ => {
                let g_prime = input.g.clone().transpose();

                let q = (g_prime.clone() * input.g.clone())
                    .try_inverse()
                    .ok_or(Error::MatrixInversionError)?;

                let p = (g_prime.clone() * w * input.g.clone())
                    .try_inverse()
                    .ok_or(Error::MatrixInversionError)?;

                let x = p.clone() * (g_prime.clone() * w * input.y.clone());
                if x[3].is_nan() {
                    return Err(Error::TimeIsNan);
                }

                Ok(Output {
                    gdop: (q[(0, 0)] + q[(1, 1)] + q[(2, 2)] + q[(3, 3)]).sqrt(),
                    pdop: (q[(0, 0)] + q[(1, 1)] + q[(2, 2)]).sqrt(),
                    tdop: q[(4, 3)].sqrt(),
                    q,
                    state: FilterState::lsq(LSQState { p, x }),
                })
            },
        }
    }
    fn kf_resolve(
        input: &Input,
        w: &DMatrix<f64>,
        p_state: Option<FilterState>,
    ) -> Result<Output, Error> {
        match p_state {
            Some(FilterState::Kf(p_state)) => {
                let x_bn = p_state.phi.clone() * p_state.x;
                let p_bn = p_state.phi.clone() * p_state.p * p_state.phi.transpose() + p_state.q;

                let p_bn_inv = p_bn.try_inverse().ok_or(Error::MatrixInversionError)?;
                let p_n = (input.g.transpose() * w * input.g.clone() + p_bn_inv.clone())
                    .try_inverse()
                    .ok_or(Error::MatrixInversionError)?;

                let w_g = input.g.transpose() * w * input.y.clone();
                let w_gy_pbn = w_g + (p_bn_inv * x_bn);
                let x_n = p_n.clone() * w_gy_pbn;

                let q_n = input.g.transpose() * input.g.clone();

                let mut q_diag = DVector::<f64>::zeros(8);
                let mut phi_diag = DVector::<f64>::zeros(8);

                for i in 0..8 {
                    phi_diag[i] = 1.0;
                }
                q_diag[7] = 1.0;

                Ok(Output {
                    gdop: (q_n[(0, 0)] + q_n[(1, 1)] + q_n[(2, 2)] + q_n[(3, 3)]).sqrt(),
                    pdop: (q_n[(0, 0)] + q_n[(1, 1)] + q_n[(2, 2)]).sqrt(),
                    tdop: q_n[(4, 3)].sqrt(),
                    q: q_n,
                    state: FilterState::kf(KFState {
                        p: p_n,
                        x: x_n,
                        q: DMatrix::<f64>::from_diagonal(&q_diag),
                        phi: DMatrix::<f64>::from_diagonal(&phi_diag),
                    }),
                })
            },
            _ => {
                let g_prime = input.g.clone().transpose();
                let q = (g_prime.clone() * input.g.clone())
                    .try_inverse()
                    .ok_or(Error::MatrixInversionError)?;

                let p = (g_prime.clone() * w * input.g.clone())
                    .try_inverse()
                    .ok_or(Error::MatrixInversionError)?;

                let x = p.clone() * (g_prime.clone() * w * input.y.clone());
                if x[3].is_nan() {
                    return Err(Error::TimeIsNan);
                }

                let mut q_diag = DVector::<f64>::zeros(8);
                let mut phi_diag = DVector::<f64>::zeros(8);

                for i in 0..8 {
                    phi_diag[i] = 1.0;
                }
                q_diag[7] = 1.0;

                Ok(Output {
                    gdop: (q[(0, 0)] + q[(1, 1)] + q[(2, 2)] + q[(3, 3)]).sqrt(),
                    pdop: (q[(0, 0)] + q[(1, 1)] + q[(2, 2)]).sqrt(),
                    tdop: q[(4, 3)].sqrt(),
                    q,
                    state: FilterState::kf(KFState {
                        p,
                        x,
                        q: DMatrix::<f64>::from_diagonal(&q_diag),
                        phi: DMatrix::<f64>::from_diagonal(&phi_diag),
                    }),
                })
            },
        }
    }
    pub fn resolve(
        &self,
        input: &Input,
        w: &DMatrix<f64>,
        p_state: Option<FilterState>,
    ) -> Result<Output, Error> {
        match self {
            Filter::Kalman => Self::kf_resolve(input, w, p_state),
            Filter::LSQ | _ => Self::lsq_resolve(input, w, p_state),
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
