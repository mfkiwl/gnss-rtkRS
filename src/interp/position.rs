use gnss::prelude::SV;
use hifitime::{Duration, Epoch};
use log::debug;
use std::collections::HashMap;

use crate::{interp::Interpolator, orbit::Orbit, prelude::AprioriPosition};

#[derive(Debug)]
pub struct OrbitInterpolator {
    order: usize,
    pub interpolators: HashMap<SV, PositionInterpolator>,
}

impl OrbitInterpolator {
    pub fn malloc(order: usize, size: usize) -> Self {
        Self {
            order,
            interpolators: HashMap::<SV, PositionInterpolator>::with_capacity(size),
        }
    }
    pub fn new_orbit(&mut self, orb: Orbit) {
        if let Some(interp) = self.interpolators.get_mut(&orb.sv) {
            interp.push((orb.epoch, (orb.position.0, orb.position.1, orb.position.2)));
        } else {
            let mut interp = PositionInterpolator::new(self.order);
            interp.push((orb.epoch, (orb.position.0, orb.position.1, orb.position.2)));
            self.interpolators.insert(orb.sv, interp);
        }
    }
    pub fn interpolate(&self, sv: SV, t_k: Epoch, apriori: &AprioriPosition) -> Option<Orbit> {
        let interp = self
            .interpolators
            .iter()
            .filter_map(|(k, v)| if *k == sv { Some(v) } else { None })
            .reduce(|k, _| k)?;
        let pos = interp.interpolate(t_k)?;
        Some(Orbit::new(sv, t_k, pos, apriori))
    }
}

/// Efficient Position Interpolator
#[derive(Debug)]
struct PositionInterpolator {
    order: usize,
    buffer: Vec<(Epoch, (f64, f64, f64))>,
}

impl Interpolator<(f64, f64, f64)> for PositionInterpolator {
    fn len(&self) -> usize {
        self.buffer.len()
    }
    fn get(&self, idx: usize) -> Option<&(Epoch, (f64, f64, f64))> {
        self.buffer.get(idx)
    }
    fn last(&self) -> Option<&(Epoch, (f64, f64, f64))> {
        self.buffer.last()
    }
    fn clear(&mut self) {
        self.buffer.clear();
    }
    fn push(&mut self, x_j: (Epoch, (f64, f64, f64))) {
        self.buffer.push(x_j);
    }
    fn interpolate(&self, t_s: Epoch) -> Option<(f64, f64, f64)> {
        if let Some(y_s) = self
            .buffer
            .iter()
            .filter_map(|(t, y)| if *t == t_s { Some(*y) } else { None })
            .reduce(|k, _| k)
        {
            Some(y_s)
        } else {
            let mut index = 0;
            for (idx, (t, _)) in self.buffer.iter().enumerate() {
                if *t > t_s {
                    break;
                }
                index = idx;
            }
            debug!("{}/{}", index, self.buffer.len());
            if index > (self.order + 1) / 2 && index < self.buffer.len() - (self.order + 1) / 2 {
                let offset = index - (self.order + 1) / 2;
                let mut polynomials = (0.0_f64, 0.0_f64, 0.0_f64);
                for i in 0..self.order + 1 {
                    let mut li = 1.0_f64;
                    let (t_i, (x_i, y_i, z_i)) = self.buffer[offset + i];
                    for j in 0..self.order + 1 {
                        let (t_j, _) = self.buffer[offset + j];
                        if j != i {
                            li *= (t_s - t_j).to_seconds();
                            li /= (t_i - t_j).to_seconds();
                        }
                    }
                    polynomials.0 += x_i * li;
                    polynomials.1 += y_i * li;
                    polynomials.2 += z_i * li;
                }
                Some(polynomials)
            } else {
                None
            }
        }
    }
    fn new(size: usize) -> Self {
        Self {
            order: size,
            buffer: Vec::<(Epoch, (f64, f64, f64))>::with_capacity(size),
        }
    }
    #[cfg(test)]
    fn snapshot(&self) -> &[(Epoch, (f64, f64, f64))] {
        &self.buffer
    }
}

#[cfg(test)]
mod test {
    use super::PositionInterpolator;
    use crate::interp::Interpolator;
    use hifitime::Epoch;
    use std::str::FromStr;
    // /*
    //  * Lagrangian interpolator theoretical limitations
    //  */
    // fn max_error(values: &Vec(Epoch, f64)>, epoch: Epoch, order: usize) -> f64 {
    //     let mut q = 1.0_f64;
    //     for (e, _) in values {
    //         q *= (epoch - e).to_seconds();
    //     }
    //     let factorial: usize = (1..=order + 1).product();
    //     q.abs() / factorial as f64 // TODO f^(n+1)[x]
    // }
    #[test]
    fn advanced() {
        for (order, max_err) in [(7, 1E-1_f64), (9, 1.0E-2_f64), (11, 0.5E-3_f64)] {
            let mut interp = PositionInterpolator::new(order);
            for (t_k, x_k, y_k, z_k) in [
                (
                    "2020-06-25T00:00:00 UTC",
                    -11562.163582,
                    14053.114306,
                    23345.128269,
                ),
                (
                    "2020-06-25T00:15:00 UTC",
                    -13618.625154,
                    13865.251337,
                    22325.739925,
                ),
                (
                    "2020-06-25T00:30:00 UTC",
                    -15578.906571,
                    13828.422191,
                    21028.690065,
                ),
                (
                    "2020-06-25T00:45:00 UTC",
                    -17408.137167,
                    13927.983030,
                    19470.096085,
                ),
                (
                    "2020-06-25T01:00:00 UTC",
                    -19074.795786,
                    14143.814957,
                    17669.329791,
                ),
                (
                    "2020-06-25T01:15:00 UTC",
                    -20551.644094,
                    14450.989573,
                    15648.777437,
                ),
                (
                    "2020-06-25T01:30:00 UTC",
                    -21816.527671,
                    14820.587810,
                    13433.562098,
                ),
                (
                    "2020-06-25T01:45:00 UTC",
                    -22853.019836,
                    15220.646336,
                    11051.231754,
                ),
                (
                    "2020-06-25T02:00:00 UTC",
                    -23650.888045,
                    15617.201843,
                    8531.416907,
                ),
                (
                    "2020-06-25T02:15:00 UTC",
                    -24206.368268,
                    15975.400418,
                    5905.461978,
                ),
                (
                    "2020-06-25T02:30:00 UTC",
                    -24522.238746,
                    16260.637103,
                    3206.035064,
                ),
                (
                    "2020-06-25T02:45:00 UTC",
                    -24607.690837,
                    16439.689804,
                    466.720941,
                ),
                (
                    "2020-06-25T03:00:00 UTC",
                    -24478.001028,
                    16481.811839,
                    -2278.397553,
                ),
                (
                    "2020-06-25T03:15:00 UTC",
                    -24154.014466,
                    16359.748711,
                    -4995.164479,
                ),
                (
                    "2020-06-25T03:30:00 UTC",
                    -23661.456259,
                    16050.647064,
                    -7649.776712,
                ),
                (
                    "2020-06-25T03:45:00 UTC",
                    -23030.092267,
                    15536.827161,
                    -10209.205433,
                ),
                (
                    "2020-06-25T04:00:00 UTC",
                    -22292.765788,
                    14806.394539,
                    -12641.607839,
                ),
                (
                    "2020-06-25T04:15:00 UTC",
                    -21484.340483,
                    13853.671544,
                    -14916.723906,
                ),
            ] {
                let t_k = Epoch::from_str(t_k).unwrap();
                interp.fill(t_k, (x_k, y_k, z_k));
            }
            for (t_s, expected) in [
                (
                    "2020-06-25T00:00:00 UTC",
                    Some((-11562.163582, 14053.114306, 23345.128269)),
                ),
                ("2020-06-25T00:00:01 UTC", None),
                ("2020-06-25T00:00:10 UTC", None),
                (
                    "2020-06-25T00:15:00 UTC",
                    Some((-13618.625154, 13865.251337, 22325.739925)),
                ),
                (
                    "2020-06-25T01:59:59 UTC",
                    Some((-23650.135944, 15616.7760377, 8534.2815443)),
                ),
            ] {
                let t_s = Epoch::from_str(t_s).unwrap();
                if let Some((x, y, z)) = expected {
                    let (x_k, y_k, z_k) = interp
                        .interpolate(t_s)
                        .expect(&format!("interpolation should be feasible @{:?}", t_s));
                    let err = ((x_k - x).abs(), (y_k - y).abs(), (z_k - z).abs());
                    assert!(err.0 < max_err, "x(err) too large @{:?}", t_s);
                    assert!(err.1 < max_err, "z(err) too large @{:?}", t_s);
                    assert!(err.2 < max_err, "y(err) too large @{:?}", t_s);
                } else {
                    assert!(
                        interp.interpolate(t_s).is_none(),
                        "unexpected interpolation results @{:?}",
                        t_s
                    );
                }
            }
        }
    }
}
