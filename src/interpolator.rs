use hifitime::{Duration, Epoch};

use log::debug;

#[derive(Debug, Copy, Clone, Default)]
pub enum Method {
    #[default]
    Linear,
    // Lagrange,
}

/// An optimized Interpolator
pub struct Interpolator {
    method: Method,
    buffer: Vec<(Epoch, f64)>,
}

impl Interpolator {
    pub fn new(method: Method, size: usize) -> Self {
        Self {
            method,
            buffer: Vec::<(Epoch, f64)>::with_capacity(size),
        }
    }
    fn linear_interp(x_s: Epoch, buffer: &Vec<(Epoch, f64)>) -> Option<f64> {
        let before = buffer.iter().filter(|(x, _)| *x <= x_s).last()?;
        let after = buffer.iter().filter(|(x, _)| *x > x_s).reduce(|k, _| k)?;
        let (before_x, before_y) = before;
        let (after_x, after_y) = after;
        let dx = (*after_x - *before_x).to_seconds();
        let mut dy = (*after_x - x_s).to_seconds() / dx * before_y;
        dy += (x_s - *before_x).to_seconds() / dx * after_y;
        Some(dy)
    }
    fn dt(&self, x_j: Epoch) -> Option<(Epoch, Duration)> {
        if self.buffer.len() > 1 {
            let (z2, _) = self.buffer.get(self.buffer.len() - 2)?;
            let (z1, _) = self.buffer.last()?;
            Some((*z1, *z1 - *z2))
        } else {
            None
        }
    }
    /// Fill in buffer, which will always be optimized in size and
    /// evenly spaced (in time). Panic on chornological order mix up.
    pub fn fill(&mut self, x_j: Epoch, y_j: f64) {
        if let Some((last, dt)) = self.dt(x_j) {
            if (x_j - last).to_seconds().is_positive() {
                if (x_j - last) > dt {
                    debug!("buffer reset on data gap");
                    self.buffer.clear();
                }
                self.buffer.push((x_j, y_j));
            } else {
                panic!("samples should be provided in chronological order");
            }
        } else {
            self.buffer.push((x_j, y_j));
        }
    }
    pub fn interpolate(&self, x_s: Epoch) -> Option<f64> {
        match self.method {
            Method::Linear => Self::linear_interp(x_s, &self.buffer),
            //Method::Lagrange => lagrange_interp(self.x_k, self.y_k),
        }
    }
    #[cfg(test)]
    pub fn snapshot(&self) -> &[(Epoch, f64)] {
        &self.buffer
    }
}

#[cfg(test)]
mod test {
    use super::{Interpolator, Method};
    use hifitime::Epoch;
    use std::str::FromStr;
    #[test]
    fn buffer_fill_in() {
        let mut interp = Interpolator::new(Method::Linear, 4);
        for (x_k, y_k, snapshot) in [
            (
                Epoch::from_str("2020-01-01T00:00:00 UTC").unwrap(),
                1.0_f64,
                vec![(Epoch::from_str("2020-01-01T00:00:00 UTC").unwrap(), 1.0_f64)],
            ),
            (
                Epoch::from_str("2020-01-01T00:00:30 UTC").unwrap(),
                2.0_f64,
                vec![
                    (Epoch::from_str("2020-01-01T00:00:00 UTC").unwrap(), 1.0_f64),
                    (Epoch::from_str("2020-01-01T00:00:30 UTC").unwrap(), 2.0_f64),
                ],
            ),
            (
                Epoch::from_str("2020-01-01T00:10:30 UTC").unwrap(),
                3.0_f64,
                vec![(Epoch::from_str("2020-01-01T00:10:30 UTC").unwrap(), 3.0_f64)],
            ),
            (
                Epoch::from_str("2020-01-01T00:10:35 UTC").unwrap(),
                4.0_f64,
                vec![
                    (Epoch::from_str("2020-01-01T00:10:30 UTC").unwrap(), 3.0_f64),
                    (Epoch::from_str("2020-01-01T00:10:35 UTC").unwrap(), 4.0_f64),
                ],
            ),
            (
                Epoch::from_str("2020-01-01T00:10:40 UTC").unwrap(),
                5.0_f64,
                vec![
                    (Epoch::from_str("2020-01-01T00:10:30 UTC").unwrap(), 3.0_f64),
                    (Epoch::from_str("2020-01-01T00:10:35 UTC").unwrap(), 4.0_f64),
                    (Epoch::from_str("2020-01-01T00:10:40 UTC").unwrap(), 5.0_f64),
                ],
            ),
            (
                Epoch::from_str("2020-01-01T00:12:00 UTC").unwrap(),
                6.0_f64,
                vec![(Epoch::from_str("2020-01-01T00:12:00 UTC").unwrap(), 6.0_f64)],
            ),
        ] {
            interp.fill(x_k, y_k);
            assert_eq!(interp.snapshot(), snapshot, "bad buffer logic");
        }
    }
    #[test]
    #[should_panic]
    fn chronological_mix_up() {
        let mut interp = Interpolator::new(Method::Linear, 4);
        for (x_k, y_k) in [
            (Epoch::from_str("2020-01-01T00:00:00 UTC").unwrap(), 1.0_f64),
            (Epoch::from_str("2020-01-01T00:00:30 UTC").unwrap(), 2.0_f64),
            (Epoch::from_str("2020-01-01T00:00:00 UTC").unwrap(), 3.0_f64),
        ] {
            interp.fill(x_k, y_k);
        }
    }
}
