use hifitime::{Duration, Epoch};
use log::debug;

trait Interpolator<T> {
    fn new(order: usize) -> Self;
    fn len(&self) -> usize;
    fn get(&self, idx: usize) -> Option<&(Epoch, T)>;
    fn last(&self) -> Option<&(Epoch, T)>;
    fn clear(&mut self);
    fn push(&mut self, x_j: (Epoch, T));
    fn interpolate(&self, x_s: Epoch) -> Option<T>;
    fn dt(&self) -> Option<(Epoch, Duration)> {
        if self.len() > 1 {
            let (z2, _) = self.get(0)?;
            let (z1, _) = self.last()?;
            Some((*z1, *z1 - *z2))
        } else {
            None
        }
    }
    /// Fill in buffer, which will always be optimized in size and
    /// evenly spaced (in time). Panic on chornological order mix up.
    fn fill(&mut self, x_j: Epoch, y_j: T) {
        if let Some((last, dt)) = self.dt() {
            if (x_j - last).to_seconds().is_sign_positive() {
                if (x_j - last) > dt {
                    debug!("buffer reset on data gap @{:?}:{}", last, x_j - last);
                    self.clear();
                }
                self.push((x_j, y_j));
            } else {
                panic!("samples should be provided in chronological order");
            }
        } else {
            self.push((x_j, y_j));
        }
    }
    #[cfg(test)]
    fn snapshot(&self) -> &[(Epoch, T)];
}

/// Efficient Time Interpolator
pub struct TimeInterpolator {
    buffer: Vec<(Epoch, f64)>,
}

impl Interpolator<f64> for TimeInterpolator {
    fn len(&self) -> usize {
        self.buffer.len()
    }
    fn get(&self, idx: usize) -> Option<&(Epoch, f64)> {
        self.buffer.get(idx)
    }
    fn last(&self) -> Option<&(Epoch, f64)> {
        self.buffer.last()
    }
    fn clear(&mut self) {
        self.buffer.clear();
    }
    fn push(&mut self, x_j: (Epoch, f64)) {
        self.buffer.push(x_j);
    }
    fn interpolate(&self, x_s: Epoch) -> Option<f64> {
        if let Some(y_s) = self
            .buffer
            .iter()
            .filter_map(|(x, y)| if *x == x_s { Some(*y) } else { None })
            .reduce(|k, _| k)
        {
            Some(y_s)
        } else {
            let before = self.buffer.iter().filter(|(x, _)| *x <= x_s).last()?;
            let after = self
                .buffer
                .iter()
                .filter(|(x, _)| *x > x_s)
                .reduce(|k, _| k)?;
            let (before_x, before_y) = before;
            let (after_x, after_y) = after;
            let dx = (*after_x - *before_x).to_seconds();
            let mut dy = (*after_x - x_s).to_seconds() / dx * before_y;
            dy += (x_s - *before_x).to_seconds() / dx * after_y;
            Some(dy)
        }
    }
    fn new(size: usize) -> Self {
        Self {
            buffer: Vec::<(Epoch, f64)>::with_capacity(size),
        }
    }
    #[cfg(test)]
    fn snapshot(&self) -> &[(Epoch, f64)] {
        &self.buffer
    }
}

/// Efficient Position Interpolator
pub struct PositionInterpolator {
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
            let before = self.buffer.iter().filter(|(t, _)| *t <= t_s).last()?;
            let after = self
                .buffer
                .iter()
                .filter(|(t, _)| *t > t_s)
                .reduce(|k, _| k)?;
            let (before_t, (before_x, before_y, before_z)) = before;
            let (after_t, (after_x, after_y, after_z)) = after;
            let dt = (*after_t - *before_t).to_seconds();
            let mut dx = (*after_t - t_s).to_seconds() / dt * before_x;
            let mut dy = (*after_t - t_s).to_seconds() / dt * before_y;
            let mut dz = (*after_t - t_s).to_seconds() / dt * before_z;
            dx += (t_s - *before_t).to_seconds() / dt * after_x;
            dy += (t_s - *before_t).to_seconds() / dt * after_y;
            dz += (t_s - *before_t).to_seconds() / dt * after_z;
            Some((dx, dy, dz))
        }
    }
    fn new(size: usize) -> Self {
        Self {
            buffer: Vec::<(Epoch, (f64, f64, f64))>::with_capacity(size),
        }
    }
    #[cfg(test)]
    fn snapshot(&self) -> &[(Epoch, (f64, f64, f64))] {
        &self.buffer
    }
}

/// Efficient Position Interpolator
#[cfg(test)]
mod test {
    use super::{Interpolator, TimeInterpolator};
    use hifitime::Epoch;
    use std::str::FromStr;
    #[test]
    fn buffer_fill_in() {
        let mut interp = TimeInterpolator::new(4);
        for (x_k, y_k, snapshot) in [
            (
                "2020-01-01T00:00:00 UTC",
                1.0_f64,
                vec![("2020-01-01T00:00:00 UTC", 1.0_f64)],
            ),
            (
                "2020-01-01T00:00:30 UTC",
                2.0_f64,
                vec![
                    ("2020-01-01T00:00:00 UTC", 1.0_f64),
                    ("2020-01-01T00:00:30 UTC", 2.0_f64),
                ],
            ),
            (
                "2020-01-01T00:10:30 UTC",
                3.0_f64,
                vec![("2020-01-01T00:10:30 UTC", 3.0_f64)],
            ),
            (
                "2020-01-01T00:10:35 UTC",
                4.0_f64,
                vec![
                    ("2020-01-01T00:10:30 UTC", 3.0_f64),
                    ("2020-01-01T00:10:35 UTC", 4.0_f64),
                ],
            ),
            (
                "2020-01-01T00:10:40 UTC",
                5.0_f64,
                vec![
                    ("2020-01-01T00:10:30 UTC", 3.0_f64),
                    ("2020-01-01T00:10:35 UTC", 4.0_f64),
                    ("2020-01-01T00:10:40 UTC", 5.0_f64),
                ],
            ),
            (
                "2020-01-01T00:12:00 UTC",
                6.0_f64,
                vec![("2020-01-01T00:12:00 UTC", 6.0_f64)],
            ),
        ] {
            let x_k = Epoch::from_str(x_k).unwrap();
            interp.fill(x_k, y_k);
            let snapshot = snapshot
                .iter()
                .map(|(t, v)| (Epoch::from_str(t).unwrap(), *v))
                .collect::<Vec<_>>();
            assert_eq!(interp.snapshot(), snapshot, "bad buffer logic");
        }
    }
    #[test]
    #[should_panic]
    fn chronological_mix_up() {
        let mut interp = TimeInterpolator::new(4);
        for (x_k, y_k) in [
            ("2020-01-01T00:00:00 UTC", 1.0_f64),
            ("2020-01-01T00:00:30 UTC", 2.0_f64),
            ("2020-01-01T00:00:00 UTC", 3.0_f64),
        ] {
            let x_k = Epoch::from_str(x_k).unwrap();
            interp.fill(x_k, y_k);
        }
    }
    #[test]
    fn time_basic() {
        let mut interp = TimeInterpolator::new(4);
        for (x_k, y_k, x_s, expected) in [
            (
                "2020-01-01T00:00:00 UTC",
                1.0_f64,
                "2020-01-01T00:00:01 UTC",
                None,
            ),
            (
                "2020-01-01T00:00:30 UTC",
                2.0_f64,
                "2020-01-01T00:00:15 UTC",
                Some(1.5),
            ),
            (
                "2020-01-01T00:01:00 UTC",
                3.0_f64,
                "2020-01-01T00:00:15 UTC",
                Some(1.5),
            ),
            (
                "2020-01-01T00:01:00 UTC",
                3.0_f64,
                "2020-01-01T00:01:00 UTC",
                Some(3.0),
            ),
        ] {
            let x_k = Epoch::from_str(x_k).unwrap();
            interp.fill(x_k, y_k);
            let x_s = Epoch::from_str(x_s).unwrap();
            assert_eq!(
                interp.interpolate(x_s),
                expected,
                "wrong interpolation results"
            );
        }
    }
    #[test]
    fn time_advanced() {
        let mut interp = TimeInterpolator::new(3);
        for (x_k, y_k) in [
            ("2019-01-08T00:00:00 UTC", 0.391711350090E-04),
            ("2019-01-08T00:00:30 UTC", 0.391710317949E-04),
            ("2019-01-08T00:01:00 UTC", 0.391708385767E-04),
            ("2019-01-08T00:01:30 UTC", 0.391709678221E-04),
            ("2019-01-08T00:02:00 UTC", 0.391708653726E-04),
            ("2019-01-08T00:02:00 UTC", 0.391708653726E-04),
            ("2019-01-08T00:02:30 UTC", 0.391709273510E-04),
            ("2019-01-08T00:03:00 UTC", 0.391708515569E-04),
            ("2019-01-08T00:03:30 UTC", 0.391706625209E-04),
        ] {
            let x_k = Epoch::from_str(x_k).unwrap();
            interp.fill(x_k, y_k);
        }
        for (x_s, y_s) in [
            ("2019-01-08T00:03:30 UTC", 0.391706625209E-04),
            (
                "2019-01-08T00:01:33 UTC",
                27.0 / 30.0 * 0.391709678221E-04 + 3.0 / 30.0 * 0.391708653726E-04,
            ),
            (
                "2019-01-08T00:01:44 UTC",
                16.0 / 30.0 * 0.391709678221E-04 + 14.0 / 30.0 * 0.391708653726E-04,
            ),
            (
                "2019-01-08T00:01:57 UTC",
                3.0 / 30.0 * 0.391709678221E-04 + 27.0 / 30.0 * 0.391708653726E-04,
            ),
        ] {
            let x_s = Epoch::from_str(x_s).unwrap();
            let y = interp.interpolate(x_s).expect(&format!(
                "interpolation should have been feasible @{:?}",
                x_s
            ));
            assert_eq!(y, y_s, "wrong interpolation results @{:?}", x_s);
        }
    }
}
