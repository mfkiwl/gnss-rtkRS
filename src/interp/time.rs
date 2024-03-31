use gnss::prelude::SV;
use hifitime::{Duration, Epoch};
use log::debug;
use std::collections::HashMap;

use crate::{clock::Clock, interp::Interpolator};

#[derive(Debug)]
pub struct ClockInterpolator {
    size: usize,
    pub interpolators: HashMap<SV, TimeInterpolator>,
}

impl ClockInterpolator {
    pub fn malloc(size: usize) -> Self {
        Self {
            size,
            interpolators: HashMap::<SV, TimeInterpolator>::with_capacity(size),
        }
    }
    pub fn new_clock(&mut self, ck: Clock) {
        if let Some(interp) = self.interpolators.get_mut(&ck.sv) {
            interp.push((ck.epoch, ck.offset));
        } else {
            let mut interp = TimeInterpolator::new(self.size);
            interp.push((ck.epoch, ck.offset));
            self.interpolators.insert(ck.sv, interp);
        }
    }
    pub fn interpolate(&self, sv: SV, t_k: Epoch) -> Option<Clock> {
        let interp = self
            .interpolators
            .iter()
            .filter_map(|(k, v)| if *k == sv { Some(v) } else { None })
            .reduce(|k, _| k)?;
        let offset = interp.interpolate(t_k)?;
        Some(Clock::new(sv, t_k, offset, None, None)) // TODO: drift + drift/r
    }
}

/// Efficient Time Interpolator
#[derive(Debug)]
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

#[cfg(test)]
mod test {
    use crate::interp::{Interpolator, TimeInterpolator};
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
    fn basic() {
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
    fn advanced() {
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
