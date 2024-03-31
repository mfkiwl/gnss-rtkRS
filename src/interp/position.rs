use crate::interp::Interpolator;
use hifitime::{Duration, Epoch};
use log::debug;

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

#[cfg(test)]
mod test {
    use crate::interp::{Interpolator, PositionInterpolator};
    use hifitime::Epoch;
    use std::str::FromStr;
    #[test]
    fn basic() {
        let mut interp = PositionInterpolator::new(4);
        for (x_k, y_k, x_s, expected) in [
            (
                "2020-01-01T00:00:00 UTC",
                (1.1_f64, 1.2_f64, 1.3_f64),
                "2020-01-01T00:00:01 UTC",
                None,
            ),
            (
                "2020-01-01T00:00:30 UTC",
                (2.1_f64, 2.2_f64, 2.3_f64),
                "2020-01-01T00:00:15 UTC",
                Some((1.6_f64, 1.7_f64, 1.8_f64)),
            ),
        ] {
            let x_k = Epoch::from_str(x_k).unwrap();
            interp.fill(x_k, y_k);
            let x_s = Epoch::from_str(x_s).unwrap();

            if let Some(expected) = expected {
                let (dx, dy, dz) = interp
                    .interpolate(x_s)
                    .expect(&format!("interpolation should have be feasible @{:?}", x_s));

                let err = (
                    (dx - expected.0).abs(),
                    (dy - expected.1).abs(),
                    (dz - expected.2).abs(),
                );

                assert!(err.0 < 1.0E-6, "x(err) too large {}@{:?}", err.0, x_s);
                assert!(err.1 < 1.0E-6, "y(err) too large {}@{:?}", err.1, x_s);
                assert!(err.2 < 1.0E-6, "z(err) too large {}@{:?}", err.2, x_s);
            } else {
                assert!(
                    interp.interpolate(x_s).is_none(),
                    "unexpected interp results @{:?}",
                    x_s
                );
            }
        }
    }
    //#[test]
    //fn advanced() {
    //    let mut interp = PositionInterpolator::new(3);
    //    for (x_k, y_k) in [
    //        ("2019-01-08T00:00:00 UTC", 0.391711350090E-04),
    //        ("2019-01-08T00:00:30 UTC", 0.391710317949E-04),
    //        ("2019-01-08T00:01:00 UTC", 0.391708385767E-04),
    //        ("2019-01-08T00:01:30 UTC", 0.391709678221E-04),
    //        ("2019-01-08T00:02:00 UTC", 0.391708653726E-04),
    //        ("2019-01-08T00:02:00 UTC", 0.391708653726E-04),
    //        ("2019-01-08T00:02:30 UTC", 0.391709273510E-04),
    //        ("2019-01-08T00:03:00 UTC", 0.391708515569E-04),
    //        ("2019-01-08T00:03:30 UTC", 0.391706625209E-04),
    //    ] {
    //        let x_k = Epoch::from_str(x_k).unwrap();
    //        interp.fill(x_k, y_k);
    //    }
    //    for (x_s, y_s) in [
    //        ("2019-01-08T00:03:30 UTC", 0.391706625209E-04),
    //        (
    //            "2019-01-08T00:01:33 UTC",
    //            27.0 / 30.0 * 0.391709678221E-04 + 3.0 / 30.0 * 0.391708653726E-04,
    //        ),
    //        (
    //            "2019-01-08T00:01:44 UTC",
    //            16.0 / 30.0 * 0.391709678221E-04 + 14.0 / 30.0 * 0.391708653726E-04,
    //        ),
    //        (
    //            "2019-01-08T00:01:57 UTC",
    //            3.0 / 30.0 * 0.391709678221E-04 + 27.0 / 30.0 * 0.391708653726E-04,
    //        ),
    //    ] {
    //        let x_s = Epoch::from_str(x_s).unwrap();
    //        let y = interp.interpolate(x_s).expect(&format!(
    //            "interpolation should have been feasible @{:?}",
    //            x_s
    //        ));
    //        assert_eq!(y, y_s, "wrong interpolation results @{:?}", x_s);
    //    }
    //}
}
