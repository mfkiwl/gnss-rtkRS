pub mod position;
pub mod time;

pub use position::PositionInterpolator;
pub use time::TimeInterpolator;

use hifitime::{Duration, Epoch};
use log::debug;

pub(crate) trait Interpolator<T> {
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
