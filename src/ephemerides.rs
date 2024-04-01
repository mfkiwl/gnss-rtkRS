use hifitime::Epoch;

/// Ephemerides data source.
pub trait EphemeridesIter {
    /// Provide Time of Issue of Ephemerides, in chronological order,
    /// with the temporary limitation to the GPST timescale.
    fn next(&mut self) -> Option<Epoch>;
}
