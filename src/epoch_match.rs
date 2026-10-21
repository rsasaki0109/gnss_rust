//! One-to-one matching of epoch streams without interpolation or silent drops.
use crate::GnssTime;
use crate::observation::ObservationEpoch;
use std::fmt;

pub const MATCH_TOLERANCE_S: f64 = 1e-7;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EpochSide {
    Rover,
    Base,
}

#[derive(Debug, Clone, PartialEq)]
pub enum EpochMatch {
    Pair {
        rover: ObservationEpoch,
        base: ObservationEpoch,
    },
    /// Includes complete observations and event flags. A caller that decimates
    /// must account for loss-of-lock/slip evidence in these epochs.
    Unmatched {
        side: EpochSide,
        epoch: ObservationEpoch,
    },
}

#[derive(Debug)]
pub enum MatchError<E> {
    Input {
        side: EpochSide,
        source: E,
    },
    OutOfOrder {
        side: EpochSide,
        previous: GnssTime,
        next: GnssTime,
    },
}
impl<E: fmt::Display> fmt::Display for MatchError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Input { side, source } => write!(f, "{side:?} epoch input: {source}"),
            Self::OutOfOrder {
                side,
                previous,
                next,
            } => write!(
                f,
                "{side:?} epoch time must increase: {previous:?} then {next:?}"
            ),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for MatchError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Input { source, .. } => Some(source),
            Self::OutOfOrder { .. } => None,
        }
    }
}

struct Stream<I> {
    input: I,
    pending: Option<ObservationEpoch>,
    previous: Option<GnssTime>,
    ended: bool,
}
impl<I> Stream<I> {
    fn new(input: I) -> Self {
        Self {
            input,
            pending: None,
            previous: None,
            ended: false,
        }
    }
    fn fill<E>(&mut self, side: EpochSide) -> Result<(), MatchError<E>>
    where
        I: Iterator<Item = Result<ObservationEpoch, E>>,
    {
        if self.pending.is_some() || self.ended {
            return Ok(());
        }
        match self.input.next() {
            Some(Ok(epoch)) => {
                if let Some(previous) = self.previous {
                    if epoch.time <= previous {
                        return Err(MatchError::OutOfOrder {
                            side,
                            previous,
                            next: epoch.time,
                        });
                    }
                }
                self.previous = Some(epoch.time);
                self.pending = Some(epoch);
            }
            Some(Err(source)) => return Err(MatchError::Input { side, source }),
            None => self.ended = true,
        }
        Ok(())
    }
}

/// Streaming merge with one lookahead epoch per input. Emits every measurement
/// epoch once, in a pair within 1e-7 s or as an unmatched epoch. Both input tails
/// are drained. Duplicate/backward times and parse errors terminate the iterator
/// with an error; the iterator is fused thereafter.
pub struct EpochMatcher<R, B> {
    rover: Stream<R>,
    base: Stream<B>,
    done: bool,
}
impl<R, B> EpochMatcher<R, B> {
    pub fn new(rover: R, base: B) -> Self {
        Self {
            rover: Stream::new(rover),
            base: Stream::new(base),
            done: false,
        }
    }
}
impl<R, B, E> Iterator for EpochMatcher<R, B>
where
    R: Iterator<Item = Result<ObservationEpoch, E>>,
    B: Iterator<Item = Result<ObservationEpoch, E>>,
{
    type Item = Result<EpochMatch, MatchError<E>>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        if let Err(error) = self
            .rover
            .fill(EpochSide::Rover)
            .and_then(|()| self.base.fill(EpochSide::Base))
        {
            self.done = true;
            return Some(Err(error));
        }
        let side = match (&self.rover.pending, &self.base.pending) {
            (Some(rover), Some(base))
                if rover.time.difference_seconds(base.time).abs() <= MATCH_TOLERANCE_S =>
            {
                return Some(Ok(EpochMatch::Pair {
                    rover: self.rover.pending.take().unwrap(),
                    base: self.base.pending.take().unwrap(),
                }));
            }
            (Some(rover), Some(base)) if rover.time < base.time => EpochSide::Rover,
            (Some(_), None) => EpochSide::Rover,
            (Some(_), Some(_)) | (None, Some(_)) => EpochSide::Base,
            (None, None) => {
                self.done = true;
                return None;
            }
        };
        let epoch = match side {
            EpochSide::Rover => self.rover.pending.take().unwrap(),
            EpochSide::Base => self.base.pending.take().unwrap(),
        };
        Some(Ok(EpochMatch::Unmatched { side, epoch }))
    }
}
impl<R, B, E> std::iter::FusedIterator for EpochMatcher<R, B>
where
    R: Iterator<Item = Result<ObservationEpoch, E>>,
    B: Iterator<Item = Result<ObservationEpoch, E>>,
{
}
