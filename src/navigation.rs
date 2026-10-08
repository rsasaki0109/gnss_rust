//! Kepler broadcast navigation for GPS, Galileo, QZSS and BeiDou (including GEO).

use crate::{EcefCoord, Error, GnssSystem, GnssTime, SatelliteId};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct KeplerOrbit {
    pub sqrt_a: f64,
    pub eccentricity: f64,
    pub inclination: f64,
    pub omega0: f64,
    pub argument_of_perigee: f64,
    pub mean_anomaly: f64,
    pub delta_n: f64,
    pub idot: f64,
    pub omega_dot: f64,
    pub cuc: f64,
    pub cus: f64,
    pub crc: f64,
    pub crs: f64,
    pub cic: f64,
    pub cis: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct BroadcastClock {
    pub af0: f64,
    pub af1: f64,
    pub af2: f64,
    /// Group delay in seconds; retained separately from state clocks.
    pub tgd: f64,
    pub secondary_tgd: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BroadcastEphemeris {
    pub satellite: SatelliteId,
    pub toe: GnssTime,
    pub toc: GnssTime,
    /// Raw toe seconds in the broadcasting constellation's time scale.
    pub toes: f64,
    pub orbit: KeplerOrbit,
    pub clock: BroadcastClock,
    pub health: u16,
    pub iode: u16,
    pub iodc: u16,
    pub data_source: u16,
    pub accuracy_m: Option<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SatelliteState {
    pub position: EcefCoord,
    pub velocity_mps: [f64; 3],
    /// Satellite clock in seconds, including relativistic correction, excluding TGD.
    pub clock_bias_s: f64,
    pub clock_drift: f64,
}

impl BroadcastEphemeris {
    pub fn validate(&self) -> Result<(), Error> {
        if !matches!(
            self.satellite.system(),
            GnssSystem::Gps | GnssSystem::Galileo | GnssSystem::Qzss | GnssSystem::BeiDou
        ) {
            return Err(Error::UnsupportedConstellation);
        }
        let o = self.orbit;
        let values = [
            self.toes,
            o.sqrt_a,
            o.eccentricity,
            o.inclination,
            o.omega0,
            o.argument_of_perigee,
            o.mean_anomaly,
            o.delta_n,
            o.idot,
            o.omega_dot,
            o.cuc,
            o.cus,
            o.crc,
            o.crs,
            o.cic,
            o.cis,
            self.clock.af0,
            self.clock.af1,
            self.clock.af2,
            self.clock.tgd,
        ];
        if !values.iter().all(|x| x.is_finite())
            || self.clock.secondary_tgd.is_some_and(|x| !x.is_finite())
            || self.accuracy_m.is_some_and(|x| !x.is_finite())
            || o.sqrt_a <= 0.0
            || !(0.0..1.0).contains(&o.eccentricity)
            || !(0.0..604_800.0).contains(&self.toes)
        {
            return Err(Error::InvalidEphemeris);
        }
        Ok(())
    }

    /// Evaluate the broadcast model. Ephemeris age/health admission belongs to
    /// NavigationData::ephemeris; this primitive can evaluate a retained record.
    pub fn satellite_state(&self, time: GnssTime) -> Result<SatelliteState, Error> {
        self.validate()?;
        let (position, clock_bias_s) = self.position_clock(time)?;
        // Match upstream's symmetric 0.5 s finite difference.
        let (forward, clk_forward) = self.position_clock(time.checked_add_seconds(0.5)?)?;
        let (backward, clk_backward) = self.position_clock(time.checked_add_seconds(-0.5)?)?;
        let velocity_mps = [
            forward.x - backward.x,
            forward.y - backward.y,
            forward.z - backward.z,
        ];
        let clock_drift = clk_forward - clk_backward;
        if !velocity_mps.iter().all(|x| x.is_finite()) || !clock_drift.is_finite() {
            return Err(Error::InvalidEphemeris);
        }
        Ok(SatelliteState {
            position,
            velocity_mps,
            clock_bias_s,
            clock_drift,
        })
    }

    /// Position and satellite clock without computing velocity (for code solvers).
    pub fn position_clock(&self, time: GnssTime) -> Result<(EcefCoord, f64), Error> {
        self.validate()?;
        let o = self.orbit;
        let a = o.sqrt_a * o.sqrt_a;
        let tk = wrap_week(time.difference_seconds(self.toe));
        let tc = wrap_week(time.difference_seconds(self.toc));
        let beidou = self.satellite.system() == GnssSystem::BeiDou;
        let mu = if beidou || self.satellite.system() == GnssSystem::Galileo {
            3.986_004_418e14
        } else {
            3.986_005e14
        };
        let earth_rotation = if beidou {
            7.292_115e-5
        } else {
            crate::constants::OMEGA_E
        };
        let mean = o.mean_anomaly + ((mu / (a * a * a)).sqrt() + o.delta_n) * tk;
        if !mean.is_finite() {
            return Err(Error::InvalidEphemeris);
        }
        let mut eccentric = mean;
        let mut converged = false;
        for _ in 0..30 {
            let previous = eccentric;
            let correction = (eccentric - o.eccentricity * eccentric.sin() - mean)
                / (1.0 - o.eccentricity * eccentric.cos());
            eccentric -= correction;
            if (eccentric - previous).abs() < 1e-14 {
                converged = true;
                break;
            }
        }
        if !converged {
            return Err(Error::NonConvergent);
        }
        let (sin_e, cos_e) = eccentric.sin_cos();
        let anomaly =
            ((1.0 - o.eccentricity * o.eccentricity).sqrt() * sin_e).atan2(cos_e - o.eccentricity);
        let phi = anomaly + o.argument_of_perigee;
        let (sin_2phi, cos_2phi) = (2.0 * phi).sin_cos();
        let u = phi + o.cus * sin_2phi + o.cuc * cos_2phi;
        let r = a * (1.0 - o.eccentricity * cos_e) + o.crs * sin_2phi + o.crc * cos_2phi;
        let i = o.inclination + o.idot * tk + o.cis * sin_2phi + o.cic * cos_2phi;
        // Unlike the legacy zero sentinel, an explicit BDT toe of zero is valid.
        let toes = self.toes;
        let x = r * u.cos();
        let y = r * u.sin();
        let (sin_i, cos_i) = i.sin_cos();
        let geo =
            beidou && (self.satellite.prn() <= 5 || (59..=63).contains(&self.satellite.prn()));
        let omega = o.omega0
            + if geo {
                o.omega_dot * tk
            } else {
                (o.omega_dot - earth_rotation) * tk
            }
            - earth_rotation * toes;
        let (sin_omega, cos_omega) = omega.sin_cos();
        let xg = x * cos_omega - y * cos_i * sin_omega;
        let yg = x * sin_omega + y * cos_i * cos_omega;
        let zg = y * sin_i;
        let position = if geo {
            let (sin_rot, cos_rot) = (earth_rotation * tk).sin_cos();
            let sin5 = -0.087_155_742_747_658_2;
            let cos5 = 0.996_194_698_091_745_6;
            EcefCoord::new(
                xg * cos_rot + yg * sin_rot * cos5 + zg * sin_rot * sin5,
                -xg * sin_rot + yg * cos_rot * cos5 + zg * cos_rot * sin5,
                -yg * sin5 + zg * cos5,
            )
        } else {
            EcefCoord::new(xg, yg, zg)
        };
        let clock = self.clock.af0 + self.clock.af1 * tc + self.clock.af2 * tc * tc
            - 4.442_807_633e-10 * o.eccentricity * o.sqrt_a * sin_e;
        if !position.is_finite() || !clock.is_finite() {
            return Err(Error::InvalidEphemeris);
        }
        Ok((position, clock))
    }
}

fn wrap_week(value: f64) -> f64 {
    let mut wrapped = value % 604_800.0;
    if wrapped > 302_400.0 {
        wrapped -= 604_800.0;
    }
    if wrapped < -302_400.0 {
        wrapped += 604_800.0;
    }
    wrapped
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct NavigationData {
    ephemerides: BTreeMap<SatelliteId, Vec<BroadcastEphemeris>>,
    pub gps_ionosphere: Option<([f64; 4], [f64; 4])>,
}
impl NavigationData {
    pub fn add_ephemeris(&mut self, eph: BroadcastEphemeris) -> Result<(), Error> {
        eph.validate()?;
        let list = self.ephemerides.entry(eph.satellite).or_default();
        list.push(eph);
        list.sort_by(|a, b| a.toe.partial_cmp(&b.toe).unwrap());
        Ok(())
    }
    pub fn len(&self) -> usize {
        self.ephemerides.values().map(Vec::len).sum()
    }
    pub fn is_empty(&self) -> bool {
        self.ephemerides.is_empty()
    }
    pub fn records(&self) -> impl Iterator<Item = &BroadcastEphemeris> {
        self.ephemerides.values().flatten()
    }
    /// Select the nearest healthy ephemeris within the caller's age budget.
    /// Galileo I/NAV/F/NAV SSR matching is not implemented in this selector yet.
    pub fn ephemeris(
        &self,
        satellite: SatelliteId,
        time: GnssTime,
        max_age_s: f64,
    ) -> Result<&BroadcastEphemeris, Error> {
        if !max_age_s.is_finite() || max_age_s < 0.0 {
            return Err(Error::InvalidConfiguration);
        }
        self.ephemerides
            .get(&satellite)
            .into_iter()
            .flatten()
            .filter(|eph| eph.health == 0 && time.difference_seconds(eph.toe).abs() <= max_age_s)
            .min_by(|a, b| {
                time.difference_seconds(a.toe)
                    .abs()
                    .total_cmp(&time.difference_seconds(b.toe).abs())
            })
            .ok_or(Error::MissingEphemeris)
    }
}
