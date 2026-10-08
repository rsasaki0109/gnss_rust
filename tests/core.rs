use gnss_rust::{Error, GnssSystem, GnssTime, SatelliteId, SignalType, constants as c};
use std::time::{Duration, UNIX_EPOCH};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn physical_constants_and_wavelengths() {
    assert_eq!(c::SPEED_OF_LIGHT, 299_792_458.0);
    near(c::WGS84_B, 6_356_752.314_245_179, 1e-9);
    near(c::GPS_L1_WAVELENGTH, 0.190_293_672_798_364_87, 1e-15);
    near(c::GPS_L2_WAVELENGTH, 0.244_210_213_424_568_25, 1e-15);
    near(c::GPS_WL_WAVELENGTH, 0.861_918_400_322_005_6, 1e-15);
    near(c::IFLC_C1 + c::IFLC_C2, 1.0, 1e-15);
    near(c::IFLC_F1 + c::IFLC_F2, 1.0, 1e-15);
    near(c::IFLC_C1, c::IFLC_F1, 1e-14);
}

#[test]
fn upstream_time_arithmetic_and_ordering() {
    let first = GnssTime::new(2000, 345_600.0).unwrap();
    let second = GnssTime::new(2000, 345_700.0).unwrap();
    assert_eq!(second.difference_seconds(first), 100.0);
    assert!(first < second);
    assert_eq!(first.checked_add_seconds(100.0).unwrap(), second);
    assert_eq!(second.checked_add_seconds(-100.0).unwrap(), first);
}

#[test]
fn time_normalizes_both_directions_across_multiple_weeks() {
    for (week, tow, expected_week, expected_tow) in [
        (2300, 604_800.25, 2301, 0.25),
        (2300, -0.25, 2299, 604_799.75),
        (0, -604_800.0, -1, 0.0),
        (2300, 1_814_400.5, 2303, 0.5),
        (2300, -1_814_400.5, 2296, 604_799.5),
    ] {
        let time = GnssTime::new(week, tow).unwrap();
        assert_eq!((time.week(), time.tow()), (expected_week, expected_tow));
    }
    assert_eq!(
        GnssTime::new(1, 0.25)
            .unwrap()
            .checked_add_seconds(-0.5)
            .unwrap(),
        GnssTime::new(0, 604_799.75).unwrap()
    );
}

#[test]
fn tiny_negative_time_remainder_and_signed_zero_stay_normalized() {
    let tiny = GnssTime::new(100, -1e-15).unwrap();
    assert!((0.0..604_800.0).contains(&tiny.tow()));
    let positive = GnssTime::new(100, 0.0).unwrap();
    let negative = GnssTime::new(100, -0.0).unwrap();
    assert_eq!(positive, negative);
    assert_eq!(
        positive.partial_cmp(&negative),
        Some(std::cmp::Ordering::Equal)
    );
}

#[test]
fn exact_equality_and_explicit_time_tolerance() {
    let first = GnssTime::new(10, 1.0).unwrap();
    let second = GnssTime::new(10, 1.000_000_5).unwrap();
    assert_ne!(first, second);
    assert!(first < second);
    assert!(second.difference_seconds(first).abs() < 1e-6);
}

#[test]
fn time_rejects_nonfinite_and_out_of_range_values() {
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(GnssTime::new(0, value), Err(Error::NonFinite));
        assert_eq!(
            GnssTime::default().checked_add_seconds(value),
            Err(Error::NonFinite)
        );
    }
    assert_eq!(GnssTime::new(0, f64::MAX), Err(Error::OutOfRange));
    assert_eq!(GnssTime::new(i32::MAX, 604_800.0), Err(Error::OutOfRange));
    assert_eq!(GnssTime::new(i32::MIN, -1.0), Err(Error::OutOfRange));
    assert_eq!(
        GnssTime::new(i32::MAX, 1.0)
            .unwrap()
            .checked_add_seconds(604_800.0),
        Err(Error::OutOfRange)
    );
}

#[test]
fn time_difference_avoids_i32_week_subtraction_overflow() {
    let high = GnssTime::new(i32::MAX, 0.0).unwrap();
    let low = GnssTime::new(i32::MIN, 0.0).unwrap();
    assert_eq!(high.difference_seconds(low), 4_294_967_295.0 * 604_800.0);
}

#[test]
fn gps_epoch_label_matches_upstream_epoch_offset() {
    let epoch = UNIX_EPOCH + Duration::from_secs(315_964_800);
    assert_eq!(
        GnssTime::default().to_system_time_gps_scale().unwrap(),
        epoch
    );
    assert_eq!(
        GnssTime::from_system_time_gps_scale(epoch).unwrap(),
        GnssTime::default()
    );
    // This is a GPS-scale label, deliberately not a modern UTC conversion.
    let later = GnssTime::new(2000, 345_600.0).unwrap();
    assert_eq!(
        later.to_system_time_gps_scale().unwrap(),
        epoch + Duration::from_secs(2000 * 604_800 + 345_600)
    );
}

#[test]
fn system_time_preserves_fractional_seconds_before_and_after_epochs() {
    let epoch = UNIX_EPOCH + Duration::from_secs(315_964_800);
    for time in [
        epoch + Duration::new(0, 123_456_789),
        epoch - Duration::new(0, 125_000_000),
        UNIX_EPOCH - Duration::new(1, 250_000_000),
        epoch + Duration::new(2300 * 604_800 + 12345, 125_000_000),
    ] {
        let gps = GnssTime::from_system_time_gps_scale(time).unwrap();
        let round_trip = gps.to_system_time_gps_scale().unwrap();
        let error = round_trip
            .duration_since(time)
            .unwrap_or_else(|e| e.duration());
        assert!(
            error <= Duration::from_nanos(1),
            "round-trip error: {error:?}"
        );
    }
}

#[test]
fn satellite_identity_format_parse_and_upstream_ordering() {
    for (system, text, prn) in [
        (GnssSystem::Gps, "G01", 1),
        (GnssSystem::Glonass, "R24", 24),
        (GnssSystem::Galileo, "E11", 11),
        (GnssSystem::BeiDou, "C06", 6),
        (GnssSystem::Qzss, "J193", 193),
        (GnssSystem::Sbas, "S120", 120),
        (GnssSystem::NavIC, "I07", 7),
    ] {
        let satellite = SatelliteId::new(system, prn).unwrap();
        assert_eq!(satellite.to_string(), text);
        assert_eq!(text.parse::<SatelliteId>().unwrap(), satellite);
        assert_eq!(satellite.system(), system);
        assert_eq!(satellite.prn(), prn);
    }
    assert!("G01".parse::<SatelliteId>().unwrap() < "G02".parse::<SatelliteId>().unwrap());
    assert!("G02".parse::<SatelliteId>().unwrap() < "R01".parse::<SatelliteId>().unwrap());
}

#[test]
fn malformed_satellite_ids_are_rejected_without_panics() {
    for value in [
        "", "G", "G1", "G00", "G256", "G-1", "G+1", "g01", "U01", "X01", "G0011", " G01", "G01 ",
        "é12", "Ｇ01", "G١",
    ] {
        assert!(value.parse::<SatelliteId>().is_err(), "accepted {value}");
    }
    assert!(SatelliteId::new(GnssSystem::Unknown, 1).is_err());
    assert!(SatelliteId::new(GnssSystem::Gps, 0).is_err());
}

#[test]
fn declared_signals_have_correct_constellation_and_frequency() {
    use SignalType::*;
    for (signal, system, frequency) in [
        (GpsL1Ca, GnssSystem::Gps, 1575.42e6),
        (GpsL1P, GnssSystem::Gps, 1575.42e6),
        (GpsL2P, GnssSystem::Gps, 1227.60e6),
        (GpsL2C, GnssSystem::Gps, 1227.60e6),
        (GpsL5, GnssSystem::Gps, 1176.45e6),
        (GalE1, GnssSystem::Galileo, 1575.42e6),
        (GalE5A, GnssSystem::Galileo, 1176.45e6),
        (GalE5B, GnssSystem::Galileo, 1207.14e6),
        (GalE6, GnssSystem::Galileo, 1278.75e6),
        (BdsB1I, GnssSystem::BeiDou, 1561.098e6),
        (BdsB2I, GnssSystem::BeiDou, 1207.14e6),
        (BdsB3I, GnssSystem::BeiDou, 1268.52e6),
        (BdsB1C, GnssSystem::BeiDou, 1575.42e6),
        (BdsB2A, GnssSystem::BeiDou, 1176.45e6),
        (QzsL1Ca, GnssSystem::Qzss, 1575.42e6),
        (QzsL2C, GnssSystem::Qzss, 1227.60e6),
        (QzsL5, GnssSystem::Qzss, 1176.45e6),
    ] {
        assert_eq!(signal.system(), system);
        assert_eq!(signal.frequency_hz(None).unwrap(), frequency);
        near(
            signal.wavelength_m(None).unwrap() * frequency,
            c::SPEED_OF_LIGHT,
            1e-7,
        );
    }
    assert_eq!(GpsL1Ca as u8, 0);
    assert_eq!(QzsL5 as u8, 20);
}

#[test]
fn glonass_fdma_channels_require_explicit_valid_channel() {
    use SignalType::*;
    for signal in [GloL1Ca, GloL1P, GloL2Ca, GloL2P] {
        assert_eq!(signal.system(), GnssSystem::Glonass);
        assert_eq!(signal.frequency_hz(None), Err(Error::MissingGlonassChannel));
        for channel in [-8, 7, i8::MIN, i8::MAX] {
            assert_eq!(
                signal.frequency_hz(Some(channel)),
                Err(Error::InvalidGlonassChannel)
            );
        }
    }
    assert_eq!(GloL1Ca.frequency_hz(Some(-7)).unwrap(), 1_598_062_500.0);
    assert_eq!(GloL2P.frequency_hz(Some(6)).unwrap(), 1_248_625_000.0);
    near(
        GloL1P.wavelength_m(Some(0)).unwrap(),
        c::SPEED_OF_LIGHT / 1_602_000_000.0,
        1e-15,
    );
}
