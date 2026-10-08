//! Explicit-body IERS 2010 kernel using native SOFA-generated fixture inputs.
//! This example does not generate ephemerides/time arguments in Rust or run PPP.
use gnss_rust::EcefCoord;
use gnss_rust::iers2010::{EarthFixedSunMoon, IersTideArguments, solid_earth_tide};

fn vector(row: &[f64]) -> EcefCoord {
    EcefCoord::new(row[0], row[1], row[2])
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    if std::env::args().len() != 1 {
        return Err("usage: iers_tide (embedded explicit-body fixtures)".into());
    }
    println!("id,mjd_utc,dx_m,dy_m,dz_m,instant_x_m,instant_y_m,instant_z_m");
    let mut count = 0;
    for line in include_str!("../tests/fixtures/upstream_iers2010_epochs.csv")
        .lines()
        .skip(1)
    {
        let row: Vec<f64> = line.split(',').map(str::parse).collect::<Result<_, _>>()?;
        let result = solid_earth_tide(
            vector(&row[7..10]),
            EarthFixedSunMoon::new(vector(&row[10..13]), vector(&row[13..16]))?,
            IersTideArguments::new(row[2], row[3])?,
        )?;
        let d = result.displacement_m;
        let p = result.instantaneous_marker_m;
        println!(
            "{:.0},{:.12},{:.15},{:.15},{:.15},{:.9},{:.9},{:.9}",
            row[0], row[1], d.x, d.y, d.z, p.x, p.y, p.z
        );
        count += 1;
    }
    eprintln!("completed {count} explicit-body IERS 2010 tide cases (component example)");
    Ok(())
}
