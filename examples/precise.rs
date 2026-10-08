//! Inspect the first midpoint epoch of a SP3 product, with optional CLK clocks.
use gnss_rust::precise::{
    PreciseProducts, precise_clock_relativistic_correction, read_clk, read_sp3,
};
use std::{env, error::Error, fs::File, io::BufReader};

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = env::args_os().skip(1);
    let path = args.next().ok_or("usage: precise SP3 [CLK]")?;
    let orbit = read_sp3(BufReader::new(File::open(path)?))?;
    let clocks = args
        .next()
        .map(|path| Ok::<_, Box<dyn Error>>(read_clk(BufReader::new(File::open(path)?))?))
        .transpose()?;
    if args.next().is_some() {
        return Err("usage: precise SP3 [CLK]".into());
    }
    let epoch = orbit
        .header()
        .start_time
        .checked_add_seconds(orbit.header().interval_s / 2.0)?;
    let products = PreciseProducts::new(orbit, clocks);
    println!(
        "satellite,week,tow,frame,reference,x_m,y_m,z_m,vx_m_s,vy_m_s,vz_m_s,clock_source,clock_bias_s,clock_drift_s_s,relativity_s"
    );
    for &satellite in &products.orbit().header().satellites {
        let state = products.interpolate(satellite, epoch, Default::default())?;
        let (source, bias, drift) = match &state.clock {
            Some(c) => (
                format!("{:?}", c.source),
                format!("{:.15e}", c.bias_s),
                format!("{:.15e}", c.drift_s_per_s),
            ),
            None => {
                eprintln!(
                    "{satellite}: clock unavailable: {}",
                    state.clock_unavailability.unwrap()
                );
                (
                    format!("{:?}", state.clock_source),
                    String::new(),
                    String::new(),
                )
            }
        };
        let p = state.position_m;
        let [vx, vy, vz] = state.velocity_m_per_s;
        let relativity = precise_clock_relativistic_correction(p, [vx, vy, vz])?;
        println!(
            "{satellite},{},{:.7},{},{:?},{:.6},{:.6},{:.6},{vx:.9},{vy:.9},{vz:.9},{source},{bias},{drift},{relativity:.15e}",
            epoch.week(),
            epoch.tow(),
            state.coordinate_system,
            state.reference_point,
            p.x,
            p.y,
            p.z
        );
    }
    Ok(())
}
