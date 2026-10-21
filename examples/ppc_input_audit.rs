//! Read-only PPC input audit; this does not output RTK or PPP solutions.
use gnss_rust::bias::{BiasKind, read_satellite_biases, read_satellite_biases_with_byte_comments};
use gnss_rust::observation::ObservationEpoch;
use gnss_rust::precise::{PreciseProducts, read_satellite_clocks, read_sp3};
use gnss_rust::rinex::{RinexObservationReader, read_navigation_for_systems};
use gnss_rust::{GnssSystem, GnssTime, SatelliteId};
use std::collections::BTreeSet;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

#[derive(Default)]
struct Counts {
    epochs: usize,
    slip_records: usize,
    first: Option<GnssTime>,
    last: Option<GnssTime>,
}
impl Counts {
    fn next(
        &mut self,
        reader: &mut RinexObservationReader<BufReader<File>>,
    ) -> Result<Option<ObservationEpoch>, Box<dyn Error>> {
        let epoch = reader.next().transpose()?;
        if let Some(epoch) = &epoch {
            if self.last.is_some_and(|last| epoch.time <= last) {
                return Err("input epochs must strictly increase".into());
            }
            self.epochs += 1;
            self.slip_records += usize::from(epoch.flag == 6);
            self.first.get_or_insert(epoch.time);
            self.last = Some(epoch.time);
        }
        Ok(epoch)
    }
}

fn gps_l1(epoch: &ObservationEpoch) -> BTreeSet<SatelliteId> {
    epoch
        .observations
        .iter()
        .filter(|o| {
            o.satellite.system() == GnssSystem::Gps
                && o.tracking_code == "1C"
                && o.pseudorange_m.is_some()
                && o.carrier_phase_cycles.is_some()
                && !o.half_cycle_ambiguity()
        })
        .map(|o| o.satellite)
        .collect()
}

fn observations(root: &Path) -> Result<(), Box<dyn Error>> {
    println!(
        "case,rover_epochs,base_epochs,coincident_pairs,rover_only,base_only,common_gps_l1_ge4_pairs,rover_phase_shift_declarations,base_phase_shift_declarations,gps_nav_records,excluded_nav_records,rover_slip_records,base_slip_records,week,first_tow,last_tow"
    );
    for city in ["nagoya", "tokyo"] {
        for run in ["run1", "run2", "run3"] {
            let directory = root.join(city).join(run);
            let open = |name| File::open(directory.join(name)).map(BufReader::new);
            let mut rover = RinexObservationReader::new(open("rover.obs")?)?;
            let mut base = RinexObservationReader::new(open("base.obs")?)?;
            let phase_counts = (
                rover.header().phase_shifts.len(),
                base.header().phase_shifts.len(),
            );
            let nav = read_navigation_for_systems(open("base.nav")?, &[GnssSystem::Gps])?;
            let (mut ro_counts, mut ba_counts) = (Counts::default(), Counts::default());
            let mut ro = ro_counts.next(&mut rover)?;
            let mut ba = ba_counts.next(&mut base)?;
            let (mut pairs, mut rover_only, mut base_only, mut sufficient) = (0, 0, 0, 0);
            loop {
                match (&ro, &ba) {
                    (Some(r), Some(b)) if r.time.difference_seconds(b.time).abs() <= 1e-7 => {
                        pairs += 1;
                        sufficient += usize::from(
                            matches!(r.flag, 0 | 1)
                                && matches!(b.flag, 0 | 1)
                                && gps_l1(r).intersection(&gps_l1(b)).count() >= 4,
                        );
                        ro = ro_counts.next(&mut rover)?;
                        ba = ba_counts.next(&mut base)?;
                    }
                    (Some(r), Some(b)) if r.time < b.time => {
                        rover_only += 1;
                        ro = ro_counts.next(&mut rover)?;
                    }
                    (Some(_), Some(_)) | (None, Some(_)) => {
                        base_only += 1;
                        ba = ba_counts.next(&mut base)?;
                    }
                    (Some(_), None) => {
                        rover_only += 1;
                        ro = ro_counts.next(&mut rover)?;
                    }
                    (None, None) => break,
                }
            }
            assert_eq!(pairs + rover_only, ro_counts.epochs);
            assert_eq!(pairs + base_only, ba_counts.epochs);
            let first = ro_counts.first.ok_or("empty rover input")?;
            println!(
                "{city}-{run},{},{},{pairs},{rover_only},{base_only},{sufficient},{},{},{},{},{},{},{},{:.7},{:.7}",
                ro_counts.epochs,
                ba_counts.epochs,
                phase_counts.0,
                phase_counts.1,
                nav.navigation.len(),
                nav.excluded_records.len(),
                ro_counts.slip_records,
                ba_counts.slip_records,
                first.week(),
                first.tow(),
                ro_counts.last.ok_or("empty rover input")?.tow()
            );
        }
    }
    Ok(())
}

fn products(root: &Path) -> Result<(), Box<dyn Error>> {
    println!(
        "date,orbit_satellites,orbit_samples,clock_satellites,clock_samples,bias_records,dsb_records,osb_records,receiver_clock_records,excluded_receiver_bias_rows,midpoint_states_with_clock,midpoint_unavailable"
    );
    let mut paths = std::fs::read_dir(root)?
        .map(|entry| entry.map(|e| e.path()))
        .collect::<Result<Vec<_>, _>>()?;
    paths.sort();
    for path in paths
        .iter()
        .filter(|path| path.extension().is_some_and(|e| e == "SP3"))
    {
        let name = path
            .file_name()
            .ok_or("missing filename")?
            .to_str()
            .ok_or("non-UTF8 filename")?;
        let stamp = name.split('_').nth(1).ok_or("missing product date")?;
        let clock = root.join(format!("IGS0OPSFIN_{stamp}_01D_30S_CLK.CLK"));
        let bias = root.join(format!("CAS0OPSRAP_{stamp}_01D_01D_DCB.BIA"));
        let orbit = read_sp3(BufReader::new(File::open(path)?))?;
        let selected_clock = read_satellite_clocks(BufReader::new(File::open(clock)?))?;
        let receiver_records = selected_clock.receiver_records.len();
        let clock = selected_clock.clocks;
        let selected_bias =
            read_satellite_biases(BufReader::new(File::open(bias)?), Default::default())?;
        let receiver_bias_rows = selected_bias.excluded_receiver_records.len();
        let bias = selected_bias.products;
        let midpoint = orbit.header().start_time.checked_add_seconds(450.0)?;
        let counts = (
            orbit.samples().len(),
            orbit.samples().values().map(Vec::len).sum::<usize>(),
            clock.samples().len(),
            clock.samples().values().map(Vec::len).sum::<usize>(),
            bias.entries().len(),
            bias.entries()
                .iter()
                .filter(|b| b.kind == BiasKind::Dsb)
                .count(),
            bias.entries()
                .iter()
                .filter(|b| b.kind == BiasKind::Osb)
                .count(),
        );
        let products = PreciseProducts::new(orbit, Some(clock));
        let (mut available, mut unavailable) = (0, 0);
        for &satellite in &products.orbit().header().satellites {
            match products.interpolate(satellite, midpoint, Default::default()) {
                Ok(state) if state.clock.is_some() => available += 1,
                _ => unavailable += 1,
            }
        }
        println!(
            "{stamp},{},{},{},{},{},{},{},{receiver_records},{receiver_bias_rows},{available},{unavailable}",
            counts.0, counts.1, counts.2, counts.3, counts.4, counts.5, counts.6
        );
    }
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    if args.len() != 2 {
        return Err(
            "usage: ppc_input_audit observations DATASET_ROOT | products PRODUCT_ROOT | code-products PRODUCT_ROOT".into(),
        );
    }
    if args[0] == "observations" {
        observations(Path::new(&args[1]))?;
    } else if args[0] == "products" {
        products(Path::new(&args[1]))?;
    } else if args[0] == "code-products" {
        code_products(Path::new(&args[1]))?;
    } else {
        return Err("expected observations, products or code-products".into());
    }
    eprintln!(
        "readiness_audit_only=true RTK_solution=false PPP_solution=false performance_accepted=false"
    );
    Ok(())
}

fn code_products(root: &Path) -> Result<(), Box<dyn Error>> {
    println!(
        "date,frame,orbit_satellites,orbit_samples,clock_satellites,clock_samples,satellite_bias_rows,excluded_receiver_bias_rows,preserved_byte_comment_lines,GPS_IF_C1W_C2W_declared,bias_APC_IGS20,GPS_sample_states_with_clock,GPS_sample_states_unavailable,code_1c2w_pairs,phase_1c2w_pairs,code_1c2l_pairs,phase_1c2l_pairs,gauge_pairs_missing,max_abs_GPS_IF_code_gauge_m"
    );
    let mut paths = std::fs::read_dir(root)?
        .map(|e| e.map(|e| e.path()))
        .collect::<Result<Vec<_>, _>>()?;
    paths.sort();
    for path in paths
        .iter()
        .filter(|p| p.extension().is_some_and(|e| e == "SP3"))
    {
        let name = path
            .file_name()
            .ok_or("missing product name")?
            .to_str()
            .ok_or("non-UTF8 product name")?;
        let stamp = name.split('_').nth(1).ok_or("missing product date")?;
        let orbit = read_sp3(BufReader::new(File::open(path)?))?;
        let clocks = read_satellite_clocks(BufReader::new(File::open(
            root.join(format!("COD0OPSRAP_{stamp}_01D_30S_CLK.CLK")),
        )?))?
        .clocks;
        let bias = read_satellite_biases_with_byte_comments(
            BufReader::new(File::open(
                root.join(format!("COD0OPSRAP_{stamp}_01D_01D_OSB.BIA")),
            )?),
            Default::default(),
        )?;
        let description = bias.selection.products.description();
        let declared = description
            .get("SATELLITE_CLOCK_REFERENCE_OBSERVABLES")
            .is_some_and(|v| v.iter().any(|v| v.as_slice() == ["G", "C1W", "C2W"]));
        let apc = description
            .get("APC_MODEL")
            .is_some_and(|v| v.len() == 1 && v[0].as_slice() == ["IGS20"]);
        let frame = orbit.header().coordinate_system.clone();
        let sample = orbit.header().start_time.checked_add_seconds(450.0)?;
        let counts = (
            orbit.samples().len(),
            orbit.samples().values().map(Vec::len).sum::<usize>(),
            clocks.samples().len(),
            clocks.samples().values().map(Vec::len).sum::<usize>(),
        );
        let products = PreciseProducts::new(orbit, Some(clocks));
        let table = &bias.selection.products;
        let (mut available, mut unavailable, mut missing, mut max_gauge) = (0, 0, 0, 0.0_f64);
        let mut pairs = [0usize; 4];
        for &sat in products
            .orbit()
            .samples()
            .keys()
            .filter(|s| s.system() == GnssSystem::Gps)
        {
            match products.interpolate(sat, sample, Default::default()) {
                Ok(s) if s.clock.is_some() => available += 1,
                _ => unavailable += 1,
            }
            match (
                table.select_osb(sat, "C1W", sample),
                table.select_osb(sat, "C2W", sample),
            ) {
                (Ok(a), Ok(b)) => {
                    let residual = gnss_rust::constants::IFLC_F1 * a.metres(None)?.0
                        + gnss_rust::constants::IFLC_F2 * b.metres(None)?.0;
                    max_gauge = max_gauge.max(residual.abs());
                }
                _ => missing += 1,
            }
            for (i, second) in ["2W", "2L"].iter().enumerate() {
                if table.select_osb(sat, "C1C", sample).is_ok()
                    && table.select_osb(sat, &format!("C{second}"), sample).is_ok()
                {
                    pairs[2 * i] += 1;
                }
                if table.select_osb(sat, "L1C", sample).is_ok()
                    && table.select_osb(sat, &format!("L{second}"), sample).is_ok()
                {
                    pairs[2 * i + 1] += 1;
                }
            }
        }
        println!(
            "{stamp},{frame},{},{},{},{},{},{},{},{declared},{apc},{available},{unavailable},{},{},{},{},{missing},{max_gauge:.17e}",
            counts.0,
            counts.1,
            counts.2,
            counts.3,
            table.entries().len(),
            bias.selection.excluded_receiver_records.len(),
            bias.preserved_free_text.len(),
            pairs[0],
            pairs[1],
            pairs[2],
            pairs[3]
        );
    }
    Ok(())
}
