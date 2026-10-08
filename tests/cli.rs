use std::{fs, process::Command};

fn fixture(name: &str) -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

#[test]
fn cli_runs_the_file_pipeline_and_reports_all_epochs() {
    let result = Command::new(env!("CARGO_BIN_EXE_gnss-rust"))
        .args(["spp", "--obs"])
        .arg(fixture("synthetic_spp.obs"))
        .arg("--nav")
        .arg(fixture("synthetic_spp.nav"))
        .arg("--no-atmosphere")
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let output = String::from_utf8(result.stdout).unwrap();
    let rows: Vec<_> = output.lines().collect();
    assert_eq!(rows.len(), 9);
    assert!(rows[0].starts_with("week,tow,x_m"));
    for row in &rows[1..] {
        assert_eq!(row.split(',').count(), 10);
        assert!(row.contains(",9,SPP,"));
    }
    assert!(
        String::from_utf8_lossy(&result.stderr).contains("processed 8 measurement epochs; valid 8")
    );
}

#[test]
fn cli_default_models_work_for_atmospheric_fixture() {
    let result = Command::new(env!("CARGO_BIN_EXE_gnss-rust"))
        .args(["spp", "--obs"])
        .arg(fixture("synthetic_spp_atmosphere.obs"))
        .arg("--nav")
        .arg(fixture("synthetic_spp.nav"))
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(String::from_utf8(result.stdout).unwrap().lines().count(), 9);
}

#[test]
fn cli_invalid_arguments_and_unimplemented_modes_do_not_claim_success() {
    for args in [
        vec!["rtk"],
        vec!["ppp"],
        vec!["spp"],
        vec!["spp", "--unknown"],
        vec!["spp", "--obs"],
    ] {
        assert!(
            !Command::new(env!("CARGO_BIN_EXE_gnss-rust"))
                .args(args)
                .output()
                .unwrap()
                .status
                .success()
        );
    }
    assert!(
        Command::new(env!("CARGO_BIN_EXE_gnss-rust"))
            .arg("--help")
            .output()
            .unwrap()
            .status
            .success()
    );
}

#[test]
fn cli_does_not_overwrite_existing_output() {
    let directory =
        std::env::temp_dir().join(format!("gnss-rust-cli-output-{}", std::process::id()));
    fs::create_dir(&directory).unwrap();
    let output = directory.join("existing.csv");
    fs::write(&output, "preserve existing user output\n").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_gnss-rust"))
        .args(["spp", "--obs"])
        .arg(fixture("synthetic_spp.obs"))
        .arg("--nav")
        .arg(fixture("synthetic_spp.nav"))
        .arg("--output")
        .arg(&output)
        .arg("--no-atmosphere")
        .output()
        .unwrap();
    let preserved = fs::read_to_string(&output).unwrap();
    fs::remove_dir_all(directory).unwrap();
    assert!(!result.status.success());
    assert_eq!(preserved, "preserve existing user output\n");
}

fn rtk_command() -> Command {
    let line = include_str!("fixtures/synthetic_rtk_truth.csv")
        .lines()
        .find(|l| !l.starts_with('#'))
        .unwrap();
    let base = line.split(',').take(3).collect::<Vec<_>>().join(",");
    let mut command = Command::new(env!("CARGO_BIN_EXE_gnss-rust"));
    command
        .args(["rtk-float", "--rover"])
        .arg(fixture("synthetic_rtk_rover.obs"))
        .arg("--base")
        .arg(fixture("synthetic_rtk_base.obs"))
        .arg("--nav")
        .arg(fixture("synthetic_rtk.nav"))
        .arg("--base-position")
        .arg(base);
    command
}

#[test]
fn cli_static_float_runs_all_pairs_and_outputs_only_float_with_candidate_diagnostics() {
    let result = rtk_command().output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let output = String::from_utf8(result.stdout).unwrap();
    let rows: Vec<_> = output.lines().collect();
    assert_eq!(rows.len(), 33);
    assert!(rows[0].ends_with("candidate_ratio"));
    for row in &rows[1..] {
        let columns: Vec<_> = row.split(',').collect();
        assert_eq!(columns.len(), 14);
        assert_eq!(columns[5], "9");
        assert_eq!(columns[6], "FLOAT");
        assert!(columns[13].parse::<f64>().unwrap().is_finite());
    }
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("processed 32 matched epoch pairs; valid FLOAT 32"));
    assert!(stderr.contains("unavailable ambiguity diagnostics 0"));
}

#[test]
fn cli_float_requires_surveyed_base_coordinates_and_checks_coordinate_options() {
    for args in [
        vec!["rtk-float"],
        vec!["rtk-float", "--base-position", "1,2"],
        vec!["rtk-float", "--base-position", "NaN,2,3"],
        vec![
            "rtk-float",
            "--base-position",
            "1,2,3",
            "--base-position",
            "4,5,6",
        ],
    ] {
        assert!(
            !Command::new(env!("CARGO_BIN_EXE_gnss-rust"))
                .args(args)
                .output()
                .unwrap()
                .status
                .success()
        );
    }
}

#[test]
fn cli_float_rejects_unmatched_epoch_sequences_and_preserves_existing_output() {
    let directory =
        std::env::temp_dir().join(format!("gnss-rust-cli-float-{}", std::process::id()));
    fs::create_dir(&directory).unwrap();
    let base = directory.join("mismatched.obs");
    let text =
        include_str!("fixtures/synthetic_rtk_base.obs").replacen("0.0000000", "0.1000000", 1);
    fs::write(&base, text).unwrap();
    // Build the same explicit arguments with a different base file.
    let mut command = rtk_command();
    let args: Vec<_> = command.get_args().map(|v| v.to_os_string()).collect();
    let mut mismatched = Command::new(command.get_program());
    let mut args = args;
    let index = args.iter().position(|v| v == "--base").unwrap();
    args[index + 1] = base.into_os_string();
    let result = mismatched.args(args).output().unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("unmatched rover/base times"));
    let out = directory.join("existing.csv");
    fs::write(&out, "preserved\n").unwrap();
    let result = command.arg("--output").arg(&out).output().unwrap();
    assert!(!result.status.success());
    assert_eq!(fs::read_to_string(&out).unwrap(), "preserved\n");
    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn cli_rtk_reports_only_validated_fixes_and_reasons_for_float() {
    for noisy in [false, true] {
        let command = rtk_command();
        let mut args: Vec<_> = command.get_args().map(|v| v.to_os_string()).collect();
        args[0] = "rtk".into();
        if noisy {
            let index = args.iter().position(|v| v == "--rover").unwrap();
            args[index + 1] = fixture("synthetic_rtk_rover_noisy.obs").into_os_string();
            let index = args.iter().position(|v| v == "--base").unwrap();
            args[index + 1] = fixture("synthetic_rtk_base_noisy.obs").into_os_string();
        }
        let result = Command::new(command.get_program())
            .args(args)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let output = String::from_utf8(result.stdout).unwrap();
        let rows: Vec<_> = output.lines().collect();
        assert_eq!(rows.len(), 33);
        assert!(rows[0].ends_with("candidate_postfit_nis_per_row"));
        for (i, row) in rows[1..].iter().enumerate() {
            let columns: Vec<_> = row.split(',').collect();
            assert_eq!(columns.len(), 18);
            assert_eq!(columns[6], if i < 5 { "FLOAT" } else { "FIX" });
            assert_eq!(
                columns[14],
                if i < 4 {
                    "insufficient_lock"
                } else if i == 4 {
                    "confirming"
                } else {
                    "accepted"
                }
            );
            if i >= 5 {
                assert!(columns[13].parse::<f64>().unwrap() >= 3.0);
                assert_eq!(columns[15], "8");
                assert!(columns[16].parse::<usize>().unwrap() >= 2);
                assert!(columns[17].parse::<f64>().unwrap() < 9.0);
            }
        }
        assert!(String::from_utf8_lossy(&result.stderr).contains("valid FLOAT 5; valid FIX 27"));
    }
}
