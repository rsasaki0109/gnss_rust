# First full-rate moving PPP on PPC

All six unmodified rover files are processed at their original 5 Hz rate:
**58,634 epochs per profile**, four predeclared profiles, **234,536 output
rows**. Profiles select 1C/2W or 1C/2L and calibrated CODE phase OSBs or explicit
raw-carrier FLOAT. Matching CODE rapid SP3/CLK/absolute OSBs and the IGS20 GPS
satellite ANTEX excerpt are used. Every invalid/event/solution row is retained.
No base observations, SPP, Doppler, IMU or reference truth enter the filter.
Only the RINEX header supplies initial position; velocity/clock start at zero.

**This baseline does not meet practical PPP acceptance.** Calibrated and raw
phase have exactly the same success times and statuses in all 12 comparisons;
their maximum position difference is 0.619 micrometres. Thus the phase OSB
choice does not explain the poor continuity in these profiles.

| Case | Raw epochs | 1C/2W FLOAT | 1C/2L FLOAT | 1C/2W 3D RMS, m | 1C/2L 3D RMS, m |
| --- | ---: | ---: | ---: | ---: | ---: |
| Nagoya run1 | 7,602 | 5 | 6 | 11.624 | 9.717 |
| Nagoya run2 | 9,451 | 1 | 1 | 7.787 | 5.824 |
| Nagoya run3 | 5,201 | 118 | 19 | 8.374 | 6.950 |
| Tokyo run1 | 11,928 | 8 | 2 | 5.382 | 5.543 |
| Tokyo run2 | 9,151 | 5 | 107 | 1.703 | 3.497 |
| Tokyo run3 | 15,301 | 4 | 17 | 4.611 | 3.619 |
| Total | **58,634** | **141** | **152** | — | — |

Errors above use the calibrated-phase profile and the antenna-origin
hypothesis, on successful epochs only. They do not describe full-route
accuracy. The JSON also retains the IMU-origin hypothesis and both conflicting
Tokyo lever arms, all unmatched-reference counts, longest successful streaks,
gaps and every terminal rejection reason. Reference timestamps are joined at
1 microsecond resolution with uniqueness guards; no truth interpolation or
solver timestamp change occurs. No successful position lacks reference data.

The solver uses its existing gate/convergence settings, constant-velocity
model and constant independent per-frequency code/phase sigmas 1 m/0.01 m.
The sigmas are the native nominal values, but native elevation/variance/QC
processing is not yet ported here. Satellite IF-PCO and nominal-yaw wind-up
are enabled. Receiver antenna calibration is unknown, satellite PCV and
solid Earth tide are absent. No correction is claimed complete by omitting it.
`arc_resets` counts reset reasons on decisions, including missing selected
satellites; it is not a count of detected cycle slips.

Most later epochs fail the unchanged innovation gate or lack enough usable
GPS rows. Nagoya run1 accepts only five 1C/2W epochs, with estimated ZTD
0.017–0.122 m and about 12 m position error despite sub-metre code residuals.
This exposes a bad initial-state/convergence/recovery regime, rather than a
successful real-data PPP result. Thresholds were not relaxed.

After the frozen baseline, an adapter calls the actual unmodified native
`PPPProcessor::initializeFilter` six times with independent RINEX-header
positions and explicit zero-clock seeds. Native precise-mode initial diagonal
variances are XYZ **100 m²**, velocity **25 (m/s)²**, clock **1e8 m²**, ZTD
**25 m²**, with climatological ZTD about **2.55 m**. Current Rust uses XYZ
1e6 m², velocity 100 (m/s)², clock 1e8 m², ZTD 100 m² and initial ZTD 2.3 m.
This component diagnostic establishes a concrete difference; it does not prove
that it causes every rejection or represent native end-to-end positioning.
Next compare a declared native initialization profile, atmosphere/noise models,
and recovery on matching inputs before any whole-solver acceptance claim.

Linux `wait4` records each individual process: 3.82–15.92 seconds including
input parsing/CSV output, peak RSS 43.29–43.74 MiB. These are descriptive
single-run measurements with no native timing control, not benchmark acceptance.
The first harness attempts failed before solver execution on reference-column
names and an unavailable GNU time command; both are retained. No numeric
threshold or solver change followed those harness failures.

The library source is unchanged from the **406-test** qualified 77b7886
milestone; new example/release/fmt/strict-Clippy checks pass. Run
`python3 -B docs/comparison/2026-10-21-ppc-ppp-first/verify.py` to recompute
statuses, comparisons and integrity; add `--verify-local-inputs` to rehash
cached products and recompute both-origin errors against reference.csv.
Source hashes are specific to this milestone: use its commit for later replay.

Whole-native PPP, missing correction models, held-out accuracy/status/runtime
and memory acceptance remain pending. Overall checklist: **18/24 (75%)**.
