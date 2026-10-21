# Explicit SPP/Doppler recovery for PVA PPP

Pinned libgnss++ motion recovery is connected to the Rust PPP filter as an
opt-in policy. Defaults remain unchanged. This is a component port and a
fixed PPC diagnostic experiment; practical RTK/PPP parity is unaccepted.

`process_epoch_with_seed` validates same-epoch SPP status, finite position,
caller-declared frame/reference and provenance. It never imports SPP clock
bias or Doppler clock drift. PVA initialization uses SPP XYZ, finite available
velocity or rest, and zero acceleration. Caller clock/ZTD initialization stays
in place; climatological ZTD is not recomputed at the SPP position. Later
recovery follows native
strict mean-XYZ-variance > 3600 m² or total position-disagreement NIS > 16.27.
It resets XYZ/V/A rows and columns, retaining the other marginal. An optional
Doppler update uses the velocity covariance marginal and full state
cross-correlations. Geometry admission previews the corrected head; the full
prior is recreated so velocity is fused once. All numeric operations remain
proposed until the GNSS update and arc commit succeed.

A same-epoch code-derived restart is a data-dependent heuristic, not an
independent position factor. On a restart, seeding velocity from the same
velocity subsequently fused gives zero initial innovation; this is not a
quality test of the seed. The source Doppler fit has its own unchanged gates.
Disabling direct initialization does not disable the separate recovery triggers,
which can also fire on the first proposed prior.
No missing velocity/covariance is replaced with a zero measurement.

## Independent components and tests

The native adapter calls actual pinned `predictState(0, seed)`,
`applyVelocityMeasurement` and `seedPositionDisagreement`: 32 cases each,
plus 32 actual initializations, per repeat. `predictState` can make additional
internal disagreement calls; the 32 disagreement calls reported here are
direct helper invocations. Dense covariance traces cover full cross terms,
strict restart boundaries, missing/invalid seeds and covariance, velocity
availability, and enabled/disabled gates. Two complete generations are byte
identical and match the existing relative 1e-10 component bound. Native source
and libraries are unchanged. The initial adapter indentation warning under
-Werror, failed command and exact source are retained in `initial-native-harness`.

Rust uses checked SPD matrices and Joseph covariance updates; native uses
LDLT/subtractive updates. Rust explicitly rejects zero/indefinite velocity
noise even if native would accept it. These validation differences are not
claimed to be bit-identical behavior on invalid inputs.

Five primitive and five PPP integration tests bring the total from 427 to
437. Tests include code-only SPP seeding on each 64-epoch clean/noisy accelerating
car epochs with the existing 1 cm/1 m bounds, seed time/frame/reference/status/
provenance failures, exactly-once velocity covariance, ignored seed clock,
missing/bad velocity covariance and optional velocity rejection, event
handling and full numeric transactionality on failed GNSS updates. Unit
velocity injections test covariance behavior and are not labeled measured
Doppler. Standard tests/fmt/clippy/release/examples pass.

## Fixed PPC experiment

Before execution, freeze native head, IndependentEpoch clock, PVA,
InterruptAllArcs and all existing GNSS gates/noises. Run six cases × two GPS
IF pairs × three modes: shadow, position, Doppler. There are 36 profiles and
351,804 rows from all 58,634 original rover epochs. The 12 shadow projections
must match the previous 35-column PPP output byte-for-byte.

Every profile computes the same independent GPS/Galileo broadcast code SPP
with native variance and default MAD/RAIM/GDOP QC, and carrier-locked Doppler
with per-system clocks. All 23 source fields must be identical across modes
and IF pairs. Position mode supplies no velocity. Shadow computes but never
uses the seed. SPP starts from the original header, independent of PPP state.
The explicit approximate broadcast WGS84/IGS20 compatibility assertion is not
a frame transform or receiver calibration. No reference/IMU enters the solver.

CSV keeps original PPP status/fields and appends source, proposed/committed
recovery diagnostics and a separate service output. Service selects PPP when
available, otherwise the current valid SPP solution. SPP is never counted as
PPP. Both reference origins and both Tokyo arms are scored, along with common,
newly accepted and lost PPP epochs and separate service comparisons.

The four-profile read-only probe at the previous 134 m Nagoya 3 epoch is
retained, including profiles still over 100 m wrong. GPS/Galileo native QC
produced 10.589 m antenna-origin error and 0.0295 m/s speed there. The currently
qualified probe reproduces its original stdout and stderr byte-for-byte.
The probe is a motivation, not an accepted recovery result.

All 1C2W modes are also compared to both previously frozen native whole-PPP
initialization profiles. These 36 comparisons separate PPP-only from combined
valid output. Features differ: native effective precise-clock random walk and
raw-phase/GPS baseline versus Rust independent clock, calibrated phase and
GPS/Galileo seeds/QC. No new whole-native PPP calls or feature-matched runtime
acceptance are claimed. Process measurements include load, seed computation
and CSV, and can overlap qualification work; they are descriptive.

## Results

PPP epochs / antenna-hypothesis 3-D RMS (m), retaining each profile's own
success set. These success-set RMS values are not common-epoch comparisons.

| Run | Pair | Shadow | Position seed | Position + Doppler |
|---|---|---:|---:|---:|
| nagoya-run1 | 1C2W | 681 / 3.331 | 1796 / 9.737 | 1603 / 9.951 |
| nagoya-run1 | 1C2L | 1300 / 3.675 | 1379 / 3.807 | 1211 / 3.613 |
| nagoya-run2 | 1C2W | 628 / 5.704 | 1631 / 13.689 | 1418 / 11.982 |
| nagoya-run2 | 1C2L | 1001 / 11.916 | 2049 / 13.088 | 1670 / 14.148 |
| nagoya-run3 | 1C2W | 573 / 12.617 | 748 / 15.624 | 681 / 14.072 |
| nagoya-run3 | 1C2L | 562 / 17.970 | 689 / 23.890 | 587 / 22.515 |
| tokyo-run1 | 1C2W | 1702 / 4.074 | 2231 / 6.976 | 2243 / 7.096 |
| tokyo-run1 | 1C2L | 1033 / 6.436 | 1366 / 6.240 | 1269 / 6.074 |
| tokyo-run2 | 1C2W | 900 / 4.254 | 2098 / 3.721 | 1884 / 4.476 |
| tokyo-run2 | 1C2L | 446 / 12.183 | 1405 / 8.296 | 1204 / 5.894 |
| tokyo-run3 | 1C2W | 4137 / 4.623 | 5174 / 5.302 | 4805 / 5.148 |
| tokyo-run3 | 1C2L | 424 / 2.233 | 3377 / 5.382 | 3072 / 5.157 |

Aggregate PPP counts are 1C2W **8,621 → 13,678 / 12,634** and 1C2L
**4,766 → 10,265 / 9,013** (position / Doppler). Fractions are 14.70% →
23.33% / 21.55% and 8.13% → 17.51% / 15.37%. IF-pair successes overlap and
must not be summed. Common Nagoya 1 1C2W Doppler points (645) worsen from
3.287 to 3.304 m; Nagoya 2 (577) worsen 5.715 to 5.916 m. Tokyo 2 1C2L
success-set RMS improves 12.183 to 5.894 m. More solutions are not sufficient
for performance acceptance.

At the former bad Nagoya 3 1C2L epoch (2325/554542.2), the fixed native-head
shadow is 134.092 m wrong. Both seed modes reject PPP with PostfitRejected
and publish a separate 10.589 m SPP fallback. Neither makes a motion restart
at this exact epoch; changed earlier trajectories/velocity affect rejection.
This is not a 10.589 m PPP result. Elsewhere Doppler PPP still reaches
126.537 m. At Nagoya 1's 69.442 m PPP error, its SPP seed is only 4.842 m
wrong: five code rows leave no code residual redundancy after XYZ/clock/ZTD.
Code RMS is 0.392 m, phase RMS 4.7e-9 m and NIS 1.359, so small residuals
alone do not establish position accuracy.

Combined valid output is 57,782/58,634 (98.55%) in each profile, but SPP
fallback has severe outliers. In 1C2W Doppler mode, combined antenna-hypothesis
RMS for Nagoya 1/2/3 and Tokyo 1/2/3 is respectively
58.672 / 54.652 / 32.868 / 43.937 / 14.377 / 16.094 m. The maximum at Nagoya 1
is **1,278.220 m**, from SPP, not PPP. These QC-admitted source failures remain
in all captures and scores. Continuity counts must not be presented as safe
or accurate coverage. Measurement quality/redundancy, ZTD/height coupling,
source SPP quality and multi-GNSS PPP remain priorities.

The original captures use the legacy header `seed_gps_clock_s`. The SPP API
actually returns the first represented system's clock, so this fixed label
was too specific. After the fixed experiment, the example header is corrected
to `seed_first_system_clock_s`, with the first system explicitly appended to
`seed_source`. Its numeric processing is unchanged. `clock-label-fix` retains
both original source files and six full shadow replays, comparing every one
of the 72 fields after only the clock-column/source-label normalization.
All 57,782 valid SPP replay rows have GPS as their first clock in these six
cases; the new explicit label also handles Galileo-only SPP. The original
36 captures and predeclared hashes remain unmodified; the
verifier checks their source snapshots and the current source/binary bridge.
Full tests concern the unchanged library; fmt/clippy/example build are repeated
after the metadata fix. The fixed-label experiment's runner accepts either
column name so its original evidence remains readable.

## Reproduction

```bash
source /workspace/.cloud-setup/activate.sh
python3 -B scripts/build_ppp_motion_fixtures.py --output-dir NEW_NATIVE_DIRECTORY
cargo build --locked --offline --release --examples
python3 -B scripts/run_ppp_motion_seeds.py --output-dir NEW_OUTPUT_DIRECTORY --native-fixtures-dir NEW_NATIVE_DIRECTORY
python3 -B scripts/analyze_ppp_motion_seeds.py --archive NEW_OUTPUT_DIRECTORY
python3 -B docs/comparison/2026-10-21-ppp-motion-seeds/verify.py --verify-local-inputs
```

The archive verifier checks capture hashes, source/input/native provenance,
repeat fixtures, qualification, all 12 old PPP projections, equal seed sources,
separate PPP/SPP status and all reference-origin/common/worst-point metrics.
Older archive verifiers use their recorded source commit; their evidence is
unchanged by this solver development.
