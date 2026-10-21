# GPS/Galileo independent-clock SPP seed

`spp_multiclock::solve_gps_galileo_epoch` adds explicit GPS 1C / Galileo 1C/1X
code SPP with separate receiver clocks. The state contains ECEF position and
one range-clock bias per represented system, with full joint covariance. No
clock prior or covariance floor repairs deficient geometry. Single GPS uses
the existing four-state solver exactly; a singleton constellation adds only
a clock, with no position information. Clock groups are derived from usable
rows at each iteration, including the final residual/covariance evaluation.

The existing broadcast transmit time, first-order Sagnac derivative,
Saastamoinen/GPS Klobuchar E1-frequency model, elevation weights and code noise
are retained. Group delay follows the pinned default legacy-primary TGD/BGD
convention; navigation source bits, primary and secondary delay fields remain
in the successful result's provenance. No unambiguous I/NAV/F/NAV source is
invented. The new API retains all raw observation rows and used/alternate/
invalid/unsupported dispositions. Its nested `spp.rejected` contains final
navigation/elevation rejections; complete row selection lives in `observations`.
Phase, Doppler and LLI values are preserved but do not enter the code fit.
Nonzero receiver epoch clock offsets are explicitly unsupported; nonfinite
offsets fail and zero/None distinctions are retained. The existing GPS API
remains unchanged.

RTK opts in with `kinematic_spp_constellations=GpsGalileo` alongside the existing
optional SPP seed. Default None/GPS seed behavior is preserved. Only a
successful position initializes the kinematic epoch: no SPP covariance is
fused, and SD ambiguities are retained. Failure uses the previous fallback
without invented clock output. Diagnostics name requested scope, satellite
count and each available clock. `ppc_rtk_float` accepts `-ggspp` after an
explicit SPP mode and appends four diagnostic columns only when enabled.

## Validation and scope

The independent native physical mixed fixture already includes distinct
2 microsecond GPS and -3 microsecond Galileo rover clocks, eight satellites,
broadcast state/clock/TGD, Sagnac and troposphere on 32 moving epochs. All
new mixed SPP positions meet the existing 3 mm truth tolerance; both clocks
meet 1e-11 s. Three GPS plus three Galileo solve without a shared-clock
assumption. Nine added tests also cover exact single-GPS equivalence,
singleton information, rank failure, tracking alternatives, missing/invalid
code, preserved phase/LLI, stale navigation, events, epoch offsets, stateless
repeat evaluation and RTK seed success/fallback.

The independent-clock QR previously used by Doppler is moved to a shared
internal module. `clock-QR-source-transfer.json` proves that its body is
unchanged except for type name, visibility and rustfmt whitespace. Existing
eleven native single-system Doppler cases, sixteen independent Eigen mixed
cases and all related regression tests pass. This does not establish whole
native SPP equivalence. The full qualification is 371 Rust tests plus examples,
fmt, all-target clippy and release build, with commands/log hashes retained.

## PPC comparison

All six cases retain GPS seed versus mixed seed, no detector / slip detector /
slips plus the existing 3 m prefit gate: 36 executions. Eighteen GPS-seed
controls reproduce every prior CSV field exactly. Thresholds and defaults are
unchanged; source/input/binary/configuration hashes precede execution. The
11,733 paired 1 Hz epochs and every unmatched/invalid/event row are retained;
these counts do not imply 5 Hz rover coverage.

| RTK FLOAT profile | GPS seed | Mixed seed | GPS-seed fallbacks → mixed |
|---|---:|---:|---:|
| no detector | 9,705 | 9,728 | 116 → 0 |
| slip detector | 9,834 | 9,856 | 118 → 0 |
| slips + prefit 3 m | 3,458 | 3,880 | 20 → 0 |

Fallback counts refer only to accepted FLOAT epochs. They are not a count of
all attempted SPP solves or original rover epochs. In the first two profiles,
every prior FLOAT epoch remains present; the mixed seed adds 23/22 epochs.
The new standalone SPP API needs only valid code, so its measurement selection
can differ from the complete code/phase rows used by the native RTK adapter.
Full native SPP weighting, FDE and QC have not been ported here.

The following antenna-assumption example uses common successful SPP-seed
epochs for the seed columns and common FLOAT epochs (including GPS fallback)
for FLOAT columns. Both origins and both Tokyo lever-arm alternatives are
retained in the machine report; the reference origin is not inferred here.
All entries are 3D RMS in metres:

| Case | GPS seed | Mixed seed | FLOAT with GPS seed | FLOAT with mixed seed |
|---|---:|---:|---:|---:|
| Nagoya run1 | 58.208 | 48.669 | 9.264 | 8.993 |
| Nagoya run2 | 102.746 | 53.155 | 19.923 | 17.301 |
| Nagoya run3 | 63.703 | 47.988 | 27.813 | 27.489 |
| Tokyo run1 | 97.189 | 60.802 | 27.843 | 22.674 |
| Tokyo run2 | 33.969 | 24.323 | 15.690 | 15.689 |
| Tokyo run3 | 20.984 | 16.377 | 23.726 | 23.829 |

This is mixed evidence for practical RTK accuracy. Without the detector,
Tokyo run1 maximum antenna-assumption FLOAT error drops 741.825→178.872 m,
while Nagoya run1 maximum worsens 145.872→146.955 m and Tokyo run3 RMS slightly
worsens. With slips plus the 3 m prefit gate, Nagoya run2 loses 43 FLOAT
epochs (682→639) and maximum error worsens 30.095→214.901 m; Tokyo run1 gate
maximum also worsens 43.836→134.492 m. A better average SPP seed and small
prefit residual do not establish reliable geometry/observation quality.
No threshold is relaxed and no new default is selected from these results.

`report.json` retains all-output and accepted-SPP-seed errors under every origin
hypothesis. `analysis.json` includes GPS/mixed status cross-tabs, added/lost
epochs, common FLOAT/common successful seed errors, and comparison with the
immutable whole-native GPS/Galileo default-SPP captures. Native is not rerun
here: it retains 10,585 RTK outputs, including 144 FIX, 878 SPP and 270 invalid
pairs. This new comparison evaluates Rust FLOAT without the FIX policy;
native AR, SPP selection, QC and other filter differences remain. Held-out
validation and runtime/memory acceptance are pending. Practical RTK/PPP parity
is not accepted and checklist progress stays 18/24 (75%). Next work is SPP
observation quality/FDE and investigation of the worsened weak-data cases.

## Replay

```bash
source /workspace/.cloud-setup/activate.sh
cargo build --locked --offline --release --example ppc_rtk_float -j2
python3 -B scripts/run_ppc_multiclock_spp.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_float \
  --output-dir /tmp/ppc-multiclock-spp-replay
python3 -B docs/comparison/2026-10-21-multiclock-spp/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
```

The verifier recomputes all retained metrics/controls and the QR transfer proof
without running positioning solvers. Source-hash verification requires this
archive's source commit; old manifests/captures are never rewritten.
`evidence-sha256.json` seals all artifacts.
