# Mixed-clock SPP MAD and RAIM refits

The opt-in QC API and `kinematic_spp_quality` RTK option connect the pinned
native MAD refit and leave-one-out RAIM selection to actual SPP solves. After
one MAD mask, a valid child is adopted immediately. Otherwise all single-SV
exclusions are fitted independently when at least two residual degrees of
freedom exist; the first strict best needs both 25% and 1 m RMS improvement.
Children never recurse. They retain the baseline's position and separate
clock initial estimates, and only its used tracking rows. Candidate failure
does not mutate the baseline or another candidate's clocks.

Each child and eventual baseline faces GDOP 50. Rust uses a clock for every
represented constellation and the actual code-fit design. **Pinned native
SPP's final DOP uses a four-column common-clock geometry.** Native's MAD
stored-seed restart and initially empty local ISB map also differ from Rust's
explicit baseline initialization. Order, noise and numeric QC thresholds are
ported; this is not a whole-native SPP equivalence claim. Full scope and
predeclared experiment are in the frozen [design](../../spp-quality-control.md).

The standalone result records original rows, final QC dispositions, all
refit exclusions/RMS/geometry/errors, adopted branch, final clocks and full
covariance. Nested `spp.rejected` retains navigation/elevation rejections.
RTK keeps a compact successful summary; terminal seed failure records its
error and follows the existing position fallback. SPP rejection does not
remove independently selected RTK measurements, and SPP covariance is not
fused. None/GPS/ElevationOnly/QC-disabled defaults and RTK gates are unchanged.

## Validation

Seven added tests cover 32 independent physical moving epochs within the
existing 3 mm / 1e-11 s truth limits, exact clean SPP/FLOAT covariance, a
100 m G01 error correctly excluded by RAIM, raw/alternate/unsupported rows,
per-clock redundancy, first failed candidates, GDOP boundary and seed fallback.
The 80 actual native variance/MAD/selector utility cases from the previous
archive are reused unchanged, with no new native build or whole-solver run.
All 383 Rust tests, release examples, fmt, all-target clippy with warnings
denied and release build pass; qualification logs and hashes are retained.

A synthetic limitation is deliberately exposed: with a 500 m G01 code error,
native-style MAD masks good G07, adopts that refit before trying RAIM, and
returns about 738 m position error with GDOP about 4.14. The regression checks
that behavior and provenance. It does not count it as accurate positioning.
Final residual-RMS/chi-square gates and robust reweighting remain separate work.

## PPC comparison

All six cases retain the prior mixed-clock/native-variance seed profiles,
no detector / slips / slips plus existing prefit 3 m gate, QC disabled/enabled:
36 executions. The 18 disabled controls match every previous CSV column.
All source/input/binary/configuration hashes precede execution. Every event,
invalid and unmatched row stays in the archive. Counts refer to 11,733 exact
paired 1 Hz epochs, not all 58,634 rover epochs at 5 Hz.

| FLOAT profile | Without QC | With QC | Accepted FLOAT seed fallbacks |
|---|---:|---:|---:|
| no detector | 9,728 | 9,728 | 0 → 4 |
| slip detector | 9,856 | 9,855 | 0 → 4 |
| slips + prefit 3 m | 4,976 | 6,353 | 0 → 0 |

The four accepted fallbacks are `SppGeometryTooWeak` (Nagoya run2 one,
Nagoya run3 two, Tokyo run1 one). These diagnostics cover accepted FLOAT
epochs, not every attempted seed or terminal processing failure. Without the
detector the FLOAT sets are identical; with it one Tokyo run1 epoch is lost.
Accepted successful seed QC branches and counts:

| Profile | Baseline | Adopted MAD | Adopted RAIM | RAIM refit attempts |
|---|---:|---:|---:|---:|
| no detector | 4,552 | 4,250 | 922 | 70,267 |
| slips | 4,589 | 4,327 | 935 | 70,915 |
| slips + prefit 3 m | 3,700 | 2,220 | 433 | 53,366 |

These attempts have no QC recursion. They are operational counts, not a runtime
or memory acceptance benchmark. With no detector, the following antenna-origin
example uses common successful seed epochs for seed columns and common FLOAT
epochs, including fallbacks, for FLOAT columns. All entries are 3D RMS metres:

| Case | Seed without QC | Seed with QC | FLOAT without QC | FLOAT with QC |
|---|---:|---:|---:|---:|
| Nagoya run1 | 33.557 | 26.796 | 8.562 | 7.956 |
| Nagoya run2 | 32.549 | 23.224 | 17.202 | 17.163 |
| Nagoya run3 | 35.384 | 30.938 | 27.432 | 27.689 |
| Tokyo run1 | 51.090 | 37.793 | 23.078 | 22.932 |
| Tokyo run2 | 16.931 | 11.399 | 15.650 | 15.631 |
| Tokyo run3 | 14.348 | 14.065 | 24.113 | 24.293 |

All six common seed RMS values fall, but FLOAT RMS worsens on Nagoya run3 and
Tokyo run3. Maximum FLOAT error worsens on several cases. The gated profile
gains outputs on every case, while losing some earlier epochs:

| Case | FLOAT without → with QC | All-output 3D RMS without → with | Common FLOAT 3D RMS without → with |
|---|---:|---:|---:|
| Nagoya run1 | 637 → 713 | 2.833 → 2.354 | 1.919 → 1.676 |
| Nagoya run2 | 775 → 1,000 | 11.146 → 5.001 | 11.106 → 2.992 |
| Nagoya run3 | 277 → 389 | 16.985 → 17.687 | 13.726 → 16.322 |
| Tokyo run1 | 1,012 → 1,435 | 7.929 → 12.181 | 6.723 → 6.221 |
| Tokyo run2 | 752 → 980 | 3.369 → 3.722 | 3.366 → 3.368 |
| Tokyo run3 | 1,523 → 1,836 | 3.438 → 4.965 | 3.024 → 2.908 |

Four gated all-output RMS values worsen. Nagoya run2 maximum falls
212.858→92.811 m, while Nagoya run3 rises 100.122→134.499 m and Tokyo run3
rises 39.734→58.170 m. The all-output sets differ, so common-set results must
also be considered. `report.json` and `analysis.json` retain horizontal,
vertical, 3D and percentile/max scores under both reference origins and both
conflicting Tokyo lever arms, plus status transitions and gained/lost counts.

## Worst-point diagnostics and limits

`outliers.json` retains each case's five worst gated QC FLOAT epochs, all raw
CSV diagnostics, all reference-origin errors and the corresponding no-QC row.
Selection is after solving and cannot feed a gate or configuration. Each
case's worst point was previously INVALID; all have only three RTK DD rows:

| Case | GPS week/TOW | Position error m | Seed error m | Final SPP RMS m | SPP GDOP | RTK GDOP |
|---|---|---:|---:|---:|---:|---:|
| Nagoya run1 | 2325/550877 | 31.466 | 18.382 | 1.894 | 8.375 | 27.334 |
| Nagoya run2 | 2323/556902 | 92.811 | 100.658 | 0.351 | 10.554 | 129.261 |
| Nagoya run3 | 2325/554421 | 134.499 | 130.962 | 0.052 | 5.089 | 494.450 |
| Tokyo run1 | 2324/188404 | 138.288 | 189.745 | 2.041 | 7.518 | 115.603 |
| Tokyo run2 | 2324/177899 | 44.127 | 25.882 | 1.676 | 3.503 | 43.421 |
| Tokyo run3 | 2324/182002 | 58.170 | 66.786 | 0.785 | 10.321 | 10.731 |

Nagoya run2/3 and Tokyo run3 use RAIM children with six rows and two clocks:
only one residual degree of freedom remains. Small residual and acceptable
SPP geometry do not guarantee an accurate seed; the independently selected
RTK set can have much worse geometry. RTK GDOP remains a diagnostic here,
without enabling its separate optional gate. Next work should address
redundancy, observation quality, seed/ambiguity initialization and integrity.

The earlier immutable whole-native baseline (10,585 RTK including 144 FIX) is
reused for common-time scores, not rerun. Standalone code-only SPP selection
still differs from native RTK's code/phase adapter. Rust FIX policy is not
evaluated here. No held-out, runtime/memory or practical RTK/PPP acceptance is
claimed. The broad checklist remains 18/24 (75%).

## Replay

```bash
source /workspace/.cloud-setup/activate.sh
cargo build --locked --offline --release --example ppc_rtk_float
python3 -B scripts/run_ppc_spp_qc.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_float --output-dir /tmp/spp-qc-replay
python3 -B docs/comparison/2026-10-21-spp-qc/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
```

The direct mode adds `-sppqc`, for example
`kinematic-galileo-spp-slips-prefit3-ggspp-nativevar-sppqc`.
The verifier recomputes all archived scores, exact controls, worst-point
diagnostics, reused native component provenance and qualification/hash checks
without running a positioning solver. Source hashes refer to this milestone;
later changes should be verified from its commit.
