# Relative GPS DSBs under an explicit clock gauge

`GpsIfDsbClock` derives additive C1W/C2W/C1C/C2L code corrections only when a
caller explicitly declares GPS IF(C1W,C2W) as the clock gauge. This is neither
an absolute OSB product nor inferred compatibility with a precise CLK file.
Phase corrections remain absent. The existing absolute-OSB API still rejects
relative products, and no default solver setting changes.

For `DSB(C1W,C2W)=d`, IF coefficients `[a,b]` impose
`a*bias_C1W+b*bias_C2W=0`. The additive corrections are `-b*d` and `a*d`.
C1C additionally subtracts DSB(C1C,C1W); C2L adds DSB(C2W,C2L). Required edges
are fixed. Reverse orientation and ns/metre units are supported; missing,
expired or ambiguous forward/reverse edges fail, without substituting graph
paths or averaging inconsistent edges. Satellite source rows, SVN, units,
validity and signed coefficients remain in provenance; receiver rows are
excluded from this gauge. See the frozen [design](../../dsb-clock-gauge.md).

Source DSB correlations are unknown. `sum(abs(coefficient)*sigma)` is retained
as a conservative standard-deviation upper bound, not an independent-error
RMS or covariance and not silently added to measurement noise. Corrections
are materialized only at the requested epoch so half-open BIA validity does
not become an inclusive end-time admission.

## Independent validation

The C++ fixture solves four simultaneous equations with Eigen full-pivot LU
on 24 cases (96 corrections), using pinned native IF constants and additive
correction sign. It does not execute an upstream DSB gauge conversion routine
or a positioning solver. Compilation is warning-free with warnings denied;
two output generations are identical and header/source/compiler hashes are
retained. The original pinned source proof is linked by hash.

Seven added Rust tests cover all fixture values and bounds, every reversed
edge, mixed ns/metre inputs, injected raw-code correction, retained absence of
phase OSBs, half-open boundary, missing/alternate/ambiguous edges, unsupported
and duplicate requests, receiver exclusion and unknown correlation bounds.
The independent equation solve caught a sign error in the initial datum term;
the corrected implementation passes without loosening tolerances. All 390
Rust tests, release examples, fmt, all-target clippy with warnings denied and
release build pass, with compressed qualification logs.

## Real product and raw observation admission

Before processing, source, input, executable and configuration hashes were
frozen. The nine already acquired CAS products are inspected at GPS-calendar
noon on each product day. Every 32 GPS identifier / four-code request is
reported, including unavailable ones: **1,152 statuses, 1,089 available and
63 missing**. Each day has 121 available / seven missing. The missing rows
are C2L on G02/G13/G16/G19/G20/G21/G22, whose required C2W-C2L difference is
absent. No replacement value is invented.

All six PPC rover files are then streamed at their original 5 Hz rate,
including every raw epoch. An eligible pair requires finite positive code
and finite phase on exact tracking identities. Zero measured phase is valid;
missing/invalid fields are not converted to zero. LLI flags and missing LLI
counts remain separate. Events/nonzero offsets have explicit categories.
The two candidate tracking pairs are audited independently; no best pair is
chosen from reference truth or accuracy. Counts can overlap:

| Case | Raw epochs | Complete 1C/2W pairs | Complete 1C/2L pairs | Epochs with ≥4 unflagged 1C/2W | Epochs with ≥4 unflagged 1C/2L |
|---|---:|---:|---:|---:|---:|
| Nagoya run1 | 7,602 | 38,339 | 37,582 | 6,197 | 5,865 |
| Nagoya run2 | 9,451 | 43,777 | 42,847 | 6,514 | 6,131 |
| Nagoya run3 | 5,201 | 22,767 | 22,326 | 3,283 | 3,213 |
| Tokyo run1 | 11,928 | 50,718 | 38,003 | 7,547 | 5,030 |
| Tokyo run2 | 9,151 | 45,614 | 34,402 | 7,125 | 5,111 |
| Tokyo run3 | 15,301 | 80,709 | 61,381 | 11,710 | 8,876 |
| Total | 58,634 | 281,924 | 236,541 | 42,376 | 34,226 |

**Every complete raw pair has its required code DSBs available**; none fails
for missing, expired, ambiguous or other DSB admission. All raw epochs here
are normal observations with supported epoch offsets, and missing LLI pair-row
count is zero. Flags are counted per candidate pair; an L1 row can participate
in both candidate pairs. These are source/admission counts, not PPP solutions,
usable geometry, converged arcs or a success rate. Four raw pairs alone do
not prove observability of all PPP nuisance states.

`report.json` recomputes per-profile categories from every compressed CSV.
Per-product CSVs retain additive values and all original term identities,
line numbers, units and signed coefficients. Per-rover CSVs retain one row
per epoch and pair-specific raw/unflagged/DSB status counts. Executions and
hashes are stored for all 15 commands; original products remain in the cloud
dataset directory with their earlier acquisition provenance.

## Limits and replay

The audit binding explicitly says clock compatibility is unverified. No
antenna calibration/reference-point compatibility is established, no phase
OSB is generated, and neither code positioning nor PPP executes. This stage
removes the missing relative-code-gauge conversion dependency while preserving
the outstanding clock/phase datum, kinematic PPP, corrections, whole-native
comparison and runtime/memory work. RTK practical quality also remains open.
The broad checklist stays 18/24 (75%).

```bash
source /workspace/.cloud-setup/activate.sh
cargo build --locked --offline --release --example ppc_dsb_clock
python3 -B scripts/run_ppc_dsb_clock.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --product-root /workspace/.cloud-setup/datasets/PPC-Precise/igs-20261021 \
  --binary target/release/examples/ppc_dsb_clock --output-dir /tmp/dsb-clock-replay
python3 -B docs/comparison/2026-10-21-dsb-clock/verify.py --verify-local-inputs
```

The verifier rechecks archived metrics, correction term sums, all 58,634 epoch
counts, independent fixture provenance, qualification logs and recursive hashes
without running positioning. `--verify-local-inputs` also checks external
input hashes; without it the archive can be checked without datasets. Source
hashes refer to this milestone; use its commit after later source changes.
