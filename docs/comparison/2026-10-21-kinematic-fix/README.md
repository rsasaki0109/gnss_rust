# Initial moving GPS/Galileo FIX policy

An opt-in full-set policy now checks moving GPS/Galileo integer candidates
against timestamped ECEF Doppler velocity. The public default remains GPS
static; the kinematic motion standard-deviation cap defaults to 1 m/s. No
FLOAT state, covariance or ambiguity is replaced by a fixed result. Partial
AR, integer holding and practical RTK/PPP acceptance remain pending.

`RtkFixPolicy::for_kinematic` requires at least two confirmations and explicit
motion evaluation. Motion covariance must be finite, symmetric and positive
definite. Missing or overly uncertain motion clears confirmation/anchor
history; invalid motion/time fails before advancing policy history. Samples
more than 2 s apart restart confirmation. Between candidate/fixed positions,
the trapezoidal mean of the endpoint velocities predicts displacement; the
existing 0.1 m history bound applies to its residual. Covariance controls
admission only and never enlarges this bound or updates FLOAT.

All existing local defaults remain: ratio 3, five lock epochs, at least four
within-system DD pairs, two confirmations, 1 m maximum FLOAT correction,
0.05 m maximum fixed coordinate standard deviation, phase RMS/absolute
0.01/0.03 m, code RMS 1 m, code-phase 0.5 m and NIS/row 9. Mixed integer labels
are normalized independently per system. The optional all-reference mode
requires the same candidate to pass every Cartesian combination of system
references, without another integer search or a favorable-reference choice.

Independent C++ physical fixtures give 27/32 mixed moving FIX results in both
reference modes and 59/64 GPS moving FIX results; all fixed positions remain
within the existing 3 mm truth tolerance. These policy tests supply
truth-derived velocity and 0.01 m/s coordinate standard deviations. They are
component/synthetic evidence, rather than whole-chain Doppler accuracy tests.
The eight added tests include independent reference changes, invalid/uncertain
motion, exact uncertainty boundaries, duplicate epochs, reset, LLI, gaps and
inconsistent velocity with recovery. All 352 Rust tests, examples, fmt,
all-target clippy and release build pass; logs and commands are retained.

## PPC results

The release example processes every paired/unmatched epoch on all six runs,
with GPS/Galileo independent-clock Doppler and the previously evaluated
finite-phase/clear-LLI admission condition. The FLOAT positions, residuals,
covariances, satellites, references and ambiguity-reset counts exactly match
the prior controls. Raw Doppler velocities and full 3x3 covariances also match.
There are 11,733 common 1 Hz pairs, drawn from 58,634 original 5 Hz rover
epochs. Unmatched or failed FLOAT epochs acquire no invented velocity/FIX.

A shape-only Nagoya run1 smoke, without reading reference truth, found all
velocity coordinate standard deviations above the default 1 m/s. Before
scoring, caps 1/2/4 m/s were declared as a sensitivity grid, alongside both
reference modes and with/without the existing CMC/Doppler slip detector. This
is 72 runs. The default cap and every local/motion-residual threshold remain
unchanged; the larger caps are experiments and were not selected as defaults.

| FLOAT/slip profile | Local reference gates | FIX at cap 1 m/s | 2 m/s | 4 m/s |
|---|---|---:|---:|---:|
| GPS/Galileo | selected | 2 | 6 | 6 |
| GPS/Galileo | all combinations | 0 | 2 | 2 |
| GPS/Galileo + slip detection | selected | 1 | 1 | 1 |
| GPS/Galileo + slip detection | all combinations | 0 | 0 | 0 |

Each configuration retains 9,727/9,855 total RTK outputs without/with slip
detection, with some FLOAT statuses replaced by FIX. At cap 2 m/s the six
selected-reference FIX points are two Nagoya run1, one Tokyo run1 and three
Tokyo run3; the all-reference intersection retains one Tokyo run1 and one
Tokyo run3. The longest FIX sequence is two 1 Hz epochs. Default cap 1 m/s
selected-reference FIX occurs only at two isolated Tokyo run3 epochs.

`report.json` retains all RTK and FIX-only horizontal/vertical/3D RMS, P50/P95
and maximum errors under the antenna, IMU and both conflicting Tokyo lever-arm
assumptions. The maximum selected-reference FIX antenna error is 0.137114 m
in Nagoya run1, 0.728253 m in Tokyo run1 and 0.063879 m in Tokyo run3. The
Tokyo run1 point also passes all-reference validation and the motion check;
its IMU-assumption errors are 0.854932/1.821790 m. This demonstrates that
repeated labels and consistent motion do not establish absolute integer
correctness. No precision or false-FIX acceptance claim follows from a local
pass, and the unknown reference origin has not been inferred from these points.

The immutable prior whole-C++ GPS/Galileo default-SPP run is compared read-only
in `analysis.json`; it was not rerun here. It yields 10,585 RTK outputs,
including 144 FIX, plus 878 SPP and 270 invalid pairs. No Rust FIX epoch
coincides with a native FIX epoch. Native/Rust status cross-tabs, common-RTK
errors for all reference hypotheses, unmatched success counts and FIX sequence
lengths are retained. Existing seed, AR, QC and filter differences remain;
this is a practical gap measurement, not solver equivalence. Runtime/memory
acceptance and held-out validation have not been performed.

At cap 2 m/s, selected-reference decisions on 9,727 valid FLOAT snapshots are
3,648 insufficient-lock, 2,577 uncertain-motion, 2,163 ratio, 825 code-phase,
149 phase-residual, 70 code-residual, 55 correction, 218 unavailable-motion,
16 confirming and 6 accepted. With slip detection, insufficient-lock rises
to 5,511 and only one FIX remains. Next dependencies are usable locked
subsets/partial AR, observation and motion uncertainty quality, and validation
against the native policy. Increasing the velocity cap alone has little effect
on FIX coverage and does not resolve the accuracy gap.

## Reproduction and retained failure

Source, input, binary and configuration hashes are in `predeclared-inputs.json`.
Every run has exact commands, return codes, raw stdout/stderr and hashes.
Input locations can be supplied again to the read-only verifier:

```bash
source /workspace/.cloud-setup/activate.sh
cargo build --locked --offline --release --example ppc_rtk_kinematic_fix -j2
python3 -B scripts/run_ppc_kinematic_fix.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_kinematic_fix \
  --output-dir /tmp/ppc-kinematic-fix-replay
python3 -B docs/comparison/2026-10-21-kinematic-fix/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
```

`initial-incomplete/` preserves the first 13 successful solver executions,
12 completed summaries, original manifest/runner and comparator failure. The
Python `x**0.5` recomputation differed from Rust `sqrt` by rounding at Nagoya
run2; replacing it with `math.sqrt` repaired the evidence check. No solver
threshold changed. A fresh complete 72-run execution produced the primary
report. The failed partial report is not counted as a complete experiment.

The verifier recomputes all controls, decisions, errors, status comparison,
published velocity predictions and accepted last-FIX displacement checks.
Candidate-only absolute positions are not reconstructed from FLOAT positions;
those history transitions are covered by the synthetic policy tests. The
9-decimal-place ECEF CSV has a documented 3 nm serialization/roundoff bound
for the independent displacement check, without changing the solver's 0.1 m
acceptance bound. `evidence-sha256.json` seals all retained artifacts.
