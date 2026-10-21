# Native relative-variance AR preselection

An opt-in `for_kinematic_variance_subset` now applies the native one-shot
relative-variance prefilter after lock qualification. The threshold is the
upper median DD variance times 10, bounded below by 1e-4 cycles squared. If
fewer than four pairs survive, keep the unfiltered lock-qualified set, matching
native resolver flow. The configured stricter policy minimum also remains
binding. No existing local/motion threshold or public default changes.

Retain the original within-system reference for surviving target pairs and
omit systems with no remaining target. The subset view marginalizes omitted
SD states using the full principal covariance, including position/SD and
cross-system terms. It never fixes omitted integers, rewrites FLOAT or adds
a covariance diagonal floor. Selection uses the input DD basis; all-reference
candidate validation does not establish basis-invariant subset search.

`progressive_variance_drop_subsets` separately ports the native nested
largest-variance-first drop path, first-input tie order, minimum size and step
limit. This component does not run integer search or declare FIX and is not
connected to policy search yet. The policy still searches exactly one subset;
adaptive candidate testing/FFRT and integer holding remain pending.

## Verification

Forty independent cases call actual pinned C++
`filterPairsByRelativeVariance` and `buildProgressiveVarianceDropSubsets`.
They cover exact/adjacent thresholds, upper median for even length, zeros,
insufficient survivors, input-order ties, zero steps and bounded drop paths.
Native compilation is warning-free, and two fixture regenerations match
exactly. Source, binary, compiler and pinned library hashes are retained.
Rust explicitly rejects nonfinite/negative variance, invalid configuration and
threshold overflow, rather than reproducing undefined native NaN ordering.

Additional policy tests preserve clean mixed results and weaken one immutable
non-reference SD variance by 1000 cycles squared. The latter excludes that
arc without fixing or mutating it and retains 27/32 moving synthetic FIX
results, all within the existing 3 mm tolerance. Policy synthetic motion is
truth-derived and is separate from whole-chain Doppler accuracy. Marginal DD
means and full covariances are also checked against native submatrices,
including omitted systems and invalid pair-index requests. The five additions
bring Rust to 362 tests; examples, fmt, all-target clippy and release pass.

## PPC result: selection works, accuracy/status gain is absent

All six runs retain caps 1 m/s (default) and 2 m/s (previous sensitivity), both
local reference modes, with/without the existing slip detector and variance
selection on/off: 96 executions. Forty-eight controls exactly reproduce all
53 original columns from the prior lock-qualified archive. All original FLOAT
positions, residuals, covariances, satellites, resets, velocity/covariance and
lock counts remain unchanged on every run. Seven new CSV fields expose source
references/targets, raw DD variances, median, threshold, retained indices and
whether the prefilter changed the problem. The verifier recomputes each mask
from these published values.

The filter acts on 13 distinct `(case, slip profile, GPS time)` combinations:
nine without slip detection in Tokyo run2, and four with it (one Tokyo run1,
three Tokyo run2). Across the two caps and two gate modes these are 52 applied
rows, each dropping one DD target. The maximum variance/upper-median ratio is
17.03315. All cap-1 applications fail motion uncertainty admission. At cap 2,
11 combinations admit motion, producing 22 rows across the two reference
modes; every one still has ratio below 3 and is rejected. Remaining applied
rows are motion-uncertain. `diagnostics.json` retains every affected time,
before/after membership, ratio and decision without reading reference truth.

Only AR membership/pair counts and, on 22 rows, candidate ratio/integers change.
Every solution status, output position, covariance, reference error and FIX
sequence stays identical. The per-condition totals therefore remain:

| Slip detector | Local gates | FIX at cap 1 m/s | FIX at cap 2 m/s |
|---|---|---:|---:|
| off | selected | 2 | 6 |
| off | all references | 0 | 2 |
| on | selected | 1 | 2 |
| on | all references | 0 | 1 |

There are 9,727/9,855 RTK outputs without/with slip detection over 11,733
paired 1 Hz epochs. Original 5 Hz rover coverage is not inferred from these
rates. Existing approximately 0.728 m antenna-assumption FIX error and maximum
two-epoch FIX streak remain. `report.json` scores all RTK/FIX outputs under
antenna and IMU assumptions, including both conflicting Tokyo lever arms.
No origin is selected from these results.

`analysis.json` reuses the frozen whole-native GPS/Galileo default-SPP run
read-only: 10,585 RTK outputs including 144 FIX, 878 SPP and 270 invalid pairs.
No whole-native rerun occurred here. Native/Rust status cross-tabs,
common-success errors and continuity are retained; no FIX times coincide.
Whole-solver settings/filter differences, held-out validation, runtime/memory
and practical RTK/PPP acceptance remain open. Checklist progress stays 18/24
(75%). This prefilter has not improved practical performance on these runs.

The next dependencies are stable partial-set confirmation, controlled
progressive/quality subset testing and motion/measurement quality. The result
does not justify changing gates or defaults to obtain additional FIX points.

## Replay

`predeclared-inputs.json` freezes source, input, binary and configuration.
All commands, return codes, raw stdout/stderr, hashes, native generator
regeneration and standard qualification logs are retained.

```bash
source /workspace/.cloud-setup/activate.sh
cargo build --locked --offline --release --example ppc_rtk_kinematic_fix -j2
python3 -B scripts/run_ppc_variance_subset.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_kinematic_fix \
  --output-dir /tmp/ppc-variance-subset-replay
python3 -B docs/comparison/2026-10-21-variance-subset-ar/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
```

Historical source-hash verifiers require their corresponding source commit.
This verifier uses old captures without rewriting their manifests and seals
all new artifacts with `evidence-sha256.json`.
