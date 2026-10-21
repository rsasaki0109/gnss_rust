# Lock-filtered partial AR on PPC

`RtkFixPolicy::for_kinematic_locked_subset` adds one opt-in AR subset of
lock-qualified satellites. It does not change the default GPS static policy
or the existing kinematic full-set constructor. Every satellite still needs
five lock epochs, each retained system needs two satellites, and their union
needs at least four within-system DD pairs. The highest-elevation qualified
reference is selected separately per system, with canonical satellite order
breaking ties. No variance/ratio-driven dropping, FFRT or holding is added.

The immutable posterior is marginalized to position and retained SD states by
taking the full principal covariance block. Position/ambiguity and cross-system
correlations remain. Omitted integers are never fixed or assigned zero
uncertainty; their observations remain in FLOAT and diagnostics. Local integer
residual checks use the explicitly declared subset. Repeated labels include
its satellite identities, so membership changes restart confirmation. Resets
of formerly qualified arcs clear history; unqualified/new arcs do not interrupt
an otherwise continuing qualified subset. Existing motion, ratio, lock,
correction, residual, covariance, minimum-pair and confirmation thresholds stay
unchanged. No fixed position/integer or velocity feeds back into FLOAT.

## Component and synthetic verification

The new independent generator calls pinned native
`buildDoubleDifferencePairsForSystem`, `extractSubset` and `solveFixedHeadState`.
Twenty-four cases include full/partial/empty selection, five-lock boundaries,
unqualified highest references, singleton/absent systems, equal elevations and
varying maturity. Full SPD covariance contains position/SD and cross-system
correlations; fixed covariance is computed separately by Eigen Schur complement.
Rust uses the existing 1e-11 × max(1,abs(expected)) component tolerance without
a covariance floor. Warning-free native compilation, library/source/binary
hashes and two exact fixture regenerations are retained in `native-generator.json`.
These checks cover the native selection and conditioning components, not its
entire subset-search policy.

Four added moving synthetic integration tests verify unchanged full-set
behavior when all arcs qualify, continued FIX while a new arc matures, fresh
confirmation when membership changes, reselection of a slipping reference,
and singleton omission without cross-system DD. Every synthetic fixed position
keeps the existing 3 mm truth tolerance and the FLOAT posterior remains
unchanged. Policy synthetic tests supply truth-derived motion, separately from
real Doppler evaluation. The five additions bring the Rust total to 357 tests;
examples, fmt, all-target clippy and release build also pass.

## All six PPC runs

Before scoring, motion caps 1 m/s (existing default) and 2 m/s (previous
sensitivity condition), selected/all-reference gates, with/without slip
detection and full/locked subsets were retained: 96 executions, including 48
full-set controls. All original 49 columns of those controls exactly match the
previous archive. On all executions, original FLOAT positions, residuals,
covariances, satellites, resets and Doppler values/covariances match. Four new
CSV columns expose AR satellites, AR references, AR DD count and all original
satellite lock counts. The verifier checks subset eligibility, singleton
exclusion, pair counts and complete reference-audit size, and confirms the
same lock ledgers across configurations.

| Slip detector | Local gates | Cap m/s | Full-set FIX | Locked-subset FIX |
|---|---|---:|---:|---:|
| off | selected | 1 | 2 | 2 |
| off | selected | 2 | 6 | 6 |
| off | all references | 1 | 0 | 0 |
| off | all references | 2 | 2 | 2 |
| on | selected | 1 | 1 | 1 |
| on | selected | 2 | 1 | 2 |
| on | all references | 1 | 0 | 0 |
| on | all references | 2 | 0 | 1 |

Every condition retains 9,727 RTK outputs without slip detection or 9,855 with
it, over 11,733 paired 1 Hz epochs. The original rover has 58,634 epochs at
5 Hz; unmatched/invalid epochs retain those statuses. The longest FIX streak
is still only two paired epochs. Across the six runs, 5,436/7,890 valid FLOAT
snapshots have fewer eligible AR satellites than FLOAT satellites without/with
slip detection, including snapshots with no usable subset. Eligibility alone
does not establish a searchable or accepted candidate.

At cap 2 m/s with slip detection, insufficient-lock decisions fall from
5,511 to 3,395, while ratio rejection rises from 1,084 to 2,029 and selected
confirming decisions rise from 5 to 31. Without slip detection, insufficient
lock falls 3,648→2,942 and ratio rejection rises 2,163→2,667. Most newly
eligible problems therefore still do not yield a valid repeated integer set.

The additional accepted point is Tokyo run1, GPS week 2324 TOW 187767. The
FLOAT update uses 11 satellites; G05/G29/E34 have four lock epochs and are
omitted from AR. The eight retained satellites give six DD pairs, ratio
7.17085 and 0.055658 m motion residual; all 16 reference combinations pass.
Its 3D reference error is 0.051338 m under the antenna assumption and
0.675148/1.543013 m under the two Tokyo IMU/lever-arm assumptions. The
reference origin remains unknown. Existing Tokyo run1 FIX with 0.728253 m
antenna error also remains, so no false-FIX or absolute accuracy acceptance
follows from this gain. `report.json` retains all RTK and FIX-only errors for
every hypothesis, rather than selecting a favorable origin.

`analysis.json` reuses the immutable prior whole-native GPS/Galileo default-SPP
run read-only; no whole native rerun occurred here. It has 10,585 RTK outputs,
including 144 FIX, plus 878 SPP and 270 invalid pairs. No Rust/native FIX times
overlap. Per-case status cross-tabs, common-RTK errors for all origin
hypotheses, unmatched success counts and FIX runs are retained. Filter/seed/QC
and full subset-search differences remain. Held-out testing and runtime/memory
acceptance have not been performed. Practical RTK/PPP parity is not accepted;
the overall checklist remains 18/24 (75%).

Next dependencies are quality/variance-based subset search, reliable integer
validation/confirmation and motion uncertainty, followed by whole-solver
comparison. Raising caps or relaxing gates is not justified by this result.

## Replay

Source/input/binary/configuration hashes are in `predeclared-inputs.json`; each
execution retains commands, return codes, raw stdout/stderr and hashes.
Native generator details and standard check logs are also retained.

```bash
source /workspace/.cloud-setup/activate.sh
cargo build --locked --offline --release --example ppc_rtk_kinematic_fix -j2
python3 -B scripts/run_ppc_locked_subset.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_kinematic_fix \
  --output-dir /tmp/ppc-locked-subset-replay
python3 -B docs/comparison/2026-10-21-locked-subset-ar/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
```

Run a historical archive's source-hash verifier at its recorded source commit.
The new verifier reads previous captured results as immutable controls without
rewriting their source manifests. `evidence-sha256.json` seals this archive.
