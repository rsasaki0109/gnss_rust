# Covariance-standardized satellite screening for moving PPP

The opt-in native sigma-4 screen is connected before the existing NIS gate.
It removes the worst row's whole satellite and rebuilds the ambiguity layout,
wind-up and proposed prior from the last committed state. Numeric state/time
never advance on failure. Direct initialization/velocity fusion is reapplied
to a fresh proposal, once per attempted update. Screened satellites cannot
survive optional rejected-epoch arc retention. Default behavior is unchanged.

The statistic is (S^-1 r)_i/sqrt((S^-1)_ii), S=HPH'+R, with strict threshold
and first-row tie behavior. For the diagonal PPP R this equals the standardized
post-fit residual. It is not raw residual / measurement sigma. Rust screens
the first linearization of each fresh prior and retains its nonlinear IEKF;
native commits one linearized update. Existing gates/noise are unchanged.
At least five code and four phase satellites are still required, so code
residual redundancy can be zero for the five XYZ/clock/ZTD unknowns.

## Component and integration qualification

Two repeated native generations call actual pinned
`worstStandardizedResidualRow` 128 times each: 64 explicit-inverse cases
matching the upstream unit-test construction, and 64 production-LDLT solves.
Both generated files are byte identical across repeats. The first 64-case
fixture remains unchanged after extending the production-path coverage.
Direct kernel decisions match exactly and all 128 statistics meet the existing
relative 1e-10 bound. Rust independently rebuilds each Gaussian S from
P=.01 I, H=I, R=S-P; all 128 reconstructed statistics meet the same bound.
Native source/libraries are unchanged; no whole-native PPP calls are made.

Initially, treating a derived numerical inverse as an independently supplied
covariance rejected small inverse asymmetry (case 8). The low-level native
kernel now checks finite shape/positive diagonals and trusts the supplied
inverse; it does not certify SPD. Actual prior/noise covariance validation in
the checked Gaussian API remains strict and unchanged. A simple Cholesky
solve also missed the statistic bound in ill-conditioned shared-state case 36
by about 1.3e-10. Two compensated residual refinements improve this new solve.
Both inverse-construction paths now pass the original bound; no tolerance was
relaxed. Initial source/logs, inverse-asymmetry measurements and a 90-digit
exact-binary64-input oracle diagnostic are retained. These arithmetic changes
apply only to the new screening calculation, not existing Kalman updates.

Nine tests bring the suite from 437 to 446: three component/boundary/invalid
input tests, five PPP integration tests and one explicit arc-exclusion test.
Both clean/noisy accelerating car traces retain identical baseline states
with no false exclusions. Single/two code-fault cases match the numeric state
of independent omitted-satellite controls and reacquire ambiguity arcs.
Exhaustion preserves full state/time, quality-rejected arcs cannot be retained,
and velocity covariance proves one fusion after a screening retry. Standard
checks and release examples pass. Unit velocity injections are covariance
checks, not measured Doppler or ground truth inputs.

## Fixed PPC experiment

Six runs × two GPS IF pairs × off/native-sigma-4 = 24 profiles / 234,536 rows,
covering all 58,634 raw rover epochs. Fix native head, IndependentEpoch clock,
PVA, GPS/Galileo native-variance-QC SPP/Doppler seeds, calibrated CODE phase,
and InterruptAllArcs. All 12 off controls must reproduce previous 72 columns
after only the previously audited clock-column/source label normalization.
Every source field must match across modes and IF pairs. No reference/IMU
enters the solver; both reference origins and both Tokyo lever arms are scored.
SPP fallback remains a separate status and never counts as PPP. Runtime/RSS
includes load, seed computation and CSV; qualification may overlap these
single-process measurements, so runtime acceptance is not claimed.

The native minimum-count screen still permits only five code rows. A new
Nagoya 3 1C2L counterexample reaches 212.252 m after removing G25 from six
code satellites: five code/four phase remain, ZTD=4.820 m, code RMS=1.381 m,
phase RMS=0.000225 m, NIS=13.194, while its SPP seed is only 7.255 m wrong.
The motion restart and removal both happened in this proposed update. This
case is retained; the native residual screen does not by itself establish
physical position accuracy. Redundancy, weak geometry and ZTD/height coupling
remain priorities. The previously observed kilometre-scale SPP source errors
are also not repaired by a PPP row screen.

All 1C2W modes are compared to both previously frozen native whole-PPP
initialization profiles, with PPP-only and combined-valid metrics separated.
Features/initialization/clock policy differ; these 24 comparisons are diagnostic,
not matched-feature native execution or practical performance acceptance.

## Results

PPP epochs / antenna-hypothesis 3-D RMS (m), each profile's own success set:

| Run | Pair | Off | Native sigma 4 |
|---|---|---:|---:|
| nagoya-run1 | 1C2W | 1603 / 9.951 | 1657 / 9.755 |
| nagoya-run1 | 1C2L | 1211 / 3.613 | 1629 / 4.942 |
| nagoya-run2 | 1C2W | 1418 / 11.982 | 1587 / 17.880 |
| nagoya-run2 | 1C2L | 1670 / 14.148 | 1843 / 19.422 |
| nagoya-run3 | 1C2W | 681 / 14.072 | 785 / 22.927 |
| nagoya-run3 | 1C2L | 587 / 22.515 | 863 / 33.769 |
| tokyo-run1 | 1C2W | 2243 / 7.096 | 2297 / 6.270 |
| tokyo-run1 | 1C2L | 1269 / 6.074 | 1441 / 6.349 |
| tokyo-run2 | 1C2W | 1884 / 4.476 | 2023 / 3.617 |
| tokyo-run2 | 1C2L | 1204 / 5.894 | 1382 / 7.146 |
| tokyo-run3 | 1C2W | 4805 / 5.148 | 4856 / 5.492 |
| tokyo-run3 | 1C2L | 3072 / 5.157 | 3296 / 5.128 |

Aggregate PPP counts are 1C2W **12,634 → 13,205** and 1C2L **9,013 → 10,454**.
There are 3,884 / 8,347 selected screen rows respectively, including attempts
where minimum counts prevent removal. Added PPP output can be inaccurate:
Nagoya 2 1C2W RMS worsens 11.982 → 17.880 m and Nagoya 3 1C2L
22.515 → 33.769 m, with the latter maximum 126.537 → 212.252 m.
Default screening remains off. Common/new/lost errors, all reference origins,
SPP fallback and frozen-native comparisons are retained rather than selecting
a favorable profile or successful subset. No practical parity is accepted.

## Reproduction

```bash
source /workspace/.cloud-setup/activate.sh
python3 -B scripts/build_ppp_screening_fixtures.py --output-dir NEW_NATIVE_DIRECTORY
cargo build --locked --offline --release --examples
python3 -B scripts/run_ppp_standardized_screening.py --output-dir NEW_OUTPUT_DIRECTORY --native-fixtures-dir NEW_NATIVE_DIRECTORY
python3 -B scripts/analyze_ppp_standardized_screening.py --archive NEW_OUTPUT_DIRECTORY
python3 -B docs/comparison/2026-10-21-ppp-standardized-screening/verify.py --verify-local-inputs
```

The verifier checks provenance, all captures, 12 prior controls, all seed fields,
screen diagnostics/counts, qualification and repeated native fixtures. With
local inputs it recomputes both origins, common/new/lost and worst-point/frozen-
native comparisons. Earlier archives retain their own recorded source commits.
