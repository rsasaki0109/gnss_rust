# Initial ZTD uncertainty only: PPC and independent car diagnostics

This parameter study changes one initial covariance entry. The library, ZTD
mean/process noise/design column, measurement noise/gates, SPP/Doppler seeds,
PVA/clock/phase/arc policies and all corrections are unchanged. ZTD remains
estimated. Example flag `-ztd-sigma SIGMA_M` is explicit; defaults are unchanged.
This is not a full native broadcast or MADOCA processing profile.

The pinned native `initialTroposphereVariance` returns 25 m² for its standard
precise path and 0.12²=0.0144 m² for per-frequency MADOCA. Its configured
broadcast default is 0.36 m². Twelve actual helper combinations are generated
twice, with identical bytes. Only these prior values are borrowed; CODE SP3/
CLK/IF processing remains in use. `prior-origin.json` identifies source headers,
pinned commit/archive, the configured native value and their hashes.

## Contract and qualification

Fix PVA, native-head climatological mean, IndependentEpoch clock, GPS/Galileo
native-variance-QC SPP/Doppler seeds, calibrated phase and InterruptAllArcs.
Standardized residual screening remains off in all modes. Compare initial
variances 25 / 0.36 / 0.0144 m² across six cases and two GPS IF pairs:
36 profiles / 351,804 rows. All 58,634 raw rover epochs are retained in each
pair/profile aggregate. All 12 baseline captures must reproduce previous
72-column CSV bytes exactly; all 23 seed fields must match across priors/pairs.
Both reference origins and both Tokyo arms, common/new/lost PPP epochs,
separate SPP fallback and all worst points are retained. No reference or IMU
enters either solver. No practical or runtime acceptance is claimed.

The library and all old numerical test inputs match base commit 0d6318f.
Current full tests still pass 446, with fmt/clippy/release and both changed
release examples passing. Twelve invalid sigma checks cover zero, negative,
nonfinite, overflow and squared-variance underflow before any dataset solve.
The flags log sigma/variance/mean/process noise explicitly.

## Independent synthetic prior-bias diagnostic

Three accelerating car traces × three priors = nine runs. The synthetic
baseline retains mean ZTD 2.3 m / variance 100 m², whereas PPC baseline has
its independently computed climatic mean / variance 25 m². Synthetic truth
is 2.4 m ZTD; it is not used to correct the initial mean. Old baseline CSVs
must match and old accuracy bounds remain unchanged. All nine emit 64 PPP
solutions, but four tighter-prior clean/unrounded profiles fail the existing
1 cm target. These are retained failures, not widened tolerances.

| Car input | Original variance 100 m² | Sigma 0.6 m | Sigma 0.12 m |
|---|---:|---:|---:|
| Clean RMS / max (m) | 0.000081 / 0.000376 | 0.014900 / 0.084423 | 0.104338 / 0.261822 |
| Noisy RMS / max (m) | 0.267524 / 0.809500 | 0.248793 / 0.753569 | 0.171618 / 0.397078 |
| Unrounded RMS / max (m) | 0.000069 / 0.000416 | 0.014903 / 0.084449 | 0.104339 / 0.261817 |

Tighter priors reduce noisy error while biasing clean data when the mean is
wrong. Unrounded observables retain the same prior-bias effect. Complete
position/velocity/acceleration metrics and all original bound failures remain
in the report. More successful epochs or a favorable noise case cannot be
used to claim that the tighter prior preserves independent synthetic accuracy.

## PPC results

PPP epochs / antenna-hypothesis 3-D RMS (m), each profile's own successful set:

| Run | Pair | Variance 25 m² | Sigma 0.6 m | Sigma 0.12 m |
|---|---|---:|---:|---:|
| nagoya-run1 | 1C2W | 1603 / 9.951 | 1588 / 9.892 | 1588 / 9.829 |
| nagoya-run1 | 1C2L | 1211 / 3.613 | 1240 / 3.654 | 1218 / 3.372 |
| nagoya-run2 | 1C2W | 1418 / 11.982 | 1875 / 17.040 | 1984 / 16.803 |
| nagoya-run2 | 1C2L | 1670 / 14.148 | 2089 / 12.237 | 2032 / 11.152 |
| nagoya-run3 | 1C2W | 681 / 14.072 | 675 / 11.393 | 635 / 11.142 |
| nagoya-run3 | 1C2L | 587 / 22.515 | 592 / 22.271 | 609 / 22.458 |
| tokyo-run1 | 1C2W | 2243 / 7.096 | 2269 / 6.780 | 2269 / 6.765 |
| tokyo-run1 | 1C2L | 1269 / 6.074 | 1257 / 5.748 | 1268 / 5.048 |
| tokyo-run2 | 1C2W | 1884 / 4.476 | 1878 / 4.407 | 1870 / 4.154 |
| tokyo-run2 | 1C2L | 1204 / 5.894 | 1204 / 5.893 | 1216 / 6.524 |
| tokyo-run3 | 1C2W | 4805 / 5.148 | 4825 / 4.948 | 4811 / 4.769 |
| tokyo-run3 | 1C2L | 3072 / 5.157 | 3093 / 4.873 | 3102 / 4.633 |

PPP counts are 1C2W **12,634 → 13,110 / 13,157** and 1C2L
**9,013 → 9,475 / 9,445** (sigma 0.6 / 0.12 m). Success sets differ; common,
new and lost epoch errors remain separate in the report. Initial uncertainty
alone gives mixed changes and does not remove the large-error problem.
Defaults remain unchanged and practical performance parity is unaccepted.

## Evaluation limits

The initial covariance alone does not constrain later estimated ZTD to the
climatic mean. Motion-only recovery retains its existing ZTD marginal. Five
code rows can still fit five XYZ/clock/ZTD unknowns without code redundancy;
phase rows with new independent ambiguities do not necessarily validate that
fit. Weak geometry, persistent code bias, ZTD/height coupling and source SPP
outliers remain priorities. Combined coverage counts include separately
labeled SPP, with its previously recorded kilometre-scale outliers.

Post-hoc diagnostics include 270 worst-origin records and 36 comparisons with
the two frozen native whole-PPP initialization profiles for 1C2W. The effective
clock policies, seeds/QC and initial priors differ; native broadcast/MADOCA
whole solvers were not executed here. Process wall/CPU/RSS measurements include
load/seed computation/CSV and may overlap qualification; they are descriptive.

## Reproduction

```bash
source /workspace/.cloud-setup/activate.sh
python3 -B scripts/build_ppp_ztd_prior_fixtures.py --output-dir NEW_NATIVE_DIRECTORY
cargo build --locked --offline --release --example ppc_ppp --example ppp_kinematic
python3 -B scripts/run_ppp_ztd_prior.py --output-dir NEW_OUTPUT_DIRECTORY --native-fixtures-dir NEW_NATIVE_DIRECTORY
python3 -B scripts/analyze_ppp_ztd_prior.py --archive NEW_OUTPUT_DIRECTORY
python3 -B docs/comparison/2026-10-21-ppp-ztd-prior/verify.py --verify-local-inputs
```

The verifier checks every capture, source/input/old-library/native provenance,
qualification, 12 exact PPC and three exact synthetic controls, sigma/mean/
process-noise diagnostics and all status/seed fields. With local inputs it
recomputes every origin/common/worst/frozen-native metric. Historical archives
use their own recorded source commits.

The first local-input verifier encountered a Python duplicate-keyword error
when merging two records for the same private native header. Its source/log
are retained in `initial-verifier`; the merge now checks matching duplicate
hashes explicitly. No solver, capture or evaluation metric changed.
