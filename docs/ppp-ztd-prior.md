# Initial ZTD uncertainty parameter study

The public PPP configuration already accepts an independent initial ZTD
variance. This experiment changes that single initial covariance entry through
an explicit example flag `-ztd-sigma SIGMA_M`; no library processing changes.
ZTD remains estimated with its existing design column, mean, process noise and
full covariance. A motion restart retains its marginal as before. Existing
measurement noise, NIS/postfit/physical gates, SPP/Doppler seeds, phase arcs,
clock/PVA models and corrections are unchanged. Defaults are unchanged.

Choose the values from pinned native `initialTroposphereVariance`, before
running: the precise-mode baseline is 25 m² (sigma 5 m), configured broadcast
variance is 0.36 m² (sigma 0.6 m), and the MADOCA per-frequency path uses
0.12²=0.0144 m². Only these initial uncertainties are borrowed; the Rust
CODE ionosphere-free model is not relabeled broadcast or MADOCA PPP. The
helper is called for 12 combinations, twice; generated tables must agree
byte-for-byte. Flags accept any finite positive sigma whose square remains
finite and positive, and diagnostics label sigma/variance/mean/process noise.

PPC predeclared grid: six cases x two IF pairs x baseline/0.6 m/0.12 m =
36 profiles / 351,804 rows, original 58,634 rover epochs per IF pair/profile.
Fix native head mean, PVA, IndependentEpoch clock, GPS/Galileo native-QC
SPP/Doppler seeds, calibrated CODE phase and InterruptAllArcs. Standardized
residual satellite screening stays off in every profile. All 12 baseline
72-column CSV controls must match the previous qualified off captures exactly.
All 23 SPP/Doppler source fields must match across priors and IF pairs.
Do not use reference or IMU in either solver. Score both reference origins,
both Tokyo arms, common/new/lost PPP, separate SPP fallback and all worst points.
No new whole-native PPP or runtime/practical parity acceptance is implied.

Run independent 64-epoch accelerating car clean/noisy/unrounded traces for all
three priors (nine runs). Synthetic baseline retains the original ZTD mean
2.3 m and variance 100 m²; the independent true ZTD is 2.4 m. New sigma flags
change only its variance and must not correct this mean using truth. Preserve
all position/velocity errors and existing 1 cm clean / 1 m noisy and 1 m/s
after-eight velocity criteria without widening them. Require the three old
baseline CSVs and bounds to match; retain any failures of the tighter priors
as results rather than declaring that the priors preserve synthetic accuracy.
PPC baseline instead keeps its independently computed climatological mean
and 25 m² variance. These distinct baseline conditions must be explicit.

Library source/tests remain unchanged from the 446-test screening qualification.
Re-use that qualified evidence with exact source hashes, then run fmt, all-target
clippy and release builds of both changed examples. Check invalid/underflow/
overflow sigma CLI inputs without running a dataset solver. Archive all source/
binary/input/native hashes before execution and all status/CSV/log captures.
A smaller prior can introduce bias if the initial climatic ZTD mean is wrong;
more PPP solutions or smaller residuals alone do not establish accuracy.
