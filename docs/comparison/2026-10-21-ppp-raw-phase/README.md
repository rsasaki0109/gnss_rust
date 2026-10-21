# Explicit raw-carrier FLOAT and PPC precise-product readiness

Four 64-epoch synthetic comparisons (clean/noisy, code-only OSBs/relative
DSBs) retain all epochs as PPP FLOAT. Maximum position change versus
calibrated phase is **0.657 micrometres**; the ambiguity shift agrees with
independent native-constant hardware-bias expectations within **0.056
micrometres**. The bound was 1 micrometre. No phase OSB is manufactured.
Strict calibrated phase remains the default. All 11 previous moving profiles
and 10 static profiles reproduce their previous **complete CSV bytes**.

The independent C++ generator reproduces both fixtures twice without changing
native source. This is a synthetic gauge check, not a whole-native PPP run.
The first comparison harness mistakenly required nine phase satellites on
every epoch; eight epochs actually have one code-only satellite. Its failure,
source and captures are retained under `initial-harness-failure/`. The harness
now counts actual supplied phase corrections. Solver and tolerance unchanged.

Matching **CODE rapid 27/27 orbit/CLK/absolute-OSB products** were downloaded
for nine central/padding days. Rust reads all nine days: GPS IF(C1W,C2W) and
IGS20 APC are explicitly declared; the sampled 32 GPS states/clock pairs are
available on each day. The maximum IF code-gauge rounding residual is under
0.064 mm. At that sample time, 1C/2W code and phase OSBs cover all 32 GPS
satellites; 1C/2L code OSBs cover 25, and phase OSBs 32. These are product
admissions, not usable rover epochs or PPP solutions. Original RINEX CLK
2.00 headers and non-UTF-8 BIA comment bytes remain unchanged in the cache.

Canonical IGS and BKG ANTEX files have identical SHA-256. The verbatim GPS
excerpt contains 116 satellite blocks with 810 explicit exclusions, and is
retained here compressed. At each of the three central dates all 32 exact
SVN/validity calibrations are available. CAS G01 says G063, conflicting with
IGS20 G049 at all three dates; CODE agrees with G049. The conflict remains
recorded and is not silently repaired. Receiver calibration is unknown,
satellite PCV is not applied, and real PPC PPP has not yet run at this milestone.

Eight meaningful added tests bring the suite to **406**, with fmt, strict
Clippy, release and all release examples passing. The technical contract is
[raw-phase-ppp.md](../../raw-phase-ppp.md). Recompute this evidence with
`python3 -B docs/comparison/2026-10-21-ppp-raw-phase/verify.py`; add
`--verify-local-inputs` to rehash cached public products. Source hashes describe
this milestone; use its commit when replaying after further development.

Practical RTK/PPP parity remains unaccepted; the broad checklist stays
**18/24 (75%)**. Next: run moving PPC PPP with the matching CODE/IGS20 inputs
and retain every success, rejection, both reference-origin hypotheses and
the missing correction assumptions before comparing the whole native solver.
