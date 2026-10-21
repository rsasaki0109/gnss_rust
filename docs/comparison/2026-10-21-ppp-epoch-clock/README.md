# Independent clock prior per epoch: positive stress test, mixed PPC results

`PppClockModel::IndependentEpoch` explicitly replaces the clock covariance
row/column at each measurement update while retaining the full remaining
position/velocity/ZTD/ambiguity marginal. Its finite prior is centred on the
last committed clock mean, with configured initial variance (1e8 m² here).
It inserts no SPP/Doppler clock measurement. Defaults, ambiguity priors and
all residual/NIS/convergence thresholds remain unchanged. See the
[model contract](../../ppp-epoch-clock.md).

Sixteen actual native `reinitializeScalarState` calls match the covariance
operation on dense SPD matrices; their fixture reproduces twice. This is a
scalar component check, not native precise-PPP clock policy or whole-solver
parity. Tests cover common receiver-clock jumps, clean/noisy moving accuracy,
and transactional rejection/event recovery. **413 tests** and all standard
checks pass.

Eight synthetic profiles retain all 64 input rows, including failures. The declared jump pattern
starts after eight epochs, repeating 30,000/−20,000/50,000/0 metres on both
code and carrier. The original random-walk model retains **22 solutions and
42 rejections**; the independent clock retains **64 solutions** for clean
and noisy data. Against its unshifted control, maximum position difference is
**0.0121 mm** and clock-offset error **0.0141 mm**, within the predeclared 1 mm
equivalence bound. Original rounded full-trace 2 mm failures remain visible;
this stress result does not waive them or establish real-data accuracy.

PPC crosses original/native initial heads with original/independent clocks,
for both calibrated CODE GPS IF pairs: **48 runs, 469,072 raw 5 Hz output
rows**. All **24 original-clock controls match their complete previous CSV
bytes**. Every case/clock/head/pair processes all original rover epochs.

| Initial head | 1C/2W FLOAT, random walk → independent | 1C/2L FLOAT, random walk → independent |
| --- | ---: | ---: |
| Original | 141 → 418 | 152 → 406 |
| Native head | 212 → 442 | 160 → 392 |

Each cell uses **58,634 epochs**, and tracking-pair counts overlap. Counts
remain below 1%; long-route continuity is unmet. Better counts do not imply
better accuracy: original-head Nagoya run2 1C/2W success-only antenna-origin
3D RMS grows **7.787 → 87.781 m**, with maximum **149.051 m**; native-head
Nagoya run2 1C/2L RMS grows **4.210 → 23.012 m**, maximum **40.374 m**.
Other conditions improve, including native-head Tokyo run3 1C/2L RMS
**2.569 → 1.734 m**, with 17 → 186 solutions. No profile is chosen as a
general improvement and the default clock remains random walk.

Both reference-origin hypotheses, both Tokyo lever arms, all terminal
rejection reasons, continuity and common/new/lost-success errors remain in
the report. `worst-points.json` adds 60 post-hoc worst-origin diagnostic
records. The worst Nagoya run2 point is newly accepted after **822.6 s**
without a commit, has **5 code/5 phase satellites**, code RMS **0.705 m**,
phase RMS about **4.4e-9 m**, covariance trace **471.21 m²**, ZTD **0.607 m**,
and antenna/IMU-hypothesis error about **149.05/150.05 m**. Estimated speed
is **1.056 m/s**, reference ENU speed **0.068 m/s** (velocity origin remains
unknown). Those facts expose a poor recovery result; they do not prove one
isolated cause. The diagnostic reads retained outputs/reference only and
does not select a solver profile or modify a threshold.

The next priority is whole-native PPP on matched inputs and correction
settings, then observation redundancy/geometry, motion prediction and
ambiguity continuity/recovery differences. Receiver calibration is unknown,
satellite PCV/tide are absent in these profiles. Native-controlled timing,
memory and practical RTK/PPP acceptance remain pending.

Verify every profile, native repeat, qualification and archive integrity with
`python3 -B docs/comparison/2026-10-21-ppp-epoch-clock/verify.py`.
Add `--verify-local-inputs` to rehash cached public inputs and recompute
all origin-hypothesis errors and worst-point diagnostics. Source hashes refer
to this milestone; use its commit to replay after further development.

Overall completion remains **18/24 (75%)**; practical parity is unaccepted.
