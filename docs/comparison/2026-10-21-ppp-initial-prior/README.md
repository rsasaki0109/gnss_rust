# Optional native initial head prior: full-rate 2×2 PPC experiment

The native compact troposphere climatology is ported with its pressure,
temperature, water-vapour and dry/wet delay outputs. All **312** actual native
latitude/height/day cases match the existing **1e-10** component bound. Six
actual native `initializeFilter` calls verify optional head states/variances.
Both native fixture generators reproduce their files twice, without modifying
upstream source. Four meaningful new tests bring the suite to **410**; fmt,
strict Clippy, release and all release examples pass.

`PppConfig::with_native_precise_head_prior(time)` explicitly opts in. It keeps
the supplied position/clock means, initializes ZTD from climatology, sets
XYZ/clock/ZTD variances to 100/1e8/25 m², and starts constant-velocity mode at
rest with 25 (m/s)². It changes neither defaults nor ambiguity priors, process
noise, observation noise, correction settings, convergence, NIS or postfit
thresholds. It is not the entire native PPP initializer or recovery policy.
See the [contract](../../ppp-initial-prior.md).

The predeclared experiment separates two factors: original/native initial head
variances and constant 2.3 m/climatological ZTD. Four combinations run both
1C/2W and 1C/2L calibrated-phase profiles on all six unmodified 5 Hz rover
files: **48 executions, 469,072 output rows**. All **12 original controls
match every CSV byte** of the previous PPC baseline. The reference never
enters the solver, and all failures/statuses are retained.

| Case | 1C/2W FLOAT, original → native head | 1C/2L FLOAT, original → native head | 1C/2W 3D RMS, m | 1C/2L 3D RMS, m |
| --- | ---: | ---: | ---: | ---: |
| nagoya-run1 | 5 → 7 | 6 → 2 | 11.624 → 10.292 | 9.717 → 2.066 |
| nagoya-run2 | 1 → 7 | 1 → 12 | 7.787 → 4.325 | 5.824 → 4.210 |
| nagoya-run3 | 118 → 118 | 19 → 19 | 8.374 → 6.518 | 6.950 → 3.844 |
| tokyo-run1 | 8 → 26 | 2 → 3 | 5.382 → 3.976 | 5.543 → 4.166 |
| tokyo-run2 | 5 → 5 | 107 → 107 | 1.703 → 1.499 | 3.497 → 3.362 |
| tokyo-run3 | 4 → 49 | 17 → 17 | 4.611 → 5.460 | 3.619 → 2.569 |
| Total | **141 → 212** | **152 → 160** | — | — |

The table uses the antenna-origin hypothesis on successful epochs only.
Both IMU-origin hypotheses, both conflicting Tokyo lever arms, all rejection
reasons and continuity/resource measurements remain in the JSON. It also
separates common-success, newly accepted and lost epochs. For Nagoya run1
1C/2W, the full success-only RMS improves, yet the **two common successes
worsen from 11.765 to 12.922 m**. Nagoya run1 1C/2L has no common successes:
the apparent improvement compares different times. Tokyo run3 1C/2W grows
from four to 49 solutions while full success-only RMS worsens from 4.611 to
5.460 m. These differences prevent a general accuracy-improvement claim.

The added initialization option does not solve poor long-route continuity.
Default behavior stays unchanged. Clock modeling, motion prediction, accepted
observation quality, and recovery after rejection remain the next comparisons;
the initial-prior difference is not claimed to explain every failure. Receiver
calibration remains unknown; satellite PCV and solid Earth tide are absent in
these profiles. Whole-native PPP has not been executed. Descriptive individual
process measurements are not native-controlled runtime/memory acceptance.

Recompute statuses, all 36 comparisons, 12 exact controls, native repeats,
qualification and integrity with
`python3 -B docs/comparison/2026-10-21-ppp-initial-prior/verify.py`.
Use `--verify-local-inputs` to rehash cached inputs and recompute every
reference-origin error. Source hashes describe this milestone; check out its
commit when replaying after later changes.

Practical RTK/PPP acceptance remains unmet. Overall checklist: **18/24 (75%)**.
