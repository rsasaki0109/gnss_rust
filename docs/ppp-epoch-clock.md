# Optional independent receiver-clock prior per epoch

`PppClockModel::IndependentEpoch` opts in to a fresh clock prior at each PPP
measurement update. Its variance is the configured initial clock variance
(`1e8 m²` in these profiles), and its centre is the last committed clock mean.
The covariance's clock row/column are cleared, and its clock diagonal is
replaced. All other means and the full position/velocity/ZTD/ambiguity marginal
covariance are retained. This is a finite weak prior, not an exactly free
clock parameter. The old mean can have a small regularizing effect.

Clock random-walk process noise is unused in this model. Motion, ambiguity,
troposphere, observation noise, gates and correction/datum/LLI handling remain
unchanged. No observed SPP or Doppler clock is inserted, and no reference truth
is fused. The old random-walk model remains the default. Output records the
clock model. Prediction operates on a copy: failed measurements and events
do not replace the committed clock or any other numerical state.

The covariance replacement is compared with 16 actual calls to the pinned
native `PPPProcessor::reinitializeScalarState` on dense SPD matrices. This
proves the scalar operation, not native precise-PPP clock policy. Native
broadcast/PVA paths can reset clocks using SPP seeds; this option applies a
separate finite prior in precise PPP without those observed seeds. Do not
label it whole-native parity or port a broadcast policy invisibly.

Synthetic controls add a common receiver-clock range term to both code and
carrier on every frequency, preserving the independent moving fixture's
geometry and raw flags. The first eight epochs are unchanged, then offsets
repeat `[30000, -20000, 50000, 0]` metres. They are declared hardware-clock
stress inputs, not measured PPC jumps. Both clean and noisy 64-epoch traces
must preserve positions and estimated clock differences versus the same
independent-clock control within 1 mm, retain existing noisy 1 m accuracy,
and recover after a separate code outlier/event. Original rounded clean
full-trace 2 mm failures remain in their original tests; no old bound is
relaxed by the new stress test's 1 cm absolute bound.

PPC comparisons predeclare original/native head priors crossed with original
random-walk/independent clocks, calibrated CODE phase, and both GPS IF tracking
pairs, on all six raw 5 Hz rover files. Retain all failures and statuses, both
reference origins, common/new/lost success errors and continuity. Original
clock controls must match every previous CSV byte. Neither clock variance
nor gates are selected against reference errors. Missing receiver calibration,
satellite PCV, tide, full native comparison and practical acceptance remain
separate dependencies.
