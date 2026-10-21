# Relative GPS code DSBs and an explicit IF clock gauge

`GpsIfDsbClock` admits relative Bias-SINEX products only under a caller's
explicit `GpsIfClockBinding`. The gauge is GPS IF(C1W,C2W), with native GPS
coefficients `a=f1²/(f1²-f2²)` and `b=-f2²/(f1²-f2²)`, so `a+b=1`.
For source differential biases `d12=b_C1W-b_C2W`, the constraint
`a*b_C1W+b*b_C2W=0` gives `b_C1W=b*d12`, `b_C2W=-a*d12`.
Raw code gets the negative bias as an additive correction:

| Requested code | Additive correction |
|---|---|
| C1W | -b*d12 |
| C2W | a*d12 |
| C1C | -b*d12 - DSB(C1C,C1W) |
| C2L | a*d12 + DSB(C2W,C2L) |

The required edges are fixed, with reversed source orientation supported.
Missing or expired edges fail explicitly; overlapping forward/reverse edges
are ambiguous even if values agree. Receiver rows and absolute OSBs are not
used. Alternate graph paths are not silently substituted or averaged.
Other systems/tracking requests remain unsupported. A per-product index avoids
scanning every DSB for each epoch. Product half-open intervals are checked
before materializing an epoch-only inclusive correction snapshot.

This is a mathematically chosen code gauge, **not an absolute OSB product**.
The binding is the caller's clock-datum assertion; neither BIA nor an arbitrary
precise CLK header proves its compatibility. Phase corrections stay absent.
In particular this API does not invent zero phase OSBs or admit PPP phase
measurements that still lack their own explicit datum. All original source
rows, units, SVN, orientation, sigma, interval labels and signed coefficients
are retained. A source standard-deviation upper bound is
`sum(abs(coefficient)*source_sigma)`, permitting unknown source correlations.
It is not an independent-error RMS or covariance, and is not fused into noise.

An independent C++ Eigen full-pivot LU solves four simultaneous bias/gauge
equations on 24 cases (96 corrections), using pinned native IF constants and
its additive correction sign. It is **not** an upstream DSB gauge-conversion
function. Rust tests reverse every edge, mix ns/metres, remove injected code
biases, preserve no-phase admission and test validity, ambiguity, graph-path,
receiver exclusion and uncertainty bounds. During development this independent
solve caught a sign error in the initial datum term; the corrected formula
above passes without loosening tolerance.

Before evaluating any positions, audit all nine already acquired CAS days at
GPS noon, all 32 GPS identifiers and four supported requests (1,152 statuses).
Also read all six PPC rover files at their original rate and count raw complete
1C/2W and 1C/2L code/phase pairs, unflagged pairs and code-DSB availability per
epoch. Retain events, unsupported epoch offsets, missing LLI counts and every
DSB failure category. No best pair is selected using reference truth; no phase
correction or positioning is executed. Clock/product/antenna compatibility,
kinematic PPP, full corrections and native PPP performance remain pending.
