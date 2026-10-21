# Lock-filtered partial AR

The opt-in `RtkFixPolicy::for_kinematic_locked_subset` searches exactly one
subset: satellites with the unchanged five-epoch lock requirement, omitting
systems with fewer than two qualified satellites. Choose the highest-elevation
qualified reference in each retained system, keeping canonical satellite order
on equal elevations. Require at least four within-system DD pairs in total.
No ratio/variance-driven subset iteration or integer holding is implemented.

Project the immutable posterior onto position and retained SD ambiguities by
extracting their full principal covariance block, including all cross-system
and position/ambiguity correlations. This is Gaussian marginalization of the
omitted states, not conditioning them to zero or discarding their influence
on the existing FLOAT posterior. Rebuild DD and residuals from retained raw
observations. Excluded raw observations remain in FLOAT and diagnostics; local
integer residual gates apply only to the declared subset. Its uncertainty
still includes the original FLOAT update's use of all observations.

Use the existing ratio, correction, covariance, residual, baseline, motion,
minimum-lock and confirmation thresholds. Repeated integer labels include
the selected satellite identities and are normalized independently per system.
When the eligible subset changes, confirmation restarts. A reset of a formerly
qualified arc clears history; a newly introduced or still unqualified arc
does not erase a continuing qualified subset. Motion unavailable/uncertain,
global reset and gaps keep the prior kinematic behavior. No partial FIX feeds
back into FLOAT, and each epoch uses fresh integer search.

Before PPC scoring, retain the previous motion cap 1 m/s plus 2 m/s sensitivity,
with selected/all-reference gates, with/without slip detection and full/locked
subsets on every run: 96 executions including 48 unchanged full-set controls.
No numeric threshold changes or truth-based subset choices are allowed. This
is the first lock-filtered partial AR path; native quality/variance fallback,
FFRT, holding and practical RTK/PPP performance acceptance remain pending.
