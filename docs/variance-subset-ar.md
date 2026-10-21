# Relative-variance preselection

The opt-in `for_kinematic_variance_subset` first forms the existing lock-qualified
set, then uses native `filterPairsByRelativeVariance` defaults: upper median of
DD variances, multiplier 10, minimum threshold 1e-4 cycles squared, minimum four
pairs (or the configured stricter policy minimum). If too few pass, retain the
unfiltered set as the native resolver does. Otherwise search one surviving set
only. No covariance diagonal floor, repeated integer search, ratio-dependent
subset choice, integer holding or feedback to FLOAT is introduced.

The original within-system references remain for surviving target pairs;
systems with no surviving target disappear. Marginalize omitted SD states with
their full principal covariance and preserve position/SD and cross-system
correlations. Existing integer labels, confirmation, motion and local gates
apply to the resulting membership. The prefilter depends on its input DD basis;
all-reference candidate validation does not make subset search basis invariant.

The native progressive variance-drop helper is also ported as a component,
including first-index tie breaking, minimum size and step limits. These paths
are not searched by the FIX policy yet. Multiple candidate testing and FFRT
remain separate work. Nonfinite/negative variance and overflowing threshold
are explicit errors; native undefined NaN ordering is not reproduced.

Before PPC scoring, retain default motion cap 1 m/s and the previous 2 m/s
sensitivity, both local reference modes, with/without slip detection and with/
without variance filtering on every case: 96 runs including 48 unchanged
lock-qualified controls. All local and motion thresholds remain unchanged.
Evaluate both origin hypotheses and Tokyo lever arms; no truth-based subset
choice, changed default or practical parity acceptance follows from this run.
