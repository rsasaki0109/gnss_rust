# Optional phase continuity across rejected PPP updates

Default `InterruptAllArcs` retains its original behavior. The explicit
`RetainContinuousArcs` policy is an experimental recovery choice, not a claim
of whole-native lifecycle equivalence. Native runs slip detection before its
update and does not globally discard all ambiguity states on every failure;
Rust still preserves transactional numeric state and uses last **accepted**
GF/MW anchors, rather than native per-observation history.

For measurement/geometry rejection only, retain an already active arc if the
current prepared proposal has the same ID and no reset reason. Existing time
gap, missing phase/code, half-cycle/LLI, power/event, GF/MW and signal/datum/
correction changes still interrupt it. New proposals cannot become accepted
through rejection. Global configuration/product/model/arithmetic errors still
break all arcs. Events keep their old behavior.

Retaining an arc does not accept this epoch: position, clock, ZTD, velocity,
full covariance, numeric state time, phase accepted time/count and GF/MW values
stay exactly at their last accepted values. Wind-up anchors remain unchanged
only for retained arcs. The next real update predicts over the full gap from
the last numeric commit. Rejected intermediate measurements are never replayed.
The usual maximum gap is measured from the last accepted phase, so repeatedly
retaining does not extend arc lifetime. IDs consumed by rejected new proposals
are never reused. Published solutions record the selected policy.

Tests must cover retention with unchanged time/count and diagnostic anchors,
intrinsic resets during a rejected update, cumulative gap expiry, failed
initialization, no pending proposal, old default interruption, numeric/wind-up
transactionality, successful recovery and event/global invalidation. Existing
synthetic accuracy bounds and all filter gates stay unchanged.

Real evaluation is predeclared on all six PPC runs at original 5 Hz, both
GPS 1C/2W and 1C/2L, original/native head priors, original/independent clock,
retention off/on: 96 runs, 938,144 rows. The 48 off controls must reproduce
their previous non-timing CSV bytes. Retention-on outputs add explicit retained
arc counts/IDs; rejected epochs remain INVALID. Score all successes and common/
new/lost epochs under both origins and both Tokyo lever arms. Do not select a
winner by truth, remove bad epochs or interpret higher counts as acceptance.
Runtime/RSS remain descriptive. No whole-native rerun is required for this
isolated policy test; the preceding whole-native benchmark remains frozen at
commit c714cd9.
