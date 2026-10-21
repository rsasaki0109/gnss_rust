# GPS/Galileo independent-clock SPP seed

The explicit mixed SPP API selects GPS 1C and Galileo E1 1C/1X code without
requiring phase or Doppler. Valid 1C precedes 1X; all raw rows and alternate/
invalid/unsupported dispositions are retained on success. The state has ECEF
position and one range clock per represented system, with full joint covariance.
Each iteration derives clock groups from usable rows; no clock prior or
covariance floor repairs deficient geometry. Single GPS uses the existing
four-state path exactly. A singleton system adds no position information.

Reuse the existing two-step broadcast transmit time, Sagnac derivative,
Saastamoinen/GPS Klobuchar E1-frequency model, elevation weights and configured
code noise. Group delay uses the pinned native default legacy-primary TGD/BGD;
message-source word and both group-delay fields are retained as provenance.
No source-specific I/NAV/F/NAV rule is guessed. This is not full native SPP
FDE/QC. Nonzero receiver epoch clock offsets are explicitly unsupported, and
nonfinite offsets fail; the existing GPS-only API is unchanged.

RTK can opt in with `kinematic_spp_constellations=GpsGalileo` alongside its
existing optional SPP seed. Default GPS seed/None behavior stays unchanged.
Only the successful position initializes the kinematic epoch; SPP covariance
is not fused and ambiguity states are retained. Failure uses the old position
fallback and reports no invented clocks. Diagnostics name requested scope,
actual satellite count and each available receiver clock.

Before PPC scoring, compare all six runs using GPS seed versus mixed seed,
with no detector, slip detection, and detection plus the existing 3 m prefit
gate: 36 executions including 18 previous GPS-seed controls. Keep every solver
threshold/default unchanged, retain unmatched/invalid epochs and score both
reference origins/Tokyo lever arms. No claim of whole-solver equality follows
from this feature; native SPP input selection and other QC differ. The shared
clock QR is also used by Doppler, whose existing native/independent regression
fixtures must remain passing.
