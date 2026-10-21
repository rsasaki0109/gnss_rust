# Mixed-clock SPP refits

The explicit QC API and RTK `kinematic_spp_quality` option connect the pinned
native one-pass MAD refit and leave-one-out RAIM. Existing defaults stay off.
The baseline has position and a clock for each actually usable constellation.
All refits start at its final position and keep its separate clock estimates.
Each child uses only baseline-used tracking rows; there is no recursive QC,
tracking substitution, clock prior, covariance floor or truth input.
This explicit stateless initialization differs from native's MAD child, which
re-enters with the processor's stored seed; native's local ISB map also starts
empty. Refitting from the baseline avoids leaking previous candidate state.

With more measurements than position-plus-clock states, apply the native
MAD mask at three sigma, with sigma floor max(1 m, configured code sigma).
Adopt a valid filtered refit immediately, as native does. If the mask is absent
or its child fails, the unchanged baseline proceeds to RAIM. At least two
residual degrees of freedom are required for RAIM. Refit every single-satellite
exclusion in deterministic satellite order. Select the first strict best RMS
only after at least 25 percent and 1 m improvement. Failed candidates never
contribute a finite selector RMS or alter the baseline/clocks.

Every child and the eventual baseline face an explicit GDOP limit of 50 by
default. **The geometry model differs from pinned native SPP:** Rust computes
unit-row geometry with all estimated constellation clocks and the same design
as the code fit. Native computes final DOP with four common-clock columns even
when it estimates inter-system biases. The numerical threshold is unchanged;
the clock geometry is explicit. This is not whole-native SPP equivalence.
Residual RMS and chi-square final gates, robust reweighting and other native
SPP correction modes remain outside this milestone.

The standalone result keeps every original observation, its tracking/missing/
LLI metadata and final disposition, full final clock covariance, pre-QC counts
and residuals, proposed MAD exclusions, every RAIM exclusion/RMS/geometry/error,
and the adopted branch. Navigation/elevation rejections remain in nested
`spp.rejected`; quality exclusions have separate dispositions/evidence. On
terminal failure the API returns an error, not a fabricated valid result.
RTK records a compact successful QC summary and clocks; a failed seed records
its error and uses the existing position fallback. The RTK measurement set
is independently selected: SPP rejection does not remove RTK code/phase rows.
SPP covariance is never fused into RTK.

Native MAD adoption is not an integrity guarantee. On the independent mixed
physical fixture, a 500 m G01 code error causes the MAD mask to remove good G07
and return about 738 m position error. This is retained as a regression showing
masking/swamping and immediate MAD adoption, not counted as accurate positioning.
RAIM-only removes a 100 m G01 error and recovers the existing 3 mm truth bound.
Full operational acceptance needs real-data evidence and further quality work.

Before PPC execution freeze all six cases, the previous native-variance mixed
seed profiles (no detector, slips, slips plus existing 3 m prefit gate), and QC
off/on: 36 runs with 18 exact old controls. Keep unchanged RTK thresholds and
None/GPS/ElevationOnly defaults, record fallback and branch diagnostics with
their observable scope, score both reference origins and both Tokyo lever arms,
and compare common epochs as well as all output sets. No FIX or practical
RTK/PPP acceptance is implied by this FLOAT experiment.
