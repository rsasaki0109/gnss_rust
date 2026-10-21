# Native SPP measurement variance and QC utilities

The explicit `NativeSnrAtmosphere` choice uses the pinned `spp_utils` variance:
code sigma/elevation, CN0 penalty below 45 dB-Hz clamped to 1..4, corrected
ionosphere sigma 0.5 times correction clamped to 0.25..8 m, and corrected
troposphere sigma 0.1 times correction clamped to 0.05..3 m. Uncorrected
ionosphere/troposphere use 2.5/sin(elevation) and 0.5/sin(elevation). The
native code sigma minimum 0.1 m and sine minimum 0.1 are preserved. These are
measurement-noise bounds, not a state-covariance floor or nuisance-clock prior.
Missing/nonpositive/nonfinite CN0 incurs no CN0 penalty, matching the utility;
raw values and missing distinctions remain in observation evidence.

Orbit, satellite clock, TGD/BGD, atmospheric correction means, signal selection
and QR are unchanged. Only row variance changes. The existing API defaults to
elevation-only; an explicit mixed SPP overload and RTK seed config opt in.
Successful model evidence exposes the applied row variance. SPP covariance is
never fused into RTK. Invalid active corrections or overflowing variance are
explicit errors instead of native fallback to an optimistic base variance.

Native robust MAD residual masks and leave-one-out RMS selection are also
ported as components. All-true mask on inadequate redundancy means no filter,
not a valid solution. With pinned caller defaults, FDE selection requires 25%/1 m improvements
and keeps the first strict best on ties. Negative candidate RMS is treated as
invalid. Neither helper is connected to repeated SPP refits in this milestone.
Geometry QC, full FDE and complete upstream SPP behavior remain separate work.

Before scoring PPC, keep the previous three mixed-seed profiles (no detector,
slips, slips plus existing 3 m prefit gate), variance on/off on all six cases:
36 runs with 18 unchanged controls. No existing gate/default or truth-based
configuration change is allowed. Report accepted seed/RTK coverage and all
reference-origin assumptions, without calling a component pass practical parity.
