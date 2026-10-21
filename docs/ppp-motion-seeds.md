# Caller-admitted SPP/Doppler motion seeds for PVA

The first primitives reproduce actual native `seedPositionDisagreement`,
PVA `predictState(dt=0, seed)` restart policy, and `applyVelocityMeasurement`.
Independent 32-case native traces per primitive cover dense SPD matrices,
cross-correlations, the strict 3600 m² mean-position variance threshold,
disabled/16.27 total-NIS disagreement gate, missing/invalid seed covariance
fallback to 100 m², unavailable velocity and absent/invalid seed. Reset XYZ,
V and A rows/columns, preserving the full other marginal. Reset variances
are 3600/100/100; the native initial precise head's XYZ variance is 100.
Repeat fixtures twice and retain the existing relative 1e-10 component bound.

Velocity updates operate on the full Gaussian state, including XYZ/clock/ZTD/
acceleration/ambiguity cross-correlations. Native uses a subtractive update;
Rust retains its checked SPD covariance and Joseph update. Normal SPD inputs
must match within component tolerance, with identical gate outcomes. Invalid
velocity/noise/dimensions fail without changing the input state. Zero or
indefinite measurement noise is rejected explicitly even if native LDLT would
accept it. Native total 3-D NIS gate maps to one third of the Rust per-row gate.

These are motion-recovery components, not an independent position-measurement
fusion. A restart based on an SPP solution from the same epoch reuses its data
as a heuristic seed; it must never be described as an independent extra factor.
An eventual PPP connection must require explicit seed time, coordinate frame,
receiver reference and source, keep clock datums/units separate, handle events,
preserve numeric transactionality and report every applied/skipped operation.
Default PPP behavior and all existing gates remain unchanged.

The read-only probe at the worst PVA Nagoya 3 point (GPS 2325/554542.2) invokes
existing GPS/Galileo code SPP and Doppler using broadcast navigation, without
truth or PPP updates. All four probe profiles are retained. Native variance/QC
with GPS/Galileo gives antenna-origin error about 10.589 m and speed 0.0295 m/s,
versus archived PPP 134.476 m/19.635 m/s and reference speed 0.003 m/s. Other
probe profiles remain over 100 m wrong. This motivates integration and quality
checks; it does not establish real-data recovery or parity. Reference origin
is unknown, and both antenna/IMU scores must be retained.

## PPP connection and experiment contract

`PppFilterSettings.motion_seed_policy` defaults to None; the explicit
`process_epoch_with_seed` API admits only valid same-epoch SPP solutions with
matching caller-declared frame/receiver reference and nonempty provenance.
Time must match within 1e-7 s. Seeds import neither clock bias nor clock drift.
Only ConstantAcceleration supports this policy. Events require no seed and
clear arcs. The optional required-seed policy rejects missing SPP solutions.
Initialization uses SPP XYZ, finite available V (else rest), zero A and the
caller/native head variances. Later restart follows the native 3600/16.27
policy and retains the clock/ZTD/ambiguity marginal. Doppler is a full-state
velocity update using the leading velocity covariance marginal, with native
ungated fusion by default. Missing/bad velocity covariance skips that update
with an explicit diagnostic. A restart can use a finite velocity mean without
a velocity covariance, matching the native heuristic; its initial innovation
is zero and therefore is not an independent quality test of that mean.
Every operation remains proposed until the GNSS update and arc commit pass.
Geometry admission previews the corrected head; full-state velocity fusion is
performed once on a fresh prior, never twice using that preview.

Before running the PPC experiment, fix native head, IndependentEpoch clock,
PVA and InterruptAllArcs, with all existing GNSS gates/noises unchanged.
Six cases x two GPS IF pairs x three seed modes = 36 runs / 351,804 raw rows:
shadow (compute unused seeds), position (no velocity supplied), Doppler.
All three compute identical independent GPS/Galileo broadcast SPP with native
variance, default MAD/RAIM/GDOP QC and carrier-locked per-system-clock Doppler.
Each SPP solve starts from the original RINEX header, independent of PPP state.
Broadcast WGS84 versus IGS20 compatibility is explicitly an approximate caller
assertion, with no coordinate transformation or receiver calibration claim.
Doppler defaults remain sigma 0.5 m/s, chi-square/dof 25 and maximum unscaled
velocity variance 400 (m/s)^2. No reference/IMU data enters either solver.
The 12 shadow PPP projections must match the previous 35-column PVA output
exactly. All source/binary/input hashes are fixed before runs; retain all rows,
common/new/lost successes, both reference origins and both Tokyo arms. Report
SPP fallback separately: it never becomes a PPP status or a carrier solution.
Reference origin uncertainty, missing receiver PCV/tides, data-dependent seed
reuse and correlated SPP/Doppler from shared satellite models remain limits.
Runtime/RSS includes seed computation/load/CSV and is descriptive; no new
whole-native PPP run or practical performance acceptance is claimed.
