# Initial GPS IF PPP with constant-velocity prediction

`PppPositionModel::ConstantVelocity(PppKinematicConfig)` explicitly opts in to
the ordinary precise native PPP motion model. Static remains the default.
The state is `[ECEF XYZ, clock_m, zenith_m, ECEF velocity_mps, IF ambiguities]`.
The first five measurement columns keep their old meaning; velocity is inferred
over time through position/velocity covariance. There are no direct Doppler,
SPP, inertial, acceleration or truth measurements. At the first epoch velocity
retains the caller's prior, with its stated uncertainty.

The native ordinary precise path uses `F(position,velocity)=dt`, propagates
the complete `F P F^T`, adds configured position/clock/trop/ambiguity diagonal
noise and adds `q_velocity*dt` to velocity diagonals. This discrete law is not
a continuous-white-acceleration model with additional dt^2/dt^3 cross terms.
Default velocity mean is zero, independent initial variance is 100 (m/s)^2,
and velocity spectral density is 1e-4 (m/s)^2/s, matching native constants.
Other finite independent priors and the existing transactional IEKF/Joseph
update remain intentional Rust choices; this is not whole-native PPP parity.

All position/clock/trop/velocity/ambiguity correlations propagate. Missing or
reset phase arcs remap after the eight-state head. Rejected updates/events
preserve the committed numeric state and velocity; the next prediction spans
the entire time since the last committed state, without double-counting a
failed prediction. Admission geometry uses the predicted position. Output
names the motion model and exposes velocity/covariance only in kinematic mode.
Existing postfit, NIS, convergence, correction and raw-data requirements remain.

The C++ adapter exposes private access in its translation unit only and calls
the actual unmodified native `PPPProcessor::predictState` on 24 dense SPD
states, marginalizing unused native clock slots into the Rust layout. It also
generates 64 moving two-frequency observation epochs with native precise
transmit/geometry/Niell physics, known clocks/zenith/ambiguities, deterministic
noise, and separate truth. Existing orbit/CLK/OSB fixtures are reused unchanged.

The original full-trace clean RINEX 2 mm target is not met: epochs 1/2/3 reach
about 2.93 mm. Retain those failures, without changing that metric. An
independent unrounded-observable control isolates carrier quantization:
restoring native unrounded phase meets 2 mm on all epochs, while restoring
only code precision does not. The raw rounded fixture remains unchanged.
Exploratory initial velocity variance 10/100/1000 and whitened QR did not
explain the effect; neither is used to hide failures. QR is not in the solver.
Their non-predeclared exploratory captures stay under `probes/`.

The final comparison freezes both rounded moving datasets, default kinematic
and static models, initial velocity variance sensitivity 10/100/1000, and
unrounded phase/code/both controls. Report all statuses, full-trace errors and
the original 2 mm target, plus velocity's initial-prior and subsequent epochs
separately. No truth enters a positioning configuration. Also build committed
static source 101240b in a temporary non-Git copy and require byte-identical
old/new outputs on all ten existing static synthetic profiles. Keep native
component parity, synthetic accuracy, and real PPC performance distinct. Real
PPP datum/phase/antenna/corrections and whole-native performance remain pending.
