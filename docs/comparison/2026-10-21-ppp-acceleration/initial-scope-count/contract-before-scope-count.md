# Optional native PVA prediction component

Add a position/velocity/acceleration option while preserving Static and
ConstantVelocity defaults and all measurement/correction/gating policies.
The Rust layout is XYZ, clock_m, ZTD_m, VXYZ, AXYZ, then real IF ambiguities.
Prediction uses F(X,V)=dt, F(X,A)=dt²/2 and F(V,A)=dt on the full covariance.
Acceleration random-walk noise is defined in local ENU at the predicted XYZ
and rotated into ECEF: horizontal sigma 1 and vertical sigma 0.1 m/s²/√s,
squared and multiplied by dt. Native PVA adds no XYZ/velocity process noise.
Clock/ZTD/ambiguity Q retain the selected policies. Failed updates keep the
last numeric commit and predict over the full gap on the next attempt.

Independent native prediction calls must cover dense SPD cross-correlations,
variable ambiguities, dt=0/0.2/1/10/120 s, hemispheres/poles, varying horizontal/
vertical sigmas and nonzero ignored XYZ noise. Use the existing relative
1e-10 component bound. Reproduce fixtures twice. Six native PVA header-seed
initializers verify native precise head and zero V/A variances 100 each.
No native solver source is modified; private access is limited to the adapter.

A separate 64-epoch 5 Hz synthetic car trajectory, generated independently
with native precise transmit/IF physics, uses accelerating and turning motion.
The trajectory lies in local ENU with zero up component, rotated to ECEF;
the fixed tangent plane introduces only its small geometric height curvature.
Both rounded RINEX and unrounded observations are retained. Predeclare full
trace position bounds 1 cm clean / 1 m noisy, and velocity error 1 m/s after
eight epochs. Keep all failures visible; do not change existing low-speed
2 mm bounds or previously retained initial rounding failures. Compare CV as
a control without requiring it to fail. No reference truth initializes either
filter; both use a perturbed RINEX header, zero V/A and weak independent priors.

This ports the native prediction component, not its whole recovery system.
SPP/Doppler velocity seeding or fusion, native covariance/disagreement-based
SPP restarts, and advancing prediction on every rejected raw epoch are not
included. No inertial/reference truth is fused. Whole native comparison and
practical acceptance remain distinct from this component and synthetic test.

Real PPC evaluation will use all six original 5 Hz files, both GPS IF pairs,
original/native head priors and both clock policies, with retention off:
48 profiles comparing CV/PVA. All 24 CV controls must match previous CSV bytes.
Native-head V variance follows each model's actual native initialization:
25 for CV, 100 for PVA; PVA A variance is 100. Original-head V variance is
100 in both models. Native-head comparisons therefore include that declared
prior difference as well as prediction/state-dimension changes.
Score all statuses, both reference origins, both Tokyo lever arms and common/
new/lost epochs. Runtime/RSS is descriptive; do not accept parity from counts.
