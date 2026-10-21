# Initial kinematic full-set FIX policy

This is an opt-in extension of the existing conservative policy. The default
remains GPS static. FLOAT state, integers and covariance are never overwritten
by FIX. Partial AR and integer holding remain separate pending work.

Before PPC evaluation, keep the existing ratio, lock, confirmation, correction,
residual, covariance, code/phase, NIS and baseline thresholds unchanged. Count
within-system DD pairs and normalize repeated integer labels separately for
each constellation. Selected-reference and complete-system-reference local
gates are explicit choices.

Kinematic evaluation requires timestamped ECEF velocity and positive-definite
3x3 velocity covariance. Default motion admission: maximum coordinate standard
deviation 1 m/s, maximum sample gap 2 s, at least two confirmations. Missing or
uncertain motion clears confirmation/anchor history and leaves FLOAT. Invalid
motion input/time is an API error and does not advance the policy epoch.

For consecutive candidate/fixed positions, predict displacement with the
trapezoidal mean of the two supplied velocities times their actual GPS-time
interval. Apply the existing 0.1 m history-jump default to the displacement
residual. Velocity covariance only controls admission; it does not widen this
bound or enter a position update. Large intervals restart confirmation. This
constant-acceleration approximation can reject valid motion during abrupt
jerk; no tolerance is changed to force acceptance.

Each FIX still requires fresh integer search, all local gates, lock and
reference-independent repeated labels. A local pass or a single candidate
cannot declare kinematic FIX. The motion check is a consistency veto and does
not establish absolute integer correctness. Doppler and carrier observations
can be correlated. Practical acceptance still requires PPC status/error,
continuity and speed/memory comparison against the whole native solvers.

An initial output-shape smoke run on Nagoya run1, without reading reference
truth, reported maximum-axis velocity standard deviations 1.026..15.750 m/s;
the default 1 m/s cap correctly withheld every FIX. Before scoring PPC, the
recorded experiment will retain caps 1, 2 and 4 m/s as a sensitivity grid, for
both local-reference modes and with/without the optional slip detector. The
policy default remains 1 m/s. This grid does not change the 0.1 m displacement
residual or any existing integer/local-validation threshold and is not a
performance acceptance or a reason to select a more permissive default.
