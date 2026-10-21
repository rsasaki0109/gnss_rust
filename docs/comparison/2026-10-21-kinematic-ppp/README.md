# Initial GPS IF kinematic PPP

`PppPositionModel::ConstantVelocity` adds explicit ECEF velocity states to the
existing GPS IF PPP. Position/velocity and all clock/trop/ambiguity correlations
propagate as `F P Fᵀ`; native ordinary precise PPP's velocity diagonal Q is
added after propagation. Static remains the default. Measurements retain the
same five head columns, exact code/phase correction admission and IEKF/Joseph
update. No Doppler, SPP, inertial, acceleration or truth input is fused.

At the first epoch, velocity is the explicit initial prior (zero, variance
100 (m/s)² by default), not a newly observed Doppler velocity. The default
velocity process density is 1e-4 (m/s)²/s. A failed update or event leaves the
committed numeric state/covariance intact; the next prediction spans the full
elapsed time. Phase loss/reacquisition remaps ambiguities after the eight-state
head, and admission uses predicted position. Successful results name the model
and expose velocity/covariance only when enabled. See the frozen
[design and scope](../../kinematic-ppp.md).

## Validation and old defaults

The C++ adapter calls the actual unmodified pinned
`PPPProcessor::predictState` on 24 dense SPD cases. Only its translation unit
exposes private access. Native unused clock slots are marginalized into the
Rust layout. Full means/covariances match at the existing 1e-10 component
tolerance for dt=0/0.2/1/30/120 s and several process-noise values. This verifies
the prediction component, without executing native PPP epoch estimation.

The same adapter uses native precise transmit, geometry and Niell physics to
generate 64 moving epochs with fixed orbit/CLK/OSB fixtures, two noise variants,
separate truth and 576 unrounded code/phase rows. All five output files are
identical on two regenerations. All four linked native library hashes match
the earlier whole-native baseline; source/compiler/build provenance is retained.

Seven integration tests and one component test cover moving positions and
velocity, first-epoch uncertainty, continuous arcs, phase outage/remapping,
rejected-state recovery, invalid motion configuration and measurement precision.
All **398 Rust tests**, release examples, fmt, all-target clippy with warnings
denied and release build pass.

The committed static source **101240b** is separately built from a temporary
copy (no Git checkout/worktree). Existing plain, wind-up, receiver antenna,
satellite antenna and legacy tide profiles, each clean/noisy, reproduce every
old CSV byte: **10 controls × 64 epochs**. Both sides' captures, source hashes,
compiler/logs and the old differing `ppp.rs` source are retained. All current
default static output/status behavior is preserved on these controls.

## Moving results

Source, executable, fixture and configuration hashes precede the final 31
executions: 11 moving-model/precision/sensitivity profiles and 20 old/new static
control runs. Every one of the 64 epoch rows remains present, including invalid
outputs. This is a calibrated synthetic GPS path at roughly 0.13–0.17 m/s and
30 s spacing. It does not validate vehicle dynamics or real PPC's 5 Hz data.

| Moving input and model | PPP FLOAT / invalid | All-output 3D RMS | Maximum 3D error |
|---|---:|---:|---:|
| rounded clean, static | 1 / 63 | 0.376 mm (1 epoch) | 0.376 mm |
| rounded noisy, static | 1 / 63 | 0.696 m (1 epoch) | 0.696 m |
| rounded clean, kinematic | 64 / 0 | 0.594 mm | 2.928 mm |
| rounded noisy, kinematic | 64 / 0 | 0.190 m | 0.696 m |
| unrounded phase only, kinematic | 64 / 0 | 0.052 mm | 0.376 mm |
| unrounded code only, kinematic | 64 / 0 | 0.596 mm | 2.942 mm |
| unrounded phase and code, kinematic | 64 / 0 | 0.057 mm | 0.416 mm |

Static's one successful point cannot be read as better trajectory accuracy.
It rejects the remaining moving epochs without reporting invented positions.
With default kinematic settings, velocity 3D RMS over all epochs is about
0.02017 m/s, including the initial zero prior with 0.1604 m/s truth difference.
From epoch 2 onward it is about 0.00222 m/s with maximum 0.00310 m/s. The
output remains a model-derived state estimate, not direct Doppler velocity.
Noisy position RMS after the first four epochs is 0.109 m; full results above
retain startup and every invalid point.

The original **clean full-trace 2 mm target is not met on rounded RINEX**.
Epochs 1/2/3 exceed it, peaking at 2.928 mm. This target and all failures stay
in `report.json`; it is not widened to claim success. After epoch 3 the clean
maximum is 0.944 mm. That diagnostic interval was chosen during investigation,
not as held-out acceptance.

Independent native unrounded observations isolate the cause: restoring phase
precision makes all 64 epochs meet 2 mm, while restoring only code precision
does not. No truth position or velocity is supplied to the solver. Initial
velocity variance 10/100/1000 leaves the rounded maximum about 2.93 mm. A
whitened-QR trial also changes it negligibly. Those exploratory, non-predeclared
captures and unused trial code are under `probes/`; QR is not retained in the
solver. The final default remains the native 100 variance and existing Joseph
update, with every old numerical gate unchanged.

This stage does not execute real PPC PPP or full native PPP epoch estimation.
Real clock/phase datum connection, product/antenna correction completeness,
urban motion, real-data convergence, solution/status comparison and runtime/
memory acceptance remain outstanding. Practical RTK/PPP parity is unaccepted;
the broad checklist stays 18/24 (75%).

## Replay

```bash
source /workspace/.cloud-setup/activate.sh
cargo build --locked --offline --release --examples
python3 -B docs/comparison/2026-10-21-kinematic-ppp/replay_static_control.py \
  --output-dir /tmp/ppp-static-control
python3 -B scripts/run_ppp_kinematic.py \
  --moving-binary target/release/examples/ppp_kinematic \
  --static-binary target/release/examples/ppp \
  --old-static-root /tmp/ppp-static-control --output-dir /tmp/kinematic-ppp-replay
python3 -B docs/comparison/2026-10-21-kinematic-ppp/verify.py
```

The replay helper reconstructs retained parent 101240b sources from their
verified hashes and archived old PPP file, without a checkout or Git history. Direct
examples accept `ppp_kinematic static|kinematic clean|noisy` and the explicit
unrounded precision variants, with optional velocity variance for sensitivity.
`verify.py` recomputes all archived errors, statuses, exact static controls,
native provenance and qualification/recursive hashes without positioning.
Source hashes refer to this milestone; use its commit after later changes.
