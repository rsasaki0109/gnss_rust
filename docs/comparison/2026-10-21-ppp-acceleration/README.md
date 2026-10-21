# PVA prediction and full-rate PPC comparison

`PppPositionModel::ConstantAcceleration` adds ECEF acceleration to the existing
position/clock/ZTD/velocity/IF-ambiguity state. Static remains the library
default, and the PPC example retains CV unless `-pva` is supplied. The
[frozen contract](../../ppp-acceleration.md) keeps measurement gates, noise,
minimum rows, bias/correction rules and ambiguity priors explicit and unchanged.
Practical RTK/PPP acceptance remains incomplete: **18/24 (75%)**.

## Component and synthetic evidence

Actual pinned `PPPProcessor::predictState` calls provide 48 dense SPD cases,
including variable ambiguities, dt=0/0.2/1/10/120 s, both hemispheres and near
poles, varying ENU noise, and ignored XYZ process noise. The complete mapped
mean/covariance agrees with the existing relative 1e-10 bound. Local horizontal/
vertical acceleration sigmas are squared, multiplied by dt and rotated at
predicted XYZ. No extra XYZ/velocity noise is added in PVA.

Six explicit PPC header-seed initializers verify native precise head means
and variances: XYZ 100, clock 1e8, ZTD 25, V/A 100 each. The generator also
initializes the native filter once for each of its 48 prediction cases, for
54 initialization calls per repeat. Two independent generations reproduce all
six files byte-for-byte. The C++ libraries/source are unchanged; private access
is confined to the adapter.

Independent native precise/IF physics generates 64 car epochs at 5 Hz, with
accelerating/turning motion in a local ENU tangent plane. The filter starts
from the perturbed RINEX header and zero V/A; no reference truth is fused.
All PVA traces pass the predeclared 1 cm clean / 1 m noisy position bound and
1 m/s velocity bound after eight epochs. CV is retained as a control:

| Trace | CV PPP / invalid | PVA PPP / invalid | PVA 3D RMS | PVA maximum | PVA max velocity error after eight epochs |
| --- | ---: | ---: | ---: | ---: | ---: |
| Rounded clean | 9 / 55 | 64 / 0 | 0.0806 mm | 0.3756 mm | 0.00205 m/s |
| Rounded noisy | 9 / 55 | 64 / 0 | 0.2675 m | 0.8095 m | 0.00757 m/s |
| Unrounded clean | 9 / 55 | 64 / 0 | 0.0689 mm | 0.4164 mm | 0.00175 m/s |

Existing low-speed 2 mm criteria and their initial rounded failures are
unchanged. Both old CV clean/noisy traces reproduce their previous CSV bytes.
Seven additional tests bring the full suite to **427**, covering actual native
prediction/initialization, ENU variance rotation and absence of extra noise,
invalid states/variances/squared noise, moving accuracy, rejection/full-gap
state equivalence, eleven-head-state ambiguity remapping and events. Full
tests, fmt, all-target clippy, release and all release examples pass.

## Real-data scope and limitations

| IF pair | Initial head | Clock | CV PPP | PVA PPP |
| --- | --- | --- | ---: | ---: |
| 1C2W | Original | Random walk | 141 | 137 |
| 1C2W | Original | Independent | 418 | 8,539 |
| 1C2W | Native | Random walk | 212 | 217 |
| 1C2W | Native | Independent | 442 | 8,621 |
| 1C2L | Original | Random walk | 152 | 233 |
| 1C2L | Original | Independent | 406 | 4,941 |
| 1C2L | Native | Random walk | 160 | 238 |
| 1C2L | Native | Independent | 392 | 4,766 |

Each profile has 58,634 input epochs. Original-head independent-clock PVA
coverage is 14.56% for 1C2W and 8.43% for 1C2L. These overlapping IF profiles
must not be added together. PVA alone does not consistently increase coverage
with the original random-walk clock. Independent-clock CV is the paired
control for the larger increases, so those counts isolate the motion-model
change under that declared clock policy.

Accuracy remains mixed. Nagoya 1 original-head independent-clock 1C2W goes
from 10 to 677 PPP epochs and own-success antenna 3D RMS 9.175 to 3.161 m;
on its five common successes, RMS slightly worsens 11.625 to 11.631 m.
Nagoya 3 1C2L increases 29 to 588 epochs but RMS worsens 6.943 to 17.209 m.
The worst PVA point is **134.476 m**, with five code/four phase rows, code
RMS 1.114 m, phase RMS 4.93e-9 m and NIS 0.778. Its estimated speed is
19.635 m/s against reference ENU speed 0.003 m/s; the last numeric commit is
only 0.4 s earlier. Low residuals do not establish correct motion or position.

All 120 worst-origin rows are retained. Forty-eight post-hoc comparisons reuse
both frozen native 1C2W PPP profiles. For original-head independent-clock PVA
versus native built-in, Nagoya 1 has 43 common PPP epochs: native RMS
15.599 m versus Rust 2.371 m. Tokyo 1 also has 43 common PPP epochs and worsens
3.711 to 4.767 m. Native SPP fallback remains outside PPP-only comparisons;
its combined valid coverage is much higher than the Rust PPP-only path.
Practical parity is unaccepted. Native SPP/Doppler initialization and recovery,
measurement quality/redundancy and multi-GNSS coverage remain priorities.

The complete grid is six runs × two GPS IF pairs × two head priors × two
clock policies × CV/PVA: **96 profiles / 938,144 rows**, at original 5 Hz.
All 48 CV controls must reproduce previous CSV bytes. Native-head V variance
is 25 for CV and 100 for PVA, following actual native initialization; the
native-head comparison includes that declared prior difference. Original-head
V variance is 100 in both models. Phase retention remains off throughout.

Both reference origins, both Tokyo lever arms and common/new/lost epochs are
scored without selecting a preferred origin or profile. Process wall/CPU/RSS
measurements include loading/CSV and remain descriptive. This is a prediction
component port: native SPP/Doppler velocity seeding/fusion, SPP-based restarts
and per-rejected-epoch prediction are absent. Failed updates preserve the
numeric state and the next attempt predicts over the full committed gap.

Post-hoc read-only diagnostics retain worst-origin rows and compare every
1C2W PVA profile against both previously frozen native PPP initialization
profiles. Native SPP outputs remain SPP. Rust uses CODE phase OSBs; the native
control uses raw-carrier FLOAT, whose constant IF ambiguity-gauge equivalence
was verified earlier. These comparisons are diagnostic, not new whole-native
runs or acceptance of matching algorithms, features or runtime.

## Reproduction

```bash
source /workspace/.cloud-setup/activate.sh
python3 -B scripts/build_ppp_acceleration_fixtures.py --output-dir NEW_NATIVE_DIRECTORY
cargo build --locked --offline --release --examples
python3 -B scripts/run_ppp_acceleration.py --output-dir NEW_OUTPUT_DIRECTORY --native-fixtures-dir NEW_NATIVE_DIRECTORY
python3 -B scripts/analyze_ppp_acceleration.py --archive NEW_OUTPUT_DIRECTORY
python3 -B docs/comparison/2026-10-21-ppp-acceleration/verify.py
python3 -B docs/comparison/2026-10-21-ppp-acceleration/verify.py --verify-local-inputs
```

The verifier recomputes synthetic bounds, validates all captures, repeated
native fixtures, qualification and old CSV controls. With local inputs it
also rehashes original files/libraries and recomputes every reference-origin,
CV/PVA and frozen-native comparison plus worst-point diagnostics.

The initial adapter build omitted the native private-header include directory;
its failed command/log and exact source are retained. The experiment's first
description mistakenly counted the explicit grid as 48 profiles / 24 controls.
That attempt was intentionally interrupted after eight PPC profiles, its full
source/contract/output evidence retained, and the complete grid rerun with the
correct 96/48 description. All 16 completed CSV captures from that attempt
must match the final run byte-for-byte. No numerical solver or gate changed
during that correction.
