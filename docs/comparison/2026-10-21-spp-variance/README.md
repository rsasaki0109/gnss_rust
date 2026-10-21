# Native SPP measurement variance and QC components

The opt-in `SppVarianceModel::NativeSnrAtmosphere` uses the pinned C++ code,
elevation, CN0 and atmospheric measurement variance. GPS/Galileo retain
independent clocks, the same orbit/clock/group-delay and atmospheric correction
means, and the same QR. Raw observations remain intact, and successful SPP
model evidence records each applied variance. RTK uses only the seed position,
without fusing its covariance. Existing None/GPS/ElevationOnly defaults and
all existing thresholds are unchanged.

Native sigma and sine lower bounds constrain measurement noise; they do not
add clock priors or floor the state covariance. Missing/nonpositive/nonfinite
CN0 gets the same lack of penalty as the native helper while retaining the
raw distinction. Rust rejects invalid active corrections and variance overflow
instead of the native optimistic fallback. See the frozen
[design and scope](../../spp-measurement-variance.md).

## Independent checks

The generator calls the actual unmodified pinned `spp_utils` functions:
32 variance, 24 MAD mask and 24 RAIM RMS-selector cases. All 80 agree in the
finite physical domain, including boundary and first-tie behavior. The two
fixture generations are identical; the four linked native library hashes
match the earlier whole-native PPC baseline. Negative candidate RMS is
explicitly invalid in Rust. MAD masks and RAIM selection remain components:
neither performs SPP refits in this milestone.

Five new tests also cover the existing independent physical 32-point moving
fixture within 3 mm and both receiver clocks within 1e-11 s, raw evidence,
covariance response and a low-CN0 injected code outlier. All 376 Rust tests,
examples, fmt, all-target clippy with warnings denied, and release build pass.
`qualification.json` retains command results and compressed log hashes.

## PPC results

All six cases run three mixed-seed profiles with variance disabled/enabled:
36 executions. The 18 disabled controls reproduce every previous CSV field.
Source, input, binary and configuration hashes precede execution. Every
invalid, event and unmatched row remains recorded. Counts below concern
11,733 paired 1 Hz epochs, not all 58,634 rover epochs at 5 Hz.

| FLOAT profile | Elevation only | Native variance |
|---|---:|---:|
| no detector | 9,728 | 9,728 |
| slip detector | 9,856 | 9,856 |
| slips + existing prefit 3 m gate | 3,880 | 4,976 |

All accepted FLOAT epochs have successful SPP seeds in these runs. This is
not a success rate for every attempted SPP solve. For the no-detector profile,
the common FLOAT and common successful seed sets coincide with each profile's
accepted set. Antenna-assumption 3D RMS, in metres:

| Case | Old seed | Native seed | Old FLOAT | Native FLOAT |
|---|---:|---:|---:|---:|
| Nagoya run1 | 49.715 | 33.557 | 8.989 | 8.562 |
| Nagoya run2 | 56.174 | 32.538 | 17.360 | 17.202 |
| Nagoya run3 | 47.844 | 35.462 | 27.480 | 27.432 |
| Tokyo run1 | 62.956 | 51.175 | 23.130 | 23.078 |
| Tokyo run2 | 24.279 | 16.931 | 15.689 | 15.650 |
| Tokyo run3 | 16.387 | 14.348 | 23.829 | 24.113 |

Normal and detector-only FLOAT RMS improve on five cases and worsen on Tokyo
run3. Maximum errors can worsen despite improved RMS. With the 3 m prefit
gate, all six cases gain output overall but also lose some previously accepted
epochs. Scores on all accepted epochs therefore describe different samples:

| Case | FLOAT old → native | All-output 3D RMS old → native | Common FLOAT 3D RMS old → native |
|---|---:|---:|---:|
| Nagoya run1 | 500 → 637 | 1.653 → 2.833 | 1.590 → 1.597 |
| Nagoya run2 | 639 → 775 | 12.117 → 11.146 | 12.213 → 12.080 |
| Nagoya run3 | 212 → 277 | 10.090 → 16.985 | 10.586 → 10.420 |
| Tokyo run1 | 808 → 1,012 | 7.724 → 7.928 | 7.340 → 7.313 |
| Tokyo run2 | 480 → 752 | 3.972 → 3.369 | 4.027 → 4.062 |
| Tokyo run3 | 1,241 → 1,523 | 2.501 → 3.438 | 2.371 → 2.380 |

In this gated profile maximum error worsens 12.592→52.760 m on Nagoya run1,
75.433→100.122 m on Nagoya run3, and 31.146→39.734 m on Tokyo run3. Better
SPP weighting and more outputs do not establish reliable positioning.

`report.json` scores every profile under antenna and IMU assumptions, including
both conflicting Tokyo lever arms. `analysis.json` retains status transitions,
common FLOAT and common successful seed scores, and counts of epochs gained
and lost. It also compares common successful epochs with the immutable prior
whole-native captures (10,585 RTK including 144 FIX), without rerunning C++.
The standalone code-only SPP selection still differs from native RTK's complete
code/phase adapter. Native SPP geometry QC and repeated MAD/RAIM refits remain
pending. Rust FIX policy was not evaluated here; no held-out or runtime/memory
acceptance is claimed. Practical RTK/PPP parity remains unaccepted at 18/24.

## Replay

```bash
source /workspace/.cloud-setup/activate.sh
cargo build --locked --offline --release --example ppc_rtk_float
python3 -B scripts/run_ppc_spp_variance.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset \
  --binary target/release/examples/ppc_rtk_float --output-dir /tmp/spp-variance-replay
python3 -B docs/comparison/2026-10-21-spp-variance/verify.py \
  --dataset-root /workspace/.cloud-setup/datasets/PPC-Dataset
```

The direct example mode adds `-nativevar` after `-ggspp`, for example
`kinematic-galileo-spp-slips-prefit3-ggspp-nativevar`. `verify.py` recomputes
the archived metrics, exact controls, native fixture provenance and all
evidence hashes without rerunning a positioning solver. Source hashes refer
to this milestone; later code changes should be checked from its commit.
