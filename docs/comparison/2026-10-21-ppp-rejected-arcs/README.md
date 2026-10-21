# Retaining continuous phase arcs after rejected PPP updates

This experiment implements the optional `RetainContinuousArcs` policy described
in the [frozen contract](../../ppp-rejected-arcs.md). The existing
`InterruptAllArcs` policy remains the default. Both retain the same numeric
transactionality, ambiguity priors, noise, minimum rows and validation gates.
Practical RTK/PPP acceptance remains incomplete: **18/24 (75%)**.

Retained arcs have an existing accepted ID, unchanged correction/signal datum
and no current reset reason. Rejection cannot create or accept an arc. It keeps
the previous accepted time/count, GF/MW and wind-up anchors; the next successful
measurement predicts over the full numeric gap. Missing phase, half-cycle/LLI,
slips, events, power failure, changed corrections and excessive accepted-time
gaps still break continuity. Configuration/product/model/arithmetic errors also
keep the old fail-closed behavior. Retained epochs remain INVALID.

This is not the whole native lifecycle: native slip detection updates its raw
GF/MW history before a filter update, whereas this Rust policy preserves the
last accepted anchors. The preceding whole-native comparison is frozen at
commit `c714cd9`; it is linked by its evidence hash, without claiming a new
native run here.

## GPS 1C/2W with the original head and clock

This isolated policy change increases total PPP output from **141 to 224**
of 58,634 epochs (0.24% to 0.38%). Coverage is still inadequate. Accuracy
does not improve uniformly, and the optional policy remains disabled by default.

| Run | Interrupt PPP | Retain PPP | Interrupt 3D RMS, m | Retain 3D RMS, m |
| --- | ---: | ---: | ---: | ---: |
| Nagoya 1 | 5 | 5 | 11.624 | 11.625 |
| Nagoya 2 | 1 | 4 | 7.787 | 7.398 |
| Nagoya 3 | 118 | 115 | 8.374 | 8.437 |
| Tokyo 1 | 8 | 18 | 5.382 | 6.216 |
| Tokyo 2 | 5 | 73 | 1.703 | 2.824 |
| Tokyo 3 | 4 | 9 | 4.611 | 5.685 |

These are antenna-origin RMS values on each profile's own success set. All
common/new/lost sets and other origin assumptions are retained separately;
additional output does not establish improved accuracy. Independent-clock
combinations also lose solutions in some cases. For example, Nagoya 1 1C/2L
original-head independent-clock output drops from 40 to 9. Tokyo 1 1C/2L
native-head independent-clock maximum error increases from 8.898 to 37.610 m.

Across the full 58,634 epochs, totals for every predeclared head/clock pair are:

| IF pair | Initial head | Clock | Interrupt PPP | Retain PPP | Rejected epochs retaining arcs |
| --- | --- | --- | ---: | ---: | ---: |
| 1C2W | Original | Random walk | 141 | 224 | 3,214 |
| 1C2W | Original | Independent | 418 | 372 | 4,391 |
| 1C2W | Native | Random walk | 212 | 293 | 3,332 |
| 1C2W | Native | Independent | 442 | 404 | 4,300 |
| 1C2L | Original | Random walk | 152 | 170 | 3,488 |
| 1C2L | Original | Independent | 406 | 325 | 4,823 |
| 1C2L | Native | Random walk | 160 | 191 | 3,436 |
| 1C2L | Native | Independent | 392 | 361 | 5,456 |

All aggregate success rates remain below 1%. The 120 retained worst-origin
records include the 37.610 m new Tokyo 1 point after a **2,113.6 s** numeric
outage. It has five code rows, four phase rows, code RMS 1.451 m and phase RMS
4.56e-9 m. Estimated speed is 0.588 m/s versus reference ENU speed 10.586 m/s.
No phase arc was retained on the preceding epoch; this is recovery after the
ordinary maximum gap, not a phase accepted through the 120-second limit.
Neither low residuals nor retained history establish a valid position.

Initialization and recovery remain concerns. In the original-head/random-walk
Nagoya 1 1C/2W retained run, the first accepted ZTD is 0.122 m and the fifth
is 0.015 m; the following update returns `NonPhysicalState`. Nagoya 2's first
accepted ZTD is 0.030 m, followed by the same error. These retained provisional
states are not realistic total atmospheric delays. Keeping their ambiguities
cannot by itself repair a poor initial geometry/troposphere estimate. Future
work should investigate initial ZTD/ambiguity coupling and motion/recovery,
alongside SPP initialization/fallback and multi-GNSS coverage, with gates kept
explicit and unchanged for controls.

## Validation and evidence

Seven meaningful additions bring the full Rust suite to **420 tests**. They
cover unchanged accepted time/count and GF anchor after rejection, cumulative
120-second expiry, 12 intrinsic discontinuities, non-acceptance of rejected
new IDs, unchanged default interruption, moving recovery, events/invalid input,
and exact full-state/covariance/wind-up equivalence to skipping rejected common
clock observations. The existing satellite-calibration expiry test now also
checks that the optional policy cannot retain arcs across that hard failure.
Existing synthetic accuracy and gate bounds are unchanged.

Fresh fmt, all-target clippy, release and all release examples pass. Qualification
logs are bound to the predeclared experiment. The real-data experiment runs all
six unmodified PPC rover files at full 5 Hz, with both GPS IF pairs, both initial
head priors, both clock policies and retention off/on: **96 runs / 938,144 rows**.
Both origin assumptions and both Tokyo lever arms are scored, including common,
new and lost success sets. Receiver calibration remains unknown. No reference
truth, SPP or Doppler enters the solver.

All source/input/binary hashes, commands, statuses, reference errors, continuity,
rejection reasons, per-process `wait4` timing/RSS and retained arc IDs/counts are
kept. Timings remain descriptive and do not accept runtime/memory parity.
The worst-origin analysis is explicitly post-hoc and read-only; it does not
select a profile, omit output positions or alter gates.

```bash
source /workspace/.cloud-setup/activate.sh
cargo build --locked --offline --release --examples
python3 -B scripts/run_ppp_rejected_arcs.py --output-dir NEW_DIRECTORY
python3 -B scripts/analyze_ppp_rejected_arcs.py --archive NEW_DIRECTORY
python3 -B docs/comparison/2026-10-21-ppp-rejected-arcs/verify.py
python3 -B docs/comparison/2026-10-21-ppp-rejected-arcs/verify.py --verify-local-inputs
```

The verifier rechecks all 96 captures and the 48 old control CSVs byte-for-byte.
With local inputs it additionally rehashes those files, recomputes all origin
and common/new/lost scores, and recreates every worst-origin diagnostic.
The historical whole-native verifier must be run using its frozen commit's
sources; its evidence is not rewritten to match later library changes.
