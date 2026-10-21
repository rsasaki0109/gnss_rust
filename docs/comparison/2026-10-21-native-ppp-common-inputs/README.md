# Whole native PPP on shared decoded PPC observations

The actual pinned C++ `PPPProcessor::processEpoch` is now compared with Rust
on GPS 1C/2W at the original 5 Hz rate: six runs, 58,634 epochs per profile,
four profiles, 234,536 output rows and 117,268 native calls. This is a controlled
GPS-only diagnostic; practical RTK/PPP parity remains **unaccepted (18/24, 75%)**.

The [comparison contract](../../ppc-native-ppp.md) describes the input and model
boundaries. The native RINEX observation parser is not run. Rust decodes the
observations into a shared ledger; native SP3/CLK/NAV/ANTEX and explicit BIA
compatibility-view loaders are run. All 894,040 selected observation groups
retain their float bits, missing values and raw flags in both native profile
readbacks. Native LLI/SSI scalars merge absent and zero for estimator use;
their source presence remains external metadata. All cases have normal flag-0
epochs and no explicit receiver offset.

Both sides use matching CODE rapid products, satellite IF-PCO and nominal-yaw
wind-up. Unknown receiver calibration, satellite PCV, tides, ocean loading
and AR are absent/disabled. The explicit BIA view repeats OBS1 into OBS2 for
64 existing GPS code OSBs because the native whitespace loader requires that
field. Source rows, values, sigmas, units and half-open validity are audited;
no phase OSBs or zero corrections are invented. All 12 native code-bias
readbacks match, with unit conversions within the existing 1e-10 m bound.

## Solution coverage

Both Rust profiles use raw-carrier FLOAT and the original random-walk clock.
`rust_native_head` selects the previously verified native head priors and
climatological ZTD. `native_builtin` retains native internal SPP initialization
and fallback. `native_header_zero` changes only the initial seed to the RINEX
header position and clock zero; future native policies are retained.

| Profile | PPP FLOAT | SPP fallback | Invalid | PPP with negative filter ZTD |
| --- | ---: | ---: | ---: | ---: |
| Rust original head | 141 | — | 58,493 | 0 |
| Rust native head | 212 | — | 58,422 | 0 |
| Native built-in | 1,959 | 54,104 | 2,571 | 852 |
| Native header/zero clock | 4,660 | 51,403 | 2,571 | 3,352 |

Native valid coverage is 56,063/58,634, mostly SPP. These SPP outputs are
never relabelled PPP. Rust currently has no equivalent fallback in this path.
Negative native ZTD is counted and retained in the accuracy scores. Extra
native solutions are not evidence that Rust's physical checks should be removed.

## Accuracy and common epochs

The following are antenna-origin 3D RMS errors, in metres. The retained report
also scores IMU origin and both Tokyo lever arms, without selecting a preferred
origin. Reference truth is read only by the scorer after solver execution.

| Run | Rust original PPP | Rust native-head PPP | Native built-in PPP | Native header PPP | Native built-in SPP |
| --- | ---: | ---: | ---: | ---: | ---: |
| Nagoya 1 | 11.624 | 10.292 | 23.629 | 22.864 | 75.421 |
| Nagoya 2 | 7.787 | 4.325 | 23.415 | 442.220 | 158.946 |
| Nagoya 3 | 8.374 | 6.518 | 455.827 | 16.532 | 57.745 |
| Tokyo 1 | 5.382 | 3.976 | 26.623 | 26.076 | 103.357 |
| Tokyo 2 | 1.703 | 1.499 | 11.620 | 17.605 | 38.450 |
| Tokyo 3 | 4.611 | 5.460 | 7.814 | 13.646 | 20.226 |

Each column uses its own successful epochs. Common/new/lost epoch sets are
separately retained for four comparisons per run. For example, Nagoya 3 has
118 common PPP epochs: native built-in RMS 2.150 m versus Rust original
8.374 m. Conversely, the native all-PPP RMS there is 455.827 m, with a
1,559.580 m maximum. Success-set averages must not be treated as paired accuracy.
These constrained GPS-only results do not characterize native multi-GNSS service
performance. Neither solver has met practical acceptance in this comparison.

## Timing, qualification and reproduction

Per-process `wait4` wall/CPU/peak-RSS and per-epoch timers are retained.
Native timers cover `processEpoch`, including SPP fallback. Rust separates
filter calls from epoch staging after decoding. Features, input staging and
success sets differ; these measurements do not accept runtime/memory parity.

The GNSS library source is unchanged from the qualified 413-test epoch-clock
commit. Fresh release examples, fmt, all-target clippy and release checks pass.
All six Rust original-head outputs, projected onto their previous 29 columns,
are byte-identical to the first real-PPC PPP archive. C++ builds with
`-Wall -Wextra -Werror` against the unchanged pinned reference libraries.

The retained 128-epoch pilot includes its own pre-timing adapter source/hash;
all its non-timing fields match the first 128 full-run outputs. Its last ZTD is
negative. The initial adapter build failure and corrected build log are retained:
misleading indentation and an uninitialized Eigen temporary were fixed in the
adapter, with no native source patch or warning suppression.

```bash
python3 -B docs/comparison/2026-10-21-native-ppp-common-inputs/verify.py
python3 -B docs/comparison/2026-10-21-native-ppp-common-inputs/verify.py --verify-local-inputs
```

The second command additionally rehashes original local inputs and recomputes
all origin/common-epoch scores. Re-run solvers into new directories with
`scripts/build_native_ppc_ppp.py` and
`scripts/run_native_ppc_ppp_comparison.py`; their CLI `--help` documents paths.
Source, inputs, binary/library hashes, commands, outputs, resources and the
recursive evidence seal are retained. Raw products and local binaries remain
in the cloud cache; the archive contains their provenance, not redistributions.

Next work should address initialization/SPP fallback and phase continuity or
recovery, then native multi-GNSS comparison. The observed differences include
native ambiguity seeding, four-row/L1 fallback, prediction during failed epochs
and permissive ZTD, whereas Rust uses strict IF pairs and transactional updates.
These are candidate explanations, not established causes of the entire gap.
