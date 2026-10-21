# Whole native PPP comparison on common decoded PPC inputs

The reference is pinned libgnss++ `72f3b7c2c4c088dc575286a6f0edf4e407bc499a`.
The adapter calls the real `PPPProcessor::processEpoch`, including internal
SPP initialization/fallback and native measurement/ambiguity/filter/recovery
policies. It does not rewrite native GNSS source. Private access exposes
diagnostics and an explicitly labelled initial-header/zero-clock control.

The controlled comparison selects GPS 1C/2W on all six unmodified PPC rover
files, at their full 5 Hz rate. Rust exports exact decoded observations into
a shared ledger, including code/phase/Doppler/SNR option distinctions, raw
LLI/SSI, epochs, receiver-offset option and RINEX-header positions. C++ imports
those values into native `ObservationData` and writes a value/flag readback.
The native RINEX **observation** parser is not executed by this benchmark;
do not label this native RINEX pipeline parity. Original RINEX bytes/hashes and
the decoded ledger remain retained. Native LLI/SSI scalar zero merges absent
and explicit zero for estimator use, while the shared ledger/readback retains
their source presence as external metadata. PPC here contains only normal
flag-0 epochs and no explicit receiver offsets; unsupported input is rejected.

Both solvers receive matching CODE rapid SP3/CLK, GPS C1C/C2W OSBs and IGS20
satellite PCO. Native's whitespace BIA loader requires OBS2 even on an OSB;
the explicit compatibility view selects 64 existing GPS code OSBs and repeats
OBS1 into OBS2. It preserves every numeric value, sigma, unit and validity
label, and retains the original standard-product hash plus source-row audit.
Rust checks each OSB's half-open validity throughout the whole case interval;
native's imported values/conversions are checked by readback. No zero bias or
phase OSB is manufactured. Raw-carrier FLOAT is explicit on the Rust side;
constant phase-bias gauge equivalence was already independently verified.

Satellite IF-PCO, nominal-yaw wind-up and relativity are enabled. Receiver
calibration remains unknown; satellite PCV, solid Earth/ocean loading and AR
are disabled/absent for these controlled profiles. Native admission, variance,
robustness, minimum rows, ambiguity initialization and SPP fallback are retained
as native behavior, not forced to equal Rust. Native can use an L1-only
fallback when L2 code is missing; Rust requires the selected IF code pair.
Those algorithmic differences belong in the comparison.

Four profiles are predeclared: native built-in initialization, native explicit
RINEX-header/zero-clock initialization, Rust original head and Rust native head.
Both Rust profiles retain the original random-walk clock and raw phase. All
native profiles keep native ordinary precise-mode kinematics and original
internal SPP settings. A pilot observes that the internal SPP uses broadcast
navigation and no DCB product; that is reported, not silently changed.
Header/zero-clock control alters initialization only, without fabricating an
observed SPP seed or changing native future policies. Native first prediction
uses dt=1 s; Rust's first committed prediction uses dt=0, as before.

Each case retains all output statuses, native raw status/validity, reported
and filter diagnostics, source-dependent clock units, decoder readback and
OSB readback. Native SPP clock output is in range metres in the actual solver
path; native PPP output divides by c into seconds. Preserve the raw value and
tag that interpretation explicitly. Native PPP does not publish a velocity
in this path; its internal filter velocity is labelled separately. Count
negative/nonfinite native ZTD and other diagnostics without filtering those
positions out of the accuracy summary.

Report all-reference-origin errors, PPP-only/native-SPP/combined validity,
continuity, common/new/lost epochs and status differences. The reference is
used only by read-only parent scoring after solver execution. Neither missing epochs nor SPP fallback
can be relabelled as PPP. Native reported RMS can refer to a different row
set from Rust postfit RMS; it is diagnostic, not an equal-residual contract.

Process wall/RSS measurements include different input staging. Timed native
`processEpoch` includes its SPP and correction/filter work; Rust timings
separate the filter and epoch staging (bias materialization/metadata/filter).
Feature coverage and staging differ, so these are descriptive timings rather
than accepted runtime/memory parity. This GPS-only comparison does not yet
establish native multi-GNSS service performance or practical RTK/PPP acceptance.
