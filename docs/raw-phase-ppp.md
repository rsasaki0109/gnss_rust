# Explicit uncalibrated FLOAT phase

The default `PhaseBiasMode::RequireSignalCorrections` still requires declared
per-signal phase corrections. `UncalibratedFloat` is an explicit alternative
for real-valued IF ambiguities when phase OSBs are unavailable. It requires
absent phase corrections (`None`), preserves raw phase metres, and estimates
the constant unknown IF phase bias in each arc's real ambiguity. It never
manufactures a zero OSB or permits integer ambiguity resolution. Supplied
phase corrections, including explicit zero, are rejected in this mode.

Code corrections and their clock datum remain mandatory. Missing phase,
half-cycle flags, LLI, receiver time offsets, correction validity, frequencies
and noise retain their existing checks. Phase mode/reference changes reset
an arc; an absent phase correction is distinct from a calibrated zero.
The solution and IF observations record the selected phase model.

An independent C++ generator uses pinned native IF constants and the existing
synthetic hardware biases. Its relative code DSB fixture explicitly declares
an IF(C1W,C2W) clock gauge. On both 64-point moving traces, calibrated phase,
code-only OSBs and relative-code/raw-phase FLOAT must agree in position to
1 micrometre, and ambiguity differences must equal the independently derived
hardware phase IF bias to the same bound. This proves gauge behavior on
synthetic data, not real PPP accuracy or native whole-solver parity.

For PPC, matching CODE rapid orbit/clock/absolute-OSB products are acquired
from BKG for all nine central/padding days. They explicitly name GPS
IF(C1W,C2W) and IGS20 APC. Their clock files use RINEX CLK 2.00; that version
shares the supported AS/AR data layout and is accepted explicitly, retaining
the original header and receiver records. Other unsupported versions fail.

Two CODE free-text lines contain non-UTF-8 bytes. The opt-in
`read_satellite_biases_with_byte_comments` preserves their original bytes,
line endings and source lines, while presenting ASCII byte escapes to the
technical parser. It does not guess an encoding or change numeric records.
Non-ASCII technical records remain errors; the old strict reader is unchanged.

Canonical IGS and BKG IGS20 ANTEX downloads are byte-identical. The selector
retains complete GPS satellite blocks verbatim and records hashes and source
ranges for every excluded block. It does not claim support for the full
receiver/azimuth/RMS ANTEX format. Admission records unique PRN/validity and
then checks the exact SVN. CAS DSB G01 names G063 in July 2024, whereas IGS20
has G049; retain that conflict. Matching CODE OSBs identify G049 correctly.
Receiver calibration is unknown and satellite PCV is not applied. Product
readiness, synthetic gauge equivalence and real PPC PPP remain separate.
