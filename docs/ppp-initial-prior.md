# Optional ordinary precise-mode initial head prior

`PppConfig::with_native_precise_head_prior(time)` explicitly changes the
initial head prior to the pinned native non-PVA precise-mode conventions.
Caller-supplied initial position and clock means are retained. ZTD comes from
the native compact latitude/season/height climatology, adapted from UNB3m;
this is not a full GPT model. GPS calendar day is explicit. Initial independent
XYZ/clock/ZTD diagonal variances are 100/1e8/25 m². Constant-velocity mode
starts at rest with 25 (m/s)² per velocity axis; static mode remains static.
Defaults are unchanged unless this method is called.

The climatology returns pressure, temperature, water-vapour pressure, dry/wet
zenith delay and their total. Negative height is clamped to zero as upstream.
Invalid calendar/coordinates or nonphysical/nonfinite values fail explicitly.
The optional prior does not silently replace a failure with constant 2.3 m.
Climatology is an initialization model, not a new slant-delay observation or
an external constraint fused into the filter.

Ambiguity initialization, process noise, observation noise, correction models,
NIS/postfit/convergence gates and transactional recovery remain unchanged.
In particular this does not port native measurement-derived ambiguity priors,
the SPP seed, native failure recovery, elevation variance or the entire native
PPP initializer. These remain separate comparison and implementation tasks.

The native compact climatology is called for 312 latitude/height/day cases,
including both hemispheres, coefficient boundaries, negative/zero/10 km
heights, seasonal extrema and leap-year day 366. Every environmental and delay
output is compared at the existing 1e-10 component bound. Six actual native
`initializeFilter` calls on independent RINEX-header seeds verify resulting
head states and variances. Original synthetic limits and default PPP behavior
are retained.

PPC diagnostics predeclare a two-by-two initialization comparison: original
versus native initial head variances, and constant 2.3 m versus climatological
ZTD. Each combination runs both 1C/2W and 1C/2L calibrated-phase profiles on
all six raw 5 Hz rover files. No truth selects a profile or changes a gate.
Compare the original profile byte-for-byte with the frozen first PPC PPP run;
retain failures, common-success versus newly accepted errors, both reference
origin hypotheses, continuity, and descriptive resource measurements.
