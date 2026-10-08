# Development goal

Read `GOALS.md` and `README.md` before changing GNSS processing. The persistent
goal is RTK **and** PPP parity with the pinned libgnss++ reference, including
accuracy, solution status, continuity, runtime and memory on matching inputs.
RINEX/SPP milestones do not complete this goal. Advance the next missing
dependency or solver feature and keep the progress/evidence in GOALS.md current.

The upstream reference is commit
`72f3b7c2c4c088dc575286a6f0edf4e407bc499a` of
`https://github.com/rsasaki0109/gnssplusplus-library`. In this cloud environment,
its source snapshot is `/workspace/.cloud-setup/reference/gnssplusplus-library`.
Use the existing isolated checkout; do not create a Git worktree unless requested.

# Validation

In this cloud environment, first source `/workspace/.cloud-setup/activate.sh`.
The repository pins Rust 1.99.0. Standard checks from the repository root:

```bash
cargo test --locked --offline
cargo fmt --check
cargo clippy --locked --offline --all-targets -- -D warnings
cargo build --locked --offline --release
```

`tests/fixtures/README.md` documents the independently generated C++ reference
data. Add meaningful boundary and end-to-end tests for each new measurement or
solver path. Preserve raw tracking codes, missing-value distinctions, LLI,
time scales, units, health and correction provenance. Do not silently discard
unsupported formats or admit stale navigation modulo one GPS week.

Keep component parity, synthetic truth tests, real-data accuracy and runtime
benchmarks distinct. Investigate numerical/status differences rather than
loosening tolerances or disabling features to declare parity. Missing external
datasets are blockers for the corresponding comparison, not passing checks.
Preserve upstream copyright notices when porting code.
