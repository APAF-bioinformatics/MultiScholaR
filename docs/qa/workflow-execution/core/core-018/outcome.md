# CORE-018 outcome: shared process floor after moving eleven app-only imports to Suggests

Decision: **accept** (stage `core-018-floor`, protocol v8.1, 2026-10-09). Both members complete, zero swap, monitor gap ≤ 70 ms, lifecycle checks passed, independent verifier accepted both. Memory alone decides; time is reported.

## Builds

| Role | Build | Source | Install |
|---|---|---|---|
| control-1 | `core_016` (series closure build) | `6d4f7c22` | `browser-smoke.Hz6YwJ` |
| candidate-1 | `core_018` (this ticket) | `906c50dc` | `browser-smoke.gQ40C7` |

Both members: `screen_p14` (8 samples, 960,000 rows), the normal 16 GiB policy, host 19.1–19.3 GiB free at start, no other managed unit.

## Measured (`core/runtime/core-018-floor/`)

| Member | Floor: anonymous peak, validation phase (GiB) | Admission phase (GiB) | Whole-slice anonymous peak (GiB) | Kernel peak (GiB) | Elapsed (s) | CPU (s) |
|---|---:|---:|---:|---:|---:|---:|
| control-1 | 2.38 | 2.42 | 2.64 | 2.72 | 47.0 | 64.1 |
| candidate-1 | 1.16 | 1.15 | 1.30 | 1.40 | 29.5 | 42.5 |
| change | **−1.22** | −1.27 | −1.34 | −1.32 | ×0.63 | ×0.66 |

The floor fell by 1.22 GiB against the 1 GiB rule. The whole-slice peak of this small member fell by the same amount, because on the small fixture the floor is most of the peak. Elapsed and CPU time also fell, by about a third; that is three processes each loading 61 fewer namespaces, and it is reported, not judged.

## What this means for the workflow

- Every managed action (controller, monitor, worker) now starts from about 1.2 GiB instead of about 2.4 GiB, on every input size, because the floor is input-independent (`core/floor/outcome.md`).
- On the full P14 input the pre-rollup floor therefore drops from about 2.4 to about 1.2 GiB. The rollup's own working memory above the floor is unchanged by this ticket, so the full-size peak is expected to fall by about 1.2 GiB, not by a ratio. **This was not measured at full size** (no full-size run in the budget) and is not claimed as a measured full-size number.
- The GUI process itself carries the same namespace, so its baseline falls by about 0.4 GiB; not measured here.

## First attempt (`core/runtime/core-018-floor.v8/`, invalid, preserved)

Under protocol v8 the stage reused the `core:lifecycle:complete` verifier selector. Both members completed with the same readings (floor 2.42 → 1.12 GiB, peaks 2.77 → 1.41 GiB, elapsed 49.0 → 30.2 s), but the independent verifier resolved the selector to the old `core-lifecycle-complete/member` folder and refused on both (`verifier launcher did not exit successfully`), so the coordinator's decision is `invalid` (`control:verifier_accepted`). After the repair (own selectors; protocol v8.1, commits `7ea3730e`, `d5e5cff8`) the v8 members could not be verified because their recorded binding carries the v8 protocol hash, so the stage was re-run as the ticket's correction (user decision; ledger row amended). Two attempts at the standalone verifier against the v8 members are in each member folder's `verification-launch-v81.log` and `verification.rds` (`failed`).

## Checks on the package change (commit `892d7fa0`)

- Fresh process loading the namespace: 127 namespaces, 349 MiB anonymous (closure build: 188, 815 MiB); none of the eleven loaded.
- `tests/testthat/test-optional-packages.R`: 10 passed (fresh-process check; classed refusal `multischolar_package_missing` from eight caller families with availability mocked absent).
- Focused files on the installed candidate: fasta 29, github 21, load-dependencies 23 passed. Heatmap-support (1 error), volcano-glimma (1 error) and enrichment (2 errors) show the same errors on the installed closure build; they expect a source load.
- Hand-edited `NAMESPACE` identical to `roxygen2::roxygenise()` output; no `export(` line changed; no S4 definition touched.
- `R CMD check --no-tests --no-manual --ignore-vignettes` on a tarball built on `/home` (`r-cmd-check.log`): the same 3 warnings and 3 notes as CORE-015.

## Rollback

`git revert 892d7fa0` restores Imports and `NAMESPACE`. The protocol v8/v8.1 stage and its evidence stay.
