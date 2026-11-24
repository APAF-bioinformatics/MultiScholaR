# KIT-184: outcome (running record; stopped at L1's stop rule)

Ticket `tickets/kit/KIT-184.md` revision 1.5. Branch `kit-184-ruv-factors-in-integration`, from `v0.9-staging`
67d38df3 (KIT-183, KIT-187 and KIT-191 merged). All computation on GitHub-hosted runners; no R on the laptop.

**State: stopped at a stop rule of ticket section 11.** The package part is built and passes its targeted runs. The
simulation's lifecycle run (L1) passed every check but one: a replicate takes 111 to 120 seconds against the 70 at
which the ticket stops and returns to the owner with a re-split (`lifecycle.md`). The simulation is therefore not
declared and not counted, the texts still hold its placeholders, and the full suite has not been run.

**After the stop (2026-11-24).** The coordinator accepted the re-split (16 jobs of 125 replicates), kept the ticket's
ceiling at 303.33 core-hours and limited the one correction S' to the affected jobs within 20.33 core-hours (ticket
section 11, revision 1.5). The declaration is still not committed and S not started: `declaration-draft.md` holds the
declaration as it would be frozen and the four changes L1's formation counts force, and waits for an independent audit.

## What was run, at which commit, what was observed

`runs.tsv` has every run with its minutes. In order:

| Run | Commit | Observed |
|---|---|---|
| R1a 45321879760 | a5229552 (tests only) | every test of new behaviour fails, as declared (`r1a/tests.tsv`) |
| R1 45325536775 (with R0 and X as its scripts) | 939a9d10 against a5229552 | the new file 432 passed, 0 failed; 17 of 25 files `same`, 5 `fixed`, 3 `message-changed` (below); check identical on both sides (0 errors, 3 warnings, 3 notes); the base-store job passes with no block skipped but the cross-platform one; the four mutations each fail the guards they are meant to fail; gate-off pages identical; the writer's cost within its ceilings (`r1/`) |
| R1b 45332543296 (instrument, VAL-54) | 621a32b5 | the two K4 fixtures given missing values: both files pass; gate G2 (`r1b/`) |
| M 45332530206 | 621a32b5 | K13 passes with MOFA2 1.22.0; the covariate result replays offline to identical terms (`m/`) |
| B 45332726284, B24 45336831513 | 83f745ee, 5e5df9b9 (the head's code plus the bundle lane) | the new file on the installed bundle (`cross-platform.tsv`) |
| L1 45333818847 | 21053957 | `lifecycle.md`: the stop |

R1's three `message-changed` files at 939a9d10: `test-workflow-integration-input.R` and
`test-integration-input-lipid-metab.R` each had one failure, the fixture guard `features < rows` of the two new K4
fixtures at two factors (the fixtures had no missing value; fixed and re-run as R1b, where both files pass);
`test-integration-input-transcriptomics.R` has three failures identical on base and head (the corpus manifest's SHA-256
recorded on the laptop). `test-kit129-run-path.R` has six failures identical on both sides.

## Acceptance criteria

| # | Criterion | State | Evidence |
|---|---|---|---|
| A1 | Gates G0 to G5; G1 exact with rank equal to k on every route's fixture | met | `feasibility.md`; R1, R1b, L1 (G3); VAL-94 (G4) |
| A2 | K1, K2 | met | R1 |
| A3 | K3, with the base store | met | R1 and its base-store job (`r1/x-guards.tsv`, variant `head`) |
| A4 | K4 per route | met, two details not built (below) | R1, R1b; `feasibility.md` |
| A5 | K5 | met | R1 |
| A6 | K6 | met | R1 |
| A7 | K7 (restated, VAL-85 (3)) | met | R1; base-store job |
| A8 | K8 | met | R1 |
| A9 | The setting: K7, K9, K12, K16; mutation run X | met | R1; `r1/x-verdict.txt` |
| A10 | K10, K11, K12 | met | R1 |
| A11 | K14, K15 | met | R1 |
| A12 | MOFA2 lane, K13 | met | M |
| A13 | Simulation | **not met: stopped at L1** (seconds per replicate above the stop). L1's other checks passed | `lifecycle.md`, `l1/` |
| A14 | Gate off | met | `r1/gate-off-verdict.txt` |
| A15 | Cross-platform record; the writer's memory and time | met on ubuntu-22.04 and windows-2025 (B) and on ubuntu-24.04 (B24): 439 passed, 0 failed, 3 skipped on each, the cross-platform block passing (`cross-platform.tsv`); the writer: 9.6 s and 740 MiB above the reading before the call at 20,000 x 500 (ceilings 60 s, 1 GiB) | B, B24; `r1/removed-factors-cost.tsv` |
| A16 | Failing first; the guards fail on a broken head | met | R1a; `r1/x-verdict.txt` |
| A17 | Full suite base against head with check | **not run**: it is the merge checkpoint and waits for the final head (the texts filled from the counted simulation) | |
| A18 | Ledger rows | rows 1 (implementation) and 2 (G1) written as `VAL-<merge>`; row 3 (the simulation) waits for it | `docs/methods-decision-ledger.md` |

## Commits

| Commit | What |
|---|---|
| a5229552 | tests first (K1 to K16; nothing under `R/`) |
| feb3d6f5, 939a9d10 | the fix: writer, reader, import, resolver, kernel, key, replay, record, Methods, tab, workbook, Rd, starter, NEWS, ledger row; a local helper renamed |
| 621a32b5 | the two K4 fixtures at two factors given missing values (instrument fix) |
| bfac4be8 | the Linux cross-platform reference (`tests/testdata/integration/kit184-cross-platform.tsv`, from R1's log) |
| 21053957 | `simulate-k184.R`, `score-k184.R` |
| 4ff76465, cc28c6ac, 74686a37 and later | records |

## Where the tests are (differs from ticket section 10 in location only)

K1 to K3, K5 to K12 and K14 to K16 are in the one new file `tests/testthat/test-integration-removed-factors.R` (with its
own fixtures), as KIT-191 kept its U1 to U13 in one file; the ticket places K6, K7, K14 and K15 in KIT-183's file and
K10 to K12 in the record, tab and workbook files. K4 is in the per-route writer files through one helper,
`expectIntegrationRemovedFactors()` (`tests/testthat/helper-integration-input.R`). K13 is in the MOFA2 lane's file.

Existing tests restated in the tests-only commit because the ticket changes their behaviour by design:
- `test-integration-enrichment-tests.R`, KIT-183 T10: a result of a test with a design names its form in its row.
- `test-integration-enrichment-tests.R`, KIT-183 T11: the RUV note beside the selector stays with the tests that have
  no design.
- `test-integration-input-matrix.R`: the file list of an input with a corrected matrix gains `removed_proteomics.parquet`.
- `test-integration-record.R` (`intRecRewrite`): its corrected matrix equals the normalised one, so the writer's classed
  warning is expected and silenced.
- `test-integration-gui.R` (the scan of a saved input): the same fixture; the classed warning is asserted.
- `test-workflow-integration-input.R`, the stand-in test of the classic route after the correlation filter: its states
  differ by a constant per protein, so it now expects `not_determined`, the classed warning and no factors file.

K4, not built as the ticket words it:
- the lipidomics declared design "with and without a missing value in the fixture": one fixture only (without); the
  principal angle to the DA's stored factors is printed for it (6.9e-12 at most);
- the transcriptomics test prints no angle to the state's RUVSeq factors: it asserts instead that the written factors
  span RUV-III's own W recomputed from `fullalpha`.

Departures recorded in the ticket (section 20.1): 6.6's "residual degrees of freedom are 16" is 15 on
metabolomics_GCMS (the corpus's metabolomics lacks one of the 20 biological samples; asserted as n - 4 and as 15);
the workbook does not send the setting with `test = "gsea"` or `"permutation"` (the ticket's words; a refusal listed
under `skipped` would be the stricter reading).

## Inferred, not observed

- That Module CI's failure on the branch's pushes (job "module-ci impacted (transcriptomics)") is not this ticket's:
  its failing expectations are the corpus manifest's SHA-256 recorded on the laptop, the same expectation that fails
  identically on base and head in R1. Module CI is the repository's own push check, outside the ticket's budget.
- That the classic proteomics route's larger residual (3e-09 against 1e-12 elsewhere) is the rounding of its level
  restoration.

## Not checked

- The simulation: nothing is declared or counted (A13). Every number of the texts that comes from it is a placeholder
  ([R3], [I3], [R5], [I5], [R5m], [R6], [I6] in `enrich_values_ruv_unadjusted_note` and `methods_enrich_limits_ruv`;
  `enrich_floor_covariates_note` says its power "has not been measured yet"). **No build may be released or merged
  with them in place.**
- The full suite base against head (A17).
- The tab in a browser (the module is driven with `shiny::testServer` and `renderTags`).
- The distance between the DA's stored factors and the removed ones on a declared design with missing values.
- That the CRAN archive's ruv 0.9.7.1 tarball is the one the lock's library was installed from (G0).
- The scorer's criteria tables on real counts: L1 exercised `score-k184.R` on 5 replicates per cell only.
- The determinism comparison: the KIT-184 part of `score-k184.R` does not make it yet (found 2026-11-24; an instrument
  gap to close before S, `declaration-draft.md` section 5).
- Cell C (not run: KIT-183's outcome is O1).
