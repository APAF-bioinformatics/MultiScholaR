# KIT-183 outcome (state at the L1 stop rule, 2026-11-23)

Ticket `tickets/kit/KIT-183.md` revision 1.4. Branch `kit-183-integration-test-of-record` from `v0.9-staging` 1e627f70.
**Stopped at the stop rule of ticket 11 after L1**: the simulation takes 55 to 70 s per replicate under S's load, above
50 s; S, the texts' numbers, R2, X, M, B and F are not run. Records: `runs.tsv`, `feasibility.md`, `lifecycle.md`.

## Observed

- F1 (limma 3.68.4 = 3.68.5 as read) passes; F3 (one seed = five) passes; F4 (MOFA2 1.22.0's enrichment) read.
- Failing first (R1a2, run 45142390361, commit cab55462 whose R/ equals 1e627f70): every new test of behaviour fails or
  errors on the base code (T1, T2, T3, T4, T6, T7, T9, T10, T11, T12, T13); T5, T8, T15 skipped (no base files then).
- Head against base (R1b, run 45142468283, b21099d6 against 1e627f70): 22 integration files same,
  `test-integration-enrichment-mofa.R` one more skip (no MOFA2 on run-tests), the new file 163 passed, 1 failed (T3:
  `gsea$test` partial-matched `tests`, a defect of the test, corrected at 11545cb3, not re-run), 4 skipped; the three
  failures of `test-integration-input-transcriptomics.R` identical on both sides (stored laptop oracle); R CMD check
  0 ERROR, 3 WARNING, 3 NOTE on both, no finding only on head; gate-off pages identical to base.
- T13's descriptive size of the change on the corpus (R1b log): 40 invented pathways, called at BH <= 0.05 by camera 0
  and by gsea 0, on proteomics (camera universe 223, gsea list 395) and on metabolomics_GCMS (133, 140).
- The base store (`kit183-base-store/store/`) was not delivered: the capture wrote it, but the runner's artifact upload
  leaves out the hidden folder `.workflow-entities`; renaming it to a visible name was refused by the permission system,
  and that step was stopped. The other base files were delivered and are committed (`gsea-deparse.txt`,
  `gsea-fixture-value.rds`, `base-terms.tsv`, `base-summary.txt`, `provenance.tsv`).
- L1: two runs, both instrumentation failures of the script's last aggregation (located in the second, corrected in
  code, not run); times per replicate and the memory case in `lifecycle.md`.

## Acceptance criteria

| # | State | Evidence |
|---|---|---|
| A1 | partly: F1, F3, F4 recorded; L1's pool counts not measured; L1's stop rule applied (fires) | `feasibility.md`, `lifecycle.md` |
| A2 | met for camera, cameraPR and the sensitivity floor against direct limma calls (T1, T2 green at b21099d6) | R1b |
| A3 | met (T3, T6, T7 green at b21099d6 except T3's test defect, corrected, not re-run) | R1b |
| A4 | met (T4 green at b21099d6; fails on base since the T4 correction) | R1b, R1a2 |
| A5 | not met: T5, T8 need the base store, not delivered; T15 not yet run with the base files | - |
| A6 | met at b21099d6 (T10) | R1b |
| A7 | met at b21099d6 (T11, T12 without the base comparison) | R1b |
| A8 | met at b21099d6 (T13) | R1b |
| A9 | not run (M) | - |
| A10 | not met: S not run (stop rule) | `lifecycle.md` |
| A11 | met (gate-off identical) | R1b scripts job |
| A12 | not run (B) | - |
| A13 | partly: failing first shown (R1a2); mutation run X not run | R1a2 |
| A14 | not run (F) | - |
| A15 | partly: ledger row (1) written; row (2) waits for S | `docs/methods-decision-ledger.md` |

## Inferred

- The package code is exercised end to end on the lane by L1: training, ranking, values, every test, the whole KEGG
  path and the comparators ran for every data set and factor before the aggregation failed (no package error).

## Not checked

- T5, T8 (no base store); T15 and T12's base comparison with the committed base files (no run after they were added).
- T3 after its test correction; the full T14 (MOFA2 lane); cross-platform (B); full suite (F); mutations (X).
- The simulation's outcome (O1 to O5), the numbers in the texts (`[R_cam]` and the others are still placeholders),
  NEWS's size of the change from S; pool counts; MOFA2 warnings at n = 6, 12, 20.
