# KIT-259 third pass (VAL-307): predictions, written before each run

Authority: VAL-307 (the owner, 2026-11-26): one further bounded pass; AUDIT-KIT-259-impl-20261126 closure check C01
(the report sentence bound to the combined QC figure), W11 (covariates notices), C03/C04 (records), C05/C06 (known
limits). I02, I04, I05, I11 are closed and not reopened. Branch `kit-259b` with origin/v0.9-staging 82a5df02 merged
(ca534f61; clean).

**Expected budget (stated before the first run):** R1 about 15 job-minutes (1.0 core-hour); the comparison run about
83 files at about 520 job-minutes with check (34.7); the render job at most 60 job-minutes (4.0). About 40 core-hours,
so KIT-259 in all about 101 of 131.33 (61.31 spent before this pass).

## R1: tests-only commit bd197791 (R/ and inst/ = 7aedc751 + staging 82a5df02); `test-kit259b-second-correction.R`; one job, ceiling 45 job-minutes

| Row | On 7aedc751's code | Reason |
|---|---|---|
| C01-E1 in process, managed | fail | no list beside the combined figure (NULL, NULL) |
| C01-E2 | fails (its pre part only) | the pre marker is not copied and no list is written; the redraw part passes (no list, marker written) |
| C01-E3 | fails | the list is not on KIT-257's removal list (TRUE, TRUE) |
| C01-E4 | fails | Copy does not carry the list |
| C01-E5 | fails | the session archive does not carry the list |
| (ii)-E3, four templates | fail | the template reads the project's markers, not the list beside the figure: "" |
| C01-E6, four templates | fail (first expectation) | markers in the project folder print a sentence; the list-without-figure and empty cases print nothing |
| I04-E1 in process, managed | fail | the covariates notice names "RUV-corrected", not "RUV-III-adjusted (for quality control)" |
| C01-G1 | passes | DIA with limpa writes no list and no combined figure |
| A3-E5, I05, I11, I04-E2, (ii)-E1, (ii)-E2, (ii)-G1 | pass | closed items, unchanged |

A G row failing on the old code, or an E row failing for another reason than stated, is a result: stop.

## R1 (run 45698227390, 14.3 job-minutes): observed

- As predicted: I04-E1 (both), C01-E3, C01-E4, (ii)-E3 (four templates), C01-E6 (four) failed by assertion for the
  stated reasons; every closed row and C01-G1 passed.
- Instrument defects (VAL-54; no result for these rows; expectations unchanged), fixed in the tests-only commit
  9331bbd2: C01-E1 (both) and C01-E2 errored at `lists[[1L]]` (the capture stored a step's NULL with `[[<-`, which drops
  it); C01-E5 skipped (the tab's export button refuses under testServer, the fixture's owner not being the test
  session's: the round trip now exports through `runProtNormExportSessionWorkflow`, as test-workflow-native-session.R
  does). Re-sent once, as the comparison run's base side (this pass's first instrument re-send of at most three).

## Comparison run: fix 40a7ae4b (+ this file) against tests-only 9331bbd2 (R/ = 7aedc751 + staging 82a5df02); with check; 5 shards per side; ceiling 10 x 80 + 2 x 90 = 980 job-minutes (expected about 546)

Scope, 83 files: R2's 77 recounted at the base, plus every test naming `drawProtNormQcFigure`, `removeProtNormOwnFiles`,
the proteomics templates, `collectDIANativeGuiFiles`, `restoreDIANativeGuiFiles`, `copyToResultsSummary` or the combined
figure's file names (adds `test-kit236-report-sentences.R`, `test-kit248-correlation-skip-record.R`,
`test-metab-05-summary-module-characterization.R`, `test-module-ci-metab-da.R`, `test-module-ci-prot-enrich.R`,
`test-workflow-checkpoint-harness.R`).

| File | Base | Head | Verdict |
|---|---|---|---|
| `test-kit259b-second-correction.R` | R1's failures, and now by assertion: C01-E1 (both: no list), C01-E2 (its pre part), C01-E5 (the restored project has no list) | every row passes | `fixed` |
| `test-kit259-pearson-histogram-figure-failure.R` | 85/0 | 85/0 | `same` |
| `test-metab-01k-s4-plot-pearson-characterization.R` | 17/0 | 17/0 | `same` |
| every other file | | | `same`, pre-existing failure lists identical |
| R CMD check | | | identical |

## Render job (closure check C.4): `docs/qa/workflow-execution/kit/kit-259/render3.R` at the same head; one scripts job, ceiling 60 job-minutes

Predictions are in the script (written before the run): V0 the head's and staging 82a5df02's templates give equal
text, drawings 3, 3, 3, 1, no placeholder sentence; V1 exactly one added sentence on LFQ, TMT and DIA ("In Figure 1,
panels a (principal component analysis after log2 transformation) and h (Pearson correlation after between-sample
normalisation) could not be drawn by the application; each shows a note giving the reason instead of the plot.", the
rest of the text equal to V0), none on DIA with limpa (no combined figure, no list); V2 (a list with no combined figure)
no sentence. Any V0 difference, or a V1 sentence with a wrong letter or stage words, is a defect: stop.
