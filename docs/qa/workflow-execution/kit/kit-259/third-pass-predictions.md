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

## Comparison run and render job: predictions added before they are sent.
