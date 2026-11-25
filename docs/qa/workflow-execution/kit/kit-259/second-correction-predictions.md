# KIT-259 second correction (VAL-287): predictions, written before the runs

Branch `kit-259b` from origin/v0.9-staging 127011429. Authority: the owner's decision VAL-287 (a second, narrow
correction: AUDIT-KIT-259-impl-20261126 I02, I04, I05, I11; and (ii) a report sentence naming the Figure 1 panels
replaced by a placeholder). New test file `tests/testthat/test-kit259b-second-correction.R`.

## R1 (tests-only commit; R/ and inst/ as staging): one job, ceiling 45 job-minutes

| Row | Predicted on the old code |
|---|---|
| A3-E5 (lipid), A3-E5 (metab) | **fail by the printed label**: the x labels contain an exponent (`e-`), and no label reads as 0 (`5.551115e-17`). **If either passes, the inferred defect I02 does not exist: stop, report, change nothing for I02.** |
| I05-E1 | fails: the three conditions are turned into a placeholder (value is a path, the file exists) |
| I04-E1 (in process), I04-E1 (managed) | fail: no warning notice after the covariates step (count 0); no "Error" notice |
| I04-E2 | fails: the classed warning escapes the delivery (1) and no notice is shown (0) |
| (ii)-E1 | fails: no marker beside the placeholder |
| (ii)-E2 | fails: the markers are not removed (TRUE, TRUE, TRUE) |
| (ii)-E3 (LFQ, TMT, DIANN, DIANN_limpa) | fail: the template has no `QC_metrics_placeholders` chunk |
| I05-G1, I11-G1 | pass |
| (ii)-G1 (four templates) | pass (no new chunk: the region knits to itself) |

Stop rule: a G row failing on the old code stops the work.

## R2 (fix against the tests-only commit; with check): predicted after R1, before it is sent (section below, added then).
