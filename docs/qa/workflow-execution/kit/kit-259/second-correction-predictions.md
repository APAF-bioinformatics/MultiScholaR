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

## R1 (run 45676566887): observed

- **I02 is real.** A3-E5 failed on the old code by the printed label, lipid and metab: the x labels were
  `-3.000000e-01 -2.000000e-01 -1.000000e-01 5.551115e-17 1.000000e-01 ... 1.000000e+00`.
- I05-E1, I04-E1 (both), (ii)-E1, (ii)-E2, (ii)-E3 (four) failed as predicted.
- Instrument defects (VAL-54; no result produced for these rows, nothing else changed): I04-E2 errored (`captured$actions`
  a closure: since KIT-257 the runner receives a function giving the actions); (ii)-G1 failed on all four templates
  because each knit used its own temporary folder and the chunk's message names it. Fixed in the tests-only commit
  417b6d7e (the comparison and every expectation unchanged). Their re-send is R2's base side (one run, below): the
  first of at most three instrument re-sends.

## R2 (fix eaab6e57 + this file against tests-only 417b6d7e; with check): one run, ceiling 10 x 80 + 2 x 90 job-minutes

Scope (77 files): the ticket's two greps at 127011429, its five named files, the audit's train additions (I07:
`test-kit244-limpa-pass-through.R` and the eight dependents), the files naming `removeProtNormOwnFiles`, the four
proteomics templates or the RUV covariates route, and the new file. 5 shards per side. No render: a report render
needs the GUI lanes' projects (KIT-236's render harness, which lands with train 12); the template rows (ii)-E3 and
(ii)-G1 knit the templates' own chunks instead.

| File | Base (tests-only on the old code) | Head | Verdict |
|---|---|---|---|
| `test-kit259b-second-correction.R` | R1's failures, except: I04-E2 fails by assertion (the warning escapes: 1; no notice: 0); (ii)-G1 passes on the four templates | every row passes | `fixed` |
| `test-kit259-pearson-histogram-figure-failure.R` | one failure: B-E1's file list lacks the marker (the declared edit) | 85 passed, 0 failed | `fixed` |
| `test-metab-01k-s4-plot-pearson-characterization.R` | 17/0 | 17/0 (the pattern still matches) | `same` |
| every other file | | | `same` |
| R CMD check | | | identical |

Any other verdict is a result: recorded, the work stops (stop rule 3).
