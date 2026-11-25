# KIT-242: outcome (branch `kit-242`)

Ticket r1.4. Branch from `v0.9-staging` e0747fc198817596e3e66397a87e06b2f869c6a7. G0 and C4b are recorded on branch
`kit-242-g0` (records at 7ea50bc1, `docs/qa/workflow-execution/kit/kit-242/`): 2.76 core-hours.

## Tests-only commit 17488ccf: the option-free rows (ticket section 7, r1.4)

`tests/testthat/test-kit242-normalisation-fallback.R`, nothing under `R/`. Rows, lipid and metab: E2 (quantile on a
sample with one value stops, `omics_normalisation_failed`, `normaliser_error`, limma's error as `parent`; and through
the declaration, no normalised state committed), E3 (an assay with no sample columns stops, `malformed_assay`), E16
(legacy RUV-III-C with k = -1 for A2 stops, `omics_ruv_failed`, `invalid_k`), E10 (the record's v2 DA entries carry
`fdr_method`, "BH" with `fdr_fallback_reason` where the manifest records a q-value fallback), NA1 to NA3 (the record
states between-sample normalisation, ITSD and RUV per assay, FALSE for an assay the step left unchanged; every assay
the RUV record names is in the RUV state). Seeds 242000 + row number. None runs the correlation filter, so they
assume nothing about KIT-255's defect.

## R1 (run 45602184961, 9.9 job-minutes, 0.66 core-hours): observed

Every one of the 14 tests fails on base (10 expectations passed, 46 failed, 2 errors), at the base values G0 observed:

| Row | Base value observed |
|---|---|
| E2 | no error is raised (G0 C2: a warning and the assay unnormalised); the error expectation fails, then the fields are NULL; `conditionMessage(NULL)` errors, which ends the test before its declaration arm (that arm runs only once the error exists) |
| E3 | no error is raised; fields NULL |
| E16 | no error is raised (G0 C12: the assay is dropped with a logger line); fields NULL |
| E10 | `fdr_method` is NULL for both entries; in this draw A1 fell back too (as in G0 C7) |
| NA1, NA2 | the per-assay outcome is absent for A1 and A2 |
| NA3 | the per-assay outcome is absent, and A2, named by the record, is absent from the RUV state |

Stop rule 1 (a row passing on base) did not fire.

## Not checked

- These rows at head: they wait for the fix, which waits for the plan auditor's re-check of r1.4.
- P0: it runs on a scratch commit carrying the fix (ticket 8.2), so it cannot run before the fix.
