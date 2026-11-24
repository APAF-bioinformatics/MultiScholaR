# KIT-214 outcome (interim: step 0, tests-only commit; G0 and R1 pending)

Ticket `tickets/kit/KIT-214.md` r1.3 (owner decision VAL-173 (1); coordinator E-1 to E-8). Branch `kit-214` from
`v0.9-staging` 74b1caee.

## Step 0 (before R1)

- Every line number and statement of the ticket's sections 2, 3, 6, 8 and 10 re-read at 74b1caee. 74b1caee differs
  from 9f3b7c49 (the ticket's reading) only in `docs/methods-decision-ledger.md`; every citation holds. The register's
  RC1-01 to RC1-03 are in r1.2 (verified, not re-applied).
- One contradiction with the code, written into r1.3: G7's base column said "same", but RC1-02 made its head
  expectation "`Samples` counts it" (the text design sample). On base `Samples` counts the numeric columns
  (`R/mod_omic_qc_tabs.R:662`), so the text sample is not counted. G7 is split into (a) `Samples` (fails on base by
  assertion; added to R1's failing-first list and A2) and (b) the rest (passes on both).
- R2 scope grep (section 10) at 74b1caee: 33 files with `workflowReplayRule`, 32 without (r1.0 recorded 32 at 91e3cdcb
  without it).

### Step-0 survey (ticket section 6)

Searched `tests/testthat` at 74b1caee for (1) mocks or injected functions with a fixed formals list of a function whose
call this ticket changes (the two duplicate kernels, `resolveDuplicateFeaturesByIntensityFn`, the ITSD builders, the
`norm$itsd_table` and `qc$itsd$metrics` lambdas), and (2) tests that run C1's declared step,
`resolveMetabDuplicateAssayData` or a QC display on an object with a numeric non-sample column against a stored value.

- `test-lipid-01g-qc-duplicates-module-characterization.R:153-157`, `:355-358` and
  `test-metab-01u-qc-duplicates-module-characterization.R:177-181` mock the kernel with `(assay_tibble, id_col,
  sample_cols)` and assert `sample_cols == c("Sample1", "Sample2")`: their designs' samples, so D, and the kernel's call
  shape is unchanged under E-7 (B). Their sample columns are double and their ids text, so the callers' conversion of
  every numeric column is a no-op on what the mock receives.
- `test-lipid-01c-qc-duplicate-module-seams.R:577` mocks the lipid kernel likewise: the file is quarantined (skip on
  line 1).
- `test-metab-01p-qc-duplicates-module-characterization.R:141-144`, `test-metab-01u-...:310`,
  `test-omics-audit-carry-overs.R:187-190`, `test-module-ci-metab-qc.R:80`, `test-omics-qc-core.R:134` call
  `resolveMetabDuplicateAssayData` without `designSamples` (NULL keeps today's rule).
- `test-omic-gui-config.R:692-696` mocks the metric kernels with `sample_columns = NULL` and calls the lambdas with four
  arguments (`:697-700`); `:933-941` mocks the two builders with `(assay_data, <id>, annotation_cols = NULL,
  is_patterns = NULL)` and calls the lambdas with three arguments. The `...` form passes nothing on those calls.
- The frozen QC displays (`test-omic-qc-module.R` `omicQcFrozen`, `omicQcItsdFrozen`; `test-omics-qc-core.R`
  `qcCoreFrozenLipidItsd`, `qcCoreFrozenMetabItsd`, `qcCoreFrozenFigures`) and `test-workflow-omics-application-reads.R`
  use the module-CI or declarations objects: text ids, every numeric column a design sample (D = N).

Result: no hit. Nothing is put to the owner (VAL-130 (4)).

## Commits

- (tests only; nothing under `R/`, `man/`, `NAMESPACE`): `tests/testthat/test-kit214-sample-columns.R` (E1 to E10,
  T6 (i) to (v), G1, G2, G3, G7 (a) and (b)), `kit-214/g0-base-values.R` (G0), `kit-214/dup-samples-control.R` (S1,
  runs in R2).

## Not checked

- No run yet: G0's values, R1's failures, and every hand-computed base value of section 10.
- The head: the fix is not written (it waits for KIT-202's merge).
- S1 has not run; it runs first in R2.
