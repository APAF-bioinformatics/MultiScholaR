# KIT-214 outcome (interim: step 0, G0 and R1 done; the fix commit waits for KIT-202's merge)

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

- cd88f95f is that commit.

## G0 and R1: run 45472096981 on cd88f95f (one request: a scripts job and one test job; base none)

`git diff 74b1caee cd88f95f -- R/ man/ NAMESPACE` is empty, so the scripts job ran the base's package code (G0 binds
the base by content). Job minutes: tests 11.6, scripts 9.1 (20.7 in all, 1.38 core-hours; ceilings 45 and 60).

Observed, G0 (`g0/g0-base-values.tsv`, `g0/console.txt`; exit 0, all 23 fixtures built):
- C1 base kept rows (`.row`), `whole` and `partitioned` identical on every fixture: E1 tables row 1 (A), E2 tables row 1
  (A), the tie row 1 (the first in input order: the oracle's tie rule holds), the code-built module-CI
  `duplicate_feature = TRUE` object row 2, through the import also row 2, the seed with a duplicated `Alignment ID`
  row 1, A2 lipid_canonical row 12 (the x 2 copy), A2 metab_lc row 5 (the x 2 copy).
- Vendor columns: the module-CI `duplicate_feature = TRUE` object through the lipid MS-DIAL import keeps only `Peak ID`
  (character), `Name` and the four samples (so S1 classes it unaffected at run time); built in code it keeps `RT (min)`,
  `Precursor m/z` and `Total Score` (numeric).
- C2 / C3: metabolite seed 5 / Samples 5, Missingness 0 (frame and shell); seed with ids 0 to 7: 5 / 5, 2.5; lipid
  MS-DIAL with `Peak ID` 1:3: 5 / 5, 13.3.
- C4 on the seed, pattern "ITSD": `ITSD_Caffeine_d9` mean 38,721.4, CV 55.893 %; `ITSD_Leucine_d10` mean 25,641.6, CV
  55.886 %; "IS with CV > 30%: 2"; long-table samples include `Alignment ID`. Objects without standards record the
  refusal "No internal standards found" (expected; no test reads them).
- C5: E9's table counts its numeric annotation (`Average Mz`) on base.
- Every hand-computed base value of section 10 equals G0's (the leucine mean, not stated in the ticket, is 25,641.6).
  No correction to the ticket's values.

Observed, R1 (`r1/tests.tsv`, `r1/failures.txt`, `r1/test-kit214-sample-columns.log`): 122 passed, 44 failed, 0 errors.
Every failure is an assertion; none is an error:
- E1 (4 cases) kept S1 = 50 (head 60); E2 (4) 10 (head 12); E3 and E9 the argument is refused (result NULL); E4 "5";
  E5 Samples 5 (both paths; Missingness 0 passes); E6 Samples 5, Missingness 2.5; E7 Samples 5, Missingness 13.3; E8
  mean 38,721.4 and CV 55.89 % (and the CV line and long table); E10 no classed error (both omics).
- T6 (i) marker row absent; (iii) and (v) no marker in the record; (iv) non-vacuity fails (both rules keep row A on
  base) and record (b) has no marker; the replays themselves report no mismatch.
- G7 (a): `Samples` 3 (head 4), both paths, both omics.
- Passing on base as the ticket expects: T6 (ii), G1 (six objects, C1 to C5 against the oracles), G2, G3, G7 (b).

Inferred: stop rule 1 does not apply (no E row passes on base); the ticket goes on to the fix commit after KIT-202's
merge, with R1b on the merged state.

## Acceptance criteria

- A1: met (G0).
- A2: base side met (R1). R1b and the head side pending.
- A3 to A6: pending (fix commit, R2 with S1, R3).

## Not checked
- The head: the fix is not written (it waits for KIT-202's merge); R1b.
- S1 has not run; it runs first in R2.
