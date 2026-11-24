# KIT-200 outcome (stopped before R3: T9 does not pass on head, for a reason outside this change)

Ticket r1.2. Branch `kit-200-prot-da-proteins-without-values` from `v0.9-staging` 67d38df3.

## Runs

| Step | Run | Commit | Job minutes | Core-hours | Result |
|---|---|---|---|---|---|
| R1 | 45321064526 | ceff82fe (tests only) | 11.7 | 0.78 | instrument defect in the fixture for seven chain tests; the rest at their stated values |
| R1-2 (instrument run 1 of 3) | 45323531564 | bd907fa2 (tests only) | 10.0 | 0.67 | every test at its stated value; every guard passes |

Runner: R 4.6.1, limma 3.68.4 (`library.tsv`), Glimma 2.15.2, reference BLAS, ubuntu-22.04. Tables: `r1/`, `r1-2/`.

**The instrument defect (R1).** `kit200Object` built the protein table with `data.frame(Protein.Ids = rownames(x), x)`,
which keeps `x`'s row names; the S4 DA's `column_to_rownames` refuses such a table ("`.data` must be a data frame
without row names."). T6, T6 (b), T6 (c) LFQ and TMT, T7, T9 and T10 stopped there, before any value was judged.
T6 (d) passed through duplicate removal, which rebuilds the table, and failed at its stated value. Fix: `row.names = NULL`
(commit bd907fa2); no expected value, threshold, comparison or stop rule changed.

## R1 (G0), observed on R1-2 at bd907fa2

| Test | Route | Base message (observed) | Stated | Verdict |
|---|---|---|---|---|
| T1 (treat NA and 0) | kernel, equal df | covariate contains NA or infinite values | same | reproduced |
| T2 (treat NA and 0) | kernel, unequal df | prior.weights contain NA values | same | reproduced |
| T3 | guard, trend off | completes, tested rows identical | passes | passes |
| T4 (a), (b) | guards | complete | pass | pass |
| T5 | all-NA input | No residual degrees of freedom in linear model fits (expect_error, wrong class) | same | reproduced |
| T6 | DIA in process, replicate filter skipped, m = g = 1 | prior.weights contain NA values (unwrapped) | same | reproduced |
| T6 (b) | groups 3, 3, 6, defaults, replicate filter skipped | prior.weights contain NA values | same | reproduced |
| T6 (c) | LFQ, TMT | prior.weights contain NA values (both) | same | reproduced |
| T6 (d) | default QC chain, normalisation RUV skip, four samples removed | prior.weights contain NA values | same | reproduced |
| T7 | managed (worker) | prior.weights contain NA values | same | reproduced |
| T8 (a) | declared design, W in the model | prior.weights contain NA values | same | reproduced |
| T8 (b) | declared none | prior.weights contain NA values | same | reproduced |
| T9 | `daAnalysisWrapperFunction` | prior.weights contain NA values | same | reproduced |
| T10 | replay (record cannot be made) | prior.weights contain NA values | same | reproduced |
| T11 | lipid and metabolomics exported S4 | purrr's "In index: 1 ... Caused by error in `fitFDistUnequalDF1()`: prior.weights contain NA values" | same message, wrapped | reproduced |
| T12 (a) | guard, trend and robust off | completes, tested rows identical | passes | passes |
| T12 (b) | trend off, robust on (measurement) | completes; equal to the reference: TRUE; largest raw_pvalue difference 0 | measurement | recorded |

Every precondition (T6 (i), (ii); T6 (d) (i), (ii); T8's R02 inputs) passed on base. Counts: 127 passed, 16 failed,
1 error (T5's `expect_error`, whose class differs; its message is the stated one).

**A2:** every route of 2.3 marked "Stops" is reproduced on base by its chain test. Stop rules 1 to 3 did not fire.

The tab's text: `runSavedOperation` passes limma's message through unwrapped (T6, T7); the tab prefixes
"Error in DA analysis: " (not run here).

## Substitutions and test design, named (section 4: "the implementer names any substitution")

- **Correlation filter (T6 to T10).** The real pairs operation runs (`computeDIAReferenceCorrelation`, threshold 0.5).
  The samples to be removed are made uncorrelated by permuting their values across proteins (seed 200002), so they fall
  below the threshold; no hand-built pairs table, so the replay (T10) recomputes the same pairs.
- **T6 (d) owner.** `protNormCoreOwner` starts from its own recipe at `protein_replicate_filtered`, so it cannot take the
  fixture through QC. T6 (d) uses its own owner and the GUI's apply steps (`runProteinIntensityFilterApplyStep`,
  `runProteinDuplicateRemovalStep`, `runProteinReplicateFilterApplyStep`), then `computeDIAReferenceNormalization` with
  RUV skip as KIT-180 T5 does. The normalisation ran; no fallback was needed.
- **T8 (a).** `P_empty` is placed in the stored normalised matrix (`args$ruv_uncorrected_quant`, the DA's input,
  decision M14) after the factors are computed on the fixture; the scatter is set in the fixture before. The precondition
  asserts the row and cells are NA in the analysed samples.
- **T8 (b).** The scatter's cells are those of the "small" recipe's 2nd to 4th proteins in runs 1, 5, 9 and 12. The head
  reference is the same DA declaration on the declared state's table without `P_empty`, in a second store.
- **T9.** The wrapper's RLE, PCA and value-count plots are replaced by empty plots: they precede the DA and are not this
  ticket's subject.
- **T1.** With a complete reference, every statistic except `AveExpr` is compared at 1e-12 (outcome 3: one-sweep and
  per-gene fits differ in the last bits, so `logFC` too); `AveExpr` is compared identical.
- **T5.** The equality with `applyLipidDaContrasts`' message holds trivially on base (both are limma's message); it is
  judged at R2, after KIT-196's merge (R03).

## Stored-reference search (6.2), at the branch point

`git grep -lI "covariate contains NA\|prior.weights contain\|No residual degrees" -- tests/ inst/ docs/qa`: KIT-124
simulation records, KIT-128 records, `parked-20261113.md`, `test-kit127-omics-ruv-factors.R`,
`test-kit128-da-assay-failure.R` (lipid and metabolomics files; not proteomics references). No proteomics stored
reference names these messages.

## After KIT-196's merge

**Merge and R1b.** The branch merged staging 0a4008a6 (KIT-196 48bfe085 and nine other tickets) at 5952c874; `git diff
--name-only 67d38df3 0a4008a6 -- R/` lists 38 files beyond KIT-196's three, so R1 repeated (R1b, 45385950606, 11.6
job-minutes): identical to R1-2, every test at its stated value.

**The fix (9779f827).** `applyProtDaContrastModel`: proteins whose `Amean` is NA (a logical index, KIT-196 audit I06)
leave between `contrasts.fit` and `eBayes`; all-empty input is refused with KIT-196's class and text. The ids are carried
as `features_without_values` (only when present) and appended as rows of NA (`appendProtDaUntestedRows`) after the
q-values and BH in `buildProtDaResultTables` and in both branches of `runTestsContrasts`, so the multiple-testing family
is the tested proteins only. Deviation from ticket 7.1, on the coordinator's instruction (KIT-196 hand-off): the rows are
appended after the adjustment, not before it with the mask excluding them; the family is the same either way.
The artifact field, like the return fields, exists only when there are such proteins (7.1's "no field" identity).

**R2 (A3, A4).**
- 45388517096 is an instrument defect: its request carried 65 files of another ticket's list, read from a scratch file
  another session had overwritten. Its tests are not this ticket's scope (all 65 `same`). Its R CMD check binds the right
  commits and is valid: **A4 met**, nothing only on head (3 WARNINGs and 3 NOTEs on both).
- 45400832062 (R2-2, section 9's 44 files): **43 `same`**, no `new-failure`, no `message-changed`, no
  `compare-error`. The new file (`only-on-head`) failed 8 tests and erred in 1, each in the fixture: no protein id table
  (the DA's annotation join), a state read from one store saved in another, the wrapper given contrasts without
  `full_format`. Fixed at b8d9377c; no expected value changed.
- 45463691826 (R2-3, the new file only, at b8d9377c): **509 passed, 1 failed**. Passing on head: T1 to T8 (with T5's
  equality to the lipidomics kernel's refusal text), T10 (the replayed DA tables identical), T11, T12 (a) and (b).
- **T9 fails on head**: `daAnalysisWrapperFunction` now passes `runTestsContrasts` and stops later, in
  `getSignificantData` ("rownames_to_column(): is.data.frame(df) is not TRUE"). The wrapper nests
  `contrasts_results$results` (a list of tables) one level too deep (`R/func_prot_da_model_wrapper.R:272-275`), so it
  stops on any data (read; not run here). Not caused by this change, and no test runs the wrapper to the end otherwise.
  Changing T9's assertion would change a declared comparison (VAL-54 (i)), so the ticket stops here for a decision.

## Not checked (yet)

- R3 (A5); A1 for T9; A6 and A7 at the end.
- The tab's displayed text after "Error in DA analysis: " (in-process wrapping only observed at the runner).
- 2.3 and 2.4's citations re-read at the branch point (A7): pending with the fix commit.

## T9 changed by the owner (VAL-180, option (a))

The owner changed T9's declared expectation (not the implementer): T9 now asserts only that the wrapper's
differential-abundance step (`runTestsContrasts`, captured inside `daAnalysisWrapperFunction`) completes and reports
`P_empty` as not tested (`features_without_values`, one untested row per table). Commit 817b026a; run R2-4.

## Unrelated findings

**`daAnalysisWrapperFunction` stops after its DA step on any data (not fixed here; its own ticket, VAL-180).**
- Read (`R/func_prot_da_model_wrapper.R`): `contrasts_results_table <- contrasts_results$results` (around `:234`) is a
  named list of tables, one per contrast. At `:272-275` it is wrapped again, `nested_list <- list(contrasts_results_table)`,
  and passed as `getSignificantData(list_of_da_tables = list(nested_list), ...)`. `getSignificantData`'s inner map
  (`R/func_prot_da_results_summary.R:250-266`) expects each element of `nested_list` to be a table; it receives the list
  of tables and stops in `rownames_to_column()`: "is.data.frame(df) is not TRUE" (observed in R2-3, run 45463691826,
  after `runTestsContrasts` had returned). The same nesting feeds `printCountDaGenesTable` at about `:306-308`.
- Not run on base data without empty proteins; the reading says it stops on any data. No test runs the wrapper to its end.
- A failing test for it: `daAnalysisWrapperFunction(obj, contrasts_tbl = data.frame(contrasts = "groupB-groupA",
  full_format = "B_vs_A=groupB-groupA"), formula_string = "~ 0 + group", treat_lfc_cutoff = 0, args_row_id = <id column>)`
  on any complete two-group `ProteinQuantitativeData` (with a `protein_id_table`; RLE, PCA and value-count plots may be
  replaced) completes, and its `significant_rows` holds one row per protein with `log2FC` equal to
  `contrasts_results$results[[1]]$logFC`. Today it fails at the error above.
- Second, smaller: when `contrasts_tbl` has no `full_format` column the wrapper builds the contrast strings with `sapply`,
  which names them by the raw contrast; `makeContrasts` then names the contrast column by that name, and `topTreat`'s
  `coef` (the full string) is out of bounds (observed in R2-2, run 45400832062).

### For the change report

The exported scripted wrapper `daAnalysisWrapperFunction` cannot complete: it nests its per-contrast result tables one
level too deep before drawing the significance table (`R/func_prot_da_model_wrapper.R:272-275`). KIT-200 only made its
DA step complete on data with proteins without values; the wrapper's own defect is left to a separate ticket (VAL-180).
