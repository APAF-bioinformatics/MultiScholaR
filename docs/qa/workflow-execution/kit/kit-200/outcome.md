# KIT-200 outcome (done; audit repairs applied; ready for the closure check)

Ticket r1.5. Branch `kit-200-prot-da-proteins-without-values` from `v0.9-staging` 67d38df3.

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
  `full_format` (a data-independent defect of the wrapper's fallback, not of the test's mechanics; the test now gives the
  documented `full_format`, no expected value changed; handed to KIT-218, VAL-188 (2)). Fixed at b8d9377c.
- 45463691826 (R2-3, the new file only, at b8d9377c): **509 passed, 1 failed**. Passing on head: T1 to T8 (with T5's
  equality to the lipidomics kernel's refusal text), T10 (the replayed DA tables identical), T11, T12 (a) and (b).
- **T9 fails on head**: `daAnalysisWrapperFunction` now passes `runTestsContrasts` and stops later, in
  `getSignificantData` ("rownames_to_column(): is.data.frame(df) is not TRUE"). The wrapper nests
  `contrasts_results$results` (a list of tables) one level too deep (`R/func_prot_da_model_wrapper.R:272-275`), so it
  stops on any data (read; not run here). Not caused by this change, and no test runs the wrapper to the end otherwise.
  Changing T9's assertion would change a declared comparison (VAL-54 (i)), so the ticket stops here for a decision.

## R3 (A5)

Run 45474317181, 149ab856 against staging 2a1543db, 20 shards per side: 467 files `same` (identical failure lists where
laptop-recorded references fail on both), the new file `only-on-head` and `ok`; no `new-failure`, `message-changed`,
`count-changed` or `compare-error`. 1,263.3 job-minutes (84.2 core-hours), longest job 67.3 minutes. **A5 met**; no stored
reference moved (6.2's verdict).

## Acceptance criteria

- A1 met (R1-2 and R1b on base; R2-4 on head, 516 passed, T9 as the owner changed it, VAL-180).
- A2 met (every route reproduced on base). A3 met (R2-2: 43 `same`; the new file passes, R2-4). A4 met (R2's check).
- A5 met (R3). A6: NEWS under the unreleased heading, ledger `VAL-211` after staging's rows; `DESCRIPTION`, lock
  files, `packaging/electron/package.json` and `NAMESPACE` unchanged; `man/runTestsContrasts.Rd` only.
- A7: 2.3 and 2.4's citations were re-read at each merge (ticket revisions 1.4 and 1.5 list each moved line); limma 3.68.4
  on the runner.
- Divergence recorded: 15 of 50 added lines in `R/` exceed R.md's 100 characters, following the file's idiom (107 such lines at base).

**Core-hours:** 126.4 of 449.4 (expected 85.8): R1 0.78, R1-2 0.67, R1b 0.77, R2 21.9, R2-2 16.2, R2-3 1.13, R2-4 0.77,
R3 84.2. R2's row (36.0) was exceeded across its instrument runs (40.0); five R2 jobs overran 60 minutes (one by 0.2).

## Not checked (yet)

- The tab's displayed text after "Error in DA analysis: " (in-process wrapping only observed at the runner).
- The amended fixtures (b8d9377c) and T9 (817b026a) were not re-run on base; by reading they stop at the same limma message before the changed lines.

## T9 changed by the owner (VAL-180, option (a))

The owner changed T9's declared expectation (not the implementer): T9 now asserts only that the wrapper's
differential-abundance step (`runTestsContrasts`, captured inside `daAnalysisWrapperFunction`) completes and reports
`P_empty` as not tested (`features_without_values`, one untested row per table). Commit 817b026a; run R2-4 (45471995943): 516 passed, 0 failed. **A1 met on head.**

**Second merge.** The branch merged staging 2a1543db at a786abd2 (ledger conflict: staging's rows first). KIT-205 (D-B) is on staging: the enrichment background now leaves rows without a p-value out, so NEWS and the ledger row say untested proteins are out of the background (c3f27bd3). R3 binds 149ab856 against 2a1543db.

## Unrelated findings

**`daAnalysisWrapperFunction` stops after its DA step on any data (not fixed here; its own ticket, VAL-180).**
- Read (`R/func_prot_da_model_wrapper.R`): `contrasts_results_table <- contrasts_results$results` (`:223`, `:232`) is a
  named list of tables, one per contrast. At `:272-275` it is wrapped again, `nested_list <- list(contrasts_results_table)`,
  and passed as `getSignificantData(list_of_da_tables = list(nested_list), ...)`. `getSignificantData`'s inner map
  (`R/func_prot_da_results_summary.R:250-266`) expects each element of `nested_list` to be a table; it receives the list
  of tables and stops in `rownames_to_column()`: "is.data.frame(df) is not TRUE" (observed in R2-3, run 45463691826,
  after `runTestsContrasts` had returned). The same nesting feeds `printCountDaGenesTable` at `:306-308`.
- Not run on base data without empty proteins; the reading says it stops on any data. No test runs the wrapper to its end.
- A failing test for it: `daAnalysisWrapperFunction(obj, contrasts_tbl = data.frame(contrasts = "groupB-groupA",
  full_format = "B_vs_A=groupB-groupA"), formula_string = "~ 0 + group", treat_lfc_cutoff = 0, args_row_id = <id column>)`
  on any complete two-group `ProteinQuantitativeData` (with a `protein_id_table`; RLE, PCA and value-count plots may be
  replaced) completes, and its `significant_rows` holds one row per protein with `log2FC` equal to
  `contrasts_results$results[[1]]$logFC`. Today it fails at the error above.
- Second, smaller: when `contrasts_tbl` has no `full_format` column the wrapper builds the contrast strings with `sapply`,
  which names them by the raw contrast; `makeContrasts` then names the contrast column by that name, and `topTreat`'s
  `coef` (the full string) is out of bounds (observed in R2-2, run 45400832062).

### The wrapper (not changed here; KIT-218)

The exported scripted wrapper `daAnalysisWrapperFunction` cannot complete: it nests its per-contrast result tables one
level too deep before drawing the significance table (`R/func_prot_da_model_wrapper.R:272-275`). KIT-200 only made its
DA step complete on data with proteins without values; the wrapper's own defect is left to a separate ticket (VAL-180).

## Audit I05: the one correction pass (a result)

The implementation audit (`AUDIT-KIT-200-impl-20261125.md`, I05) found that on a matrix without row names the tested rows
after an empty protein were labelled one position too low (direct scripted calls only). This is the ticket's one
correction pass for a result. Tests-only commit a51b230f adds the register's test; run C1 shows it failing on that head.
Fix dc23416c: the register's lines in `applyProtDaContrastModel` name an unnamed fit's rows by position before the
removal. Observed: C1 (45493507825, at a51b230f): the I05 test fails 18 times (both trend settings, three contrasts): the tested
rows after position 3 carry the label one lower, "3" holds a tested protein, and the expected labels 1 to 61 are not all
present. R2b (below) passes it on the fix.

## R2b on the merge with staging 0cf0e440 (audit I01)

The branch merged staging 0cf0e440 at 48c1b100 (ledger conflict only: staging's rows first, KIT-200's after VAL-190).
R2b (45493851112): 58fbe5cc against 0cf0e440, section 9 re-grepped on the merged tree (47 files: R2-2's 44 plus
`test-kit205-enrichment-background.R`, `test-kit208-gprofiler-identifiers.R`, `test-reference-mapping-enrichment.R`), with
R CMD check. **46 `same`**; the new file `only-on-head` and ok (538 passed, the I05 test included); no `new-failure`,
`message-changed` or `compare-error`. Check: 3 WARNINGs and 3 NOTEs on both sides, identical, none only on head
(**A4 on the merged tree**). 304.2 job-minutes (20.3 core-hours, inside the row's 36.0); one job 65.7 minutes. Every R2b
file was in R3's full suite, and only NEWS and the ledger were hand-resolved, so R3 (45474317181) stands (6.4).

**Core-hours, total:** 147.7 of 449.4 (C1 1.08, R2b 20.28 added to 126.4).

## Hand-offs (audit notes I09 to I11, as the register words them)

- **I09 (VAL-162's wording; for the owner).** T4 (b) (EMPIRICAL, passing on base and head): a protein with one value per group (zero residual df) has a finite `logFC` and `raw_pvalue`. By limma's code (`R/squeezeVar.R:17`, `R/ebayes.R:57-66`) its variance is the prior's, s0² at its mean, on d0 degrees of freedom, and its p enters the q-value family. VAL-162 says such a feature "is not tested but stays in the estimation of limma's variance prior". The second half holds; the first does not (KIT-204 r1.1 line 93 says it correctly). VAL-162 is scoped to lipidomics and metabolomics; proteomics does the same with no row. Unchanged by KIT-200.
- **I10 (limpa; INFERRED, for the owner).** `limpa::dpcQuant` returns complete protein values, and the EList is subset to the design's samples at the DA (`R/func_prot_da_model_methods.R:326`). A protein whose observations were all in samples the correlation filter later removed is therefore tested on model-based values with their precision weights, where the limma routes now report it untested. This is limpa's design, not a defect of this ticket. limpa's source was not read.
- **I11 (counts; for the owner).** `R/func_prot_da_model_methods.R:597-606` and `printCountDaGenesTable` put rows with a missing q under "Not significant", and the summary's total counts every row (`R/utils_omic_gui_config.R:2190`; `R/mod_omic_da_panels.R:290` counts `logFC`). After this merge an untested protein is shown as "Not significant" until KIT-204 (VAL-129 (1)) merges. This is pre-existing for proteins missing a compared group, and KIT-204 is already ordered before the release commit.

## For the change report

**What was wrong, who was affected.** The proteomics differential-abundance test (limma with the variance trend and the robust option, the settings of every application route) stopped with limma's message 'prior.weights contain NA values' ('covariate contains NA or infinite values' when every protein had the same residual degrees of freedom) as soon as one protein had no value in any analysed sample: limma gives such a protein no mean, and the trend takes that mean as its covariate. No result was produced. The DA tab showed "Error in DA analysis: …", and the managed worker and a replay failed the same way. In the application such a protein arises only when samples are removed after protein QC by the correlation filter: on the default QC chain the filter must remove every sample holding the protein's values, at least four (two in each of two groups); with the replicate filter skipped, as few as one. Also affected: scripted analyses with missing-value cutoffs of 100 %; the DIA workbooks (they call `protein_deAnalysisWrapperFunction` with the trend on from `config.ini`); direct calls of `runTestsContrasts` with the trend on; and the exported `differentialAbundanceAnalysis()` on lipidomics and metabolomics objects, which runs `runTestsContrasts`. DIA with limpa was not affected (its quantified values are complete; inferred, limpa's source not read).

**Code.** Base 2a1543db (as at the branch point 67d38df3): `R/func_prot_da_model_stats.R` `applyProtDaContrastModel` `:291-321`, with `contrasts.fit` at `:303` and `eBayes` at `:306-307`; `buildProtDaResultTables` `:384-419`; `runTestsContrasts` `:494-831` (topTreat `:593`, topTable `:715`, return `:824`). Fix 9779f827 (unchanged at 3f8203c2): removal and refusal `:304-318`; `appendProtDaUntestedRows` `:348-353`; rows appended after the q-values and BH at `:442-446`, `:730`, `:853`; `features_without_values` at `:336-337`, `:453-456`, `:868-869`; `man/runTestsContrasts.Rd`. Audit I05 (dc23416c) adds lines after `:314`: an unnamed fit's rows are named by position before the removal.

**What changed in results.**
- Analyses without such a protein are byte-identical: the new step is entered only when one exists (R2-2 45400832062: 43 files `same`; R3 45474317181: 467 files `same`, no stored reference moved).
- Analyses that stopped now complete. The protein appears once in each comparison table with every statistic missing: `AveExpr` missing where a table carries it; in the long table and the export it looks like a protein with no value in either compared group. It is not in `fit.eb`, it is listed in `features_without_values`, and it is in no q-value or Benjamini-Hochberg family. Since KIT-205 it is also out of the enrichment lists and background.
- Every other protein's result equals the analysis of the same data without it: identical where a missing value remains, to 1e-12 where the reduced data are complete.
- Input in which no protein has a value is refused with 'no feature has a value in any sample, so nothing can be tested' (`omics_da_no_observed_features`).
- Scripted calls with the trend off, on data with such proteins: their rows move to the end of each table, with `AveExpr` missing instead of NaN, and they leave `fit.eb`. Other results are unchanged, except with the robust option on and unequal residual degrees of freedom, where they can change slightly (T12 (b) on base: equal, largest difference 0).
- Data of one or two proteins, one of them empty: the B statistic changes.

**Evidence.**
- Failing first on the tests-only commits: R1-2 45323531564 (bd907fa2) and R1b 45385950606 (5952c874). Every route stopped at its stated limma message (T1 to T11).
- Passing on head: R2-4 45471995943 (817b026a, 516 passed) and R3 (the new file ok on the merge with 2a1543db).
- Check: R2 45388517096 (9779f827 against 0a4008a6, findings identical).
- Data: simulated only (60 proteins, three groups of four, seed 200001; chain fixtures of 12 and 18 samples). Reference BLAS, limma 3.68.4, R 4.6.1, ubuntu-22.04.
- Not measured: how often such proteins occur in real data.
