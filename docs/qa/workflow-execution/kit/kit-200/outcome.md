# KIT-200 outcome (in progress: R1 done; the fix waits for KIT-196's merge, R03)

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

## Not checked (yet)

- The fix, R2, R3 and A3 to A7: after KIT-196's merge (R03).
- The tab's displayed text after "Error in DA analysis: " (in-process wrapping only observed at the runner).
- 2.3 and 2.4's citations re-read at the branch point (A7): pending with the fix commit.
