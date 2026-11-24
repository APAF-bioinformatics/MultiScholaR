# KIT-182 outcome: proteomics normalisation and the correlation filter use the current quality-control result

Ticket `tickets/kit/KIT-182.md` revision 1.4 (the third re-check's minors C01 to C11 applied; W1, W3, W4 answered by
VAL-99). Authority VAL-48, VAL-91 (D1, D2), VAL-99. Branch `kit-182-normalisation-rerun`.

## Commits

| Commit | What |
|---|---|
| af0b5772 | Tests first: `tests/testthat/test-prot-norm-input-current.R` (F1, T1 to T10) and T8's managed block in `test-prot-norm-worker-figures.R`; nothing under `R/` |
| e0974603 | The fix (`R/utils_workflow_dia_products.R`, `R/mod_prot_norm_workflow_helpers.R`, `R/mod_prot_norm_observer_helpers.R`, the six `R/mod_prot_qc_protein_*.R`), NEWS, ledger row `VAL-<merge>` |
| 922681df, d4e863fc | Instrument repairs of test mechanics (VAL-54), no expected value changed |
| 633e3cbf, 868c907e | Merges of `v0.9-staging` (ledger-only conflicts; the second brings KIT-187) |
| run branch `run-k182-r1-base2` (4a1de896, 2c03214f) | af0b5772 with the same test repairs: the base side of the instrument reruns (its `R/` equals the base's) |

## Runs (`runs.tsv`)

Correction pass: R1b (45326674826) on the merge 5add430b against staging 0b926295 (KIT-177 merged first); it stopped at F1 (h) (mechanics), so R2 was not repeated.


R1 (45294678011): fix against tests-only on the regression scope (31 files) with R CMD check. R1-2 (45302634386), R1-3 (45306612548) and R1-3b (45306797628): instrument runs 1, 2 and 3 of 3 (VAL-54 (iv)) on the two affected files. R1-3 stopped at plan (empty base field) and ran no test; R1-3b repeated it with the base filled. R1's instrument budget is spent: a further instrument defect on R1's scope stops for the owner. R2 (45310866571): full suite on the merge candidate against `v0.9-staging`.

## Observed

- **Failing first (A1).** On the tests-only side (R1-3b base, 2c03214f) the new file has 140 failed expectations and 0
  errors, all in the failing-first tests at their stated values: T1 sources PRF1 and 64 rows normalised (DIA, LFQ, TMT;
  managed input PRF1); T1 (e) PRF1; T3 PRF1 twice; T4 (a) to (g) no refusal, states committed; T5 (b), T6 (b) PRF1; T7 (a)
  the completion notice only; T8 `pre_qc` 0 calls and `pre_norm_pca.png` unchanged; T9 (a) to (c) the steps run and
  commit, managed `apply` called; T10 (a) to (g) correlation committed, managed `pearson-pairs` dispatched; the managed
  T8 block in the second file fails 4 expectations (actions begin with `normalization-ruv`, inputs PRF1, PCA unchanged).
  T2, T5 (a), T6 (a), T7 (b), T9 (d), (e), T10 (e) first part, T10 (f) allowed parts and every F1 case pass on both. On
  the fix every expectation of both files passes (317 and 71).
- **F1 to F3 (A2):** `feasibility.md`. Every F1 expectation met (clause per case recorded); F2 k = 6; F3 all routes ran.
- **Regression scope (A3).** R1: 29 files `same`, the new file `fixed`, `test-prot-norm-worker-figures.R` `same` in R1
  (its T8 block errored on both sides then) and `fixed` in R1-3b. Failures present on both sides with identical messages:
  `test-prot-norm-core.R` (272), `test-kit129-run-path.R` (6), `test-workflow-saved-da.R` (171),
  `test-gui-entry-proteomics.R` (10), `test-workflow-ruv-pearson-diagnostic.R` (2 errors): stored laptop references, as
  recorded before (KIT-156).
- **Check (A4).** R1 and R2: 0 ERROR, 3 WARNING, 3 NOTE on both sides; nothing only on head.
- **Full suite (A5).** R2: 456 files `same`, the new file `only-on-head` (317 passed), `test-prot-norm-worker-figures.R`
  `count-changed` (63 to 71 passed, 0 failed on both; the four existing blocks unchanged), no `new-failure`, no
  `message-changed`; `test-gui-entry-proteomics.R` fails identically on both (10). One head shard of attempt 1 was lost in
  the library install step; `gh run rerun --failed` repeated it.
- **Effect (A6).** `effect.tsv`: per route, the stale re-run normalised 64 proteins on the tests-only commit and 58 on the
  fix (k = 6); managed: the table sent only; the managed action's input was PRF1 and then the current QC result. A mechanism check on a synthetic
  fixture, not an effect size on real data.
- On the merge with KIT-187 (R2, 868c907e) the effect lines are unchanged: k = 6, 64 to 58 proteins (R2's head log of
  `test-prot-norm-input-current.R`).
- **A7.** NEWS under `## MultiScholaR (unreleased; becomes 0.9.1, the first alpha for testing on real data)`; ledger row
  `VAL-<merge>` after the last VAL row; `DESCRIPTION`, `NAMESPACE`, `man/`, lock files and `packaging/` unchanged.

## Inferred (not run directly)

- Replay is unchanged: no guard is on its path (the guards sit in the module steps, `diaNormalizationRequest`,
  `diaCorrelationSource` and the tab), and `test-workflow-analysis-record.R` is `same`.

## Unrelated discoveries (candidates; none changed here)

- (a) Test harness: `omicsManagedChainMocks` mocks `readWorkflowExecutionComponents(manager, reference, fields = NULL,
  ...)` (`tests/testthat/helper-workflow-omics-worker.R:166-168`) while the package default is the three assay tables, so a
  caller relying on the default reads zero-row prototypes under the mock. Test fidelity, not a product defect.
- (b) Product, read only: the native export declares `<ruv_>normalised_results_cln_with_replicates.tsv` and
  `filtered_session_summary.txt` unconditionally (`R/utils_workflow_dia_sessions.R:1630-1632`, `:1647`), the stub following
  the current RUV radio; a missing dependency fails with the misleading "saved table changed during read".
- (c) Helper limit: F1 (g)'s "selective tables require readSavedWorkflowTable" comes from the helper's display read
  (`helper-prot-norm-core.R:115`, `readManagedOmicsResult(..., NULL)`), not from the product's managed QC path.

## Divergences from the ticket text (recorded)

1. `diaNormalisedLineage(manager, spec, head = FALSE)`: the ticket's unused `reference` argument is not taken.
2. `registerProtNormServerObservers`: the generator reaches the tab's run through the default of
   `runNormalizationWorkflowFn` (a closure that adds `generatePreNormalizationQcFn`), not through the call, because
   `test-prot-05b-norm-module-contracts.R:3100` passes a delegate without that formal (it would fail with "unused
   argument"); the observer's call is unchanged.
3. T4 (c): the replicate filter re-commits `protein_replicate_filtered`, which overwrites PRF1 by name; the step-name walk
   then names the step before it, not "Protein Replicate Filter". The test checks the class only in (c). The (c) message is shown by reading only.
4. T8's managed block also needed one optional return field (`inputs`) in `protNormK122Drive`; its existing blocks are
   unchanged.
5. F1 (h) checks the lineage after `restoreProtNormTabState` without a session export (see `feasibility.md`).
6. `rlang::abort` (not `cli::cli_abort`) for the four refusals, so the message text is not reformatted; the texts are the
   ticket's 7.4 texts as proposed (VAL-41; Revert, not Reset, VAL-99 (2)).
7. R2 ran with `check: true` (the merged `R/` holds KIT-187, which R1's check did not cover), which the request did not price: R2 used 1,115.2 job-minutes (74.35 core-hours) against its row ceiling of 1,080 (72.0), 35.2 over: 29.1 for the check, 5.6 for the shard lost in library install and the failed results job of attempt 1, 0.5 in the test jobs. The ticket's total stays within its ceiling (105.6 of 324.0 core-hours).

## Acceptance criteria

| # | Status | Evidence |
|---|---|---|
| A0 | met | ticket sections 14, 21 |
| A1 | met | R1-3b `r1-3b-compare.tsv`, `r1-3b-failures-base.txt` |
| A2 | partly met: F1 (h) as declared (export and reopen) not answered: in R1b the export and the archive's opening ran, the adoption into a new owner stopped on the fixture's missing contrasts table | `feasibility.md`, R1b |
| A3 | met | `r1-compare.tsv`, `r1-3b-compare.tsv` |
| A4 | met | `r1-check-compare.txt`, `r2-check-compare.txt` |
| A5 | met | `r2-full-compare.tsv`, `r2-full-compare-failures.txt` |
| A6 | met | `effect.tsv` |
| A7 | met | the branch diff |
| A8 | met | this file, `runs.tsv` |

## Not checked

- A session export and its reopen into a new owner, then the lineage (F1 (h) as first written): the export stops on this
  fixture on both commits.
- The running application (no Save Design after normalisation pressed in a browser; S9 is tested by `saveState`).
- The other omics (KIT-194, VAL-99 (1)).
- Real data; any effect size.
- Whether KIT-177 changes anything here: it had not merged by R2; R1b covers it if KIT-177 merges first (section 5 of the implementation audit).
- If the pre-normalisation figures' redraw fails (in process `pre_qc()` raises, or on the worker the leading `normalization-qc-figures` action fails or is cancelled), `qc_sources$pre_normalization` already names the new source, so the next run does not redraw and its publication copies the earlier table's `pre_norm_*.png` (read, not run; the tab entry has the same order, `R/mod_prot_norm_qc_generation_helpers.R:8-11`).
