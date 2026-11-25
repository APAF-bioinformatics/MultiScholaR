# KIT-221 R1: purpose, binding, ceiling, predictions (written before the request is sent)

**Purpose.** On the tests-only commit (nothing under `R/` differs from base 176f9d3ff) show that each failing-first
row of the ticket (section 10: W1, W1b, W2, W3, W4, W5, W7, W8, W9) FAILS (not errors) at its stated base value,
that the guard G4 (`test-kit218-da-wrapper-tail.R`) passes, and take the base-side digests of G1 and G3 (`g1.R`).
Record M1 (W1b's base message). Decides stop rules 1, 2 and 6.

**Binding.** package_commit = the tests-only commit (the tip of branch `kit-221` when the request is written);
base_commit "none". Files: `test-kit221-da-several-comparisons.R`, `test-kit218-da-wrapper-tail.R`; one shard.
Script `docs/qa/workflow-execution/kit/kit-221/g1.R` with arguments `kit221-g1 base`, output `kit221-g1`.
`blas: reference` (default). No R CMD check.

**Core-hours (4-core runners, job-minutes x 4 / 60).** Expected: tests job about 15 min + script job about 15 min =
30 min = 2.0 core-hours (ticket row R1: 1.0 for one job; the script job is the base side of G1 and G3, which the
ticket's R1 row does not price). Ceiling: `file_timeout_minutes` 20 and script `timeout_minutes` 30, so at most
about 60 + 40 = 100 min = 6.67 core-hours; the ticket's total ceiling (473.33) is unchanged.

## Predictions, `test-kit221-da-several-comparisons.R` (9 tests; status failed; 0 errors, 0 skipped)

Counts are expectations (a `expect_null` on the caught message counts as one). Predicted: 42 passed, 39 failed,
0 errors, 0 skipped.

| Test | Predicted base passed / failed | The failing assertion(s), in order |
|---|---|---|
| W1 limma, two-row data.frame | 2 / 8 | `unique(significant_rows$comparison)` is `B_vs_A` only (ticket's first assertion); `nrow(significant_rows)` 60 not 120; `nrow(contrasts_results_table)` 60; `num_sig_da_molecules` and its first-go table list `B_vs_A` only; `nrow(list_of_volcano_plots)` 1; no `log2FC:C_vs_A` column; `da_proteins_long` `B_vs_A` only. Passes: the call completed; `da_proteins_wide` has 60 rows. No error. |
| W1b tibble | 0 / 1 | the call stops: `expect_null(run$error)` fails and prints the message; M1 = limma's 'number of items to replace is not a multiple of replacement length' (by reading; recorded, not judged) |
| W2 equal per comparison (limma) | 8 / 2 | `comparisons` of `significant_rows` is `B_vs_A` only; C_vs_A rows 0 against 60. B_vs_A's rows, ids, `log2FC`, `raw_pvalue`, `fdr_qvalue` are identical to the one-comparison call (pass) |
| W3 sibling, two-row data.frame | 8 / 5 | `length(contrasts_results$results)` 1 not 2 (ticket's base value); comparisons `B_vs_A` only; 60 rows not 120; per-comparison: comparisons set, C_vs_A rows 0 against 60. B_vs_A identical to its tibble-slice call |
| W4 limpa (fake `dpcDE`) | 10 / 10 | as W1 (8 failures), plus the comparisons set and C_vs_A rows 0 against 60 of W2's form. Passes: three calls completed, `da_model$model` is `dpcDE`, wide table 60 rows, B_vs_A identical |
| W5 RUV factors in the model, three groups | 9 / 5 | `comparisons` `B_vs_A` only; `length(results)` 1; `nrow(significant_rows)` n not 2n; the comparisons set; C_vs_A rows 0 against n. B_vs_A identical. If the route cannot be reached on the three-group fixture the test errors or fails earlier: that is E-4's case (record; no new machinery) |
| W7 factor column, limma | 5 / 5 | comparisons set (`C_vs_A` only); B_vs_A rows 0 against 60; the rows labelled `C_vs_A` carry B-A's statistics: `log2FC`, `raw_pvalue`, `fdr_qvalue` of C_vs_A differ from the one-comparison C_vs_A call (3 failures); row count and ids of C_vs_A pass. Mean `log2FC` of P11..P16 under `C_vs_A` printed (about 1) |
| W8 sibling, factor column | 0 / 1 | the call stops with 'subscript out of bounds' (`expect_null(run$error)` fails) |
| W9 character vector, both functions | 0 / 2 | each call stops with an error that is not of class `prot_da_contrasts_invalid` ('incorrect number of dimensions'); the class assertion fails |

**Stop rule 6 (register R01).** If W7's `log2FC`, `raw_pvalue` and `fdr_qvalue` for C_vs_A pass on base, or W8's call
completes, `purrr::map` passed labels, not factor codes: the mislabelling premise is false; recorded; W7 and W8 stay
as guards with their head expectations unchanged. A different outcome in a row (a passing row predicted to fail, an
error where a failure is predicted) is a missed prediction: stop and report.

## Predictions, `test-kit218-da-wrapper-tail.R` (G4)

0 failed, 0 errors; passed and skipped as KIT-218's record (its G1 and M2 skip because the fix is present). G4 passes.

## Predictions, `g1.R` (script)

Exit 0, writes `g1-base.tsv` with one row per digest for the cases `limma`, `limpa`, `ruv`, `kernel`, `alias`; each case has
at least ten rows (significant_rows, contrasts_results_table, wide and long tables, count tables, volcano data, the
comparison table of `contrasts_results`, the model record). The digests are the base side of G1 and G3; they decide
nothing alone.

## Head (after the fix; judged in the same file)

All 124 expectations pass (W1 10, W1b 10, W2 14, W3 17, W4 24, W5 18, W7 14, W8 15, W9 2), 0 failed. G4 passes.

## Correction of 2026-11-26 (after run 45713440865)

Predicted for the new file: 42 passed / 39 failed. Observed: 43 / 38. One expectation differs: W5's
`length(contrasts_results$results) == 2` passes on base, because `differentialAbundanceAnalysis()` keeps every
comparison in `contrasts_results` (ticket 2.1); the prediction mis-read 2.1. No expectation was changed after the
result (only a comment marking it a guard, which passes on base and head); the ticket's correction pass is not used.
Coordinator's ruling: R1's result stands and is not re-sent.

## g1.R rerun (instrument rerun 1 of 3; written before it is sent)

Purpose: the BASE digests of G1 and G3. R1's script job stopped at its first call (`ProteinQuantitativeData` is not
exported): mechanics, no digest produced. Fix: the namespace accessor in `make_object`. Binding: the commit that
carries the fix and nothing under `R/` (R/ identical to base 176f9d3ff); script only, no tests; base none. Predicted:
exit 0, `g1-base.tsv` with the cases limma, limpa, ruv, kernel, alias, at least ten rows each. Expected about 1
core-hour (one job, about 15 min); ceiling 30 min script timeout plus set-up, about 40 min = 2.67 core-hours. No
expected value, threshold or comparison changes (condition i).
