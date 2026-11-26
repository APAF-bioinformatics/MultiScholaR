# KIT-280 step 0: what was read at the branch point, and what moved

Branch `kit-280` from `origin/v0.9-staging` e4a6d4aea (merge train 37). Read with `git grep` and the files; nothing run.

## Statements of the ticket that hold at e4a6d4aea (re-located)

| Statement | Now at |
|---|---|
| the line that makes a threshold-level value missing | `R/func_prot_qc_filtering_helpers.R:243-247` (`is_missing = if_else(!is.na(x) & x > threshold, 0, 1)`); threshold :226-233 (type 7); KIT-187 marker :258-271 |
| the normalisation declaration's `fixed` | `R/utils_workflow_dia_products.R:536` (ticket: :526): now `list(missing_value_comparison = "counts", protein_normalisation_fallback = "remove")`; both entries kept, one added |
| the filter declaration's `fixed` | `R/utils_workflow_prot_declarations.R:464` (unchanged) |
| marker table | `R/utils_workflow_analysis_record.R:703-` ; last row on staging is `ruv_fallback` (KIT-242's four rows `normalisation_fallback`, `protein_normalisation_fallback`, `itsd_outcome`, `ruv_fallback` follow `correlation_exclusion`); the new row is appended after `ruv_fallback` |
| record map | `workflowRecordReportSources`, `proteins_intensity_cutoff_percentile` entry at :1494 region; a new `threshold_values_rule` entry is added |
| display column lists | `R/utils_workflow_dia_da.R:378` (page), :416-417 (display; ticket :416), `R/utils_workflow_dia_da_views.R:249` (fallback) |
| labels | `R/utils_workflow_operations.R:575-577` |
| DA label back-mapping | `R/func_prot_da_model_methods.R:150-200, 447-453`: the helper renames digit-leading labels to `grp_...` in the design and the contrasts and maps back in `comparison` only; `createDaResultsLongFormat` (`R/func_prot_da_results_long_format.R:82-83`) takes `left_group`/`right_group` from the expression column of the safe-labelled table, so they read `grp_1h`. The flag's match therefore strips a leading `grp_` before a digit (T3d states the outcome). The reading holds; stop rule 1 is not triggered |
| `createDaResultsLongFormat` callers | `func_prot_da_model_methods.R:557`, `func_prot_s4_da_methods.R:357`, `func_prot_da_model_wrapper.R:384`, `func_prot_da_compat.R:218` |

## KIT-242's note for this ticket

The commit function's `config$outcomes <- list(samples_removed = ...)` (`utils_workflow_dia_products.R:1074-1077`) replaces the list. This ticket writes nothing into any `outcomes`: the flag is a column of `protein_id_table`, the marked count a line of the Apply text. There is nothing to merge; the line is not touched.

## What moved since the ticket was written

- KIT-242 and KIT-260 are on staging (train 37). KIT-279 is still tickets only: no Revert guard on a repeated Apply, so T5b needs none (it applies the filter once per owner).
- Regression scope: the ticket's pattern gives **104 files** at e4a6d4aea (ticket: 102 at a5de1c2eb), plus `test-omic-qc-module.R`, plus the new file = **106**. `test-kit266-limpa-combined-qc-figure.R` (1,714.9 s) is in scope.
- Timings (train 37, `merge-train-20261126aa-tests-head.tsv`): the 105 files sum to 20,707.6 s; `test-kit247-design-resave.R` 4,119.5 s alone; the other 104 plus the new file at 300 s: 16,888 s = 281.5 min per side. Five shards dealt by position (balanced by ordering the list: 56.3 test-minutes each, 67.3 with install) and one dedicated job per side for `test-kit247-design-resave.R` (timeout 120 min).
- Rule B is built exactly as section 4.1 (`is.na(x)` in the current-token branch, nothing else changed): KIT-281's gate result (VAL-378) transfers.
