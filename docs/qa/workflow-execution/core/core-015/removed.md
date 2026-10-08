# CORE-015 removal record (written before any deletion, 2026-10-08)

Scope (user, 2026-10-08, "trimmed"):
- remove proven-dead rollup code, the experimental rollup profiles and arms, and the `multischolar.perf108.*` options, with their tests and harness selector patterns;
- add a preconditions document with classed refusals for the two opaque errors the review found;
- write a short design note.

Centralising shared constants (ticket item 3) is not done.

Reachability rule: a function or branch is removed only when every product entry point (exported functions, the GUI modules in `R/mod_*.R`, the execution registry in `R/utils_workflow_execution.R`) is shown not to reach it. Otherwise it stays.

## Removed

| Item | Where | Reachability proof | Replacement |
|---|---|---|---|
| `persistWorkflowPrecursorRollup` | `utils_workflow_dia_qc_rollup.R` | Its only call (line 152) follows `if (!is.null(persist)) return(persistIndexedWorkflowPrecursorRollup(...))` at line 94, so it would need a non-NULL `persist` that has already returned. No other call in `R/`, `tests/` or `tools/`. | none |
| `readWorkflowRollupAnnotations` | same | Only call (line 314) is inside `persistWorkflowPrecursorRollup`. | none |
| Profiles P0–P13, P15, P16 | `diaRollupProfile` and their branches in `utils_workflow_dia_relational_rollup.R` | Product callers: `computeDIAReferencePrecursorRollup` passes `"P14"`; `computeSavedWorkflowPrecursorRollup` passes `getOption("multischolar.dia_rollup_profile")`, which no product code sets (only tests and `tools/profiling`). | Any non-P14 name keeps the existing `dia_rollup_unavailable` refusal ("unsupported DIA rollup profile"). |
| Arms C and D: the DBI DuckDB connection, the loaded `dia_rsum` extension and the pinned duckplyr relational frontend | `newDIARelationalScope`, `prepareDIARelationalAggregate`, `persistDIARelationalRollup`, `diaRollupExtensionInfo` | Both product callers pass `arm = "F"`. Arms C and D are selected only in `tests/testthat` (harness-gated preflights) and `tools/profiling`. | `arm` other than `"F"` refuses with `dia_rollup_unavailable`. `validateWorkflowActionParameters` checks `diaRollupNativeInfo()`, the provider P14 actually uses, in place of `diaRollupExtensionInfo()`. |
| P5/P6/P7 intermediates, P8–P10 sample-aligned staging (`prepareDIARelationalLayout`, `stageDIARelationalCohort`), P11/P12 non-native result reads | relational rollup | Branches keyed only on those profile names. | none |
| `multischolar.perf108.ordered_groups` and `R/utils_workflow_dia_ordered_reduce.R` (`prepareDIAOrderedAggregate`, `assertDIAOrderedSchema`) | relational rollup, ordered-reduce file, `DESCRIPTION` Collate | Default `FALSE`; set only by `tools/profiling` and the `test-perf108-ordered-reduce.R` gate. | none |
| `multischolar.perf108.native_output` and `streamWorkflowQueryToWriter` | relational rollup, `persistDIAFilteredSamples` and its execution fields (qc_rollup), codec | Default `FALSE`; set only by `tools/profiling` and perf108 tests. | The product route stages output, as today. |
| `multischolar.perf108.paired_state` | `prepareDIARelationalAggregate` | Default `FALSE`; set only by `tools/profiling` and perf108 tests. The C provider keeps registering `dia_paired_r_sum`; only the R branch goes. | none |
| `multischolar.perf108.partition_validation` (per-sample visit branches) | `validateDIARelationalInput`, `estimateDIARollupWorkingMemory` | Default `FALSE`; set only by `tools/profiling` (structural V) and `test-perf108-partition-validation.R`. | The per-batch checks (CORE-017, re-applied in CORE-010) stay. `visitWorkflowCompletePartitions(batches = NULL)` stays as a generic helper. |
| `multischolar.perf108.persistent_reader` | codec reader branch, `streamDIARelationalResultFile` cursor branch | Default `FALSE`; set only by `tools/profiling` and perf108 tests. | none |
| `streamDIARelationalQuery` non-fast branches (DBI Arrow, `nanoarrow` streaming) | relational rollup | Reached only with a DBI connection (arms C and D) or the P2 profile. | fastduckdb branch only. |
| Tests for removed paths | `test-perf108-ordered-reduce.R` and `test-perf108-cursor-integration.R` (deleted); the non-P14 and arm C/D parts and the `MULTISCHOLAR_PERF108_QUERY_REPAIR` plan comparison of `test-workflow-bulk-relational.R` | They exercise only the removed options, profiles or arms. | Kept and updated: `test-perf108-output-cache.R` (the staged reader P14 still uses; now ungated, new signature), `test-perf108-partition-validation.R` (generic per-sample helper; option lines dropped), `test-perf-105.R` (now runs P14 on the native route). Kept unchanged: `test-perf108-native-null.R` and `test-perf108-paired-state.R` (they test the C provider, which still registers both sums). |
| Harness selector patterns for removed profiles | `run_workflow_checkpoint_proof.R`, `.sh` | Those selectors launch only removed profiles. | P14 selectors stay. |

## Kept, with reason

- **The non-exact DuckDB `sum()` arm** in `computeSavedWorkflowPrecursorRollup`. It is reachable from the GUI: `mod_prot_qc_peptide_rollup.R:112`, then `computeDIACapacityRollup` (`utils_workflow_dia_state.R:1626`), when `multischolar.dia_capacity_trial` is set. CORE-012 retained the ordered sum, but this experimental GUI route is app-owned and goes to the PERF-107 hand-off (CORE-016).
- **`multischolar.dia_rollup_batches`** (an explicit fixed P14 count). The installed adapters member's launch expression sets it; removing it changes member launch, which would require a lifecycle re-proof.
- **Mandatory `worker_cap` for arm F** (deferred from CORE-010) stays deferred. The harness-gated tests and unmanaged callers still call without a cap. Handed to CORE-016.
- **`multischolar.dia_rollup_profile`**, the test and harness selector for `computeSavedWorkflowPrecursorRollup`. It now accepts only `"P14"`.
- **The `dia_rsum.duckdb_extension` binary and its manifest** under `inst/duckdb`. The manifest also carries the native provider binding that `diaRollupNativeInfo` reads. Removing package files is out of this trimmed scope.
- **The perf108 `structural_options` launch branch in `test-workflow-execution-adapters.R` and the `worker_options` data in `test-perf108-structural-decision.R`** set the removed options for historical members. They are kept so member launch is unchanged; the options are no longer read.
- **`helper-workflow-bulk-rollup.R`** now maps profile-named runs to arm F and drops the arm C local switch. This is a harness file, so the bindings are re-frozen.
- **Historical PERF tooling in `tools/profiling`** (`workflow_perf106.R`, `workflow_perf108_ordered.R`, `workflow_perf108_structural.R`, `perf108-fastduckdb/freeze-structural.R`, and the PERF-108 stages of `workflow_checkpoint_dia_driver.R` and `helper-workflow-bulk-rollup.R`). It still names removed options and functions, and fails if run at HEAD. Accepted PERF-108 evidence is reproduced by checking out its commit (`9df060f0` or earlier), not by keeping dead routes in the product. This deviates from the ticket's "no reference remains in `tools/profiling`" criterion, by user scope.

## Added during implementation

- **The platform check** (Linux x86_64 with 80-bit long double) moved from the removed `diaRollupExtensionInfo` into `diaRollupNativeInfo` with the same `dia_rollup_unavailable` class. A test simulates another platform.
- **Classed refusals** were added for `workflowHostAvailable` (no `MemAvailable`) and `openWorkflowPeakWindow` (non-writable `memory.peak`), each with a test simulating its absence.
- **Exact P14 cases** were run in-process from source with the writer-r preflight gate, before and after: all five exact (ordinary 2,304 rows, empty, missing factors and annotations 2,304, boundary, threshold).
- **`duckplyr` moved from Imports to Suggests** (`R CMD check`: "Namespace in Imports field not imported from: duckplyr" after the arm D removal). It is still used by `test-workflow-memory-mechanism-member.R` and the historical driver. The version pin is unchanged.

## Verification (2026-10-08)

- **Exact P14 cases** (in-process, writer-r gate): before 153 checks passed, after 173 (checks that were conditional on perf108 environment variables now always run). All five cases are exact and identical before and after.
- **Affected suites from source:** `test-workflow-bulk-relational.R` 45, `test-perf-105.R` 62 (now on the native route), `test-perf108-output-cache.R` 26, `test-workflow-saved-qc-rollup.R` 76, `test-workflow-saved-resources.R` 302, `test-workflow-execution.R` 64, `test-perf-103.R` 85. No failures.
- **Full `tests/testthat` in one process** (`testthat::test_dir`, source, control library first). Before, at `20eae34e` in a clean worktree: 334 files, 26,533 passed, 233 failed, 107 errors. After: 332 files, 26,574 passed, 233 failed, 105 errors. No file has more failures or errors after the change; two have fewer (`test-workflow-perf101-diagnostic.R`, `test-workflow-ruv-pearson-diagnostic.R`). The failures common to both runs come from running every file in one process; the files concerned pass when run alone.
- **`R CMD check --no-manual --no-tests`** (`r-cmd-check.log`): no errors; 3 warnings and 3 notes, all present before this ticket (the `cli` `::` use, which is not in Imports at `20eae34e` either; Makevars flags and GNU make extensions in unchanged files; hidden files, long paths and the source-package note). The `duckplyr` unused-import warning this ticket introduced is fixed by moving it to Suggests.
- **Remaining references** (`grep-remaining.txt`) are only in the kept historical files listed above.
