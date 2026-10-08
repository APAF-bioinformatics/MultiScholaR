# CORE series hand-off to PERF-107 (app and GUI work)

Written by CORE-016. Line numbers are at the CORE closure source (`6d4f7c22`); PERF-107 rebases on the closure commit before resuming; since the post-closure CORE-018 (2026-10-09) changed `DESCRIPTION`, `NAMESPACE` and the plotting files that use Glimma, ComplexHeatmap, circlize, ggraph, GGally, ggpubr, seqinr, git2r, Biobase, cowplot and devtools (now Suggests behind `requireOptionalPackage()`), the rebase target is the `performance-v2` head that includes CORE-018, recorded in the series plan. Each item was checked against the code at that source.

## App findings excluded from CORE (review, 2026-10-07)

1. **`.GlobalEnv` assignments reachable from Shiny server code.**
   - `R/utils_workflow_saved_exports.R` lines 824, 828, 875: `filtering_progress`.
   - `R/utils_workflow_dia_sessions.R` lines 1561 and 1571: the checkpoint metadata globals restored on Open.
   - The review also named `utils_workflow_dia_qc_rollup.R` and `utils_workflow_execution.R`; neither assigns into `.GlobalEnv` at the closure source.
2. **Checkpoint Open ordering.** In `openDIAPeptideCheckpoint` (`R/utils_workflow_dia_sessions.R`), `restoreDIANativeGuiFiles(imported, paths)` runs before the stale-project check (`identical(workflow_data$state_manager, previous)`) and before the rollback handler is registered. A stale or failed Open can therefore leave restored GUI files behind.
3. **Checkpoint states and validator disagree.** `isDIAPeptideCheckpoint` accepts `raw_data_s4` and `qvalue_filtered`, which come before the design step, but `validateDIAPeptideCheckpointMetadata` requires a valid design matrix and design list. A checkpoint taken at those states cannot be reopened.
4. **Entry-fixture code in production.** `runWorkflowGUIEntryFixture` in `R/mod_prot_qc_peptide_sample.R` (from line 94) is a qualification fixture compiled into the production module. It is gated by `application$config$entry_fixture_test`.

## Found during the CORE series

5. **Store waste** (CORE-014, `core/core-014/outcome.md`).
   - Stores keep 39–58% of their bytes outside what Save archives.
   - Each Open leaves the previous session's private folder (`.dia-native-import-*` or `.dia-gui-adoption-*`, about 731 MB in the measured session) beside the project.
   - The safe fix: after a successful Open swap, once the previous manager is released and its rollback handler has exited, remove the replaced folder. A live-store clean-up becomes safe only after the app declares the references it holds (export metadata in `context$owned_export`, adoption working objects from `prepareDIAGUIReferences`).
6. **Peptide QC cost is the per-group kernel** (CORE-008, `core/core-008/outcome.md`). In `runDIAReferencePeptideGroups`, the per-protein-group scientific kernel (`dplyr::summarise` per group) took 17–33% of an imputation step and grew with input. Table appends took 1–3%.
7. **Experimental capacity trial** (CORE-015, `core/core-015/removed.md`). `multischolar.dia_capacity_trial` routes the GUI rollup module (`mod_prot_qc_peptide_rollup.R`, then `computeDIACapacityRollup`) to the non-exact DuckDB `sum()` arm. It is the last caller of that arm; CORE-012 decided the product keeps the exact ordered sum.
8. **Mandatory worker cap** (deferred from CORE-010 and CORE-015). `computeDIARelationalRollup` still accepts no `worker_cap` and then skips admission. Managed actions always pass the dispatched cap; unmanaged callers do not.
9. **Admission now follows free host memory** (CORE-010, `workflowAdmissionCap`). On a memory-short machine the P14 rollup takes more batches instead of failing. GUI messages about memory refusals may need to mention free memory as well as the 16 GiB policy.
