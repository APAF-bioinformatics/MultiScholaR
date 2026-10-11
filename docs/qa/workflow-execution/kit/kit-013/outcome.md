# KIT-013: app defects from the CORE hand-off

Status 2026-10-12: **implemented and verified in-process**; the installed-build regression (the `entry-*` and `checkpoint-*` harness scenarios through the generic entry hook) waits for the next recorded install and the user's go. Authority: `tickets/kit/00-KIT-PLAN.md` (run ledger KIT-013: no counted member; the CORE-014 read-only store probe); protocol v1.3. Control: `kit_011` (`core/browser-smoke.AFaqC1`). Predecessor cases: PERF-107 hand-off items 1 to 5 (`../../core/perf-107-handoff.md`); audit findings A4, B29, B30, B31 (`../KIT-SERIES-AUDIT-20261009.md`).

## What was built

1. **Global-environment pairs replaced by owned payloads on `workflow_data`** (hand-off item 1, audit B30, A4).
   - `filtering_progress`: `updateProteinFiltering` gains `progress_owner` (default `workflowFilteringProgressOwner()`: the bus the proteomics server registers as `session$userData$workflow_data`, else the global environment so scripts and workbooks keep their behaviour); `readWorkflowFilteringProgress`/`writeWorkflowFilteringProgress` in `R/func_prot_qc_filtering_progress_helpers.R`. `publishDIAQvalueOutput` seeds, renders and acknowledges on that owner; the Q-value, peptide-group and protein-output product captures read `workflow$filtering_progress`; the worker publication (`publishWorkflowPeptideResult`) owns a private environment for its render instead of a global binding.
   - `config_list` / `contrasts_tbl` / `uniprot_dat_cln` / `aa_seq_tbl_final`: every proteomics writer stops copying to the global environment (`mod_proteomics.R`, design state/import/builder helpers, import config/processing helpers, DA load-session handler, `openDIAPeptideCheckpoint`, normalisation helpers, summary session helper) and every reader takes the bus (`mod_prot_da_handlers_run_analysis.R`, `mod_prot_da_handlers_init.R`, `mod_prot_summary_report_helpers.R`, `mod_prot_qc_protein_cleanup.R`, the enrichment analysis/input/registration/setup helpers through `resolveProtEnrichContrastsTable`, `runProtNormBetweenSamplesStep`/`recordProtNormMethod` and `updateProtNormRuvAuditTrail` through `workflowData`, `createWorkflowArgsFromConfig` through `workflow_data$config_list`). The enrichment helpers keep their `globalEnv` parameter as a script/test seam with default `NULL`; the application never passes one. The normalisation method and RUV audit trail are now recorded on the bus configuration the study parameters are written from (before, they went into a global copy that the summary overwrote).
2. **Open ordering** (item 2): `openDIAPeptideCheckpoint` registers the rollback before the first live assignment and before `restoreDIANativeGuiFiles` runs; the stale-project check precedes the restore. `restoreDIANativeGuiFiles` returns the files it created and removes them itself when a copy fails; `discardDIANativeGuiFiles` removes them from the rollback, and the failed Open also removes the private folder it extracted the archive into (found by the new test: before, every failed Open left a `.dia-native-import-*` folder beside the project), so a stale or failed Open leaves nothing behind.
3. **States and validator** (item 3, audit B31): reproduced first (`test-workflow-native-session.R`, "a peptide checkpoint taken at every accepted state reopens"): a checkpoint saved at each of the eight accepted peptide states, `raw_data_s4` and `qvalue_filtered` included, reopens to that state with identical peptide data. The hand-off premise does not hold: `raw_data_s4` is committed by the design module after the design matrix exists (`buildProtDesignStateCheckpoint`), so the design the validator requires is present at every accepted state. **Direction recorded: no change** to `isDIAPeptideCheckpoint` or `validateDIAPeptideCheckpointMetadata`.
4. **Entry hook** (item 4, audit B29): `WorkflowExecutionApplication` carries `entry_hook`; `workflowManagedMain` installs it from the launch record through `workflowEntryHook(config)` (today the one hook is the PERF-107 small entry fixture, selected by `entry_fixture_test`); `runWorkflowGUIEntryFixture` moved from the sample module to `R/utils_workflow_execution.R`. The sample module's observer (`input$entry_hook`) runs whatever hook the application carries, validates only that the command is a plain identifier, writes the `entry-<command>.rds` receipt once and records a hook failure in it; without a hook the input is ignored. No qualification name remains in the module. The checkpoint driver's `requestPerf107EntryFixture` sets the renamed input.
5. **Replaced folder clean-up** (item 5, CORE-014): `removeDIAReplacedSessionFolder(root, keep, paths)` removes the replaced manager's private folder (`.dia-native-*`, `.dia-native-import-*`, `.dia-gui-adoption-*`, directly under the project's source directory, never the adopted owner's root and never any other folder) after a successful swap: in `openDIAPeptideCheckpoint` from an exit handler registered after the rollback (so it runs after the previous owner is closed and the rollback handler has exited), and in `adoptDIANativeSession` after adoption.

## Results

| Check | Evidence | Result |
|---|---|---|
| Open probe (CORE-014 classification before and after two Opens, screen fixture) | [open-probe.R](open-probe.R), [open-probe.json](open-probe.json) | before: one private folder (`.dia-native-…`, 325.6 MB, 266 entities, 132 outside the closure, 178.7 MB); after Open 1: only `.dia-native-import-…` (233.4 MB; 123 entities, 1 outside the closure); after Open 2: only the second import folder; the replaced folder removed each time (`replaced_folder_removed_each_time` TRUE) |
| Every accepted state reopens (item 3) | `test-workflow-native-session.R` | eight states saved and reopened, data identical |
| Failed and stale Opens leave no restored file; second Open removes the first folder; the rule refuses user folders, folders elsewhere and the adopted root (items 2, 5) | `test-workflow-native-session.R` | passes |
| Entry hook: installed from the record, run by the module, receipt once, failure recorded, bad command ignored (item 4) | `test-workflow-execution-launch-record.R`, `test-perf107-entry.R` | passes |
| Progress owned by the bus; global environment untouched; script default unchanged (item 1) | `test-prot-02r-qc-reporting-helpers-shared.R` | passes |
| Grep: proteomics `.GlobalEnv` writes/reads of the five pairs | `grep -rn "\.GlobalEnv" R/mod_prot*.R R/mod_proteomics.R R/utils_workflow_*.R R/func_prot_qc_filtering_progress_helpers.R` | none for `filtering_progress`, `config_list`, `contrasts_tbl`, `uniprot_dat_cln`, `aa_seq_tbl_final`; remaining hits are `project_dirs`/`experiment_paths` (not in scope, below), the worker-process `local_bindings` of `config_list`/`project_dirs` in `publishWorkflowPeptideResult` (kernel inputs inside the worker, no Shiny server), and the documented script default `globalenv()` |

Test files run after the change (working tree, proof library on `R_LIBS`):

| File | Pass | Fail | Error | Skip |
|---|---|---|---|---|
| test-workflow-native-session.R (two KIT-013 tests added) | 481 | 1 (baseline artefact) | 1 (same) | 1 |
| test-perf107-checkpoints.R | 49 | 0 | 0 | 0 |
| test-perf107-entry.R (hook dispatch test added) | 34 | 0 | 0 | 0 |
| test-workflow-execution-launch-record.R (hook test added) | 22 | 0 | 0 | 0 |
| test-perf107-consumers.R / test-perf107-audit.R / test-perf-107.R | 24 / 74 / 21 | 0 | 0 | 0 |
| test-prot-02r-qc-reporting-helpers-shared.R (owner test added) | 52 | 0 | 0 | 0 |
| test-prot-02e-qc-protein-cleanup-shared.R | 53 | 0 | 0 | 0 |
| test-prot-02b-qc-peptide-module-contracts.R / 02i / 02q | 772 / (02i, 02q green) | 0 | 0 | 0 |
| test-dia-baseline-compatibility.R | 42 | 0 | 0 | 0 |
| test-prot-01c / 01d import, test-module-ci-prot-import.R | 410 / 79 / 90 | 0 | 0 | 0 |
| test-prot-04a design, test-module-ci-prot-design.R | 230 / 96 | 0 | 0 | 0 |
| test-prot-05b / 05d normalisation | 562 / 199 | 0 | 0 | 0 |
| test-prot-07b-da-handlers-compat.R | 55 | 0 | 0 | 0 |
| test-prot-11 / 11a enrichment | 1253 / 143 | 0 | 0 | 0 |
| test-prot-12 / 12a summary, test-module-ci-prot-summary.R | 446 / 62 / 110 | 0 | 0 | 0 |
| test-workflow-gui-memory-lifetimes.R, test-workflow-execution.R, test-workflow-kit-harness.R | 252 / 64 / 203 | 0 | 0 | 3 / 10 / 0 |
| test-module-ci-prot-norm.R, test-module-ci-prot-enrich.R | 99 / 96 | 2 / 0 | 3 / 1 | pre-existing (below) |


Pre-existing failures, identical at the committed baseline (`32499f43`, checked in a worktree): `test-module-ci-prot-norm.R` MCI-008.2, MCI-008.4, MCI-008.6; `test-module-ci-prot-enrich.R` MCI-010.2; `test-workflow-native-session.R` "native publication recovers after process loss" (needs the installed package; recorded at KIT-001).

## Claim limits and unrelated discoveries

- The installed-build check (ticket regression: the `entry-*` and `checkpoint-*` driver scenarios through the generic hook on an installed build) has not run; it needs an install of this tree and the proof service, recorded not counted, after the user's go. Until then item 4 is verified in-process only (hook installation, module dispatch, receipts).
- Out of scope, left as found: `project_dirs` written to the global environment by `setupDirectories` and read by `updateProteinFiltering`, `publishDIAProteinProgress` and `mod_prot_summary_publication_helpers.R`; `experiment_paths` read from the global environment in `mod_prot_qc_protein_cleanup.R`; the `func_general_filemgmt_*` script helpers' global fallbacks for `contrasts_tbl`/`design_matrix`/`config_list` (the application passes these explicitly; `createWorkflowArgsFromConfig` now reads the bus before that fallback); lipidomics and metabolomics globals (stage 3).
- `updateProteinFiltering` is exported; its script behaviour (global `filtering_progress`) is unchanged by the default owner. In the application the owner is found through the Shiny session (`session$userData$workflow_data`, registered by the proteomics server), the same route `getDIAActiveManager()` uses for the state owner.
- The probe exercised the DA tab's native-session Open (`adoptDIANativeSession`) on the screen member's exported session; the peptide-checkpoint Open is covered by the test with the synthetic small fixture.

## Files

`R/func_prot_qc_filtering_progress_helpers.R`, `R/utils_workflow_saved_exports.R`, `R/utils_workflow_dia_state.R`, `R/utils_workflow_dia_qc_rollup.R`, `R/utils_workflow_execution.R`, `R/utils_workflow_dia_sessions.R`, `R/mod_proteomics.R`, `R/mod_prot_qc_peptide_sample.R`, `R/mod_prot_qc_protein_cleanup.R`, `R/mod_prot_design_state_helpers.R`, `R/mod_prot_design_import_helpers.R`, `R/mod_prot_design_builder_helpers.R`, `R/mod_prot_da_handlers_load_session.R`, `R/mod_prot_da_handlers_run_analysis.R`, `R/mod_prot_da_handlers_init.R`, `R/mod_prot_summary_report_helpers.R`, `R/mod_prot_summary_session_helpers.R`, `R/func_general_filemgmt_study_parameter_writer.R`, `R/mod_prot_enrich_analysis_helpers.R`, `R/mod_prot_enrich_input_helpers.R`, `R/mod_prot_enrich_registration_helpers.R`, `R/mod_prot_enrich_setup_helpers.R`, `R/mod_prot_norm_ruv_helpers.R`, `R/mod_prot_norm_workflow_helpers.R`, `R/mod_prot_import_config_helpers.R`, `R/mod_prot_import_processing_helpers.R`, `R/mod_prot_import_orchestration_helpers.R`, `R/mod_prot_import_server_helpers.R`, `tools/profiling/workflow_checkpoint_dia_driver.R`, `man/updateProteinFiltering.Rd`; tests as listed in the ticket.
