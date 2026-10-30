# KIT-070 outcome: the owner for lipidomics and metabolomics

Status 2026-10-30: **completed**. Code d7daa43e (item 1, the owner layer general) + a2c40f2a (item 2, both modules on the owner) + 1c5eb450 (bindings, the lipid import lane's assertion) on `performance-v2`; base 12367700 (KIT-069 audited). Install `kit_070` = `core/browser-smoke.XXcwKH` (a2c40f2a) over the kit_030 library. No counted member. Plan: [plan.md](plan.md) (the planning gate, authored at ticket start).

## What changed

**Item 1, the owner layer general (shared code, gated by the proteomics replay and records).**

| Where | Change |
|---|---|
| `workflowDeclaredClassNames(kind)` (new, class declarations) | the declared class names by kind |
| owner `setWorkflowType` | accepts `LIPIDOMICS` and `METABOLOMICS` beside the four proteomics types |
| `readWorkflowStateMetadata` | its class list is every declared assay class (was the two proteomics names) |
| `attachDIAStateSession(..., omic_type)` | the session slot `dia_capacity_manager` (what `getWorkflowActiveOwner()` and the kit readers resolve) is written and cleared by the proteomics owner only; a lipid or metab owner's flush hook checks identity against its module's bus (`workflow_data$state_manager`). Needed because one session runs one module server per selected omic (`app_server.R:589`) and the slot is single |
| `createDIASessionOwner(..., omic_type = "proteomics")` | identity `<OMIC>-pending`; the omic on the owner's context and the attach; `bindDIASessionOwnerType` carries it to the bound owner |
| `bindWorkflowOwnerTypeAtDesign(workflow_data, type, format)` (new) | binds a pending lipid/metab owner at the design step (the import of these omics writes files only); the same type again passes (a second design pass); another type aborts; a non-owner manager (a test double of the module helpers) passes through |
| `workflowRoutesByType` | `LIPIDOMICS = "lipidomics"`, `METABOLOMICS = "metabolomics"` (the kit's route names) |
| `workflowNativeStoreLabel` | `<TYPE>-native` when no vendor format is known |
| analysis record | carries `omic_type` (the owner's, `proteomics` when absent); the replay's `setupDirectories` and `<omic>_dia_baseline_smoke` alias by it; input roles `lipid_assay`, `metab_assay` accepted. Lipid and metab steps still fall to "cannot be replayed" until KIT-031/032 add their branches |
| worker publication (`publishWorkflow{Da,Peptide,Protein}Result`) | the directories indexed by `entry$omic_type` instead of the literal `$proteomics`; `workflowRequestExperimentPaths` (no entry in scope) keeps its `$proteomics` fallback (recorded) |
| harness `kitCurrentInputComponentId` | the component through the class declaration's assay; a `lists` assay gives the sorted element ids. `kitFindStore`, `kitReadHead`, `kitFingerprintState` unchanged: the `.dia-native-` prefix and `dia-state-head.rds` are kept for compatibility and the schema gate is class-free |

**Item 2, both modules on the owner.**

| Where | Change |
|---|---|
| `mod_lipidomics.R`, `mod_metabolomics.R` | `state_manager = NULL` in the bus; `createDIASessionOwner(workflow_data, experiment_paths, session, omic_type = omic_type)` once the paths are known |
| the imports | the type recorded on the pending owner is `LIPIDOMICS` / `METABOLOMICS` (the owner's names; the in-memory class's `*_standard` spellings retire with it) |
| the four design sites (builder and import, both omics) | the `WorkflowState$new()` fallbacks and their `workflowStateClass` / `stateManagerFactory` seams removed; `bindWorkflowOwnerTypeAtDesign(workflowData, <TYPE>, workflowData$data_format)` before the `*_raw_data_s4` save |
| the 26 `saveState` sites, the 4 `revertToState` sites, every `getState`/`getHistory`/`current_state` read | **unchanged in text**: on the owner `saveState` is `checkpointWorkflowExecutionResult` + `commitWorkflowReferenceState` of the declared class (the entry point proteomics' `saveDIAOwnerState` forwards to), `revertToState` reverts by the stored reference and hydrates the reverted object, `getState` hydrates the S4. The lipid and metab modules therefore read whole objects back from the store on every observer run (the memory and time cost of this stage; KIT-031/032's declarations move the kernels onto references, the KIT-048 pattern) |

Divergences from the ticket's wording, recorded: the 26 sites were not rewritten to `commitWorkflowReferenceState` calls because the owner's `saveState` *is* that call (the smallest sufficient change; the round-trip test proves every state a component checkpoint); the reverts stay on `revertToState` (by reference inside the owner; the hydrated return is ignored at all four sites, as before); the install is over the kit_030 library (the kit_069 library is the same chain). The `WorkflowState` class itself is untouched (still the embedded bookkeeping record of every owner).

## Tests

- New `tests/testthat/test-workflow-omics-owner.R` (74 expectations): the class list and routes; a lipidomics owner with its omic (identity `LIPIDOMICS-pending`, the proteomics slot untouched, the import's type on the pending owner, the design-step binding to `LIPIDSEARCH-native`, the second pass and the other-type refusal); a real owner in a temporary project for each omic: `saveState` → a component checkpoint with one table component per assay, `getState` round trip identical in every assay table and the design, a second commit, `revertToState` by name (current state, history, the reverted object), the analysis record with `omic_type` validating; both module servers under `shiny::testServer` create the owner at start (identity, omic, paths, the proteomics slot `NULL`, the KIT-030 bus registered); no `R/mod_lipid*`, `R/mod_metab*` file constructs `WorkflowState$new` or keeps the removed seams.
- Updated to the owner contract, none deleted: `test-workflow-native-session.R` untouched (proteomics owner unchanged); `test-workflow-server-shared.R` (the owner mock takes `...`); `test-metab-04` (the factory seam is gone: a missing owner is refused, a double passes through); `test-lipid-13a`, `test-module-ci-lipid-design` (the class seam removed); the module-CI design helpers hold a recording double instead of `NULL`; `test-lipid-00b/00c`, `test-metab-00b/00d`, `test-module-ci-lipid-import` (the type names).

## Results

| Check | Evidence | Result |
|---|---|---|
| Item 1 gates before item 2 (Simulation Setup) | [replay-check-kit070-item1.json](replay-check-kit070-item1.json); [compare-lfq-tmt-item1.txt](compare-lfq-tmt-item1.txt), [compare-limpa-item1.txt](compare-limpa-item1.txt) ([compare-oracles.R](compare-oracles.R) against `kit-065/oracles-*.rds`); [item1-shared-tests-summary.txt](item1-shared-tests-summary.txt) | KIT-028 replay **exact, 16 states, 0 mismatches, 291 s**; LFQ (11 steps), TMT (11) and DIA_LIMPA (14) records **identical** to the kit_065 oracles (0 reproduce mismatches, 0 per-state mismatches); 27 shared test files 4,888 pass, 0 fail |
| Closure gates | [replay-check-kit070-final.json](replay-check-kit070-final.json); [compare-records-final.txt](compare-records-final.txt) | replay **exact, 16 states, 286 s**; the three records **identical** again |
| Lipid/metab suite (125 files incl. the new one), each in its own process | [sweep-lipid-metab-final.txt](sweep-lipid-metab-final.txt) | **8,134 pass, 0 fail, 0 errors**, 19 skips |
| Module-CI 12 files before (12367700) / after | [module-ci-comparison.txt](module-ci-comparison.txt) | identical counts (one assertion updated to the type name); 230 paired paths 193 byte-identical; the others the `pkgload` object, the timestamp-bearing session RDS, workbooks and `study_parameters.txt`, all identical at content level except the provenance block (detached worktree vs working tree); 9 before-only paths are a `devtools_install_*` temp directory of the before process. The lanes drive doubles, so this proves the module code around the doubles, not the owner (the members do) |
| The lipid lifecycle member on kit_070 | `kit-lifecycle-lipid-stage.json` | **`lifecycle_complete_observed`**, 242 s, peak 2.46 GiB anonymous, max sampler interval 0.07 s, no gaps. **A store exists**: `scripts/lipidomics/.dia-native-*/.workflow-entities`, found by `kitFindStore`; head `lipid_correlation_filtered`; history `initial > lipid_raw_data_s4 > lipid_intensity_filtered > lipid_duplicates_resolved > lipid_qc_complete > lipid_itsd_norm > lipid_log2 > lipid_normalized > lipid_ruv_corrected > lipid_correlation_filtered` (rows 432 > 392 > 377 > 377 > 365 …); the current state's fingerprint class `LipidomicsAssayData`, components `lipid_data/LCMS_Pos`, `lipid_data/LCMS_Neg`, `design_matrix`; `kitCurrentInputComponentId` returns the two assay components' ids. The DA tab's session restore re-commits under the exported state's name, so the history has no `loaded_for_de` (as the in-memory manager behaved) |
| The metab lifecycle member on kit_070 | `kit-lifecycle-metab-stage.json` | **`lifecycle_complete_observed`**, 234 s, peak 2.38 GiB anonymous, max sampler interval 0.06 s, no gaps. A store exists under `scripts/metabolomics/.dia-native-*`; head `metab_correlation_filtered`; history `initial > metab_raw_data_s4 > metab_intensity_filtered > metab_duplicates_resolved > metab_qc_complete > metab_itsd_norm > metab_log2 > metab_normalized > metab_ruv_corrected > metab_correlation_filtered` (rows 612 > 576 > 576 > 576 > 564 …); fingerprint class `MetaboliteAssayData`, components `metabolite_data/LCMS_Pos`, `metabolite_data/GCMS`, `design_matrix`. The first launch of this member was refused on a dirty tree (an evidence file filed mid-run; the same slip as KIT-069's), committed and relaunched; no application attempt was lost |
| Member projects against the kit_069 members' (same corpora, same journey) | [compare-content.R](../kit-030/compare-content.R), [compare-groups.R](../kit-069/compare-groups.R) | lipid: every DA long table, publication workbook, NumSigDeMolecules table and `study_parameters.txt` identical; the exported session RDS identical with its `timestamp` fields removed; the `.docx`, PDFs, debug log by embedded times, the volcano PNGs by `ggrepel`, the session summary by the dated file name; the store directory is new. Metab: the same, 91 paired files, every data product identical, the store directory new ([member-project-comparison-kit069-vs-kit070.txt](member-project-comparison-kit069-vs-kit070.txt)). |
| R CMD check (`--no-tests --no-manual --no-build-vignettes`, `_R_CHECK_FORCE_SUGGESTS_=false`) on the a2c40f2a export | `<proof>/kit/scratch-070-check` | **at baseline**: 3 WARNINGs, 3 NOTEs |
| Harness test on the regenerated bindings | `test-workflow-kit-harness.R` | 483 pass |
| Backward compatibility | read + the gates | proteomics behaviour unchanged (replay and records identical); the head file and `.dia-native-` prefix unchanged; public API and S4 classes unchanged; older lipid/metab RDS sessions are still read by the DA tab's Open (it re-saves the object onto the owner; adoption proper is KIT-071) |
| Rollback | per-item `git revert` (d7daa43e, a2c40f2a, 1c5eb450) | stores written by kit_070's lipid/metab modules are not read by the in-memory modules of older builds (recorded) |

## Claim limits

- The lifecycle members are in-memory-rule members (`exact = NA`, no oracle): they prove the journeys complete through the report on the owner with a store the harness reads, and that the data products match the kit_069 members'; they do not prove state-by-state exactness against a record (KIT-031/032's recorders).
- Memory: the modules hydrate whole S4 objects from the store on every read; the members' peaks (2.46 GiB lipid) are reported, not claimed as an improvement (KIT-069's lipid member peaked at 2.21 GiB in memory).
- `getWorkflowActiveOwner()` still resolves the proteomics owner only; lipid and metab kit readers reach their owner through the per-omic bus. KIT-031/032 decide whether the active-owner resolution takes an omic.

## Hand-off

KIT-071 takes the owner in both modules (owned exports, adoption of legacy RDS sessions, the summary's final state from the head). KIT-031 meets a general owner: every lipid state a component checkpoint of the declared class with the store the harness snapshots.
