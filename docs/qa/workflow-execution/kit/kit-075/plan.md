# KIT-075 implementation plan (planning gate)

Written 2026-11-02 against the KIT-074 tree (a0001f4a). The ticket's line numbers predate KIT-073; the sites below were re-read on this tree.

## Order

1. **The walk first, so it fails before the change.** `tests/testthat/test-workflow-omics-residency.R`:
   - **Setup.** Per omic, a small screen corpus from `tools/profiling/kit_omics_corpus.R` (seeded recipe, in a temporary directory). Import and design go through the observers' helpers, as `kit-031/oracle-record.R` and `kit-032/oracle-record.R` call them, on the server's real `workflow_data` and owner.
   - **Journey.** Then the QC, norm, DA and summary tab servers under `testServer`, driven by `session$setInputs`. The input ids are the kit driver's (`runKitOmicsStep`).
   - **The walk, after each action.** It follows environments, lists, attributes, closures and R6 objects, stopping at namespaces and the global environment. Roots: each module's server environment, `workflow_data` and `session$userData`. The test fails on any `LipidomicsAssayData` / `MetaboliteAssayData` it reaches.
   - **Owner-lookup tests (item 1)** in the same file:
     - (a) a lipid value with no proteomics owner;
     - (b) with a proteomics owner attached.
2. **Item 1.** `R/utils_workflow_dia_state.R`: `manager = NULL` on `getWorkflowOwnedInput`, `getWorkflowOwnedPayload`, `setWorkflowOwnedPayload` and `savedReactiveVal`. NULL keeps today's lookup. The receipt probe also counts a list-of-tables value and `design_matrix`.
3. **Items 2 and 6, wrappers and revert.**
   - **Wrappers.** `runWorkflowOmicsReferenceOperation` and the core wrappers get `hydrate = FALSE`. The lipid and metab wrappers default to FALSE too, keeping the argument, because the KIT-031/032 recorders pass it. The DA wrapper drops `theObject`.
   - **Module `$object` readers** (`mod_lipid_qc_intensity_server_helpers`, `mod_lipid_qc_duplicates_server_helpers`, `mod_lipid_norm_workflow_helpers`, `mod_metab_qc_intensity`, `mod_metab_qc_duplicates_server_helpers`, `mod_metab_norm_pipeline_helpers`, `mod_metab_norm_session_helpers`): use the committed reference, or a per-action read for the QC plots.
   - **Revert.** The four GUI `revertToState` calls move to `selectWorkflowStateReference`.
4. **QC (items 4 and 5).**
   - **Metadata reads:**
     - `mod_{lipid,metab}_qc.R`: the tab's render and observe;
     - `mod_lipid_qc_s4.R` and `mod_metab_qc_s4_server_helpers.R`: the finalise, summary and stats;
     - the intensity pre-counts;
     - the duplicates owner path.
   - **Per-action reads, inside `withWorkflowStateScope`:** duplicate detection, ITSD candidate detection, and the QC plots of each action.
5. **Norm (items 3 to 5).**
   - **`normData` slots** become state references (the proteomics `correlation_filtered_obj` pattern). Set sites: `mod_lipid_norm_workflow_helpers`, `mod_metab_norm_pipeline_helpers`, `mod_metab_norm_session_helpers`.
   - **Readers** resolve per action (`mod_lipid_norm_support_helpers:264,325`, `mod_metab_norm_render_helpers:348,379`, `mod_metab_norm_session_helpers:104,260,387`).
   - **Metadata reads:** the observers' assay names.
   - **Per-action reads:** the pre-QC, the pipeline entry and the export.
   - **ITSD tables:** the closures capture the candidate table, not the S4.
   - **Norm session list:** `current_s4_object` becomes a reference.
6. **DA (item 3).**
   - **Session data:** `current_s4_object` becomes the state reference (`resolveWorkflowOmicsOwnedSession`, the readers in `mod_lipid_da_session_helpers` and `mod_metab_da_observer_helpers`).
   - **Metadata reads:** `lipidDaCurrentS4` and `metabDaCurrentS4`, and their status and preflight callers.
   - **Plots and exports:** the Glimma volcano, the heatmap and the exports' heatmaps read the displayed assay's components per action.
7. **Summary (items 4 and 5):**
   - **Metadata read:** `@args`.
   - **Per-action read:** the integration RDS.
8. **Import and design bus values (item 3).**
   - **What changes:** `workflow_data$data_tbl`, `data_cln` and `design_matrix`, the builder `reactiveVal`s and the import `localData`. They become owned values with the owner passed, on the proteomics `setDIAOwnedPayload` / `savedReactiveVal` pattern.
   - **Readers** resolve with the owner.
   - **Why it's last:** these are tables, not S4 assay objects. They are not part of the walk's assertion, so they are tested on their own.
9. **Checks.**
   - the walk and lookup tests, which fail first and then pass;
   - the affected tests, the dependents grep and R CMD check;
   - `kit_075`, the exactness set, both lifecycle members (stop if either is above kit_032), then the full suite.

## Not changed

Kernels; the declarations' ids, parameters and identities; `revertToState`'s contract; proteomics holders; saved formats.
