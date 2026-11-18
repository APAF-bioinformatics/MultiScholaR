# KIT-146 debug-output inventory

Base: performance-v2 `b35d8089`, 2026-11-19. Scope: `R/` and `inst/` files ending `.R`/`.Rmd`/`.rmd`.

## Patterns

| Key | Pattern | How counted |
|---|---|---|
| A | `DEBUG66` / `debug66` | every occurrence (also comments, a roxygen line and one default argument); "A in code" = non-comment lines |
| B | other debug tags: `DEBUG_66`, `[D66]`, `d66_log(...)`, `D66 DEBUG`, `DEBUG_STR`, `DEBUG:` inside `message`/`cat`/`print` | non-comment lines without A |
| C | entering/exiting traces: `--- Entering`, `--- Exiting`, `=== ENTER`, `=== EXIT` inside `message`/`cat` | non-comment lines without A |
| D | console dumps of internal objects: statement `str(...)` / `utils::str(...)`, `print(head|str|names|colnames|summary|dim|glimpse(...))` | non-comment lines |

Searched and absent: other numbered tags (`DEBUG[0-9]+` other than 66: none), `print("---` (none), `browser()` (none).

Counted but not part of this pass (needs a scope decision; tests assert some of them, e.g. `test-prot-04-design.R` expects the `mod_prot_design_server Arg:`/`Step:` lines):

- `Step:` / `Arg:` / `Data State` trace messages without a debug tag: 306 lines in 24 `R/` files.
- untagged `print(...)` of other values (e.g. `print(rowinfo_vector)`, `print("Starting enrichment analysis")`): not inventoried line by line.
- commented-out debug code (`# message("DEBUG: ...")` etc.): 19 lines, out of scope.

## Open branches checked

`git diff --name-only performance-v2...<branch>` for `kit-129-integration` (33cfc99b, 169 files), `kit-052-win-paths` (463abfb9, 8), `kit-145-loess-minimum` (5ce114f9, 18), `kit-128` (ace1f4de, 31), `kit-126-audit` (95f9d1d3, 10), `kit-127-audit` (b3bd74c3, 9), `kit-134` (77c6e79d, 58). A file any of them changes is held for the second part.

## Before counts per file

| File | A DEBUG66 (all) | A in code | B other tags | C enter/exit | D dumps | Open branches touching it | Pass |
|---|---:|---:|---:|---:|---:|---|---|
| `R/app_server.R` | 0 | 0 | 0 | 2 | 0 | none | first part |
| `R/func_general_enrichment_core.R` | 0 | 0 | 0 | 0 | 1 | none | first part |
| `R/func_general_enrichment_engine_helpers.R` | 5 | 4 | 0 | 4 | 5 | none | first part |
| `R/func_general_enrichment_orchestration.R` | 14 | 9 | 7 | 0 | 12 | none | first part |
| `R/func_general_filemgmt_config_io_helpers.R` | 23 | 23 | 0 | 0 | 0 | none | first part |
| `R/func_general_filemgmt_path_helpers.R` | 12 | 12 | 0 | 0 | 0 | none | first part |
| `R/func_general_filemgmt_report_render_helpers.R` | 34 | 34 | 0 | 0 | 0 | none | first part |
| `R/func_general_plotting_pca_rle_helpers.R` | 19 | 18 | 0 | 0 | 0 | none | first part |
| `R/func_general_plotting_volcano_glimma_helpers.R` | 0 | 0 | 0 | 2 | 0 | none | first part |
| `R/func_general_runtime_helpers.R` | 2 | 1 | 0 | 0 | 0 | none | first part |
| `R/func_github_management.R` | 0 | 0 | 1 | 0 | 0 | none | first part |
| `R/func_lipid_da_model.R` | 0 | 0 | 23 | 0 | 0 | kit-127-audit, kit-128, kit-129-integration | second part, after merge |
| `R/func_lipid_da_plotting.R` | 0 | 0 | 30 | 0 | 0 | none | first part |
| `R/func_lipid_norm_ruv_optimization_helpers.R` | 17 | 17 | 0 | 0 | 0 | none | first part |
| `R/func_lipid_s4_correlation_methods.R` | 28 | 28 | 0 | 0 | 0 | none | first part |
| `R/func_lipid_s4_da_analysis_methods.R` | 0 | 0 | 0 | 3 | 0 | none | first part |
| `R/func_lipid_s4_plotting_methods.R` | 0 | 0 | 0 | 1 | 0 | none | first part |
| `R/func_lipid_s4_results_accessors.R` | 0 | 0 | 0 | 0 | 1 | none | first part |
| `R/func_lipid_s4_ruv_control_methods.R` | 46 | 46 | 0 | 0 | 0 | none | first part |
| `R/func_lipid_s4_ruv_methods.R` | 93 | 93 | 0 | 0 | 0 | none | first part |
| `R/func_metab_da_core_engine.R` | 0 | 0 | 9 | 0 | 0 | kit-128, kit-129-integration | second part, after merge |
| `R/func_metab_da_orchestration.R` | 0 | 0 | 21 | 1 | 0 | kit-127-audit, kit-128, kit-129-integration | second part, after merge |
| `R/func_metab_da_results_long_format.R` | 0 | 0 | 16 | 1 | 0 | none | first part |
| `R/func_metab_da_volcano_glimma.R` | 0 | 0 | 39 | 0 | 0 | none | first part |
| `R/func_metab_norm.R` | 17 | 17 | 0 | 0 | 0 | none | first part |
| `R/func_metab_s4_da_methods.R` | 0 | 0 | 0 | 3 | 0 | none | first part |
| `R/func_metab_s4_da_results.R` | 0 | 0 | 0 | 0 | 1 | none | first part |
| `R/func_metab_s4_design_control_methods.R` | 46 | 46 | 0 | 0 | 0 | none | first part |
| `R/func_metab_s4_plotting_methods.R` | 0 | 0 | 0 | 1 | 0 | none | first part |
| `R/func_metab_s4_qc_methods.R` | 28 | 28 | 0 | 0 | 0 | none | first part |
| `R/func_metab_s4_ruv_methods.R` | 93 | 93 | 0 | 0 | 0 | none | first part |
| `R/func_multiomics_enrich_kegg_helpers.R` | 0 | 0 | 20 | 3 | 17 | none | first part |
| `R/func_multiomics_enrich_metabolomics_helpers.R` | 0 | 0 | 0 | 4 | 5 | none | first part |
| `R/func_multiomics_enrich_reactome_helpers.R` | 0 | 0 | 7 | 2 | 7 | none | first part |
| `R/func_multiomics_enrich_stringdb_reporting_helpers.R` | 0 | 0 | 2 | 0 | 0 | none | first part |
| `R/func_omics_da_list_methods.R` | 0 | 0 | 0 | 0 | 1 | none | first part |
| `R/func_omics_design_shared.R` | 1 | 0 | 0 | 0 | 0 | kit-128, kit-129-integration | second part, after merge |
| `R/func_pept_s4_qc_methods.R` | 2 | 2 | 0 | 0 | 0 | none | first part |
| `R/func_prot_annotation_accession_ranking_helpers.R` | 4 | 4 | 0 | 0 | 0 | none | first part |
| `R/func_prot_annotation_ensembl_helpers.R` | 0 | 0 | 0 | 0 | 1 | none | first part |
| `R/func_prot_annotation_uniprot_batch_helpers.R` | 0 | 0 | 0 | 0 | 2 | none | first part |
| `R/func_prot_annotation_uniprot_full_helpers.R` | 6 | 6 | 0 | 0 | 4 | none | first part |
| `R/func_prot_annotation_uniprot_helpers.R` | 27 | 27 | 0 | 0 | 5 | none | first part |
| `R/func_prot_da_model_methods.R` | 38 | 38 | 0 | 4 | 6 | kit-126-audit, kit-128, kit-129-integration | second part, after merge |
| `R/func_prot_da_model_stats.R` | 0 | 0 | 0 | 2 | 0 | kit-128, kit-129-integration | second part, after merge |
| `R/func_prot_da_model_wrapper.R` | 1 | 1 | 22 | 0 | 5 | kit-126-audit, kit-129-integration | second part, after merge |
| `R/func_prot_da_results_io.R` | 0 | 0 | 29 | 0 | 0 | none | first part |
| `R/func_prot_da_results_long_format.R` | 8 | 7 | 0 | 0 | 0 | none | first part |
| `R/func_prot_da_results_methods.R` | 0 | 0 | 0 | 2 | 0 | none | first part |
| `R/func_prot_da_results_summary.R` | 6 | 6 | 0 | 0 | 4 | none | first part |
| `R/func_prot_import_tmt.R` | 0 | 0 | 0 | 4 | 0 | none | first part |
| `R/func_prot_limpa_qc_helpers.R` | 0 | 0 | 0 | 2 | 0 | none | first part |
| `R/func_prot_norm_optimization_helpers.R` | 31 | 31 | 0 | 0 | 0 | kit-126-audit, kit-128, kit-129-integration | second part, after merge |
| `R/func_prot_qc_correlation_helpers.R` | 2 | 2 | 0 | 0 | 0 | none | first part |
| `R/func_prot_qc_filtering_helpers.R` | 7 | 7 | 0 | 0 | 0 | none | first part |
| `R/func_prot_qc_filtering_progress_helpers.R` | 2 | 1 | 0 | 0 | 0 | kit-052-win-paths | second part, after merge |
| `R/func_prot_qc_peptide_group_filters.R` | 6 | 6 | 0 | 0 | 0 | none | first part |
| `R/func_prot_qc_peptide_methods.R` | 0 | 0 | 0 | 2 | 0 | none | first part |
| `R/func_prot_s4_accession_methods.R` | 3 | 3 | 0 | 0 | 0 | none | first part |
| `R/func_prot_s4_da_methods.R` | 37 | 37 | 0 | 2 | 6 | kit-126-audit, kit-129-integration | second part, after merge |
| `R/func_prot_s4_da_results.R` | 0 | 0 | 0 | 2 | 0 | none | first part |
| `R/func_prot_s4_norm_controls.R` | 7 | 7 | 0 | 0 | 0 | none | first part |
| `R/func_prot_s4_norm_methods.R` | 7 | 7 | 0 | 0 | 0 | none | first part |
| `R/func_prot_s4_qc_correlation_methods.R` | 3 | 3 | 0 | 0 | 0 | none | first part |
| `R/func_prot_s4_qc_filter_methods.R` | 63 | 62 | 0 | 0 | 0 | none | first part |
| `R/func_prot_s4_qc_plot_methods.R` | 0 | 0 | 0 | 0 | 1 | none | first part |
| `R/mod_omic_design_builder.R` | 1 | 0 | 0 | 0 | 0 | kit-128, kit-129-integration | second part, after merge |
| `R/mod_omic_enrich.R` | 0 | 0 | 0 | 3 | 0 | none | first part |
| `R/mod_prot_da_handlers_init.R` | 0 | 0 | 0 | 4 | 4 | none | first part |
| `R/mod_prot_design.R` | 0 | 0 | 0 | 1 | 0 | none | first part |
| `R/mod_prot_design_state_helpers.R` | 6 | 6 | 0 | 0 | 0 | none | first part |
| `R/mod_prot_import_server.R` | 0 | 0 | 0 | 1 | 0 | none | first part |
| `R/mod_prot_norm_correlation_helpers.R` | 1 | 1 | 0 | 0 | 0 | none | first part |
| `R/mod_prot_norm_render_helpers.R` | 1 | 1 | 0 | 0 | 0 | none | first part |
| `R/mod_prot_norm_ruv_helpers.R` | 0 | 0 | 0 | 0 | 1 | none | first part |
| `R/utils_workflow_omics_design.R` | 1 | 0 | 0 | 0 | 0 | none | first part |
| `inst/reports/lipidomics/lipidomics_report.rmd` | 0 | 0 | 39 | 0 | 0 | none | held: writes D66_debug.log, stored gui-entry oracles expect it |
| `inst/reports/metabolomics/metabolomics_report.rmd` | 0 | 0 | 39 | 0 | 0 | none | held: writes D66_debug.log, stored gui-entry oracles expect it |
| `inst/workbooks/metabolomics/metabolomics_workflow_experienced.rmd` | 0 | 0 | 0 | 0 | 3 | none | out of scope (workbook narrative output) |
| `inst/workbooks/metabolomics/metabolomics_workflow_starter.rmd` | 0 | 0 | 0 | 0 | 1 | none | out of scope (workbook narrative output) |
| `inst/workbooks/proteomics/DIA_workflow_experienced.rmd` | 0 | 0 | 0 | 0 | 1 | none | out of scope (workbook narrative output) |
| `inst/workbooks/proteomics/DIA_workflow_limpa.rmd` | 0 | 0 | 0 | 0 | 5 | none | out of scope (workbook narrative output) |
| `inst/workbooks/proteomics/DIA_workflow_limpa_experienced.rmd` | 0 | 0 | 0 | 0 | 1 | none | out of scope (workbook narrative output) |
| `inst/workbooks/proteomics/DIA_workflow_limpa_starter.rmd` | 0 | 0 | 0 | 0 | 1 | none | out of scope (workbook narrative output) |
| `inst/workbooks/proteomics/DIA_workflow_starter.rmd` | 0 | 0 | 0 | 0 | 1 | none | out of scope (workbook narrative output) |
| **Total (85 files)** | 748 | 734 | 304 | 56 | 102 | | |

## Outcome, first part (2026-11-19, branch kit-146-debug-cleanup)

Removed in the 54 `R/` files no open branch changes: whole debug statements only (A, B, C), found with the R parser so multi-line calls go as a unit; the boxed banners around DEBUG66 headers; the emptied `[D66:START]`/`[D66:END]` blocks and the local `d66_log` helpers once nothing called them; `if`/`for`/`else` bodies that held nothing but removed debug lines (conditions without side effects); eight `DEBUG66` label comments over removed code. A statement was removed only when it was not the value of its block, or the block's value is unused (loop bodies, `if` statements followed by more code, observer bodies). 1,182 lines removed in total, in 52 files (after restoring two test-expected lines in two files).

Not done in this pass: D (console dumps) is left whole, together with the untagged `Step:`/`Arg:`/`Data State` traces: most dumps sit under an untagged label line of that family (e.g. `message("Data State (up_matrix) Head:")`), so removing only the dump leaves the label. Needs a scope decision.

### Counts before -> after

| Group | A all | A in code | B | C | D |
|---|---|---|---|---|---|
| First part (cleaned now) | 637 -> 16 | 626 -> 13 | 151 -> 5 | 47 -> 2 | 72 -> 72 |
| Second part (held, open branch) | 111 -> 111 | 108 -> 108 | 75 -> 75 | 9 -> 9 | 17 -> 17 |
| inst/reports (held) | 0 -> 0 | 0 -> 0 | 78 -> 78 | 0 -> 0 | 0 -> 0 |
| inst/workbooks (out of scope) | 0 -> 0 | 0 -> 0 | 0 -> 0 | 0 -> 0 | 13 -> 13 |
| **All** | 748 -> 127 | 734 -> 121 | 304 -> 158 | 56 -> 11 | 102 -> 102 |

### Kept in first-part files, with the reason

| File:line (after) | Line | Reason |
|---|---|---|
| `R/func_prot_annotation_uniprot_helpers.R:49` | `DEBUG66: Cache hit, loading from cache` | `helper-e2e-browser.R` (`e2e_assert_no_uniprot_network`) reads this text from the app log |
| `R/func_prot_annotation_uniprot_helpers.R:55` | `DEBUG66: Cache miss, calling directUniprotDownload...` | same test pair (`test-e2e-browser-harness.R`) |
| `R/func_metab_s4_plotting_methods.R:197` | `--- Entering plotRle for MetaboliteAssayData ---` | expected by `test-metab-01h-s4-plot-rle-characterization.R`, `test-metab-s4-plotting-direct-shared.R` (the lipid twin was removed: no test expects it) |
| `R/mod_prot_design.R:142` | `--- Entering mod_prot_design_server ---` | expected by `test-prot-04-design.R` (call log) |
| `R/mod_prot_norm_render_helpers.R:112` | `messageFn("--- DEBUG66 [final_qc_plot]: Rendering ---")` | restored after runner run 44366929002: `test-prot-05b-norm-module-contracts.R` expects a message containing "Rendering" |
| `R/mod_prot_norm_correlation_helpers.R:139` | `messageFn(sprintf("--- DEBUG66 [mod_prot_norm]: %s ---", startMessage))` | restored after the same run: `test-prot-05b-norm-module-contracts.R` expects the caller's `startMessage` ("Starting Correlation Filter Flow") in the messages |
| `R/func_prot_qc_filtering_helpers.R:175,203,259` | `DEBUG66: ... Memory: %.1f MB` | the message arguments call `gc()`; removing them changes garbage-collection timing and so memory peaks in the KIT memory measurements; needs a decision (drop, or keep an explicit `gc()`) |
| `R/func_general_runtime_helpers.R:105,108` | roxygen line and default `log_level = "DEBUG66"` of exported `checkMemoryBoth()` | public signature; second part decides |
| `R/func_prot_da_results_long_format.R:76` | comment over commented-out prints | commented-out debug block, out of scope |
| `R/utils_workflow_omics_design.R:10` | provenance comment naming `DEBUG66` | records a KIT-103 divergence; reword in the second part if the acceptance grep must be empty |
| `R/func_github_management.R:108` | `stop("DEBUG: Failed to retrieve GitHub token ...")` | a user-facing error, not output; wording left |

### Needs a proper message (left as is; the debug line is the only failure signal)

| File:line (after) | Line |
|---|---|
| `R/func_general_filemgmt_config_io_helpers.R:493` | `DEBUG66: Error converting '%s': %s` (handler while formatting the config) |
| `R/func_general_plotting_pca_rle_helpers.R:181-184` | `DEBUG66 [plotPcaHelper] mixOmics::pca() failed (%s); falling back to stats::prcomp()` |
| `R/func_multiomics_enrich_kegg_helpers.R:611` | `runKeggEnrichment DEBUG: First join attempt failed, trying second format...` |
| `R/func_multiomics_enrich_kegg_helpers.R:620` | `runKeggEnrichment DEBUG: Second join attempt failed, trying numeric part matching...` |
| `R/func_prot_annotation_uniprot_full_helpers.R:325` | `=== DEBUG66: ERROR in getUniprotAnnotationsFull ===` (with the `print(str(e))` after it) |
| `R/func_prot_da_results_io.R:259` | `DEBUG: ERROR in wide_annot creation:` |
| `R/func_prot_da_results_io.R:288` | `DEBUG: da_proteins_long is NULL - cannot create long_annot` (before `return(NULL)`) |
| `R/func_prot_da_results_io.R:329` | `DEBUG: long_annot is empty or null - annotation pipeline failed` |
| `R/func_prot_s4_qc_filter_methods.R:213` | `DEBUG66 WARNING: theObject@args is NULL!` |
| `R/func_prot_s4_qc_filter_methods.R:230` | `DEBUG66 WARNING: ruv_grouping_variable is NULL!` |

Unrelated, seen and left: in `R/func_prot_da_results_io.R` the `tryCatch` error handler assigns the fallback `da_proteins_wide_annot` inside the handler function, so the fallback never reaches the caller's frame (pre-existing).

Held for the second part (an open branch changes the file): `func_lipid_da_model.R` (127-audit, 128, 129), `func_metab_da_core_engine.R` (128, 129), `func_metab_da_orchestration.R` (127-audit, 128, 129), `func_omics_design_shared.R` (128, 129), `func_prot_da_model_methods.R` (126-audit, 128, 129), `func_prot_da_model_stats.R` (128, 129), `func_prot_da_model_wrapper.R` (126-audit, 129), `func_prot_norm_optimization_helpers.R` (126-audit, 128, 129), `func_prot_qc_filtering_progress_helpers.R` (052-win-paths), `func_prot_s4_da_methods.R` (126-audit, 129), `mod_omic_design_builder.R` (128, 129). Held separately: `inst/reports/{lipidomics,metabolomics}/*_report.rmd` write `D66_debug.log` into the report source directory and the stored gui-entry oracles (`tests/testdata/gui-entry/{metab_gc,metab_lc,metab_combined,lipid_canonical}.json`) list that file as present, so removing it changes a stored reference; needs a decision and a re-record.

### Verification

Local: every touched file parses (`parse()`, one at a time). No local tests, build, check or install were run (user order 2026-11-19: compute on GitHub runners). Required on a runner, against performance-v2 `b35d8089` as the base for any failure:

1. Build, `R CMD check` with `_R_CHECK_FORCE_SUGGESTS_=false`, install into a scratch library; compare findings with the base.
2. The 187 test files below (union of the files that name a function from a touched file), one R process per file; every file passing on the base must pass here. `func_prot_s4_da_results.R` has no test that names its functions.

Per touched file:

| Touched file | Test files |
|---|---|
| `R/app_server.R` | test-app-shell-shared.R test-module-ci-trx-app.R |
| `R/func_general_enrichment_engine_helpers.R` | test-general-enrichment-shared.R test-omic-enrich-module.R test-omics-enrich-core.R test-workflow-class-declarations.R |
| `R/func_general_enrichment_orchestration.R` | test-general-enrichment-shared.R test-module-ci-prot-enrich.R test-omic-enrich-module.R test-prot-11a-enrichment-module-characterization.R test-prot-11-enrichment-module-contracts.R |
| `R/func_general_filemgmt_config_io_helpers.R` | test-general-filemgmt-config-contracts.R test-general-filemgmt-config-direct-shared.R test-general-filemgmt-config-roundtrip-contracts.R test-module-ci-prot-peptide-qc.R test-prot-01d-import-module-characterization.R test-prot-02c-qc-peptide-intensity-shared.R test-prot-02d-qc-protein-intensity-shared.R test-prot-02g-qc-peptide-qvalue-shared.R test-prot-02h-peptide-audit-shared.R test-prot-02i-qc-peptide-sample-shared.R test-prot-02j-qc-peptide-protein-shared.R test-prot-02k-qc-peptide-impute-shared.R test-prot-02r-qc-reporting-helpers-shared.R test-workflow-audit-memory-member.R test-workflow-checkpoint-dia-state.R test-workflow-checkpoint-harness.R test-workflow-checkpoints.R test-workflow-gui-capacity-verifier.R test-workflow-gui-diagnostic.R test-workflow-gui-memory-lifetimes.R test-workflow-gui-memory-repairs.R test-workflow-memory-mechanism-member.R test-workflow-native-session.R test-workflow-prot-declarations.R |
| `R/func_general_filemgmt_path_helpers.R` | test-general-filemgmt-path-contracts.R test-general-filemgmt-report-helpers-shared.R test-multiomics-enrich-shared.R test-omics-summary-core.R test-prot-07d-da-model-stats-shared.R test-prot-da-compat-shared.R |
| `R/func_general_filemgmt_report_render_helpers.R` | test-general-enrichment-shared.R test-general-filemgmt-export-contracts.R test-general-filemgmt-report-helpers-shared.R test-lipid-14a-summary-module-characterization.R test-lipid-15-summary-module-seams.R test-metab-03b-norm-module-helper-characterization.R test-metab-03d-norm-module-characterization.R test-metab-05-summary-module-characterization.R test-omic-summary-module.R test-prot-07d-da-model-stats-shared.R test-prot-12-summary-module-contracts.R test-workflow-omics-figures.R test-workflow-saved-qc.R |
| `R/func_general_plotting_pca_rle_helpers.R` | test-general-pca-rle-direct-shared.R test-general-plotting-pca-helper-characterization.R test-lipid-10-s4-plotting-methods.R test-metab-01g-s4-plot-pca-characterization.R test-metab-01h-s4-plot-rle-characterization.R test-metab-s4-plotting-direct-shared.R test-optional-packages.R test-prot-04d-peptide-s4-normalization-characterization.R test-prot-s4-qc-direct-shared.R |
| `R/func_general_plotting_volcano_glimma_helpers.R` | test-general-volcano-glimma-direct-shared.R test-lipid-07-s4-da-plot-methods.R test-metab-02l-da-volcano-plot-characterization.R test-module-ci-prot-da.R test-optional-packages.R test-prot-07c-da-results-output-shared.R test-prot-07-da-analysis.R test-prot-07d-da-model-stats-shared.R test-prot-08-volcano.R test-prot-da-compat-shared.R |
| `R/func_github_management.R` | test-github-management-direct-shared.R test-omic-summary-module.R test-optional-packages.R |
| `R/func_lipid_da_plotting.R` | test-glimma-plot.R test-lipid-02ab-da-module-characterization.R test-lipid-da-export-shared.R test-module-ci-lipid-da.R test-omic-da-module.R test-omic-gui-config.R test-omics-audit-carry-overs.R test-omics-da-core.R test-workflow-omics-application-reads.R test-workflow-omics-figures.R |
| `R/func_lipid_norm_ruv_optimization_helpers.R` | test-lipid-12a-norm-module-characterization.R test-lipid-12b-norm-ruv-helpers-shared.R test-module-ci-lipid-norm.R test-omic-norm-module.R test-omics-norm-core.R test-workflow-omics-ruv-cancor-once.R test-workflow-omics-ruv-complete-controls.R |
| `R/func_lipid_s4_correlation_methods.R` | test-lipid-05-s4-correlation-helpers.R test-lipid-12a-norm-module-characterization.R test-metab-01j-s4-pearson-correlation-characterization.R test-metab-01k-s4-plot-pearson-characterization.R test-metab-03d-norm-module-characterization.R test-metab-norm-standalone-shared.R test-metab-s4-qc-shared.R test-module-ci-lipid-norm.R test-module-ci-lipid-qc.R test-module-ci-metab-norm.R test-module-ci-metab-qc.R test-prot-04-design.R test-prot-04d-peptide-s4-normalization-characterization.R test-prot-05d-norm-module-characterization.R test-prot-s4-qc-direct-shared.R test-workflow-native-session.R test-workflow-omics-figures.R test-workflow-omics-pearson-matrix.R test-workflow-ruv-pearson-diagnostic.R |
| `R/func_lipid_s4_da_analysis_methods.R` | test-lipid-08-s4-da-analysis-methods.R test-metab-02p-da-analysis-characterization.R test-metab-02q-da-analysis-helper-characterization.R test-metab-s4-da-methods-direct-shared.R test-workflow-saved-da.R |
| `R/func_lipid_s4_plotting_methods.R` | test-general-filemgmt-config-contracts.R test-lipid-10-s4-plotting-methods.R test-metab-01g-s4-plot-pca-characterization.R test-metab-01h-s4-plot-rle-characterization.R test-metab-01i-s4-plot-density-characterization.R test-metab-03d-norm-module-characterization.R test-metab-norm-standalone-shared.R test-metab-s4-plotting-direct-shared.R test-module-ci-lipid-qc.R test-module-ci-metab-qc.R test-module-ci-prot-norm.R test-module-ci-prot-peptide-qc.R test-module-ci-prot-protein-qc.R test-omic-norm-module.R test-omics-norm-core.R test-prot-04c-peptide-plot-density-ggplot-shared.R test-prot-04-design.R test-prot-04d-peptide-s4-normalization-characterization.R test-prot-05d-norm-module-characterization.R test-prot-07-da-analysis.R test-prot-07d-da-model-stats-shared.R test-prot-da-compat-shared.R test-prot-s4-qc-direct-shared.R test-workflow-gui-memory-repairs.R test-workflow-native-session.R test-workflow-omics-figures.R test-workflow-saved-qc.R |
| `R/func_lipid_s4_ruv_control_methods.R` | test-lipid-09-s4-normalization-methods.R test-lipid-12b-norm-ruv-helpers-shared.R test-lipid-12c-ruv-s4-characterization.R test-metab-02y-get-neg-ctrl-metab-anova-characterization.R test-metab-04-ruv-optimization.R test-metab-norm-standalone-shared.R |
| `R/func_lipid_s4_ruv_methods.R` | test-general-filemgmt-config-contracts.R test-general-filemgmt-config-direct-shared.R test-lipid-09-s4-normalization-methods.R test-lipid-12a-norm-module-characterization.R test-lipid-12b-norm-ruv-helpers-shared.R test-lipid-12c-ruv-s4-characterization.R test-metab-02z-ruv-cancor-characterization.R test-metab-03a-ruviii-c-varying-characterization.R test-metab-04-ruv-optimization.R test-metab-norm-standalone-shared.R test-module-ci-metab-norm.R test-omic-norm-module.R test-omics-norm-core.R test-prot-04d-peptide-s4-normalization-characterization.R test-prot-06-ruv.R test-prot-ruv-levels.R test-prot-ruv-replicate-structure.R test-prot-ruv-sweep-imputation.R test-prot-s4-norm-methods-direct-shared.R test-prot-s4-qc-direct-shared.R test-workflow-lipid-declarations.R test-workflow-metab-declarations.R test-workflow-native-session.R test-workflow-omics-ruv-cancor-once.R test-workflow-omics-ruv-complete-controls.R test-workflow-ruv-pearson-diagnostic.R test-workflow-saved-qc.R |
| `R/func_metab_da_results_long_format.R` | test-metab-02-da-results-long-format-characterization.R test-metab-02h-da-orchestration-characterization.R test-metab-da-orchestration-direct-shared.R test-metab-import-metrics-organisms-direct-shared.R |
| `R/func_metab_da_volcano_glimma.R` | test-glimma-plot.R test-metab-02ab-da-module-characterization.R test-metab-02e-da-volcano-glimma-characterization.R test-omics-audit-carry-overs.R test-workflow-omics-application-reads.R |
| `R/func_metab_norm.R` | test-metab-03d-norm-module-characterization.R test-metab-04-ruv-optimization.R test-metab-norm-standalone-shared.R test-omic-norm-module.R test-omics-norm-core.R test-omics-norm-itsd-selection.R test-workflow-omics-figures.R test-workflow-omics-ruv-cancor-once.R test-workflow-omics-ruv-complete-controls.R |
| `R/func_metab_s4_da_methods.R` | test-lipid-08-s4-da-analysis-methods.R test-metab-02p-da-analysis-characterization.R test-metab-02q-da-analysis-helper-characterization.R test-metab-s4-da-methods-direct-shared.R test-workflow-saved-da.R |
| `R/func_metab_s4_design_control_methods.R` | test-lipid-09-s4-normalization-methods.R test-lipid-12b-norm-ruv-helpers-shared.R test-lipid-12c-ruv-s4-characterization.R test-metab-02w-normalise-between-samples-characterization.R test-metab-02x-clean-design-matrix-characterization.R test-metab-02y-get-neg-ctrl-metab-anova-characterization.R test-metab-03a-ruviii-c-varying-characterization.R test-metab-04-ruv-optimization.R test-metab-norm-standalone-shared.R test-metab-s4-norm-core-shared.R test-prot-s4-norm-methods-direct-shared.R test-prot-s4-qc-direct-shared.R |
| `R/func_metab_s4_qc_methods.R` | test-lipid-05-s4-correlation-helpers.R test-lipid-12a-norm-module-characterization.R test-metab-01j-s4-pearson-correlation-characterization.R test-metab-01k-s4-plot-pearson-characterization.R test-metab-01n-s4-find-duplicate-feature-ids-characterization.R test-metab-01p-qc-duplicates-module-characterization.R test-metab-01q-qc-intensity-module-characterization.R test-metab-01u-qc-duplicates-module-characterization.R test-metab-02t-metabolite-intensity-filtering-characterization.R test-metab-02v-sample-correlation-filter-characterization.R test-metab-03d-norm-module-characterization.R test-metab-norm-standalone-shared.R test-metab-qc-shared.R test-metab-s4-qc-shared.R test-module-ci-lipid-norm.R test-module-ci-lipid-qc.R test-module-ci-metab-norm.R test-module-ci-metab-qc.R test-omics-audit-carry-overs.R test-omics-qc-core.R test-prot-04-design.R test-prot-04d-peptide-s4-normalization-characterization.R test-prot-05d-norm-module-characterization.R test-prot-s4-qc-direct-shared.R test-workflow-metab-declarations.R test-workflow-native-session.R test-workflow-omics-application-reads.R test-workflow-omics-figures.R test-workflow-omics-pearson-matrix.R test-workflow-ruv-pearson-diagnostic.R |
| `R/func_metab_s4_ruv_methods.R` | test-general-filemgmt-config-contracts.R test-general-filemgmt-config-direct-shared.R test-lipid-09-s4-normalization-methods.R test-lipid-12a-norm-module-characterization.R test-lipid-12b-norm-ruv-helpers-shared.R test-lipid-12c-ruv-s4-characterization.R test-metab-02z-ruv-cancor-characterization.R test-metab-03a-ruviii-c-varying-characterization.R test-metab-04-ruv-optimization.R test-metab-norm-standalone-shared.R test-module-ci-metab-norm.R test-omic-norm-module.R test-omics-norm-core.R test-prot-04d-peptide-s4-normalization-characterization.R test-prot-06-ruv.R test-prot-ruv-levels.R test-prot-ruv-replicate-structure.R test-prot-ruv-sweep-imputation.R test-prot-s4-norm-methods-direct-shared.R test-prot-s4-qc-direct-shared.R test-workflow-lipid-declarations.R test-workflow-metab-declarations.R test-workflow-native-session.R test-workflow-omics-ruv-cancor-once.R test-workflow-omics-ruv-complete-controls.R test-workflow-ruv-pearson-diagnostic.R test-workflow-saved-qc.R |
| `R/func_multiomics_enrich_kegg_helpers.R` | test-multiomics-enrich-shared.R |
| `R/func_multiomics_enrich_metabolomics_helpers.R` | test-multiomics-enrich-shared.R |
| `R/func_multiomics_enrich_reactome_helpers.R` | test-multiomics-enrich-shared.R |
| `R/func_multiomics_enrich_stringdb_reporting_helpers.R` | test-multiomics-enrich-shared.R |
| `R/func_pept_s4_qc_methods.R` | test-general-filemgmt-config-contracts.R test-lipid-05-s4-correlation-helpers.R test-lipid-10-s4-plotting-methods.R test-lipid-12a-norm-module-characterization.R test-metab-01g-s4-plot-pca-characterization.R test-metab-01h-s4-plot-rle-characterization.R test-metab-01i-s4-plot-density-characterization.R test-metab-01j-s4-pearson-correlation-characterization.R test-metab-01k-s4-plot-pearson-characterization.R test-metab-03d-norm-module-characterization.R test-metab-norm-standalone-shared.R test-metab-s4-plotting-direct-shared.R test-metab-s4-qc-shared.R test-module-ci-lipid-norm.R test-module-ci-lipid-qc.R test-module-ci-metab-norm.R test-module-ci-metab-qc.R test-module-ci-prot-norm.R test-module-ci-prot-peptide-qc.R test-module-ci-prot-protein-qc.R test-omic-norm-module.R test-omics-norm-core.R test-prot-04c-peptide-plot-density-ggplot-shared.R test-prot-04-design.R test-prot-04d-peptide-s4-normalization-characterization.R test-prot-05d-norm-module-characterization.R test-prot-07-da-analysis.R test-prot-07d-da-model-stats-shared.R test-prot-da-compat-shared.R test-prot-s4-qc-direct-shared.R test-workflow-gui-memory-repairs.R test-workflow-native-session.R test-workflow-omics-figures.R test-workflow-omics-pearson-matrix.R test-workflow-ruv-pearson-diagnostic.R test-workflow-saved-qc.R |
| `R/func_prot_annotation_accession_ranking_helpers.R` | test-module-ci-prot-protein-qc.R test-prot-04d-peptide-s4-normalization-characterization.R test-prot-05d-norm-module-characterization.R test-prot-annotation-accession-ranking-direct-shared.R test-prot-annotation-processing-cleaning-direct-shared.R |
| `R/func_prot_annotation_uniprot_full_helpers.R` | test-omic-gui-config.R test-prot-04a-design-module-characterization.R test-prot-annotation-ensembl-full-direct-shared.R test-prot-design-core.R |
| `R/func_prot_annotation_uniprot_helpers.R` | test-e2e-browser-harness.R test-prot-annotation-ensembl-full-direct-shared.R test-prot-annotation-uniprot-direct-shared.R test-prot-design-core.R |
| `R/func_prot_da_results_io.R` | test-prot-07d-da-model-stats-shared.R |
| `R/func_prot_da_results_long_format.R` | test-prot-07c-da-results-characterization.R test-prot-07c-da-results-output-shared.R test-prot-07d-da-model-stats-shared.R test-prot-da-compat-shared.R |
| `R/func_prot_da_results_methods.R` | test-module-ci-prot-da.R test-prot-07b-da-handlers-compat.R test-prot-07c-da-results-output-shared.R |
| `R/func_prot_da_results_summary.R` | test-lipid-07-s4-da-plot-methods.R test-metab-02k-da-numsig-barplot-characterization.R test-module-ci-prot-da.R test-prot-07c-da-results-characterization.R test-prot-07c-da-results-output-shared.R test-prot-07-da-analysis.R test-prot-07d-da-model-stats-shared.R test-prot-da-compat-shared.R |
| `R/func_prot_import_tmt.R` | test-module-ci-prot-import.R test-prot-import-core.R test-prot-import-quantity-text.R test-prot-import-tmt-direct-shared.R test-workflow-native-session.R |
| `R/func_prot_limpa_qc_helpers.R` | test-prot-05f-limpa-qc-helpers-shared.R |
| `R/func_prot_qc_correlation_helpers.R` | test-metab-01j-s4-pearson-correlation-characterization.R test-module-ci-prot-norm.R test-prot-qc-correlation-direct-shared.R |
| `R/func_prot_qc_filtering_helpers.R` | test-prot-02-qc-filtering.R test-prot-s4-qc-direct-shared.R test-workflow-gui-memory-repairs.R |
| `R/func_prot_qc_peptide_group_filters.R` | test-prot-02-qc-peptide-groupaware.R test-prot-qc-peptide-intensity-group-batches.R test-prot-qc-peptide-methods-direct-shared.R test-workflow-partitions.R |
| `R/func_prot_qc_peptide_methods.R` | test-general-filemgmt-config-contracts.R test-general-filemgmt-config-roundtrip-contracts.R test-general-filemgmt-duplicate-resolution-contracts.R test-perf-003.R test-perf-105.R test-perf-106.R test-perf107-consumers.R test-perf107-synthetic.R test-prot-02b-qc-peptide-module-contracts.R test-prot-02c-qc-peptide-intensity-shared.R test-prot-02g-qc-peptide-qvalue-shared.R test-prot-02h-peptide-audit-shared.R test-prot-02i-qc-peptide-sample-shared.R test-prot-02j-qc-peptide-protein-shared.R test-prot-02m-qc-peptide-replicate-shared.R test-prot-02r-qc-reporting-helpers-shared.R test-prot-02s-identification-evidence.R test-prot-02t-qc-peptide-exclusions.R test-prot-qc-core.R test-prot-qc-peptide-methods-direct-shared.R test-prot-qc-peptide-replicate-filters-direct-shared.R test-prot-s4-qc-direct-shared.R test-workflow-analysis-record.R test-workflow-audit-memory-member.R test-workflow-checkpoint-dia-state.R test-workflow-checkpoint-harness.R test-workflow-checkpoints.R test-workflow-execution-adapters.R test-workflow-execution-fixtures.R test-workflow-gui-memory-lifetimes.R test-workflow-kit-harness.R test-workflow-native-session.R test-workflow-partitions.R test-workflow-prot-declarations.R test-workflow-saved-exports.R |
| `R/func_prot_s4_accession_methods.R` | test-prot-02e-qc-protein-cleanup-shared.R test-prot-04-design.R test-prot-04d-peptide-s4-normalization-characterization.R test-prot-05d-norm-module-characterization.R test-workflow-analysis-record.R test-workflow-lipid-declarations.R test-workflow-metab-declarations.R test-workflow-operation-registry.R |
| `R/func_prot_s4_da_results.R` | none found |
| `R/func_prot_s4_norm_controls.R` | test-prot-norm-optimization-direct-shared.R test-prot-ruv-sweep-imputation.R |
| `R/func_prot_s4_norm_methods.R` | test-general-filemgmt-config-contracts.R test-general-filemgmt-config-direct-shared.R test-lipid-09-s4-normalization-methods.R test-lipid-12a-norm-module-characterization.R test-lipid-12b-norm-ruv-helpers-shared.R test-lipid-12c-ruv-s4-characterization.R test-metab-02w-normalise-between-samples-characterization.R test-metab-02z-ruv-cancor-characterization.R test-metab-03a-ruviii-c-varying-characterization.R test-metab-04-ruv-optimization.R test-metab-norm-standalone-shared.R test-metab-s4-norm-core-shared.R test-module-ci-lipid-norm.R test-module-ci-metab-norm.R test-omic-norm-module.R test-omics-norm-core.R test-prot-04-design.R test-prot-04d-peptide-s4-normalization-characterization.R test-prot-05b-norm-module-contracts.R test-prot-05d-norm-module-characterization.R test-prot-06-ruv.R test-prot-norm-optimization-direct-shared.R test-prot-ruv-levels.R test-prot-ruv-replicate-structure.R test-prot-ruv-sweep-imputation.R test-prot-s4-norm-methods-direct-shared.R test-prot-s4-qc-direct-shared.R test-workflow-lipid-declarations.R test-workflow-metab-declarations.R test-workflow-native-session.R test-workflow-omics-ruv-cancor-once.R test-workflow-omics-ruv-complete-controls.R test-workflow-prot-declarations.R test-workflow-ruv-pearson-diagnostic.R test-workflow-saved-consumers.R test-workflow-saved-exports.R test-workflow-saved-installed.R test-workflow-saved-qc.R |
| `R/func_prot_s4_qc_correlation_methods.R` | test-lipid-05-s4-correlation-helpers.R test-lipid-12a-norm-module-characterization.R test-metab-01j-s4-pearson-correlation-characterization.R test-metab-01k-s4-plot-pearson-characterization.R test-metab-03d-norm-module-characterization.R test-metab-norm-standalone-shared.R test-metab-s4-qc-shared.R test-module-ci-lipid-norm.R test-module-ci-lipid-qc.R test-module-ci-metab-norm.R test-module-ci-metab-qc.R test-prot-04-design.R test-prot-04d-peptide-s4-normalization-characterization.R test-prot-05d-norm-module-characterization.R test-prot-s4-qc-direct-shared.R test-workflow-native-session.R test-workflow-omics-figures.R test-workflow-omics-pearson-matrix.R test-workflow-ruv-pearson-diagnostic.R |
| `R/func_prot_s4_qc_filter_methods.R` | test-general-filemgmt-config-contracts.R test-general-filemgmt-duplicate-resolution-contracts.R test-general-filemgmt-update-missing-value-shared.R test-lipid-09-s4-normalization-methods.R test-lipid-12c-ruv-s4-characterization.R test-metab-02w-normalise-between-samples-characterization.R test-metab-02x-clean-design-matrix-characterization.R test-metab-03a-ruviii-c-varying-characterization.R test-metab-s4-norm-core-shared.R test-module-ci-prot-protein-qc.R test-prot-02b-qc-peptide-module-contracts.R test-prot-02c-qc-protein-module-contracts.R test-prot-02d-qc-protein-intensity-shared.R test-prot-02h-qc-protein-replicate-shared.R test-prot-02j-qc-peptide-protein-shared.R test-prot-02s-identification-evidence.R test-prot-05d-norm-module-characterization.R test-prot-qc-peptide-methods-direct-shared.R test-prot-s4-norm-methods-direct-shared.R test-prot-s4-qc-direct-shared.R test-workflow-analysis-record.R test-workflow-checkpoint-dia-state.R test-workflow-checkpoint-harness.R test-workflow-dia-products.R test-workflow-gui-diagnostic.R test-workflow-gui-memory-lifetimes.R test-workflow-gui-memory-repairs.R test-workflow-memory-mechanism-member.R test-workflow-native-session.R test-workflow-partitions.R test-workflow-prot-declarations.R |
| `R/mod_omic_enrich.R` | test-module-ci-prot-da.R test-module-ci-prot-enrich.R test-prot-11a-enrichment-module-characterization.R |
| `R/mod_prot_da_handlers_init.R` | test-omic-gui-config.R test-prot-07b-da-handlers-compat.R |
| `R/mod_prot_design_state_helpers.R` | test-prot-04a-design-module-characterization.R test-prot-design-core.R |
| `R/mod_prot_import_server.R` | test-module-ci-prot-import.R test-prot-01c-import-module-contracts.R test-prot-01d-import-module-characterization.R test-workflow-server-shared.R |
| `R/mod_prot_norm_correlation_helpers.R` | test-dia-baseline-compatibility.R test-module-ci-prot-norm.R test-prot-05b-norm-module-contracts.R test-prot-05d-norm-module-characterization.R test-workflow-native-session.R |
| `R/mod_prot_norm_render_helpers.R` | test-omic-norm-module.R test-prot-05b-norm-module-contracts.R test-workflow-dia-products.R |
| `R/mod_prot_norm_ruv_helpers.R` | test-dia-baseline-compatibility.R test-module-ci-prot-dia-limpa.R test-module-ci-prot-norm.R test-prot-05b-norm-module-contracts.R test-prot-05d-norm-module-characterization.R test-prot-norm-core.R |
