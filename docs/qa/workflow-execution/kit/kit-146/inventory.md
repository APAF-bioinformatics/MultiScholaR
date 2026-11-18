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
