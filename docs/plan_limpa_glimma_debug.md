# Orchestration Plan: Limpa Script, Glimma Plot, and Test Data Capture Fixes

This document outlines the precise plan to fix three issues reported in the current branch (`bugfix/volcano-plot-glimma`), intended to be executed via `\orchestrate debug`.

## Issue 1: Missing `deAnalysisWrapperFunction` in Limpa Script
**Context:** The `limpa` script fails because `deAnalysisWrapperFunction` is missing in the current branch. It previously existed in the `BookChapter` branch.
**Plan:**
1. **Locate Component: (Completed)** Extract the `deAnalysisWrapperFunction` and related interactive volcano plotting wrappers (`writeInteractiveVolcanoPlotProteomics`, `writeInteractiveVolcanoPlotProteomicsWidget`) from `R/de_analysis_function_wrapper.R` in the `BookChapter` branch.
2. **Implement in Current Branch:(Completed)** 
   - Create or update the relevant script (e.g., `R/da_analysis_function_wrapper.R` or `R/func_prot_da.R`) to include `deAnalysisWrapperFunction`.
   - Ensure the function is properly exported via Roxygen tags (`#' @export`).
   - Run `devtools::document()` to update the `NAMESPACE`.
3. **Verify Connection:** Ensure that the input `theObject` and `contrasts_tbl` conform to what `deAnalysisWrapperFunction` expects and that the limpa DPC analysis executes without "function not found" errors.

## Issue 1.5: Missing `protein_deAnalysisWrapperFunction` Parameter Errors in Wrapper (Completed)
**Context:** It appears the function `protein_deAnalysisWrapperFunction` in the `R/da_analysis_function_wrapper.R` script, an error is logged stating: `unused argument (list_of_de_tables = list(contrasts_results_table))` when calling `getSignificantData()`. This happens because `getSignificantData` was updated to accept `list_of_da_tables` instead.
**Plan:**
1. **Fix Parameter Names:** Update instances of `list_of_de_tables` to `list_of_da_tables` in `protein_deAnalysisWrapperFunction` for `getSignificantData` and potentially `printCountDeGenesTable` where applicable.
2. **Re-Document:** Run `devtools::document()` to make the exported function available to the RMD script.

## Issue 2: Glimma Interactive Volcano Plot Missing Table/Expression Plot (Completed)
**Context:** The R/Shiny Glimma volcano plot is missing the interactive table and the expression plot. This is a known issue likely related to column naming, missing `counts` data, or incorrect parameter mapping in `glimmaXY`.
**Plan:**
1. **Check Data Types and Column Names:** In `generateProtDAVolcanoPlotGlimma` (`R/func_prot_da.R`):
   - Ensure `display_df` (the `anno` parameter) has NO full stops (`.`) in column names. They should be strictly formatted (e.g., replacing spaces and dots with underscores), as Glimma drops the table if invalid column names exist.
   - Verify that `display.columns` passed to `glimmaXY` exactly matches the names in `display_df`.
2. **Fix Counts and Groups Mapping:** 
   - Investigate the block where `counts_mat` and `groups_vec` are generated (around line 371 in `R/func_prot_da.R`). If case-sensitive or whitespace issues exist between `colnames(counts_mat)` and `this_design[[sample_id_col]]`, `counts_mat` becomes NULL, which suppresses the expression plot.
   - Refactor the sample matching to be completely robust, ensuring `counts_mat` correctly preserves rownames (`Protein.Ids`) and matches `groups_vec`.
3. **Validate Widget UI Output:** Ensure the Shiny UI uses the correct widget rendering (e.g., `htmlwidgets::shinyRenderWidget` or `renderGlimma`) instead of a generic HTML output, if required for Glimma 2.0+.

## Issue 3: Capturing Test Data to Analysis Folder (Completed)
**Context:** The `.capture_checkpoint` currently saves test data to the `tests/testdata/...` folder, but it needs to save directly to the analysis output folder.
**Plan:**
1. **Modify `.capture_checkpoint` in `R/utils_shiny_logging.R`:**
   - Update the `base_dir` resolution logic. Investigate global options like `multischolar.analysis_dir` or extract the path dynamically to point to the current analysis output folder rather than `tests/testdata/`.
   - Ensure that if the targeted analysis directory doesn't exist, it correctly creates the required subdirectories.
2. **Update Integration:** Ensure the `differentialAbundanceAnalysisHelper` and other modules correctly configure the output path via `options()` before capturing test data.

## Execution
Run the following YOLO command in terminal to commence fixing:
`/orchestrate debug @docs/plan_limpa_glimma_debug.md`
