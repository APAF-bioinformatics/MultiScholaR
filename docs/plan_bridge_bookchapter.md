# Plan: Bridge BookChapter Branch and Fix Missing Functions

This plan outlines the steps to resolve the missing `formatDIANN` function and other branch-related discrepancies to ensure the `DIA_workflow_limpa_starter.rmd` works correctly with the local `MultiScholaR` installation.

## Background

The user encountered an error: `Error in formatDIANN(data_tbl_parquet_filt) : could not find function "formatDIANN"`. This happened because the `bugfix/volcano-plot-glimma` branch lacks the export for `formatDIANN`, even though a version of it was recently added to `R/func_prot_import.R`.

## Identified Issues

1. **`formatDIANN` Export Missing**: The function is defined in `R/func_prot_import.R` but lacks Roxygen2 headers and the `@export` tag.
2. **Parquet Support**: The `formatDIANNParquet` function exists in the current branch but not in `BookChapter`. The Rmd currently uses `formatDIANN`.
3. **Branch Divergence**: Several functions and structures (like `R/file_management.R`) exist in `BookChapter` but have been reorganized or renamed in the current branch.

## Proposed Changes

### [MultiScholaR Package]

#### [MODIFY] [R/func_prot_import.R](file:///Users/ignatiuspang/Workings/2025/MultiScholaR-volcano-fix/R/func_prot_import.R)
- Add Roxygen2 headers and `@export` to `formatDIANN`.
- Ensure `formatDIANN` logic matches the requirement for processing parquet-originated `EList` objects and is compatible with the `EList` output from `limpa::readDIANN`.
- Verify `formatDIANNParquet` is correctly implemented and exported.
- Add necessary `@importFrom` tags for `dplyr`, `tidyr`, and `tibble` to avoid conflicts.

#### [MODIFY] [NAMESPACE](file:///Users/ignatiuspang/Workings/2025/MultiScholaR-volcano-fix/NAMESPACE)
- Update by running `devtools::document()`.

## Verification Plan

### Automated Tests
1. **Regenerate Documentation**:
   Run `devtools::document(pkg = "/Users/ignatiuspang/Workings/2025/MultiScholaR-volcano-fix")`.
2. **Verify Exports**:
   Check that `export(formatDIANN)` and `export(formatDIANNParquet)` appear in the `NAMESPACE` file.
3. **Run Existing Tests**:
   Run `devtools::test(filter = "prot-01-import")` to ensure `importDIANNData` still works.
4. **New Unit Test**:
   Create a temporary test script to verify `formatDIANN` correctly transforms a mock `EList` into a DIA-NN-like data frame.

### Manual Verification
1. **Rmd Validation**:
   Run the following chunk in [DIA_workflow_limpa_starter.rmd](file:///Users/ignatiuspang/Workings/2026/neurolincs_book_chapter_analysis/scripts/proteomics/DIA_workflow_limpa_starter.rmd) and verify it no longer errors:
   ```R
   data_tbl_parquet_filt <- limpa::readDIANN(...)
   data_tbl <- formatDIANN(data_tbl_parquet_filt)
   ```
2. **Project Setup**:
   Ensure `setupDirectories` and parquet reading work together by verifying `data_tbl` is correctly populated.
3. **Console Verification**:
   Confirm that `formatDIANN` can be called from the console after `devtools::load_all()` or re-installing the package.

<!-- APAF Bioinformatics | plan_bridge_bookchapter.md | Approved | 2026-03-14 -->
