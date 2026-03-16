# Implementation Plan: Universal "Sanitize Sample Names" Feature

## 1. Problem Description
Sample names in experimental datasets often contain characters that are not "syntactically valid" in R. Common issues include:
- **Spaces:** `Sample 1`
- **Special Characters:** `Control-A!`, `Batch@2`
- **Leading Numbers:** `123_Test`

When these names are processed by R packages like `iq`, `limma`, or `vroom`, they are often silently transformed into valid names (e.g., `123_Test` becomes `X123_Test`). If this transformation happens to the quantitative data columns but NOT to the design matrix `Run` column, the S4 object's validity check fails:
`"Samples in protein data and design matrix must be the same"`

## 2. Proposed Solution
We will implement an optional, user-controlled "Sanitize Sample Names" feature at the **Data Import** stage for all three omics pipelines (Proteomics, Lipidomics, Metabolomics). This ensures that sample IDs are "R-safe" across the entire downstream analysis pipeline.

## 3. Implementation Details

### Phase A: UI Updates (Import Modules)
Add a "Sanitize Sample Names" checkbox and help text to the following files:
- `R/mod_prot_data_import.R`
- `R/mod_lipid_data_import.R`
- `R/mod_metab_data_import.R`

**UI Code Example:**
```r
shiny::checkboxInput(ns("sanitize_names"), "Sanitize Sample Names", value = TRUE)
shiny::helpText("Clean sample IDs (e.g., '123-Sample!' -> 'x123_sample') for better compatibility with downstream analysis.")
```

### Phase B: Logic Updates (Server Side)
Modify the data processing logic in each module to apply `janitor::make_clean_names()` if the option is enabled:

1.  **Proteomics:**
    - Clean `design_matrix$Run`.
    - Clean the `Run` column in the long-format peptide data.
2.  **Lipidomics / Metabolomics:**
    - Clean `design_matrix$Run`.
    - Clean the column names of the wide-format quantitative data (excluding the primary ID column).

### Phase C: Logging and Notifications
- **Traceability:** Log the transformation (Original -> Cleaned) using `logger::log_info`.
- **Transparency:** Show a `shiny::showNotification` to inform the user that their sample names were adjusted for safety.

## 4. Targeted Modules

| Pipeline | UI File | S4 Constructor Point |
| :--- | :--- | :--- |
| **Proteomics** | `R/mod_prot_data_import.R` | Before `PeptideQuantitativeData()` |
| **Lipidomics** | `R/mod_lipid_data_import.R` | Before `createLipidomicsAssayData()` |
| **Metabolomics** | `R/mod_metab_data_import.R` | Before `createMetaboliteAssayData()` |

## 5. Verification Plan

### Automated Tests
1.  **Mock Data Test:** Create a design matrix with names `c("1st Sample", "Sample-A!", "Run 3")`.
2.  **Transformation Verification:** Assert that `janitor::make_clean_names` correctly yields `c("x1st_sample", "sample_a", "run_3")`.
3.  **S4 Integrity:** Verify the S4 object constructor accepts the sanitized data and passes its internal validity check.

### Manual Verification
1.  Launch the MultiScholaR Shiny App.
2.  Import a dataset with "bad" sample names (e.g., the Neurolincs parquet dataset).
3.  Enable "Sanitize Sample Names".
4.  Confirm the data loads successfully and proceed to Normalization to ensure no downstream mismatches occur.

<!-- APAF Bioinformatics | plan_sanitize_sample_names_20260316.md | Approved | 2026-03-16 -->
