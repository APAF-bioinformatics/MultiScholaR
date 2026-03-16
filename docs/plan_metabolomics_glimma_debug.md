# Debug Plan: Metabolomics & Lipidomics Glimma Plot Fix

**Objective:** Fix the missing interactive table and expression plots in the Metabolomics (and simultaneously Lipidomics) Glimma volcano plots. 

## Root Cause Analysis
Based on the successful fixes applied to the Proteomics Glimma plot earlier today, the Metabolomics and Lipidomics plots are suffering from the exact same underlying issues:

1. **Row Name Synchronization:** `Glimma::glimmaVolcano` (v2) requires the `rownames(anno)` to perfectly match `rownames(counts)`. Currently, the `anno` dataframe is created inline during the function call without explicitly setting its row names to the metabolite/lipid IDs.
2. **JSON Serialization Failures:** If the annotation data contains `NA` values or is treated as non-character (e.g., factors/logicals), it can silently break the JSON serialization when constructing the HTML widget, resulting in the interactive table and expression plot disappearing.

## Affected Files & Functions
1. **File:** `R/func_metab_da.R`
   - **Function:** `generateMetabDAVolcanoPlotGlimma`
2. **File:** `R/func_lipid_da.R` (Bonus: Proactive fix to prevent the same issue for lipids)
   - **Function:** `generateLipidDAVolcanoPlotGlimma`

## Step-by-Step Implementation Plan

### 1. Fix Metabolomics (`R/func_metab_da.R`)

Locate the `generateMetabDAVolcanoPlotGlimma` function. Instead of creating the `anno` dataframe inline inside the `glimmaVolcano()` call, build and sanitize it directly **before** the call.

**Change from this:**
```R
glimma_widget <- Glimma::glimmaVolcano(
    x = fit_obj,
    coef = coef_index,
    counts = counts_mat,
    groups = groups,
    anno = data.frame(
        ID = volcano_tab$metabolite_id,
        Name = volcano_tab$display_name,
        Assay = volcano_tab$assay
    ),
    # ...
```

**To this:**
```R
# [FIX]: Prepare annotation dataframe perfectly matching row names and char types without NAs
anno_df <- data.frame(
    ID = volcano_tab$metabolite_id,
    Name = volcano_tab$display_name,
    Assay = volcano_tab$assay,
    stringsAsFactors = FALSE
)

# Clean NAs to empty strings and cast everything to character 
anno_df <- anno_df |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ tidyr::replace_na(., ""))) |>
    as.data.frame()

# CRITICAL: Rownames MUST match the fit object / counts matrix
rownames(anno_df) <- as.character(volcano_tab$metabolite_id)

glimma_widget <- Glimma::glimmaVolcano(
    x = fit_obj,
    coef = coef_index,
    counts = counts_mat,
    groups = groups,
    anno = anno_df,
    # ...
```

### 2. Fix Lipidomics (`R/func_lipid_da.R`)

Perform the exact same structural fix in `generateLipidDAVolcanoPlotGlimma`. 

**Make building `anno_df` explicit before the `glimmaVolcano` call:**
```R
# [FIX]: Prepare annotation dataframe perfectly matching row names and char types without NAs
anno_df <- data.frame(
    ID = volcano_tab$lipid_id,
    Name = volcano_tab$display_name,
    Assay = volcano_tab$assay,
    stringsAsFactors = FALSE
)

# Clean NAs to empty strings and cast everything to character 
anno_df <- anno_df |>
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ tidyr::replace_na(., ""))) |>
    as.data.frame()

# CRITICAL: Rownames MUST match the fit object / counts matrix
rownames(anno_df) <- as.character(volcano_tab$lipid_id)

glimma_widget <- Glimma::glimmaVolcano(
    x = fit_obj,
    coef = coef_index,
    counts = counts_mat,
    groups = groups,
    anno = anno_df,
    # ...
```

### 3. Verification

Once these files are updated, the pipeline should be rerun (e.g. `DIA_workflow_limpa_starter.rmd` or Shiny App). Checking the generated `.html` Glimma volcano plots will instantly confirm that the tables and expression point plots are visible alongside the main scatter plot.

<!-- APAF Bioinformatics | R_is_for_Robot | Approved -->
