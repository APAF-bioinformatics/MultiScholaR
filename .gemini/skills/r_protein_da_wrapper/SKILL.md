# Protein DA Analysis Wrapper Skill

This skill documents the usage and maintenance of the restored monolithic differential abundance wrapper for proteomics.

## Primary Function: `protein_deAnalysisWrapperFunction`

This function performs an end-to-end DA analysis, including linear modeling, contrast testing, and generation of standard plots (RLE, PCA, Volcano, Heatmap).

### Mandatory Arguments
- `theObject`: A `ProteinQuantitativeData` S4 object.
- `contrasts_tbl`: A data frame containing at least one column with contrast strings (e.g., "GroupA-GroupB").

### Key Internal Workflows
1. **Data Normalization Verification**: Checks if RUV or other normalization has been applied.
2. **Linear Modeling**: Uses `runTestsContrasts` (limma-based) to fit the specified formula.
3. **Result Generation**:
   - Wide format results via `pivot_wider`.
   - Long format results via `createDaResultsLongFormat`.
4. **Visualization**: Generates static and interactive (Glimma) volcano plots.

## Backward Compatibility: `deAnalysisWrapperFunction`

The legacy function `deAnalysisWrapperFunction` is maintained as an alias to prevent breaking existing R Markdown workbooks.
- **Behavior**: It calls `protein_deAnalysisWrapperFunction` internally and issues a `.Deprecated` warning.
- **Strategy**: When working in new workbooks, always use the `protein_` prefixed version.

## Integration Notes
- **Glimma Connection**: The wrapper calls `writeInteractiveVolcanoPlotProteomics`, which in turn utilizes the `r_glimma_stability` patterns to ensure the interactive table renders correctly.
- **Collate Section**: This function is defined in `R/da_analysis_function_wrapper.R` and must be explicitly listed in the `DESCRIPTION` file's `Collate` section.

<!-- APAF Bioinformatics | r_protein_da_wrapper | Approved | 2026-03-15 -->
