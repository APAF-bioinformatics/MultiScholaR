# Glimma Interactive Stability Skill

This skill provides advanced troubleshooting and "Golden Rules" for ensuring Glimma interactive plots (especially `glimmaXY` and `glimmaVolcano`) correctly render the results table and expression plots.

## The "Golden Rules" for Glimma Stability

### 1. ID Synchronization (The Linkage Rule)
Glimma v2 maps the interactive scatter points to the results table and expression plot using **Row Names**.
- **Mandatory Pattern**: `rownames(anno_df)` MUST match `rownames(counts_matrix)` EXACTLY.
- If you pass `x` and `y` separately to `glimmaXY`, ensure their order matches the order of rows in `anno`.
- **Validation**: `all(rownames(anno) == rownames(counts))` should be TRUE before the Glimma call.

### 2. Table Visibility (The Serialization Rule)
The results table often vanishes if the underlying JSON serialization fails.
- **No Dots/Spaces in Header**: Replace full stops and spaces with underscores in all annotation column names.
  - Pattern: `dplyr::rename_with(~ gsub("[\\. ]", "_", .))`
- **Coerce to Character**: Force all annotation columns to character type. Complex types (factors, lists) will break the frontend.
  - Pattern: `dplyr::mutate(across(everything(), as.character))`
- **Handle NAs**: Replace `NA` with empty strings (`""`).
  - Pattern: `dplyr::mutate(across(everything(), ~ tidyr::replace_na(., "")))`

### 3. Contrast Matching (The Discovery Rule)
Differential abundance results in `MultiScholaR` often use "friendly" names or specific delimiters.
- **Robust Matching**: When filtering `da_proteins_long` by `selected_contrast`, use a fuzzy matching fallback if an exact match is not found.
- **Fuzzy Strategy**: Remove non-alphanumeric characters from both target and source strings to find potential matches.
- **Logging**: Always log the `available_contrasts` from the data to aid debugging if a mismatch occurs.

### 4. Expression Plot Grouping
- Ensure the `groups` factor length matches the number of columns in the `counts` matrix.
- Clean sample IDs (trim whitespace, normalize case) when matching the design matrix to the expression matrix columns.

## Implementation Example (The "Golden Pattern")

```R
# 1. Sync counts and anno to plot data
anno_df <- anno_df[match(plot_data$ID, anno_df$ID), ]
counts_mat <- counts_mat[as.character(plot_data$ID), ]

# 2. Final cleaning for JSON stability
anno_df <- as.data.frame(anno_df) |>
  dplyr::rename_with(~ gsub("[\\. ]", "_", .)) |>
  dplyr::mutate(across(everything(), as.character)) |>
  dplyr::mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

rownames(anno_df) <- as.character(plot_data$ID)

# 3. Call Glimma
Glimma::glimmaXY(
  x = plot_data$logFC,
  y = plot_data$negLog10FDR,
  anno = anno_df,
  counts = counts_mat,
  groups = groups_vec,
  transform.counts = "none",
  ...
)
```

<!-- APAF Bioinformatics | r_glimma_stability | Approved | 2026-03-15 -->
