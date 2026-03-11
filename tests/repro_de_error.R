
library(dplyr)
library(tidyr)
library(tibble)
library(rlang)

# Mock data
lfc_qval_tbl <- tibble(
  metabolite_id = c("M1", "M2"),
  logFC = c(1.2, -0.5),
  P.Value = c(0.01, 0.04),
  fdr_qvalue = c(0.05, 0.1),
  fdr_value_bh = c(0.05, 0.1),
  raw_pvalue = c(0.01, 0.04),
  total_proteins_quantified = 2, # Just to have some columns
  assay = "LCMS_Pos",
  comparison = "groupTreatment-groupControl",
  friendly_name = "Treatment_vs_Control",
  significant = "Up"
)

expr_matrix <- matrix(
  c(10, 11, 12, 13, 11, 10),
  nrow = 2,
  dimnames = list(c("M1", "M2"), c("S1", "S2", "S3"))
)

design_matrix <- data.frame(
  Run = c("S1", "S2", "S3"),
  group = c("Control", "Treatment", "Treatment"),
  stringsAsFactors = FALSE
)

# Test function
test_repro <- function() {
  sample_id_col <- "Run"
  group_id_col <- "group"
  metabolite_id_col <- "metabolite_id"
  
  message("Testing createMetabDeResultsLongFormat logic...")
  
  # 1. Pivot longer
  intensity_long <- expr_matrix |>
    as.data.frame() |>
    tibble::rownames_to_column(metabolite_id_col) |>
    tidyr::pivot_longer(
      cols = -!!rlang::sym(metabolite_id_col),
      names_to = sample_id_col,
      values_to = "intensity"
    ) |>
    dplyr::left_join(design_matrix, by = sample_id_col)
  
  # 2. Pivot wider
  intensity_wide <- intensity_long |>
    dplyr::mutate(
      col_name = paste0("intensity.", !!rlang::sym(sample_id_col), ".", !!rlang::sym(group_id_col))
    ) |>
    dplyr::select(!!rlang::sym(metabolite_id_col), col_name, intensity) |>
    tidyr::pivot_wider(
      id_cols = !!rlang::sym(metabolite_id_col),
      names_from = col_name,
      values_from = intensity
    )
  
  # 3. Separate wider
  message("Testing separate_wider_delim...")
  # Use a local copy
  current_tbl <- lfc_qval_tbl
  
  tryCatch({
    current_tbl <- current_tbl |>
      tidyr::separate_wider_delim(
        comparison,
        delim = "-",
        names = c("numerator", "denominator"),
        cols_remove = FALSE
      )
  }, error = function(e) {
    message("Caught expected error in separate_wider_delim: ", e$message)
  })
  
  message("Success!")
}

# Run cases
message("\n--- Case 1: Valid ---")
test_repro()

message("\n--- Case 2: Invalid comparison (no delimiter) ---")
lfc_qval_tbl$comparison <- "groupTreatment"
tryCatch({
  lfc_qval_tbl |>
    tidyr::separate_wider_delim(
      comparison,
      delim = "-",
      names = c("numerator", "denominator"),
      cols_remove = FALSE
    )
}, error = function(e) {
  message("Error message: ", e$message)
})

message("\n--- Case 3: Missing comparison column ---")
lfc_qval_tbl_missing <- lfc_qval_tbl |> select(-comparison)
tryCatch({
  lfc_qval_tbl_missing |>
    tidyr::separate_wider_delim(
      comparison,
      delim = "-",
      names = c("numerator", "denominator"),
      cols_remove = FALSE
    )
}, error = function(e) {
  message("Error message: ", e$message)
})
