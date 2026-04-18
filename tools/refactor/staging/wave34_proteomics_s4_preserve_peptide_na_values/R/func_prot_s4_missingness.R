#' @export
setMethod(
  f = "preservePeptideNaValues",
  signature = c("PeptideQuantitativeData", "ProteinQuantitativeData"),
  definition = function(peptide_obj, protein_obj) {
    preservePeptideNaValuesHelper(peptide_obj, protein_obj)
  }
)

#' @export
preservePeptideNaValuesHelper <- function(peptide_obj, protein_obj) {
  sample_id_column <- peptide_obj@sample_id
  protein_id_column <- peptide_obj@protein_id_column

  check_peptide_value <- peptide_obj@peptide_data |>
    group_by(!!sym(sample_id_column), !!sym(protein_id_column)) |>
    summarise(
      Peptide.Normalised = sum(Peptide.Normalised, na.rm = TRUE),
      is_na = sum(is.na(Peptide.Normalised)),
      num_values = n()
    ) |>
    mutate(Peptide.Normalised = if_else(is_na == num_values, NA_real_, Peptide.Normalised)) |>
    ungroup() |>
    arrange(!!sym(sample_id_column)) |>
    pivot_wider(
      id_cols = !!sym(protein_id_column),
      names_from = !!sym(sample_id_column),
      values_from = Peptide.Normalised,
      values_fill = NA_real_
    )

  check_peptide_value_cln <- check_peptide_value[
    rownames(protein_obj@protein_quant_table),
    colnames(protein_obj@protein_quant_table)
  ]

  if (length(which(rownames(protein_obj@protein_quant_table) == rownames(check_peptide_value_cln))) != nrow(check_peptide_value_cln)) {
    stop("The rows in the protein object and the peptide object do not match")
  }

  if (length(which(colnames(protein_obj@protein_quant_table) == colnames(check_peptide_value_cln))) != ncol(check_peptide_value_cln)) {
    stop("The columns in the protein object and the peptide object do not match")
  }

  protein_obj@protein_quant_table[is.na(check_peptide_value_cln)] <- NA

  protein_obj
}

