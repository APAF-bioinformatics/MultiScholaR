# ----------------------------------------------------------------------------
# prepareDataForVolcanoPlot
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Prepare data for volcano plot
#' @param input_table The input table with the log fold-change and q-value columns.
#' @param protein_id_column The name of the column representing the protein ID (tidyverse style).
#' @param uniprot_table The uniprot table with the imporatnt info on each protein
#' @param uniprot_protein_id_column The name of the column representing the protein ID in the uniprot table (tidyverse style).
#' @param number_of_genes The number of genes to show in the volcano plot.
#' @param fdr_threshold The FDR threshold for the volcano plot.
#' @param fdr_column The name of the column representing the FDR value (tidyverse style).
#' @param log2FC_column The name of the column representing the log fold-change (tidyverse style).
#' @return A table with the following columns:
#' label  The label of the significant proteins.
#' log2FC The log2 fold-change of the significant proteins.
#' lqm  The -log10 of the q-value.
#' colour The colour of the significant proteins.
#' rank_positive: The rank of the positive fold-change values.
#' rank_negative: The rank of the negative fold-change values.
#' gene_name_significant  The gene name of the significant proteins.
#'
#' @export
prepareDataForVolcanoPlot <- function(
  input_table,
  protein_id_column = uniprot_acc,
  uniprot_table,
  uniprot_protein_id_column = uniprot_acc_first,
  gene_name_column = gene_name,
  number_of_genes = 3000,
  fdr_threshold = 0.05,
  fdr_column = q.mod,
  log2FC_column = log2FC
) {
  temp_col_name <- as_string(as_name(enquo(protein_id_column)))

  # Check if we're dealing with metabolite data (no protein IDs)
  # This handles the case where the input table is for metabolites and doesn't have Protein.Ids or uniprot_acc
  is_metabolite_data <- !any(grepl("Protein.Ids|uniprot_acc", names(input_table)))

  if (is_metabolite_data) {
    message("   Processing metabolite data (no UniProt ID mapping)...")
    proteomics_volcano_tbl <- input_table |>
      dplyr::select({{ protein_id_column }}, {{ fdr_column }}, {{ log2FC_column }}) |>
      mutate(
        colour = case_when(
          {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} > 0 ~ "red",
          {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} < 0 ~ "blue",
          TRUE ~ "grey"
        ),
        lqm = -log10({{ fdr_column }}),
        label = case_when(
          {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} > 0 ~ "Significant Increase",
          {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} < 0 ~ "Significant Decrease",
          TRUE ~ "Not significant"
        ),
        label = factor(label, levels = c(
          "Significant Increase",
          "Significant Decrease",
          "Not significant"
        )),
        rank_positive = case_when(
          {{ log2FC_column }} > 0 ~ {{ fdr_column }},
          TRUE ~ NA_real_
        ) |> rank(),
        rank_negative = case_when(
          {{ log2FC_column }} < 0 ~ {{ fdr_column }},
          TRUE ~ NA_real_
        ) |> rank(),
        # For metabolites, we use the protein_id_column (usually Name) as the "gene name"
        gene_name_significant = case_when(
          {{ fdr_column }} < fdr_threshold &
            (rank_positive <= number_of_genes |
              rank_negative <= number_of_genes) ~ as.character({{ protein_id_column }}),
          TRUE ~ NA_character_
        )
      )
  } else {
    message("   Processing protein data (performing UniProt ID mapping)...")
    proteomics_volcano_tbl <- input_table |>
      dplyr::mutate(uniprot_acc_first = purrr::map_chr({{ protein_id_column }}, \(x) {
        str_split(x, ":")[[1]][1]
      })) |>
      dplyr::relocate(uniprot_acc_first, .after = temp_col_name) |>
      dplyr::select(uniprot_acc_first, {{ fdr_column }}, {{ log2FC_column }}) |>
      left_join(uniprot_table,
        by = join_by(uniprot_acc_first == {{ uniprot_protein_id_column }})
      ) |>
      mutate(colour = case_when(
        {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} > 0 ~ "red",
        {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} < 0 ~ "blue",
        TRUE ~ "grey"
      )) |>
      mutate(lqm = -log10({{ fdr_column }})) |>
      mutate(label = case_when(
        {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} > 0 ~ "Significant Increase",
        {{ fdr_column }} < fdr_threshold & {{ log2FC_column }} < 0 ~ "Significant Decrease",
        TRUE ~ "Not significant"
      )) |>
      mutate(label = factor(label, levels = c(
        "Significant Increase",
        "Significant Decrease",
        "Not significant"
      ))) |>
      mutate(rank_positive = case_when(
        {{ log2FC_column }} > 0 ~ {{ fdr_column }},
        TRUE ~ NA_real_
      ) |> rank()) |>
      mutate(rank_negative = case_when(
        {{ log2FC_column }} < 0 ~ {{ fdr_column }},
        TRUE ~ NA_real_
      ) |> rank()) |>
      mutate({{ gene_name_column }} := purrr::map_chr({{ gene_name_column }}, \(x) {
        str_split(x, " ")[[1]][1]
      })) |>
      mutate(gene_name_significant = case_when(
        {{ fdr_column }} < fdr_threshold &
          (rank_positive <= number_of_genes |
            rank_negative <= number_of_genes) ~ {{ gene_name_column }},
        TRUE ~ NA
      ))
  }

  proteomics_volcano_tbl
}

