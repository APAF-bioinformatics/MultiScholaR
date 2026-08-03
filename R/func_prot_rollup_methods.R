# ----------------------------------------------------------------------------
# rollUpPrecursorToPeptideHelper
# ----------------------------------------------------------------------------
#' @title Rollup Precursors to Peptides
#' @description Charge states and modified forms are summed to stripped-peptide
#'   quantities while preserving frozen experiment-wide identification counts.
#' @export
rollUpPrecursorToPeptideHelper <- function( input_table
                                      , sample_id_column = Run
                                      , protein_id_column = Protein.Ids
                                      , peptide_sequence_column = Stripped.Sequence
                                      , modified_peptide_sequence_column = Modified.Sequence
                                      , precursor_quantity_column = Precursor.Quantity
                                      , precursor_normalised_column = Precursor.Normalised
                                      , core_utilisation) {

  peptide_normalised_tbl <- NA
  use_cluster <- inherits(core_utilisation, "multidplyr_cluster")
  if (!use_cluster) {

    peptide_normalised_tbl <- input_table  |>
      group_by( {{sample_id_column}}, {{protein_id_column}}, {{peptide_sequence_column}}, {{modified_peptide_sequence_column}} ) |>
      summarise( Peptide.RawQuantity = sum( {{precursor_quantity_column}} )
                 ,  Peptide.Normalised = sum( {{precursor_normalised_column}} ) ) |>
      ungroup() |>
      group_by( {{sample_id_column}}, {{protein_id_column}}, {{peptide_sequence_column}} ) |>
      summarise( Peptide.RawQuantity = sum( Peptide.RawQuantity )
                 ,  Peptide.Normalised = sum( Peptide.Normalised )
                 ,  peptidoform_count = n()) |>
      ungroup()

  } else {
    peptide_normalised_tbl <- input_table  |>

      group_by( {{sample_id_column}}, {{protein_id_column}}, {{peptide_sequence_column}}, {{modified_peptide_sequence_column}} ) |>
      partition(core_utilisation) |>
      summarise( Peptide.RawQuantity = sum( {{precursor_quantity_column}} )
                 ,  Peptide.Normalised = sum( {{precursor_normalised_column}} ) ) |>
      collect() |>
      ungroup() |>

      group_by( {{sample_id_column}}, {{protein_id_column}}, {{peptide_sequence_column}} ) |>
      partition(core_utilisation) |>
      summarise( Peptide.RawQuantity = sum( Peptide.RawQuantity )
                 ,  Peptide.Normalised = sum( Peptide.Normalised )
                 , peptidoform_count = n() ) |>
      collect() |>
      ungroup()

  }

  protein_id_str <- rlang::as_name(rlang::enquo(protein_id_column))
  evidence_columns <- intersect(
    c(
      "identification_peptide_count",
      "identification_peptidoform_count",
      "peptides_for_protein_count",
      "peptidoforms_for_protein_count"
    ),
    names(input_table)
  )

  if (length(evidence_columns) > 0L) {
    frozen_evidence <- input_table |>
      dplyr::select(
        dplyr::all_of(protein_id_str),
        dplyr::all_of(evidence_columns)
      ) |>
      dplyr::distinct()

    conflicting_evidence <- frozen_evidence |>
      dplyr::count(.data[[protein_id_str]], name = ".evidence_rows") |>
      dplyr::filter(.data$.evidence_rows != 1L)

    if (nrow(conflicting_evidence) > 0L) {
      stop(
        "rollUpPrecursorToPeptideHelper: conflicting frozen identification evidence for one or more protein groups.",
        call. = FALSE
      )
    }

    peptide_normalised_tbl <- peptide_normalised_tbl |>
      dplyr::left_join(frozen_evidence, by = protein_id_str)
  }

  # DIA-NN Protein.Group is the quantitative key, while Protein.Ids remains
  # useful accession provenance for annotation and existing reporting code.
  if (protein_id_str != "Protein.Ids" && "Protein.Ids" %in% names(input_table)) {
    collapse_protein_ids <- function(values) {
      tokens <- unlist(strsplit(as.character(values), ";", fixed = TRUE))
      tokens <- trimws(tokens)
      tokens <- tokens[!is.na(tokens) & nzchar(tokens)]
      paste(sort(unique(tokens)), collapse = ";")
    }

    protein_id_annotations <- input_table |>
      dplyr::group_by(.data[[protein_id_str]]) |>
      dplyr::summarise(
        Protein.Ids = collapse_protein_ids(.data$Protein.Ids),
        .groups = "drop"
      )

    peptide_normalised_tbl <- peptide_normalised_tbl |>
      dplyr::left_join(protein_id_annotations, by = protein_id_str)
  }

  peptide_normalised_tbl
}

# ----------------------------------------------------------------------------
# rollUpPrecursorToPeptide
# ----------------------------------------------------------------------------
#'@export
setMethod(f="rollUpPrecursorToPeptide"
          , signature="PeptideQuantitativeData"
          , definition=function (theObject, core_utilisation = NULL) {

            peptide_data <- theObject@peptide_data
            protein_id_column <- theObject@protein_id_column
            peptide_sequence_column <- theObject@peptide_sequence_column
            q_value_column <- theObject@q_value_column
            global_q_value_column <- theObject@global_q_value_column
            proteotypic_peptide_sequence_column <- theObject@proteotypic_peptide_sequence_column
            raw_quantity_column <- theObject@raw_quantity_column
            norm_quantity_column <- theObject@norm_quantity_column

            is_logged_data <- theObject@is_logged_data

            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id
            group_id <- theObject@group_id
            technical_replicate_id <- theObject@technical_replicate_id

            core_utilisation <- checkParamsObjectFunctionSimplify( theObject, "core_utilisation", NA)
            theObject <- updateParamInObject(theObject, "core_utilisation")

            theObject@peptide_data <- rollUpPrecursorToPeptideHelper(input_table = peptide_data
                                                               , sample_id_column = !!sym(sample_id)
                                                               , protein_id_column = !!sym(protein_id_column)
                                                               , peptide_sequence_column = !!sym(peptide_sequence_column)
                                                               , precursor_quantity_column = !!sym(raw_quantity_column)
                                                               , precursor_normalised_column = !!sym(norm_quantity_column)
                                                               , core_utilisation = core_utilisation)

             theObject@raw_quantity_column   <- "Peptide.RawQuantity"
             theObject@norm_quantity_column <- "Peptide.Normalised"

             theObject <- cleanDesignMatrixPeptide(theObject)

            return(theObject)
          })
