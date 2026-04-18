# ----------------------------------------------------------------------------
# rollUpPrecursorToPeptideHelper
# ----------------------------------------------------------------------------
#' @title Rollup Precursors to Peptides
#' @description  Peptides of with charges and modifications are rolled up (summed) together
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
  if( length(which(is.na(core_utilisation))) == 0 ) {

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

