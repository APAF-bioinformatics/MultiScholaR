# ============================================================================
# func_prot_rollup.R
# ============================================================================
# Purpose: Peptide-to-protein rollup and counting functions
# 
# This file contains functions for rolling up precursor/peptide data to
# protein level and counting peptides/proteins across samples. Functions in
# this file are used by mod_prot_qc_peptide_rollup.R, mod_prot_qc_protein_rollup.R,
# and related QC modules.
#
# NOTE: Rollup is part of the QC workflow but is a distinct operation,
# hence the separate file name without "qc" prefix.
#
# Functions to extract here:
# - rollUpPrecursorToPeptide(): S4 method for precursor-to-peptide rollup
# - rollUpPrecursorToPeptideHelper(): Helper for rollup operations
# - calcPeptidesPerProtein(): Calculate peptides per protein
# - calcTotalPeptides(): Calculate total unique peptides
# - countPeptidesPerRun(): Count peptides per run/sample
# - countProteinsPerRun(): Count proteins per run/sample
# - Additional counting and rollup helper functions
#
# Dependencies:
# - dplyr, tidyr
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: rollUpPrecursorToPeptide()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Rolls up precursor data to peptide level
# setMethod(f = "rollUpPrecursorToPeptide", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 2: rollUpPrecursorToPeptideHelper()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Helper function for precursor-to-peptide rollup
# rollUpPrecursorToPeptideHelper <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 3: calcPeptidesPerProtein()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Calculates number of peptides per protein
# calcPeptidesPerProtein <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 4: calcTotalPeptides()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Calculates total unique peptides
# calcTotalPeptides <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 5: countPeptidesPerRun()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Counts peptides per run/sample
# countPeptidesPerRun <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 6: countProteinsPerRun()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Counts proteins per run/sample
# countProteinsPerRun <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 7: countUniqueProteins()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Counts unique proteins
# countUniqueProteins <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 8: count_num_peptides()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Counts number of peptides (helper)
# count_num_peptides <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 9: count_num_proteins()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Counts number of proteins (helper)
# count_num_proteins <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 10: count_num_samples()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Counts number of samples (helper)
# count_num_samples <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }


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
# calcPeptidesPerProtein
# ----------------------------------------------------------------------------
#' Calculate peptides per protein
#' 
#' @param data A data frame or S4 object containing protein/peptide data
#' @return Data frame with protein IDs and peptide counts
#' @export
calcPeptidesPerProtein <- function(data) {
  # For protein quantification data, return empty data frame
  if (isS4(data)) {
    if ("protein_quant_table" %in% slotNames(data)) {
      return(data.frame(Protein.Ids = character(), 
                       n_peptides = integer()))
    }
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |>
             group_by(Protein.Ids) |>
             summarise(n_peptides = n_distinct(Stripped.Sequence), 
                      .groups = "drop"))
    }
  }
  
  # For regular dataframes, check if it's protein quantification data
  if ("Protein.Ids" %in% names(data)) {
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      return(data.frame(Protein.Ids = character(), 
                       n_peptides = integer()))
    }
    
    if ("Stripped.Sequence" %in% names(data)) {
      return(data |>
             group_by(Protein.Ids) |>
             summarise(n_peptides = n_distinct(Stripped.Sequence), 
                      .groups = "drop"))
    }
  }
  stop("Required columns not found")
}


# ----------------------------------------------------------------------------
# calcTotalPeptides
# ----------------------------------------------------------------------------
#' Calculate total unique peptides
#' 
#' @param data A data frame or S4 object containing protein/peptide data
#' @return Integer count of unique peptide-protein combinations
#' @export
calcTotalPeptides <- function(data) {
  # For protein quantification data, return NA
  if (isS4(data)) {
    if ("protein_quant_table" %in% slotNames(data)) {
      return(NA_integer_)
    }
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |>
             distinct(Protein.Ids, Stripped.Sequence) |>
             nrow())
    }
  }
  
  # For regular dataframes, check if it's protein quantification data
  if ("Protein.Ids" %in% names(data)) {
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      return(NA_integer_)
    }
    
    if ("Stripped.Sequence" %in% names(data)) {
      return(distinct(data, Protein.Ids, Stripped.Sequence) |> nrow())
    }
  }
  stop("Required columns not found")
}


# ----------------------------------------------------------------------------
# countPeptidesPerRun
# ----------------------------------------------------------------------------
#' Count peptides per run
#' 
#' @param data A data frame or S4 object containing protein/peptide data
#' @return Data frame with run IDs and peptide counts
#' @export
countPeptidesPerRun <- function(data) {
  # For protein quantification data, return empty data frame
  if (isS4(data)) {
    if ("protein_quant_table" %in% slotNames(data)) {
      return(data.frame(Run = character(), 
                       n_peptides = integer()))
    }
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |>
             group_by(Run) |>
             summarise(n_peptides = n_distinct(Stripped.Sequence), 
                      .groups = "drop") |>
             arrange(Run))
    }
  }
  
  # For regular dataframes, check if it's protein quantification data
  if ("Protein.Ids" %in% names(data)) {
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      return(data.frame(Run = character(), 
                       n_peptides = integer()))
    }
    
    if (all(c("Run", "Stripped.Sequence") %in% names(data))) {
      return(data |>
             group_by(Run) |>
             summarise(n_peptides = n_distinct(Stripped.Sequence), 
                      .groups = "drop") |>
             arrange(Run))
    }
  }
  stop("Required columns not found")
}


# ----------------------------------------------------------------------------
# count_num_peptides
# ----------------------------------------------------------------------------
# Count the number of peptides in the input table
#' @export
count_num_peptides <- function( input_table
                                , protein_id_column = Protein.Ids
                                , peptide_sequence_column = Stripped.Sequence ) {
  num_peptides <- input_table |>
    distinct( {{protein_id_column}}, {{peptide_sequence_column}}) |>
    count()

  num_peptides[[1,1]]
}


# ----------------------------------------------------------------------------
# countProteinsPerRun
# ----------------------------------------------------------------------------
#' Count proteins per run
#' 
#' @param data A data frame or S4 object containing protein/peptide data
#' @return Data frame with run IDs and protein counts
#' @export
countProteinsPerRun <- function(data) {
  if (isS4(data)) {
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |>
             group_by(Run) |>
             summarise(n_proteins = n_distinct(Protein.Ids), 
                      .groups = "drop") |>
             arrange(Run))
    }
    if ("protein_quant_table" %in% slotNames(data)) {
      data <- data@protein_quant_table
      run_cols <- setdiff(names(data), "Protein.Ids")
      
      # For each run (column), count non-NA values
      result <- data.frame(
        Run = run_cols,
        n_proteins = sapply(run_cols, function(col) {
          sum(!is.na(data[[col]]))
        })
      ) |> arrange(Run)
      
      return(result)
    }
  }
  
  # For regular dataframes
  if ("Protein.Ids" %in% names(data)) {
    # Check if it's a protein quantification table
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      run_cols <- setdiff(names(data), "Protein.Ids")
      
      # For each run (column), count non-NA values
      result <- data.frame(
        Run = run_cols,
        n_proteins = sapply(run_cols, function(col) {
          sum(!is.na(data[[col]]))
        })
      ) |> arrange(Run)
      
      return(result)
    }
    
    # For peptide data
    if ("Run" %in% names(data)) {
      return(data |>
             group_by(Run) |>
             summarise(n_proteins = n_distinct(Protein.Ids), 
                      .groups = "drop") |>
             arrange(Run))
    }
  }
  stop("Required columns not found")
}


# ----------------------------------------------------------------------------
# countUniqueProteins
# ----------------------------------------------------------------------------
#' Count unique proteins in peptide or protein data
#' 
#' @param data A data frame or S4 object containing protein/peptide data
#' @return Integer count of unique proteins
#' @export
countUniqueProteins <- function(data) {
  if (isS4(data)) {
    if ("peptide_data" %in% slotNames(data)) {
      return(data@peptide_data |> 
             distinct(Protein.Ids) |> 
             nrow())
    }
    if ("protein_quant_table" %in% slotNames(data)) {
      return(nrow(data@protein_quant_table))
    }
  }
  
  # For regular dataframes
  if ("Protein.Ids" %in% names(data)) {
    # Check if it's a protein quantification table
    if (all(sapply(data[setdiff(names(data), "Protein.Ids")], is.numeric))) {
      return(nrow(data))  # Each row is a unique protein
    }
    return(distinct(data, Protein.Ids) |> nrow())
  }
  stop("No Protein.Ids column found")
}


# ----------------------------------------------------------------------------
# count_num_proteins
# ----------------------------------------------------------------------------
# Count the number of peptides in the input table
#' @export
count_num_proteins <- function( input_table
                                , protein_id_column = Protein.Ids) {
  num_proteins <- input_table |>
    distinct( {{protein_id_column}}) |>
    count()

  num_proteins[[1,1]]
}


# ----------------------------------------------------------------------------
# count_num_samples
# ----------------------------------------------------------------------------
# Count the number of samples in the input table
#' @export
count_num_samples <- function( input_table
                               , sample_id_column = Run) {
  num_samples <- input_table |>
    distinct( {{sample_id_column}}) |>
    count()

  num_samples[[1,1]]
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

