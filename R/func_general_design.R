# ============================================================================
# func_general_design.R
# ============================================================================
# Purpose: Design matrix cleaning and manipulation functions
# 
# This file contains functions for cleaning and manipulating experimental
# design matrices used across all omics types. Functions in this file are
# used by design builder modules and DE analysis workflows.
#
# Functions to extract here:
# - cleanDesignMatrix(): S4 method for cleaning design matrix
# - cleanDesignMatrixPeptide(): S4 method for cleaning peptide design matrix
# - cleanDesignMatrixCleanCategories(): Cleans design matrix categories
# - cleanDesignMatrixCleanCategoriesMap(): Maps cleaned categories
# - cleanDesignMatrixCreateEachVersusAllColumns(): Creates each vs all columns
# - Additional design matrix helper functions
#
# Dependencies:
# - dplyr, tidyr
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: cleanDesignMatrix()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Cleans design matrix for protein data
# setMethod(f = "cleanDesignMatrix", signature = "ProteinQuantitativeData", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 2: cleanDesignMatrixPeptide()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Cleans design matrix for peptide data
# setMethod(f = "cleanDesignMatrixPeptide", signature = "PeptideQuantitativeData", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 3: cleanDesignMatrixCleanCategories()
# Current location: R/helper_functions.R
# Description: Cleans design matrix categories
# cleanDesignMatrixCleanCategories <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 4: cleanDesignMatrixCleanCategoriesMap()
# Current location: R/helper_functions.R
# Description: Maps cleaned categories in design matrix
# cleanDesignMatrixCleanCategoriesMap <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 5: cleanDesignMatrixCreateEachVersusAllColumns()
# Current location: R/helper_functions.R
# Description: Creates each versus all columns in design matrix
# cleanDesignMatrixCreateEachVersusAllColumns <- function(...) {
#   # Extract from R/helper_functions.R
# }


# ----------------------------------------------------------------------------
# cleanDesignMatrixCleanCategories
# ----------------------------------------------------------------------------
#' @export
cleanDesignMatrixCleanCategories <- function(x ) {

  str_replace_all(x, ">=", "ge") |>
    str_replace_all( "<=", "le") |>
    str_replace_all( ">", "gt") |>
    str_replace_all( "<", "lt") |>
    str_replace_all( "\\+", ".POS") |>
    str_replace_all( "\\-", ".NEG") |>
    str_replace_all( " ", "\\.") |>
    str_replace_all("&", "and") |>
    str_replace_all("/", "_") |>
    str_replace_all("\\:", ".")
}


# ----------------------------------------------------------------------------
# cleanDesignMatrixCleanCategoriesMap
# ----------------------------------------------------------------------------
#' @export
cleanDesignMatrixCleanCategoriesMap <- function( input_table, column ) {
  input_table |>
    mutate( {{column}} := purrr::map_chr( {{column}}, cleanDesignMatrixCleanCategories ) )

}


# ----------------------------------------------------------------------------
# cleanDesignMatrixCreateEachVersusAllColumns
# ----------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
# Codes to format experimental design table for pairwise comparison of groups
#' @export
cleanDesignMatrixCreateEachVersusAllColumns <- function(input_table, id_cols, column ) {

  id_col_name <-  as_string(as_name(enquo(id_cols)))
  column_string <- as_string(as_name(enquo(column)))

  new_columns_tab <-  input_table |>
    mutate( my_value = TRUE ) |>
    pivot_wider( id_cols = {{id_cols}}
                 , names_from = {{column}}
                 , values_from = my_value
                 , values_fill = FALSE
                 , names_prefix = paste0(column_string, ".")  )

  new_column_names <- base::setdiff( colnames( new_columns_tab), id_col_name )

  # print (id_col_name)
  # print(new_column_names)

  return_table <- input_table |>
    left_join( new_columns_tab
               , by = join_by( {{id_cols}} )) |>
    relocate( all_of(new_column_names), .after={{column}} )

  return_table

}


# ----------------------------------------------------------------------------
# cleanDesignMatrixPeptide
# ----------------------------------------------------------------------------
#'@exportMethod cleanDesignMatrixPeptide
setMethod( f ="cleanDesignMatrixPeptide"
           , signature = "PeptideQuantitativeData"
           , definition=function( theObject ) {

             samples_id_vector <- theObject@peptide_data |> distinct(!!sym(theObject@sample_id)) |> dplyr::pull(!!sym(theObject@sample_id))

             theObject@design_matrix <- data.frame( temp_sample_id = samples_id_vector )  |>
               inner_join( theObject@design_matrix
                           , by = join_by ( temp_sample_id == !!sym(theObject@sample_id)) ) |>
               dplyr::rename( !!sym(theObject@sample_id) := "temp_sample_id" ) |>
               dplyr::filter( !!sym( theObject@sample_id) %in% samples_id_vector )


             return(theObject)
           })

