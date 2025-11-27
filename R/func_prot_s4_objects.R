# ============================================================================
# func_prot_s4_objects.R
# ============================================================================
# Purpose: Proteomics S4 class definitions and methods
# 
# This file contains S4 class definitions and methods specific to proteomics,
# including PeptideQuantitativeData, ProteinQuantitativeData classes,
# their constructors, and proteomics-specific S4 methods.
#
# Consolidated from:
# - peptideVsSamplesS4Objects.R (PeptideQuantitativeData class + 23 methods)
# - proteinVsSamplesS4Objects.R (ProteinQuantitativeData method)
# - protein_de_analysis_wrapper.R (DE analysis methods)
#
# Dependencies:
# - methods package
# - func_general_s4_generics.R (for generic definitions)
# ============================================================================

# ==========================================
# ProteinQuantitativeData S4 Class
# ==========================================

#' ProteinQuantitativeData S4 Class
#'
#' @description
#' An S4 class to store and manage protein-level quantitative data
#' from proteomics experiments along with experimental design metadata.
#'
#' @slot protein_quant_table A data.frame containing protein quantification data
#' @slot protein_id_column Character string naming the protein ID column
#' @slot design_matrix A data.frame containing experimental design
#' @slot protein_id_table A data.frame containing protein ID information
#' @slot sample_id Character string naming the sample ID column
#' @slot group_id Character string naming the group column
#' @slot technical_replicate_id Character string naming the replicate column
#' @slot args A list of additional arguments
#'
#' @exportClass ProteinQuantitativeData
ProteinQuantitativeData <- setClass("ProteinQuantitativeData"

         , slots = c(

                      # Protein vs Sample quantitative data

                      protein_quant_table = "data.frame"

                      , protein_id_column = "character"



                      # Design Matrix Information

                      , design_matrix = "data.frame"

                      , protein_id_table = "data.frame"

                      , sample_id="character"

                      , group_id="character"

                      , technical_replicate_id="character"

                      , args = "list")



         , prototype = list(

           # Protein vs Sample quantitative data

           protein_id_column = "Protein.Ids"

          , protein_id_table = data.frame()

           # Design Matrix Information

           , sample_id="Sample_id"

           , group_id="group"

           , technical_replicate_id="replicates"

           , args = NULL

           )



         , validity = function(object) {

           if( !is.data.frame(object@protein_quant_table) ) {

             stop("protein_quant_table must be a data.frame")

           }

           if( !is.character(object@protein_id_column) ) {

             stop("protein_id_column must be a character")

           }

           if( !is.data.frame(object@design_matrix) ) {

             stop("design_matrix must be a data.frame")

           }

           if( !is.character(object@sample_id) ) {

             stop("sample_id must be a character")

           }

           if( !is.character(object@group_id) ) {

             stop("group_id must be a character")

           }

           if( !is.character(object@technical_replicate_id) ) {

             stop("technical_replicate_id must be a character")

           }



            if( ! object@protein_id_column %in% colnames(object@protein_quant_table) ) {

                stop("Protein ID column must be in the protein data table")

            }



           if( ! object@sample_id %in% colnames(object@design_matrix) ) {

             stop("Sample ID column must be in the design matrix")

           }





           #Need to check the rows names in design matrix and the column names of the data table

           samples_in_protein_quant_table <- setdiff(colnames( object@protein_quant_table), object@protein_id_column)

           samples_in_design_matrix <- object@design_matrix |> dplyr::pull( !! sym( object@sample_id ) )



           if( length( which( sort(samples_in_protein_quant_table) != sort(samples_in_design_matrix) )) > 0 ) {

             stop("Samples in protein data and design matrix must be the same" )

           }



         }



)



# Initialize method to ensure design_matrix's sample_id column is character

setMethod("initialize", "ProteinQuantitativeData",

  function(.Object, ...) {

    # Capture all arguments passed to the constructor

    args_list <- list(...)



    # If design_matrix and sample_id are provided in the arguments

    if (!is.null(args_list$design_matrix) && !is.null(args_list$sample_id)) {

      if (args_list$sample_id %in% names(args_list$design_matrix)) {

        args_list$design_matrix[[args_list$sample_id]] <- as.character(args_list$design_matrix[[args_list$sample_id]])

        logger::log_debug(

          "Initialize ProteinQuantitativeData: Converted column '{args_list$sample_id}' in design_matrix to character."

        )

      } else {

        # This case should ideally be caught by validity, but good to be aware

        logger::log_warn(

          "Initialize ProteinQuantitativeData: Sample ID column '{args_list$sample_id}' not found in the provided design_matrix during initialization."

        )

      }

    }

    

    # Reconstruct the call to the default initializer with potentially modified args_list

    # This uses do.call to pass arguments correctly after modification

    .Object <- do.call(callNextMethod, c(list(.Object), args_list))

    

    # Explicitly return the object

    return(.Object)

  }

)

# ==========================================
# Content from peptideVsSamplesS4Objects.R
# ==========================================
## Create S4 class for proteomics protein level abundance data
#'@exportClass PeptideQuantitativeData
PeptideQuantitativeData <- setClass("PeptideQuantitativeData"

                                     , slots = c(
                                       # Protein vs Sample quantitative data
                                       peptide_data = "data.frame"
                                       , peptide_matrix = "matrix"
                                       , protein_id_column = "character"
                                       , peptide_sequence_column = "character"
                                       , q_value_column = "character"
                                       , global_q_value_column = "character"
                                       , proteotypic_peptide_sequence_column = "character"
                                       , raw_quantity_column = "character"
                                       , norm_quantity_column = "character"
                                       , is_logged_data = "logical"

                                       # Design Matrix Information
                                       , design_matrix = "data.frame"
                                       , sample_id="character"
                                       , group_id="character"
                                       , technical_replicate_id="character"
                                       , args = "list"
                                     )

                                     , prototype = list(
                                       # Protein vs Sample quantitative data
                                       peptide_matrix = matrix(numeric(0), nrow = 0, ncol = 0)
                                       , protein_id_column = "Protein_Ids"
                                       , peptide_sequence_column = "Stripped.Sequence"
                                       , q_value_column = "Q.Value"
                                       , global_q_value_column = "Global.Q.Value"
                                       , proteotypic_peptide_sequence_column = "Proteotypic"
                                       , raw_quantity_column = "Precursor.Quantity"
                                       , norm_quantity_column = "Precursor.Normalised"
                                       , is_logged_data = FALSE

                                       # Design Matrix Information
                                       , sample_id="Sample_id"
                                       , group_id="group"
                                       , technical_replicate_id="replicates"

                                       # Parameters for methods and functions
                                       , args = NULL

                                     )

                                     , validity = function(object) {
                                       if( !is.data.frame(object@peptide_data) ) {
                                         stop("peptide_data must be a data.frame")
                                       }
                                       if( !is.character(object@protein_id_column) ) {
                                         stop("protein_id_column must be a character")
                                       }

                                       if( !is.character(object@peptide_sequence_column) ) {
                                         stop("peptide_sequence_column must be a character")
                                       }
                                       if( !is.character(object@q_value_column) ) {
                                         stop("q_value_column must be a character")
                                       }
                                       if(!is.character(object@global_q_value_column) ) {
                                         stop("global_q_value_column must be a character")
                                       }
                                       if(!is.character(object@proteotypic_peptide_sequence_column) ) {
                                         stop("proteotypic_peptide_sequence_column must be a character")
                                       }
                                       if(!is.character(object@raw_quantity_column)  ) {
                                         stop("raw_quantity_column must be a character")
                                       }

                                       if(!is.character(object@norm_quantity_column) ) {
                                         stop("norm_quantity_column must be a character")
                                       }

                                       if( !is.data.frame(object@design_matrix) ) {
                                         stop("design_matrix must be a data.frame")
                                       }

                                       if( !is.character(object@sample_id) ) {
                                         stop("sample_id must be a character")
                                       }

                                       if( !is.character(object@group_id) ) {
                                         stop("group_id must be a character")
                                       }

                                       if( !is.character(object@technical_replicate_id) ) {
                                         stop("technical_replicate_id must be a character")
                                       }

                                       if( ! object@protein_id_column %in% colnames(object@peptide_data) ) {
                                         print(protein_id_column)
                                         print( colnames(object@peptide_data) )
                                         stop("Protein ID column must be in the peptide data table")
                                       }

                                       if(!object@peptide_sequence_column %in% colnames(object@peptide_data) ) {
                                         stop("Peptide sequence column must be in the peptide data table")
                                       }

                                       if(!object@q_value_column %in% colnames(object@peptide_data) ) {
                                         stop("Q value column must be in the peptide data table")
                                       }

                                       if(!object@raw_quantity_column %in% colnames(object@peptide_data) &
                                          !object@norm_quantity_column %in% colnames(object@peptide_data) ) {
                                         stop("Precursor raw quantity or normalised quantity column must be in the peptide data table")
                                       }

                                       if( ! object@sample_id %in% colnames(object@design_matrix) ) {
                                         stop("Sample ID column must be in the design matrix")
                                       }


                                       #Need to check the rows names in design matrix and the column names of the data table
                                       samples_in_peptide_data <-  object@peptide_data |> distinct(!!sym(object@sample_id)) |> dplyr::pull(!!sym(object@sample_id))
                                       samples_in_design_matrix <- object@design_matrix |> dplyr::pull( !! sym( object@sample_id ) )

                                       if( length( which( sort(samples_in_peptide_data) != sort(samples_in_design_matrix) )) > 0 ) {
                                         stop("Samples in peptide data and design matrix must be the same" )
                                       }

                                     }

)


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @export
PeptideQuantitativeDataDiann <- function( peptide_data
                                          , design_matrix
                                          , sample_id = "Run"
                                          , group_id = "group"
                                          , technical_replicate_id = "replicates"
                                          , args = NA) {



  peptide_data <- new( "PeptideQuantitativeData"

    # Protein vs Sample quantitative data
    , peptide_data = peptide_data
    , protein_id_column = "Protein.Ids"
    , peptide_sequence_column = "Stripped.Sequence"
    , q_value_column = "Q.Value"
    , global_q_value_column = "Global.Q.Value"
    , proteotypic_peptide_sequence_column = "Proteotypic"
    , raw_quantity_column = "Precursor.Quantity"
    , norm_quantity_column = "Precursor.Normalised"
    , is_logged_data = FALSE

    # Design Matrix Information
    , design_matrix = design_matrix
    , sample_id= sample_id
    , group_id= group_id
    , technical_replicate_id= technical_replicate_id
    , args = args
  )

  peptide_data

}



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Format the design matrix so that only metadata for samples in the protein data are retained, and also
# sort the sample IDs in the same order as the data matrix


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
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------



#'@export
setMethod( f ="srlQvalueProteotypicPeptideClean"
           , signature="PeptideQuantitativeData"
           , definition=function ( theObject
                                  , qvalue_threshold = NULL
                                  , global_qvalue_threshold = NULL
                                  , choose_only_proteotypic_peptide = NULL
                                  , input_matrix_column_ids =  NULL
                                  ) {
             peptide_data <- theObject@peptide_data
             protein_id_column <- theObject@protein_id_column
             q_value_column <- theObject@q_value_column
             global_q_value_column <- theObject@global_q_value_column
             peptide_sequence_column <- theObject@peptide_sequence_column
             proteotypic_peptide_sequence_column <- theObject@proteotypic_peptide_sequence_column
             raw_quantity_column <- theObject@raw_quantity_column
             norm_quantity_column <- theObject@norm_quantity_column

             qvalue_threshold <- checkParamsObjectFunctionSimplify( theObject, "qvalue_threshold", 0.01)

             global_qvalue_threshold <- checkParamsObjectFunctionSimplify( theObject, "global_qvalue_threshold", 0.01)

             choose_only_proteotypic_peptide <- checkParamsObjectFunctionSimplify( theObject
                                                                                   , "choose_only_proteotypic_peptide"
                                                                                   , 1 )

             input_matrix_column_ids <- checkParamsObjectFunctionSimplify( theObject
                                                                           , "input_matrix_column_ids" )

             theObject <- updateParamInObject(theObject, "qvalue_threshold")
             theObject <- updateParamInObject(theObject, "global_qvalue_threshold")
             theObject <- updateParamInObject(theObject, "choose_only_proteotypic_peptide")

             dia_nn_default_columns <- c("Protein.Ids"
                                        , "Stripped.Sequence"
                                        , "Q.Value"
                                        , "Global.Q.Value"
                                        , "Precursor.Quantity"
                                        , "Precursor.Normalised")

             theObject <- updateParamInObject(theObject, "input_matrix_column_ids")

             # print( paste("qvalue_threshold: ", qvalue_threshold))
             search_srl_quant_cln <- srlQvalueProteotypicPeptideCleanHelper( input_table = peptide_data
                                                                       , input_matrix_column_ids = unique(c(input_matrix_column_ids
                                                                                                      , protein_id_column
                                                                                                      , peptide_sequence_column
                                                                                                      , peptide_sequence_column))
                                                                       , protein_id_column = !!sym(protein_id_column)
                                                                       , q_value_column = !!sym(q_value_column)
                                                                       , global_q_value_column = !!sym(global_q_value_column)
                                                                       , global_qvalue_threshold = global_qvalue_threshold
                                                                       , qvalue_threshold = qvalue_threshold
                                                                       , choose_only_proteotypic_peptide = choose_only_proteotypic_peptide)

             theObject@peptide_data <- search_srl_quant_cln

             theObject <- cleanDesignMatrixPeptide(theObject)

             return(theObject)
           })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


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

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setMethod( f="peptideIntensityFiltering"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject, peptides_intensity_cutoff_percentile = NULL, peptides_proportion_of_samples_below_cutoff = NULL, core_utilisation = NULL) {
             print("--- Entering peptideIntensityFiltering S4 Method ---")
             
             peptide_data <- theObject@peptide_data
             raw_quantity_column <- theObject@raw_quantity_column
             norm_quantity_column <- theObject@norm_quantity_column

             print("   peptideIntensityFiltering: Extracting input data...")
             print(sprintf("      Arg: raw_quantity_column = %s", raw_quantity_column))
             print(sprintf("      Arg: norm_quantity_column = %s", norm_quantity_column))
             print(sprintf("      Data State (peptide_data): Dims = %d rows, %d cols", nrow(peptide_data), ncol(peptide_data)))

             print("   peptideIntensityFiltering: Resolving parameters with checkParamsObjectFunctionSimplify...")
             peptides_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify( theObject
                                                                                    , "peptides_intensity_cutoff_percentile")
             print(sprintf("      Resolved peptides_intensity_cutoff_percentile = %g", peptides_intensity_cutoff_percentile))

             peptides_proportion_of_samples_below_cutoff <- checkParamsObjectFunctionSimplify( theObject
                                                                                                , "peptides_proportion_of_samples_below_cutoff")
             print(sprintf("      Resolved peptides_proportion_of_samples_below_cutoff = %g", peptides_proportion_of_samples_below_cutoff))

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject, "core_utilisation", NA)
             print(sprintf("      Resolved core_utilisation = %s", ifelse(is.na(core_utilisation), "NA", as.character(core_utilisation))))

             print("   peptideIntensityFiltering: Updating parameters in S4 object...")
             theObject <- updateParamInObject(theObject, "peptides_intensity_cutoff_percentile")
             theObject <- updateParamInObject(theObject, "peptides_proportion_of_samples_below_cutoff")
             theObject <- updateParamInObject(theObject, "core_utilisation")

             print("   peptideIntensityFiltering: Calculating intensity threshold...")
             # Get non-missing values for threshold calculation
             valid_values <- peptide_data |> dplyr::pull(!!sym(raw_quantity_column))
             valid_values <- valid_values[!is.na(valid_values) & !is.nan(valid_values) & !is.infinite(valid_values)]
             
             print(sprintf("      peptideIntensityFiltering: Found %d valid intensity values", length(valid_values)))
             print(sprintf("      peptideIntensityFiltering: Valid values range: %g to %g", min(valid_values, na.rm=TRUE), max(valid_values, na.rm=TRUE)))
             
             min_peptide_intensity_threshold <- ceiling( quantile( peptide_data |> dplyr::pull(!!sym(raw_quantity_column)), na.rm=TRUE, probs = c(peptides_intensity_cutoff_percentile/100) ))[1]
             print(sprintf("      peptideIntensityFiltering: Calculated min_peptide_intensity_threshold = %g (percentile %g%%)", 
                          min_peptide_intensity_threshold, peptides_intensity_cutoff_percentile))

             print("   peptideIntensityFiltering: About to call helper function...")
             print(sprintf("      Helper Args: min_peptide_intensity_threshold = %g", min_peptide_intensity_threshold))
             print(sprintf("      Helper Args: peptides_proportion_of_samples_below_cutoff = %g", peptides_proportion_of_samples_below_cutoff))
             print(sprintf("      Helper Args: protein_id_column = %s", theObject@protein_id_column))
             print(sprintf("      Helper Args: peptide_sequence_column = %s", theObject@peptide_sequence_column))
             print(sprintf("      Helper Args: peptide_quantity_column = %s", raw_quantity_column))

             peptide_normalised_pif_cln <- peptideIntensityFilteringHelper( peptide_data
                                                                      , min_peptide_intensity_threshold = min_peptide_intensity_threshold
                                                                      , peptides_proportion_of_samples_below_cutoff = peptides_proportion_of_samples_below_cutoff
                                                                      , protein_id_column = !!sym( theObject@protein_id_column)
                                                                      , peptide_sequence_column = !!sym(theObject@peptide_sequence_column)
                                                                      , peptide_quantity_column = !!sym(raw_quantity_column)
                                                                      , core_utilisation = core_utilisation)

             print(sprintf("   peptideIntensityFiltering: Helper function returned. New dims = %d rows, %d cols", 
                          nrow(peptide_normalised_pif_cln), ncol(peptide_normalised_pif_cln)))

             theObject@peptide_data <- peptide_normalised_pif_cln

             print("   peptideIntensityFiltering: Cleaning design matrix...")
             theObject <- cleanDesignMatrixPeptide(theObject)

             print("--- Exiting peptideIntensityFiltering S4 Method ---")
             return(theObject)
           })

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setMethod( f = "removePeptidesWithMissingValuesPercent"
           , signature="PeptideQuantitativeData"
           , definition=function( theObject
                                  , grouping_variable = NULL
                                  , groupwise_percentage_cutoff = NULL
                                  , max_groups_percentage_cutoff = NULL
                                  , peptides_intensity_cutoff_percentile = NULL) {

             peptide_data <- theObject@peptide_data
             protein_id_column <- theObject@protein_id_column
             peptide_sequence_column <- theObject@peptide_sequence_column
             raw_quantity_column <- theObject@raw_quantity_column
             norm_quantity_column <- theObject@norm_quantity_column
             sample_id <- theObject@sample_id

             design_matrix <- theObject@design_matrix

             grouping_variable <- checkParamsObjectFunctionSimplify( theObject
                                                                   , "grouping_variable"
                                                                   , NULL)
             groupwise_percentage_cutoff <- checkParamsObjectFunctionSimplify( theObject
                                                                                   , "groupwise_percentage_cutoff"
                                                                                   , 50)
             max_groups_percentage_cutoff <- checkParamsObjectFunctionSimplify( theObject
                                                                                   , "max_groups_percentage_cutoff"
                                                                                   , 50)
             peptides_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify( theObject
                                                                                    , "peptides_intensity_cutoff_percentile"
                                                                                    , 50)

             theObject <- updateParamInObject(theObject, "grouping_variable")
             theObject <- updateParamInObject(theObject, "groupwise_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "max_groups_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "peptides_intensity_cutoff_percentile")

             min_protein_intensity_threshold <- ceiling( quantile( peptide_data |>
                                                                     dplyr::filter( !is.nan(!!sym(norm_quantity_column)) & !is.infinite(!!sym(norm_quantity_column))) |>
                                                                     dplyr::pull(!!sym(norm_quantity_column))
                                                                   , na.rm=TRUE
                                                                   , probs = c(peptides_intensity_cutoff_percentile/100) ))[1]

             # print(min_protein_intensity_threshold )

             theObject@peptide_data <- removePeptidesWithMissingValuesPercentHelper( peptide_data
                                                                               , design_matrix = design_matrix
                                                                               , sample_id = !!sym(sample_id)
                                                                               , protein_id_column = !!sym(protein_id_column)
                                                                               , peptide_sequence_column = !!sym(peptide_sequence_column)
                                                                               , grouping_variable = !!sym(grouping_variable)
                                                                               , groupwise_percentage_cutoff = groupwise_percentage_cutoff
                                                                               , max_groups_percentage_cutoff = max_groups_percentage_cutoff
                                                                               , abundance_threshold = peptides_intensity_cutoff_percentile
                                                                               , abundance_column =  norm_quantity_column )


             theObject <- cleanDesignMatrixPeptide(theObject)

             return(theObject)

           })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------




#'@title Filter the proteins based on the number of peptides and peptidoforms
#'@description Keep the proteins only if they have two or more peptides.
#'@param theObject Object of class PeptideQuantitativeData
#'@param num_peptides_per_protein_thresh Minimum number of peptides per protein
#'@param num_peptidoforms_per_protein_thresh Minimum number of peptidoforms per protein
#'@param core_utilisation core_utilisation to use for parallel processing
#'@export
setMethod( f="filterMinNumPeptidesPerProtein"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject, ... ) {
             
             # Extract specific parameters from ...
             args <- list(...)
             num_peptides_per_protein_thresh <- args$num_peptides_per_protein_thresh
             num_peptidoforms_per_protein_thresh <- args$num_peptidoforms_per_protein_thresh
             core_utilisation <- args$core_utilisation
             peptide_data <- theObject@peptide_data
             protein_id_column <- theObject@protein_id_column

             num_peptides_per_protein_thresh <- checkParamsObjectFunctionSimplify( theObject
                                                                                   , "num_peptides_per_protein_thresh"
                                                                                   , 1)

             num_peptidoforms_per_protein_thresh <- checkParamsObjectFunctionSimplify( theObject
                                                                                       , "num_peptidoforms_per_protein_thresh"
                                                                                       , 2)

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject, "core_utilisation", NA)


             theObject <- updateParamInObject(theObject, "num_peptides_per_protein_thresh")
             theObject <- updateParamInObject(theObject, "num_peptidoforms_per_protein_thresh")
             theObject <- updateParamInObject(theObject, "core_utilisation")

             theObject@peptide_data <- filterMinNumPeptidesPerProteinHelper ( input_table = peptide_data
                                                                        , num_peptides_per_protein_thresh = num_peptides_per_protein_thresh
                                                                        , num_peptidoforms_per_protein_thresh = num_peptidoforms_per_protein_thresh
                                                                        , protein_id_column = !!sym(protein_id_column)
                                                                        , core_utilisation = core_utilisation)

             theObject <- cleanDesignMatrixPeptide(theObject)

             theObject
           })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setMethod( f="filterMinNumPeptidesPerSample"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject
                                    , peptides_per_sample_cutoff = NULL
                                    , core_utilisation = NULL
                                    , inclusion_list = NULL) {

             peptide_data <- theObject@peptide_data
             sample_id_column <- theObject@sample_id

             peptides_per_sample_cutoff <- checkParamsObjectFunctionSimplify( theObject
                                                                              , "peptides_per_sample_cutoff"
                                                                              , 5000)

             inclusion_list <- checkParamsObjectFunctionSimplifyAcceptNull( theObject
                                                                            , "inclusion_list"
                                                                            , NULL)

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject, "core_utilisation", NA)

             theObject <- updateParamInObject(theObject, "peptides_per_sample_cutoff")
             theObject <- updateParamInObject(theObject, "inclusion_list")
             theObject <- updateParamInObject(theObject, "core_utilisation")

             theObject@peptide_data <- filterMinNumPeptidesPerSampleHelper( peptide_data
                                            , peptides_per_sample_cutoff = peptides_per_sample_cutoff
                                            , sample_id_column = !!sym(sample_id_column)
                                            , core_utilisation
                                            , inclusion_list = inclusion_list )

             theObject <- cleanDesignMatrixPeptide(theObject)

             theObject
           })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setMethod( f="removePeptidesWithOnlyOneReplicate"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject, replicate_group_column = NULL, core_utilisation = NULL) {

             peptide_data <- theObject@peptide_data
             sample_id_column <- theObject@sample_id
             design_matrix <- theObject@design_matrix


             grouping_variable <- checkParamsObjectFunctionSimplifyAcceptNull( theObject
                                                                       , "replicate_group_column"
                                                                       , NULL)

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject
                                                           , "core_utilisation"
                                                           , NA)

             theObject <- updateParamInObject(theObject, "replicate_group_column")
             theObject <- updateParamInObject(theObject, "core_utilisation")

             theObject@peptide_data <- removePeptidesWithOnlyOneReplicateHelper( input_table = peptide_data
                                                                                             , samples_id_tbl = design_matrix
                                                                                             , input_table_sample_id_column = !!sym(sample_id_column)
                                                                                             , sample_id_tbl_sample_id_column  = !!sym(sample_id_column)
                                                                                             , replicate_group_column = !!sym(replicate_group_column)
                                                                                             , core_utilisation = core_utilisation)
             theObject <- cleanDesignMatrixPeptide(theObject)

             theObject
           })



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------




#'@export
setMethod( f="plotPeptidesProteinsCountsPerSample"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject ) {

             plotPeptidesProteinsCountsPerSampleHelper( theObject@peptide_data
                                                  , intensity_column =  !!sym( theObject@norm_quantity_column)
                                                  , protein_id_column = !!sym(theObject@protein_id_column)
                                                  , peptide_id_column = !!sym(theObject@peptide_sequence_column)
                                                  , sample_id_column = !!sym( theObject@sample_id )
                                                  , peptide_sequence_column = !!sym( theObject@peptide_sequence_column) )


           })
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#'@export
setMethod( f="peptideMissingValueImputation"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject,  imputed_value_column = NULL, proportion_missing_values = NULL, core_utilisation = NULL) {
             peptide_data <- theObject@peptide_data
             raw_quantity_column <- theObject@raw_quantity_column
             sample_id_column <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id
             design_matrix <- theObject@design_matrix


             imputed_value_column <- checkParamsObjectFunctionSimplifyAcceptNull( theObject
                                                           , "imputed_value_column"
                                                           , NULL)

             proportion_missing_values <- checkParamsObjectFunctionSimplifyAcceptNull( theObject
                                                           , "proportion_missing_values"
                                                           , NULL)

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject
                                                           , "core_utilisation"
                                                           , NA)

             theObject <- updateParamInObject(theObject, "imputed_value_column")
             theObject <- updateParamInObject(theObject, "proportion_missing_values")
             theObject <- updateParamInObject(theObject, "core_utilisation")

             peptide_values_imputed <- peptideMissingValueImputationHelper( input_table = peptide_data
                                                                      , metadata_table = design_matrix
                                                                      , quantity_to_impute_column = !!sym( raw_quantity_column )
                                                                      , imputed_value_column = !!sym(imputed_value_column)
                                                                      , core_utilisation = core_utilisation
                                                                      , input_table_sample_id_column = !!sym( sample_id_column)
                                                                      , sample_id_tbl_sample_id_column = !!sym( sample_id_column)
                                                                      , replicate_group_column = !!sym( replicate_group_column)
                                                                      , proportion_missing_values = proportion_missing_values )

             theObject@peptide_data <- peptide_values_imputed

             theObject <- cleanDesignMatrixPeptide(theObject)

             theObject
           })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# I want to input two peptide data objects and compare them,
# to see how the number of proteins and peptides changes and how the number of samples changed
# Use set diff or set intersect to compare the peptides, proteins, samples in the two objects
#'@export
compareTwoPeptideDataObjects <- function( object_a, object_b) {

  object_a_peptides <- object_a@peptide_data |>
    distinct(!!sym(object_a@protein_id_column), !!sym(object_a@peptide_sequence_column))

  object_b_peptides <- object_b@peptide_data |>
    distinct(!!sym(object_b@protein_id_column), !!sym(object_b@peptide_sequence_column))

  object_a_proteins <- object_a@peptide_data |>
    distinct(!!sym(object_a@protein_id_column)) |>
    dplyr::pull(!!sym(object_a@protein_id_column))

  object_b_proteins <- object_b@peptide_data |>
    distinct(!!sym(object_b@protein_id_column)) |>
    dplyr::pull(!!sym(object_b@protein_id_column))

  object_a_samples <- object_a@design_matrix |>
    distinct(!!sym(object_a@sample_id)) |>
    dplyr::pull(!!sym(object_a@sample_id))

  object_b_samples <- object_b@design_matrix |>
    distinct(!!sym(object_b@sample_id)) |>
    dplyr::pull(!!sym(object_b@sample_id))


  peptides_in_a_not_b <- nrow( dplyr::setdiff( object_a_peptides, object_b_peptides) )
  peptides_intersect_a_and_b <- nrow( dplyr::intersect( object_a_peptides, object_b_peptides) )
  peptides_in_b_not_a <- nrow(  dplyr::setdiff( object_b_peptides, object_a_peptides) )

  proteins_in_a_not_b <- length( setdiff( object_a_proteins, object_b_proteins) )
  proteins_intersect_a_and_b <- length( intersect( object_a_proteins, object_b_proteins) )
  proteins_in_b_not_a <- length( setdiff( object_b_proteins, object_a_proteins) )


  samples_in_a_not_b <- length( setdiff( object_a_samples, object_b_samples) )
  samples_intersect_a_and_b <- length( intersect( object_a_samples, object_b_samples) )
  samples_in_b_not_a <- length( setdiff( object_b_samples, object_a_samples) )

  comparisons_list <- list(  peptides = list( in_a_not_b = peptides_in_a_not_b
                                             , intersect_a_and_b = peptides_intersect_a_and_b
                                             , in_b_not_a = peptides_in_b_not_a)
                            , proteins = list( in_a_not_b = proteins_in_a_not_b
                                               , intersect_a_and_b = proteins_intersect_a_and_b
                                               , in_b_not_a = proteins_in_b_not_a)
                            , samples = list( in_a_not_b = samples_in_a_not_b
                                              , intersect_a_and_b = samples_intersect_a_and_b
                                              , in_b_not_a = samples_in_b_not_a)
  )

  comparison_tibble <- comparisons_list |>
    purrr::map_df( tibble::as_tibble) |>
    add_column( Levels = c("peptides", "proteins", "samples")) |>
    relocate( Levels, .before="in_a_not_b")

  comparison_tibble


}


#'@export
summarisePeptideObject <- function(theObject) {

  num_peptides <- theObject@peptide_data |>
    distinct(!!sym(theObject@protein_id_column), !!sym(theObject@peptide_sequence_column))

  num_proteins <- theObject@peptide_data |>
    distinct(!!sym(theObject@protein_id_column)) |>
    dplyr::pull(!!sym(theObject@protein_id_column))

  num_samples <- theObject@design_matrix |>
    distinct(!!sym(theObject@sample_id)) |>
    dplyr::pull(!!sym(theObject@sample_id))

  summary_list <- list( num_peptides = nrow(num_peptides)
                       , num_proteins = length(num_proteins)
                       , num_samples = length(num_samples))

  summary_list


}

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export

#'@export
setMethod(f="calcPeptideMatrix"
          , signature="PeptideQuantitativeData"
          , definition=function( theObject ) {
            
            peptide_data <- theObject@peptide_data
            # Use the NORMALIZED quantity column for the matrix, not the raw one.
            # This is critical for all downstream stats (limpa, etc.)
            quantity_column_to_use <- theObject@norm_quantity_column 
            sample_id_column <- theObject@sample_id
            protein_id_column <- theObject@protein_id_column
            peptide_sequence_column <- theObject@peptide_sequence_column

            peptide_quant_table <- peptide_data |>
              dplyr::select(!!sym(sample_id_column)
                            , !!sym(theObject@protein_id_column)
                            , !!sym(theObject@peptide_sequence_column)
                            , !!sym(quantity_column_to_use)) |>
              mutate( !!sym(quantity_column_to_use)  := purrr::map_dbl(!!sym(quantity_column_to_use), as.numeric )) |>
              tidyr::pivot_wider(names_from = !!sym(sample_id_column)
                                 , values_from = !!sym(quantity_column_to_use)) |>
              dplyr::rename(Protein.Ids = !!sym(theObject@protein_id_column)
                            , Stripped.Sequence = !!sym(theObject@peptide_sequence_column))
            
            
            normalised_frozen_peptide_matrix_filt <- peptide_quant_table |>
              mutate(Protein.Ids = as.character(Protein.Ids)
                     , Stripped.Sequence = as.character(Stripped.Sequence) ) |>
              mutate(peptide_ids = paste(Protein.Ids, Stripped.Sequence, sep = "%")) |>
              dplyr::select(-Protein.Ids,-Stripped.Sequence) |>
              column_to_rownames("peptide_ids") |>
              as.matrix()           
            
            theObject@peptide_matrix <- normalised_frozen_peptide_matrix_filt
                                                   
            # print(head( normalised_frozen_peptide_matrix_filt))                                                   
            
            return(theObject)
            
          } )           


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export

#'@export
setMethod(f="getNegCtrlProtAnovaPeptides"
          , signature="PeptideQuantitativeData"
          , definition=function( theObject
                                 , ruv_grouping_variable = NULL
                                 , percentage_as_neg_ctrl = NULL
                                 , num_neg_ctrl = NULL
                                 , ruv_qval_cutoff = NULL
                                 , ruv_fdr_method = NULL ) {
            
            peptide_data <- theObject@peptide_data
            raw_quantity_column <- theObject@raw_quantity_column
            sample_id_column <- theObject@sample_id
            replicate_group_column <- theObject@technical_replicate_id
            design_matrix <- theObject@design_matrix
            protein_id_column <- theObject@protein_id_column
            peptide_sequence_column <- theObject@peptide_sequence_column
            group_id <- theObject@group_id
            
            normalised_frozen_peptide_matrix_filt <- theObject@peptide_matrix
            
            ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", "replicates")
            percentage_as_neg_ctrl <- checkParamsObjectFunctionSimplify( theObject, "percentage_as_neg_ctrl", 10)
            num_neg_ctrl <- checkParamsObjectFunctionSimplify( theObject
                                                               , "num_neg_ctrl"
                                                               , round(nrow( normalised_frozen_peptide_matrix_filt) * percentage_as_neg_ctrl / 100, 0))
            ruv_qval_cutoff <- checkParamsObjectFunctionSimplify( theObject, "ruv_qval_cutoff", 0.05)
            ruv_fdr_method <- checkParamsObjectFunctionSimplify( theObject, "ruv_fdr_method", "BH")
            
            theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
            theObject <- updateParamInObject(theObject, "percentage_as_neg_ctrl")
            theObject <- updateParamInObject(theObject, "num_neg_ctrl")
            theObject <- updateParamInObject(theObject, "ruv_qval_cutoff")
            theObject <- updateParamInObject(theObject, "ruv_fdr_method")
            
            control_genes_index <- getNegCtrlProtAnovaHelper( normalised_frozen_peptide_matrix_filt[,design_matrix |> dplyr::pull(!!sym(sample_id_column)) ]
                                                              , design_matrix = design_matrix |>
                                                                column_to_rownames(sample_id_column) |>
                                                                dplyr::select( -!!sym(group_id))
                                                              , grouping_variable = ruv_grouping_variable
                                                              , percentage_as_neg_ctrl = percentage_as_neg_ctrl
                                                              , num_neg_ctrl = num_neg_ctrl
                                                              , ruv_qval_cutoff = ruv_qval_cutoff
                                                              , ruv_fdr_method = ruv_fdr_method )
            
            return(control_genes_index)
          })



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' RUV Canonical Correlation Analysis for Peptide Data
#' 
#' Performs Remove Unwanted Variation (RUV) using canonical correlation analysis
#' with DPC imputation for missing values.
#' 
#' @param theObject A PeptideQuantitativeData object
#' @param ctrl Vector of control gene indices for RUV
#' @param num_components_to_impute Number of principal components for imputation
#' @param ruv_grouping_variable Column name in design matrix for RUV grouping
#' 
#' @export
#' @title RUV Canonical Correlation for PeptideQuantitativeData
#' @name ruvCancor,PeptideQuantitativeData-method
#' @export
setMethod( f = "ruvCancor"
           , signature="PeptideQuantitativeData"
           , definition=function( theObject, ctrl= NULL, num_components_to_impute=NULL, ruv_grouping_variable = NULL) {
             
             
             peptide_matrix <- theObject@peptide_matrix
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             sample_id <- theObject@sample_id
             
             ctrl <- checkParamsObjectFunctionSimplify( theObject, "ctrl", NULL)
             num_components_to_impute <- checkParamsObjectFunctionSimplify( theObject, "num_components_to_impute", 2)
             ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", NULL)
             
             theObject <- updateParamInObject(theObject, "ctrl")
             theObject <- updateParamInObject(theObject, "num_components_to_impute")
             theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
             
             if(! ruv_grouping_variable %in% colnames(design_matrix)) {
               stop( paste0("The 'ruv_grouping_variable = "
                            , ruv_grouping_variable
                            , "' is not a column in the design matrix.") )
             }
             
             if( is.na(num_components_to_impute) || num_components_to_impute < 1) {
               stop(paste0("The num_components_to_impute = ", num_components_to_impute, " value is invalid."))
             }
             
             if( length( ctrl) < 5 ) {
               stop(paste0( "The number of negative control molecules entered is less than 5. Please check the 'ctl' parameter."))
             }

             # Remove or winsorize extreme values
             peptide_matrix_clean <- apply(log2(peptide_matrix), 2, function(x) {
               q <- quantile(x, probs = c(0.01, 0.99), na.rm=TRUE)
               x[x < q[1]] <- q[1]
               x[x > q[2]] <- q[2]
               return(x)
             })
             
             
             dpcfit <- dpc(peptide_matrix_clean)
             
             peptide_matrix_complete <- dpcImpute(peptide_matrix_clean, dpc=dpcfit)$E
             
             
             Y <-  t( peptide_matrix_complete[,design_matrix |> dplyr::pull(!!sym(sample_id))]) 

             print("steps")
             
             cancorplot_r2 <- ruv_cancorplot( Y ,
                                              X = design_matrix |>
                                                dplyr::pull(!!sym(ruv_grouping_variable)),
                                              ctl = ctrl)
             cancorplot_r2
             
             
           })

#' @title Fast version of ruvCancor for optimization
#' @name ruvCancorFast,PeptideQuantitativeData-method
#' @description This function provides a lightweight alternative to ruvCancor that skips
#' the expensive DPC imputation step, making it suitable for use during
#' optimization loops where speed is critical.
#' @param theObject A PeptideQuantitativeData object
#' @param ctrl Control features for RUV
#' @param num_components_to_impute Number of components to impute
#' @param ruv_grouping_variable Grouping variable for RUV
#' @param simple_imputation_method Method for simple missing value handling.
#'   Options: "none" (default), "mean", "median", "min"
#' @export
setMethod( f = "ruvCancorFast"
           , signature="PeptideQuantitativeData"
           , definition=function( theObject, ctrl= NULL, num_components_to_impute=NULL, 
                                  ruv_grouping_variable = NULL, simple_imputation_method = "none") {
             
             peptide_matrix <- theObject@peptide_matrix
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             sample_id <- theObject@sample_id
             
             ctrl <- checkParamsObjectFunctionSimplify( theObject, "ctrl", NULL)
             num_components_to_impute <- checkParamsObjectFunctionSimplify( theObject, "num_components_to_impute", 2)
             ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", NULL)
             
             theObject <- updateParamInObject(theObject, "ctrl")
             theObject <- updateParamInObject(theObject, "num_components_to_impute")
             theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
             
             if(! ruv_grouping_variable %in% colnames(design_matrix)) {
               stop( paste0("The 'ruv_grouping_variable = "
                            , ruv_grouping_variable
                            , "' is not a column in the design matrix.") )
             }
             
             if( is.na(num_components_to_impute) || num_components_to_impute < 1) {
               stop(paste0("The num_components_to_impute = ", num_components_to_impute, " value is invalid."))
             }
             
             if( length( ctrl) < 5 ) {
               stop(paste0( "The number of negative control molecules entered is less than 5. Please check the 'ctl' parameter."))
             }
             
             # Skip expensive DPC imputation - use simple approach for optimization
             peptide_matrix_working <- log2(peptide_matrix + 1)  # Add small constant to avoid log(0)
             
             # Handle missing values with simple method if requested
             if (simple_imputation_method != "none" && anyNA(peptide_matrix_working)) {
               if (simple_imputation_method == "mean") {
                 peptide_matrix_working <- apply(peptide_matrix_working, 2, function(x) {
                   x[is.na(x)] <- mean(x, na.rm = TRUE)
                   return(x)
                 })
               } else if (simple_imputation_method == "median") {
                 peptide_matrix_working <- apply(peptide_matrix_working, 2, function(x) {
                   x[is.na(x)] <- median(x, na.rm = TRUE)
                   return(x)
                 })
               } else if (simple_imputation_method == "min") {
                 peptide_matrix_working <- apply(peptide_matrix_working, 2, function(x) {
                   x[is.na(x)] <- min(x, na.rm = TRUE)
                   return(x)
                 })
               }
             }
             
             # Remove or winsorize extreme values (lightweight version)
             peptide_matrix_clean <- apply(peptide_matrix_working, 2, function(x) {
               q <- quantile(x, probs = c(0.01, 0.99), na.rm=TRUE)
               x[x < q[1]] <- q[1]
               x[x > q[2]] <- q[2]
               return(x)
             })
             
             # Use the matrix directly without expensive DPC imputation
             Y <- t(peptide_matrix_clean[, design_matrix |> dplyr::pull(!!sym(sample_id))])
             
             # Generate canonical correlation plot (this is what we actually need for optimization)
             cancorplot_r2 <- ruv_cancorplot( Y ,
                                              X = design_matrix |>
                                                dplyr::pull(!!sym(ruv_grouping_variable)),
                                              ctl = ctrl)
             
             return(cancorplot_r2)
           })

#'@export
setMethod( f = "ruvIII_C_Varying"
           , signature="PeptideQuantitativeData"
           , definition=function( theObject, ruv_grouping_variable = NULL, ruv_number_k = NULL, ctrl = NULL) {
             
             peptide_data <- theObject@peptide_data
             protein_id_column <- theObject@protein_id_column
             design_matrix <- theObject@design_matrix
             group_id <- theObject@group_id
             sample_id <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id
             
             ruv_grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "ruv_grouping_variable", NULL)
             k <- checkParamsObjectFunctionSimplify( theObject, "ruv_number_k", NULL)
             ctrl <- checkParamsObjectFunctionSimplify( theObject, "ctrl", NULL)
             
             theObject <- updateParamInObject(theObject, "ruv_grouping_variable")
             theObject <- updateParamInObject(theObject, "ruv_number_k")
             theObject <- updateParamInObject(theObject, "ctrl")

             # Create matrix from peptide_data (exactly like protein version)
             current_quant_column <- theObject@norm_quantity_column
             
             # Create unique peptide identifier for matrix rownames
             peptide_unique_id <- paste(peptide_data[[theObject@protein_id_column]], 
                                       peptide_data[[theObject@peptide_sequence_column]], 
                                       sep = "%")
             
             # Create wide format data frame then convert to matrix (exactly like protein version)
             temp_peptide_wide <- peptide_data |>
               mutate(peptide_id = peptide_unique_id) |>
               select(peptide_id, !!sym(sample_id), !!sym(current_quant_column)) |>
               pivot_wider(names_from = !!sym(sample_id), 
                          values_from = !!sym(current_quant_column),
                          values_fn = mean)

             # Convert to matrix exactly like protein version: column_to_rownames() then as.matrix()
             normalised_frozen_peptide_matrix_filt <- temp_peptide_wide |>
               column_to_rownames("peptide_id") |>
               as.matrix()
             
             Y <-  t( normalised_frozen_peptide_matrix_filt[,design_matrix |> dplyr::pull(!!sym(sample_id))])
             
             M <- getRuvIIIReplicateMatrixHelper( design_matrix
                                                  , !!sym(sample_id)
                                                  , !!sym(ruv_grouping_variable))
             
             cln_mat <- RUVIII_C_Varying( k = ruv_number_k
                                          , Y = Y
                                          , M = M
                                          , toCorrect = colnames(Y)
                                          , potentialControls = names( ctrl[which(ctrl)] ) )
             
             # Remove samples with no values
             cln_mat_2 <- cln_mat[rowSums(is.na(cln_mat) | is.nan(cln_mat)) != ncol(cln_mat),]
             
             # Remove proteins with no values
             cln_mat_3 <- t(cln_mat_2)
             cln_mat_4 <- cln_mat_3[rowSums(is.na(cln_mat_3) | is.nan(cln_mat_3)) != ncol(cln_mat_3),]
             
             # Convert back to data frame (exactly like protein version)
             ruv_normalised_results_cln <- cln_mat_4 |>
               as.data.frame() |>
               rownames_to_column("peptide_id")

             # Update peptide_data by joining with RUV-corrected values and replacing the quant column
             updated_peptide_data <- peptide_data |>
               mutate(peptide_id = peptide_unique_id) |>
               select(-!!sym(current_quant_column)) |>
               left_join(ruv_normalised_results_cln |> 
                        pivot_longer(cols = -peptide_id, 
                                    names_to = sample_id, 
                                    values_to = current_quant_column),
                        by = c("peptide_id", sample_id)) |>
               select(-peptide_id)

             # Update both slots
             theObject@peptide_data <- updated_peptide_data
             theObject@peptide_matrix <- cln_mat_4
             
             theObject <- cleanDesignMatrixPeptide(theObject)
             
             return( theObject )
           })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Find the Best Negative Control Percentage for RUV-III Analysis on Peptide Data
#'
#' This function automatically determines the optimal percentage of peptides to use
#' as negative controls for RUV-III analysis by testing different percentages and
#' evaluating the separation quality between "All" and "Control" groups in canonical
#' correlation plots.
#'
#' @param peptide_matrix_obj A PeptideQuantitativeData object containing
#'   the peptide quantification data
#' @param percentage_range A numeric vector specifying the range of percentages to test.
#'   Default is seq(1, 20, by = 1) for testing 1% to 20% in 1% increments
#' @param num_components_to_impute Number of components to use for imputation in ruvCancor.
#'   Default is 5
#' @param ruv_grouping_variable The grouping variable to use for RUV analysis.
#'   Default is "group"
#' @param ruv_qval_cutoff The FDR threshold for negative control selection.
#'   Default is 0.05
#' @param ruv_fdr_method The FDR calculation method. Default is "qvalue"
#' @param separation_metric The metric to use for evaluating separation quality.
#'   Options: "max_difference" (default), "mean_difference", "auc", "weighted_difference"
#' @param k_penalty_weight Weight for penalizing high k values in composite score.
#'   Default is 0.5. Higher values penalize high k more strongly
#' @param max_acceptable_k Maximum acceptable k value. k values above this get heavy penalty.
#'   Default is 3
#' @param adaptive_k_penalty Whether to automatically adjust max_acceptable_k based on sample size.
#'   Default is TRUE (recommended). Set to FALSE only if you need exact reproducibility with previous results
#' @param verbose Whether to print progress messages. Default is TRUE
#' @param ensure_matrix Whether to ensure peptide matrix is calculated. Default is TRUE
#'
#' @return A list containing:
#'   \itemize{
#'     \item best_percentage: The optimal percentage as a numeric value
#'     \item best_k: The optimal k value from findBestK() for the best percentage
#'     \item best_control_genes_index: The control genes index for the best percentage
#'     \item best_separation_score: The separation score for the best percentage
#'     \item best_composite_score: The composite score (separation penalized by k value)
#'     \item optimization_results: A data frame with all tested percentages and their scores
#'     \item best_cancor_plot: The canonical correlation plot for the best percentage
#'     \item separation_metric_used: The separation metric that was used
#'     \item k_penalty_weight: The k penalty weight that was used
#'     \item max_acceptable_k: The maximum acceptable k value that was used
#'   }
#'
#' @importFrom logger log_info log_warn
#' @importFrom purrr imap map_dfr map_dbl
#' @export

#' @export
setMethod(f="findBestNegCtrlPercentagePeptides"
          , signature="PeptideQuantitativeData"
          , definition=function(theObject,
                                percentage_range = seq(1, 20, by = 1),
                                num_components_to_impute = 5,
                                ruv_grouping_variable = "group",
                                ruv_qval_cutoff = 0.05,
                                ruv_fdr_method = "qvalue",
                                separation_metric = "max_difference",
                                k_penalty_weight = 0.5,
                                max_acceptable_k = 3,
                                adaptive_k_penalty = TRUE,
                                verbose = TRUE,
                                ensure_matrix = TRUE) {
  
  # Input validation
  if (!inherits(theObject, "PeptideQuantitativeData")) {
    stop("theObject must be a PeptideQuantitativeData object")
  }
  
  # Ensure peptide matrix is calculated if requested
  if (ensure_matrix && (!"peptide_matrix" %in% slotNames(theObject) || is.null(theObject@peptide_matrix) || length(theObject@peptide_matrix) == 0)) {
    if (verbose) {
      log_info("Peptide matrix not found. Calculating peptide matrix...")
    }
    theObject <- calcPeptideMatrix(theObject)
  }
  
  if (length(percentage_range) == 0 || any(percentage_range <= 0) || any(percentage_range > 100)) {
    stop("percentage_range must contain values between 0 and 100")
  }
  
  if (!separation_metric %in% c("max_difference", "mean_difference", "auc", "weighted_difference")) {
    stop("separation_metric must be one of: 'max_difference', 'mean_difference', 'auc', 'weighted_difference'")
  }
  
  if (k_penalty_weight < 0 || k_penalty_weight > 1) {
    stop("k_penalty_weight must be between 0 and 1")
  }
  
  if (max_acceptable_k < 1 || !is.numeric(max_acceptable_k)) {
    stop("max_acceptable_k must be a positive number >= 1")
  }
  
  if (adaptive_k_penalty && max_acceptable_k < 3) {
    stop("max_acceptable_k must be at least 3 when adaptive_k_penalty is TRUE")
  }
  
  # Calculate adaptive max_acceptable_k if requested
  if (adaptive_k_penalty) {
    # Get sample size from the peptide matrix
    sample_size <- ncol(theObject@peptide_matrix)
    
    # Calculate adaptive max_acceptable_k based on sample size
    adaptive_max_k <- .peptide_calculateAdaptiveMaxK(sample_size)
    
    if (verbose) {
      log_info("Adaptive penalty enabled: Sample size = {sample_size}, Adaptive max_acceptable_k = {adaptive_max_k} (original = {max_acceptable_k})")
    }
    
    max_acceptable_k <- adaptive_max_k
  }
  
  # Detect small datasets and warn/adjust percentage range if needed
  sample_size <- ncol(theObject@peptide_matrix)
  if (sample_size < 15 && max(percentage_range) < 30) {
    if (verbose) {
      log_warn("Small dataset detected (n={sample_size}). Consider testing higher percentages (up to 30-50%) for better negative control identification.")
      log_warn("Current range: {paste(range(percentage_range), collapse = '-')}%. May need wider range for optimal results.")
    }
  }
  
  if (verbose) {
    log_info("Starting optimization of negative control percentage for PEPTIDE data with k value consideration...")
    log_info("Testing {length(percentage_range)} different percentages: {paste(range(percentage_range), collapse = '-')}%")
    log_info("K penalty weight: {k_penalty_weight}, Max acceptable k: {max_acceptable_k}")
    if (adaptive_k_penalty) {
      log_info("Using adaptive k penalty based on sample size")
    }
  }
  
  # Process all percentages using functional programming
  if (verbose) {
    log_info("Processing {length(percentage_range)} percentages using vectorized operations...")
  }
  
  # Create a function to process a single percentage
  process_percentage <- function(current_percentage, index) {
    if (verbose && index %% 5 == 0) {
      log_info("Testing percentage {index}/{length(percentage_range)}: {current_percentage}%")
    }
    
    tryCatch({
      # Get negative control peptides for current percentage
      control_genes_index <- getNegCtrlProtAnovaPeptides(
        theObject,
        ruv_grouping_variable = ruv_grouping_variable,
        percentage_as_neg_ctrl = current_percentage,
        ruv_qval_cutoff = ruv_qval_cutoff,
        ruv_fdr_method = ruv_fdr_method
      )
      
      # Check if we have enough control peptides
      num_controls <- sum(control_genes_index, na.rm = TRUE)
      if (num_controls < 5) {
        if (verbose) {
          log_warn("Percentage {current_percentage}%: Only {num_controls} control peptides found (minimum 5 required). Skipping.")
        }
        return(list(
          percentage = current_percentage,
          separation_score = NA_real_,
          best_k = NA_real_,
          composite_score = NA_real_,
          num_controls = num_controls,
          valid_plot = FALSE,
          control_genes_index = NULL,
          cancor_plot = NULL
        ))
      }
      
      # Generate canonical correlation plot using FAST version (skips expensive DPC imputation)
      cancorplot <- ruvCancorFast(
        theObject,
        ctrl = control_genes_index,
        num_components_to_impute = num_components_to_impute,
        ruv_grouping_variable = ruv_grouping_variable,
        simple_imputation_method = "mean"  # Use simple mean imputation for speed
      )
      
      # Calculate separation score
      separation_score <- .peptide_calculateSeparationScore(cancorplot, separation_metric)
      
      # Calculate the best k using the existing findBestK function
      best_k <- tryCatch({
        findBestK(cancorplot)
      }, error = function(e) {
        if (verbose) {
          log_warn("Percentage {current_percentage}%: Error calculating best k: {e$message}")
        }
        return(NA_real_)
      })
      
      # Calculate composite score that considers both separation and k value
      composite_score <- .peptide_calculateCompositeScore(
        separation_score, 
        best_k, 
        k_penalty_weight, 
        max_acceptable_k
      )
      
      return(list(
        percentage = current_percentage,
        separation_score = separation_score,
        best_k = best_k,
        composite_score = composite_score,
        num_controls = num_controls,
        valid_plot = TRUE,
        control_genes_index = control_genes_index,
        cancor_plot = cancorplot
      ))
      
    }, error = function(e) {
      if (verbose) {
        log_warn("Percentage {current_percentage}%: Error occurred - {e$message}")
      }
      return(list(
        percentage = current_percentage,
        separation_score = NA_real_,
        best_k = NA_real_,
        composite_score = NA_real_,
        num_controls = NA_integer_,
        valid_plot = FALSE,
        control_genes_index = NULL,
        cancor_plot = NULL
      ))
    })
  }
  
  # Use purrr::imap() for functional processing
  all_results <- percentage_range |>
    purrr::imap(process_percentage)
  
  # Extract results into proper data frame
  results <- all_results |>
    purrr::map_dfr(~ data.frame(
      percentage = .x$percentage,
      separation_score = .x$separation_score,
      best_k = .x$best_k,
      composite_score = .x$composite_score,
      num_controls = .x$num_controls,
      valid_plot = .x$valid_plot
    ))
  
  # Find the best result using composite score (considers both separation and k value)
  valid_results <- all_results[!is.na(purrr::map_dbl(all_results, "composite_score"))]
  
  if (length(valid_results) == 0) {
    stop("No valid percentage found. Please check your data and parameters.")
  }
  
  best_index <- which.max(purrr::map_dbl(valid_results, "composite_score"))
  best_result <- valid_results[[best_index]]
  
  best_percentage <- best_result$percentage
  best_control_genes_index <- best_result$control_genes_index
  best_cancor_plot <- best_result$cancor_plot
  best_separation_score <- best_result$separation_score
  best_composite_score <- best_result$composite_score
  best_k <- best_result$best_k
  
  # Final validation and logging
  if (verbose) {
    log_info("Optimization complete!")
    log_info("Best percentage: {best_percentage}% (composite score: {round(best_composite_score, 4)})")
    log_info("  - Separation score: {round(best_separation_score, 4)}")
    log_info("  - Best k value: {best_k}")
    log_info("  - Number of control peptides: {sum(best_control_genes_index, na.rm = TRUE)}")
  }
  
  # Return comprehensive results
  return(list(
    best_percentage = best_percentage,
    best_k = best_k,
    best_control_genes_index = best_control_genes_index,
    best_separation_score = best_separation_score,
    best_composite_score = best_composite_score,
    optimization_results = results,
    best_cancor_plot = best_cancor_plot,
    separation_metric_used = separation_metric,
    k_penalty_weight = k_penalty_weight,
    max_acceptable_k = max_acceptable_k,
    adaptive_k_penalty_used = adaptive_k_penalty,
    sample_size = if(adaptive_k_penalty) ncol(theObject@peptide_matrix) else NA
  ))
})

#' Calculate Separation Score for Canonical Correlation Plot (Peptide version)
#'
#' Internal helper function to calculate separation quality between "All" and "Control"
#' groups in a canonical correlation plot.
#'
#' @param cancorplot A ggplot object from ruvCancor
#' @param metric The separation metric to calculate
#'
#' @return A numeric separation score (higher is better)
#'
#' @keywords internal
.peptide_calculateSeparationScore <- function(cancorplot, metric = "max_difference") {
  
  # Extract data from the plot
  if (!inherits(cancorplot, "ggplot") || is.null(cancorplot$data)) {
    return(NA_real_)
  }
  
  plot_data <- cancorplot$data
  
  # Check required columns exist
  if (!all(c("featureset", "cc", "K") %in% colnames(plot_data))) {
    return(NA_real_)
  }
  
  # Get indices for Control and All groups
  controls_idx <- which(plot_data$featureset == "Control")
  all_idx <- which(plot_data$featureset == "All")
  
  if (length(controls_idx) == 0 || length(all_idx) == 0) {
    return(NA_real_)
  }
  
  # Calculate differences between All and Control
  difference_between_all_ctrl <- plot_data$cc[all_idx] - plot_data$cc[controls_idx]
  
  # Remove any NA or infinite values
  valid_diffs <- difference_between_all_ctrl[is.finite(difference_between_all_ctrl)]
  
  if (length(valid_diffs) == 0) {
    return(NA_real_)
  }
  
  # Calculate score based on specified metric
  score <- switch(metric,
    "max_difference" = max(valid_diffs, na.rm = TRUE),
    "mean_difference" = mean(valid_diffs, na.rm = TRUE),
    "auc" = {
      # Area under the curve (trapezoidal rule approximation)
      k_values <- plot_data$K[all_idx][is.finite(difference_between_all_ctrl)]
      if (length(k_values) < 2) return(NA_real_)
      
      # Sort by K value
      sorted_idx <- order(k_values)
      k_sorted <- k_values[sorted_idx]
      diff_sorted <- valid_diffs[sorted_idx]
      
      # Calculate AUC using trapezoidal rule
      sum(diff(k_sorted) * (head(diff_sorted, -1) + tail(diff_sorted, -1)) / 2)
    },
    "weighted_difference" = {
      # Weight differences by their K value (higher K gets more weight)
      k_values <- plot_data$K[all_idx][is.finite(difference_between_all_ctrl)]
      if (length(k_values) == 0) return(NA_real_)
      
      weights <- k_values / max(k_values, na.rm = TRUE)
      sum(valid_diffs * weights, na.rm = TRUE) / sum(weights, na.rm = TRUE)
    },
    NA_real_
  )
  
  return(as.numeric(score))
}

#' Calculate Composite Score for Percentage Optimization (Peptide version)
#'
#' Internal helper function to calculate a composite score that considers both
#' separation quality and the resulting k value from findBestK(). This prevents
#' over-optimization towards percentages that give good separation but unreasonably
#' high k values that would remove biological signal.
#'
#' @param separation_score The separation score from calculateSeparationScore()
#' @param best_k The best k value from findBestK()
#' @param k_penalty_weight Weight for k penalty (0-1). Higher values penalize high k more
#' @param max_acceptable_k Maximum acceptable k value. k values above this get heavy penalty
#'
#' @return A numeric composite score (higher is better)
#'
#' @keywords internal
.peptide_calculateCompositeScore <- function(separation_score, best_k, k_penalty_weight, max_acceptable_k) {
  
  # Handle NA cases
  if (is.na(separation_score) || is.na(best_k)) {
    return(NA_real_)
  }
  
  # Ensure positive values
  if (separation_score <= 0) {
    return(0)
  }
  
  # Calculate k penalty
  if (best_k <= max_acceptable_k) {
    # Linear penalty within acceptable range: penalty = 0 at k=1, penalty = k_penalty_weight at k=max_acceptable_k
    k_penalty <- k_penalty_weight * (best_k - 1) / (max_acceptable_k - 1)
  } else {
    # Heavy penalty for k values above max_acceptable_k
    # Exponential penalty: starts at k_penalty_weight and increases rapidly
    excess_k <- best_k - max_acceptable_k
    k_penalty <- k_penalty_weight + (1 - k_penalty_weight) * (1 - exp(-excess_k))
  }
  
  # Ensure k_penalty is between 0 and 1
  k_penalty <- pmax(0, pmin(1, k_penalty))
  
  # Calculate composite score: separation_score * (1 - k_penalty)
  # This means:
  # - k=1: no penalty (multiply by 1)
  # - k=max_acceptable_k: penalty = k_penalty_weight (multiply by 1-k_penalty_weight)
  # - k>max_acceptable_k: heavy penalty (multiply by value approaching 0)
  composite_score <- separation_score * (1 - k_penalty)
  
  return(as.numeric(composite_score))
}

#' Calculate Adaptive Maximum Acceptable K Based on Sample Size (Peptide version)
#'
#' Internal helper function to determine an appropriate max_acceptable_k value
#' based on the number of samples in the dataset. This prevents over-correction
#' in small datasets and allows more flexibility in large datasets.
#'
#' @param sample_size Number of samples in the dataset
#'
#' @return An integer representing the adaptive max_acceptable_k value
#'
#' @details
#' The adaptive calculation follows these principles:
#' - Small datasets (n < 15): Conservative approach, max_k = 2
#' - Medium datasets (n = 15-40): Standard approach, max_k = 3  
#' - Large datasets (n = 40-80): Moderate approach, max_k = 4
#' - Very large datasets (n > 80): Permissive approach, max_k = 5
#' 
#' The rationale is that with more samples, you have more degrees of freedom
#' and statistical power, making higher k values less problematic.
#'
#' @keywords internal
.peptide_calculateAdaptiveMaxK <- function(sample_size) {
  
  if (sample_size < 15) {
    # Small datasets: be very conservative
    # Each k factor consumes significant degrees of freedom
    return(2L)
  } else if (sample_size < 40) {
    # Medium datasets: standard approach
    # This is the typical proteomics experiment size
    return(3L)
  } else if (sample_size < 80) {
    # Large datasets: can afford one extra k factor
    # Sufficient statistical power to handle k=4
    return(4L)
  } else {
    # Very large datasets: most permissive
    # Abundant statistical power allows k=5 if separation justifies it
    return(5L)
  }
}

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Peptide Missing Value Imputation using limpa Package
#'
#' This function uses the limpa package's detection probability curve (DPC) approach
#' for sophisticated missing value imputation specifically designed for proteomics data.
#' This method is more robust than traditional imputation as it models the missing value
#' mechanism based on detection probabilities.
#'
#' @param theObject A PeptideQuantitativeData object
#' @param imputed_value_column Name for the new column containing imputed values.
#'   Default is "Peptide.Imputed.Limpa"
#' @param use_log2_transform Whether to log2 transform the data before imputation.
#'   Default is TRUE (recommended by limpa)
#' @param verbose Whether to print progress messages. Default is TRUE
#' @param ensure_matrix Whether to ensure peptide matrix is calculated. Default is TRUE
#'
#' @details
#' The limpa package uses a detection probability curve (DPC) to model the relationship
#' between peptide intensity and the probability of detection. This allows for more
#' sophisticated imputation that accounts for the intensity-dependent nature of missing
#' values in proteomics data, rather than assuming they are missing at random.
#'
#' The process follows these steps:
#' 1. Estimate the detection probability curve using dpc()
#' 2. Perform row-wise imputation using dpcImpute()
#' 3. Transform results back to original scale if needed
#'
#' @return Updated PeptideQuantitativeData object with imputed values
#'
#' @importFrom limpa dpc dpcImpute
#' @export

#' @export
setMethod(f="peptideMissingValueImputationLimpa"
          , signature="PeptideQuantitativeData"
          , definition = function(theObject, 
                                  imputed_value_column = NULL, 
                                  use_log2_transform = TRUE,
                                  verbose = TRUE,
                                  ensure_matrix = TRUE) {
            
            # Load required packages
            if (!requireNamespace("limpa", quietly = TRUE)) {
              stop("limpa package is required but not installed. Please install it using: BiocManager::install('limpa')")
            }
            
            # Parameter validation and defaults
            imputed_value_column <- checkParamsObjectFunctionSimplifyAcceptNull(
              theObject, "imputed_value_column", "Peptide.Imputed.Limpa"
            )
            
            use_log2_transform <- checkParamsObjectFunctionSimplify(
              theObject, "use_log2_transform", TRUE
            )
            
            verbose <- checkParamsObjectFunctionSimplify(
              theObject, "verbose", TRUE
            )
            
            # Update parameters in object
            theObject <- updateParamInObject(theObject, "imputed_value_column")
            theObject <- updateParamInObject(theObject, "use_log2_transform")
            theObject <- updateParamInObject(theObject, "verbose")
            
            # Ensure peptide matrix is calculated if requested
            if (ensure_matrix && (!"peptide_matrix" %in% slotNames(theObject) || 
                                  is.null(theObject@peptide_matrix) || 
                                  length(theObject@peptide_matrix) == 0)) {
              if (verbose) {
                log_info("Peptide matrix not found. Calculating peptide matrix...")
              }
              theObject <- calcPeptideMatrix(theObject)
            }
            
            # Extract data
            peptide_data <- theObject@peptide_data
            peptide_matrix <- theObject@peptide_matrix
            raw_quantity_column <- theObject@raw_quantity_column
            sample_id_column <- theObject@sample_id
            design_matrix <- theObject@design_matrix
            
            if (verbose) {
              log_info("Starting limpa-based missing value imputation...")
              log_info("Data dimensions: {nrow(peptide_matrix)} peptides x {ncol(peptide_matrix)} samples")
              log_info("Missing value percentage: {round(100 * mean(is.na(peptide_matrix)), 1)}%")
            }
            
            # Prepare data for limpa (peptides as rows, samples as columns)
            # limpa expects log2-transformed data
            y_peptide <- peptide_matrix
            
            # Transform to log2 if requested and data is not already log-transformed
            if (use_log2_transform && !theObject@is_logged_data) {
              if (verbose) {
                log_info("Applying log2 transformation...")
              }
              # Add small constant to avoid log(0)
              y_peptide <- log2(y_peptide + 1)
            } else if (use_log2_transform && theObject@is_logged_data) {
              if (verbose) {
                log_warn("Data already log2 transformed, skipping additional transformation")
              }
              # Data already log2, use as-is
            } else if (!use_log2_transform && !theObject@is_logged_data) {
              if (verbose) {
                log_info("Converting raw intensities to log2 scale for limpa...")
              }
              # limpa expects log2 data, so transform raw data
              y_peptide <- log2(y_peptide + 1)
            } else {
              # !use_log2_transform && theObject@is_logged_data
              if (verbose) {
                log_info("Using existing log2 transformed data (no additional transformation)")
              }
              # Data already log2, use as-is - this is the correct case!
            }
            
            # Check for infinite or NaN values
            if (any(is.infinite(y_peptide) | is.nan(y_peptide), na.rm = TRUE)) {
              if (verbose) {
                log_warn("Infinite or NaN values detected. Replacing with NA...")
              }
              y_peptide[is.infinite(y_peptide) | is.nan(y_peptide)] <- NA
            }
            
            # Estimate Detection Probability Curve
            if (verbose) {
              log_info("Estimating detection probability curve...")
            }
            
            tryCatch({
              dpcfit <- limpa::dpc(y_peptide)
              
              if (verbose) {
                log_info("DPC parameters estimated:")
                log_info("  beta0 (intercept): {round(dpcfit$dpc[1], 4)}")
                log_info("  beta1 (slope): {round(dpcfit$dpc[2], 4)}")
                
                # Interpret the slope
                slope_interpretation <- if (dpcfit$dpc[2] < 0.3) {
                  "nearly random missing"
                } else if (dpcfit$dpc[2] < 0.7) {
                  "moderate intensity-dependent missing"
                } else if (dpcfit$dpc[2] < 1.2) {
                  "strong intensity-dependent missing"
                } else {
                  "very strong intensity-dependent missing (approaching left-censoring)"
                }
                log_info("  Interpretation: {slope_interpretation}")
              }
              
              # Perform row-wise imputation using limpa
              if (verbose) {
                log_info("Performing row-wise imputation using DPC model...")
              }
              
              y_imputed <- limpa::dpcImpute(y_peptide, dpc = dpcfit)
              
              if (verbose) {
                log_info("Imputation completed successfully")
                log_info("No missing values remaining: {!any(is.na(y_imputed$E))}")
              }
              
              # Extract the imputed matrix
              imputed_matrix <- y_imputed$E
              
              # Transform back to original scale if necessary
              if (use_log2_transform && !theObject@is_logged_data) {
                if (verbose) {
                  log_info("Converting back from log2 scale...")
                }
                imputed_matrix <- 2^imputed_matrix - 1
                # Ensure no negative values
                imputed_matrix[imputed_matrix < 0] <- 0
              }
              
              # Convert back to long format and merge with original data
              if (verbose) {
                log_info("Converting imputed data back to original format...")
              }
              
              # Create peptide IDs that match the matrix rownames
              peptide_ids <- rownames(imputed_matrix)
              
              # Convert imputed matrix to long format
              imputed_long <- imputed_matrix |>
                as.data.frame() |>
                tibble::rownames_to_column("peptide_id") |>
                tidyr::pivot_longer(cols = -peptide_id, 
                                   names_to = sample_id_column, 
                                   values_to = imputed_value_column) |>
                tidyr::separate(peptide_id, 
                               into = c(theObject@protein_id_column, theObject@peptide_sequence_column), 
                               sep = "%")
              
              # Merge with original peptide data
              updated_peptide_data <- peptide_data |>
                dplyr::left_join(imputed_long, 
                                by = c(theObject@protein_id_column, 
                                      theObject@peptide_sequence_column,
                                      sample_id_column))
              
              # Update the object
              theObject@peptide_data <- updated_peptide_data
              theObject@peptide_matrix <- imputed_matrix
              
              # Update norm_quantity_column to point to the new imputed column
              # This ensures plotting functions use the final imputed data
              theObject@norm_quantity_column <- imputed_value_column
              
              # Store DPC results in the object for future reference
              if (is.null(theObject@args)) {
                theObject@args <- list()
              }
              theObject@args$limpa_dpc_results <- list(
                dpc_parameters = dpcfit$dpc,  # Numeric vector c(intercept, slope)
                dpc_object = dpcfit,          # Full DPC object (preferred for dpcQuant)
                missing_percentage_before = round(100 * mean(is.na(y_peptide)), 1),
                missing_percentage_after = round(100 * mean(is.na(imputed_matrix)), 1),
                slope_interpretation = slope_interpretation,
                dpc_method = "limpa_dpc",
                # Store the original y_peptide data for recreating DPC plot
                y_peptide_for_dpc = y_peptide
              )
              
              # Clean design matrix
              theObject <- cleanDesignMatrixPeptide(theObject)
              
              if (verbose) {
                log_info("limpa-based imputation completed successfully!")
                log_info("New imputed column: {imputed_value_column}")
                log_info("DPC parameters stored in object@args$limpa_dpc_results")
              }
              
              return(theObject)
              
            }, error = function(e) {
              log_error("Error during limpa imputation: {e$message}")
              stop(paste("limpa imputation failed:", e$message))
            })
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting Methods for PeptideQuantitativeData
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export
setMethod(f="plotRle"
          , signature="PeptideQuantitativeData"
          , definition=function(theObject, grouping_variable, yaxis_limit = c(), sample_label = NULL) {
            peptide_matrix <- theObject@peptide_matrix
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id

            design_matrix <- as.data.frame(design_matrix)

            if(!is.null(sample_label)) {
              if (sample_label %in% colnames(design_matrix)) {
                rownames(design_matrix) <- design_matrix[,sample_label]
                colnames(peptide_matrix) <- design_matrix[,sample_label]
              } 
            } else {
              rownames(design_matrix) <- design_matrix[,sample_id]
            }

            rowinfo_vector <- NA
            if(!is.na(grouping_variable)){
              rowinfo_vector <- design_matrix[colnames(peptide_matrix), grouping_variable]
            }

            rle_plot <- plotRleHelper(t(peptide_matrix)
                                     , rowinfo = rowinfo_vector
                                     , yaxis_limit = yaxis_limit)

            return(rle_plot)
          })

plotPca <- function(theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size=8, cv_percentile = 0.90) {
  # Defensive checks
  if (!is.character(grouping_variable) || length(grouping_variable) != 1) {
    stop("grouping_variable must be a single character string")
  }
  
  if (!is.null(shape_variable) && (!is.character(shape_variable) || length(shape_variable) != 1)) {
    stop("shape_variable must be NULL or a single character string")
  }
  
  if (!grouping_variable %in% colnames(theObject@design_matrix)) {
    stop(sprintf("grouping_variable '%s' not found in design matrix", grouping_variable))
  }
  
  if (!is.null(shape_variable) && !shape_variable %in% colnames(theObject@design_matrix)) {
    stop(sprintf("shape_variable '%s' not found in design matrix", shape_variable))
  }
  
  data_matrix <- NULL
  ## I want to check the class of theObject here
 if( class(theObject) == "PeptideQuantitativeData") {
   data_matrix <- theObject@peptide_matrix
   
 } else if( class(theObject) == "ProteinQuantitativeData") {
   data_matrix <- theObject@protein_quant_table |>
     column_to_rownames(var = "Protein.Ids") |>
     as.matrix()
 }
  
  design_matrix <- theObject@design_matrix
  sample_id <- theObject@sample_id
  
  # Prepare matrix for PCA (data should already be log2 transformed)
  data_matrix_pca <- data_matrix
  data_matrix_pca[!is.finite(data_matrix_pca)] <- NA
  
  if(is.na(label_column) || label_column == "") {
    label_column <- ""
  }
  
  pca_plot <- plotPcaHelper(data_matrix_pca
                            , design_matrix = design_matrix
                            , sample_id_column = sample_id
                            , grouping_variable = grouping_variable
                            , shape_variable = shape_variable
                            , label_column = label_column
                            , title = title
                            , geom.text.size = font_size
                            , cv_percentile = cv_percentile)
  
  return(pca_plot)
}

# ##'@export
# #setMethod(f="plotPca"
#           , signature="PeptideQuantitativeData"
#           , definition=function(theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size=8) {
#             # Defensive checks
#             if (!is.character(grouping_variable) || length(grouping_variable) != 1) {
#               stop("grouping_variable must be a single character string")
#             }
            
#             if (!is.null(shape_variable) && (!is.character(shape_variable) || length(shape_variable) != 1)) {
#               stop("shape_variable must be NULL or a single character string")
#             }
            
#             if (!grouping_variable %in% colnames(theObject@design_matrix)) {
#               stop(sprintf("grouping_variable '%s' not found in design matrix", grouping_variable))
#             }
            
#             if (!is.null(shape_variable) && !shape_variable %in% colnames(theObject@design_matrix)) {
#               stop(sprintf("shape_variable '%s' not found in design matrix", shape_variable))
#             }

#             peptide_matrix <- theObject@peptide_matrix
#             design_matrix <- theObject@design_matrix
#             sample_id <- theObject@sample_id

#             # Prepare matrix for PCA (data should already be log2 transformed)
#             peptide_matrix_pca <- peptide_matrix
#             peptide_matrix_pca[!is.finite(peptide_matrix_pca)] <- NA

#             if(is.na(label_column) || label_column == "") {
#               label_column <- ""
#             }

#             pca_plot <- plotPcaHelper(peptide_matrix_pca
#                                      , design_matrix = design_matrix
#                                      , sample_id_column = sample_id
#                                      , grouping_variable = grouping_variable
#                                      , shape_variable = shape_variable
#                                      , label_column = label_column
#                                      , title = title
#                                      , geom.text.size = font_size)

#             return(pca_plot)
#           })

# Set up S4 class definitions for ggplot objects
setOldClass(c("gg", "ggplot"))
setOldClass("ggplot2::ggplot")

#' Create Density Plots from PCA data
#'
#' This function takes a ggplot object (presumably a PCA plot) and creates
#' density plots for the first two principal components.

setMethod(f="plotDensity"
          , signature="ggplot2::ggplot"
          , definition=function(theObject, grouping_variable, title = "", font_size = 8) {
            # First try to get data directly from the ggplot object's data element
            if (!is.null(theObject$data) && is.data.frame(theObject$data)) {
              pca_data <- as_tibble(theObject$data)
            } else {
              # Fall back to other extraction methods
              pca_data <- as_tibble(ggplot_build(theObject)$data[[1]])
              
              # If the data doesn't have PC1/PC2, try to extract from the plot's environment
              if (!("PC1" %in% colnames(pca_data) && "PC2" %in% colnames(pca_data))) {
                # Try to get the data from the plot's environment
                if (exists("data", envir = environment(theObject$mapping$x))) {
                  pca_data <- as_tibble(get("data", envir = environment(theObject$mapping$x)))
                } else {
                  stop("Could not extract PCA data from the ggplot object")
                }
              }
            }
            
            # Check if grouping variable exists in the data
            if (!grouping_variable %in% colnames(pca_data)) {
              stop(sprintf("grouping_variable '%s' not found in the data", grouping_variable))
            }

            # Create PC1 density plot
            pc1_density <- ggplot(pca_data, aes(x = PC1, fill = !!sym(grouping_variable), color = !!sym(grouping_variable))) +
              geom_density(alpha = 0.3) +
              theme_bw() +
              labs(title = title,
                   x = "PC1",
                   y = "Density") +
              theme(
                text = element_text(size = font_size),
                plot.margin = margin(b = 0, t = 5, l = 5, r = 5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
              )

            # Create PC2 density plot
            pc2_density <- ggplot(pca_data, aes(x = PC2, fill = !!sym(grouping_variable), color = !!sym(grouping_variable))) +
              geom_density(alpha = 0.3) +
              theme_bw() +
              labs(x = "PC2",
                   y = "Density") +
              theme(
                text = element_text(size = font_size),
                plot.margin = margin(t = 0, b = 5, l = 5, r = 5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
              )
            
            # Combine plots with minimal spacing
            combined_plot <- pc1_density / pc2_density +
              plot_layout(heights = c(1, 1)) +
              plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))
            
            return(combined_plot)
          }) 

#'@export
setMethod(f="plotPearson"
          , signature="PeptideQuantitativeData"
          , definition=function(theObject, tech_rep_remove_regex, correlation_group = NA) {
            peptide_matrix <- theObject@peptide_matrix
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id
            
            correlation_group_to_use <- correlation_group
            
            if(is.na(correlation_group)) {
              correlation_group_to_use <- theObject@technical_replicate_id
            }

            # Create a temporary ProteinQuantitativeData-like structure to use existing correlation function
            # We'll adapt the pearsonCorForSamplePairs function logic for peptides
            
            # Convert peptide matrix to data frame format similar to protein_quant_table
            peptide_data_for_corr <- peptide_matrix |>
              as.data.frame() |>
              rownames_to_column("peptide_id")
            
            # Create a temporary object structure for correlation calculation
            temp_obj <- list(
              data_table = peptide_data_for_corr,
              id_column = "peptide_id",
              design_matrix = design_matrix,
              sample_id = sample_id
            )
            
            # Calculate correlations between sample pairs
            correlation_vec <- calculatePeptidePearsonCorrelation(temp_obj, 
                                                                 tech_rep_remove_regex,
                                                                 correlation_group_to_use)

            pearson_plot <- correlation_vec |>
              ggplot(aes(pearson_correlation)) +
              geom_histogram(breaks = seq(min(round(correlation_vec$pearson_correlation - 0.5, 2), na.rm = TRUE), 1, 0.001)) +
              scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4)) +
              xlab("Pearson Correlation") +
              ylab("Counts") +
              theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank())

            return(pearson_plot)
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Helper function for peptide Pearson correlation calculation
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

calculatePeptidePearsonCorrelation <- function(temp_obj, tech_rep_remove_regex, correlation_group) {
  data_table <- temp_obj$data_table
  id_column <- temp_obj$id_column
  design_matrix <- temp_obj$design_matrix
  sample_id <- temp_obj$sample_id
  
  # Get sample columns (exclude ID column)
  sample_columns <- setdiff(colnames(data_table), id_column)
  
  # Filter out technical replicates if regex provided
  if(!is.null(tech_rep_remove_regex) && tech_rep_remove_regex != "") {
    sample_columns <- sample_columns[!grepl(tech_rep_remove_regex, sample_columns)]
  }
  
  # Create correlation matrix
  peptide_matrix_for_corr <- data_table |>
    column_to_rownames(id_column) |>
    select(all_of(sample_columns)) |>
    as.matrix()
  
  # Calculate correlations between all sample pairs
  sample_correlations <- cor(peptide_matrix_for_corr, use = "pairwise.complete.obs")
  
  # Extract upper triangle (avoid duplicate pairs and self-correlations)
  upper_tri_indices <- which(upper.tri(sample_correlations), arr.ind = TRUE)
  
  correlation_results <- data.frame(
    sample1 = rownames(sample_correlations)[upper_tri_indices[,1]],
    sample2 = colnames(sample_correlations)[upper_tri_indices[,2]],
    pearson_correlation = sample_correlations[upper_tri_indices],
    stringsAsFactors = FALSE
  )
  
  # Remove NA correlations
  correlation_results <- correlation_results[!is.na(correlation_results$pearson_correlation), ]
  
  return(correlation_results)
}

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Normalization Methods for PeptideQuantitativeData
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

## normalise between Arrays
#'@export
#'@param theObject Object of class PeptideQuantitativeData
#'@param normalisation_method Method to use for normalisation. Options are cyclicloess, quantile, scale, none
#' @title Normalize Between Samples for PeptideQuantitativeData
#' @name normaliseBetweenSamples,PeptideQuantitativeData-method
#' @export
setMethod(f="normaliseBetweenSamples"
          , signature="PeptideQuantitativeData"
          , definition=function(theObject, normalisation_method = NULL) {
            peptide_data <- theObject@peptide_data
            design_matrix <- theObject@design_matrix
            sample_id <- theObject@sample_id
            current_quant_column <- theObject@norm_quantity_column

            normalisation_method <- checkParamsObjectFunctionSimplify(theObject
                                                                      , "normalisation_method"
                                                                      , "cyclicloess")

            theObject <- updateParamInObject(theObject, "normalisation_method")

            # Create matrix from peptide_data (like protein version from protein_quant_table)
            # Create unique peptide identifier for matrix rownames
            peptide_unique_id <- paste(peptide_data[[theObject@protein_id_column]], 
                                      peptide_data[[theObject@peptide_sequence_column]], 
                                      sep = "%")
            
            # Create wide format data frame then convert to matrix (exactly like protein version)
            temp_peptide_wide <- peptide_data |>
              mutate(peptide_id = peptide_unique_id) |>
              select(peptide_id, !!sym(sample_id), !!sym(current_quant_column)) |>
              pivot_wider(names_from = !!sym(sample_id), 
                         values_from = !!sym(current_quant_column),
                         values_fn = mean)

            # Convert to matrix exactly like protein version: column_to_rownames() then as.matrix()
            frozen_peptide_matrix <- temp_peptide_wide |>
              column_to_rownames("peptide_id") |>
              as.matrix()

            frozen_peptide_matrix[!is.finite(frozen_peptide_matrix)] <- NA

            normalised_frozen_peptide_matrix <- frozen_peptide_matrix

            print(paste0("normalisation_method = ", normalisation_method))

            switch(normalisation_method
                   , cyclicloess = {
                     normalised_frozen_peptide_matrix <- limma::normalizeCyclicLoess(frozen_peptide_matrix)
                   }
                   , quantile = {
                     normalised_frozen_peptide_matrix <- limma::normalizeQuantiles(frozen_peptide_matrix)
                   }
                   , scale = {
                     normalised_frozen_peptide_matrix <- limma::normalizeMedianAbsValues(frozen_peptide_matrix)
                   }
                   , none = {
                     normalised_frozen_peptide_matrix <- frozen_peptide_matrix
                   }
            )

            normalised_frozen_peptide_matrix[!is.finite(normalised_frozen_peptide_matrix)] <- NA

            # Convert back to data frame (exactly like protein version)
            normalised_peptide_table <- normalised_frozen_peptide_matrix |>
              as.data.frame() |>
              rownames_to_column("peptide_id")

            # Update peptide_data by joining with normalized values and replacing the quant column
            updated_peptide_data <- peptide_data |>
              mutate(peptide_id = peptide_unique_id) |>
              select(-!!sym(current_quant_column)) |>
              left_join(normalised_peptide_table |> 
                       pivot_longer(cols = -peptide_id, 
                                   names_to = sample_id, 
                                   values_to = current_quant_column),
                       by = c("peptide_id", sample_id)) |>
              select(-peptide_id)

            # Update both slots
            theObject@peptide_data <- updated_peptide_data
            theObject@peptide_matrix <- normalised_frozen_peptide_matrix

            theObject <- cleanDesignMatrixPeptide(theObject)

            return(theObject)
          })

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Log2 Transform Peptide Matrix
#'
#' Transforms raw peptide intensity values to log2 scale for downstream normalization and RUV analysis.
#' This should be called after calcPeptideMatrix() and before normaliseBetweenSamples().
#' @export
#'
#' @param theObject A PeptideQuantitativeData object
#' @return PeptideQuantitativeData object with log2 transformed data
#' @export
setMethod(f="log2TransformPeptideMatrix"
          , signature="PeptideQuantitativeData"
          , definition=function(theObject) {
            
            if (theObject@is_logged_data) {
              warning("Data appears to already be log-transformed (is_logged_data = TRUE). Skipping transformation.")
              return(theObject)
            }
            
            peptide_data <- theObject@peptide_data
            peptide_matrix <- theObject@peptide_matrix
            current_quant_column <- theObject@norm_quantity_column
            
            # Log2 transform the matrix
            log2_peptide_matrix <- peptide_matrix
            
            # Handle zeros and negative values
            log2_peptide_matrix[log2_peptide_matrix <= 0] <- NA
            log2_peptide_matrix <- log2(log2_peptide_matrix)
            
            # Update matrix
            theObject@peptide_matrix <- log2_peptide_matrix
            
            # Also update peptide_data to maintain consistency
            # Create long format from log2 matrix
            log2_long <- log2_peptide_matrix |>
              as.data.frame() |>
              rownames_to_column("peptide_row_id") |>
              pivot_longer(cols = -peptide_row_id, 
                          names_to = theObject@sample_id, 
                          values_to = "log2_value")

            # Update peptide_data: match by protein%peptide ID and sample
            updated_peptide_data <- peptide_data |>
              mutate(peptide_row_id = paste(!!sym(theObject@protein_id_column), 
                                           !!sym(theObject@peptide_sequence_column), 
                                           sep = "%")) |>
              left_join(log2_long, by = c("peptide_row_id", theObject@sample_id)) |>
              mutate(!!sym(current_quant_column) := ifelse(!is.na(log2_value), 
                                                           log2_value, 
                                                           !!sym(current_quant_column))) |>
              select(-peptide_row_id, -log2_value)

            theObject@peptide_data <- updated_peptide_data
            
            # Mark as logged
            theObject@is_logged_data <- TRUE
            
            theObject <- cleanDesignMatrixPeptide(theObject)
            
            message("Peptide data successfully log2 transformed. Raw intensities converted to log2 scale.")
            message(paste("is_logged_data flag set to:", theObject@is_logged_data))
            
            return(theObject)
          })


#'@export
setMethod(f = "chooseBestProteinAccession"
          , signature="PeptideQuantitativeData"
          , definition=function(theObject, delim=NULL, seqinr_obj=NULL
                              , seqinr_accession_column=NULL
                              , replace_zero_with_na = NULL # Note: this param is kept for signature consistency but not used at peptide level
                              , aggregation_method = NULL) {

            peptide_data <- theObject@peptide_data
            protein_id_column <- theObject@protein_id_column
            is_logged <- theObject@is_logged_data
            verbose <- TRUE # Assuming verbose is desired, can be made a parameter

            # --- Parameter Handling ---
            delim <- checkParamsObjectFunctionSimplify(theObject, "delim",  default_value =  ";")
            seqinr_obj <- checkParamsObjectFunctionSimplify(theObject, "seqinr_obj",  default_value = NULL)
            seqinr_accession_column <- checkParamsObjectFunctionSimplify(theObject
                                                                       , "seqinr_accession_column"
                                                                       , default_value = "uniprot_acc")
            aggregation_method <- checkParamsObjectFunctionSimplify(theObject
                                                                  , "aggregation_method"
                                                                  , default_value = "mean")

            if (!aggregation_method %in% c("sum", "mean", "median")) {
              stop("aggregation_method must be one of: 'sum', 'mean', 'median'")
            }

            theObject <- updateParamInObject(theObject, "delim")
            theObject <- updateParamInObject(theObject, "seqinr_obj")
            theObject <- updateParamInObject(theObject, "seqinr_accession_column")
            theObject <- updateParamInObject(theObject, "aggregation_method")
            
            if (verbose) {
              log_info("Choosing best protein accession at the peptide level...")
              log_info("Aggregation method for duplicate peptides will be: '{aggregation_method}'")
            }

            # --- Map Old to New IDs ---
            accession_mapping <- chooseBestProteinAccessionHelper(input_tbl = peptide_data,
                                                                  acc_detail_tab = seqinr_obj,
                                                                  accessions_column = !!sym(protein_id_column),
                                                                  row_id_column = seqinr_accession_column,
                                                                  group_id = !!sym(protein_id_column),
                                                                  delim = delim)

            # --- Update Protein IDs ---
            updated_peptide_data <- peptide_data |>
              dplyr::left_join(accession_mapping |> dplyr::select(!!sym(protein_id_column), !!sym(seqinr_accession_column)), 
                               by = setNames(protein_id_column, protein_id_column)) |>
              dplyr::select(-!!sym(protein_id_column)) |>
              dplyr::rename(!!sym(protein_id_column) := !!sym(seqinr_accession_column))

            # --- Aggregate Duplicates ---
            if (verbose) {
              log_info("Aggregating peptide quantities for newly resolved protein groups...")
            }

            grouping_cols <- c(theObject@sample_id, protein_id_column, theObject@peptide_sequence_column)
            quant_cols <- c(theObject@raw_quantity_column, theObject@norm_quantity_column)
            meta_cols <- setdiff(colnames(updated_peptide_data), c(grouping_cols, quant_cols))

            aggregated_data <- updated_peptide_data |>
              dplyr::group_by(across(all_of(grouping_cols))) |>
              dplyr::summarise(
                # Aggregate quantitative columns based on the chosen method
                across(all_of(quant_cols), function(x) {
                  # Remove NAs for calculation
                  x_clean <- x[!is.na(x)]
                  if (length(x_clean) == 0) return(NA_real_)
                  
                  # Perform aggregation
                  if (aggregation_method == "mean") {
                    mean(x_clean)
                  } else if (aggregation_method == "median") {
                    median(x_clean)
                  } else if (aggregation_method == "sum") {
                    if (is_logged) {
                      log2(sum(2^x_clean)) # Sum on linear scale, then log2 transform back
                    } else {
                      sum(x_clean)
                    }
                  }
                }),
                # Keep the first value for all metadata columns
                across(all_of(meta_cols), dplyr::first),
                .groups = "drop"
              )

            if (verbose) {
              log_info("Aggregation complete. Before: {nrow(updated_peptide_data)} rows. After: {nrow(aggregated_data)} rows.")
            }

            theObject@peptide_data <- aggregated_data

            # --- Recalculate Matrix ---
            if (verbose) {
              log_info("Recalculating peptide matrix with updated protein accessions...")
            }
            theObject <- calcPeptideMatrix(theObject)

            return(theObject)
          })
# ==========================================
# Content from proteinVsSamplesS4Objects.R
# ==========================================
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Protein-Level Missing Value Imputation using limpa Package
#'
#' This function applies limpa's DPC-based missing value imputation directly to 
#' protein-level quantification data. This is useful when you already have protein
#' quantification data with missing values that need to be handled.
#'
#' @param theObject A ProteinQuantitativeData object with protein-level data
#' @param dpc_results DPC results to use. If NULL, will estimate using dpc_slope
#' @param dpc_slope Default DPC slope to use if no DPC results available (default: 0.8)
#' @param quantified_protein_column Name for the column containing quantified protein values
#' @param verbose Whether to print progress messages. Default is TRUE
#' @param chunk When verbose=TRUE, how often to output progress information (default: 1000)
#'
#' @details
#' This method treats each protein as a separate "feature" and applies DPC-based
#' imputation using limpa's dpcImpute function. This is appropriate when you have
#' protein-level data with missing values that follow intensity-dependent patterns.
#'
#' @return Updated ProteinQuantitativeData object with imputed protein values
#'
#' @importFrom limpa dpc dpcImpute
#' @export
setMethod(f="proteinMissingValueImputationLimpa"
          , signature="ProteinQuantitativeData"
          , definition = function(theObject, 
                                  dpc_results = NULL,
                                  dpc_slope = 0.8,
                                  quantified_protein_column = NULL,
                                  verbose = TRUE,
                                  chunk = 1000) {
            
            # Load required packages
            if (!requireNamespace("limpa", quietly = TRUE)) {
              stop("limpa package is required but not installed. Please install it using: BiocManager::install('limpa')")
            }
            
            # Parameter validation and defaults
            quantified_protein_column <- if (is.null(quantified_protein_column)) {
              "Protein.Imputed.Limpa"
            } else {
              quantified_protein_column
            }
            
            # Extract data from protein object
            protein_quant_table <- theObject@protein_quant_table
            protein_id_column <- theObject@protein_id_column
            sample_id_column <- theObject@sample_id
            design_matrix <- theObject@design_matrix
            
            if (verbose) {
              log_info("Starting limpa-based protein-level missing value imputation...")
            }
            
            # Convert to matrix format (proteins as rows, samples as columns)
            # First, identify sample columns (exclude protein ID and other metadata)
            sample_columns <- setdiff(colnames(protein_quant_table), 
                                     c(protein_id_column, "description", "gene_name", 
                                       "protein_name", "organism", "length"))
            
            if (verbose) {
              log_info("Converting protein data to matrix format...")
              log_info("Found {length(sample_columns)} sample columns")
            }
            
            # Create protein matrix
            protein_matrix <- protein_quant_table |>
              dplyr::select(all_of(c(protein_id_column, sample_columns))) |>
              tibble::column_to_rownames(protein_id_column) |>
              as.matrix()

            if (verbose) {
              log_info("Protein matrix dimensions: {nrow(protein_matrix)} proteins x {ncol(protein_matrix)} samples")
              log_info("Missing value percentage: {round(100 * mean(is.na(protein_matrix)), 1)}%")
            }
            
            # Check if we need log2 transformation
            # Assume if max value > 50, data is not log2 transformed
            max_val <- max(protein_matrix, na.rm = TRUE)
            needs_log_transform <- max_val > 50
            
            if (needs_log_transform) {
              if (verbose) {
                log_info("Converting to log2 scale for limpa (max value: {round(max_val, 2)})...")
              }
              protein_matrix <- log2(protein_matrix + 1)
            } else {
              if (verbose) {
                log_info("Data appears to be log2-scale already (max value: {round(max_val, 2)})")
              }
            }
            
            # Handle infinite or NaN values
            if (any(is.infinite(protein_matrix) | is.nan(protein_matrix), na.rm = TRUE)) {
              if (verbose) {
                log_warn("Infinite or NaN values detected. Replacing with NA...")
              }
              protein_matrix[is.infinite(protein_matrix) | is.nan(protein_matrix)] <- NA
            }
            
            # Get or estimate DPC parameters
            dpc_params <- NULL
            if (!is.null(dpc_results)) {
              dpc_params <- dpc_results
              if (verbose) {
                log_info("Using provided DPC results")
              }
            } else {
              if (verbose) {
                log_info("Estimating DPC from protein data...")
              }
            tryCatch({
                dpcfit <- limpa::dpc(protein_matrix)
                dpc_params <- dpcfit
                if (verbose) {
                  log_info("DPC parameters estimated:")
                  log_info("  beta0 (intercept): {round(dpcfit$dpc[1], 4)}")
                  log_info("  beta1 (slope): {round(dpcfit$dpc[2], 4)}")
                }
            }, error = function(e) {
                if (verbose) {
                  log_warn("DPC estimation failed, using default slope: {dpc_slope}")
                }
                dpc_params <- NULL
              })
            }
            
            # Apply missing value imputation
            if (verbose) {
              log_info("Applying protein-level missing value imputation...")
            }
            
              tryCatch({
              # Apply dpcImpute to protein matrix
              if (!is.null(dpc_params)) {
                imputed_result <- limpa::dpcImpute(protein_matrix, dpc = dpc_params, verbose = verbose, chunk = chunk)
            } else {
                imputed_result <- limpa::dpcImpute(protein_matrix, dpc.slope = dpc_slope, verbose = verbose, chunk = chunk)
              }
              
              if (verbose) {
                log_info("Protein-level imputation completed successfully")
                log_info("No missing values remaining: {!any(is.na(imputed_result$E))}")
              }
              
              # Extract imputed matrix
              imputed_matrix <- imputed_result$E
              
              # Transform back to original scale if necessary
              if (needs_log_transform) {
                if (verbose) {
                  log_info("Converting back from log2 scale...")
                }
                imputed_matrix <- 2^imputed_matrix - 1
                # Ensure no negative values
                imputed_matrix[imputed_matrix < 0] <- 0
              }
              
              # Convert back to long format
              if (verbose) {
                log_info("Converting imputed data back to original format...")
              }
              
              imputed_long <- imputed_matrix |>
                      as.data.frame() |>
                tibble::rownames_to_column(protein_id_column) |>
                tidyr::pivot_longer(cols = -all_of(protein_id_column), 
                                   names_to = sample_id_column, 
                                   values_to = quantified_protein_column)
              
              # Merge with original protein data
              updated_protein_data <- protein_quant_table |>
                dplyr::left_join(imputed_long, by = c(protein_id_column, sample_id_column))
              
              # Update the object
              theObject@protein_quant_table <- updated_protein_data
              
              # Store DPC results for future reference
              if (is.null(theObject@args)) {
                theObject@args <- list()
              }
              
              theObject@args$limpa_protein_imputation_results <- list(
                dpc_parameters_used = if (!is.null(dpc_params)) {
                  if (is.list(dpc_params) && !is.null(dpc_params$dpc)) {
                    dpc_params$dpc  # Extract parameters from DPC object
                  } else if (is.numeric(dpc_params)) {
                    dpc_params  # Already numeric parameters  
             } else {
                    c(NA, dpc_slope)
                  }
                } else {
                  c(NA, dpc_slope)
                },
                dpc_object_used = if (is.list(dpc_params) && !is.null(dpc_params$dpc)) dpc_params else NULL,
                quantified_protein_column = quantified_protein_column,
                missing_percentage_before = round(100 * mean(is.na(protein_matrix)), 1),
                missing_percentage_after = 0,  # DPC imputation produces complete data
                imputation_method = "limpa_dpc_protein_imputation",
                total_proteins_imputed = nrow(imputed_matrix)
              )
              
              if (verbose) {
                log_info("limpa protein-level imputation completed successfully!")
                log_info("New imputed column: {quantified_protein_column}")
                log_info("DPC results stored in object@args$limpa_protein_imputation_results")
              }
            
            return(theObject)
              
              }, error = function(e) {
              log_error(paste("Error during limpa protein imputation:", e$message))
              stop(paste("limpa protein imputation failed:", e$message))
            })
          })

# ==========================================
# Content from protein_de_analysis_wrapper.R
# ==========================================
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Differential Expression Analysis for Proteomics
#' 
#' This file contains modular functions for proteomics differential expression analysis,
#' breaking up the monolithic deAnalysisWrapperFunction into focused components.
#' 
#' Functions:
#' - differentialExpressionAnalysis: Main S4 generic for DE analysis
#' - differentialExpressionAnalysisHelper: Core limma-based analysis
#' - generateVolcanoPlotGlimma: Interactive volcano plot generation
#' - generateDEHeatmap: Heatmap visualization with clustering
#' 
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#' @export
setMethod( f ="differentialExpressionAnalysis"
           , signature = "ProteinQuantitativeData"
           , definition=function( theObject
                                  , contrasts_tbl = NULL
                                  , formula_string = NULL
                                  , group_id = NULL
                                  , de_q_val_thresh = NULL
                                  , treat_lfc_cutoff = NULL
                                  , eBayes_trend = NULL
                                  , eBayes_robust = NULL
                                  , args_group_pattern = NULL
                                  , args_row_id = NULL
                                  , qvalue_column = "fdr_qvalue"
                                  , raw_pvalue_column = "raw_pvalue" ) {

  # IMMEDIATE ERROR CATCH - Check if we even get here
  message("*** ENTERING differentialExpressionAnalysis METHOD ***")
  message(sprintf("*** METHOD SIGNATURE MATCHED: ProteinQuantitativeData ***"))
  
  # Try to catch the index error immediately
  tryCatch({
    message("*** Testing parameter access ***")
    if (!is.null(contrasts_tbl)) {
      test <- contrasts_tbl[[1]]
      message("*** Parameter access successful ***")
    }
  }, error = function(e) {
    message(sprintf("*** IMMEDIATE ERROR: %s ***", e$message))
    message(sprintf("*** ERROR CLASS: %s ***", class(e)))
    stop(e)
  })

  message("--- Entering differentialExpressionAnalysis ---")
  message(sprintf("   differentialExpressionAnalysis: theObject class = %s", class(theObject)))
  message(sprintf("   differentialExpressionAnalysis: contrasts_tbl provided = %s", !is.null(contrasts_tbl)))
  if(!is.null(contrasts_tbl)) {
    message(sprintf("   differentialExpressionAnalysis: contrasts_tbl dims = %d x %d", nrow(contrasts_tbl), ncol(contrasts_tbl)))
    message(sprintf("   differentialExpressionAnalysis: contrasts_tbl content = %s", paste(contrasts_tbl[[1]], collapse=", ")))
  }

  # Wrap the helper function call in tryCatch to get better error info
  message("   differentialExpressionAnalysis: About to call differentialExpressionAnalysisHelper...")
  
  results_list <- tryCatch({
    differentialExpressionAnalysisHelper(  theObject
                                          , contrasts_tbl = contrasts_tbl
                                          , formula_string = formula_string
                                          , group_id = group_id
                                          , de_q_val_thresh = de_q_val_thresh
                                          , treat_lfc_cutoff = treat_lfc_cutoff
                                          , eBayes_trend = eBayes_trend
                                          , eBayes_robust = eBayes_robust
                                          , args_group_pattern = args_group_pattern
                                          , args_row_id = args_row_id
                                          , qvalue_column = qvalue_column
                                          , raw_pvalue_column = raw_pvalue_column
                                          )
  }, error = function(e) {
    # CRITICAL FIX: Use paste() for logger calls in error handlers to avoid interpolation bug
    message(paste("   differentialExpressionAnalysis ERROR in helper function:", e$message))
    message(paste("   differentialExpressionAnalysis ERROR call stack:", capture.output(traceback())))
    stop(e)
  })

  message("   differentialExpressionAnalysis: Helper function completed successfully!")
  message("--- Exiting differentialExpressionAnalysis ---")
  return(results_list)

})

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @export
setMethod( f ="differentialExpressionAnalysisHelper"
           , signature = "ProteinQuantitativeData"
           , definition=function( theObject
                                  , contrasts_tbl = NULL
                                  , formula_string = NULL
                                  , group_id = NULL
                                  , de_q_val_thresh = NULL
                                  , treat_lfc_cutoff = NULL
                                  , eBayes_trend = NULL
                                  , eBayes_robust = NULL
                                  , args_group_pattern = NULL
                                  , args_row_id = NULL
                                  , qvalue_column = "fdr_qvalue"
                                  , raw_pvalue_column = "raw_pvalue" ) {

  message("--- Entering differentialExpressionAnalysisHelper ---")
  
  # IMMEDIATE PARAMETER VALIDATION TO CATCH INDEX ERROR
  message("   PARAMETER VALIDATION Step: Checking all input parameters...")
  message(sprintf("      theObject is NULL: %s", is.null(theObject)))
  message(sprintf("      contrasts_tbl is NULL: %s", is.null(contrasts_tbl)))
  if (!is.null(contrasts_tbl)) {
    message(sprintf("      contrasts_tbl class: %s", class(contrasts_tbl)))
    message(sprintf("      contrasts_tbl type: %s", typeof(contrasts_tbl)))
    message("      contrasts_tbl print:")
    print(contrasts_tbl)
    message("      Trying to access contrasts_tbl[[1]]...")
    tryCatch({
      first_col <- contrasts_tbl[[1]]
      message(sprintf("      SUCCESS: contrasts_tbl[[1]] class: %s", class(first_col)))
      message(sprintf("      SUCCESS: contrasts_tbl[[1]] length: %d", length(first_col)))
      message(sprintf("      SUCCESS: contrasts_tbl[[1]] content: %s", paste(first_col, collapse=", ")))
    }, error = function(e) {
      message(sprintf("      ERROR accessing contrasts_tbl[[1]]: %s", e$message))
    })
  }
  
  message("   DEBUG66: Initial parameter inspection")
  message(sprintf("      DEBUG66: theObject class = %s", class(theObject)))
  message(sprintf("      DEBUG66: contrasts_tbl param is.null = %s", is.null(contrasts_tbl)))
  if (!is.null(contrasts_tbl)) {
    message(sprintf("      DEBUG66: contrasts_tbl param class = %s", class(contrasts_tbl)))
    message("      DEBUG66: contrasts_tbl param structure:")
    str(contrasts_tbl)
  }

  message("   differentialExpressionAnalysisHelper Step: Extracting parameters...")
  # Extract parameters from S4 object with fallbacks
  contrasts_tbl <- checkParamsObjectFunctionSimplify( theObject, "contrasts_tbl", contrasts_tbl)
  formula_string <- checkParamsObjectFunctionSimplify( theObject, "formula_string", "~ 0 + group")
  group_id <- checkParamsObjectFunctionSimplify( theObject, "group_id", "group")
  de_q_val_thresh <- checkParamsObjectFunctionSimplify( theObject, "de_q_val_thresh", 0.05)
  treat_lfc_cutoff <- checkParamsObjectFunctionSimplify( theObject, "treat_lfc_cutoff", 0)
  eBayes_trend <- checkParamsObjectFunctionSimplify( theObject, "eBayes_trend", TRUE)
  eBayes_robust <- checkParamsObjectFunctionSimplify( theObject, "eBayes_robust", TRUE)
  args_group_pattern <- checkParamsObjectFunctionSimplify( theObject, "args_group_pattern", "(\\d+)")
  args_row_id <- checkParamsObjectFunctionSimplify( theObject, "args_row_id", "uniprot_acc")
  
  message("   DEBUG66: After parameter extraction")
  message(sprintf("      DEBUG66: contrasts_tbl class = %s", class(contrasts_tbl)))
  message("      DEBUG66: contrasts_tbl structure:")
  str(contrasts_tbl)

  message(sprintf("   differentialExpressionAnalysisHelper: formula_string = %s", formula_string))
  message(sprintf("   differentialExpressionAnalysisHelper: group_id = %s", group_id))
  message(sprintf("   differentialExpressionAnalysisHelper: de_q_val_thresh = %f", de_q_val_thresh))

  message("   differentialExpressionAnalysisHelper Step: Processing group names...")
  # Handle group names that start with numbers (same pattern as original wrapper)
  design_matrix <- theObject@design_matrix
  group_col <- design_matrix[[group_id]]
  message(sprintf("   differentialExpressionAnalysisHelper: group_col length = %d", length(group_col)))
  message(sprintf("   differentialExpressionAnalysisHelper: group_col content = %s", paste(head(group_col, 10), collapse=", ")))
  
  # Check if any group names start with numbers and create mapping
  starts_with_number <- grepl("^[0-9]", group_col)
  message(sprintf("   differentialExpressionAnalysisHelper: any start with number? %s", any(starts_with_number)))
  
  if(any(starts_with_number)) {
    message("   differentialExpressionAnalysisHelper Step: Fixing group names that start with numbers...")
    
    original_groups <- unique(group_col)
    message(sprintf("   differentialExpressionAnalysisHelper: original_groups = %s", paste(original_groups, collapse=", ")))
    message(sprintf("   differentialExpressionAnalysisHelper: About to process %d original groups with purrr::map_chr", length(original_groups)))
    
    tryCatch({
      message("      DEBUG66: Before safe_groups purrr::map_chr")
      message(sprintf("      DEBUG66: original_groups class = %s", class(original_groups)))
      message(sprintf("      DEBUG66: original_groups length = %d", length(original_groups)))
      message("      DEBUG66: original_groups content:")
      print(original_groups)
      
      safe_groups <- purrr::map_chr(original_groups, \(x) {
        message(sprintf("         DEBUG66: Processing group item: '%s' (class: %s)", x, class(x)))
        result <- if(grepl("^[0-9]", x)) paste0("grp_", x) else x
        message(sprintf("         DEBUG66: Result for '%s' -> '%s'", x, result))
        result
      })
      message("   differentialExpressionAnalysisHelper: safe_groups processing SUCCESS")
      message("      DEBUG66: safe_groups result:")
      print(safe_groups)
    }, error = function(e) {
      message(sprintf("   differentialExpressionAnalysisHelper ERROR in safe_groups purrr::map_chr: %s", e$message))
      message("      DEBUG66: Error details:")
      message(sprintf("      DEBUG66: Error class: %s", class(e)))
      print(str(e))
      stop(e)
    })
    
    group_mapping <- setNames(original_groups, safe_groups)
    message("   differentialExpressionAnalysisHelper: group_mapping created")
    
    # Update design matrix with safe names
    message("   differentialExpressionAnalysisHelper: About to update design matrix with purrr::map_chr")
    tryCatch({
      design_matrix[[group_id]] <- purrr::map_chr(group_col, \(x) {
        if(grepl("^[0-9]", x)) paste0("grp_", x) else x
      })
      message("   differentialExpressionAnalysisHelper: design matrix update SUCCESS")
    }, error = function(e) {
      message(sprintf("   differentialExpressionAnalysisHelper ERROR in design matrix purrr::map_chr: %s", e$message))
      stop(e)
    })
    
    # Update contrasts table if it exists
    if(!is.null(contrasts_tbl)) {
      message("   differentialExpressionAnalysisHelper: About to update contrasts table with purrr::map_chr")
      message(sprintf("   differentialExpressionAnalysisHelper: contrasts_tbl[[1]] = %s", paste(contrasts_tbl[[1]], collapse=", ")))
      
      message("      DEBUG66: Contrasts table inspection before purrr::map_chr")
      message(sprintf("      DEBUG66: contrasts_tbl class = %s", class(contrasts_tbl)))
      message(sprintf("      DEBUG66: contrasts_tbl dim = %d x %d", nrow(contrasts_tbl), ncol(contrasts_tbl)))
      message("      DEBUG66: contrasts_tbl structure:")
      str(contrasts_tbl)
      message(sprintf("      DEBUG66: contrasts_tbl[[1]] class = %s", class(contrasts_tbl[[1]])))
      message(sprintf("      DEBUG66: contrasts_tbl[[1]] length = %d", length(contrasts_tbl[[1]])))
      message("      DEBUG66: group_mapping content:")
      print(group_mapping)
      
      tryCatch({
        contrasts_tbl[[1]] <- purrr::map_chr(contrasts_tbl[[1]], \(x) {
          message(sprintf("         DEBUG66: Processing contrast: '%s' (class: %s)", x, class(x)))
          message(sprintf("         DEBUG66: group_mapping names: %s", paste(names(group_mapping), collapse=", ")))
          result <- x
          for(orig in names(group_mapping)) {
            message(sprintf("            DEBUG66: Checking gsub: '%s' -> '%s'", group_mapping[orig], orig))
            result <- gsub(group_mapping[orig], orig, result, fixed = TRUE)
          }
          message(sprintf("         DEBUG66: Final result: '%s' -> '%s'", x, result))
          result
        })
        message("   differentialExpressionAnalysisHelper: contrasts table update SUCCESS")
        message("      DEBUG66: Updated contrasts_tbl[[1]]:")
        print(contrasts_tbl[[1]])
      }, error = function(e) {
        message(sprintf("   differentialExpressionAnalysisHelper ERROR in contrasts table purrr::map_chr: %s", e$message))
        message("      DEBUG66: Error details:")
        message(sprintf("      DEBUG66: Error class: %s", class(e)))
        print(str(e))
        stop(e)
      })
    }
    
    theObject@design_matrix <- design_matrix
  }

  message("   differentialExpressionAnalysisHelper Step: Updating S4 object parameters...")
  # Update object with validated parameters
  theObject <- updateParamInObject(theObject, "contrasts_tbl")
  theObject <- updateParamInObject(theObject, "formula_string")
  theObject <- updateParamInObject(theObject, "group_id")
  theObject <- updateParamInObject(theObject, "de_q_val_thresh")
  theObject <- updateParamInObject(theObject, "treat_lfc_cutoff")
  theObject <- updateParamInObject(theObject, "eBayes_trend")
  theObject <- updateParamInObject(theObject, "eBayes_robust")
  theObject <- updateParamInObject(theObject, "args_group_pattern")
  theObject <- updateParamInObject(theObject, "args_row_id")

  return_list <- list()
  return_list$theObject <- theObject

  message("   differentialExpressionAnalysisHelper Step: Generating QC plots...")

  # Generate QC plots (RLE, PCA, density)
  rle_plot <- plotRle(theObject = theObject, theObject@group_id) +
    theme(axis.text.x = element_text(size = 13)) +
    theme(axis.text.y = element_text(size = 13)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12)) +
    xlab("Samples")

  return_list$rle_plot <- rle_plot

  pca_plot <- plotPca( theObject
                      , grouping_variable = theObject@group_id
                      , label_column = ""
                      , title = ""
                      , font_size = 8) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))

  return_list$pca_plot <- pca_plot

  pca_plot_with_labels <- plotPca( theObject
                                  , grouping_variable = theObject@group_id
                                  , label_column = theObject@sample_id
                                  , title = ""
                                  , font_size = 8) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(axis.title.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 12))

  return_list$pca_plot_with_labels <- pca_plot_with_labels

  # Count the number of values
  return_list$plot_num_of_values <- plotNumOfValuesNoLog(theObject@protein_quant_table)

  message("   differentialExpressionAnalysisHelper Step: Running limma contrasts analysis...")

  # Prepare design matrix for limma
  message("   DEBUG66: About to prepare design matrix rownames")
  message(paste("      DEBUG66: theObject@sample_id =", theObject@sample_id))
  message(paste("      DEBUG66: design_matrix dims before =", nrow(theObject@design_matrix), "x", ncol(theObject@design_matrix)))
  message("      DEBUG66: design_matrix column names:")
  print(colnames(theObject@design_matrix))
  
  rownames( theObject@design_matrix ) <- theObject@design_matrix |> dplyr::pull( one_of(theObject@sample_id ))
  
  message("      DEBUG66: design_matrix rownames set successfully")

  # Run the core limma analysis using existing function
  protein_quant_matrix <- as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column))
  contrast_strings_to_use <- contrasts_tbl[, 1]
  message(paste("   differentialExpressionAnalysisHelper: About to call runTestsContrasts with", length(contrast_strings_to_use), "contrasts"))
  message(paste("   differentialExpressionAnalysisHelper: protein_quant_matrix dims =", nrow(protein_quant_matrix), "x", ncol(protein_quant_matrix)))
  
  contrasts_results <- runTestsContrasts(protein_quant_matrix,
                                         contrast_strings = contrast_strings_to_use,
                                         design_matrix = theObject@design_matrix,
                                         formula_string = formula_string,
                                         weights = NA,
                                         treat_lfc_cutoff = as.double(treat_lfc_cutoff),
                                         eBayes_trend = as.logical(eBayes_trend),
                                         eBayes_robust = as.logical(eBayes_robust))

  message("   differentialExpressionAnalysisHelper Step: runTestsContrasts completed successfully!")

  # The result from runTestsContrasts is a LIST with a 'results'
  # element, which ITSELF is a list of tables. For a single contrast run, we 
  # need to extract the first table from the nested list.
  
  # Check if the expected structure exists
  if (is.null(contrasts_results) || is.null(contrasts_results$results) || length(contrasts_results$results) == 0) {
      stop("Error: DE analysis function did not return results in the expected format.")
  }
  
  # Extract the results table (it's the first element in the nested list)
  de_results_table <- contrasts_results$results[[1]]
  
  # Get the name of the contrast from the list element name
  contrast_name <- names(contrasts_results$results)[1]
  
  # Ensure it's a data frame before proceeding
  if (!is.data.frame(de_results_table)) {
    stop("Error: DE analysis results are not in the expected data frame format.")
  }

  # CRITICAL FIX 3.0: The 'topTreat' table needs a 'comparison' column added to it,
  # containing the name of the contrast. It also needs the protein IDs from rownames.
  de_results_table <- de_results_table |>
    tibble::rownames_to_column(var = args_row_id) |>
    dplyr::mutate(comparison = contrast_name)

  # Map back to original group names in results if needed
  if(exists("group_mapping")) {
    contrasts_results_table <- de_results_table |>
      dplyr::mutate(comparison = purrr::map_chr(comparison, \(x) {
        result <- x
        for(safe_name in names(group_mapping)) {
          result <- gsub(safe_name, group_mapping[safe_name], result, fixed = TRUE)
        }
        result
      }))
  } else {
    contrasts_results_table <- de_results_table
  }

  return_list$contrasts_results <- contrasts_results
  return_list$contrasts_results_table <- contrasts_results_table
  
  # Extract qvalue warnings if present
  if (!is.null(contrasts_results$qvalue_warnings) && length(contrasts_results$qvalue_warnings) > 0) {
    return_list$qvalue_warnings <- contrasts_results$qvalue_warnings
    message(sprintf("   differentialExpressionAnalysisHelper Step: qvalue() failed for %d contrast(s)", length(contrasts_results$qvalue_warnings)))
  }

  message("   differentialExpressionAnalysisHelper Step: Preparing data for visualization...")
  message(paste("   DEBUG66: contrast_name for list naming =", contrast_name))

  # Prepare data for volcano plots
  # CRITICAL FIX: getSignificantData expects list_of_de_tables to be a list where each element 
  # is itself a list of data.frames. Since we're processing one contrast at a time, we need to 
  # wrap the single data.frame in an additional list layer.
  # CRITICAL FIX 2: The list element name MUST contain "=" delimiter because countStatDeGenesHelper
  # expects to split on "=" to separate comparison name from expression. Use contrast_name which
  # contains the full format like "H4_vs_WT=groupH4-groupWT"
  nested_list <- list(contrasts_results_table)
  names(nested_list) <- contrast_name  # Use actual contrast name with "=" delimiter
  significant_rows <- getSignificantData(list_of_de_tables = list(nested_list),
                                         list_of_descriptions = list("RUV applied"),
                                         row_id = !!sym(args_row_id),
                                         p_value_column = !!sym(raw_pvalue_column),
                                         q_value_column = !!sym(qvalue_column),
                                         fdr_value_column = fdr_value_bh_adjustment,
                                         log_q_value_column = lqm,
                                         log_fc_column = logFC,
                                         comparison_column = "comparison",
                                         expression_column = "log_intensity",
                                         facet_column = analysis_type,
                                         q_val_thresh = de_q_val_thresh) |>
    dplyr::rename(log2FC = "logFC")

  return_list$significant_rows <- significant_rows

  # Generate volcano plots
  volplot_plot <- plotVolcano(significant_rows,
                              log_q_value_column = lqm,
                              log_fc_column = log2FC,
                              q_val_thresh = de_q_val_thresh,
                              formula_string = "analysis_type ~ comparison")

  return_list$volplot_plot <- volplot_plot

  # Count significant molecules
  # CRITICAL FIX: Same as above - wrap in additional list layer and use contrast_name
  nested_list_for_count <- list(contrasts_results_table)
  names(nested_list_for_count) <- contrast_name  # Use actual contrast name with "=" delimiter
  num_sig_de_molecules <- printCountDeGenesTable(list_of_de_tables = list(nested_list_for_count),
                                                 list_of_descriptions = list("RUV applied"),
                                                 formula_string = "analysis_type ~ comparison")

  return_list$num_sig_de_molecules_first_go <- num_sig_de_molecules

  # P-values distribution plot
  pvalhist <- printPValuesDistribution(significant_rows,
                                       p_value_column = !!sym(raw_pvalue_column),
                                       formula_string = "analysis_type ~ comparison")

  return_list$pvalhist <- pvalhist

  message("   differentialExpressionAnalysisHelper Step: Creating results tables...")

  # Create wide format output
  counts_table_to_use <- theObject@protein_quant_table

  norm_counts <- counts_table_to_use |>
    as.data.frame() |>
    column_to_rownames(args_row_id) |>
    set_colnames(paste0(colnames(counts_table_to_use[-1]), ".log2norm")) |>
    rownames_to_column(args_row_id)

  return_list$norm_counts <- norm_counts

  de_proteins_wide <- significant_rows |>
    dplyr::filter(analysis_type == "RUV applied") |>
    dplyr::select(-lqm, -colour, -analysis_type) |>
    pivot_wider(id_cols = c(!!sym(args_row_id)),
                names_from = c(comparison),
                names_sep = ":",
                values_from = c(log2FC, !!sym(qvalue_column), !!sym(raw_pvalue_column))) |>
    left_join(counts_table_to_use, by = join_by( !!sym(args_row_id) == !!sym(theObject@protein_id_column))) |>
    left_join(theObject@protein_id_table, by = join_by( !!sym(args_row_id) == !!sym(theObject@protein_id_column))) |>
    dplyr::arrange(across(matches(paste0("!!sym(", qvalue_column, ")")))) |>
    distinct()

  return_list$de_proteins_wide <- de_proteins_wide

  # Create long format output
  de_proteins_long <- createDeResultsLongFormat( lfc_qval_tbl = significant_rows |>
                                                   dplyr::filter(analysis_type == "RUV applied"),
                                                 norm_counts_input_tbl = as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column)),
                                                 raw_counts_input_tbl = as.matrix(column_to_rownames(theObject@protein_quant_table, theObject@protein_id_column)),
                                                 row_id = args_row_id,
                                                 sample_id = theObject@sample_id,
                                                 group_id = group_id,
                                                 group_pattern = args_group_pattern,
                                                 design_matrix_norm = theObject@design_matrix,
                                                 design_matrix_raw = theObject@design_matrix,
                                                 protein_id_table = theObject@protein_id_table)

  return_list$de_proteins_long <- de_proteins_long

  # Static volcano plots with gene names
  static_volcano_plot_data <- de_proteins_long |>
    mutate( lqm = -log10(!!sym(qvalue_column))) |>
    dplyr::mutate(label = case_when( !!sym(qvalue_column) < de_q_val_thresh ~ "Significant",
                                     TRUE ~ "Not sig.")) |>
    dplyr::mutate(colour = case_when( !!sym(qvalue_column) < de_q_val_thresh ~ "purple",
                                      TRUE ~ "black")) |>
    dplyr::mutate(colour = factor(colour, levels = c("black", "purple")))

  list_of_volcano_plots <- static_volcano_plot_data %>%
    group_by( comparison) %>%
    nest() %>%
    ungroup() %>%
    mutate( title = paste( comparison)) %>%
    mutate( plot = purrr::map2( data, title, \(x,y) { plotOneVolcanoNoVerticalLines(x, y
                                                                                     , log_q_value_column = lqm
                                                                                     , log_fc_column = log2FC) } ) )

  return_list$list_of_volcano_plots <- list_of_volcano_plots

  # Additional summary statistics
  num_sig_de_molecules <- significant_rows %>%
    dplyr::mutate(status = case_when(!!sym(qvalue_column) >= de_q_val_thresh ~ "Not significant",
                                     log2FC >= 0 & !!sym(qvalue_column) < de_q_val_thresh ~ "Significant and Up",
                                     log2FC < 0 & !!sym(qvalue_column) < de_q_val_thresh ~ "Significant and Down",
                                     TRUE ~ "Not significant")) %>%
    group_by( comparison, status) %>%
    summarise(counts = n()) %>%
    ungroup()

  return_list$num_sig_de_molecules <- num_sig_de_molecules

  # Create barplots if significant results exist
  if (num_sig_de_molecules %>%
      dplyr::filter(status != "Not significant") |>
      nrow() > 0 ) {

    num_sig_de_genes_barplot_only_significant <- num_sig_de_molecules %>%
      dplyr::filter(status != "Not significant") %>%
      ggplot(aes(x = status, y = counts)) +
      geom_bar(stat = "identity") +
      geom_text(stat = 'identity', aes(label = counts), vjust = -0.5) +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(~ comparison)

    return_list$num_sig_de_genes_barplot_only_significant <- num_sig_de_genes_barplot_only_significant

    num_sig_de_genes_barplot_with_not_significant <- num_sig_de_molecules %>%
      ggplot(aes(x = status, y = counts)) +
      geom_bar(stat = "identity") +
      geom_text(stat = 'identity', aes(label = counts), vjust = -0.5) +
      theme(axis.text.x = element_text(angle = 90)) +
      facet_wrap(~ comparison)

    return_list$num_sig_de_genes_barplot_with_not_significant <- num_sig_de_genes_barplot_with_not_significant
  }

  message("--- Exiting differentialExpressionAnalysisHelper ---")
  return(return_list)

})

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Generate Interactive Volcano Plot using Glimma
#' 
#' @export
generateVolcanoPlotGlimma <- function( de_results_list
                                       , selected_contrast = NULL
                                       , uniprot_tbl = NULL
                                       , args_row_id = "uniprot_acc"
                                       , fdr_column = "fdr_qvalue"
                                       , raw_p_value_column = "raw_pvalue"
                                       , log2fc_column = "log2FC"
                                       , de_q_val_thresh = 0.05
                                       , uniprot_id_column = "Entry"
                                       , gene_names_column = "gene_names"
                                       , display_columns = c( "best_uniprot_acc" )) {

  message("--- Entering generateVolcanoPlotGlimma ---")
  message(paste("   generateVolcanoPlotGlimma Arg: selected_contrast =", selected_contrast))

  if (is.null(de_results_list) || is.null(de_results_list$de_proteins_long)) {
    message("   generateVolcanoPlotGlimma: No DE results available")
    return(NULL)
  }

  if (is.null(selected_contrast)) {
    message("   generateVolcanoPlotGlimma: No contrast selected")
    return(NULL)
  }

  # Get the contrast-specific results
  de_proteins_long <- de_results_list$de_proteins_long
  contrasts_results <- de_results_list$contrasts_results

  # CRITICAL FIX: Extract comparison name from selected_contrast if it contains "="
  # The selected_contrast might be the full name like "GA_Elevated.minus.GA_Control=groupGA_Elevated-groupGA_Control"
  # But the comparison column in de_proteins_long contains only "GA_Elevated.minus.GA_Control"
  comparison_to_search <- stringr::str_extract(selected_contrast, "^[^=]+")
  message(paste("   generateVolcanoPlotGlimma: Searching for comparison =", comparison_to_search))
  message(paste("   generateVolcanoPlotGlimma: Original selected_contrast =", selected_contrast))

  # Filter for selected contrast using the extracted comparison name
  contrast_data <- de_proteins_long |>
    dplyr::filter(comparison == comparison_to_search)

  if (nrow(contrast_data) == 0) {
    message(paste("   generateVolcanoPlotGlimma: No data found for contrast", comparison_to_search))
    # DEBUG: Show what comparisons are actually available
    available_comparisons <- unique(de_proteins_long$comparison)
    message(paste("   generateVolcanoPlotGlimma: Available comparisons:", paste(available_comparisons, collapse = ", ")))
    return(NULL)
  }

  message(paste("   generateVolcanoPlotGlimma: Found", nrow(contrast_data), "rows for contrast", comparison_to_search))
  message("   generateVolcanoPlotGlimma Step: Preparing volcano plot data...")

  # Prepare volcano plot table
  volcano_plot_tab <- contrast_data |>
    dplyr::mutate(best_uniprot_acc = str_split(!!sym(args_row_id), ":") |> purrr::map_chr(1))

  # Add UniProt annotations if available
  if (!is.null(uniprot_tbl)) {
    volcano_plot_tab <- volcano_plot_tab |>
      left_join(uniprot_tbl, by = join_by( best_uniprot_acc == !!sym( uniprot_id_column) )) |>
      dplyr::rename( UNIPROT_GENENAME = !!sym(gene_names_column) ) |>
      mutate( UNIPROT_GENENAME = purrr::map_chr( UNIPROT_GENENAME, \(x){
        tryCatch({
          if(is.na(x) || is.null(x) || x == "") {
            ""
          } else {
            split_result <- str_split(x, " |:")[[1]]
            if(length(split_result) > 0) split_result[1] else ""
          }
        }, error = function(e) "")
      })) |>
      dplyr::mutate(gene_name = UNIPROT_GENENAME )
  } else {
    volcano_plot_tab <- volcano_plot_tab |>
      dplyr::mutate(gene_name = best_uniprot_acc)
  }

  # Add visualization columns
  volcano_plot_tab <- volcano_plot_tab |>
    mutate( lqm = -log10(!!sym(fdr_column))) |>
    dplyr::mutate(label = case_when(abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= de_q_val_thresh ~ "Not sig., logFC >= 1",
                                    abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < de_q_val_thresh ~ "Sig., logFC >= 1",
                                    abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < de_q_val_thresh ~ "Sig., logFC < 1",
                                    TRUE ~ "Not sig.")) |>
    dplyr::mutate(colour = case_when(abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) >= de_q_val_thresh ~ "orange",
                                     abs(!!sym(log2fc_column)) >= 1 & !!sym(fdr_column) < de_q_val_thresh ~ "purple",
                                     abs(!!sym(log2fc_column)) < 1 & !!sym(fdr_column) < de_q_val_thresh ~ "blue",
                                     TRUE ~ "black")) |>
    dplyr::mutate(analysis_type = comparison) |>
    dplyr::select( best_uniprot_acc, lqm, !!sym(fdr_column), !!sym(raw_p_value_column), !!sym(log2fc_column), comparison, label, colour, gene_name
                   , any_of( display_columns )) |>
    dplyr::mutate( my_alpha = case_when ( gene_name != "" ~ 1
                                          , TRUE ~ 0.5))

  message("   generateVolcanoPlotGlimma Step: Generating interactive plot widget...")

  # Find the coefficient index for this contrast
  # CRITICAL FIX: Use grep for reliable pattern matching instead of stringr
  coef_names <- colnames(contrasts_results$fit.eb$coefficients)
  message(paste("   generateVolcanoPlotGlimma: Available coefficient names:", paste(coef_names, collapse = ", ")))
  message(paste("   generateVolcanoPlotGlimma: Looking for pattern:", paste0("^", comparison_to_search, "=")))
  
  # Use grep to find coefficient that starts with friendly name followed by "="
  coef_index <- grep(paste0("^", comparison_to_search, "="), coef_names)
  message(paste("   generateVolcanoPlotGlimma: Pattern match found indices:", paste(coef_index, collapse = ", ")))
  
  # Fallback: try exact match with the friendly name (for backwards compatibility)
  if (length(coef_index) == 0) {
    message("   generateVolcanoPlotGlimma: Pattern match failed, trying exact match with friendly name")
    coef_index <- which(coef_names == comparison_to_search)
  }
  
  # Last resort: try exact match with selected_contrast
  if (length(coef_index) == 0) {
    message("   generateVolcanoPlotGlimma: Exact friendly name match failed, trying selected_contrast")
    coef_index <- which(coef_names == selected_contrast)
  }

  if (length(coef_index) == 0) {
    message(paste("   generateVolcanoPlotGlimma: FINAL FAILURE - No coefficient found for:", selected_contrast))
    message(paste("   generateVolcanoPlotGlimma: Also tried:", comparison_to_search))
    message("   generateVolcanoPlotGlimma: Available coefficient patterns:")
    for (i in seq_along(coef_names)) {
      message(sprintf("     [%d] %s", i, coef_names[i]))
    }
    message(paste("   generateVolcanoPlotGlimma: Pattern attempted:", paste0("^", comparison_to_search, "=")))
    return(NULL)
  }

  message(paste("   generateVolcanoPlotGlimma: Using coefficient index", coef_index, "for contrast", coef_names[coef_index]))

  # Prepare counts matrix and groups
  counts_mat <- de_results_list$theObject@protein_quant_table |>
    column_to_rownames(de_results_list$theObject@protein_id_column) |>
    as.matrix()

  this_design_matrix <- de_results_list$theObject@design_matrix
  rownames( this_design_matrix ) <- this_design_matrix[, de_results_list$theObject@sample_id]
  this_groups <- this_design_matrix[colnames( counts_mat), "group"]

  # Generate the Glimma widget
  glimma_widget <- getGlimmaVolcanoProteomicsWidget( contrasts_results$fit.eb
                                                    , coef = coef_index
                                                    , volcano_plot_tab = volcano_plot_tab
                                                    , uniprot_column = best_uniprot_acc
                                                    , gene_name_column = gene_name
                                                    , display_columns = display_columns
                                                    , counts_tbl = counts_mat
                                                    , groups = this_groups )

  message("--- Exiting generateVolcanoPlotGlimma ---")
  return(glimma_widget)
}

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Generate DE Results Heatmap with Advanced Clustering
#' 
#' @export
generateDEHeatmap <- function( de_results_list
                              , selected_contrast = NULL
                              , top_n_genes = 50
                              , clustering_method = "ward.D2"
                              , distance_method = "euclidean"
                              , cluster_rows = TRUE
                              , cluster_cols = TRUE
                              , scale_data = "row"
                              , color_scheme = "RdBu"
                              , show_gene_names = FALSE
                              , de_q_val_thresh = 0.05
                              , qvalue_column = "fdr_qvalue"
                              , log2fc_column = "log2FC") {

  message("--- Entering generateDEHeatmap ---")
  message(sprintf("   generateDEHeatmap Arg: selected_contrast = %s", selected_contrast))
  message(sprintf("   generateDEHeatmap Arg: top_n_genes = %d", top_n_genes))

  if (is.null(de_results_list) || is.null(de_results_list$de_proteins_long)) {
    message("   generateDEHeatmap: No DE results available")
    return(NULL)
  }

  if (is.null(selected_contrast)) {
    message("   generateDEHeatmap: No contrast selected")
    return(NULL)
  }

  message("   generateDEHeatmap Step: Filtering data for selected contrast...")

  # Get the contrast-specific results
  de_proteins_long <- de_results_list$de_proteins_long
  
  # Filter for selected contrast and significant results
  contrast_data <- de_proteins_long |>
    dplyr::filter(comparison == selected_contrast) |>
    dplyr::filter(!!sym(qvalue_column) < de_q_val_thresh)

  if (nrow(contrast_data) == 0) {
    message(sprintf("   generateDEHeatmap: No significant genes found for contrast %s", selected_contrast))
    return(NULL)
  }

  message(sprintf("   generateDEHeatmap Step: Found %d significant genes", nrow(contrast_data)))

  # Order by absolute fold change and take top N
  top_genes_data <- contrast_data |>
    dplyr::arrange(desc(abs(!!sym(log2fc_column)))) |>
    head(top_n_genes)

  message(sprintf("   generateDEHeatmap Step: Selected top %d genes for heatmap", nrow(top_genes_data)))

  # Extract expression matrix for these genes
  protein_ids <- top_genes_data |> dplyr::pull(!!sym(names(top_genes_data)[1])) # First column should be protein ID

  # Get expression data from the S4 object
  expr_matrix <- de_results_list$theObject@protein_quant_table |>
    dplyr::filter(!!sym(de_results_list$theObject@protein_id_column) %in% protein_ids) |>
    column_to_rownames(de_results_list$theObject@protein_id_column) |>
    as.matrix()

  if (nrow(expr_matrix) == 0) {
    message("   generateDEHeatmap: No expression data found for selected genes")
    return(NULL)
  }

  message(sprintf("   generateDEHeatmap Step: Expression matrix dims = %d rows, %d cols", nrow(expr_matrix), ncol(expr_matrix)))

  # Apply scaling
  if (scale_data == "row") {
    message("   generateDEHeatmap Step: Applying row scaling...")
    expr_matrix <- t(scale(t(expr_matrix)))
  } else if (scale_data == "column") {
    message("   generateDEHeatmap Step: Applying column scaling...")
    expr_matrix <- scale(expr_matrix)
  } else if (scale_data == "both") {
    message("   generateDEHeatmap Step: Applying both row and column scaling...")
    expr_matrix <- scale(t(scale(t(expr_matrix))))
  } else {
    message("   generateDEHeatmap Step: No scaling applied")
  }

  message("   generateDEHeatmap Step: Setting up clustering parameters...")

  # Set up clustering
  row_clust <- NULL
  col_clust <- NULL

  if (cluster_rows || cluster_cols) {
    # Calculate distance matrices
    if (distance_method %in% c("pearson", "spearman")) {
      if (cluster_rows) {
        if (distance_method == "pearson") {
          row_dist <- as.dist(1 - cor(t(expr_matrix), method = "pearson", use = "complete.obs"))
        } else {
          row_dist <- as.dist(1 - cor(t(expr_matrix), method = "spearman", use = "complete.obs"))
        }
        row_clust <- hclust(row_dist, method = clustering_method)
      }
      
      if (cluster_cols) {
        if (distance_method == "pearson") {
          col_dist <- as.dist(1 - cor(expr_matrix, method = "pearson", use = "complete.obs"))
        } else {
          col_dist <- as.dist(1 - cor(expr_matrix, method = "spearman", use = "complete.obs"))
        }
        col_clust <- hclust(col_dist, method = clustering_method)
      }
    } else {
      # Standard distance metrics
      if (cluster_rows) {
        row_dist <- dist(expr_matrix, method = distance_method)
        row_clust <- hclust(row_dist, method = clustering_method)
      }
      
      if (cluster_cols) {
        col_dist <- dist(t(expr_matrix), method = distance_method)
        col_clust <- hclust(col_dist, method = clustering_method)
      }
    }
  }

  message("   generateDEHeatmap Step: Setting up color scheme...")

  # Choose color palette
  colors <- switch(color_scheme,
    "RdBu" = colorRampPalette(c("red", "white", "blue"))(100),
    "RdYlBu" = colorRampPalette(c("red", "yellow", "blue"))(100),
    "coolwarm" = colorRampPalette(c("blue", "white", "red"))(100),
    "viridis" = viridis::viridis(100),
    "plasma" = viridis::plasma(100),
    "inferno" = viridis::inferno(100),
    colorRampPalette(c("red", "white", "blue"))(100)
  )

  message("   generateDEHeatmap Step: Generating heatmap...")

  # Generate the heatmap
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    message("   generateDEHeatmap: Using ComplexHeatmap for advanced visualization")
    
    heatmap_plot <- ComplexHeatmap::Heatmap(
      expr_matrix,
      name = "Expression",
      col = colors,
      cluster_rows = if (cluster_rows) row_clust else FALSE,
      cluster_columns = if (cluster_cols) col_clust else FALSE,
      show_row_names = show_gene_names,
      show_column_names = TRUE,
      show_row_dend = cluster_rows,
      show_column_dend = cluster_cols,
      column_title = paste("Heatmap - Top", nrow(expr_matrix), "genes\nContrast:", selected_contrast),
      heatmap_legend_param = list(title = "Scaled\nExpression")
    )
  } else {
    message("   generateDEHeatmap: Using base R heatmap")
    
    # Create a plot function that returns the heatmap
    heatmap_plot <- function() {
      heatmap(
        expr_matrix,
        main = paste("Heatmap - Top", nrow(expr_matrix), "genes\nContrast:", selected_contrast),
        col = colors,
        Rowv = if (cluster_rows) as.dendrogram(row_clust) else NA,
        Colv = if (cluster_cols) as.dendrogram(col_clust) else NA,
        labRow = if (show_gene_names) rownames(expr_matrix) else rep("", nrow(expr_matrix)),
        cexRow = if (show_gene_names) 0.8 else 0.1
      )
    }
  }

  message("--- Exiting generateDEHeatmap ---")
  return(heatmap_plot)
}

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @export
setMethod(f = "outputDeResultsAllContrasts",
          signature = "ProteinQuantitativeData",
          definition = function(theObject,
                               de_results_list_all_contrasts = NULL,
                               uniprot_tbl = NULL,
                               de_output_dir = NULL,
                               publication_graphs_dir = NULL,
                               file_prefix = "de_proteins",
                               args_row_id = NULL,
                               gene_names_column = "gene_names",
                               uniprot_id_column = "Entry") {

  message("--- Entering outputDeResultsAllContrasts ---")
  message(sprintf("   outputDeResultsAllContrasts: de_output_dir = %s", de_output_dir))
  message(sprintf("   outputDeResultsAllContrasts: file_prefix = %s", file_prefix))
  
  # Extract parameters from S4 object with fallbacks
  uniprot_tbl <- checkParamsObjectFunctionSimplify(theObject, "uniprot_tbl", uniprot_tbl)
  de_output_dir <- checkParamsObjectFunctionSimplify(theObject, "de_output_dir", de_output_dir)
  publication_graphs_dir <- checkParamsObjectFunctionSimplify(theObject, "publication_graphs_dir", publication_graphs_dir)
  args_row_id <- checkParamsObjectFunctionSimplify(theObject, "args_row_id", args_row_id)
  gene_names_column <- checkParamsObjectFunctionSimplify(theObject, "gene_names_column", gene_names_column)
  uniprot_id_column <- checkParamsObjectFunctionSimplify(theObject, "uniprot_id_column", uniprot_id_column)
  
  # Ensure output directory exists (CRITICAL FIX!)
  if (!dir.exists(de_output_dir)) {
    dir.create(de_output_dir, recursive = TRUE, showWarnings = FALSE)
    message(sprintf("   outputDeResultsAllContrasts: Created output directory: %s", de_output_dir))
  }
  
  message(sprintf("   outputDeResultsAllContrasts: Processing %d contrasts", length(de_results_list_all_contrasts)))
  
  # Write results for each contrast with UNIQUE filenames
  for (contrast_name in names(de_results_list_all_contrasts)) {
    message(sprintf("   outputDeResultsAllContrasts: Writing files for contrast: %s", contrast_name))
    
    contrast_result <- de_results_list_all_contrasts[[contrast_name]]
    
    if (!is.null(contrast_result$de_proteins_long)) {
      # Clean contrast name for safe filename
      safe_contrast_name <- gsub("[^A-Za-z0-9_-]", "_", contrast_name)
      
             # Create annotated version (FIXED: proper conditional logic)
       de_proteins_long_annot <- contrast_result$de_proteins_long |>
         mutate(uniprot_acc_cleaned = str_split(!!sym(args_row_id), "-") |>
                  purrr::map_chr(1))
       
       # Add UniProt annotations if available
       if (!is.null(uniprot_tbl) && nrow(uniprot_tbl) > 0) {
         de_proteins_long_annot <- de_proteins_long_annot |>
           left_join(uniprot_tbl, by = join_by(uniprot_acc_cleaned == !!sym(uniprot_id_column))) |>
           dplyr::select(-uniprot_acc_cleaned)
       } else {
         de_proteins_long_annot <- de_proteins_long_annot |>
           dplyr::select(-uniprot_acc_cleaned)
       }
       
       # Add gene_name column with proper conditional logic
       if ("gene_names" %in% names(de_proteins_long_annot)) {
         de_proteins_long_annot <- de_proteins_long_annot |>
           mutate(gene_name = purrr::map_chr(gene_names, \(x){
             tryCatch({
               if(is.na(x) || is.null(x) || x == "") {
                 ""
               } else {
                 split_result <- str_split(x, " |:")[[1]]
                 if(length(split_result) > 0) split_result[1] else ""
               }
             }, error = function(e) "")
           }))
       } else if (gene_names_column %in% names(de_proteins_long_annot)) {
         de_proteins_long_annot <- de_proteins_long_annot |>
           mutate(gene_name = purrr::map_chr(.data[[gene_names_column]], \(x){
             tryCatch({
               if(is.na(x) || is.null(x) || x == "") {
                 ""
               } else {
                 split_result <- str_split(x, " |:")[[1]]
                 if(length(split_result) > 0) split_result[1] else ""
               }
             }, error = function(e) "")
           }))
       } else {
         de_proteins_long_annot <- de_proteins_long_annot |>
           mutate(gene_name = "")
       }
       
       # Relocate gene_name column
       de_proteins_long_annot <- de_proteins_long_annot |>
         relocate(gene_name, .after = !!sym(args_row_id))
      
      # CRITICAL FIX: Use contrast-specific filename!
      contrast_filename <- paste0(file_prefix, "_", safe_contrast_name, "_long_annot.tsv")
      output_path <- file.path(de_output_dir, contrast_filename)
      
      message(sprintf("   outputDeResultsAllContrasts: Writing %s", output_path))
      
      # Write the file
      vroom::vroom_write(de_proteins_long_annot, output_path)
      
      # Verify file was written
      if (file.exists(output_path)) {
        file_size <- file.size(output_path)
        message(sprintf("   outputDeResultsAllContrasts: SUCCESS - File written: %s (%d bytes)", 
                       contrast_filename, file_size))
      } else {
        message(sprintf("   outputDeResultsAllContrasts: ERROR - File NOT written: %s", output_path))
      }
      
      # Also write Excel version
      excel_filename <- paste0(file_prefix, "_", safe_contrast_name, "_long_annot.xlsx")
      excel_path <- file.path(de_output_dir, excel_filename)
      writexl::write_xlsx(de_proteins_long_annot, excel_path)
      
      message(sprintf("   outputDeResultsAllContrasts: Also wrote Excel: %s", excel_filename))
      
      # âœ… NEW: Generate volcano plot for this contrast
      message(sprintf("   outputDeResultsAllContrasts: Generating volcano plot for contrast: %s", contrast_name))
      
      # Extract parameters - try direct access first, then use checkParams
      de_q_val_thresh <- if (!is.null(theObject@args$outputDeResultsAllContrasts$de_q_val_thresh)) {
        theObject@args$outputDeResultsAllContrasts$de_q_val_thresh
      } else {
        checkParamsObjectFunctionSimplify(theObject, "de_q_val_thresh", 0.05)
      }
      
      fdr_column <- if (!is.null(theObject@args$outputDeResultsAllContrasts$fdr_column)) {
        theObject@args$outputDeResultsAllContrasts$fdr_column
      } else {
        checkParamsObjectFunctionSimplify(theObject, "fdr_column", "fdr_qvalue")
      }
      
      log2fc_column <- if (!is.null(theObject@args$outputDeResultsAllContrasts$log2fc_column)) {
        theObject@args$outputDeResultsAllContrasts$log2fc_column
      } else {
        checkParamsObjectFunctionSimplify(theObject, "log2fc_column", "log2FC")
      }
      
      # Prepare data for volcano plot (same logic as in differentialExpressionAnalysisHelper)
      volcano_data <- contrast_result$de_proteins_long |>
        mutate(lqm = -log10(!!sym(fdr_column))) |>
        dplyr::mutate(label = case_when(!!sym(fdr_column) < de_q_val_thresh ~ "Significant",
                                       TRUE ~ "Not sig.")) |>
        dplyr::mutate(colour = case_when(!!sym(fdr_column) < de_q_val_thresh ~ "purple",
                                        TRUE ~ "black")) |>
        dplyr::mutate(colour = factor(colour, levels = c("black", "purple")))
      
      # Generate the volcano plot
      volcano_plot <- plotOneVolcanoNoVerticalLines(
        volcano_data, 
        contrast_name,
        log_q_value_column = lqm,
        log_fc_column = !!sym(log2fc_column)
      )
      
      # Create Volcano_Plots directory if it doesn't exist
      volcano_dir <- file.path(publication_graphs_dir, "Volcano_Plots")
      if (!dir.exists(volcano_dir)) {
        dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
        message(sprintf("   outputDeResultsAllContrasts: Created volcano plots directory: %s", volcano_dir))
      }
      
      # Save volcano plot with contrast-specific filename
      volcano_png_path <- file.path(volcano_dir, paste0(safe_contrast_name, ".png"))
      volcano_pdf_path <- file.path(volcano_dir, paste0(safe_contrast_name, ".pdf"))
      
      # Save as PNG
      tryCatch({
        ggplot2::ggsave(volcano_png_path, volcano_plot, width = 7, height = 7, dpi = 300)
        message(sprintf("   outputDeResultsAllContrasts: Saved volcano plot PNG: %s", basename(volcano_png_path)))
      }, error = function(e) {
        message(sprintf("   outputDeResultsAllContrasts: ERROR saving PNG: %s", e$message))
      })
      
      # Save as PDF
      tryCatch({
        ggplot2::ggsave(volcano_pdf_path, volcano_plot, width = 7, height = 7)
        message(sprintf("   outputDeResultsAllContrasts: Saved volcano plot PDF: %s", basename(volcano_pdf_path)))
      }, error = function(e) {
        message(sprintf("   outputDeResultsAllContrasts: ERROR saving PDF: %s", e$message))
      })
    }
  }
  
  # âœ… NEW: Generate combined multi-page PDF with all volcano plots
  message("   outputDeResultsAllContrasts: Creating combined volcano plots PDF...")
  
  volcano_dir <- file.path(publication_graphs_dir, "Volcano_Plots")
  combined_pdf_path <- file.path(volcano_dir, "all_volcano_plots_combined.pdf")
  
  tryCatch({
    # Re-extract parameters for combined PDF generation (in case they were modified)
    de_q_val_thresh <- if (!is.null(theObject@args$outputDeResultsAllContrasts$de_q_val_thresh)) {
      theObject@args$outputDeResultsAllContrasts$de_q_val_thresh
    } else {
      0.05
    }
    
    fdr_column <- if (!is.null(theObject@args$outputDeResultsAllContrasts$fdr_column)) {
      theObject@args$outputDeResultsAllContrasts$fdr_column
    } else {
      "fdr_qvalue"
    }
    
    log2fc_column <- if (!is.null(theObject@args$outputDeResultsAllContrasts$log2fc_column)) {
      theObject@args$outputDeResultsAllContrasts$log2fc_column
    } else {
      "log2FC"
    }
    
    # Collect all volcano plots
    all_volcano_plots <- list()
    
    for (contrast_name in names(de_results_list_all_contrasts)) {
      contrast_result <- de_results_list_all_contrasts[[contrast_name]]
      
      if (!is.null(contrast_result$de_proteins_long)) {
        # Recreate the volcano plot
        volcano_data <- contrast_result$de_proteins_long |>
          mutate(lqm = -log10(!!sym(fdr_column))) |>
          dplyr::mutate(label = case_when(!!sym(fdr_column) < de_q_val_thresh ~ "Significant",
                                         TRUE ~ "Not sig.")) |>
          dplyr::mutate(colour = case_when(!!sym(fdr_column) < de_q_val_thresh ~ "purple",
                                          TRUE ~ "black")) |>
          dplyr::mutate(colour = factor(colour, levels = c("black", "purple")))
        
        volcano_plot <- plotOneVolcanoNoVerticalLines(
          volcano_data, 
          contrast_name,
          log_q_value_column = lqm,
          log_fc_column = !!sym(log2fc_column)
        )
        
        all_volcano_plots[[contrast_name]] <- volcano_plot
      }
    }
    
    # Generate multi-page PDF
    if (length(all_volcano_plots) > 0) {
      pdf(file = combined_pdf_path, width = 7, height = 7, onefile = TRUE)
      purrr::walk(all_volcano_plots, print)
      invisible(dev.off())
      message(sprintf("   outputDeResultsAllContrasts: Created combined PDF with %d volcano plots: %s", 
                     length(all_volcano_plots), basename(combined_pdf_path)))
    }
    
  }, error = function(e) {
    message(sprintf("   outputDeResultsAllContrasts: ERROR creating combined PDF: %s", e$message))
  })
  
  message("--- Exiting outputDeResultsAllContrasts ---")
  return(TRUE)
}) 
