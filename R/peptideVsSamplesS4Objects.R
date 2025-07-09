


## Create S4 class for proteomics protein level abundance data
#'@exportClass PeptideQuantitativeData
PeptideQuantitativeData <- setClass("PeptideQuantitativeData"

                                     , slots = c(
                                       # Protein vs Sample quantitative data
                                       peptide_data = "data.frame"
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
                                       protein_id_column = "Protein_Ids"
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
#' @export PeptideQuantitativeData

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
setGeneric(name="srlQvalueProteotypicPeptideClean"
           , def=function( theObject, qvalue_threshold = NULL, global_qvalue_threshold = NULL, choose_only_proteotypic_peptide = NULL, input_matrix_column_ids = NULL) {
             standardGeneric("srlQvalueProteotypicPeptideClean")
           }
           , signature = c("theObject", "qvalue_threshold", "global_qvalue_threshold", "choose_only_proteotypic_peptide", "input_matrix_column_ids") )


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
setGeneric(name="rollUpPrecursorToPeptide"
           , def=function( theObject, core_utilisation = NULL) {
             standardGeneric("rollUpPrecursorToPeptide")
           }
           , signature=c("theObject", "core_utilisation"))

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
setGeneric(name="peptideIntensityFiltering"
           , def=function( theObject, peptides_intensity_cutoff_percentile = NULL, peptides_proportion_of_samples_below_cutoff = NULL, core_utilisation = NULL) {
             standardGeneric("peptideIntensityFiltering")
           }
           , signature=c("theObject", "peptides_intensity_cutoff_percentile", "peptides_proportion_of_samples_below_cutoff", "core_utilisation"))

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
setGeneric(name="removePeptidesWithMissingValuesPercent"
           , def=function( theObject
                           , grouping_variable = NULL
                           , groupwise_percentage_cutoff = NULL
                           , max_groups_percentage_cutoff = NULL
                           , peptides_intensity_cutoff_percentile = NULL) {
             standardGeneric("removePeptidesWithMissingValuesPercent")
           }
           , signature=c("theObject"
                         , "grouping_variable"
                         , "groupwise_percentage_cutoff"
                         , "max_groups_percentage_cutoff"
                         , "peptides_intensity_cutoff_percentile" ))

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

#'@export
setGeneric(name="filterMinNumPeptidesPerProtein"
           , def=function( theObject
                           , num_peptides_per_protein_thresh = NULL
                           , num_peptidoforms_per_protein_thresh = NULL
                           , core_utilisation = NULL ) {
             standardGeneric("filterMinNumPeptidesPerProtein")
           }
           , signature=c("theObject"
                         , "num_peptides_per_protein_thresh"
                         , "num_peptidoforms_per_protein_thresh"
                         , "core_utilisation"))

#'@export
#'@description
#' Keep the proteins only if they have two or more peptides.
#'@param theObject Object of class PeptideQuantitativeData
#'@param num_peptides_per_protein_thresh Minimum number of peptides per protein
#'@param num_peptidoforms_per_protein_thresh Minimum number of peptidoforms per protein
#'@param core_utilisation core_utilisation to use for parallel processing
setMethod( f="filterMinNumPeptidesPerProtein"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject
                                    , num_peptides_per_protein_thresh = NULL
                                    , num_peptidoforms_per_protein_thresh = NULL
                                    , core_utilisation = NULL

                                    ) {
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
setGeneric( name="filterMinNumPeptidesPerSample"
            , def=function( theObject, peptides_per_sample_cutoff = NULL, core_utilisation = NULL, inclusion_list = NULL) {
              standardGeneric("filterMinNumPeptidesPerSample")
           }
           , signature=c("theObject", "peptides_per_sample_cutoff", "core_utilisation", "inclusion_list" ))

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

#' @export
setGeneric( name="removePeptidesWithOnlyOneReplicate"
            , def=function( theObject, replicate_group_column = NULL, core_utilisation = NULL) {
              standardGeneric("removePeptidesWithOnlyOneReplicate")
            }
            , signature=c("theObject", "replicate_group_column", "core_utilisation" ))

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

setGeneric( name="plotPeptidesProteinsCountsPerSample"
            , def=function( theObject ) {
              standardGeneric("plotPeptidesProteinsCountsPerSample")
            }
            , signature=c("theObject" ))



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
setGeneric( name="peptideMissingValueImputation"
            , def=function( theObject,  imputed_value_column = NULL, proportion_missing_values = NULL, core_utilisation = NULL) {
              standardGeneric("peptideMissingValueImputation")
            }
            , signature=c("theObject", "imputed_value_column", "proportion_missing_values", "core_utilisation" ))

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
setGeneric(name="calcPeptideMatrix"
           , def=function( theObject ) {
             standardGeneric("calcPeptideMatrix")
           }
           , signature=c("theObject"))

#'@export
setMethod(f="calcPeptideMatrix"
          , signature="PeptideQuantitativeData"
          , definition=function( theObject ) {
            
            peptide_data <- theObject@peptide_data
            raw_quantity_column <- theObject@raw_quantity_column
            sample_id_column <- theObject@sample_id
            protein_id_column <- theObject@protein_id_column
            peptide_sequence_column <- theObject@peptide_sequence_column

            peptide_quant_table <- peptide_data |>
              dplyr::select(!!sym(sample_id_column)
                            , !!sym(theObject@protein_id_column)
                            , !!sym(theObject@peptide_sequence_column)
                            , !!sym(raw_quantity_column)) |>
              mutate( !!sym(raw_quantity_column)  := purrr::map_dbl(!!sym(raw_quantity_column), as.numeric )) |>
              tidyr::pivot_wider(names_from = !!sym(sample_id_column)
                                 , values_from = !!sym(raw_quantity_column)) |>
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
setGeneric(name="getNegCtrlProtAnovaPeptides"
           , def=function( theObject
                           , ruv_grouping_variable  = NULL
                           , percentage_as_neg_ctrl  = NULL
                           , num_neg_ctrl  = NULL
                           , ruv_qval_cutoff = NULL
                           , ruv_fdr_method = NULL ) {
             standardGeneric("getNegCtrlProtAnovaPeptides")
           }
           , signature=c("theObject", "ruv_grouping_variable", "num_neg_ctrl", "ruv_qval_cutoff", "ruv_fdr_method"))

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
#'@export
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

#'@export
setMethod( f = "ruvIII_C_Varying"
           , signature="PeptideQuantitativeData"
           , definition=function( theObject, ruv_grouping_variable = NULL, ruv_number_k = NULL, ctrl = NULL) {
             
             peptide_matrix <- theObject@peptide_matrix
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

             
             Y <-  t( peptide_matrix[,design_matrix |> dplyr::pull(!!sym(sample_id))])
             
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
setGeneric(name="findBestNegCtrlPercentagePeptides"
           , def=function(theObject,
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
             standardGeneric("findBestNegCtrlPercentagePeptides")
           }
           , signature=c("theObject"))

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
      
      # Generate canonical correlation plot
      cancorplot <- ruvCancor(
        theObject,
        ctrl = control_genes_index,
        num_components_to_impute = num_components_to_impute,
        ruv_grouping_variable = ruv_grouping_variable
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
#' @param dpc_method Method for DPC estimation. Default is "eb" (empirical Bayes)
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
setGeneric(name="peptideMissingValueImputationLimpa"
           , def=function(theObject, 
                          imputed_value_column = NULL, 
                          use_log2_transform = TRUE,
                          dpc_method = "eb",
                          verbose = TRUE,
                          ensure_matrix = TRUE) {
             standardGeneric("peptideMissingValueImputationLimpa")
           }
           , signature=c("theObject"))

#' @export
setMethod(f="peptideMissingValueImputationLimpa"
          , signature="PeptideQuantitativeData"
          , definition = function(theObject, 
                                  imputed_value_column = NULL, 
                                  use_log2_transform = TRUE,
                                  dpc_method = "eb",
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
            
            dpc_method <- checkParamsObjectFunctionSimplify(
              theObject, "dpc_method", "eb"
            )
            
            verbose <- checkParamsObjectFunctionSimplify(
              theObject, "verbose", TRUE
            )
            
            # Update parameters in object
            theObject <- updateParamInObject(theObject, "imputed_value_column")
            theObject <- updateParamInObject(theObject, "use_log2_transform")
            theObject <- updateParamInObject(theObject, "dpc_method")
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
            } else if (!use_log2_transform && theObject@is_logged_data) {
              if (verbose) {
                log_info("Converting from log2 scale...")
              }
              y_peptide <- 2^y_peptide - 1
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
              log_info("Estimating detection probability curve using method: {dpc_method}")
            }
            
            tryCatch({
              dpcfit <- limpa::dpc(y_peptide, method = dpc_method)
              
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
              
              # Store DPC results in the object for future reference
              if (is.null(theObject@args)) {
                theObject@args <- list()
              }
              theObject@args$limpa_dpc_results <- list(
                dpc_parameters = dpcfit$dpc,
                dpc_method = dpc_method,
                missing_percentage_before = round(100 * mean(is.na(y_peptide)), 1),
                slope_interpretation = slope_interpretation
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