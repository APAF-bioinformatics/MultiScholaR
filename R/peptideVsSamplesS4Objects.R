


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
                                         stop("Precursor raw quantity or normalized quantity column must be in the peptide data table")
                                       }

                                       if( ! object@sample_id %in% colnames(object@design_matrix) ) {
                                         stop("Sample ID column must be in the design matrix")
                                       }


                                       #Need to check the rows names in design matrix and the column names of the data table
                                       samples_in_peptide_data <-  object@peptide_data |> distinct(!!sym(object@sample_id)) |> pull(!!sym(object@sample_id))
                                       samples_in_design_matrix <- object@design_matrix |> pull( !! sym( object@sample_id ) )

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
#'@exportGeneric cleanDesignMatrixPeptide
setGeneric(name ="cleanDesignMatrixPeptide"
           , def=function( theObject) {
             standardGeneric("cleanDesignMatrixPeptide")
           })

#'@exportMethod cleanDesignMatrixPeptide
setMethod( f ="cleanDesignMatrixPeptide"
           , signature = "PeptideQuantitativeData"
           , definition=function( theObject ) {

             samples_id_vector <- theObject@peptide_data |> distinct(!!sym(theObject@sample_id)) |> pull(!!sym(theObject@sample_id))

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

             qvalue_threshold <- checkParamsObjectFunctionSimplify( theObject, "qvalue_threshold", qvalue_threshold, 0.01)

             global_qvalue_threshold <- checkParamsObjectFunctionSimplify( theObject, "global_qvalue_threshold", global_qvalue_threshold, 0.01)

             choose_only_proteotypic_peptide <- checkParamsObjectFunctionSimplify( theObject
                                                                                   , "choose_only_proteotypic_peptide"
                                                                                   , choose_only_proteotypic_peptide
                                                                                   , 1)

             theObject <- updateParamInObject(theObject, "qvalue_threshold", qvalue_threshold)
             theObject <- updateParamInObject(theObject, "global_qvalue_threshold", global_qvalue_threshold)
             theObject <- updateParamInObject(theObject, "choose_only_proteotypic_peptide", choose_only_proteotypic_peptide)

             dia_nn_default_columns <- c("Protein.Ids"
                                        , "Stripped.Sequence"
                                        , "Q.Value"
                                        , "Global.Q.Value"
                                        , "Precursor.Quantity"
                                        , "Precursor.Normalised")

             input_matrix_column_ids <- checkParamsObjectFunctionSimplify( theObject
                                                                     , "input_matrix_column_ids"
                                                                     , input_matrix_column_ids
                                                                     , diann_default_columns )

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

            core_utilisation <- checkParamsObjectFunctionSimplify( theObject, "core_utilisation", core_utilisation, NA)
            theObject <- updateParamInObject(theObject, "core_utilisation", core_utilisation)

            theObject@peptide_data <- rollUpPrecursorToPeptideHelper(input_table = peptide_data
                                                               , sample_id_column = !!sym(sample_id)
                                                               , protein_id_column = !!sym(protein_id_column)
                                                               , peptide_sequence_column = !!sym(peptide_sequence_column)
                                                               , precursor_quantity_column = !!sym(raw_quantity_column)
                                                               , precursor_normalized_column = !!sym(norm_quantity_column)
                                                               , core_utilisation = core_utilisation)

             theObject@raw_quantity_column   <- "Peptide.RawQuantity"
             theObject@norm_quantity_column <- "Peptide.Normalized"

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
             peptide_data <- theObject@peptide_data
             raw_quantity_column <- theObject@raw_quantity_column
             norm_quantity_column <- theObject@norm_quantity_column

             peptides_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify( theObject
                                                                                    , "peptides_intensity_cutoff_percentile"
                                                                                    , peptides_intensity_cutoff_percentile)

             peptides_proportion_of_samples_below_cutoff <- checkParamsObjectFunctionSimplify( theObject
                                                                                                , "peptides_proportion_of_samples_below_cutoff"
                                                                                                , peptides_proportion_of_samples_below_cutoff)

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject, "core_utilisation", core_utilisation, NA)

             min_peptide_intensity_threshold <- ceiling( quantile( peptide_data |> pull(!!sym(raw_quantity_column)), na.rm=TRUE, probs = c(peptides_intensity_cutoff_percentile) ))[1]

             theObject <- updateParamInObject(theObject, "peptides_intensity_cutoff_percentile", peptides_intensity_cutoff_percentile)
             theObject <- updateParamInObject(theObject, "peptides_proportion_of_samples_below_cutoff", peptides_proportion_of_samples_below_cutoff)
             theObject <- updateParamInObject(theObject, "core_utilisation", core_utilisation)

             peptide_normalized_pif_cln <- peptideIntensityFilteringHelper( peptide_data
                                                                      , min_peptide_intensity_threshold = min_peptide_intensity_threshold
                                                                      , peptides_proportion_of_samples_below_cutoff = peptides_proportion_of_samples_below_cutoff
                                                                      , protein_id_column = !!sym( theObject@protein_id_column)
                                                                      , peptide_sequence_column = !!sym(theObject@peptide_sequence_column)
                                                                      , peptide_quantity_column = !!sym(raw_quantity_column)
                                                                      , core_utilisation = core_utilisation)

             theObject@peptide_data <- peptide_normalized_pif_cln

             theObject <- cleanDesignMatrixPeptide(theObject)

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
                                                                   , grouping_variable
                                                                   , NULL)
             groupwise_percentage_cutoff <- checkParamsObjectFunctionSimplify( theObject
                                                                                   , "groupwise_percentage_cutoff"
                                                                                   , groupwise_percentage_cutoff
                                                                                   , 50)
             max_groups_percentage_cutoff <- checkParamsObjectFunctionSimplify( theObject
                                                                                   , "max_groups_percentage_cutoff"
                                                                                   , max_groups_percentage_cutoff
                                                                                   , 50)
             peptides_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify( theObject
                                                                                    , "peptides_intensity_cutoff_percentile"
                                                                                    , peptides_intensity_cutoff_percentile
                                                                                    , 50)

             theObject <- updateParamInObject(theObject, "grouping_variable", grouping_variable)
             theObject <- updateParamInObject(theObject, "groupwise_percentage_cutoff", groupwise_percentage_cutoff)
             theObject <- updateParamInObject(theObject, "max_groups_percentage_cutoff", max_groups_percentage_cutoff)
             theObject <- updateParamInObject(theObject, "peptides_intensity_cutoff_percentile", peptides_intensity_cutoff_percentile)

             min_protein_intensity_threshold <- ceiling( quantile( peptide_data |>
                                                                     dplyr::filter( !is.nan(!!sym(norm_quantity_column)) & !is.infinite(!!sym(norm_quantity_column))) |>
                                                                     pull(!!sym(norm_quantity_column))
                                                                   , na.rm=TRUE
                                                                   , probs = c(peptides_intensity_cutoff_percentile/100) ))[1]

             # print(min_protein_intensity_threshold )

             theObject@peptide_data <- removePeptidesWithMissingValuesPercentHelper( peptide_data
                                                                               , cols= !matches("row_id")
                                                                               , design_matrix = design_matrix
                                                                               , sample_id = !!sym(sample_id)
                                                                               , protein_id_column = !!sym(protein_id_column)
                                                                               , peptide_sequence_column = !!sym(peptide_sequence_column)
                                                                               , group_column = !!sym(grouping_variable)
                                                                               , groupwise_percentage_cutoff = groupwise_percentage_cutoff
                                                                               , max_groups_percentage_cutoff = max_groups_percentage_cutoff
                                                                               , abundance_threshold = min_protein_intensity_threshold
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
                                                                                   , num_peptides_per_protein_thresh
                                                                                   , 1)

             num_peptidoforms_per_protein_thresh <- checkParamsObjectFunctionSimplify( theObject
                                                                                       , "num_peptidoforms_per_protein_thresh"
                                                                                       , num_peptidoforms_per_protein_thresh
                                                                                       , 2)

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject, "core_utilisation", core_utilisation, NA)


             theObject <- updateParamInObject(theObject, "num_peptides_per_protein_thresh", num_peptides_per_protein_thresh)
             theObject <- updateParamInObject(theObject, "num_peptidoforms_per_protein_thresh", num_peptidoforms_per_protein_thresh)
             theObject <- updateParamInObject(theObject, "core_utilisation", core_utilisation)


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
                                                                              , peptides_per_sample_cutoff
                                                                              , 5000)

             inclusion_list <- checkParamsObjectFunctionSimplifyAcceptNull( theObject
                                                                            , "inclusion_list"
                                                                            , inclusion_list
                                                                            , NULL)

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject, "core_utilisation", core_utilisation, NA)

             theObject <- updateParamInObject(theObject, "peptides_per_sample_cutoff", peptides_per_sample_cutoff)
             theObject <- updateParamInObject(theObject, "inclusion_list", inclusion_list)
             theObject <- updateParamInObject(theObject, "core_utilisation", core_utilisation)

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
setGeneric( name="removePeptidesWithOnlyOneReplicate"
            , def=function( theObject, grouping_variable = NULL, core_utilisation = NULL, inclusion_list = NULL) {
              standardGeneric("removePeptidesWithOnlyOneReplicate")
            }
            , signature=c("theObject", "grouping_variable", "core_utilisation", "inclusion_list" ))

#'@export
setMethod( f="removePeptidesWithOnlyOneReplicate"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject, grouping_variable = NULL, core_utilisation = NULL, inclusion_list = c()) {

             peptide_data <- theObject@peptide_data
             sample_id_column <- theObject@sample_id
             replicate_group_column <- theObject@technical_replicate_id
             design_matrix <- theObject@design_matrix


             grouping_variable <- checkParamsObjectFunctionSimplifyAcceptNull( theObject
                                                                       , "grouping_variable"
                                                                       , grouping_variable
                                                                       , NULL)

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject
                                                           , "core_utilisation"
                                                           , core_utilisation
                                                           , NA)

             inclusion_list <- checkParamsObjectFunctionSimplifyAcceptNull( theObject
                                                                       , "inclusion_list"
                                                                       , inclusion_list
                                                                       , NULL)

             theObject <- updateParamInObject(theObject, "grouping_variable", grouping_variable)
             theObject <- updateParamInObject(theObject, "inclusion_list", inclusion_list)
             theObject <- updateParamInObject(theObject, "core_utilisation", core_utilisation)

             theObject@peptide_data <- removePeptidesWithOnlyOneReplicateHelper( input_table = peptide_data
                                                                                             , samples_id_tbl = design_matrix
                                                                                             , input_table_sample_id_column = !!sym(sample_id_column)
                                                                                             , sample_id_tbl_sample_id_column  = !!sym(sample_id_column)
                                                                                             , replicate_group_column = !!sym(grouping_variable)
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
                                                           , imputed_value_column
                                                           , NULL)

             proportion_missing_values <- checkParamsObjectFunctionSimplifyAcceptNull( theObject
                                                           , "proportion_missing_values"
                                                           , proportion_missing_values
                                                           , NULL)

             core_utilisation <- checkParamsObjectFunctionSimplify( theObject
                                                           , "core_utilisation"
                                                           , core_utilisation
                                                           , NA)

             theObject <- updateParamInObject(theObject, "imputed_value_column", imputed_value_column)
             theObject <- updateParamInObject(theObject, "proportion_missing_values", proportion_missing_values)
             theObject <- updateParamInObject(theObject, "core_utilisation", core_utilisation)

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


