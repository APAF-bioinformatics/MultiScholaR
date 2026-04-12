# ----------------------------------------------------------------------------
# peptideIntensityFiltering
# ----------------------------------------------------------------------------
#'@export
setMethod( f="peptideIntensityFiltering"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject, 
                                    grouping_variable = NULL, 
                                    groupwise_percentage_cutoff = NULL, 
                                    max_groups_percentage_cutoff = NULL, 
                                    peptides_intensity_cutoff_percentile = NULL, 
                                    core_utilisation = NULL) {
             message("--- Entering peptideIntensityFiltering S4 Method ---")
             
             peptide_data <- theObject@peptide_data
             raw_quantity_column <- theObject@raw_quantity_column
             norm_quantity_column <- theObject@norm_quantity_column

             message("   peptideIntensityFiltering: Extracting input data...")
             message(sprintf("      Arg: raw_quantity_column = %s", raw_quantity_column))
             message(sprintf("      Arg: norm_quantity_column = %s", norm_quantity_column))
             message(sprintf("      Data State (peptide_data): Dims = %d rows, %d cols", nrow(peptide_data), ncol(peptide_data)))
             message(sprintf("      Columns: %s", paste(colnames(peptide_data), collapse = ", ")))

             message("   peptideIntensityFiltering: Resolving parameters...")
             grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "grouping_variable", "group")
             groupwise_percentage_cutoff <- checkParamsObjectFunctionSimplify( theObject, "groupwise_percentage_cutoff", 1)
             max_groups_percentage_cutoff <- checkParamsObjectFunctionSimplify( theObject, "max_groups_percentage_cutoff", 50)
             peptides_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify( theObject, "peptides_intensity_cutoff_percentile", 1)
             core_utilisation <- checkParamsObjectFunctionSimplify( theObject, "core_utilisation", NA)

             message(sprintf("      Resolved: grouping_variable = %s", grouping_variable))
             message(sprintf("      Resolved: groupwise_cutoff = %g%%, max_groups_fail = %g%%", groupwise_percentage_cutoff, max_groups_percentage_cutoff))
             message(sprintf("      Resolved: intensity_percentile = %g%%", peptides_intensity_cutoff_percentile))

             message("   peptideIntensityFiltering: Updating parameters in S4 object...")
             theObject <- updateParamInObject(theObject, "grouping_variable")
             theObject <- updateParamInObject(theObject, "groupwise_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "max_groups_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "peptides_intensity_cutoff_percentile")
             theObject <- updateParamInObject(theObject, "core_utilisation")

             message("   peptideIntensityFiltering: Calculating intensity threshold...")
             # Get non-missing values for threshold calculation
             valid_values <- peptide_data |> dplyr::pull(!!sym(raw_quantity_column))
             valid_values <- valid_values[!is.na(valid_values) & !is.nan(valid_values) & !is.infinite(valid_values)]
             
             if (length(valid_values) == 0) {
               message("      WARNING: No valid intensity values found for threshold calculation.")
               min_peptide_intensity_threshold <- 0
             } else {
               min_peptide_intensity_threshold <- ceiling( quantile( valid_values, na.rm=TRUE, probs = c(peptides_intensity_cutoff_percentile/100) ))[1]
             }
             message(sprintf("      Calculated min_peptide_intensity_threshold = %g", min_peptide_intensity_threshold))

             message("   peptideIntensityFiltering: Calling helper function...")
             peptide_normalised_pif_cln <- peptideIntensityFilteringHelper( 
                                              input_table = peptide_data
                                              , design_matrix = theObject@design_matrix
                                              , min_peptide_intensity_threshold = min_peptide_intensity_threshold
                                              , sample_id_column = theObject@sample_id
                                              , grouping_variable = grouping_variable
                                              , groupwise_percentage_cutoff = groupwise_percentage_cutoff
                                              , max_groups_percentage_cutoff = max_groups_percentage_cutoff
                                              , protein_id_column = theObject@protein_id_column
                                              , peptide_sequence_column = theObject@peptide_sequence_column
                                              , peptide_quantity_column = raw_quantity_column
                                              , core_utilisation = core_utilisation)

             message(sprintf("   peptideIntensityFiltering: Helper returned %d rows", nrow(peptide_normalised_pif_cln)))

             theObject@peptide_data <- peptide_normalised_pif_cln

             message("   peptideIntensityFiltering: Cleaning design matrix...")
             theObject <- cleanDesignMatrixPeptide(theObject)

             message("--- Exiting peptideIntensityFiltering S4 Method ---")
             return(theObject)
           })

# ----------------------------------------------------------------------------
# removePeptidesWithMissingValuesPercent
# ----------------------------------------------------------------------------
#'@export
setMethod( f = "removePeptidesWithMissingValuesPercent"
           , signature="PeptideQuantitativeData"
           , definition=function( theObject
                                  , grouping_variable = NULL
                                  , groupwise_percentage_cutoff = NULL
                                  , max_groups_percentage_cutoff = NULL
                                  , peptides_intensity_cutoff_percentile = NULL) {
             
             message("--- Entering removePeptidesWithMissingValuesPercent S4 Method ---")
             
             peptide_data <- theObject@peptide_data
             protein_id_column <- theObject@protein_id_column
             peptide_sequence_column <- theObject@peptide_sequence_column
             raw_quantity_column <- theObject@raw_quantity_column
             norm_quantity_column <- theObject@norm_quantity_column
             sample_id <- theObject@sample_id
             design_matrix <- theObject@design_matrix

             message("   removePeptidesWithMissingValuesPercent: Resolving parameters...")
             grouping_variable <- checkParamsObjectFunctionSimplify( theObject, "grouping_variable", "group")
             groupwise_percentage_cutoff <- checkParamsObjectFunctionSimplify( theObject, "groupwise_percentage_cutoff", 50)
             max_groups_percentage_cutoff <- checkParamsObjectFunctionSimplify( theObject, "max_groups_percentage_cutoff", 50)
             peptides_intensity_cutoff_percentile <- checkParamsObjectFunctionSimplify( theObject, "peptides_intensity_cutoff_percentile", 1)

             message(sprintf("      Resolved: grouping_variable = %s", grouping_variable))
             message(sprintf("      Resolved: groupwise_cutoff = %g%%, max_groups_fail = %g%%", groupwise_percentage_cutoff, max_groups_percentage_cutoff))
             message(sprintf("      Resolved: intensity_percentile = %g%%", peptides_intensity_cutoff_percentile))

             theObject <- updateParamInObject(theObject, "grouping_variable")
             theObject <- updateParamInObject(theObject, "groupwise_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "max_groups_percentage_cutoff")
             theObject <- updateParamInObject(theObject, "peptides_intensity_cutoff_percentile")

             message("   removePeptidesWithMissingValuesPercent: Calculating intensity threshold...")
             # Filter out non-numeric/invalid values for threshold calculation
             valid_values <- peptide_data |> 
               dplyr::pull(!!sym(norm_quantity_column))
             valid_values <- valid_values[!is.na(valid_values) & !is.nan(valid_values) & !is.infinite(valid_values)]
             
             if (length(valid_values) == 0) {
               message("      WARNING: No valid intensity values found for threshold calculation.")
               min_peptide_intensity_threshold <- 0
             } else {
               min_peptide_intensity_threshold <- ceiling( quantile( valid_values, na.rm=TRUE, probs = c(peptides_intensity_cutoff_percentile/100) ))[1]
             }
             message(sprintf("      Calculated min_peptide_intensity_threshold = %g", min_peptide_intensity_threshold))

             message("   removePeptidesWithMissingValuesPercent: Calling helper function...")
             theObject@peptide_data <- removePeptidesWithMissingValuesPercentHelper( 
                                                 input_table = peptide_data
                                               , design_matrix = design_matrix
                                               , sample_id = !!sym(sample_id)
                                               , protein_id_column = !!sym(protein_id_column)
                                               , peptide_sequence_column = !!sym(peptide_sequence_column)
                                               , grouping_variable = !!sym(grouping_variable)
                                               , groupwise_percentage_cutoff = groupwise_percentage_cutoff
                                               , max_groups_percentage_cutoff = max_groups_percentage_cutoff
                                               , abundance_threshold = min_peptide_intensity_threshold
                                               , abundance_column = !!sym(norm_quantity_column) )


             message("   removePeptidesWithMissingValuesPercent: Cleaning design matrix...")
             theObject <- cleanDesignMatrixPeptide(theObject)

             message("--- Exiting removePeptidesWithMissingValuesPercent S4 Method ---")
             return(theObject)

           })

# ----------------------------------------------------------------------------
# removePeptidesWithOnlyOneReplicate
# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# filterMinNumPeptidesPerProtein
# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# filterMinNumPeptidesPerSample
# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# srlQvalueProteotypicPeptideClean
# ----------------------------------------------------------------------------
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

