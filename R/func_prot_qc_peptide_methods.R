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
             valid_values <- peptide_data[[raw_quantity_column]]
             valid_values <- valid_values[!is.na(valid_values) & !is.nan(valid_values) & !is.infinite(valid_values)]
             
             if (length(valid_values) == 0) {
               message("      WARNING: No valid intensity values found for threshold calculation.")
               min_peptide_intensity_threshold <- 0
             } else {
               min_peptide_intensity_threshold <- unname(ceiling( quantile( valid_values, na.rm=TRUE, probs = c(peptides_intensity_cutoff_percentile/100) ))[1])
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
             valid_values <- peptide_data[[norm_quantity_column]]
             valid_values <- valid_values[!is.na(valid_values) & !is.nan(valid_values) & !is.infinite(valid_values)]
             
             if (length(valid_values) == 0) {
               message("      WARNING: No valid intensity values found for threshold calculation.")
               min_peptide_intensity_threshold <- 0
             } else {
               min_peptide_intensity_threshold <- unname(ceiling( quantile( valid_values, na.rm=TRUE, probs = c(peptides_intensity_cutoff_percentile/100) ))[1])
             }
             message(sprintf("      Calculated min_peptide_intensity_threshold = %g", min_peptide_intensity_threshold))

             message("   removePeptidesWithMissingValuesPercent: Calling helper function...")
             theObject@peptide_data <- removePeptidesWithMissingValuesPercentHelper( 
                                                 input_table = peptide_data
                                               , design_matrix = design_matrix
                                               , sample_id = sample_id
                                               , protein_id_column = protein_id_column
                                               , peptide_sequence_column = peptide_sequence_column
                                               , grouping_variable = grouping_variable
                                               , groupwise_percentage_cutoff = groupwise_percentage_cutoff
                                               , max_groups_percentage_cutoff = max_groups_percentage_cutoff
                                               , abundance_threshold = min_peptide_intensity_threshold
                                               , abundance_column = norm_quantity_column )


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
                                                                                             , input_table_sample_id_column = sample_id_column
                                                                                             , sample_id_tbl_sample_id_column  = sample_id_column
                                                                                             , replicate_group_column = replicate_group_column
                                                                                             , protein_id_column = theObject@protein_id_column
                                                                                             , peptide_sequence_column = theObject@peptide_sequence_column
                                                                                             , core_utilisation = core_utilisation)
             theObject <- cleanDesignMatrixPeptide(theObject)

             theObject
           })

# ----------------------------------------------------------------------------
# filterMinNumPeptidesPerProtein
# ----------------------------------------------------------------------------
#'@title Filter the proteins based on the number of peptides and peptidoforms
#'@description Keep protein groups with at least one distinct stripped peptide
#'  and at least two distinct peptidoforms by default.
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

             parameter_section <- theObject@args$filterMinNumPeptidesPerProtein

             if (is.null(num_peptides_per_protein_thresh) && is.list(parameter_section)) {
               num_peptides_per_protein_thresh <- parameter_section$num_peptides_per_protein_thresh
               if (is.null(num_peptides_per_protein_thresh)) {
                 num_peptides_per_protein_thresh <- parameter_section$peptides_per_protein_cutoff
               }
             }

             if (is.null(num_peptidoforms_per_protein_thresh) && is.list(parameter_section)) {
               num_peptidoforms_per_protein_thresh <- parameter_section$num_peptidoforms_per_protein_thresh
               if (is.null(num_peptidoforms_per_protein_thresh)) {
                 num_peptidoforms_per_protein_thresh <- parameter_section$peptidoforms_per_protein_cutoff
               }
             }

             num_peptides_per_protein_thresh <- .resolvePeptideQcMethodParam(
               theObject = theObject,
               section_name = "filterMinNumPeptidesPerProtein",
               param_name = "num_peptides_per_protein_thresh",
               explicit_value = num_peptides_per_protein_thresh,
               default_value = 1
             )

             num_peptidoforms_per_protein_thresh <- .resolvePeptideQcMethodParam(
               theObject = theObject,
               section_name = "filterMinNumPeptidesPerProtein",
               param_name = "num_peptidoforms_per_protein_thresh",
               explicit_value = num_peptidoforms_per_protein_thresh,
               default_value = 2
             )

             core_utilisation <- .resolvePeptideQcMethodParam(
               theObject = theObject,
               section_name = "filterMinNumPeptidesPerProtein",
               param_name = "core_utilisation",
               explicit_value = core_utilisation,
               default_value = NA
             )


             if (is.null(theObject@args$filterMinNumPeptidesPerProtein)) {
               theObject@args$filterMinNumPeptidesPerProtein <- list()
             }
             theObject@args$filterMinNumPeptidesPerProtein$num_peptides_per_protein_thresh <- num_peptides_per_protein_thresh
             theObject@args$filterMinNumPeptidesPerProtein$num_peptidoforms_per_protein_thresh <- num_peptidoforms_per_protein_thresh
             theObject@args$filterMinNumPeptidesPerProtein$core_utilisation <- core_utilisation

             theObject@peptide_data <- filterMinNumPeptidesPerProteinHelper ( input_table = peptide_data
                                                                        , num_peptides_per_protein_thresh = num_peptides_per_protein_thresh
                                                                        , num_peptidoforms_per_protein_thresh = num_peptidoforms_per_protein_thresh
                                                                        , protein_id_column = protein_id_column
                                                                        , peptide_sequence_column = theObject@peptide_sequence_column
                                                                        , modified_peptide_sequence_column = "Modified.Sequence"
                                                                        , core_utilisation = core_utilisation)

             theObject <- cleanDesignMatrixPeptide(theObject)

             theObject
           })

# ----------------------------------------------------------------------------
# filterMinNumPeptidesPerSample
# ----------------------------------------------------------------------------
.resolvePeptideQcMethodParam <- function(theObject,
                                         section_name,
                                         param_name,
                                         explicit_value = NULL,
                                         default_value = NULL,
                                         accept_null = FALSE,
                                         check_fn = checkParamsObjectFunctionSimplify) {
  if (!is.null(explicit_value)) {
    return(explicit_value)
  }

  section <- theObject@args[[section_name]]
  has_object_value <- is.list(section) && param_name %in% names(section)
  object_value <- if (has_object_value) {
    section[[param_name]]
  } else {
    NULL
  }
  if (has_object_value && (!is.null(object_value) || accept_null)) {
    return(object_value)
  }

  checked_value <- tryCatch(
    suppressWarnings(check_fn(theObject, param_name, default_value)),
    error = function(e) NULL
  )
  if (!is.null(checked_value) || accept_null) {
    return(checked_value)
  }

  default_value
}

#'@export
setMethod( f="filterMinNumPeptidesPerSample"
           , signature="PeptideQuantitativeData"
           , definition = function( theObject
                                    , peptides_per_sample_cutoff = NULL
                                    , core_utilisation = NULL
                                    , inclusion_list = NULL) {

             peptide_data <- theObject@peptide_data
             sample_id_column <- theObject@sample_id

             peptides_per_sample_cutoff <- .resolvePeptideQcMethodParam(
               theObject = theObject,
               section_name = "filterMinNumPeptidesPerSample",
               param_name = "peptides_per_sample_cutoff",
               explicit_value = peptides_per_sample_cutoff,
               default_value = 5000
             )

             inclusion_list <- .resolvePeptideQcMethodParam(
               theObject = theObject,
               section_name = "filterMinNumPeptidesPerSample",
               param_name = "inclusion_list",
               explicit_value = inclusion_list,
               default_value = NULL,
               accept_null = TRUE,
               check_fn = checkParamsObjectFunctionSimplifyAcceptNull
             )
             if (is.null(inclusion_list) || identical(inclusion_list, "")) {
               inclusion_list <- character(0)
             }

             core_utilisation <- .resolvePeptideQcMethodParam(
               theObject = theObject,
               section_name = "filterMinNumPeptidesPerSample",
               param_name = "core_utilisation",
               explicit_value = core_utilisation,
               default_value = NA
             )

             theObject <- updateParamInObject(theObject, "peptides_per_sample_cutoff")
             theObject <- updateParamInObject(theObject, "inclusion_list")
             theObject <- updateParamInObject(theObject, "core_utilisation")

             theObject@peptide_data <- filterMinNumPeptidesPerSampleHelper( peptide_data
                                            , peptides_per_sample_cutoff = peptides_per_sample_cutoff
                                            , sample_id_column = sample_id_column
                                            , core_utilisation
                                            , inclusion_list = inclusion_list )

             theObject <- cleanDesignMatrixPeptide(theObject)

             theObject
           })

# ----------------------------------------------------------------------------
# srlQvalueProteotypicPeptideClean
# ----------------------------------------------------------------------------
.resolveSrlQvaluePeptideParam <- function(theObject,
                                          param_name,
                                          explicit_value = NULL,
                                          default_value = NULL,
                                          section_name = "srlQvalueProteotypicPeptideClean") {
  if (!is.null(explicit_value)) {
    return(explicit_value)
  }

  checked_value <- tryCatch(
    checkParamsObjectFunctionSimplify(theObject, param_name),
    error = function(e) NULL
  )
  if (!is.null(checked_value)) {
    return(checked_value)
  }

  section <- theObject@args[[section_name]]
  object_value <- if (is.list(section)) section[[param_name]] else NULL
  if (!is.null(object_value)) {
    return(object_value)
  }

  if (!is.null(default_value)) {
    return(default_value)
  }

  stop(sprintf("%s: '%s' is not defined.", section_name, param_name), call. = FALSE)
}

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

             dia_nn_default_columns <- c("Protein.Ids"
                                        , "Stripped.Sequence"
                                        , "Modified.Sequence"
                                        , "Q.Value"
                                        , "Global.Q.Value"
                                        , "Precursor.Quantity"
                                        , "Precursor.Normalised")

             qvalue_threshold <- .resolveSrlQvaluePeptideParam(
               theObject,
               "qvalue_threshold",
               explicit_value = qvalue_threshold,
               default_value = 0.01
             )

             global_qvalue_threshold <- .resolveSrlQvaluePeptideParam(
               theObject,
               "global_qvalue_threshold",
               explicit_value = global_qvalue_threshold,
               default_value = 0.01
             )

             choose_only_proteotypic_peptide <- .resolveSrlQvaluePeptideParam(
               theObject,
               "choose_only_proteotypic_peptide",
               explicit_value = choose_only_proteotypic_peptide,
               default_value = 1
             )

             input_matrix_column_ids <- .resolveSrlQvaluePeptideParam(
               theObject,
               "input_matrix_column_ids",
               explicit_value = input_matrix_column_ids,
               default_value = dia_nn_default_columns
             )

             if (!norm_quantity_column %in% names(peptide_data) &&
                 raw_quantity_column %in% names(peptide_data)) {
               peptide_data[[norm_quantity_column]] <- peptide_data[[raw_quantity_column]]
             }

             if (!global_q_value_column %in% names(peptide_data) &&
                 "PG.Q.Value" %in% names(peptide_data)) {
               peptide_data[[global_q_value_column]] <- peptide_data[["PG.Q.Value"]]
             }

             use_proteotypic_filter <- isTRUE(suppressWarnings(as.numeric(choose_only_proteotypic_peptide)) == 1)
             if (use_proteotypic_filter &&
                 !proteotypic_peptide_sequence_column %in% names(peptide_data)) {
               peptide_data[[proteotypic_peptide_sequence_column]] <- 1
             }

             if (is.null(theObject@args$srlQvalueProteotypicPeptideClean)) {
               theObject@args$srlQvalueProteotypicPeptideClean <- list()
             }
             theObject@args$srlQvalueProteotypicPeptideClean$qvalue_threshold <- qvalue_threshold
             theObject@args$srlQvalueProteotypicPeptideClean$global_qvalue_threshold <- global_qvalue_threshold
             theObject@args$srlQvalueProteotypicPeptideClean$choose_only_proteotypic_peptide <- choose_only_proteotypic_peptide
             theObject@args$srlQvalueProteotypicPeptideClean$input_matrix_column_ids <- input_matrix_column_ids

             # print( paste("qvalue_threshold: ", qvalue_threshold))
             search_srl_quant_cln <- srlQvalueProteotypicPeptideCleanHelper( input_table = peptide_data
                                                                       , input_matrix_column_ids = unique(c(input_matrix_column_ids
                                                                                                      , protein_id_column
                                                                                                      , peptide_sequence_column
                                                                                                      , "Modified.Sequence"))
                                                                       , protein_id_column = protein_id_column
                                                                       , peptide_sequence_column = peptide_sequence_column
                                                                       , modified_peptide_sequence_column = "Modified.Sequence"
                                                                       , q_value_column = q_value_column
                                                                       , global_q_value_column = global_q_value_column
                                                                       , proteotypic_peptide_sequence_column = proteotypic_peptide_sequence_column
                                                                       , global_qvalue_threshold = global_qvalue_threshold
                                                                       , qvalue_threshold = qvalue_threshold
                                                                       , choose_only_proteotypic_peptide = choose_only_proteotypic_peptide)

             theObject@peptide_data <- search_srl_quant_cln

             theObject <- cleanDesignMatrixPeptide(theObject)

             return(theObject)
           })
