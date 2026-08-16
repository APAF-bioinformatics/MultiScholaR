# ----------------------------------------------------------------------------
# peptideIntensityFiltering
# ----------------------------------------------------------------------------
#'@export
setMethod( f="peptideIntensityFiltering"
           , signature="PeptideQuantitativeData"
           , definition = function(theObject,
                                   grouping_variable = NULL,
                                   groupwise_percentage_cutoff = NULL,
                                   max_groups_percentage_cutoff = NULL,
                                   peptides_intensity_cutoff_percentile = NULL,
                                   core_utilisation = NULL,
                                   min_reps_per_group = NULL,
                                   min_groups = NULL,
                                   strict_mode = NULL,
                                   peptide_quantity_column = NULL) {
             section <- theObject@args$peptideIntensityFiltering
             if (!is.list(section)) {
               section <- list()
             }
             resolve_section_value <- function(explicit_value, parameter_name, default_value = NULL) {
               if (!is.null(explicit_value)) {
                 return(explicit_value)
               }
               if (parameter_name %in% names(section) && !is.null(section[[parameter_name]])) {
                 return(section[[parameter_name]])
               }
               default_value
             }
             as_flag <- function(value, parameter_name) {
               if (is.logical(value) && length(value) == 1L && !is.na(value)) {
                 return(value)
               }
               if (is.numeric(value) && length(value) == 1L && value %in% c(0, 1)) {
                 return(as.logical(value))
               }
               if (is.character(value) && length(value) == 1L) {
                 normalised <- tolower(trimws(value))
                 if (normalised %in% c("true", "t", "yes", "y", "1")) return(TRUE)
                 if (normalised %in% c("false", "f", "no", "n", "0")) return(FALSE)
               }
               stop(
                 sprintf("peptideIntensityFiltering: `%s` must be TRUE or FALSE.", parameter_name),
                 call. = FALSE
               )
             }

             grouping_variable <- resolve_section_value(grouping_variable, "grouping_variable", "group")
             min_reps_per_group <- resolve_section_value(min_reps_per_group, "min_reps_per_group", NULL)
             min_groups <- resolve_section_value(min_groups, "min_groups", NULL)
             strict_mode <- as_flag(
               resolve_section_value(strict_mode, "strict_mode", FALSE),
               "strict_mode"
             )
             groupwise_percentage_cutoff <- resolve_section_value(
               groupwise_percentage_cutoff,
               "groupwise_percentage_cutoff",
               1
             )
             max_groups_percentage_cutoff <- resolve_section_value(
               max_groups_percentage_cutoff,
               "max_groups_percentage_cutoff",
               50
             )
             peptides_intensity_cutoff_percentile <- resolve_section_value(
               peptides_intensity_cutoff_percentile,
               "peptides_intensity_cutoff_percentile",
               1
             )
             core_utilisation <- resolve_section_value(core_utilisation, "core_utilisation", NA)
             peptide_quantity_column <- resolve_section_value(
               peptide_quantity_column,
               "peptide_quantity_column",
               theObject@norm_quantity_column
             )
             if (is.null(peptide_quantity_column) || !nzchar(trimws(peptide_quantity_column))) {
               peptide_quantity_column <- theObject@norm_quantity_column
             }

             percentile <- suppressWarnings(as.numeric(peptides_intensity_cutoff_percentile))
             if (length(percentile) != 1L || is.na(percentile) ||
                 !is.finite(percentile) || percentile < 0 || percentile > 100) {
               stop(
                 "peptideIntensityFiltering: `peptides_intensity_cutoff_percentile` must be between 0 and 100.",
                 call. = FALSE
               )
             }
             if (!peptide_quantity_column %in% names(theObject@peptide_data)) {
               stop(
                 sprintf(
                   "peptideIntensityFiltering: quantity column `%s` is absent from peptide_data.",
                   peptide_quantity_column
                 ),
                 call. = FALSE
               )
             }
             quantity_values <- theObject@peptide_data[[peptide_quantity_column]]
             if (!is.numeric(quantity_values)) {
               stop(
                 sprintf(
                   "peptideIntensityFiltering: quantity column `%s` must be numeric.",
                   peptide_quantity_column
                 ),
                 call. = FALSE
               )
             }
             valid_values <- quantity_values[is.finite(quantity_values)]
             min_peptide_intensity_threshold <- if (length(valid_values) == 0L) {
               0
             } else {
               unname(ceiling(stats::quantile(
                 valid_values,
                 na.rm = TRUE,
                 probs = percentile / 100
               ))[[1L]])
             }

             filter_result <- peptideIntensityFilteringHelper(
               input_table = theObject@peptide_data,
               design_matrix = theObject@design_matrix,
               min_peptide_intensity_threshold = min_peptide_intensity_threshold,
               sample_id_column = theObject@sample_id,
               grouping_variable = grouping_variable,
               groupwise_percentage_cutoff = groupwise_percentage_cutoff,
               max_groups_percentage_cutoff = max_groups_percentage_cutoff,
               protein_id_column = theObject@protein_id_column,
               peptide_sequence_column = theObject@peptide_sequence_column,
               peptide_quantity_column = peptide_quantity_column,
               core_utilisation = core_utilisation,
               min_reps_per_group = min_reps_per_group,
               min_groups = min_groups,
               strict_mode = strict_mode,
               return_filter_result = TRUE
             )
             if (is.data.frame(filter_result)) {
               filter_result <- list(
                 data = filter_result,
                 support_table = NULL,
                 removal_ledger = NULL,
                 summary = list(
                   rule_mode = "legacy_helper_result",
                   intensity_threshold = min_peptide_intensity_threshold,
                   intensity_quantity_column = peptide_quantity_column
                 )
               )
             }

             section$grouping_variable <- grouping_variable
             section$min_reps_per_group <- min_reps_per_group
             section$min_groups <- min_groups
             section$strict_mode <- strict_mode
             section$groupwise_percentage_cutoff <- groupwise_percentage_cutoff
             section$max_groups_percentage_cutoff <- max_groups_percentage_cutoff
             section$peptides_intensity_cutoff_percentile <- percentile
             section$core_utilisation <- core_utilisation
             section$peptide_quantity_column <- peptide_quantity_column
             section$min_peptide_intensity_threshold <- min_peptide_intensity_threshold
             section$support_table <- filter_result$support_table
             section$removal_ledger <- filter_result$removal_ledger
             section$filter_summary <- filter_result$summary
             theObject@args$peptideIntensityFiltering <- section
             theObject@peptide_data <- filter_result$data
             theObject <- cleanDesignMatrixPeptide(theObject)
             theObject
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
             section <- theObject@args$removePeptidesWithOnlyOneReplicate
             if (!is.list(section)) {
               section <- list()
             }
             if (is.null(replicate_group_column)) {
               replicate_group_column <- section$replicate_group_column
             }
             if (is.null(replicate_group_column)) {
               replicate_group_column <- section$grouping_variable
             }
             if (is.null(replicate_group_column) || !nzchar(trimws(replicate_group_column))) {
               replicate_group_column <- "replicates"
             }
             if (is.null(core_utilisation)) {
               core_utilisation <- section$core_utilisation
             }
             if (is.null(core_utilisation)) {
               core_utilisation <- NA
             }

             filter_result <- removePeptidesWithOnlyOneReplicateHelper(
               input_table = theObject@peptide_data,
               samples_id_tbl = theObject@design_matrix,
               input_table_sample_id_column = theObject@sample_id,
               sample_id_tbl_sample_id_column = theObject@sample_id,
               replicate_group_column = replicate_group_column,
               protein_id_column = theObject@protein_id_column,
               peptide_sequence_column = theObject@peptide_sequence_column,
               core_utilisation = core_utilisation,
               minimum_distinct_runs = 2L,
               retention_policy = "supported_in_any_group",
               return_filter_result = TRUE
             )
             if (is.data.frame(filter_result)) {
               filter_result <- list(
                 data = filter_result,
                 support_table = NULL,
                 removal_ledger = NULL,
                 summary = list(retention_policy = "legacy_helper_result")
               )
             }
             section$replicate_group_column <- replicate_group_column
             section$grouping_variable <- replicate_group_column
             section$core_utilisation <- core_utilisation
             section$minimum_distinct_runs <- 2L
             section$retention_policy <- "supported_in_any_group"
             section$support_table <- filter_result$support_table
             section$removal_ledger <- filter_result$removal_ledger
             section$filter_summary <- filter_result$summary
             theObject@args$removePeptidesWithOnlyOneReplicate <- section
             theObject@peptide_data <- filter_result$data
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
#' @param ... Filtering overrides, including `num_peptides_per_protein_thresh`,
#'   `num_peptidoforms_per_protein_thresh`, and `core_utilisation`.
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
               default_value = 500
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

             filter_result <- filterMinNumPeptidesPerSampleHelper(
               input_table = peptide_data,
               peptides_per_sample_cutoff = peptides_per_sample_cutoff,
               sample_id_column = sample_id_column,
               core_utilisation = core_utilisation,
               inclusion_list = inclusion_list,
               protein_id_column = theObject@protein_id_column,
               peptide_sequence_column = theObject@peptide_sequence_column,
               design_matrix = theObject@design_matrix,
               design_sample_id_column = sample_id_column,
               return_filter_result = TRUE
             )
             if (is.data.frame(filter_result)) {
               filter_result <- list(
                 data = filter_result,
                 support_table = NULL,
                 removal_ledger = NULL,
                 summary = list(count_definition = "legacy_helper_result")
               )
             }

             section <- theObject@args$filterMinNumPeptidesPerSample
             if (!is.list(section)) {
               section <- list()
             }
             section$peptides_per_sample_cutoff <- peptides_per_sample_cutoff
             section$inclusion_list <- inclusion_list
             section$core_utilisation <- core_utilisation
             section$support_table <- filter_result$support_table
             section$removal_ledger <- filter_result$removal_ledger
             section$filter_summary <- filter_result$summary
             theObject@args$filterMinNumPeptidesPerSample <- section
             theObject@peptide_data <- filter_result$data

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
                                  , global_pg_qvalue_threshold = NULL
                                  , confidence_provenance_mode = NULL
                                  , exclude_decoys = NULL
                                  , exclude_contaminants = NULL
                                  , contaminant_manifest = NULL
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
                                        , "Global.PG.Q.Value"
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

             global_pg_qvalue_threshold <- .resolveSrlQvaluePeptideParam(
               theObject,
               "global_pg_qvalue_threshold",
               explicit_value = global_pg_qvalue_threshold,
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

             confidence_provenance_mode <- .resolveSrlQvaluePeptideParam(
               theObject,
               "confidence_provenance_mode",
               explicit_value = confidence_provenance_mode,
               default_value = "diann_main_report"
             )

             exclude_decoys <- .resolveSrlQvaluePeptideParam(
               theObject,
               "exclude_decoys",
               explicit_value = exclude_decoys,
               default_value = TRUE
             )

             exclude_contaminants <- .resolveSrlQvaluePeptideParam(
               theObject,
               "exclude_contaminants",
               explicit_value = exclude_contaminants,
               default_value = TRUE
             )

             contaminant_manifest <- .resolveSrlQvaluePeptideParam(
               theObject,
               "contaminant_manifest",
               explicit_value = contaminant_manifest,
               default_value = ""
             )

             .validateDiaNnIdentificationContract(
               qvalue_threshold = qvalue_threshold,
               global_qvalue_threshold = global_qvalue_threshold,
               global_pg_qvalue_threshold = global_pg_qvalue_threshold,
               choose_only_proteotypic_peptide = choose_only_proteotypic_peptide
             )

             confidence_provenance_mode <- match.arg(
               confidence_provenance_mode,
               c("diann_main_report", "upstream_prefiltered")
             )
             global_pg_q_value_column <- "Global.PG.Q.Value"

             if (!norm_quantity_column %in% names(peptide_data) &&
                 raw_quantity_column %in% names(peptide_data)) {
               peptide_data[[norm_quantity_column]] <- peptide_data[[raw_quantity_column]]
             }

             if (is.null(theObject@args$srlQvalueProteotypicPeptideClean)) {
               theObject@args$srlQvalueProteotypicPeptideClean <- list()
             }
             theObject@args$srlQvalueProteotypicPeptideClean$qvalue_threshold <- qvalue_threshold
             theObject@args$srlQvalueProteotypicPeptideClean$global_qvalue_threshold <- global_qvalue_threshold
             theObject@args$srlQvalueProteotypicPeptideClean$global_pg_qvalue_threshold <- global_pg_qvalue_threshold
             theObject@args$srlQvalueProteotypicPeptideClean$choose_only_proteotypic_peptide <- choose_only_proteotypic_peptide
             theObject@args$srlQvalueProteotypicPeptideClean$input_matrix_column_ids <- input_matrix_column_ids
             theObject@args$srlQvalueProteotypicPeptideClean$confidence_provenance_mode <- confidence_provenance_mode
             theObject@args$srlQvalueProteotypicPeptideClean$exclude_decoys <- .asPeptideExclusionFlag(
               exclude_decoys,
               "exclude_decoys"
             )
             theObject@args$srlQvalueProteotypicPeptideClean$exclude_contaminants <- .asPeptideExclusionFlag(
               exclude_contaminants,
               "exclude_contaminants"
             )
             theObject@args$srlQvalueProteotypicPeptideClean$contaminant_manifest <- contaminant_manifest
             theObject@args$srlQvalueProteotypicPeptideClean$confidence_metric_columns <- list(
               q_value = q_value_column,
               global_q_value = global_q_value_column,
               global_pg_q_value = global_pg_q_value_column,
               proteotypic = proteotypic_peptide_sequence_column
             )
             confidence_metric_columns <- c(
               q_value = q_value_column,
               global_q_value = global_q_value_column,
               global_pg_q_value = global_pg_q_value_column,
               proteotypic = proteotypic_peptide_sequence_column
             )
             theObject@args$srlQvalueProteotypicPeptideClean$confidence_metrics_present <- stats::setNames(
               unname(confidence_metric_columns) %in% names(peptide_data),
               names(confidence_metric_columns)
             )

             # print( paste("qvalue_threshold: ", qvalue_threshold))
             filter_result <- srlQvalueProteotypicPeptideCleanHelper( input_table = peptide_data
                                                                       , input_matrix_column_ids = unique(c(input_matrix_column_ids
                                                                                                      , protein_id_column
                                                                                                      , peptide_sequence_column
                                                                                                      , "Modified.Sequence"))
                                                                       , protein_id_column = protein_id_column
                                                                       , peptide_sequence_column = peptide_sequence_column
                                                                       , modified_peptide_sequence_column = "Modified.Sequence"
                                                                       , q_value_column = q_value_column
                                                                       , global_q_value_column = global_q_value_column
                                                                       , global_pg_q_value_column = global_pg_q_value_column
                                                                       , proteotypic_peptide_sequence_column = proteotypic_peptide_sequence_column
                                                                       , global_qvalue_threshold = global_qvalue_threshold
                                                                       , global_pg_qvalue_threshold = global_pg_qvalue_threshold
                                                                       , qvalue_threshold = qvalue_threshold
                                                                       , choose_only_proteotypic_peptide = choose_only_proteotypic_peptide
                                                                       , confidence_provenance_mode = confidence_provenance_mode
                                                                       , exclude_decoys = exclude_decoys
                                                                       , exclude_contaminants = exclude_contaminants
                                                                       , contaminant_manifest = contaminant_manifest
                                                                       , protein_ids_column = "Protein.Ids"
                                                                       , return_exclusion_result = TRUE)

             if (is.data.frame(filter_result)) {
               search_srl_quant_cln <- filter_result
               exclusion_summary <- data.frame(
                 input_rows = nrow(filter_result),
                 excluded_rows = 0L,
                 biological_rows = nrow(filter_result),
                 classification_status = "legacy_helper_without_classification",
                 stringsAsFactors = FALSE
               )
               exclusion_ledger <- filter_result[0, , drop = FALSE]
               manifest_provenance <- readPeptideContaminantManifest(NULL)
             } else {
               search_srl_quant_cln <- filter_result$data
               exclusion_summary <- filter_result$exclusion_summary
               exclusion_ledger <- filter_result$exclusion_ledger
               manifest_provenance <- filter_result$contaminant_manifest
             }

             theObject@args$srlQvalueProteotypicPeptideClean$biological_exclusion_summary <- exclusion_summary
             theObject@args$srlQvalueProteotypicPeptideClean$biological_exclusion_ledger <- exclusion_ledger
             theObject@args$srlQvalueProteotypicPeptideClean$contaminant_manifest_provenance <- manifest_provenance

             theObject@peptide_data <- search_srl_quant_cln

             theObject <- cleanDesignMatrixPeptide(theObject)

             return(theObject)
           })
