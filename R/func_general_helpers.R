# ============================================================================
# func_general_helpers.R
# ============================================================================
# Purpose: General helper and utility functions
# 
# This file contains general-purpose helper functions used across all
# omics types and workflows, including parameter checking, parsing,
# data manipulation, correlation calculations, and result extraction.
# Functions in this file are used by multiple modules across the package.
#
# Functions to extract here:
# - Parameter checking functions (checkParamsObjectFunctionSimplify, etc.)
# - Parsing functions (parseList, parseNumList, parseString, parseType)
# - Data manipulation functions (list2df, list2graph, etc.)
# - Correlation functions (calculatePearsonCorrelationOptimized, etc.)
# - Result extraction functions (extractResults, extractRuvResults, etc.)
# - Filtering progress functions (updateProteinFiltering, etc.)
# - Additional utility helper functions
#
# Dependencies:
# - dplyr, tidyr, purrr
# - func_general_filemgmt.R (for file utilities)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# === Parameter Checking Functions ===

# Function 1: checkParamsObjectFunctionSimplify()
# Current location: R/helper_functions.R
# Description: Checks and simplifies parameters from object
# checkParamsObjectFunctionSimplify <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 2: checkParamsObjectFunctionSimplifyAcceptNull()
# Current location: R/helper_functions.R
# Description: Checks parameters accepting NULL values
# checkParamsObjectFunctionSimplifyAcceptNull <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 3: testRequiredArguments()
# Current location: R/helper_functions.R
# Description: Tests if required arguments are present
# testRequiredArguments <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 4: testRequiredFiles()
# Current location: R/helper_functions.R
# Description: Tests if required files exist
# testRequiredFiles <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 5: testRequiredFilesWarning()
# Current location: R/helper_functions.R
# Description: Tests required files with warnings
# testRequiredFilesWarning <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 6: isArgumentDefined()
# Current location: R/helper_functions.R
# Description: Checks if argument is defined
# isArgumentDefined <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 7: getFunctionName()
# Current location: R/helper_functions.R
# Description: Gets current function name
# getFunctionName <- function() {
#   # Extract from R/helper_functions.R
# }

# Function 8: getFunctionNameSecondLevel()
# Current location: R/helper_functions.R
# Description: Gets function name from second level
# getFunctionNameSecondLevel <- function() {
#   # Extract from R/helper_functions.R
# }

# === Parsing Functions ===

# Function 9: parseList()
# Current location: R/helper_functions.R
# Description: Parses list from arguments
# parseList <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 10: parseNumList()
# Current location: R/enrichment_functions.R
# Description: Parses numeric list from text
# parseNumList <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 11: parseString()
# Current location: R/helper_functions.R
# Description: Parses string from arguments
# parseString <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 12: parseType()
# Current location: R/helper_functions.R
# Description: Parses type from arguments
# parseType <- function(...) {
#   # Extract from R/helper_functions.R
# }

# === Data Manipulation Functions ===

# Function 13: list2df()
# Current location: R/enrichment_functions.R
# Description: Converts list to data frame
# list2df <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 14: list2graph()
# Current location: R/enrichment_functions.R
# Description: Converts list to graph
# list2graph <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 15: listifyTableByColumn()
# Current location: R/enrichment_functions.R
# Description: Converts table to list by column
# listifyTableByColumn <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 16: changeToCategorical()
# Current location: R/enrichment_functions.R
# Description: Changes variable to categorical
# changeToCategorical <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 17: getTypeOfGrouping()
# Current location: R/enrichment_functions.R
# Description: Gets type of grouping variable
# getTypeOfGrouping <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === Averaging and Replicate Functions ===

# Function 18: averageValuesFromReplicates()
# Current location: R/enrichment_functions.R
# Description: Averages values from replicates
# averageValuesFromReplicates <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 19: avgReplicateProteinIntensity()
# Current location: R/enrichment_functions.R
# Description: Averages replicate protein intensities
# avgReplicateProteinIntensity <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === Missing Value Functions ===

# Function 20: removeRowsWithMissingValues()
# Current location: R/helper_functions.R
# Description: Removes rows with missing values
# removeRowsWithMissingValues <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 21: removeRowsWithMissingValuesPercentHelper()
# Current location: R/helper_functions.R
# Description: Helper for removing rows with missing values
# removeRowsWithMissingValuesPercentHelper <- function(...) {
#   # Extract from R/helper_functions.R
# }

# === Correlation Functions ===

# Function 22: calculatePearsonCorrelationOptimized()
# Current location: R/enrichment_functions.R
# Description: Calculates Pearson correlation (optimized)
# calculatePearsonCorrelationOptimized <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 23: calulatePearsonCorrelation()
# Current location: R/enrichment_functions.R
# Description: Calculates Pearson correlation
# calulatePearsonCorrelation <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 24: calulatePearsonCorrelationForSamplePairsHelper()
# Current location: R/enrichment_functions.R
# Description: Helper for calculating correlation for sample pairs
# calulatePearsonCorrelationForSamplePairsHelper <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 25: filterSamplesByPeptideCorrelationThreshold()
# Current location: R/enrichment_functions.R
# Description: Filters samples by peptide correlation threshold
# filterSamplesByPeptideCorrelationThreshold <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 26: filterSamplesByProteinCorrelationThresholdHelper()
# Current location: R/enrichment_functions.R
# Description: Helper for filtering by protein correlation
# filterSamplesByProteinCorrelationThresholdHelper <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 27: findSamplesPairBelowPeptideCorrelationThreshold()
# Current location: R/enrichment_functions.R
# Description: Finds sample pairs below correlation threshold
# findSamplesPairBelowPeptideCorrelationThreshold <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 28: getPairsOfSamplesTable()
# Current location: R/enrichment_functions.R
# Description: Gets table of sample pairs
# getPairsOfSamplesTable <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 29: compareTwoPeptideDataObjects()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Compares two peptide data objects
# compareTwoPeptideDataObjects <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# === Comparison and UMAP Functions ===

# Function 30: compareUmapComponentsPairs()
# Current location: R/de_proteins_functions.R
# Description: Compares UMAP component pairs
# compareUmapComponentsPairs <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 31: umap_factor_plot()
# Current location: R/de_proteins_functions.R
# Description: Creates UMAP factor plot
# umap_factor_plot <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# === Data Extraction Functions ===

# Function 32: getRowsToKeepList()
# Current location: R/enrichment_functions.R
# Description: Gets list of rows to keep
# getRowsToKeepList <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 33: getSignificantData()
# Current location: R/enrichment_functions.R
# Description: Gets significant data
# getSignificantData <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 34: extractResults()
# Current location: R/enrichment_functions.R
# Description: Extracts results
# extractResults <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 35: extractRuvResults()
# Current location: R/helper_functions.R
# Description: Extracts RUV results
# extractRuvResults <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 36: extract_experiment()
# Current location: R/helper_functions.R
# Description: Extracts experiment data
# extract_experiment <- function(...) {
#   # Extract from R/helper_functions.R
# }

# === Enrichment Result Functions ===

# Function 37: createDEResultsForEnrichment()
# Current location: R/enrichment_functions.R
# Description: Creates DE results for enrichment
# createDEResultsForEnrichment <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 38: createEnrichmentResults()
# Current location: R/enrichment_functions.R
# Description: Creates enrichment results
# createEnrichmentResults <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 39: getEnrichmentResult()
# Current location: R/enrichment_functions.R
# Description: Gets enrichment result
# getEnrichmentResult <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 40: getEnrichmentSummary()
# Current location: R/enrichment_functions.R
# Description: Gets enrichment summary
# getEnrichmentSummary <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 41: summarize_enrichment()
# Current location: R/enrichment_functions.R
# Description: Summarizes enrichment results
# summarize_enrichment <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 42: processEnrichments()
# Current location: R/enrichment_functions.R
# Description: Processes enrichment results
# processEnrichments <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 43: perform_enrichment()
# Current location: R/enrichment_functions.R
# Description: Performs enrichment analysis
# perform_enrichment <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 44: generate_enrichment_plots()
# Current location: R/enrichment_functions.R
# Description: Generates enrichment plots
# generate_enrichment_plots <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === DE Statistics Functions ===

# Function 45: countStatDeGenes()
# Current location: R/enrichment_functions.R
# Description: Counts statistically significant DE genes
# countStatDeGenes <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 46: countStatDeGenesHelper()
# Current location: R/enrichment_functions.R
# Description: Helper for counting DE genes
# countStatDeGenesHelper <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 47: printCountDeGenesTable()
# Current location: R/enrichment_functions.R
# Description: Prints table of DE gene counts
# printCountDeGenesTable <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 48: printPValuesDistribution()
# Current location: R/de_proteins_functions.R
# Description: Prints p-value distribution
# printPValuesDistribution <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# === Parameter Update Functions ===

# Function 49: get_param_change_message()
# Current location: R/enrichment_functions.R
# Description: Gets parameter change message
# get_param_change_message <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 50: updateParamInObject()
# Current location: R/helper_functions.R
# Description: Updates parameter in object
# updateParamInObject <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 51: update_n()
# Current location: R/enrichment_functions.R
# Description: Updates n value
# update_n <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 52: setArgsDefault()
# Current location: R/helper_functions.R
# Description: Sets default arguments
# setArgsDefault <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 53: updateMissingValueParameters()
# Current location: R/helper_functions.R
# Description: Updates missing value parameters
# updateMissingValueParameters <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 54: updateRuvParameters()
# Current location: R/helper_functions.R
# Description: Updates RUV parameters
# updateRuvParameters <- function(...) {
#   # Extract from R/helper_functions.R
# }

# === Filtering Progress Functions ===

# Function 55: updateProteinFiltering()
# Current location: R/QC_visualisation.R
# Description: Updates and visualizes protein filtering progress
# updateProteinFiltering <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 56: updateFilteringProgressMetabolomics()
# Current location: R/QC_visualisation.R
# Description: Updates metabolomics filtering progress
# updateFilteringProgressMetabolomics <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 57: filtering_progress()
# Current location: R/enrichment_functions.R
# Description: Gets filtering progress
# filtering_progress <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 58: getFilteringProgressMetabolomics()
# Current location: R/QC_visualisation.R
# Description: Gets metabolomics filtering progress
# getFilteringProgressMetabolomics <- function() {
#   # Extract from R/QC_visualisation.R
# }

# === Summary Functions ===

# Function 59: summarizeQCPlot()
# Current location: R/helper_functions.R
# Description: Summarizes QC plot
# summarizeQCPlot <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 60: summarisePeptideObject()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Summarizes peptide object
# summarisePeptideObject <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# === Utility Functions ===

# Function 61: calcHtSize()
# Current location: R/enrichment_functions.R
# Description: Calculates heatmap size
# calcHtSize <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 62: get_ggrepel_segsize()
# Current location: R/enrichment_functions.R
# Description: Gets ggrepel segment size
# get_ggrepel_segsize <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 63: createIdToAttributeHash()
# Current location: R/helper_functions.R
# Description: Creates ID to attribute hash
# createIdToAttributeHash <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 64: convertKeyToAttribute()
# Current location: R/helper_functions.R
# Description: Converts key to attribute using hash
# convertKeyToAttribute <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 65: loadDependencies()
# Current location: R/helper_functions.R
# Description: Loads package dependencies
# loadDependencies <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 66: processAndFilterData()
# Current location: R/enrichment_functions.R
# Description: Processes and filters data
# processAndFilterData <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 67: peptidesIntensityMatrixPivotLonger()
# Current location: R/enrichment_functions.R
# Description: Pivots peptide intensity matrix to long format
# peptidesIntensityMatrixPivotLonger <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 68: proteinIntensityMatrixPivotLonger()
# Current location: R/enrichment_functions.R
# Description: Pivots protein intensity matrix to long format
# proteinIntensityMatrixPivotLonger <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 69: proteinTechRepCorrelationHelper()
# Current location: R/enrichment_functions.R
# Description: Helper for protein technical replicate correlation
# proteinTechRepCorrelationHelper <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 70: createWordCloudDataFrame()
# Current location: R/enrichment_functions.R
# Description: Creates word cloud data frame
# createWordCloudDataFrame <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 71: saveTimeRecord()
# Current location: R/enrichment_functions.R
# Description: Saves time record
# saveTimeRecord <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 72: appender_shiny()
# Current location: R/enrichment_functions.R
# Description: Appender for Shiny logging
# appender_shiny <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 73: taxonIdToGprofilerOrganism()
# Current location: R/enrichment_functions.R
# Description: Converts taxon ID to gprofiler organism
# taxonIdToGprofilerOrganism <- function(...) {
#   # Extract from R/enrichment_functions.R
# }


# ----------------------------------------------------------------------------
# saveTimeRecord
# ----------------------------------------------------------------------------
#---------------------------------------------------------------------------------------
#' @title Save Time Record
#' @description Save timing information for compute steps
#' @param compute_time_record Data frame to store timing records
#' @param step_name Name of the step being timed
#' @param toc_output Output from tictoc::toc()
#' @importFrom grid convertX convertY unit
#' @export
saveTimeRecord <- function ( compute_time_record, step_name, toc_output ) {

  if(length(which(compute_time_record[,"step_name"] == "")) == 0) {
    stop("saveTimeRecord: No more room on compute time record table.")
  }
  # find_first_empty_slot
  first_empty_slot <- which(compute_time_record[,"step_name"] == "")[1]

  print(first_empty_slot)

  compute_time_record[first_empty_slot, "step_name"] <- step_name
  compute_time_record[first_empty_slot, "time_elapsed"] <- toc_output$callback_msg

  return(compute_time_record)
}


# ----------------------------------------------------------------------------
# changeToCategorical
# ----------------------------------------------------------------------------
#' getOneContinousPalette
#' @export
changeToCategorical <- function(metadata_tbl, column_name, num_colours=9) {

  list_of_values <-  metadata_tbl |>
    dplyr::select( all_of(column_name)  ) |>
    #dplyr::filter(!is.na(!!sym(column_name))) |>
    pull()

  min_value <- min(list_of_values, na.rm =TRUE)
  max_value <-  max(list_of_values, na.rm =TRUE)

  if(min_value > 1) {
    min_value <- floor(min_value)
  }

  if(max_value > 1) {
    max_value <- ceiling(max_value)
  }

  formatted_list_of_values <- cut(list_of_values, breaks=seq( min_value, max_value, length.out=num_colours) )

  formatted_list_of_values
}


# ----------------------------------------------------------------------------
# peptidesIntensityMatrixPivotLonger
# ----------------------------------------------------------------------------
#' @title Pivot Peptide Matrix Long
#' @description
#'  Pivot peptide intensity matrix into long format table.
#' @export
peptidesIntensityMatrixPivotLonger <- function( input_matrix
                                                , sample_id_column
                                                , sequence_column
                                                , protein_id_column
                                                , quantity_column
                                                , unlog_data = TRUE) {

  output_matrix <- input_matrix |>
    as.data.frame() |>
    rownames_to_column(protein_id_column) |>
    pivot_longer(cols=!contains(protein_id_column)
                 , names_to = sample_id_column
                 , values_to =  quantity_column ) |>
    separate( col = protein_id_column
              , into=c(protein_id_column, "Stripped.Sequence"), sep="_")

  if ( unlog_data == TRUE) {
    output_matrix <- output_matrix|>
      mutate( {{quantity_column}} := 2^(!!sym(quantity_column)))

  }

  output_matrix
}


# ----------------------------------------------------------------------------
# proteinIntensityMatrixPivotLonger
# ----------------------------------------------------------------------------
#' @title Pivot Protein Matrix Long
#' @description
#' Pivot protein intensity matrix into long format
#' @export
proteinIntensityMatrixPivotLonger <- function( input_matrix
                                               , sample_id_column
                                               , protein_id_column
                                               , quantity_column) {

  output_matrix <- input_matrix |>
    as.data.frame() |>
    rownames_to_column(protein_id_column) |>
    pivot_longer(cols=!contains(protein_id_column)
                 , names_to = sample_id_column
                 , values_to =  quantity_column )

  output_matrix
}


# ----------------------------------------------------------------------------
# createIdToAttributeHash
# ----------------------------------------------------------------------------
## Output:
## An environment that act as a hash to convert keys to attributes.
#' @export
createIdToAttributeHash <- function(keys, attributes) {
	keys <- as.character( as.vector(keys))
	attribute <- as.vector(attributes)

	hash <- new.env(hash = TRUE, parent = parent.frame())

	if ( length(keys) != length(attributes))  {
		warning('Length of keys != Length of attributes list.')
		return(1)
	}

	mapply(function(k, a) base::assign(k, a, envir = hash), 
		   keys, attributes)

	return(hash)
}


# ----------------------------------------------------------------------------
# convertKeyToAttribute
# ----------------------------------------------------------------------------
## Output:
## A value that correspond to the query key value.
#' @export
convertKeyToAttribute <- function(key, hash) {

  if ( base::exists(key, hash) ) {
    return ( base::get(key, hash))
  }
  else {
    return (NA)
  }
}


# ----------------------------------------------------------------------------
# setArgsDefault
# ----------------------------------------------------------------------------
#' @export
setArgsDefault <- function(args, value_name, as_func, default_val=NA ) {

  if(isArgumentDefined(args, value_name))
  {
    args<-parseType(args,
                    c(value_name)
                    ,as_func)
  }else {
    logwarn(paste0( value_name, " is undefined, default value set to ", as.character(default_val), "."))
    args[[ value_name ]] <- default_val
  }

  return(args)
}


# ----------------------------------------------------------------------------
# getFunctionName
# ----------------------------------------------------------------------------
#' @export
getFunctionName <- function() {
  calls <- sys.calls()
  current_call <- calls[[length(calls) - 1]]
  as.character(current_call[1])
}


# ----------------------------------------------------------------------------
# getFunctionNameSecondLevel
# ----------------------------------------------------------------------------
#' @export
getFunctionNameSecondLevel <- function() {
  calls <- sys.calls()
  current_call <- calls[[length(calls) - 2]]
  as.character(current_call[1])
}


# ----------------------------------------------------------------------------
# checkParamsObjectFunctionSimplify
# ----------------------------------------------------------------------------
#' Check the parameters in the arguments list and the function parameters to see what param applies
#' @export
checkParamsObjectFunctionSimplify <- function(theObject, param_name_string, default_value = NULL) {

  function_name <- getFunctionNameSecondLevel()

  # print(function_name)
  param_value <- dynGet(param_name_string)

  # Fix: Safely access nested list to avoid index errors
  object_value <- NULL
  if (!is.null(theObject@args) && 
      !is.null(theObject@args[[function_name]]) && 
      is.list(theObject@args[[function_name]])) {
    object_value <- theObject@args[[function_name]][[param_name_string]]
  }

  # print(paste0("param_value = ", param_value))

  error <- paste0(function_name,  paste0(": '", param_name_string, "' is not defined.\n") )

  if( !is.null(param_value) ) {
    return( param_value)
  } else if( !is.null(object_value) ) {
    # print("use object value")
    return( object_value)
  } else if( !is.null(default_value) ) {
    return( default_value)
  } else {
    stop( error )
  }

}


# ----------------------------------------------------------------------------
# checkParamsObjectFunctionSimplifyAcceptNull
# ----------------------------------------------------------------------------
#' Check the parameters in the arguments list and the function parameters to see what param applies
#' @export
checkParamsObjectFunctionSimplifyAcceptNull <- function(theObject, param_name_string, default_value = NULL) {

  function_name <- getFunctionNameSecondLevel()

  # print(function_name)
  param_value <- dynGet(param_name_string)
  object_value <- theObject@args[[function_name]][[param_name_string]]

  # print(paste0("param_value = ", param_value))

  error <- paste0(function_name, ": '", param_name_string, "' is not defined.\n")

  if( !is.null(param_value) ) {
    return( param_value)
  } else if( !is.null(object_value) ) {
    return( object_value)
  } else if ( !is.null(default_value) ) {
    return( default_value)
  } else {
    warning(error)
    return( NULL )
  }

}


# ----------------------------------------------------------------------------
# updateParamInObject
# ----------------------------------------------------------------------------
#' Update the parameter in the object
#'@export
updateParamInObject <- function(theObject, param_name_string) {

  function_name <- getFunctionNameSecondLevel()

  theObject@args[[function_name]][[param_name_string]] <- dynGet(param_name_string)

  theObject

}


# ----------------------------------------------------------------------------
# updateConfigParameter
# ----------------------------------------------------------------------------
#' @title Update Parameter in S4 Object Args and Global Config List
#' @description Modifies a specific parameter within an S4 object's @args slot
#'              and also updates the corresponding value in a global list named
#'              'config_list'.
#'
#' @param theObject The S4 object whose @args slot needs updating.
#' @param function_name The name identifying the parameter section (character string,
#'                      e.g., "peptideIntensityFiltering"). Corresponds to the
#'                      first-level key in both @args and config_list.
#' @param parameter_name The specific parameter name to update (character string,
#'                       e.g., "peptides_proportion_of_samples_below_cutoff").
#'                       Corresponds to the second-level key.
#' @param new_value The new value to assign to the parameter.
#' @param config_list_name The name of the global list variable holding the
#'                         configuration (defaults to "config_list").
#' @param env The environment where the global config list resides (defaults to
#'            .GlobalEnv).
#'
#' @return The modified S4 object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assume 'myPeptideData' is a PeptideQuantitativeData object
#' # Assume 'config_list' exists in the global environment
#'
#' # Check initial values (example)
#' # print(myPeptideData@args$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff)
#' # print(config_list$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff)
#'
#' # Update the parameter to 0.7
#' myPeptideData <- updateConfigParameter(
#'   theObject = myPeptideData,
#'   function_name = "peptideIntensityFiltering",
#'   parameter_name = "peptides_proportion_of_samples_below_cutoff",
#'   new_value = 0.7
#' )
#'
#' # Verify changes (example)
#' # print(myPeptideData@args$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff) # Should be 0.7
#' # print(config_list$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff) # Should be 0.7
#' }
updateConfigParameter <- function(theObject,
                                function_name,
                                parameter_name,
                                new_value,
                                config_list_name = "config_list",
                                env = .GlobalEnv) {

  # --- Input Validation ---
  if (!isS4(theObject)) {
    stop("'theObject' must be an S4 object.")
  }
  if (!"args" %in% methods::slotNames(theObject)) {
      stop("'theObject' must have an '@args' slot.")
  }
  if (!is.character(function_name) || length(function_name) != 1) {
    stop("'function_name' must be a single character string.")
  }
  if (!is.character(parameter_name) || length(parameter_name) != 1) {
    stop("'parameter_name' must be a single character string.")
  }
  if (!exists(config_list_name, envir = env)) {
      stop("Global config list '", config_list_name, "' not found in the specified environment.")
  }

  # Retrieve the global list safely
  current_config_list <- get(config_list_name, envir = env)

  if (!is.list(current_config_list)) {
      stop("Global variable '", config_list_name, "' is not a list.")
  }
   if (!function_name %in% names(current_config_list)) {
    warning("Function name '", function_name, "' not found in global config list '", config_list_name, "'. Adding it.")
    current_config_list[[function_name]] <- list()
  }
  if (!parameter_name %in% names(current_config_list[[function_name]])) {
       warning("Parameter '", parameter_name, "' not found under '", function_name, "' in global config list '", config_list_name, "'. Adding it.")
  }


  # --- Update S4 Object @args ---
  if (is.null(theObject@args)) {
       theObject@args <- list() # Initialize args if it's NULL
   }
  if (!is.list(theObject@args[[function_name]])) {
     # Initialize the sub-list if it doesn't exist or isn't a list
     theObject@args[[function_name]] <- list()
  }
  theObject@args[[function_name]][[parameter_name]] <- new_value
  message("Updated @args$", function_name, "$", parameter_name, " in S4 object.")


  # --- Update Global Config List ---
  current_config_list[[function_name]][[parameter_name]] <- new_value
  # Assign the modified list back to the global environment
  assign(config_list_name, current_config_list, envir = env)
  message("Updated ", config_list_name, "$", function_name, "$", parameter_name, " in global environment.")

  # --- Return the modified object ---
  return(theObject)
}


# ----------------------------------------------------------------------------
# extract_experiment
# ----------------------------------------------------------------------------
##################################################################################################################
#' @title Extract Substrings from Underscore-Separated Strings
#'
#' @description
#' Extracts substrings from underscore-separated strings using different modes:
#' range (between positions), start (first element), or end (last element).
#'
#' @param x Character vector containing the strings to process
#' @param mode Character string specifying extraction mode:
#'   * "range": Extract elements between two underscore positions
#'   * "start": Extract from start to first underscore
#'   * "end": Extract from last underscore to end
#' @param start Integer: Starting position for range mode (1-based)
#' @param end Integer: Ending position for range mode (1-based, required for range mode)
#'
#' @return Character vector with extracted strings. Returns NA for strings where
#'   requested positions are out of bounds.
#'
#' @examples
#' x <- "20140602_ffs_expt1_r1_junk"
#' extract_experiment(x, mode = "range", start = 1, end = 3)  # "20140602_ffs_expt1"
#' extract_experiment(x, mode = "start")  # "20140602"
#' extract_experiment(x, mode = "end")    # "junk"
#'
#' # Multiple strings
#' x <- c("20140602_ffs_expt1_r1_junk", "20140603_ffs_expt2_r2_test")
#' extract_experiment(x, mode = "range", start = 2, end = 3)  # c("ffs_expt1", "ffs_expt2")
#'
#' @export
extract_experiment <- function(x, mode = "range", start = 1, end = NULL) {
  if (!mode %in% c("range", "start", "end")) {
    stop("Mode must be one of: 'range', 'start', 'end'")
  }

  process_string <- function(str) {
    parts <- unlist(strsplit(str, "_"))

    if (mode == "range") {
      if (is.null(end)) stop("End position required for range mode")
      if (start > length(parts) || end > length(parts)) {
        warning("Position out of bounds for string: ", str)

        return(NA_character_)
      }
      return(paste(parts[start:end], collapse = "_"))
    }

    else if (mode == "start") {
      return(parts[1])
    }

    else if (mode == "end") {
      return(parts[length(parts)])
    }
  }

  sapply(x, process_string)
}


# ----------------------------------------------------------------------------
# formatConfigList
# ----------------------------------------------------------------------------
#' @title Format Configuration List
#' @param config_list List of configuration parameters
#' @param indent Number of spaces for indentation
#' @export
formatConfigList <- function(config_list, indent = 0) {
    message(sprintf("--- DEBUG66: Entering formatConfigList with %d items, indent=%d ---", length(config_list), indent))
    output <- character()

    # Exclude internal_workflow_source_dir from printing
    names_to_process <- names(config_list)
    names_to_process <- names_to_process[names_to_process != "internal_workflow_source_dir"]
    message(sprintf("   DEBUG66: Processing %d items after exclusions", length(names_to_process)))

    # FUNCTIONAL APPROACH - NO FOR LOOPS that cause hanging - Works in Shiny AND .rmd
    output <- if (requireNamespace("purrr", quietly = TRUE)) {
        purrr::map_chr(names_to_process, function(name) {
        message(sprintf("   DEBUG66: Processing config item '%s'", name))
        value <- config_list[[name]]
        message(sprintf("   DEBUG66: Item '%s' class: %s", name, paste(class(value), collapse=", ")))
        
        # Skip core_utilisation, seqinr_obj, and complex objects from display
        if (name == "core_utilisation" || name == "seqinr_obj" ||
            any(class(value) %in% c("process", "R6", "multidplyr_cluster", "cluster", "SOCKcluster", "tbl_df", "tbl", "data.frame"))) {
            message(sprintf("   DEBUG66: Skipping '%s' due to complex class or large data frame", name))
            return("")  # Return empty string instead of next
        }

        # Format the name
        name_formatted <- gsub("\\.", " ", name)
        name_formatted <- gsub("_", " ", name_formatted)
        name_formatted <- tools::toTitleCase(name_formatted)

        # Handle different value types
        if (is.list(value)) {
            message(sprintf("   DEBUG66: '%s' is a list with %d elements", name, length(value)))
            if (length(value) > 0 && !is.null(names(value))) {
                output_lines <- paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ":")
                message(sprintf("   DEBUG66: Recursing into list '%s'", name))
                recursive_result <- formatConfigList(value, indent + 2)
                return(paste(c(output_lines, recursive_result), collapse = "\n"))
            } else if (length(value) > 0) { # Unnamed list, process elements functionally
                message(sprintf("   DEBUG66: '%s' is unnamed list, processing elements", name))
                header <- paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ":")
                
                                 # FUNCTIONAL processing of list elements - NO FOR LOOP
                 element_lines <- if (requireNamespace("purrr", quietly = TRUE)) {
                     purrr::imap_chr(value, function(item_val, item_idx) {
                         message(sprintf("   DEBUG66: Processing unnamed list item %d, class: %s", item_idx, paste(class(item_val), collapse=", ")))
                         if (is.atomic(item_val) && length(item_val) == 1) {
                             paste0(paste(rep(" ", indent + 2), collapse=""), "- ", as.character(item_val))
                         } else {
                             paste0(paste(rep(" ", indent + 2), collapse=""), "- [Complex List Element]")
                         }
                     })
                 } else {
                     # Base R fallback for list processing
                     sapply(seq_along(value), function(item_idx) {
                         item_val <- value[[item_idx]]
                         if (is.atomic(item_val) && length(item_val) == 1) {
                             paste0(paste(rep(" ", indent + 2), collapse=""), "- ", as.character(item_val))
                         } else {
                             paste0(paste(rep(" ", indent + 2), collapse=""), "- [Complex List Element]")
                         }
                     })
                 }
                return(paste(c(header, element_lines), collapse = "\n"))
            } else { # Empty list
                message(sprintf("   DEBUG66: '%s' is empty list", name))
                return(paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ": [Empty List]"))
            }
        } else {
            message(sprintf("   DEBUG66: '%s' is atomic/non-list, attempting conversion", name))
            
            # SAFE value display formatting
            value_display <- tryCatch({
                if (is.character(value) && length(value) > 5) {
                    paste(c(utils::head(value, 5), "..."), collapse = ", ")
                } else if (is.character(value) && length(value) > 1) {
                    paste(value, collapse = ", ")
                } else {
                    message(sprintf("   DEBUG66: About to convert '%s' to character", name))
                    result <- as.character(value)
                    message(sprintf("   DEBUG66: Conversion successful for '%s'", name))
                    if (length(result) == 0) "[Empty/NULL]" else result
                }
            }, error = function(e) {
                message(sprintf("   DEBUG66: Error converting '%s': %s", name, e$message))
                "[CONVERSION ERROR]"
            })
            
                         return(paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ": ", value_display))
         }
     }) |> 
     (\(x) x[x != ""])()  # Remove empty strings from skipped items
     } else {
         # Fallback using base R lapply if purrr not available
         result_list <- lapply(names_to_process, function(name) {
             value <- config_list[[name]]
             
             # Skip complex objects
             if (name == "core_utilisation" ||
                 any(class(value) %in% c("process", "R6", "multidplyr_cluster", "cluster", "SOCKcluster"))) {
                 return("")
             }
             
             # Format name
             name_formatted <- gsub("\\.", " ", name)
             name_formatted <- gsub("_", " ", name_formatted)
             name_formatted <- tools::toTitleCase(name_formatted)
             
             # Simple formatting for fallback
             value_display <- tryCatch({
                 if (is.list(value)) {
                     "[List Object]"
                 } else {
                     as.character(value)[1]  # Take first element only for safety
                 }
             }, error = function(e) "[CONVERSION ERROR]")
             
             paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ": ", value_display)
         })
         unlist(result_list[result_list != ""])
     }
     
     # Flatten the nested strings  
     output <- unlist(strsplit(output, "\n"))
    message(sprintf("--- DEBUG66: Exiting formatConfigList, returning %d lines ---", length(output)))
    return(output)
}


# ----------------------------------------------------------------------------
# updateMissingValueParameters
# ----------------------------------------------------------------------------
#' Update Missing Value Parameters in Configuration List and S4 Object
#'
#' @description
#' Automatically calculates and updates the missing value filtering parameters in the configuration list
#' and S4 object @args based on the experimental design matrix. The function ensures at least a specified 
#' number of groups have sufficient quantifiable values for analysis.
#'
#' @param theObject An S4 object containing the experimental data and design matrix.
#' @param min_reps_per_group Integer specifying the minimum number of replicates required in each passing group.
#'                          If a group has fewer total replicates than this value, the minimum is adjusted.
#' @param min_groups Integer specifying the minimum number of groups required to have sufficient
#'                  quantifiable values. Default is 2.
#' @param config_list_name The name of the global config list variable (defaults to "config_list").
#' @param env The environment where the global config list resides (defaults to .GlobalEnv).
#'
#' @return Updated S4 object with synchronized @args and global config_list
#'
#' @details
#' The function calculates:
#' - groupwise_percentage_cutoff: Based on minimum required replicates per group
#' - max_groups_percentage_cutoff: Based on minimum required groups
#' 
#' Both the S4 object's @args slot and the global config_list are updated to maintain synchronization.
#'
#' @examples
#' \dontrun{
#' protein_log2_quant_cln <- updateMissingValueParameters(
#'   theObject = protein_log2_quant_cln, 
#'   min_reps_per_group = 2, 
#'   min_groups = 2
#' )
#' }
#'
#' @export
updateMissingValueParameters <- function(theObject, 
                                       min_reps_per_group = 2, 
                                       min_groups = 2,
                                       config_list_name = "config_list",
                                       env = .GlobalEnv) {
    
    # --- Input Validation ---
    if (!isS4(theObject)) {
        stop("'theObject' must be an S4 object.")
    }
    if (!"design_matrix" %in% methods::slotNames(theObject)) {
        stop("'theObject' must have a '@design_matrix' slot.")
    }
    if (!"args" %in% methods::slotNames(theObject)) {
        stop("'theObject' must have an '@args' slot.")
    }
    if (!exists(config_list_name, envir = env)) {
        stop("Global config list '", config_list_name, "' not found in the specified environment.")
    }
    
    # Extract design matrix from S4 object
    design_matrix <- theObject@design_matrix
    
    # Retrieve the global config list
    config_list <- get(config_list_name, envir = env)
    
    # Get number of replicates per group
    reps_per_group_tbl <- design_matrix |>
        group_by(group) |>
        summarise(n_reps = n()) |>
        ungroup()
    
    # Get total number of groups
    total_groups <- nrow(reps_per_group_tbl)
    
    if (min_groups > total_groups) {
        stop("min_groups cannot be larger than total number of groups")
    }
    
    # Calculate percentage missing allowed for each group
    group_thresholds <- reps_per_group_tbl |>
        mutate(
            adjusted_min_reps = pmin(n_reps, min_reps_per_group),
            max_missing = n_reps - adjusted_min_reps,
            missing_percent = round((max_missing / n_reps) * 100, 3)
        )
    
    # Use a consistent percentage threshold across all groups
    # Take the maximum percentage to ensure all groups meet minimum requirements
    groupwise_cutoff <- max(group_thresholds$missing_percent)
    
    # Calculate maximum failing groups allowed
    max_failing_groups <- total_groups - min_groups
    max_groups_cutoff <- round((max_failing_groups / total_groups) * 100, 3)
    
    # Update global config_list
    if (is.null(config_list$removeRowsWithMissingValuesPercent)) {
        config_list$removeRowsWithMissingValuesPercent <- list()
    }
    config_list$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff <- groupwise_cutoff
    config_list$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff <- max_groups_cutoff
    config_list$removeRowsWithMissingValuesPercent$proteins_intensity_cutoff_percentile <- 1
    
    # Assign updated config_list back to global environment
    assign(config_list_name, config_list, envir = env)
    
    # Update S4 object @args to match config_list
    if (is.null(theObject@args)) {
        theObject@args <- list()
    }
    if (is.null(theObject@args$removeRowsWithMissingValuesPercent)) {
        theObject@args$removeRowsWithMissingValuesPercent <- list()
    }
    
    theObject@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff <- groupwise_cutoff
    theObject@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff <- max_groups_cutoff
    theObject@args$removeRowsWithMissingValuesPercent$proteins_intensity_cutoff_percentile <- 1
    
    # Create informative message
    basic_msg <- sprintf(
        "Updated missing value parameters in both config_list and S4 object @args:
    - Requiring at least %d replicates in each passing group (varies by group size)
    - Requiring at least %d out of %d groups to pass (%.3f%% failing groups allowed)
    - groupwise_percentage_cutoff set to %.3f%%
    - max_groups_percentage_cutoff set to %.3f%%",
        min_reps_per_group,
        min_groups,
        total_groups,
        max_groups_cutoff,
        groupwise_cutoff,
        max_groups_cutoff
    )
    
    # Add details for each group
    group_detail_strings <- group_thresholds |>
        mutate(
            detail = sprintf("    Group %s: %d out of %d replicates required (%.3f%% missing allowed)",
                             group, adjusted_min_reps, n_reps, missing_percent)
        ) |>
        dplyr::pull(detail)
        
    group_details <- paste(group_detail_strings, collapse = "\n")
    
    # Print the message
    message(paste(basic_msg, "\n\nGroup details:", group_details, sep = "\n"))
    message("✅ S4 object @args and global config_list are now synchronized")
    
    return(theObject)
}


# ----------------------------------------------------------------------------
# chooseBestProteinAccession_s3
# ----------------------------------------------------------------------------
#' Choose Best Protein Accession for Data Frame (S3 Version)
#'
#' @description A simplified S3 version of chooseBestProteinAccession that works
#' directly on data frames. Resolves protein groups by selecting the best
#' representative protein based on FASTA sequence data and aggregates
#' quantitative values.
#'
#' @param data_tbl Data frame containing protein quantitative data
#' @param protein_id_column Character string specifying the protein ID column name
#' @param seqinr_obj Data frame with FASTA sequence information (aa_seq_tbl_final)
#' @param seqinr_accession_column Character string specifying the accession column in seqinr_obj
#' @param delim Character string delimiter used to separate protein groups (default: ";")
#' @param replace_zero_with_na Logical, whether to replace zeros with NA (default: TRUE)
#' @param aggregation_method Character string aggregation method ("mean", "median", "sum") (default: "mean")
#'
#' @return Data frame with cleaned protein IDs and aggregated values
#'
#' @details
#' This function processes protein groups (e.g., "P12345;P67890;Q11111") by:
#' \itemize{
#'   \item Splitting groups using the specified delimiter
#'   \item Finding the best representative protein using FASTA sequence data
#'   \item Aggregating quantitative values for the chosen protein
#'   \item Returning cleaned data with single protein IDs
#' }
#'
#' The selection of "best" protein is based on:
#' \itemize{
#'   \item Presence in FASTA sequence database
#'   \item Protein length (longer sequences preferred)
#'   \item Alphabetical order as tiebreaker
#' }
#'
#' @examples
#' \dontrun{
#' cleaned_data <- chooseBestProteinAccession_s3(
#'   data_tbl = protein_data,
#'   protein_id_column = "Protein.Group",
#'   seqinr_obj = aa_seq_tbl_final,
#'   seqinr_accession_column = "uniprot_acc"
#' )
#' }
#'
#' @export
chooseBestProteinAccession_s3 <- function(data_tbl,
                                          protein_id_column,
                                          seqinr_obj,
                                          seqinr_accession_column = "uniprot_acc",
                                          delim = ";",
                                          replace_zero_with_na = TRUE,
                                          aggregation_method = "mean") {
  
  cat("=== Starting chooseBestProteinAccession_s3 ===\n")
  
  tryCatch({
    # Validate inputs with detailed error messages
    cat("*** S3 CLEANUP: Validating inputs ***\n")
    
    if (is.null(data_tbl) || !is.data.frame(data_tbl)) {
      stop("data_tbl must be a non-null data frame")
    }
    
    if (nrow(data_tbl) == 0) {
      stop("data_tbl cannot be empty")
    }
    
    cat(sprintf("*** S3 CLEANUP: data_tbl dimensions: %d rows x %d cols ***\n", nrow(data_tbl), ncol(data_tbl)))
    cat(sprintf("*** S3 CLEANUP: data_tbl column names: %s ***\n", paste(head(names(data_tbl), 10), collapse = ", ")))
    
    if (!protein_id_column %in% names(data_tbl)) {
      stop(sprintf("Protein ID column '%s' not found in data. Available columns: %s", 
                   protein_id_column, paste(names(data_tbl), collapse = ", ")))
    }
    
    if (is.null(seqinr_obj) || !is.data.frame(seqinr_obj)) {
      stop("seqinr_obj must be a non-null data frame")
    }
    
    if (nrow(seqinr_obj) == 0) {
      stop("seqinr_obj cannot be empty")
    }
    
    cat(sprintf("*** S3 CLEANUP: seqinr_obj dimensions: %d rows x %d cols ***\n", nrow(seqinr_obj), ncol(seqinr_obj)))
    cat(sprintf("*** S3 CLEANUP: seqinr_obj column names: %s ***\n", paste(head(names(seqinr_obj), 10), collapse = ", ")))
    
    if (!seqinr_accession_column %in% names(seqinr_obj)) {
      stop(sprintf("Accession column '%s' not found in sequence data. Available columns: %s", 
                   seqinr_accession_column, paste(names(seqinr_obj), collapse = ", ")))
    }
    
    cat("*** S3 CLEANUP: Input validation passed ***\n")
    
    # Get non-quantitative columns (metadata columns) with error handling
    cat("*** S3 CLEANUP: Identifying column types ***\n")
    
    numeric_cols <- tryCatch({
      sapply(data_tbl, is.numeric)
    }, error = function(e) {
      stop(sprintf("Error checking column types: %s", e$message))
    })
    
    quant_cols <- names(data_tbl)[numeric_cols]
    meta_cols <- names(data_tbl)[!numeric_cols]
    
    cat(sprintf("*** S3 CLEANUP: Found %d quantitative columns and %d metadata columns ***\n", 
                length(quant_cols), length(meta_cols)))
    cat(sprintf("*** S3 CLEANUP: Quantitative columns: %s ***\n", paste(head(quant_cols, 5), collapse = ", ")))
    cat(sprintf("*** S3 CLEANUP: Metadata columns: %s ***\n", paste(head(meta_cols, 5), collapse = ", ")))
  
    # Replace zeros with NA if requested
    if (replace_zero_with_na && length(quant_cols) > 0) {
      cat("*** S3 CLEANUP: Replacing zeros with NA in quantitative columns ***\n")
      
      tryCatch({
        data_tbl[quant_cols] <- lapply(data_tbl[quant_cols], function(x) {
          if (is.numeric(x)) {
            x[x == 0] <- NA
          }
          x
        })
        cat("*** S3 CLEANUP: Zero replacement completed ***\n")
      }, error = function(e) {
        stop(sprintf("Error replacing zeros with NA: %s", e$message))
      })
    }
    
    # Get unique protein groups
    cat("*** S3 CLEANUP: Getting unique protein groups ***\n")
    
    protein_groups <- tryCatch({
      unique(data_tbl[[protein_id_column]])
    }, error = function(e) {
      stop(sprintf("Error getting unique protein groups: %s", e$message))
    })
    
    cat(sprintf("*** S3 CLEANUP: Processing %d unique protein groups ***\n", length(protein_groups)))
    cat(sprintf("*** S3 CLEANUP: First 5 protein groups: %s ***\n", paste(head(protein_groups, 5), collapse = ", ")))
    
    # Create lookup table for sequence information
    cat("*** S3 CLEANUP: Creating sequence lookup table ***\n")
    
    seq_lookup <- tryCatch({
      if ("length" %in% names(seqinr_obj)) {
        cat("*** S3 CLEANUP: Using existing length column ***\n")
        seqinr_obj |>
          dplyr::select(dplyr::all_of(c(seqinr_accession_column, "length"))) |>
          dplyr::distinct()
      } else {
        cat("*** S3 CLEANUP: Creating default length column ***\n")
        # If no length column, create one or use a default
        seqinr_obj |>
          dplyr::select(dplyr::all_of(seqinr_accession_column)) |>
          dplyr::distinct() |>
          dplyr::mutate(length = 1000)  # Default length
      }
    }, error = function(e) {
      stop(sprintf("Error creating sequence lookup: %s", e$message))
    })
    
    cat(sprintf("*** S3 CLEANUP: Created sequence lookup with %d entries ***\n", nrow(seq_lookup)))
  
    # Function to choose best protein from a group
    cat("*** S3 CLEANUP: Defining chooseBestFromGroup function ***\n")
    
    chooseBestFromGroup <- function(group_string) {
      tryCatch({
        if (is.na(group_string) || group_string == "") {
          return(NA_character_)
        }
        
        # Split the group
        proteins <- stringr::str_split(group_string, delim)[[1]]
        proteins <- stringr::str_trim(proteins)  # Remove whitespace
        proteins <- proteins[proteins != ""]     # Remove empty strings
        
        if (length(proteins) == 1) {
          return(proteins[1])
        }
        
        # Find proteins that exist in sequence database
        proteins_in_db <- proteins[proteins %in% seq_lookup[[seqinr_accession_column]]]
        
        if (length(proteins_in_db) == 0) {
          # No proteins found in database, return first one
          return(proteins[1])
        }
        
        if (length(proteins_in_db) == 1) {
          return(proteins_in_db[1])
        }
        
        # Multiple proteins in database - choose by length (longest first)
        protein_info <- seq_lookup |>
          dplyr::filter(!!sym(seqinr_accession_column) %in% proteins_in_db) |>
          dplyr::arrange(desc(length), !!sym(seqinr_accession_column))
        
        return(protein_info[[seqinr_accession_column]][1])
        
      }, error = function(e) {
        warning(sprintf("Error processing protein group '%s': %s", group_string, e$message))
        return(group_string)  # Return original if processing fails
      })
    }
    
    # Create mapping of original groups to best proteins
    cat("*** S3 CLEANUP: Creating protein group mapping ***\n")
    
    protein_mapping <- tryCatch({
      data.frame(
        original_group = protein_groups,
        best_protein = sapply(protein_groups, chooseBestFromGroup, USE.NAMES = FALSE),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      stop(sprintf("Error creating protein mapping: %s", e$message))
    })
    
    # Count reductions
    groups_with_multiple <- tryCatch({
      sum(stringr::str_detect(protein_mapping$original_group, delim), na.rm = TRUE)
    }, error = function(e) {
      cat(sprintf("*** S3 CLEANUP: Warning - could not count multi-protein groups: %s ***\n", e$message))
      0
    })
    
    cat(sprintf("*** S3 CLEANUP: Found %d protein groups with multiple proteins ***\n", groups_with_multiple))
  
    # Add mapping to data
    cat("*** S3 CLEANUP: Applying protein mapping to data ***\n")
    
    data_with_mapping <- tryCatch({
      data_tbl |>
        dplyr::left_join(
          protein_mapping,
          by = setNames("original_group", protein_id_column)
        )
    }, error = function(e) {
      stop(sprintf("Error applying protein mapping: %s", e$message))
    })
    
    cat(sprintf("*** S3 CLEANUP: Data mapping completed. Mapped data has %d rows ***\n", nrow(data_with_mapping)))
    
    # Group by best protein and aggregate
    cat(sprintf("*** S3 CLEANUP: Aggregating data using method: %s ***\n", aggregation_method))
    
    # Define aggregation function
    agg_func <- switch(aggregation_method,
      "mean" = function(x) mean(x, na.rm = TRUE),
      "median" = function(x) median(x, na.rm = TRUE),
      "sum" = function(x) sum(x, na.rm = TRUE),
      function(x) mean(x, na.rm = TRUE)  # Default to mean
    )
    
    # Separate quantitative and metadata columns for different aggregation
    cat("*** S3 CLEANUP: Aggregating quantitative data ***\n")
    
    quant_data <- tryCatch({
      if (length(quant_cols) > 0) {
        data_with_mapping |>
          dplyr::select(best_protein, dplyr::all_of(quant_cols)) |>
          dplyr::group_by(best_protein) |>
          dplyr::summarise(
            dplyr::across(dplyr::all_of(quant_cols), agg_func),
            .groups = "drop"
          )
      } else {
        data.frame(best_protein = unique(data_with_mapping$best_protein))
      }
    }, error = function(e) {
      stop(sprintf("Error aggregating quantitative data: %s", e$message))
    })
    
    cat("*** S3 CLEANUP: Aggregating metadata ***\n")
    
    # For metadata columns, take the first non-NA value for each best protein
    meta_data <- tryCatch({
      if (length(meta_cols) > 0) {
        data_with_mapping |>
          dplyr::select(best_protein, dplyr::all_of(meta_cols)) |>
          dplyr::group_by(best_protein) |>
          dplyr::summarise(
            dplyr::across(dplyr::all_of(meta_cols), ~ dplyr::first(na.omit(.))[1]),
            .groups = "drop"
          )
      } else {
        data.frame(best_protein = unique(data_with_mapping$best_protein))
      }
    }, error = function(e) {
      stop(sprintf("Error aggregating metadata: %s", e$message))
    })
    
    cat("*** S3 CLEANUP: Combining results ***\n")
    
    # Combine results
    result_data <- tryCatch({
      quant_data |>
        dplyr::left_join(meta_data, by = "best_protein") |>
        dplyr::rename(!!sym(protein_id_column) := best_protein)
    }, error = function(e) {
      stop(sprintf("Error combining results: %s", e$message))
    })
    
    cat("*** S3 CLEANUP: Reordering columns ***\n")
    
    # Reorder columns to match original
    result_data <- tryCatch({
      result_data |>
        dplyr::select(dplyr::all_of(names(data_tbl)))
    }, error = function(e) {
      stop(sprintf("Error reordering columns: %s", e$message))
    })
    
    # Final statistics
    original_proteins <- length(unique(data_tbl[[protein_id_column]]))
    final_proteins <- length(unique(result_data[[protein_id_column]]))
    reduction_pct <- round(((original_proteins - final_proteins) / original_proteins) * 100, 1)
    
    cat(sprintf("*** S3 CLEANUP: Protein cleanup completed: %d -> %d proteins (%.1f%% reduction) ***\n",
                     original_proteins, final_proteins, reduction_pct))
    
    return(result_data)
    
  }, error = function(e) {
    cat(sprintf("*** S3 CLEANUP FATAL ERROR: %s ***\n", e$message))
    cat("*** S3 CLEANUP: Returning original data ***\n")
    return(data_tbl)  # Return original data if everything fails
  })
}


# ----------------------------------------------------------------------------
# calcHtSize
# ----------------------------------------------------------------------------
#' @export
calcHtSize = function(ht, unit = "inch") {
  pdf(NULL)
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()

  c(w, h)
}


# ============================================================================
# MEMORY DEBUGGING UTILITIES
# ============================================================================
# These functions track ACTUAL process memory (what Task Manager shows),
# not just R's internal heap reporting which misses environment capture.
# ============================================================================

#' Get Actual Process Memory (Task Manager View)
#' 
#' Returns the actual memory used by the R process as reported by the OS.
#' This is what you see in Task Manager / Activity Monitor, and includes:
#' - R heap memory
#' - Environment captures in ggplot/tidyverse objects
#' - Lazy evaluation promises
#' - Memory-mapped files
#' 
#' @return Memory usage in MB (numeric)
#' @export
getProcessMemoryMB <- function() {
  tryCatch({
    if (.Platform$OS.type == "windows") {
      # Windows: Use tasklist command (more reliable than wmic on modern Windows)
      pid <- Sys.getpid()
      cmd_output <- system(sprintf('tasklist /FI "PID eq %d" /FO CSV /NH', pid), intern = TRUE)
      # Parse CSV output: "process.exe","PID","Session","Session#","Mem Usage"
      # Mem Usage format: "1,234,567 K"
      if (length(cmd_output) > 0 && nchar(cmd_output[1]) > 0) {
        mem_str <- gsub('"', '', strsplit(cmd_output[1], ',')[[1]][5])
        mem_kb <- as.numeric(gsub('[^0-9]', '', mem_str))
        return(mem_kb / 1024)  # Convert KB to MB
      }
    } else {
      # Unix/Mac: Use ps command
      pid <- Sys.getpid()
      mem_kb <- as.numeric(system(sprintf("ps -o rss= -p %d", pid), intern = TRUE))
      return(mem_kb / 1024)  # Convert KB to MB
    }
    return(NA_real_)
  }, error = function(e) {
    return(NA_real_)
  })
}

#' Get R Heap Memory
#' 
#' Returns R's internal heap memory (what pryr::mem_used() shows).
#' NOTE: This does NOT include environment captures!
#' 
#' @return Memory usage in MB (numeric)
#' @export
getRHeapMemoryMB <- function() {
  sum(gc()[, 2])
}

#' Comprehensive Memory Checkpoint
#' 
#' Logs BOTH R heap memory AND actual process memory to help identify

#' environment capture issues where R heap stays flat but process memory spikes.
#' 
#' @param step_name Character string describing the current step
#' @param context Optional context string (e.g., function name)
#' @param log_level One of "DEBUG66", "INFO", "WARN"
#' @return Invisibly returns a list with r_heap_mb and process_mb
#' @export
checkMemoryBoth <- function(step_name, context = "", log_level = "DEBUG66") {
  # Disabled for performance to avoid slow system calls
  invisible(list(r_heap_mb = 0, process_mb = 0, hidden_mb = 0))
}

#' Memory Delta Reporter
#' 
#' Compares current memory to a baseline and reports the change.
#' Useful for bracketing operations to see exactly how much memory they consume.
#' 
#' @param baseline_mem List from a previous checkMemoryBoth() call
#' @param step_name Character string describing the completed step
#' @param context Optional context string
#' @return Invisibly returns the delta values
#' @export
reportMemoryDelta <- function(baseline_mem, step_name, context = "") {
  # Disabled for performance to avoid slow system calls
  invisible(list(delta_r_heap = 0, delta_process = 0))
}
