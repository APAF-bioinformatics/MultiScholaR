# ============================================================================
# func_metab_de.R
# ============================================================================
# Purpose: Metabolomics differential abundance analysis functions
# 
# This file contains functions for metabolomics differential abundance
# analysis, including limma-based analysis and result formatting. Functions
# in this file are used by metabolomics DE modules and related workflows.
#
# Functions to extract here:
# - differentialAbundanceAnalysis(): S4 method for DE analysis (metabolite)
# - differentialAbundanceAnalysisHelper(): Helper for DE analysis
# - getCountsTable(): Gets counts table from object
# - Additional metabolomics DE helper functions
#
# Dependencies:
# - limma, edgeR
# - func_general_plotting.R (for visualization)
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: differentialAbundanceAnalysis() (metabolite method)
# Current location: R/metaboliteVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Performs differential abundance analysis on metabolites
# setMethod(f = "differentialAbundanceAnalysis", signature = "MetaboliteAssayData", ...) {
#   # Extract from R/metaboliteVsSamplesS4Objects.R
# }

# Function 2: differentialAbundanceAnalysisHelper()
# Current location: R/metaboliteVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Helper function for metabolite DE analysis
# setMethod(f = "differentialAbundanceAnalysisHelper", ...) {
#   # Extract from R/metaboliteVsSamplesS4Objects.R
# }

# Function 3: getCountsTable()
# Current location: R/metabolite_de_analysis_wrapper.R, R/metabolites_de_analysis_wrapper.R
# Description: Gets counts table from object
# getCountsTable <- function(obj) {
#   # Extract from R/metabolite_de_analysis_wrapper.R or R/metabolites_de_analysis_wrapper.R
# }


# ----------------------------------------------------------------------------
# getCountsTable
# ----------------------------------------------------------------------------
# Helper function to get counts table
getCountsTable <- function(obj) {
  if (inherits(obj, "MetaboliteAssayData")) {
    message(sprintf("   Getting counts table for object of class: %s", class(obj)[1]))
    message(sprintf("   Returning metabolite_data with dimensions: %d rows, %d cols",
                    nrow(obj@metabolite_data), ncol(obj@metabolite_data)))
    obj@metabolite_data
  } else if (inherits(obj, "ProteinQuantitativeData")) {
    message(sprintf("   Returning protein_quant_table with dimensions: %d rows, %d cols",
                    nrow(obj@protein_quant_table), ncol(obj@protein_quant_table)))
    obj@protein_quant_table
  } else {
    message(sprintf("   ERROR: Unsupported object type: %s", class(obj)[1]))
    stop("Unsupported object type")
  }
}

