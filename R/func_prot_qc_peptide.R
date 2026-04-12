# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ============================================================================
# func_prot_qc_peptide.R
# ============================================================================
# Purpose: Peptide-level quality control and filtering functions
# 
# This file contains functions for peptide-level QC filtering, including
# intensity filtering, missing value filtering, replicate filtering, and
# q-value filtering. Functions in this file are used by mod_prot_qc_peptide.R
# and related QC modules.
#
# NOTE: Peptides are part of the proteomics workflow, hence "prot" prefix.
# This file contains peptide-specific QC functions within the proteomics context.
#
# Functions to extract here:
# - peptideIntensityFiltering(): S4 method for peptide intensity filtering
# - peptideIntensityFilteringHelper(): Helper for peptide intensity filtering
# - removePeptidesWithMissingValuesPercent(): S4 method for missing value filtering
# - removePeptidesWithMissingValuesPercentHelper(): Helper for missing value filtering
# - removePeptidesWithOnlyOneReplicate(): S4 method for replicate filtering
# - removePeptidesWithOnlyOneReplicateHelper(): Helper for replicate filtering
# - filterMinNumPeptidesPerProtein(): S4 method for filtering by peptides per protein
# - filterMinNumPeptidesPerSample(): S4 method for filtering by peptides per sample
# - srlQvalueProteotypicPeptideClean(): S4 method for q-value filtering
# - Additional peptide QC helper functions
#
# Dependencies:
# - dplyr, tidyr
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: peptideIntensityFiltering()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Filters peptides based on intensity thresholds
# setMethod(f = "peptideIntensityFiltering", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 2: peptideIntensityFilteringHelper()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Helper function for peptide intensity filtering
# peptideIntensityFilteringHelper <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 3: removePeptidesWithMissingValuesPercent()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Removes peptides with excessive missing values
# setMethod(f = "removePeptidesWithMissingValuesPercent", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 4: removePeptidesWithMissingValuesPercentHelper()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Helper for removing peptides with missing values
# removePeptidesWithMissingValuesPercentHelper <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 5: removePeptidesWithOnlyOneReplicate()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Removes peptides present in only one replicate
# setMethod(f = "removePeptidesWithOnlyOneReplicate", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 6: removePeptidesWithOnlyOneReplicateHelper()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Helper for removing peptides with only one replicate
# removePeptidesWithOnlyOneReplicateHelper <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 7: filterMinNumPeptidesPerProtein()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Filters proteins with insufficient peptides
# setMethod(f = "filterMinNumPeptidesPerProtein", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 8: filterMinNumPeptidesPerProteinHelper()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Helper for filtering by peptides per protein
# filterMinNumPeptidesPerProteinHelper <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 9: filterMinNumPeptidesPerSample()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Filters samples with insufficient peptides
# setMethod(f = "filterMinNumPeptidesPerSample", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 10: filterMinNumPeptidesPerSampleHelper()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Helper for filtering by peptides per sample
# filterMinNumPeptidesPerSampleHelper <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 11: srlQvalueProteotypicPeptideClean()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Filters peptides based on q-value and proteotypic status
# setMethod(f = "srlQvalueProteotypicPeptideClean", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 12: srlQvalueProteotypicPeptideCleanHelper()
# Current location: R/peptideVsSamplesS4Objects.R
# Description: Helper for q-value and proteotypic filtering
# srlQvalueProteotypicPeptideCleanHelper <- function(...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 13: checkPeptideNAPercentages()
# Current location: R/da_proteins_functions.R
# Description: Checks peptide NA percentages
# checkPeptideNAPercentages <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 14: removePeptidesOnlyInHek293()
# Current location: R/da_proteins_functions.R
# Description: Removes peptides only found in HEK293 cells
# removePeptidesOnlyInHek293 <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 15: removePeptidesWithoutAbundances()
# Current location: R/da_proteins_functions.R
# Description: Removes peptides without abundance values
# removePeptidesWithoutAbundances <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 16: filterByScoreAndGetSimilarPeptides()
# Current location: R/enrichment_functions.R
# Description: Filters peptides by score and finds similar peptides
# filterByScoreAndGetSimilarPeptides <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 17: filterPeptideAndExtractProbabilities()
# Current location: R/enrichment_functions.R
# Description: Filters peptides and extracts probabilities
# filterPeptideAndExtractProbabilities <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 18: groupParalogPeptides()
# Current location: R/enrichment_functions.R
# Description: Groups paralog peptides
# groupParalogPeptides <- function(...) {
#   # Extract from R/enrichment_functions.R
# }



# ----------------------------------------------------------------------------
# resolvePeptideQcColumnName
# ----------------------------------------------------------------------------
#' @keywords internal
resolvePeptideQcColumnName <- function(column_expr, env = parent.frame()) {
  if (is.null(column_expr)) {
    return(NULL)
  }

  if (is.character(column_expr) && length(column_expr) == 1) {
    return(column_expr)
  }

  column_value <- tryCatch(
    eval(column_expr, envir = env),
    error = function(e) NULL
  )
  if (is.character(column_value) && length(column_value) == 1) {
    return(column_value)
  }
  if (is.symbol(column_value)) {
    return(as.character(column_value))
  }

  if (is.symbol(column_expr)) {
    return(as.character(column_expr))
  }

  rlang::expr_text(column_expr)
}

































