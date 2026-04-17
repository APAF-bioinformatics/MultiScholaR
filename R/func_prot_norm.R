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
# func_prot_norm.R
# ============================================================================
# Purpose: Protein normalization functions
# 
# This file contains functions for protein-level normalization, including
# RUV-III-C normalization, negative control selection, and normalization
# parameter optimization. Functions in this file are used by mod_prot_norm.R
# and related normalization modules.
#
# Functions to extract here:
# - normaliseBetweenSamples(): S4 method for between-sample normalization
# - normaliseUntransformedData(): S4 method for untransformed data normalization
# - ruvIII_C_Varying(): S4 method for RUV-III-C normalization (protein)
# - ruvCancor(): S4 method for RUV canonical correlation analysis
# - ruvCancorFast(): S4 method for fast RUV canonical correlation
# - findBestNegCtrlPercentage(): Find optimal negative control percentage
# - findBestK(): Find optimal k value for RUV
# - getNegCtrlProtAnova(): Get negative control proteins using ANOVA
# - Additional normalization helper functions
#
# Dependencies:
# - limma, RUVSeq (or custom RUV implementation)
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: normaliseBetweenSamples()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Normalizes data between samples using various methods
# setMethod(f = "normaliseBetweenSamples", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 2: normaliseUntransformedData()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Normalizes untransformed data
# setMethod(f = "normaliseUntransformedData", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 3: ruvIII_C_Varying() (protein method)
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Applies RUV-III-C normalization with varying k
# setMethod(f = "ruvIII_C_Varying", signature = "ProteinQuantitativeData", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 4: ruvCancor()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Performs RUV canonical correlation analysis
# setMethod(f = "ruvCancor", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 5: ruvCancorFast()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Fast version of RUV canonical correlation analysis
# setMethod(f = "ruvCancorFast", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 6: findBestNegCtrlPercentage()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Finds optimal percentage of negative control proteins
# findBestNegCtrlPercentage <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 7: findBestK()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Finds optimal k value for RUV normalization
# findBestK <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 8: findBestKForAssayList()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Finds best k for a list of assays
# findBestKForAssayList <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 9: getNegCtrlProtAnova()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Gets negative control proteins using ANOVA
# setMethod(f = "getNegCtrlProtAnova", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 10: getNegCtrlProtAnovaHelper()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Helper for getting negative control proteins
# getNegCtrlProtAnovaHelper <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 11: getRuvIIIReplicateMatrixHelper()
# Current location: R/proteinVsSamplesS4Objects.R
# Description: Helper for creating RUV-III replicate matrix
# getRuvIIIReplicateMatrixHelper <- function(...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 12: extractRuvResults()
# Current location: R/helper_functions.R
# Description: Extracts RUV normalization results
# extractRuvResults <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 13: updateRuvParameters()
# Current location: R/helper_functions.R
# Description: Updates RUV parameters in config
# updateRuvParameters <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 14: scaleCenterAndFillMissing()
# Current location: R/helper_functions.R
# Description: Scales, centers, and fills missing values
# scaleCenterAndFillMissing <- function(...) {
#   # Extract from R/helper_functions.R
# }























