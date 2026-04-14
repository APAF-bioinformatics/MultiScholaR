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
# func_lipid_da.R
# ============================================================================
# Purpose: Lipidomics differential abundance analysis functions
#
# This file contains functions for lipidomics differential abundance
# analysis, including limma-based analysis and result formatting. Functions
# in this file are used by lipidomics DA modules and related workflows.
#
# Functions to extract here:
# - differentialAbundanceAnalysis(): S4 method for DA analysis (lipid)
# - differentialAbundanceAnalysisHelper(): Helper for DA analysis
# - getCountsTable(): Gets counts table from object
# - Additional lipidomics DE helper functions
#
# Dependencies:
# - limma, edgeR
# - func_general_plotting.R (for visualization)
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: differentialAbundanceAnalysis() (lipid method)
# Current location: R/lipidVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Performs differential abundance analysis on lipids
# setMethod(f = "differentialAbundanceAnalysis", signature = "LipidomicsAssayData", ...) {
#   # Extract from R/lipidVsSamplesS4Objects.R
# }

# Function 2: differentialAbundanceAnalysisHelper()
# Current location: R/lipidVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Helper function for lipid DA analysis
# setMethod(f = "differentialAbundanceAnalysisHelper", ...) {
#   # Extract from R/lipidVsSamplesS4Objects.R
# }

# Function 3: getCountsTable()
# Current location: R/lipid_da_analysis_wrapper.R, R/lipids_da_analysis_wrapper.R
# Description: Gets counts table from object
# getCountsTable <- function(obj) {
#   # Extract from R/lipid_da_analysis_wrapper.R or R/lipids_da_analysis_wrapper.R
# }


# ----------------------------------------------------------------------------
# getCountsTable
# ----------------------------------------------------------------------------
# Helper function to get counts table
getCountsTable <- function(obj) {
    if (inherits(obj, "LipidomicsAssayData")) {
        message(sprintf("   Getting counts table for object of class: %s", class(obj)[1]))
        message(sprintf(
            "   Returning lipid_data with dimensions: %d rows, %d cols",
            nrow(obj@lipid_data), ncol(obj@lipid_data)
        ))
        obj@lipid_data
    } else if (inherits(obj, "ProteinQuantitativeData")) {
        message(sprintf(
            "   Returning protein_quant_table with dimensions: %d rows, %d cols",
            nrow(obj@protein_quant_table), ncol(obj@protein_quant_table)
        ))
        obj@protein_quant_table
    } else {
        message(sprintf("   ERROR: Unsupported object type: %s", class(obj)[1]))
        stop("Unsupported object type")
    }
}


# ============================================================================
# METABOLOMICS DIFFERENTIAL EXPRESSION ANALYSIS FUNCTIONS
# ============================================================================
# These functions provide differential expression analysis for lipidomics
# data stored as multiple assays (e.g., LCMS_Pos, LCMS_Neg) in the
# LipidomicsAssayData S4 object.
# ============================================================================








# ----------------------------------------------------------------------------
# getLipidQuantData
# ----------------------------------------------------------------------------
#' Extract quantitative data columns from an assay data frame
#'
#' @description Helper function to separate lipid quantitative data
#'   (sample columns) from metadata columns (ID, annotation, etc.).
#'
#' @param assay_df Data frame containing lipid data for one assay.
#' @param lipid_id_col Name of the lipid ID column.
#' @param annotation_col Name of the annotation column.
#' @param additional_meta_cols Additional columns to exclude from quant data.
#'
#' @return A list with:
#'   - quant_data: Data frame with only sample columns
#'   - meta_data: Data frame with only metadata columns
#'   - sample_cols: Names of sample columns
#'
#' @export
getLipidQuantData <- function(
  assay_df,
  lipid_id_col = "Alignment ID",
  annotation_col = "Lipid name",
  additional_meta_cols = NULL
) {
    # Identify metadata columns to exclude
    meta_cols <- c(lipid_id_col, annotation_col)
    if (!is.null(additional_meta_cols)) {
        meta_cols <- c(meta_cols, additional_meta_cols)
    }

    # Get sample columns (everything that's not metadata)
    all_cols <- colnames(assay_df)
    sample_cols <- setdiff(all_cols, meta_cols)

    # Also exclude any obviously non-numeric columns
    sample_cols <- sample_cols[sapply(assay_df[, sample_cols, drop = FALSE], is.numeric)]

    list(
        quant_data = assay_df[, sample_cols, drop = FALSE],
        meta_data = assay_df[, intersect(meta_cols, all_cols), drop = FALSE],
        sample_cols = sample_cols
    )
}








