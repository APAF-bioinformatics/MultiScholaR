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
# func_prot_s4_objects.R
# ============================================================================
# Purpose: Proteomics S4 class definitions and methods
#
# This file contains S4 class definitions and methods specific to proteomics,
# including PeptideQuantitativeData, ProteinQuantitativeData classes,
# their constructors, and proteomics-specific S4 methods.
#
# Consolidated from:
# - peptideVsSamplesS4Objects.R (PeptideQuantitativeData class + 23 methods)
# - proteinVsSamplesS4Objects.R (ProteinQuantitativeData method)
# - protein_da_analysis_wrapper.R (DE analysis methods)
#
# Dependencies:
# - methods package
# - func_general_s4_generics.R (for generic definitions)
# ============================================================================

# ==========================================
# ProteinQuantitativeData S4 Class
# ==========================================



# Initialize method to ensure design_matrix's sample_id column is character


# ==========================================
# Content from proteinVsSamplesS4Objects.R
# ==========================================
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


# ==========================================
# Content from protein_da_analysis_wrapper.R
# ==========================================
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Differential Expression Analysis for Proteomics
#'
#' This file contains modular functions for proteomics differential expression analysis,
#' breaking up the monolithic deAnalysisWrapperFunction into focused components.
#'
#' Functions:
#' - differentialAbundanceAnalysis: Main S4 generic for DE analysis
#' - differentialAbundanceAnalysisHelper: Core limma-based analysis
#' - generateProtDAVolcanoPlotGlimma: Interactive volcano plot generation
#' - generateDEHeatmap: Heatmap visualization with clustering
#'
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Note: generateProtDAVolcanoPlotGlimma has been moved to R/func_prot_da.R

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# NOTE: generateDEHeatmap is defined in func_prot_da.R
# This duplicate has been removed to avoid conflicts.
# The canonical version in func_prot_da.R includes:
# - circlize::colorRamp2 for proper color scaling
# - ComplexHeatmap::HeatmapAnnotation for group annotations
# - logger calls instead of message()
# - tryCatch error handling
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


# ==========================================
# Missing S4 Methods recovered from GUI branch
# ==========================================




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
setGeneric(
  name = "removeProteinsWithOnlyOneReplicate",
  def = function(theObject, core_utilisation = NULL, grouping_variable = NULL) {
    standardGeneric("removeProteinsWithOnlyOneReplicate")
  },
  signature = c("theObject", "core_utilisation", "grouping_variable")
)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
setGeneric(
  name = "removeRowsWithMissingValuesPercent",
  def = function(
    theObject,
    ruv_grouping_variable = NULL,
    groupwise_percentage_cutoff = NULL,
    max_groups_percentage_cutoff = NULL,
    proteins_intensity_cutoff_percentile = NULL
  ) {
    standardGeneric("removeRowsWithMissingValuesPercent")
  },
  signature = c(
    "theObject",
    "ruv_grouping_variable",
    "groupwise_percentage_cutoff",
    "max_groups_percentage_cutoff",
    "proteins_intensity_cutoff_percentile"
  )
)


# NOTE: setGeneric for chooseBestProteinAccession is defined in allGenerics.R
# Do NOT redefine it here as it invalidates methods defined in func_pept_s4_objects.R


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#' @export
# NOTE: setGeneric for chooseBestProteinAccessionSumDuplicates is defined in allGenerics.R
# The signature in allGenerics.R needs to match this method's parameters



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @export
setGeneric(
  name = "filterSamplesByProteinCorrelationThreshold",
  def = function(theObject, pearson_correlation_per_pair = NULL, min_pearson_correlation_threshold = NULL) {
    standardGeneric("filterSamplesByProteinCorrelationThreshold")
  },
  signature = c("theObject", "pearson_correlation_per_pair", "min_pearson_correlation_threshold")
)


# Format the design matrix so that only metadata for samples in the protein data are retained, and also
# sort the sample IDs in the same order as the data matrix

## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## MISSING METHODS COPIED FROM proteins4.R - Added during golem audit fix
## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------





## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------






#' @export
setGeneric("InitialiseGrid", function(dummy = NULL) {
  standardGeneric("InitialiseGrid")
})



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create a QC composite figure

#' @export
setGeneric(
  name = "createGridQC",
  def = function(theObject, pca_titles = NULL, density_titles = NULL, rle_titles = NULL, pearson_titles = NULL, cancor_titles = NULL, limpa_titles = NULL, ncol = 3, save_path = NULL, file_name = "pca_density_rle_pearson_corr_plots_merged", workflow_name = NULL) {
    standardGeneric("createGridQC")
  },
  signature = c("theObject")
)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @export
setGeneric(
  name = "getNegCtrlProtAnova",
  def = function(
    theObject,
    ruv_grouping_variable = NULL,
    percentage_as_neg_ctrl = NULL,
    num_neg_ctrl = NULL,
    ruv_qval_cutoff = NULL,
    ruv_fdr_method = NULL,
    exclude_pool_samples = TRUE
  ) {
    standardGeneric("getNegCtrlProtAnova")
  },
  signature = c("theObject", "ruv_grouping_variable", "num_neg_ctrl", "ruv_qval_cutoff", "ruv_fdr_method")
)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


#' @export
setGeneric(
  name = "getRuvIIIReplicateMatrix",
  def = function(theObject, ruv_grouping_variable = NULL) {
    standardGeneric("getRuvIIIReplicateMatrix")
  },
  signature = c("theObject", "ruv_grouping_variable")
)




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @export
setGeneric(
  name = "preservePeptideNaValues",
  def = function(peptide_obj, protein_obj) {
    standardGeneric("preservePeptideNaValues")
  },
  signature = c("peptide_obj", "protein_obj")
)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------


