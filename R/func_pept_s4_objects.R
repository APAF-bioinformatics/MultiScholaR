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



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------




##----------------------------------------------------------------------------------------------------------------------------------------------------------------------



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#'@export







##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Find the Best Negative Control Percentage for RUV-III Analysis on Peptide Data
#'
#' This function automatically determines the optimal percentage of peptides to use
#' as negative controls for RUV-III analysis by testing different percentages and
#' evaluating the separation quality between "All" and "Control" groups in canonical
#' correlation plots.
#'
#' @param peptide_matrix_obj A PeptideQuantitativeData object containing
#'   the peptide quantification data
#' @param percentage_range A numeric vector specifying the range of percentages to test.
#'   Default is seq(1, 20, by = 1) for testing 1% to 20% in 1% increments
#' @param num_components_to_impute Number of components to use for imputation in ruvCancor.
#'   Default is 5
#' @param ruv_grouping_variable The grouping variable to use for RUV analysis.
#'   Default is "group"
#' @param ruv_qval_cutoff The FDR threshold for negative control selection.
#'   Default is 0.05
#' @param ruv_fdr_method The FDR calculation method. Default is "qvalue"
#' @param separation_metric The metric to use for evaluating separation quality.
#'   Options: "max_difference" (default), "mean_difference", "auc".
#'   "weighted_difference" is deprecated; use "max_difference" or "mean_difference" instead.
#' @param k_penalty_weight Weight for penalizing high k values in composite score.
#'   Default is 0.5. Higher values penalize high k more strongly
#' @param max_acceptable_k Maximum acceptable k value. k values above this get heavy penalty.
#'   Default is 3
#' @param adaptive_k_penalty Whether to automatically adjust max_acceptable_k based on sample size.
#'   Default is TRUE (recommended). Set to FALSE only if you need exact reproducibility with previous results
#' @param verbose Whether to print progress messages. Default is TRUE
#' @param ensure_matrix Whether to ensure peptide matrix is calculated. Default is TRUE
#'
#' @return A list containing:
#'   \itemize{
#'     \item best_percentage: The optimal percentage as a numeric value
#'     \item best_k: The optimal k value from findBestKElbow() for the best percentage
#'     \item best_control_genes_index: The control genes index for the best percentage
#'     \item best_separation_score: The separation score for the best percentage
#'     \item best_composite_score: The composite score (separation penalized by k value)
#'     \item optimization_results: A data frame with all tested percentages and their scores
#'     \item best_cancor_plot: The canonical correlation plot for the best percentage
#'     \item separation_metric_used: The separation metric that was used
#'     \item k_penalty_weight: The k penalty weight that was used
#'     \item max_acceptable_k: The maximum acceptable k value that was used
#'   }
#'
#' @importFrom logger log_info log_warn
#' @importFrom purrr imap map_dfr map_dbl
#' @export





##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' Peptide Missing Value Imputation using limpa Package
#'
#' This function uses the limpa package's detection probability curve (DPC) approach
#' for sophisticated missing value imputation specifically designed for proteomics data.
#' This method is more robust than traditional imputation as it models the missing value
#' mechanism based on detection probabilities.
#'
#' @param theObject A PeptideQuantitativeData object
#' @param imputed_value_column Name for the new column containing imputed values.
#'   Default is "Peptide.Imputed.Limpa"
#' @param use_log2_transform Whether to log2 transform the data before imputation.
#'   Default is TRUE (recommended by limpa)
#' @param verbose Whether to print progress messages. Default is TRUE
#' @param ensure_matrix Whether to ensure peptide matrix is calculated. Default is TRUE
#'
#' @details
#' The limpa package uses a detection probability curve (DPC) to model the relationship
#' between peptide intensity and the probability of detection. This allows for more
#' sophisticated imputation that accounts for the intensity-dependent nature of missing
#' values in proteomics data, rather than assuming they are missing at random.
#'
#' The process follows these steps:
#' 1. Estimate the detection probability curve using dpc()
#' 2. Perform row-wise imputation using dpcImpute()
#' 3. Transform results back to original scale if needed
#'
#' @return Updated PeptideQuantitativeData object with imputed values
#'
#' @export


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plotting Methods for PeptideQuantitativeData
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Set up S4 class definitions for ggplot objects
setOldClass(c("gg", "ggplot"))
setOldClass("ggplot2::ggplot")

#' Create Density Plots from PCA data
#'
#' This function takes a ggplot object (presumably a PCA plot) and creates
#' density plots for the first two principal components.



##----------------------------------------------------------------------------------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Normalization Methods for PeptideQuantitativeData
##----------------------------------------------------------------------------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------------------------------------------------------------------------------




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Pearson correlation plot for PeptideQuantitativeData


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pearsonCorForSamplePairs for PeptideQuantitativeData

