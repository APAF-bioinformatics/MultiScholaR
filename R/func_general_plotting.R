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
# func_general_plotting.R
# ============================================================================
# Purpose: General plotting utilities and visualization functions
#
# This file contains general-purpose plotting functions used across all
# omics types, including volcano plots, PCA plots, RLE plots, correlation
# plots, enrichment plots, and various QC visualizations. Functions in
# this file are used by multiple modules across proteomics, metabolomics,
# and multiomics workflows.
#
# Functions to extract here:
# - Volcano plot functions (plotVolcano, plotOneVolcano, etc.)
# - PCA plot functions (plotPca, plotPcaHelper, etc.)
# - RLE plot functions (plotRle, plotRleHelper, etc.)
# - Density plot functions (plotDensity, etc.)
# - Correlation plot functions (plotPearson, etc.)
# - Enrichment plot functions (plotEnrichmentBarplot, etc.)
# - Heatmap functions (getProteinsHeatMap, etc.)
# - Interactive plot functions (writeInteractiveVolcanoPlotProteomics, etc.)
# - Plot saving functions (savePlot, save_plot, etc.)
# - Theme and color palette functions (apafTheme, get_color_palette, etc.)
#
# Dependencies:
# - ggplot2, plotly, Glimma, gridExtra, patchwork
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# === Volcano Plot Functions ===

# Function 1: plotVolcano()
# Current location: R/metaboliteVsSamplesS4Objects.R
# Type: S4 generic (exportMethods)
# Description: Creates volcano plots
# setGeneric(name="plotVolcano", ...) {
#   # Extract from R/metaboliteVsSamplesS4Objects.R
# }

# Function 2: plotOneVolcano()
# Current location: R/da_proteins_functions.R
# Description: Creates a single volcano plot
# plotOneVolcano <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 3: plotOneVolcanoNoVerticalLines()
# Current location: R/da_proteins_functions.R
# Description: Creates volcano plot without vertical lines
# plotOneVolcanoNoVerticalLines <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 4: plotVolcanoS4()
# Current location: R/metaboliteVsSamplesS4Objects.R
# Type: S4 generic (exportMethods)
# Description: Creates volcano plots for S4 objects
# setGeneric(name="plotVolcanoS4", ...) {
#   # Extract from R/metaboliteVsSamplesS4Objects.R
# }

# Function 5: plotInteractiveVolcano()
# Current location: R/metaboliteVsSamplesS4Objects.R
# Type: S4 generic (exportMethods)
# Description: Creates interactive volcano plots
# setGeneric(name="plotInteractiveVolcano", ...) {
#   # Extract from R/metaboliteVsSamplesS4Objects.R
# }

# Function 6: prepareDataForVolcanoPlot()
# Current location: R/da_proteins_functions.R
# Description: Prepares data for volcano plot
# prepareDataForVolcanoPlot <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 7: writeInteractiveVolcanoPlotProteomics()
# Current location: R/da_proteins_functions.R
# Description: Writes interactive volcano plot for proteomics
# writeInteractiveVolcanoPlotProteomics <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 8: writeInteractiveVolcanoPlotProteomicsMain()
# Current location: R/da_proteins_functions.R
# Description: Main function for interactive volcano plots
# writeInteractiveVolcanoPlotProteomicsMain <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 9: writeInteractiveVolcanoPlotProteomicsWidget()
# Current location: R/da_proteins_functions.R
# Description: Creates interactive volcano plot widget
# writeInteractiveVolcanoPlotProteomicsWidget <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 10: getGlimmaVolcanoProteomics()
# Current location: R/da_proteins_functions.R
# Description: Gets Glimma volcano plot for proteomics
# getGlimmaVolcanoProteomics <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 11: getGlimmaVolcanoProteomicsWidget()
# Current location: R/da_proteins_functions.R
# Description: Gets Glimma volcano plot widget
# getGlimmaVolcanoProteomicsWidget <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 12: getGlimmaVolcanoPhosphoproteomics()
# Current location: R/da_proteins_functions.R
# Description: Gets Glimma volcano plot for phosphoproteomics
# getGlimmaVolcanoPhosphoproteomics <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# === PCA Plot Functions ===

# Function 13: plotPca()
# Current location: R/peptideVsSamplesS4Objects.R, R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Creates PCA plots
# setMethod(f = "plotPca", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R or R/proteinVsSamplesS4Objects.R
# }

# Function 14: plotPcaHelper()
# Current location: R/da_proteins_functions.R
# Description: Helper function for PCA plots
# plotPcaHelper <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 15: plotPcaListHelper()
# Current location: R/da_proteins_functions.R
# Description: Helper for creating list of PCA plots
# plotPcaListHelper <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 16: plotPcaGgpairs()
# Current location: R/da_proteins_functions.R
# Description: Creates PCA plots using ggpairs
# plotPcaGgpairs <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 17: rlePcaPlotList()
# Current location: R/da_proteins_functions.R
# Description: Creates list of RLE and PCA plots
# rlePcaPlotList <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# === RLE Plot Functions ===

# Function 18: plotRle()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Creates RLE plots
# setMethod(f = "plotRle", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 19: plotRleHelper()
# Current location: R/da_proteins_functions.R
# Description: Helper function for RLE plots
# plotRleHelper <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 20: plotRleQc()
# Current location: R/da_proteins_functions.R
# Description: Creates RLE QC plots
# plotRleQc <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 21: getOneRlePlotData()
# Current location: R/da_proteins_functions.R
# Description: Gets data for one RLE plot
# getOneRlePlotData <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# === Density Plot Functions ===

# Function 22: plotDensity()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Creates density plots
# setMethod(f = "plotDensity", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 23: plotDensityOfProteinIntensityPerSample()
# Current location: R/da_proteins_functions.R
# Description: Creates density plots of protein intensity per sample
# plotDensityOfProteinIntensityPerSample <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# === Missing Value Plot Functions ===

# Function 24: plotNumMissingValues()
# Current location: R/da_proteins_functions.R
# Description: Plots number of missing values
# plotNumMissingValues <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 25: plotNumOfValues()
# Current location: R/da_proteins_functions.R
# Description: Plots number of values
# plotNumOfValues <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 26: plotNumOfValuesNoLog()
# Current location: R/da_proteins_functions.R
# Description: Plots number of values without log transformation
# plotNumOfValuesNoLog <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 27: plotHistogramOfPercentMissingPerIndvidual()
# Current location: R/da_proteins_functions.R
# Description: Plots histogram of percent missing per individual
# plotHistogramOfPercentMissingPerIndvidual <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# === Count Plot Functions ===

# Function 28: plotPeptidesProteinsCountsPerSample()
# Current location: R/peptideVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Plots peptide and protein counts per sample
# setMethod(f = "plotPeptidesProteinsCountsPerSample", ...) {
#   # Extract from R/peptideVsSamplesS4Objects.R
# }

# Function 29: plotPeptidesProteinsCountsPerSampleHelper()
# Current location: R/da_proteins_functions.R
# Description: Helper for plotting counts per sample
# plotPeptidesProteinsCountsPerSampleHelper <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 30: plotPercentSamplesVsProteinQuantified()
# Current location: R/da_proteins_functions.R
# Description: Plots percent samples vs proteins quantified
# plotPercentSamplesVsProteinQuantified <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 31: plotNumSigDiffExpBarPlot()
# Current location: R/metaboliteVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Plots bar plot of number of significant DE features
# setMethod(f = "plotNumSigDiffExpBarPlot", ...) {
#   # Extract from R/metaboliteVsSamplesS4Objects.R
# }

# === Correlation Plot Functions ===

# Function 32: plotPearson()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Creates Pearson correlation plots
# setMethod(f = "plotPearson", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 33: pearsonCorForSamplePairs()
# Current location: R/proteinVsSamplesS4Objects.R
# Type: S4 method (exportMethods)
# Description: Calculates Pearson correlation for sample pairs
# setMethod(f = "pearsonCorForSamplePairs", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 34: getSamplesCorrelationHeatMap()
# Current location: R/da_proteins_functions.R
# Description: Creates correlation heatmap for samples
# getSamplesCorrelationHeatMap <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 35: getSamplesCorrelationMatrix()
# Current location: R/da_proteins_functions.R
# Description: Gets correlation matrix for samples
# getSamplesCorrelationMatrix <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# === Heatmap Functions ===

# Function 36: getProteinsHeatMap()
# Current location: R/da_proteins_functions.R
# Description: Creates heatmap of proteins
# getProteinsHeatMap <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 37: generateDEHeatmap()
# Current location: R/protein_da_analysis_wrapper.R
# Description: Generates heatmap for DE results
# generateDEHeatmap <- function(...) {
#   # Extract from R/protein_da_analysis_wrapper.R
# }

# Function 38: getEnrichmentHeatmap()
# Current location: R/enrichment_functions.R
# Description: Creates enrichment heatmap
# getEnrichmentHeatmap <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 39: getEnrichmentPlotly()
# Current location: R/enrichment_functions.R
# Description: Creates interactive enrichment plot
# getEnrichmentPlotly <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 40: drawListOfFunctionalEnrichmentHeatmaps()
# Current location: R/enrichment_functions.R
# Description: Draws list of functional enrichment heatmaps
# drawListOfFunctionalEnrichmentHeatmaps <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 41: drawListOfFunctionalEnrichmentHeatmapsScholar()
# Current location: R/enrichment_functions.R
# Description: Draws list of enrichment heatmaps (Scholar version)
# drawListOfFunctionalEnrichmentHeatmapsScholar <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 42: saveListOfFunctionalEnrichmentHeatmaps()
# Current location: R/enrichment_functions.R
# Description: Saves list of functional enrichment heatmaps
# saveListOfFunctionalEnrichmentHeatmaps <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === Enrichment Plot Functions ===

# Function 43: plotEnrichmentBarplot()
# Current location: R/enrichment_functions.R
# Description: Creates enrichment bar plot
# plotEnrichmentBarplot <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 44: enrichedGoTermBarPlot()
# Current location: R/enrichment_functions.R
# Description: Creates enriched GO term bar plot
# enrichedGoTermBarPlot <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 45: enrichedPathwayBarPlot()
# Current location: R/enrichment_functions.R
# Description: Creates enriched pathway bar plot
# enrichedPathwayBarPlot <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 46: plotStringDbEnrichmentResults()
# Current location: R/string_enrichment_functions_refactored.R
# Description: Plots StringDB enrichment results
# plotStringDbEnrichmentResults <- function(...) {
#   # Extract from R/string_enrichment_functions_refactored.R
# }

# Function 47: printStringDbFunctionalEnrichmentBarGraph()
# Current location: R/multiomics_enrichment_functions.R
# Description: Prints StringDB functional enrichment bar graph
# printStringDbFunctionalEnrichmentBarGraph <- function(...) {
#   # Extract from R/multiomics_enrichment_functions.R
# }

# Function 48: generate_enrichment_plots()
# Current location: R/enrichment_functions.R
# Description: Generates enrichment plots
# generate_enrichment_plots <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === Specialized Plot Functions ===

# Function 49: generateLimpaQCPlots()
# Current location: R/limpa_functions.R
# Description: Generates QC plots for limpa imputation
# generateLimpaQCPlots <- function(...) {
#   # Extract from R/limpa_functions.R
# }

# Function 50: generateMetaboliteFilteringPlots()
# Current location: R/QC_visualisation.R
# Description: Generates filtering progress plots for metabolomics
# generateMetaboliteFilteringPlots <- function(...) {
#   # Extract from R/QC_visualisation.R
# }

# Function 51: plotMofaWeights()
# Current location: R/multiomics_functions_MOFA.R
# Description: Plots MOFA factor weights
# plotMofaWeights <- function(...) {
#   # Extract from R/multiomics_functions_MOFA.R
# }

# Function 52: compareUmapComponentsPairs()
# Current location: R/da_proteins_functions.R
# Description: Compares UMAP component pairs
# compareUmapComponentsPairs <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 53: umap_factor_plot()
# Current location: R/da_proteins_functions.R
# Description: Creates UMAP factor plot
# umap_factor_plot <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 54: getMaxMinBoxplot()
# Current location: R/da_proteins_functions.R
# Description: Creates max/min boxplot
# getMaxMinBoxplot <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 55: printPValuesDistribution()
# Current location: R/da_proteins_functions.R
# Description: Prints p-value distribution
# printPValuesDistribution <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# === Plot Saving Functions ===

# Function 56: savePlot()
# Current location: R/helper_functions.R
# Description: Saves plots to file
# savePlot <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 57: save_plot()
# Current location: R/helper_functions.R
# Description: Saves plot (alternative version)
# save_plot <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 58: gg_save_logging()
# Current location: R/helper_functions.R
# Description: Saves plot with logging
# gg_save_logging <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 59: saveListOfPdfs()
# Current location: R/enrichment_functions.R
# Description: Saves list of PDFs
# saveListOfPdfs <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# === Theme and Color Functions ===

# Function 60: apafTheme()
# Current location: R/helper_functions.R
# Description: APAF custom theme for plots
# apafTheme <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 61: get_color_palette()
# Current location: R/helper_functions.R
# Description: Gets color palette
# get_color_palette <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 62: get_enrichplot_color()
# Current location: R/enrichment_functions.R
# Description: Gets enrichplot color
# get_enrichplot_color <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 63: set_enrichplot_color()
# Current location: R/enrichment_functions.R
# Description: Sets enrichplot color
# set_enrichplot_color <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 64: getCategoricalColourPalette()
# Current location: R/da_proteins_functions.R
# Description: Gets categorical color palette
# getCategoricalColourPalette <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 65: getCategoricalAndContinuousColourRules()
# Current location: R/da_proteins_functions.R
# Description: Gets categorical and continuous color rules
# getCategoricalAndContinuousColourRules <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 66: getContinousColourRules()
# Current location: R/da_proteins_functions.R
# Description: Gets continuous color rules
# getContinousColourRules <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 67: getOneContinousPalette()
# Current location: R/da_proteins_functions.R
# Description: Gets one continuous palette
# getOneContinousPalette <- function(...) {
#   # Extract from R/da_proteins_functions.R
# }

# Function 68: get_ggrepel_segsize()
# Current location: R/enrichment_functions.R
# Description: Gets ggrepel segment size
# get_ggrepel_segsize <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 69: calcHtSize()
# Current location: R/enrichment_functions.R
# Description: Calculates heatmap size
# calcHtSize <- function(...) {
#   # Extract from R/enrichment_functions.R
# }











































































