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
# Current location: R/de_proteins_functions.R
# Description: Creates a single volcano plot
# plotOneVolcano <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 3: plotOneVolcanoNoVerticalLines()
# Current location: R/de_proteins_functions.R
# Description: Creates volcano plot without vertical lines
# plotOneVolcanoNoVerticalLines <- function(...) {
#   # Extract from R/de_proteins_functions.R
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
# Current location: R/de_proteins_functions.R
# Description: Prepares data for volcano plot
# prepareDataForVolcanoPlot <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 7: writeInteractiveVolcanoPlotProteomics()
# Current location: R/de_proteins_functions.R
# Description: Writes interactive volcano plot for proteomics
# writeInteractiveVolcanoPlotProteomics <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 8: writeInteractiveVolcanoPlotProteomicsMain()
# Current location: R/de_proteins_functions.R
# Description: Main function for interactive volcano plots
# writeInteractiveVolcanoPlotProteomicsMain <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 9: writeInteractiveVolcanoPlotProteomicsWidget()
# Current location: R/de_proteins_functions.R
# Description: Creates interactive volcano plot widget
# writeInteractiveVolcanoPlotProteomicsWidget <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 10: getGlimmaVolcanoProteomics()
# Current location: R/de_proteins_functions.R
# Description: Gets Glimma volcano plot for proteomics
# getGlimmaVolcanoProteomics <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 11: getGlimmaVolcanoProteomicsWidget()
# Current location: R/de_proteins_functions.R
# Description: Gets Glimma volcano plot widget
# getGlimmaVolcanoProteomicsWidget <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 12: getGlimmaVolcanoPhosphoproteomics()
# Current location: R/de_proteins_functions.R
# Description: Gets Glimma volcano plot for phosphoproteomics
# getGlimmaVolcanoPhosphoproteomics <- function(...) {
#   # Extract from R/de_proteins_functions.R
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
# Current location: R/de_proteins_functions.R
# Description: Helper function for PCA plots
# plotPcaHelper <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 15: plotPcaListHelper()
# Current location: R/de_proteins_functions.R
# Description: Helper for creating list of PCA plots
# plotPcaListHelper <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 16: plotPcaGgpairs()
# Current location: R/de_proteins_functions.R
# Description: Creates PCA plots using ggpairs
# plotPcaGgpairs <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 17: rlePcaPlotList()
# Current location: R/de_proteins_functions.R
# Description: Creates list of RLE and PCA plots
# rlePcaPlotList <- function(...) {
#   # Extract from R/de_proteins_functions.R
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
# Current location: R/de_proteins_functions.R
# Description: Helper function for RLE plots
# plotRleHelper <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 20: plotRleQc()
# Current location: R/de_proteins_functions.R
# Description: Creates RLE QC plots
# plotRleQc <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 21: getOneRlePlotData()
# Current location: R/de_proteins_functions.R
# Description: Gets data for one RLE plot
# getOneRlePlotData <- function(...) {
#   # Extract from R/de_proteins_functions.R
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
# Current location: R/de_proteins_functions.R
# Description: Creates density plots of protein intensity per sample
# plotDensityOfProteinIntensityPerSample <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# === Missing Value Plot Functions ===

# Function 24: plotNumMissingValues()
# Current location: R/de_proteins_functions.R
# Description: Plots number of missing values
# plotNumMissingValues <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 25: plotNumOfValues()
# Current location: R/de_proteins_functions.R
# Description: Plots number of values
# plotNumOfValues <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 26: plotNumOfValuesNoLog()
# Current location: R/de_proteins_functions.R
# Description: Plots number of values without log transformation
# plotNumOfValuesNoLog <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 27: plotHistogramOfPercentMissingPerIndvidual()
# Current location: R/de_proteins_functions.R
# Description: Plots histogram of percent missing per individual
# plotHistogramOfPercentMissingPerIndvidual <- function(...) {
#   # Extract from R/de_proteins_functions.R
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
# Current location: R/de_proteins_functions.R
# Description: Helper for plotting counts per sample
# plotPeptidesProteinsCountsPerSampleHelper <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 30: plotPercentSamplesVsProteinQuantified()
# Current location: R/de_proteins_functions.R
# Description: Plots percent samples vs proteins quantified
# plotPercentSamplesVsProteinQuantified <- function(...) {
#   # Extract from R/de_proteins_functions.R
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
# Current location: R/de_proteins_functions.R
# Description: Creates correlation heatmap for samples
# getSamplesCorrelationHeatMap <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 35: getSamplesCorrelationMatrix()
# Current location: R/de_proteins_functions.R
# Description: Gets correlation matrix for samples
# getSamplesCorrelationMatrix <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# === Heatmap Functions ===

# Function 36: getProteinsHeatMap()
# Current location: R/de_proteins_functions.R
# Description: Creates heatmap of proteins
# getProteinsHeatMap <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 37: generateDEHeatmap()
# Current location: R/protein_de_analysis_wrapper.R
# Description: Generates heatmap for DE results
# generateDEHeatmap <- function(...) {
#   # Extract from R/protein_de_analysis_wrapper.R
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
# Current location: R/de_proteins_functions.R
# Description: Compares UMAP component pairs
# compareUmapComponentsPairs <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 53: umap_factor_plot()
# Current location: R/de_proteins_functions.R
# Description: Creates UMAP factor plot
# umap_factor_plot <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 54: getMaxMinBoxplot()
# Current location: R/de_proteins_functions.R
# Description: Creates max/min boxplot
# getMaxMinBoxplot <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 55: printPValuesDistribution()
# Current location: R/de_proteins_functions.R
# Description: Prints p-value distribution
# printPValuesDistribution <- function(...) {
#   # Extract from R/de_proteins_functions.R
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
# Current location: R/de_proteins_functions.R
# Description: Gets categorical color palette
# getCategoricalColourPalette <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 65: getCategoricalAndContinuousColourRules()
# Current location: R/de_proteins_functions.R
# Description: Gets categorical and continuous color rules
# getCategoricalAndContinuousColourRules <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 66: getContinousColourRules()
# Current location: R/de_proteins_functions.R
# Description: Gets continuous color rules
# getContinousColourRules <- function(...) {
#   # Extract from R/de_proteins_functions.R
# }

# Function 67: getOneContinousPalette()
# Current location: R/de_proteins_functions.R
# Description: Gets one continuous palette
# getOneContinousPalette <- function(...) {
#   # Extract from R/de_proteins_functions.R
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


# ----------------------------------------------------------------------------
# plotPeptidesProteinsCountsPerSampleHelper
# ----------------------------------------------------------------------------
#' plotPeptidesProteinsCountsPerSampleHelper
#' @description Plot the number of proteins and peptides identified per sample
#' @export
plotPeptidesProteinsCountsPerSampleHelper <- function( input_table
                                                 , intensity_column = Peptide.RawQuantity
                                                 , protein_id_column = Protein.Ids
                                                 , peptide_id_column = Stripped.Sequence
                                                 , sample_id_column = Run
                                                 , peptide_sequence_column = Stripped.Sequence ) {

  num_proteins_per_sample <- input_table |>
    dplyr::filter( !is.na(  {{intensity_column}} )) |>
    distinct( {{sample_id_column}}, {{protein_id_column}} ) |>
    group_by( {{sample_id_column}} ) |>
    summarise( count = n()  ) |>
    ungroup()

  num_peptides_per_sample <- input_table |>
    dplyr::filter( !is.na(  {{intensity_column}} )) |>
    distinct( {{sample_id_column}}, {{protein_id_column}}, {{peptide_id_column}} ) |>
    group_by( {{sample_id_column}}) |>
    summarise( count = n()  ) |>
    ungroup()

  combined_counts <- num_proteins_per_sample |>
    mutate( type = "Protein" ) |>
    bind_rows( num_peptides_per_sample |>
                 mutate( type = "Peptide")) |>
    pivot_wider( id_cols = {{sample_id_column}}
                 , names_from = type
                 , values_from = count)

  output_plot <- combined_counts |>
    ggplot( aes( reorder({{sample_id_column}}, Peptide) )) +
    geom_point(aes(y = Peptide/10, shape="Peptide" ),  show.legend = TRUE) +
    geom_point(aes(y = Protein, shape="Protein"  ),  show.legend = TRUE) +
    scale_y_continuous(name = "Protein",
                       sec.axis = sec_axis(\(x) { x*10 }, name =  "Peptide")) +
    apafTheme() +
    theme( axis.text.x=element_blank()
           , axis.ticks.x=element_blank()
           , panel.grid.major.x = element_blank() ) +
    xlab("Samples") +
    scale_shape_manual ( values = c("Peptide" = 1
                                    , "Protein" = 2) ) +
    labs( shape = "Category")

  output_plot
}


# ----------------------------------------------------------------------------
# plotHistogramOfPercentMissingPerIndvidual
# ----------------------------------------------------------------------------
#' @export
plotHistogramOfPercentMissingPerIndvidual <- function( percent_missing_table
                                                       , percent_missing_column = percent_missing) {
  percent_missing_table |>
    ggplot( aes( {{percent_missing_column}})) +
    geom_histogram() +
    apafTheme() +
    xlab("Percent Missing") +
    ylab("Count")


}


# ----------------------------------------------------------------------------
# getOneRlePlotData
# ----------------------------------------------------------------------------
#' @export
getOneRlePlotData <- function( input_matrix ) {


  #if(!( length(which( is.na(input_matrix[, 1]) | is.nan(input_matrix[, 1]) | is.infinite(input_matrix[, 1]) )) > 0 )){

  input_matrix[is.infinite(input_matrix)  | is.nan(input_matrix) ] <- NA

  deviations <- input_matrix - Biobase::rowMedians(input_matrix, na.rm=TRUE)

  stats <-  graphics::boxplot(
    deviations,
    outcol="lightgray",
    cex=0.1,
    cex.axis=0.7,
    las=2,
    outline=FALSE)

  rownames(stats$stats) <- c("lower whisker", "lower hinge", "median", "upper hinge", "upper whisker")
  colnames(stats$stats ) <- colnames(deviations )

  results <- stats$stats |>
    as.data.frame() |>
    rownames_to_column("Quantiles") |>
    pivot_longer( cols = !contains("Quantiles")) |>
    mutate( Quantiles = factor( Quantiles, levels=rev(c( "lower whisker", "lower hinge", "median", "upper hinge", "upper whisker" ))))

  # print(head( results) )



  return(results)
  #}

}


# ----------------------------------------------------------------------------
# plotRleQc
# ----------------------------------------------------------------------------
#' @export
plotRleQc <- function( input_table
                       , x_value = name
                       , y_value = value
                       , quantiles_column = Quantiles ) {

  rle_results <- input_table |>
    ggplot(aes( x={{x_value}}, y={{y_value}},  group={{quantiles_column}}, col={{quantiles_column}})) +
    geom_line()  +
    apafTheme() +
    theme(axis.text.x = element_blank()) +
    xlab("Samples") +
    ylab("Relative log expression") +
    labs(col = "Boxplot features")

  rle_results
}


# ----------------------------------------------------------------------------
# compareUmapComponentsPairs
# ----------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#' @export
compareUmapComponentsPairs <- function(input_table, columns = c("V1", "V2","V3","V4"), covariate) {

  pm <- umap_data |>
    ggpairs( columns = columns, ggplot2::aes(colour = {{covariate}}), legend = 1)  +
    apafTheme()

  pm

}


# ----------------------------------------------------------------------------
# umap_factor_plot
# ----------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#' @export
umap_factor_plot <- function(input_data, header, legend_label, x = V1, y = V2, colour_rule) {


  input_data |>
    mutate( !!sym( {{header}}) := factor( !!sym( {{header}})) ) |>
    ggplot(aes( {{x}}, {{y}}, color = !!sym( {{header}}) )) +
    geom_point() +
    scale_colour_manual( name = legend_label
                         , values=colour_rule
                         , breaks=names( colour_rule)) +
    apafTheme()


}


# ----------------------------------------------------------------------------
# getCategoricalColourPalette
# ----------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------
#' @export
getCategoricalColourPalette <- function() {
  set1_colour <- brewer.pal(9,'Set1')
  set2_colour <- brewer.pal(8,'Set2')
  set3_colour <- brewer.pal(12,'Set3')
  pastel1_colour <- brewer.pal(9,'Pastel1')
  pastel2_colour <- brewer.pal(8,'Pastel2')
  dark2_colour <- brewer.pal(8,'Dark2')
  accent_colour <- brewer.pal(8,'Accent')
  paired_colour <- brewer.pal(12,'Paired')

  set1_2_3_colour <- c( set1_colour, set2_colour, set3_colour
                        , pastel1_colour, pastel2_colour, dark2_colour
                        , accent_colour, paired_colour )

  return(set1_2_3_colour)
}


# ----------------------------------------------------------------------------
# getOneContinousPalette
# ----------------------------------------------------------------------------
#' @export
getOneContinousPalette <- function(metadata_tbl, column_name, palette_name, na_colour = "white") {
  number_of_values <- metadata_tbl |>
    dplyr::select( all_of(column_name)  ) |>
    dplyr::filter(!is.na(!!sym(column_name))) |>
    distinct()  |>
    arrange( !!sym( column_name)  ) |>
    pull() |>
    length()

  list_of_names <- metadata_tbl |>
    dplyr::select( all_of(column_name)  )  |>
    dplyr::filter(!is.na(!!sym(column_name))) |>
    distinct()  |>
    arrange( !!sym( column_name)  )  |>
    pull()

  na_name <- metadata_tbl |>
    dplyr::select( all_of(column_name)  ) |>
    distinct()  |>
    dplyr::filter(is.na(!!sym(column_name))) |>
    pull()

  list_of_colours <- brewer.pal(number_of_values, palette_name)
  names(list_of_colours) <- list_of_names

  if(length(na_name) == 1) {
    new_list_of_colours <- c(list_of_colours, na_colour )
    names(new_list_of_colours) <- c(names(list_of_colours), "NA")
    return( new_list_of_colours )
  }

  return( list_of_colours )
}


# ----------------------------------------------------------------------------
# getContinousColourRules
# ----------------------------------------------------------------------------
#' getContinousColourRules
#' @export
getContinousColourRules <- function( metadata_tbl
                                     , metadata_column_labels
                                     , metadata_column_selected
                                     , continous_scale_columns
                                     , na_colour = "white" ) {

  metadata_column_labels_copy <- metadata_column_labels
  names( metadata_column_labels_copy) <- metadata_column_selected

  list_of_continuous_colour_palette <- c( "Greys", "Blues", "Greens", "Purples", "Reds", "Oranges", "BuGn"
                                          , "BuPu", "GnBu", "OrRd", "PuBu", "PuBuGn", "PuRd"
                                          ,  "RdPu", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd" )

  if( length(continous_scale_columns) > length( list_of_continuous_colour_palette )) {
    list_of_continuous_colour_palette <- rep( list_of_continuous_colour_palette, length.out = length(continous_scale_columns))
  }

  list_of_continous_colour_rules <- purrr::map2(continous_scale_columns
                                                ,  list_of_continuous_colour_palette[seq_along(continous_scale_columns)]
                                                , \(column, palette_name) { getOneContinousPalette(metadata_tbl
                                                                                                   , column
                                                                                                   , palette_name
                                                                                                   , na_colour = na_colour) } )

  names( list_of_continous_colour_rules) <- metadata_column_labels_copy[continous_scale_columns]

  return(list_of_continous_colour_rules)

}


# ----------------------------------------------------------------------------
# getCategoricalAndContinuousColourRules
# ----------------------------------------------------------------------------
#' getCategoricalAndContinuousColourRules
#' @param metadata_tbl This is the table containing sample ID and other columns containing clinical variables / metadata
#' @param metadata_column_labels This is the nice
#' @export
getCategoricalAndContinuousColourRules <- function( metadata_tbl
                                                    , metadata_column_labels
                                                    , metadata_column_selected
                                                    , categorical_columns
                                                    , continous_scale_columns
                                                    , ms_machine_column
                                                    , sample_id_column = Run
                                                    , columns_to_exclude
                                                    , na_colour = "white" ) {

  metadata_column_labels_copy <- metadata_column_labels
  names( metadata_column_labels_copy) <- metadata_column_selected

  if( ms_machine_column %in% columns_to_exclude ) {
    metadata_column_selected <- setdiff( metadata_column_selected, ms_machine_column)
  }

  cln_meatadata_tbl <- metadata_tbl |>
    column_to_rownames(as_name( enquo(sample_id_column ))) |>
    dplyr::select( all_of( c(metadata_column_selected) ) )

  colour_rules <- getCategoricalColourRules( metadata_tbl =  cln_meatadata_tbl
                                             , metadata_column_labels = metadata_column_labels
                                             , metadata_column_selected = metadata_column_selected
                                             , categorical_columns = categorical_columns
                                             , ms_machine_column = ms_machine_column
                                             , columns_to_exclude = columns_to_exclude
                                             , na_colour = na_colour)

  print("Add column annotation")
  colnames(cln_meatadata_tbl) <-  metadata_column_labels_copy[metadata_column_selected]


  continous_colour_list <- getContinousColourRules( metadata_tbl
                                                    , metadata_column_labels
                                                    , metadata_column_selected
                                                    , continous_scale_columns
                                                    , na_colour = na_colour)

  categorical_and_continuous_colour_rules <- c( colour_rules, continous_colour_list)

  columns_to_use <- setdiff(names(categorical_and_continuous_colour_rules), metadata_column_labels_copy[columns_to_exclude])
  categorical_and_continuous_colour_rules_filt <- categorical_and_continuous_colour_rules[columns_to_use]

  return(categorical_and_continuous_colour_rules_filt)
}


# ----------------------------------------------------------------------------
# getSamplesCorrelationHeatMap
# ----------------------------------------------------------------------------
#' getSamplesCorrelationHeatMap
#' @description get the
#' @param correlation_matrix Output from the `getSamplesCorrelationMatrix` function
#' @param metadata_tbl This is the table containing sample ID and other columns containing clinical variables / metadata
#' @param is_HEK_column A logical column in the metadata table that indicates if the sample is a HEK sample
#' @param metadata_column_selected A list of column names in string selected from the metadata tbl
#' @param metadata_column_labels A list of column names in string to rename each of the columns selected in the param `metadata_column_selected`
#' @param categorical_columns A vector of string with all the names of the categorical data column  present in the `metadata_tbl` table
#' @param continous_scale_columns  A vector of string with all the names of the continuous data column  present in the `metadata_tbl` table
#' @param ms_machine_column A string of the column name describing the mass spectrometer machine used to analyze each sample
#' @param sample_id_column A string describing the column name of the sample ID column
#' @export
getSamplesCorrelationHeatMap <- function(correlation_matrix
                                         , metadata_tbl
                                         , is_HEK_column = is_HEK
                                         , metadata_column_labels
                                         , metadata_column_selected
                                         , colour_rules
                                         , columns_to_exclude
                                         , sample_id_column = Run
                                         , use_raster = TRUE
                                         , raster_device = "CairoPDF"
                                         , heatmap_legend_param = list(title = "Correlation")
                                         , heatmap_width = ncol(correlation_matrix)*unit(0.05, "cm")
                                         , heatmap_height = nrow(correlation_matrix)*unit(0.05, "cm")
) {

  names( metadata_column_labels) <- metadata_column_selected

  without_hek_samples <- metadata_tbl |>
    dplyr::filter( {{is_HEK_column}} == FALSE) |>
    pull({{sample_id_column}})

  correlation_samples_to_use <- intersect( colnames(correlation_matrix), without_hek_samples) |> sort()

  cln_meatadata_orig_col_name <- metadata_tbl |>
    dplyr::filter( {{is_HEK_column}} == FALSE) |>
    dplyr::filter( {{sample_id_column}} %in% correlation_samples_to_use) |>
    arrange( {{sample_id_column}}) |>
    column_to_rownames( as_name( enquo( sample_id_column)  ))|>
    dplyr::select( all_of( setdiff( metadata_column_selected, columns_to_exclude) ) )

  columns_to_use <- setdiff(names(colour_rules), metadata_column_labels[columns_to_exclude])
  colour_rules_filt <- colour_rules[columns_to_use]

  print("Add column annotation")
  cln_meatadata_tbl <- cln_meatadata_orig_col_name
  colnames(cln_meatadata_tbl) <-  metadata_column_labels[colnames(cln_meatadata_orig_col_name)]

  top_annotation <- HeatmapAnnotation(df = cln_meatadata_tbl |>
                                        dplyr::select( - any_of(metadata_column_labels[columns_to_exclude] ))
                                      , col = colour_rules_filt
                                      , show_legend = FALSE )

  print("Add row annotation")
  row_ha <- rowAnnotation( df = cln_meatadata_tbl |>
                             dplyr::select( - any_of(metadata_column_labels[columns_to_exclude] ))
                           , col = colour_rules_filt
                           , show_legend = FALSE)

  output_heatmap <- Heatmap(correlation_matrix[correlation_samples_to_use, correlation_samples_to_use]
                            , name="Correlation"
                            , left_annotation = row_ha
                            , top_annotation = top_annotation
                            , show_row_names=FALSE
                            , show_column_names = FALSE
                            , use_raster = use_raster
                            , raster_device = raster_device
                            , row_title_gp = gpar(fontsize = 13.2)
                            , column_title_gp = gpar(fontsize = 13.2)
                            , row_names_gp = gpar(fontsize = 12, fontfamily = "sans")
                            , column_names_gp = gpar(fontsize = 12)
                            , heatmap_legend_param = heatmap_legend_param
                            , heatmap_height = heatmap_height
                            , heatmap_width = heatmap_width
  )

  output_legends <- purrr::map2 (colour_rules_filt
                                 , names(colour_rules_filt)
                                 , \(rule, title) { Legend(labels = names(rule), title = title,
                                                           legend_gp = gpar(fill = rule))} )

  # output_legends <- packLegend(list = list_of_legends)

  return( list( heatmap = output_heatmap
                , legend = output_legends ))


}


# ----------------------------------------------------------------------------
# plotDensityOfProteinIntensityPerSample
# ----------------------------------------------------------------------------
#' @export
plotDensityOfProteinIntensityPerSample <- function( protein_intensity_long_tbl
                                                    , number_of_peptides_per_protein_per_sample
                                                    , protein_id_column = Protein.Ids
                                                    , sample_id_column = Run
                                                    , num_peptides_column = num_peptides_after_impute
                                                    , protein_intensity_column = Log2.Protein.Imputed) {

  protein_intensity_vs_num_peptides_for_replicates <- protein_intensity_long_tbl |>
    left_join( number_of_peptides_per_protein_per_sample
               , by = join_by( {{protein_id_column}}, {{sample_id_column}} )) |>
    mutate( peptides_status = ifelse( {{num_peptides_column}} == 1, "Multiple Peptides"
                                      , "Single Peptide")) |>
    ggplot(aes( {{protein_intensity_column}}, group=peptides_status, fill= peptides_status, alpha=0.5 )) +
    geom_density()  +
    scale_alpha(guide = 'none') +
    apafTheme()  +
    xlab("log2 Protein Intensity") +
    ylab("Density") +
    labs( fill = "Peptide") +
    scale_y_continuous( expand = expansion(  mult=c(0, 0.1)))

}


# ----------------------------------------------------------------------------
# plotPercentSamplesVsProteinQuantified
# ----------------------------------------------------------------------------
#' @export
plotPercentSamplesVsProteinQuantified <- function ( protein_intensity_long_tbl = frozen_protein_table
                                                    , number_of_peptides_per_protein_per_sample = number_of_peptides_per_protein_per_sample
                                                    , protein_id_column = Protein.Ids
                                                    , sample_id_column = Run
                                                    , num_peptides_column = num_peptides_after_impute
                                                    , protein_intensity_column = Log2.Protein.Imputed) {

  samples_vs_intensity <-  protein_intensity_long_tbl  |>
    left_join( number_of_peptides_per_protein_per_sample
               , by = join_by(  {{protein_id_column}}, {{sample_id_column}} )) |>
    mutate( peptides_status = ifelse( num_peptides_after_impute == 1, "Multiple Peptides"
                                      , "Single Peptide"))

  total_num_samples <- samples_vs_intensity |>
    distinct( {{sample_id_column}}) |>
    nrow()

  summarise_peptide_status <- function ( input_vector) {

    if( "Multiple Peptides" %in% input_vector  ) {
      return ( "Multiple Peptides" )
    } else {
      return ( "Single Peptide")
    }
  }

  num_samples_per_protein <- samples_vs_intensity |>
    dplyr::filter(!is.na({{protein_intensity_column}})) |>
    group_by( {{protein_id_column}} ) |>
    summarise( num_values = n()
               , peptides_status = summarise_peptide_status(peptides_status )) |>
    ungroup()  |>
    mutate ( percentage = num_values / total_num_samples * 100 )

  num_samples_per_protein |>
    mutate( percentage_bin = cut( percentage, breaks = c(0, 20, 40, 60, 80, 100) )) |>
    ggplot( aes( percentage_bin, fill = peptides_status, group = peptides_status)) +
    geom_bar( position = "dodge" ) +
    apafTheme()  +
    xlab("Percentage of Samples") +
    ylab("Num. Quantified Proteins ") +
    labs( fill = "Peptide") +
    scale_y_continuous( expand = expansion(  mult=c(0, 0.1)))

}


# ----------------------------------------------------------------------------
# getProteinsHeatMap
# ----------------------------------------------------------------------------
#' @title get protein intensity heatmap
#' @description Generates a heatmap of protein intensities
#'@export
getProteinsHeatMap <- function( protein_matrix
                                , metadata_tbl
                                , is_HEK_column = is_HEK
                                , metadata_column_selected
                                , metadata_column_labels
                                # , categorical_columns
                                # , continous_scale_columns
                                # , ms_machine_column
                                , colour_rules
                                , columns_to_exclude
                                , core_utilisation_samples = TRUE
                                , sort_by_sample_id = TRUE
                                , sample_id_column = Run
                                , use_raster = TRUE
                                , raster_device = "CairoTIFF"
                                , heatmap_legend_param = list(title = "Intensity")) {

  metadata_column_labels_copy <- metadata_column_labels
  names( metadata_column_labels_copy) <- metadata_column_selected

  print("Without HEK samples")
  without_hek_samples <- metadata_tbl |>
    dplyr::filter( {{is_HEK_column}} == FALSE) |>
    pull({{sample_id_column}})

  samples_to_use <- intersect( colnames(protein_matrix), without_hek_samples)

  if( sort_by_sample_id == TRUE) {
    samples_to_use <- samples_to_use |>
      sort()
  }

  cln_meatadata_tbl_orig_col_names <- metadata_tbl |>
    dplyr::filter( {{is_HEK_column}} == FALSE)  |>
    dplyr::filter( {{sample_id_column}} %in% samples_to_use) |>
    arrange( {{sample_id_column}}) |>
    dplyr::select( {{sample_id_column}}, all_of( setdiff(metadata_column_selected, columns_to_exclude) ) ) |>
    distinct()  |>
    column_to_rownames( as_label( enquo(sample_id_column) ) )

  print("Add column annotation")
  cln_meatadata_tbl <- cln_meatadata_tbl_orig_col_names[ samples_to_use
                                                         , setdiff(metadata_column_selected, columns_to_exclude)]
  colnames(cln_meatadata_tbl) <- metadata_column_labels_copy[setdiff(metadata_column_selected, columns_to_exclude)]
  colour_rules_filt <- colour_rules[metadata_column_labels_copy[setdiff(metadata_column_selected, columns_to_exclude)] ]

  # print(colour_rules)
  # print(metadata_column_labels_copy[setdiff(metadata_column_selected, columns_to_exclude)] )
  # print(colour_rules)
  # print(colour_rules_filt)

  # print(colour_rules_filt)

  print("Set top located annotation")
  top_annotation <- HeatmapAnnotation( df = cln_meatadata_tbl
                                       , col = colour_rules_filt
                                       , show_legend = FALSE
                                       , annotation_name_side = "left" )

  print("Print Heatmap")

  heatmap <- Heatmap( protein_matrix[, samples_to_use]
                      , name="Intensity"
                      , top_annotation = top_annotation
                      , show_row_names=FALSE
                      , show_column_names = FALSE
                      , use_raster = use_raster
                      , raster_device = raster_device
                      , row_title_gp = gpar(fontsize = 13.2)
                      , column_title_gp = gpar(fontsize = 13.2)
                      , row_names_gp = gpar(fontsize = 12, fontfamily = "sans")
                      , column_names_gp = gpar(fontsize = 12)
                      , heatmap_legend_param = heatmap_legend_param
                      , core_utilisation_columns = core_utilisation_samples )

  output_legends <- purrr::map2 (colour_rules_filt
                                 , names(colour_rules_filt)
                                 , \(rule, title) { Legend(labels = names(rule), title = title,
                                                           legend_gp = gpar(fill = rule))} )

  return(list( heatmap = heatmap
               , legend = output_legends ))

}


# ----------------------------------------------------------------------------
# apafTheme
# ----------------------------------------------------------------------------
#' @title ProCan ggplot2 theme. Rectangle box around each plot.
#' @description Standard ggplot2 theme for ProCan plots
#' @export
apafTheme <- function() {
  theme(
    # Set font family and size
    text = element_text(family = "Arial", size = 12),
    # Add rectangular box around the plot
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    # Add grid lines
    panel.grid.major = element_line(color = "gray", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    # Set plot background color
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    # Set axis line and tick colors
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    # Set axis label colors and sizes
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    # Set legend title and label colors and sizes
    legend.title = element_text(color = "black", size = 12),
    legend.text = element_text(color = "black", size = 10),
    # Set plot title and subtitle colors and sizes
    plot.title = element_text(color = "black", size = 14),
    plot.subtitle = element_text(color = "black", size = 12),
    # Set plot margin sizes
    plot.margin = unit(c(1, 1, 1, 1), "cm")

  )
}


# ----------------------------------------------------------------------------
# get_color_palette
# ----------------------------------------------------------------------------
#' Generate a color palette
#' 
#' @param n Number of colors needed
#' @param base_color Base color to use
#' @return Vector of colors
#' @export 
get_color_palette <- function(n, base_color) {
  colorRampPalette(c(base_color, "black"))(n)
}


# ----------------------------------------------------------------------------
# plotNumMissingValues
# ----------------------------------------------------------------------------
#' Plot the number of missing values in each sample
#'@param input_table  Data matrix with each row as a protein and each column a sample.
#'@return A ggplot2 bar plot showing the number of missing values per column.
#'@export
plotNumMissingValues <- function(input_table) {

  plot_num_missing_values <- apply(data.matrix(log2(input_table)), 2,
                                   function(x) { length(which(!is.finite(x))) }) |>
    t() |>
    t() |>
    set_colnames("No. of Missing Values") |>
    as.data.frame() |>
    rownames_to_column("Samples ID") |>
    ggplot(aes(x = `Samples ID`, y = `No. of Missing Values`)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90))

  plot_num_missing_values
}


# ----------------------------------------------------------------------------
# plotNumOfValues
# ----------------------------------------------------------------------------
#' Plot the number of values in each sample
#'@param input_table  Data matrix with each row as a protein and each column a sample.
#'@return A ggplot2 bar plot showing the number of missing values per column.
#'@export
plotNumOfValues <- function(input_table) {

  plot_num_missing_values <- apply(data.matrix(log2(input_table)), 2,
                                   function(x) { length(which(!is.na(x))) }) |>
    t() |>
    t() |>
    set_colnames("No. of Values") |>
    as.data.frame() |>
    rownames_to_column("Samples ID") |>
    ggplot(aes(x = `Samples ID`, y = `No. of Values`)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90))

  plot_num_missing_values
}


# ----------------------------------------------------------------------------
# plotNumOfValuesNoLog
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Plot the number of values in each sample
#'@param input_table  Data matrix with each row as a protein and each column a sample.
#'@return A ggplot2 bar plot showing the number of missing values per column.
#'@export
plotNumOfValuesNoLog <- function(input_table) {

  plot_num_missing_values <- apply(data.matrix(input_table), 2,
                                   function(x) { length(which(!is.na(x))) }) |>
    t() |>
    t() |>
    set_colnames("No. of Values") |>
    as.data.frame() |>
    rownames_to_column("Samples ID") |>
    ggplot(aes(x = `Samples ID`, y = `No. of Values`)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90))

  plot_num_missing_values
}


# ----------------------------------------------------------------------------
# plotPcaHelper
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Plot PCA Helper
#'
#' @description
#' This function performs a Principal Component Analysis (PCA) on the given data and generates a ggplot object for visualization.
#' It includes filtering steps to select the most variable features before running PCA.
#'
#' @param data A matrix or data frame of quantitative data, with features in rows and samples in columns.
#' @param design_matrix A data frame containing sample metadata.
#' @param sample_id_column The name of the column in `design_matrix` that identifies samples. Defaults to "Sample_ID".
#' @param grouping_variable The column in `design_matrix` used for coloring points in the PCA plot. Defaults to "group".
#' @param shape_variable An optional column in `design_matrix` to be used for the shape aesthetic of points.
#' @param label_column An optional column in `design_matrix` to be used for labeling points.
#' @param title The title of the plot.
#' @param geom.text.size The size of the text labels if `label_column` is provided. Defaults to 11.
#' @param ncomp The number of principal components to compute. Defaults to 2.
#' @param cv_percentile The percentile threshold for selecting features based on coefficient of variation. Defaults to 0.90 (top 10% most variable features).
#' @param ... Additional arguments passed to other methods (not currently used).
#'
#' @return A `ggplot` object representing the PCA plot.
#'
#' @importFrom mixOmics pca
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab labs theme element_blank scale_x_continuous scale_y_continuous coord_cartesian
#' @importFrom ggrepel geom_text_repel
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom rlang sym
#'
#' @examples
#' # Create dummy data
#' data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' colnames(data) <- paste0("Sample", 1:10)
#' rownames(data) <- paste0("Protein", 1:100)
#'
#' # Create dummy design matrix
#' design_matrix <- data.frame(
#'   Sample_ID = paste0("Sample", 1:10),
#'   group = rep(c("A", "B"), each = 5),
#'   batch = rep(c("X", "Y"), times = 5)
#' )
#'
#' # Generate PCA plot
#' plotPcaHelper(
#'   data = data,
#'   design_matrix = design_matrix,
#'   title = "PCA of Dummy Data",
#'   grouping_variable = "group",
#'   shape_variable = "batch"
#' )
#'
#' @export
plotPcaHelper <- function(data,
                          design_matrix,
                          sample_id_column = "Sample_ID",
                          grouping_variable = "group",
                          shape_variable = NULL,
                          label_column = NULL,
                          title, geom.text.size = 11, ncomp = 2,
                          cv_percentile = 0.90,
                          ...) {
  
  # Ensure design_matrix is a data frame
  design_matrix <- as.data.frame(design_matrix)
  
  # --- START of modification ---
  if(nrow(data) > 1) {
    # 1. Filter by presence in samples
    # Keep features that have values in at least 5% of the samples.
    num_samples <- ncol(data)
    min_samples_present <- floor(0.05 * num_samples)
    data_abundant <- data[rowSums(!is.na(data)) >= min_samples_present, ]
    
    # 2. Filter by Coefficient of Variation (CV) on the already filtered data
    if(nrow(data_abundant) > 1) {
      cvs <- apply(data_abundant, 1, function(x) {
        if (all(is.na(x))) return(NA)
        s <- sd(x, na.rm = TRUE)
        m <- mean(x, na.rm = TRUE)
        if (is.na(s) || is.na(m) || abs(m) < 1e-6) return(0)
        s / m
      })
      
      # Find the CV threshold based on the specified percentile
      cv_threshold <- quantile(cvs, cv_percentile, na.rm = TRUE)
      
      if(!is.na(cv_threshold) && cv_threshold > 0) {
        data_filtered <- data_abundant[which(cvs >= cv_threshold), ]
        print(nrow(data))
        print(nrow(data_filtered))
      } else {
        data_filtered <- data_abundant
      }
    } else {
      data_filtered <- data_abundant
    }
  } else {
    data_filtered <- data
  }
  # --- END of modification ---
  
  pca.res <- mixOmics::pca(t(as.matrix(data_filtered)), ncomp = ncomp)
  proportion_explained <- pca.res$prop_expl_var
  
  temp_tbl <- pca.res$variates$X |>
    as.data.frame() |>
    rownames_to_column(var = sample_id_column) |>
    left_join(design_matrix, by = sample_id_column)
  
  # --- START of fix for shape aesthetic ---
  # Ensure the shape variable is treated as a discrete factor for plotting
  if (!is.null(shape_variable) && shape_variable %in% colnames(temp_tbl)) {
    temp_tbl[[shape_variable]] <- as.factor(temp_tbl[[shape_variable]])
  }
  # --- END of fix for shape aesthetic ---
  
  # More defensive check for grouping variables
  if (!grouping_variable %in% colnames(temp_tbl)) {
    stop(sprintf("Grouping variable '%s' not found in the data", grouping_variable))
  }
  
  if (!is.null(shape_variable) && !shape_variable %in% colnames(temp_tbl)) {
    stop(sprintf("Shape variable '%s' not found in the data", shape_variable))
  }
  
  # Create base plot with appropriate aesthetics based on whether shape_variable is NULL
  if (is.null(label_column) || label_column == "") {
    if (is.null(shape_variable)) {
      # No shape variation, only color
      base_plot <- temp_tbl |>
        ggplot(aes(PC1, PC2, color = !!sym(grouping_variable)))
    } else {
      # Both color and shape
      base_plot <- temp_tbl |>
        ggplot(aes(PC1, PC2, color = !!sym(grouping_variable), shape = !!sym(shape_variable)))
    }
  } else {
    if (!label_column %in% colnames(temp_tbl)) {
      stop(sprintf("Label column '%s' not found in the data", label_column))
    }
    
    if (is.null(shape_variable)) {
      # No shape variation, only color, with labels
      base_plot <- temp_tbl |>
        ggplot(aes(PC1, PC2, color = !!sym(grouping_variable), label = !!sym(label_column)))
    } else {
      # Both color and shape, with labels
      base_plot <- temp_tbl |>
        ggplot(aes(PC1, PC2, color = !!sym(grouping_variable), shape = !!sym(shape_variable), 
                   label = !!sym(label_column)))
    }
  }
  
  # Calculate the percentage label for axis (e.g., 0.90 -> "top 10%")
  cv_percent_label <- paste0("top ", round((1 - cv_percentile) * 100, 0), "%")
  
  output <- base_plot +
    geom_point(size = 3) +
    xlab(paste("PC1 (", round(proportion_explained$X[["PC1"]] * 100, 0), "% of ", cv_percent_label, " CV)", sep = "")) +
    ylab(paste("PC2 (", round(proportion_explained$X[["PC2"]] * 100, 0), "% of ", cv_percent_label, " CV)", sep = "")) +
    labs(title = title) +
    theme(legend.title = element_blank()) +
    scale_x_continuous(labels = function(x) format(x, scientific = FALSE, digits = 3)) +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE, digits = 3))
  
  # Add explicit color and shape scales to support >6 discrete levels
  # Get color palette (supports many levels)
  categorical_colors <- getCategoricalColourPalette()
  
  # Add explicit color scale
  output <- output + scale_color_manual(values = categorical_colors)
  
  # Add explicit shape scale (15 distinct shapes - first 6 match ggplot2 defaults)
  if (!is.null(shape_variable)) {
    shape_values <- c(16, 17, 15, 3, 7, 8, 0, 1, 2, 4, 5, 6, 9, 10, 11)
    output <- output + scale_shape_manual(values = shape_values)
  }

  # Calculate axis limits based on the full range of data
  pc1_range <- range(temp_tbl$PC1, na.rm = TRUE)
  pc2_range <- range(temp_tbl$PC2, na.rm = TRUE)
  buffer_pc1 <- (pc1_range[2] - pc1_range[1]) * 0.05 # 5% buffer
  buffer_pc2 <- (pc2_range[2] - pc2_range[1]) * 0.05 # 5% buffer
  
  output <- output + coord_cartesian(
    xlim = c(pc1_range[1] - buffer_pc1, pc1_range[2] + buffer_pc1),
    ylim = c(pc2_range[1] - buffer_pc2, pc2_range[2] + buffer_pc2)
  )
  
  if (!is.null(label_column) && label_column != "") {
    output <- output + geom_text_repel(size = geom.text.size, show.legend = FALSE)
  }
  
  #class(output) <- "ggplot"
  output
}


# ----------------------------------------------------------------------------
# plotPcaListHelper
# ----------------------------------------------------------------------------
#'@export
plotPcaListHelper <- function(data,
                              design_matrix,
                              sample_id_column = "Sample_ID",
                              grouping_variables_list = c("group"),
                              label_column = NULL,
                              title, geom.text.size = 11, ncomp = 2,
                              cv_percentile = 0.90,
                              ...) {
  
  
  # --- START of modification ---
  if(nrow(data) > 1) {
    # 1. Filter by presence in samples
    # Keep features that have values in at least 5% of the samples.
    num_samples <- ncol(data)
    min_samples_present <- floor(0.05 * num_samples)
    data_abundant <- data[rowSums(!is.na(data)) >= min_samples_present, ]
    
    # 2. Filter by Coefficient of Variation (CV) on the already filtered data
    if(nrow(data_abundant) > 1) {
      cvs <- apply(data_abundant, 1, function(x) {
        if (all(is.na(x))) return(NA)
        s <- sd(x, na.rm = TRUE)
        m <- mean(x, na.rm = TRUE)
        if (is.na(s) || is.na(m) || abs(m) < 1e-6) return(0)
        s / m
      })
      
      # Find the CV threshold based on the specified percentile
      cv_threshold <- quantile(cvs, cv_percentile, na.rm = TRUE)
      
      if(!is.na(cv_threshold) && cv_threshold > 0) {
        data_filtered <- data_abundant[which(cvs >= cv_threshold), ]
        print(nrow(data))
        print(nrow(data_filtered))
      } else {
        data_filtered <- data_abundant
      }
    } else {
      data_filtered <- data_abundant
    }
  } else {
    data_filtered <- data
  }
  # --- END of modification ---
  

  pca.res <- mixOmics::pca(t(as.matrix(data_filtered)), ncomp = ncomp)
  proportion_explained <- pca.res$prop_expl_var

  temp_tbl <- pca.res$variates$X |>
    as.data.frame() |>
    rownames_to_column(var = sample_id_column) |>
    left_join(design_matrix, by = sample_id_column)


  plotOneGgplotPca <- function( grouping_variable ) {
    unique_groups <- temp_tbl |> distinct(!!sym(grouping_variable)) |> dplyr::pull(!!sym(grouping_variable))

    if (is.null(label_column) || label_column == "") {
      output <- temp_tbl |>
        ggplot(aes(PC1, PC2, col = !!sym(grouping_variable))) +
        geom_point() +
        xlab(paste("PC1 (", round(proportion_explained$X[["PC1"]] * 100, 0), "%)", sep = "")) +
        ylab(paste("PC2 (", round(proportion_explained$X[["PC2"]] * 100, 0), "%)", sep = "")) +
        labs(title = title) +
        theme(legend.title = element_blank())
    } else {
      output <- temp_tbl |>
        ggplot(aes(PC1, PC2, col = !!sym(grouping_variable), label = !!sym(label_column))) +
        geom_point() +
        geom_text_repel(size = geom.text.size, show.legend = FALSE) +
        xlab(paste("PC1 (", round(proportion_explained$X[["PC1"]] * 100, 0), "%)", sep = "")) +
        ylab(paste("PC2 (", round(proportion_explained$X[["PC2"]] * 100, 0), "%)", sep = "")) +
        labs(title = title) +
        theme(legend.title = element_blank())
    }

    return(output)
  }

  output_list <- purrr::map( grouping_variables_list, plotOneGgplotPca)

  output_list
}


# ----------------------------------------------------------------------------
# plotPcaGgpairs
# ----------------------------------------------------------------------------
#'@export
plotPcaGgpairs <- function( data_matrix
                            , design_matrix
                            , grouping_variable
                            , sample_id_column
                            , ncomp=2 ) {

  pca.res <- mixOmics::pca(t(as.matrix(data_matrix)), ncomp=ncomp)


  pca_prop_explained_helper <- function( pca_obj, comp_idx ) {
    proportion_explained <- pca.res$prop_expl_var

    pc_label <- paste0("PC", comp_idx)

    perc_label <- paste( paste0(pc_label, " ("), round(proportion_explained$X[[pc_label]]*100, 0),"%)", sep="")

    perc_label
  }

  pc_list <- purrr::map_chr( seq_len(ncomp)
                             , \(comp_idx){ pca_prop_explained_helper(pca.res, comp_idx)})

  pca_variates_x <- pca.res$variates$X

  colnames(pca_variates_x) <- pc_list

  pca_plot_ggpairs <- pca_variates_x |>
    as.data.frame() |>
    rownames_to_column( sample_id_column ) |>
    left_join( design_matrix
               , by = join_by( !!sym(sample_id_column) ==  !!sym(sample_id_column)) ) |>
    ggpairs( columns=pc_list, aes( colour = !!sym(grouping_variable), fill= !!sym(grouping_variable), alpha=0.4)
             , legend = 1)

  pca_plot_ggpairs
}


# ----------------------------------------------------------------------------
# plotRleHelper
# ----------------------------------------------------------------------------
#'@title Plot RLE
#'@export
#'@param Y  Rows = Samples, Columns = Proteins or Peptides
plotRleHelper <- function(Y, rowinfo = NULL, probs = c(0.05, 0.25, 0.5, 0.75,
                                                       0.95), yaxis_limit = c(-0.5, 0.5))
{
  #  checks = check.ggplot()
  # if (checks) {
  rle <- t(apply(t(Y) - apply(Y, 2, function(x){median(x, na.rm=TRUE)}), 2, function(x){quantile(x, probs = probs, na.rm=TRUE)}))
  colnames(rle) <- c("min", "lower", "middle", "upper",
                     "max")
  df <- cbind(data.frame(rle.x.factor = rownames(rle)), data.frame(rle))

  if (!is.null(rowinfo)) {
    rowinfo <- data.frame(rowinfo = rowinfo)
    df_temp <- cbind(df, rowinfo)

    my.x.factor.levels <- df_temp |>
      arrange(rowinfo) |>
      distinct(rle.x.factor) |>
      dplyr::pull(rle.x.factor)

    df <- df_temp |>
      mutate(rle.x.factor = factor(rle.x.factor,
                                   levels = my.x.factor.levels)) |>
      arrange(rowinfo)
  }

  rleplot <- ggplot(df, aes(x = .data[["rle.x.factor"]])) +
    geom_boxplot(aes(lower = .data[["lower"]]
                     , middle = .data[["middle"]]
                     , upper = .data[["upper"]]
                     , max = .data[["max"]]
                     , min = .data[["min"]]),
                 stat = "identity") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90) #, axis.ticks.x = element_blank()
    ) +
    theme(axis.title.y = element_blank(), axis.text.y = element_text(size = rel(1.5))) +
    geom_hline(yintercept = 0)


  if( length( yaxis_limit ) ==2 ) {

    rleplot <- rleplot +
      coord_cartesian(ylim = yaxis_limit)

  }


  if (!is.null(rowinfo)) {
    if (ncol(rowinfo) == 1) {
      rleplot <- rleplot + aes(fill = rowinfo) + labs(fill = "")
    }
  }

  return(rleplot)
  # }
  # else return(FALSE)
}


# ----------------------------------------------------------------------------
# getMaxMinBoxplot
# ----------------------------------------------------------------------------
#'@title Get Max and Min Boxplot
#' @export
#' @description Input a ggplot2 boxplot, return the maximum and minimum data point adjusted by the adjust_factor.
#' @param input_boxplot A ggplot2 boxplot object.
#' @param adjust_factor A numeric value to adjust the maximum and minimum data point.
getMaxMinBoxplot <- function( input_boxplot, adjust_factor = 0.05) {

  df_min <- min( input_boxplot$data$min, na.rm=TRUE)

  df_max <- max( input_boxplot$data$max, na.rm=TRUE )

  if( df_min > 0 ) {
    df_min <- df_min*(1-adjust_factor)
  } else {
    df_min <- df_min*(1+adjust_factor)

  }

  if( df_max > 0 ) {
    df_max <- df_max*(1+adjust_factor)
  } else {
    df_max <- df_max*(1-adjust_factor)

  }

  return( c(df_min, df_max))
}


# ----------------------------------------------------------------------------
# rlePcaPlotList
# ----------------------------------------------------------------------------
#'@export
rlePcaPlotList <- function(list_of_data_matrix, list_of_design_matrix,
                           sample_id_column = Sample_ID, grouping_variable = group, list_of_descriptions) {

  rle_list <- purrr::pmap( list( data_matrix=list_of_data_matrix, description=list_of_descriptions, design_matrix=list_of_design_matrix),
                           function( data_matrix, description, design_matrix) { plotRleHelper(t(as.matrix(data_matrix)),
                                                                                              rowinfo = design_matrix[colnames(data_matrix), as_name(enquo(grouping_variable))]  )  +
                               labs(title = description)} )

  pca_list <- purrr::pmap(list( data_matrix=list_of_data_matrix, description=list_of_descriptions, design_matrix=list_of_design_matrix),
                          function( data_matrix, description, design_matrix) { plotPcaHelper(data_matrix,
                                                                                             design_matrix = design_matrix,
                                                                                             sample_id_column = sample_id_column ,
                                                                                             grouping_variable =  grouping_variable ,
                                                                                             title = description, cex = 7) })

  list_of_plots <- c(rle_list, pca_list)

  rle_pca_plots_arranged <- ggarrange(plotlist = list_of_plots, nrow = 2, ncol = length(list_of_descriptions),
                                      common.legend = FALSE, legend = "bottom", widths = 10, heights = 10)

  rle_pca_plots_arranged
}


# ----------------------------------------------------------------------------
# plotOneVolcano
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Draw the volcano plot, used in publication graphs
#' @param input_data Input data with the `log_q_value_column`, the `log_fc_column`, the `points_type_label` and the `points_color` columns.
#' @param log_q_value_column The name of the column representing the log q-value.
#' @param log_fc_column The name of the column representing the log fold-change.
#' @param points_type_label A column in input table with the type of points based on log fold-change and q-value (e.g. "Not sig., logFC >= 1" = "orange" , "Sig., logFC >= 1" = "purple" , "Sig., logFC < 1" = "blue" , "Not sig." )
#' @param points_color A column in input table with the colour of the points corresponding to each type of points (e.g. orange, purple, blue black, )
#' @param q_val_thresh A numerical value specifying the q-value threshold for statistically significant proteins.
#' @param log2FC_thresh A numerical value specifying the log fold-change threshold to draw a vertical line
#' @param formula_string The formula string used in the facet_grid command for the ggplot scatter plot.
#' @export
plotOneVolcano <- function( input_data, input_title,
                            log_q_value_column = lqm,
                            log_fc_column = logFC,
                            points_type_label = label,
                            points_color = colour,
                            q_val_thresh=0.05 ,
                            log2FC_thresh = 1) {

  colour_tbl <- input_data |>
    distinct( {{points_type_label}}, {{points_color}} )

  # print(colour_tbl)

  colour_map <- colour_tbl |>
    dplyr::pull({{points_color}} ) |>
    as.vector()

  names( colour_map ) <- colour_tbl |>
    dplyr::pull({{points_type_label}} )

  avail_labels <- input_data |>
    distinct({{points_type_label}}) |>
    dplyr::pull({{points_type_label}})

  avail_colours <- colour_map[avail_labels]

  # print(avail_labels)
  # print(avail_colours)

  volcano_plot <-  input_data |>
    ggplot(aes(y = {{log_q_value_column}},
               x = {{log_fc_column}} )) +
    geom_point(aes(col = label)) +
    scale_colour_manual(values = avail_colours)  +
    geom_hline(yintercept = -log10(q_val_thresh)) +
    theme_bw() +
    xlab(expression(Log[2](`fold-change`))) +
    ylab(expression(-log[10](`q-value`))) +
    labs(title = input_title)+  # Remove legend title
    theme(legend.title = element_blank()) +
    # theme(legend.position = "none")  +
    theme(axis.text.x = element_text(size = 13))   +
    theme(axis.text.y = element_text(size = 13))  +
    theme(axis.title.x = element_text(size = 12))  +
    theme(axis.title.y = element_text(size = 12))  +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) # +
  # theme(legend.title = element_text(size = 12))

  if( !is.na(log2FC_thresh) ) {
  volcano_plot <- volcano_plot+
    geom_vline(xintercept = log2FC_thresh, colour = "black", linewidth = 0.2) +
    geom_vline(xintercept = -log2FC_thresh, colour = "black", linewidth = 0.2)

  }

  volcano_plot
}


# ----------------------------------------------------------------------------
# plotOneVolcanoNoVerticalLines
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Draw the volcano plot, used in publication graphs
#' @param input_data Input data with the `log_q_value_column`, the `log_fc_column`, the `points_type_label` and the `points_color` columns.
#' @param log_q_value_column The name of the column representing the log q-value.
#' @param log_fc_column The name of the column representing the log fold-change.
#' @param points_type_label A column in input table with the type of points based on log fold-change and q-value (e.g. "Not sig., logFC >= 1" = "orange" , "Sig., logFC >= 1" = "purple" , "Sig., logFC < 1" = "blue" , "Not sig." )
#' @param points_color A column in input table with the colour of the points corresponding to each type of points (e.g. orange, purple, blue black, )
#' @param q_val_thresh A numerical value specifying the q-value threshold for statistically significant proteins.
#' @param gene_name The column representing the gene name
#' @param formula_string The formula string used in the facet_grid command for the ggplot scatter plot.
#' @export
plotOneVolcanoNoVerticalLines <- function( input_data, input_title,
                            log_q_value_column = lqm,
                            log_fc_column = logFC,
                            points_type_label = label,
                            points_color = colour,
                            q_val_thresh=0.05) {

  colour_tbl <- input_data |>
    distinct( {{points_type_label}}, {{points_color}} )

  colour_map <- colour_tbl |>
    dplyr::pull({{points_color}} ) |>
    as.vector()

  names( colour_map ) <- colour_tbl |>
    dplyr::pull({{points_type_label}} )

  avail_labels <- input_data |>
    distinct({{points_type_label}}) |>
    dplyr::pull({{points_type_label}})

  avail_colours <- colour_map[avail_labels]

  volcano_plot <-  input_data |>
    ggplot(aes(y = {{log_q_value_column}},
               x = {{log_fc_column}},
               col={{points_type_label}})) +
    geom_point()

  volcano_plot <-   volcano_plot +
    scale_colour_manual(values = avail_colours) +
    # geom_vline(xintercept = 1, colour = "black", size = 0.2) +
    # geom_vline(xintercept = -1, colour = "black", size = 0.2) +
    geom_hline(yintercept = -log10(q_val_thresh)) +
    theme_bw() +
    xlab(expression(Log[2](`fold-change`))) +
    ylab(expression(-log[10](FDR))) +
    labs(title = input_title)+  # Remove legend title
    theme(legend.title = element_blank()) +
    # theme(legend.position = "none")  +
    theme(axis.text.x = element_text(size = 13))   +
    theme(axis.text.y = element_text(size = 13))  +
    theme(axis.title.x = element_text(size = 12))  +
    theme(axis.title.y = element_text(size = 12))  +
    theme(plot.title = element_text(size = 12)) +
    theme(legend.text = element_text(size = 12)) # +
  # theme(legend.title = element_text(size = 12))

  volcano_plot
}


# ----------------------------------------------------------------------------
# printOneVolcanoPlotWithProteinLabel
# ----------------------------------------------------------------------------
#'  This function creates a volcano plot with protein labels
#' @param input_table The input table to be used for the volcano plot, contains the protein_id_column, fdr_column and log2FC_column
#' @param uniprot_table The uniprot table to be used for the volcano plot, contains the uniprot_protein_id_column and gene_name_column
#' @param protein_id_column The column name in the input_table that contains the protein ids (tidyverse format)
#' @param uniprot_protein_id_column The column name in the uniprot_table that contains the uniprot protein ids (tidyverse format)
#' @param gene_name_column The column name in the uniprot_table that contains the gene names (tidyverse format)
#' @param number_of_genes Increasing P-value rank for the number of proteins to display on the volcano plot, default is 100`
#' @param fdr_threshold The FDR threshold to use for the volcano plot, default is 0.05
#' @param fdr_column The column name in the input_table that contains the FDR values (tidyverse format)
#' @param log2FC_column The column name in the input_table that contains the log2FC values (tidyverse format)
#' @param input_title The title to use for the volcano plot
#' @param max.overlaps The maximum number of overlaps to allow for the protein labels using ggrepel (default is 20)
#' @return A ggplot object with the volcano plot and protein labels
printOneVolcanoPlotWithProteinLabel <- function( input_table
                                                 , uniprot_table
                                                 , protein_id_column = Protein.Ids
                                                 , uniprot_protein_id_column = Entry
                                                 , gene_name_column = gene_name
                                                 , number_of_genes = 100
                                                 , fdr_threshold = 0.05
                                                 , fdr_column = fdr_qvalue
                                                 , log2FC_column = log2FC
                                                 , input_title = "Proteomics"
                                                 , include_protein_label = TRUE
                                                 , max.overlaps = 20) {
  proteomics_volcano_tbl <- prepareDataForVolcanoPlot( input_table
                                                       , protein_id_column = {{protein_id_column}}
                                                       , uniprot_table = uniprot_table
                                                       , uniprot_protein_id_column = {{uniprot_protein_id_column}}
                                                       , gene_name_column = {{gene_name_column}}
                                                       , number_of_genes = number_of_genes
                                                       , fdr_threshold = fdr_threshold
                                                       , fdr_column = {{fdr_column}}
                                                       , log2FC_column = {{log2FC_column}})

  proteomics_volcano_plot <- plotOneVolcanoNoVerticalLines ( proteomics_volcano_tbl,
                                                                                 input_title = input_title,
                                                                                 log_q_value_column = lqm,
                                                                                 log_fc_column = log2FC,
                                                                                 points_type_label = label,
                                                                                 points_color = colour,
                                                                                 q_val_thresh= fdr_threshold)

  proteomics_volcano_plot_with_proteins_label <- proteomics_volcano_plot

  if( include_protein_label == TRUE) {

    proteomics_volcano_plot_with_proteins_label <- proteomics_volcano_plot +
      ggrepel::geom_text_repel( aes( label = gene_name_significant)
                                , show.legend = FALSE
                                , max.overlaps = max.overlaps)

  }

  proteomics_volcano_plot_with_proteins_label

}


# ----------------------------------------------------------------------------
# getGlimmaVolcanoProteomics
# ----------------------------------------------------------------------------
getGlimmaVolcanoProteomics <- function( r_obj
                                        , coef
                                        , volcano_plot_tab
                                        , uniprot_column = best_uniprot_acc
                                        , gene_name_column = gene_name
                                        , display_columns = c(  "PROTEIN_NAMES"   )
                                        , additional_annotations = NULL
                                        , additional_annotations_join_column = NULL
                                        , counts_tbl = NULL
                                        , groups = NULL
                                        , output_dir) {

  if( coef <= ncol(r_obj$coefficients )) {

    best_uniprot_acc <- str_split(rownames(r_obj@.Data[[1]]), " |:" ) |>
      purrr::map_chr(1)

    volcano_plot_tab_cln <- volcano_plot_tab |>
      # dplyr::rename( best_uniprot_acc =  {{uniprot_column}}
      #                , gene_name = {{gene_name_column}}   ) |>
      dplyr::select( {{uniprot_column}}
                     , {{gene_name_column}}, any_of( display_columns) ) |>
      distinct()


    if( !is.null( additional_annotations )
        & !is.null( additional_annotations_join_column ) ) {

      volcano_plot_tab_cln <- volcano_plot_tab_cln |>
        left_join( additional_annotations
                   , by = join_by( {{uniprot_column}} == {{additional_annotations_join_column}} ) ) |>
        dplyr::select( {{uniprot_column}}
                       , {{gene_name_column}}
                       ,any_of( display_columns))

    }

    anno_tbl <- data.frame( uniprot_acc = rownames(r_obj@.Data[[1]])
                            , temp_column = best_uniprot_acc ) |>
      dplyr::rename( {{uniprot_column}} := temp_column)

      anno_tbl <- anno_tbl |>
        left_join( volcano_plot_tab_cln
                   , by = join_by({{uniprot_column}} == {{uniprot_column}}) )  |>
        mutate( gene_name = case_when( is.na( gene_name) ~ {{uniprot_column}},
                                       TRUE ~ gene_name) )

      gene_names <- anno_tbl |>
        dplyr::pull(gene_name)

      rownames( r_obj@.Data[[1]] ) <- gene_names

      r_obj$p.value[,coef] <- qvalue( r_obj$p.value[,coef])$qvalues

      htmlwidgets::saveWidget( widget = Glimma::glimmaVolcano(r_obj
                                                      , coef=coef
                                                      , anno=anno_tbl
                                                      , counts = counts_tbl
                                                      , groups = groups
                                                      , display.columns = colnames(anno_tbl )
                                                      , status=decideTests(r_obj, adjust.method="none")
                                                      , p.adj.method = "none"
                                                      , transform.counts='none'
                                                      ) #the plotly object
                               , file = file.path( output_dir
                                                   , paste0(colnames(r_obj$coefficients)[coef], ".html"))  #the path & file name
                               , selfcontained = TRUE #creates a single html file
      )

    }

}


# ----------------------------------------------------------------------------
# getGlimmaVolcanoProteomicsWidget
# ----------------------------------------------------------------------------
#' @export
getGlimmaVolcanoProteomicsWidget <- function( r_obj
                                        , coef
                                        , volcano_plot_tab
                                        , uniprot_column = best_uniprot_acc
                                        , gene_name_column = gene_name
                                        , display_columns = c(  "PROTEIN_NAMES"   )
                                        , additional_annotations = NULL
                                        , additional_annotations_join_column = NULL
                                        , counts_tbl = NULL
                                        , groups = NULL) {

  message("--- Entering getGlimmaVolcanoProteomicsWidget ---")
  message(sprintf("   getGlimmaVolcanoProteomicsWidget Arg: coef = %d", coef))
  message(sprintf("   getGlimmaVolcanoProteomicsWidget Arg: r_obj class = %s", class(r_obj)))
  message(sprintf("   getGlimmaVolcanoProteomicsWidget Arg: volcano_plot_tab dims = %d rows, %d cols", nrow(volcano_plot_tab), ncol(volcano_plot_tab)))
  message("   getGlimmaVolcanoProteomicsWidget Arg: volcano_plot_tab structure:")
  utils::str(volcano_plot_tab)
  message("   getGlimmaVolcanoProteomicsWidget Arg: volcano_plot_tab head:")
  print(head(volcano_plot_tab))

  if( coef <= ncol(r_obj$coefficients )) {
    message("   getGlimmaVolcanoProteomicsWidget Step: Coefficient validation passed...")
    message(sprintf("      Data State: r_obj$coefficients has %d columns", ncol(r_obj$coefficients)))

    message("   getGlimmaVolcanoProteomicsWidget Step: Extracting best_uniprot_acc from rownames...")
    message("      Data State: r_obj@.Data[[1]] structure before extraction:")
    utils::str(r_obj@.Data[[1]])
    message(sprintf("      Data State: r_obj@.Data[[1]] dims = %d rows, %d cols", nrow(r_obj@.Data[[1]]), ncol(r_obj@.Data[[1]])))
    message("      Data State: rownames(r_obj@.Data[[1]]) head:")
    print(head(rownames(r_obj@.Data[[1]])))

    best_uniprot_acc <- str_split(rownames(r_obj@.Data[[1]]), " |:" ) |>
      purrr::map_chr(1)
    
    message("   getGlimmaVolcanoProteomicsWidget Step: best_uniprot_acc extraction completed.")
    message(sprintf("      Data State: best_uniprot_acc length = %d", length(best_uniprot_acc)))
    message("      Data State: best_uniprot_acc head:")
    print(head(best_uniprot_acc))

    # print(paste("nrow = ", nrow(r_obj@.Data[[1]])))
    # print(head(best_uniprot_acc))

    message("   getGlimmaVolcanoProteomicsWidget Step: Cleaning volcano_plot_tab...")
    message("      Data State: volcano_plot_tab before cleaning:")
    utils::str(volcano_plot_tab)
    
    volcano_plot_tab_cln <- volcano_plot_tab  |>
      dplyr::select ( {{uniprot_column}}
                       , {{gene_name_column}}, any_of( display_columns)) |>
      distinct()
    
    message("   getGlimmaVolcanoProteomicsWidget Step: volcano_plot_tab cleaning completed.")
    message(sprintf("      Data State: volcano_plot_tab_cln dims = %d rows, %d cols", nrow(volcano_plot_tab_cln), ncol(volcano_plot_tab_cln)))
    message("      Data State: volcano_plot_tab_cln head:")
    print(head(volcano_plot_tab_cln))

    # print (head( volcano_plot_tab_cln))

    if( !is.null( additional_annotations )
        & !is.null( additional_annotations_join_column ) ) {
      
      message("   getGlimmaVolcanoProteomicsWidget Step: Processing additional annotations...")
      message("      Data State: additional_annotations structure:")
      utils::str(additional_annotations)

      volcano_plot_tab_cln <- volcano_plot_tab_cln |>
        left_join( additional_annotations
                   , by = join_by( {{uniprot_column}} == {{additional_annotations_join_column}} ) ) |>
        dplyr::select( {{uniprot_column}}
                       , {{gene_name_column}}
                       , any_of( display_columns))
                       
      message("   getGlimmaVolcanoProteomicsWidget Step: Additional annotations processing completed.")
      message(sprintf("      Data State: volcano_plot_tab_cln after annotation dims = %d rows, %d cols", nrow(volcano_plot_tab_cln), ncol(volcano_plot_tab_cln)))
    } else {
      message("   getGlimmaVolcanoProteomicsWidget Step: No additional annotations to process.")
    }

    message("   getGlimmaVolcanoProteomicsWidget Step: Creating annotation table...")
    
    anno_tbl <- data.frame( uniprot_acc = rownames(r_obj@.Data[[1]]) # This uniprot_acc does not matter, only shows in glimma Volcano table
                            , temp_column = best_uniprot_acc ) |>
      dplyr::rename( {{uniprot_column}} := temp_column) |>
      left_join( volcano_plot_tab_cln
                 , by = join_by({{uniprot_column}} == {{uniprot_column}}) )  |>
      mutate( gene_name = case_when( is.na( gene_name) ~ {{uniprot_column}},
                                     TRUE ~ gene_name) )
                                     
    message("   getGlimmaVolcanoProteomicsWidget Step: Annotation table creation completed.")
    message(sprintf("      Data State: anno_tbl dims = %d rows, %d cols", nrow(anno_tbl), ncol(anno_tbl)))
    message("      Data State: anno_tbl structure:")
    utils::str(anno_tbl)
    message("      Data State: anno_tbl head:")
    print(head(anno_tbl))

    message("   getGlimmaVolcanoProteomicsWidget Step: Extracting gene names...")
    gene_names <- anno_tbl |>
      dplyr::pull(gene_name)
      
    message("   getGlimmaVolcanoProteomicsWidget Step: Gene names extraction completed.")
    message(sprintf("      Data State: gene_names length = %d", length(gene_names)))
    message("      Data State: gene_names head:")
    print(head(gene_names))

    message("   getGlimmaVolcanoProteomicsWidget Step: Updating rownames in r_obj...")
    message("      Data State: r_obj@.Data[[1]] rownames before update:")
    print(head(rownames(r_obj@.Data[[1]])))
    
    rownames( r_obj@.Data[[1]] ) <- gene_names
    
    message("   getGlimmaVolcanoProteomicsWidget Step: Rownames update completed.")
    message("      Data State: r_obj@.Data[[1]] rownames after update:")
    print(head(rownames(r_obj@.Data[[1]])))

    message("   getGlimmaVolcanoProteomicsWidget Step: Updating p-values with qvalue...")
    message(sprintf("      Data State: r_obj$p.value dimensions = %d rows, %d cols", nrow(r_obj$p.value), ncol(r_obj$p.value)))
    message(sprintf("      Data State: coef = %d, ncol(r_obj$p.value) = %d", coef, ncol(r_obj$p.value)))
    message("      Data State: r_obj$p.value[,coef] before qvalue transformation:")
    print(head(r_obj$p.value[,coef]))

    r_obj$p.value[,coef] <- qvalue( r_obj$p.value[,coef])$qvalues
    
    message("   getGlimmaVolcanoProteomicsWidget Step: P-value qvalue transformation completed.")
    message("      Data State: r_obj$p.value[,coef] after qvalue transformation:")
    print(head(r_obj$p.value[,coef]))

    message("   getGlimmaVolcanoProteomicsWidget Step: Calling glimmaVolcano...")
    message("      glimmaVolcano parameters:")
    message(sprintf("        coef = %d", coef))
    message(sprintf("        counts_tbl is.null = %s", is.null(counts_tbl)))
    if (!is.null(counts_tbl)) {
      message(sprintf("        counts_tbl dims = %d rows, %d cols", nrow(counts_tbl), ncol(counts_tbl)))
    }
    message(sprintf("        groups is.null = %s", is.null(groups)))
    if (!is.null(groups)) {
      message(sprintf("        groups length = %d", length(groups)))
      message("        groups head:")
      print(head(groups))
    }
    message(sprintf("        display_columns = %s", paste(display_columns, collapse = ", ")))

    result <- tryCatch({
      message("      About to call glimmaVolcano with parameters...")
      
      # Validate key parameters before calling
      message(sprintf("      Validation: r_obj class = %s", class(r_obj)))
      message(sprintf("      Validation: coef = %d", coef))
      message(sprintf("      Validation: ncol(r_obj$coefficients) = %d", ncol(r_obj$coefficients)))
      message(sprintf("      Validation: anno_tbl nrow = %d", nrow(anno_tbl)))
      message(sprintf("      Validation: display_columns = %s", paste(display_columns, collapse = ", ")))
      
      # Check if decideTests works and clean up any issues
      message("      Testing decideTests...")
      status_result <- decideTests(r_obj, adjust.method="none")
      message(sprintf("      decideTests completed: %d x %d", nrow(status_result), ncol(status_result)))
      
      # CRITICAL FIX: Handle NA values that cause dimension mismatches in glimmaVolcano
      message("      Checking for NA values in status_result...")
      na_count <- sum(is.na(status_result))
      message(sprintf("      Found %d NA values in status_result", na_count))
      
      if (na_count > 0) {
        message("      Cleaning NA values from status_result (setting to 0 = nonDE)...")
        status_result[is.na(status_result)] <- 0
        message("      NA values cleaned")
      }
      
      # Additional safety: Check for dimension alignment
      message("      Verifying dimensions before glimmaVolcano call...")
      message(sprintf("      r_obj coefficients: %d x %d", nrow(r_obj$coefficients), ncol(r_obj$coefficients)))
      message(sprintf("      r_obj p.value: %d x %d", nrow(r_obj$p.value), ncol(r_obj$p.value)))
      message(sprintf("      status_result: %d x %d", nrow(status_result), ncol(status_result)))
      message(sprintf("      anno_tbl: %d rows", nrow(anno_tbl)))
      message(sprintf("      counts_tbl: %d x %d", nrow(counts_tbl), ncol(counts_tbl)))
      
      message("      Now calling glimmaVolcano...")
      widget <- glimmaVolcano(r_obj
                     , coef=coef
                     , counts = counts_tbl
                     , groups = groups
                     , anno=anno_tbl
                     , display.columns = display_columns
                     , status=status_result
                     , p.adj.method="none"
                     , transform.counts='none') #the plotly object
      
      message("      glimmaVolcano call returned successfully!")
      message(sprintf("      Widget class: %s", class(widget)))
      message(sprintf("      Widget is.null: %s", is.null(widget)))
      
      return(widget)
      
    }, error = function(e) {
      message(sprintf("      ERROR in glimmaVolcano call: %s", e$message))
      message("      Full error details:")
      message(capture.output(print(e)))
      message("      Traceback:")
      message(capture.output(traceback()))
      return(NULL)
    })
                   
    message("   getGlimmaVolcanoProteomicsWidget Step: glimmaVolcano call completed successfully.")
    message("      Data State: result class:")
    print(class(result))
    
    message("--- Exiting getGlimmaVolcanoProteomicsWidget ---")
    return(result)

  } else {
    message(sprintf("   getGlimmaVolcanoProteomicsWidget Condition FALSE: coef (%d) > ncol(r_obj$coefficients) (%d)", coef, ncol(r_obj$coefficients)))
    message("--- Exiting getGlimmaVolcanoProteomicsWidget (early exit) ---")
    return(NULL)
  }

}


# ----------------------------------------------------------------------------
# getGlimmaVolcanoPhosphoproteomics
# ----------------------------------------------------------------------------
getGlimmaVolcanoPhosphoproteomics <- function( r_obj
                                        , coef
                                        , volcano_plot_tab
                                        , sites_id_column = sites_id
                                        , sites_id_display_column = sites_id_short
                                        , display_columns = c(  "sequence", "PROTEIN_NAMES"   )
                                        , additional_annotations = NULL
                                        , additional_annotations_join_column = NULL
                                        , counts_tbl = NULL
                                        , output_dir) {

  if( coef <= ncol(r_obj$coefficients )) {

    volcano_plot_tab_cln <- volcano_plot_tab |>
      dplyr::distinct( {{sites_id_column}}
                       , {{sites_id_display_column}}
                       , any_of( display_columns) )

    if( !is.null( additional_annotations )
        & !is.null( additional_annotations_join_column ) ) {

      volcano_plot_tab_cln <- volcano_plot_tab_cln |>
        left_join( additional_annotations
                   , by = join_by( {{sites_id_column}} == {{additional_annotations_join_column}} ) ) |>
        dplyr::select( any_of( display_columns))
    }

    anno_tbl <-  data.frame(  sites_id = rownames(r_obj@.Data[[1]])) |>
      left_join( volcano_plot_tab_cln
                 , by = join_by(sites_id == {{sites_id_column}} ) )

    sites_id_short_list <- anno_tbl |>
                   dplyr::pull(sites_id_short)

    rownames( r_obj@.Data[[1]] ) <- sites_id_short_list

    #coef <- seq_len( ncol(r_obj$coefficients))[1]

    r_obj$p.value[,coef] <- qvalue( r_obj$p.value[,coef])$qvalues


    htmlwidgets::saveWidget( widget = glimmaVolcano(r_obj, coef=coef
                                                    , counts = counts_tbl
                                                    , anno=anno_tbl
                                                    , display.columns=display_columns
                                                    , p.adj.method = "none"
                                                    , transform.counts='none' ) #the plotly object
                             , file = file.path( output_dir
                                                 , paste0(colnames(r_obj$coefficients)[coef], ".html"))  #the path & file name
                             , selfcontained = TRUE #creates a single html file
    )
  }

}


# ----------------------------------------------------------------------------
# printPValuesDistribution
# ----------------------------------------------------------------------------
#' Draw the p-values distribution plot.
#' @param selected_data A table that is generated by running the function \code{\link{get_significant_data}}.
#' @param log_p_value_column The name of the column representing the p-value.
#' @param formula_string The formula string used in the facet_grid command for the ggplot scatter plot.
#'@export
printPValuesDistribution <- function(selected_data, p_value_column = raw_pvalue, formula_string = "is_ruv_applied ~ comparison") {

  breaks <- c(0, 0.001, 0.01, 0.05,
              seq(0.1, 1, by = 0.1))

# after_stat(density)
  pvalhist <- ggplot(selected_data, aes({ { p_value_column } })) +
    theme(axis.title.y = element_blank()) +
    xlab("P-value") +
    geom_histogram(aes(y = after_stat(density)),
                   breaks = breaks,
                   position = "identity",
                   color = "black") +
    geom_histogram(aes(y = after_stat(density)),
                   breaks = breaks,
                   position = "identity")

  if (!is.na(formula_string)) {

    pvalhist <- pvalhist +
      facet_grid(as.formula(formula_string))
  }


  pvalhist

}


# ----------------------------------------------------------------------------
# gg_save_logging
# ----------------------------------------------------------------------------
#' @export
gg_save_logging <- function( input_plot
                             , file_name_part
                             , plots_format
                             , width=7
                             , height=7) {
  for( format_ext in plots_format) {
    file_name <- paste0(file_name_part, format_ext)
    captured_output<-capture.output(
      ggsave(plot=input_plot
             , filename = file_name
             , width=width
             , height=height )
      ,type = "message"
    )
    logdebug(captured_output)
  }
}


# ----------------------------------------------------------------------------
# summarizeQCPlot
# ----------------------------------------------------------------------------
summarizeQCPlot <- function(qc_figure) {
            cat("RLE Plots:\n")
            for (plot_name in names(qc_figure@rle_plots)) {
              cat(paste(" -", plot_name, "\n"))
              print(qc_figure@rle_plots[[plot_name]])
            }

            cat("\nPCA Plots:\n")
            for (plot_name in names(qc_figure@pca_plots)) {
              cat(paste(" -", plot_name, "\n"))
              print(qc_figure@pca_plots[[plot_name]])
            }

            cat("\nDensity Plots:\n")
            for (plot_name in names(qc_figure@density_plots)) {
              cat(paste(" -", plot_name, "\n"))
              print(qc_figure@density_plots[[plot_name]])
            }

            cat("\nPearson Correlation Plots:\n")
            for (plot_name in names(qc_figure@pearson_plots)) {
              cat(paste(" -", plot_name, "\n"))
              print(qc_figure@pearson_plots[[plot_name]])
            }
}


# ----------------------------------------------------------------------------
# plotVolcano
# ----------------------------------------------------------------------------
#' Draw the volcano plot.
#' @param selected_data A table that is generated by running the function \code{\link{get_significant_data}}.
#' @param log_q_value_column The name of the column representing the log q-value.
#' @param log_fc_column The name of the column representing the log fold-change.
#' @param q_val_thresh A numerical value specifying the q-value threshold for statistically significant proteins.
#' @param formula_string The formula string used in the facet_grid command for the ggplot scatter plot.
#' @export
plotVolcano <- function( selected_data
                         , log_q_value_column = lqm
                         , log_fc_column = logFC
                         , q_val_thresh = 0.05
                         , formula_string = "analysis_type ~ comparison" ) {

  volplot_gg.all <- selected_data |>
    ggplot(aes(y = { { log_q_value_column } }, x = { { log_fc_column } })) +
    geom_point(aes(col = colour)) +
    scale_colour_manual(values = c(levels(selected_data$colour)),
                        labels = c(paste0("Not significant, logFC > ",
                                          1),
                                   paste0("Significant, logFC >= ",
                                          1),
                                   paste0("Significant, logFC <",
                                          1),
                                   "Not Significant")) +
    geom_vline(xintercept = 1, colour = "black", linewidth = 0.2) +
    geom_vline(xintercept = -1, colour = "black", linewidth = 0.2) +
    geom_hline(yintercept = -log10(q_val_thresh)) +
    theme_bw() +
    xlab("Log fold changes") +
    ylab("-log10 q-value") +
    theme(legend.position = "none")

  volplot_gg.plot <- volplot_gg.all
  if( !is.na(formula_string) | formula_string != "" ) {
    volplot_gg.plot <- volplot_gg.all +
      facet_grid( as.formula(formula_string),
                  labeller = labeller(facet_category = label_wrap_gen(width = 10)))
  }

  volplot_gg.plot
}


# ----------------------------------------------------------------------------
# plotPca
# ----------------------------------------------------------------------------
plotPca <- function(theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size=8, cv_percentile = 0.90) {
  # Defensive checks
  if (!is.character(grouping_variable) || length(grouping_variable) != 1) {
    stop("grouping_variable must be a single character string")
  }
  
  if (!is.null(shape_variable) && (!is.character(shape_variable) || length(shape_variable) != 1)) {
    stop("shape_variable must be NULL or a single character string")
  }
  
  if (!grouping_variable %in% colnames(theObject@design_matrix)) {
    stop(sprintf("grouping_variable '%s' not found in design matrix", grouping_variable))
  }
  
  if (!is.null(shape_variable) && !shape_variable %in% colnames(theObject@design_matrix)) {
    stop(sprintf("shape_variable '%s' not found in design matrix", shape_variable))
  }
  
  data_matrix <- NULL
  ## I want to check the class of theObject here
 if( class(theObject) == "PeptideQuantitativeData") {
   data_matrix <- theObject@peptide_matrix
   
 } else if( class(theObject) == "ProteinQuantitativeData") {
   data_matrix <- theObject@protein_quant_table |>
     column_to_rownames(var = "Protein.Ids") |>
     as.matrix()
 }
  
  design_matrix <- theObject@design_matrix
  sample_id <- theObject@sample_id
  
  # Prepare matrix for PCA (data should already be log2 transformed)
  data_matrix_pca <- data_matrix
  data_matrix_pca[!is.finite(data_matrix_pca)] <- NA
  
  if(is.na(label_column) || label_column == "") {
    label_column <- ""
  }
  
  pca_plot <- plotPcaHelper(data_matrix_pca
                            , design_matrix = design_matrix
                            , sample_id_column = sample_id
                            , grouping_variable = grouping_variable
                            , shape_variable = shape_variable
                            , label_column = label_column
                            , title = title
                            , geom.text.size = font_size
                            , cv_percentile = cv_percentile)
  
  return(pca_plot)
}

