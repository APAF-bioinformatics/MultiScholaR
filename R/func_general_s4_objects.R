# ============================================================================
# func_general_s4_objects.R
# ============================================================================
# Purpose: General/shared S4 object classes and utilities
# 
# This file contains S4 class definitions and utilities that are shared across
# multiple omics types or are general-purpose. This includes FilteringProgress
# classes and other shared S4 utilities.
#
# NOTE: Omic-specific S4 classes should go in:
# - func_prot_s4_objects.R (ProteinQuantitativeData, PeptideQuantitativeData)
# - func_metab_s4_objects.R (MetaboliteAssayData)
# - func_lipid_s4_objects.R (LipidAssayData - placeholder)
#
# Functions to extract here:
# - FilteringProgress class definition (shared across omics)
# - FilteringProgressMetabolomics class definition
# - Shared S4 utility methods
# - S4 object validation helpers (if shared)
#
# Dependencies:
# - methods package
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# === Shared S4 Class Definitions ===

# Function 1: FilteringProgress class definition
# Current location: R/qc_and_rollup.R
# Description: S4 class for filtering progress tracking (used by proteomics)
# setClass("FilteringProgress", ...) {
#   # Extract from R/qc_and_rollup.R
# }

# Function 2: FilteringProgressMetabolomics class definition
# Current location: R/qc_and_rollup.R
# Description: S4 class for metabolomics filtering progress
# setClass("FilteringProgressMetabolomics", ...) {
#   # Extract from R/qc_and_rollup.R
# }

# === Shared S4 Utilities ===

# Function 3: Shared S4 utility methods
# Current location: Various files
# Description: Any S4 methods that are shared across multiple omics types
# # Extract shared methods if any exist


# ----------------------------------------------------------------------------
# FilteringProgress
# ----------------------------------------------------------------------------
#' FilteringProgress Class
#' 
#' @description
#' An S4 class to track and store the progress of protein filtering steps in
#' proteomics data analysis. This class maintains records of protein and peptide
#' counts at each filtering stage.
#' 
#' @slot steps Character vector storing names of filtering steps
#' @slot proteins Numeric vector storing protein counts for each step
#' @slot total_peptides Numeric vector storing total peptide counts for each step
#' @slot peptides_per_protein List storing peptides per protein distributions for each step
#' @slot proteins_per_run List storing proteins per run counts for each step
#' @slot peptides_per_run List storing peptides per run counts for each step
#' 
#' @export
setClass("FilteringProgress",
  slots = list(
    steps = "character",
    proteins = "numeric",
    total_peptides = "numeric",
    peptides_per_protein = "list",
    proteins_per_run = "list",
    peptides_per_run = "list"
  )
)

##################################################################################################################

#' Initialize a new FilteringProgress object
#' 
#' @description
#' Creates a new FilteringProgress object with empty slots to track protein
#' filtering progress.
#' 
#' @return A new FilteringProgress object
#' 
#' @examples
#' filtering_progress <- new("FilteringProgress",
#'   steps = character(),
#'   proteins = numeric(),
#'   total_peptides = numeric(),
#'   peptides_per_protein = list(),
#'   proteins_per_run = list(),
#'   peptides_per_run = list()
#' )
#' 
#' @export
filtering_progress <- new("FilteringProgress",
  steps = character(),
  proteins = numeric(),
  total_peptides = numeric(),
  peptides_per_protein = list(),
  proteins_per_run = list(),
  peptides_per_run = list()
)

##################################################################################################################

#' Generate a color palette
#' 
#' @param n Number of colors needed
#' @param base_color Base color to use
#' @return Vector of colors
#' @export 
get_color_palette <- function(n, base_color) {
  colorRampPalette(c(base_color, "black"))(n)
}

# ==========================================
# DirectoryManager Class (from helper_functions.R)
# ==========================================

#' DirectoryManager S4 Class
#' 
#' @description
#' An S4 class to manage project directory paths for MultiScholaR analyses.
#' This class stores paths to various output directories used during analysis.
#' 
#' @slot base_dir Character string for base project directory
#' @slot results_dir Character string for results directory
#' @slot data_dir Character string for data directory
#' @slot source_dir Character string for source files directory
#' @slot de_output_dir Character string for differential expression output
#' @slot publication_graphs_dir Character string for publication-ready graphs
#' @slot timestamp Character string for timestamp used in file naming
#' @slot qc_dir Character string for QC output directory
#' @slot time_dir Character string for time-stamped directory
#' @slot results_summary_dir Character string for results summary directory
#' @slot pathway_dir Character string for pathway analysis output
#' 
#' @export
setClass("DirectoryManager",
    slots = c(
        base_dir = "character",
        results_dir = "character",
        data_dir = "character",
        source_dir = "character",
        de_output_dir = "character",
        publication_graphs_dir = "character",
        timestamp = "character",
        qc_dir = "character",
        time_dir = "character",
        results_summary_dir = "character",
        pathway_dir = "character"
    )
)

# ==========================================
# Enrichment S4 Classes (from functional_enrichment.R)
# ==========================================

#' de_results_for_enrichment S4 Class
#' 
#' @description
#' An S4 class to store differential expression results formatted for
#' enrichment analysis. Contains contrast definitions, DE data, and
#' experimental design information.
#' 
#' @slot contrasts A tibble containing contrast information
#' @slot de_data A list of DE results data frames
#' @slot design_matrix A data frame containing the design matrix
#' 
#' @export
setClass("de_results_for_enrichment",
         slots = list(
           contrasts = "tbl_df",
           de_data = "list",
           design_matrix = "data.frame"
         ))

#' Create DE Results For Enrichment
#'
#' @param contrasts_tbl A tibble containing contrast information
#' @param design_matrix A data frame containing the design matrix
#' @param de_output_dir Directory containing DE results files
#' @return An S4 object of class de_results_for_enrichment
#' @export
createDEResultsForEnrichment <- function(contrasts_tbl, design_matrix, de_output_dir) {
  # Helper function to format contrast filename
  format_contrast_filename <- function(contrast_string) {
    contrast_name <- stringr::str_split(contrast_string, "=")[[1]][1] |>
      stringr::str_replace_all("\\.", "_")

    paste0("de_proteins_", contrast_name, "_long_annot.tsv")
  }

  # Create new S4 object
  de_results <- new("de_results_for_enrichment")

  # Convert contrasts_tbl to tibble if it isn't already
  contrasts_tbl <- tibble::as_tibble(contrasts_tbl)

  # Fill slots
  de_results@contrasts <- contrasts_tbl
  de_results@design_matrix <- design_matrix
  de_results@de_data <- contrasts_tbl$contrasts |>
    purrr::set_names() |>
    purrr::map(function(contrast) {
      filename <- format_contrast_filename(contrast)
      filepath <- file.path(de_output_dir, filename)

      if (!file.exists(filepath)) {
        warning("File not found: ", filepath)
        return(NULL)
      }

      readr::read_tsv(filepath, show_col_types = FALSE)
    })

  return(de_results)
}

#' EnrichmentResults S4 Class
#' 
#' @description
#' An S4 class to store enrichment analysis results including
#' enrichment data, plots, and summaries.
#' 
#' @slot contrasts A tibble containing contrast information
#' @slot enrichment_data A list of enrichment results
#' @slot enrichment_plots A list of gostplot objects
#' @slot enrichment_plotly A list of interactive plotly objects
#' @slot enrichment_summaries A list of summary data
#' 
#' @export
setClass("EnrichmentResults",
         slots = list(
           contrasts = "tbl_df",
           enrichment_data = "list",
           enrichment_plots = "list",
           enrichment_plotly = "list",
           enrichment_summaries = "list"
         ))

#' Create EnrichmentResults Object
#' 
#' @param contrasts_tbl A tibble containing contrast information
#' @return An S4 object of class EnrichmentResults
#' @export
createEnrichmentResults <- function(contrasts_tbl) {
  new("EnrichmentResults",
      contrasts = contrasts_tbl,
      enrichment_data = list(),
      enrichment_plots = list(),
      enrichment_plotly = list(),
      enrichment_summaries = list())
}

