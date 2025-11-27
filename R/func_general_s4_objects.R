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

