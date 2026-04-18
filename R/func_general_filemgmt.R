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
# func_general_filemgmt.R
# ============================================================================
# Purpose: File management and directory setup functions
#
# This file contains functions for file management, directory creation,
# configuration file handling, and project path management. Functions in
# this file are used across all omics types and workflows.
#
# Functions to extract here:
# - setupDirectories(): Main function for setting up project directories
# - createDirIfNotExists(): Creates directory if it doesn't exist
# - createDirectoryIfNotExists(): Alternative directory creation function
# - createOutputDir(): Creates output directory
# - setupAndShowDirectories(): Sets up and shows directory structure
# - getProjectPaths(): Gets project paths
# - readConfigFile(): Reads configuration file
# - readConfigFileSection(): Reads section from config file
# - formatConfigList(): Formats config list
# - updateConfigParameter(): Updates config parameter
# - createStudyParametersFile(): Creates study parameters file
# - createWorkflowArgsFromConfig(): Creates workflow args from config
# - getDefaultProteomicsConfig(): Gets default proteomics config
# - copyToResultsSummary(): Copies files to results summary
# - write_results(): Writes results to file
# - pushProjectToGithub(): Pushes project to GitHub
# - pushProjectToGithubFromDirs(): Pushes project from directories
# - sourceRmdFile(): Sources R Markdown file
# - sourceRmdFileSimple(): Simple R Markdown sourcing
# - Additional file management helper functions
#
# Dependencies:
# - fs, here, configr, ini
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: setupDirectories()
# Current location: R/file_management.R
# Description: Sets up project directories for omics analysis
# setupDirectories <- function(...) {
#   # Extract from R/file_management.R
# }

# Function 2: createDirIfNotExists()
# Current location: R/helper_functions.R
# Description: Creates directory if it doesn't exist
# createDirIfNotExists <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 3: createDirectoryIfNotExists()
# Current location: R/helper_functions.R
# Description: Alternative directory creation function
# createDirectoryIfNotExists <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 4: createOutputDir()
# Current location: R/helper_functions.R
# Description: Creates output directory
# createOutputDir <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 5: setupAndShowDirectories()
# Current location: R/helper_functions.R
# Description: Sets up and shows directory structure
# setupAndShowDirectories <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 6: getProjectPaths()
# Current location: R/helper_functions.R
# Description: Gets project paths
# getProjectPaths <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 7: readConfigFile()
# Current location: R/helper_functions.R
# Description: Reads configuration file
# readConfigFile <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 8: readConfigFileSection()
# Current location: R/helper_functions.R
# Description: Reads section from config file
# readConfigFileSection <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 9: formatConfigList()
# Current location: R/helper_functions.R
# Description: Formats config list
# formatConfigList <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 10: updateConfigParameter()
# Current location: R/helper_functions.R
# Description: Updates config parameter
# updateConfigParameter <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 11: createStudyParametersFile()
# Current location: R/helper_functions.R
# Description: Creates study parameters file
# createStudyParametersFile <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 12: createWorkflowArgsFromConfig()
# Current location: R/helper_functions.R
# Description: Creates workflow args from config
# createWorkflowArgsFromConfig <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 13: getDefaultProteomicsConfig()
# Current location: R/mod_prot_import.R
# Description: Gets default proteomics configuration
# getDefaultProteomicsConfig <- function() {
#   # Extract from R/mod_prot_import.R
# }

# Function 14: copyToResultsSummary()
# Current location: R/helper_functions.R
# Description: Copies files to results summary directory
# copyToResultsSummary <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 15: write_results()
# Current location: R/helper_functions.R
# Description: Writes results to file
# write_results <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 16: pushProjectToGithub()
# Current location: R/github_managment.R
# Description: Pushes project to GitHub
# pushProjectToGithub <- function(...) {
#   # Extract from R/github_managment.R
# }

# Function 17: pushProjectToGithubFromDirs()
# Current location: R/github_managment.R
# Description: Pushes project from directories
# pushProjectToGithubFromDirs <- function(...) {
#   # Extract from R/github_managment.R
# }

# Function 18: sourceRmdFile()
# Current location: R/helper_functions.R
# Description: Sources R Markdown file
# sourceRmdFile <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 19: sourceRmdFileSimple()
# Current location: R/helper_functions.R
# Description: Simple R Markdown sourcing
# sourceRmdFileSimple <- function(...) {
#   # Extract from R/helper_functions.R
# }












# ----------------------------------------------------------------------------
# setupDirectories
# ----------------------------------------------------------------------------
#' @title Setup directories
#' @description Creates and manages project directories with version control for multiple omic types.
#' @param base_dir Base directory path (optional, defaults to here::here())
#' @param omic_types Character vector or comma/space-separated string indicating the types of omics data (e.g., c("proteomics", "metabolomics"), "proteomics, transcriptomics").
#' @param label Optional label to append to each omic type directory name (e.g., "proteomics_MyLabel").
#' @param force Logical; if TRUE, skips user confirmation for overwriting existing directories (default: FALSE).
#' @param reuse_existing Logical; if TRUE, uses an existing directory structure without prompting or overwriting (default: FALSE). For non-interactive use.
#' @return A named list where each name is an omic type. Each element is a list containing the relevant directory paths for that omic type.
#' @importFrom logger log_info log_warn log_error log_debug
#' @importFrom stringr str_split str_trim str_detect str_remove str_replace str_ends
#' @importFrom here here
#' @importFrom purrr map map_chr map_dbl imap imap_chr iwalk walk walk2
#' @importFrom dplyr filter select mutate pull first distinct summarise group_by arrange
#' @importFrom readr read_tsv write_tsv write_lines
#' @importFrom vroom vroom vroom_write
#' @importFrom tools file_path_as_absolute file_path_sans_ext toTitleCase
#' @importFrom rlang abort
#' @importFrom tibble is_tibble
#' @importFrom knitr purl
#' @importFrom rmarkdown render
#' @importFrom openxlsx createWorkbook addWorksheet read.xlsx writeData setColWidths addStyle createStyle writeDataTable saveWorkbook
#' @importFrom gh gh
#' @importFrom rstudioapi isAvailable getActiveDocumentContext documentSave
#' @importFrom methods slotNames
#' @importFrom utils head
#' @importFrom rappdirs user_cache_dir
#' @export
setupDirectories <- function(base_dir = here::here(), omic_types, label = NULL, force = FALSE, reuse_existing = FALSE) {
    parsed_omic_types <- parseSetupDirectoriesOmicTypes(omic_types)

    # --- Initialization ---
    all_created_paths <- list() # List to store paths for each omic type
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S") # Single timestamp for the run

    # --- Loop Through Each Omic Type ---
    for (current_omic_type in parsed_omic_types) {
        logger::log_info("--- Processing setup for omic type: {current_omic_type} ---")

        # Fetch Configuration for the current omic type
        omic_config <- getSetupDirectoriesOmicConfig(current_omic_type)

        # Construct Directory Names for the current omic type
        omic_label_dirname <- paste0(current_omic_type, if (!is.null(label)) paste0("_", substr(label, 1, 30)) else "")

        results_path <- file.path(base_dir, "results", omic_label_dirname)
        results_summary_path <- file.path(base_dir, "results_summary", omic_label_dirname)
        scripts_path <- file.path(base_dir, "scripts", omic_label_dirname) # Destination scripts path for this omic

        directory_handling <- handleSetupDirectoriesExistingDirs(
            current_omic_type = current_omic_type,
            omic_label_dirname = omic_label_dirname,
            results_path = results_path,
            results_summary_path = results_summary_path,
            scripts_path = scripts_path,
            force = force,
            reuse_existing = reuse_existing
        )
        process_current_omic <- directory_handling$process_current_omic
        reuse_current_omic_dirs <- directory_handling$reuse_current_omic_dirs

        if (!process_current_omic) {
            next # Skip to the next omic type in the loop
        }

        current_omic_paths_def <- initializeSetupDirectoriesOmicPaths(
            base_dir = base_dir,
            current_omic_type = current_omic_type,
            results_path = results_path,
            results_summary_path = results_summary_path,
            scripts_path = scripts_path,
            omic_config = omic_config,
            timestamp = timestamp
        )

        materializeSetupDirectoriesStructure(
            current_omic_paths_def = current_omic_paths_def,
            omic_config = omic_config,
            current_omic_type = current_omic_type,
            omic_label_dirname = omic_label_dirname,
            reuse_current_omic_dirs = reuse_current_omic_dirs
        )


        # --- Build and Store Path Variables for the current omic type ---
        omic_specific_paths_list <- buildSetupDirectoriesPathList(
            base_dir = base_dir,
            current_omic_type = current_omic_type,
            omic_label_dirname = omic_label_dirname,
            timestamp = timestamp,
            current_omic_paths_def = current_omic_paths_def,
            omic_config = omic_config
        )

        all_created_paths[[omic_label_dirname]] <- omic_specific_paths_list # Use omic_label_dirname as key
        logger::log_info("Stored paths for {omic_label_dirname}.")
    } # End loop over omic types

    printSetupDirectoriesSummary(all_created_paths)

    # Return the collected paths invisibly
    invisible(all_created_paths)
}








































##################################################################################################################
#' @title Load MultiScholaR Dependencies
#'
#' @description
#' Installs and loads all required packages for MultiScholaR This includes packages from CRAN, Bioconductor, and GitHub.
#'
#' @param verbose logical; if TRUE (default), displays progress messages during
#'   package installation and loading
#'
#' @details
#' Checks for and installs missing packages, then loads all required
#' dependencies. It handles special cases for GitHub packages (RUVIIIC and
#' MultiScholaR) and ensures all necessary packages are available for the
#' workflows.
#'
#' @return None (called for side effects)
#'
#' @examples
#' \dontrun{
#' # Load with default verbose messaging
#' loadDependencies()
#'
#' # Load silently
#' loadDependencies(verbose = FALSE)
#' }
#'
#' @importFrom devtools install_github
#' @importFrom pacman p_load
#'
#' @export
loadDependencies <- function(verbose = TRUE) {
    # --- Install Core Managers ---
    if (!requireNamespace("pacman", quietly = TRUE)) {
        if (verbose) message("Installing pacman...")
        utils::install.packages("pacman")
    }
    library(pacman)

    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        if (verbose) message("Installing BiocManager...")
        utils::install.packages("BiocManager")
    }
    # Ensure BiocManager is loaded for installation checks, but don't need library()
    # library(BiocManager) # Not strictly needed just for install()

    if (!requireNamespace("devtools", quietly = TRUE)) {
        if (verbose) message("Installing devtools...")
        utils::install.packages("devtools")
    }
    # library(devtools) # Not strictly needed just for install_github()

    # --- Define Packages by Source ---
    cran_packages <- c(
        "tidyverse", "seqinr", "lazyeval", "rlang", "glue", "GGally",
        "here", "tibble", "magrittr", "future.apply", "tictoc",
        "beepr", "furrr", "readxl", "writexl", "RColorBrewer",
        "multidplyr", "RSpectra", "progress", "Rcpp", "RcppEigen",
        "ruv", "iq", "ggrepel", "patchwork", "dplyr", "gtools",
        "shiny", "DT", "gh", "openxlsx", "plotly", "vroom",
        "gplots", "iheatmapr", "UpSetR", "gt", "gprofiler2",
        "htmltools", "rstudioapi", "flextable", "viridis", "here",
        "git2r", "fs", "logger", "httr", "jsonlite",
        "configr", "webshot2", "shiny", "shinyjs", "shinyWidgets",
        "shinydashboard", "shinythemes", "shinycssloaders", "golem",

        # Added for BookChapter
        # Note: "conflicted" removed from auto-loading - it causes hard errors on namespace
        # conflicts during pacman::p_load(), creating an install loop. Users can still
        # use conflicted manually in RMarkdown workflows with explicit conflict_prefer() calls.
        "tidytext",
        # Added from Suggests:
        "testthat", "ggplot2", "ggpubr", "svglite",
        "ggraph", "reticulate", "shinyFiles", "arrow"
    )

    bioc_packages <- c(
        "UniProt.ws", "mixOmics", "limma", "qvalue",
        "clusterProfiler", "GO.db", # GO.db is often a dependency, ensure it's listed
        "basilisk", "limpa", "ComplexHeatmap"
    )

    github_packages <- list(
        RUVIIIC = "cran/RUVIIIC", # Hosted on CRAN's GitHub mirror? Check source if issues
        Glimma = "APAF-bioinformatics/GlimmaV2"
    )

    # --- Installation and Loading Logic ---

    # Helper function to install/load
    install_and_load <- function(pkg, installer_func, source_name, verbose) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            if (verbose) message(sprintf("Installing %s from %s...", pkg, source_name))
            tryCatch(
                {
                    installer_func(pkg)
                    # After install, load it
                    pacman::p_load(char = pkg, character.only = TRUE)
                },
                error = function(e) {
                    warning(sprintf("Failed to install %s from %s: %s", pkg, source_name, e$message))
                }
            )
        } else {
            if (verbose) message(sprintf("%s is already installed, loading...", pkg))
            pacman::p_load(char = pkg, character.only = TRUE)
        }
    }

    # CRAN Packages
    if (verbose) message("\n--- Processing CRAN Packages ---")
    purrr::walk(cran_packages, ~ install_and_load(
        pkg = .x,
        # Use base R install.packages directly
        installer_func = function(p) utils::install.packages(p, dependencies = TRUE),
        source_name = "CRAN",
        verbose = verbose
    ))

    # Bioconductor Packages
    if (verbose) message("\n--- Processing Bioconductor Packages ---")
    purrr::walk(bioc_packages, ~ install_and_load(
        pkg = .x,
        installer_func = function(p) BiocManager::install(p, update = FALSE, ask = FALSE), # Use BiocManager
        source_name = "Bioconductor",
        verbose = verbose
    ))

    # GitHub Packages
    if (verbose) message("\n--- Processing GitHub Packages ---")
    purrr::iwalk(github_packages, ~ {
        pkg_name <- .y # Name of the package (e.g., "RUVIIIC")
        repo <- .x # Repository path (e.g., "cran/RUVIIIC")
        if (!requireNamespace(pkg_name, quietly = TRUE)) {
            if (verbose) message(sprintf("Installing %s from GitHub (%s)...", pkg_name, repo))
            tryCatch(
                {
                    # Force installation to handle potentially corrupt states
                    devtools::install_github(repo, force = TRUE)
                    pacman::p_load(char = pkg_name, character.only = TRUE)
                },
                error = function(e) {
                    warning(sprintf("Failed to install %s from GitHub (%s): %s", pkg_name, repo, e$message))
                }
            )
        } else {
            if (verbose) message(sprintf("%s is already installed, loading...", pkg_name))
            pacman::p_load(char = pkg_name, character.only = TRUE)
        }
    })

    if (verbose) message("\nAll dependencies processed successfully!")
}
