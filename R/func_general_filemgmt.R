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
    # --- Input Parsing and Validation ---

    # Handle comma/space separated string or vector input
    if (is.character(omic_types) && length(omic_types) == 1) {
        # Split by comma or space, trim whitespace, remove empty strings
        parsed_omic_types <- omic_types |>
            stringr::str_split("[,\\s]+") |> # Split by comma or one or more spaces
            unlist() |>
            stringr::str_trim()
        parsed_omic_types <- parsed_omic_types[parsed_omic_types != ""]
    } else if (is.character(omic_types)) {
        parsed_omic_types <- stringr::str_trim(omic_types)
        parsed_omic_types <- parsed_omic_types[parsed_omic_types != ""]
    } else {
        logger::log_error("`omic_types` must be a character vector or a single comma/space-separated string.")
        stop("`omic_types` must be a character vector or a single comma/space-separated string.")
    }

    if (length(parsed_omic_types) == 0) {
        logger::log_error("No valid `omic_types` provided after parsing.")
        stop("No valid `omic_types` provided.")
    }

    # Define Omics-Specific Configuration
    valid_omic_types <- c("proteomics", "metabolomics", "transcriptomics", "lipidomics", "integration") # Expanded

    # Validate all provided types
    invalid_types <- setdiff(parsed_omic_types, valid_omic_types)
    if (length(invalid_types) > 0) {
        err_msg <- sprintf("Invalid omic_type(s) specified: %s. Choose from: %s.",
                           paste(invalid_types, collapse = ", "),
                           paste(valid_omic_types, collapse = ", "))
        logger::log_error(err_msg)
        stop(err_msg)
    }

    # --- Initialization ---
    all_created_paths <- list() # List to store paths for each omic type
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S") # Single timestamp for the run

    # --- Loop Through Each Omic Type ---
    for (current_omic_type in parsed_omic_types) {
        logger::log_info("--- Processing setup for omic type: {current_omic_type} ---")

        # Fetch Configuration for the current omic type
        omic_config <- switch(current_omic_type,
            proteomics = list(
                results_subdirs = c(
                    "protein_qc", "peptide_qc", "clean_proteins", "de_proteins",
                    "publication_graphs", "pathway_enrichment",
                    file.path("publication_graphs", "filtering_qc")
                ),
                results_summary_subdirs = c("QC_figures", "Publication_figures", "Publication_tables", "Study_report"),
                scripts_source_leaf = "proteomics",
                global_vars = list(
                    de_output_leaf = "de_proteins",
                    pathway_leaf = "pathway_enrichment",
                    feature_qc_leaf = "protein_qc",
                    subfeature_qc_leaf = "peptide_qc",
                    clean_features_leaf = "clean_proteins"
                )
            ),
            metabolomics = list(
                results_subdirs = c(
                    "metabolite_qc", "feature_qc", "clean_metabolites", "de_metabolites",
                    "publication_graphs", "pathway_enrichment",
                    file.path("publication_graphs", "filtering_qc")
                ),
                results_summary_subdirs = c("QC_figures", "Publication_figures", "Publication_tables", "Study_report"),
                scripts_source_leaf = "metabolomics",
                global_vars = list(
                    de_output_leaf = "de_metabolites",
                    pathway_leaf = "pathway_enrichment",
                    feature_qc_leaf = "metabolite_qc",
                    subfeature_qc_leaf = NULL, # No peptide equivalent for metabolomics typically
                    clean_features_leaf = "clean_metabolites"
                )
            ),
            transcriptomics = list(
                results_subdirs = c(
                    "gene_qc", "count_data", "normalized_counts", "de_genes",
                    "publication_graphs", "pathway_enrichment",
                    file.path("publication_graphs", "filtering_qc")
                ),
                results_summary_subdirs = c("QC_figures", "Publication_figures", "Publication_tables", "Study_report"),
                scripts_source_leaf = "transcriptomics",
                global_vars = list(
                    de_output_leaf = "de_genes",
                    pathway_leaf = "pathway_enrichment",
                    feature_qc_leaf = "gene_qc",
                    subfeature_qc_leaf = NULL,
                    clean_features_leaf = "normalized_counts", # For normalized data, distinct from raw counts
                    raw_counts_leaf = "count_data" # Specific for transcriptomics
                )
            ),
            lipidomics = list( # New configuration for lipidomics
                results_subdirs = c(
                    "lipid_qc", "feature_qc", "clean_lipids", "de_lipids",
                    "publication_graphs", "pathway_enrichment",
                    file.path("publication_graphs", "filtering_qc")
                ),
                results_summary_subdirs = c("QC_figures", "Publication_figures", "Publication_tables", "Study_report"),
                scripts_source_leaf = "lipidomics",
                global_vars = list(
                    de_output_leaf = "de_lipids",
                    pathway_leaf = "pathway_enrichment",
                    feature_qc_leaf = "lipid_qc",
                    subfeature_qc_leaf = NULL,
                    clean_features_leaf = "clean_lipids"
                )
            ),
            integration = list( # New configuration for integration
                results_subdirs = c(
                    "multiomic_inputs",
                    file.path("mofa", "inputs"), file.path("mofa", "plots"), file.path("mofa", "tables"),
                    file.path("mixomics", "inputs"), file.path("mixomics", "plots"), file.path("mixomics", "tables"),
                    file.path("integration_enrichment", "inputs"), file.path("integration_enrichment", "plots"), file.path("integration_enrichment", "tables")
                ),
                results_summary_subdirs = c("QC_figures", "Publication_figures", "Publication_tables", "Study_report"), # Standard for now
                scripts_source_leaf = "integration",
                global_vars = list(
                    multiomic_inputs_leaf = "multiomic_inputs",
                    mofa_inputs_leaf = file.path("mofa", "inputs"),
                    mofa_plots_leaf = file.path("mofa", "plots"),
                    mofa_tables_leaf = file.path("mofa", "tables"),
                    mixomics_inputs_leaf = file.path("mixomics", "inputs"),
                    mixomics_plots_leaf = file.path("mixomics", "plots"),
                    mixomics_tables_leaf = file.path("mixomics", "tables"),
                    integration_enrichment_inputs_leaf = file.path("integration_enrichment", "inputs"),
                    integration_enrichment_plots_leaf = file.path("integration_enrichment", "plots"),
                    integration_enrichment_tables_leaf = file.path("integration_enrichment", "tables")
                    # No generic feature_qc, de_output etc. for integration
                )
            ),
            { # Default case for the switch, though validation should prevent this
                err_msg <- sprintf("Internal error: Configuration not found for validated omic_type: %s", current_omic_type)
                logger::log_error(err_msg)
                stop(err_msg)
            }
        )

        # Construct Directory Names for the current omic type
        omic_label_dirname <- paste0(current_omic_type, if (!is.null(label)) paste0("_", substr(label, 1, 30)) else "")

        results_path <- file.path(base_dir, "results", omic_label_dirname)
        results_summary_path <- file.path(base_dir, "results_summary", omic_label_dirname)
        scripts_path <- file.path(base_dir, "scripts", omic_label_dirname) # Destination scripts path for this omic

        # --- User Interaction for Directory Handling ---
        process_current_omic <- TRUE
        reuse_current_omic_dirs <- FALSE

        # Check Existing Directories and User Confirmation
        existing_dirs_found <- dir.exists(results_path) || dir.exists(results_summary_path) || dir.exists(scripts_path)

        if (reuse_existing && existing_dirs_found) {
            logger::log_info("Reusing existing directory structure for '{omic_label_dirname}' as requested by non-interactive call.")
            reuse_current_omic_dirs <- TRUE
        } else if (!force && existing_dirs_found) {
            logger::log_warn("Key directory(ies) for '{current_omic_type}' (label: '{omic_label_dirname}') already exist:")
            if (dir.exists(results_path)) logger::log_warn("- {results_path}")
            if (dir.exists(results_summary_path)) logger::log_warn("- {results_summary_path}")
            if (dir.exists(scripts_path)) logger::log_warn("- {scripts_path}")

            response_overwrite <- readline(prompt = sprintf("Overwrite directories for '%s'? (y/n): ", omic_label_dirname))

            if (tolower(substr(response_overwrite, 1, 1)) == "y") {
                logger::log_info("Overwriting existing directories for '{omic_label_dirname}' as requested.")
                if (dir.exists(results_path)) unlink(results_path, recursive = TRUE, force = TRUE)
                if (dir.exists(results_summary_path)) unlink(results_summary_path, recursive = TRUE, force = TRUE)
                if (dir.exists(scripts_path)) unlink(scripts_path, recursive = TRUE, force = TRUE)
                # reuse_current_omic_dirs remains FALSE, proceed to create new structure
            } else {
                response_reuse <- readline(prompt = sprintf("Use existing directory structure for '%s' this session? (y/n): ", omic_label_dirname))
                if (tolower(substr(response_reuse, 1, 1)) == "y") {
                    logger::log_info("Reusing existing directory structure for '{omic_label_dirname}'.")
                    reuse_current_omic_dirs <- TRUE
                } else {
                    logger::log_info("Setup for '{omic_label_dirname}' cancelled by user. Directory variables will not be set for this omic type.")
                    process_current_omic <- FALSE
                }
            }
        } else if (force && existing_dirs_found) {
            logger::log_info("Force=TRUE: Overwriting existing directories for '{omic_label_dirname}' if they exist.")
            if (dir.exists(results_path)) unlink(results_path, recursive = TRUE, force = TRUE)
            if (dir.exists(results_summary_path)) unlink(results_summary_path, recursive = TRUE, force = TRUE)
            if (dir.exists(scripts_path)) unlink(scripts_path, recursive = TRUE, force = TRUE)
        }

        if (!process_current_omic) {
            next # Skip to the next omic type in the loop
        }

        # Define Paths for the current omic type
        publication_graphs_dir_base <- file.path(results_path, "publication_graphs")
        # For omics types like 'integration', 'publication_graphs' might not be a primary subdir.
        # We still define qc_dir_base relative to where 'filtering_qc' would be if 'publication_graphs' existed.
        # If 'publication_graphs' is not in omic_config$results_subdirs, it will be created if time_dir needs it.
        qc_dir_base <- file.path(results_path, "publication_graphs", "filtering_qc")

        # Create Omic-Specific Data Directory (e.g., base_dir/data/proteomics)
        # This is different from the previous common_data_dir
        omic_specific_data_dir <- file.path(base_dir, "data", current_omic_type) # Uses current_omic_type, no label
        dir.create(omic_specific_data_dir, recursive = TRUE, showWarnings = FALSE) # Idempotent
        logger::log_info("Ensured omic-specific data directory exists: {omic_specific_data_dir}")

        current_omic_paths_def <- list( # Using a temporary name to avoid confusion with return list
            results_base = results_path,
            # results_subdirs paths will be fully qualified below
            results_summary_base = results_summary_path,
            # results_summary_subdirs paths will be fully qualified below
            data_dir = omic_specific_data_dir, # Corrected to omic-specific data directory
            scripts_dest_dir = scripts_path,
            scripts_source_dir = file.path(base_dir, "scripts", omic_config$scripts_source_leaf),
            publication_graphs_dir = publication_graphs_dir_base, # May not exist for all omics
            qc_dir = qc_dir_base, # Base for timestamped, may not exist for all omics
            time_dir = file.path(qc_dir_base, timestamp) # Timestamped dir, ensures parent qc_dir and pub_graphs exists if used
        )
        
        # Ensure the timestamped directory for THIS session exists if we are processing this omic.
        # This includes its parent directories like qc_dir and publication_graphs_dir if they are part of the path.
        dir.create(current_omic_paths_def$time_dir, recursive = TRUE, showWarnings = FALSE)
        logger::log_info("Ensured timestamped directory exists: {current_omic_paths_def$time_dir}")


        # --- Create Directories and Copy Scripts (if not reusing) ---
        if (!reuse_current_omic_dirs) {
            logger::log_info("Creating directory structure for '{omic_label_dirname}'...")

            # Create results base and defined subdirs
            dir.create(current_omic_paths_def$results_base, recursive = TRUE, showWarnings = FALSE)
            qualified_results_subdirs <- file.path(current_omic_paths_def$results_base, omic_config$results_subdirs)
            invisible(sapply(qualified_results_subdirs, dir.create, recursive = TRUE, showWarnings = FALSE))

            # Create results_summary base and defined subdirs
            dir.create(current_omic_paths_def$results_summary_base, recursive = TRUE, showWarnings = FALSE)
            qualified_results_summary_subdirs <- file.path(current_omic_paths_def$results_summary_base, omic_config$results_summary_subdirs)
            invisible(sapply(qualified_results_summary_subdirs, dir.create, recursive = TRUE, showWarnings = FALSE))
            
            # Note: common_data_dir and time_dir (and its parents) are already created by this point.

            # --- Handle Scripts Directory Copying ---
            if (dir.exists(current_omic_paths_def$scripts_source_dir)) {
                dir.create(current_omic_paths_def$scripts_dest_dir, recursive = TRUE, showWarnings = FALSE) # Ensure destination exists

                script_files <- list.files(current_omic_paths_def$scripts_source_dir, full.names = TRUE, recursive = TRUE)
                script_files <- script_files[!grepl("renv/|renv\\.lock|\\.Rmd$", script_files, ignore.case = TRUE)]

                if (length(script_files) > 0) {
                    logger::log_info("Copying scripts for {current_omic_type} from {current_omic_paths_def$scripts_source_dir} to {current_omic_paths_def$scripts_dest_dir}...")
                    copied_scripts_ok <- TRUE
                    # Normalize source path to forward slashes for Windows regex compatibility
                    source_abs_normalized <- gsub("\\\\", "/", tools::file_path_as_absolute(current_omic_paths_def$scripts_source_dir))
                    for (f in script_files) {
                        # Normalize file path to forward slashes for consistent regex matching on Windows
                        f_abs_normalized <- gsub("\\\\", "/", tools::file_path_as_absolute(f))
                        rel_path <- sub(paste0("^", source_abs_normalized, "/?"), "", f_abs_normalized)
                        dest_file <- file.path(current_omic_paths_def$scripts_dest_dir, rel_path)
                        dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
                        if (!file.copy(f, dest_file, overwrite = TRUE)) {
                            logger::log_warn("Failed to copy script: {f} to {dest_file}")
                            copied_scripts_ok <- FALSE
                        }
                    }
                    if (!copied_scripts_ok) {
                        logger::log_warn("Some script files failed to copy for {current_omic_type}.")
                    }
                } else {
                    logger::log_info("No scripts (excluding Rmd/renv) found in {current_omic_paths_def$scripts_source_dir} to copy.")
                }
            } else {
                dir.create(current_omic_paths_def$scripts_dest_dir, recursive = TRUE, showWarnings = FALSE) # Still create dest dir
                logger::log_warn("Source script directory not found for {current_omic_type}: {current_omic_paths_def$scripts_source_dir}. No scripts copied.")
            }
        } else {
            logger::log_info("Skipping directory creation and script copying for '{omic_label_dirname}' as existing structure is being reused.")
            # Timestamped dir (current_omic_paths_def$time_dir) is already created.
        }


        # Create Integration Directory (Shared across all omics)
        integration_dir <- file.path(base_dir, "integration")
        dir.create(integration_dir, recursive = TRUE, showWarnings = FALSE)

        # --- Build and Store Path Variables for the current omic type ---
        # Start with common/generic names relative to the current omic's paths
        omic_specific_paths_list <- list(
            base_dir = base_dir,
            omic_type = current_omic_type, # Add omic_type for context
            omic_label = omic_label_dirname, # Add omic_label for context
            results_dir = current_omic_paths_def$results_base,
            data_dir = current_omic_paths_def$data_dir,
            source_dir = current_omic_paths_def$scripts_dest_dir,
            timestamp = timestamp,
            integration_dir = integration_dir,
            # qc_dir and publication_graphs_dir might not be primary for all (e.g. integration)
            # but time_dir uses them in its path, so they are defined.
            # Their direct relevance as top-level named paths depends on the omic type.
            publication_graphs_dir = current_omic_paths_def$publication_graphs_dir,
            qc_dir = current_omic_paths_def$qc_dir, # Base of timestamped outputs
            time_dir = current_omic_paths_def$time_dir, # Session-specific output dir
            results_summary_dir = current_omic_paths_def$results_summary_base,
            
            # Standard summary subdirectories, created if results_summary_subdirs is not empty
            qc_figures_dir = file.path(current_omic_paths_def$results_summary_base, "QC_figures"),
            publication_figures_dir = file.path(current_omic_paths_def$results_summary_base, "Publication_figures"),
            publication_tables_dir = file.path(current_omic_paths_def$results_summary_base, "Publication_tables"),
            study_report_dir = file.path(current_omic_paths_def$results_summary_base, "Study_report")
        )

        # Add UniProt Annotation Directory specifically for proteomics
        if (current_omic_type == "proteomics") {
            uniprot_annotation_path <- file.path(base_dir, "data", "UniProt")
            dir.create(uniprot_annotation_path, recursive = TRUE, showWarnings = FALSE)
            logger::log_info("Ensured UniProt annotation directory exists: {uniprot_annotation_path}")
            omic_specific_paths_list$uniprot_annotation_dir <- uniprot_annotation_path
        }

        # Add omics-category specific paths from global_vars config
        # These are more specific than the general subdirs and provide convenient named variables
        if (!is.null(omic_config$global_vars)) {
            for (var_name_suffix in names(omic_config$global_vars)) {
                leaf_path <- omic_config$global_vars[[var_name_suffix]]
                if (!is.null(leaf_path)) {
                    # Construct the full path relative to the results_dir for this omic type
                    full_path <- file.path(current_omic_paths_def$results_base, leaf_path)
                    # The actual variable name will be like 'de_output_dir', 'feature_qc_dir'
                    # If var_name_suffix is "de_output_leaf", actual name is "de_output_dir"
                    actual_var_name <- sub("_leaf$", "_dir", var_name_suffix)
                    omic_specific_paths_list[[actual_var_name]] <- full_path
                }
            }
        }
        
        # Add specific variable names that might map to the more generic ones, for convenience or legacy.
        # This part also ensures specific directories (like transcriptomics count_data_dir) are explicitly listed if defined in results_subdirs.
        specific_mapped_names <- list()
        if (current_omic_type == "proteomics") {
            if(!is.null(omic_specific_paths_list$feature_qc_dir)) specific_mapped_names$protein_qc_dir = omic_specific_paths_list$feature_qc_dir
            if(!is.null(omic_specific_paths_list$subfeature_qc_dir)) specific_mapped_names$peptide_qc_dir = omic_specific_paths_list$subfeature_qc_dir
            if(!is.null(omic_specific_paths_list$clean_features_dir)) specific_mapped_names$clean_proteins_dir = omic_specific_paths_list$clean_features_dir
        } else if (current_omic_type == "metabolomics") {
            if(!is.null(omic_specific_paths_list$feature_qc_dir)) specific_mapped_names$metabolite_qc_dir = omic_specific_paths_list$feature_qc_dir
            if(!is.null(omic_specific_paths_list$clean_features_dir)) specific_mapped_names$clean_metabolites_dir = omic_specific_paths_list$clean_features_dir
        } else if (current_omic_type == "transcriptomics") {
            if(!is.null(omic_specific_paths_list$feature_qc_dir)) specific_mapped_names$gene_qc_dir = omic_specific_paths_list$feature_qc_dir
            # raw_counts_dir is already added from global_vars if raw_counts_leaf is defined
            if(!is.null(omic_specific_paths_list$clean_features_dir)) specific_mapped_names$normalized_counts_dir = omic_specific_paths_list$clean_features_dir
        } else if (current_omic_type == "lipidomics") {
            if(!is.null(omic_specific_paths_list$feature_qc_dir)) specific_mapped_names$lipid_qc_dir = omic_specific_paths_list$feature_qc_dir
            if(!is.null(omic_specific_paths_list$clean_features_dir)) specific_mapped_names$clean_lipids_dir = omic_specific_paths_list$clean_features_dir
        }
        # For 'integration', specific paths like mofa_inputs_dir, etc., are already directly added from global_vars.

        omic_specific_paths_list <- c(omic_specific_paths_list, specific_mapped_names[!sapply(specific_mapped_names, is.null)])

        # Ensure all defined directories exist (belt-and-suspenders, esp. for reused structures or complex paths)
        # This primarily targets directories that are values in omic_specific_paths_list
        for(path_entry in omic_specific_paths_list) {
            if(is.character(path_entry) && (endsWith(path_entry, "_dir") || endsWith(path_entry, "_base") || grepl("mofa|mixomics|integration_enrichment", path_entry))) {
                 # Avoid trying to create dir for base_dir, timestamp, omic_type, omic_label etc.
                 # A bit heuristic here; might need refinement if path names are ambiguous.
                 # Check if it's likely a directory path constructed within the project.
                if (startsWith(path_entry, base_dir) && path_entry != base_dir && !is.null(path_entry) && path_entry != "") {
                    suppressWarnings(dir.create(path_entry, recursive = TRUE, showWarnings = FALSE))
                }
            }
        }

        all_created_paths[[omic_label_dirname]] <- omic_specific_paths_list # Use omic_label_dirname as key
        logger::log_info("Stored paths for {omic_label_dirname}.")
    } # End loop over omic types

    # --- Print Structure ---
    logger::log_info("--- Directory Structure Setup Complete ---")
    if (length(all_created_paths) > 0) {
        for (omic_key_name in names(all_created_paths)) {
            cat(sprintf("\nPaths for Omics Label: '%s' (Type: %s)\n", omic_key_name, all_created_paths[[omic_key_name]]$omic_type))
            current_paths_to_print <- all_created_paths[[omic_key_name]]
            
            # Remove non-path-like informational fields for printing
            current_paths_to_print$omic_type <- NULL
            current_paths_to_print$omic_label <- NULL

            print_order <- c("base_dir", "data_dir", "source_dir", "results_dir", "results_summary_dir", "timestamp", "time_dir")
            remaining_names <- setdiff(names(current_paths_to_print), print_order)
            sorted_remaining_names <- sort(remaining_names)
            final_print_order <- c(print_order, sorted_remaining_names)
            
            # Filter out any names that might not be in current_paths_to_print (e.g. if an optional path was NULL)
            final_print_order <- final_print_order[final_print_order %in% names(current_paths_to_print)]


            for (path_name in final_print_order) {
                p_val <- current_paths_to_print[[path_name]]
                if (is.null(p_val)) next # Skip if path was not defined (e.g. optional subfeature for some omics)

                if (is.character(p_val) && p_val != "") {
                    # Heuristic to decide if it's a directory to count items in
                    is_dir_to_count <-endsWith(path_name, "_dir") || endsWith(path_name, "_base") || path_name %in% c("qc_dir", "time_dir")
                    
                    if (is_dir_to_count && dir.exists(p_val)) {
                        tryCatch(
                            {
                                files_in_dir <- list.files(p_val, recursive = FALSE, all.files = TRUE, no.. = TRUE) # Non-recursive for cleaner count
                                item_count <- length(files_in_dir)
                                cat(sprintf("  %-30s : %s (%d items)\n", path_name, p_val, item_count))
                            },
                            error = function(e) {
                                cat(sprintf("  %-30s : %s (Error accessing content)\n", path_name, p_val))
                            }
                        )
                    } else if (is_dir_to_count && !dir.exists(p_val)) {
                         cat(sprintf("  %-30s : %s (Directory not found/created!)\n", path_name, p_val))
                    } else {
                        # For non-directory strings like timestamp
                        cat(sprintf("  %-30s : %s\n", path_name, p_val))
                    }
                }
            }
        }
    } else {
        logger::log_info("No directories were set up for any omic type (e.g., setup cancelled for all).")
    }
    cat("------------------------------------------\n")


    # Return the collected paths invisibly
    invisible(all_created_paths)
}


# ----------------------------------------------------------------------------
# saveListOfPdfs
# ----------------------------------------------------------------------------
#' @export
saveListOfPdfs <- function(list, filename) {
  #start pdf
  cairo_pdf(filename)

  #loop
  #purrr::walk( list, print)
  for (p in list) {
    print(p)
  }

  #end pdf
  dev.off()

  invisible(NULL)
}


# ----------------------------------------------------------------------------
# getProjectPaths
# ----------------------------------------------------------------------------
#' Get Project Paths with Fallback Logic
#' 
#' @description
#' Safely retrieves project paths from project_dirs with automatic fallback between
#' different key formats. Handles both Shiny (keys: "proteomics") and RMarkdown 
#' (keys: "proteomics_MyLabel") workflows.
#' 
#' @param omic_type Character string indicating the type of omics data (e.g., "proteomics")
#' @param experiment_label Character string with experiment label (optional for Shiny, required for RMarkdown)
#' @param project_dirs_object_name Name of the global project_dirs variable (default: "project_dirs")
#' @param env Environment where project_dirs resides (default: .GlobalEnv)
#' 
#' @return List containing the directory paths for the specified omic type
#' 
#' @export
getProjectPaths <- function(omic_type, 
                           experiment_label = NULL, 
                           project_dirs_object_name = "project_dirs",
                           env = .GlobalEnv) {
  
  message("--- DEBUG66: Entering getProjectPaths ---")
  message(sprintf("   DEBUG66: omic_type = '%s'", omic_type))
  message(sprintf("   DEBUG66: experiment_label = '%s'", ifelse(is.null(experiment_label), "NULL", experiment_label)))
  message(sprintf("   DEBUG66: project_dirs_object_name = '%s'", project_dirs_object_name))
  
  # Validate project_dirs exists
  if (!exists(project_dirs_object_name, envir = env)) {
    rlang::abort(paste0("Global object ", sQuote(project_dirs_object_name), 
                       " not found. Run setupDirectories() first or ensure project is properly initialized."))
  }
  
  project_dirs_global <- get(project_dirs_object_name, envir = env)
  
  message("   DEBUG66: project_dirs object found in global environment")
  message(sprintf("      DEBUG66: project_dirs has %d keys", length(names(project_dirs_global))))
  message(sprintf("      DEBUG66: Available keys: %s", paste(names(project_dirs_global), collapse=", ")))
  
  # Try multiple key formats in order of preference
  potential_keys <- c()
  
  # 1. RMarkdown format with label: "omic_type_experiment_label"
  if (!is.null(experiment_label) && nzchar(experiment_label)) {
    potential_keys <- c(potential_keys, paste0(omic_type, "_", experiment_label))
  }
  
  # 2. Shiny format (no label): "omic_type"
  potential_keys <- c(potential_keys, omic_type)
  
  # 3. Try to find any key that starts with omic_type (most permissive fallback)
  matching_keys <- names(project_dirs_global)[grepl(paste0("^", omic_type, "(_|$)"), names(project_dirs_global))]
  if (length(matching_keys) > 0) {
    potential_keys <- c(potential_keys, matching_keys[1])  # Take first match
  }
  
  message(sprintf("   DEBUG66: Trying %d potential keys: %s", length(potential_keys), paste(potential_keys, collapse=", ")))
  
  # Try each potential key
  for (key in potential_keys) {
    message(sprintf("   DEBUG66: Trying key: '%s'", key))
    if (key %in% names(project_dirs_global)) {
      message(sprintf("   DEBUG66: SUCCESS - Found match with key: '%s'", key))
      logger::log_debug("Found project paths using key: {key}")
      return(project_dirs_global[[key]])
    } else {
      message(sprintf("      DEBUG66: Key '%s' not found, continuing...", key))
    }
  }
  
  # If we get here, no valid key was found
  message("   DEBUG66: No valid key found - ABORTING")
  available_keys <- paste(names(project_dirs_global), collapse = ", ")
  rlang::abort(paste0(
    "Could not find project paths for omic_type='", omic_type, "'",
    if (!is.null(experiment_label)) paste0(", experiment_label='", experiment_label, "'") else "",
    ".\nAvailable keys in ", sQuote(project_dirs_object_name), ": ", available_keys,
    "\nEnsure setupDirectories() was called with the correct omic_type and label."
  ))
}


# ----------------------------------------------------------------------------
# createDirectoryIfNotExists
# ----------------------------------------------------------------------------
#' @export
createDirectoryIfNotExists <- function(file_path, mode = "0777") {

  ### Create directory recursively if it doesn't exist
  if (! file.exists(file_path)){

    dir.create(file_path, showWarnings = TRUE, recursive = TRUE, mode = mode)

  }

}


# ----------------------------------------------------------------------------
# createDirIfNotExists
# ----------------------------------------------------------------------------
#' @export
createDirIfNotExists  <- function(file_path, mode = "0777") {
  createDirectoryIfNotExists(file_path, mode = "0777")
}


# ----------------------------------------------------------------------------
# sourceRmdFileSimple
# ----------------------------------------------------------------------------
## Function to source Rmd files
# https://stackoverflow.com/questions/10966109/how-to-source-r-markdown-file-like-sourcemyfile-r
#' @export
sourceRmdFileSimple <- function(x, ...) {
  source(purl(x, output = tempfile()), ...)
}


# ----------------------------------------------------------------------------
# sourceRmdFile
# ----------------------------------------------------------------------------
#' https://gist.github.com/noamross/a549ee50e8a4fd68b8b1
#' Source the R code from an knitr file, optionally skipping plots
#'
#' @param file the knitr file to source
#' @param skip_plots whether to make plots. If TRUE (default) sets a null graphics device
#'
#' @return This function is called for its side effects
#' @export
sourceRmdFile <- function(file, skip_plots = TRUE) {
  temp = tempfile(fileext=".R")
  knitr::purl(file, output=temp)

  if(skip_plots) {
    old_dev = getOption('device')
    options(device = function(...) {
      .Call("R_GD_nullDevice", PACKAGE = "grDevices")
    })
  }
  source(temp)
  if(skip_plots) {
    options(device = old_dev)
  }
}


# ----------------------------------------------------------------------------
# createOutputDir
# ----------------------------------------------------------------------------
##################################################################################################################
#=====================================================================================================
#' @export
createOutputDir <- function(output_dir, no_backup) {
  if (output_dir == "") {
    logerror("output_dir is an empty string")
    q()
  }
  if (dir.exists(output_dir)) {
    if (no_backup) {
      unlink(output_dir, recursive = TRUE)
    }
    else {
      backup_name <- paste(output_dir, "_prev", sep = "")
      if (dir.exists(backup_name)) { unlink(backup_name, recursive = TRUE) }
      system(paste("mv", output_dir, backup_name)) }
  }
  dir.create(output_dir, recursive = TRUE)
}


# ----------------------------------------------------------------------------
# testRequiredFiles
# ----------------------------------------------------------------------------
#' @export
testRequiredFiles <- function(files) {
  missing_files <- !file.exists(files)
  invisible(sapply(files[missing_files], function(file) {
    logerror("Missing required file: %s", file)
    q()
  }))
}


# ----------------------------------------------------------------------------
# testRequiredFilesWarning
# ----------------------------------------------------------------------------
#' @export
testRequiredFilesWarning <- function(files) {
  missing_files <- !file.exists(files)
  invisible(sapply(files[missing_files], function(file) {
    logwarn("Missing required file: %s", file)
  }))
}


# ----------------------------------------------------------------------------
# testRequiredArguments
# ----------------------------------------------------------------------------
#' @export
testRequiredArguments <- function(arg_list, parameters) {
  invisible(sapply(parameters, function(par) {
    if (!par %in% names(arg_list)) {
      logerror("Missing required argument: %s", par)
      q()
    }
  }))
}


# ----------------------------------------------------------------------------
# parseType
# ----------------------------------------------------------------------------
#' @export
parseType<-function (arg_list,parameters,functType){
  invisible(sapply(parameters, function(key) {
    arg_list[key]<-functType(arg_list[key])
  }))
  return (arg_list)
}


# ----------------------------------------------------------------------------
# parseString
# ----------------------------------------------------------------------------
#' @export
parseString<-function (arg_list,parameters){
  invisible(sapply(parameters, function(key) {
    if(key %in% names(arg_list)) {
      arg_list[key] <- str_replace_all(arg_list[key], c("\"" = "", "\'" = ""))
    }
  }))
  return (arg_list)
}


# ----------------------------------------------------------------------------
# parseList
# ----------------------------------------------------------------------------
#' @export
parseList<-function (arg_list,parameters){
  invisible(sapply(parameters, function(key) {
    items<-str_replace_all(as.character(arg_list[key])," ","")
    arg_list[key]<- base::strsplit(items,split=",")
  }))
  return (arg_list)
}


# ----------------------------------------------------------------------------
# isArgumentDefined
# ----------------------------------------------------------------------------
#' @export
isArgumentDefined<-function(arg_list,parameter){
  return (!is.null(arg_list[parameter]) & (parameter %in% names(arg_list)) & as.character(arg_list[parameter]) != "")
}


# ----------------------------------------------------------------------------
# savePlot
# ----------------------------------------------------------------------------
#' Save a plot in multiple formats
#'
#' This function saves a given plot in multiple specified formats and also save the ggplot object as a rds file.
#'
#' @param plot The plot object to be saved
#' @param base_path The base directory path where the plot will be saved
#' @param plot_name The name to be used for the saved plot files
#' @param formats A vector of file formats to save the plot in (default: c("pdf", "png"))
#' @param width The width of the plot (default: 7)
#' @param height The height of the plot (default: 7)
#' @param ... Additional arguments to be passed to ggsave
#'
#' @return This function is called for its side effects (saving files)
#' @export
#'
#'
savePlot <- function(plot, base_path, plot_name, formats = c("pdf", "png"), width=7, height=7, ... ) {
  # Always save the RDS (works for both single plots and lists)
  saveRDS( plot, file.path(base_path, paste0(plot_name, ".rds")))
  
  # Check if plot is a list of plots
  if (is.list(plot) && !inherits(plot, "gg")) {
    # It's a list of plots - save each one individually
    plot_names <- names(plot)
    if (is.null(plot_names)) {
      plot_names <- paste0("plot_", seq_along(plot))
    }
    
    purrr::walk2(plot, plot_names, function(p, pname) {
      if (inherits(p, "gg")) {
        purrr::walk(formats, function(format) {
          file_path <- file.path(base_path, paste0(plot_name, "_", pname, ".", format))
          ggsave(filename = file_path, plot = p, device = format, width=width, height=height, ...)
        })
      }
    })
  } else {
    # Single plot - original behavior
    purrr::walk( formats, \(format){
      file_path <- file.path(base_path, paste0(plot_name, ".", format))
      ggsave(filename = file_path, plot = plot, device = format, width=width, height=height, ...)
    })
  }
}


# ----------------------------------------------------------------------------
# save_plot
# ----------------------------------------------------------------------------
#' Save a plot in multiple formats
#'
#' This function saves a given plot in multiple specified formats and also save the ggplot object as a rds file.
#'
#' @param plot The plot object to be saved
#' @param base_path The base directory path where the plot will be saved
#' @param plot_name The name to be used for the saved plot files
#' @param formats A vector of file formats to save the plot in (default: c("pdf", "png"))
##' @param width The width of the plot (default: 7)
#' @param height The height of the plot (default: 7)
#' @param ... Additional arguments to be passed to ggsave
#'
#' @return This function is called for its side effects (saving files)
#' @export
#'
#'
save_plot <- function(plot, base_path, plot_name, formats = c("pdf", "png"), width=7, height=7, ... ) {
  savePlot(plot, base_path, plot_name, formats, width, height, ...)
}


# ----------------------------------------------------------------------------
# write_results
# ----------------------------------------------------------------------------
write_results <- function(data, filename) {
  vroom::vroom_write(data, file.path(results_dir, "protein_qc", filename))
}


# ----------------------------------------------------------------------------
# readConfigFile
# ----------------------------------------------------------------------------
#' @title Read the config file and return the list of parameters
#' @description Read the config file and return the list of parameters
#' @param file The file path to the config file
#' @param file_type The type of the file (default: "ini")
#' @export
readConfigFile <- function( file=file.path(source_dir, "config.ini")) {

  config_list <- read.config(file=file, file.type = "ini" )

  # to set the number of cores to be used in the parallel processing
  if("globalParameters" %in% names(config_list)) {
    if ( "number_of_cpus" %in% names( config_list[["globalParameters"]])  ) {

      print(paste0("Read globalParameters: number_of_cpus = "
                   , config_list$globalParameters$number_of_cpus))
      core_utilisation <- new_cluster(config_list$globalParameters$number_of_cpus)
      cluster_library(core_utilisation, c("tidyverse", "glue", "rlang", "lazyeval"))

      list_of_multithreaded_functions <- c("rollUpPrecursorToPeptide"
                                           , "peptideIntensityFiltering"
                                           , "filterMinNumPeptidesPerProtein"
                                           , "filterMinNumPeptidesPerSample"
                                           , "removePeptidesWithOnlyOneReplicate"
                                           , "peptideMissingValueImputation"
                                           , "removeProteinsWithOnlyOneReplicate")

      setCoreUtilisation <- function(config_list, function_name) {
        if (!function_name %in% names(config_list)) {
          config_list[[function_name]] <- list()
        }
        config_list[[function_name]][["core_utilisation"]] <- core_utilisation

        config_list
      }

      for( x in list_of_multithreaded_functions) {
        config_list <- setCoreUtilisation(config_list, x)
      }

      config_list[["globalParameters"]][["plots_format"]] <- str_split(config_list[["globalParameters"]][["plots_format"]], ",")[[1]]

    }}

  getConfigValue <- function (config_list, section, value) {
    config_list[[section]][[value]]
  }

  setConfigValueAsNumeric <- function (config_list, section, value) {
    config_list[[section]][[value]] <- as.numeric(config_list[[section]][[value]])
    config_list
  }

  if("srlQvalueProteotypicPeptideClean" %in% names(config_list)) {
    # Parse input_matrix_column_ids - handle both formats:
    # 1. Already a vector (from fresh workflow)
    # 2. String like 'c("Run", "Precursor.Id", ...)' (from saved config.ini)
    raw_value <- config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]]
    
    if (is.character(raw_value) && length(raw_value) == 1 && grepl("^c\\(", raw_value)) {
      # Format: 'c("Run", "Precursor.Id", ...)' - parse as R code
      config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]] <- 
        eval(parse(text = raw_value))
    } else if (is.character(raw_value) && length(raw_value) == 1) {
      # Format: 'Run, Precursor.Id, ...' - split by comma
      raw_split <- str_split(raw_value, ",")[[1]]
      config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]] <- trimws(raw_split)
    }
    # else: already a vector, keep as is
    
    # Remove any empty strings
    config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]] <- 
      config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]][
        config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]] != ""
      ]

    print(paste0("Read srlQvalueProteotypicPeptideClean: input_matrix_column_ids = "
                 , paste0(config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]]
                          , collapse=", ")))

    config_list <- setConfigValueAsNumeric(config_list
                                           , "srlQvalueProteotypicPeptideClean"
                                           , "qvalue_threshold")
    config_list <- setConfigValueAsNumeric(config_list
                                           , "srlQvalueProteotypicPeptideClean"
                                           , "global_qvalue_threshold")
    config_list <- setConfigValueAsNumeric(config_list
                                           , "srlQvalueProteotypicPeptideClean"
                                           , "choose_only_proteotypic_peptide")

  }


  if("peptideIntensityFiltering" %in% names(config_list)) {
    config_list <- setConfigValueAsNumeric(config_list
                                           , "peptideIntensityFiltering"
                                           , "peptides_intensity_cutoff_percentile")
    config_list <- setConfigValueAsNumeric(config_list
                                           , "peptideIntensityFiltering"
                                           , "peptides_proportion_of_samples_below_cutoff")
  }


  if("filterMinNumPeptidesPerProtein" %in% names(config_list)) {
    config_list <- setConfigValueAsNumeric(config_list
                                           , "filterMinNumPeptidesPerProtein"
                                           , "peptides_per_protein_cutoff")
    config_list <- setConfigValueAsNumeric(config_list
                                           , "filterMinNumPeptidesPerProtein"
                                           , "peptidoforms_per_protein_cutoff")
    # config_list <- setConfigValueAsNumeric(config_list
    #                                        , ""
    #                                        , "")
  }

  if("filterMinNumPeptidesPerSample" %in% names(config_list)) {
    config_list <- setConfigValueAsNumeric(config_list
                                           , "filterMinNumPeptidesPerSample"
                                           , "peptides_per_sample_cutoff")

    if(!"inclusion_list" %in% names( config_list[["filterMinNumPeptidesPerSample"]])) {
      config_list[["filterMinNumPeptidesPerSample"]][["inclusion_list"]] <- ""
    }

    config_list[["filterMinNumPeptidesPerSample"]][["inclusion_list"]] <- str_split(config_list[["filterMinNumPeptidesPerSample"]][["inclusion_list"]], ",")[[1]]

  }

  if("peptideMissingValueImputation" %in% names(config_list)) {
    config_list <- setConfigValueAsNumeric(config_list
                                           , "peptideMissingValueImputation"
                                           , "proportion_missing_values")
  }

  if("removeRowsWithMissingValuesPercent" %in% names(config_list)) {
    config_list <- setConfigValueAsNumeric(config_list
                                           , "removeRowsWithMissingValuesPercent"
                                           , "groupwise_percentage_cutoff")

    config_list <- setConfigValueAsNumeric(config_list
                                           , "removeRowsWithMissingValuesPercent"
                                           , "max_groups_percentage_cutoff")

    config_list <- setConfigValueAsNumeric(config_list
                                           , "removeRowsWithMissingValuesPercent"
                                           , "proteins_intensity_cutoff_percentile")

  }


  if("ruvIII_C_Varying" %in% names(config_list)) {
    config_list <- setConfigValueAsNumeric(config_list
                                           , "ruvIII_C_Varying"
                                           , "ruv_number_k")
  }

  if("plotRle" %in% names(config_list)) {
    config_list[["plotRle"]][["yaxis_limit"]] <- str_split(config_list[["plotRle"]][["yaxis_limit"]], ",")[[1]] |>
      purrr::map_dbl( \(x) as.numeric(x) )

    print(paste0("Read plotRle: yaxis_limit = "
                 , paste0(config_list[["plotRle"]][["yaxis_limit"]], collapse=", ")))
  }

   if("deAnalysisParameters" %in% names(config_list)) {
    # Handle plots_format as array
    config_list[["deAnalysisParameters"]][["plots_format"]] <-
      str_split(config_list[["deAnalysisParameters"]][["plots_format"]], ",")[[1]]

    # Add new lfc_cutoff parameter
    config_list[["deAnalysisParameters"]][["lfc_cutoff"]] <- FALSE

    # Modify treat_lfc_cutoff to use ifelse
    config_list[["deAnalysisParameters"]][["treat_lfc_cutoff"]] <-
      ifelse(config_list[["deAnalysisParameters"]][["lfc_cutoff"]], log2(1.5), 0)

    # Handle args_group_pattern - remove quotes and fix escaping
    if("args_group_pattern" %in% names(config_list[["deAnalysisParameters"]])) {
      config_list[["deAnalysisParameters"]][["args_group_pattern"]] <-
        gsub('^"|"$', '', config_list[["deAnalysisParameters"]][["args_group_pattern"]]) |>
        gsub(pattern = "\\\\", replacement = "\\")
    }

    # Convert numeric parameters
    config_list <- setConfigValueAsNumeric(config_list,
                                         "deAnalysisParameters",
                                         "de_q_val_thresh")

    # Convert boolean parameters
    config_list[["deAnalysisParameters"]][["eBayes_trend"]] <-
      tolower(config_list[["deAnalysisParameters"]][["eBayes_trend"]]) == "true"
    config_list[["deAnalysisParameters"]][["eBayes_robust"]] <-
      tolower(config_list[["deAnalysisParameters"]][["eBayes_robust"]]) == "true"

    print(paste0("Read deAnalysisParameters: formula_string = ",
                config_list[["deAnalysisParameters"]][["formula_string"]]))
}

  config_list
}






#' @title Read the config file and specify the section and or parameter to update the object
#' @description Read the config file and specify the section and or parameter to update the object
#' @param theObject The object to be updated
#' @param file The file path to the config file
#' @param section The section to be updated
#' @param value The parameter value to be updated
#' @export
readConfigFileSection <- function( theObject
                            , file=file.path(source_dir, "config.ini")
                            , function_name
                            , parameter_name = NULL ) {

  config_list <- readConfigFile( file=file
                              , file_type = "ini" )

  if ( is.null(parameter_name) ) {
    theObject@args[[function_name]] <- config_list[[function_name]]
  } else {
    theObject@args[[function_name]][[parameter_name]] <- config_list[[function_name]][[parameter_name]]
  }

  theObject
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
            tryCatch({
                installer_func(pkg)
                # After install, load it
                pacman::p_load(char = pkg, character.only = TRUE)
            }, error = function(e) {
                warning(sprintf("Failed to install %s from %s: %s", pkg, source_name, e$message))
            })
        } else {
            if (verbose) message(sprintf("%s is already installed, loading...", pkg))
            pacman::p_load(char = pkg, character.only = TRUE)
        }
    }

    # CRAN Packages
    if (verbose) message("\n--- Processing CRAN Packages ---")
    purrr::walk(cran_packages, ~install_and_load(
        pkg = .x,
        # Use base R install.packages directly
        installer_func = function(p) utils::install.packages(p, dependencies = TRUE),
        source_name = "CRAN",
        verbose = verbose
    ))

    # Bioconductor Packages
    if (verbose) message("\n--- Processing Bioconductor Packages ---")
    purrr::walk(bioc_packages, ~install_and_load(
        pkg = .x,
        installer_func = function(p) BiocManager::install(p, update = FALSE, ask = FALSE), # Use BiocManager
        source_name = "Bioconductor",
        verbose = verbose
    ))

    # GitHub Packages
    if (verbose) message("\n--- Processing GitHub Packages ---")
    purrr::iwalk(github_packages, ~{
        pkg_name <- .y # Name of the package (e.g., "RUVIIIC")
        repo <- .x     # Repository path (e.g., "cran/RUVIIIC")
        if (!requireNamespace(pkg_name, quietly = TRUE)) {
             if (verbose) message(sprintf("Installing %s from GitHub (%s)...", pkg_name, repo))
             tryCatch({
                 # Force installation to handle potentially corrupt states
                 devtools::install_github(repo, force = TRUE)
                 pacman::p_load(char = pkg_name, character.only = TRUE)
             }, error = function(e) {
                 warning(sprintf("Failed to install %s from GitHub (%s): %s", pkg_name, repo, e$message))
             })
         } else {
             if (verbose) message(sprintf("%s is already installed, loading...", pkg_name))
             pacman::p_load(char = pkg_name, character.only = TRUE)
         }
    })

    if (verbose) message("\nAll dependencies processed successfully!")
}

##################################################################################################################
#' @title Extract Substrings from Underscore-Separated Strings
#'
#' @description
#' Extracts substrings from underscore-separated strings using different modes:
#' range (between positions), start (first element), or end (last element).
#'
#' @param x Character vector containing the strings to process
#' @param mode Character string specifying extraction mode:
#'   * "range": Extract elements between two underscore positions
#'   * "start": Extract from start to first underscore
#'   * "end": Extract from last underscore to end
#' @param start Integer: Starting position for range mode (1-based)
#' @param end Integer: Ending position for range mode (1-based, required for range mode)
#'
#' @return Character vector with extracted strings. Returns NA for strings where
#'   requested positions are out of bounds.
#'
#' @examples
#' x <- "20140602_ffs_expt1_r1_junk"
#' extract_experiment(x, mode = "range", start = 1, end = 3)  # "20140602_ffs_expt1"
#' extract_experiment(x, mode = "start")  # "20140602"
#' extract_experiment(x, mode = "end")    # "junk"
#'
#' # Multiple strings
#' x <- c("20140602_ffs_expt1_r1_junk", "20140603_ffs_expt2_r2_test")
#' extract_experiment(x, mode = "range", start = 2, end = 3)  # c("ffs_expt1", "ffs_expt2")
#'
#' @export
extract_experiment <- function(x, mode = "range", start = 1, end = NULL) {
  if (!mode %in% c("range", "start", "end")) {
    stop("Mode must be one of: 'range', 'start', 'end'")
  }

  process_string <- function(str) {
    parts <- unlist(strsplit(str, "_"))

    if (mode == "range") {
      if (is.null(end)) stop("End position required for range mode")
      if (start > length(parts) || end > length(parts)) {
        warning("Position out of bounds for string: ", str)

        return(NA_character_)
      }
      return(paste(parts[start:end], collapse = "_"))
    }

    else if (mode == "start") {
      return(parts[1])
    }

    else if (mode == "end") {
      return(parts[length(parts)])
    }
  }

  sapply(x, process_string)
}

#' @title Setup Project Directories
#' @description Creates and manages project directories with version control. If directories exist, prompts user to overwrite or reuse.
#' @param base_dir Base directory path (optional, defaults to here::here())
#' @param label Optional label to append to proteomics directory name (e.g., "proteomics_MyLabel")
#' @param force Logical; if TRUE, skips user confirmation and overwrites existing directories (default: FALSE)
#' @return List of directory paths assigned to the global environment.
#' @export
setupAndShowDirectories <- function(base_dir = here::here(), label = NULL, force = FALSE) {
    # --- Define Base Paths and Names ---
    proteomics_dirname <- if (!is.null(label)) paste0("proteomics_", substr(label, 1, 30)) else "proteomics"
    
    # Define all expected paths in one structure for easier management
    paths <- list(
        results = list(
            base = file.path(base_dir, "results", proteomics_dirname),
            subdirs = c("protein_qc", "peptide_qc", "clean_proteins", "de_proteins", 
                       "publication_graphs", "pathway_enrichment",
                       file.path("publication_graphs", "filtering_qc"))
        ),
        results_summary = list(
            base = file.path(base_dir, "results_summary", proteomics_dirname),
            subdirs = c("QC_figures", "Publication_figures", "Publication_tables", "Study_report")
        ),
        special = list(
            data = file.path(base_dir, "data"),
            scripts = file.path(base_dir, "scripts", proteomics_dirname), # Project-specific scripts
            qc_base = file.path(base_dir, "results", proteomics_dirname, "publication_graphs", "filtering_qc")
            # 'time' directory path defined later using timestamp
        )
    )
    
    results_path <- paths$results$base
    results_summary_path <- paths$results_summary$base
    scripts_path <- paths$special$scripts
    
    # Flag to track if we should skip creation/copying and just set vars
    reuse_existing <- FALSE 

    # --- Handle Existing Directories ---
    if (!force && (dir.exists(results_path) || dir.exists(results_summary_path) || dir.exists(scripts_path))) {
        cat(sprintf("\nWarning: Key directory(ies) already exist:\n"))
        if (dir.exists(results_path)) cat(sprintf("- %s\n", results_path))
        if (dir.exists(results_summary_path)) cat(sprintf("- %s\n", results_summary_path))
        if (dir.exists(scripts_path)) cat(sprintf("- %s\n", scripts_path))
        
        response_overwrite <- readline(prompt = "Do you want to overwrite these directories? (y/n): ")
        
        if (tolower(substr(response_overwrite, 1, 1)) == "y") {
            # User chose to overwrite: Delete existing directories
            message("Overwriting existing directories as requested...")
            if (dir.exists(results_path)) unlink(results_path, recursive = TRUE, force = TRUE)
            if (dir.exists(results_summary_path)) unlink(results_summary_path, recursive = TRUE, force = TRUE)
            if (dir.exists(scripts_path)) unlink(scripts_path, recursive = TRUE, force = TRUE)
            # reuse_existing remains FALSE, proceed to create new structure
        } else {
            # User chose NOT to overwrite
            response_reuse <- readline(prompt = "Do you wish to use the existing directory structure for this session? (y/n): ")
            if (tolower(substr(response_reuse, 1, 1)) == "y") {
                message("Reusing existing directory structure and setting environment variables.")
                reuse_existing <- TRUE # Set flag to skip creation/copying
            } else {
                message("Setup cancelled by user. Directory variables not set.")
                return(invisible(NULL)) # Abort if user neither overwrites nor reuses
            }
        }
    } else if (force) {
        message("Force=TRUE: Overwriting existing directories if they exist.")
        if (dir.exists(results_path)) unlink(results_path, recursive = TRUE, force = TRUE)
        if (dir.exists(results_summary_path)) unlink(results_summary_path, recursive = TRUE, force = TRUE)
        if (dir.exists(scripts_path)) unlink(scripts_path, recursive = TRUE, force = TRUE)
    }
    
    # --- Timestamp and Final Path Definitions ---
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    # Define timestamp-specific path (needed regardless of reuse)
    paths$special$time <- file.path(paths$special$qc_base, timestamp)

    # --- Directory Creation and Script Copying ---
    # Only run if NOT reusing existing structure
    if (!reuse_existing) {
        message("Creating directory structure...")
        # Create results and results_summary base and subdirs
    lapply(c("results", "results_summary"), function(type) {
        dir.create(paths[[type]]$base, recursive = TRUE, showWarnings = FALSE)
        invisible(sapply(
            file.path(paths[[type]]$base, paths[[type]]$subdirs),
            dir.create, recursive = TRUE, showWarnings = FALSE
        ))
    })
    
        # Create special directories (ensure qc_base exists for time dir)
        dir.create(paths$special$qc_base, recursive = TRUE, showWarnings = FALSE)
        dir.create(paths$special$data, recursive = TRUE, showWarnings = FALSE)
        dir.create(paths$special$scripts, recursive = TRUE, showWarnings = FALSE)
        dir.create(paths$special$time, recursive = TRUE, showWarnings = FALSE) # Timestamped dir

        # Handle scripts directory copying from template
        scripts_template_source <- file.path(base_dir, "scripts", "proteomics") # Base template location
        if (dir.exists(scripts_template_source)) {
            message("Copying script files (excluding .Rmd) from template...")
            script_files <- list.files(scripts_template_source, full.names = TRUE, recursive = TRUE)
            script_files <- script_files[!grepl("\\.[rR][mM][dD]$", script_files)] # Filter out .Rmd
        
        if (length(script_files) > 0) {
            # Normalize source path to forward slashes for Windows regex compatibility
            source_abs_normalized <- gsub("\\\\", "/", tools::file_path_as_absolute(scripts_template_source))
            sapply(script_files, function(f) {
                    # Calculate relative path within the source template
                    # Normalize file path to forward slashes for consistent regex matching on Windows
                    f_abs_normalized <- gsub("\\\\", "/", tools::file_path_as_absolute(f))
                    rel_path <- sub(paste0("^", source_abs_normalized, "/?"), "", f_abs_normalized)
                    dest_file <- file.path(paths$special$scripts, rel_path) # Destination in project scripts dir
                dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
                file.copy(f, dest_file, overwrite = TRUE)
            })
            } else {
                 message("No non-Rmd script files found in template directory.")
            }
        } else {
            message(paste("Script template directory not found at:", scripts_template_source, "- skipping script copy."))
        }
    } else {
         message("Skipping directory creation and script copying as existing structure is being reused.")
         # IMPORTANT: Ensure the timestamped directory for THIS session exists even when reusing.
         dir.create(paths$special$time, recursive = TRUE, showWarnings = FALSE)
    }

    # --- Define Final Directory Paths List for Environment ---
    # This part runs whether creating new or reusing existing structure.
    publication_graphs_dir <- file.path(paths$results$base, "publication_graphs")
    # qc_dir now refers to the timestamped directory base
    qc_dir <- paths$special$qc_base 
    
    # Create integration directory
    integration_dir <- file.path(base_dir, "integration")
    dir.create(integration_dir, recursive = TRUE, showWarnings = FALSE)

    dir_paths <- list(
        base_dir = base_dir,
        integration_dir = integration_dir,
        results_dir = paths$results$base,
        data_dir = paths$special$data,
        source_dir = paths$special$scripts, # This is the *project-specific* scripts dir
        de_output_dir = file.path(paths$results$base, "de_proteins"),
        publication_graphs_dir = publication_graphs_dir,
        timestamp = timestamp,
        qc_dir = qc_dir, # Base QC dir
        time_dir = paths$special$time, # Timestamped dir for current run output
        results_summary_dir = paths$results_summary$base,
        pathway_dir = file.path(paths$results$base, "pathway_enrichment"),
        protein_qc_dir = file.path(paths$results$base, "protein_qc"),
        peptide_qc_dir = file.path(paths$results$base, "peptide_qc"),
        clean_proteins_dir = file.path(paths$results$base, "clean_proteins"),
        qc_figures_dir = file.path(paths$results_summary$base, "QC_figures"),
        publication_figures_dir = file.path(paths$results_summary$base, "Publication_figures"),
        publication_tables_dir = file.path(paths$results_summary$base, "Publication_tables"),
        study_report_dir = file.path(paths$results_summary$base, "Study_report")
    )
    
    # --- Assign to Global Environment ---
    message("Assigning directory paths to global environment...")
    list2env(dir_paths, envir = .GlobalEnv)
    
    # --- Print Structure Summary ---
    cat("\nFinal Directory Structure Set in Global Environment:\n")
    print_paths <- sort(names(dir_paths)) # Sort for consistent order
    invisible(lapply(print_paths, function(name) {
        p <- dir_paths[[name]]
        # Check existence only for paths expected to be directories within the project
        is_project_dir <- startsWith(p, base_dir) && name != "base_dir" && name != "timestamp"
        
        if (is_project_dir && dir.exists(p)) {
            file_count <- length(list.files(p, recursive = TRUE, all.files = TRUE))
            cat(sprintf("%s = %s (%d files/dirs)\n", name, p, file_count))
        } else if (is.character(p)) {
             # Print non-dirs (like timestamp) or dirs outside base_dir (shouldn't happen here)
            cat(sprintf("%s = %s\n", name, p))
        } else {
            cat(sprintf("%s = %s [Non-character path]\n", name, p)) # Should not happen
        }
    }))
    
    # Return the list of paths invisibly
    invisible(dir_paths)
}

##################################################################################################################
# S4 class WorkflowArgs removed - replaced with createStudyParametersFile function

# Old complex createWorkflowArgsFromConfig function removed - replaced with simple wrapper below

#' @title Format Configuration List
#' @param config_list List of configuration parameters
#' @param indent Number of spaces for indentation
#' @export
formatConfigList <- function(config_list, indent = 0) {
    message(sprintf("--- DEBUG66: Entering formatConfigList with %d items, indent=%d ---", length(config_list), indent))
    output <- character()

    # Exclude internal_workflow_source_dir from printing
    names_to_process <- names(config_list)
    names_to_process <- names_to_process[names_to_process != "internal_workflow_source_dir"]
    message(sprintf("   DEBUG66: Processing %d items after exclusions", length(names_to_process)))

    # FUNCTIONAL APPROACH - NO FOR LOOPS that cause hanging - Works in Shiny AND .rmd
    output <- if (requireNamespace("purrr", quietly = TRUE)) {
        purrr::map_chr(names_to_process, function(name) {
        message(sprintf("   DEBUG66: Processing config item '%s'", name))
        value <- config_list[[name]]
        message(sprintf("   DEBUG66: Item '%s' class: %s", name, paste(class(value), collapse=", ")))
        
        # Skip core_utilisation, seqinr_obj, and complex objects from display
        if (name == "core_utilisation" || name == "seqinr_obj" ||
            any(class(value) %in% c("process", "R6", "multidplyr_cluster", "cluster", "SOCKcluster", "tbl_df", "tbl", "data.frame"))) {
            message(sprintf("   DEBUG66: Skipping '%s' due to complex class or large data frame", name))
            return("")  # Return empty string instead of next
        }

        # Format the name
        name_formatted <- gsub("\\.", " ", name)
        name_formatted <- gsub("_", " ", name_formatted)
        name_formatted <- tools::toTitleCase(name_formatted)

        # Handle different value types
        if (is.list(value)) {
            message(sprintf("   DEBUG66: '%s' is a list with %d elements", name, length(value)))
            if (length(value) > 0 && !is.null(names(value))) {
                output_lines <- paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ":")
                message(sprintf("   DEBUG66: Recursing into list '%s'", name))
                recursive_result <- formatConfigList(value, indent + 2)
                return(paste(c(output_lines, recursive_result), collapse = "\n"))
            } else if (length(value) > 0) { # Unnamed list, process elements functionally
                message(sprintf("   DEBUG66: '%s' is unnamed list, processing elements", name))
                header <- paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ":")
                
                                 # FUNCTIONAL processing of list elements - NO FOR LOOP
                 element_lines <- if (requireNamespace("purrr", quietly = TRUE)) {
                     purrr::imap_chr(value, function(item_val, item_idx) {
                         message(sprintf("   DEBUG66: Processing unnamed list item %d, class: %s", item_idx, paste(class(item_val), collapse=", ")))
                         if (is.atomic(item_val) && length(item_val) == 1) {
                             paste0(paste(rep(" ", indent + 2), collapse=""), "- ", as.character(item_val))
                         } else {
                             paste0(paste(rep(" ", indent + 2), collapse=""), "- [Complex List Element]")
                         }
                     })
                 } else {
                     # Base R fallback for list processing
                     sapply(seq_along(value), function(item_idx) {
                         item_val <- value[[item_idx]]
                         if (is.atomic(item_val) && length(item_val) == 1) {
                             paste0(paste(rep(" ", indent + 2), collapse=""), "- ", as.character(item_val))
                         } else {
                             paste0(paste(rep(" ", indent + 2), collapse=""), "- [Complex List Element]")
                         }
                     })
                 }
                return(paste(c(header, element_lines), collapse = "\n"))
            } else { # Empty list
                message(sprintf("   DEBUG66: '%s' is empty list", name))
                return(paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ": [Empty List]"))
            }
        } else {
            message(sprintf("   DEBUG66: '%s' is atomic/non-list, attempting conversion", name))
            
            # SAFE value display formatting
            value_display <- tryCatch({
                if (is.character(value) && length(value) > 5) {
                    paste(c(utils::head(value, 5), "..."), collapse = ", ")
                } else if (is.character(value) && length(value) > 1) {
                    paste(value, collapse = ", ")
                } else {
                    message(sprintf("   DEBUG66: About to convert '%s' to character", name))
                    result <- as.character(value)
                    message(sprintf("   DEBUG66: Conversion successful for '%s'", name))
                    if (length(result) == 0) "[Empty/NULL]" else result
                }
            }, error = function(e) {
                message(sprintf("   DEBUG66: Error converting '%s': %s", name, e$message))
                "[CONVERSION ERROR]"
            })
            
                         return(paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ": ", value_display))
         }
     }) |> 
     (\(x) x[x != ""])()  # Remove empty strings from skipped items
     } else {
         # Fallback using base R lapply if purrr not available
         result_list <- lapply(names_to_process, function(name) {
             value <- config_list[[name]]
             
             # Skip complex objects
             if (name == "core_utilisation" ||
                 any(class(value) %in% c("process", "R6", "multidplyr_cluster", "cluster", "SOCKcluster"))) {
                 return("")
             }
             
             # Format name
             name_formatted <- gsub("\\.", " ", name)
             name_formatted <- gsub("_", " ", name_formatted)
             name_formatted <- tools::toTitleCase(name_formatted)
             
             # Simple formatting for fallback
             value_display <- tryCatch({
                 if (is.list(value)) {
                     "[List Object]"
                 } else {
                     as.character(value)[1]  # Take first element only for safety
                 }
             }, error = function(e) "[CONVERSION ERROR]")
             
             paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ": ", value_display)
         })
         unlist(result_list[result_list != ""])
     }
     
     # Flatten the nested strings  
     output <- unlist(strsplit(output, "\n"))
    message(sprintf("--- DEBUG66: Exiting formatConfigList, returning %d lines ---", length(output)))
    return(output)
}

# S4 show method for WorkflowArgs removed - replaced with createStudyParametersFile function


##################################################################################################################

#' @title Copy Files to Results Summary and Show Status
#' @description Copies specified files to results summary directory and displays copy status.
#'              It now derives paths from a global `project_dirs` object using `omic_type` and `experiment_label`.
#' @param omic_type Character string specifying the omics type (e.g., "proteomics", "metabolomics").
#' @param experiment_label Character string specifying the experiment label.
#' @param contrasts_tbl A tibble containing contrast information (typically from the global environment if not passed directly).
#' @param design_matrix A data frame for the design matrix (typically from the global environment if not passed directly).
#' @param force Logical; if TRUE, skips user confirmation for backup (default: FALSE).
#' @param current_rmd Optional path to the current .Rmd file being worked on; if NULL (default),
#'                    the function will attempt to detect and save the currently active .Rmd file in RStudio.
#' @param project_dirs_object_name The name of the list object in the global environment that holds the directory structures (typically the output of setupDirectories). Defaults to "project_dirs".
#' @return Invisible list of failed copies for error handling.
#' @importFrom rlang abort
#' @importFrom tools file_path_as_absolute
#' @export
copyToResultsSummary <- function(omic_type,
                                 experiment_label,
                                 contrasts_tbl = NULL,
                                 design_matrix = NULL,
                                 force = FALSE, 
                                 current_rmd = NULL,
                                 project_dirs_object_name = "project_dirs") {

    # --- Start: Path Derivation and Validation --- #
    if (missing(omic_type) || !is.character(omic_type) || length(omic_type) != 1 || omic_type == "") {
        rlang::abort("`omic_type` must be a single non-empty character string.")
    }
    if (missing(experiment_label) || !is.character(experiment_label) || length(experiment_label) != 1 || experiment_label == "") {
        rlang::abort("`experiment_label` must be a single non-empty character string.")
    }
    
    # Use the new helper function with automatic fallback
    current_paths <- tryCatch({
        getProjectPaths(
            omic_type = omic_type,
            experiment_label = experiment_label,
            project_dirs_object_name = project_dirs_object_name
        )
    }, error = function(e) {
        rlang::abort(paste0("Failed to get project paths: ", e$message))
    })
    
    # Use the new helper function with automatic fallback
    current_paths <- tryCatch({
        getProjectPaths(
            omic_type = omic_type,
            experiment_label = experiment_label,
            project_dirs_object_name = project_dirs_object_name
        )
    }, error = function(e) {
        rlang::abort(paste0("Failed to get project paths: ", e$message))
    })
    
    # Validate that current_paths is a list and contains essential directory paths
    required_paths_in_current <- c("results_dir", "results_summary_dir", "publication_graphs_dir", 
                                   "time_dir", "qc_dir", "de_output_dir", "pathway_dir", "source_dir", "feature_qc_dir")
    if (!is.list(current_paths) || !all(required_paths_in_current %in% names(current_paths))) {
        missing_req <- setdiff(required_paths_in_current, names(current_paths))
        rlang::abort(paste0("Essential paths missing from project_dirs: ", paste(missing_req, collapse=", ")))
        rlang::abort(paste0("Essential paths missing from project_dirs: ", paste(missing_req, collapse=", ")))
    }
    # --- End: Path Derivation and Validation ---

    # Create descriptive label for logging (optional, for user-friendly messages)
    omic_label <- if (!is.null(experiment_label) && nzchar(experiment_label)) {
        paste0(omic_type, "_", experiment_label)
    } else {
        omic_type
    }

    # Create descriptive label for logging (optional, for user-friendly messages)
    omic_label <- if (!is.null(experiment_label) && nzchar(experiment_label)) {
        paste0(omic_type, "_", experiment_label)
    } else {
        omic_type
    }

    # Track failed copies
    failed_copies <- list()
    
    cat("\nRelevant directory paths:\n")
    cat("\nRelevant directory paths:\n")
    cat(sprintf("Results Dir: %s\n", current_paths$results_dir))
    cat(sprintf("Results Summary Dir: %s\n", current_paths$results_summary_dir))
    cat(sprintf("Publication Graphs Dir: %s\n", current_paths$publication_graphs_dir))
    cat(sprintf("Time Dir (current run): %s\n", current_paths$time_dir))
    cat(sprintf("DE Output Dir: %s\n", current_paths$de_output_dir))
    cat(sprintf("Pathway Dir: %s\n", current_paths$pathway_dir))
    cat(sprintf("Source (Scripts) Dir: %s\n", current_paths$source_dir))
    cat(sprintf("Feature QC Dir: %s\n", current_paths$feature_qc_dir))
    if(!is.null(current_paths$subfeature_qc_dir)) cat(sprintf("Sub-feature QC Dir: %s\n", current_paths$subfeature_qc_dir))
    cat("\n")

    # ROBUST: Try to get contrasts_tbl and design_matrix from environment first, then from files
    cat("Checking for required objects...\n")
    
    # Try contrasts_tbl from environment first
    if (is.null(contrasts_tbl)) {
        if (exists("contrasts_tbl", envir = parent.frame())) {
            contrasts_tbl <- get("contrasts_tbl", envir = parent.frame())
            cat("✓ Using 'contrasts_tbl' from calling environment\n")
        } else if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
            cat("✓ Using 'contrasts_tbl' from global environment\n")
        } else {
            # Fallback: try to load from file
            contrasts_file_options <- c(
                file.path(current_paths$source_dir, "contrasts_tbl.tab"),
                file.path(current_paths$source_dir, "contrast_strings.tab"),
                file.path(current_paths$source_dir, "contrasts.tab")
            )
            
            contrasts_loaded <- FALSE
            for (contrasts_file in contrasts_file_options) {
                if (file.exists(contrasts_file)) {
                    tryCatch({
                        contrasts_tbl <- readr::read_tsv(contrasts_file, show_col_types = FALSE)
                        cat(sprintf("✓ Loaded 'contrasts_tbl' from file: %s\n", basename(contrasts_file)))
                        contrasts_loaded <- TRUE
                        break
                    }, error = function(e) {
                        cat(sprintf("⚠ Failed to load contrasts from %s: %s\n", basename(contrasts_file), e$message))
                    })
                }
            }
            
            if (!contrasts_loaded) {
                cat("✗ 'contrasts_tbl' not found in environment or files\n")
            }
        }
    } else {
        cat("✓ Using provided 'contrasts_tbl' parameter\n")
    }
    
    # Try design_matrix from environment first
    if (is.null(design_matrix)) {
        if (exists("design_matrix", envir = parent.frame())) {
            design_matrix <- get("design_matrix", envir = parent.frame())
            cat("✓ Using 'design_matrix' from calling environment\n")
        } else if (exists("design_matrix", envir = .GlobalEnv)) {
            design_matrix <- get("design_matrix", envir = .GlobalEnv)
            cat("✓ Using 'design_matrix' from global environment\n")
        } else {
            # Fallback: try to load from file
            design_matrix_file <- file.path(current_paths$source_dir, "design_matrix.tab")
            
            if (file.exists(design_matrix_file)) {
                tryCatch({
                    design_matrix <- readr::read_tsv(design_matrix_file, show_col_types = FALSE)
                    cat(sprintf("✓ Loaded 'design_matrix' from file: %s\n", basename(design_matrix_file)))
                }, error = function(e) {
                    cat(sprintf("✗ Failed to load design_matrix from %s: %s\n", basename(design_matrix_file), e$message))
                    design_matrix <- NULL
                })
            } else {
                cat(sprintf("✗ 'design_matrix' not found in environment or at expected file location: %s\n", design_matrix_file))
            }
        }
    } else {
        cat("✓ Using provided 'design_matrix' parameter\n")
    }
    
    # Handle current Rmd file
    if (is.null(current_rmd) && exists("rstudioapi") && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
        context <- rstudioapi::getActiveDocumentContext()
        if (!is.null(context) && !is.null(context$path) && grepl("\\.rmd$", context$path, ignore.case = TRUE)) {
            current_rmd <- context$path
            cat(sprintf("Detected active .Rmd file: %s\n", current_rmd))
        }
    }
    
    if (!is.null(current_rmd) && file.exists(current_rmd)) {
        tryCatch({
            if (exists("rstudioapi") && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
                context <- rstudioapi::getActiveDocumentContext()
                if (!is.null(context) && !is.null(context$path) && normalizePath(tools::file_path_as_absolute(context$path)) == normalizePath(tools::file_path_as_absolute(current_rmd))) {
                    rstudioapi::documentSave(context$id)
                    cat(sprintf("Saved current Rmd file: %s\n", current_rmd))
                } else {
                    docs <- rstudioapi::getSourceEditorContexts()
                    for (doc in docs) {
                        if (!is.null(doc$path) && normalizePath(tools::file_path_as_absolute(doc$path)) == normalizePath(tools::file_path_as_absolute(current_rmd))) {
                            rstudioapi::documentSave(doc$id)
                            cat(sprintf("Saved Rmd file: %s\n", current_rmd))
                            break
                        }
                    }
                }
            }
            dest_file <- file.path(current_paths$source_dir, basename(current_rmd))
            file.copy(current_rmd, dest_file, overwrite = TRUE)
            cat(sprintf("Copied Rmd file to project scripts directory: %s\n", dest_file))
        }, error = function(e) {
            warning(sprintf("Failed to save/copy Rmd file: %s", e$message))
            failed_copies[[length(failed_copies) + 1]] <- list(type = "rmd_copy", source = current_rmd, destination = current_paths$source_dir, display_name = "Current Rmd File", error = e$message)
        })
    }
    
    # Define target directories using derived paths
    pub_tables_dir <- file.path(current_paths$results_summary_dir, "Publication_tables")
    
    # Ensure results_summary_dir exists before checking its contents (it should, from setupDirectories)
    if (!dir.exists(current_paths$results_summary_dir)) {
        dir.create(current_paths$results_summary_dir, recursive = TRUE, showWarnings = FALSE)
        logger::log_info("Created results_summary_dir as it was missing: {current_paths$results_summary_dir}")
    }
    
    # Check if the results_summary_dir has any content
    contents_of_summary_dir <- list.files(current_paths$results_summary_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    
    if (length(contents_of_summary_dir) > 0) {
        logger::log_info("Results summary directory for {omic_label} ({current_paths$results_summary_dir}) contains existing files/folders.")
        backup_dirname <- paste0(omic_label, "_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        logger::log_info("Results summary directory for {omic_label} ({current_paths$results_summary_dir}) contains existing files/folders.")
        backup_dirname <- paste0(omic_label, "_backup_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        backup_dir <- file.path(dirname(current_paths$results_summary_dir), backup_dirname)
        
        should_proceed_with_backup <- if (!force) {
            cat(sprintf("\\nResults summary directory for %s contains content:\\n- %s\\n", omic_label, current_paths$results_summary_dir))
            cat(sprintf("\\nResults summary directory for %s contains content:\\n- %s\\n", omic_label, current_paths$results_summary_dir))
            repeat {
                response <- readline(prompt = "Do you want to backup existing directory and proceed by overwriting? (y/n): ")
                response <- tolower(substr(response, 1, 1))
                if (response %in% c("y", "n")) break
                cat("Please enter 'y' or 'n'\\n")
            }
            response == "y"
        } else {
            logger::log_info("Force mode enabled - backing up and proceeding with overwrite for {omic_label}...")
            logger::log_info("Force mode enabled - backing up and proceeding with overwrite for {omic_label}...")
            TRUE
        }
        
        if (!should_proceed_with_backup) {
            logger::log_info("Overwrite of {current_paths$results_summary_dir} for {omic_label} cancelled by user. No backup made, original files untouched.")
            logger::log_info("Overwrite of {current_paths$results_summary_dir} for {omic_label} cancelled by user. No backup made, original files untouched.")
            # Return a list indicating cancellation, which can be checked by the caller.
            return(invisible(list(status="cancelled", omic_key=omic_label, message=paste0("Backup and overwrite for ", omic_label, " cancelled by user."))))
            return(invisible(list(status="cancelled", omic_key=omic_label, message=paste0("Backup and overwrite for ", omic_label, " cancelled by user."))))
        }
        
        # Proceed with backup and clearing
        if (!dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(backup_dir)) {
             logger::log_warn("Failed to create backup directory: {backup_dir} for {omic_label}. Original directory will not be cleared.")
             logger::log_warn("Failed to create backup directory: {backup_dir} for {omic_label}. Original directory will not be cleared.")
             failed_copies[[length(failed_copies) + 1]] <- list(type = "backup_dir_creation", path = backup_dir, error = "Failed to create backup directory")
             # Do not proceed with unlink if backup dir creation fails
        } else {
            items_to_backup <- list.files(current_paths$results_summary_dir, full.names = TRUE, all.files = TRUE, no.. = TRUE)
            backup_copy_successful <- TRUE
            if (length(items_to_backup) > 0) {
                dir.create(backup_dir, showWarnings = FALSE, recursive = TRUE) # Ensure backup_dir itself exists
                copy_results <- file.copy(from = items_to_backup, to = backup_dir, overwrite = TRUE, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE)
                if (!all(copy_results)) {
                    backup_copy_successful <- FALSE
                    logger::log_warn("Not all items were successfully copied from {current_paths$results_summary_dir} to backup directory {backup_dir}.")
                }
            }

            backup_has_content <- length(list.files(backup_dir, recursive=TRUE, all.files=TRUE, no..=TRUE)) > 0
            
            if (backup_copy_successful && (backup_has_content || length(contents_of_summary_dir) == 0) ) {
                logger::log_info("Successfully backed up content of {current_paths$results_summary_dir} to: {backup_dir}")
                backup_info <- data.frame(original_dir = current_paths$results_summary_dir, backup_time = Sys.time(), omic_key = omic_label, stringsAsFactors = FALSE)
                backup_info <- data.frame(original_dir = current_paths$results_summary_dir, backup_time = Sys.time(), omic_key = omic_label, stringsAsFactors = FALSE)
                tryCatch(
                    write.table(backup_info, file = file.path(backup_dir, "backup_info.txt"), sep = "\\t", row.names = FALSE, quote = FALSE),
                    error = function(e) logger::log_warn("Failed to write backup_info.txt: {e$message}")
                )
                
                logger::log_info("Clearing original results_summary_dir: {current_paths$results_summary_dir}")
                unlink_success <- tryCatch({
                    unlink(current_paths$results_summary_dir, recursive = TRUE, force = TRUE)
                    TRUE # Return TRUE on success
                }, error = function(e) {
                    logger::log_error("Error unlinking {current_paths$results_summary_dir}: {e$message}")
                    FALSE # Return FALSE on error
                })

                if(unlink_success) {
                    if (!dir.create(current_paths$results_summary_dir, recursive = TRUE, showWarnings = FALSE) && !dir.exists(current_paths$results_summary_dir)) {
                        warning_msg <- sprintf("CRITICAL: Failed to recreate results_summary_dir %s after backup and unlink. Subsequent operations will likely fail.", current_paths$results_summary_dir)
                        logger::log_error(warning_msg)
                        failed_copies[[length(failed_copies) + 1]] <- list(type = "critical_dir_recreation", path = current_paths$results_summary_dir, error = warning_msg)
                        # Potentially stop execution or return an error status here
                    } else {
                         logger::log_info("Successfully cleared and recreated results_summary_dir: {current_paths$results_summary_dir}")
                    }
                } else {
                     warning_msg <- sprintf("Failed to clear original results_summary_dir %s after backup. Subsequent operations may overwrite or mix files.", current_paths$results_summary_dir)
                     logger::log_warn(warning_msg)
                     failed_copies[[length(failed_copies) + 1]] <- list(type = "dir_clear_failure", path = current_paths$results_summary_dir, error = warning_msg)
                }
            } else {
                logger::log_warn("Failed to copy all items to backup for {omic_label}, or backup is unexpectedly empty. Original directory {current_paths$results_summary_dir} was NOT cleared.")
                logger::log_warn("Failed to copy all items to backup for {omic_label}, or backup is unexpectedly empty. Original directory {current_paths$results_summary_dir} was NOT cleared.")
                failed_copies[[length(failed_copies) + 1]] <- list(type = "backup_content_copy", source = current_paths$results_summary_dir, destination = backup_dir, error = "Failed to copy items to backup or backup empty; original not cleared")
            }
        }
    } else {
        logger::log_info("Results summary directory for {omic_label} ({current_paths$results_summary_dir}) is empty. No backup needed. Proceeding to create subdirectories.")
        logger::log_info("Results summary directory for {omic_label} ({current_paths$results_summary_dir}) is empty. No backup needed. Proceeding to create subdirectories.")
        # Ensure the main directory exists (it should if we got here and it was empty, or it was just created if missing)
        if (!dir.exists(current_paths$results_summary_dir)) {
             if(!dir.create(current_paths$results_summary_dir, recursive = TRUE, showWarnings = FALSE)) {
                 warning_msg <- sprintf("CRITICAL: Failed to create initially empty results_summary_dir %s. Subsequent operations may fail.", current_paths$results_summary_dir)
                 logger::log_error(warning_msg)
                 failed_copies[[length(failed_copies) + 1]] <- list(type = "critical_empty_dir_creation", path = current_paths$results_summary_dir, error = warning_msg)
             }
        }
    }
    
    summary_subdirs <- c("QC_figures", "Publication_figures", "Publication_tables", "Study_report")
    sapply(summary_subdirs, \(subdir) {
        dir_path <- file.path(current_paths$results_summary_dir, subdir)
        if (!dir.create(dir_path, recursive = TRUE, showWarnings = FALSE) && !dir.exists(dir_path)) {
            warning(sprintf("Failed to create directory: %s", dir_path))
            failed_copies[[length(failed_copies) + 1]] <- list(type = "directory_creation", path = dir_path, error = "Failed to create directory")
        }
    })

    # Define files to copy - make paths relative to current_paths elements
    files_to_copy <- list(
        list(source = file.path(current_paths$time_dir, "12_correlation_filtered_combined_plots.png"), dest = "QC_figures", is_dir = FALSE, display_name = "Correlation Filtered Plots"),
        list(source = file.path(current_paths$feature_qc_dir, "composite_QC_figure.pdf"), dest = "QC_figures", is_dir = FALSE, display_name = "Composite QC (PDF)", new_name = paste0(omic_type, "_composite_QC_figure.pdf")),
        list(source = file.path(current_paths$feature_qc_dir, "composite_QC_figure.png"), dest = "QC_figures", is_dir = FALSE, display_name = "Composite QC (PNG)", new_name = paste0(omic_type, "_composite_QC_figure.png")),
        list(source = file.path(current_paths$publication_graphs_dir, "Interactive_Volcano_Plots"), dest = "Publication_figures/Interactive_Volcano_Plots", is_dir = TRUE, display_name = "Interactive Volcano Plots"),
        list(source = file.path(current_paths$publication_graphs_dir, "NumSigDeMolecules"), dest = "Publication_figures/NumSigDeMolecules", is_dir = TRUE, display_name = "Num Sig DE Molecules"),
        list(source = file.path(current_paths$publication_graphs_dir, "Volcano_Plots"), dest = "Publication_figures/Volcano_Plots", is_dir = TRUE, display_name = "Volcano Plots"),
        list(source = current_paths$pathway_dir, dest = "Publication_figures/Enrichment_Plots", is_dir = TRUE, display_name = "Pathway Enrichment Plots"),
        list(source = "contrasts_tbl", dest = "Study_report", type = "object", save_as = "contrasts_tbl.tab", display_name = "Contrasts Table"),
        list(source = "design_matrix", dest = "Study_report", type = "object", save_as = "design_matrix.tab", display_name = "Design Matrix"),
        list(source = file.path(current_paths$source_dir, "study_parameters.txt"), dest = "Study_report", is_dir = FALSE, display_name = "Study Parameters")
    )

    # Logic to find and copy the num_sig_de_molecules tab file to Publication_tables
    num_sig_dir <- file.path(current_paths$publication_graphs_dir, "NumSigDeMolecules")
    if (dir.exists(num_sig_dir)) {
        # Look for the tab file
        sig_files <- list.files(num_sig_dir, pattern = "_num_sig_de_molecules\\.tab$", full.names = TRUE)
        if (length(sig_files) > 0) {
            # Use the first match
            files_to_copy <- c(files_to_copy, list(
                list(source = sig_files[1], dest = "Publication_tables", is_dir = FALSE, display_name = "Num Sig DE Molecules Tab", new_name = paste0("de_", omic_type, "_num_sig_de_molecules.tab"))
            ))
        }
    }

    if (omic_type == "proteomics") {
        # Check for both RUV and non-RUV filenames (conditional on whether RUV was applied)
        ruv_tsv_path <- file.path(current_paths$feature_qc_dir, "ruv_normalised_results_cln_with_replicates.tsv")
        norm_tsv_path <- file.path(current_paths$feature_qc_dir, "normalised_results_cln_with_replicates.tsv")
        ruv_rds_path <- file.path(current_paths$feature_qc_dir, "ruv_normalised_results_cln_with_replicates.RDS")
        norm_rds_path <- file.path(current_paths$feature_qc_dir, "normalised_results_cln_with_replicates.RDS")
        
        # Determine which files exist and use appropriate naming
        if (file.exists(ruv_tsv_path)) {
            # RUV was applied - use RUV filenames
            files_to_copy <- c(files_to_copy, list(
                list(source = ruv_tsv_path, dest = "Publication_tables", is_dir = FALSE, display_name = "RUV Normalized Results TSV", new_name = "RUV_normalised_results.tsv"),
                list(source = ruv_rds_path, dest = "Publication_tables", is_dir = FALSE, display_name = "RUV Normalized Results RDS", new_name = "ruv_normalised_results.RDS")
            ))
            cat("COPY: Using RUV-normalized file names (RUV was applied)\n")
        } else if (file.exists(norm_tsv_path)) {
            # RUV was skipped - use normalized filenames
            files_to_copy <- c(files_to_copy, list(
                list(source = norm_tsv_path, dest = "Publication_tables", is_dir = FALSE, display_name = "Normalized Results TSV", new_name = "normalised_results.tsv"),
                list(source = norm_rds_path, dest = "Publication_tables", is_dir = FALSE, display_name = "Normalized Results RDS", new_name = "normalised_results.RDS")
            ))
            cat("COPY: Using normalized file names (RUV was skipped)\n")
        } else {
            # Neither exists - try RUV path anyway (will fail gracefully)
            files_to_copy <- c(files_to_copy, list(
                list(source = ruv_tsv_path, dest = "Publication_tables", is_dir = FALSE, display_name = "Normalized Results TSV", new_name = "normalised_results.tsv"),
                list(source = ruv_rds_path, dest = "Publication_tables", is_dir = FALSE, display_name = "Normalized Results RDS", new_name = "normalised_results.RDS")
            ))
            cat("COPY: Warning - neither RUV nor normalized files found, attempting RUV paths\n")
        }
    } # Add else if for other omic-specific files if necessary
    
    # Excel files paths
    de_results_excel_path <- file.path(pub_tables_dir, paste0("DE_results_", omic_type, ".xlsx"))
    enrichment_excel_path <- file.path(pub_tables_dir, paste0("Pathway_enrichment_results_", omic_type, ".xlsx"))

    # Create combined DE workbook
    de_wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(de_wb, "DE_Results_Index")
    de_index_data <- data.frame(Sheet = character(), Description = character(), stringsAsFactors = FALSE)
    de_files <- list.files(path = current_paths$de_output_dir, pattern = paste0("de_.+_long_annot\\.xlsx$"), full.names = TRUE) # Changed from \\w+ to .+ to allow hyphens
    
    purrr::imap(de_files, \(file, idx) {
        sheet_name <- sprintf("DE_Sheet%d", idx)
        base_name <- basename(file) |> stringr::str_remove("^de_") |> stringr::str_remove("_long_annot\\.xlsx$")
        de_index_data <<- rbind(de_index_data, data.frame(Sheet = sheet_name, Description = base_name, stringsAsFactors = FALSE))
        data_content <- tryCatch(openxlsx::read.xlsx(file), error = function(e) NULL)
        if (!is.null(data_content)) {
             openxlsx::addWorksheet(de_wb, sheet_name)
             openxlsx::writeData(de_wb, sheet_name, data_content)
        } else {
            warning(paste0("Failed to read DE Excel file: ", file))
            failed_copies[[length(failed_copies) + 1]] <- list(type = "de_excel_read", source = file, error = "Failed to read")
        }
    })
    openxlsx::writeData(de_wb, "DE_Results_Index", de_index_data)
    openxlsx::setColWidths(de_wb, "DE_Results_Index", cols = 1:2, widths = c(15, 50))
    openxlsx::addStyle(de_wb, "DE_Results_Index", style = openxlsx::createStyle(textDecoration = "bold"), rows = 1, cols = 1:2)

    # Create combined Enrichment workbook
    enrichment_wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(enrichment_wb, "Enrichment_Index")
    enrichment_index_data <- data.frame(Sheet = character(), Contrast = character(), Direction = character(), stringsAsFactors = FALSE)
    enrichment_files <- list.files(path = current_paths$pathway_dir, pattern = "_enrichment_results\\.tsv$", full.names = TRUE)
    
    purrr::imap(enrichment_files, \(file, idx) {
        base_name <- basename(file) |> stringr::str_remove("_enrichment_results\\.tsv$")
        direction <- ifelse(stringr::str_ends(base_name, "_up"), "up", ifelse(stringr::str_ends(base_name, "_down"), "down", "unknown"))
        contrast_label <- stringr::str_replace(base_name, "_(up|down)$", "")
        sheet_name <- sprintf("Enrich_Sheet%d", idx)
        enrichment_index_data <<- rbind(enrichment_index_data, data.frame(Sheet = sheet_name, Contrast = contrast_label, Direction = direction, stringsAsFactors = FALSE))
        data_content <- tryCatch(readr::read_tsv(file, show_col_types = FALSE), error = function(e) NULL)
         if (!is.null(data_content)) {
            openxlsx::addWorksheet(enrichment_wb, sheet_name)
            openxlsx::writeData(enrichment_wb, sheet_name, data_content)
        } else {
            warning(paste0("Failed to read Enrichment TSV file: ", file))
            failed_copies[[length(failed_copies) + 1]] <- list(type = "enrichment_tsv_read", source = file, error = "Failed to read")
        }
    })
    openxlsx::writeDataTable(enrichment_wb, "Enrichment_Index", enrichment_index_data, tableStyle = "TableStyleLight9", headerStyle = openxlsx::createStyle(textDecoration = "bold"), withFilter = TRUE)
    openxlsx::writeData(enrichment_wb, "Enrichment_Index", data.frame(Note = "Contrast represents the comparison (e.g., Group1_minus_Group2). Direction shows up-regulated or down-regulated genes."), startRow = nrow(enrichment_index_data) + 3)
    
    dir.create(pub_tables_dir, recursive = TRUE, showWarnings = FALSE)
    tryCatch({ 
      openxlsx::saveWorkbook(de_wb, de_results_excel_path, overwrite = TRUE)
      cat(paste("Successfully saved DE results to:", de_results_excel_path, "\n"))
    }, error = function(e) { 
      failed_copies[[length(failed_copies) + 1]] <- list(type = "workbook_save", path = de_results_excel_path, error = e$message)
    })
    
    tryCatch({ 
      openxlsx::saveWorkbook(enrichment_wb, enrichment_excel_path, overwrite = TRUE)
      cat(paste("Successfully saved Enrichment results to:", enrichment_excel_path, "\n"))
    }, error = function(e) { 
      failed_copies[[length(failed_copies) + 1]] <- list(type = "workbook_save", path = enrichment_excel_path, error = e$message)
    })

    cat("\nCopying individual files/folders to Results Summary for ", omic_label, "...\n")
    cat("\nCopying individual files/folders to Results Summary for ", omic_label, "...\n")
    cat("==============================================================\n\n")

    files_to_copy |>
        lapply(\(file_spec) {
            dest_dir_final <- gsub("//", "/", file.path(current_paths$results_summary_dir, file_spec$dest))
            file_spec$source <- gsub("//", "/", file_spec$source)
            source_display <- file_spec$source # For display purposes

            copy_success <- TRUE
            error_msg <- NULL

            if (!is.null(file_spec$type) && file_spec$type == "object") {
                # ENHANCED: Use robust object sourcing - check environment first, then file
                obj <- NULL
                source_exists <- FALSE
                
                # Try to get object from environments first
                if (!is.null(get0(file_spec$source, envir = parent.frame()))) {
                    obj <- get(file_spec$source, envir = parent.frame())
                    source_exists <- TRUE
                    cat(sprintf("    ✓ Found '%s' in calling environment\n", file_spec$source))
                } else if (!is.null(get0(file_spec$source, envir = .GlobalEnv))) {
                    obj <- get(file_spec$source, envir = .GlobalEnv)
                    source_exists <- TRUE
                    cat(sprintf("    ✓ Found '%s' in global environment\n", file_spec$source))
                } else {
                    # Fallback: try to load from file
                    if (file_spec$source == "design_matrix") {
                        design_matrix_file <- file.path(current_paths$source_dir, "design_matrix.tab")
                        if (file.exists(design_matrix_file)) {
                            tryCatch({
                                obj <- readr::read_tsv(design_matrix_file, show_col_types = FALSE)
                                source_exists <- TRUE
                                cat(sprintf("    ✓ Loaded '%s' from file: %s\n", file_spec$source, basename(design_matrix_file)))
                            }, error = function(e) {
                                error_msg <- sprintf("Failed to load %s from file %s: %s", file_spec$source, basename(design_matrix_file), e$message)
                                cat(sprintf("    ✗ %s\n", error_msg))
                            })
                        } else {
                            error_msg <- sprintf("Object '%s' not found in environment and file not found: %s", file_spec$source, design_matrix_file)
                        }
                    } else if (file_spec$source == "contrasts_tbl") {
                        contrasts_file_options <- c(
                            file.path(current_paths$source_dir, "contrasts_tbl.tab"),
                            file.path(current_paths$source_dir, "contrast_strings.tab"),
                            file.path(current_paths$source_dir, "contrasts.tab")
                        )
                        
                        for (contrasts_file in contrasts_file_options) {
                            if (file.exists(contrasts_file)) {
                                tryCatch({
                                    obj <- readr::read_tsv(contrasts_file, show_col_types = FALSE)
                                    source_exists <- TRUE
                                    cat(sprintf("    ✓ Loaded '%s' from file: %s\n", file_spec$source, basename(contrasts_file)))
                                    break
                                }, error = function(e) {
                                    cat(sprintf("    ⚠ Failed to load contrasts from %s: %s\n", basename(contrasts_file), e$message))
                                })
                            }
                        }
                        
                        if (!source_exists) {
                            error_msg <- sprintf("Object '%s' not found in environment and no readable contrasts files found", file_spec$source)
                        }
                    } else {
                        # For other objects, keep original behavior
                        error_msg <- sprintf("Object '%s' not found in parent/global environment", file_spec$source)
                    }
                }
                
                if (!source_exists) {
                    if (is.null(error_msg)) {
                        error_msg <- sprintf("Object '%s' not found in environment or files", file_spec$source)
                    }
                }
            } else {
                source_exists <- if (file_spec$is_dir) dir.exists(file_spec$source) else file.exists(file_spec$source)
                if (!source_exists) error_msg <- sprintf("Source %s not found: %s", if(file_spec$is_dir) "directory" else "file", file_spec$source)
            }

            if (source_exists) {
                dir.create(dest_dir_final, recursive = TRUE, showWarnings = FALSE)
                if (!is.null(file_spec$type) && file_spec$type == "object") {
                    tryCatch({
                        # Use the already-loaded object instead of getting from environment again
                        if (is.null(obj)) {
                            obj <- get(file_spec$source, envir = parent.frame())  # Fallback for other objects
                        }
                        dest_path <- file.path(dest_dir_final, file_spec$save_as)
                        write.table(obj, file = dest_path, sep = "\t", row.names = FALSE, quote = FALSE)
                        if (!file.exists(dest_path) || (file.exists(dest_path) && file.size(dest_path) == 0 && nrow(obj) > 0) ) {
                            copy_success <- FALSE
                            error_msg <- "Failed to write object or file is empty"
                        }
                    }, error = function(e) { copy_success <<- FALSE; error_msg <<- sprintf("Error writing object: %s", e$message) })
                } else if (file_spec$is_dir) {
                    # REWRITTEN: Robust directory copy logic using relative paths
                    all_source_files_rel <- list.files(file_spec$source, recursive = TRUE, full.names = FALSE)
                    
                    if (length(all_source_files_rel) > 0) {
                        failed_in_dir <- 0
                        
                        copy_results <- sapply(all_source_files_rel, function(rel_path) {
                            src_file <- file.path(file_spec$source, rel_path)
                            dest_file <- file.path(dest_dir_final, rel_path)
                            
                            # Ensure destination directory exists
                            dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
                            
                            # Copy the file
                            success <- file.copy(src_file, dest_file, overwrite = TRUE)
                            if (!success) failed_in_dir <<- failed_in_dir + 1
                            return(success)
                        })
                        
                        if (!all(copy_results)) {
                            copy_success <- FALSE
                            error_msg <- sprintf("Failed to copy %d/%d files from %s", failed_in_dir, length(all_source_files_rel), file_spec$source)
                        } else {
                            # Verify counts match roughly
                             dest_files_count <- length(list.files(dest_dir_final, recursive = TRUE))
                             if (dest_files_count < length(all_source_files_rel)) {
                                 # This might happen if overwrite failed silently or something odd, but file.copy returned true?
                                 # We rely on file.copy return value mostly.
                             }
                        }
                    } else {
                        message(sprintf("Source directory %s is empty. Nothing to copy.", file_spec$source))
                    }
                } else {
                    dest_path <- file.path(dest_dir_final, if (!is.null(file_spec$new_name)) file_spec$new_name else basename(file_spec$source))
                    if (!file.copy(from = file_spec$source, to = dest_path, overwrite = TRUE)) {
                        copy_success <- FALSE
                        error_msg <- "Failed to copy file"
                    } else if (file.exists(file_spec$source) && file.exists(dest_path) && file.size(file_spec$source) != file.size(dest_path)) {
                        copy_success <- FALSE
                        error_msg <- sprintf("File size mismatch: source=%d, dest=%d bytes", file.size(file_spec$source), file.size(dest_path))
                    }
                }
            }

            if (!source_exists || !copy_success) {
                failed_copies[[length(failed_copies) + 1]] <- list(type = if(!is.null(file_spec$type) && file_spec$type == "object") "object" else if(file_spec$is_dir) "directory" else "file", source = source_display, destination = dest_dir_final, display_name = file_spec$display_name, error = error_msg)
            }
            cat(sprintf("%-35s [%s -> %s] %s\n", file_spec$display_name, if(source_exists) "✓" else "✗", if(copy_success && source_exists) "✓" else "✗", if(!is.null(file_spec$type) && file_spec$type == "object") "Object" else if(file_spec$is_dir) "Directory" else "File"))
            if (!is.null(error_msg)) cat(sprintf("%35s Error: %s\n", "", error_msg))
        })

    cat("\nLegend: ✓ = exists/success, ✗ = missing/failed\n")
    cat("Arrow (->) shows source -> destination status\n")

    if (length(failed_copies) > 0) {
        cat("\nFailed Copies Summary:\n")
        cat("\nFailed Copies Summary:\n")
        cat("=====================================\n")
        lapply(failed_copies, function(failure) {
            cat(sprintf("\n%s: %s\n", failure$display_name, failure$error))
            cat(sprintf("  Source: %s\n", as.character(failure$source))) # Ensure source is char
            cat(sprintf("  Destination Attempted: %s\n", failure$destination))
        })
        warning(sprintf("%d files/objects/directories failed to copy correctly", length(failed_copies)))
        warning(sprintf("%d files/objects/directories failed to copy correctly", length(failed_copies)))
    }
    cat("--- End of copyToResultsSummary ---\n\n")
    cat("--- End of copyToResultsSummary ---\n\n")
    invisible(failed_copies)
}


#' Update Missing Value Parameters in Configuration List and S4 Object
#'
#' @description
#' Automatically calculates and updates the missing value filtering parameters in the configuration list
#' and S4 object @args based on the experimental design matrix. The function ensures at least a specified 
#' number of groups have sufficient quantifiable values for analysis.
#'
#' @param theObject An S4 object containing the experimental data and design matrix.
#' @param min_reps_per_group Integer specifying the minimum number of replicates required in each passing group.
#'                          If a group has fewer total replicates than this value, the minimum is adjusted.
#' @param min_groups Integer specifying the minimum number of groups required to have sufficient
#'                  quantifiable values. Default is 2.
#' @param config_list_name The name of the global config list variable (defaults to "config_list").
#' @param env The environment where the global config list resides (defaults to .GlobalEnv).
#'
#' @return Updated S4 object with synchronized @args and global config_list
#'
#' @details
#' The function calculates:
#' - groupwise_percentage_cutoff: Based on minimum required replicates per group
#' - max_groups_percentage_cutoff: Based on minimum required groups
#' 
#' Both the S4 object's @args slot and the global config_list are updated to maintain synchronization.
#'
#' @examples
#' \dontrun{
#' protein_log2_quant_cln <- updateMissingValueParameters(
#'   theObject = protein_log2_quant_cln, 
#'   min_reps_per_group = 2, 
#'   min_groups = 2
#' )
#' }
#'
#' @export
updateMissingValueParameters <- function(theObject, 
                                       min_reps_per_group = 2, 
                                       min_groups = 2,
                                       config_list_name = "config_list",
                                       env = .GlobalEnv) {
    
    # --- Input Validation ---
    if (!isS4(theObject)) {
        stop("'theObject' must be an S4 object.")
    }
    if (!"design_matrix" %in% methods::slotNames(theObject)) {
        stop("'theObject' must have a '@design_matrix' slot.")
    }
    if (!"args" %in% methods::slotNames(theObject)) {
        stop("'theObject' must have an '@args' slot.")
    }
    if (!exists(config_list_name, envir = env)) {
        stop("Global config list '", config_list_name, "' not found in the specified environment.")
    }
    
    # Extract design matrix from S4 object
    design_matrix <- theObject@design_matrix
    
    # Retrieve the global config list
    config_list <- get(config_list_name, envir = env)
    
    # Get number of replicates per group
    reps_per_group_tbl <- design_matrix |>
        group_by(group) |>
        summarise(n_reps = n()) |>
        ungroup()
    
    # Get total number of groups
    total_groups <- nrow(reps_per_group_tbl)
    
    if (min_groups > total_groups) {
        stop("min_groups cannot be larger than total number of groups")
    }
    
    # Calculate percentage missing allowed for each group
    group_thresholds <- reps_per_group_tbl |>
        mutate(
            adjusted_min_reps = pmin(n_reps, min_reps_per_group),
            max_missing = n_reps - adjusted_min_reps,
            missing_percent = round((max_missing / n_reps) * 100, 3)
        )
    
    # Use a consistent percentage threshold across all groups
    # Take the maximum percentage to ensure all groups meet minimum requirements
    groupwise_cutoff <- max(group_thresholds$missing_percent)
    
    # Calculate maximum failing groups allowed
    max_failing_groups <- total_groups - min_groups
    max_groups_cutoff <- round((max_failing_groups / total_groups) * 100, 3)
    
    # Update global config_list
    if (is.null(config_list$removeRowsWithMissingValuesPercent)) {
        config_list$removeRowsWithMissingValuesPercent <- list()
    }
    config_list$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff <- groupwise_cutoff
    config_list$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff <- max_groups_cutoff
    config_list$removeRowsWithMissingValuesPercent$proteins_intensity_cutoff_percentile <- 1
    
    # Assign updated config_list back to global environment
    assign(config_list_name, config_list, envir = env)
    
    # Update S4 object @args to match config_list
    if (is.null(theObject@args)) {
        theObject@args <- list()
    }
    if (is.null(theObject@args$removeRowsWithMissingValuesPercent)) {
        theObject@args$removeRowsWithMissingValuesPercent <- list()
    }
    
    theObject@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff <- groupwise_cutoff
    theObject@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff <- max_groups_cutoff
    theObject@args$removeRowsWithMissingValuesPercent$proteins_intensity_cutoff_percentile <- 1
    
    # Create informative message
    basic_msg <- sprintf(
        "Updated missing value parameters in both config_list and S4 object @args:
    - Requiring at least %d replicates in each passing group (varies by group size)
    - Requiring at least %d out of %d groups to pass (%.3f%% failing groups allowed)
    - groupwise_percentage_cutoff set to %.3f%%
    - max_groups_percentage_cutoff set to %.3f%%",
        min_reps_per_group,
        min_groups,
        total_groups,
        max_groups_cutoff,
        groupwise_cutoff,
        max_groups_cutoff
    )
    
    # Add details for each group
    group_detail_strings <- group_thresholds |>
        mutate(
            detail = sprintf("    Group %s: %d out of %d replicates required (%.3f%% missing allowed)",
                             group, adjusted_min_reps, n_reps, missing_percent)
        ) |>
        dplyr::pull(detail)
        
    group_details <- paste(group_detail_strings, collapse = "\n")
    
    # Print the message
    message(paste(basic_msg, "\n\nGroup details:", group_details, sep = "\n"))
    message("✅ S4 object @args and global config_list are now synchronized")
    
    return(theObject)
}

##################################################################################################################
#' @export
updateRuvParameters <- function(config_list, best_k, control_genes_index, percentage_as_neg_ctrl) {
  config_list$ruvParameters$best_k <- best_k
  config_list$ruvParameters$num_neg_ctrl <- length(control_genes_index)
  config_list$ruvParameters$percentage_as_neg_ctrl <- percentage_as_neg_ctrl
  
  # Print the number of negative controls (as in the original code)
  config_list$ruvParameters$num_neg_ctrl
  
  # Return the updated config list
  return(config_list)
}

#' @title Download Report Template from GitHub
#' @description Downloads a report template from the MultiScholaR GitHub repository
#'              and caches it locally for future use.
#' @param omic_type Character string, the omic type (e.g., "proteomics")
#' @param rmd_filename Character string, the template filename (e.g., "DIANN_limpa_report.rmd")
#' @return Character string, path to the downloaded/cached template file
#' @keywords internal
downloadReportTemplate <- function(omic_type, rmd_filename) {
    # Create cache directory using tools (base R, no extra dependencies)
    if (requireNamespace("rappdirs", quietly = TRUE)) {
        cache_base <- rappdirs::user_cache_dir("MultiScholaR")
    } else {
        # Fallback to temp directory if rappdirs not available
        cache_base <- file.path(tempdir(), "MultiScholaR_cache")
    }
    
    cache_dir <- file.path(cache_base, "report_templates", omic_type, "report")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Local cache file path
    cached_file <- file.path(cache_dir, rmd_filename)
    
    # If already cached and less than 7 days old, return it
    if (file.exists(cached_file)) {
        file_age_days <- as.numeric(difftime(Sys.time(), file.info(cached_file)$mtime, units = "days"))
        if (file_age_days < 7) {
            logger::log_info("Using cached template: {cached_file}")
            return(cached_file)
        } else {
            logger::log_info("Cached template is older than 7 days, re-downloading...")
        }
    }
    
    # GitHub URL (raw content from main branch)
    github_url <- sprintf(
        "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/Workbooks/%s/report/%s",
        omic_type, rmd_filename
    )
    
    logger::log_info("Downloading template from GitHub: {github_url}")
    
    # Download
    tryCatch({
        download.file(github_url, cached_file, mode = "wb", quiet = TRUE)
        logger::log_info("Successfully downloaded template to: {cached_file}")
        return(cached_file)
    }, error = function(e) {
        rlang::abort(paste0(
            "Failed to download template '", rmd_filename, "' from GitHub.\n",
            "URL: ", github_url, "\n",
            "Error: ", e$message
        ))
    })
}

##################################################################################################################
#' @export
#' @importFrom rlang abort
#' @importFrom tools file_path_sans_ext
RenderReport <- function(omic_type,
                         experiment_label,
                         rmd_filename = "DIANN_report.rmd",
                         project_dirs_object_name = "project_dirs",
                         output_format = NULL) {

    message("--- DEBUG66: Entering RenderReport ---")
    message(sprintf("   DEBUG66: omic_type = '%s'", omic_type))
    message(sprintf("   DEBUG66: experiment_label = '%s'", experiment_label))
    message(sprintf("   DEBUG66: rmd_filename = '%s'", rmd_filename))
    message(sprintf("   DEBUG66: project_dirs_object_name = '%s'", project_dirs_object_name))

    # --- Validate Inputs ---
    if (missing(omic_type) || !is.character(omic_type) || length(omic_type) != 1 || omic_type == "") {
        rlang::abort("`omic_type` must be a single non-empty character string.")
    }
    if (missing(experiment_label) || !is.character(experiment_label) || length(experiment_label) != 1 || experiment_label == "") {
        rlang::abort("`experiment_label` must be a single non-empty character string.")
    }
    if (!is.character(rmd_filename) || length(rmd_filename) != 1 || rmd_filename == "") {
        rlang::abort("`rmd_filename` must be a single non-empty character string.")
    }
    
    message("   DEBUG66: Input validation passed")

    # --- Retrieve Paths from Global Project Directories Object ---
    # Use the new helper function with automatic fallback
    message("   DEBUG66: Calling getProjectPaths...")
    current_paths <- tryCatch({
        getProjectPaths(
            omic_type = omic_type,
            experiment_label = experiment_label,
            project_dirs_object_name = project_dirs_object_name
        )
    }, error = function(e) {
        rlang::abort(paste0("Failed to get project paths: ", e$message))
    })

    message("   DEBUG66: getProjectPaths completed")
    if (!is.null(current_paths)) {
      message(sprintf("      DEBUG66: current_paths is list: %s", is.list(current_paths)))
      message(sprintf("      DEBUG66: current_paths names: %s", paste(names(current_paths), collapse=", ")))
      if ("base_dir" %in% names(current_paths)) {
        message(sprintf("      DEBUG66: base_dir = %s", current_paths$base_dir))
      }
      if ("results_summary_dir" %in% names(current_paths)) {
        message(sprintf("      DEBUG66: results_summary_dir = %s", current_paths$results_summary_dir))
      }
      if ("source_dir" %in% names(current_paths)) {
        message(sprintf("      DEBUG66: source_dir = %s", current_paths$source_dir))
      }
    } else {
      message("      DEBUG66: current_paths is NULL!")
    }

    if (!is.list(current_paths) || 
        is.null(current_paths$base_dir) || # Need base_dir to find the template Rmd
        is.null(current_paths$results_summary_dir)) {
        rlang::abort(paste0("Essential paths (base_dir, results_summary_dir) missing from project_dirs"))
        rlang::abort(paste0("Essential paths (base_dir, results_summary_dir) missing from project_dirs"))
    }
    
    message("   DEBUG66: Path validation passed")

    # --- Determine the source Rmd template directory (e.g., scripts/proteomics) ---
    # This logic mirrors part of setupDirectories to find the correct unlabelled script source leaf.
    omic_script_template_leaf <- switch(omic_type,
        proteomics = "proteomics",
        metabolomics = "metabolomics",
        transcriptomics = "transcriptomics",
        lipidomics = "lipidomics",
        integration = "integration",
        {
            rlang::abort(paste0("Internal error: Unrecognized omic_type ", sQuote(omic_type), " for Rmd template path construction."))
        }
    )
    
    rmd_template_dir <- file.path(current_paths$base_dir, "scripts", omic_script_template_leaf)
    rmd_input_path <- file.path(rmd_template_dir, rmd_filename)
    
    message(sprintf("   DEBUG66: rmd_input_path = %s", rmd_input_path))
    message(sprintf("   DEBUG66: Template file exists: %s", file.exists(rmd_input_path)))
    
    if (!file.exists(rmd_input_path)) {
        rlang::abort(paste0("R Markdown template file not found at the expected location: ", sQuote(rmd_input_path),
                           ". This should be in the general scripts/<omic_type> directory (e.g., scripts/proteomics)."))
    }

    # --- Construct Output Path (in the labelled results_summary directory) ---
    message("   DEBUG66: Constructing output file path...")
    output_file_basename <- paste0(tools::file_path_sans_ext(rmd_filename), 
                                   "_", omic_type, 
                                   "_", experiment_label)
                                   
    output_ext <- ".docx" # Default
    if (!is.null(output_format)){
        if(output_format == "word_document" || grepl("word", output_format, ignore.case=TRUE)){
            output_ext <- ".docx"
        } else if (output_format == "html_document" || grepl("html", output_format, ignore.case=TRUE)){
            output_ext <- ".html"
        } else if (grepl("pdf", output_format, ignore.case=TRUE)){
            output_ext <- ".pdf"
        }
    }

    output_file_path <- file.path(current_paths$results_summary_dir, paste0(output_file_basename, output_ext))

    message(sprintf("   DEBUG66: output_file_path = %s", output_file_path))
    message(sprintf("   DEBUG66: output_ext = %s", output_ext))

    logger::log_info("Attempting to render report:")
    logger::log_info("- Rmd Source (Template): {rmd_input_path}")
    logger::log_info("- Output File: {output_file_path}")
    logger::log_info("- Params: omic_type=\'{omic_type}\', experiment_label=\'{experiment_label}\'")

    # Read study_parameters.txt to extract workflow_name and timestamp
    message("   DEBUG66: Reading study_parameters.txt...")
    params_path <- file.path(current_paths$source_dir, "study_parameters.txt")
    message(sprintf("      DEBUG66: params_path = %s", params_path))
    message(sprintf("      DEBUG66: params file exists: %s", file.exists(params_path)))
    if (file.exists(params_path)) {
        lines <- readLines(params_path)
        workflow_name_line <- grep("Workflow Name:", lines, value = TRUE)
        timestamp_line <- grep("Timestamp:", lines, value = TRUE)
        workflow_name <- if (length(workflow_name_line) > 0) trimws(sub("Workflow Name:", "", workflow_name_line[1])) else "Unknown Workflow"
        timestamp <- if (length(timestamp_line) > 0) trimws(sub("Timestamp:", "", timestamp_line[1])) else format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    } else {
        workflow_name <- "Unknown Workflow"
        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    }
    
    message(sprintf("   DEBUG66: workflow_name = '%s'", workflow_name))
    message(sprintf("   DEBUG66: timestamp = '%s'", timestamp))

    # --- Render the Report ---
    message("   DEBUG66: About to call rmarkdown::render...")
    message("      DEBUG66: Render parameters:")
    message(sprintf("         input = %s", rmd_input_path))
    message(sprintf("         output_file = %s", output_file_path))
    message(sprintf("         output_format = %s", ifelse(is.null(output_format), "NULL (using Rmd default)", output_format)))
    
    rendered_path <- tryCatch({
        rmarkdown::render(
            input = rmd_input_path,
            params = list(
                omic_type = omic_type,
                experiment_label = experiment_label,
                workflow_name = workflow_name,
                timestamp = timestamp
            ),
            output_file = output_file_path,
            output_format = output_format, # Pass this along; if NULL, Rmd default is used
            envir = new.env(parent = globalenv()) # Render in a clean environment
        )
    }, error = function(e) {
        message("   DEBUG66: rmarkdown::render FAILED")
        message(sprintf("      DEBUG66: Error message: %s", e$message))
        message("      DEBUG66: Error traceback:")
        print(e)
        logger::log_error("Failed to render R Markdown report: {e$message}")
        logger::log_error("Input path: {rmd_input_path}")
        logger::log_error("Output path: {output_file_path}")
        NULL # Return NULL on failure
    })

    message(sprintf("   DEBUG66: rmarkdown::render returned: %s", ifelse(is.null(rendered_path), "NULL", rendered_path)))
    if (!is.null(rendered_path)) {
        message(sprintf("   DEBUG66: Output file exists: %s", file.exists(rendered_path)))
    }

    if (!is.null(rendered_path) && file.exists(rendered_path)) {
        logger::log_info("Report successfully rendered to: {rendered_path}")
    } else {
        logger::log_warn("Report rendering failed or output file not found at expected location.")
    }
    
    message("--- DEBUG66: Exiting RenderReport ---")
    message(sprintf("   DEBUG66: Returning: %s", ifelse(is.null(rendered_path), "NULL", rendered_path)))
    
    invisible(rendered_path)
}

#' @title Update Parameter in S4 Object Args and Global Config List
#' @description Modifies a specific parameter within an S4 object's @args slot
#'              and also updates the corresponding value in a global list named
#'              'config_list'.
#'
#' @param theObject The S4 object whose @args slot needs updating.
#' @param function_name The name identifying the parameter section (character string,
#'                      e.g., "peptideIntensityFiltering"). Corresponds to the
#'                      first-level key in both @args and config_list.
#' @param parameter_name The specific parameter name to update (character string,
#'                       e.g., "peptides_proportion_of_samples_below_cutoff").
#'                       Corresponds to the second-level key.
#' @param new_value The new value to assign to the parameter.
#' @param config_list_name The name of the global list variable holding the
#'                         configuration (defaults to "config_list").
#' @param env The environment where the global config list resides (defaults to
#'            .GlobalEnv).
#'
#' @return The modified S4 object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assume 'myPeptideData' is a PeptideQuantitativeData object
#' # Assume 'config_list' exists in the global environment
#'
#' # Check initial values (example)
#' # print(myPeptideData@args$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff)
#' # print(config_list$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff)
#'
#' # Update the parameter to 0.7
#' myPeptideData <- updateConfigParameter(
#'   theObject = myPeptideData,
#'   function_name = "peptideIntensityFiltering",
#'   parameter_name = "peptides_proportion_of_samples_below_cutoff",
#'   new_value = 0.7
#' )
#'
#' # Verify changes (example)
#' # print(myPeptideData@args$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff) # Should be 0.7
#' # print(config_list$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff) # Should be 0.7
#' }
updateConfigParameter <- function(theObject,
                                function_name,
                                parameter_name,
                                new_value,
                                config_list_name = "config_list",
                                env = .GlobalEnv) {

  # --- Input Validation ---
  if (!isS4(theObject)) {
    stop("'theObject' must be an S4 object.")
  }
  if (!"args" %in% methods::slotNames(theObject)) {
      stop("'theObject' must have an '@args' slot.")
  }
  if (!is.character(function_name) || length(function_name) != 1) {
    stop("'function_name' must be a single character string.")
  }
  if (!is.character(parameter_name) || length(parameter_name) != 1) {
    stop("'parameter_name' must be a single character string.")
  }
  if (!exists(config_list_name, envir = env)) {
      stop("Global config list '", config_list_name, "' not found in the specified environment.")
  }

  # Retrieve the global list safely
  current_config_list <- get(config_list_name, envir = env)

  if (!is.list(current_config_list)) {
      stop("Global variable '", config_list_name, "' is not a list.")
  }
   if (!function_name %in% names(current_config_list)) {
    warning("Function name '", function_name, "' not found in global config list '", config_list_name, "'. Adding it.")
    current_config_list[[function_name]] <- list()
  }
  if (!parameter_name %in% names(current_config_list[[function_name]])) {
       warning("Parameter '", parameter_name, "' not found under '", function_name, "' in global config list '", config_list_name, "'. Adding it.")
  }


  # --- Update S4 Object @args ---
  if (is.null(theObject@args)) {
       theObject@args <- list() # Initialize args if it's NULL
   }
  if (!is.list(theObject@args[[function_name]])) {
     # Initialize the sub-list if it doesn't exist or isn't a list
     theObject@args[[function_name]] <- list()
  }
  theObject@args[[function_name]][[parameter_name]] <- new_value
  message("Updated @args$", function_name, "$", parameter_name, " in S4 object.")


  # --- Update Global Config List ---
  current_config_list[[function_name]][[parameter_name]] <- new_value
  # Assign the modified list back to the global environment
  assign(config_list_name, current_config_list, envir = env)
  message("Updated ", config_list_name, "$", function_name, "$", parameter_name, " in global environment.")

  # --- Return the modified object ---
  return(theObject)
}

#' Choose Best Protein Accession for Data Frame (S3 Version)
#'
#' @description A simplified S3 version of chooseBestProteinAccession that works
#' directly on data frames. Resolves protein groups by selecting the best
#' representative protein based on FASTA sequence data and aggregates
#' quantitative values.
#'
#' @param data_tbl Data frame containing protein quantitative data
#' @param protein_id_column Character string specifying the protein ID column name
#' @param seqinr_obj Data frame with FASTA sequence information (aa_seq_tbl_final)
#' @param seqinr_accession_column Character string specifying the accession column in seqinr_obj
#' @param delim Character string delimiter used to separate protein groups (default: ";")
#' @param replace_zero_with_na Logical, whether to replace zeros with NA (default: TRUE)
#' @param aggregation_method Character string aggregation method ("mean", "median", "sum") (default: "mean")
#'
#' @return Data frame with cleaned protein IDs and aggregated values
#'
#' @details
#' This function processes protein groups (e.g., "P12345;P67890;Q11111") by:
#' \itemize{
#'   \item Splitting groups using the specified delimiter
#'   \item Finding the best representative protein using FASTA sequence data
#'   \item Aggregating quantitative values for the chosen protein
#'   \item Returning cleaned data with single protein IDs
#' }
#'
#' The selection of "best" protein is based on:
#' \itemize{
#'   \item Presence in FASTA sequence database
#'   \item Protein length (longer sequences preferred)
#'   \item Alphabetical order as tiebreaker
#' }
#'
#' @examples
#' \dontrun{
#' cleaned_data <- chooseBestProteinAccession_s3(
#'   data_tbl = protein_data,
#'   protein_id_column = "Protein.Group",
#'   seqinr_obj = aa_seq_tbl_final,
#'   seqinr_accession_column = "uniprot_acc"
#' )
#' }
#'
#' @export
chooseBestProteinAccession_s3 <- function(data_tbl,
                                          protein_id_column,
                                          seqinr_obj,
                                          seqinr_accession_column = "uniprot_acc",
                                          delim = ";",
                                          replace_zero_with_na = TRUE,
                                          aggregation_method = "mean") {
  
  cat("=== Starting chooseBestProteinAccession_s3 ===\n")
  
  tryCatch({
    # Validate inputs with detailed error messages
    cat("*** S3 CLEANUP: Validating inputs ***\n")
    
    if (is.null(data_tbl) || !is.data.frame(data_tbl)) {
      stop("data_tbl must be a non-null data frame")
    }
    
    if (nrow(data_tbl) == 0) {
      stop("data_tbl cannot be empty")
    }
    
    cat(sprintf("*** S3 CLEANUP: data_tbl dimensions: %d rows x %d cols ***\n", nrow(data_tbl), ncol(data_tbl)))
    cat(sprintf("*** S3 CLEANUP: data_tbl column names: %s ***\n", paste(head(names(data_tbl), 10), collapse = ", ")))
    
    if (!protein_id_column %in% names(data_tbl)) {
      stop(sprintf("Protein ID column '%s' not found in data. Available columns: %s", 
                   protein_id_column, paste(names(data_tbl), collapse = ", ")))
    }
    
    if (is.null(seqinr_obj) || !is.data.frame(seqinr_obj)) {
      stop("seqinr_obj must be a non-null data frame")
    }
    
    if (nrow(seqinr_obj) == 0) {
      stop("seqinr_obj cannot be empty")
    }
    
    cat(sprintf("*** S3 CLEANUP: seqinr_obj dimensions: %d rows x %d cols ***\n", nrow(seqinr_obj), ncol(seqinr_obj)))
    cat(sprintf("*** S3 CLEANUP: seqinr_obj column names: %s ***\n", paste(head(names(seqinr_obj), 10), collapse = ", ")))
    
    if (!seqinr_accession_column %in% names(seqinr_obj)) {
      stop(sprintf("Accession column '%s' not found in sequence data. Available columns: %s", 
                   seqinr_accession_column, paste(names(seqinr_obj), collapse = ", ")))
    }
    
    cat("*** S3 CLEANUP: Input validation passed ***\n")
    
    # Get non-quantitative columns (metadata columns) with error handling
    cat("*** S3 CLEANUP: Identifying column types ***\n")
    
    numeric_cols <- tryCatch({
      sapply(data_tbl, is.numeric)
    }, error = function(e) {
      stop(sprintf("Error checking column types: %s", e$message))
    })
    
    quant_cols <- names(data_tbl)[numeric_cols]
    meta_cols <- names(data_tbl)[!numeric_cols]
    
    cat(sprintf("*** S3 CLEANUP: Found %d quantitative columns and %d metadata columns ***\n", 
                length(quant_cols), length(meta_cols)))
    cat(sprintf("*** S3 CLEANUP: Quantitative columns: %s ***\n", paste(head(quant_cols, 5), collapse = ", ")))
    cat(sprintf("*** S3 CLEANUP: Metadata columns: %s ***\n", paste(head(meta_cols, 5), collapse = ", ")))
  
    # Replace zeros with NA if requested
    if (replace_zero_with_na && length(quant_cols) > 0) {
      cat("*** S3 CLEANUP: Replacing zeros with NA in quantitative columns ***\n")
      
      tryCatch({
        data_tbl[quant_cols] <- lapply(data_tbl[quant_cols], function(x) {
          if (is.numeric(x)) {
            x[x == 0] <- NA
          }
          x
        })
        cat("*** S3 CLEANUP: Zero replacement completed ***\n")
      }, error = function(e) {
        stop(sprintf("Error replacing zeros with NA: %s", e$message))
      })
    }
    
    # Get unique protein groups
    cat("*** S3 CLEANUP: Getting unique protein groups ***\n")
    
    protein_groups <- tryCatch({
      unique(data_tbl[[protein_id_column]])
    }, error = function(e) {
      stop(sprintf("Error getting unique protein groups: %s", e$message))
    })
    
    cat(sprintf("*** S3 CLEANUP: Processing %d unique protein groups ***\n", length(protein_groups)))
    cat(sprintf("*** S3 CLEANUP: First 5 protein groups: %s ***\n", paste(head(protein_groups, 5), collapse = ", ")))
    
    # Create lookup table for sequence information
    cat("*** S3 CLEANUP: Creating sequence lookup table ***\n")
    
    seq_lookup <- tryCatch({
      if ("length" %in% names(seqinr_obj)) {
        cat("*** S3 CLEANUP: Using existing length column ***\n")
        seqinr_obj |>
          dplyr::select(dplyr::all_of(c(seqinr_accession_column, "length"))) |>
          dplyr::distinct()
      } else {
        cat("*** S3 CLEANUP: Creating default length column ***\n")
        # If no length column, create one or use a default
        seqinr_obj |>
          dplyr::select(dplyr::all_of(seqinr_accession_column)) |>
          dplyr::distinct() |>
          dplyr::mutate(length = 1000)  # Default length
      }
    }, error = function(e) {
      stop(sprintf("Error creating sequence lookup: %s", e$message))
    })
    
    cat(sprintf("*** S3 CLEANUP: Created sequence lookup with %d entries ***\n", nrow(seq_lookup)))
  
    # Function to choose best protein from a group
    cat("*** S3 CLEANUP: Defining chooseBestFromGroup function ***\n")
    
    chooseBestFromGroup <- function(group_string) {
      tryCatch({
        if (is.na(group_string) || group_string == "") {
          return(NA_character_)
        }
        
        # Split the group
        proteins <- stringr::str_split(group_string, delim)[[1]]
        proteins <- stringr::str_trim(proteins)  # Remove whitespace
        proteins <- proteins[proteins != ""]     # Remove empty strings
        
        if (length(proteins) == 1) {
          return(proteins[1])
        }
        
        # Find proteins that exist in sequence database
        proteins_in_db <- proteins[proteins %in% seq_lookup[[seqinr_accession_column]]]
        
        if (length(proteins_in_db) == 0) {
          # No proteins found in database, return first one
          return(proteins[1])
        }
        
        if (length(proteins_in_db) == 1) {
          return(proteins_in_db[1])
        }
        
        # Multiple proteins in database - choose by length (longest first)
        protein_info <- seq_lookup |>
          dplyr::filter(!!sym(seqinr_accession_column) %in% proteins_in_db) |>
          dplyr::arrange(desc(length), !!sym(seqinr_accession_column))
        
        return(protein_info[[seqinr_accession_column]][1])
        
      }, error = function(e) {
        warning(sprintf("Error processing protein group '%s': %s", group_string, e$message))
        return(group_string)  # Return original if processing fails
      })
    }
    
    # Create mapping of original groups to best proteins
    cat("*** S3 CLEANUP: Creating protein group mapping ***\n")
    
    protein_mapping <- tryCatch({
      data.frame(
        original_group = protein_groups,
        best_protein = sapply(protein_groups, chooseBestFromGroup, USE.NAMES = FALSE),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      stop(sprintf("Error creating protein mapping: %s", e$message))
    })
    
    # Count reductions
    groups_with_multiple <- tryCatch({
      sum(stringr::str_detect(protein_mapping$original_group, delim), na.rm = TRUE)
    }, error = function(e) {
      cat(sprintf("*** S3 CLEANUP: Warning - could not count multi-protein groups: %s ***\n", e$message))
      0
    })
    
    cat(sprintf("*** S3 CLEANUP: Found %d protein groups with multiple proteins ***\n", groups_with_multiple))
  
    # Add mapping to data
    cat("*** S3 CLEANUP: Applying protein mapping to data ***\n")
    
    data_with_mapping <- tryCatch({
      data_tbl |>
        dplyr::left_join(
          protein_mapping,
          by = setNames("original_group", protein_id_column)
        )
    }, error = function(e) {
      stop(sprintf("Error applying protein mapping: %s", e$message))
    })
    
    cat(sprintf("*** S3 CLEANUP: Data mapping completed. Mapped data has %d rows ***\n", nrow(data_with_mapping)))
    
    # Group by best protein and aggregate
    cat(sprintf("*** S3 CLEANUP: Aggregating data using method: %s ***\n", aggregation_method))
    
    # Define aggregation function
    agg_func <- switch(aggregation_method,
      "mean" = function(x) mean(x, na.rm = TRUE),
      "median" = function(x) median(x, na.rm = TRUE),
      "sum" = function(x) sum(x, na.rm = TRUE),
      function(x) mean(x, na.rm = TRUE)  # Default to mean
    )
    
    # Separate quantitative and metadata columns for different aggregation
    cat("*** S3 CLEANUP: Aggregating quantitative data ***\n")
    
    quant_data <- tryCatch({
      if (length(quant_cols) > 0) {
        data_with_mapping |>
          dplyr::select(best_protein, dplyr::all_of(quant_cols)) |>
          dplyr::group_by(best_protein) |>
          dplyr::summarise(
            dplyr::across(dplyr::all_of(quant_cols), agg_func),
            .groups = "drop"
          )
      } else {
        data.frame(best_protein = unique(data_with_mapping$best_protein))
      }
    }, error = function(e) {
      stop(sprintf("Error aggregating quantitative data: %s", e$message))
    })
    
    cat("*** S3 CLEANUP: Aggregating metadata ***\n")
    
    # For metadata columns, take the first non-NA value for each best protein
    meta_data <- tryCatch({
      if (length(meta_cols) > 0) {
        data_with_mapping |>
          dplyr::select(best_protein, dplyr::all_of(meta_cols)) |>
          dplyr::group_by(best_protein) |>
          dplyr::summarise(
            dplyr::across(dplyr::all_of(meta_cols), ~ dplyr::first(na.omit(.))[1]),
            .groups = "drop"
          )
      } else {
        data.frame(best_protein = unique(data_with_mapping$best_protein))
      }
    }, error = function(e) {
      stop(sprintf("Error aggregating metadata: %s", e$message))
    })
    
    cat("*** S3 CLEANUP: Combining results ***\n")
    
    # Combine results
    result_data <- tryCatch({
      quant_data |>
        dplyr::left_join(meta_data, by = "best_protein") |>
        dplyr::rename(!!sym(protein_id_column) := best_protein)
    }, error = function(e) {
      stop(sprintf("Error combining results: %s", e$message))
    })
    
    cat("*** S3 CLEANUP: Reordering columns ***\n")
    
    # Reorder columns to match original
    result_data <- tryCatch({
      result_data |>
        dplyr::select(dplyr::all_of(names(data_tbl)))
    }, error = function(e) {
      stop(sprintf("Error reordering columns: %s", e$message))
    })
    
    # Final statistics
    original_proteins <- length(unique(data_tbl[[protein_id_column]]))
    final_proteins <- length(unique(result_data[[protein_id_column]]))
    reduction_pct <- round(((original_proteins - final_proteins) / original_proteins) * 100, 1)
    
    cat(sprintf("*** S3 CLEANUP: Protein cleanup completed: %d -> %d proteins (%.1f%% reduction) ***\n",
                     original_proteins, final_proteins, reduction_pct))
    
    return(result_data)
    
  }, error = function(e) {
    cat(sprintf("*** S3 CLEANUP FATAL ERROR: %s ***\n", e$message))
    cat("*** S3 CLEANUP: Returning original data ***\n")
    return(data_tbl)  # Return original data if everything fails
  })
}

#' Create Study Parameters File
#'
#' @description
#' Creates a study parameters text file directly without using S4 objects.
#' This replaces the overly complex createWorkflowArgsFromConfig + WorkflowArgs show() approach.
#'
#' @param workflow_name Character string, name of the workflow
#' @param description Character string, description of the analysis
#' @param organism_name Character string, organism name (optional)
#' @param taxon_id Character string or numeric, taxon ID (optional)
#' @param source_dir_path Character string, path to save the study_parameters.txt file
#' @param contrasts_tbl Data frame, contrasts table (optional)
#' @param config_list_name Character string, name of the global config list (default: "config_list")
#' @param env Environment where the config list resides (default: .GlobalEnv)
#'
#' @return Character string path to the created study_parameters.txt file
#'
#' @examples
#' \dontrun{
#' file_path <- createStudyParametersFile(
#'   workflow_name = "proteomics_analysis",
#'   description = "DIA proteomics analysis",
#'   source_dir_path = "/path/to/scripts"
#' )
#' }
#'
#' @export
createStudyParametersFile <- function(workflow_name, 
                                      description = "", 
                                      organism_name = NULL, 
                                      taxon_id = NULL,
                                      source_dir_path = NULL,
                                      contrasts_tbl = NULL,
                                      config_list_name = "config_list",
                                      env = .GlobalEnv) {
    
    # Validate required inputs
    if (missing(workflow_name) || !is.character(workflow_name) || length(workflow_name) != 1) {
        stop("workflow_name must be a single character string")
    }
    
    if (is.null(source_dir_path) || !is.character(source_dir_path) || length(source_dir_path) != 1) {
        stop("source_dir_path must be a single character string")
    }
    
    if (!dir.exists(source_dir_path)) {
        stop("source_dir_path directory does not exist: ", source_dir_path)
    }
    
    # Get git information
    git_info <- tryCatch({
        if (requireNamespace("gh", quietly = TRUE)) {
            # Make the API call (works for public repos without token)
            branch_info <- gh::gh("/repos/APAF-BIOINFORMATICS/MultiScholaR/branches/main")
            list(
                commit_sha = branch_info$commit$sha,
                branch = "main",
                repo = "MultiScholaR", 
                timestamp = branch_info$commit$commit$author$date
            )
        } else {
            list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
        }
    }, error = function(e) {
        message("Error fetching GitHub info: ", e$message)
        list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
    })
    
    # Get organism information from session if not provided
    if (is.null(organism_name) && exists("organism_name", envir = .GlobalEnv)) {
        organism_name <- get("organism_name", envir = .GlobalEnv)
    }
    
    if (is.null(taxon_id) && exists("taxon_id", envir = .GlobalEnv)) {
        taxon_id <- get("taxon_id", envir = .GlobalEnv)
    }
    
    # Get config list
    config_list <- if (exists(config_list_name, envir = env)) {
        get(config_list_name, envir = env)
    } else {
        list()
    }
    
    # Get contrasts_tbl if not provided
    if (is.null(contrasts_tbl)) {
        if (exists("contrasts_tbl", envir = parent.frame())) {
            contrasts_tbl <- get("contrasts_tbl", envir = parent.frame())
        } else if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
        }
    }
    
    # Build the output content
    output_lines <- c(
        "Study Parameters",
        "================",
        "",
        "Basic Information:",
        "-----------------",
        paste("Workflow Name:", workflow_name),
        paste("Description:", if(nzchar(description)) description else "N/A"),
        paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        ""
    )
    
    # Add git information
    if (!is.null(git_info) && is.list(git_info)) {
        git_lines <- c(
            "Git Information:",
            "---------------",
            paste("Repository:", ifelse(!is.null(git_info$repo) && !is.na(git_info$repo), git_info$repo, "N/A")),
            paste("Branch:", ifelse(!is.null(git_info$branch) && !is.na(git_info$branch), git_info$branch, "N/A")),
            paste("Commit:", ifelse(!is.null(git_info$commit_sha) && !is.na(git_info$commit_sha), substr(git_info$commit_sha, 1, 7), "N/A")),
            paste("Git Timestamp:", ifelse(!is.null(git_info$timestamp) && !is.na(git_info$timestamp), git_info$timestamp, "N/A")),
            ""
        )
        output_lines <- c(output_lines, git_lines)
    }
    
    # Add organism information
    if (!is.null(organism_name) && !is.na(organism_name) && nzchar(organism_name)) {
        organism_lines <- c(
            "Organism Information:",
            "---------------------",
            paste("Organism Name:", organism_name)
        )
        if (!is.null(taxon_id) && !is.na(taxon_id) && nzchar(taxon_id)) {
            organism_lines <- c(organism_lines, paste("Taxon ID:", taxon_id))
        }
        organism_lines <- c(organism_lines, "")
        output_lines <- c(output_lines, organism_lines)
    }
    
    # Add configuration parameters
    config_lines <- c(
        "Configuration Parameters:",
        "-------------------------"
    )
    
    # Clean the config list (remove problematic objects)
    clean_config <- config_list
    # Remove the internal source dir if it exists
    if (!is.null(clean_config$internal_workflow_source_dir)) {
        clean_config$internal_workflow_source_dir <- NULL
    }
    
    # Format the config list
    if (length(clean_config) > 0) {
        config_params <- formatConfigList(clean_config)
        config_lines <- c(config_lines, config_params)
    } else {
        config_lines <- c(config_lines, "No configuration parameters available")
    }
    
    output_lines <- c(output_lines, config_lines)
    
    # Add contrasts information
    if (!is.null(contrasts_tbl) && (is.data.frame(contrasts_tbl) || tibble::is_tibble(contrasts_tbl)) && nrow(contrasts_tbl) > 0) {
        contrasts_lines <- c(
            "",
            "Contrasts:",
            "----------"
        )
        
        if ("contrasts" %in% colnames(contrasts_tbl)) {
            contrasts_info <- tryCatch({
                contrasts_col <- contrasts_tbl[["contrasts"]]
                paste("  ", as.character(contrasts_col))
            }, error = function(e) {
                paste("  [Error extracting contrasts:", e$message, "]")
            })
            contrasts_lines <- c(contrasts_lines, contrasts_info)
        } else {
            contrasts_lines <- c(contrasts_lines, "  [Column 'contrasts' not found in contrasts_tbl]")
        }
        
        output_lines <- c(output_lines, contrasts_lines)
    }
    
    # Write to file
    output_file <- file.path(source_dir_path, "study_parameters.txt")
    
    tryCatch({
        writeLines(output_lines, output_file)
        message("Study parameters saved to: ", output_file)
        return(output_file)
    }, error = function(e) {
        stop("Failed to write study parameters file: ", e$message)
    })
}

# Updated function to extract parameters from S4 @args slot
#' @title Create workflow arguments from config
#' @param workflow_name Character string, name of the workflow
#' @param description Character string, description of the analysis
#' @param organism_name Character string, organism name (optional)
#' @param taxon_id Character string or numeric, taxon ID (optional)
#' @param source_dir_path Character string, path to save the study_parameters.txt file
#' @param final_s4_object S4 object containing workflow parameters in @args slot (optional)
#' @param contrasts_tbl Data frame, contrasts table (optional)
#' @param workflow_data List containing workflow state and optimization results (optional)
#' @export
createWorkflowArgsFromConfig <- function(workflow_name, description = "", 
                                        organism_name = NULL, taxon_id = NULL,
                                        source_dir_path = NULL, 
                                        final_s4_object = NULL,
                                        contrasts_tbl = NULL,
                                        workflow_data = NULL) {
    # Validate required inputs
    if (missing(workflow_name) || !is.character(workflow_name) || length(workflow_name) != 1) {
        stop("workflow_name must be a single character string")
    }
    
    if (is.null(source_dir_path) || !is.character(source_dir_path) || length(source_dir_path) != 1) {
        stop("source_dir_path must be a single character string")
    }
    
    if (!dir.exists(source_dir_path)) {
        stop("source_dir_path directory does not exist: ", source_dir_path)
    }
    
    cat("WORKFLOW ARGS: Starting parameter extraction\n")
    
    # Extract parameters from S4 object @args if provided
    s4_params <- list()
         # Check if final_s4_object has @args slot
     s4_has_args <- tryCatch({
       !is.null(final_s4_object) && isS4(final_s4_object) && !is.null(final_s4_object@args)
     }, error = function(e) {
       FALSE
     })
     
     if (s4_has_args) {
        cat("WORKFLOW ARGS: Extracting parameters from S4 @args slot\n")
        
        tryCatch({
            s4_params <- final_s4_object@args
            cat(sprintf("WORKFLOW ARGS: Found %d function groups in S4 @args\n", length(s4_params)))
            
            # Log the function groups for debugging
            for (func_name in names(s4_params)) {
                param_count <- length(s4_params[[func_name]])
                cat(sprintf("WORKFLOW ARGS: Function '%s' has %d parameters\n", func_name, param_count))
            }
        }, error = function(e) {
            cat(sprintf("WORKFLOW ARGS: Error extracting S4 params: %s\n", e$message))
            s4_params <- list()
        })
    } else {
        cat("WORKFLOW ARGS: No S4 object provided or S4 @args is NULL\n")
    }
    
    # Get fallback config_list from global environment
    config_list <- if (exists("config_list", envir = .GlobalEnv)) {
        get("config_list", envir = .GlobalEnv)
    } else {
        list()
    }
    
         cat(sprintf("WORKFLOW ARGS: Global config_list has %d sections\n", length(config_list)))
     
     # Extract RUV optimization results from workflow_data if available
     ruv_optimization_result <- NULL
     if (!is.null(workflow_data)) {
         cat("WORKFLOW ARGS: Checking for RUV optimization results in workflow_data\n")
         
         # Check multiple possible locations for RUV optimization results
         if (!is.null(workflow_data$ruv_optimization_result)) {
             ruv_optimization_result <- workflow_data$ruv_optimization_result
             cat("WORKFLOW ARGS: Found RUV optimization results in workflow_data$ruv_optimization_result\n")
         } else if (!is.null(workflow_data$state_manager)) {
             # Try to get RUV results from state manager config
             tryCatch({
                 current_state_config <- workflow_data$state_manager$getStateConfig(workflow_data$state_manager$current_state)
                 if (!is.null(current_state_config$ruv_optimization_result)) {
                     ruv_optimization_result <- current_state_config$ruv_optimization_result
                     cat("WORKFLOW ARGS: Found RUV optimization results in state manager config\n")
                 }
             }, error = function(e) {
                 cat(sprintf("WORKFLOW ARGS: Could not extract RUV results from state manager: %s\n", e$message))
             })
         }
         
         # Also check if there's a saved RUV file
         if (is.null(ruv_optimization_result) && !is.null(source_dir_path)) {
             ruv_file <- file.path(source_dir_path, "ruv_optimization_results.RDS")
             if (file.exists(ruv_file)) {
                 tryCatch({
                     loaded_result <- readRDS(ruv_file)
                     # Check if loaded result indicates RUV was skipped
                     if (isTRUE(loaded_result$ruv_skipped)) {
                         cat(sprintf("WORKFLOW ARGS: Loaded RUV skip result from file: %s\n", ruv_file))
                         ruv_optimization_result <- loaded_result
                     } else {
                         # Only use file if it's NOT a skip result (backwards compatibility)
                         ruv_optimization_result <- loaded_result
                         cat(sprintf("WORKFLOW ARGS: Loaded RUV optimization results from file: %s\n", ruv_file))
                     }
                 }, error = function(e) {
                     cat(sprintf("WORKFLOW ARGS: Could not load RUV file: %s\n", e$message))
                 })
             }
         }
     } 
     # Merge S4 parameters with config_list (S4 takes precedence)
     merged_config <- config_list
     
     if (length(s4_params) > 0) {
         cat("WORKFLOW ARGS: Merging S4 parameters with config_list\n")
         
         # S4 parameters override config_list parameters
         for (func_name in names(s4_params)) {
             if (is.list(s4_params[[func_name]]) && length(s4_params[[func_name]]) > 0) {
                 merged_config[[func_name]] <- s4_params[[func_name]]
                 cat(sprintf("WORKFLOW ARGS: Updated '%s' section with %d S4 parameters\n", 
                            func_name, length(s4_params[[func_name]])))
             }
         }
     }
    
    # Get git information
    git_info <- tryCatch({
        if (requireNamespace("gh", quietly = TRUE)) {
            # Make the API call (works for public repos without token)
            branch_info <- gh::gh("/repos/APAF-BIOINFORMATICS/MultiScholaR/branches/main")
            list(
                commit_sha = branch_info$commit$sha,
                branch = "main",
                repo = "MultiScholaR", 
                timestamp = branch_info$commit$commit$author$date
            )
        } else {
            list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
        }
    }, error = function(e) {
        message("Error fetching GitHub info: ", e$message)
        list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
    })
    
    # Get organism information from session if not provided
    if (is.null(organism_name) && exists("organism_name", envir = .GlobalEnv)) {
        organism_name <- get("organism_name", envir = .GlobalEnv)
    }
    
    if (is.null(taxon_id) && exists("taxon_id", envir = .GlobalEnv)) {
        taxon_id <- get("taxon_id", envir = .GlobalEnv)
    }
    
    # Get contrasts_tbl if not provided
    if (is.null(contrasts_tbl)) {
        if (exists("contrasts_tbl", envir = parent.frame())) {
            contrasts_tbl <- get("contrasts_tbl", envir = parent.frame())
        } else if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
        }
    }
    
    # Build the output content
    output_lines <- c(
        "Study Parameters",
        "================",
        "",
        "Basic Information:",
        "-----------------",
        paste("Workflow Name:", workflow_name),
        paste("Description:", if(nzchar(description)) description else "N/A"),
        paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        ""
    )
    
    # Add git information
    if (!is.null(git_info) && is.list(git_info)) {
        git_lines <- c(
            "Git Information:",
            "---------------",
            paste("Repository:", ifelse(!is.null(git_info$repo) && !is.na(git_info$repo), git_info$repo, "N/A")),
            paste("Branch:", ifelse(!is.null(git_info$branch) && !is.na(git_info$branch), git_info$branch, "N/A")),
            paste("Commit:", ifelse(!is.null(git_info$commit_sha) && !is.na(git_info$commit_sha), substr(git_info$commit_sha, 1, 7), "N/A")),
            paste("Git Timestamp:", ifelse(!is.null(git_info$timestamp) && !is.na(git_info$timestamp), git_info$timestamp, "N/A")),
            ""
        )
        output_lines <- c(output_lines, git_lines)
    }
    
         # Add organism information
     if (!is.null(organism_name) && !is.na(organism_name) && nzchar(organism_name)) {
         organism_lines <- c(
             "Organism Information:",
             "---------------------",
             paste("Organism Name:", organism_name)
         )
         if (!is.null(taxon_id) && !is.na(taxon_id) && nzchar(taxon_id)) {
             organism_lines <- c(organism_lines, paste("Taxon ID:", taxon_id))
         }
         organism_lines <- c(organism_lines, "")
         output_lines <- c(output_lines, organism_lines)
     }
     
    # Add RUV optimization results if available
    if (!is.null(ruv_optimization_result) && is.list(ruv_optimization_result)) {
        
        # Check if RUV was skipped
        ruv_was_skipped <- isTRUE(ruv_optimization_result$ruv_skipped)
        
        if (ruv_was_skipped) {
            cat("WORKFLOW ARGS: RUV was skipped - writing skip section\n")
            ruv_lines <- c(
                "RUV-III Batch Correction:",
                "-------------------------",
                "• Status: Not Applied",
                "• Reason: User determined RUV was not appropriate due to dataset constraints",
                ""
            )
            output_lines <- c(output_lines, ruv_lines)
        } else {
            cat("WORKFLOW ARGS: Formatting RUV optimization results\n")
            
            ruv_lines <- c(
                "Automatic RUV Optimization Results:",
                "-----------------------------------"
            )
         
         # Extract and format RUV optimization values
         tryCatch({
             best_percentage <- if (!is.null(ruv_optimization_result$best_percentage)) {
                 sprintf("%.1f%%", ruv_optimization_result$best_percentage)
             } else { "N/A" }
             
             best_k <- if (!is.null(ruv_optimization_result$best_k)) {
                 as.character(ruv_optimization_result$best_k)
             } else { "N/A" }
             
             separation_score <- if (!is.null(ruv_optimization_result$best_separation_score)) {
                 sprintf("%.4f", ruv_optimization_result$best_separation_score)
             } else { "N/A" }
             
             composite_score <- if (!is.null(ruv_optimization_result$best_composite_score)) {
                 sprintf("%.4f", ruv_optimization_result$best_composite_score)
             } else { "N/A" }
             
             control_genes_count <- if (!is.null(ruv_optimization_result$best_control_genes_index)) {
                 as.character(sum(ruv_optimization_result$best_control_genes_index, na.rm = TRUE))
             } else { "N/A" }
             
             separation_metric <- if (!is.null(ruv_optimization_result$separation_metric_used)) {
                 as.character(ruv_optimization_result$separation_metric_used)
             } else { "N/A" }
             
             k_penalty_weight <- if (!is.null(ruv_optimization_result$k_penalty_weight)) {
                 sprintf("%.1f", ruv_optimization_result$k_penalty_weight)
             } else { "N/A" }
             
             adaptive_penalty <- if (!is.null(ruv_optimization_result$adaptive_k_penalty_used)) {
                 ifelse(ruv_optimization_result$adaptive_k_penalty_used, "TRUE", "FALSE")
             } else { "N/A" }
             
             sample_size <- if (!is.null(ruv_optimization_result$sample_size)) {
                 as.character(ruv_optimization_result$sample_size)
             } else { "N/A" }
             
             # Also try to get RUV grouping variable from S4 parameters
             ruv_grouping_variable <- "N/A"
             if (!is.null(s4_params$ruvIII_C_Varying$ruv_grouping_variable)) {
                 ruv_grouping_variable <- s4_params$ruvIII_C_Varying$ruv_grouping_variable
             } else if (!is.null(s4_params$getNegCtrlProtAnova$ruv_grouping_variable)) {
                 ruv_grouping_variable <- s4_params$getNegCtrlProtAnova$ruv_grouping_variable
             }
             
             ruv_lines <- c(ruv_lines,
                 paste("• Best percentage:", best_percentage),
                 paste("• Best k value:", best_k),
                 paste("• Separation score:", separation_score),
                 paste("• Composite score:", composite_score),
                 paste("• Control genes:", control_genes_count),
                 paste("• RUV grouping variable:", ruv_grouping_variable),
                 paste("• Separation metric:", separation_metric),
                 paste("• K penalty weight:", k_penalty_weight),
                 paste("• Adaptive penalty:", adaptive_penalty),
                 paste("• Sample size:", sample_size),
                 ""
             )
             
             cat("WORKFLOW ARGS: Successfully formatted RUV optimization results\n")
             
        }, error = function(e) {
            cat(sprintf("WORKFLOW ARGS: Error formatting RUV results: %s\n", e$message))
            ruv_lines <- c(ruv_lines,
                paste("• [Error formatting RUV optimization results:", e$message, "]"),
                ""
            )
        })
        
        output_lines <- c(output_lines, ruv_lines)
        } # End of else block for RUV applied
   } else {
       cat("WORKFLOW ARGS: No RUV optimization results available\n")
   }
    
    # Add FASTA Processing Information
    if (!is.null(workflow_data) && !is.null(workflow_data$fasta_metadata)) {
        cat("WORKFLOW ARGS: Adding FASTA processing information\n")
        fasta_lines <- c(
            "FASTA File Processing:",
            "---------------------",
            paste("• Format:", workflow_data$fasta_metadata$fasta_format),
            paste("• Sequences:", workflow_data$fasta_metadata$num_sequences),
            paste("• Protein Evidence Available:", workflow_data$fasta_metadata$has_protein_evidence),
            paste("• Gene Names Available:", workflow_data$fasta_metadata$has_gene_names),
            paste("• Isoform Information Available:", workflow_data$fasta_metadata$has_isoform_info),
            paste("• Status Information Available:", workflow_data$fasta_metadata$has_status_info),
            ""
        )
        output_lines <- c(output_lines, fasta_lines)
    } else {
        cat("WORKFLOW ARGS: No FASTA metadata available\n")
    }
    
    # ✅ NEW: Add Mixed Species FASTA Analysis section
    if (!is.null(workflow_data) && !is.null(workflow_data$mixed_species_analysis)) {
        cat("WORKFLOW ARGS: Adding mixed species analysis information\n")
        mixed_species_info <- workflow_data$mixed_species_analysis
        
        mixed_species_lines <- c(
            "Mixed Species FASTA Analysis:",
            "-----------------------------",
            paste("• Multi-Species FASTA Used:", ifelse(isTRUE(mixed_species_info$enabled), "Yes", "No"))
        )
        
        if (isTRUE(mixed_species_info$enabled)) {
            mixed_species_lines <- c(mixed_species_lines,
                paste("• Selected Primary Organism:", mixed_species_info$selected_organism %||% "N/A"),
                paste("• Selected Taxon ID:", mixed_species_info$selected_taxon_id %||% "N/A"),
                paste("• Filtered at Import:", ifelse(isTRUE(mixed_species_info$filter_applied_at_import), "Yes", "No"))
            )
            
            # Add organism distribution summary if available
            if (!is.null(mixed_species_info$organism_distribution) && 
                is.data.frame(mixed_species_info$organism_distribution) &&
                nrow(mixed_species_info$organism_distribution) > 0) {
                mixed_species_lines <- c(mixed_species_lines, 
                    "",
                    "  Organism Distribution in FASTA:"
                )
                
                # Add top organisms (limit to top 5)
                top_orgs <- utils::head(mixed_species_info$organism_distribution, 5)
                for (i in seq_len(nrow(top_orgs))) {
                    org_line <- sprintf("    - %s (Taxon %s): %d proteins (%.1f%%)",
                        top_orgs$organism_name[i] %||% "Unknown",
                        top_orgs$taxon_id[i] %||% "N/A",
                        top_orgs$protein_count[i] %||% 0,
                        top_orgs$percentage[i] %||% 0
                    )
                    mixed_species_lines <- c(mixed_species_lines, org_line)
                }
            }
        }
        
        mixed_species_lines <- c(mixed_species_lines, "")
        output_lines <- c(output_lines, mixed_species_lines)
    }
    
    # ✅ NEW: Add Enrichment Organism Filtering section
    if (!is.null(workflow_data) && !is.null(workflow_data$enrichment_organism_filter)) {
        cat("WORKFLOW ARGS: Adding enrichment organism filtering information\n")
        filter_info <- workflow_data$enrichment_organism_filter
        
        filter_lines <- c(
            "Enrichment Analysis - Organism Filtering:",
            "-----------------------------------------",
            paste("• Organism Filter Enabled:", ifelse(isTRUE(filter_info$enabled), "Yes", "No"))
        )
        
        if (isTRUE(filter_info$filter_applied)) {
            filter_lines <- c(filter_lines,
                paste("• Filter Applied:", "Yes"),
                paste("• Target Taxon ID:", filter_info$target_taxon_id %||% "N/A"),
                paste("• Proteins Before Filtering:", filter_info$proteins_before %||% "N/A"),
                paste("• Proteins After Filtering:", filter_info$proteins_after %||% "N/A"),
                paste("• Proteins Removed:", filter_info$proteins_removed %||% "N/A")
            )
            
            if (!is.null(filter_info$proteins_before) && filter_info$proteins_before > 0) {
                retention_pct <- round((filter_info$proteins_after / filter_info$proteins_before) * 100, 1)
                filter_lines <- c(filter_lines, 
                    paste("• Retention Rate:", paste0(retention_pct, "%"))
                )
            }
        }
        
        filter_lines <- c(filter_lines, "")
        output_lines <- c(output_lines, filter_lines)
    }
    
    # Add Accession Cleanup Results
    if (!is.null(workflow_data) && !is.null(workflow_data$accession_cleanup_results)) {
        cat("WORKFLOW ARGS: Adding accession cleanup results\n")
        cleanup_lines <- c(
            "Protein Accession Cleanup:",
            "-------------------------",
            paste("• Cleanup Applied:", workflow_data$accession_cleanup_results$cleanup_applied),
            paste("• Aggregation Method:", workflow_data$accession_cleanup_results$aggregation_method),
            paste("• Delimiter Used:", workflow_data$accession_cleanup_results$delimiter_used),
            paste("• Proteins Before Cleanup:", workflow_data$accession_cleanup_results$proteins_before),
            paste("• Proteins After Cleanup:", workflow_data$accession_cleanup_results$proteins_after),
            paste("• Full UniProt Metadata:", workflow_data$accession_cleanup_results$had_full_metadata),
            ""
        )
        output_lines <- c(output_lines, cleanup_lines)
    } else {
        cat("WORKFLOW ARGS: No accession cleanup results available\n")
    }
    
    # Add Protein Filtering Summary
    if (!is.null(workflow_data) && !is.null(workflow_data$protein_counts)) {
        cat("WORKFLOW ARGS: Adding protein filtering summary\n")
        counts_lines <- c(
            "Protein Filtering Summary:",
            "-------------------------",
            paste("• Proteins after QC filtering:", workflow_data$protein_counts$after_qc_filtering %||% "N/A"),
            paste("• Proteins after RUV filtering:", workflow_data$protein_counts$after_ruv_filtering %||% "N/A"),
            paste("• Final proteins for DE analysis:", workflow_data$protein_counts$final_for_de %||% "N/A"),
            ""
        )
        output_lines <- c(output_lines, counts_lines)
    } else {
        cat("WORKFLOW ARGS: No protein counts available\n")
    }
    
    # Add configuration parameters (now includes S4 parameters)
    config_lines <- c(
        "Workflow Parameters:",
        "-------------------",
        ""
    )
    
    # Clean the merged config (remove problematic objects)
    clean_config <- merged_config
    # Remove the internal source dir if it exists
    if (!is.null(clean_config$internal_workflow_source_dir)) {
        clean_config$internal_workflow_source_dir <- NULL
    }
    
    # Add special section for S4-derived parameters  
    if (length(s4_params) > 0) {
        cat("WORKFLOW ARGS: About to format S4 parameters using functional approach\n")
        
        config_lines <- c(config_lines, 
                         "  Parameters from Final S4 Object:",
                         "  --------------------------------")
        
        # SAFE parameter formatting function
        formatParameterValue <- function(param_value) {
            tryCatch({
                if (is.null(param_value)) {
                    "NULL"
                } else if (is.data.frame(param_value)) {
                    # Skip data frames (like seqinr_obj) - too large to serialize
                    sprintf("[Data frame: %d rows x %d cols - omitted for brevity]", 
                            nrow(param_value), ncol(param_value))
                } else if (is.logical(param_value)) {
                    if (length(param_value) == 1) {
                        ifelse(param_value, "TRUE", "FALSE")
                    } else if (length(param_value) > 50) {
                        # Handle large logical vectors (like control genes index)
                        true_count <- sum(param_value, na.rm = TRUE)
                        total_count <- length(param_value)
                        sprintf("logical vector [%d TRUE, %d FALSE out of %d total]", 
                               true_count, total_count - true_count, total_count)
                    } else {
                        # Show first few values for smaller vectors
                        preview <- ifelse(utils::head(param_value, 5), "TRUE", "FALSE")
                        if (length(param_value) > 5) {
                            paste0("c(", paste(preview, collapse = ", "), ", ...)")
                        } else {
                            paste0("c(", paste(preview, collapse = ", "), ")")
                        }
                    }
                } else if (is.numeric(param_value)) {
                    if (length(param_value) == 1) {
                        as.character(param_value)
                    } else if (length(param_value) > 5) {
                        sprintf("numeric vector [%d values: %s, ...]", 
                               length(param_value), 
                               paste(as.character(utils::head(param_value, 3)), collapse = ", "))
                    } else {
                        paste0("c(", paste(as.character(param_value), collapse = ", "), ")")
                    }
                } else if (is.character(param_value)) {
                    if (length(param_value) == 1) {
                        param_value
                    } else if (length(param_value) > 5) {
                        sprintf("character vector [%d values: %s, ...]", 
                               length(param_value), 
                               paste(shQuote(utils::head(param_value, 3)), collapse = ", "))
                    } else {
                        paste0("c(", paste(shQuote(param_value), collapse = ", "), ")")
                    }
                } else {
                    # SAFE fallback - no dput() which was causing the hang
                    paste0("[", class(param_value)[1], " object]")
                }
            }, error = function(e) {
                "[SERIALIZATION ERROR]"
            })
        }
        
         s4_sections <- if (requireNamespace("purrr", quietly = TRUE)) {
             purrr::imap(s4_params, function(func_params, func_name) {
                 tryCatch({
                     if (!is.list(func_params) || length(func_params) == 0) return("")
                     
                     cat(sprintf("WORKFLOW ARGS: Processing S4 function group '%s'\n", func_name))
                     header <- sprintf("[%s]", func_name)
                     
                     # Use map instead of imap_chr for safer handling
                     param_lines <- purrr::imap(func_params, function(param_value, param_name) {
                         tryCatch({
                             param_str <- formatParameterValue(param_value)
                             sprintf("  %s = %s", param_name, param_str)
                         }, error = function(e) {
                             cat(sprintf("WORKFLOW ARGS: Error formatting parameter '%s' in '%s': %s\n", param_name, func_name, e$message))
                             sprintf("  %s = [ERROR: %s]", param_name, e$message)
                         })
                     })
                     
                     # Convert list to character vector safely
                     param_lines_char <- unlist(param_lines)
                     paste(c(header, param_lines_char, ""), collapse = "\n")
                 }, error = function(e) {
                     cat(sprintf("WORKFLOW ARGS: Error processing S4 function group '%s': %s\n", func_name, e$message))
                     sprintf("[%s]\n  [ERROR: %s]\n", func_name, e$message)
                 })
             })
         } else {
             # Fallback using base R if purrr not available
             lapply(names(s4_params), function(func_name) {
                 func_params <- s4_params[[func_name]]
                 if (!is.list(func_params) || length(func_params) == 0) return("")
                 
                 header <- sprintf("[%s]", func_name)
                 
                 param_lines <- lapply(names(func_params), function(param_name) {
                     param_value <- func_params[[param_name]]
                     param_str <- formatParameterValue(param_value)
                     sprintf("  %s = %s", param_name, param_str)
                 })
                 
                 paste(c(header, unlist(param_lines), ""), collapse = "\n")
             })
         }
        
        # Add all S4 sections at once - safely handle list from purrr::imap()
        s4_sections_char <- tryCatch({
            if (is.list(s4_sections)) {
                # Convert list to character vector
                unlist(s4_sections)
            } else {
                s4_sections
            }
        }, error = function(e) {
            cat(sprintf("WORKFLOW ARGS: Error converting S4 sections: %s\n", e$message))
            "[ERROR: Could not process S4 sections]"
        })
        
        config_lines <- c(config_lines, unlist(strsplit(paste(s4_sections_char, collapse = ""), "\n")))
        
        cat("WORKFLOW ARGS: S4 parameters formatted successfully\n")
    }
    
    # Add UI parameters from workflow_data if available (DE and Enrichment UI inputs)
    if (!is.null(workflow_data)) {
        cat("WORKFLOW ARGS: Checking for UI parameters in workflow_data\n")
        ui_sections <- c()
        
        # Check for DE UI parameters
        if (!is.null(workflow_data$de_ui_params)) {
            cat("WORKFLOW ARGS: Found DE UI parameters in workflow_data\n")
            de_ui_lines <- c(
                "[Differential Expression UI Parameters]",
                sprintf("  q_value_threshold = %s", ifelse(!is.null(workflow_data$de_ui_params$q_value_threshold), workflow_data$de_ui_params$q_value_threshold, "N/A")),
                sprintf("  log_fold_change_cutoff = %s", ifelse(!is.null(workflow_data$de_ui_params$log_fold_change_cutoff), workflow_data$de_ui_params$log_fold_change_cutoff, "N/A")),
                sprintf("  treat_enabled = %s", ifelse(!is.null(workflow_data$de_ui_params$treat_enabled), workflow_data$de_ui_params$treat_enabled, "N/A")),
                ""
            )
            ui_sections <- c(ui_sections, de_ui_lines)
        }
        
        # Check for Enrichment UI parameters  
        if (!is.null(workflow_data$enrichment_ui_params)) {
            cat("WORKFLOW ARGS: Found Enrichment UI parameters in workflow_data\n")
            enrichment_ui_lines <- c(
                "[Enrichment Analysis UI Parameters]",
                sprintf("  up_log2fc_cutoff = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$up_log2fc_cutoff), workflow_data$enrichment_ui_params$up_log2fc_cutoff, "N/A")),
                sprintf("  down_log2fc_cutoff = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$down_log2fc_cutoff), workflow_data$enrichment_ui_params$down_log2fc_cutoff, "N/A")),
                sprintf("  q_value_cutoff = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$q_value_cutoff), workflow_data$enrichment_ui_params$q_value_cutoff, "N/A")),
                sprintf("  organism_selected = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$organism_selected), workflow_data$enrichment_ui_params$organism_selected, "N/A")),
                sprintf("  database_source = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$database_source), workflow_data$enrichment_ui_params$database_source, "N/A")),
                ""
            )
            ui_sections <- c(ui_sections, enrichment_ui_lines)
        }
        
        # Add UI sections to config_lines if any were found
        if (length(ui_sections) > 0) {
            config_lines <- c(config_lines, 
                             "  User Interface Parameters:",
                             "  -------------------------", 
                             ui_sections)
            cat("WORKFLOW ARGS: Added UI parameters to output\n")
        } else {
            cat("WORKFLOW ARGS: No UI parameters found in workflow_data\n")
        }
    }
    
         # Format the remaining config list
     cat("WORKFLOW ARGS: About to format remaining config list\n")
     if (length(clean_config) > 0) {
         config_lines <- c(config_lines, 
                          "Additional Configuration Parameters:",
                          "-----------------------------------")
         
         cat("WORKFLOW ARGS: Calling formatConfigList...\n")
         tryCatch({
           # Check if formatConfigList exists before calling
           if (exists("formatConfigList", mode = "function")) {
               config_params <- formatConfigList(clean_config)
               config_lines <- c(config_lines, config_params)
               cat("WORKFLOW ARGS: formatConfigList completed successfully\n")
           } else {
               cat("WORKFLOW ARGS: formatConfigList function not found, using basic formatting\n")
               # Basic fallback formatting without for loops
               basic_config <- unlist(clean_config, recursive = TRUE)
               config_lines <- c(config_lines, names(basic_config), " = ", as.character(basic_config))
           }
         }, error = function(e) {
           cat(sprintf("WORKFLOW ARGS: formatConfigList failed: %s\n", e$message))
           config_lines <- c(config_lines, paste("Error formatting config:", e$message))
         })
     } else {
         config_lines <- c(config_lines, "No additional configuration parameters available")
     }
    
         cat("WORKFLOW ARGS: Adding config lines to output\n")
     output_lines <- c(output_lines, config_lines)
     
     # Add contrasts information
     cat("WORKFLOW ARGS: Processing contrasts information\n")
     if (!is.null(contrasts_tbl) && (is.data.frame(contrasts_tbl) || tibble::is_tibble(contrasts_tbl)) && nrow(contrasts_tbl) > 0) {
         cat("WORKFLOW ARGS: Adding contrasts to output\n")
         contrasts_lines <- c(
             "",
             "Contrasts:",
             "----------"
         )
         
         if ("contrasts" %in% colnames(contrasts_tbl)) {
             contrasts_info <- tryCatch({
                 contrasts_col <- contrasts_tbl[["contrasts"]]
                 paste("  ", as.character(contrasts_col))
             }, error = function(e) {
                 paste("  [Error extracting contrasts:", e$message, "]")
             })
             contrasts_lines <- c(contrasts_lines, contrasts_info)
         } else {
             contrasts_lines <- c(contrasts_lines, "  [Column 'contrasts' not found in contrasts_tbl]")
         }
         
         output_lines <- c(output_lines, contrasts_lines)
         cat("WORKFLOW ARGS: Contrasts added successfully\n")
     } else {
         cat("WORKFLOW ARGS: No contrasts to add\n")
     }
     
     # Write to file
     cat("WORKFLOW ARGS: About to write file\n")
     output_file <- file.path(source_dir_path, "study_parameters.txt")
     cat(sprintf("WORKFLOW ARGS: Target file path: %s\n", output_file))
     
     tryCatch({
         cat("WORKFLOW ARGS: Calling writeLines...\n")
         writeLines(output_lines, output_file)
         cat(sprintf("WORKFLOW ARGS: Study parameters saved to: %s\n", output_file))
         return(output_file)
     }, error = function(e) {
         cat(sprintf("WORKFLOW ARGS: Error writing file: %s\n", e$message))
         stop("Failed to write study parameters file: ", e$message)
     })
}

#' Check Missing Value Percentages in Peptide Data
#' 
#' @description Calculate and report the percentage of missing values (NAs) in peptide data
#' at different levels: total dataset, per sample, and per group.
#' 
#' @param peptide_obj A PeptideQuantitativeData S4 object
#' @param verbose Logical, whether to print detailed results (default: TRUE)
#' 
#' @return A list containing:
#' \itemize{
#'   \item total_na_percent: Overall percentage of NAs in the dataset
#'   \item per_sample_na: Data frame with NA percentages per sample
#'   \item per_group_na: Data frame with NA percentages per group
#'   \item summary_stats: Summary statistics of NA distribution
#' }
#' 
#' @export
checkPeptideNAPercentages <- function(peptide_obj, verbose = TRUE) {
  
  # Validate input
  if (!is(peptide_obj, "PeptideQuantitativeData")) {
    stop("Input must be a PeptideQuantitativeData S4 object")
  }
  
  # Extract data from S4 object
  peptide_matrix <- peptide_obj@peptide_matrix
  design_matrix <- peptide_obj@design_matrix
  sample_id_col <- peptide_obj@sample_id
  group_id_col <- peptide_obj@group_id
  
  # Validate that matrix and design matrix are compatible
  if (ncol(peptide_matrix) != nrow(design_matrix)) {
    stop("Number of samples in peptide_matrix doesn't match design_matrix rows")
  }
  
  # Calculate total NA percentage
  total_values <- length(peptide_matrix)
  total_nas <- sum(is.na(peptide_matrix))
  total_na_percent <- (total_nas / total_values) * 100
  
  # Calculate per-sample NA percentages
  sample_na_counts <- apply(peptide_matrix, 2, function(x) sum(is.na(x)))
  sample_na_percentages <- (sample_na_counts / nrow(peptide_matrix)) * 100
  
  per_sample_na <- data.frame(
    sample = colnames(peptide_matrix),
    na_count = sample_na_counts,
    na_percentage = sample_na_percentages,
    stringsAsFactors = FALSE
  )
  
  # Add group information to per-sample results
  per_sample_na <- merge(per_sample_na, design_matrix, 
                        by.x = "sample", by.y = sample_id_col, all.x = TRUE)
  
  # Calculate per-group NA percentages
  per_group_na <- per_sample_na %>%
    group_by(!!sym(group_id_col)) %>%
    summarise(
      num_samples = n(),
      mean_na_percentage = mean(na_percentage, na.rm = TRUE),
      median_na_percentage = median(na_percentage, na.rm = TRUE),
      min_na_percentage = min(na_percentage, na.rm = TRUE),
      max_na_percentage = max(na_percentage, na.rm = TRUE),
      sd_na_percentage = sd(na_percentage, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(mean_na_percentage)
  
  # Calculate summary statistics
  summary_stats <- list(
    total_peptides = nrow(peptide_matrix),
    total_samples = ncol(peptide_matrix),
    total_groups = length(unique(design_matrix[[group_id_col]])),
    total_values = total_values,
    total_nas = total_nas,
    mean_na_per_sample = mean(sample_na_percentages),
    median_na_per_sample = median(sample_na_percentages),
    min_na_per_sample = min(sample_na_percentages),
    max_na_per_sample = max(sample_na_percentages)
  )
  
  # Print results if verbose
  if (verbose) {
    cat("\n=== Peptide Data Missing Value Analysis ===\n")
    cat(sprintf("Dataset dimensions: %d peptides × %d samples\n", 
                nrow(peptide_matrix), ncol(peptide_matrix)))
    cat(sprintf("Number of groups: %d\n", summary_stats$total_groups))
    cat(sprintf("Total missing values: %s out of %s (%.2f%%)\n", 
                format(total_nas, big.mark = ","),
                format(total_values, big.mark = ","),
                total_na_percent))
    
    cat("\n--- Per-Sample Missing Value Summary ---\n")
    cat(sprintf("Mean NA%% per sample: %.2f%%\n", summary_stats$mean_na_per_sample))
    cat(sprintf("Median NA%% per sample: %.2f%%\n", summary_stats$median_na_per_sample))
    cat(sprintf("Range: %.2f%% - %.2f%%\n", 
                summary_stats$min_na_per_sample, summary_stats$max_na_per_sample))
    
    cat("\n--- Per-Group Missing Value Summary ---\n")
    print(per_group_na)
    
    cat("\n--- Samples with Highest Missing Values ---\n")
    top_missing_samples <- per_sample_na %>%
      arrange(desc(na_percentage)) %>%
      head(min(5, nrow(per_sample_na)))
    print(top_missing_samples[, c("sample", group_id_col, "na_percentage")])
    
    cat("\n--- Samples with Lowest Missing Values ---\n")
    bottom_missing_samples <- per_sample_na %>%
      arrange(na_percentage) %>%
      head(min(5, nrow(per_sample_na)))
    print(bottom_missing_samples[, c("sample", group_id_col, "na_percentage")])
  }
  
  # Return results
  results <- list(
    total_na_percent = total_na_percent,
    per_sample_na = per_sample_na,
    per_group_na = per_group_na,
    summary_stats = summary_stats
  )
  
  return(invisible(results))
}

#' Validate Post-Imputation Peptide Data
#' 
#' @description A simple wrapper to validate peptide data after imputation,
#' specifically checking if imputation was successful (should show 0% NAs).
#' 
#' @param peptide_obj A PeptideQuantitativeData S4 object (post-imputation)
#' @param expected_na_percent Expected NA percentage (default: 0 for post-imputation)
#' @param tolerance Tolerance for expected percentage (default: 0.1%)
#' 
#' @return Logical indicating if validation passed, with detailed output
#' 
#' @export
validatePostImputationData <- function(peptide_obj, expected_na_percent = 0, tolerance = 0.1) {
  
  cat("\n=== POST-IMPUTATION VALIDATION ===\n")
  
  # Run the full NA analysis
  na_results <- checkPeptideNAPercentages(peptide_obj, verbose = TRUE)
  
  # Check if imputation was successful
  actual_na_percent <- na_results$total_na_percent
  is_valid <- abs(actual_na_percent - expected_na_percent) <= tolerance
  
  cat("\n--- VALIDATION RESULT ---\n")
  cat(sprintf("Expected NA%%: %.2f%% (± %.2f%%)\n", expected_na_percent, tolerance))
  cat(sprintf("Actual NA%%: %.2f%%\n", actual_na_percent))
  
  if (is_valid) {
    cat("✓ VALIDATION PASSED: Imputation appears successful!\n")
  } else {
    cat("✗ VALIDATION FAILED: Unexpected NA percentage detected!\n")
    if (actual_na_percent > expected_na_percent + tolerance) {
      cat("  → Issue: More NAs than expected. Imputation may have failed.\n")
    } else {
      cat("  → Issue: Fewer NAs than expected. Check data integrity.\n")
    }
  }
  
  # Additional warnings for common issues
  if (actual_na_percent > 10) {
    cat("⚠ WARNING: High NA percentage suggests imputation problems!\n")
  }
  
  if (na_results$summary_stats$max_na_per_sample > actual_na_percent + 5) {
    cat("⚠ WARNING: Large variation in NA% between samples detected!\n")
  }
  
  cat("\n")
  return(invisible(list(
    is_valid = is_valid,
    actual_na_percent = actual_na_percent,
    expected_na_percent = expected_na_percent,
    full_results = na_results
  )))
}

#' Get Recommendations for Handling Protein-Level Missing Values
#' 
#' @description Provides specific recommendations for dealing with missing values
#' in protein data based on the percentage and distribution of NAs.
#' 
#' @param protein_obj A ProteinQuantitativeData S4 object
#' @param include_code Logical, whether to include example R code (default: TRUE)
#' 
#' @return Prints recommendations and invisibly returns a list of strategies
#' 
#' @export
getProteinNARecommendations <- function(protein_obj, include_code = TRUE) {
  
  # Get NA analysis
  na_results <- checkProteinNAPercentages(protein_obj, verbose = FALSE)
  na_percent <- na_results$total_na_percent
  
  cat("\n=== PROTEIN NA HANDLING RECOMMENDATIONS ===\n")
  cat(sprintf("Your data: %.1f%% NAs across %d proteins\n\n", 
              na_percent, na_results$summary_stats$total_proteins))
  
  if (na_percent < 15) {
    cat("🎯 RECOMMENDATION: Complete Case Analysis\n")
    cat("• Your data has excellent protein coverage\n")
    cat("• Can proceed with standard analysis on proteins with complete data\n")
    if (include_code) {
      cat("\n📝 Example code:\n")
      cat("complete_proteins <- protein_obj@protein_quant_table[complete.cases(protein_obj@protein_quant_table), ]\n")
    }
    
  } else if (na_percent >= 15 && na_percent < 40) {
    cat("🎯 RECOMMENDATION: Consider Protein-Level Imputation\n")
    cat("• Moderate missing values - imputation could be beneficial\n")
    cat("• Options: KNN, minimum value, or mixed imputation strategies\n")
    cat("• Alternative: Filter to proteins detected in ≥X samples per group\n")
    if (include_code) {
      cat("\n📝 Example filtering code:\n")
      cat("# Keep proteins detected in ≥50% of samples per group\n")
      cat("filtered_proteins <- filterProteinsByGroupDetection(protein_obj, min_detection_rate = 0.5)\n")
    }
    
  } else if (na_percent >= 40 && na_percent < 60) {
    cat("🎯 RECOMMENDATION: Strict Filtering + Targeted Imputation\n")
    cat("• High missing values suggest challenging sample/detection conditions\n")
    cat("• Focus on well-detected proteins (present in majority of samples)\n")
    cat("• Consider group-wise detection requirements\n")
    if (include_code) {
      cat("\n📝 Example approach:\n")
      cat("# Keep proteins detected in ≥70% of samples in at least one group\n")
      cat("robust_proteins <- filterProteinsByGroupwise(protein_obj, min_group_detection = 0.7)\n")
    }
    
  } else {
    cat("⚠️  RECOMMENDATION: Review Data Quality\n")
    cat("• Very high missing values (>60%) suggest potential issues\n")
    cat("• Check: sample quality, peptide identification, rollup parameters\n")
    cat("• Consider more stringent protein identification criteria\n")
    cat("• May need to focus only on highly abundant/well-detected proteins\n")
  }
  
  cat("\n📚 STRATEGIES SUMMARY:\n")
  cat("1. Complete Case: Use only proteins with no NAs\n")
  cat("2. Filtering: Remove proteins with >X% missing values\n")
  cat("3. Group-wise: Require detection in ≥Y% samples per group\n")
  cat("4. Imputation: Fill NAs with estimated values (KNN, minimum, etc.)\n")
  cat("5. Hybrid: Combine filtering + imputation\n")
  
  cat("\n💡 TIP: Protein NAs ≠ Data Quality Issues\n")
  cat("Missing proteins often reflect:\n")
  cat("• Low abundance proteins below detection limit\n")
  cat("• Sample-specific biology (some proteins not expressed)\n")
  cat("• Normal variation in complex proteomes\n\n")
  
  strategies <- list(
    na_percent = na_percent,
    primary_recommendation = if (na_percent < 15) "complete_case" 
                            else if (na_percent < 40) "imputation_or_filtering"
                            else if (na_percent < 60) "strict_filtering"
                            else "data_quality_review",
    alternative_strategies = c("complete_case", "group_wise_filtering", "imputation", "hybrid")
  )
  
  return(invisible(strategies))
}

#' Check Missing Value Percentages in Protein Data
#' 
#' @description Calculate and report the percentage of missing values (NAs) in protein data
#' at different levels: total dataset, per sample, and per group.
#' 
#' @param protein_obj A ProteinQuantitativeData S4 object
#' @param verbose Logical, whether to print detailed results (default: TRUE)
#' 
#' @return A list containing:
#' \itemize{
#'   \item total_na_percent: Overall percentage of NAs in the dataset
#'   \item per_sample_na: Data frame with NA percentages per sample
#'   \item per_group_na: Data frame with NA percentages per group
#'   \item summary_stats: Summary statistics of NA distribution
#' }
#' 
#' @export
checkProteinNAPercentages <- function(protein_obj, verbose = TRUE) {
  
  # Validate input
  if (!is(protein_obj, "ProteinQuantitativeData")) {
    stop("Input must be a ProteinQuantitativeData S4 object")
  }
  
  # Extract data from S4 object
  protein_quant_table <- protein_obj@protein_quant_table
  design_matrix <- protein_obj@design_matrix
  sample_id_col <- protein_obj@sample_id
  group_id_col <- protein_obj@group_id
  protein_id_col <- protein_obj@protein_id_column
  
  # Identify sample columns (exclude protein ID column)
  sample_columns <- setdiff(colnames(protein_quant_table), protein_id_col)
  
  # Validate that sample columns match design matrix
  if (length(sample_columns) != nrow(design_matrix)) {
    stop("Number of sample columns doesn't match design_matrix rows")
  }
  
  # Extract quantitative data matrix (samples only)
  protein_matrix <- as.matrix(protein_quant_table[, sample_columns])
  rownames(protein_matrix) <- protein_quant_table[[protein_id_col]]
  
  # Calculate total NA percentage
  total_values <- length(protein_matrix)
  total_nas <- sum(is.na(protein_matrix))
  total_na_percent <- (total_nas / total_values) * 100
  
  # Calculate per-sample NA percentages
  sample_na_counts <- apply(protein_matrix, 2, function(x) sum(is.na(x)))
  sample_na_percentages <- (sample_na_counts / nrow(protein_matrix)) * 100
  
  per_sample_na <- data.frame(
    sample = names(sample_na_counts),
    na_count = sample_na_counts,
    na_percentage = sample_na_percentages,
    stringsAsFactors = FALSE
  )
  
  # Add group information to per-sample results
  per_sample_na <- merge(per_sample_na, design_matrix, 
                        by.x = "sample", by.y = sample_id_col, all.x = TRUE)
  
  # Calculate per-group NA percentages
  per_group_na <- per_sample_na %>%
    group_by(!!sym(group_id_col)) %>%
    summarise(
      num_samples = n(),
      mean_na_percentage = mean(na_percentage, na.rm = TRUE),
      median_na_percentage = median(na_percentage, na.rm = TRUE),
      min_na_percentage = min(na_percentage, na.rm = TRUE),
      max_na_percentage = max(na_percentage, na.rm = TRUE),
      sd_na_percentage = sd(na_percentage, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(mean_na_percentage)
  
  # Calculate summary statistics
  summary_stats <- list(
    total_proteins = nrow(protein_matrix),
    total_samples = ncol(protein_matrix),
    total_groups = length(unique(design_matrix[[group_id_col]])),
    total_values = total_values,
    total_nas = total_nas,
    mean_na_per_sample = mean(sample_na_percentages),
    median_na_per_sample = median(sample_na_percentages),
    min_na_per_sample = min(sample_na_percentages),
    max_na_per_sample = max(sample_na_percentages)
  )
  
  # Print results if verbose
  if (verbose) {
    cat("\n=== Protein Data Missing Value Analysis ===\n")
    cat(sprintf("Dataset dimensions: %d proteins × %d samples\n", 
                nrow(protein_matrix), ncol(protein_matrix)))
    cat(sprintf("Number of groups: %d\n", summary_stats$total_groups))
    cat(sprintf("Total missing values: %s out of %s (%.2f%%)\n", 
                format(total_nas, big.mark = ","),
                format(total_values, big.mark = ","),
                total_na_percent))
    
    cat("\n--- Per-Sample Missing Value Summary ---\n")
    cat(sprintf("Mean NA%% per sample: %.2f%%\n", summary_stats$mean_na_per_sample))
    cat(sprintf("Median NA%% per sample: %.2f%%\n", summary_stats$median_na_per_sample))
    cat(sprintf("Range: %.2f%% - %.2f%%\n", 
                summary_stats$min_na_per_sample, summary_stats$max_na_per_sample))
    
    cat("\n--- Per-Group Missing Value Summary ---\n")
    print(per_group_na)
    
    cat("\n--- Samples with Highest Missing Values ---\n")
    top_missing_samples <- per_sample_na %>%
      arrange(desc(na_percentage)) %>%
      head(min(5, nrow(per_sample_na)))
    print(top_missing_samples[, c("sample", group_id_col, "na_percentage")])
    
    cat("\n--- Samples with Lowest Missing Values ---\n")
    bottom_missing_samples <- per_sample_na %>%
      arrange(na_percentage) %>%
      head(min(5, nrow(per_sample_na)))
    print(bottom_missing_samples[, c("sample", group_id_col, "na_percentage")])
  }
  
  # Return results
  results <- list(
    total_na_percent = total_na_percent,
    per_sample_na = per_sample_na,
    per_group_na = per_group_na,
    summary_stats = summary_stats
  )
  
  return(invisible(results))
}

#' Validate Post-Imputation Protein Data
#' 
#' @description A simple wrapper to validate protein data after imputation,
#' specifically checking if imputation was successful.
#' 
#' @param protein_obj A ProteinQuantitativeData S4 object (post-imputation)
#' @param expected_na_percent Expected NA percentage (default: varies based on protein data)
#' @param tolerance Tolerance for expected percentage (default: 10%)
#' 
#' @return Logical indicating if validation passed, with detailed output
#' 
#' @export
validatePostImputationProteinData <- function(protein_obj, expected_na_percent = NULL, tolerance = 10) {
  
  cat("\n=== POST-IMPUTATION PROTEIN DATA VALIDATION ===\n")
  cat("Note: Protein-level NAs occur even after peptide imputation because:\n")
  cat("• Proteins need ≥1 detected peptide to get a quantification\n")
  cat("• Some proteins detected only in subset of samples\n")
  cat("• This is normal proteomics data behavior!\n\n")
  
  # Run the full NA analysis
  na_results <- checkProteinNAPercentages(protein_obj, verbose = TRUE)
  
  # Set expected NA percentage if not provided (proteins often have some NAs)
  if (is.null(expected_na_percent)) {
    # For protein data, NAs are very common due to missing peptides/proteins
    # Typical ranges: 20-50% depending on sample complexity and detection method
    expected_na_percent <- 35  # Realistic expectation for protein data
    cat(sprintf("Note: Using default expected NA%% of %.1f%% for protein data\n", expected_na_percent))
    cat("(Protein-level NAs are normal due to incomplete protein detection across samples)\n")
  }
  
  # Check if validation passes
  actual_na_percent <- na_results$total_na_percent
  is_valid <- abs(actual_na_percent - expected_na_percent) <= tolerance
  
  cat("\n--- VALIDATION RESULT ---\n")
  cat(sprintf("Expected NA%%: %.2f%% (± %.2f%%)\n", expected_na_percent, tolerance))
  cat(sprintf("Actual NA%%: %.2f%%\n", actual_na_percent))
  
  if (is_valid) {
    cat("✓ VALIDATION PASSED: Protein data NA levels are within expected range!\n")
  } else {
    cat("✗ VALIDATION FAILED: Unexpected NA percentage detected!\n")
    if (actual_na_percent > expected_na_percent + tolerance) {
      cat("  → Issue: More NAs than expected. Check for missing proteins/peptides.\n")
    } else {
      cat("  → Issue: Fewer NAs than expected. Possible over-imputation.\n")
    }
  }
  
  # Additional warnings for common issues
  if (actual_na_percent > 50) {
    cat("⚠ WARNING: Very high NA percentage (>50%) suggests data quality issues!\n")
  }
  
  if (actual_na_percent < 10) {
    cat("ℹ INFO: Very low NA percentage (<10%) - excellent protein coverage!\n")
  }
  
  # Educational information about protein NAs
  if (actual_na_percent > 20 && actual_na_percent < 50) {
    cat("ℹ INFO: NA percentage is typical for protein-level data\n")
    cat("  → This reflects biological reality: not all proteins detected in all samples\n")
    cat("  → Consider: protein-level imputation OR complete-case analysis\n")
  }
  
  if (na_results$summary_stats$max_na_per_sample > actual_na_percent + 10) {
    cat("⚠ WARNING: Large variation in NA% between samples detected!\n")
    cat("  → Some samples may have much lower protein coverage.\n")
  }
  
  # Check for problematic samples (>80% missing)
  high_missing_samples <- na_results$per_sample_na[na_results$per_sample_na$na_percentage > 80, ]
  if (nrow(high_missing_samples) > 0) {
    cat("⚠ WARNING: Samples with >80% missing proteins detected:\n")
    print(high_missing_samples[, c("sample", "na_percentage")])
  }
  
  cat("\n")
  return(invisible(list(
    is_valid = is_valid,
    actual_na_percent = actual_na_percent,
    expected_na_percent = expected_na_percent,
    full_results = na_results
  )))
}

