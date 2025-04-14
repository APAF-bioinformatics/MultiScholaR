#' @description Creates and manages project directories with version control for multiple omic types.
#' @param base_dir Base directory path (optional, defaults to here::here())
#' @param omic_types Character vector or comma/space-separated string indicating the types of omics data (e.g., c("proteomics", "metabolomics"), "proteomics, transcriptomics").
#' @param label Optional label to append to each omic type directory name (e.g., "proteomics_MyLabel").
#' @param force Logical; if TRUE, skips user confirmation for overwriting existing directories (default: FALSE).
#' @return A named list where each name is an omic type. Each element is a list containing the relevant directory paths for that omic type.
#' @importFrom logger log_info log_warn
#' @importFrom stringr str_split str_trim
#' @export
multisetupAndShowDirectories <- function(base_dir = here::here(), omic_types, label = NULL, force = FALSE) {
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
        stop("`omic_types` must be a character vector or a single comma/space-separated string.")
    }

    if (length(parsed_omic_types) == 0) {
        stop("No valid `omic_types` provided.")
    }

    # Define Omics-Specific Configuration
    valid_omic_types <- c("proteomics", "metabolomics", "transcriptomics") # Expand as needed

    # Validate all provided types
    invalid_types <- setdiff(parsed_omic_types, valid_omic_types)
    if (length(invalid_types) > 0) {
        stop(
            "Invalid omic_type(s) specified: ", paste(invalid_types, collapse = ", "),
            ". Choose from: ", paste(valid_omic_types, collapse = ", ")
        )
    }

    # --- Initialization ---
    all_created_paths <- list() # List to store paths for each omic type
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S") # Single timestamp for the run

    # --- Loop Through Each Omic Type ---
    for (current_omic_type in parsed_omic_types) {
        logger::log_info("Processing setup for omic type: {current_omic_type}")

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
                    "metabolite_qc", "feature_qc", "de_metabolites", # Note: feature_qc is generic here
                    "publication_graphs", "pathway_enrichment",
                    file.path("publication_graphs", "filtering_qc")
                ),
                results_summary_subdirs = c("QC_figures", "Publication_figures", "Publication_tables", "Study_report"),
                scripts_source_leaf = "metabolomics",
                global_vars = list(
                    de_output_leaf = "de_metabolites",
                    pathway_leaf = "pathway_enrichment",
                    feature_qc_leaf = "metabolite_qc", # Specific leaf for path construction
                    subfeature_qc_leaf = NULL,
                    clean_features_leaf = "clean_metabolites"
                )
            ),
            transcriptomics = list(
                results_subdirs = c(
                    "gene_qc", "count_data", "de_genes",
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
                    clean_features_leaf = "normalized_counts"
                )
            ),
            stop("Internal error: Configuration not found for validated omic_type: ", current_omic_type) # Should not happen
        )

        # Construct Directory Names for the current omic type
        omic_label_dirname <- paste0(current_omic_type, if (!is.null(label)) paste0("_", substr(label, 1, 30)) else "")

        results_path <- file.path(base_dir, "results", omic_label_dirname)
        results_summary_path <- file.path(base_dir, "results_summary", omic_label_dirname)
        scripts_path <- file.path(base_dir, "scripts", omic_label_dirname) # Destination scripts path for this omic

        # Check Existing Directories and User Confirmation (only if force = FALSE)
        if (!force && (dir.exists(results_path) || dir.exists(results_summary_path) || dir.exists(scripts_path))) {
            cat(sprintf("\nWarning: Directory(ies) for '%s' already exist:\n", omic_label_dirname))
            if (dir.exists(results_path)) cat(sprintf("- %s\n", results_path))
            if (dir.exists(results_summary_path)) cat(sprintf("- %s\n", results_summary_path))
            if (dir.exists(scripts_path)) cat(sprintf("- %s\n", scripts_path))

            response <- readline(prompt = sprintf("Do you want to overwrite directories for '%s'? (y/n): ", current_omic_type))

            if (tolower(substr(response, 1, 1)) != "y") {
                logger::log_info("Skipping setup for {current_omic_type} due to user cancellation.")
                next # Skip to the next omic type in the loop
            } else {
                logger::log_info("Overwriting existing directories for {current_omic_type} as requested.")
                # Remove existing directories if user confirmed overwrite
                if (dir.exists(results_path)) unlink(results_path, recursive = TRUE, force = TRUE)
                if (dir.exists(results_summary_path)) unlink(results_summary_path, recursive = TRUE, force = TRUE)
                if (dir.exists(scripts_path)) unlink(scripts_path, recursive = TRUE, force = TRUE)
            }
        } else if (force && (dir.exists(results_path) || dir.exists(results_summary_path) || dir.exists(scripts_path))) {
            logger::log_info("Force is TRUE. Overwriting existing directories for {current_omic_type}.")
            # Remove existing directories if force is TRUE
            if (dir.exists(results_path)) unlink(results_path, recursive = TRUE, force = TRUE)
            if (dir.exists(results_summary_path)) unlink(results_summary_path, recursive = TRUE, force = TRUE)
            if (dir.exists(scripts_path)) unlink(scripts_path, recursive = TRUE, force = TRUE)
        }


        # Define Paths for the current omic type
        publication_graphs_dir_base <- file.path(results_path, "publication_graphs")
        qc_dir_base <- file.path(publication_graphs_dir_base, "filtering_qc")
        # Ensure the common data directory exists (create if not, idempotent)
        common_data_dir <- file.path(base_dir, "data")
        dir.create(common_data_dir, recursive = TRUE, showWarnings = FALSE)

        # List to hold paths specific to the current omic type
        current_omic_paths <- list(
            results_base = results_path,
            results_subdirs = file.path(results_path, omic_config$results_subdirs),
            results_summary_base = results_summary_path,
            results_summary_subdirs = file.path(results_summary_path, omic_config$results_summary_subdirs),
            data_dir = common_data_dir, # Reference the common data directory
            scripts_dest_dir = scripts_path,
            scripts_source_dir = file.path(base_dir, "scripts", omic_config$scripts_source_leaf),
            publication_graphs_dir = publication_graphs_dir_base,
            qc_dir = qc_dir_base,
            # Note: Timestamped dir is specific to this omic's results
            time_dir = file.path(qc_dir_base, timestamp)
        )

        # --- Create Directories for the current omic type ---
        # Create results base and subdirs
        dir.create(current_omic_paths$results_base, recursive = TRUE, showWarnings = FALSE)
        invisible(sapply(current_omic_paths$results_subdirs, dir.create, recursive = TRUE, showWarnings = FALSE))

        # Create results_summary base and subdirs
        dir.create(current_omic_paths$results_summary_base, recursive = TRUE, showWarnings = FALSE)
        invisible(sapply(current_omic_paths$results_summary_subdirs, dir.create, recursive = TRUE, showWarnings = FALSE))

        # Create timestamped qc dir (includes parents like qc_dir and pub_graphs)
        dir.create(current_omic_paths$time_dir, recursive = TRUE, showWarnings = FALSE)

        # --- Handle Scripts Directory Copying for the current omic type ---
        if (dir.exists(current_omic_paths$scripts_source_dir)) {
            dir.create(current_omic_paths$scripts_dest_dir, recursive = TRUE, showWarnings = FALSE) # Create destination

            script_files <- list.files(current_omic_paths$scripts_source_dir, full.names = TRUE, recursive = TRUE)
            # Exclude Rmd files and renv directory/lockfile
            script_files <- script_files[!grepl("renv/|renv\\.lock|\\.Rmd$", script_files, ignore.case = TRUE)]

            if (length(script_files) > 0) {
                copied_scripts_ok <- TRUE
                for (f in script_files) {
                    # Calculate relative path from source base to file f
                    rel_path <- sub(paste0("^", tools::file_path_as_absolute(current_omic_paths$scripts_source_dir), "/?"), "", tools::file_path_as_absolute(f))
                    dest_file <- file.path(current_omic_paths$scripts_dest_dir, rel_path)
                    # Ensure destination subdirectory exists before copying
                    dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
                    if (!file.copy(f, dest_file, overwrite = TRUE)) {
                        logger::log_warn("Failed to copy script: {f} to {dest_file}")
                        copied_scripts_ok <- FALSE
                    }
                }
                if (copied_scripts_ok) {
                    logger::log_info("Copied scripts for {current_omic_type} from {current_omic_paths$scripts_source_dir} to {current_omic_paths$scripts_dest_dir}")
                } else {
                    logger::log_warn("Some script files failed to copy for {current_omic_type}.")
                }
            } else {
                logger::log_info("No scripts (excluding Rmd/renv) found in {current_omic_paths$scripts_source_dir} to copy.")
            }
        } else {
            # Still create the destination directory even if source doesn't exist
            dir.create(current_omic_paths$scripts_dest_dir, recursive = TRUE, showWarnings = FALSE)
            logger::log_warn("Source script directory not found for {current_omic_type}: {current_omic_paths$scripts_source_dir}. No scripts copied.")
        }


        # --- Build and Store Path Variables for the current omic type ---
        # Start with generic names relative to the current omic's paths
        omic_specific_paths_list <- list(
            base_dir = base_dir, # Overall project base
            results_dir = current_omic_paths$results_base,
            data_dir = current_omic_paths$data_dir, # Common data dir
            source_dir = current_omic_paths$scripts_dest_dir, # Destination scripts dir for this omic
            de_output_dir = file.path(current_omic_paths$results_base, omic_config$global_vars$de_output_leaf),
            publication_graphs_dir = current_omic_paths$publication_graphs_dir,
            timestamp = timestamp, # Same timestamp across all omics in this run
            qc_dir = current_omic_paths$qc_dir,
            time_dir = current_omic_paths$time_dir,
            results_summary_dir = current_omic_paths$results_summary_base,
            pathway_dir = file.path(current_omic_paths$results_base, omic_config$global_vars$pathway_leaf),
            feature_qc_dir = file.path(current_omic_paths$results_base, omic_config$global_vars$feature_qc_leaf), # Generic name
            clean_features_dir = file.path(current_omic_paths$results_base, omic_config$global_vars$clean_features_leaf), # Generic name
            qc_figures_dir = file.path(current_omic_paths$results_summary_base, "QC_figures"),
            publication_figures_dir = file.path(current_omic_paths$results_summary_base, "Publication_figures"),
            publication_tables_dir = file.path(current_omic_paths$results_summary_base, "Publication_tables"),
            study_report_dir = file.path(current_omic_paths$results_summary_base, "Study_report")
        )

        # Add optional sub-feature QC dir if defined in config
        if (!is.null(omic_config$global_vars$subfeature_qc_leaf)) {
            subfeature_qc_dir_path <- file.path(current_omic_paths$results_base, omic_config$global_vars$subfeature_qc_leaf)
            omic_specific_paths_list$subfeature_qc_dir <- subfeature_qc_dir_path # Add generic name
            # Ensure the directory exists (should be covered by results_subdirs creation, but belt-and-suspenders)
            dir.create(subfeature_qc_dir_path, recursive = TRUE, showWarnings = FALSE)
        }

        # Dynamically add specific variable names based on omic_type (e.g., protein_qc_dir)
        # These map the generic names (like feature_qc_dir) to the specific names for convenience
        specific_global_names <- switch(current_omic_type,
            proteomics = list(
                protein_qc_dir = omic_specific_paths_list$feature_qc_dir,
                peptide_qc_dir = omic_specific_paths_list$subfeature_qc_dir, # Will be NULL if subfeature not defined
                clean_proteins_dir = omic_specific_paths_list$clean_features_dir
            ),
            metabolomics = list(
                metabolite_qc_dir = omic_specific_paths_list$feature_qc_dir,
                # no subfeature defined for metabolomics in this config
                clean_metabolites_dir = omic_specific_paths_list$clean_features_dir
            ),
            transcriptomics = list(
                gene_qc_dir = omic_specific_paths_list$feature_qc_dir,
                count_data_dir = file.path(current_omic_paths$results_base, "count_data"), # Specific dir from results_subdirs
                normalized_counts_dir = omic_specific_paths_list$clean_features_dir
                # no subfeature defined for transcriptomics
            ),
            list() # Default empty list if type not matched (shouldn't happen)
        )

        # Ensure specific directories that might not be generic exist
        # Example: transcriptomics count_data_dir
        if (current_omic_type == "transcriptomics" && !is.null(specific_global_names$count_data_dir)) {
            dir.create(specific_global_names$count_data_dir, recursive = TRUE, showWarnings = FALSE)
        }


        # Add non-NULL specific names to the list for this omic type
        # Filter out NULLs which occur if e.g. subfeature_qc_dir wasn't defined
        omic_specific_paths_list <- c(omic_specific_paths_list, specific_global_names[!sapply(specific_global_names, is.null)])

        # Store the collected paths for this omic type in the main list
        all_created_paths[[current_omic_type]] <- omic_specific_paths_list
        logger::log_info("Stored paths for {current_omic_type}.")
    } # End loop over omic types

    # --- Print Structure ---
    cat("\n--- Directory Structure Created/Verified ---\n")
    if (length(all_created_paths) > 0) {
        for (omic_name in names(all_created_paths)) {
            cat(sprintf("\nPaths for Omics Type: '%s'\n", omic_name))
            current_paths <- all_created_paths[[omic_name]]
            # Filter to print only character paths, exclude the timestamp itself if desired
            print_paths <- current_paths[sapply(current_paths, is.character)]

            # Order paths for slightly nicer printing (optional)
            print_paths <- print_paths[order(names(print_paths))]

            invisible(lapply(names(print_paths), function(path_name) {
                p <- print_paths[[path_name]]
                # Check if it looks like a directory path intended to exist
                is_dir_path <- endsWith(path_name, "_dir") || endsWith(path_name, "_base") || path_name %in% c("qc_dir", "time_dir", "data_dir")

                if (is_dir_path && dir.exists(p)) {
                    # Use tryCatch for potential permission issues when listing files
                    tryCatch(
                        {
                            # List files, excluding '.' and '..' for count
                            # Note: This can be slow on very large directories. Consider removing file count if performance is an issue.
                            # all.files=TRUE includes hidden files (like .gitkeep if used)
                            # no.. = TRUE excludes '.' and '..' from the list itself
                            files_in_dir <- list.files(p, recursive = TRUE, all.files = TRUE, no.. = TRUE)
                            # Filter out directory entries themselves if recursive, count only files
                            # This interpretation might need refinement based on desired count semantics
                            # file_count <- sum(!dir.exists(file.path(p, files_in_dir))) # Count only files - potentially slow
                            file_count <- length(files_in_dir) # Count files and subdirs

                            cat(sprintf("  %-25s : %s (%d items)\n", path_name, p, file_count))
                        },
                        error = function(e) {
                            cat(sprintf("  %-25s : %s (Error accessing content: %s)\n", path_name, p, e$message))
                        }
                    )
                } else if (is_dir_path) {
                    # Path is intended as a dir, but doesn't exist (shouldn't happen if creation worked)
                    cat(sprintf("  %-25s : %s (Directory not found!)\n", path_name, p))
                } else {
                    # Likely not a directory path (e.g., timestamp, base_dir if it's just the root)
                    cat(sprintf("  %-25s : %s\n", path_name, p))
                }
            }))
        }
    } else {
        cat("No directories were set up (possibly due to user cancellation for all types).\n")
    }
    cat("------------------------------------------\n")


    # Return the collected paths invisibly
    invisible(all_created_paths)
}
