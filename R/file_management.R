#' @description Creates and manages project directories with version control for multiple omic types.
#' @param base_dir Base directory path (optional, defaults to here::here())
#' @param omic_types Character vector or comma/space-separated string indicating the types of omics data (e.g., c("proteomics", "metabolomics"), "proteomics, transcriptomics").
#' @param label Optional label to append to each omic type directory name (e.g., "proteomics_MyLabel").
#' @param force Logical; if TRUE, skips user confirmation for overwriting existing directories (default: FALSE).
#' @return A named list where each name is an omic type. Each element is a list containing the relevant directory paths for that omic type.
#' @importFrom stringr str_split str_trim
#' @export
setupDirectories <- function(base_dir = here::here(), omic_types, label = NULL, force = FALSE) {
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
        stop("No valid `omic_types` provided after parsing.")
    }

    # Define Omics-Specific Configuration
    valid_omic_types <- c("proteomics", "metabolomics", "transcriptomics", "lipidomics", "integration") # Expanded

    # Validate all provided types
    invalid_types <- setdiff(parsed_omic_types, valid_omic_types)
    if (length(invalid_types) > 0) {
        err_msg <- sprintf("Invalid omic_type(s) specified: %s. Choose from: %s.",
                           paste(invalid_types, collapse = ", "),
                           paste(valid_omic_types, collapse = ", "))
        stop(err_msg)
    }

    # --- Initialization ---
    all_created_paths <- list() # List to store paths for each omic type
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S") # Single timestamp for the run

    # --- Loop Through Each Omic Type ---
    for (current_omic_type in parsed_omic_types) {
        message("--- Processing setup for omic type: ", current_omic_type, " ---")

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

        # Check Existing Directories and User Confirmation (only if force = FALSE)
        existing_dirs_found <- dir.exists(results_path) || dir.exists(results_summary_path) || dir.exists(scripts_path)
        
        # SAFETY: Prevent accidental deletion of template directories (those without labels)
        template_scripts_path <- file.path(base_dir, "scripts", current_omic_type)
        if (is.null(label) && scripts_path == template_scripts_path) {
            message("ERROR: Cannot overwrite template directory: ", template_scripts_path, ". This contains original templates needed for analysis directories.")
            message("Please use setupDirectories() with a label parameter (e.g., label = 'your_analysis') to create analysis-specific directories.")
            stop("Attempted to overwrite template directory without a label. Use a label to create analysis directories.")
        }

        if (!force && existing_dirs_found) {
            message("WARNING: Analysis directories for '", current_omic_type, "' with label '", omic_label_dirname, "' already exist:")
            if (dir.exists(results_path)) message("- ", results_path)
            if (dir.exists(results_summary_path)) message("- ", results_summary_path)
            if (dir.exists(scripts_path)) message("- ", scripts_path)
            
            if (!is.null(label)) {
                message("Note: Template directory ", template_scripts_path, " will be preserved and used as source.")
            }

            response_overwrite <- readline(prompt = sprintf("Overwrite ANALYSIS directories for '%s'? (y/n): ", omic_label_dirname))

            if (tolower(substr(response_overwrite, 1, 1)) == "y") {
                message("Overwriting existing directories for '", omic_label_dirname, "' as requested.")
                if (dir.exists(results_path)) unlink(results_path, recursive = TRUE, force = TRUE)
                if (dir.exists(results_summary_path)) unlink(results_summary_path, recursive = TRUE, force = TRUE)
                if (dir.exists(scripts_path)) unlink(scripts_path, recursive = TRUE, force = TRUE)
                # reuse_current_omic_dirs remains FALSE, proceed to create new structure
            } else {
                response_reuse <- readline(prompt = sprintf("Use existing directory structure for '%s' this session? (y/n): ", omic_label_dirname))
                if (tolower(substr(response_reuse, 1, 1)) == "y") {
                    message("Reusing existing directory structure for '", omic_label_dirname, "'.")
                    reuse_current_omic_dirs <- TRUE
                } else {
                    message("Setup for '", omic_label_dirname, "' cancelled by user. Directory variables will not be set for this omic type.")
                    process_current_omic <- FALSE
                }
            }
        } else if (force && existing_dirs_found) {
            # Additional safety check for force mode
            if (is.null(label) && scripts_path == template_scripts_path) {
                message("ERROR: Force mode cannot overwrite template directory: ", template_scripts_path)
                stop("Force mode attempted to overwrite template directory without a label.")
            }
            message("Force=TRUE: Overwriting existing ANALYSIS directories for '", omic_label_dirname, "'.")
            if (!is.null(label)) {
                message("Template directory ", template_scripts_path, " will be preserved.")
            }
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
        message("Ensured omic-specific data directory exists: ", omic_specific_data_dir)

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
        message("Ensured timestamped directory exists: ", current_omic_paths_def$time_dir)


        # --- Create Directories and Copy Scripts (if not reusing) ---
        if (!reuse_current_omic_dirs) {
            message("Creating directory structure for '", omic_label_dirname, "'...")

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

                # Get all files including .Rmd templates but exclude renv files
                script_files <- list.files(current_omic_paths_def$scripts_source_dir, full.names = TRUE, recursive = TRUE)
                script_files <- script_files[!grepl("renv/|renv\\.lock$", script_files, ignore.case = TRUE)]

                if (length(script_files) > 0) {
                    message("Copying scripts and templates for ", current_omic_type, " from ", current_omic_paths_def$scripts_source_dir, " to ", current_omic_paths_def$scripts_dest_dir, "...")
                    copied_scripts_ok <- TRUE
                    for (f in script_files) {
                        rel_path <- sub(paste0("^", tools::file_path_as_absolute(current_omic_paths_def$scripts_source_dir), "/?"), "", tools::file_path_as_absolute(f))
                        dest_file <- file.path(current_omic_paths_def$scripts_dest_dir, rel_path)
                        dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
                        if (!file.copy(f, dest_file, overwrite = TRUE)) {
                            message("WARNING: Failed to copy script: ", f, " to ", dest_file)
                            copied_scripts_ok <- FALSE
                        }
                    }
                    if (!copied_scripts_ok) {
                        message("WARNING: Some script files failed to copy for ", current_omic_type, ".")
                    } else {
                        message("Successfully copied ", length(script_files), " files including any report templates for ", current_omic_type, ".")
                    }
                } else {
                    message("No scripts found in ", current_omic_paths_def$scripts_source_dir, " to copy for ", current_omic_type, ".")
                }
            } else {
                dir.create(current_omic_paths_def$scripts_dest_dir, recursive = TRUE, showWarnings = FALSE) # Still create dest dir
                message("WARNING: Source script directory not found for ", current_omic_type, ": ", current_omic_paths_def$scripts_source_dir, ". No scripts copied.")
            }
        } else {
            message("Skipping directory creation and script copying for '", omic_label_dirname, "' as existing structure is being reused.")
            # Timestamped dir (current_omic_paths_def$time_dir) is already created.
        }


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
            message("Ensured UniProt annotation directory exists: ", uniprot_annotation_path)
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
        message("Stored paths for ", omic_label_dirname, ".")
    } # End loop over omic types

    # --- Print Structure ---
    message("--- Directory Structure Setup Complete ---")
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
        message("No directories were set up for any omic type (e.g., setup cancelled for all).")
    }
    cat("------------------------------------------\n")


    # Return the collected paths invisibly
    invisible(all_created_paths)
}
