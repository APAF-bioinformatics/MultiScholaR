# ----------------------------------------------------------------------------
# parseSetupDirectoriesOmicTypes
# ----------------------------------------------------------------------------
parseSetupDirectoriesOmicTypes <- function(omic_types) {
    # Handle comma/space separated string or vector input.
    if (is.character(omic_types) && length(omic_types) == 1) {
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

    valid_omic_types <- c("proteomics", "metabolomics", "transcriptomics", "lipidomics", "integration")
    invalid_types <- setdiff(parsed_omic_types, valid_omic_types)

    if (length(invalid_types) > 0) {
        err_msg <- sprintf(
            "Invalid omic_type(s) specified: %s. Choose from: %s.",
            paste(invalid_types, collapse = ", "),
            paste(valid_omic_types, collapse = ", ")
        )
        logger::log_error(err_msg)
        stop(err_msg)
    }

    parsed_omic_types
}

# ----------------------------------------------------------------------------
# getSetupDirectoriesOmicConfig
# ----------------------------------------------------------------------------
getSetupDirectoriesOmicConfig <- function(current_omic_type) {
    switch(current_omic_type,
        proteomics = list(
            results_subdirs = c(
                "protein_qc", "peptide_qc", "clean_proteins", "da_proteins",
                "publication_graphs", "pathway_enrichment",
                file.path("publication_graphs", "filtering_qc")
            ),
            results_summary_subdirs = c("QC_figures", "Publication_figures", "Publication_tables", "Study_report"),
            scripts_source_leaf = "proteomics",
            global_vars = list(
                da_output_leaf = "da_proteins",
                pathway_leaf = "pathway_enrichment",
                feature_qc_leaf = "protein_qc",
                subfeature_qc_leaf = "peptide_qc",
                clean_features_leaf = "clean_proteins"
            )
        ),
        metabolomics = list(
            results_subdirs = c(
                "metabolite_qc", "feature_qc", "clean_metabolites", "da_metabolites",
                "publication_graphs", "pathway_enrichment",
                file.path("publication_graphs", "filtering_qc")
            ),
            results_summary_subdirs = c("QC_figures", "Publication_figures", "Publication_tables", "Study_report"),
            scripts_source_leaf = "metabolomics",
            global_vars = list(
                da_output_leaf = "da_metabolites",
                pathway_leaf = "pathway_enrichment",
                feature_qc_leaf = "metabolite_qc",
                subfeature_qc_leaf = NULL, # No peptide equivalent for metabolomics typically
                clean_features_leaf = "clean_metabolites"
            )
        ),
        transcriptomics = list(
            results_subdirs = c(
                "gene_qc", "count_data", "normalized_counts", "da_genes",
                "publication_graphs", "pathway_enrichment",
                file.path("publication_graphs", "filtering_qc")
            ),
            results_summary_subdirs = c("QC_figures", "Publication_figures", "Publication_tables", "Study_report"),
            scripts_source_leaf = "transcriptomics",
            global_vars = list(
                da_output_leaf = "da_genes",
                pathway_leaf = "pathway_enrichment",
                feature_qc_leaf = "gene_qc",
                subfeature_qc_leaf = NULL,
                clean_features_leaf = "normalized_counts", # For normalized data, distinct from raw counts
                raw_counts_leaf = "count_data" # Specific for transcriptomics
            )
        ),
        lipidomics = list( # New configuration for lipidomics
            results_subdirs = c(
                "lipid_qc", "feature_qc", "clean_lipids", "da_lipids",
                "publication_graphs", "pathway_enrichment",
                file.path("publication_graphs", "filtering_qc")
            ),
            results_summary_subdirs = c("QC_figures", "Publication_figures", "Publication_tables", "Study_report"),
            scripts_source_leaf = "lipidomics",
            global_vars = list(
                da_output_leaf = "da_lipids",
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
}

# ----------------------------------------------------------------------------
# buildSetupDirectoriesPathList
# ----------------------------------------------------------------------------
buildSetupDirectoriesPathList <- function(
    base_dir,
    current_omic_type,
    omic_label_dirname,
    timestamp,
    current_omic_paths_def,
    omic_config
) {
    # Create Integration Directory (Shared across all omics)
    integration_dir <- file.path(base_dir, "integration")
    dir.create(integration_dir, recursive = TRUE, showWarnings = FALSE)

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
                # The actual variable name will be like 'da_output_dir', 'feature_qc_dir'
                # If var_name_suffix is "de_output_leaf", actual name is "da_output_dir"
                actual_var_name <- sub("_leaf$", "_dir", var_name_suffix)
                omic_specific_paths_list[[actual_var_name]] <- full_path
            }
        }
    }

    # Add specific variable names that might map to the more generic ones, for convenience or legacy.
    # This part also ensures specific directories (like transcriptomics count_data_dir) are explicitly listed if defined in results_subdirs.
    specific_mapped_names <- list()
    if (current_omic_type == "proteomics") {
        if (!is.null(omic_specific_paths_list$feature_qc_dir)) specific_mapped_names$protein_qc_dir <- omic_specific_paths_list$feature_qc_dir
        if (!is.null(omic_specific_paths_list$subfeature_qc_dir)) specific_mapped_names$peptide_qc_dir <- omic_specific_paths_list$subfeature_qc_dir
        if (!is.null(omic_specific_paths_list$clean_features_dir)) specific_mapped_names$clean_proteins_dir <- omic_specific_paths_list$clean_features_dir
    } else if (current_omic_type == "metabolomics") {
        if (!is.null(omic_specific_paths_list$feature_qc_dir)) specific_mapped_names$metabolite_qc_dir <- omic_specific_paths_list$feature_qc_dir
        if (!is.null(omic_specific_paths_list$clean_features_dir)) specific_mapped_names$clean_metabolites_dir <- omic_specific_paths_list$clean_features_dir
    } else if (current_omic_type == "transcriptomics") {
        if (!is.null(omic_specific_paths_list$feature_qc_dir)) specific_mapped_names$gene_qc_dir <- omic_specific_paths_list$feature_qc_dir
        # raw_counts_dir is already added from global_vars if raw_counts_leaf is defined
        if (!is.null(omic_specific_paths_list$clean_features_dir)) specific_mapped_names$normalized_counts_dir <- omic_specific_paths_list$clean_features_dir
    } else if (current_omic_type == "lipidomics") {
        if (!is.null(omic_specific_paths_list$feature_qc_dir)) specific_mapped_names$lipid_qc_dir <- omic_specific_paths_list$feature_qc_dir
        if (!is.null(omic_specific_paths_list$clean_features_dir)) specific_mapped_names$clean_lipids_dir <- omic_specific_paths_list$clean_features_dir
    }
    # For 'integration', specific paths like mofa_inputs_dir, etc., are already directly added from global_vars.

    omic_specific_paths_list <- c(omic_specific_paths_list, specific_mapped_names[!sapply(specific_mapped_names, is.null)])

    # Ensure all defined directories exist (belt-and-suspenders, esp. for reused structures or complex paths)
    # This primarily targets directories that are values in omic_specific_paths_list
    for (path_entry in omic_specific_paths_list) {
        if (is.character(path_entry) && (endsWith(path_entry, "_dir") || endsWith(path_entry, "_base") || grepl("mofa|mixomics|integration_enrichment", path_entry))) {
            # Avoid trying to create dir for base_dir, timestamp, omic_type, omic_label etc.
            # A bit heuristic here; might need refinement if path names are ambiguous.
            # Check if it's likely a directory path constructed within the project.
            if (startsWith(path_entry, base_dir) && path_entry != base_dir && !is.null(path_entry) && path_entry != "") {
                suppressWarnings(dir.create(path_entry, recursive = TRUE, showWarnings = FALSE))
            }
        }
    }

    omic_specific_paths_list
}

# ----------------------------------------------------------------------------
# handleSetupDirectoriesExistingDirs
# ----------------------------------------------------------------------------
handleSetupDirectoriesExistingDirs <- function(
    current_omic_type,
    omic_label_dirname,
    results_path,
    results_summary_path,
    scripts_path,
    force,
    reuse_existing
) {
    process_current_omic <- TRUE
    reuse_current_omic_dirs <- FALSE
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

    list(
        process_current_omic = process_current_omic,
        reuse_current_omic_dirs = reuse_current_omic_dirs
    )
}

# ----------------------------------------------------------------------------
# materializeSetupDirectoriesStructure
# ----------------------------------------------------------------------------
materializeSetupDirectoriesStructure <- function(
    current_omic_paths_def,
    omic_config,
    current_omic_type,
    omic_label_dirname,
    reuse_current_omic_dirs
) {
    if (!reuse_current_omic_dirs || (reuse_current_omic_dirs && !dir.exists(current_omic_paths_def$results_base))) {
        if (!reuse_current_omic_dirs) {
            logger::log_info("Creating directory structure for '{omic_label_dirname}'...")
        } else {
            logger::log_info("Ensuring directory structure exists for '{omic_label_dirname}' (Reuse Mode)...")
        }

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
        # Only copy scripts if NOT reusing (preserve existing scripts) or if forced?
        # For dynamic add, we might want to copy if missing.
        # Let's keep script copying strict to non-reuse mode to avoid overwriting user edits.
        if (!reuse_current_omic_dirs) {
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
        }
    } else {
        # If reusing AND base results dir exists, we assume structure is fine.
        # But maybe we should double check subdirs?
        # For robustness, let's just ensure subdirs exist even in reuse mode.

        logger::log_info("Reusing existing directory structure for '{omic_label_dirname}'. Checking subdirectories...")

        # Ensure subdirs exist
        dir.create(current_omic_paths_def$results_base, recursive = TRUE, showWarnings = FALSE)
        qualified_results_subdirs <- file.path(current_omic_paths_def$results_base, omic_config$results_subdirs)
        invisible(sapply(qualified_results_subdirs, dir.create, recursive = TRUE, showWarnings = FALSE))

        dir.create(current_omic_paths_def$results_summary_base, recursive = TRUE, showWarnings = FALSE)
        qualified_results_summary_subdirs <- file.path(current_omic_paths_def$results_summary_base, omic_config$results_summary_subdirs)
        invisible(sapply(qualified_results_summary_subdirs, dir.create, recursive = TRUE, showWarnings = FALSE))
    }
}

# ----------------------------------------------------------------------------
# printSetupDirectoriesSummary
# ----------------------------------------------------------------------------
printSetupDirectoriesSummary <- function(all_created_paths) {
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
                    is_dir_to_count <- endsWith(path_name, "_dir") || endsWith(path_name, "_base") || path_name %in% c("qc_dir", "time_dir")

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
}

# ----------------------------------------------------------------------------
# initializeSetupDirectoriesOmicPaths
# ----------------------------------------------------------------------------
initializeSetupDirectoriesOmicPaths <- function(
    base_dir,
    current_omic_type,
    results_path,
    results_summary_path,
    scripts_path,
    omic_config,
    timestamp
) {
    publication_graphs_dir_base <- file.path(results_path, "publication_graphs")
    # For omics types like 'integration', 'publication_graphs' might not be a
    # primary subdir, but the timestamped QC path still nests underneath it.
    qc_dir_base <- file.path(results_path, "publication_graphs", "filtering_qc")

    # Create the omic-specific data directory (for example data/proteomics).
    omic_specific_data_dir <- file.path(base_dir, "data", current_omic_type)
    dir.create(omic_specific_data_dir, recursive = TRUE, showWarnings = FALSE)
    logger::log_info("Ensured omic-specific data directory exists: {omic_specific_data_dir}")

    current_omic_paths_def <- list(
        results_base = results_path,
        results_summary_base = results_summary_path,
        data_dir = omic_specific_data_dir,
        scripts_dest_dir = scripts_path,
        scripts_source_dir = file.path(base_dir, "scripts", omic_config$scripts_source_leaf),
        publication_graphs_dir = publication_graphs_dir_base,
        qc_dir = qc_dir_base,
        time_dir = file.path(qc_dir_base, timestamp)
    )

    # Ensure the session timestamp directory exists before downstream
    # materialization references it.
    dir.create(current_omic_paths_def$time_dir, recursive = TRUE, showWarnings = FALSE)
    logger::log_info("Ensured timestamped directory exists: {current_omic_paths_def$time_dir}")

    current_omic_paths_def
}

