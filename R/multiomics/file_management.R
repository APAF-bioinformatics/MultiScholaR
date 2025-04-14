#' @description Creates and manages project directories with version control for a specific omic type.
#' @param base_dir Base directory path (optional, defaults to here::here())
#' @param omic_type Character string indicating the type of omics data (e.g., "proteomics", "metabolomics", "transcriptomics"). Defaults to "proteomics".
#' @param label Optional label to append to the omic type directory name (e.g., "proteomics_MyLabel").
#' @param force Logical; if TRUE, skips user confirmation (default: FALSE).
#' @return List of directory paths assigned to the global environment.
#' @export
multisetupAndShowDirectories <- function(base_dir = here::here(), omic_type = "proteomics", label = NULL, force = FALSE) {
    # Define Omics-Specific Configuration
    valid_omic_types <- c("proteomics", "metabolomics", "transcriptomics") # Added transcriptomics
    if (!omic_type %in% valid_omic_types) {
        stop("Invalid omic_type specified. Choose from: ", paste(valid_omic_types, collapse = ", "))
    }

    omic_config <- switch(omic_type,
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
                "metabolite_qc", "feature_qc", "de_metabolites",
                "publication_graphs", "pathway_enrichment",
                file.path("publication_graphs", "filtering_qc")
            ),
            results_summary_subdirs = c("QC_figures", "Publication_figures", "Publication_tables", "Study_report"),
            scripts_source_leaf = "metabolomics",
            global_vars = list(
                de_output_leaf = "de_metabolites",
                pathway_leaf = "pathway_enrichment",
                feature_qc_leaf = "metabolite_qc",
                subfeature_qc_leaf = NULL,
                clean_features_leaf = "clean_metabolites"
            )
        ),
        transcriptomics = list(
            results_subdirs = c(
                "gene_qc", "count_data", "de_genes", # Specific dirs for transcriptomics
                "publication_graphs", "pathway_enrichment",
                file.path("publication_graphs", "filtering_qc")
            ),
            results_summary_subdirs = c("QC_figures", "Publication_figures", "Publication_tables", "Study_report"),
            scripts_source_leaf = "transcriptomics", # Assumes scripts/transcriptomics exists or will exist
            global_vars = list(
                de_output_leaf = "de_genes",
                pathway_leaf = "pathway_enrichment",
                feature_qc_leaf = "gene_qc", # QC related to gene-level data
                subfeature_qc_leaf = NULL, # Or potentially "transcript_qc" if needed
                clean_features_leaf = "normalized_counts" # Or "clean_genes"
            )
        ),
        # Add other omics types here
        stop("Configuration not defined for omic_type: ", omic_type) # Fallback
    )

    # Construct Directory Names
    omic_label_dirname <- paste0(omic_type, if (!is.null(label)) paste0("_", substr(label, 1, 30)) else "")

    results_path <- file.path(base_dir, "results", omic_label_dirname)
    results_summary_path <- file.path(base_dir, "results_summary", omic_label_dirname)
    scripts_path <- file.path(base_dir, "scripts", omic_label_dirname) # Destination scripts path

    # Check Existing Directories and User Confirmation
    if (!force && (dir.exists(results_path) || dir.exists(results_summary_path) || dir.exists(scripts_path))) {
        cat(sprintf("\nWarning: Directory(ies) for '%s' already exist:\n", omic_label_dirname))
        if (dir.exists(results_path)) cat(sprintf("- %s\n", results_path))
        if (dir.exists(results_summary_path)) cat(sprintf("- %s\n", results_summary_path))
        if (dir.exists(scripts_path)) cat(sprintf("- %s\n", scripts_path))

        response <- readline(prompt = "Do you want to overwrite? (y/n): ")

        if (tolower(substr(response, 1, 1)) != "y") {
            message("Setup cancelled by user")
            return(invisible(NULL))
        }

        # Remove existing directories if user confirmed
        if (dir.exists(results_path)) unlink(results_path, recursive = TRUE)
        if (dir.exists(results_summary_path)) unlink(results_summary_path, recursive = TRUE)
        if (dir.exists(scripts_path)) unlink(scripts_path, recursive = TRUE)
    }

    # Define Paths
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    publication_graphs_dir_base <- file.path(results_path, "publication_graphs")
    qc_dir_base <- file.path(publication_graphs_dir_base, "filtering_qc")

    paths <- list(
        results_base = results_path,
        results_subdirs = file.path(results_path, omic_config$results_subdirs),
        results_summary_base = results_summary_path,
        results_summary_subdirs = file.path(results_summary_path, omic_config$results_summary_subdirs),
        data_dir = file.path(base_dir, "data"), # Common data directory
        scripts_dest_dir = scripts_path,
        scripts_source_dir = file.path(base_dir, "scripts", omic_config$scripts_source_leaf),
        publication_graphs_dir = publication_graphs_dir_base,
        qc_dir = qc_dir_base,
        time_dir = file.path(qc_dir_base, timestamp)
    )

    # Create Directories
    # Create results base and subdirs
    dir.create(paths$results_base, recursive = TRUE, showWarnings = FALSE)
    invisible(sapply(paths$results_subdirs, dir.create, recursive = TRUE, showWarnings = FALSE))

    # Create results_summary base and subdirs
    dir.create(paths$results_summary_base, recursive = TRUE, showWarnings = FALSE)
    invisible(sapply(paths$results_summary_subdirs, dir.create, recursive = TRUE, showWarnings = FALSE))

    # Create data and timestamped qc dir
    dir.create(paths$data_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$time_dir, recursive = TRUE, showWarnings = FALSE) # Includes qc_dir and pub_graphs parent

    # Handle Scripts Directory Copying
    if (dir.exists(paths$scripts_source_dir)) {
        dir.create(paths$scripts_dest_dir, recursive = TRUE, showWarnings = FALSE) # Create destination

        script_files <- list.files(paths$scripts_source_dir, full.names = TRUE, recursive = TRUE)
        script_files <- script_files[!grepl("\\.[rR][mM][dD]$", script_files)] # Filter out Rmd

        if (length(script_files) > 0) {
            copied_scripts <- sapply(script_files, function(f) {
                # Calculate relative path from source base to file f
                rel_path <- sub(paste0("^", paths$scripts_source_dir, "/?"), "", f)
                dest_file <- file.path(paths$scripts_dest_dir, rel_path)
                dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)
                file.copy(f, dest_file, overwrite = TRUE)
            })
            if (any(!copied_scripts)) warning("Some script files failed to copy.")
        }
        log_info("Copied scripts (excluding Rmd) from {paths$scripts_source_dir} to {paths$scripts_dest_dir}")
    } else {
        dir.create(paths$scripts_dest_dir, recursive = TRUE, showWarnings = FALSE) # Still create dest dir
        log_warn("Source script directory not found: {paths$scripts_source_dir}. No scripts copied.")
    }

    # Build and Assign Global Variables
    dir_paths_global <- list(
        base_dir = base_dir,
        results_dir = paths$results_base,
        data_dir = paths$data_dir,
        source_dir = paths$scripts_dest_dir, # Use the destination scripts dir for consistency
        de_output_dir = file.path(paths$results_base, omic_config$global_vars$de_output_leaf),
        publication_graphs_dir = paths$publication_graphs_dir,
        timestamp = timestamp,
        qc_dir = paths$qc_dir,
        time_dir = paths$time_dir,
        results_summary_dir = paths$results_summary_base,
        pathway_dir = file.path(paths$results_base, omic_config$global_vars$pathway_leaf),
        feature_qc_dir = file.path(paths$results_base, omic_config$global_vars$feature_qc_leaf),
        clean_features_dir = file.path(paths$results_base, omic_config$global_vars$clean_features_leaf),
        qc_figures_dir = file.path(paths$results_summary_base, "QC_figures"),
        publication_figures_dir = file.path(paths$results_summary_base, "Publication_figures"),
        publication_tables_dir = file.path(paths$results_summary_base, "Publication_tables"),
        study_report_dir = file.path(paths$results_summary_base, "Study_report")
    )

    # Add optional sub-feature QC dir and associated global variable if defined
    if (!is.null(omic_config$global_vars$subfeature_qc_leaf)) {
        subfeature_qc_dir_path <- file.path(paths$results_base, omic_config$global_vars$subfeature_qc_leaf)
        dir_paths_global$subfeature_qc_dir <- subfeature_qc_dir_path # Assign path to the list
        # Ensure the directory exists
        dir.create(subfeature_qc_dir_path, recursive = TRUE, showWarnings = FALSE)
    }

    # Dynamically add specific global variable names based on omic_type (e.g., protein_qc_dir)
    # This maps the generic 'feature_qc_leaf' etc. to the specific global variable name
    specific_global_names <- list(
        protein_qc_dir = if (omic_type == "proteomics") dir_paths_global$feature_qc_dir else NULL,
        peptide_qc_dir = if (omic_type == "proteomics") dir_paths_global$subfeature_qc_dir else NULL,
        clean_proteins_dir = if (omic_type == "proteomics") dir_paths_global$clean_features_dir else NULL,
        metabolite_qc_dir = if (omic_type == "metabolomics") dir_paths_global$feature_qc_dir else NULL,
        clean_metabolites_dir = if (omic_type == "metabolomics") dir_paths_global$clean_features_dir else NULL,
        gene_qc_dir = if (omic_type == "transcriptomics") dir_paths_global$feature_qc_dir else NULL,
        count_data_dir = if (omic_type == "transcriptomics") file.path(paths$results_base, "count_data") else NULL, # Add specific dir
        normalized_counts_dir = if (omic_type == "transcriptomics") dir_paths_global$clean_features_dir else NULL
        # Add other specific names as needed
    )

    # Add non-NULL specific names to the main list
    dir_paths_global <- c(dir_paths_global, specific_global_names[!sapply(specific_global_names, is.null)])

    # Assign to global environment
    list2env(dir_paths_global, envir = .GlobalEnv)
    log_info("Assigned directory paths to global environment.")

    # Print Structure
    cat("\nDirectory Structure Created/Verified:\n")
    print_paths <- dir_paths_global[sapply(dir_paths_global, is.character)] # Only print paths
    invisible(lapply(print_paths, function(p) {
        if (dir.exists(p)) {
            # Use tryCatch for permission issues, though list.files might be enough
            tryCatch(
                {
                    file_count <- length(list.files(p, recursive = TRUE, all.files = TRUE, no.. = TRUE))
                    cat(sprintf("  %s (%d files/dirs)\n", p, file_count))
                },
                error = function(e) {
                    cat(sprintf("  %s (Error accessing content)\n", p))
                }
            )
        } else {
            cat(sprintf("  %s (Path only, not a directory)\n", p)) # e.g., timestamp
        }
    }))

    invisible(dir_paths_global) # Return the list assigned to global env
}
