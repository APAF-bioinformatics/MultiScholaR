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
            subdirs = c(
                "protein_qc", "peptide_qc", "clean_proteins", "da_proteins",
                "publication_graphs", "pathway_enrichment",
                file.path("publication_graphs", "filtering_qc")
            )
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
                dir.create,
                recursive = TRUE, showWarnings = FALSE
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
        da_output_dir = file.path(paths$results$base, "da_proteins"),
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

