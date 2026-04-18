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
    current_paths <- tryCatch(
        {
            getProjectPaths(
                omic_type = omic_type,
                experiment_label = experiment_label,
                project_dirs_object_name = project_dirs_object_name
            )
        },
        error = function(e) {
            rlang::abort(paste0("Failed to get project paths: ", e$message))
        }
    )

    # Validate that current_paths is a list and contains essential directory paths
    required_paths_in_current <- c(
        "results_dir", "results_summary_dir", "publication_graphs_dir",
        "time_dir", "qc_dir", "da_output_dir", "pathway_dir", "source_dir", "feature_qc_dir"
    )
    if (!is.list(current_paths) || !all(required_paths_in_current %in% names(current_paths))) {
        missing_req <- setdiff(required_paths_in_current, names(current_paths))
        rlang::abort(paste0("Essential paths missing from project_dirs: ", paste(missing_req, collapse = ", ")))
        rlang::abort(paste0("Essential paths missing from project_dirs: ", paste(missing_req, collapse = ", ")))
    }
    # --- End: Path Derivation and Validation ---

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
    cat(sprintf("DA Output Dir: %s\n", current_paths$da_output_dir))
    cat(sprintf("Pathway Dir: %s\n", current_paths$pathway_dir))
    cat(sprintf("Source (Scripts) Dir: %s\n", current_paths$source_dir))
    cat(sprintf("Feature QC Dir: %s\n", current_paths$feature_qc_dir))
    if (!is.null(current_paths$subfeature_qc_dir)) cat(sprintf("Sub-feature QC Dir: %s\n", current_paths$subfeature_qc_dir))
    cat("\n")

    # ROBUST: Try to get contrasts_tbl and design_matrix from environment first, then from files
    cat("Checking for required objects...\n")

    # Try contrasts_tbl from environment first
    if (is.null(contrasts_tbl)) {
        if (exists("contrasts_tbl", envir = parent.frame())) {
            contrasts_tbl <- get("contrasts_tbl", envir = parent.frame())
            cat("[OK] Using 'contrasts_tbl' from calling environment\n")
        } else if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
            cat("[OK] Using 'contrasts_tbl' from global environment\n")
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
                    tryCatch(
                        {
                            contrasts_tbl <- readr::read_tsv(contrasts_file, show_col_types = FALSE)
                            cat(sprintf("[OK] Loaded 'contrasts_tbl' from file: %s\n", basename(contrasts_file)))
                            contrasts_loaded <- TRUE
                            break
                        },
                        error = function(e) {
                            cat(sprintf("[WARNING] Failed to load contrasts from %s: %s\n", basename(contrasts_file), e$message))
                        }
                    )
                }
            }

            if (!contrasts_loaded) {
                cat("[FAIL] 'contrasts_tbl' not found in environment or files\n")
            }
        }
    } else {
        cat("[OK] Using provided 'contrasts_tbl' parameter\n")
    }

    # Try design_matrix from environment first
    if (is.null(design_matrix)) {
        if (exists("design_matrix", envir = parent.frame())) {
            design_matrix <- get("design_matrix", envir = parent.frame())
            cat("[OK] Using 'design_matrix' from calling environment\n")
        } else if (exists("design_matrix", envir = .GlobalEnv)) {
            design_matrix <- get("design_matrix", envir = .GlobalEnv)
            cat("[OK] Using 'design_matrix' from global environment\n")
        } else {
            # Fallback: try to load from file
            design_matrix_file <- file.path(current_paths$source_dir, "design_matrix.tab")

            if (file.exists(design_matrix_file)) {
                tryCatch(
                    {
                        design_matrix <- readr::read_tsv(design_matrix_file, show_col_types = FALSE)
                        cat(sprintf("[OK] Loaded 'design_matrix' from file: %s\n", basename(design_matrix_file)))
                    },
                    error = function(e) {
                        cat(sprintf("[FAIL] Failed to load design_matrix from %s: %s\n", basename(design_matrix_file), e$message))
                        design_matrix <- NULL
                    }
                )
            } else {
                cat(sprintf("[FAIL] 'design_matrix' not found in environment or at expected file location: %s\n", design_matrix_file))
            }
        }
    } else {
        cat("[OK] Using provided 'design_matrix' parameter\n")
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
        tryCatch(
            {
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
            },
            error = function(e) {
                warning(sprintf("Failed to save/copy Rmd file: %s", e$message))
                failed_copies[[length(failed_copies) + 1]] <- list(type = "rmd_copy", source = current_rmd, destination = current_paths$source_dir, display_name = "Current Rmd File", error = e$message)
            }
        )
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
            return(invisible(list(status = "cancelled", omic_key = omic_label, message = paste0("Backup and overwrite for ", omic_label, " cancelled by user."))))
            return(invisible(list(status = "cancelled", omic_key = omic_label, message = paste0("Backup and overwrite for ", omic_label, " cancelled by user."))))
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

            backup_has_content <- length(list.files(backup_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)) > 0

            if (backup_copy_successful && (backup_has_content || length(contents_of_summary_dir) == 0)) {
                logger::log_info("Successfully backed up content of {current_paths$results_summary_dir} to: {backup_dir}")
                backup_info <- data.frame(original_dir = current_paths$results_summary_dir, backup_time = Sys.time(), omic_key = omic_label, stringsAsFactors = FALSE)
                backup_info <- data.frame(original_dir = current_paths$results_summary_dir, backup_time = Sys.time(), omic_key = omic_label, stringsAsFactors = FALSE)
                tryCatch(
                    write.table(backup_info, file = file.path(backup_dir, "backup_info.txt"), sep = "\\t", row.names = FALSE, quote = FALSE),
                    error = function(e) logger::log_warn("Failed to write backup_info.txt: {e$message}")
                )

                logger::log_info("Clearing original results_summary_dir: {current_paths$results_summary_dir}")
                unlink_success <- tryCatch(
                    {
                        unlink(current_paths$results_summary_dir, recursive = TRUE, force = TRUE)
                        TRUE # Return TRUE on success
                    },
                    error = function(e) {
                        logger::log_error("Error unlinking {current_paths$results_summary_dir}: {e$message}")
                        FALSE # Return FALSE on error
                    }
                )

                if (unlink_success) {
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
            if (!dir.create(current_paths$results_summary_dir, recursive = TRUE, showWarnings = FALSE)) {
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
        list(source = file.path(current_paths$publication_graphs_dir, "NumSigDaMolecules"), dest = "Publication_figures/NumSigDaMolecules", is_dir = TRUE, display_name = "Num Sig DA Molecules"),
        list(source = file.path(current_paths$publication_graphs_dir, "Volcano_Plots"), dest = "Publication_figures/Volcano_Plots", is_dir = TRUE, display_name = "Volcano Plots"),
        list(source = file.path(current_paths$publication_graphs_dir, "Heatmap"), dest = "Publication_figures/Heatmap", is_dir = TRUE, display_name = "Interactive Heatmaps"),
        list(source = current_paths$pathway_dir, dest = "Publication_figures/Enrichment_Plots", is_dir = TRUE, display_name = "Pathway Enrichment Plots"),
        list(source = "contrasts_tbl", dest = "Study_report", type = "object", save_as = "contrasts_tbl.tab", display_name = "Contrasts Table"),
        list(source = "design_matrix", dest = "Study_report", type = "object", save_as = "design_matrix.tab", display_name = "Design Matrix"),
        list(source = file.path(current_paths$source_dir, "study_parameters.txt"), dest = "Study_report", is_dir = FALSE, display_name = "Study Parameters")
    )

    # Logic to find and copy the num_sig_da_molecules tab file to Publication_tables
    num_sig_dir <- file.path(current_paths$publication_graphs_dir, "NumSigDaMolecules")
    if (dir.exists(num_sig_dir)) {
        # Look for the tab file
        sig_files <- list.files(num_sig_dir, pattern = "(_num_sig_da_molecules|_num_significant_differentially_abundant_all)\\.tab$", full.names = TRUE)
        if (length(sig_files) > 0) {
            # Use the first match
            files_to_copy <- c(files_to_copy, list(
                list(source = sig_files[1], dest = "Publication_tables", is_dir = FALSE, display_name = "Num Sig DA Molecules Tab", new_name = paste0("da_", omic_type, "_num_sig_da_molecules.tab"))
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
    } else if (omic_type == "metabolomics") {
        cat("COPY: Adding metabolomics-specific files\n")

        # Composite QC figure (generated by mod_metab_norm.R)
        composite_file <- file.path(current_paths$feature_qc_dir, "composite_QC_figure.png")
        if (file.exists(composite_file)) {
            files_to_copy <- c(files_to_copy, list(
                list(source = composite_file, dest = "QC_figures", is_dir = FALSE, display_name = "Composite QC (PNG)", new_name = "metabolomics_composite_QC_figure.png")
            ))
            cat("COPY: Added metabolomics composite QC figure\n")
        } else {
            cat("COPY: Warning - Composite QC figure not found at", composite_file, "\n")
        }

        # Per-assay stage-based QC plots (pattern: *_pre_norm_*.png, *_post_norm_*.png, *_ruv_corrected_*.png)
        # Matches: lcms_neg_pre_norm_pca.png, lcms_pos_post_norm_density.png, etc.
        stage_qc_files <- list.files(current_paths$feature_qc_dir, pattern = "_(pre_norm|post_norm|ruv_corrected)_(pca|density|rle|correlation)\\.png$", full.names = TRUE)
        if (length(stage_qc_files) > 0) {
            for (qc_file in stage_qc_files) {
                files_to_copy <- c(files_to_copy, list(
                    list(source = qc_file, dest = "QC_figures", is_dir = FALSE, display_name = sprintf("Metab QC: %s", basename(qc_file)))
                ))
            }
            cat(sprintf("COPY: Added %d metabolomics stage-based QC files\n", length(stage_qc_files)))
        }

        # Legacy patterns: Per-assay QC composites (pattern: *_assay_metrics.png, *_combined_plots.png)
        metab_qc_files <- list.files(current_paths$feature_qc_dir, pattern = "(assay_metrics|combined_plots)\\.png$", full.names = TRUE)
        if (length(metab_qc_files) > 0) {
            for (qc_file in metab_qc_files) {
                files_to_copy <- c(files_to_copy, list(
                    list(source = qc_file, dest = "QC_figures", is_dir = FALSE, display_name = sprintf("Metab QC: %s", basename(qc_file)))
                ))
            }
            cat(sprintf("COPY: Added %d metabolomics legacy QC files\n", length(metab_qc_files)))
        }

        # ITSD normalization figures (pattern: itsd_*.png)
        itsd_files <- list.files(current_paths$feature_qc_dir, pattern = "itsd.*\\.png$", full.names = TRUE)
        if (length(itsd_files) > 0) {
            for (itsd_file in itsd_files) {
                files_to_copy <- c(files_to_copy, list(
                    list(source = itsd_file, dest = "QC_figures", is_dir = FALSE, display_name = sprintf("ITSD: %s", basename(itsd_file)))
                ))
            }
            cat(sprintf("COPY: Added %d ITSD normalization figures\n", length(itsd_files)))
        }

        # Correlation heatmaps per assay (legacy pattern for any correlation files not caught above)
        corr_files <- list.files(current_paths$feature_qc_dir, pattern = "^[^_]*_?correlation[^/]*\\.png$", full.names = TRUE)
        # Filter out files already added by stage_qc_files pattern
        corr_files <- setdiff(corr_files, stage_qc_files)
        if (length(corr_files) > 0) {
            for (corr_file in corr_files) {
                files_to_copy <- c(files_to_copy, list(
                    list(source = corr_file, dest = "QC_figures", is_dir = FALSE, display_name = sprintf("Correlation: %s", basename(corr_file)))
                ))
            }
            cat(sprintf("COPY: Added %d additional correlation heatmap files\n", length(corr_files)))
        }

        # Heatmap directory for DA heatmaps
        heatmap_dir <- file.path(current_paths$publication_graphs_dir, "Heatmaps")
        if (dir.exists(heatmap_dir)) {
            files_to_copy <- c(files_to_copy, list(
                list(source = heatmap_dir, dest = "Publication_figures/Heatmaps", is_dir = TRUE, display_name = "DA Heatmaps")
            ))
            cat("COPY: Added Heatmaps directory\n")
        }

        # Normalized results (check for RUV vs non-RUV like proteomics)
        ruv_rds_path <- file.path(current_paths$feature_qc_dir, "ruv_normalised_results.RDS")
        norm_rds_path <- file.path(current_paths$feature_qc_dir, "normalised_results.RDS")

        if (file.exists(ruv_rds_path)) {
            files_to_copy <- c(files_to_copy, list(
                list(source = ruv_rds_path, dest = "Publication_tables", is_dir = FALSE, display_name = "RUV Normalized Results RDS (Metab)", new_name = "ruv_normalised_results.RDS")
            ))
            cat("COPY: Using RUV-normalized file for metabolomics\n")
        } else if (file.exists(norm_rds_path)) {
            files_to_copy <- c(files_to_copy, list(
                list(source = norm_rds_path, dest = "Publication_tables", is_dir = FALSE, display_name = "Normalized Results RDS (Metab)", new_name = "normalised_results.RDS")
            ))
            cat("COPY: Using normalized file for metabolomics (RUV skipped)\n")
        } else {
            cat("COPY: Warning - no normalized RDS files found for metabolomics\n")
        }
    }

    # Excel files paths
    da_results_excel_path <- file.path(pub_tables_dir, paste0("DA_results_", omic_type, ".xlsx"))
    enrichment_excel_path <- file.path(pub_tables_dir, paste0("Pathway_enrichment_results_", omic_type, ".xlsx"))

    # Create combined DA workbook
    da_wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(da_wb, "DA_Results_Index")
    da_index_data <- data.frame(Sheet = character(), Description = character(), stringsAsFactors = FALSE)
    da_files <- list.files(path = current_paths$da_output_dir, pattern = paste0("(da|de)_.+_long_annot\\.xlsx$"), full.names = TRUE)

    purrr::imap(da_files, \(file, idx) {
        sheet_name <- sprintf("DA_Sheet%d", idx)
        base_name <- basename(file) |>
            stringr::str_remove("^da_") |>
            stringr::str_remove("^de_") |>
            stringr::str_remove("_long_annot\\.xlsx$")
        da_index_data <<- rbind(da_index_data, data.frame(Sheet = sheet_name, Description = base_name, stringsAsFactors = FALSE))
        data_content <- tryCatch(openxlsx::read.xlsx(file), error = function(e) NULL)
        if (!is.null(data_content)) {
            openxlsx::addWorksheet(da_wb, sheet_name)
            openxlsx::writeData(da_wb, sheet_name, data_content)
        } else {
            warning(paste0("Failed to read DA Excel file: ", file))
            failed_copies[[length(failed_copies) + 1]] <- list(type = "da_excel_read", source = file, error = "Failed to read")
        }
    })
    openxlsx::writeData(da_wb, "DA_Results_Index", da_index_data)
    openxlsx::setColWidths(da_wb, "DA_Results_Index", cols = 1:2, widths = c(15, 50))
    openxlsx::addStyle(da_wb, "DA_Results_Index", style = openxlsx::createStyle(textDecoration = "bold"), rows = 1, cols = 1:2)

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
    tryCatch(
        {
            openxlsx::saveWorkbook(da_wb, da_results_excel_path, overwrite = TRUE)
            cat(paste("Successfully saved DA results to:", da_results_excel_path, "\n"))
        },
        error = function(e) {
            failed_copies[[length(failed_copies) + 1]] <- list(type = "workbook_save", path = da_results_excel_path, error = e$message)
        }
    )

    tryCatch(
        {
            openxlsx::saveWorkbook(enrichment_wb, enrichment_excel_path, overwrite = TRUE)
            cat(paste("Successfully saved Enrichment results to:", enrichment_excel_path, "\n"))
        },
        error = function(e) {
            failed_copies[[length(failed_copies) + 1]] <- list(type = "workbook_save", path = enrichment_excel_path, error = e$message)
        }
    )

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
                    cat(sprintf("    [OK] Found '%s' in calling environment\n", file_spec$source))
                } else if (!is.null(get0(file_spec$source, envir = .GlobalEnv))) {
                    obj <- get(file_spec$source, envir = .GlobalEnv)
                    source_exists <- TRUE
                    cat(sprintf("    [OK] Found '%s' in global environment\n", file_spec$source))
                } else {
                    # Fallback: try to load from file
                    if (file_spec$source == "design_matrix") {
                        design_matrix_file <- file.path(current_paths$source_dir, "design_matrix.tab")
                        if (file.exists(design_matrix_file)) {
                            tryCatch(
                                {
                                    obj <- readr::read_tsv(design_matrix_file, show_col_types = FALSE)
                                    source_exists <- TRUE
                                    cat(sprintf("    [OK] Loaded '%s' from file: %s\n", file_spec$source, basename(design_matrix_file)))
                                },
                                error = function(e) {
                                    error_msg <- sprintf("Failed to load %s from file %s: %s", file_spec$source, basename(design_matrix_file), e$message)
                                    cat(sprintf("    [FAIL] %s\n", error_msg))
                                }
                            )
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
                                tryCatch(
                                    {
                                        obj <- readr::read_tsv(contrasts_file, show_col_types = FALSE)
                                        source_exists <- TRUE
                                        cat(sprintf("    [OK] Loaded '%s' from file: %s\n", file_spec$source, basename(contrasts_file)))
                                        break
                                    },
                                    error = function(e) {
                                        cat(sprintf("    [WARNING] Failed to load contrasts from %s: %s\n", basename(contrasts_file), e$message))
                                    }
                                )
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
                if (!source_exists) error_msg <- sprintf("Source %s not found: %s", if (file_spec$is_dir) "directory" else "file", file_spec$source)
            }

            if (source_exists) {
                dir.create(dest_dir_final, recursive = TRUE, showWarnings = FALSE)
                if (!is.null(file_spec$type) && file_spec$type == "object") {
                    tryCatch(
                        {
                            # Use the already-loaded object instead of getting from environment again
                            if (is.null(obj)) {
                                obj <- get(file_spec$source, envir = parent.frame()) # Fallback for other objects
                            }
                            dest_path <- file.path(dest_dir_final, file_spec$save_as)
                            write.table(obj, file = dest_path, sep = "\t", row.names = FALSE, quote = FALSE)
                            if (!file.exists(dest_path) || (file.exists(dest_path) && file.size(dest_path) == 0 && nrow(obj) > 0)) {
                                copy_success <- FALSE
                                error_msg <- "Failed to write object or file is empty"
                            }
                        },
                        error = function(e) {
                            copy_success <<- FALSE
                            error_msg <<- sprintf("Error writing object: %s", e$message)
                        }
                    )
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
                failed_copies[[length(failed_copies) + 1]] <- list(type = if (!is.null(file_spec$type) && file_spec$type == "object") "object" else if (file_spec$is_dir) "directory" else "file", source = source_display, destination = dest_dir_final, display_name = file_spec$display_name, error = error_msg)
            }
            cat(sprintf("%-35s [%s -> %s] %s\n", file_spec$display_name, if (source_exists) "[OK]" else "[FAIL]", if (copy_success && source_exists) "[OK]" else "[FAIL]", if (!is.null(file_spec$type) && file_spec$type == "object") "Object" else if (file_spec$is_dir) "Directory" else "File"))
            if (!is.null(error_msg)) cat(sprintf("%35s Error: %s\n", "", error_msg))
        })

    cat("\nLegend: [OK] = exists/success, [FAIL] = missing/failed\n")
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
    tryCatch(
        {
            download.file(github_url, cached_file, mode = "wb", quiet = TRUE)
            logger::log_info("Successfully downloaded template to: {cached_file}")
            return(cached_file)
        },
        error = function(e) {
            rlang::abort(paste0(
                "Failed to download template '", rmd_filename, "' from GitHub.\n",
                "URL: ", github_url, "\n",
                "Error: ", e$message
            ))
        }
    )
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
    current_paths <- tryCatch(
        {
            getProjectPaths(
                omic_type = omic_type,
                experiment_label = experiment_label,
                project_dirs_object_name = project_dirs_object_name
            )
        },
        error = function(e) {
            rlang::abort(paste0("Failed to get project paths: ", e$message))
        }
    )

    message("   DEBUG66: getProjectPaths completed")
    if (!is.null(current_paths)) {
        message(sprintf("      DEBUG66: current_paths is list: %s", is.list(current_paths)))
        message(sprintf("      DEBUG66: current_paths names: %s", paste(names(current_paths), collapse = ", ")))
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
        rlang::abort(paste0(
            "R Markdown template file not found at the expected location: ", sQuote(rmd_input_path),
            ". This should be in the general scripts/<omic_type> directory (e.g., scripts/proteomics)."
        ))
    }

    # --- Construct Output Path (in the labelled results_summary directory) ---
    message("   DEBUG66: Constructing output file path...")
    output_file_basename <- paste0(
        tools::file_path_sans_ext(rmd_filename),
        "_", omic_type,
        "_", experiment_label
    )

    output_ext <- ".docx" # Default
    if (!is.null(output_format)) {
        if (output_format == "word_document" || grepl("word", output_format, ignore.case = TRUE)) {
            output_ext <- ".docx"
        } else if (output_format == "html_document" || grepl("html", output_format, ignore.case = TRUE)) {
            output_ext <- ".html"
        } else if (grepl("pdf", output_format, ignore.case = TRUE)) {
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

    rendered_path <- tryCatch(
        {
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
        },
        error = function(e) {
            message("   DEBUG66: rmarkdown::render FAILED")
            message(sprintf("      DEBUG66: Error message: %s", e$message))
            message("      DEBUG66: Error traceback:")
            print(e)
            logger::log_error("Failed to render R Markdown report: {e$message}")
            logger::log_error("Input path: {rmd_input_path}")
            logger::log_error("Output path: {output_file_path}")
            NULL # Return NULL on failure
        }
    )

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

# ----------------------------------------------------------------------------
# saveListOfPdfs
# ----------------------------------------------------------------------------
#' @export
saveListOfPdfs <- function(list, filename) {
    # start pdf
    cairo_pdf(filename)

    # loop
    # purrr::walk( list, print)
    for (p in list) {
        print(p)
    }

    # end pdf
    dev.off()

    invisible(NULL)
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
    temp <- tempfile(fileext = ".R")
    knitr::purl(file, output = temp)

    if (skip_plots) {
        old_dev <- getOption("device")
        options(device = function(...) {
            .Call("R_GD_nullDevice", PACKAGE = "grDevices")
        })
    }
    source(temp)
    if (skip_plots) {
        options(device = old_dev)
    }
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
savePlot <- function(plot, base_path, plot_name, formats = c("pdf", "png"), width = 7, height = 7, ...) {
    # Always save the RDS (works for both single plots and lists)
    saveRDS(plot, file.path(base_path, paste0(plot_name, ".rds")))

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
                    
                    # Use cairo_pdf for PDF format to avoid font issues on macOS
                    save_device <- format
                    if (format == "pdf") {
                        save_device <- grDevices::cairo_pdf
                    }
                    
                    ggsave(filename = file_path, plot = p, device = save_device, width = width, height = height, ...)
                })
            }
        })
    } else {
        # Single plot - original behavior
        purrr::walk(formats, \(format){
            file_path <- file.path(base_path, paste0(plot_name, ".", format))
            
            # Use cairo_pdf for PDF format to avoid font issues on macOS
            save_device <- format
            if (format == "pdf") {
                save_device <- grDevices::cairo_pdf
            }
            
            ggsave(filename = file_path, plot = plot, device = save_device, width = width, height = height, ...)
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
save_plot <- function(plot, base_path, plot_name, formats = c("pdf", "png"), width = 7, height = 7, ...) {
    savePlot(plot, base_path, plot_name, formats, width, height, ...)
}

# ----------------------------------------------------------------------------
# write_results
# ----------------------------------------------------------------------------
write_results <- function(data, filename) {
    vroom::vroom_write(data, file.path(results_dir, "protein_qc", filename))
}

