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

