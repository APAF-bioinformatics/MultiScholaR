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

