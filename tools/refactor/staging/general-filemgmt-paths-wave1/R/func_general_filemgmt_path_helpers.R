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
        rlang::abort(paste0(
            "Global object ", sQuote(project_dirs_object_name),
            " not found. Run setupDirectories() first or ensure project is properly initialized."
        ))
    }

    project_dirs_global <- get(project_dirs_object_name, envir = env)

    message("   DEBUG66: project_dirs object found in global environment")
    message(sprintf("      DEBUG66: project_dirs has %d keys", length(names(project_dirs_global))))
    message(sprintf("      DEBUG66: Available keys: %s", paste(names(project_dirs_global), collapse = ", ")))

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
        potential_keys <- c(potential_keys, matching_keys[1]) # Take first match
    }

    message(sprintf("   DEBUG66: Trying %d potential keys: %s", length(potential_keys), paste(potential_keys, collapse = ", ")))

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
    if (!file.exists(file_path)) {
        dir.create(file_path, showWarnings = TRUE, recursive = TRUE, mode = mode)
    }
}

# ----------------------------------------------------------------------------
# createDirIfNotExists
# ----------------------------------------------------------------------------
#' @export
createDirIfNotExists <- function(file_path, mode = "0777") {
    createDirectoryIfNotExists(file_path, mode = "0777")
}

