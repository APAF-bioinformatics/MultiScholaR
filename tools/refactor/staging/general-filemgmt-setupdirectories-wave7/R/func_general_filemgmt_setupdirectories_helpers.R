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

