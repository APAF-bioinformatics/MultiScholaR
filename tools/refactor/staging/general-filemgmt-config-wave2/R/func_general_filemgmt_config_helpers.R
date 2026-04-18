#' @title Format Configuration List
#' @param config_list List of configuration parameters
#' @param indent Number of spaces for indentation
#' @export
formatConfigList <- function(config_list, indent = 0) {
    message(sprintf("--- DEBUG66: Entering formatConfigList with %d items, indent=%d ---", length(config_list), indent))
    output <- character()

    # Exclude internal_workflow_source_dir from printing
    names_to_process <- names(config_list)
    names_to_process <- names_to_process[names_to_process != "internal_workflow_source_dir"]
    message(sprintf("   DEBUG66: Processing %d items after exclusions", length(names_to_process)))

    # FUNCTIONAL APPROACH - NO FOR LOOPS that cause hanging - Works in Shiny AND .rmd
    output <- if (requireNamespace("purrr", quietly = TRUE)) {
        purrr::map_chr(names_to_process, function(name) {
            message(sprintf("   DEBUG66: Processing config item '%s'", name))
            value <- config_list[[name]]
            message(sprintf("   DEBUG66: Item '%s' class: %s", name, paste(class(value), collapse = ", ")))

            # Skip core_utilisation, seqinr_obj, and complex objects from display
            if (name == "core_utilisation" || name == "seqinr_obj" ||
                any(class(value) %in% c("process", "R6", "multidplyr_cluster", "cluster", "SOCKcluster", "tbl_df", "tbl", "data.frame"))) {
                message(sprintf("   DEBUG66: Skipping '%s' due to complex class or large data frame", name))
                return("") # Return empty string instead of next
            }

            # Format the name
            name_formatted <- gsub("\\.", " ", name)
            name_formatted <- gsub("_", " ", name_formatted)
            name_formatted <- tools::toTitleCase(name_formatted)

            # Handle different value types
            if (is.list(value)) {
                message(sprintf("   DEBUG66: '%s' is a list with %d elements", name, length(value)))
                if (length(value) > 0 && !is.null(names(value))) {
                    output_lines <- paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ":")
                    message(sprintf("   DEBUG66: Recursing into list '%s'", name))
                    recursive_result <- formatConfigList(value, indent + 2)
                    return(paste(c(output_lines, recursive_result), collapse = "\n"))
                } else if (length(value) > 0) { # Unnamed list, process elements functionally
                    message(sprintf("   DEBUG66: '%s' is unnamed list, processing elements", name))
                    header <- paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ":")

                    # FUNCTIONAL processing of list elements - NO FOR LOOP
                    element_lines <- if (requireNamespace("purrr", quietly = TRUE)) {
                        purrr::imap_chr(value, function(item_val, item_idx) {
                            message(sprintf("   DEBUG66: Processing unnamed list item %d, class: %s", item_idx, paste(class(item_val), collapse = ", ")))
                            if (is.atomic(item_val) && length(item_val) == 1) {
                                paste0(paste(rep(" ", indent + 2), collapse = ""), "- ", as.character(item_val))
                            } else {
                                paste0(paste(rep(" ", indent + 2), collapse = ""), "- [Complex List Element]")
                            }
                        })
                    } else {
                        # Base R fallback for list processing
                        sapply(seq_along(value), function(item_idx) {
                            item_val <- value[[item_idx]]
                            if (is.atomic(item_val) && length(item_val) == 1) {
                                paste0(paste(rep(" ", indent + 2), collapse = ""), "- ", as.character(item_val))
                            } else {
                                paste0(paste(rep(" ", indent + 2), collapse = ""), "- [Complex List Element]")
                            }
                        })
                    }
                    return(paste(c(header, element_lines), collapse = "\n"))
                } else { # Empty list
                    message(sprintf("   DEBUG66: '%s' is empty list", name))
                    return(paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ": [Empty List]"))
                }
            } else {
                message(sprintf("   DEBUG66: '%s' is atomic/non-list, attempting conversion", name))

                # SAFE value display formatting
                value_display <- tryCatch(
                    {
                        if (is.character(value) && length(value) > 5) {
                            paste(c(utils::head(value, 5), "..."), collapse = ", ")
                        } else if (is.character(value) && length(value) > 1) {
                            paste(value, collapse = ", ")
                        } else {
                            message(sprintf("   DEBUG66: About to convert '%s' to character", name))
                            result <- as.character(value)
                            message(sprintf("   DEBUG66: Conversion successful for '%s'", name))
                            if (length(result) == 0) "[Empty/NULL]" else result
                        }
                    },
                    error = function(e) {
                        message(sprintf("   DEBUG66: Error converting '%s': %s", name, e$message))
                        "[CONVERSION ERROR]"
                    }
                )

                return(paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ": ", value_display))
            }
        }) |>
            (\(x) x[x != ""])() # Remove empty strings from skipped items
    } else {
        # Fallback using base R lapply if purrr not available
        result_list <- lapply(names_to_process, function(name) {
            value <- config_list[[name]]

            # Skip complex objects
            if (name == "core_utilisation" ||
                any(class(value) %in% c("process", "R6", "multidplyr_cluster", "cluster", "SOCKcluster"))) {
                return("")
            }

            # Format name
            name_formatted <- gsub("\\.", " ", name)
            name_formatted <- gsub("_", " ", name_formatted)
            name_formatted <- tools::toTitleCase(name_formatted)

            # Simple formatting for fallback
            value_display <- tryCatch(
                {
                    if (is.list(value)) {
                        "[List Object]"
                    } else {
                        as.character(value)[1] # Take first element only for safety
                    }
                },
                error = function(e) "[CONVERSION ERROR]"
            )

            paste0(paste(rep(" ", indent), collapse = ""), name_formatted, ": ", value_display)
        })
        unlist(result_list[result_list != ""])
    }

    # Flatten the nested strings
    output <- unlist(strsplit(output, "\n"))
    message(sprintf("--- DEBUG66: Exiting formatConfigList, returning %d lines ---", length(output)))
    return(output)
}

#' @title Update Parameter in S4 Object Args and Global Config List
#' @description Modifies a specific parameter within an S4 object's @args slot
#'              and also updates the corresponding value in a global list named
#'              'config_list'.
#'
#' @param theObject The S4 object whose @args slot needs updating.
#' @param function_name The name identifying the parameter section (character string,
#'                      e.g., "peptideIntensityFiltering"). Corresponds to the
#'                      first-level key in both @args and config_list.
#' @param parameter_name The specific parameter name to update (character string,
#'                       e.g., "peptides_proportion_of_samples_below_cutoff").
#'                       Corresponds to the second-level key.
#' @param new_value The new value to assign to the parameter.
#' @param config_list_name The name of the global list variable holding the
#'                         configuration (defaults to "config_list").
#' @param env The environment where the global config list resides (defaults to
#'            .GlobalEnv).
#'
#' @return The modified S4 object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assume 'myPeptideData' is a PeptideQuantitativeData object
#' # Assume 'config_list' exists in the global environment
#'
#' # Check initial values (example)
#' # print(myPeptideData@args$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff)
#' # print(config_list$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff)
#'
#' # Update the parameter to 0.7
#' myPeptideData <- updateConfigParameter(
#'     theObject = myPeptideData,
#'     function_name = "peptideIntensityFiltering",
#'     parameter_name = "peptides_proportion_of_samples_below_cutoff",
#'     new_value = 0.7
#' )
#'
#' # Verify changes (example)
#' # print(myPeptideData@args$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff) # Should be 0.7
#' # print(config_list$peptideIntensityFiltering$peptides_proportion_of_samples_below_cutoff) # Should be 0.7
#' }
updateConfigParameter <- function(theObject,
                                  function_name,
                                  parameter_name,
                                  new_value,
                                  config_list_name = "config_list",
                                  env = .GlobalEnv) {
    # --- Input Validation ---
    if (!isS4(theObject)) {
        stop("'theObject' must be an S4 object.")
    }
    if (!"args" %in% methods::slotNames(theObject)) {
        stop("'theObject' must have an '@args' slot.")
    }
    if (!is.character(function_name) || length(function_name) != 1) {
        stop("'function_name' must be a single character string.")
    }
    if (!is.character(parameter_name) || length(parameter_name) != 1) {
        stop("'parameter_name' must be a single character string.")
    }
    if (!exists(config_list_name, envir = env)) {
        stop("Global config list '", config_list_name, "' not found in the specified environment.")
    }

    # Retrieve the global list safely
    current_config_list <- get(config_list_name, envir = env)

    if (!is.list(current_config_list)) {
        stop("Global variable '", config_list_name, "' is not a list.")
    }
    if (!function_name %in% names(current_config_list)) {
        warning("Function name '", function_name, "' not found in global config list '", config_list_name, "'. Adding it.")
        current_config_list[[function_name]] <- list()
    }
    if (!parameter_name %in% names(current_config_list[[function_name]])) {
        warning("Parameter '", parameter_name, "' not found under '", function_name, "' in global config list '", config_list_name, "'. Adding it.")
    }


    # --- Update S4 Object @args ---
    if (is.null(theObject@args)) {
        theObject@args <- list() # Initialize args if it's NULL
    }
    if (!is.list(theObject@args[[function_name]])) {
        # Initialize the sub-list if it doesn't exist or isn't a list
        theObject@args[[function_name]] <- list()
    }
    theObject@args[[function_name]][[parameter_name]] <- new_value
    message("Updated @args$", function_name, "$", parameter_name, " in S4 object.")


    # --- Update Global Config List ---
    current_config_list[[function_name]][[parameter_name]] <- new_value
    # Assign the modified list back to the global environment
    assign(config_list_name, current_config_list, envir = env)
    message("Updated ", config_list_name, "$", function_name, "$", parameter_name, " in global environment.")

    # --- Return the modified object ---
    return(theObject)
}

#' Create Study Parameters File
#'
#' @description
#' Creates a study parameters text file directly without using S4 objects.
#' This replaces the overly complex createWorkflowArgsFromConfig + WorkflowArgs show() approach.
#'
#' @param workflow_name Character string, name of the workflow
#' @param description Character string, description of the analysis
#' @param organism_name Character string, organism name (optional)
#' @param taxon_id Character string or numeric, taxon ID (optional)
#' @param source_dir_path Character string, path to save the study_parameters.txt file
#' @param contrasts_tbl Data frame, contrasts table (optional)
#' @param config_list_name Character string, name of the global config list (default: "config_list")
#' @param env Environment where the config list resides (default: .GlobalEnv)
#'
#' @return Character string path to the created study_parameters.txt file
#'
#' @examples
#' \dontrun{
#' file_path <- createStudyParametersFile(
#'     workflow_name = "proteomics_analysis",
#'     description = "DIA proteomics analysis",
#'     source_dir_path = "/path/to/scripts"
#' )
#' }
#'
#' @export
createStudyParametersFile <- function(workflow_name,
                                      description = "",
                                      organism_name = NULL,
                                      taxon_id = NULL,
                                      source_dir_path = NULL,
                                      contrasts_tbl = NULL,
                                      config_list_name = "config_list",
                                      env = .GlobalEnv) {
    # Validate required inputs
    if (missing(workflow_name) || !is.character(workflow_name) || length(workflow_name) != 1) {
        stop("workflow_name must be a single character string")
    }

    if (is.null(source_dir_path) || !is.character(source_dir_path) || length(source_dir_path) != 1) {
        stop("source_dir_path must be a single character string")
    }

    if (!dir.exists(source_dir_path)) {
        stop("source_dir_path directory does not exist: ", source_dir_path)
    }

    # Get git information
    git_info <- tryCatch(
        {
            if (requireNamespace("gh", quietly = TRUE)) {
                # Make the API call (works for public repos without token)
                branch_info <- gh::gh("/repos/APAF-BIOINFORMATICS/MultiScholaR/branches/main")
                list(
                    commit_sha = branch_info$commit$sha,
                    branch = "main",
                    repo = "MultiScholaR",
                    timestamp = branch_info$commit$commit$author$date
                )
            } else {
                list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
            }
        },
        error = function(e) {
            message("Error fetching GitHub info: ", e$message)
            list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
        }
    )

    # Get organism information from session if not provided
    if (is.null(organism_name) && exists("organism_name", envir = .GlobalEnv)) {
        organism_name <- get("organism_name", envir = .GlobalEnv)
    }

    if (is.null(taxon_id) && exists("taxon_id", envir = .GlobalEnv)) {
        taxon_id <- get("taxon_id", envir = .GlobalEnv)
    }

    # Get config list
    config_list <- if (exists(config_list_name, envir = env)) {
        get(config_list_name, envir = env)
    } else {
        list()
    }

    # Get contrasts_tbl if not provided
    if (is.null(contrasts_tbl)) {
        if (exists("contrasts_tbl", envir = parent.frame())) {
            contrasts_tbl <- get("contrasts_tbl", envir = parent.frame())
        } else if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
        }
    }

    # Build the output content
    output_lines <- c(
        "Study Parameters",
        "================",
        "",
        "Basic Information:",
        "-----------------",
        paste("Workflow Name:", workflow_name),
        paste("Description:", if (nzchar(description)) description else "N/A"),
        paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        ""
    )

    # Add git information
    if (!is.null(git_info) && is.list(git_info)) {
        git_lines <- c(
            "Git Information:",
            "---------------",
            paste("Repository:", ifelse(!is.null(git_info$repo) && !is.na(git_info$repo), git_info$repo, "N/A")),
            paste("Branch:", ifelse(!is.null(git_info$branch) && !is.na(git_info$branch), git_info$branch, "N/A")),
            paste("Commit:", ifelse(!is.null(git_info$commit_sha) && !is.na(git_info$commit_sha), substr(git_info$commit_sha, 1, 7), "N/A")),
            paste("Git Timestamp:", ifelse(!is.null(git_info$timestamp) && !is.na(git_info$timestamp), git_info$timestamp, "N/A")),
            ""
        )
        output_lines <- c(output_lines, git_lines)
    }

    # Add organism information
    if (!is.null(organism_name) && !is.na(organism_name) && nzchar(organism_name)) {
        organism_lines <- c(
            "Organism Information:",
            "---------------------",
            paste("Organism Name:", organism_name)
        )
        if (!is.null(taxon_id) && !is.na(taxon_id) && nzchar(taxon_id)) {
            organism_lines <- c(organism_lines, paste("Taxon ID:", taxon_id))
        }
        organism_lines <- c(organism_lines, "")
        output_lines <- c(output_lines, organism_lines)
    }

    # Add configuration parameters
    config_lines <- c(
        "Configuration Parameters:",
        "-------------------------"
    )

    # Clean the config list (remove problematic objects)
    clean_config <- config_list
    # Remove the internal source dir if it exists
    if (!is.null(clean_config$internal_workflow_source_dir)) {
        clean_config$internal_workflow_source_dir <- NULL
    }

    # Format the config list
    if (length(clean_config) > 0) {
        config_params <- formatConfigList(clean_config)
        config_lines <- c(config_lines, config_params)
    } else {
        config_lines <- c(config_lines, "No configuration parameters available")
    }

    output_lines <- c(output_lines, config_lines)

    # Add contrasts information
    if (!is.null(contrasts_tbl) && (is.data.frame(contrasts_tbl) || tibble::is_tibble(contrasts_tbl)) && nrow(contrasts_tbl) > 0) {
        contrasts_lines <- c(
            "",
            "Contrasts:",
            "----------"
        )

        if ("contrasts" %in% colnames(contrasts_tbl)) {
            contrasts_info <- tryCatch(
                {
                    contrasts_col <- contrasts_tbl[["contrasts"]]
                    paste("  ", as.character(contrasts_col))
                },
                error = function(e) {
                    paste("  [Error extracting contrasts:", e$message, "]")
                }
            )
            contrasts_lines <- c(contrasts_lines, contrasts_info)
        } else {
            contrasts_lines <- c(contrasts_lines, "  [Column 'contrasts' not found in contrasts_tbl]")
        }

        output_lines <- c(output_lines, contrasts_lines)
    }

    # Write to file
    output_file <- file.path(source_dir_path, "study_parameters.txt")

    tryCatch(
        {
            writeLines(output_lines, output_file)
            message("Study parameters saved to: ", output_file)
            return(output_file)
        },
        error = function(e) {
            stop("Failed to write study parameters file: ", e$message)
        }
    )
}

# Updated function to extract parameters from S4 @args slot
#' @title Create workflow arguments from config
#' @param workflow_name Character string, name of the workflow
#' @param description Character string, description of the analysis
#' @param organism_name Character string, organism name (optional)
#' @param taxon_id Character string or numeric, taxon ID (optional)
#' @param source_dir_path Character string, path to save the study_parameters.txt file
#' @param final_s4_object S4 object containing workflow parameters in @args slot (optional)
#' @param contrasts_tbl Data frame, contrasts table (optional)
#' @param workflow_data List containing workflow state and optimization results (optional)
#' @export
createWorkflowArgsFromConfig <- function(workflow_name, description = "",
                                         organism_name = NULL, taxon_id = NULL,
                                         source_dir_path = NULL,
                                         final_s4_object = NULL,
                                         contrasts_tbl = NULL,
                                         workflow_data = NULL) {
    # Validate required inputs
    if (missing(workflow_name) || !is.character(workflow_name) || length(workflow_name) != 1) {
        stop("workflow_name must be a single character string")
    }

    if (is.null(source_dir_path) || !is.character(source_dir_path) || length(source_dir_path) != 1) {
        stop("source_dir_path must be a single character string")
    }

    if (!dir.exists(source_dir_path)) {
        stop("source_dir_path directory does not exist: ", source_dir_path)
    }

    cat("WORKFLOW ARGS: Starting parameter extraction\n")

    # Extract parameters from S4 object @args if provided
    s4_params <- list()
    # Check if final_s4_object has @args slot
    s4_has_args <- tryCatch(
        {
            !is.null(final_s4_object) && isS4(final_s4_object) && !is.null(final_s4_object@args)
        },
        error = function(e) {
            FALSE
        }
    )

    if (s4_has_args) {
        cat("WORKFLOW ARGS: Extracting parameters from S4 @args slot\n")

        tryCatch(
            {
                s4_params <- final_s4_object@args
                cat(sprintf("WORKFLOW ARGS: Found %d function groups in S4 @args\n", length(s4_params)))

                # Log the function groups for debugging
                for (func_name in names(s4_params)) {
                    param_count <- length(s4_params[[func_name]])
                    cat(sprintf("WORKFLOW ARGS: Function '%s' has %d parameters\n", func_name, param_count))
                }
            },
            error = function(e) {
                cat(sprintf("WORKFLOW ARGS: Error extracting S4 params: %s\n", e$message))
                s4_params <- list()
            }
        )
    } else {
        cat("WORKFLOW ARGS: No S4 object provided or S4 @args is NULL\n")
    }

    # Get fallback config_list from global environment
    config_list <- if (exists("config_list", envir = .GlobalEnv)) {
        get("config_list", envir = .GlobalEnv)
    } else {
        list()
    }

    cat(sprintf("WORKFLOW ARGS: Global config_list has %d sections\n", length(config_list)))

    # Extract RUV optimization results from workflow_data if available
    ruv_optimization_result <- NULL
    if (!is.null(workflow_data)) {
        cat("WORKFLOW ARGS: Checking for RUV optimization results in workflow_data\n")

        # Check multiple possible locations for RUV optimization results
        # First check for plural (metabolomics per-assay format)
        if (!is.null(workflow_data$ruv_optimization_results)) {
            ruv_optimization_result <- workflow_data$ruv_optimization_results
            cat("WORKFLOW ARGS: Found RUV optimization results in workflow_data$ruv_optimization_results (plural)\n")
        } else if (!is.null(workflow_data$ruv_optimization_result)) {
            ruv_optimization_result <- workflow_data$ruv_optimization_result
            cat("WORKFLOW ARGS: Found RUV optimization results in workflow_data$ruv_optimization_result\n")
        } else if (!is.null(workflow_data$state_manager)) {
            # Try to get RUV results from state manager config
            tryCatch(
                {
                    current_state_config <- workflow_data$state_manager$getStateConfig(workflow_data$state_manager$current_state)
                    if (!is.null(current_state_config$ruv_optimization_result)) {
                        ruv_optimization_result <- current_state_config$ruv_optimization_result
                        cat("WORKFLOW ARGS: Found RUV optimization results in state manager config\n")
                    }
                },
                error = function(e) {
                    cat(sprintf("WORKFLOW ARGS: Could not extract RUV results from state manager: %s\n", e$message))
                }
            )
        }

        # Also check if there's a saved RUV file
        if (is.null(ruv_optimization_result) && !is.null(source_dir_path)) {
            ruv_file <- file.path(source_dir_path, "ruv_optimization_results.RDS")
            if (file.exists(ruv_file)) {
                tryCatch(
                    {
                        loaded_result <- readRDS(ruv_file)
                        # Check if loaded result indicates RUV was skipped
                        if (isTRUE(loaded_result$ruv_skipped)) {
                            cat(sprintf("WORKFLOW ARGS: Loaded RUV skip result from file: %s\n", ruv_file))
                            ruv_optimization_result <- loaded_result
                        } else {
                            # Only use file if it's NOT a skip result (backwards compatibility)
                            ruv_optimization_result <- loaded_result
                            cat(sprintf("WORKFLOW ARGS: Loaded RUV optimization results from file: %s\n", ruv_file))
                        }
                    },
                    error = function(e) {
                        cat(sprintf("WORKFLOW ARGS: Could not load RUV file: %s\n", e$message))
                    }
                )
            }
        }
    }
    # Merge S4 parameters with config_list (S4 takes precedence)
    merged_config <- config_list

    if (length(s4_params) > 0) {
        cat("WORKFLOW ARGS: Merging S4 parameters with config_list\n")

        # S4 parameters override config_list parameters
        for (func_name in names(s4_params)) {
            if (is.list(s4_params[[func_name]]) && length(s4_params[[func_name]]) > 0) {
                merged_config[[func_name]] <- s4_params[[func_name]]
                cat(sprintf(
                    "WORKFLOW ARGS: Updated '%s' section with %d S4 parameters\n",
                    func_name, length(s4_params[[func_name]])
                ))
            }
        }
    }

    # Get git information - prefer local repo info over remote API
    git_info <- tryCatch(
        {
            # Try to find the package source directory (works with devtools::load_all())
            pkg_path <- tryCatch(
                {
                    # First try: find via package namespace path
                    pkg_ns_path <- system.file(package = "MultiScholaR")
                    if (nzchar(pkg_ns_path) && dir.exists(file.path(dirname(pkg_ns_path), ".git"))) {
                        dirname(pkg_ns_path)
                    } else {
                        # Second try: check common development locations
                        dev_paths <- c(
                            getwd(),
                            Sys.getenv("MULTISCHOLAR_PATH", unset = NA),
                            file.path(Sys.getenv("HOME"), "Documents", "MultiScholaR")
                        )
                        found_path <- NA_character_
                        for (p in dev_paths) {
                            if (!is.na(p) && nzchar(p) && dir.exists(file.path(p, ".git"))) {
                                found_path <- p
                                break
                            }
                        }
                        found_path
                    }
                },
                error = function(e) NA_character_
            )

            # If we found a git repo, use local git commands
            if (!is.na(pkg_path) && dir.exists(file.path(pkg_path, ".git"))) {
                cat(sprintf("WORKFLOW ARGS: Found local git repo at %s\n", pkg_path))

                # Get current branch
                branch <- tryCatch(
                    {
                        result <- system2("git",
                            args = c("-C", pkg_path, "branch", "--show-current"),
                            stdout = TRUE, stderr = FALSE
                        )
                        if (length(result) > 0 && nzchar(result[1])) result[1] else NA_character_
                    },
                    error = function(e) NA_character_
                )

                # Get current commit SHA
                commit_sha <- tryCatch(
                    {
                        result <- system2("git",
                            args = c("-C", pkg_path, "rev-parse", "HEAD"),
                            stdout = TRUE, stderr = FALSE
                        )
                        if (length(result) > 0 && nzchar(result[1])) result[1] else NA_character_
                    },
                    error = function(e) NA_character_
                )

                # Get commit timestamp
                timestamp <- tryCatch(
                    {
                        result <- system2("git",
                            args = c("-C", pkg_path, "log", "-1", "--format=%cI"),
                            stdout = TRUE, stderr = FALSE
                        )
                        if (length(result) > 0 && nzchar(result[1])) result[1] else NA_character_
                    },
                    error = function(e) NA_character_
                )

                # Get repo name from remote or directory name
                repo_name <- tryCatch(
                    {
                        remote <- system2("git",
                            args = c("-C", pkg_path, "remote", "get-url", "origin"),
                            stdout = TRUE, stderr = FALSE
                        )
                        if (length(remote) > 0 && nzchar(remote[1])) {
                            # Extract repo name from URL
                            basename(gsub("\\.git$", "", remote[1]))
                        } else {
                            basename(pkg_path)
                        }
                    },
                    error = function(e) basename(pkg_path)
                )

                list(
                    commit_sha = commit_sha,
                    branch = branch,
                    repo = repo_name,
                    timestamp = timestamp,
                    source = "local_git"
                )
            } else {
                # Fallback: No local git, return package version info
                cat("WORKFLOW ARGS: No local git repo found, using package info\n")
                pkg_version <- tryCatch(
                    {
                        as.character(utils::packageVersion("MultiScholaR"))
                    },
                    error = function(e) "unknown"
                )

                list(
                    commit_sha = NA_character_,
                    branch = NA_character_,
                    repo = "MultiScholaR",
                    timestamp = NA_character_,
                    package_version = pkg_version,
                    source = "package_info"
                )
            }
        },
        error = function(e) {
            cat(sprintf("WORKFLOW ARGS: Error getting git info: %s\n", e$message))
            list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
        }
    )

    # Get organism information from session if not provided
    if (is.null(organism_name) && exists("organism_name", envir = .GlobalEnv)) {
        organism_name <- get("organism_name", envir = .GlobalEnv)
    }

    if (is.null(taxon_id) && exists("taxon_id", envir = .GlobalEnv)) {
        taxon_id <- get("taxon_id", envir = .GlobalEnv)
    }

    # Get contrasts_tbl if not provided
    if (is.null(contrasts_tbl)) {
        if (exists("contrasts_tbl", envir = parent.frame())) {
            contrasts_tbl <- get("contrasts_tbl", envir = parent.frame())
        } else if (exists("contrasts_tbl", envir = .GlobalEnv)) {
            contrasts_tbl <- get("contrasts_tbl", envir = .GlobalEnv)
        }
    }

    # Build the output content
    output_lines <- c(
        "Study Parameters",
        "================",
        "",
        "Basic Information:",
        "-----------------",
        paste("Workflow Name:", workflow_name),
        paste("Description:", if (nzchar(description)) description else "N/A"),
        paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        ""
    )

    # Add git/version information
    if (!is.null(git_info) && is.list(git_info)) {
        git_lines <- c(
            "Version Information:",
            "--------------------",
            paste("Repository:", ifelse(!is.null(git_info$repo) && !is.na(git_info$repo), git_info$repo, "N/A")),
            paste("Branch:", ifelse(!is.null(git_info$branch) && !is.na(git_info$branch), git_info$branch, "N/A")),
            paste("Commit:", ifelse(!is.null(git_info$commit_sha) && !is.na(git_info$commit_sha), substr(git_info$commit_sha, 1, 7), "N/A")),
            paste("Commit Timestamp:", ifelse(!is.null(git_info$timestamp) && !is.na(git_info$timestamp), git_info$timestamp, "N/A"))
        )

        # Add package version if available
        if (!is.null(git_info$package_version)) {
            git_lines <- c(git_lines, paste("Package Version:", git_info$package_version))
        }

        # Add source info (local_git or package_info)
        if (!is.null(git_info$source)) {
            source_text <- if (git_info$source == "local_git") {
                "Local Git Repository"
            } else {
                "Installed Package"
            }
            git_lines <- c(git_lines, paste("Source:", source_text))
        }

        git_lines <- c(git_lines, "")
        output_lines <- c(output_lines, git_lines)
    }

    # Add organism information
    if (!is.null(organism_name) && !is.na(organism_name) && nzchar(organism_name)) {
        organism_lines <- c(
            "Organism Information:",
            "---------------------",
            paste("Organism Name:", organism_name)
        )
        if (!is.null(taxon_id) && !is.na(taxon_id) && nzchar(taxon_id)) {
            organism_lines <- c(organism_lines, paste("Taxon ID:", taxon_id))
        }
        organism_lines <- c(organism_lines, "")
        output_lines <- c(output_lines, organism_lines)
    }

    # -- Metabolomics-Specific Sections --
    metab_params <- NULL
    if (!is.null(final_s4_object) && inherits(final_s4_object, "MetaboliteAssayData")) {
        cat("WORKFLOW ARGS: Detected MetaboliteAssayData - extracting metabolomics parameters\n")
        metab_params <- extractMetabS4Params(final_s4_object)

        # Add Assay Information section
        if (!is.null(metab_params$assay_info)) {
            assay_info <- metab_params$assay_info
            assay_lines <- c(
                "Assay Information:",
                "-------------------",
                paste("* Number of Assays:", assay_info$num_assays %||% "N/A"),
                paste("* Assay Names:", assay_info$assay_names %||% "N/A"),
                paste("* Total Metabolites:", assay_info$total_metabolites %||% "N/A")
            )

            # Add per-assay metabolite counts
            if (!is.null(assay_info$metabolites_per_assay) && length(assay_info$metabolites_per_assay) > 0) {
                assay_lines <- c(assay_lines, "* Metabolites per Assay:")
                for (assay_name in names(assay_info$metabolites_per_assay)) {
                    count <- assay_info$metabolites_per_assay[[assay_name]]
                    assay_lines <- c(assay_lines, sprintf("    - %s: %s", assay_name, count))
                }
            }

            assay_lines <- c(assay_lines, "")
            output_lines <- c(output_lines, assay_lines)
            cat("WORKFLOW ARGS: Added Assay Information section\n")
        }

        # Add ITSD Normalization section
        if (!is.null(metab_params$itsd_normalization)) {
            itsd <- metab_params$itsd_normalization
            itsd_lines <- c(
                "ITSD Normalization:",
                "-------------------",
                paste("* Applied:", ifelse(isTRUE(itsd$applied), "Yes", "No"))
            )

            if (isTRUE(itsd$applied)) {
                itsd_lines <- c(
                    itsd_lines,
                    paste("* Method Type:", itsd$method_type %||% "N/A"),
                    paste("* Aggregation:", itsd$aggregation %||% "N/A"),
                    paste("* Pattern Columns:", itsd$pattern_columns %||% "N/A"),
                    paste("* ITSD Removed After Normalization:", ifelse(isTRUE(itsd$removed_after_norm), "Yes", "No"))
                )
                if (!is.null(itsd$timestamp) && !is.na(itsd$timestamp)) {
                    itsd_lines <- c(itsd_lines, paste("* Timestamp:", itsd$timestamp))
                }

                # Add per-assay ITSD feature names
                if (!is.null(itsd$features_per_assay) && is.list(itsd$features_per_assay) && length(itsd$features_per_assay) > 0) {
                    itsd_lines <- c(itsd_lines, "* ITSD Features Per Assay:")
                    for (assay_name in names(itsd$features_per_assay)) {
                        features <- itsd$features_per_assay[[assay_name]]
                        count <- length(features)
                        # Show first few feature names, truncate if many
                        preview <- if (count <= 5) {
                            paste(features, collapse = ", ")
                        } else {
                            paste(c(head(features, 5), sprintf("... +%d more", count - 5)), collapse = ", ")
                        }
                        itsd_lines <- c(itsd_lines, sprintf("    - %s (%d features): %s", assay_name, count, preview))
                    }
                    cat(sprintf("WORKFLOW ARGS: Added ITSD features for %d assays\n", length(itsd$features_per_assay)))
                }
            }

            itsd_lines <- c(itsd_lines, "")
            output_lines <- c(output_lines, itsd_lines)
            cat("WORKFLOW ARGS: Added ITSD Normalization section\n")
        }

        # Add Log Transformation section
        if (!is.null(metab_params$log_transformation)) {
            log_trans <- metab_params$log_transformation
            log_lines <- c(
                "Log Transformation:",
                "------------------",
                paste("* Applied:", ifelse(isTRUE(log_trans$applied), "Yes", "No"))
            )

            if (isTRUE(log_trans$applied)) {
                log_lines <- c(
                    log_lines,
                    paste("* Offset:", log_trans$offset %||% "N/A")
                )
            }

            log_lines <- c(log_lines, "")
            output_lines <- c(output_lines, log_lines)
            cat("WORKFLOW ARGS: Added Log Transformation section\n")
        }

        # Add Between-Sample Normalization section
        if (!is.null(metab_params$normalisation_method) && !is.na(metab_params$normalisation_method)) {
            norm_lines <- c(
                "Between-Sample Normalization:",
                "-----------------------------",
                paste("* Method:", metab_params$normalisation_method),
                ""
            )
            output_lines <- c(output_lines, norm_lines)
            cat("WORKFLOW ARGS: Added Normalization Method section\n")
        }

        # Add Metadata section
        if (!is.null(metab_params$metadata) && length(metab_params$metadata) > 0) {
            meta <- metab_params$metadata
            meta_lines <- c(
                "Metabolomics Metadata:",
                "----------------------",
                paste("* Metabolite ID Column:", meta$metabolite_id_column %||% "N/A"),
                paste("* Annotation ID Column:", meta$annotation_id_column %||% "N/A"),
                paste("* Database Identifier Type:", meta$database_identifier_type %||% "N/A"),
                paste("* Sample ID Column:", meta$sample_id %||% "N/A"),
                paste("* Group ID Column:", meta$group_id %||% "N/A")
            )

            if (!is.null(meta$n_samples)) {
                meta_lines <- c(
                    meta_lines,
                    paste("* Number of Samples:", meta$n_samples),
                    paste("* Number of Groups:", meta$n_groups %||% "N/A"),
                    paste("* Group Names:", meta$group_names %||% "N/A")
                )
            }

            if (!is.null(meta$internal_standard_regex) && !is.na(meta$internal_standard_regex)) {
                meta_lines <- c(
                    meta_lines,
                    paste("* Internal Standard Regex:", meta$internal_standard_regex)
                )
            }

            meta_lines <- c(meta_lines, "")
            output_lines <- c(output_lines, meta_lines)
            cat("WORKFLOW ARGS: Added Metabolomics Metadata section\n")
        }

        # Add RUV-III Batch Correction section from S4 params (metabolomics-specific)
        if (!is.null(metab_params$ruv_params) && isTRUE(metab_params$ruv_params$applied)) {
            ruv_section <- c(
                "RUV-III Batch Correction:",
                "-------------------------",
                "* Status: Applied"
            )

            if (isTRUE(metab_params$ruv_params$per_assay)) {
                # Per-assay RUV parameters
                ruv_section <- c(
                    ruv_section,
                    paste("* Grouping Variable:", metab_params$ruv_params$grouping_variable %||% "N/A"),
                    "* Per-Assay Parameters:"
                )
                for (assay_name in names(metab_params$ruv_params$k_per_assay)) {
                    k_val <- metab_params$ruv_params$k_per_assay[[assay_name]]
                    ctrl_count <- if (!is.null(metab_params$ruv_params$ctrl_per_assay[[assay_name]])) {
                        metab_params$ruv_params$ctrl_per_assay[[assay_name]]
                    } else {
                        "N/A"
                    }
                    ruv_section <- c(
                        ruv_section,
                        sprintf("    - %s: k=%d, controls=%s", assay_name, k_val, ctrl_count)
                    )
                }
                cat("WORKFLOW ARGS: Added per-assay RUV-III section\n")
            } else {
                # Single RUV value
                ruv_section <- c(
                    ruv_section,
                    paste("* Grouping Variable:", metab_params$ruv_params$grouping_variable %||% "N/A"),
                    paste("* k value:", metab_params$ruv_params$number_k %||% "N/A"),
                    paste("* Control features:", metab_params$ruv_params$ctrl_count %||% "N/A")
                )
                cat("WORKFLOW ARGS: Added single RUV-III section\n")
            }

            ruv_section <- c(ruv_section, "")
            output_lines <- c(output_lines, ruv_section)
        } else if (!is.null(metab_params) && (is.null(metab_params$ruv_params) || !isTRUE(metab_params$ruv_params$applied))) {
            # Metabolomics without RUV applied
            ruv_section <- c(
                "RUV-III Batch Correction:",
                "-------------------------",
                "* Status: Not Applied",
                ""
            )
            output_lines <- c(output_lines, ruv_section)
            cat("WORKFLOW ARGS: RUV not applied for metabolomics\n")
        }
    }

    # -- Lipidomics-Specific Sections --
    lipid_params <- NULL
    if (!is.null(final_s4_object) && inherits(final_s4_object, "LipidomicsAssayData")) {
        cat("WORKFLOW ARGS: Detected LipidomicsAssayData - extracting lipidomics parameters\n")
        lipid_params <- extractLipidS4Params(final_s4_object)

        # Add Assay Information section
        if (!is.null(lipid_params$assay_info)) {
            assay_info <- lipid_params$assay_info
            assay_lines <- c(
                "Assay Information:",
                "-------------------",
                paste("* Number of Assays:", assay_info$num_assays %||% "N/A"),
                paste("* Assay Names:", assay_info$assay_names %||% "N/A"),
                paste("* Total Lipids:", assay_info$total_lipids %||% "N/A")
            )

            # Add per-assay lipid counts
            if (!is.null(assay_info$lipids_per_assay) && length(assay_info$lipids_per_assay) > 0) {
                assay_lines <- c(assay_lines, "* Lipids per Assay:")
                for (assay_name in names(assay_info$lipids_per_assay)) {
                    count <- assay_info$lipids_per_assay[[assay_name]]
                    assay_lines <- c(assay_lines, sprintf("    - %s: %s", assay_name, count))
                }
            }

            assay_lines <- c(assay_lines, "")
            output_lines <- c(output_lines, assay_lines)
            cat("WORKFLOW ARGS: Added Lipidomics Assay Information section\n")
        }

        # Add ITSD Normalization section
        if (!is.null(lipid_params$itsd_normalization)) {
            itsd <- lipid_params$itsd_normalization
            itsd_lines <- c(
                "ITSD Normalization:",
                "-------------------",
                paste("* Applied:", ifelse(isTRUE(itsd$applied), "Yes", "No"))
            )

            if (isTRUE(itsd$applied)) {
                itsd_lines <- c(
                    itsd_lines,
                    paste("* Method Type:", itsd$method_type %||% "N/A"),
                    paste("* Aggregation:", itsd$aggregation %||% "N/A"),
                    paste("* Pattern Columns:", itsd$pattern_columns %||% "N/A"),
                    paste("* ITSD Removed After Normalization:", ifelse(isTRUE(itsd$removed_after_norm), "Yes", "No"))
                )
                if (!is.null(itsd$timestamp) && !is.na(itsd$timestamp)) {
                    itsd_lines <- c(itsd_lines, paste("* Timestamp:", itsd$timestamp))
                }

                # Add per-assay ITSD feature names
                if (!is.null(itsd$features_per_assay) && is.list(itsd$features_per_assay) && length(itsd$features_per_assay) > 0) {
                    itsd_lines <- c(itsd_lines, "* ITSD Features Per Assay:")
                    for (assay_name in names(itsd$features_per_assay)) {
                        features <- itsd$features_per_assay[[assay_name]]
                        count <- length(features)
                        # Show first few feature names, truncate if many
                        preview <- if (count <= 5) {
                            paste(features, collapse = ", ")
                        } else {
                            paste(c(head(features, 5), sprintf("... +%d more", count - 5)), collapse = ", ")
                        }
                        itsd_lines <- c(itsd_lines, sprintf("    - %s (%d features): %s", assay_name, count, preview))
                    }
                    cat(sprintf("WORKFLOW ARGS: Added ITSD features for %d assays\n", length(itsd$features_per_assay)))
                }
            }

            itsd_lines <- c(itsd_lines, "")
            output_lines <- c(output_lines, itsd_lines)
            cat("WORKFLOW ARGS: Added ITSD Normalization section\n")
        }

        # Add Log Transformation section
        if (!is.null(lipid_params$log_transformation)) {
            log_trans <- lipid_params$log_transformation
            log_lines <- c(
                "Log Transformation:",
                "------------------",
                paste("* Applied:", ifelse(isTRUE(log_trans$applied), "Yes", "No"))
            )

            if (isTRUE(log_trans$applied)) {
                log_lines <- c(
                    log_lines,
                    paste("* Offset:", log_trans$offset %||% "N/A")
                )
            }

            log_lines <- c(log_lines, "")
            output_lines <- c(output_lines, log_lines)
            cat("WORKFLOW ARGS: Added Log Transformation section\n")
        }

        # Add Between-Sample Normalization section
        if (!is.null(lipid_params$normalisation_method) && !is.na(lipid_params$normalisation_method)) {
            norm_lines <- c(
                "Between-Sample Normalization:",
                "-----------------------------",
                paste("* Method:", lipid_params$normalisation_method),
                ""
            )
            output_lines <- c(output_lines, norm_lines)
            cat("WORKFLOW ARGS: Added Normalization Method section\n")
        }

        # Add Metadata section
        if (!is.null(lipid_params$metadata) && length(lipid_params$metadata) > 0) {
            meta <- lipid_params$metadata
            meta_lines <- c(
                "Lipidomics Metadata:",
                "----------------------",
                paste("* Lipid ID Column:", meta$lipid_id_column %||% "N/A"),
                paste("* Annotation ID Column:", meta$annotation_id_column %||% "N/A"),
                paste("* Database Identifier Type:", meta$database_identifier_type %||% "N/A"),
                paste("* Sample ID Column:", meta$sample_id %||% "N/A"),
                paste("* Group ID Column:", meta$group_id %||% "N/A")
            )

            if (!is.null(meta$n_samples)) {
                meta_lines <- c(
                    meta_lines,
                    paste("* Number of Samples:", meta$n_samples),
                    paste("* Number of Groups:", meta$n_groups %||% "N/A"),
                    paste("* Group Names:", meta$group_names %||% "N/A")
                )
            }

            if (!is.null(meta$internal_standard_regex) && !is.na(meta$internal_standard_regex)) {
                meta_lines <- c(
                    meta_lines,
                    paste("* Internal Standard Regex:", meta$internal_standard_regex)
                )
            }

            meta_lines <- c(meta_lines, "")
            output_lines <- c(output_lines, meta_lines)
            cat("WORKFLOW ARGS: Added Lipidomics Metadata section\n")
        }

        # Add RUV-III Batch Correction section from S4 params (lipidomics-specific)
        if (!is.null(lipid_params$ruv_params) && isTRUE(lipid_params$ruv_params$applied)) {
            ruv_section <- c(
                "RUV-III Batch Correction:",
                "-------------------------",
                "* Status: Applied"
            )

            if (isTRUE(lipid_params$ruv_params$per_assay)) {
                # Per-assay RUV parameters
                ruv_section <- c(
                    ruv_section,
                    paste("* Grouping Variable:", lipid_params$ruv_params$grouping_variable %||% "N/A"),
                    "* Per-Assay Parameters:"
                )
                for (assay_name in names(lipid_params$ruv_params$k_per_assay)) {
                    k_val <- lipid_params$ruv_params$k_per_assay[[assay_name]]
                    ctrl_count <- if (!is.null(lipid_params$ruv_params$ctrl_per_assay[[assay_name]])) {
                        lipid_params$ruv_params$ctrl_per_assay[[assay_name]]
                    } else {
                        "N/A"
                    }
                    ruv_section <- c(
                        ruv_section,
                        sprintf("    - %s: k=%d, controls=%s", assay_name, k_val, ctrl_count)
                    )
                }
                cat("WORKFLOW ARGS: Added per-assay RUV-III section\n")
            } else {
                # Single RUV value
                ruv_section <- c(
                    ruv_section,
                    paste("* Grouping Variable:", lipid_params$ruv_params$grouping_variable %||% "N/A"),
                    paste("* k value:", lipid_params$ruv_params$number_k %||% "N/A"),
                    paste("* Control features:", lipid_params$ruv_params$ctrl_count %||% "N/A")
                )
                cat("WORKFLOW ARGS: Added single RUV-III section\n")
            }

            ruv_section <- c(ruv_section, "")
            output_lines <- c(output_lines, ruv_section)
        } else if (!is.null(lipid_params) && (is.null(lipid_params$ruv_params) || !isTRUE(lipid_params$ruv_params$applied))) {
            # Lipidomics without RUV applied
            ruv_section <- c(
                "RUV-III Batch Correction:",
                "-------------------------",
                "* Status: Not Applied",
                ""
            )
            output_lines <- c(output_lines, ruv_section)
            cat("WORKFLOW ARGS: RUV not applied for lipidomics\n")
        }
    }

    # Add RUV optimization results if available
    if (!is.null(ruv_optimization_result) && is.list(ruv_optimization_result)) {
        # Check if RUV was skipped
        ruv_was_skipped <- isTRUE(ruv_optimization_result$ruv_skipped)

        if (ruv_was_skipped) {
            cat("WORKFLOW ARGS: RUV was skipped - writing skip section\n")
            ruv_lines <- c(
                "RUV-III Batch Correction:",
                "-------------------------",
                "* Status: Not Applied",
                "* Reason: User determined RUV was not appropriate due to dataset constraints",
                ""
            )
            output_lines <- c(output_lines, ruv_lines)
        } else {
            cat("WORKFLOW ARGS: Formatting RUV optimization results\n")

            ruv_lines <- c(
                "Automatic RUV Optimization Results:",
                "-----------------------------------"
            )

            # Extract and format RUV optimization values
            tryCatch(
                {
                    best_percentage <- if (!is.null(ruv_optimization_result$best_percentage)) {
                        sprintf("%.1f%%", ruv_optimization_result$best_percentage)
                    } else {
                        "N/A"
                    }

                    best_k <- if (!is.null(ruv_optimization_result$best_k)) {
                        as.character(ruv_optimization_result$best_k)
                    } else {
                        "N/A"
                    }

                    separation_score <- if (!is.null(ruv_optimization_result$best_separation_score)) {
                        sprintf("%.4f", ruv_optimization_result$best_separation_score)
                    } else {
                        "N/A"
                    }

                    composite_score <- if (!is.null(ruv_optimization_result$best_composite_score)) {
                        sprintf("%.4f", ruv_optimization_result$best_composite_score)
                    } else {
                        "N/A"
                    }

                    control_genes_count <- if (!is.null(ruv_optimization_result$best_control_genes_index)) {
                        as.character(sum(ruv_optimization_result$best_control_genes_index, na.rm = TRUE))
                    } else {
                        "N/A"
                    }

                    separation_metric <- if (!is.null(ruv_optimization_result$separation_metric_used)) {
                        as.character(ruv_optimization_result$separation_metric_used)
                    } else {
                        "N/A"
                    }

                    k_penalty_weight <- if (!is.null(ruv_optimization_result$k_penalty_weight)) {
                        sprintf("%.1f", ruv_optimization_result$k_penalty_weight)
                    } else {
                        "N/A"
                    }

                    adaptive_penalty <- if (!is.null(ruv_optimization_result$adaptive_k_penalty_used)) {
                        ifelse(ruv_optimization_result$adaptive_k_penalty_used, "TRUE", "FALSE")
                    } else {
                        "N/A"
                    }

                    sample_size <- if (!is.null(ruv_optimization_result$sample_size)) {
                        as.character(ruv_optimization_result$sample_size)
                    } else {
                        "N/A"
                    }

                    # Also try to get RUV grouping variable from S4 parameters
                    ruv_grouping_variable <- "N/A"
                    if (!is.null(s4_params$ruvIII_C_Varying$ruv_grouping_variable)) {
                        ruv_grouping_variable <- s4_params$ruvIII_C_Varying$ruv_grouping_variable
                    } else if (!is.null(s4_params$getNegCtrlProtAnova$ruv_grouping_variable)) {
                        ruv_grouping_variable <- s4_params$getNegCtrlProtAnova$ruv_grouping_variable
                    }

                    ruv_lines <- c(
                        ruv_lines,
                        paste("* Best percentage:", best_percentage),
                        paste("* Best k value:", best_k),
                        paste("* Separation score:", separation_score),
                        paste("* Composite score:", composite_score),
                        paste("* Control genes:", control_genes_count),
                        paste("* RUV grouping variable:", ruv_grouping_variable),
                        paste("* Separation metric:", separation_metric),
                        paste("* K penalty weight:", k_penalty_weight),
                        paste("* Adaptive penalty:", adaptive_penalty),
                        paste("* Sample size:", sample_size),
                        ""
                    )

                    cat("WORKFLOW ARGS: Successfully formatted RUV optimization results\n")
                },
                error = function(e) {
                    cat(sprintf("WORKFLOW ARGS: Error formatting RUV results: %s\n", e$message))
                    ruv_lines <- c(
                        ruv_lines,
                        paste("* [Error formatting RUV optimization results:", e$message, "]"),
                        ""
                    )
                }
            )

            output_lines <- c(output_lines, ruv_lines)
        } # End of else block for RUV applied
    } else {
        cat("WORKFLOW ARGS: No RUV optimization results available\n")
    }

    # Add FASTA Processing Information
    if (!is.null(workflow_data) && !is.null(workflow_data$fasta_metadata)) {
        cat("WORKFLOW ARGS: Adding FASTA processing information\n")
        fasta_lines <- c(
            "FASTA File Processing:",
            "---------------------",
            paste("* Format:", workflow_data$fasta_metadata$fasta_format),
            paste("* Sequences:", workflow_data$fasta_metadata$num_sequences),
            paste("* Protein Evidence Available:", workflow_data$fasta_metadata$has_protein_evidence),
            paste("* Gene Names Available:", workflow_data$fasta_metadata$has_gene_names),
            paste("* Isoform Information Available:", workflow_data$fasta_metadata$has_isoform_info),
            paste("* Status Information Available:", workflow_data$fasta_metadata$has_status_info),
            ""
        )
        output_lines <- c(output_lines, fasta_lines)
    } else {
        cat("WORKFLOW ARGS: No FASTA metadata available\n")
    }

    # [OK] NEW: Add Mixed Species FASTA Analysis section
    if (!is.null(workflow_data) && !is.null(workflow_data$mixed_species_analysis)) {
        cat("WORKFLOW ARGS: Adding mixed species analysis information\n")
        mixed_species_info <- workflow_data$mixed_species_analysis

        mixed_species_lines <- c(
            "Mixed Species FASTA Analysis:",
            "-----------------------------",
            paste("* Multi-Species FASTA Used:", ifelse(isTRUE(mixed_species_info$enabled), "Yes", "No"))
        )

        if (isTRUE(mixed_species_info$enabled)) {
            mixed_species_lines <- c(
                mixed_species_lines,
                paste("* Selected Primary Organism:", mixed_species_info$selected_organism %||% "N/A"),
                paste("* Selected Taxon ID:", mixed_species_info$selected_taxon_id %||% "N/A"),
                paste("* Filtered at Import:", ifelse(isTRUE(mixed_species_info$filter_applied_at_import), "Yes", "No"))
            )

            # Add organism distribution summary if available
            if (!is.null(mixed_species_info$organism_distribution) &&
                is.data.frame(mixed_species_info$organism_distribution) &&
                nrow(mixed_species_info$organism_distribution) > 0) {
                mixed_species_lines <- c(
                    mixed_species_lines,
                    "",
                    "  Organism Distribution in FASTA:"
                )

                # Add top organisms (limit to top 5)
                top_orgs <- utils::head(mixed_species_info$organism_distribution, 5)
                for (i in seq_len(nrow(top_orgs))) {
                    org_line <- sprintf(
                        "    - %s (Taxon %s): %d proteins (%.1f%%)",
                        top_orgs$organism_name[i] %||% "Unknown",
                        top_orgs$taxon_id[i] %||% "N/A",
                        top_orgs$protein_count[i] %||% 0,
                        top_orgs$percentage[i] %||% 0
                    )
                    mixed_species_lines <- c(mixed_species_lines, org_line)
                }
            }
        }

        mixed_species_lines <- c(mixed_species_lines, "")
        output_lines <- c(output_lines, mixed_species_lines)
    }

    # [OK] NEW: Add Enrichment Organism Filtering section
    if (!is.null(workflow_data) && !is.null(workflow_data$enrichment_organism_filter)) {
        cat("WORKFLOW ARGS: Adding enrichment organism filtering information\n")
        filter_info <- workflow_data$enrichment_organism_filter

        filter_lines <- c(
            "Enrichment Analysis - Organism Filtering:",
            "-----------------------------------------",
            paste("* Organism Filter Enabled:", ifelse(isTRUE(filter_info$enabled), "Yes", "No"))
        )

        if (isTRUE(filter_info$filter_applied)) {
            filter_lines <- c(
                filter_lines,
                paste("* Filter Applied:", "Yes"),
                paste("* Target Taxon ID:", filter_info$target_taxon_id %||% "N/A"),
                paste("* Proteins Before Filtering:", filter_info$proteins_before %||% "N/A"),
                paste("* Proteins After Filtering:", filter_info$proteins_after %||% "N/A"),
                paste("* Proteins Removed:", filter_info$proteins_removed %||% "N/A")
            )

            if (!is.null(filter_info$proteins_before) && filter_info$proteins_before > 0) {
                retention_pct <- round((filter_info$proteins_after / filter_info$proteins_before) * 100, 1)
                filter_lines <- c(
                    filter_lines,
                    paste("* Retention Rate:", paste0(retention_pct, "%"))
                )
            }
        }

        filter_lines <- c(filter_lines, "")
        output_lines <- c(output_lines, filter_lines)
    }

    # Add Accession Cleanup Results
    if (!is.null(workflow_data) && !is.null(workflow_data$accession_cleanup_results)) {
        cat("WORKFLOW ARGS: Adding accession cleanup results\n")
        cleanup_lines <- c(
            "Protein Accession Cleanup:",
            "-------------------------",
            paste("* Cleanup Applied:", workflow_data$accession_cleanup_results$cleanup_applied),
            paste("* Aggregation Method:", workflow_data$accession_cleanup_results$aggregation_method),
            paste("* Delimiter Used:", workflow_data$accession_cleanup_results$delimiter_used),
            paste("* Proteins Before Cleanup:", workflow_data$accession_cleanup_results$proteins_before),
            paste("* Proteins After Cleanup:", workflow_data$accession_cleanup_results$proteins_after),
            paste("* Full UniProt Metadata:", workflow_data$accession_cleanup_results$had_full_metadata),
            ""
        )
        output_lines <- c(output_lines, cleanup_lines)
    } else {
        cat("WORKFLOW ARGS: No accession cleanup results available\n")
    }

    # Add Protein Filtering Summary
    if (!is.null(workflow_data) && !is.null(workflow_data$protein_counts)) {
        cat("WORKFLOW ARGS: Adding protein filtering summary\n")
        counts_lines <- c(
            "Protein Filtering Summary:",
            "-------------------------",
            paste("* Proteins after QC filtering:", workflow_data$protein_counts$after_qc_filtering %||% "N/A"),
            paste("* Proteins after RUV filtering:", workflow_data$protein_counts$after_ruv_filtering %||% "N/A"),
            paste("* Final proteins for DE analysis:", workflow_data$protein_counts$final_for_de %||% "N/A"),
            ""
        )
        output_lines <- c(output_lines, counts_lines)
    } else {
        cat("WORKFLOW ARGS: No protein counts available\n")
    }

    # Add configuration parameters (now includes S4 parameters)
    config_lines <- c(
        "Workflow Parameters:",
        "-------------------",
        ""
    )

    # Clean the merged config (remove problematic objects)
    clean_config <- merged_config
    # Remove the internal source dir if it exists
    if (!is.null(clean_config$internal_workflow_source_dir)) {
        clean_config$internal_workflow_source_dir <- NULL
    }

    # Add special section for S4-derived parameters
    if (length(s4_params) > 0) {
        cat("WORKFLOW ARGS: About to format S4 parameters using functional approach\n")

        config_lines <- c(
            config_lines,
            "  Parameters from Final S4 Object:",
            "  --------------------------------"
        )

        # SAFE parameter formatting function
        formatParameterValue <- function(param_value) {
            tryCatch(
                {
                    if (is.null(param_value)) {
                        "NULL"
                    } else if (is.data.frame(param_value)) {
                        # Skip data frames (like seqinr_obj) - too large to serialize
                        sprintf(
                            "[Data frame: %d rows x %d cols - omitted for brevity]",
                            nrow(param_value), ncol(param_value)
                        )
                    } else if (is.logical(param_value)) {
                        if (length(param_value) == 1) {
                            ifelse(param_value, "TRUE", "FALSE")
                        } else if (length(param_value) > 50) {
                            # Handle large logical vectors (like control genes index)
                            true_count <- sum(param_value, na.rm = TRUE)
                            total_count <- length(param_value)
                            sprintf(
                                "logical vector [%d TRUE, %d FALSE out of %d total]",
                                true_count, total_count - true_count, total_count
                            )
                        } else {
                            # Show first few values for smaller vectors
                            preview <- ifelse(utils::head(param_value, 5), "TRUE", "FALSE")
                            if (length(param_value) > 5) {
                                paste0("c(", paste(preview, collapse = ", "), ", ...)")
                            } else {
                                paste0("c(", paste(preview, collapse = ", "), ")")
                            }
                        }
                    } else if (is.numeric(param_value)) {
                        if (length(param_value) == 1) {
                            as.character(param_value)
                        } else if (length(param_value) > 5) {
                            sprintf(
                                "numeric vector [%d values: %s, ...]",
                                length(param_value),
                                paste(as.character(utils::head(param_value, 3)), collapse = ", ")
                            )
                        } else {
                            paste0("c(", paste(as.character(param_value), collapse = ", "), ")")
                        }
                    } else if (is.character(param_value)) {
                        if (length(param_value) == 1) {
                            param_value
                        } else if (length(param_value) > 5) {
                            sprintf(
                                "character vector [%d values: %s, ...]",
                                length(param_value),
                                paste(shQuote(utils::head(param_value, 3)), collapse = ", ")
                            )
                        } else {
                            paste0("c(", paste(shQuote(param_value), collapse = ", "), ")")
                        }
                    } else {
                        # SAFE fallback - no dput() which was causing the hang
                        paste0("[", class(param_value)[1], " object]")
                    }
                },
                error = function(e) {
                    "[SERIALIZATION ERROR]"
                }
            )
        }

        s4_sections <- if (requireNamespace("purrr", quietly = TRUE)) {
            purrr::imap(s4_params, function(func_params, func_name) {
                tryCatch(
                    {
                        if (!is.list(func_params) || length(func_params) == 0) {
                            return("")
                        }

                        cat(sprintf("WORKFLOW ARGS: Processing S4 function group '%s'\n", func_name))
                        header <- sprintf("[%s]", func_name)

                        # Use map instead of imap_chr for safer handling
                        param_lines <- purrr::imap(func_params, function(param_value, param_name) {
                            tryCatch(
                                {
                                    param_str <- formatParameterValue(param_value)
                                    sprintf("  %s = %s", param_name, param_str)
                                },
                                error = function(e) {
                                    cat(sprintf("WORKFLOW ARGS: Error formatting parameter '%s' in '%s': %s\n", param_name, func_name, e$message))
                                    sprintf("  %s = [ERROR: %s]", param_name, e$message)
                                }
                            )
                        })

                        # Convert list to character vector safely
                        param_lines_char <- unlist(param_lines)
                        paste(c(header, param_lines_char, ""), collapse = "\n")
                    },
                    error = function(e) {
                        cat(sprintf("WORKFLOW ARGS: Error processing S4 function group '%s': %s\n", func_name, e$message))
                        sprintf("[%s]\n  [ERROR: %s]\n", func_name, e$message)
                    }
                )
            })
        } else {
            # Fallback using base R if purrr not available
            lapply(names(s4_params), function(func_name) {
                func_params <- s4_params[[func_name]]
                if (!is.list(func_params) || length(func_params) == 0) {
                    return("")
                }

                header <- sprintf("[%s]", func_name)

                param_lines <- lapply(names(func_params), function(param_name) {
                    param_value <- func_params[[param_name]]
                    param_str <- formatParameterValue(param_value)
                    sprintf("  %s = %s", param_name, param_str)
                })

                paste(c(header, unlist(param_lines), ""), collapse = "\n")
            })
        }

        # Add all S4 sections at once - safely handle list from purrr::imap()
        s4_sections_char <- tryCatch(
            {
                if (is.list(s4_sections)) {
                    # Convert list to character vector
                    unlist(s4_sections)
                } else {
                    s4_sections
                }
            },
            error = function(e) {
                cat(sprintf("WORKFLOW ARGS: Error converting S4 sections: %s\n", e$message))
                "[ERROR: Could not process S4 sections]"
            }
        )

        config_lines <- c(config_lines, unlist(strsplit(paste(s4_sections_char, collapse = ""), "\n")))

        cat("WORKFLOW ARGS: S4 parameters formatted successfully\n")
    }

    # Add UI parameters from workflow_data if available (DE and Enrichment UI inputs)
    if (!is.null(workflow_data)) {
        cat("WORKFLOW ARGS: Checking for UI parameters in workflow_data\n")
        ui_sections <- c()

        # Check for DE UI parameters
        if (!is.null(workflow_data$da_ui_params)) {
            cat("WORKFLOW ARGS: Found DE UI parameters in workflow_data\n")
            da_ui_lines <- c(
                "[Differential Expression UI Parameters]",
                sprintf("  q_value_threshold = %s", ifelse(!is.null(workflow_data$da_ui_params$q_value_threshold), workflow_data$da_ui_params$q_value_threshold, "N/A")),
                sprintf("  log_fold_change_cutoff = %s", ifelse(!is.null(workflow_data$da_ui_params$log_fold_change_cutoff), workflow_data$da_ui_params$log_fold_change_cutoff, "N/A")),
                sprintf("  treat_enabled = %s", ifelse(!is.null(workflow_data$da_ui_params$treat_enabled), workflow_data$da_ui_params$treat_enabled, "N/A")),
                ""
            )
            ui_sections <- c(ui_sections, da_ui_lines)
        }

        # Check for Enrichment UI parameters
        if (!is.null(workflow_data$enrichment_ui_params)) {
            cat("WORKFLOW ARGS: Found Enrichment UI parameters in workflow_data\n")
            enrichment_ui_lines <- c(
                "[Enrichment Analysis UI Parameters]",
                sprintf("  up_log2fc_cutoff = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$up_log2fc_cutoff), workflow_data$enrichment_ui_params$up_log2fc_cutoff, "N/A")),
                sprintf("  down_log2fc_cutoff = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$down_log2fc_cutoff), workflow_data$enrichment_ui_params$down_log2fc_cutoff, "N/A")),
                sprintf("  q_value_cutoff = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$q_value_cutoff), workflow_data$enrichment_ui_params$q_value_cutoff, "N/A")),
                sprintf("  organism_selected = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$organism_selected), workflow_data$enrichment_ui_params$organism_selected, "N/A")),
                sprintf("  database_source = %s", ifelse(!is.null(workflow_data$enrichment_ui_params$database_source), workflow_data$enrichment_ui_params$database_source, "N/A")),
                ""
            )
            ui_sections <- c(ui_sections, enrichment_ui_lines)
        }

        # Add UI sections to config_lines if any were found
        if (length(ui_sections) > 0) {
            config_lines <- c(
                config_lines,
                "  User Interface Parameters:",
                "  -------------------------",
                ui_sections
            )
            cat("WORKFLOW ARGS: Added UI parameters to output\n")
        } else {
            cat("WORKFLOW ARGS: No UI parameters found in workflow_data\n")
        }
    }

    # Format the remaining config list
    cat("WORKFLOW ARGS: About to format remaining config list\n")
    if (length(clean_config) > 0) {
        config_lines <- c(
            config_lines,
            "Additional Configuration Parameters:",
            "-----------------------------------"
        )

        cat("WORKFLOW ARGS: Calling formatConfigList...\n")
        tryCatch(
            {
                # Check if formatConfigList exists before calling
                if (exists("formatConfigList", mode = "function")) {
                    config_params <- formatConfigList(clean_config)
                    config_lines <- c(config_lines, config_params)
                    cat("WORKFLOW ARGS: formatConfigList completed successfully\n")
                } else {
                    cat("WORKFLOW ARGS: formatConfigList function not found, using basic formatting\n")
                    # Basic fallback formatting without for loops
                    basic_config <- unlist(clean_config, recursive = TRUE)
                    config_lines <- c(config_lines, names(basic_config), " = ", as.character(basic_config))
                }
            },
            error = function(e) {
                cat(sprintf("WORKFLOW ARGS: formatConfigList failed: %s\n", e$message))
                config_lines <- c(config_lines, paste("Error formatting config:", e$message))
            }
        )
    } else {
        config_lines <- c(config_lines, "No additional configuration parameters available")
    }

    cat("WORKFLOW ARGS: Adding config lines to output\n")
    output_lines <- c(output_lines, config_lines)

    # Add contrasts information
    cat("WORKFLOW ARGS: Processing contrasts information\n")
    if (!is.null(contrasts_tbl) && (is.data.frame(contrasts_tbl) || tibble::is_tibble(contrasts_tbl)) && nrow(contrasts_tbl) > 0) {
        cat("WORKFLOW ARGS: Adding contrasts to output\n")
        contrasts_lines <- c(
            "",
            "Contrasts:",
            "----------"
        )

        if ("contrasts" %in% colnames(contrasts_tbl)) {
            contrasts_info <- tryCatch(
                {
                    contrasts_col <- contrasts_tbl[["contrasts"]]
                    paste("  ", as.character(contrasts_col))
                },
                error = function(e) {
                    paste("  [Error extracting contrasts:", e$message, "]")
                }
            )
            contrasts_lines <- c(contrasts_lines, contrasts_info)
        } else {
            contrasts_lines <- c(contrasts_lines, "  [Column 'contrasts' not found in contrasts_tbl]")
        }

        output_lines <- c(output_lines, contrasts_lines)
        cat("WORKFLOW ARGS: Contrasts added successfully\n")
    } else {
        cat("WORKFLOW ARGS: No contrasts to add\n")
    }

    # Write to file
    cat("WORKFLOW ARGS: About to write file\n")
    output_file <- file.path(source_dir_path, "study_parameters.txt")
    cat(sprintf("WORKFLOW ARGS: Target file path: %s\n", output_file))

    tryCatch(
        {
            cat("WORKFLOW ARGS: Calling writeLines...\n")
            writeLines(output_lines, output_file)
            cat(sprintf("WORKFLOW ARGS: Study parameters saved to: %s\n", output_file))
            return(output_file)
        },
        error = function(e) {
            cat(sprintf("WORKFLOW ARGS: Error writing file: %s\n", e$message))
            stop("Failed to write study parameters file: ", e$message)
        }
    )
}

