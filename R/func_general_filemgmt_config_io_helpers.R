# ----------------------------------------------------------------------------
# readConfigFile
# ----------------------------------------------------------------------------
#' @title Read the config file and return the list of parameters
#' @description Read the config file and return the list of parameters
#' @param file The file path to the config file
#' @param file_type The type of the file (default: "ini")
#' @export
readConfigFile <- function(file = file.path(source_dir, "config.ini"), file_type = "ini") {
    message(sprintf("DEBUG66: readConfigFile() called with file = %s", file))
    message(sprintf("DEBUG66: file exists = %s", file.exists(file)))

    config_list <- read.config(file = file, file.type = file_type)

    message(sprintf(
        "DEBUG66: read.config() completed, config_list names = %s",
        paste(names(config_list), collapse = ", ")
    ))
    message(sprintf(
        "DEBUG66: config_list class = %s, length = %d",
        class(config_list)[1], length(config_list)
    ))

    # to set the number of cores to be used in the parallel processing
    if ("globalParameters" %in% names(config_list)) {
        if ("number_of_cpus" %in% names(config_list[["globalParameters"]])) {
            print(paste0(
                "Read globalParameters: number_of_cpus = ",
                config_list$globalParameters$number_of_cpus
            ))
            core_utilisation <- new_cluster(config_list$globalParameters$number_of_cpus)
            cluster_library(core_utilisation, c("tidyverse", "glue", "rlang", "lazyeval"))

            list_of_multithreaded_functions <- c(
                "rollUpPrecursorToPeptide",
                "peptideIntensityFiltering",
                "filterMinNumPeptidesPerProtein",
                "filterMinNumPeptidesPerSample",
                "removePeptidesWithOnlyOneReplicate",
                "peptideMissingValueImputation",
                "removeProteinsWithOnlyOneReplicate"
            )

            setCoreUtilisation <- function(config_list, function_name) {
                if (!function_name %in% names(config_list)) {
                    config_list[[function_name]] <- list()
                }
                config_list[[function_name]][["core_utilisation"]] <- core_utilisation

                config_list
            }

            for (x in list_of_multithreaded_functions) {
                config_list <- setCoreUtilisation(config_list, x)
            }

            config_list[["globalParameters"]][["plots_format"]] <- str_split(config_list[["globalParameters"]][["plots_format"]], ",")[[1]]
        }
    }

    getConfigValue <- function(config_list, section, value) {
        config_list[[section]][[value]]
    }

    setConfigValueAsNumeric <- function(config_list, section, value) {
        config_list[[section]][[value]] <- as.numeric(config_list[[section]][[value]])
        config_list
    }

    if ("srlQvalueProteotypicPeptideClean" %in% names(config_list)) {
        # Parse input_matrix_column_ids - handle both formats:
        # 1. Already a vector (from fresh workflow)
        # 2. String like 'c("Run", "Precursor.Id", ...)' (from saved config.ini)
        raw_value <- config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]]

        if (is.character(raw_value) && length(raw_value) == 1 && grepl("^c\\(", raw_value)) {
            # Format: 'c("Run", "Precursor.Id", ...)' - parse as R code
            config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]] <-
                eval(parse(text = raw_value))
        } else if (is.character(raw_value) && length(raw_value) == 1) {
            # Format: 'Run, Precursor.Id, ...' - split by comma
            raw_split <- str_split(raw_value, ",")[[1]]
            config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]] <- trimws(raw_split)
        }
        # else: already a vector, keep as is

        # Remove any empty strings
        config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]] <-
            config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]][
                config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]] != ""
            ]

        print(paste0(
            "Read srlQvalueProteotypicPeptideClean: input_matrix_column_ids = ",
            paste0(config_list[["srlQvalueProteotypicPeptideClean"]][["input_matrix_column_ids"]],
                collapse = ", "
            )
        ))

        config_list <- setConfigValueAsNumeric(
            config_list,
            "srlQvalueProteotypicPeptideClean",
            "qvalue_threshold"
        )
        config_list <- setConfigValueAsNumeric(
            config_list,
            "srlQvalueProteotypicPeptideClean",
            "global_qvalue_threshold"
        )
        if ("global_pg_qvalue_threshold" %in%
            names(config_list[["srlQvalueProteotypicPeptideClean"]])) {
            config_list <- setConfigValueAsNumeric(
                config_list,
                "srlQvalueProteotypicPeptideClean",
                "global_pg_qvalue_threshold"
            )
        } else {
            config_list[["srlQvalueProteotypicPeptideClean"]][["global_pg_qvalue_threshold"]] <- 0.01
        }
        config_list <- setConfigValueAsNumeric(
            config_list,
            "srlQvalueProteotypicPeptideClean",
            "choose_only_proteotypic_peptide"
        )
        for (exclusion_flag in c("exclude_decoys", "exclude_contaminants")) {
            if (exclusion_flag %in%
                names(config_list[["srlQvalueProteotypicPeptideClean"]])) {
                config_list <- setConfigValueAsNumeric(
                    config_list,
                    "srlQvalueProteotypicPeptideClean",
                    exclusion_flag
                )
            }
        }
    }


    if ("peptideIntensityFiltering" %in% names(config_list)) {
        numeric_peptide_intensity_fields <- c(
            "peptides_intensity_cutoff_percentile",
            "peptides_proportion_of_samples_below_cutoff",
            "min_reps_per_group",
            "min_groups"
        )
        for (field_name in numeric_peptide_intensity_fields) {
            if (field_name %in% names(config_list[["peptideIntensityFiltering"]])) {
                config_list <- setConfigValueAsNumeric(
                    config_list,
                    "peptideIntensityFiltering",
                    field_name
                )
            }
        }
        if ("strict_mode" %in% names(config_list[["peptideIntensityFiltering"]])) {
            strict_value <- config_list[["peptideIntensityFiltering"]][["strict_mode"]]
            config_list[["peptideIntensityFiltering"]][["strict_mode"]] <-
                tolower(trimws(as.character(strict_value))) %in%
                    c("1", "true", "t", "yes", "y")
        }
    }


    if ("filterMinNumPeptidesPerProtein" %in% names(config_list)) {
        config_list <- setConfigValueAsNumeric(
            config_list,
            "filterMinNumPeptidesPerProtein",
            "peptides_per_protein_cutoff"
        )
        config_list <- setConfigValueAsNumeric(
            config_list,
            "filterMinNumPeptidesPerProtein",
            "peptidoforms_per_protein_cutoff"
        )
        # config_list <- setConfigValueAsNumeric(config_list
        #                                        , ""
        #                                        , "")
    }

    if ("filterMinNumPeptidesPerSample" %in% names(config_list)) {
        config_list <- setConfigValueAsNumeric(
            config_list,
            "filterMinNumPeptidesPerSample",
            "peptides_per_sample_cutoff"
        )

        if (!"inclusion_list" %in% names(config_list[["filterMinNumPeptidesPerSample"]])) {
            config_list[["filterMinNumPeptidesPerSample"]][["inclusion_list"]] <- ""
        }

        config_list[["filterMinNumPeptidesPerSample"]][["inclusion_list"]] <- str_split(config_list[["filterMinNumPeptidesPerSample"]][["inclusion_list"]], ",")[[1]]
    }

    if ("peptideMissingValueImputation" %in% names(config_list)) {
        config_list <- setConfigValueAsNumeric(
            config_list,
            "peptideMissingValueImputation",
            "proportion_missing_values"
        )
    }

    if ("removeRowsWithMissingValuesPercent" %in% names(config_list)) {
        config_list <- setConfigValueAsNumeric(
            config_list,
            "removeRowsWithMissingValuesPercent",
            "groupwise_percentage_cutoff"
        )

        config_list <- setConfigValueAsNumeric(
            config_list,
            "removeRowsWithMissingValuesPercent",
            "max_groups_percentage_cutoff"
        )

        config_list <- setConfigValueAsNumeric(
            config_list,
            "removeRowsWithMissingValuesPercent",
            "proteins_intensity_cutoff_percentile"
        )
    }


    if ("ruvIII_C_Varying" %in% names(config_list)) {
        config_list <- setConfigValueAsNumeric(
            config_list,
            "ruvIII_C_Varying",
            "ruv_number_k"
        )
    }

    if ("plotRle" %in% names(config_list)) {
        config_list[["plotRle"]][["yaxis_limit"]] <- str_split(config_list[["plotRle"]][["yaxis_limit"]], ",")[[1]] |>
            purrr::map_dbl(\(x) as.numeric(x))

        print(paste0(
            "Read plotRle: yaxis_limit = ",
            paste0(config_list[["plotRle"]][["yaxis_limit"]], collapse = ", ")
        ))
    }

    if ("deAnalysisParameters" %in% names(config_list)) {
        message("DEBUG66: Processing deAnalysisParameters section")
        message(sprintf(
            "DEBUG66: deAnalysisParameters keys = %s",
            paste(names(config_list[["deAnalysisParameters"]]), collapse = ", ")
        ))

        # Handle plots_format as array (only if it exists)
        if ("plots_format" %in% names(config_list[["deAnalysisParameters"]])) {
            plots_val <- config_list[["deAnalysisParameters"]][["plots_format"]]
            if (!is.null(plots_val) && nzchar(plots_val)) {
                config_list[["deAnalysisParameters"]][["plots_format"]] <-
                    str_split(plots_val, ",")[[1]]
            }
        }

        # Add new lfc_cutoff parameter
        config_list[["deAnalysisParameters"]][["lfc_cutoff"]] <- FALSE

        # Modify treat_lfc_cutoff to use ifelse
        config_list[["deAnalysisParameters"]][["treat_lfc_cutoff"]] <-
            ifelse(config_list[["deAnalysisParameters"]][["lfc_cutoff"]], log2(1.5), 0)

        # Handle args_group_pattern - remove quotes and fix escaping
        if ("args_group_pattern" %in% names(config_list[["deAnalysisParameters"]])) {
            config_list[["deAnalysisParameters"]][["args_group_pattern"]] <-
                gsub('^"|"$', "", config_list[["deAnalysisParameters"]][["args_group_pattern"]]) |>
                gsub(pattern = "\\\\", replacement = "\\")
        }

        # Convert numeric parameters (only if they exist)
        if ("da_q_val_thresh" %in% names(config_list[["deAnalysisParameters"]])) {
            config_list <- setConfigValueAsNumeric(
                config_list,
                "deAnalysisParameters",
                "da_q_val_thresh"
            )
        }

        # Convert boolean parameters (only if they exist)
        if ("eBayes_trend" %in% names(config_list[["deAnalysisParameters"]])) {
            config_list[["deAnalysisParameters"]][["eBayes_trend"]] <-
                tolower(config_list[["deAnalysisParameters"]][["eBayes_trend"]]) == "true"
        }
        if ("eBayes_robust" %in% names(config_list[["deAnalysisParameters"]])) {
            config_list[["deAnalysisParameters"]][["eBayes_robust"]] <-
                tolower(config_list[["deAnalysisParameters"]][["eBayes_robust"]]) == "true"
        }

        message("DEBUG66: deAnalysisParameters processing complete")
        print(paste0(
            "Read deAnalysisParameters: formula_string = ",
            config_list[["deAnalysisParameters"]][["formula_string"]]
        ))
    }

    message("DEBUG66: readConfigFile() returning successfully")
    config_list
}

#' @title Read the config file and specify the section and or parameter to update the object
#' @description Read the config file and specify the section and or parameter to update the object
#' @param theObject The object to be updated
#' @param file The file path to the config file
#' @param section The section to be updated
#' @param value The parameter value to be updated
#' @export
readConfigFileSection <- function(
  theObject,
  file = file.path(source_dir, "config.ini"),
  function_name,
  parameter_name = NULL
) {
    config_list <- readConfigFile(
        file = file,
        file_type = "ini"
    )

    if (is.null(parameter_name)) {
        theObject@args[[function_name]] <- config_list[[function_name]]
    } else {
        theObject@args[[function_name]][[parameter_name]] <- config_list[[function_name]][[parameter_name]]
    }

    theObject
}

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

formatStudyParameterContrasts <- function(contrasts_tbl) {
    if (is.null(contrasts_tbl) ||
        !(is.data.frame(contrasts_tbl) || tibble::is_tibble(contrasts_tbl)) ||
        nrow(contrasts_tbl) == 0) {
        return(character())
    }

    raw_col <- intersect(c("contrasts", "contrast_string", "contrast"), colnames(contrasts_tbl))
    raw_values <- if (length(raw_col) > 0) {
        as.character(contrasts_tbl[[raw_col[[1]]]])
    } else {
        rep(NA_character_, nrow(contrasts_tbl))
    }

    friendly_values <- if ("friendly_names" %in% colnames(contrasts_tbl)) {
        as.character(contrasts_tbl[["friendly_names"]])
    } else {
        rep(NA_character_, nrow(contrasts_tbl))
    }

    full_values <- if ("full_format" %in% colnames(contrasts_tbl)) {
        as.character(contrasts_tbl[["full_format"]])
    } else {
        rep(NA_character_, nrow(contrasts_tbl))
    }

    vapply(seq_len(nrow(contrasts_tbl)), function(i) {
        full_value <- full_values[[i]]
        raw_value <- raw_values[[i]]
        friendly_value <- friendly_values[[i]]

        if (!is.na(full_value) && nzchar(full_value)) {
            return(paste("  ", full_value))
        }

        if (!is.na(raw_value) && nzchar(raw_value)) {
            if (!is.na(friendly_value) && nzchar(friendly_value) && !identical(friendly_value, raw_value)) {
                return(paste0("  ", friendly_value, " = ", raw_value))
            }
            return(paste("  ", raw_value))
        }

        "  [No contrast expression available]"
    }, character(1))
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

        contrasts_info <- tryCatch(
            formatStudyParameterContrasts(contrasts_tbl),
            error = function(e) {
                paste("  [Error extracting contrasts:", e$message, "]")
            }
        )
        contrasts_lines <- c(contrasts_lines, contrasts_info)

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

