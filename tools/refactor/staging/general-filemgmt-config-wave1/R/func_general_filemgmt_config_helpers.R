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
        config_list <- setConfigValueAsNumeric(
            config_list,
            "srlQvalueProteotypicPeptideClean",
            "choose_only_proteotypic_peptide"
        )
    }


    if ("peptideIntensityFiltering" %in% names(config_list)) {
        config_list <- setConfigValueAsNumeric(
            config_list,
            "peptideIntensityFiltering",
            "peptides_intensity_cutoff_percentile"
        )
        config_list <- setConfigValueAsNumeric(
            config_list,
            "peptideIntensityFiltering",
            "peptides_proportion_of_samples_below_cutoff"
        )
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

