#' Extract Metabolomics-Specific Parameters from MetaboliteAssayData S4 Object
#'
#' @description Helper function to extract metabolomics-specific workflow parameters
#'              from the MetaboliteAssayData S4 object's @args slot and other slots.
#'
#' @param metab_s4 A MetaboliteAssayData S4 object
#'
#' @return A list containing:
#' \itemize{
#'   \item assay_info: List with assay names, counts, metabolite totals
#'   \item itsd_normalization: ITSD normalization parameters (if applied)
#'   \item log_transformation: Log transformation parameters
#'   \item normalisation_method: Between-sample normalization method
#'   \item ruv_params: RUV-III parameters (if applied)
#'   \item metadata: Slot metadata (metabolite_id_column, sample_id, group_id, etc.)
#' }
#'
#' @keywords internal
extractMetabS4Params <- function(metab_s4) {
    if (!inherits(metab_s4, "MetaboliteAssayData")) {
        warning("extractMetabS4Params: Object is not MetaboliteAssayData, returning empty list")
        return(list())
    }

    cat("METAB PARAMS: Extracting parameters from MetaboliteAssayData S4 object\n")

    result <- list()

    # -- Assay Information --
    tryCatch(
        {
            assay_data <- metab_s4@metabolite_data
            assay_names <- names(assay_data)

            # Get metabolite counts per assay
            metabolite_counts <- sapply(assay_names, function(name) {
                assay_df <- assay_data[[name]]
                if (is.data.frame(assay_df) && metab_s4@metabolite_id_column %in% names(assay_df)) {
                    length(unique(assay_df[[metab_s4@metabolite_id_column]]))
                } else if (is.data.frame(assay_df)) {
                    nrow(assay_df)
                } else {
                    NA_integer_
                }
            })

            result$assay_info <- list(
                num_assays = length(assay_names),
                assay_names = paste(assay_names, collapse = ", "),
                metabolites_per_assay = metabolite_counts,
                total_metabolites = sum(metabolite_counts, na.rm = TRUE)
            )
            cat(sprintf(
                "METAB PARAMS: Found %d assays: %s\n",
                length(assay_names), paste(assay_names, collapse = ", ")
            ))
        },
        error = function(e) {
            cat(sprintf("METAB PARAMS: Error extracting assay info: %s\n", e$message))
            result$assay_info <- list(num_assays = NA, assay_names = NA)
        }
    )

    # -- Extract @args parameters --
    s4_args <- tryCatch(
        {
            if (!is.null(metab_s4@args)) metab_s4@args else list()
        },
        error = function(e) {
            list()
        }
    )

    # -- ITSD Normalization --
    tryCatch(
        {
            if (!is.null(s4_args$ITSDNormalization)) {
                itsd <- s4_args$ITSDNormalization
                result$itsd_normalization <- list(
                    applied = isTRUE(itsd$applied),
                    method_type = itsd$method_type %||% NA_character_,
                    aggregation = itsd$itsd_aggregation %||% NA_character_,
                    pattern_columns = if (!is.null(itsd$itsd_pattern_columns)) {
                        paste(itsd$itsd_pattern_columns, collapse = ", ")
                    } else {
                        NA_character_
                    },
                    removed_after_norm = isTRUE(itsd$removed_itsd),
                    timestamp = if (!is.null(itsd$timestamp)) {
                        format(itsd$timestamp, "%Y-%m-%d %H:%M:%S")
                    } else {
                        NA_character_
                    }
                )

                # Extract per-assay ITSD feature names (new)
                if (!is.null(itsd$itsd_features_per_assay) && is.list(itsd$itsd_features_per_assay)) {
                    result$itsd_normalization$features_per_assay <- itsd$itsd_features_per_assay
                    result$itsd_normalization$counts_per_assay <- itsd$itsd_counts_per_assay
                    cat(sprintf("METAB PARAMS: Found ITSD features for %d assays\n", length(itsd$itsd_features_per_assay)))
                }

                cat("METAB PARAMS: Extracted ITSD normalization parameters\n")
            } else {
                result$itsd_normalization <- list(applied = FALSE)
                cat("METAB PARAMS: No ITSD normalization found in @args\n")
            }
        },
        error = function(e) {
            cat(sprintf("METAB PARAMS: Error extracting ITSD params: %s\n", e$message))
            result$itsd_normalization <- list(applied = FALSE, error = e$message)
        }
    )

    # -- Log Transformation --
    tryCatch(
        {
            result$log_transformation <- list(
                applied = isTRUE(s4_args$log_transformed),
                offset = s4_args$log_transform_offset %||% NA_real_
            )
            cat(sprintf(
                "METAB PARAMS: Log transformation applied: %s\n",
                isTRUE(s4_args$log_transformed)
            ))
        },
        error = function(e) {
            result$log_transformation <- list(applied = FALSE)
        }
    )

    # -- Between-Sample Normalization --
    tryCatch(
        {
            result$normalisation_method <- s4_args$normalisation_method %||% NA_character_
            cat(sprintf(
                "METAB PARAMS: Normalization method: %s\n",
                result$normalisation_method %||% "N/A"
            ))
        },
        error = function(e) {
            result$normalisation_method <- NA_character_
        }
    )

    # -- RUV-III Parameters --
    tryCatch(
        {
            if (!is.null(s4_args$ruv_number_k) || !is.null(s4_args$ruv_grouping_variable)) {
                # Check if ruv_number_k is per-assay (list) or single value
                if (is.list(s4_args$ruv_number_k)) {
                    # Per-assay RUV (metabolomics with multiple assays)
                    result$ruv_params <- list(
                        applied = TRUE,
                        per_assay = TRUE,
                        grouping_variable = s4_args$ruv_grouping_variable %||% NA_character_,
                        k_per_assay = s4_args$ruv_number_k,
                        ctrl_per_assay = if (!is.null(s4_args$ctrl) && is.list(s4_args$ctrl)) {
                            lapply(s4_args$ctrl, function(x) sum(x, na.rm = TRUE))
                        } else {
                            NULL
                        }
                    )
                    cat(sprintf(
                        "METAB PARAMS: RUV-III applied per-assay with k values: %s\n",
                        paste(names(s4_args$ruv_number_k), "=", unlist(s4_args$ruv_number_k), collapse = ", ")
                    ))
                } else {
                    # Single RUV value (proteomics or single-assay metabolomics)
                    result$ruv_params <- list(
                        applied = TRUE,
                        per_assay = FALSE,
                        grouping_variable = s4_args$ruv_grouping_variable %||% NA_character_,
                        number_k = s4_args$ruv_number_k %||% NA_integer_,
                        ctrl_count = if (!is.null(s4_args$ctrl)) sum(s4_args$ctrl, na.rm = TRUE) else NA_integer_
                    )
                    cat(sprintf(
                        "METAB PARAMS: RUV-III applied with k=%s\n",
                        s4_args$ruv_number_k %||% "N/A"
                    ))
                }
            } else {
                result$ruv_params <- list(applied = FALSE)
                cat("METAB PARAMS: No RUV-III parameters found\n")
            }
        },
        error = function(e) {
            result$ruv_params <- list(applied = FALSE)
        }
    )

    # -- Slot Metadata --
    tryCatch(
        {
            result$metadata <- list(
                metabolite_id_column = metab_s4@metabolite_id_column,
                annotation_id_column = metab_s4@annotation_id_column,
                database_identifier_type = metab_s4@database_identifier_type,
                sample_id = metab_s4@sample_id,
                group_id = metab_s4@group_id,
                internal_standard_regex = if (!is.na(metab_s4@internal_standard_regex)) {
                    metab_s4@internal_standard_regex
                } else {
                    NA_character_
                }
            )

            # Design matrix info
            if (!is.null(metab_s4@design_matrix) && is.data.frame(metab_s4@design_matrix)) {
                dm <- metab_s4@design_matrix
                result$metadata$n_samples <- nrow(dm)
                result$metadata$n_groups <- length(unique(dm[[metab_s4@group_id]]))
                result$metadata$group_names <- paste(unique(dm[[metab_s4@group_id]]), collapse = ", ")
            }
            cat("METAB PARAMS: Extracted slot metadata\n")
        },
        error = function(e) {
            cat(sprintf("METAB PARAMS: Error extracting metadata: %s\n", e$message))
            result$metadata <- list()
        }
    )

    # -- Pass through any other @args sections --
    result$raw_args <- s4_args

    cat("METAB PARAMS: Extraction complete. Found %d parameter groups\n", length(result))
    return(result)
}

#' Extract Lipidomics-Specific Parameters from LipidomicsAssayData S4 Object
#'
#' @description Helper function to extract lipidomics-specific workflow parameters
#'              from the LipidomicsAssayData S4 object's @args slot and other slots.
#'
#' @param lipid_s4 A LipidomicsAssayData S4 object
#'
#' @return A list containing:
#' \itemize{
#'   \item assay_info: List with assay names, counts, lipid totals
#'   \item itsd_normalization: ITSD normalization parameters (if applied)
#'   \item log_transformation: Log transformation parameters
#'   \item normalisation_method: Between-sample normalization method
#'   \item ruv_params: RUV-III parameters (if applied)
#'   \item metadata: Slot metadata (lipid_id_column, sample_id, group_id, etc.)
#' }
#'
#' @keywords internal
extractLipidS4Params <- function(lipid_s4) {
    if (!inherits(lipid_s4, "LipidomicsAssayData")) {
        warning("extractLipidS4Params: Object is not LipidomicsAssayData, returning empty list")
        return(list())
    }

    cat("LIPID PARAMS: Extracting parameters from LipidomicsAssayData S4 object\n")

    result <- list()

    # -- Assay Information --
    tryCatch(
        {
            assay_data <- lipid_s4@lipid_data
            assay_names <- names(assay_data)

            # Get lipid counts per assay
            lipid_counts <- sapply(assay_names, function(name) {
                assay_df <- assay_data[[name]]
                if (is.data.frame(assay_df) && lipid_s4@lipid_id_column %in% names(assay_df)) {
                    length(unique(assay_df[[lipid_s4@lipid_id_column]]))
                } else if (is.data.frame(assay_df)) {
                    nrow(assay_df)
                } else {
                    NA_integer_
                }
            })

            result$assay_info <- list(
                num_assays = length(assay_names),
                assay_names = paste(assay_names, collapse = ", "),
                lipids_per_assay = lipid_counts,
                total_lipids = sum(lipid_counts, na.rm = TRUE)
            )
            cat(sprintf(
                "LIPID PARAMS: Found %d assays: %s\n",
                length(assay_names), paste(assay_names, collapse = ", ")
            ))
        },
        error = function(e) {
            cat(sprintf("LIPID PARAMS: Error extracting assay info: %s\n", e$message))
            result$assay_info <- list(num_assays = NA, assay_names = NA)
        }
    )

    # -- Extract @args parameters --
    s4_args <- tryCatch(
        {
            if (!is.null(lipid_s4@args)) lipid_s4@args else list()
        },
        error = function(e) {
            list()
        }
    )

    # -- ITSD Normalization --
    tryCatch(
        {
            if (!is.null(s4_args$ITSDNormalization)) {
                itsd <- s4_args$ITSDNormalization
                result$itsd_normalization <- list(
                    applied = isTRUE(itsd$applied),
                    method_type = itsd$method_type %||% NA_character_,
                    aggregation = itsd$itsd_aggregation %||% NA_character_,
                    pattern_columns = if (!is.null(itsd$itsd_pattern_columns)) {
                        paste(itsd$itsd_pattern_columns, collapse = ", ")
                    } else {
                        NA_character_
                    },
                    removed_after_norm = isTRUE(itsd$removed_itsd),
                    timestamp = if (!is.null(itsd$timestamp)) {
                        format(itsd$timestamp, "%Y-%m-%d %H:%M:%S")
                    } else {
                        NA_character_
                    }
                )

                # Extract per-assay ITSD feature names
                if (!is.null(itsd$itsd_features_per_assay) && is.list(itsd$itsd_features_per_assay)) {
                    result$itsd_normalization$features_per_assay <- itsd$itsd_features_per_assay
                    result$itsd_normalization$counts_per_assay <- itsd$itsd_counts_per_assay
                    cat(sprintf("LIPID PARAMS: Found ITSD features for %d assays\n", length(itsd$itsd_features_per_assay)))
                }

                cat("LIPID PARAMS: Extracted ITSD normalization parameters\n")
            } else {
                result$itsd_normalization <- list(applied = FALSE)
                cat("LIPID PARAMS: No ITSD normalization found in @args\n")
            }
        },
        error = function(e) {
            cat(sprintf("LIPID PARAMS: Error extracting ITSD params: %s\n", e$message))
            result$itsd_normalization <- list(applied = FALSE, error = e$message)
        }
    )

    # -- Log Transformation --
    tryCatch(
        {
            result$log_transformation <- list(
                applied = isTRUE(s4_args$log_transformed),
                offset = s4_args$log_transform_offset %||% NA_real_
            )
            cat(sprintf(
                "LIPID PARAMS: Log transformation applied: %s\n",
                isTRUE(s4_args$log_transformed)
            ))
        },
        error = function(e) {
            result$log_transformation <- list(applied = FALSE)
        }
    )

    # -- Between-Sample Normalization --
    tryCatch(
        {
            result$normalisation_method <- s4_args$normalisation_method %||% NA_character_
            cat(sprintf(
                "LIPID PARAMS: Normalization method: %s\n",
                result$normalisation_method %||% "N/A"
            ))
        },
        error = function(e) {
            result$normalisation_method <- NA_character_
        }
    )

    # -- RUV-III Parameters --
    tryCatch(
        {
            if (!is.null(s4_args$ruv_number_k) || !is.null(s4_args$ruv_grouping_variable)) {
                # Check if ruv_number_k is per-assay (list) or single value
                if (is.list(s4_args$ruv_number_k)) {
                    # Per-assay RUV (lipidomics with multiple assays)
                    result$ruv_params <- list(
                        applied = TRUE,
                        per_assay = TRUE,
                        grouping_variable = s4_args$ruv_grouping_variable %||% NA_character_,
                        k_per_assay = s4_args$ruv_number_k,
                        ctrl_per_assay = if (!is.null(s4_args$ctrl) && is.list(s4_args$ctrl)) {
                            lapply(s4_args$ctrl, function(x) sum(x, na.rm = TRUE))
                        } else {
                            NULL
                        }
                    )
                    cat(sprintf(
                        "LIPID PARAMS: RUV-III applied per-assay with k values: %s\n",
                        paste(names(s4_args$ruv_number_k), "=", unlist(s4_args$ruv_number_k), collapse = ", ")
                    ))
                } else {
                    # Single RUV value
                    result$ruv_params <- list(
                        applied = TRUE,
                        per_assay = FALSE,
                        grouping_variable = s4_args$ruv_grouping_variable %||% NA_character_,
                        number_k = s4_args$ruv_number_k %||% NA_integer_,
                        ctrl_count = if (!is.null(s4_args$ctrl)) sum(s4_args$ctrl, na.rm = TRUE) else NA_integer_
                    )
                    cat(sprintf(
                        "LIPID PARAMS: RUV-III applied with k=%s\n",
                        s4_args$ruv_number_k %||% "N/A"
                    ))
                }
            } else {
                result$ruv_params <- list(applied = FALSE)
                cat("LIPID PARAMS: No RUV-III parameters found\n")
            }
        },
        error = function(e) {
            result$ruv_params <- list(applied = FALSE)
        }
    )

    # -- Slot Metadata --
    tryCatch(
        {
            result$metadata <- list(
                lipid_id_column = lipid_s4@lipid_id_column,
                annotation_id_column = lipid_s4@annotation_id_column,
                database_identifier_type = lipid_s4@database_identifier_type,
                sample_id = lipid_s4@sample_id,
                group_id = lipid_s4@group_id,
                internal_standard_regex = if (!is.na(lipid_s4@internal_standard_regex)) {
                    lipid_s4@internal_standard_regex
                } else {
                    NA_character_
                }
            )

            # Design matrix info
            if (!is.null(lipid_s4@design_matrix) && is.data.frame(lipid_s4@design_matrix)) {
                dm <- lipid_s4@design_matrix
                result$metadata$n_samples <- nrow(dm)
                result$metadata$n_groups <- length(unique(dm[[lipid_s4@group_id]]))
                result$metadata$group_names <- paste(unique(dm[[lipid_s4@group_id]]), collapse = ", ")
            }
            cat("LIPID PARAMS: Extracted slot metadata\n")
        },
        error = function(e) {
            cat(sprintf("LIPID PARAMS: Error extracting metadata: %s\n", e$message))
            result$metadata <- list()
        }
    )

    # -- Pass through any other @args sections --
    result$raw_args <- s4_args

    cat(sprintf("LIPID PARAMS: Extraction complete. Found %d parameter groups\n", length(result)))
    return(result)
}

