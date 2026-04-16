# ==========================================
# Content from lipidVsSamplesS4Objects.R
# ==========================================
#' Lipid Assay Data S4 Class
#'
#' An S4 class to store and manage multiple lipidomics quantitative datasets
#' derived from different platforms or ionization modes, along with associated
#' experimental design and metadata.
#'
#' @slot lipid_data A named list of data frames. Each data frame contains
#'   quantitative data for one assay (e.g., LCMS-Pos). Rows should represent
#'   lipids/features, and columns should represent samples. Each data frame
#'   MUST contain the column specified by `lipid_id_column`. The names of
#'   the list elements should describe the assay source (e.g., "LCMS_Pos").
#' @slot lipid_id_column Character string. The name of the column within
#'   each assay data frame that contains the **primary feature identifier** (used
#'   for uniquely identifying rows in quantitative analysis).
#' @slot annotation_id_column Character string. The name of the column (often
#'   in feature metadata, but could be in assay tables) that contains the
#'   **biological or database annotation identifier** (e.g., HMDB ID, KEGG ID,
#'   chemical name, CHEBI ID). This is the ID typically used for downstream
#'   biological mapping.
#' @slot database_identifier_type Character string. Describes the type or format
#'   of identifiers found in the column specified by `annotation_id_column`
#'   (e.g., "HMDB", "KEGG", "CHEBI", "Mixed_CHEBI_Unknown", "InternalName").
#' @slot internal_standard_regex Character string. A regular expression used to
#'   identify features that are internal standards based on their identifier in
#'   the `lipid_id_column` or `annotation_id_column`. Set to `NA_character_` or `""`
#'   if not applicable or no internal standards are used.
#' @slot design_matrix A data frame containing the experimental design. Must
#'   include the column specified by `sample_id`.
#' @slot sample_id Character string. The name of the column in `design_matrix`
#'   that contains unique sample identifiers. These identifiers must correspond
#'   to the sample column names in the assay data frames.
#' @slot group_id Character string. The name of the column in `design_matrix`
#'   that defines the primary experimental groups or conditions.
#' @slot technical_replicate_id Character string. The name of the column in
#'   `design_matrix` that identifies technical replicates, if applicable. Use
#'   `NA_character_` if there are no technical replicates.
#' @slot args A list, typically populated from a configuration file, holding
#'   parameters used during processing.
#'
#' @importFrom methods setClass slotNames slot new
#' @importFrom dplyr pull distinct
#' @importFrom rlang sym !!
#' @exportClass LipidomicsAssayData
LipidomicsAssayData <- setClass("LipidomicsAssayData",
    slots = c(
        lipid_data = "list",
        lipid_id_column = "character",
        annotation_id_column = "character",
        database_identifier_type = "character",
        internal_standard_regex = "character",
        design_matrix = "data.frame",
        sample_id = "character",
        group_id = "character",
        technical_replicate_id = "character",
        args = "list"
    ),
    prototype = list(
        lipid_data = list(),
        lipid_id_column = "database_identifier",
        annotation_id_column = "lipid_identification",
        database_identifier_type = "Unknown",
        internal_standard_regex = NA_character_,
        design_matrix = data.frame(),
        sample_id = "Sample_ID",
        group_id = "group",
        technical_replicate_id = NA_character_,
        args = list()
    ),
    validity = function(object) {
        errors <- character()
        # --- Get required info ---
        sample_id_col <- object@sample_id
        lipid_id_col <- object@lipid_id_column
        design_matrix <- object@design_matrix
        lipid_data <- object@lipid_data

        # --- Basic slot type checks (as before) ---
        if (!is.list(lipid_data)) {
            errors <- c(errors, "`lipid_data` must be a list.")
        } else if (length(lipid_data) > 0 && !all(sapply(lipid_data, is.data.frame))) {
            errors <- c(errors, "All elements in `lipid_data` must be data frames.")
        }
        if (!is.character(object@lipid_id_column) || length(object@lipid_id_column) != 1) {
            errors <- c(errors, "`lipid_id_column` must be a single character string.")
        }
        if (!is.character(object@annotation_id_column) || length(object@annotation_id_column) != 1) {
            errors <- c(errors, "`annotation_id_column` must be a single character string.")
        }
        if (!is.character(object@database_identifier_type) || length(object@database_identifier_type) != 1) {
            errors <- c(errors, "`database_identifier_type` must be a single character string.")
        }
        if (!is.character(object@internal_standard_regex) || length(object@internal_standard_regex) != 1) {
            errors <- c(errors, "`internal_standard_regex` must be a single character string (can be NA_character_).")
        }
        if (!is.data.frame(object@design_matrix)) {
            errors <- c(errors, "`design_matrix` must be a data frame.")
        }
        if (!is.character(object@sample_id) || length(object@sample_id) != 1) {
            errors <- c(errors, "`sample_id` must be a single character string.")
        }
        if (!is.character(object@group_id) || length(object@group_id) != 1) {
            errors <- c(errors, "`group_id` must be a single character string.")
        }
        if (!is.character(object@technical_replicate_id) || length(object@technical_replicate_id) != 1) {
            errors <- c(errors, "`technical_replicate_id` must be a single character string (can be NA_character_).")
        }
        if (!is.list(object@args)) {
            errors <- c(errors, "`args` must be a list.")
        }

        # --- Content Checks ---
        # Check design matrix first
        if (!is.data.frame(design_matrix)) {
            # Error already added by basic checks, but prevent further processing
        } else if (!sample_id_col %in% colnames(design_matrix)) {
            errors <- c(errors, paste0("`sample_id` column ('", sample_id_col, "') not found in `design_matrix`."))
        } else {
            # Get unique, sorted sample IDs from design matrix (ensure character)
            samples_in_design <- tryCatch(
                design_matrix[[sample_id_col]] |> as.character() |> unique() |> sort(),
                error = function(e) {
                    errors <- c(errors, "Error extracting sample IDs from design matrix.")
                    character(0)
                }
            )

            if (length(samples_in_design) == 0 && length(errors) == 0) {
                errors <- c(errors, "No valid sample IDs found in design matrix.")
            }

            # Proceed with assay checks only if design matrix looks okay so far
            if (length(lipid_data) > 0 && length(errors) == 0) {
                assay_names_vec <- names(lipid_data)
                if (is.null(assay_names_vec)) assay_names_vec <- paste0("Assay_", seq_along(lipid_data))
                names(lipid_data) <- assay_names_vec # Ensure the list is named for lapply output

                # Use lapply to check each assay and collect results/errors
                assay_check_results <- lapply(assay_names_vec, function(assay_name) {
                    assay_df <- lipid_data[[assay_name]]
                    assay_errors <- character()

                    # Check lipid ID column exists
                    if (!lipid_id_col %in% colnames(assay_df)) {
                        assay_errors <- c(assay_errors, paste0("Assay '", assay_name, "': `lipid_id_column` ('", lipid_id_col, "') not found."))
                    }

                    # Identify actual sample columns in the assay
                    assay_colnames <- colnames(assay_df)
                    actual_sample_cols_in_assay <- intersect(assay_colnames, samples_in_design) |> sort()

                    # Store results for later checks
                    list(
                        errors = assay_errors,
                        sample_cols = actual_sample_cols_in_assay
                    )
                })

                # Aggregate errors from individual assay checks
                all_assay_errors <- unlist(lapply(assay_check_results, `[[`, "errors"))
                errors <- c(errors, all_assay_errors)

                # Perform cross-assay consistency checks if no individual errors found yet
                if (length(errors) == 0 && length(assay_check_results) > 1) {
                    first_assay_samples <- assay_check_results[[1]]$sample_cols
                    consistency_check <- sapply(assay_check_results[-1], function(res) {
                        identical(res$sample_cols, first_assay_samples)
                    })
                    if (!all(consistency_check)) {
                        mismatched_assays <- assay_names_vec[c(FALSE, !consistency_check)] # Get names of inconsistent assays
                        errors <- c(errors, paste0("Actual sample columns differ between assays. First mismatch found in: ", mismatched_assays[1]))
                    }
                }

                # Perform comparison with design matrix if no errors found yet
                if (length(errors) == 0 && length(assay_check_results) >= 1) {
                    first_assay_samples <- assay_check_results[[1]]$sample_cols # Get samples from first (or only) assay
                    if (!identical(first_assay_samples, samples_in_design)) {
                        errors <- c(errors, paste0("Sample columns in assays do not exactly match unique sample IDs ('", sample_id_col, "') in `design_matrix`."))
                        # Add more detail:
                        missing_in_assay <- setdiff(samples_in_design, first_assay_samples)
                        extra_in_assay <- setdiff(first_assay_samples, samples_in_design)
                        if (length(missing_in_assay) > 0) errors <- c(errors, paste0("   Samples in design_matrix missing from assay columns: ", paste(utils::head(missing_in_assay, 10), collapse = ", "), ifelse(length(missing_in_assay) > 10, "...", "")))
                        if (length(extra_in_assay) > 0) errors <- c(errors, paste0("   Sample columns in assay not found in design_matrix: ", paste(utils::head(extra_in_assay, 10), collapse = ", "), ifelse(length(extra_in_assay) > 10, "...", "")))
                    }
                }
            }
        }

        # --- Final Check ---
        if (length(errors) == 0) TRUE else errors
    }
)

