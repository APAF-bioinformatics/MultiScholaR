#' Metabolite Assay Data S4 Class
#'
#' An S4 class to store and manage multiple metabolomics quantitative datasets
#' derived from different platforms or ionization modes, along with associated
#' experimental design and metadata.
#'
#' @slot metabolite_data A named list of data frames. Each data frame contains
#'   quantitative data for one assay (e.g., LCMS-Pos). Rows should represent
#'   metabolites/features, and columns should represent samples. Each data frame
#'   MUST contain the column specified by `metabolite_id_column`. The names of
#'   the list elements should describe the assay source (e.g., "LCMS_Pos").
#' @slot metabolite_id_column Character string. The name of the column within
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
#'   the `metabolite_id_column` or `annotation_id_column`. Set to `NA_character_` or `""`
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
#' @exportClass MetaboliteAssayData
MetaboliteAssayData <- setClass("MetaboliteAssayData",
    slots = c(
        metabolite_data = "list",
        metabolite_id_column = "character",
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
        metabolite_data = list(),
        metabolite_id_column = "database_identifier",
        annotation_id_column = "metabolite_identification",
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

        # --- Basic Type Checks ---
        if (!is.list(object@metabolite_data)) {
            errors <- c(errors, "`metabolite_data` must be a list.")
        }
        if (length(object@metabolite_data) > 0 && !all(sapply(object@metabolite_data, is.data.frame))) {
            errors <- c(errors, "All elements in `metabolite_data` must be data frames.")
        }
        if (!is.character(object@metabolite_id_column) || length(object@metabolite_id_column) != 1) {
            errors <- c(errors, "`metabolite_id_column` must be a single character string.")
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
        if (length(object@metabolite_data) > 0) {
            # Check primary metabolite ID column exists in all assays
            metabolite_col_exists <- sapply(object@metabolite_data, function(df) {
                object@metabolite_id_column %in% colnames(df)
            })
            if (!all(metabolite_col_exists)) {
                missing_in <- names(object@metabolite_data)[!metabolite_col_exists]
                errors <- c(errors, paste0("Primary `metabolite_id_column` ('", object@metabolite_id_column, "') not found in assay(s): ", paste(missing_in, collapse = ", ")))
            }

            # Check sample ID column exists in design matrix
            if (!(object@sample_id %in% colnames(object@design_matrix))) {
                errors <- c(errors, paste0("`sample_id` column ('", object@sample_id, "') not found in `design_matrix`."))
            } else {
                # Check consistency of samples between assays and design matrix
                samples_in_design <- tryCatch(
                    object@design_matrix |> dplyr::pull(!!rlang::sym(object@sample_id)) |> as.character() |> unique() |> sort(),
                    error = function(e) character(0) # Handle case where column doesn't exist
                )

                if (length(samples_in_design) == 0 && !(object@sample_id %in% colnames(object@design_matrix))) {
                    # Error already captured above
                } else {
                    first_assay_samples <- NULL
                    for (assay_name in names(object@metabolite_data)) {
                        assay_df <- object@metabolite_data[[assay_name]]
                        assay_samples <- setdiff(colnames(assay_df), object@metabolite_id_column) |> sort()

                        if (is.null(first_assay_samples)) {
                            first_assay_samples <- assay_samples
                        }

                        # Check all assays have the same sample columns
                        if (!identical(assay_samples, first_assay_samples)) {
                            errors <- c(errors, paste0("Sample columns differ between assays. Check assay '", assay_name, "'."))
                            # Stop checking further assays for sample consistency after first mismatch
                            break
                        }
                        # Check assay samples match design matrix samples (only need to check once if all assays are identical)
                        if (assay_name == names(object@metabolite_data)[1]) {
                            if (!identical(assay_samples, samples_in_design)) {
                                errors <- c(errors, paste0("Sample columns in assays do not exactly match unique sample IDs ('", object@sample_id, "') in `design_matrix`."))
                                # Stop checking further assays for sample consistency after first mismatch
                                break
                            }
                        }
                    }
                }
            }
        }

        # --- Final Check ---
        if (length(errors) == 0) TRUE else errors
    }
)
#' @export MetaboliteAssayData


#' Create MetaboliteAssayData Object
#'
#' Constructor function for the MetaboliteAssayData class.
#'
#' @param metabolite_data Named list of data frames (assays).
#' @param design_matrix Experimental design data frame.
#' @param metabolite_id_column Name of the **primary feature ID** column within assays (e.g., `"database_identifier"`).
#' @param annotation_id_column Name of the **annotation ID** column (e.g., `"metabolite_identification"`).
#' @param sample_id Name of the sample ID column in design_matrix and assays.
#' @param group_id Name of the group column in design_matrix.
#' @param technical_replicate_id Name of the technical replicate column in design_matrix (use NA_character_ if none).
#' @param database_identifier_type Type of identifier in the `annotation_id_column` (e.g., `"Mixed_CHEBI_Unknown"`).
#' @param internal_standard_regex Regex to identify internal standards. Use `NA_character_` or `""` if none.
#' @param args List of arguments (e.g., from config).
#'
#' @return A MetaboliteAssayData object.
#' @export
#' @examples
#' \dontrun{
#' # Assuming lcms_pos_df, lcms_neg_df, gcms_df are data frames
#' # with 'Metabolite' as ID column and samples as other columns
#' # Assuming design_df has 'SampleID', 'Group', 'Replicate' columns
#' assays_list <- list(
#'     LCMS_Pos = lcms_pos_df,
#'     LCMS_Neg = lcms_neg_df,
#'     GCMS = gcms_df
#' )
#' config <- list(...) # Your config list
#'
#' met_assay_obj <- createMetaboliteAssayData(
#'     metabolite_data = assays_list,
#'     design_matrix = design_df,
#'     metabolite_id_column = "Metabolite",
#'     sample_id = "SampleID",
#'     group_id = "Group",
#'     technical_replicate_id = "Replicate",
#'     database_identifier_type = "InternalName",
#'     internal_standard_regex = "^IS_",
#'     args = config
#' )
#' }
createMetaboliteAssayData <- function(
    metabolite_data,
    design_matrix,
    metabolite_id_column = "database_identifier",
    annotation_id_column = "metabolite_identification",
    sample_id = "Sample_ID",
    group_id = "group",
    technical_replicate_id = NA_character_,
    database_identifier_type = "Unknown",
    internal_standard_regex = NA_character_,
    args = list()) {
    # Perform basic checks before creating the object
    stopifnot(is.list(metabolite_data))
    stopifnot(all(sapply(metabolite_data, is.data.frame)))
    stopifnot(is.data.frame(design_matrix))
    # Add more checks as needed...

    obj <- new("MetaboliteAssayData",
        metabolite_data = metabolite_data,
        metabolite_id_column = metabolite_id_column,
        annotation_id_column = annotation_id_column,
        database_identifier_type = database_identifier_type,
        internal_standard_regex = internal_standard_regex,
        design_matrix = design_matrix,
        sample_id = sample_id,
        group_id = group_id,
        technical_replicate_id = technical_replicate_id,
        args = args
    )
    # Validity check is automatically called by 'new'
    return(obj)
}
