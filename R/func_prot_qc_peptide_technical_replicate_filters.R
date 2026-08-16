# ----------------------------------------------------------------------------
# removePeptidesWithOnlyOneReplicateHelper
# ----------------------------------------------------------------------------
.validatePeptideRunDesign <- function(input_table,
                                      design_matrix,
                                      input_sample_column,
                                      design_sample_column,
                                      design_group_column = NULL,
                                      caller = "peptide QC") {
  input_runs <- as.character(input_table[[input_sample_column]])
  design_runs <- as.character(design_matrix[[design_sample_column]])
  missing_input <- is.na(input_runs) | !nzchar(trimws(input_runs))
  missing_design <- is.na(design_runs) | !nzchar(trimws(design_runs))
  if (any(missing_input)) {
    stop(
      sprintf(
        "%s: input run identity is missing at row(s) %s.",
        caller,
        paste(utils::head(which(missing_input), 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (any(missing_design)) {
    stop(
      sprintf(
        "%s: design run identity is missing at row(s) %s.",
        caller,
        paste(utils::head(which(missing_design), 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  duplicate_design_runs <- unique(design_runs[
    duplicated(design_runs) | duplicated(design_runs, fromLast = TRUE)
  ])
  if (length(duplicate_design_runs) > 0L) {
    stop(
      sprintf(
        "%s: design must contain exactly one row per run; duplicate run ID(s): %s.",
        caller,
        paste(utils::head(duplicate_design_runs, 10L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  unmatched_input <- setdiff(unique(input_runs), design_runs)
  if (length(unmatched_input) > 0L) {
    stop(
      paste0(
        caller, ": input/design run identities do not match. ",
        "Unmapped input run(s): ",
        paste(utils::head(unmatched_input, 10L), collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  mapping <- data.frame(
    .qc_run = design_runs,
    stringsAsFactors = FALSE
  )
  if (!is.null(design_group_column)) {
    design_groups <- as.character(design_matrix[[design_group_column]])
    missing_group <- is.na(design_groups) | !nzchar(trimws(design_groups))
    if (any(missing_group)) {
      stop(
        sprintf(
          "%s: design replicate-group identity is missing at row(s) %s.",
          caller,
          paste(utils::head(which(missing_group), 5L), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    mapping$.qc_replicate_group <- design_groups
  }
  mapping
}

#---------------------------------------------------------------------------------------
#' @title Remove Peptides with Single Replicate
#' @description Retain peptide identities supported by at least two distinct
#'   runs in any declared replicate group. Per-group support is completed
#'   against the validated design, while singleton observations in other groups
#'   remain in the globally retained feature.
#' @export
#' @param input_table,samples_id_tbl,input_table_sample_id_column,sample_id_tbl_sample_id_column,replicate_group_column,protein_id_column,peptide_sequence_column,core_utilisation,minimum_distinct_runs,retention_policy,return_filter_result Runtime inputs used by this function; see the usage section for accepted values.
removePeptidesWithOnlyOneReplicateHelper <- function(input_table,
                                                      samples_id_tbl,
                                                      input_table_sample_id_column = Run,
                                                      sample_id_tbl_sample_id_column = ms_filename,
                                                      replicate_group_column = general_sample_info,
                                                      protein_id_column = Protein.Ids,
                                                      peptide_sequence_column = Stripped.Sequence,
                                                      core_utilisation,
                                                      minimum_distinct_runs = 2L,
                                                      retention_policy = "supported_in_any_group",
                                                      return_filter_result = FALSE) {

  input_sample_str <- resolvePeptideQcColumnArgument(substitute(input_table_sample_id_column), input_table_sample_id_column, names(input_table), environment())
  sample_tbl_sample_str <- resolvePeptideQcColumnArgument(substitute(sample_id_tbl_sample_id_column), sample_id_tbl_sample_id_column, names(samples_id_tbl), environment())
  replicate_group_str <- resolvePeptideQcColumnArgument(substitute(replicate_group_column), replicate_group_column, names(samples_id_tbl), environment())
  protein_id_str <- resolvePeptideQcColumnArgument(substitute(protein_id_column), protein_id_column, names(input_table), environment())
  peptide_seq_str <- resolvePeptideQcColumnArgument(substitute(peptide_sequence_column), peptide_sequence_column, names(input_table), environment())
  required_columns <- c(input_sample_str, protein_id_str, peptide_seq_str)
  missing_input_columns <- setdiff(required_columns, names(input_table))
  missing_design_columns <- setdiff(
    c(sample_tbl_sample_str, replicate_group_str),
    names(samples_id_tbl)
  )
  if (length(missing_input_columns) > 0L || length(missing_design_columns) > 0L) {
    stop(
      paste0(
        "removePeptidesWithOnlyOneReplicateHelper: missing required column(s): ",
        paste(c(missing_input_columns, missing_design_columns), collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }
  minimum_distinct_runs <- suppressWarnings(as.numeric(minimum_distinct_runs))
  if (length(minimum_distinct_runs) != 1L || is.na(minimum_distinct_runs) ||
      !is.finite(minimum_distinct_runs) || minimum_distinct_runs < 2L ||
      minimum_distinct_runs != floor(minimum_distinct_runs)) {
    stop(
      "removePeptidesWithOnlyOneReplicateHelper: `minimum_distinct_runs` must be an integer >= 2.",
      call. = FALSE
    )
  }
  minimum_distinct_runs <- as.integer(minimum_distinct_runs)
  if (!identical(retention_policy, "supported_in_any_group")) {
    stop(
      "removePeptidesWithOnlyOneReplicateHelper: only `supported_in_any_group` is currently supported.",
      call. = FALSE
    )
  }
  if (inherits(core_utilisation, "multidplyr_cluster")) {
    warning(
      paste0(
        "removePeptidesWithOnlyOneReplicateHelper: distinct-run support is ",
        "calculated deterministically in-process; the supplied cluster is not used."
      ),
      call. = FALSE
    )
  }

  mapping <- .validatePeptideRunDesign(
    input_table = input_table,
    design_matrix = samples_id_tbl,
    input_sample_column = input_sample_str,
    design_sample_column = sample_tbl_sample_str,
    design_group_column = replicate_group_str,
    caller = "removePeptidesWithOnlyOneReplicateHelper"
  )
  missing_feature <- is.na(input_table[[protein_id_str]]) |
    !nzchar(trimws(as.character(input_table[[protein_id_str]]))) |
    is.na(input_table[[peptide_seq_str]]) |
    !nzchar(trimws(as.character(input_table[[peptide_seq_str]])))
  if (any(missing_feature)) {
    stop(
      sprintf(
        "removePeptidesWithOnlyOneReplicateHelper: protein/peptide identity is missing at row(s) %s.",
        paste(utils::head(which(missing_feature), 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  observations <- input_table |>
    dplyr::transmute(
      .qc_run = as.character(.data[[input_sample_str]]),
      .qc_protein = as.character(.data[[protein_id_str]]),
      .qc_peptide = as.character(.data[[peptide_seq_str]])
    ) |>
    dplyr::distinct() |>
    dplyr::left_join(mapping, by = ".qc_run")
  features <- observations |>
    dplyr::distinct(.data$.qc_protein, .data$.qc_peptide)
  replicate_groups <- mapping |>
    dplyr::distinct(.data$.qc_replicate_group)
  observed_support <- observations |>
    dplyr::group_by(
      .data$.qc_replicate_group,
      .data$.qc_protein,
      .data$.qc_peptide
    ) |>
    dplyr::summarise(
      distinct_run_count = dplyr::n_distinct(.data$.qc_run),
      .groups = "drop"
    )

  support_table <- tidyr::expand_grid(features, replicate_groups) |>
    dplyr::left_join(
      observed_support,
      by = c(".qc_protein", ".qc_peptide", ".qc_replicate_group")
    ) |>
    dplyr::mutate(
      distinct_run_count = dplyr::coalesce(.data$distinct_run_count, 0L),
      required_distinct_runs = minimum_distinct_runs,
      group_supports_feature = .data$distinct_run_count >= minimum_distinct_runs
    ) |>
    dplyr::arrange(
      .data$.qc_protein,
      .data$.qc_peptide,
      .data$.qc_replicate_group
    )
  feature_decisions <- support_table |>
    dplyr::group_by(.data$.qc_protein, .data$.qc_peptide) |>
    dplyr::summarise(
      supporting_group_count = sum(.data$group_supports_feature),
      feature_retained = any(.data$group_supports_feature),
      .groups = "drop"
    )
  retained_features <- feature_decisions |>
    dplyr::filter(.data$feature_retained) |>
    dplyr::transmute(
      !!protein_id_str := .data$.qc_protein,
      !!peptide_seq_str := .data$.qc_peptide
    )
  filtered_data <- input_table |>
    dplyr::inner_join(
      retained_features,
      by = c(protein_id_str, peptide_seq_str)
    )
  removal_ledger <- feature_decisions |>
    dplyr::filter(!.data$feature_retained) |>
    dplyr::transmute(
      !!protein_id_str := .data$.qc_protein,
      !!peptide_seq_str := .data$.qc_peptide,
      failure_reason = "no_group_with_minimum_distinct_runs",
      supporting_group_count = .data$supporting_group_count,
      required_distinct_runs = minimum_distinct_runs
    )

  names(support_table)[names(support_table) == ".qc_protein"] <- protein_id_str
  names(support_table)[names(support_table) == ".qc_peptide"] <- peptide_seq_str
  names(support_table)[names(support_table) == ".qc_replicate_group"] <- replicate_group_str
  summary <- list(
    identity_definition = paste(protein_id_str, peptide_seq_str, sep = " + "),
    run_count_definition = "distinct_run_ids",
    retention_policy = retention_policy,
    replicate_group_column = replicate_group_str,
    minimum_distinct_runs = minimum_distinct_runs,
    input_feature_count = nrow(features),
    retained_feature_count = sum(feature_decisions$feature_retained),
    removed_feature_count = sum(!feature_decisions$feature_retained),
    replicate_group_count = nrow(replicate_groups)
  )

  if (isTRUE(return_filter_result)) {
    return(list(
      data = filtered_data,
      support_table = support_table,
      removal_ledger = removal_ledger,
      summary = summary
    ))
  }
  filtered_data
}

