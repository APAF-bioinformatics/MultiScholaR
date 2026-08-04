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

# ----------------------------------------------------------------------------
# filterMinNumPeptidesPerProteinHelper
# ----------------------------------------------------------------------------
# Identification evidence must be calculated while the original stripped and
# modified sequence identities are still available.  In particular, neither a
# row count nor a sum of per-run peptidoform counts is an experiment-wide
# identity count.
.proteinIdentificationEvidenceColumns <- c(
  "identification_peptide_count",
  "identification_peptidoform_count"
)

.calculateProteinIdentificationEvidence <- function(
    input_table,
    protein_id_column = Protein.Group,
    peptide_sequence_column = Stripped.Sequence,
    modified_peptide_sequence_column = Modified.Sequence) {
  input_columns <- names(input_table)
  protein_id_str <- resolvePeptideQcColumnArgument(
    substitute(protein_id_column),
    protein_id_column,
    input_columns,
    environment()
  )
  peptide_sequence_str <- resolvePeptideQcColumnArgument(
    substitute(peptide_sequence_column),
    peptide_sequence_column,
    input_columns,
    environment()
  )
  modified_sequence_str <- resolvePeptideQcColumnArgument(
    substitute(modified_peptide_sequence_column),
    modified_peptide_sequence_column,
    input_columns,
    environment()
  )

  input_table |>
    dplyr::filter(
      !is.na(.data[[protein_id_str]]),
      trimws(as.character(.data[[protein_id_str]])) != ""
    ) |>
    dplyr::mutate(
      .identification_peptide_identity = trimws(
        as.character(.data[[peptide_sequence_str]])
      ),
      .identification_peptidoform_identity = trimws(
        as.character(.data[[modified_sequence_str]])
      )
    ) |>
    dplyr::group_by(.data[[protein_id_str]]) |>
    dplyr::summarise(
      identification_peptide_count = dplyr::n_distinct(
        .data$.identification_peptide_identity[
          !is.na(.data$.identification_peptide_identity) &
            nzchar(.data$.identification_peptide_identity)
        ]
      ),
      identification_peptidoform_count = dplyr::n_distinct(
        .data$.identification_peptidoform_identity[
          !is.na(.data$.identification_peptidoform_identity) &
            nzchar(.data$.identification_peptidoform_identity)
        ]
      ),
      .groups = "drop"
    )
}

.annotateProteinIdentificationEvidence <- function(
    input_table,
    protein_id_column = Protein.Group,
    peptide_sequence_column = Stripped.Sequence,
    modified_peptide_sequence_column = Modified.Sequence) {
  protein_id_str <- resolvePeptideQcColumnArgument(
    substitute(protein_id_column),
    protein_id_column,
    names(input_table),
    environment()
  )
  peptide_sequence_str <- resolvePeptideQcColumnArgument(
    substitute(peptide_sequence_column),
    peptide_sequence_column,
    names(input_table),
    environment()
  )
  modified_sequence_str <- resolvePeptideQcColumnArgument(
    substitute(modified_peptide_sequence_column),
    modified_peptide_sequence_column,
    names(input_table),
    environment()
  )

  evidence <- .calculateProteinIdentificationEvidence(
    input_table = input_table,
    protein_id_column = protein_id_str,
    peptide_sequence_column = peptide_sequence_str,
    modified_peptide_sequence_column = modified_sequence_str
  )

  input_table |>
    dplyr::select(-dplyr::any_of(.proteinIdentificationEvidenceColumns)) |>
    dplyr::left_join(evidence, by = protein_id_str)
}

#' @export
#' @title Filter Proteins by Minimum Distinct Peptide Evidence
#' @description Keep protein groups only when they meet the configured numbers
#'   of distinct stripped peptide sequences and distinct modified peptidoforms.
#'   Frozen identification counts created at the q-value filtering stage are
#'   preferred so later quantitative filtering cannot change identification
#'   support.
#' @param input_table Peptide quantities table in long format
#' @param num_peptides_per_protein_thresh Minimum number of peptides per protein
#' @param num_peptidoforms_per_protein_thresh Minimum number of peptidoforms per protein
#' @param protein_id_column Protein ID column name as string
#' @param peptide_sequence_column Peptide sequence column name
#' @param modified_peptide_sequence_column Modified peptide sequence column name
#' @param peptidoform_ids_column List-column containing the modified peptide
#'   identities retained by an older precursor-to-peptide rollup. This is kept
#'   only for compatibility with already-created in-memory objects.
#' @param core_utilisation core_utilisation to use for parallel processing
filterMinNumPeptidesPerProteinHelper <- function( input_table
          , num_peptides_per_protein_thresh = 1
          , num_peptidoforms_per_protein_thresh = 2
          , protein_id_column = Protein.Ids
          , peptide_sequence_column = Stripped.Sequence
          , modified_peptide_sequence_column = Modified.Sequence
          , peptidoform_ids_column = "peptidoform_ids"
          , core_utilisation) {

  protein_id_str <- resolvePeptideQcColumnArgument(substitute(protein_id_column), protein_id_column, names(input_table), environment())
  peptide_sequence_str <- resolvePeptideQcColumnArgument(substitute(peptide_sequence_column), peptide_sequence_column, names(input_table), environment())
  modified_sequence_str <- resolvePeptideQcColumnArgument(substitute(modified_peptide_sequence_column), modified_peptide_sequence_column, names(input_table), environment())

  frozen_columns_present <- .proteinIdentificationEvidenceColumns %in% names(input_table)
  if (any(frozen_columns_present) && !all(frozen_columns_present)) {
    stop(
      "filterMinNumPeptidesPerProtein: frozen identification evidence is incomplete.",
      call. = FALSE
    )
  }

  if (all(frozen_columns_present)) {
    frozen_evidence <- input_table |>
      dplyr::select(
        dplyr::all_of(protein_id_str),
        dplyr::all_of(.proteinIdentificationEvidenceColumns)
      ) |>
      dplyr::distinct()

    conflicting_evidence <- frozen_evidence |>
      dplyr::count(.data[[protein_id_str]], name = ".evidence_rows") |>
      dplyr::filter(.data$.evidence_rows != 1L)

    if (nrow(conflicting_evidence) > 0L) {
      stop(
        "filterMinNumPeptidesPerProtein: conflicting frozen identification evidence for one or more protein groups.",
        call. = FALSE
      )
    }

    num_peptides_per_protein <- frozen_evidence |>
      dplyr::transmute(
        !!protein_id_str := .data[[protein_id_str]],
        peptides_for_protein_count = as.integer(
          .data$identification_peptide_count
        ),
        peptidoforms_for_protein_count = as.integer(
          .data$identification_peptidoform_count
        )
      )
  } else if (modified_sequence_str %in% names(input_table)) {
    exact_evidence <- .calculateProteinIdentificationEvidence(
      input_table = input_table,
      protein_id_column = protein_id_str,
      peptide_sequence_column = peptide_sequence_str,
      modified_peptide_sequence_column = modified_sequence_str
    )

    num_peptides_per_protein <- exact_evidence |>
      dplyr::transmute(
        !!protein_id_str := .data[[protein_id_str]],
        peptides_for_protein_count = as.integer(
          .data$identification_peptide_count
        ),
        peptidoforms_for_protein_count = as.integer(
          .data$identification_peptidoform_count
        )
      )
  } else if (peptidoform_ids_column %in% names(input_table)) {
    warning(
      "Using legacy peptidoform_ids provenance. Re-run q-value filtering to create frozen identification evidence.",
      call. = FALSE
    )
    peptide_counts <- input_table |>
      dplyr::distinct(.data[[protein_id_str]], .data[[peptide_sequence_str]]) |>
      dplyr::count(
        .data[[protein_id_str]],
        name = "peptides_for_protein_count"
      )

    peptidoform_counts <- input_table |>
      dplyr::select(
        dplyr::all_of(protein_id_str),
        dplyr::all_of(peptidoform_ids_column)
      ) |>
      tidyr::unnest_longer(
        dplyr::all_of(peptidoform_ids_column),
        values_to = ".peptidoform_id",
        keep_empty = FALSE
      ) |>
      dplyr::filter(
        !is.na(.data$.peptidoform_id),
        trimws(as.character(.data$.peptidoform_id)) != ""
      ) |>
      dplyr::distinct(.data[[protein_id_str]], .data$.peptidoform_id) |>
      dplyr::count(
        .data[[protein_id_str]],
        name = "peptidoforms_for_protein_count"
      )

    num_peptides_per_protein <- peptide_counts |>
      dplyr::left_join(peptidoform_counts, by = protein_id_str)
  } else {
    stop(
      paste0(
        "filterMinNumPeptidesPerProtein: exact experiment-wide peptide evidence is unavailable. ",
        "Run q-value filtering before precursor rollup, or supply Modified.Sequence."
      ),
      call. = FALSE
    )
  }

  num_peptides_per_protein <- num_peptides_per_protein |>
    dplyr::mutate(
      peptidoforms_for_protein_count = dplyr::coalesce(
        .data$peptidoforms_for_protein_count,
        0L
      )
    )

  input_table_for_join <- input_table |>
    dplyr::select(
      -dplyr::any_of(c(
        "peptides_for_protein_count",
        "peptidoforms_for_protein_count"
      ))
    )

  protein_peptide_cln <- NA
  valid_thresholds <- length(num_peptides_per_protein_thresh) == 1L &&
    length(num_peptidoforms_per_protein_thresh) == 1L &&
    is.finite(num_peptides_per_protein_thresh) &&
    is.finite(num_peptidoforms_per_protein_thresh) &&
    num_peptides_per_protein_thresh >= 1 &&
    num_peptidoforms_per_protein_thresh >= 1

  if (valid_thresholds) {

    protein_peptide_cln <- input_table_for_join |>
      inner_join( num_peptides_per_protein
                  , by = protein_id_str) |>
      dplyr::filter(   peptidoforms_for_protein_count >= num_peptidoforms_per_protein_thresh
                      ,
                      peptides_for_protein_count >= num_peptides_per_protein_thresh
                     )
  } else {
    stop(
      paste0(
        "filterMinNumPeptidesPerProtein: num_peptides_per_protein_thresh and ",
        "num_peptidoforms_per_protein_thresh must be provided as finite values >= 1."
      )
    )
  }

  protein_peptide_cln
}

# ----------------------------------------------------------------------------
# filterMinNumPeptidesPerSampleHelper
# ----------------------------------------------------------------------------
#' @export
#' @title Filter Samples by Minimum Number of Peptides
#' @description Retain samples containing at least the configured number of
#'   distinct `(protein group, stripped peptide)` identities. Repeated rows do
#'   not add evidence, and inclusion-list entries never manufacture samples.
#' @param inclusion_list Samples to retain regardless of their observed
#'   distinct peptide-identity count.
filterMinNumPeptidesPerSampleHelper <- function(input_table,
                                                 peptides_per_sample_cutoff = 500,
                                                 sample_id_column = Run,
                                                 core_utilisation,
                                                 inclusion_list = c(),
                                                 protein_id_column = Protein.Ids,
                                                 peptide_sequence_column = Stripped.Sequence,
                                                 design_matrix = NULL,
                                                 design_sample_id_column = sample_id_column,
                                                 return_filter_result = FALSE) {
  input_names <- names(input_table)
  sample_id_str <- resolvePeptideQcColumnArgument(
    substitute(sample_id_column), sample_id_column, input_names, environment()
  )
  protein_id_str <- resolvePeptideQcColumnArgument(
    substitute(protein_id_column), protein_id_column, input_names, environment()
  )
  peptide_sequence_str <- resolvePeptideQcColumnArgument(
    substitute(peptide_sequence_column), peptide_sequence_column,
    input_names, environment()
  )
  missing_columns <- setdiff(
    c(sample_id_str, protein_id_str, peptide_sequence_str),
    input_names
  )
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "filterMinNumPeptidesPerSampleHelper: missing required column(s): %s.",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  threshold <- suppressWarnings(as.numeric(peptides_per_sample_cutoff))
  if (length(threshold) != 1L || is.na(threshold) || !is.finite(threshold) ||
      threshold < 1L || threshold != floor(threshold)) {
    stop(
      "filterMinNumPeptidesPerSampleHelper: `peptides_per_sample_cutoff` must be a positive integer.",
      call. = FALSE
    )
  }
  threshold <- as.integer(threshold)
  inclusion_list <- unique(trimws(as.character(inclusion_list)))
  inclusion_list <- inclusion_list[!is.na(inclusion_list) & nzchar(inclusion_list)]
  if (inherits(core_utilisation, "multidplyr_cluster")) {
    warning(
      paste0(
        "filterMinNumPeptidesPerSampleHelper: distinct peptide identities are ",
        "counted deterministically in-process; the supplied cluster is not used."
      ),
      call. = FALSE
    )
  }

  missing_identity <- is.na(input_table[[sample_id_str]]) |
    !nzchar(trimws(as.character(input_table[[sample_id_str]]))) |
    is.na(input_table[[protein_id_str]]) |
    !nzchar(trimws(as.character(input_table[[protein_id_str]]))) |
    is.na(input_table[[peptide_sequence_str]]) |
    !nzchar(trimws(as.character(input_table[[peptide_sequence_str]])))
  if (any(missing_identity)) {
    stop(
      sprintf(
        "filterMinNumPeptidesPerSampleHelper: sample/protein/peptide identity is missing at row(s) %s.",
        paste(utils::head(which(missing_identity), 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  design_status <- "not_supplied_legacy"
  if (!is.null(design_matrix)) {
    design_sample_str <- resolvePeptideQcColumnArgument(
      substitute(design_sample_id_column), design_sample_id_column,
      names(design_matrix), environment()
    )
    if (!design_sample_str %in% names(design_matrix)) {
      stop(
        sprintf(
          "filterMinNumPeptidesPerSampleHelper: design sample column `%s` is missing.",
          design_sample_str
        ),
        call. = FALSE
      )
    }
    .validatePeptideRunDesign(
      input_table = input_table,
      design_matrix = design_matrix,
      input_sample_column = sample_id_str,
      design_sample_column = design_sample_str,
      caller = "filterMinNumPeptidesPerSampleHelper"
    )
    design_status <- "validated_one_row_per_run"
  }

  sample_support <- input_table |>
    dplyr::transmute(
      .qc_sample = as.character(.data[[sample_id_str]]),
      .qc_protein = as.character(.data[[protein_id_str]]),
      .qc_peptide = as.character(.data[[peptide_sequence_str]])
    ) |>
    dplyr::distinct() |>
    dplyr::count(.data$.qc_sample, name = "distinct_peptide_identity_count") |>
    dplyr::mutate(
      inclusion_override = .data$.qc_sample %in% inclusion_list,
      sample_passes = .data$distinct_peptide_identity_count >= threshold |
        .data$inclusion_override,
      decision_reason = dplyr::case_when(
        .data$distinct_peptide_identity_count >= threshold ~ "meets_distinct_identity_threshold",
        .data$inclusion_override ~ "explicit_inclusion_override",
        TRUE ~ "below_distinct_identity_threshold"
      )
    ) |>
    dplyr::arrange(.data$.qc_sample)
  passing_samples <- sample_support |>
    dplyr::filter(.data$sample_passes) |>
    dplyr::transmute(!!sample_id_str := .data$.qc_sample)
  filtered_table <- input_table |>
    dplyr::inner_join(passing_samples, by = sample_id_str)
  removal_ledger <- sample_support |>
    dplyr::filter(!.data$sample_passes) |>
    dplyr::transmute(
      !!sample_id_str := .data$.qc_sample,
      distinct_peptide_identity_count = .data$distinct_peptide_identity_count,
      required_distinct_peptide_identities = threshold,
      failure_reason = .data$decision_reason
    )
  names(sample_support)[names(sample_support) == ".qc_sample"] <- sample_id_str
  input_samples <- unique(as.character(input_table[[sample_id_str]]))
  summary <- list(
    identity_definition = paste(protein_id_str, peptide_sequence_str, sep = " + "),
    count_definition = "distinct_protein_group_stripped_peptide_pairs_per_run",
    peptides_per_sample_cutoff = threshold,
    design_validation = design_status,
    input_sample_count = length(input_samples),
    retained_sample_count = sum(sample_support$sample_passes),
    removed_sample_count = sum(!sample_support$sample_passes),
    inclusion_list = inclusion_list,
    inclusion_entries_absent_from_data = setdiff(inclusion_list, input_samples)
  )

  if (isTRUE(return_filter_result)) {
    return(list(
      data = filtered_table,
      support_table = sample_support,
      removal_ledger = removal_ledger,
      summary = summary
    ))
  }
  filtered_table
}

# ----------------------------------------------------------------------------
# srlQvalueProteotypicPeptideCleanHelper
# ----------------------------------------------------------------------------
.diaNnIdentificationQvalueCutoff <- 0.01

.validateDiaNnIdentificationContract <- function(
    qvalue_threshold,
    global_qvalue_threshold,
    global_pg_qvalue_threshold,
    choose_only_proteotypic_peptide) {
  threshold_values <- suppressWarnings(as.numeric(c(
    qvalue_threshold,
    global_qvalue_threshold,
    global_pg_qvalue_threshold
  )))
  scalar_inputs <- all(vapply(
    list(
      qvalue_threshold,
      global_qvalue_threshold,
      global_pg_qvalue_threshold
    ),
    length,
    integer(1)
  ) == 1L)
  fixed_thresholds <- scalar_inputs &&
    length(threshold_values) == 3L &&
    all(is.finite(threshold_values)) &&
    all(abs(threshold_values - .diaNnIdentificationQvalueCutoff) < .Machine$double.eps^0.5)

  if (!fixed_thresholds) {
    stop(
      paste0(
        "DIA-NN identification q-value thresholds are fixed at 0.01 for ",
        "Q.Value, Global.Q.Value, and Global.PG.Q.Value."
      ),
      call. = FALSE
    )
  }

  proteotypic_value <- suppressWarnings(as.numeric(choose_only_proteotypic_peptide))
  if (length(proteotypic_value) != 1L ||
      !is.finite(proteotypic_value) ||
      proteotypic_value != 1) {
    stop(
      "DIA-NN identification evidence requires Proteotypic == 1.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' @title Filter Peptides by Q-Value and Proteotypic Status
#' @description Keep spectrum-peptide matches that are within the q-value
#'   thresholds and proteotypic, then freeze experiment-wide peptide and
#'   peptidoform identity counts for each protein group.
#' @param input_table Peptide quantities table in long format
#' @param qvalue_threshold Maximum precursor q-value
#' @param global_qvalue_threshold Maximum global q-value
#' @param global_pg_qvalue_threshold Maximum global protein-group q-value
#' @param choose_only_proteotypic_peptide Whether to retain only proteotypic
#'   identifications
#' @param confidence_provenance_mode Either `diann_main_report`, which requires
#'   all confidence fields, or the explicit `upstream_prefiltered` compatibility
#'   mode, which records any unavailable fields as unverified.
#' @param input_matrix_column_ids Columns retained in the filtered table
#' @param protein_id_column Protein-group identity column
#' @param peptide_sequence_column Stripped peptide identity column
#' @param modified_peptide_sequence_column Modified peptidoform identity column
#' @param q_value_column Precursor q-value column
#' @param global_q_value_column Global precursor q-value column
#' @param global_pg_q_value_column Global protein-group q-value column
#' @param proteotypic_peptide_sequence_column Proteotypic indicator column
#' @param exclude_decoys Exclude confidently classified decoys from the
#'   biological analysis view.
#' @param exclude_contaminants Exclude confidently classified contaminants from
#'   the biological analysis view.
#' @param contaminant_manifest Optional user-supplied contaminant accessions.
#' @param protein_ids_column Accession-provenance column used for classification.
#' @param return_exclusion_result Return classification/ledger metadata together
#'   with the filtered data. The default preserves the historical data-frame API.
#' @export
srlQvalueProteotypicPeptideCleanHelper <- function(input_table
                                             , qvalue_threshold = 0.01
                                             , global_qvalue_threshold = 0.01
                                             , choose_only_proteotypic_peptide = 1
                                             ,   input_matrix_column_ids = c("Run"
                                                                       , "Precursor.Id"
                                                                       , "Protein.Ids"
                                                                       , "Stripped.Sequence"
                                                                       , "Modified.Sequence"
                                                                       , "Precursor.Charge"
                                                                       , "Precursor.Quantity"
                                                                       , "Precursor.Normalised")
                                             , protein_id_column = Protein.Ids
                                             , peptide_sequence_column = Stripped.Sequence
                                             , modified_peptide_sequence_column = Modified.Sequence
                                             , q_value_column = Q.Value
                                             , global_q_value_column = Global.Q.Value
                                             , proteotypic_peptide_sequence_column = Proteotypic
                                             , global_pg_qvalue_threshold = 0.01
                                             , confidence_provenance_mode = "diann_main_report"
                                             , global_pg_q_value_column = Global.PG.Q.Value
                                             , exclude_decoys = TRUE
                                             , exclude_contaminants = TRUE
                                             , contaminant_manifest = NULL
                                             , protein_ids_column = "Protein.Ids"
                                             , return_exclusion_result = FALSE) {

  confidence_provenance_mode <- match.arg(
    confidence_provenance_mode,
    c("diann_main_report", "upstream_prefiltered")
  )
  .validateDiaNnIdentificationContract(
    qvalue_threshold = qvalue_threshold,
    global_qvalue_threshold = global_qvalue_threshold,
    global_pg_qvalue_threshold = global_pg_qvalue_threshold,
    choose_only_proteotypic_peptide = choose_only_proteotypic_peptide
  )

  qvalue_threshold <- .diaNnIdentificationQvalueCutoff
  global_qvalue_threshold <- .diaNnIdentificationQvalueCutoff
  global_pg_qvalue_threshold <- .diaNnIdentificationQvalueCutoff

  protein_id_name <- resolvePeptideQcColumnArgument(
    substitute(protein_id_column),
    protein_id_column,
    names(input_table),
    environment()
  )
  peptide_sequence_name <- resolvePeptideQcColumnArgument(
    substitute(peptide_sequence_column),
    peptide_sequence_column,
    names(input_table),
    environment()
  )
  modified_sequence_name <- resolvePeptideQcColumnArgument(
    substitute(modified_peptide_sequence_column),
    modified_peptide_sequence_column,
    names(input_table),
    environment()
  )

  q_val_name <- resolvePeptideQcColumnArgument(
    substitute(q_value_column),
    q_value_column,
    names(input_table),
    environment()
  )
  global_q_val_name <- resolvePeptideQcColumnArgument(
    substitute(global_q_value_column),
    global_q_value_column,
    names(input_table),
    environment()
  )
  global_pg_q_val_name <- resolvePeptideQcColumnArgument(
    substitute(global_pg_q_value_column),
    global_pg_q_value_column,
    names(input_table),
    environment()
  )
  proteotypic_name <- resolvePeptideQcColumnArgument(
    substitute(proteotypic_peptide_sequence_column),
    proteotypic_peptide_sequence_column,
    names(input_table),
    environment()
  )
  confidence_columns <- c(
    q_value = q_val_name,
    global_q_value = global_q_val_name,
    global_pg_q_value = global_pg_q_val_name,
    proteotypic = proteotypic_name
  )

  # [OK] DIAGNOSTIC + DEFENSIVE: Check output column availability
  missing_cols <- input_matrix_column_ids[!input_matrix_column_ids %in% names(input_table)]
  if (identical(confidence_provenance_mode, "upstream_prefiltered")) {
    missing_cols <- setdiff(missing_cols, unname(confidence_columns))
  }
  
  if (length(missing_cols) > 0) {
    error_msg <- paste0(
      "Q-value filter error: Required output columns not found in data.\n",
      "Missing columns: ", paste(missing_cols, collapse = ", "), "\n",
      "Required columns: ", paste(input_matrix_column_ids, collapse = ", "), "\n",
      "Available columns: ", paste(names(input_table), collapse = ", "), "\n\n",
      "This may be caused by:\n",
      "1. Whitespace in column names from config.ini parsing\n",
      "2. Column names with special characters or encoding issues\n",
      "3. Importing data from a different workflow stage"
    )
    logger::log_error(error_msg)
    stop(error_msg)
  }
  
  # [OK] ALSO CHECK: Filter columns exist
  filter_cols <- unname(confidence_columns)
  missing_filter_cols <- filter_cols[!filter_cols %in% names(input_table)]

  if (identical(confidence_provenance_mode, "diann_main_report") &&
      length(missing_filter_cols) > 0) {
    error_msg <- paste0(
      "Q-value filter error: Required filter columns not found in data.\n",
      "Missing filter columns: ", paste(missing_filter_cols, collapse = ", "), "\n",
      "Available columns: ", paste(names(input_table), collapse = ", ")
    )
    logger::log_error(error_msg)
    stop(error_msg)
  }

  if (identical(confidence_provenance_mode, "upstream_prefiltered") &&
      length(missing_filter_cols) > 0) {
    warning(
      paste0(
        "Explicit upstream_prefiltered mode: the following confidence fields ",
        "are unavailable and were not verified: ",
        paste(missing_filter_cols, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  qvalue_filter <- rep(TRUE, nrow(input_table))
  confidence_thresholds <- c(
    q_value = qvalue_threshold,
    global_q_value = global_qvalue_threshold,
    global_pg_q_value = global_pg_qvalue_threshold
  )
  for (metric_name in names(confidence_thresholds)) {
    metric_column <- confidence_columns[[metric_name]]
    if (metric_column %in% names(input_table)) {
      qvalue_filter <- qvalue_filter &
        input_table[[metric_column]] <= confidence_thresholds[[metric_name]]
    }
  }
  if (proteotypic_name %in% names(input_table)) {
    qvalue_filter <- qvalue_filter & input_table[[proteotypic_name]] == 1
  }

  q_valid_rows <- input_table |>
    dplyr::mutate(.source_row_id = dplyr::row_number()) |>
    dplyr::filter(qvalue_filter) |>
    dplyr::select(dplyr::all_of(unique(c(
      ".source_row_id",
      input_matrix_column_ids[input_matrix_column_ids %in% names(input_table)],
      filter_cols[filter_cols %in% names(input_table)]
    ))))

  exclusion_result <- classifyPeptideBiologicalExclusions(
    input_table = q_valid_rows,
    protein_id_column = protein_id_name,
    protein_ids_column = protein_ids_column,
    contaminant_manifest = contaminant_manifest,
    exclude_decoys = .asPeptideExclusionFlag(exclude_decoys, "exclude_decoys"),
    exclude_contaminants = .asPeptideExclusionFlag(
      exclude_contaminants,
      "exclude_contaminants"
    )
  )

  search_srl_quant_cln <- .annotateProteinIdentificationEvidence(
    input_table = exclusion_result$analysis_data,
    protein_id_column = protein_id_name,
    peptide_sequence_column = peptide_sequence_name,
    modified_peptide_sequence_column = modified_sequence_name
  )

  if (isTRUE(return_exclusion_result)) {
    return(list(
      data = search_srl_quant_cln,
      q_valid_classified_data = exclusion_result$classified_data,
      exclusion_ledger = exclusion_result$exclusion_ledger,
      exclusion_summary = exclusion_result$summary,
      contaminant_manifest = exclusion_result$manifest
    ))
  }

  search_srl_quant_cln

}
