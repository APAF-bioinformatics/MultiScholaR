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
#' @param input_table,peptides_per_sample_cutoff,sample_id_column,core_utilisation,protein_id_column,peptide_sequence_column,design_matrix,design_sample_id_column,return_filter_result Runtime inputs used by this function; see the usage section for accepted values.
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

