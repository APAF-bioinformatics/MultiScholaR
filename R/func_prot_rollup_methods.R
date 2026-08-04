# ----------------------------------------------------------------------------
# Precursor-to-peptide rollup internals
# ----------------------------------------------------------------------------
.rollupSumObservedOrNA <- function(values) {
  if (all(is.na(values))) {
    return(NA_real_)
  }
  sum(values, na.rm = TRUE)
}

.rollupAssertPopulatedKey <- function(input_table, column_name) {
  values <- input_table[[column_name]]
  missing_key <- is.na(values) |
    (is.character(values) & !nzchar(trimws(values)))
  if (any(missing_key)) {
    example_rows <- paste(utils::head(which(missing_key), 5L), collapse = ", ")
    stop(
      sprintf(
        "rollUpPrecursorToPeptideHelper: `%s` contains missing/blank rollup keys at row(s) %s. Supply complete precursor identities before rollup.",
        column_name,
        example_rows
      ),
      call. = FALSE
    )
  }
}

.rollupAssertLinearQuantity <- function(input_table, column_name) {
  values <- input_table[[column_name]]
  if (!is.numeric(values)) {
    stop(
      sprintf(
        "rollUpPrecursorToPeptideHelper: `%s` must be numeric linear-scale abundance, not %s.",
        column_name,
        paste(class(values), collapse = "/")
      ),
      call. = FALSE
    )
  }

  bad <- is.nan(values) | (!is.na(values) & (!is.finite(values) | values < 0))
  if (any(bad)) {
    example_rows <- utils::head(which(bad), 5L)
    example_values <- paste(as.character(values[example_rows]), collapse = ", ")
    stop(
      sprintf(
        paste0(
          "rollUpPrecursorToPeptideHelper: `%s` contains negative or non-finite ",
          "abundance at row(s) %s (value(s): %s). Rollup requires non-negative, ",
          "finite linear-scale values; replace invalid values with NA only when they are truly missing."
        ),
        column_name,
        paste(example_rows, collapse = ", "),
        example_values
      ),
      call. = FALSE
    )
  }
}

.rollupDuplicateExamples <- function(input_table, identity_columns, maximum = 5L) {
  duplicated_rows <- duplicated(input_table[identity_columns]) |
    duplicated(input_table[identity_columns], fromLast = TRUE)
  duplicate_keys <- unique(input_table[duplicated_rows, identity_columns, drop = FALSE])
  duplicate_keys <- utils::head(duplicate_keys, maximum)

  apply(duplicate_keys, 1L, function(values) {
    paste(sprintf("%s=%s", identity_columns, as.character(values)), collapse = ", ")
  })
}

.rollupSummarisePrepared <- function(prepared, core_utilisation) {
  grouped <- prepared |>
    dplyr::group_by(
      .data$.rollup_sample,
      .data$.rollup_protein,
      .data$.rollup_peptide
    )

  if (inherits(core_utilisation, "multidplyr_cluster")) {
    multidplyr::cluster_assign(
      core_utilisation,
      .rollupSumObservedOrNA = .rollupSumObservedOrNA
    )
    grouped <- grouped |>
      multidplyr::partition(core_utilisation)
  }

  output <- grouped |>
    dplyr::summarise(
      Peptide.RawQuantity = .rollupSumObservedOrNA(.data$.rollup_raw),
      Peptide.Normalised = .rollupSumObservedOrNA(.data$.rollup_normalised),
      peptidoform_count = length(unique(stats::na.omit(.data$.rollup_modified)))
    )

  if (inherits(core_utilisation, "multidplyr_cluster")) {
    output <- dplyr::collect(output)
  }

  output |>
    dplyr::ungroup() |>
    dplyr::arrange(
      .data$.rollup_sample,
      .data$.rollup_protein,
      .data$.rollup_peptide
    )
}

# ----------------------------------------------------------------------------
# rollUpPrecursorToPeptideHelper
# ----------------------------------------------------------------------------
#' @title Rollup Precursors to Peptides
#' @description Validates and sums linear precursor quantities across charge
#'   states and modified/unmodified forms into stripped-peptide quantities. An
#'   all-missing component set remains `NA`; otherwise observed components are
#'   summed. Frozen experiment-wide identification counts are preserved.
#' @param input_table Long precursor-level input table.
#' @param sample_id_column Run/sample identifier column.
#' @param protein_id_column Active protein grouping column, normally
#'   `Protein.Group` for DIA-NN data.
#' @param peptide_sequence_column Stripped peptide sequence column.
#' @param modified_peptide_sequence_column Modified peptide sequence column.
#' @param precursor_quantity_column Raw, linear precursor quantity column.
#' @param precursor_normalised_column Normalised, linear precursor quantity
#'   column.
#' @param core_utilisation Optional `multidplyr_cluster`; other values use the
#'   sequential path.
#' @param precursor_id_column Declared precursor identity column. Legacy inputs
#'   without this column use modified sequence and charge as a derived identity.
#' @param precursor_charge_column Precursor charge column used only for the
#'   legacy derived identity.
#' @param is_logged_data Whether abundance columns are log transformed. Logged
#'   inputs are rejected because summing them is not quantitatively valid.
#' @param return_rollup_result If `TRUE`, return a list containing `data` and an
#'   audit `summary`; otherwise return the historical data-frame result.
#' @export
rollUpPrecursorToPeptideHelper <- function(input_table,
                                           sample_id_column = Run,
                                           protein_id_column = Protein.Ids,
                                           peptide_sequence_column = Stripped.Sequence,
                                           modified_peptide_sequence_column = Modified.Sequence,
                                           precursor_quantity_column = Precursor.Quantity,
                                           precursor_normalised_column = Precursor.Normalised,
                                           core_utilisation,
                                           precursor_id_column = Precursor.Id,
                                           precursor_charge_column = Precursor.Charge,
                                           is_logged_data = FALSE,
                                           return_rollup_result = FALSE) {
  if (!is.data.frame(input_table)) {
    stop("rollUpPrecursorToPeptideHelper: `input_table` must be a data.frame.", call. = FALSE)
  }
  if (isTRUE(is_logged_data)) {
    stop(
      paste0(
        "rollUpPrecursorToPeptideHelper: logged abundance cannot be summed. ",
        "Roll up non-negative linear precursor quantities before log transformation."
      ),
      call. = FALSE
    )
  }

  input_names <- names(input_table)
  sample_id_str <- resolvePeptideQcColumnArgument(
    substitute(sample_id_column), sample_id_column, input_names, environment()
  )
  protein_id_str <- resolvePeptideQcColumnArgument(
    substitute(protein_id_column), protein_id_column, input_names, environment()
  )
  peptide_sequence_str <- resolvePeptideQcColumnArgument(
    substitute(peptide_sequence_column), peptide_sequence_column, input_names, environment()
  )
  modified_sequence_str <- resolvePeptideQcColumnArgument(
    substitute(modified_peptide_sequence_column), modified_peptide_sequence_column, input_names, environment()
  )
  precursor_quantity_str <- resolvePeptideQcColumnArgument(
    substitute(precursor_quantity_column), precursor_quantity_column, input_names, environment()
  )
  precursor_normalised_str <- resolvePeptideQcColumnArgument(
    substitute(precursor_normalised_column), precursor_normalised_column, input_names, environment()
  )
  precursor_id_str <- resolvePeptideQcColumnArgument(
    substitute(precursor_id_column), precursor_id_column, input_names, environment()
  )
  precursor_charge_str <- resolvePeptideQcColumnArgument(
    substitute(precursor_charge_column), precursor_charge_column, input_names, environment()
  )

  required_columns <- c(
    sample_id_str,
    protein_id_str,
    peptide_sequence_str,
    modified_sequence_str,
    precursor_quantity_str,
    precursor_normalised_str
  )
  missing_columns <- setdiff(required_columns, input_names)
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "rollUpPrecursorToPeptideHelper: missing required column(s): %s.",
        paste(sprintf("`%s`", missing_columns), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  for (key_column in c(sample_id_str, protein_id_str, peptide_sequence_str, modified_sequence_str)) {
    .rollupAssertPopulatedKey(input_table, key_column)
  }
  .rollupAssertLinearQuantity(input_table, precursor_quantity_str)
  .rollupAssertLinearQuantity(input_table, precursor_normalised_str)

  if (precursor_id_str %in% input_names) {
    .rollupAssertPopulatedKey(input_table, precursor_id_str)
    identity_columns <- c(sample_id_str, protein_id_str, precursor_id_str)
    identity_mode <- "declared_precursor_id"
  } else {
    identity_columns <- c(sample_id_str, protein_id_str, modified_sequence_str)
    if (precursor_charge_str %in% input_names) {
      .rollupAssertPopulatedKey(input_table, precursor_charge_str)
      identity_columns <- c(identity_columns, precursor_charge_str)
    }
    identity_mode <- "legacy_derived_modified_sequence_and_charge"
    warning(
      paste0(
        "rollUpPrecursorToPeptideHelper: `", precursor_id_str,
        "` is absent; validating uniqueness with modified sequence",
        if (precursor_charge_str %in% input_names) " and precursor charge" else "",
        ". Add a declared precursor identity column for fully auditable rollup."
      ),
      call. = FALSE
    )
  }

  duplicate_rows <- duplicated(input_table[identity_columns]) |
    duplicated(input_table[identity_columns], fromLast = TRUE)
  if (any(duplicate_rows)) {
    duplicate_examples <- .rollupDuplicateExamples(input_table, identity_columns)
    stop(
      paste0(
        "rollUpPrecursorToPeptideHelper: duplicate precursor identities would inflate abundance. ",
        "Deduplicate or resolve these keys before rollup. Example duplicate key(s): ",
        paste(duplicate_examples, collapse = " | ")
      ),
      call. = FALSE
    )
  }

  modified_values <- as.character(input_table[[modified_sequence_str]])
  modified_values[!nzchar(trimws(modified_values))] <- NA_character_
  prepared <- data.frame(
    .rollup_sample = input_table[[sample_id_str]],
    .rollup_protein = input_table[[protein_id_str]],
    .rollup_peptide = input_table[[peptide_sequence_str]],
    .rollup_modified = modified_values,
    .rollup_raw = input_table[[precursor_quantity_str]],
    .rollup_normalised = input_table[[precursor_normalised_str]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  missingness <- prepared |>
    dplyr::group_by(
      .data$.rollup_sample,
      .data$.rollup_protein,
      .data$.rollup_peptide
    ) |>
    dplyr::summarise(
      .raw_observed = sum(!is.na(.data$.rollup_raw)),
      .raw_missing = sum(is.na(.data$.rollup_raw)),
      .normalised_observed = sum(!is.na(.data$.rollup_normalised)),
      .normalised_missing = sum(is.na(.data$.rollup_normalised)),
      .groups = "drop"
    )

  peptide_normalised_tbl <- .rollupSummarisePrepared(prepared, core_utilisation)
  names(peptide_normalised_tbl)[names(peptide_normalised_tbl) == ".rollup_sample"] <- sample_id_str
  names(peptide_normalised_tbl)[names(peptide_normalised_tbl) == ".rollup_protein"] <- protein_id_str
  names(peptide_normalised_tbl)[names(peptide_normalised_tbl) == ".rollup_peptide"] <- peptide_sequence_str

  evidence_columns <- intersect(
    c(
      "identification_peptide_count",
      "identification_peptidoform_count",
      "peptides_for_protein_count",
      "peptidoforms_for_protein_count"
    ),
    input_names
  )

  if (length(evidence_columns) > 0L) {
    frozen_evidence <- input_table |>
      dplyr::select(
        dplyr::all_of(protein_id_str),
        dplyr::all_of(evidence_columns)
      ) |>
      dplyr::distinct()

    conflicting_evidence <- frozen_evidence |>
      dplyr::count(.data[[protein_id_str]], name = ".evidence_rows") |>
      dplyr::filter(.data$.evidence_rows != 1L)

    if (nrow(conflicting_evidence) > 0L) {
      stop(
        "rollUpPrecursorToPeptideHelper: conflicting frozen identification evidence for one or more protein groups.",
        call. = FALSE
      )
    }

    peptide_normalised_tbl <- peptide_normalised_tbl |>
      dplyr::left_join(frozen_evidence, by = protein_id_str)
  }

  # DIA-NN Protein.Group is the quantitative key. Protein.Ids is retained as
  # accession provenance for annotation and reporting, but never substituted
  # for that active grouping key.
  if (protein_id_str != "Protein.Ids" && "Protein.Ids" %in% input_names) {
    collapse_protein_ids <- function(values) {
      tokens <- unlist(strsplit(as.character(values), ";", fixed = TRUE))
      tokens <- trimws(tokens)
      tokens <- tokens[!is.na(tokens) & nzchar(tokens)]
      paste(sort(unique(tokens)), collapse = ";")
    }

    protein_id_annotations <- input_table |>
      dplyr::group_by(.data[[protein_id_str]]) |>
      dplyr::summarise(
        Protein.Ids = collapse_protein_ids(.data$Protein.Ids),
        .groups = "drop"
      )

    peptide_normalised_tbl <- peptide_normalised_tbl |>
      dplyr::left_join(protein_id_annotations, by = protein_id_str)
  }

  peptide_normalised_tbl <- peptide_normalised_tbl |>
    dplyr::arrange(
      .data[[sample_id_str]],
      .data[[protein_id_str]],
      .data[[peptide_sequence_str]]
    )

  rollup_summary <- list(
    aggregation_rule = "linear_sum_modified_unmodified_and_charge_states",
    active_protein_key = protein_id_str,
    precursor_identity_mode = identity_mode,
    input_precursor_rows = nrow(input_table),
    unique_precursor_identities = nrow(unique(input_table[identity_columns])),
    output_stripped_peptide_identities = nrow(peptide_normalised_tbl),
    raw_partial_missing_outputs = sum(
      missingness$.raw_observed > 0L & missingness$.raw_missing > 0L
    ),
    raw_all_missing_outputs = sum(missingness$.raw_observed == 0L),
    normalised_partial_missing_outputs = sum(
      missingness$.normalised_observed > 0L & missingness$.normalised_missing > 0L
    ),
    normalised_all_missing_outputs = sum(missingness$.normalised_observed == 0L),
    frozen_identification_evidence = if (length(evidence_columns) > 0L) "preserved" else "absent",
    protein_ids_provenance = if (protein_id_str != "Protein.Ids" && "Protein.Ids" %in% input_names) {
      "preserved"
    } else if (protein_id_str == "Protein.Ids") {
      "active_key"
    } else {
      "absent"
    }
  )

  if (isTRUE(return_rollup_result)) {
    return(list(data = peptide_normalised_tbl, summary = rollup_summary))
  }
  peptide_normalised_tbl
}

# ----------------------------------------------------------------------------
# rollUpPrecursorToPeptide
# ----------------------------------------------------------------------------
#' @export
setMethod(
  f = "rollUpPrecursorToPeptide",
  signature = "PeptideQuantitativeData",
  definition = function(theObject, core_utilisation = NULL) {
    peptide_data <- theObject@peptide_data
    protein_id_column <- theObject@protein_id_column
    peptide_sequence_column <- theObject@peptide_sequence_column
    raw_quantity_column <- theObject@raw_quantity_column
    norm_quantity_column <- theObject@norm_quantity_column
    sample_id <- theObject@sample_id
    is_logged_data <- theObject@is_logged_data

    core_utilisation <- checkParamsObjectFunctionSimplify(
      theObject,
      "core_utilisation",
      NA
    )
    theObject <- updateParamInObject(theObject, "core_utilisation")

    rollup_result <- rollUpPrecursorToPeptideHelper(
      input_table = peptide_data,
      sample_id_column = sample_id,
      protein_id_column = protein_id_column,
      peptide_sequence_column = peptide_sequence_column,
      modified_peptide_sequence_column = "Modified.Sequence",
      precursor_quantity_column = raw_quantity_column,
      precursor_normalised_column = norm_quantity_column,
      core_utilisation = core_utilisation,
      precursor_id_column = "Precursor.Id",
      precursor_charge_column = "Precursor.Charge",
      is_logged_data = is_logged_data,
      return_rollup_result = TRUE
    )

    theObject@peptide_data <- rollup_result$data
    theObject@raw_quantity_column <- "Peptide.RawQuantity"
    theObject@norm_quantity_column <- "Peptide.Normalised"
    theObject@args$rollUpPrecursorToPeptide$rollup_summary <- rollup_result$summary
    theObject <- cleanDesignMatrixPeptide(theObject)

    theObject
  }
)
