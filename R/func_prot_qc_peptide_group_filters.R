# ----------------------------------------------------------------------------
# check_case_collision_columns
# ----------------------------------------------------------------------------
#' @title Check for Case-Collision Column Names
#' @description Detects columns in a data frame that differ only by capitalisation
#'   (e.g., both "Group" and "group"). Emits a warning for each collision found.
#' @param df A data frame to check
#' @param df_label A label for the data frame (used in warning messages)
#' @return Invisible NULL. Called for its side-effect (warnings).
#' @keywords internal
#' @export
check_case_collision_columns <- function(df, df_label = "data frame") {
  col_names <- colnames(df)
  lower_names <- tolower(col_names)
  
  # Find groups of columns that have the same lowercase name
  dup_groups <- split(col_names, lower_names)
  collisions <- dup_groups[lengths(dup_groups) > 1]
  
  if (length(collisions) > 0) {
    lapply(collisions, function(collision_set) {
      warning(sprintf(
        "Case-collision detected in %s: columns %s differ only by capitalisation. This may cause unexpected behaviour.",
        df_label, paste(sprintf("'%s'", collision_set), collapse = " and ")
      ), call. = FALSE)
    })
  }
  
  invisible(NULL)
}

# ----------------------------------------------------------------------------
# peptideIntensityFilteringHelper
# ----------------------------------------------------------------------------
#' @export
#' @title Filter Peptides by Direct Observed-Run Counts
#' @description Retain peptide identities observed at or above an intensity
#'   threshold in the required number of distinct design runs and biological
#'   groups. Historical percentage parameters remain available as an explicit,
#'   deprecated compatibility path.
peptideIntensityFilteringHelper <- function(input_table = NULL,
                                             design_matrix = NULL,
                                             min_peptide_intensity_threshold = 15,
                                             sample_id_column = "Run",
                                             grouping_variable = "group",
                                             groupwise_percentage_cutoff = 1,
                                             max_groups_percentage_cutoff = 50,
                                             protein_id_column = Protein.Ids,
                                             peptide_sequence_column = Stripped.Sequence,
                                             peptide_quantity_column = Peptide.Normalised,
                                             core_utilisation = NA,
                                             min_reps_per_group = NULL,
                                             min_groups = NULL,
                                             strict_mode = FALSE,
                                             return_filter_result = FALSE,
                                             ...) {
  extra_args <- list(...)
  if (is.null(input_table)) {
    input_table <- extra_args$input_data
  }
  if (is.null(input_table)) {
    stop(
      "peptideIntensityFilteringHelper: input_table or input_data must be provided.",
      call. = FALSE
    )
  }
  if (!is.data.frame(input_table) || !is.data.frame(design_matrix)) {
    stop(
      "peptideIntensityFilteringHelper: `input_table` and `design_matrix` must be data.frames.",
      call. = FALSE
    )
  }
  if (!is.null(extra_args$raw_quantity_column)) {
    peptide_quantity_column <- extra_args$raw_quantity_column
  }
  if (inherits(core_utilisation, "multidplyr_cluster")) {
    warning(
      paste0(
        "peptideIntensityFilteringHelper: the direct distinct-run count path ",
        "executes deterministically in-process; the supplied multidplyr cluster is not used."
      ),
      call. = FALSE
    )
  }

  input_cols <- names(input_table)
  design_cols <- names(design_matrix)
  sample_id_str <- resolvePeptideQcColumnArgument(
    substitute(sample_id_column), sample_id_column,
    candidates = union(input_cols, design_cols), environment()
  )
  group_var_str <- resolvePeptideQcColumnArgument(
    substitute(grouping_variable), grouping_variable,
    candidates = design_cols, environment()
  )
  peptide_quant_str <- resolvePeptideQcColumnArgument(
    substitute(peptide_quantity_column), peptide_quantity_column,
    candidates = input_cols, environment()
  )
  protein_id_str <- resolvePeptideQcColumnArgument(
    substitute(protein_id_column), protein_id_column,
    candidates = input_cols, environment()
  )
  peptide_seq_str <- resolvePeptideQcColumnArgument(
    substitute(peptide_sequence_column), peptide_sequence_column,
    candidates = input_cols, environment()
  )

  check_case_collision_columns(design_matrix, "design_matrix")
  resolve_design_column <- function(requested) {
    if (requested %in% design_cols) {
      return(requested)
    }
    matches <- which(tolower(design_cols) == tolower(requested))
    if (length(matches) != 1L) {
      stop(
        sprintf(
          "peptideIntensityFilteringHelper: column `%s` not found unambiguously in design_matrix. Available: %s.",
          requested,
          paste(design_cols, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    design_cols[[matches]]
  }
  sample_id_str <- resolve_design_column(sample_id_str)
  group_var_str <- resolve_design_column(group_var_str)

  required_input_keys <- c(protein_id_str, peptide_seq_str)
  missing_input_keys <- setdiff(required_input_keys, input_cols)
  if (length(missing_input_keys) > 0L) {
    stop(
      sprintf(
        "peptideIntensityFilteringHelper: missing input identity column(s): %s.",
        paste(sprintf("`%s`", missing_input_keys), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  design_minimal <- design_matrix |>
    dplyr::transmute(
      .filter_run = as.character(.data[[sample_id_str]]),
      .filter_group = as.character(.data[[group_var_str]])
    )
  missing_design_key <- is.na(design_minimal$.filter_run) |
    !nzchar(trimws(design_minimal$.filter_run)) |
    is.na(design_minimal$.filter_group) |
    !nzchar(trimws(design_minimal$.filter_group))
  if (any(missing_design_key)) {
    stop(
      sprintf(
        "peptideIntensityFilteringHelper: design run/group identity is missing at row(s) %s.",
        paste(utils::head(which(missing_design_key), 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  duplicated_design_runs <- unique(design_minimal$.filter_run[
    duplicated(design_minimal$.filter_run) |
      duplicated(design_minimal$.filter_run, fromLast = TRUE)
  ])
  if (length(duplicated_design_runs) > 0L) {
    stop(
      sprintf(
        "peptideIntensityFilteringHelper: design_matrix must contain one row per run; duplicate run ID(s): %s.",
        paste(utils::head(duplicated_design_runs, 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  feature_missing <- is.na(input_table[[protein_id_str]]) |
    !nzchar(trimws(as.character(input_table[[protein_id_str]]))) |
    is.na(input_table[[peptide_seq_str]]) |
    !nzchar(trimws(as.character(input_table[[peptide_seq_str]])))
  if (any(feature_missing)) {
    stop(
      sprintf(
        "peptideIntensityFilteringHelper: protein/peptide identity is missing at row(s) %s.",
        paste(utils::head(which(feature_missing), 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  is_long <- sample_id_str %in% input_cols
  if (is_long) {
    if (!peptide_quant_str %in% input_cols) {
      stop(
        sprintf(
          "peptideIntensityFilteringHelper: quantity column `%s` is missing from long input.",
          peptide_quant_str
        ),
        call. = FALSE
      )
    }
    if (!is.numeric(input_table[[peptide_quant_str]])) {
      stop(
        sprintf(
          "peptideIntensityFilteringHelper: quantity column `%s` must be numeric.",
          peptide_quant_str
        ),
        call. = FALSE
      )
    }
    abundance_long <- input_table |>
      dplyr::transmute(
        .filter_protein = as.character(.data[[protein_id_str]]),
        .filter_peptide = as.character(.data[[peptide_seq_str]]),
        .filter_run = as.character(.data[[sample_id_str]]),
        .filter_quantity = .data[[peptide_quant_str]]
      )

    duplicate_quantitative_identity <- abundance_long |>
      dplyr::count(
        .data$.filter_protein,
        .data$.filter_peptide,
        .data$.filter_run,
        name = ".identity_rows"
      ) |>
      dplyr::filter(.data$.identity_rows > 1L)
    if (nrow(duplicate_quantitative_identity) > 0L) {
      examples <- duplicate_quantitative_identity |>
        dplyr::slice_head(n = 5L) |>
        dplyr::transmute(
          .example = sprintf(
            "%s/%s/%s",
            .data$.filter_protein,
            .data$.filter_peptide,
            .data$.filter_run
          )
        ) |>
        dplyr::pull(.data$.example)
      stop(
        paste0(
          "peptideIntensityFilteringHelper: duplicate (protein, peptide, run) ",
          "quantitative identities detected: ", paste(examples, collapse = ", "),
          ". Roll up or deduplicate precursor rows before intensity filtering."
        ),
        call. = FALSE
      )
    }
  } else {
    design_runs <- design_minimal$.filter_run
    missing_run_columns <- setdiff(design_runs, input_cols)
    if (length(missing_run_columns) > 0L) {
      stop(
        sprintf(
          "peptideIntensityFilteringHelper: wide input is missing design run column(s): %s.",
          paste(utils::head(missing_run_columns, 10L), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    non_numeric_runs <- design_runs[!vapply(
      input_table[design_runs], is.numeric, logical(1)
    )]
    if (length(non_numeric_runs) > 0L) {
      stop(
        sprintf(
          "peptideIntensityFilteringHelper: wide run column(s) must be numeric: %s.",
          paste(non_numeric_runs, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    abundance_long <- input_table |>
      dplyr::transmute(
        .filter_protein = as.character(.data[[protein_id_str]]),
        .filter_peptide = as.character(.data[[peptide_seq_str]]),
        dplyr::across(dplyr::all_of(design_runs))
      ) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(design_runs),
        names_to = ".filter_run",
        values_to = ".filter_quantity"
      )
  }

  input_runs <- unique(abundance_long$.filter_run)
  design_runs <- design_minimal$.filter_run
  missing_input_run <- is.na(abundance_long$.filter_run) |
    !nzchar(trimws(abundance_long$.filter_run))
  if (any(missing_input_run)) {
    stop(
      sprintf(
        "peptideIntensityFilteringHelper: input run identity is missing at row(s) %s.",
        paste(utils::head(which(missing_input_run), 5L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  unmatched_input_runs <- setdiff(input_runs, design_runs)
  absent_design_runs <- setdiff(design_runs, input_runs)
  if (length(unmatched_input_runs) > 0L || length(absent_design_runs) > 0L) {
    stop(
      paste0(
        "peptideIntensityFilteringHelper: input/design run identities do not match. ",
        if (length(unmatched_input_runs) > 0L) {
          paste0("Unmapped input run(s): ", paste(utils::head(unmatched_input_runs, 10L), collapse = ", "), ". ")
        } else {
          ""
        },
        if (length(absent_design_runs) > 0L) {
          paste0("Design run(s) absent from input: ", paste(utils::head(absent_design_runs, 10L), collapse = ", "), ".")
        } else {
          ""
        }
      ),
      call. = FALSE
    )
  }

  threshold <- suppressWarnings(as.numeric(min_peptide_intensity_threshold))
  if (length(threshold) != 1L || is.na(threshold) || !is.finite(threshold)) {
    stop(
      "peptideIntensityFilteringHelper: `min_peptide_intensity_threshold` must be one finite number.",
      call. = FALSE
    )
  }
  strict_mode <- isTRUE(strict_mode)
  direct_count_mode <- !is.null(min_reps_per_group) || !is.null(min_groups)
  if (direct_count_mode && (is.null(min_reps_per_group) || is.null(min_groups))) {
    stop(
      "peptideIntensityFilteringHelper: provide both `min_reps_per_group` and `min_groups` for direct-count filtering.",
      call. = FALSE
    )
  }

  group_sizes <- design_minimal |>
    dplyr::count(.data$.filter_group, name = "design_run_count") |>
    dplyr::arrange(.data$.filter_group)
  total_groups <- nrow(group_sizes)
  if (direct_count_mode) {
    validate_count <- function(value, name, maximum = Inf) {
      numeric_value <- suppressWarnings(as.numeric(value))
      if (length(numeric_value) != 1L || is.na(numeric_value) ||
          !is.finite(numeric_value) || numeric_value < 1 ||
          numeric_value != floor(numeric_value) || numeric_value > maximum) {
        stop(
          sprintf(
            "peptideIntensityFilteringHelper: `%s` must be an integer from 1 to %s.",
            name,
            if (is.finite(maximum)) maximum else "Inf"
          ),
          call. = FALSE
        )
      }
      as.integer(numeric_value)
    }
    min_reps_per_group <- validate_count(min_reps_per_group, "min_reps_per_group")
    min_groups <- validate_count(
      min_groups,
      "min_groups",
      if (strict_mode) Inf else total_groups
    )
    group_sizes$required_observed_runs <- min_reps_per_group
    rule_mode <- if (strict_mode) "strict_all_design_runs" else "direct_counts"
  } else {
    percentages <- c(groupwise_percentage_cutoff, max_groups_percentage_cutoff)
    if (length(percentages) != 2L || any(!is.finite(percentages)) ||
        any(percentages < 0 | percentages > 100)) {
      stop(
        "peptideIntensityFilteringHelper: legacy percentage cutoffs must be finite values from 0 to 100.",
        call. = FALSE
      )
    }
    warning(
      paste0(
        "peptideIntensityFilteringHelper: percentage-only missingness arguments are deprecated. ",
        "Use `min_reps_per_group` and `min_groups`; this call uses an exact deterministic adapter."
      ),
      call. = FALSE
    )
    group_sizes$required_observed_runs <- ceiling(
      group_sizes$design_run_count * (1 - groupwise_percentage_cutoff / 100) - 1e-12
    )
    allowed_failed_groups <- floor(
      total_groups * max_groups_percentage_cutoff / 100 + 1e-12
    )
    min_groups <- total_groups - allowed_failed_groups
    min_reps_per_group <- NA_integer_
    rule_mode <- if (strict_mode) "strict_all_design_runs" else "legacy_percentage_adapter"
  }

  if (strict_mode) {
    group_sizes$required_observed_runs <- group_sizes$design_run_count
    min_groups <- total_groups
  }

  abundance_long <- abundance_long |>
    dplyr::left_join(design_minimal, by = ".filter_run") |>
    dplyr::mutate(
      .is_observed = is.finite(.data$.filter_quantity) &
        .data$.filter_quantity >= threshold
    )

  features <- abundance_long |>
    dplyr::distinct(.data$.filter_protein, .data$.filter_peptide)
  observed_support <- abundance_long |>
    dplyr::filter(.data$.is_observed) |>
    dplyr::group_by(
      .data$.filter_protein,
      .data$.filter_peptide,
      .data$.filter_group
    ) |>
    dplyr::summarise(
      observed_run_count = dplyr::n_distinct(.data$.filter_run),
      .groups = "drop"
    )

  support_table <- tidyr::expand_grid(
    features,
    .filter_group = group_sizes$.filter_group
  ) |>
    dplyr::left_join(
      observed_support,
      by = c(".filter_protein", ".filter_peptide", ".filter_group")
    ) |>
    dplyr::mutate(observed_run_count = dplyr::coalesce(.data$observed_run_count, 0L)) |>
    dplyr::left_join(group_sizes, by = ".filter_group") |>
    dplyr::mutate(
      group_passes = .data$observed_run_count >= .data$required_observed_runs,
      group_failure_reason = dplyr::if_else(
        .data$group_passes,
        NA_character_,
        if (strict_mode) "strict_missing_design_run" else "insufficient_observed_runs"
      )
    ) |>
    dplyr::arrange(
      .data$.filter_protein,
      .data$.filter_peptide,
      .data$.filter_group
    )

  feature_decisions <- support_table |>
    dplyr::group_by(.data$.filter_protein, .data$.filter_peptide) |>
    dplyr::summarise(
      passing_group_count = sum(.data$group_passes),
      total_group_count = dplyr::n(),
      feature_passes = .data$passing_group_count >= min_groups,
      failure_reason = dplyr::if_else(
        .data$feature_passes,
        NA_character_,
        if (strict_mode) "strict_missing_design_run" else "insufficient_passing_groups"
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$.filter_protein, .data$.filter_peptide)

  removal_ledger <- feature_decisions |>
    dplyr::filter(!.data$feature_passes) |>
    dplyr::transmute(
      !!protein_id_str := .data$.filter_protein,
      !!peptide_seq_str := .data$.filter_peptide,
      failure_reason = .data$failure_reason,
      passing_group_count = .data$passing_group_count,
      required_passing_groups = min_groups
    )

  filtered_data <- input_table
  if (nrow(removal_ledger) > 0L) {
    filtered_data <- input_table |>
      dplyr::anti_join(
        removal_ledger |>
          dplyr::select(dplyr::all_of(c(protein_id_str, peptide_seq_str))),
        by = c(protein_id_str, peptide_seq_str)
      )
  }

  names(support_table)[names(support_table) == ".filter_protein"] <- protein_id_str
  names(support_table)[names(support_table) == ".filter_peptide"] <- peptide_seq_str
  names(support_table)[names(support_table) == ".filter_group"] <- group_var_str

  filter_summary <- list(
    rule_mode = rule_mode,
    strict_mode = strict_mode,
    min_reps_per_group = min_reps_per_group,
    min_groups = as.integer(min_groups),
    intensity_threshold = threshold,
    intensity_quantity_column = peptide_quant_str,
    input_feature_count = nrow(features),
    retained_feature_count = sum(feature_decisions$feature_passes),
    removed_feature_count = sum(!feature_decisions$feature_passes),
    design_run_count = nrow(design_minimal),
    design_group_count = total_groups,
    non_finite_measurement_count = sum(!is.finite(abundance_long$.filter_quantity))
  )

  if (isTRUE(return_filter_result)) {
    return(list(
      data = filtered_data,
      support_table = support_table,
      removal_ledger = removal_ledger,
      summary = filter_summary
    ))
  }
  filtered_data
}

# ----------------------------------------------------------------------------
# removePeptidesWithMissingValuesPercentHelper
# ----------------------------------------------------------------------------
#' @title Remove Peptides with Missing Values
#'@param input_table An input table with a column containing the row ID and the rest of the columns representing abundance values for each sample.
#'@param cols A tidyselect command to select the columns. This includes the functions starts_with(), ends_with(), contains(), matches(), and num_range()
#'@param design_matrix A data frame with a column containing the sample ID (as per the sample_id param) and the experimental group (as per the group param). Each row as the sample ID as row name in the data frame.
#'@param sample_id The name of the column in design_matrix table that has the sample ID.
#'@param protein_id_column Protein ID column name
#'@param peptide_sequence_column Peptide sequence column name
#'@param grouping_variable The name of the column in design_matrix table that has the experimental group.
#'@param groupwise_percentage_cutoff The maximum percentage of values below threshold allow in each group for a peptide .
#'@param max_groups_percentage_cutoff The maximum percentage of groups allowed with too many samples with peptide abundance values below threshold.
#'@param abundance_threshold Abundance threshold in which the protein in the sample must be above for it to be considered for inclusion into data analysis.
#'@param abundance_column The name of the column containing the abundance values.
#'@return A filtered data.frame
#'@export
removePeptidesWithMissingValuesPercentHelper <- function(input_table
                                               , design_matrix
                                               , sample_id
                                               , protein_id_column
                                               , peptide_sequence_column
                                               , grouping_variable
                                               , groupwise_percentage_cutoff = 1
                                               , max_groups_percentage_cutoff = 50
                                               , abundance_threshold = 0
                                               , abundance_column = "Abundance") {

  message("+===========================================================================+")
  message("|  DEBUG66: Entering removePeptidesWithMissingValuesPercentHelper (OPTIMIZED)|")
  message("+===========================================================================+")

  input_cols <- colnames(input_table)
  design_cols <- colnames(design_matrix)
  sample_id_str <- resolvePeptideQcColumnArgument(
    substitute(sample_id),
    sample_id,
    candidates = union(input_cols, design_cols),
    environment()
  )
  group_var_str <- resolvePeptideQcColumnArgument(
    substitute(grouping_variable),
    grouping_variable,
    candidates = design_cols,
    environment()
  )
  protein_id_str <- resolvePeptideQcColumnArgument(
    substitute(protein_id_column),
    protein_id_column,
    candidates = input_cols,
    environment()
  )
  peptide_seq_str <- resolvePeptideQcColumnArgument(
    substitute(peptide_sequence_column),
    peptide_sequence_column,
    candidates = input_cols,
    environment()
  )
  abundance_col_str <- resolvePeptideQcColumnArgument(
    substitute(abundance_column),
    abundance_column,
    candidates = input_cols,
    environment()
  )

  message(sprintf("   DEBUG66: Args: threshold=%g, group_cutoff=%g%%, max_groups_fail=%g%%", 
                  abundance_threshold, groupwise_percentage_cutoff, max_groups_percentage_cutoff))
  check_case_collision_columns(design_matrix, "design_matrix")

  # Case-insensitive column resolution for design matrix
  dm_cols <- colnames(design_matrix)
  if (!group_var_str %in% dm_cols) {
    match_idx <- which(tolower(dm_cols) == tolower(group_var_str))
    if (length(match_idx) == 1) {
      message(sprintf("   DEBUG66: Column '%s' not found in design matrix, using '%s' (case-insensitive match)", 
                      group_var_str, dm_cols[match_idx]))
      group_var_str <- dm_cols[match_idx]
    } else {
      stop(sprintf("Column '%s' not found in design matrix. Available: %s", 
                    group_var_str, paste(dm_cols, collapse = ", ")))
    }
  }
  if (!sample_id_str %in% dm_cols) {
    match_idx <- which(tolower(dm_cols) == tolower(sample_id_str))
    if (length(match_idx) == 1) {
      message(sprintf("   DEBUG66: Column '%s' not found in design matrix, using '%s' (case-insensitive match)", 
                      sample_id_str, dm_cols[match_idx]))
      sample_id_str <- dm_cols[match_idx]
    }
  }

  # Prepare minimal design matrix
  design_matrix_minimal <- design_matrix |>
    dplyr::select(!!rlang::sym(sample_id_str), !!rlang::sym(group_var_str)) |>
    dplyr::mutate(!!rlang::sym(sample_id_str) := as.character(!!rlang::sym(sample_id_str)))

  # Handle long vs wide input
  is_long <- sample_id_str %in% colnames(input_table)
  
  if (is_long) {
    message("   DEBUG66: Input table appears to be in LONG format. Skipping pivot.")
    abundance_long <- input_table |>
      dplyr::mutate(
        !!rlang::sym(sample_id_str) := as.character(!!rlang::sym(sample_id_str)),
        !!rlang::sym(abundance_col_str) := dplyr::if_else(
          is.finite(!!rlang::sym(abundance_col_str)),
          !!rlang::sym(abundance_col_str),
          NA_real_
        )
      ) |>
      dplyr::left_join(design_matrix_minimal, by = sample_id_str)
  } else {
    message("   DEBUG66: Input table appears to be in WIDE format. Pivoting long.")
    abundance_long <- input_table |>
      tidyr::pivot_longer(
        cols = !c(!!rlang::sym(protein_id_str), !!rlang::sym(peptide_seq_str)),
        names_to = sample_id_str,
        values_to = abundance_col_str
      ) |>
      dplyr::mutate(
        !!rlang::sym(sample_id_str) := as.character(!!rlang::sym(sample_id_str)),
        !!rlang::sym(abundance_col_str) := dplyr::if_else(
          is.finite(!!rlang::sym(abundance_col_str)),
          !!rlang::sym(abundance_col_str),
          NA_real_
        )
      ) |>
      dplyr::left_join(design_matrix_minimal, by = sample_id_str)
  }

  gc()

  # Count samples per group
  count_values_per_group <- abundance_long |>
    dplyr::distinct(!!rlang::sym(sample_id_str), !!rlang::sym(group_var_str)) |>
    dplyr::count(!!rlang::sym(group_var_str), name = "num_per_group")

  # Calculate missing stats per group
  count_percent_missing_per_group <- abundance_long |>
    dplyr::mutate(is_below = dplyr::if_else(
      !is.na(!!rlang::sym(abundance_col_str)) & 
      !!rlang::sym(abundance_col_str) >= abundance_threshold, 
      0, 1
    )) |>
    dplyr::group_by(!!rlang::sym(protein_id_str), !!rlang::sym(peptide_seq_str), !!rlang::sym(group_var_str)) |>
    dplyr::summarise(num_observed_above = sum(1 - is_below, na.rm = TRUE), .groups = "drop") |>
    # CRITICAL: Ensure all peptide-group combinations exist. 
    # If a peptide has 0 rows for a group, complete() will add it with num_observed_above = 0.
    tidyr::complete(
      tidyr::nesting(!!rlang::sym(protein_id_str), !!rlang::sym(peptide_seq_str)),
      !!rlang::sym(group_var_str),
      fill = list(num_observed_above = 0)
    ) |>
    dplyr::left_join(count_values_per_group, by = group_var_str) |>
    dplyr::mutate(num_below_per_group = num_per_group - num_observed_above) |>
    dplyr::mutate(perc_below_per_group = num_below_per_group / num_per_group * 100)

  # Identify peptides to remove
  total_num_of_groups <- nrow(count_values_per_group)
  
  peptides_to_remove <- count_percent_missing_per_group |>
    dplyr::group_by(!!rlang::sym(protein_id_str), !!rlang::sym(peptide_seq_str)) |>
    dplyr::summarise(
      num_groups_failed = sum(perc_below_per_group > groupwise_percentage_cutoff, na.rm = TRUE),
      percent_groups_failed = (num_groups_failed / total_num_of_groups) * 100,
      .groups = "drop"
    ) |>
    dplyr::filter(percent_groups_failed > max_groups_percentage_cutoff)
    
  num_removed <- nrow(peptides_to_remove)
  message(sprintf("   Results: Identified %d peptides to remove.", num_removed))

  # Filter original table
  if (num_removed > 0) {
    if (is_long) {
      filtered_tbl <- input_table |>
        dplyr::anti_join(peptides_to_remove, by = c(protein_id_str, peptide_seq_str))
    } else {
      # For wide format, convert peptides_to_remove back to IDs to filter
      filtered_tbl <- input_table |>
        dplyr::anti_join(peptides_to_remove, by = c(protein_id_str, peptide_seq_str))
    }
  } else {
    filtered_tbl <- input_table
  }

  message(sprintf("   Results: Peptides: %d -> %d", nrow(input_table), nrow(filtered_tbl)))
  
  message("+===========================================================================+")
  message("|  Exiting removePeptidesWithMissingValuesPercentHelper                     |")
  message("+===========================================================================+")

  gc()
  return(filtered_tbl)

}
