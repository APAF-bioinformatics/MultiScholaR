module_ci_assert_file_nonempty <- function(path) {
  testthat::expect_true(file.exists(path), info = sprintf("File does not exist: %s", path))
  testthat::expect_true(file.info(path)$size > 0, info = sprintf("File is empty: %s", path))
  invisible(path)
}

module_ci_assert_sample_identity <- function(before, after, expected_dropped = character()) {
  before <- as.character(before)
  after <- as.character(after)
  expected_dropped <- as.character(expected_dropped)
  unexpected_missing <- setdiff(before, c(after, expected_dropped))
  unexpected_added <- setdiff(after, before)

  testthat::expect_true(
    length(unexpected_missing) == 0L && length(unexpected_added) == 0L,
    info = sprintf(
      "Unexpected sample drift; missing=[%s], added=[%s]",
      paste(unexpected_missing, collapse = ", "),
      paste(unexpected_added, collapse = ", ")
    )
  )
  invisible(TRUE)
}

module_ci_assert_feature_identity <- function(before, after, expected_dropped = character()) {
  before <- as.character(before)
  after <- as.character(after)
  expected_dropped <- as.character(expected_dropped)
  unexpected_missing <- setdiff(before, c(after, expected_dropped))
  unexpected_added <- setdiff(after, before)

  testthat::expect_true(
    length(unexpected_missing) == 0L && length(unexpected_added) == 0L,
    info = sprintf(
      "Unexpected feature drift; missing=[%s], added=[%s]",
      paste(unexpected_missing, collapse = ", "),
      paste(unexpected_added, collapse = ", ")
    )
  )
  invisible(TRUE)
}

module_ci_extract_assay_names <- function(object) {
  if (is.null(object)) {
    return(character())
  }
  if (is.data.frame(object) && "assay" %in% names(object)) {
    return(unique(as.character(object$assay)))
  }
  if (is.data.frame(object) && "Assay" %in% names(object)) {
    return(unique(as.character(object$Assay)))
  }
  if (is.list(object) && !is.null(object$assay_names)) {
    return(as.character(unlist(object$assay_names, use.names = FALSE)))
  }
  if (is.list(object) && !is.null(names(object)) && all(nzchar(names(object)))) {
    return(names(object))
  }
  if (methods::is(object, "LipidomicsAssayData") || methods::is(object, "MetaboliteAssayData")) {
    if ("lipid_data" %in% methods::slotNames(object)) {
      return(names(object@lipid_data))
    }
    if ("metabolite_data" %in% methods::slotNames(object)) {
      return(names(object@metabolite_data))
    }
  }
  character()
}

module_ci_assert_assay_provenance <- function(object, expected_assays) {
  expected_assays <- as.character(expected_assays)
  actual <- module_ci_extract_assay_names(object)
  testthat::expect_setequal(actual, expected_assays)
  invisible(actual)
}

module_ci_assert_design_alignment <- function(
    design,
    sample_columns,
    sample_col = "sample"
) {
  if (!sample_col %in% names(design)) {
    stop(sprintf("sample_col not found in design: %s", sample_col), call. = FALSE)
  }
  design_samples <- as.character(design[[sample_col]])
  sample_columns <- as.character(sample_columns)
  testthat::expect_false(anyDuplicated(design_samples) > 0L, info = "Design has duplicate samples")
  testthat::expect_setequal(design_samples, sample_columns)
  invisible(TRUE)
}

module_ci_assert_numeric_finite <- function(
    data,
    columns = NULL,
    allow_na = TRUE
) {
  if (is.null(columns)) {
    columns <- names(data)[vapply(data, is.numeric, logical(1))]
  }
  for (column in columns) {
    values <- data[[column]]
    testthat::expect_true(is.numeric(values), info = sprintf("Column is not numeric: %s", column))
    finite_or_na <- is.finite(values) | (isTRUE(allow_na) & is.na(values))
    testthat::expect_true(all(finite_or_na), info = sprintf("Column has non-finite values: %s", column))
  }
  invisible(TRUE)
}

module_ci_assert_no_duplicate_keys <- function(data, keys) {
  missing_keys <- setdiff(keys, names(data))
  if (length(missing_keys) > 0L) {
    stop(sprintf("Missing key columns: %s", paste(missing_keys, collapse = ", ")), call. = FALSE)
  }
  key_values <- do.call(paste, c(data[keys], sep = "\r"))
  testthat::expect_false(anyDuplicated(key_values) > 0L, info = sprintf("Duplicate keys: %s", paste(keys, collapse = ", ")))
  invisible(TRUE)
}

module_ci_assert_table_schema <- function(
    path,
    required_columns,
    min_rows = 1L,
    sep = "\t"
) {
  module_ci_assert_file_nonempty(path)
  table <- utils::read.delim(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
  testthat::expect_true(nrow(table) >= min_rows)
  testthat::expect_true(all(required_columns %in% names(table)))
  invisible(table)
}

module_ci_table_checksum <- function(data, columns = NULL) {
  data <- as.data.frame(data, stringsAsFactors = FALSE, check.names = FALSE)
  if (is.null(columns)) {
    columns <- names(data)
  }

  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf("Cannot checksum missing columns: %s", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  temp <- tempfile(fileext = ".tsv")
  on.exit(unlink(temp), add = TRUE)
  utils::write.table(
    data[, columns, drop = FALSE],
    file = temp,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
  unname(tools::md5sum(temp))
}

module_ci_assert_checksum_equal <- function(data, expected_checksum, columns = NULL) {
  actual <- module_ci_table_checksum(data, columns = columns)
  label_columns <- if (is.null(columns)) names(data) else columns
  testthat::expect_identical(
    actual,
    as.character(expected_checksum),
    info = sprintf("Fixture checksum mismatch for columns: %s", paste(label_columns, collapse = ", "))
  )
  invisible(actual)
}

module_ci_assert_checksum_unchanged <- function(before, after, columns = NULL) {
  before_checksum <- module_ci_table_checksum(before, columns = columns)
  after_checksum <- module_ci_table_checksum(after, columns = columns)
  label_columns <- if (is.null(columns)) names(before) else columns
  testthat::expect_identical(
    after_checksum,
    before_checksum,
    info = sprintf("Selected input columns changed unexpectedly: %s", paste(label_columns, collapse = ", "))
  )
  invisible(after_checksum)
}

module_ci_missingness_by_column <- function(data, columns = NULL) {
  if (is.null(columns)) {
    columns <- names(data)
  }
  stats::setNames(
    vapply(columns, function(column) sum(is.na(data[[column]])), integer(1)),
    columns
  )
}

module_ci_assert_expected_missingness <- function(data, expected_na_by_column) {
  expected <- unlist(expected_na_by_column, use.names = TRUE)
  actual <- module_ci_missingness_by_column(data, names(expected))
  testthat::expect_identical(as.integer(actual), as.integer(expected))
  invisible(actual)
}

module_ci_nonfinite_columns <- function(data, columns = NULL) {
  if (is.null(columns)) {
    columns <- names(data)[vapply(data, is.numeric, logical(1))]
  }
  columns[vapply(columns, function(column) {
    any(!is.finite(data[[column]]) & !is.na(data[[column]]))
  }, logical(1))]
}

module_ci_assert_nonfinite_columns <- function(data, expected_columns, columns = NULL) {
  actual <- module_ci_nonfinite_columns(data, columns = columns)
  testthat::expect_setequal(actual, as.character(unlist(expected_columns, use.names = FALSE)))
  invisible(actual)
}
