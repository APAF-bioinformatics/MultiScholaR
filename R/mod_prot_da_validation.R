protDaScalarString <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) && nzchar(trimws(value))
}

normaliseProtDaContrastsTable <- function(contrasts_tbl) {
  if (!is.data.frame(contrasts_tbl) || nrow(contrasts_tbl) == 0L) {
    stop("DA contrasts table must be a non-empty data frame.", call. = FALSE)
  }

  raw_contrasts <- if ("contrasts" %in% names(contrasts_tbl)) {
    contrasts_tbl$contrasts
  } else {
    contrasts_tbl[[1L]]
  }
  raw_contrasts <- as.character(raw_contrasts)

  if (any(is.na(raw_contrasts) | !nzchar(trimws(raw_contrasts)))) {
    stop("DA contrasts table contains an empty contrast.", call. = FALSE)
  }
  if (anyDuplicated(raw_contrasts) > 0L) {
    stop("DA contrasts table contains duplicate contrast definitions.", call. = FALSE)
  }

  full_format <- NULL
  if ("full_format" %in% names(contrasts_tbl)) {
    full_format <- as.character(contrasts_tbl$full_format)
  }
  if (is.null(full_format) || any(is.na(full_format) | !nzchar(trimws(full_format)))) {
    full_format <- vapply(raw_contrasts, function(contrast_string) {
      clean_string <- gsub("^group", "", contrast_string)
      clean_string <- gsub("-group", "-", clean_string)
      friendly_name <- gsub("-", "_vs_", clean_string)
      paste0(friendly_name, "=", contrast_string)
    }, character(1), USE.NAMES = FALSE)
  }

  friendly_names <- NULL
  if ("friendly_names" %in% names(contrasts_tbl)) {
    friendly_names <- as.character(contrasts_tbl$friendly_names)
  }
  if (is.null(friendly_names) || any(is.na(friendly_names) | !nzchar(trimws(friendly_names)))) {
    friendly_names <- sub("=.*$", "", full_format)
  }

  if (anyDuplicated(full_format) > 0L || anyDuplicated(friendly_names) > 0L) {
    stop("DA contrasts table contains duplicate contrast labels.", call. = FALSE)
  }

  contrasts_tbl$contrasts <- raw_contrasts
  contrasts_tbl$full_format <- full_format
  contrasts_tbl$friendly_names <- friendly_names
  contrasts_tbl
}

protDaAnalysisContrastStrings <- function(contrasts_tbl) {
  normaliseProtDaContrastsTable(contrasts_tbl)$full_format
}

protDaContrastExpression <- function(contrast_string) {
  sub("^[^=]+=", "", as.character(contrast_string))
}

protDaContrastTerms <- function(contrast_string) {
  expression <- protDaContrastExpression(contrast_string)
  expression <- gsub("`", "", expression, fixed = TRUE)
  expression <- gsub("[()+*/^\\-]", " ", expression)
  terms <- unlist(strsplit(expression, "\\s+"), use.names = FALSE)
  terms <- terms[nzchar(terms)]
  terms[!grepl("^[0-9.]+$", terms)]
}

validateProtDaFilteredSession <- function(session_data, source_dir = NULL) {
  if (!is.null(source_dir) && !dir.exists(source_dir)) {
    stop("Could not find the source directory to load session data.", call. = FALSE)
  }
  if (!is.list(session_data)) {
    stop("Filtered session data must be an R list.", call. = FALSE)
  }

  object <- session_data$current_s4_object
  if (is.null(object) || !isS4(object) || !methods::is(object, "ProteinQuantitativeData")) {
    stop("Filtered session data does not contain a ProteinQuantitativeData object.", call. = FALSE)
  }
  if (!is.data.frame(session_data$design_matrix) || nrow(session_data$design_matrix) == 0L) {
    stop("Filtered session data does not contain a valid design matrix.", call. = FALSE)
  }
  if (!is.data.frame(session_data$contrasts_tbl) || nrow(session_data$contrasts_tbl) == 0L) {
    stop("Filtered session data does not contain a valid contrasts table.", call. = FALSE)
  }
  normaliseProtDaContrastsTable(session_data$contrasts_tbl)

  final_protein_count <- suppressWarnings(as.integer(session_data$final_protein_count))
  final_sample_count <- suppressWarnings(as.integer(session_data$final_sample_count))
  if (length(final_protein_count) != 1L || is.na(final_protein_count) || final_protein_count < 1L) {
    stop("Filtered session data has an invalid protein count.", call. = FALSE)
  }
  if (length(final_sample_count) != 1L || is.na(final_sample_count) || final_sample_count < 1L) {
    stop("Filtered session data has an invalid sample count.", call. = FALSE)
  }

  invisible(TRUE)
}

validateProtDaModelAndContrasts <- function(theObject, formula_string, contrasts_tbl) {
  if (is.null(theObject) || !isS4(theObject) || !methods::is(theObject, "ProteinQuantitativeData")) {
    stop("DA analysis requires a ProteinQuantitativeData object.", call. = FALSE)
  }
  if (!protDaScalarString(formula_string)) {
    stop("DA formula must be a non-empty string.", call. = FALSE)
  }

  design_matrix <- theObject@design_matrix
  if (!is.data.frame(design_matrix) || nrow(design_matrix) == 0L) {
    stop("DA analysis requires a non-empty design matrix.", call. = FALSE)
  }
  if (!theObject@sample_id %in% names(design_matrix)) {
    stop(sprintf("DA design matrix is missing sample column '%s'.", theObject@sample_id), call. = FALSE)
  }

  sample_columns <- setdiff(names(theObject@protein_quant_table), theObject@protein_id_column)
  missing_samples <- setdiff(as.character(design_matrix[[theObject@sample_id]]), sample_columns)
  if (length(missing_samples) > 0L) {
    stop(
      sprintf("DA design matrix contains samples absent from protein data: %s", paste(missing_samples, collapse = ", ")),
      call. = FALSE
    )
  }

  formula_obj <- tryCatch(
    stats::as.formula(formula_string),
    error = function(e) stop(sprintf("Invalid DA formula: %s", e$message), call. = FALSE)
  )
  design_model <- tryCatch(
    {
      model_frame <- stats::model.frame(formula_obj, design_matrix, na.action = stats::na.pass)
      stats::model.matrix(formula_obj, model_frame)
    },
    error = function(e) stop(sprintf("Invalid DA formula: %s", e$message), call. = FALSE)
  )
  if (nrow(design_model) == 0L || ncol(design_model) == 0L) {
    stop("DA formula produced an empty design matrix.", call. = FALSE)
  }

  contrasts_tbl <- normaliseProtDaContrastsTable(contrasts_tbl)
  contrast_strings <- contrasts_tbl$full_format
  contrast_expressions <- protDaContrastExpression(contrast_strings)
  if (any(!nzchar(trimws(contrast_expressions)))) {
    stop("DA contrasts table contains an empty contrast expression.", call. = FALSE)
  }
  if (any(grepl("(^\\s*[+*/^\\-])|([+*/^\\-]\\s*$)", contrast_expressions))) {
    stop("DA contrasts table contains a malformed contrast expression.", call. = FALSE)
  }

  design_terms <- colnames(design_model)
  for (idx in seq_along(contrast_strings)) {
    terms <- protDaContrastTerms(contrast_strings[[idx]])
    if (length(terms) == 0L) {
      stop("DA contrasts table contains an empty contrast expression.", call. = FALSE)
    }
    missing_terms <- setdiff(terms, design_terms)
    if (length(missing_terms) > 0L) {
      stop(
        sprintf(
          "DA contrast '%s' references terms absent from the model matrix: %s",
          contrasts_tbl$contrasts[[idx]],
          paste(missing_terms, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  list(
    contrasts_tbl = contrasts_tbl,
    contrast_strings = contrast_strings,
    design_model = design_model
  )
}
