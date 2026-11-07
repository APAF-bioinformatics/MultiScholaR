# KIT-102 W2: the feature-table import core of lipid, metab and trx (scoping I1 to I13). Functions
# over a request and the omic's GUI configuration (`omicGuiConfig()`): the header read, the preview,
# sample and column resolution, the assay-name rules, assembly with the second-assay reader,
# sanitising, the numeric-like conversion, validation, the column mapping, the artifacts and the
# processing log. `importOmicsTables()` runs them in order and returns the payload, or signals a
# classed condition with today's message. The GUI keeps the notifications and the commit to the bus.
# Moved from `R/mod_lipid_import_{selection,processing}_helpers.R`,
# `R/mod_metab_import_{processing,column}_helpers.R` and `R/mod_transcriptomics.R`; the trx count
# policy and preview are in `R/func_trx_counts.R`. Recorded divergence from conventions/R.md: the
# per-omic validation messages are held here by omic type (`omicsImportRules()`), as the
# configuration does not carry them yet.

#' The columns the lipid import's "exclude normalised" option drops from the sample columns
#' (`resolveLipidImportSampleColumns`, mod_lipid_import_processing_helpers.R:107).
#' @noRd
OMICS_IMPORT_NORMALISED_COLUMNS <- "_norm(ali[sz]ed|allized)?$|normalized"

# ------------------------------------------------------------------------------------------------
# Small helpers and the per-omic rules.
# ------------------------------------------------------------------------------------------------

#' Whether a value is one non-empty, non-missing string.
#' @noRd
omicsImportScalarString <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

#' A value's first element as a string; `""` for NULL or empty.
#' @noRd
omicsImportScalarValue <- function(value) {
    if (is.null(value) || length(value) == 0L) return("")
    as.character(value[[1]])
}

#' Refuse an import with a message, classed `omics_import_invalid` (and any extra classes first).
#' @noRd
omicsImportAbort <- function(message, class = character()) {
    rlang::abort(message, class = c(class, "omics_import_invalid"))
}

#' The per-omic validation rules and messages (each omic's today). Lipid and trx run the lipid
#' validation (`validateLipidImportProcessingState`, mod_lipid_import_processing_helpers.R:479);
#' trx refuses it with its count-policy class (`trxValidateImport`, mod_transcriptomics.R:115) and
#' takes every column but the id and name columns as samples (mod_transcriptomics.R:134). Metab
#' runs `validateMetabImportWorkflowPayload` (mod_metab_import_processing_helpers.R:114).
#' @param omic_type `lipidomics`, `metabolomics` or `transcriptomics`.
#' @return A list: `messages`, `assay_prefix` (mapping errors and warnings prefixed by the assay),
#'   `unique` (errors and warnings de-duplicated), `numeric_like` (the numeric check is
#'   `omicsImportValuesAreNumericLike`, else `is.numeric`), `samples` (`request` or
#'   `all_but_id`), `refusal_class`.
#' @noRd
omicsImportRules <- function(omic_type) {
    lipid <- list(
        messages = c(no_assays = "No lipidomics assays were assembled",
            unnamed = "All lipidomics assays must have non-empty names",
            duplicate = "Duplicate lipidomics assay names are not allowed",
            no_samples = "No sample columns specified",
            not_data_frame = "Assay '%s' is not a data frame",
            non_numeric = "Assay '%s' has non-numeric sample columns: %s"),
        assay_prefix = TRUE, unique = FALSE, numeric_like = FALSE, samples = "request",
        refusal_class = character())
    switch(omic_type,
        lipidomics = lipid,
        transcriptomics = utils::modifyList(lipid, list(samples = "all_but_id",
            refusal_class = c("trx_counts_invalid", "trx_error"))),
        metabolomics = list(
            messages = c(no_assays = "No assay data provided.",
                unnamed = "All assays must have non-empty names.",
                duplicate = "Duplicate assay names are not allowed: %s",
                no_id = "Metabolite ID column is not specified.",
                no_samples = "No sample columns specified.",
                not_data_frame = "Assay '%s' is not a data frame.",
                non_numeric = "Assay '%s' has non-numeric sample column(s): %s"),
            assay_prefix = FALSE, unique = TRUE, numeric_like = TRUE, samples = "request",
            refusal_class = character()),
        cli::cli_abort("no import rules for omic {.val {omic_type}}.",
            class = "omics_import_invalid"))
}

#' The reader for a format: the configured reader of that name, else the default.
#' @param readers A named list of reader functions with a `default`.
#' @noRd
omicsImportReader <- function(readers, format) {
    key <- tolower(omicsImportScalarValue(format))
    if (nzchar(key) && key %in% names(readers)) readers[[key]] else readers$default
}

# ------------------------------------------------------------------------------------------------
# I1 to I4: headers, preview, sample and column resolution.
# ------------------------------------------------------------------------------------------------

#' I1. A delimited file's header names: the first line split on commas (when it has one) or tabs,
#' with surrounding quotes removed. `character(0)` when the file cannot be read.
#' @param path The file.
#' @noRd
readOmicsImportHeaders <- function(path) {
    tryCatch({
        con <- file(path, "r")
        on.exit(close(con), add = TRUE)
        first_line <- readLines(con, n = 1)
        raw_headers <- if (grepl(",", first_line)) {
            strsplit(first_line, ",")[[1]]
        } else {
            strsplit(first_line, "\t")[[1]]
        }
        gsub('^"|"$', "", raw_headers)
    }, error = \(e) character(0))
}

#' I2. The feature-table preview of a first assay: headers, the detected format, the format's reader
#' and the columns it proposes (`loadLipidImportAssayPreview`,
#' `prepareMetabImportAssaySelectionState`).
#' @param path The assay file.
#' @param import The omic's `config$import` (`detect`, `readers`, `id_mapping_key`).
#' @return A list: `headers`, `detected_format`, `format_confidence`, `import_result` (the
#'   reader's), `assay_data`, and `updates` with the `id` and `annotation` choices and
#'   selections and the proposed `is_pattern` (NULL when none).
#' @noRd
previewOmicsImport <- function(path, import) {
    headers <- readOmicsImportHeaders(path)
    if (length(headers) == 0L) omicsImportAbort("Could not read headers from file")
    format_info <- import$detect(headers = headers, filename = basename(path))
    result <- omicsImportReader(import$readers, format_info$format)(path)
    # The readers propose the id column as `lipid_id` or `metabolite_id`: the mapping key's stem.
    id_key <- sub("_col$", "", import$id_mapping_key)
    is_pattern <- result$is_pattern
    list(headers = headers, detected_format = format_info$format,
        format_confidence = format_info$confidence, import_result = result,
        assay_data = result$data,
        updates = list(
            id = list(choices = headers, selected = result$detected_columns[[id_key]]),
            annotation = list(choices = c("(None)" = "", headers),
                selected = result$detected_columns$annotation),
            is_pattern = if (!is.null(is_pattern) && !is.na(is_pattern)) is_pattern else NULL))
}

#' I3. The sample columns: for a custom format, the columns matching the pattern (when any match);
#' else the reader's sample columns, or the numeric columns; then, given `exclude_pattern`, those
#' not matching it (`resolveLipidImportSampleColumns`, `resolveMetabImportSampleColumns`).
#' @param exclude_pattern NULL, or `OMICS_IMPORT_NORMALISED_COLUMNS` for lipid's option.
#' @noRd
resolveOmicsImportSampleColumns <- function(assay_data, vendor_format, pattern = NULL,
                                            import_result = NULL, exclude_pattern = NULL) {
    if (identical(vendor_format, "custom") && !is.null(pattern) && nzchar(pattern)) {
        matched <- names(assay_data)[grepl(pattern, names(assay_data), ignore.case = TRUE)]
        if (length(matched) > 0) return(matched)
    }
    sample_columns <- if (!is.null(import_result)) {
        import_result$sample_columns
    } else {
        names(assay_data)[vapply(assay_data, is.numeric, logical(1))]
    }
    if (!is.null(exclude_pattern)) {
        excluded <- grepl(exclude_pattern, sample_columns, ignore.case = TRUE)
        sample_columns <- sample_columns[!excluded]
    }
    sample_columns
}

#' I4. The id or annotation column: the selected one, or for a custom format the typed name matched
#' case-insensitively to a header (kept as typed when none matches).
#' @noRd
resolveOmicsImportColumn <- function(headers, vendor_format, selected, custom) {
    if (!identical(vendor_format, "custom")) return(selected)
    if (is.null(custom) || !nzchar(custom)) return(custom)
    index <- match(tolower(custom), tolower(headers))
    if (is.na(index)) custom else headers[[index]]
}

# ------------------------------------------------------------------------------------------------
# I5 to I7: assay names, assembly, sanitising.
# ------------------------------------------------------------------------------------------------

#' Whether the second assay counts: a file and a non-empty name.
#' @noRd
omicsImportSecondAssayCounts <- function(assay2_file, assay2_name) {
    !is.null(assay2_file) && omicsImportScalarString(assay2_name)
}

#' I5. The assay-name rule (decision D3, from `config$import$assay_names`). `require_pair` (lipid,
#' trx): a primary name, no path separator, a second name exactly when a second file, two distinct
#' names. `skip_incomplete` (metab): a second assay counts only with a file and a name; the names
#' are non-empty and distinct; the path rule runs in the validation.
#' @return TRUE, invisibly; otherwise an `omics_import_invalid` error with the configured message.
#' @noRd
validateOmicsImportAssayNames <- function(assay1_name, assay2_name, assay2_file, rules) {
    messages <- rules$messages
    assay1_name <- omicsImportScalarValue(assay1_name)
    assay2_name <- omicsImportScalarValue(assay2_name)
    if (identical(rules$rule, "skip_incomplete")) {
        counted <- c(assay1_name,
            if (omicsImportSecondAssayCounts(assay2_file, assay2_name)) assay2_name)
        if (any(!nzchar(counted))) omicsImportAbort(messages[["required"]])
        if (anyDuplicated(counted)) {
            omicsImportAbort(paste0(messages[["duplicate"]],
                paste(unique(counted[duplicated(counted)]), collapse = ", ")))
        }
        return(invisible(TRUE))
    }
    assay2_file <- omicsImportScalarValue(assay2_file)
    has_file <- omicsImportScalarString(assay2_file)
    has_name <- omicsImportScalarString(assay2_name)
    if (!omicsImportScalarString(assay1_name)) omicsImportAbort(messages[["required"]])
    if (any(grepl(rules$path_pattern, c(assay1_name, if (has_name) assay2_name)))) {
        omicsImportAbort(messages[["path"]])
    }
    if (has_file != has_name) omicsImportAbort(messages[["pair"]])
    if (has_file && identical(assay1_name, assay2_name)) omicsImportAbort(messages[["duplicate"]])
    invisible(TRUE)
}

#' I6. The named assay list: the first assay, and the second read by the configured reader for the
#' format (given the id and annotation columns when `second_assay_columns` and the reader takes
#' them).
#' @param import The omic's `config$import`.
#' @param format The format the second assay's reader is chosen by.
#' @noRd
assembleOmicsImportAssays <- function(assay1_name, assay1_data, assay2_file, assay2_name, import,
                                      format, id_col = NULL, annotation_col = NULL) {
    assays <- list()
    assays[[assay1_name]] <- assay1_data
    if (!isTRUE(import$second_assay) || !omicsImportSecondAssayCounts(assay2_file, assay2_name)) {
        return(assays)
    }
    reader <- omicsImportReader(import$readers, format)
    args <- list(assay2_file)
    if (isTRUE(import$second_assay_columns)) {
        reader_args <- names(formals(reader))
        if ("lipid_id_column" %in% reader_args) args$lipid_id_column <- id_col
        if ("annotation_column" %in% reader_args) args$annotation_column <- annotation_col
    }
    assays[[assay2_name]] <- do.call(reader, args)$data
    assays
}

#' I7. Sample names made R-safe (`janitor::make_clean_names`) in every assay, when asked.
#' @return A list: `assays`, `sample_columns` and `sanitized`.
#' @noRd
sanitizeOmicsImportSampleNames <- function(assays, sample_columns, sanitize) {
    if (!isTRUE(sanitize)) {
        return(list(assays = assays, sample_columns = sample_columns, sanitized = FALSE))
    }
    cleaned <- janitor::make_clean_names(sample_columns)
    assays <- lapply(assays, \(assay) {
        index <- match(sample_columns, names(assay))
        found <- !is.na(index)
        if (any(found)) names(assay)[index[found]] <- cleaned[found]
        assay
    })
    list(assays = assays, sample_columns = cleaned, sanitized = TRUE)
}

# ------------------------------------------------------------------------------------------------
# I8, I9: the numeric-like rule, the conversion and the validation.
# ------------------------------------------------------------------------------------------------

#' Whether values are numeric-like: numeric and finite or NA; or text (or factor) whose every value,
#' trimmed, is empty, NA or a finite number (`metabImportValuesAreNumericLike`).
#' @noRd
omicsImportValuesAreNumericLike <- function(values) {
    if (is.numeric(values)) return(all(is.na(values) | is.finite(values)))
    if (is.factor(values)) values <- as.character(values)
    if (!is.character(values)) return(FALSE)
    trimmed <- trimws(values)
    missing_like <- is.na(trimmed) | !nzchar(trimmed)
    coerced <- suppressWarnings(as.numeric(trimmed))
    all(missing_like | (!is.na(coerced) & is.finite(coerced)))
}

#' I9. The shared parity conversion (KIT-102): each named column that is text or factor and
#' numeric-like is converted with `as.numeric` (empty text becomes NA). Other columns are not
#' touched; the omic's validation refuses them. A table with nothing to convert is returned as is.
#' @param table A data frame.
#' @param columns The column names to consider; absent names are ignored.
#' @return `table`, with the numeric-like text columns as doubles.
#' @noRd
coerceOmicsNumericLikeColumns <- function(table, columns) {
    for (column in intersect(columns, names(table))) {
        values <- table[[column]]
        is_text <- is.character(values) || is.factor(values)
        if (is_text && omicsImportValuesAreNumericLike(values)) {
            table[[column]] <- suppressWarnings(as.numeric(trimws(as.character(values))))
        }
    }
    table
}

#' I8. Validate the assembled assays with the omic's checks and messages (`omicsImportRules()`):
#' the list and its names, the id and sample columns, each assay through the configured mapping
#' validator, and the numeric check of the sample columns. Metab's assay-name path rule runs here.
#' @param validate_mapping `config$import$validate_mapping`, called as (data, id column, samples).
#' @param assay_names `config$import$assay_names`.
#' @return A list: `valid`, `errors`, `warnings` and `summary` (each assay's mapping summary).
#' @noRd
validateOmicsImportTables <- function(assays, id_col, sample_columns, validate_mapping, rules,
                                      assay_names) {
    messages <- rules$messages
    assay_names_seen <- names(assays)
    errors <- character()
    warnings <- character()
    summary <- list()
    if (!is.list(assays) || length(assays) == 0L) errors <- c(errors, messages[["no_assays"]])
    if (is.null(assay_names_seen) || any(!nzchar(assay_names_seen))) {
        errors <- c(errors, messages[["unnamed"]])
    }
    if (anyDuplicated(assay_names_seen) > 0L) {
        repeated <- unique(assay_names_seen[duplicated(assay_names_seen)])
        errors <- c(errors, if (grepl("%s", messages[["duplicate"]], fixed = TRUE)) {
            sprintf(messages[["duplicate"]], paste(repeated, collapse = ", "))
        } else {
            messages[["duplicate"]]
        })
    }
    if (identical(assay_names$rule, "skip_incomplete") &&
        any(grepl(assay_names$path_pattern, assay_names_seen))) {
        errors <- c(errors, assay_names$messages[["path"]])
    }
    if ("no_id" %in% names(messages) && (is.null(id_col) || !nzchar(id_col))) {
        errors <- c(errors, messages[["no_id"]])
    }
    if (is.null(sample_columns) || length(sample_columns) == 0L) {
        errors <- c(errors, messages[["no_samples"]])
    }
    is_numeric <- if (rules$numeric_like) omicsImportValuesAreNumericLike else is.numeric
    for (name in assay_names_seen) {
        assay <- assays[[name]]
        if (!is.data.frame(assay)) {
            errors <- c(errors, sprintf(messages[["not_data_frame"]], name))
            next
        }
        mapping <- validate_mapping(assay, id_col, sample_columns)
        if (rules$assay_prefix) {
            if (!isTRUE(mapping$valid)) {
                errors <- c(errors, sprintf("Assay '%s': %s", name,
                    paste(mapping$errors, collapse = "; ")))
            }
            if (length(mapping$warnings) > 0L) {
                warnings <- c(warnings, sprintf("Assay '%s': %s", name, mapping$warnings))
            }
        } else {
            errors <- c(errors, mapping$errors)
            warnings <- c(warnings, mapping$warnings)
        }
        summary[[name]] <- mapping$summary
        present <- intersect(sample_columns, names(assay))
        non_numeric <- present[!vapply(assay[present], is_numeric, logical(1))]
        if (length(non_numeric) > 0L) {
            errors <- c(errors, sprintf(messages[["non_numeric"]], name,
                paste(non_numeric, collapse = ", ")))
        }
    }
    if (rules$unique) {
        errors <- unique(errors)
        warnings <- unique(warnings)
    }
    list(valid = length(errors) == 0L, errors = errors, warnings = warnings, summary = summary)
}

# ------------------------------------------------------------------------------------------------
# I10 to I12: the column mapping, the artifacts and the processing log.
# ------------------------------------------------------------------------------------------------

#' I10. The column mapping, under the omic's id key (`lipid_id_col` for lipid and trx,
#' `metabolite_id_col` for metab).
#' @noRd
buildOmicsImportColumnMapping <- function(id_key, id_col, annotation_col, sample_columns,
                                          is_pattern) {
    stats::setNames(list(id_col,
        if (omicsImportScalarString(annotation_col)) annotation_col else NULL,
        sample_columns,
        if (omicsImportScalarString(is_pattern)) is_pattern else NA_character_),
        c(id_key, "annotation_col", "sample_columns", "is_pattern"))
}

#' I11. The import's source artifacts in `experiment_paths$source_dir`: each assay as
#' `data_cln_<assay>.tab`, the manifest, the column mapping, copies of the source files and the
#' omic's summary file (`writeLipidImportSourceArtifacts`, `writeMetabImportSourceArtifacts`).
#' @param source_files A named list of source paths; NULL entries are skipped.
#' @param summary_file `config$import$summary_file`.
#' @return A list: `written`, and `paths` and `source_copies`; or `written = FALSE` with a reason.
#' @noRd
writeOmicsImportArtifacts <- function(assays, column_mapping, experiment_paths, source_files,
                                      summary_file) {
    source_dir <- experiment_paths$source_dir
    if (is.null(source_dir) || !nzchar(source_dir) || !dir.exists(source_dir)) {
        return(list(written = FALSE, reason = "source_dir unavailable"))
    }
    assay_names <- names(assays)
    paths <- list()
    for (name in assay_names) {
        path <- file.path(source_dir, paste0("data_cln_", name, ".tab"))
        utils::write.table(assays[[name]], file = path, sep = "\t", row.names = FALSE,
            quote = FALSE)
        paths[[paste0("data_cln_", name)]] <- path
    }
    paths$assay_manifest <- file.path(source_dir, "assay_manifest.txt")
    writeLines(assay_names, paths$assay_manifest)
    paths$column_mapping <- file.path(source_dir, "column_mapping.json")
    jsonlite::write_json(column_mapping, paths$column_mapping, auto_unbox = TRUE, pretty = TRUE,
        null = "null")
    copy_dir <- file.path(source_dir, "import_sources")
    dir.create(copy_dir, recursive = TRUE, showWarnings = FALSE)
    copies <- character()
    source_files <- Filter(Negate(is.null), as.list(source_files))
    for (name in names(source_files)) {
        source <- source_files[[name]]
        if (omicsImportScalarString(source) && file.exists(source)) {
            destination <- file.path(copy_dir, paste0(name, "_", basename(source)))
            file.copy(source, destination, overwrite = TRUE)
            copies <- c(copies, destination)
        }
    }
    summary_rows <- data.frame(assay_name = assay_names,
        feature_count = vapply(assays, nrow, integer(1)),
        sample_count = length(column_mapping$sample_columns),
        sample_columns = paste(column_mapping$sample_columns, collapse = ","),
        stringsAsFactors = FALSE)
    paths$import_summary <- file.path(source_dir, summary_file)
    utils::write.table(summary_rows, file = paths$import_summary, sep = "\t", row.names = FALSE,
        quote = FALSE)
    list(written = TRUE, paths = paths, source_copies = copies)
}

#' I12. The import's processing log, with the feature count per assay under the omic's key
#' (`n_lipids`, `n_metabolites`): the unique ids when the id column is present, else the rows.
#' @noRd
describeOmicsImport <- function(assays, data_format, id_col, sample_columns, count_key,
                                artifacts, timestamp = Sys.time()) {
    counts <- vapply(assays, \(assay) {
        if (id_col %in% names(assay)) length(unique(assay[[id_col]])) else nrow(assay)
    }, integer(1))
    log <- list(timestamp = timestamp, n_assays = length(assays), assay_names = names(assays),
        detected_format = data_format)
    log[[count_key]] <- counts
    c(log, list(n_samples = length(sample_columns), artifacts = artifacts))
}

# ------------------------------------------------------------------------------------------------
# I13: the orchestration.
# ------------------------------------------------------------------------------------------------

#' I13. Import a feature table (lipid, metab, trx): the assay-name rule, assembly, sanitising, the
#' numeric-like conversion of the sample columns, the omic's `check` plugin (trx's count policy),
#' validation, the column mapping, the artifacts and the processing log. Nothing is committed to the
#' workflow: the GUI binds the owner, holds the tables and sets the status.
#' @param request A list: `assay1_name`, `assay1_data`, `assay1_file`, `assay2_file`,
#'   `assay2_name`, `vendor_format`, `detected_format`, `id_col`, `annotation_col`,
#'   `sample_columns` (ignored for trx, whose samples are every column but the id and name),
#'   `sanitize_names`, `is_pattern`, `experiment_paths`.
#' @param config The omic's `omicGuiConfig()`.
#' @return A list: `assays`, `sample_columns`, `data_format`, `column_mapping`, `processing_log`,
#'   `validation`, `artifacts` and `sanitized`. A refusal is an `omics_import_invalid` error (trx:
#'   also `trx_counts_invalid`) with today's message; the check plugin's and readers' errors pass
#'   through unchanged.
#' @noRd
importOmicsTables <- function(request, config) {
    import <- config$import
    rules <- omicsImportRules(config$shell$omic_type)
    id_col <- request$id_col
    annotation_col <- request$annotation_col
    second <- isTRUE(import$second_assay)
    assay2_file <- if (second) request$assay2_file
    assay2_name <- if (second) request$assay2_name
    validateOmicsImportAssayNames(request$assay1_name, assay2_name, assay2_file, import$assay_names)
    custom <- identical(request$vendor_format, "custom")
    reader_format <- if (custom) {
        request$detected_format
    } else {
        request$detected_format %||% request$vendor_format
    }
    assays <- assembleOmicsImportAssays(request$assay1_name, request$assay1_data, assay2_file,
        assay2_name, import, reader_format, id_col, annotation_col)
    sample_columns <- if (identical(rules$samples, "all_but_id")) {
        setdiff(names(request$assay1_data), c(id_col, annotation_col))
    } else {
        request$sample_columns
    }
    sanitized <- sanitizeOmicsImportSampleNames(assays, sample_columns, request$sanitize_names)
    sample_columns <- sanitized$sample_columns
    assays <- lapply(sanitized$assays, coerceOmicsNumericLikeColumns, columns = sample_columns)
    if (is.function(import$check)) import$check(assays, id_col, annotation_col)
    validation <- validateOmicsImportTables(assays, id_col, sample_columns, import$validate_mapping,
        rules, import$assay_names)
    if (!isTRUE(validation$valid)) {
        suffix <- if (length(rules$refusal_class)) {
            sprintf(" (%s)", rules$refusal_class[[1]])
        } else {
            ""
        }
        omicsImportAbort(paste0(import$refusal_prefix, ": ",
            paste(validation$errors, collapse = "; "), suffix), rules$refusal_class)
    }
    data_format <- if (custom) "custom" else request$detected_format
    column_mapping <- buildOmicsImportColumnMapping(import$id_mapping_key, id_col, annotation_col,
        sample_columns, request$is_pattern)
    source_files <- list(assay1 = request$assay1_file,
        assay2 = if (omicsImportSecondAssayCounts(assay2_file, assay2_name)) assay2_file)
    artifacts <- writeOmicsImportArtifacts(assays, column_mapping, request$experiment_paths,
        source_files, import$summary_file)
    list(assays = assays, sample_columns = sample_columns, data_format = data_format,
        column_mapping = column_mapping,
        processing_log = describeOmicsImport(assays, data_format, id_col, sample_columns,
            import$count_key, artifacts),
        validation = validation, artifacts = artifacts, sanitized = sanitized$sanitized)
}
