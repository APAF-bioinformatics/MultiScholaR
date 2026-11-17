# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# KIT-061: the per-omic integration input (KIT-060 audit section 4; the user's decisions D1, D2). Each omic's Summary
# tab writes one folder, `<project>/integration/<omic>/`, of plain tables under fixed names with a JSON manifest; the
# integration reads it through `readIntegrationInput()` and never an omic's store or S4 class.
#
#   integration_input.json      the manifest (INTEGRATION_INPUT_SCHEMA)
#   matrix_<view>.tsv           features x samples; first column `feature_id`, then one column per sample id
#   features_<view>.tsv         one row per feature, in the matrix's order: `feature_id`, `display_name`, annotation
#   samples.tsv                 the design of the matrix's samples, in the matrix's column order
#
# The manifest carries the omic and its route, the matrix's role (`corrected` or `normalised`), the RUV design that
# produced a corrected matrix, the source state's entity id and hash, every file's sha256, the package version and the
# time of writing.

INTEGRATION_INPUT_SCHEMA <- "multischolar.integration_input.v1"
INTEGRATION_INPUT_MANIFEST <- "integration_input.json"
INTEGRATION_INPUT_SAMPLES <- "samples.tsv"
INTEGRATION_INPUT_ROLES <- c("corrected", "normalised")

#' Refuse an integration input with the contract's classed error.
#' @param message What failed, lower case.
#' @noRd
abortIntegrationInput <- function(message) {
    cli::cli_abort(message, class = "integration_input_invalid", call = NULL)
}

#' A view's two file names (fixed by the view name).
#' @param view The view name.
#' @return `list(matrix, features)`.
#' @noRd
integrationInputViewFiles <- function(view) {
    list(matrix = sprintf("matrix_%s.tsv", view), features = sprintf("features_%s.tsv", view))
}

#' Validate the contract's tables against the manifest's column names and views: the matrix's first column is
#' `feature_id` with unique, non-missing ids and numeric sample columns; the features file lists the same ids in the
#' same order with a display name; the samples file lists exactly the matrix's samples, once each, with a group.
#' @param manifest The manifest (its `views`, `sample_id_column`, `group_column`, `feature_id_column`,
#'   `display_name_column`).
#' @param samples The samples table.
#' @param matrices,features Named lists of tables by view.
#' @return TRUE invisibly, or an `integration_input_invalid` error.
#' @noRd
validateIntegrationInputTables <- function(manifest, samples, matrices, features) {
    sample_column <- manifest$sample_id_column
    group_column <- manifest$group_column
    feature_column <- manifest$feature_id_column
    display_column <- manifest$display_name_column
    if (!is.data.frame(samples) || !all(c(sample_column, group_column) %in% names(samples))) {
        abortIntegrationInput(sprintf("the samples table needs the columns %s and %s", sample_column, group_column))
    }
    sample_ids <- as.character(samples[[sample_column]])
    if (!length(sample_ids) || anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
        abortIntegrationInput("the sample ids must be present, non-empty and unique")
    }
    views <- vapply(manifest$views, function(view) as.character(view$name), character(1))
    if (!length(views) || !setequal(views, names(matrices)) || !setequal(views, names(features))) {
        abortIntegrationInput("every view needs one matrix and one features table")
    }
    for (view in views) {
        matrix <- matrices[[view]]
        if (!is.data.frame(matrix) || ncol(matrix) < 2L || !identical(names(matrix)[[1L]], feature_column)) {
            abortIntegrationInput(sprintf("the matrix of view %s must start with the column %s", view, feature_column))
        }
        ids <- as.character(matrix[[feature_column]])
        if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
            abortIntegrationInput(sprintf("the feature ids of view %s must be present, non-empty and unique", view))
        }
        columns <- names(matrix)[-1L]
        if (anyDuplicated(columns) || !setequal(columns, sample_ids) || length(columns) != length(sample_ids)) {
            abortIntegrationInput(sprintf("the samples of view %s must be the samples table's, once each", view))
        }
        if (!all(vapply(matrix[columns], is.numeric, logical(1)))) {
            abortIntegrationInput(sprintf("the sample columns of view %s must be numeric", view))
        }
        annotation <- features[[view]]
        if (!is.data.frame(annotation) || !all(c(feature_column, display_column) %in% names(annotation)) ||
            !identical(as.character(annotation[[feature_column]]), ids)) {
            abortIntegrationInput(sprintf("the features of view %s must list the matrix's ids in order with %s",
                view, display_column))
        }
    }
    invisible(TRUE)
}

#' Write one omic's integration input (KIT-061): the tables under their fixed names, then the manifest with every
#' file's sha256, in a staging folder that replaces `dir` once complete.
#' @param dir The omic's folder, `<project>/integration/<omic>/`.
#' @param manifest The omic's fields: `omic_type`, `workflow`, `matrix_role`, `scale`, `correction` (NULL when
#'   `normalised`), `source_state`, `annotation`, `sample_id_column`, `group_column`, `technical_replicate_column`.
#' @param samples The design of the matrix's samples.
#' @param matrices Named list by view of features x samples tables whose first column is `feature_id`.
#' @param features Named list by view of feature tables (`feature_id`, `display_name`, annotation columns).
#' @param now The time of writing.
#' @return The manifest's path, invisibly.
#' @noRd
writeIntegrationInput <- function(dir, manifest, samples, matrices, features, now = Sys.time()) {
    views <- names(matrices)
    if (!length(views) || any(!grepl("^[A-Za-z0-9_-]+$", views))) {
        abortIntegrationInput("view names must be letters, digits, underscores or hyphens")
    }
    if (!isTRUE(manifest$matrix_role %in% INTEGRATION_INPUT_ROLES)) {
        abortIntegrationInput("the matrix role must be corrected or normalised")
    }
    manifest$schema <- INTEGRATION_INPUT_SCHEMA
    manifest$feature_id_column <- "feature_id"
    manifest$display_name_column <- "display_name"
    manifest$views <- lapply(views, function(view) {
        files <- integrationInputViewFiles(view)
        list(name = view, matrix_file = files$matrix, features_file = files$features,
            rows = nrow(matrices[[view]]), columns = ncol(matrices[[view]]) - 1L)
    })
    sample_ids <- as.character(samples[[manifest$sample_id_column]])
    validateIntegrationInputTables(manifest, samples, matrices, features)
    parent <- dirname(dir)
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    staging <- tempfile(paste0(".", basename(dir), "-"), tmpdir = parent)
    dir.create(staging)
    on.exit(unlink(staging, recursive = TRUE), add = TRUE)
    write <- function(table, name) {
        vroom::vroom_write(table, file.path(staging, name), delim = "\t", na = "NA", progress = FALSE)
        name
    }
    written <- write(samples, INTEGRATION_INPUT_SAMPLES)
    for (view in views) {
        files <- integrationInputViewFiles(view)
        # The matrix's columns in the samples table's order.
        written <- c(written, write(matrices[[view]][c("feature_id", sample_ids)], files$matrix),
            write(features[[view]], files$features))
    }
    manifest$samples_file <- INTEGRATION_INPUT_SAMPLES
    manifest$files <- stats::setNames(lapply(written, function(name) {
        digest::digest(file = file.path(staging, name), algo = "sha256")
    }), written)
    manifest$software <- list(package = "MultiScholaR", version = as.character(utils::packageVersion("MultiScholaR")))
    manifest$written_at <- format(now, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    jsonlite::write_json(manifest, file.path(staging, INTEGRATION_INPUT_MANIFEST), auto_unbox = TRUE, pretty = TRUE,
        digits = NA, null = "null", na = "null")
    unlink(dir, recursive = TRUE)
    if (!file.rename(staging, dir)) abortIntegrationInput(sprintf("the integration input could not be published at %s", dir))
    invisible(file.path(dir, INTEGRATION_INPUT_MANIFEST))
}

#' Read and validate one omic's integration input (KIT-061): the manifest, then every file it lists against its
#' sha256, then the tables' shapes and ids against the manifest and each other.
#' @param dir The omic's folder, `<project>/integration/<omic>/`.
#' @return `list(manifest, samples, views)`; `views` is named by view, each `list(matrix, features)`. An input that
#'   fails any check is refused with an `integration_input_invalid` error.
#' @noRd
readIntegrationInput <- function(dir) {
    path <- file.path(dir, INTEGRATION_INPUT_MANIFEST)
    if (!file.exists(path)) abortIntegrationInput(sprintf("no integration input manifest at %s", path))
    manifest <- tryCatch(jsonlite::read_json(path, simplifyVector = FALSE), error = function(error) {
        abortIntegrationInput(sprintf("the integration input manifest cannot be read: %s", conditionMessage(error)))
    })
    if (!identical(manifest$schema, INTEGRATION_INPUT_SCHEMA)) {
        abortIntegrationInput(sprintf("the manifest's schema must be %s", INTEGRATION_INPUT_SCHEMA))
    }
    required <- c("omic_type", "workflow", "matrix_role", "source_state", "sample_id_column", "group_column",
        "feature_id_column", "display_name_column", "views", "samples_file", "files", "software", "written_at")
    missing <- setdiff(required, names(manifest))
    if (length(missing)) abortIntegrationInput(sprintf("the manifest lacks %s", paste(missing, collapse = ", ")))
    if (!isTRUE(manifest$matrix_role %in% INTEGRATION_INPUT_ROLES)) {
        abortIntegrationInput("the matrix role must be corrected or normalised")
    }
    if (identical(manifest$matrix_role, "corrected") && !is.list(manifest$correction)) {
        abortIntegrationInput("a corrected matrix needs the correction that produced it")
    }
    files <- manifest$files
    if (!is.list(files) || !length(files) || any(!grepl("^[A-Za-z0-9_.-]+$", names(files)))) {
        abortIntegrationInput("the manifest's files must be plain file names with their sha256")
    }
    for (name in names(files)) {
        file <- file.path(dir, name)
        if (!file.exists(file)) abortIntegrationInput(sprintf("the file %s is missing", name))
        if (!identical(digest::digest(file = file, algo = "sha256"), files[[name]])) {
            abortIntegrationInput(sprintf("the file %s does not match its sha256", name))
        }
    }
    listed <- c(manifest$samples_file, unlist(lapply(manifest$views, function(view) c(view$matrix_file, view$features_file))))
    if (!all(listed %in% names(files))) abortIntegrationInput("every table must be listed with its sha256")
    read <- function(name, types) {
        table <- suppressWarnings(readr::read_tsv(file.path(dir, name), col_types = types, na = "NA",
            name_repair = "minimal", progress = FALSE, show_col_types = FALSE))
        if (nrow(readr::problems(table))) abortIntegrationInput(sprintf("the file %s cannot be parsed", name))
        as.data.frame(table, stringsAsFactors = FALSE, optional = TRUE)
    }
    text <- readr::cols(.default = readr::col_character())
    samples <- read(manifest$samples_file, text)
    matrices <- list()
    features <- list()
    for (view in manifest$views) {
        matrix <- read(view$matrix_file, text)
        if (!identical(nrow(matrix), as.integer(view$rows)) || !identical(ncol(matrix) - 1L, as.integer(view$columns))) {
            abortIntegrationInput(sprintf("the matrix of view %s is not %s x %s", view$name, view$rows, view$columns))
        }
        # Numbers are parsed by R's own conversion, which reads back the written digits exactly (readr's parser can
        # differ in the last bit).
        for (column in seq_along(matrix)[-1L]) {
            values <- suppressWarnings(as.numeric(matrix[[column]]))
            if (any(is.na(values) & !is.nan(values) & !is.na(matrix[[column]]))) {
                abortIntegrationInput(sprintf("the matrix of view %s has a value that is not a number", view$name))
            }
            matrix[[column]] <- values
        }
        matrices[[view$name]] <- matrix
        features[[view$name]] <- read(view$features_file, text)
    }
    validateIntegrationInputTables(manifest, samples, matrices, features)
    list(manifest = manifest, samples = samples, views = stats::setNames(lapply(names(matrices), function(view) {
        list(matrix = matrices[[view]], features = features[[view]])
    }), names(matrices)))
}
