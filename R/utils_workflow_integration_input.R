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
# tab writes one folder, `<project>/integration/<omic>/`, of tables under fixed names with a JSON manifest; the
# integration reads it through `readIntegrationInput()` and never an omic's store or S4 class.
#
# KIT-147 (plan item A, contract v2): the tables are Parquet files.
#
#   integration_input.json      the manifest (INTEGRATION_INPUT_SCHEMA)
#   matrix_<view>.parquet       features x samples, wide: `feature_id` (text), then one double column per sample id in
#                               the samples table's order; row groups of 10,000 rows
#   features_<view>.parquet     one row per feature, in the matrix's order: `feature_id`, `display_name`, annotation
#   samples.parquet             the design of the matrix's samples, in the matrix's column order
#
# Every Parquet file is written with format version 2.6, snappy compression and row groups of at most 10,000 rows; the
# samples and features tables hold text columns (as the v1 reader returned them), the matrix text ids and doubles. The
# manifest carries the omic and its route, the matrix's role (`corrected` or `normalised`), the RUV design that produced
# a corrected matrix, the source state's entity id and hash, every file's sha256, the Parquet layout (`parquet`) and
# each view's row-group count, the package version and the time of writing. A v1 folder (TSV tables, schema v1) was
# only ever written by development builds and is refused (`integration_input_outdated`).
#
# KIT-153 (INT-23): each data set's `matrix_<view>.parquet` is the normalised matrix before correction for unwanted
# variation (RUV), stage `normalised` (`matrix_stage` of the view; the manifest's `matrix_role` is `normalised`, the
# role of that file); where the omic applied RUV to the data set, the RUV-corrected matrix is written beside it as
#
#   corrected_<view>.parquet    the same features in the same order and the same samples, after RUV (the manifest's
#                               `correction` names it; per view `corrected_matrix_file` and `corrected_row_groups`)
#
# and the manifest's `normalisation` names the normalisation and the state or product that held the uncorrected
# matrix. A v2 reader (KIT-147) reads such a folder as what it is: a normalised matrix per data set, the corrected file
# one more file with its sha256. A v2 folder of the earlier build has no `matrix_stage`: it held one matrix per data set,
# the RUV-corrected one where the omic applied RUV, so it is refused in plain words (`integration_input_outdated`).
#
# KIT-184 (VAL-53 (1), VAL-58; an addition inside contract v2): where a data set has a corrected matrix, the unwanted
# factors the correction removed are written beside it as
#
#   removed_<view>.parquet      `feature_id` (`ruv_W1` to `ruv_Wk`), then one double column per run in the samples
#                               table's order: an orthonormal basis, orthogonal to the constant, of the row space of
#                               (normalised minus corrected) over the features measured in every run
#
# with, in the data set's manifest entry, `removed_factors_file` and `removed_factors` (`method`, `k_declared`, `rank`,
# `k`, `share`, `features`, `tolerance`, `status`, `reason`; `integrationInputRemovedFactors()`). No stored W of any RUV
# step is read: the factors are taken from what was removed. A data set without a corrected matrix has neither; one
# whose factors cannot be determined has the record with its reason and no file. An earlier build reads such a folder as
# before (one more listed file, two fields it does not name); this build reads a folder saved before it as a data set
# without removed factors.

INTEGRATION_INPUT_SCHEMA <- "multischolar.integration_input.v2"
INTEGRATION_INPUT_SCHEMA_V1 <- "multischolar.integration_input.v1"
INTEGRATION_INPUT_MANIFEST <- "integration_input.json"
INTEGRATION_INPUT_SAMPLES <- "samples.parquet"
INTEGRATION_INPUT_ROLES <- c("corrected", "normalised")
# KIT-153: the stage of each data set's `matrix_file` (the normalised matrix before RUV).
INTEGRATION_INPUT_MATRIX_STAGE <- "normalised"
# The Parquet layout every v2 file follows, recorded in the manifest as `parquet` and required as written.
INTEGRATION_INPUT_PARQUET <- list(format_version = "2.6", compression = "snappy", row_group_rows = 10000L,
    column_types = list(feature_id = "string", runs = "double", samples = "string", features = "string"))
# KIT-184: the relative bound below which a singular value of the removed part, and the removed part itself, count as
# zero: far above double-precision rounding and far below any variation in data (ticket 5.4).
INTEGRATION_INPUT_REMOVED_TOLERANCE <- 1e-8

#' Refuse an integration input with the contract's classed error.
#' @param message What failed, lower case.
#' @param class A more specific class before `integration_input_invalid`, or NULL.
#' @noRd
abortIntegrationInput <- function(message, class = NULL) {
    cli::cli_abort(message, class = c(class, "integration_input_invalid"), call = NULL)
}

#' A view's file names (fixed by the view name). KIT-153: `corrected`, the RUV-corrected matrix where the omic has one.
#' KIT-184: `removed`, the unwanted factors the correction removed.
#' @param view The view name.
#' @return `list(matrix, features, corrected, removed)`.
#' @noRd
integrationInputViewFiles <- function(view) {
    list(matrix = sprintf("matrix_%s.parquet", view), features = sprintf("features_%s.parquet", view),
        corrected = sprintf("corrected_%s.parquet", view), removed = sprintf("removed_%s.parquet", view))
}

#' The number of factors a correction declares for a data set (KIT-184 5.4): `correction$views[[view]]$k` where the
#' correction lists data sets (lipidomics, metabolomics), else `correction$k`.
#' @param correction The manifest's `correction`, or NULL.
#' @param view The data set's name.
#' @return One whole number of at least 1, or `NA_integer_` when the correction declares none.
#' @noRd
integrationInputDeclaredK <- function(correction, view) {
    k <- if (is.list(correction)) {
        per_view <- if (is.list(correction$views)) correction$views[[view]]
        if (is.list(per_view) && !is.null(per_view$k)) per_view$k else correction$k
    }
    if (!is.numeric(k) || length(k) != 1L || is.na(k) || k != round(k) || k < 1) return(NA_integer_)
    as.integer(k)
}

#' The unwanted factors a correction removed from one data set (KIT-184 5.4; VAL-53 (1), VAL-58): an orthonormal basis
#' of the row space of (normalised minus corrected), each row centred on its mean, over the features with a finite
#' value in every run of both matrices. On every RUV route the correction of such a feature subtracts a vector in the
#' span of one k-column W per data set (plus a constant where the feature's level is restored), so that row space is the
#' centred W's span and has dimension k or less. The numerical rank r counts the singular values above
#' `INTEGRATION_INPUT_REMOVED_TOLERANCE` times the largest; min(k, r) factors are written, each with its entry of
#' largest absolute value positive. `exact`: r is at most k (everything removed from those features lies in the
#' factors). `approximate`: r is above k (a run removed after RUV left features that were corrected one by one among the
#' complete ones); `share` is then the part of the removed sum of squares the factors hold. `not_determined`: no
#' declared k, fewer than k + 1 complete features, or a difference that is zero or a constant per feature.
#' @param uncorrected,corrected The normalised and the corrected matrix: `feature_id`, then one column per run (any
#'   column order; the same features in the same order).
#' @param k The correction's declared number of factors (`integrationInputDeclaredK()`).
#' @param sample_ids The samples table's run ids, in its order.
#' @return `list(table, record)`: `table` the factors in the contract's matrix layout (`feature_id` `ruv_W1` to
#'   `ruv_Wk`, then the runs in `sample_ids`' order), NULL when not determined; `record` `list(method, k_declared, rank,
#'   k, share, features, tolerance, status, reason)`.
#' @noRd
integrationInputRemovedFactors <- function(uncorrected, corrected, k, sample_ids) {
    tolerance <- INTEGRATION_INPUT_REMOVED_TOLERANCE
    declared <- if (is.numeric(k) && length(k) == 1L && !is.na(k) && k == round(k) && k >= 1) as.integer(k) else NA_integer_
    u <- as.matrix(as.data.frame(uncorrected, optional = TRUE)[sample_ids])
    v <- as.matrix(as.data.frame(corrected, optional = TRUE)[sample_ids])
    storage.mode(u) <- "double"
    storage.mode(v) <- "double"
    complete <- rowSums(!is.finite(u)) == 0L & rowSums(!is.finite(v)) == 0L
    features <- sum(complete)
    undetermined <- function(reason) list(table = NULL, record = list(method = "difference_svd", k_declared = declared,
        rank = NA_integer_, k = 0L, share = NA_real_, features = features, tolerance = tolerance, status = "not_determined",
        reason = reason))
    if (is.na(declared)) return(undetermined("the correction records no number of removed factors"))
    if (features < declared + 1L) {
        # KIT-255: a run with no finite value (set missing by the sample-correlation filter, or empty before) is named.
        empty <- sample_ids[colSums(is.finite(u)) == 0L | colSums(is.finite(v)) == 0L]
        if (length(empty)) {
            return(undetermined(sprintf("%s no value in this data set, and no feature is then measured in every run",
                if (length(empty) == 1L) paste("run", empty, "has") else paste("runs", paste(empty, collapse = ", "), "have"))))
        }
        return(undetermined(sprintf("%d features are measured in every run, fewer than the %d needed for %d factors",
            features, declared + 1L, declared)))
    }
    difference <- u[complete, , drop = FALSE] - v[complete, , drop = FALSE]
    difference <- difference - rowMeans(difference)
    if (max(abs(difference)) <= tolerance * max(abs(u[complete, , drop = FALSE]))) {
        return(undetermined(paste("the normalised and the corrected values are equal or differ by a constant per feature,",
            "so the correction removed no direction across the runs")))
    }
    decomposition <- svd(difference, nu = 0L)
    d <- decomposition$d
    rank <- sum(d > tolerance * d[[1L]])
    kept <- min(declared, rank)
    factors <- decomposition$v[, seq_len(kept), drop = FALSE]
    for (j in seq_len(kept)) if (factors[which.max(abs(factors[, j])), j] < 0) factors[, j] <- -factors[, j]
    table <- data.frame(feature_id = paste0("ruv_W", seq_len(kept)), t(factors), stringsAsFactors = FALSE)
    names(table) <- c("feature_id", sample_ids)
    list(table = table, record = list(method = "difference_svd", k_declared = declared, rank = as.integer(rank),
        k = as.integer(kept), share = sum(d[seq_len(kept)]^2) / sum(d^2), features = features, tolerance = tolerance,
        status = if (rank <= declared) "exact" else "approximate", reason = NULL))
}

#' The normalised matrix that entered RUV, on the corrected matrix's features and samples (KIT-153): its rows of the
#' corrected matrix's feature ids in that order, its columns `feature_id` and the corrected matrix's samples. Refused
#' when it lacks one of them (the matrix that entered RUV is then not the one found).
#' @param uncorrected The normalised matrix (`feature_id`, then one column per run).
#' @param corrected The corrected matrix, the same form.
#' @param view The data set's name (for the message).
#' @return The uncorrected matrix in the corrected one's shape.
#' @noRd
integrationInputUncorrectedOn <- function(uncorrected, corrected, view) {
    at <- match(as.character(corrected$feature_id), as.character(uncorrected$feature_id))
    runs <- names(corrected)[-1L]
    absent <- setdiff(runs, names(uncorrected))
    if (anyNA(at) || length(absent)) {
        abortIntegrationInput(sprintf(paste("the normalised matrix that entered RUV in data set %s lacks %d of the corrected",
            "matrix's features and %d of its runs"), view, sum(is.na(at)), length(absent)))
    }
    out <- as.data.frame(uncorrected, stringsAsFactors = FALSE, optional = TRUE)[at, c("feature_id", runs), drop = FALSE]
    out$feature_id <- as.character(out$feature_id)
    rownames(out) <- NULL
    out
}

#' A view's name for an omic with several assays (INT-4, user 2026-11-19): `<omic>_<assay>` (for example
#' `lipidomics_LCMS_Pos`), since lipidomics and metabolomics both have an assay `LCMS_Pos`; a character the contract's
#' view names do not take becomes an underscore.
#' @param omic The omic (`lipidomics`, `metabolomics`).
#' @param assay The assay's name.
#' @return The view name.
#' @noRd
integrationInputViewName <- function(omic, assay) {
    paste(omic, gsub("[^A-Za-z0-9_-]", "_", assay), sep = "_")
}

#' The design column that names each run's biological sample (INT-2, user 2026-11-19: technical replicates are averaged
#' per biological sample before the model, which is the model's step; the writer names the column): the design's
#' declared technical-replicate column, else `tech_rep_group` (the group and replicate number, the biological unit since
#' KIT-128), else the run's own id (each run its own biological sample).
#' KIT-222 (owner's OD-3): a design whose biological-sample names hold runs of more than one technical-replicate marking
#' action is refused (`refuseOmicsDesignReplicateKeyCollision`), as the integration would average two biological
#' samples and pair the average by name.
#' @param design The input's samples table.
#' @param technical The design's technical-replicate column, or NA.
#' @param sample The run id column.
#' @return A column name of `design`.
#' @noRd
integrationInputBiologicalSampleColumn <- function(design, technical, sample) {
    column <- if (is.character(technical) && length(technical) == 1L && isTRUE(technical %in% names(design))) {
        technical
    } else if ("tech_rep_group" %in% names(design)) {
        "tech_rep_group"
    } else {
        sample
    }
    if (!identical(column, sample)) {
        frame <- as.data.frame(design)
        if ("group" %in% names(frame)) {
            resolveProtTechnicalReplicateBlock(frame, "group", column)
        } else {
            refuseOmicsDesignReplicateKeyCollision(frame, key = column, run = sample)
        }
    }
    column
}

#' Validate the contract's tables against the manifest's column names and views: the matrix's first column is
#' `feature_id` with unique, non-missing ids and numeric sample columns; the features file lists the same ids in the
#' same order with a display name; the samples file lists exactly the matrix's samples, once each, with a group.
#' @param manifest The manifest (its `views`, `sample_id_column`, `group_column`, `feature_id_column`,
#'   `display_name_column`).
#' @param samples The samples table.
#' @param matrices,features Named lists of tables by view.
#' @param corrected KIT-153: named list by view of the RUV-corrected matrices (a subset of the views); each has its
#'   view's feature ids in the same order and its samples.
#' @return TRUE invisibly, or an `integration_input_invalid` error.
#' @noRd
validateIntegrationInputTables <- function(manifest, samples, matrices, features, corrected = list()) {
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
    if (length(corrected) && (is.null(names(corrected)) || !all(names(corrected) %in% views) || anyDuplicated(names(corrected)))) {
        abortIntegrationInput("every corrected matrix must be one of a view's")
    }
    for (view in names(corrected)) {
        matrix <- corrected[[view]]
        if (!is.data.frame(matrix) || !identical(names(matrix)[[1L]], feature_column) ||
            !identical(as.character(matrix[[feature_column]]), as.character(matrices[[view]][[feature_column]]))) {
            abortIntegrationInput(sprintf("the corrected matrix of view %s must list the view's feature ids in order", view))
        }
        columns <- names(matrix)[-1L]
        if (anyDuplicated(columns) || !setequal(columns, sample_ids) || length(columns) != length(sample_ids)) {
            abortIntegrationInput(sprintf("the samples of the corrected matrix of view %s must be the samples table's, once each",
                view))
        }
        if (!all(vapply(matrix[columns], is.numeric, logical(1)))) {
            abortIntegrationInput(sprintf("the sample columns of the corrected matrix of view %s must be numeric", view))
        }
    }
    invisible(TRUE)
}

# ---------------------------------------------------------------------------------------------------------------
# KIT-147: the Parquet files of contract v2.
# ---------------------------------------------------------------------------------------------------------------

#' A table as the contract's samples and features files hold it: a plain data frame of text columns (each value as
#' `as.character` gives it, missing values kept), as the v1 reader returned these tables.
#' @param table A data frame.
#' @noRd
integrationInputTextTable <- function(table) {
    table <- as.data.frame(table, stringsAsFactors = FALSE, optional = TRUE)
    table[] <- lapply(table, function(column) if (is.character(column) && is.null(attributes(column))) column else
        as.character(column))
    rownames(table) <- NULL
    table
}

#' A matrix as the contract's matrix file holds it: `feature_id` as text, then one double column per sample id in the
#' given order.
#' @param matrix A features x samples table whose first column is `feature_id`.
#' @param sample_ids The samples table's ids, in its order.
#' @noRd
integrationInputMatrixTable <- function(matrix, sample_ids) {
    table <- as.data.frame(matrix, stringsAsFactors = FALSE, optional = TRUE)[c("feature_id", sample_ids)]
    table$feature_id <- as.character(table$feature_id)
    table[sample_ids] <- lapply(table[sample_ids], as.double)
    rownames(table) <- NULL
    table
}

#' The zero-row matrix of a data set: `feature_id` as text, then one double column per run.
#' @param columns `feature_id`, then the run ids.
#' @noRd
integrationInputMatrixPrototype <- function(columns) {
    structure(stats::setNames(c(list(character()), rep(list(double()), length(columns) - 1L)), columns),
        class = "data.frame", row.names = integer())
}

#' Write one table as a contract Parquet file: format version 2.6, snappy, row groups of 10,000 rows, no R metadata
#' (the same table gives the same bytes on the same platform's Arrow build and the same content on every platform;
#' a file's sha256 identifies it on the machine that wrote it and travels with the project).
#' @param table A plain data frame of text or double columns.
#' @param path The file to write.
#' @noRd
writeIntegrationInputParquet <- function(table, path) {
    layout <- INTEGRATION_INPUT_PARQUET
    arrow_table <- arrow::Table$create(table)
    arrow_table <- arrow_table$ReplaceSchemaMetadata(list(multischolar_integration_input = INTEGRATION_INPUT_SCHEMA))
    sink <- arrow::FileOutputStream$create(path)
    on.exit(try(sink$close(), silent = TRUE), add = TRUE)
    writer <- arrow::ParquetFileWriter$create(arrow_table$schema, sink,
        properties = arrow::ParquetWriterProperties$create(names(table), version = layout$format_version,
            compression = layout$compression))
    writer$WriteTable(arrow_table, chunk_size = layout$row_group_rows)
    writer$Close()
    sink$close()
    invisible(path)
}

#' Open a contract Parquet file, run `read` on its reader and close it (a plain stream, closed on exit: KIT-131).
#' @param path The file.
#' @param name The file's name in messages.
#' @param read `function(reader)`.
#' @return `read`'s value.
#' @noRd
withIntegrationInputParquet <- function(path, name, read) {
    unreadable <- function(error) abortIntegrationInput(sprintf("the file %s cannot be read as Parquet: %s", name,
        conditionMessage(error)))
    stream <- tryCatch(arrow::ReadableFile$create(path), error = unreadable)
    on.exit(try(stream$close(), silent = TRUE), add = TRUE)
    reader <- tryCatch(arrow::ParquetFileReader$create(stream), error = unreadable)
    on.exit(try(reader$.unsafe_delete(), silent = TRUE), add = TRUE, after = FALSE)
    read(reader)
}

#' An Arrow table as a plain data frame.
#' @noRd
integrationInputFrame <- function(table) {
    frame <- as.data.frame(table)
    class(frame) <- "data.frame"
    frame
}

#' A Parquet file's column names and Arrow types.
#' @noRd
integrationInputSchema <- function(reader) {
    schema <- reader$GetSchema()
    list(names = schema$names, types = vapply(schema$fields, function(field) field$type$ToString(), character(1)))
}

#' Read a small contract table whole (the samples or a features file): every column text.
#' @param path The file.
#' @param name The file's name in messages.
#' @return A plain data frame.
#' @noRd
readIntegrationInputTextTable <- function(path, name) {
    withIntegrationInputParquet(path, name, function(reader) {
        schema <- integrationInputSchema(reader)
        if (!length(schema$names) || anyDuplicated(schema$names) || !all(schema$types == "string")) {
            abortIntegrationInput(sprintf("the file %s must hold uniquely named text columns", name))
        }
        integrationInputFrame(reader$ReadTable())
    })
}

#' Read one row group of a data set's matrix file (KIT-147): the rows of that group, every column.
#' @param path The matrix file.
#' @param name The file's name in messages.
#' @param index The row group, from 1.
#' @param columns The columns the file must have, in order (`feature_id`, then the runs).
#' @return A plain data frame of at most 10,000 rows.
#' @noRd
readIntegrationInputRowGroup <- function(path, name, index, columns) {
    withIntegrationInputParquet(path, name, function(reader) {
        if (!identical(integrationInputSchema(reader)$names, columns)) {
            abortIntegrationInput(sprintf("the file %s no longer has the columns it was validated with", name))
        }
        if (!(index >= 1L && index <= reader$num_row_groups)) {
            abortIntegrationInput(sprintf("the file %s has no row group %d", name, index))
        }
        page <- integrationInputFrame(reader$ReadRowGroup(as.integer(index) - 1L))
        if (nrow(page) > INTEGRATION_INPUT_PARQUET$row_group_rows) {
            abortIntegrationInput(sprintf("the file %s has a row group of more than %d rows", name,
                INTEGRATION_INPUT_PARQUET$row_group_rows))
        }
        page
    })
}

#' Read a data set's whole matrix (KIT-147): for small inputs and tests; the import streams row groups instead.
#' @inheritParams readIntegrationInputRowGroup
#' @noRd
readIntegrationInputMatrix <- function(path, name, columns) {
    withIntegrationInputParquet(path, name, function(reader) {
        if (!identical(integrationInputSchema(reader)$names, columns)) {
            abortIntegrationInput(sprintf("the file %s no longer has the columns it was validated with", name))
        }
        integrationInputFrame(reader$ReadTable())
    })
}

#' Check a data set's matrix file against the manifest without reading its values (KIT-147): the schema (`feature_id`
#' as text first, then the samples table's ids once each in its order as doubles), the shape and row-group count the
#' manifest states, each row group at most 10,000 rows, and the feature ids (read from the id column alone) present,
#' non-empty and unique.
#' @param dir The omic's folder.
#' @param view The manifest's view entry.
#' @param feature_column The manifest's feature id column.
#' @param sample_ids The samples table's ids, in its order.
#' @param file,row_groups KIT-153: the matrix file and its row-group count (the view's `matrix_file` and `row_groups`, or
#'   its `corrected_matrix_file` and `corrected_row_groups`).
#' @param rows KIT-184: the rows the file must have (the view's features; for the removed factors' file their number).
#' @return The feature ids, in file order.
#' @noRd
checkIntegrationInputMatrix <- function(dir, view, feature_column, sample_ids, file = view$matrix_file,
    row_groups = view$row_groups, rows = view$rows) {
    name <- file
    expected_rows <- rows
    withIntegrationInputParquet(file.path(dir, name), name, function(reader) {
        schema <- integrationInputSchema(reader)
        rows <- as.integer(reader$num_rows)
        if (!identical(rows, as.integer(expected_rows)) || !identical(length(schema$names) - 1L, as.integer(view$columns))) {
            abortIntegrationInput(sprintf("the matrix %s of view %s is not %s x %s", name, view$name, expected_rows,
                view$columns))
        }
        if (!identical(schema$names[[1L]], feature_column) || !identical(schema$types[[1L]], "string")) {
            abortIntegrationInput(sprintf("the matrix of view %s must start with the text column %s", view$name, feature_column))
        }
        if (!identical(schema$names[-1L], sample_ids)) {
            abortIntegrationInput(sprintf("the samples of view %s must be the samples table's, once each, in its order",
                view$name))
        }
        if (!all(schema$types[-1L] == INTEGRATION_INPUT_PARQUET$column_types$runs)) {
            abortIntegrationInput(sprintf("the sample columns of view %s must be numeric (double)", view$name))
        }
        groups <- as.integer(reader$num_row_groups)
        if (!identical(groups, as.integer(row_groups))) {
            abortIntegrationInput(sprintf("the matrix %s of view %s does not have the %s row groups its manifest states",
                name, view$name, row_groups))
        }
        ids <- character(rows)
        at <- 0L
        for (group in seq_len(groups)) {
            part <- integrationInputFrame(reader$ReadRowGroup(group - 1L, 0L))[[1L]]
            if (length(part) > INTEGRATION_INPUT_PARQUET$row_group_rows) {
                abortIntegrationInput(sprintf("the matrix of view %s has a row group of more than %d rows", view$name,
                    INTEGRATION_INPUT_PARQUET$row_group_rows))
            }
            ids[at + seq_along(part)] <- part
            at <- at + length(part)
        }
        if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
            abortIntegrationInput(sprintf("the feature ids of view %s must be present, non-empty and unique", view$name))
        }
        ids
    })
}

#' Write one omic's integration input (KIT-061; KIT-147: Parquet, contract v2): the tables under their fixed names,
#' then the manifest with every file's sha256 and the Parquet layout, in a staging folder that replaces `dir` once
#' complete.
#' @param dir The omic's folder, `<project>/integration/<omic>/`.
#' @param manifest The omic's fields: `omic_type`, `workflow`, `matrix_role` (KIT-153: `normalised`, the role of each
#'   data set's matrix file), `scale`, `normalisation` (the normalisation and the state or product that held the
#'   uncorrected matrix), `correction` (the RUV that produced the corrected matrices; NULL when there are none),
#'   `source_state`, `annotation`, `sample_id_column`, `group_column`, `technical_replicate_column`.
#' @param samples The design of the matrix's samples.
#' @param matrices Named list by view of features x samples tables whose first column is `feature_id`: the normalised
#'   matrices before RUV (KIT-153).
#' @param features Named list by view of feature tables (`feature_id`, `display_name`, annotation columns).
#' @param now The time of writing.
#' @param renameFn The folder rename (`file.rename`; a test injects a failure).
#' @param corrected KIT-153: named list by view of the RUV-corrected matrices of the views the omic corrected (the same
#'   features in the same order and the same samples as the view's matrix); empty when the omic applied no RUV.
#' @return The manifest's path, invisibly.
#' @noRd
writeIntegrationInput <- function(dir, manifest, samples, matrices, features, now = Sys.time(), renameFn = file.rename,
    corrected = list()) {
    views <- names(matrices)
    if (!length(views) || any(!grepl("^[A-Za-z0-9_-]+$", views))) {
        abortIntegrationInput("view names must be letters, digits, underscores or hyphens")
    }
    # KIT-153: each data set's matrix file is the normalised matrix before RUV; a corrected matrix is written beside it.
    if (!identical(manifest$matrix_role, INTEGRATION_INPUT_MATRIX_STAGE)) {
        abortIntegrationInput(paste("the matrix role must be normalised: each data set's matrix is the normalised matrix",
            "before RUV, and a corrected matrix is given beside it"))
    }
    if (length(corrected) && !is.list(manifest$correction)) {
        abortIntegrationInput("a corrected matrix needs the correction that produced it")
    }
    manifest$schema <- INTEGRATION_INPUT_SCHEMA
    manifest$feature_id_column <- "feature_id"
    manifest$display_name_column <- "display_name"
    manifest$views <- lapply(views, function(view) {
        files <- integrationInputViewFiles(view)
        list(name = view, matrix_file = files$matrix, features_file = files$features,
            rows = nrow(matrices[[view]]), columns = ncol(matrices[[view]]) - 1L, matrix_stage = INTEGRATION_INPUT_MATRIX_STAGE)
    })
    sample_ids <- as.character(samples[[manifest$sample_id_column]])
    validateIntegrationInputTables(manifest, samples, matrices, features, corrected)
    parent <- dirname(dir)
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    staging <- tempfile(paste0(".", basename(dir), "-"), tmpdir = parent)
    if (!dir.create(staging)) abortIntegrationInput(sprintf("the staging folder %s could not be created", staging))
    on.exit(unlink(staging, recursive = TRUE), add = TRUE)
    write <- function(table, name) {
        writeIntegrationInputParquet(table, file.path(staging, name))
        name
    }
    written <- write(integrationInputTextTable(samples), INTEGRATION_INPUT_SAMPLES)
    for (i in seq_along(views)) {
        files <- integrationInputViewFiles(views[[i]])
        # The matrix's columns in the samples table's order.
        written <- c(written, write(integrationInputMatrixTable(matrices[[views[[i]]]], sample_ids), files$matrix),
            write(integrationInputTextTable(features[[views[[i]]]]), files$features))
        manifest$views[[i]]$row_groups <- withIntegrationInputParquet(file.path(staging, files$matrix), files$matrix,
            function(reader) as.integer(reader$num_row_groups))
        if (!is.null(corrected[[views[[i]]]])) {
            written <- c(written, write(integrationInputMatrixTable(corrected[[views[[i]]]], sample_ids), files$corrected))
            manifest$views[[i]]$corrected_matrix_file <- files$corrected
            manifest$views[[i]]$corrected_row_groups <- withIntegrationInputParquet(file.path(staging, files$corrected),
                files$corrected, function(reader) as.integer(reader$num_row_groups))
            # KIT-184: the unwanted factors the correction removed, from the two matrices themselves.
            removed <- integrationInputRemovedFactors(matrices[[views[[i]]]], corrected[[views[[i]]]],
                integrationInputDeclaredK(manifest$correction, views[[i]]), sample_ids)
            manifest$views[[i]]$removed_factors <- removed$record
            if (is.null(removed$table)) {
                cli::cli_warn(sprintf(paste("the unwanted factors removed from data set %s could not be determined (%s);",
                    "the integration input is written without them"), views[[i]], removed$record$reason),
                    class = "integration_input_removed_factors_undetermined")
            } else {
                written <- c(written, write(removed$table, files$removed))
                manifest$views[[i]]$removed_factors_file <- files$removed
            }
        }
    }
    manifest$samples_file <- INTEGRATION_INPUT_SAMPLES
    manifest$parquet <- INTEGRATION_INPUT_PARQUET
    manifest$files <- stats::setNames(lapply(written, function(name) {
        digest::digest(file = file.path(staging, name), algo = "sha256")
    }), written)
    manifest$software <- list(package = "MultiScholaR", version = as.character(utils::packageVersion("MultiScholaR")))
    manifest$written_at <- format(now, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
    jsonlite::write_json(manifest, file.path(staging, INTEGRATION_INPUT_MANIFEST), auto_unbox = TRUE, pretty = TRUE,
        digits = NA, null = "null", na = "null")
    # K061-2: the previous input is moved aside, then the staging folder renamed into place; on a failed publish the
    # previous input is renamed back. Only a published input removes the one it replaced.
    previous <- NULL
    if (file.exists(dir)) {
        previous <- tempfile(paste0(".", basename(dir), "-previous-"), tmpdir = parent)
        if (!isTRUE(renameFn(dir, previous))) {
            abortIntegrationInput(sprintf("the previous integration input at %s could not be moved aside; it is kept", dir))
        }
    }
    if (!isTRUE(renameFn(staging, dir))) {
        if (!is.null(previous) && !isTRUE(renameFn(previous, dir))) {
            abortIntegrationInput(sprintf("the integration input could not be published at %s, and the previous input could not be restored from %s",
                dir, previous))
        }
        # KIT-228 (Re-check 1 R01): a previous input is kept only when there was one; a Save removes it at its start.
        abortIntegrationInput(if (is.null(previous)) sprintf("the integration input could not be published at %s", dir) else
            sprintf("the integration input could not be published at %s; the previous input is kept", dir))
    }
    if (!is.null(previous) && (unlink(previous, recursive = TRUE) != 0L || file.exists(previous))) {
        cli::cli_warn(sprintf("the replaced integration input could not be removed from %s", previous))
    }
    invisible(file.path(dir, INTEGRATION_INPUT_MANIFEST))
}

#' Remove an omic's integration input (K061-3): when a save writes none, an earlier save's folder must not validate.
#' @param dir The omic's folder, `<project>/integration/<omic>/`.
#' @return TRUE invisibly, or an `integration_input_invalid` error when the folder cannot be removed.
#' @noRd
removeIntegrationInput <- function(dir) {
    if (file.exists(dir) && (unlink(dir, recursive = TRUE) != 0L || file.exists(dir))) {
        abortIntegrationInput(sprintf("the earlier integration input at %s could not be removed", dir))
    }
    invisible(TRUE)
}

#' Read and validate one omic's integration input (KIT-061; KIT-147: contract v2, no matrix loaded): the manifest,
#' then every file it lists against its sha256, then the samples table, then per view the matrix file's schema, shape
#' and row groups with its feature ids (the id column alone) and the features file against those ids. A v1 folder is
#' refused (`integration_input_outdated`): it was written by a development build and its omic's Summary must be saved
#' again.
#' @param dir The omic's folder, `<project>/integration/<omic>/`.
#' @return `list(manifest, samples, views)`; `views` is named by view, each `list(features, rows, columns,
#'   row_groups)` (the matrix is read through `openIntegrationInputView()`). An input that fails any check is refused
#'   with an `integration_input_invalid` error.
#' @noRd
readIntegrationInput <- function(dir) {
    path <- file.path(dir, INTEGRATION_INPUT_MANIFEST)
    if (!file.exists(path)) abortIntegrationInput(sprintf("no integration input manifest at %s", path))
    manifest <- tryCatch(jsonlite::read_json(path, simplifyVector = FALSE), error = function(error) {
        abortIntegrationInput(sprintf("the integration input manifest cannot be read: %s", conditionMessage(error)))
    })
    if (identical(manifest$schema, INTEGRATION_INPUT_SCHEMA_V1)) {
        abortIntegrationInput(sprintf(paste("the integration input in %s was saved by an earlier development version of",
            "MultiScholaR and cannot be read; save the %s Summary again"), dir, basename(dir)),
            class = "integration_input_outdated")
    }
    if (!identical(manifest$schema, INTEGRATION_INPUT_SCHEMA)) {
        abortIntegrationInput(sprintf("the manifest's schema must be %s", INTEGRATION_INPUT_SCHEMA))
    }
    required <- c("omic_type", "workflow", "matrix_role", "source_state", "sample_id_column", "group_column",
        "feature_id_column", "display_name_column", "views", "samples_file", "parquet", "files", "software", "written_at")
    missing <- setdiff(required, names(manifest))
    if (length(missing)) abortIntegrationInput(sprintf("the manifest lacks %s", paste(missing, collapse = ", ")))
    # KIT-153: every data set names its matrix's stage. The earlier build named none: it kept one matrix per data set,
    # RUV-corrected where the omic applied RUV, which cannot be read as the normalised matrix before RUV.
    stages <- vapply(if (is.list(manifest$views)) manifest$views else list(), function(view) {
        stage <- if (is.list(view)) view$matrix_stage
        if (is.character(stage) && length(stage) == 1L) stage else NA_character_
    }, character(1))
    if (length(stages) && all(is.na(stages))) {
        abortIntegrationInput(sprintf(paste("the integration input in %s was saved by an earlier development version of",
            "MultiScholaR, which kept one matrix per data set (corrected for unwanted variation where the omic applied RUV)",
            "and did not say which; save the %s Summary again"), dir, basename(dir)), class = "integration_input_outdated")
    }
    if (anyNA(stages) || any(stages != INTEGRATION_INPUT_MATRIX_STAGE)) {
        abortIntegrationInput("every data set's matrix must be the normalised matrix before RUV (matrix_stage normalised)")
    }
    # K061-6: the omic is its folder's; the views are named once each, as the writer names them, and their files are
    # the view's fixed names.
    if (!identical(manifest$omic_type, basename(dir))) {
        abortIntegrationInput(sprintf("the manifest's omic must be its folder's, %s", basename(dir)))
    }
    view_names <- vapply(if (is.list(manifest$views)) manifest$views else list(), function(view) {
        if (is.character(view$name) && length(view$name) == 1L && !is.na(view$name)) view$name else NA_character_
    }, character(1))
    if (!length(view_names) || anyNA(view_names) || any(!grepl("^[A-Za-z0-9_-]+$", view_names))) {
        abortIntegrationInput("view names must be letters, digits, underscores or hyphens")
    }
    if (anyDuplicated(view_names)) abortIntegrationInput("the view names must be unique")
    for (view in manifest$views) {
        fixed <- integrationInputViewFiles(view$name)
        if (!identical(view$matrix_file, fixed$matrix) || !identical(view$features_file, fixed$features)) {
            abortIntegrationInput(sprintf("the files of view %s must be %s and %s", view$name, fixed$matrix, fixed$features))
        }
        if (!is.null(view$corrected_matrix_file) && !identical(view$corrected_matrix_file, fixed$corrected)) {
            abortIntegrationInput(sprintf("the corrected matrix of view %s must be %s", view$name, fixed$corrected))
        }
        # KIT-184: the removed factors' file has its fixed name and belongs to a data set with a corrected matrix.
        if (!is.null(view$removed_factors_file) && (!identical(view$removed_factors_file, fixed$removed) ||
            is.null(view$corrected_matrix_file))) {
            abortIntegrationInput(sprintf("the removed factors of view %s must be %s, beside its corrected matrix", view$name,
                fixed$removed))
        }
    }
    corrected_views <- Filter(function(view) !is.null(view$corrected_matrix_file), manifest$views)
    if (!identical(manifest$samples_file, INTEGRATION_INPUT_SAMPLES)) {
        abortIntegrationInput(sprintf("the samples file must be %s", INTEGRATION_INPUT_SAMPLES))
    }
    if (!identical(manifest$parquet, INTEGRATION_INPUT_PARQUET)) {
        abortIntegrationInput(paste("the manifest's Parquet layout must be the contract's: format 2.6, snappy, row groups",
            "of 10,000 rows, text ids and annotation, double values"))
    }
    if (!identical(manifest$matrix_role, INTEGRATION_INPUT_MATRIX_STAGE)) {
        abortIntegrationInput("the matrix role must be normalised: the role of each data set's matrix file")
    }
    if (length(corrected_views) && !is.list(manifest$correction)) {
        abortIntegrationInput("a corrected matrix needs the correction that produced it")
    }
    files <- manifest$files
    if (!is.list(files) || !length(files) || any(!grepl("^[A-Za-z0-9_.-]+$", names(files))) ||
        any(names(files) %in% c(".", ".."))) {
        abortIntegrationInput("the manifest's files must be plain file names with their sha256")
    }
    for (name in names(files)) {
        file <- file.path(dir, name)
        if (!file.exists(file)) abortIntegrationInput(sprintf("the file %s is missing", name))
        if (!identical(digest::digest(file = file, algo = "sha256"), files[[name]])) {
            abortIntegrationInput(sprintf("the file %s does not match its sha256", name))
        }
    }
    listed <- c(manifest$samples_file, unlist(lapply(manifest$views, function(view) c(view$matrix_file, view$features_file,
        view$corrected_matrix_file, view$removed_factors_file))))
    if (!all(listed %in% names(files))) abortIntegrationInput("every table must be listed with its sha256")
    sample_column <- manifest$sample_id_column
    samples <- readIntegrationInputTextTable(file.path(dir, manifest$samples_file), manifest$samples_file)
    if (!all(c(sample_column, manifest$group_column) %in% names(samples))) {
        abortIntegrationInput(sprintf("the samples table needs the columns %s and %s", sample_column, manifest$group_column))
    }
    sample_ids <- samples[[sample_column]]
    if (!length(sample_ids) || anyNA(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
        abortIntegrationInput("the sample ids must be present, non-empty and unique")
    }
    feature_column <- manifest$feature_id_column
    display_column <- manifest$display_name_column
    views <- list()
    for (view in manifest$views) {
        ids <- checkIntegrationInputMatrix(dir, view, feature_column, sample_ids)
        annotation <- readIntegrationInputTextTable(file.path(dir, view$features_file), view$features_file)
        if (!all(c(feature_column, display_column) %in% names(annotation)) ||
            !identical(annotation[[feature_column]], ids)) {
            abortIntegrationInput(sprintf("the features of view %s must list the matrix's ids in order with %s",
                view$name, display_column))
        }
        # KIT-153: the corrected matrix, where there is one, has the same features in the same order and the same samples.
        corrected <- !is.null(view$corrected_matrix_file)
        if (corrected && !identical(checkIntegrationInputMatrix(dir, view, feature_column, sample_ids,
            view$corrected_matrix_file, view$corrected_row_groups), ids)) {
            abortIntegrationInput(sprintf("the corrected matrix of view %s must list the matrix's feature ids in order",
                view$name))
        }
        views[[view$name]] <- list(features = annotation, rows = length(ids), columns = length(sample_ids),
            row_groups = as.integer(view$row_groups), corrected = corrected,
            corrected_row_groups = if (corrected) as.integer(view$corrected_row_groups) else NA_integer_)
        # KIT-184: the removed factors, where the data set has their file: checked here, kept with the file's sha256.
        # A data set without the file adds nothing (its entry is the one the build before KIT-184 gave).
        if (!is.null(view$removed_factors_file)) {
            checkIntegrationInputRemovedFactors(dir, view, feature_column, sample_ids)
            views[[view$name]]$removed <- c(view$removed_factors, list(sha256 = files[[view$removed_factors_file]]))
        }
    }
    list(manifest = manifest, samples = samples, views = views)
}

#' Check a data set's removed-factors file against its manifest entry (KIT-184 7.1): the record names a whole number
#' of factors of at least 1; the file has that many rows named `ruv_W1` to `ruv_Wk`, the runs as double columns in the
#' samples table's order, one row group and no missing value.
#' @param dir The omic's folder.
#' @param view The manifest's view entry (with `removed_factors_file` and `removed_factors`).
#' @param feature_column The manifest's feature id column.
#' @param sample_ids The samples table's ids, in its order.
#' @return TRUE invisibly, or an `integration_input_invalid` error.
#' @noRd
checkIntegrationInputRemovedFactors <- function(dir, view, feature_column, sample_ids) {
    k <- if (is.list(view$removed_factors)) view$removed_factors$k
    if (!is.numeric(k) || length(k) != 1L || is.na(k) || k != round(k) || k < 1) {
        abortIntegrationInput(sprintf("the removed factors of view %s must state their number", view$name))
    }
    ids <- checkIntegrationInputMatrix(dir, view, feature_column, sample_ids, view$removed_factors_file, 1L, as.integer(k))
    if (!identical(ids, paste0("ruv_W", seq_len(k)))) {
        abortIntegrationInput(sprintf("the removed factors of view %s must be named ruv_W1 to ruv_W%d", view$name, as.integer(k)))
    }
    values <- readIntegrationInputMatrix(file.path(dir, view$removed_factors_file), view$removed_factors_file,
        c(feature_column, sample_ids))
    if (!all(vapply(values[sample_ids], function(column) all(is.finite(column)), logical(1)))) {
        abortIntegrationInput(sprintf("the removed factors of view %s must have a value in every run", view$name))
    }
    invisible(TRUE)
}

#' The removed factors of one data set of a validated integration input (KIT-184 8.3).
#' @param dir The omic's folder.
#' @param view The data set's name.
#' @param input `readIntegrationInput(dir)`'s value, or NULL to read and validate the folder here.
#' @return A numeric matrix, runs (in the samples table's order) by factors (`ruv_W1` to `ruv_Wk`), with dimnames;
#'   refused (`integration_input_invalid`) when the data set has no removed factors.
#' @noRd
readIntegrationInputRemovedFactors <- function(dir, view, input = NULL) {
    if (is.null(input)) input <- readIntegrationInput(dir)
    entry <- Filter(function(item) identical(item$name, view), input$manifest$views)
    if (length(entry) != 1L || is.null(entry[[1L]]$removed_factors_file)) {
        abortIntegrationInput(sprintf("the integration input in %s has no removed factors for data set %s", dir, view))
    }
    name <- entry[[1L]]$removed_factors_file
    sample_ids <- input$samples[[input$manifest$sample_id_column]]
    table <- readIntegrationInputMatrix(file.path(dir, name), name, c(input$manifest$feature_id_column, sample_ids))
    factors <- t(as.matrix(table[sample_ids]))
    dimnames(factors) <- list(sample_ids, as.character(table[[1L]]))
    factors
}

#' A reader over one data set of a validated integration input (KIT-147): its size, its zero-row prototype and its
#' matrix by row group (each read opens and closes the file; no row group is kept).
#' @param dir The omic's folder.
#' @param view The data set's name.
#' @param input `readIntegrationInput(dir)`'s value, or NULL to read and validate the folder here.
#' @param matrix KIT-153: `uncorrected` (the normalised matrix before RUV, the view's `matrix_file`) or `ruv_corrected`
#'   (its `corrected_matrix_file`; refused when the data set has none).
#' @return `list(view, path, sha256, rows, columns, row_groups, prototype, read_row_group(index), read())`: `columns`
#'   is `feature_id` then the run ids in the samples table's order; `read_row_group` gives the rows of one row group
#'   (from 1), `read` the whole matrix (small inputs and tests).
#' @noRd
openIntegrationInputView <- function(dir, view, input = NULL, matrix = c("uncorrected", "ruv_corrected")) {
    matrix <- match.arg(matrix)
    if (is.null(input)) input <- readIntegrationInput(dir)
    entry <- Filter(function(item) identical(item$name, view), input$manifest$views)
    if (length(entry) != 1L) abortIntegrationInput(sprintf("the integration input in %s has no data set %s", dir, view))
    entry <- entry[[1L]]
    name <- entry$matrix_file
    row_groups <- entry$row_groups
    if (identical(matrix, "ruv_corrected")) {
        if (is.null(entry$corrected_matrix_file)) {
            abortIntegrationInput(sprintf("the integration input in %s has no RUV-corrected matrix for data set %s", dir, view))
        }
        name <- entry$corrected_matrix_file
        row_groups <- entry$corrected_row_groups
    }
    path <- file.path(dir, name)
    columns <- c(input$manifest$feature_id_column, input$samples[[input$manifest$sample_id_column]])
    list(view = view, path = path, sha256 = input$manifest$files[[name]], rows = as.integer(entry$rows),
        columns = columns, row_groups = as.integer(row_groups),
        prototype = integrationInputMatrixPrototype(columns),
        read_row_group = function(index) readIntegrationInputRowGroup(path, name, index, columns),
        read = function() readIntegrationInputMatrix(path, name, columns))
}
