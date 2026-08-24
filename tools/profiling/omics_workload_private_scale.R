privateScaleAbort <- function(message, class) {
    condition <- structure(
        list(message = message, call = NULL),
        class = c(class, "omics_private_scale_error", "error", "condition")
    )
    stop(condition)
}

privateScaleRequireSalt <- function(salt) {
    if (!is.character(salt) || length(salt) != 1L || is.na(salt) ||
        nchar(salt, type = "bytes") < 16L) {
        privateScaleAbort(
            "private scale inspection requires a non-empty secret salt",
            "omics_private_scale_unsalted"
        )
    }
    invisible(salt)
}

privateScaleFileState <- function(path) {
    info <- file.info(path)
    if (!file.exists(path) || dir.exists(path) || !isTRUE(info$isdir == FALSE)) {
        privateScaleAbort(
            "private scale source must be one readable regular file",
            "omics_private_scale_unreadable"
        )
    }
    list(
        byte_size = unname(as.numeric(info$size)),
        modified = as.numeric(info$mtime),
        digest = digest::digest(file = path, algo = "sha256", serialize = FALSE)
    )
}

privateScaleInspectRectangularTsv <- function(path, chunk_size = 10000L) {
    connection <- file(path, open = "rt", encoding = "UTF-8")
    on.exit(close(connection), add = TRUE)
    header <- readLines(connection, n = 1L, warn = FALSE)
    if (length(header) != 1L || !nzchar(header)) {
        privateScaleAbort(
            "private scale TSV has no header",
            "omics_private_scale_nonrectangular"
        )
    }
    column_count <- lengths(regmatches(header, gregexpr("\t", header))) + 1L
    row_count <- 0L
    repeat {
        lines <- readLines(connection, n = chunk_size, warn = FALSE)
        if (!length(lines)) break
        observed <- lengths(regmatches(lines, gregexpr("\t", lines))) + 1L
        if (any(observed != column_count)) {
            privateScaleAbort(
                "private scale TSV is not rectangular",
                "omics_private_scale_nonrectangular"
            )
        }
        row_count <- row_count + length(lines)
    }
    if (row_count < 1L || column_count < 2L) {
        privateScaleAbort(
            "private scale TSV has insufficient dimensions",
            "omics_private_scale_nonrectangular"
        )
    }
    list(row_count = as.integer(row_count), column_count = as.integer(column_count))
}

inspectPrivateScaleTsv <- function(path, salt) {
    privateScaleRequireSalt(salt)
    before <- privateScaleFileState(path)
    shape <- privateScaleInspectRectangularTsv(path)
    after <- privateScaleFileState(path)
    if (!identical(before, after)) {
        privateScaleAbort(
            "private scale source changed during inspection",
            "omics_private_scale_changed"
        )
    }
    result <- list(
        row_count = shape$row_count,
        column_count = shape$column_count,
        byte_size = before$byte_size,
        salted_source_fingerprint = digest::digest(
            paste0(salt, ":", before$digest),
            algo = "sha256",
            serialize = FALSE
        )
    )
    if (!identical(
        names(result),
        c(
            "row_count", "column_count", "byte_size",
            "salted_source_fingerprint"
        )
    )) {
        privateScaleAbort(
            "private scale manifest contains unexpected fields",
            "omics_private_scale_unsafe_manifest"
        )
    }
    result
}

mapPrivateScaleToMetabolomics <- function(scale_metadata) {
    required <- c(
        "row_count", "column_count", "byte_size",
        "salted_source_fingerprint"
    )
    if (!is.list(scale_metadata) || !identical(names(scale_metadata), required)) {
        privateScaleAbort(
            "private scale mapping requires the sanitized manifest",
            "omics_private_scale_unsafe_manifest"
        )
    }
    rows <- as.integer(scale_metadata$row_count)
    report_columns <- as.integer(scale_metadata$column_count)
    bytes <- as.numeric(scale_metadata$byte_size)
    if (!is.finite(rows) || rows < 1L || !is.finite(report_columns) ||
        report_columns < 2L || !is.finite(bytes) || bytes < 1) {
        privateScaleAbort(
            "private scale metadata is malformed",
            "omics_private_scale_unsafe_manifest"
        )
    }
    list(
        feature_count = max(10000L, min(25000L, as.integer(round(rows / 4L)))),
        sample_count = 12L,
        assay_mix = list(LCMS_Pos = 1L, LCMS_Neg = 1L, GCMS = 1L),
        source_report_column_count = report_columns,
        source_byte_class = if (bytes >= 64 * 1024^2) "64_to_128_mib" else "under_64_mib",
        mapping_policy = "rows_to_bounded_features_columns_io_only_v1"
    )
}
