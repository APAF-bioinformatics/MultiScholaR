metabPublicationOutputMembers <- function(route, profile_id) {
    profile <- metabPublicationAssayProfiles()[[profile_id]]
    if (is.null(profile)) metabPublicationAbort("assay profile is unsupported")
    assays <- unlist(profile$assays, use.names = FALSE)
    if (identical(route, "custom")) return(list("custom.tsv"))
    if (!identical(route, "msdial")) {
        metabPublicationAbort("serializer route is unsupported")
    }
    as.list(paste0(tolower(assays), ".csv"))
}

metabPublicationBlockRanges <- function(feature_count, chunk_rows) {
    starts <- seq.int(1L, feature_count, by = chunk_rows)
    lapply(starts, \(start) {
        seq.int(start, min(feature_count, start + chunk_rows - 1L))
    })
}

metabPublicationFormatNumeric <- function(value) {
    output <- rep.int("", length(value))
    finite <- is.finite(value)
    output[finite] <- formatC(
        value[finite],
        format = "f",
        digits = 6L,
        decimal.mark = "."
    )
    output[is.infinite(value) & value > 0] <- "Inf"
    output[is.infinite(value) & value < 0] <- "-Inf"
    output
}

metabPublicationWriteTable <- function(data, path, delimiter, append) {
    numeric <- vapply(data, is.double, logical(1))
    data[numeric] <- lapply(data[numeric], metabPublicationFormatNumeric)
    utils::write.table(
        data,
        file = path,
        append = append,
        quote = FALSE,
        sep = delimiter,
        eol = "\n",
        na = "",
        row.names = FALSE,
        col.names = !append,
        fileEncoding = "UTF-8"
    )
    invisible(path)
}

metabPublicationAnnotation <- function(index, width) {
    prefix <- sprintf("Synthetic metabolite %09d ", index)
    padding <- pmax(0L, width - nchar(prefix, type = "bytes"))
    paste0(prefix, vapply(padding, \(count) {
        paste(rep("X", count), collapse = "")
    }, character(1)))
}

metabPublicationCustomBlock <- function(plan, block, feature_index) {
    width <- as.integer(metabPublicationParameter(
        plan$parameters,
        "annotation_width"
    ))
    data <- data.frame(
        assay = block$features$assay,
        feature_id = block$features$feature_id,
        annotation = metabPublicationAnnotation(feature_index, width),
        internal_standard = block$features$internal_standard,
        correlated_group_id = block$features$correlated_group_id,
        mz = round(75 + feature_index * 0.001, digits = 6L),
        retention_time = round(0.1 + feature_index * 0.0001, digits = 6L),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    quantities <- as.data.frame(block$values, check.names = FALSE)
    names(quantities) <- plan$design$sample_id
    cbind(data, quantities)
}

metabPublicationMsdialBlock <- function(plan, block, feature_index, assay) {
    selected <- block$features$assay == assay
    indices <- feature_index[selected]
    features <- block$features[selected, , drop = FALSE]
    values <- block$values[selected, , drop = FALSE]
    width <- as.integer(metabPublicationParameter(
        plan$parameters,
        "annotation_width"
    ))
    data <- data.frame(
        `Alignment ID` = indices,
        `Average Rt(min)` = round(0.1 + indices * 0.0001, digits = 6L),
        `Average Mz` = round(75 + indices * 0.001, digits = 6L),
        Name = metabPublicationAnnotation(indices, width),
        `Adduct type` = if (identical(assay, "LCMS_Neg")) "[M-H]-" else {
            if (identical(assay, "LCMS_Pos")) "[M+H]+" else ""
        },
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    quantities <- as.data.frame(values, check.names = FALSE)
    names(quantities) <- plan$design$sample_id
    cbind(data, quantities)
}

metabPublicationPayloadBinding <- function(paths) {
    records <- lapply(sort(paths, method = "radix"), \(path) {
        list(
            filename = basename(path),
            sha256 = publicationFileDigest(path),
            size_bytes = as.numeric(file.info(path)$size)
        )
    })
    list(
        members = records,
        member_count = as.integer(length(records)),
        payload_set_sha256 = publicationObjectDigest(records),
        size_bytes = sum(vapply(records, `[[`, numeric(1), "size_bytes"))
    )
}

metabPublicationSerialize <- function(
    plan,
    route,
    profile_id,
    output_root,
    chunk_rows,
    observer = NULL
) {
    if (file.exists(output_root) || dir.exists(output_root)) {
        metabPublicationAbort("metabolomics payload root already exists")
    }
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    members <- unlist(
        metabPublicationOutputMembers(route, profile_id),
        use.names = FALSE
    )
    paths <- file.path(output_root, members)
    assays <- unlist(
        metabPublicationAssayProfiles()[[profile_id]]$assays,
        use.names = FALSE
    )
    written <- stats::setNames(rep(FALSE, length(paths)), members)
    ranges <- metabPublicationBlockRanges(nrow(plan$features), chunk_rows)
    for (block_index in seq_along(ranges)) {
        indices <- ranges[[block_index]]
        block <- metabPublicationGenerateBlock(plan, indices)
        if (identical(route, "custom")) {
            payload <- metabPublicationCustomBlock(plan, block, indices)
            metabPublicationWriteTable(
                payload,
                paths[[1L]],
                "\t",
                append = written[[1L]]
            )
            written[[1L]] <- TRUE
        } else {
            for (assay_index in seq_along(assays)) {
                assay <- assays[[assay_index]]
                if (!any(block$features$assay == assay)) next
                payload <- metabPublicationMsdialBlock(
                    plan,
                    block,
                    indices,
                    assay
                )
                metabPublicationWriteTable(
                    payload,
                    paths[[assay_index]],
                    ",",
                    append = written[[assay_index]]
                )
                written[[assay_index]] <- TRUE
            }
        }
        if (is.function(observer)) observer(block, block_index)
        rm(block, payload)
        if (block_index %% 10L == 0L) gc(verbose = FALSE)
    }
    if (!all(written) || !all(file.exists(paths))) {
        metabPublicationAbort("metabolomics payload members are incomplete")
    }
    metabPublicationPayloadBinding(paths)
}
