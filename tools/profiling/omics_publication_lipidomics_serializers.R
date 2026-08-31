lipidPublicationOutputMembers <- function(route, profile_id) {
    profile <- lipidPublicationAssayProfiles()[[profile_id]]
    if (is.null(profile)) lipidPublicationAbort("assay profile is unsupported")
    if (!route %in% names(lipidPublicationCapabilities())) {
        lipidPublicationAbort("serializer route is unsupported")
    }
    assays <- unlist(profile$assays, use.names = FALSE)
    extension <- switch(route, lipidsearch = "txt", msdial = "csv", "tsv")
    as.list(paste0(route, "_", tolower(assays), ".", extension))
}

lipidPublicationBlockRanges <- function(feature_count, chunk_rows) {
    starts <- seq.int(1L, feature_count, by = chunk_rows)
    lapply(starts, \(start) {
        seq.int(start, min(feature_count, start + chunk_rows - 1L))
    })
}

lipidPublicationFormatNumeric <- function(value) {
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

lipidPublicationWriteTable <- function(data, path, delimiter, append) {
    numeric <- vapply(data, is.double, logical(1))
    data[numeric] <- lapply(data[numeric], lipidPublicationFormatNumeric)
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

lipidPublicationAnnotation <- function(index, width) {
    prefix <- sprintf("Synthetic lipid %09d ", index)
    padding <- pmax(0L, width - nchar(prefix, type = "bytes"))
    paste0(prefix, vapply(padding, \(count) {
        paste(rep("X", count), collapse = "")
    }, character(1)))
}

lipidPublicationCoordinates <- function(features, parameters) {
    mz <- 150 + features$family_index * 0.01
    rt <- 0.5 + features$family_index * 0.001
    paired <- !is.na(features$isomer_pair_id)
    if (any(paired)) {
        pair_member <- ave(
            seq_len(sum(paired)),
            features$isomer_pair_id[paired],
            FUN = seq_along
        )
        direction <- ifelse(pair_member %% 2L == 0L, 1, -1)
        ppm <- lipidPublicationParameter(
            parameters,
            "isomer_mass_delta_ppm_max"
        )
        rt_delta <- lipidPublicationParameter(
            parameters,
            "isomer_rt_delta_minutes_max"
        )
        mz[paired] <- mz[paired] * (1 + direction * ppm * 0.9 * 1e-6 / 2)
        rt[paired] <- rt[paired] + direction * rt_delta / 2
    }
    list(
        mz = round(mz, digits = 6L),
        retention_time = round(rt, digits = 6L)
    )
}

lipidPublicationSyntheticName <- function(features, feature_index, width) {
    prefix <- paste(
        features$lipid_class,
        features$composition_family_id,
        sprintf("%09d", feature_index)
    )
    padding <- pmax(0L, width - nchar(prefix, type = "bytes"))
    paste0(prefix, vapply(padding, function(count) {
        paste(rep("X", count), collapse = "")
    }, character(1)))
}

lipidPublicationCustomBlock <- function(plan, block, feature_index, assay) {
    selected <- block$features$assay == assay
    features <- block$features[selected, , drop = FALSE]
    indices <- feature_index[selected]
    values <- block$values[selected, , drop = FALSE]
    width <- as.integer(lipidPublicationParameter(
        plan$parameters,
        "annotation_width"
    ))
    coordinates <- lipidPublicationCoordinates(features, plan$parameters)
    data <- data.frame(
        lipid_id = features$feature_id,
        annotation = lipidPublicationSyntheticName(
            features,
            indices,
            width
        ),
        lipid_class = features$lipid_class,
        adduct = features$adduct,
        ion_mode = features$ion_mode,
        composition_family_id = features$composition_family_id,
        isomer_pair_id = features$isomer_pair_id,
        internal_standard = features$internal_standard,
        mz = coordinates$mz,
        retention_time = coordinates$retention_time,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    quantities <- as.data.frame(values, check.names = FALSE)
    names(quantities) <- plan$design$sample_id
    cbind(data, quantities)
}

lipidPublicationMsdialBlock <- function(plan, block, feature_index, assay) {
    selected <- block$features$assay == assay
    indices <- feature_index[selected]
    features <- block$features[selected, , drop = FALSE]
    values <- block$values[selected, , drop = FALSE]
    width <- as.integer(lipidPublicationParameter(
        plan$parameters,
        "annotation_width"
    ))
    coordinates <- lipidPublicationCoordinates(features, plan$parameters)
    data <- data.frame(
        `Alignment ID` = indices,
        `Average Rt(min)` = coordinates$retention_time,
        `Average Mz` = coordinates$mz,
        Name = lipidPublicationSyntheticName(features, indices, width),
        Ontology = features$lipid_class,
        `Adduct type` = features$adduct,
        Comment = paste(
            features$composition_family_id,
            ifelse(
                is.na(features$isomer_pair_id),
                "none",
                features$isomer_pair_id
            ),
            sep = ";"
        ),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    quantities <- as.data.frame(values, check.names = FALSE)
    names(quantities) <- plan$design$sample_id
    cbind(data, quantities)
}

lipidPublicationLipidSearchBlock <- function(
    plan,
    block,
    feature_index,
    assay
) {
    selected <- block$features$assay == assay
    indices <- feature_index[selected]
    features <- block$features[selected, , drop = FALSE]
    values <- block$values[selected, , drop = FALSE]
    width <- as.integer(lipidPublicationParameter(
        plan$parameters,
        "annotation_width"
    ))
    coordinates <- lipidPublicationCoordinates(features, plan$parameters)
    names <- lipidPublicationSyntheticName(features, indices, width)
    names[features$internal_standard] <- paste0(
        "IS_",
        names[features$internal_standard],
        "_d7"
    )
    data <- data.frame(
        Idx = indices,
        LipidName = names,
        LipidClass = features$lipid_class,
        FattyAcid = features$composition_family_id,
        LipidGroup = ifelse(
            is.na(features$isomer_pair_id),
            features$composition_family_id,
            features$isomer_pair_id
        ),
        IonType = features$adduct,
        BaseRt = coordinates$retention_time,
        CalcMz = coordinates$mz,
        Grade = "Synthetic",
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
    quantities <- as.data.frame(values, check.names = FALSE)
    names(quantities) <- plan$design$sample_id
    cbind(data, quantities)
}

lipidPublicationPayloadBinding <- function(paths) {
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

lipidPublicationSerialize <- function(
    plan,
    route,
    profile_id,
    output_root,
    chunk_rows,
    observer = NULL
) {
    if (file.exists(output_root) || dir.exists(output_root)) {
        lipidPublicationAbort("lipidomics payload root already exists")
    }
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    members <- unlist(
        lipidPublicationOutputMembers(route, profile_id),
        use.names = FALSE
    )
    paths <- file.path(output_root, members)
    assays <- unlist(
        lipidPublicationAssayProfiles()[[profile_id]]$assays,
        use.names = FALSE
    )
    written <- stats::setNames(rep(FALSE, length(paths)), members)
    ranges <- lipidPublicationBlockRanges(nrow(plan$features), chunk_rows)
    for (block_index in seq_along(ranges)) {
        indices <- ranges[[block_index]]
        block <- lipidPublicationGenerateBlock(plan, indices)
        for (assay_index in seq_along(assays)) {
            assay <- assays[[assay_index]]
            if (!any(block$features$assay == assay)) next
            payload <- switch(
                route,
                lipidsearch = lipidPublicationLipidSearchBlock(
                    plan,
                    block,
                    indices,
                    assay
                ),
                msdial = lipidPublicationMsdialBlock(
                    plan,
                    block,
                    indices,
                    assay
                ),
                custom = lipidPublicationCustomBlock(
                    plan,
                    block,
                    indices,
                    assay
                )
            )
            delimiter <- if (identical(route, "msdial")) "," else "\t"
            lipidPublicationWriteTable(
                payload,
                paths[[assay_index]],
                delimiter,
                append = written[[assay_index]]
            )
            written[[assay_index]] <- TRUE
        }
        if (is.function(observer)) observer(block, block_index)
        rm(block, payload)
        if (block_index %% 10L == 0L) gc(verbose = FALSE)
    }
    if (!all(written) || !all(file.exists(paths))) {
        lipidPublicationAbort("lipidomics payload members are incomplete")
    }
    lipidPublicationPayloadBinding(paths)
}
