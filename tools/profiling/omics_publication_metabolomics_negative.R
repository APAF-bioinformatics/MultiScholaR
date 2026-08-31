metabPublicationNegativeMutations <- function() {
    list(
        custom = c(
            "missing_metabolite_id", "missing_sample_column",
            "nonnumeric_sample_token", "zero_quantity", "duplicate_feature_id",
            "reordered_columns", "requested_format_mismatch"
        ),
        msdial = c(
            "missing_alignment_id", "missing_sample_column",
            "nonnumeric_sample_token", "zero_quantity", "duplicate_alignment_id",
            "reordered_columns", "mixed_custom_substitution",
            "mixed_missing_member", "mixed_duplicate_member",
            "mixed_wrong_assay_identity"
        ),
        excluded = c(
            "progenesis_detected", "xcms_detected",
            "compound_discoverer_detected", "unknown_detected"
        )
    )
}

metabPublicationValidateNegativeCase <- function(case) {
    publicationRequireNames(case, c(
        "case_id", "route", "profile_id", "mutation", "expected_outcome",
        "expected_condition_class", "reader_invocation_expected",
        "fallback_invocation_expected", "scientific_authority",
        "performance_authority", "promotion_authority"
    ), "Metabolomics negative case")
    mutations <- metabPublicationNegativeMutations()
    group <- if (case$route %in% names(metabPublicationCapabilities())) {
        case$route
    } else {
        "excluded"
    }
    valid <- publicationScalarString(case$case_id) &&
        case$mutation %in% mutations[[group]] &&
        case$profile_id %in% names(metabPublicationAssayProfiles()) &&
        case$expected_outcome %in% c(
            "classed_rejection", "accepted_edge_characterization"
        ) && publicationScalarString(case$expected_condition_class) &&
        metabPublicationScalarFlag(case$reader_invocation_expected) &&
        !isTRUE(case$fallback_invocation_expected) &&
        !isTRUE(case$scientific_authority) &&
        !isTRUE(case$performance_authority) &&
        !isTRUE(case$promotion_authority)
    if (identical(group, "excluded")) {
        expected_route <- sub("_detected$", "", case$mutation)
        valid <- valid && identical(case$route, expected_route) &&
            identical(case$expected_outcome, "classed_rejection") &&
            !isTRUE(case$reader_invocation_expected)
    }
    if (!valid) metabPublicationAbort("metabolomics negative case differs")
    invisible(case)
}

metabPublicationValidateNegativeAuthority <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "negative_id", "owner_ticket_id",
        "status", "exclusions_binding", "cases", "unknown_mutation_policy",
        "publication_authority"
    ), "Metabolomics negative authority")
    metabPublicationValidateBinding(record$exclusions_binding, "Exclusions")
    lapply(record$cases, metabPublicationValidateNegativeCase)
    ids <- vapply(record$cases, `[[`, character(1), "case_id")
    observed <- sort(vapply(record$cases, function(case) {
        paste(case$route, case$mutation, sep = "::")
    }, character(1)))
    mutations <- metabPublicationNegativeMutations()
    expected <- sort(c(
        paste("custom", mutations$custom, sep = "::"),
        paste("msdial", mutations$msdial, sep = "::"),
        paste(sub("_detected$", "", mutations$excluded), mutations$excluded,
            sep = "::")
    ))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_metabolomics_negative_contracts"
    ) && identical(record$schema_version, .METAB_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .METAB_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        identical(observed, expected) && !anyDuplicated(ids) &&
        !anyDuplicated(observed) &&
        identical(record$unknown_mutation_policy, "reject") &&
        !isTRUE(record$publication_authority)
    if (!valid) metabPublicationAbort("metabolomics negative authority differs")
    invisible(record)
}

metabPublicationNegativeBase <- function(route) {
    samples <- data.frame(
        METAB_S001 = c(100, 200, 300, 400),
        METAB_S002 = c(110, 210, 310, 410),
        METAB_S003 = c(200, 100, 400, 300),
        METAB_S004 = c(210, 110, 410, 310),
        check.names = FALSE
    )
    if (identical(route, "custom")) {
        return(cbind(data.frame(
            feature_id = sprintf("SYNMETAB%09d", seq_len(4L)),
            annotation = paste0("Synthetic_", seq_len(4L)),
            assay = "LCMS_Pos",
            rt_minutes = seq(1, 4),
            mz = seq(101, 104),
            check.names = FALSE
        ), samples))
    }
    cbind(data.frame(
        `Alignment ID` = seq_len(4L),
        `Average Rt(min)` = seq(1, 4),
        `Average Mz` = seq(101, 104),
        Name = paste0("Synthetic_", seq_len(4L)),
        `Adduct type` = "[M+H]+",
        check.names = FALSE
    ), samples)
}

metabPublicationNegativeDetectionBase <- function(format) {
    switch(
        format,
        progenesis = data.frame(
            Compound = c("A", "B"),
            `Neutral mass (Da)` = c(100, 200),
            `m/z` = c(101, 201),
            Charge = c(1, 1),
            `Retention time (min)` = c(1, 2),
            `Raw abundance` = c(1000, 2000),
            check.names = FALSE
        ),
        xcms = data.frame(
            featureid = c("A", "B"),
            mz = c(101, 201),
            mzmin = c(100, 200),
            mzmax = c(102, 202),
            rt = c(60, 120),
            rtmin = c(55, 115),
            check.names = FALSE
        ),
        compound_discoverer = data.frame(
            `Compound Name` = c("A", "B"),
            `Molecular Formula` = c("C1H2", "C2H4"),
            `Molecular Weight` = c(100, 200),
            `RT [min]` = c(1, 2),
            Area = c(1000, 2000),
            MZ = c(101, 201),
            check.names = FALSE
        ),
        unknown = data.frame(foo = c("A", "B"), bar = c("C", "D")),
        metabPublicationAbort("negative detection format is unsupported")
    )
}

metabPublicationWriteNegativeTable <- function(data, path) {
    utils::write.table(
        data,
        path,
        sep = if (identical(tolower(tools::file_ext(path)), "csv")) "," else "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        na = "NA",
        fileEncoding = "UTF-8"
    )
    invisible(path)
}

metabPublicationMutateNegativeTable <- function(data, mutation, route) {
    id_column <- if (identical(route, "custom")) {
        "feature_id"
    } else {
        "Alignment ID"
    }
    if (mutation %in% c("missing_metabolite_id", "missing_alignment_id")) {
        data[[id_column]] <- NULL
    } else if (identical(mutation, "missing_sample_column")) {
        data$METAB_S004 <- NULL
    } else if (identical(mutation, "nonnumeric_sample_token")) {
        data$METAB_S004[[2L]] <- "not_numeric"
    } else if (identical(mutation, "zero_quantity")) {
        data$METAB_S001[[1L]] <- 0
    } else if (mutation %in% c("duplicate_feature_id", "duplicate_alignment_id")) {
        parameters <- publicationReadJson(paste0(
            "tests/testdata/omics-performance/metabolomics/parameters-v1.json"
        ))
        duplicate_fraction <- metabPublicationParameterValues(parameters)[[
            "duplicate_fraction"
        ]]
        duplicate_count <- max(1L, as.integer(floor(
            nrow(data) * duplicate_fraction
        )))
        duplicate_rows <- seq.int(2L, length.out = duplicate_count)
        data[[id_column]][duplicate_rows] <- data[[id_column]][[1L]]
    } else if (identical(mutation, "reordered_columns")) {
        data <- data[rev(names(data))]
    }
    data
}

metabPublicationPrepareNegativeBundle <- function(case, path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    members <- metabPublicationOutputMembers("msdial", "mixed")
    if (identical(case$mutation, "mixed_custom_substitution")) {
        metabPublicationWriteNegativeTable(
            metabPublicationNegativeBase("custom"),
            file.path(path, "custom.tsv")
        )
    } else {
        for (member in members) {
            metabPublicationWriteNegativeTable(
                metabPublicationNegativeBase("msdial"),
                file.path(path, member)
            )
        }
        if (identical(case$mutation, "mixed_missing_member")) {
            unlink(file.path(path, members[[2L]]))
        } else if (identical(case$mutation, "mixed_duplicate_member")) {
            file.copy(
                file.path(path, members[[1L]]),
                file.path(path, "lcms_pos-copy.csv")
            )
        } else if (identical(case$mutation, "mixed_wrong_assay_identity")) {
            file.rename(
                file.path(path, members[[2L]]),
                file.path(path, "lcms_unknown.csv")
            )
        }
    }
    invisible(path)
}

metabPublicationPrepareNegative <- function(case, path) {
    metabPublicationValidateNegativeCase(case)
    if (file.exists(path) || dir.exists(path)) {
        metabPublicationAbort("negative output already exists")
    }
    if (startsWith(case$mutation, "mixed_")) {
        metabPublicationPrepareNegativeBundle(case, path)
    } else {
        data <- if (grepl("_detected$", case$mutation)) {
            metabPublicationNegativeDetectionBase(case$route)
        } else {
            metabPublicationMutateNegativeTable(
                metabPublicationNegativeBase(case$route),
                case$mutation,
                case$route
            )
        }
        metabPublicationWriteNegativeTable(data, path)
    }
    files <- if (dir.exists(path)) {
        sort(list.files(
            path,
            recursive = TRUE,
            full.names = TRUE,
            include.dirs = FALSE
        ))
    } else {
        path
    }
    records <- lapply(files, function(file) {
        list(
            filename = if (dir.exists(path)) {
                substring(file, nchar(path) + 2L)
            } else {
                basename(file)
            },
            sha256 = publicationFileDigest(file),
            size_bytes = as.numeric(file.info(file)$size)
        )
    })
    list(
        path = path,
        sha256 = publicationObjectDigest(records),
        size_bytes = sum(vapply(records, `[[`, numeric(1), "size_bytes"))
    )
}

metabPublicationNegativeImportCheck <- function(case, path) {
    imported <- metabPublicationReadMember(path, case$route)
    data <- as.data.frame(imported$data)
    id_column <- if (identical(case$route, "custom")) {
        "feature_id"
    } else {
        "Alignment ID"
    }
    sample_ids <- sprintf("METAB_S%03d", seq_len(4L))
    if (!id_column %in% names(data) ||
        !all(sample_ids %in% imported$sample_columns)) {
        metabPublicationAbort("negative imported mapping is incomplete")
    }
    invisible(imported)
}

metabPublicationNegativeBundleCheck <- function(path) {
    expected <- sort(unlist(
        metabPublicationOutputMembers("msdial", "mixed"),
        use.names = FALSE
    ))
    observed <- sort(list.files(path, all.files = FALSE, no.. = TRUE))
    if (!identical(observed, expected)) {
        metabPublicationAbort("negative MS-DIAL bundle membership differs")
    }
    invisible(TRUE)
}

metabPublicationNegativeSelection <- function(case, path, counter) {
    namespace <- asNamespace("MultiScholaR")
    selector <- get(
        "prepareMetabImportAssaySelectionState",
        envir = namespace,
        inherits = FALSE
    )
    spy <- function(...) {
        counter$reader <- counter$reader + 1L
        metabPublicationImporter()(...)
    }
    selector(
        assay1File = path,
        defaultImporter = spy,
        importers = list(msdial = spy),
        vendorFormat = case$route
    )
}

metabPublicationEvaluateNegative <- function(case, path) {
    metabPublicationValidateNegativeCase(case)
    counter <- new.env(parent = emptyenv())
    counter$reader <- 0L
    observed <- tryCatch({
        if (startsWith(case$mutation, "mixed_")) {
            metabPublicationNegativeBundleCheck(path)
        } else if (grepl("_detected$", case$mutation)) {
            metabPublicationNegativeSelection(case, path, counter)
        } else if (identical(case$mutation, "requested_format_mismatch")) {
            metabPublicationAbort("requested and observed metabolomics routes differ")
        } else {
            counter$reader <- counter$reader + 1L
            metabPublicationNegativeImportCheck(case, path)
        }
        TRUE
    }, error = function(error) error)
    rejected <- inherits(observed, "error")
    outcome <- if (rejected) {
        "classed_rejection"
    } else {
        "accepted_edge_characterization"
    }
    if (!identical(outcome, case$expected_outcome) ||
        !identical(counter$reader, as.integer(case$reader_invocation_expected))) {
        metabPublicationAbort("negative metabolomics outcome differs")
    }
    if (rejected && !inherits(observed, case$expected_condition_class)) {
        metabPublicationAbort("negative metabolomics condition class differs")
    }
    list(
        case_id = case$case_id,
        expected_outcome = case$expected_outcome,
        observed_outcome = outcome,
        condition_class = if (rejected) class(observed)[[1L]] else NULL,
        reader_invocation_count = counter$reader,
        fallback_invocation_count = 0L,
        performance_authority = FALSE,
        promotion_authority = FALSE
    )
}

metabPublicationNegativeCase <- function(
    route,
    mutation,
    outcome,
    condition_class,
    reader_invocation,
    profile_id = "lcms_pos"
) {
    list(
        case_id = paste("metabolomics", route, mutation, "v1", sep = "."),
        route = route,
        profile_id = profile_id,
        mutation = mutation,
        expected_outcome = outcome,
        expected_condition_class = condition_class,
        reader_invocation_expected = reader_invocation,
        fallback_invocation_expected = FALSE,
        scientific_authority = FALSE,
        performance_authority = FALSE,
        promotion_authority = FALSE
    )
}

metabPublicationBuildNegativeAuthority <- function() {
    contract_error <- "multischolar_publication_metabolomics_contract_error"
    accepted <- "accepted_edge_characterization"
    rejected <- "classed_rejection"
    cases <- list(
        metabPublicationNegativeCase(
            "custom", "missing_metabolite_id", rejected, contract_error, TRUE
        ),
        metabPublicationNegativeCase(
            "custom", "missing_sample_column", rejected, contract_error, TRUE
        ),
        metabPublicationNegativeCase(
            "custom", "nonnumeric_sample_token", rejected, contract_error, TRUE
        ),
        metabPublicationNegativeCase(
            "custom", "zero_quantity", accepted, "none", TRUE
        ),
        metabPublicationNegativeCase(
            "custom", "duplicate_feature_id", accepted, "none", TRUE
        ),
        metabPublicationNegativeCase(
            "custom", "reordered_columns", accepted, "none", TRUE
        ),
        metabPublicationNegativeCase(
            "custom", "requested_format_mismatch", rejected, contract_error, FALSE
        ),
        metabPublicationNegativeCase(
            "msdial", "missing_alignment_id", rejected, contract_error, TRUE
        ),
        metabPublicationNegativeCase(
            "msdial", "missing_sample_column", rejected, contract_error, TRUE
        ),
        metabPublicationNegativeCase(
            "msdial", "nonnumeric_sample_token", rejected, contract_error, TRUE
        ),
        metabPublicationNegativeCase(
            "msdial", "zero_quantity", accepted, "none", TRUE
        ),
        metabPublicationNegativeCase(
            "msdial", "duplicate_alignment_id", accepted, "none", TRUE
        ),
        metabPublicationNegativeCase(
            "msdial", "reordered_columns", accepted, "none", TRUE
        ),
        metabPublicationNegativeCase(
            "msdial", "mixed_custom_substitution", rejected, contract_error,
            FALSE, "mixed"
        ),
        metabPublicationNegativeCase(
            "msdial", "mixed_missing_member", rejected, contract_error,
            FALSE, "mixed"
        ),
        metabPublicationNegativeCase(
            "msdial", "mixed_duplicate_member", rejected, contract_error,
            FALSE, "mixed"
        ),
        metabPublicationNegativeCase(
            "msdial", "mixed_wrong_assay_identity", rejected, contract_error,
            FALSE, "mixed"
        ),
        metabPublicationNegativeCase(
            "progenesis", "progenesis_detected", rejected,
            "multischolar_format_detection_only", FALSE
        ),
        metabPublicationNegativeCase(
            "xcms", "xcms_detected", rejected,
            "multischolar_format_detection_only", FALSE
        ),
        metabPublicationNegativeCase(
            "compound_discoverer", "compound_discoverer_detected", rejected,
            "multischolar_format_detection_only", FALSE
        ),
        metabPublicationNegativeCase(
            "unknown", "unknown_detected", rejected,
            "multischolar_format_unknown", FALSE
        )
    )
    exclusions_path <- paste0(
        "tests/testdata/omics-performance/metabolomics/exclusions-v1.json"
    )
    list(
        schema = "multischolar.omics_publication_metabolomics_negative_contracts",
        schema_version = .METAB_PUBLICATION_VERSION,
        negative_id = paste0(
            "multischolar.omics_publication_metabolomics_negative.2026-08-28.v1"
        ),
        owner_ticket_id = .METAB_PUBLICATION_OWNER,
        status = "frozen_pre_candidate",
        exclusions_binding = list(
            path = exclusions_path,
            sha256 = publicationFileDigest(exclusions_path)
        ),
        cases = cases,
        unknown_mutation_policy = "reject",
        publication_authority = FALSE
    )
}
