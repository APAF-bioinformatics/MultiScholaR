lipidPublicationNegativeBaseMutations <- function(route) {
    id_mutation <- switch(
        route,
        lipidsearch = "missing_lipid_name",
        msdial = "missing_alignment_id",
        custom = "missing_lipid_id"
    )
    c(
        id_mutation,
        "missing_sample_column",
        "nonnumeric_sample_token",
        "zero_quantity",
        "duplicate_identifier",
        "reordered_columns",
        "requested_format_mismatch"
    )
}

lipidPublicationNegativeBundleMutations <- function() {
    c(
        "mixed_single_table_substitution",
        "mixed_missing_member",
        "mixed_duplicate_member",
        "mixed_wrong_assay_identity",
        "mixed_three_member_bundle"
    )
}

lipidPublicationNegativeExcludedMutations <- function() {
    c(
        "progenesis_detected",
        "xcms_detected",
        "compound_discoverer_detected",
        "unknown_detected"
    )
}

lipidPublicationValidateNegativeCase <- function(case) {
    publicationRequireNames(case, c(
        "case_id", "route", "profile_id", "mutation", "expected_outcome",
        "expected_condition_class", "reader_invocation_expected",
        "fallback_invocation_expected", "scientific_authority",
        "performance_authority", "promotion_authority"
    ), "Lipidomics negative case")
    supported <- case$route %in% names(lipidPublicationCapabilities())
    allowed <- if (supported) {
        c(
            lipidPublicationNegativeBaseMutations(case$route),
            lipidPublicationNegativeBundleMutations()
        )
    } else {
        lipidPublicationNegativeExcludedMutations()
    }
    valid <- publicationScalarString(case$case_id) &&
        case$mutation %in% allowed &&
        case$profile_id %in% names(lipidPublicationAssayProfiles()) &&
        case$expected_outcome %in% c(
            "classed_rejection", "accepted_edge_characterization"
        ) && publicationScalarString(case$expected_condition_class) &&
        lipidPublicationScalarFlag(case$reader_invocation_expected) &&
        !isTRUE(case$fallback_invocation_expected) &&
        !isTRUE(case$scientific_authority) &&
        !isTRUE(case$performance_authority) &&
        !isTRUE(case$promotion_authority)
    if (!supported) {
        expected_route <- sub("_detected$", "", case$mutation)
        valid <- valid && identical(case$route, expected_route) &&
            identical(case$expected_outcome, "classed_rejection") &&
            !isTRUE(case$reader_invocation_expected)
    }
    if (!valid) lipidPublicationAbort("lipidomics negative case differs")
    invisible(case)
}

lipidPublicationValidateNegativeAuthority <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "negative_id", "owner_ticket_id",
        "status", "exclusions_binding", "cases", "unknown_mutation_policy",
        "publication_authority"
    ), "Lipidomics negative authority")
    lipidPublicationValidateBinding(record$exclusions_binding, "Exclusions")
    lapply(record$cases, lipidPublicationValidateNegativeCase)
    ids <- vapply(record$cases, `[[`, character(1), "case_id")
    observed <- sort(vapply(record$cases, function(case) {
        paste(case$route, case$mutation, sep = "::")
    }, character(1)))
    expected <- character()
    for (route in names(lipidPublicationCapabilities())) {
        expected <- c(
            expected,
            paste(
                route,
                c(
                    lipidPublicationNegativeBaseMutations(route),
                    lipidPublicationNegativeBundleMutations()
                ),
                sep = "::"
            )
        )
    }
    excluded <- lipidPublicationNegativeExcludedMutations()
    expected <- sort(c(
        expected,
        paste(sub("_detected$", "", excluded), excluded, sep = "::")
    ))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_lipidomics_negative_contracts"
    ) && identical(record$schema_version, .LIPID_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .LIPID_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        identical(observed, expected) && !anyDuplicated(ids) &&
        !anyDuplicated(observed) &&
        identical(record$unknown_mutation_policy, "reject") &&
        !isTRUE(record$publication_authority)
    if (!valid) lipidPublicationAbort("lipidomics negative authority differs")
    invisible(record)
}

lipidPublicationNegativeSamples <- function() {
    data.frame(
        LIPID_S001 = c(100, 200, 300, 400),
        LIPID_S002 = c(110, 210, 310, 410),
        LIPID_S003 = c(105, 205, 305, 405),
        LIPID_S004 = c(200, 100, 400, 300),
        LIPID_S005 = c(210, 110, 410, 310),
        LIPID_S006 = c(205, 105, 405, 305),
        check.names = FALSE
    )
}

lipidPublicationNegativeBase <- function(route) {
    samples <- lipidPublicationNegativeSamples()
    if (identical(route, "custom")) {
        return(cbind(data.frame(
            lipid_id = sprintf("SYNLIPID%09d", seq_len(4L)),
            annotation = paste0("Synthetic_lipid_", seq_len(4L)),
            lipid_class = "SYN_PC",
            adduct = "[M+H]+",
            ion_mode = "positive",
            composition_family_id = sprintf("SYNCOMP%08d", seq_len(4L)),
            isomer_pair_id = NA_character_,
            internal_standard = FALSE,
            mz = seq(101, 104),
            retention_time = seq(1, 4),
            check.names = FALSE
        ), samples))
    }
    if (identical(route, "msdial")) {
        return(cbind(data.frame(
            `Alignment ID` = seq_len(4L),
            `Average Rt(min)` = seq(1, 4),
            `Average Mz` = seq(101, 104),
            Name = paste0("Synthetic_lipid_", seq_len(4L)),
            Ontology = "SYN_PC",
            `Adduct type` = "[M+H]+",
            Comment = paste0("SYNCOMP0000000", seq_len(4L), ";none"),
            check.names = FALSE
        ), samples))
    }
    cbind(data.frame(
        Idx = seq_len(4L),
        LipidName = paste0("SYN_PC_", seq_len(4L)),
        LipidClass = "SYN_PC",
        FattyAcid = paste0("SYNCOMP0000000", seq_len(4L)),
        LipidGroup = paste0("SYNCOMP0000000", seq_len(4L)),
        IonType = "[M+H]+",
        BaseRt = seq(1, 4),
        CalcMz = seq(101, 104),
        Grade = "Synthetic",
        check.names = FALSE
    ), samples)
}

lipidPublicationNegativeDetectionBase <- function(format) {
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
        lipidPublicationAbort("negative detection format is unsupported")
    )
}

lipidPublicationWriteNegativeTable <- function(data, path) {
    delimiter <- if (identical(tolower(tools::file_ext(path)), "csv")) {
        ","
    } else {
        "\t"
    }
    lipidPublicationWriteTable(data, path, delimiter, append = FALSE)
}

lipidPublicationNegativeIdColumn <- function(route) {
    switch(
        route,
        lipidsearch = "LipidName",
        msdial = "Alignment ID",
        custom = "lipid_id"
    )
}

lipidPublicationMutateNegativeTable <- function(data, mutation, route) {
    id_column <- lipidPublicationNegativeIdColumn(route)
    if (startsWith(mutation, "missing_") && grepl(
        "(lipid_name|alignment_id|lipid_id)$",
        mutation
    )) {
        data[[id_column]] <- NULL
    } else if (identical(mutation, "missing_sample_column")) {
        data$LIPID_S006 <- NULL
    } else if (identical(mutation, "nonnumeric_sample_token")) {
        data$LIPID_S006[[2L]] <- "not_numeric"
    } else if (identical(mutation, "zero_quantity")) {
        data$LIPID_S001[[1L]] <- 0
    } else if (identical(mutation, "duplicate_identifier")) {
        parameters <- publicationReadJson(paste0(
            "tests/testdata/omics-performance/lipidomics/parameters-v1.json"
        ))
        fraction <- lipidPublicationParameterValues(parameters)[[
            "duplicate_fraction"
        ]]
        count <- max(1L, as.integer(floor(nrow(data) * fraction)))
        rows <- seq.int(2L, length.out = count)
        data[[id_column]][rows] <- data[[id_column]][[1L]]
    } else if (identical(mutation, "reordered_columns")) {
        data <- data[rev(names(data))]
    }
    data
}

lipidPublicationPrepareNegativeBundle <- function(case, path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    members <- unlist(lipidPublicationOutputMembers(case$route, "mixed_lc"))
    if (identical(case$mutation, "mixed_single_table_substitution")) {
        lipidPublicationWriteNegativeTable(
            lipidPublicationNegativeBase("custom"),
            file.path(path, "combined.tsv")
        )
        return(invisible(path))
    }
    for (member in members) {
        lipidPublicationWriteNegativeTable(
            lipidPublicationNegativeBase(case$route),
            file.path(path, member)
        )
    }
    if (identical(case$mutation, "mixed_missing_member")) {
        unlink(file.path(path, members[[2L]]))
    } else if (identical(case$mutation, "mixed_duplicate_member")) {
        file.copy(
            file.path(path, members[[1L]]),
            file.path(path, paste0("copy-", members[[1L]]))
        )
    } else if (identical(case$mutation, "mixed_wrong_assay_identity")) {
        file.rename(
            file.path(path, members[[2L]]),
            file.path(path, paste0(case$route, "_unknown.tsv"))
        )
    } else if (identical(case$mutation, "mixed_three_member_bundle")) {
        lipidPublicationWriteNegativeTable(
            lipidPublicationNegativeBase(case$route),
            file.path(path, paste0(case$route, "_gcms.tsv"))
        )
    }
    invisible(path)
}

lipidPublicationPrepareNegative <- function(case, path) {
    lipidPublicationValidateNegativeCase(case)
    if (file.exists(path) || dir.exists(path)) {
        lipidPublicationAbort("negative output already exists")
    }
    if (startsWith(case$mutation, "mixed_")) {
        lipidPublicationPrepareNegativeBundle(case, path)
    } else {
        data <- if (grepl("_detected$", case$mutation)) {
            lipidPublicationNegativeDetectionBase(case$route)
        } else {
            lipidPublicationMutateNegativeTable(
                lipidPublicationNegativeBase(case$route),
                case$mutation,
                case$route
            )
        }
        lipidPublicationWriteNegativeTable(data, path)
    }
    files <- if (dir.exists(path)) {
        sort(list.files(path, full.names = TRUE))
    } else {
        path
    }
    records <- lapply(files, function(file) {
        list(
            filename = basename(file),
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

lipidPublicationNegativeImportCheck <- function(case, path) {
    imported <- lipidPublicationReadMember(path, case$route)
    data <- as.data.frame(imported$data)
    sample_ids <- sprintf("LIPID_S%03d", seq_len(6L))
    id_column <- lipidPublicationNegativeIdColumn(case$route)
    if (!id_column %in% names(data) ||
        !all(sample_ids %in% imported$sample_columns)) {
        lipidPublicationAbort("negative imported mapping is incomplete")
    }
    invisible(imported)
}

lipidPublicationNegativeBundleCheck <- function(case, path) {
    expected <- sort(unlist(
        lipidPublicationOutputMembers(case$route, "mixed_lc"),
        use.names = FALSE
    ))
    observed <- sort(list.files(path, all.files = FALSE, no.. = TRUE))
    if (!identical(observed, expected)) {
        lipidPublicationAbort("negative lipidomics bundle membership differs")
    }
    invisible(TRUE)
}

lipidPublicationNegativeSelection <- function(case, path, counter) {
    selector <- get(
        "loadLipidImportAssayPreview",
        envir = asNamespace("MultiScholaR"),
        inherits = FALSE
    )
    msdial_spy <- function(...) {
        counter$reader <- counter$reader + 1L
        lipidPublicationImporter("msdial")(...)
    }
    lipidsearch_spy <- function(...) {
        counter$reader <- counter$reader + 1L
        lipidPublicationImporter("lipidsearch")(...)
    }
    selector(
        assay1File = path,
        importMsdial = msdial_spy,
        importLipidSearch = lipidsearch_spy,
        logInfo = function(...) NULL,
        vendorFormat = case$route
    )
}

lipidPublicationEvaluateNegative <- function(case, path) {
    lipidPublicationValidateNegativeCase(case)
    counter <- new.env(parent = emptyenv())
    counter$reader <- 0L
    observed <- tryCatch({
        if (startsWith(case$mutation, "mixed_")) {
            lipidPublicationNegativeBundleCheck(case, path)
        } else if (grepl("_detected$", case$mutation)) {
            lipidPublicationNegativeSelection(case, path, counter)
        } else if (identical(case$mutation, "requested_format_mismatch")) {
            lipidPublicationAbort("requested and observed lipidomics routes differ")
        } else {
            counter$reader <- counter$reader + 1L
            lipidPublicationNegativeImportCheck(case, path)
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
        lipidPublicationAbort("negative lipidomics outcome differs")
    }
    if (rejected && !inherits(observed, case$expected_condition_class)) {
        lipidPublicationAbort("negative lipidomics condition class differs")
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

lipidPublicationNegativeCase <- function(
    route,
    mutation,
    outcome,
    condition_class,
    reader_invocation,
    profile_id = "lcms_pos"
) {
    list(
        case_id = paste("lipidomics", route, mutation, "v1", sep = "."),
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

lipidPublicationBuildNegativeAuthority <- function() {
    contract_error <- "multischolar_publication_lipidomics_contract_error"
    accepted <- "accepted_edge_characterization"
    rejected <- "classed_rejection"
    cases <- list()
    for (route in names(lipidPublicationCapabilities())) {
        for (mutation in lipidPublicationNegativeBaseMutations(route)) {
            accepted_mutation <- mutation %in% c(
                "zero_quantity", "duplicate_identifier", "reordered_columns"
            )
            cases[[length(cases) + 1L]] <- lipidPublicationNegativeCase(
                route,
                mutation,
                if (accepted_mutation) accepted else rejected,
                if (accepted_mutation) "none" else contract_error,
                !identical(mutation, "requested_format_mismatch")
            )
        }
        for (mutation in lipidPublicationNegativeBundleMutations()) {
            cases[[length(cases) + 1L]] <- lipidPublicationNegativeCase(
                route,
                mutation,
                rejected,
                contract_error,
                FALSE,
                "mixed_lc"
            )
        }
    }
    for (mutation in lipidPublicationNegativeExcludedMutations()) {
        route <- sub("_detected$", "", mutation)
        condition <- if (identical(route, "unknown")) {
            "multischolar_format_unknown"
        } else {
            "multischolar_format_detection_only"
        }
        cases[[length(cases) + 1L]] <- lipidPublicationNegativeCase(
            route,
            mutation,
            rejected,
            condition,
            FALSE
        )
    }
    exclusions_path <- paste0(
        "tests/testdata/omics-performance/lipidomics/exclusions-v1.json"
    )
    list(
        schema = "multischolar.omics_publication_lipidomics_negative_contracts",
        schema_version = .LIPID_PUBLICATION_VERSION,
        negative_id = paste0(
            "multischolar.omics_publication_lipidomics_negative.2026-08-28.v1"
        ),
        owner_ticket_id = .LIPID_PUBLICATION_OWNER,
        status = "frozen_pre_candidate",
        exclusions_binding = lipidPublicationAuthorityBinding(exclusions_path),
        cases = cases,
        unknown_mutation_policy = "reject",
        publication_authority = FALSE
    )
}
