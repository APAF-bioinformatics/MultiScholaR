lipidPublicationFixtureSourceMap <- function() {
    list(
        LCMS_Pos = paste0(
            "tests/testdata/e2e/lipid_canonical/",
            "lipidsearch_lcms_pos.txt"
        ),
        LCMS_Neg = paste0(
            "tests/testdata/e2e/lipid_canonical/",
            "lipidsearch_lcms_neg.txt"
        ),
        GCMS = paste0(
            "tests/testdata/e2e/lipid_canonical/",
            "lipidsearch_gcms.txt"
        )
    )
}

lipidPublicationFixtureAssays <- function(profile_id) {
    unlist(
        lipidPublicationAssayProfiles()[[profile_id]]$assays,
        use.names = FALSE
    )
}

lipidPublicationFixtureSourcePaths <- function(route, profile_id) {
    if (!identical(route, "lipidsearch")) return(character())
    assays <- unlist(
        lipidPublicationAssayProfiles()[[profile_id]]$assays,
        use.names = FALSE
    )
    unname(unlist(lipidPublicationFixtureSourceMap()[assays]))
}

lipidPublicationFixtureBindings <- function(paths) {
    lapply(as.list(paths), function(path) {
        list(path = path, sha256 = publicationFileDigest(path))
    })
}

lipidPublicationReadFixtureTable <- function(path) {
    delimiter <- if (identical(tolower(tools::file_ext(path)), "csv")) {
        ","
    } else {
        "\t"
    }
    utils::read.table(
        publicationPath(path),
        header = TRUE,
        sep = delimiter,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = ""
    )
}

lipidPublicationCanonicalSamples <- function(data) {
    columns <- grep("^(WT|KO)_[0-9]+$", names(data), value = TRUE)
    columns <- c(
        columns[startsWith(columns, "WT_")],
        columns[startsWith(columns, "KO_")]
    )
    values <- data[columns]
    names(values) <- sprintf("LIPID_S%03d", seq_along(columns))
    values
}

lipidPublicationCanonicalLipidSearch <- function(data) {
    metadata <- c(
        "LipidName", "LipidClass", "FattyAcid", "IonType", "BaseRt",
        "CalcMz", "Grade"
    )
    if (!all(metadata %in% names(data))) {
        lipidPublicationAbort("reviewed LipidSearch fixture schema differs")
    }
    cbind(
        data[metadata],
        lipidPublicationCanonicalSamples(data)
    )
}

lipidPublicationFixtureStreams <- function(route, profile_id) {
    route_index <- match(route, names(lipidPublicationCapabilities()))
    profile_index <- match(profile_id, names(lipidPublicationAssayProfiles()))
    seed <- 510000L + route_index * 10000L + profile_index * 100L
    as.list(stats::setNames(
        as.integer(seed + seq_len(8L) * 100L),
        c(
            "hierarchy", "feature_offsets", "class_residual", "batch",
            "residual", "mcar", "mar", "mnar"
        )
    ))
}

lipidPublicationBuildGeneratedReaderFixture <- function(
    route,
    profile_id,
    root
) {
    parameters_path <- paste0(
        "tests/testdata/omics-performance/lipidomics/parameters-v1.json"
    )
    parameters <- publicationReadJson(parameters_path)
    lipidPublicationValidateParameters(parameters)
    assays <- unlist(
        lipidPublicationAssayProfiles()[[profile_id]]$assays,
        use.names = FALSE
    )
    counts <- stats::setNames(as.list(rep(12L, length(assays))), assays)
    plan <- lipidPublicationModelPlan(
        counts,
        sample_count = 6L,
        parameter_authority = parameters,
        streams = lipidPublicationFixtureStreams(route, profile_id),
        chunk_rows = 100L
    )
    payload_root <- file.path(root, "payload")
    payload <- lipidPublicationSerialize(
        plan,
        route,
        profile_id,
        payload_root,
        chunk_rows = 100L
    )
    members <- unlist(lipidPublicationOutputMembers(route, profile_id))
    source_paths <- c(
        parameters_path,
        "tools/profiling/omics_publication_lipidomics_model.R",
        "tools/profiling/omics_publication_lipidomics_serializers.R"
    )
    list(
        payload = payload,
        paths = file.path(payload_root, members),
        source_authority_bindings = lipidPublicationFixtureBindings(source_paths)
    )
}

lipidPublicationBuildReviewedFixture <- function(route, profile_id, root) {
    source_paths <- lipidPublicationFixtureSourcePaths(route, profile_id)
    assays <- unlist(
        lipidPublicationAssayProfiles()[[profile_id]]$assays,
        use.names = FALSE
    )
    members <- unlist(lipidPublicationOutputMembers(route, profile_id))
    payload_root <- file.path(root, "payload")
    dir.create(payload_root, recursive = TRUE, showWarnings = FALSE)
    paths <- character(length(members))
    for (index in seq_along(members)) {
        source <- lipidPublicationReadFixtureTable(source_paths[[index]])
        canonical <- lipidPublicationCanonicalLipidSearch(source)
        paths[[index]] <- file.path(payload_root, members[[index]])
        lipidPublicationWriteTable(
            canonical,
            paths[[index]],
            delimiter = "\t",
            append = FALSE
        )
    }
    authorities <- c(
        source_paths,
        "tests/testdata/e2e/lipid_canonical/README.md",
        "tools/profiling/omics_publication_lipidomics_fixtures.R"
    )
    list(
        payload = lipidPublicationPayloadBinding(paths),
        paths = paths,
        source_authority_bindings = lipidPublicationFixtureBindings(
            unique(authorities)
        )
    )
}

lipidPublicationBuildFixturePayload <- function(route, profile_id, root) {
    if (identical(route, "lipidsearch")) {
        lipidPublicationBuildReviewedFixture(route, profile_id, root)
    } else {
        lipidPublicationBuildGeneratedReaderFixture(route, profile_id, root)
    }
}

lipidPublicationFixtureAssay <- function(filename, profile_id) {
    assays <- unlist(
        lipidPublicationAssayProfiles()[[profile_id]]$assays,
        use.names = FALSE
    )
    if (length(assays) == 1L) return(assays[[1L]])
    found <- assays[vapply(assays, function(assay) {
        grepl(paste0("_", tolower(assay), "[.]"), filename)
    }, logical(1))]
    if (length(found) != 1L) {
        lipidPublicationAbort("fixture member assay identity differs")
    }
    found[[1L]]
}

lipidPublicationFixtureColumns <- function(route, data) {
    if (identical(route, "custom")) {
        return(list(
            id = as.integer(sub("^SYNLIPID", "", data$lipid_id)),
            lipid_class = as.character(data$lipid_class),
            family = as.character(data$composition_family_id),
            pair = as.character(data$isomer_pair_id),
            annotation = as.character(data$annotation)
        ))
    }
    if (identical(route, "msdial")) {
        comment <- as.character(data$Comment)
        pair <- sub("^[^;]*;", "", comment)
        pair[pair == "none"] <- NA_character_
        return(list(
            id = as.integer(data[["Alignment ID"]]),
            lipid_class = as.character(data$Ontology),
            family = sub(";.*$", "", comment),
            pair = pair,
            annotation = as.character(data$Name)
        ))
    }
    list(
        id = NULL,
        lipid_class = as.character(data$LipidClass),
        family = as.character(data$FattyAcid),
        pair = rep(NA_character_, nrow(data)),
        annotation = as.character(data$LipidName)
    )
}

lipidPublicationFixtureDirectTruth <- function(route, profile_id, prepared) {
    sample_ids <- sprintf("LIPID_S%03d", seq_len(6L))
    values <- numeric()
    ids <- integer()
    assays <- character()
    classes <- character()
    families <- character()
    pairs <- character()
    annotations <- character()
    log2fc <- numeric()
    offset <- 0L
    for (path in prepared$paths) {
        data <- lipidPublicationReadFixtureTable(path)
        matrix <- as.matrix(data[sample_ids])
        storage.mode(matrix) <- "double"
        columns <- lipidPublicationFixtureColumns(route, data)
        member_ids <- columns$id
        if (is.null(member_ids)) member_ids <- offset + seq_len(nrow(data))
        control <- seq_len(3L)
        treatment <- seq.int(4L, 6L)
        values <- c(values, as.numeric(matrix))
        ids <- c(ids, member_ids)
        assays <- c(assays, rep(
            lipidPublicationFixtureAssay(basename(path), profile_id),
            nrow(data)
        ))
        classes <- c(classes, columns$lipid_class)
        families <- c(families, columns$family)
        pairs <- c(pairs, columns$pair)
        annotations <- c(annotations, columns$annotation)
        log2fc <- c(log2fc, rowMeans(
            log2(matrix[, treatment, drop = FALSE] + 1),
            na.rm = TRUE
        ) - rowMeans(
            log2(matrix[, control, drop = FALSE] + 1),
            na.rm = TRUE
        ))
        offset <- offset + nrow(data)
    }
    up <- ids[is.finite(log2fc) & log2fc > 0.25]
    down <- ids[is.finite(log2fc) & log2fc < -0.25]
    if (!length(c(up, down)) ||
        !all(is.finite(values[!is.na(values)])) || any(!nzchar(annotations))) {
        lipidPublicationAbort("fixture arithmetic authority is incomplete")
    }
    counts <- table(factor(
        assays,
        levels = unlist(
            lipidPublicationAssayProfiles()[[profile_id]]$assays,
            use.names = FALSE
        )
    ))
    list(
        sample_ids = sample_ids,
        values = values,
        ids = ids,
        assays = assays,
        classes = classes,
        families = families,
        pairs = pairs,
        log2fc = log2fc,
        up = up,
        down = down,
        assay_feature_counts = as.list(counts)
    )
}

lipidPublicationFixtureEvidenceClass <- function(route, profile_id) {
    if (!identical(route, "lipidsearch")) {
        return("generated_reader_characterization")
    }
    if (profile_id %in% c("gcms", "mixed_lc_gcms")) {
        "application_assay_name_smoke"
    } else {
        "reviewed_fixture_correctness"
    }
}

lipidPublicationFixtureOracleMethod <- function(route, profile_id) {
    if (!identical(route, "lipidsearch")) {
        return("generated_reader_schema_truth")
    }
    if (profile_id %in% c("gcms", "mixed_lc_gcms")) {
        "direct_arithmetic_application_assay_name_smoke"
    } else {
        "direct_raw_table_arithmetic_and_reviewed_e2e_authority"
    }
}

lipidPublicationBuildFixtureTruth <- function(route, profile_id, prepared) {
    direct <- lipidPublicationFixtureDirectTruth(route, profile_id, prepared)
    capability <- lipidPublicationExpectedCapability(route)
    member_count <- length(prepared$paths)
    dimensions <- list(
        aggregate_feature_count = as.integer(length(direct$ids)),
        assay_feature_counts = direct$assay_feature_counts,
        sample_count = 6L,
        assay_count = as.integer(length(unique(direct$assays))),
        quantity_count = as.numeric(length(direct$values)),
        payload_member_count = as.integer(member_count)
    )
    profile <- list(
        profile_id = profile_id,
        assays = lipidPublicationAssayProfiles()[[profile_id]]$assays,
        payload_mode = if (member_count == 2L) "bundle" else "single",
        member_count = as.integer(member_count)
    )
    workload_id <- lipidPublicationWorkloadId(
        route,
        profile_id,
        "fixture_correctness"
    )
    pair_count <- length(unique(stats::na.omit(direct$pairs)))
    list(
        schema = "multischolar.omics_publication_lipidomics_fixture_truth",
        schema_version = .LIPID_PUBLICATION_VERSION,
        owner_ticket_id = .LIPID_PUBLICATION_OWNER,
        truth_id = paste0(workload_id, ".truth"),
        workload_id = workload_id,
        workload_class = "fixture_correctness",
        capability = capability,
        assay_profile = profile,
        payload = prepared$payload,
        source_authority_bindings = prepared$source_authority_bindings,
        dimensions = dimensions,
        expected_import = list(
            member_count = as.integer(member_count),
            aggregate_feature_count = dimensions$aggregate_feature_count,
            assay_feature_counts = dimensions$assay_feature_counts,
            sample_count = 6L,
            quantity_count = dimensions$quantity_count,
            quantity_na_count = as.numeric(sum(is.na(direct$values))),
            quantity_sum = sum(direct$values, na.rm = TRUE),
            lipid_class_count = as.integer(length(unique(direct$classes))),
            composition_family_count = as.integer(length(unique(direct$families))),
            isomer_like_pair_count = as.integer(pair_count),
            numerical_tolerance = max(
                1,
                abs(sum(direct$values, na.rm = TRUE))
            ) * 1e-12
        ),
        effects = list(
            effect_log2 = min(abs(direct$log2fc[
                is.finite(direct$log2fc) & abs(direct$log2fc) > 0.25
            ])),
            minimum_sign_agreement = 1,
            median_margin_fraction = 0,
            up_count = as.integer(length(direct$up)),
            down_count = as.integer(length(direct$down)),
            unaffected_count = as.integer(
                length(direct$ids) - length(direct$up) - length(direct$down)
            ),
            internal_standard_effect_count = 0L,
            up_feature_ids = as.list(direct$up),
            down_feature_ids = as.list(direct$down),
            assignment_sha256 = publicationObjectDigest(data.frame(
                feature_id = direct$ids,
                log2fc = direct$log2fc
            ))
        ),
        design = list(
            sample_ids = as.list(direct$sample_ids),
            control_ids = as.list(direct$sample_ids[1:3]),
            treatment_ids = as.list(direct$sample_ids[4:6])
        ),
        mapping = list(
            route = route,
            source = switch(
                route,
                lipidsearch = "lipidsearch_schema_autodetection",
                msdial = "msdial_schema_autodetection",
                custom = "explicit_user_mapping_contract"
            ),
            fallback_allowed = FALSE
        ),
        oracle_method = lipidPublicationFixtureOracleMethod(route, profile_id),
        verified_stages = capability$verified_stages,
        fixture_evidence_class = lipidPublicationFixtureEvidenceClass(
            route,
            profile_id
        ),
        gc_vendor_authority = FALSE,
        chemical_identity_claim = FALSE,
        limitations = list(
            "Fixture values establish correctness, not cross-project scale.",
            "Synthetic labels and isomer-like pairs are not chemical identities.",
            "GCMS-named LipidSearch fixtures are application naming smoke only."
        ),
        publication_authority = FALSE
    )
}

lipidPublicationFreezeFixture <- function(route, profile_id, output_root) {
    if (file.exists(output_root) || dir.exists(output_root)) {
        lipidPublicationAbort("fixture freeze output already exists")
    }
    prepared <- lipidPublicationBuildFixturePayload(route, profile_id, output_root)
    truth <- lipidPublicationBuildFixtureTruth(route, profile_id, prepared)
    truth_path <- file.path(output_root, "truth.json")
    publicationWriteJson(truth, truth_path)
    list(
        route = route,
        profile_id = profile_id,
        payload = prepared$payload,
        payload_paths = prepared$paths,
        source_authority_bindings = prepared$source_authority_bindings,
        truth = list(
            path = truth_path,
            sha256 = publicationFileDigest(truth_path),
            size_bytes = as.numeric(file.info(truth_path)$size),
            record = truth
        )
    )
}
