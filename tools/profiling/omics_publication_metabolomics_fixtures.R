metabPublicationFixtureSourcePaths <- function(route, profile_id) {
    custom <- list(
        lcms_pos = "tests/testdata/e2e/metab_lc/lcms_pos_features.tsv",
        lcms_neg = "tests/testdata/e2e/metab_lc/lcms_neg_features.tsv",
        gcms = "tests/testdata/e2e/metab_gc/gcms_features.tsv"
    )
    msdial <- list(
        lcms_pos = "tests/testdata/e2e/metab_lc/seed_lcms_pos.csv",
        lcms_neg = "tests/testdata/e2e/metab_lc/seed_lcms_neg.csv",
        gcms = "tests/testdata/e2e/metab_gc/seed_gcms.csv"
    )
    sources <- if (identical(route, "custom")) custom else msdial
    if (identical(profile_id, "mixed")) {
        return(unname(sources[c("lcms_pos", "lcms_neg", "gcms")]))
    }
    sources[[profile_id]]
}

metabPublicationFixtureReviewPaths <- function(profile_id) {
    if (identical(profile_id, "gcms")) {
        return("tests/testdata/e2e/metab_gc/README.md")
    }
    if (identical(profile_id, "mixed")) {
        return(c(
            "tests/testdata/e2e/metab_lc/README.md",
            "tests/testdata/e2e/metab_gc/README.md",
            "tests/testdata/e2e/metab_combined/README.md"
        ))
    }
    "tests/testdata/e2e/metab_lc/README.md"
}

metabPublicationFixtureBindings <- function(paths) {
    lapply(as.list(paths), function(path) {
        list(path = path, sha256 = publicationFileDigest(path))
    })
}

metabPublicationReadFixtureSource <- function(path) {
    reader <- if (identical(tolower(tools::file_ext(path)), "csv")) {
        utils::read.csv
    } else {
        function(file, ...) utils::read.delim(file, ...)
    }
    reader(
        publicationPath(path),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}

metabPublicationFixtureSampleColumns <- function(data) {
    columns <- grep(
        "^(WT|KO|ctrl|treat)_[0-9]+$",
        names(data),
        value = TRUE
    )
    controls <- columns[grepl("^(WT|ctrl)_", columns)]
    treatments <- columns[grepl("^(KO|treat)_", columns)]
    c(controls, treatments)
}

metabPublicationCanonicalSampleData <- function(data) {
    columns <- metabPublicationFixtureSampleColumns(data)
    values <- data[columns]
    names(values) <- sprintf("METAB_S%03d", seq_along(columns))
    values
}

metabPublicationCanonicalCustom <- function(data, assay, first_id) {
    feature_column <- if ("Feature.Name" %in% names(data)) {
        "Feature.Name"
    } else {
        metabPublicationAbort("custom fixture feature column is absent")
    }
    rt_column <- intersect(c("Retention.Time", "Average Rt(min)"), names(data))
    mz_column <- intersect(c("m.z", "Average Mz"), names(data))
    count <- nrow(data)
    cbind(data.frame(
        feature_id = sprintf(
            "SYNMETAB%09d",
            seq.int(first_id, length.out = count)
        ),
        annotation = as.character(data[[feature_column]]),
        assay = assay,
        rt_minutes = if (length(rt_column)) data[[rt_column[[1L]]]] else NA_real_,
        mz = if (length(mz_column)) data[[mz_column[[1L]]]] else NA_real_,
        check.names = FALSE,
        stringsAsFactors = FALSE
    ), metabPublicationCanonicalSampleData(data))
}

metabPublicationCanonicalMsdial <- function(data, first_id) {
    rt_column <- intersect(c("Average Rt(min)", "RT (min)"), names(data))
    mz_column <- intersect(c("Average Mz", "Precursor m/z"), names(data))
    name_column <- intersect(c("Name", "Metabolite name"), names(data))
    adduct_column <- intersect(c("Adduct type", "Adduct"), names(data))
    count <- nrow(data)
    cbind(data.frame(
        `Alignment ID` = seq.int(first_id, length.out = count),
        `Average Rt(min)` = if (length(rt_column)) {
            data[[rt_column[[1L]]]]
        } else {
            NA_real_
        },
        `Average Mz` = if (length(mz_column)) {
            data[[mz_column[[1L]]]]
        } else {
            NA_real_
        },
        Name = if (length(name_column)) {
            as.character(data[[name_column[[1L]]]])
        } else {
            sprintf("Reviewed_feature_%03d", seq_len(count))
        },
        `Adduct type` = if (length(adduct_column)) {
            as.character(data[[adduct_column[[1L]]]])
        } else {
            NA_character_
        },
        check.names = FALSE,
        stringsAsFactors = FALSE
    ), metabPublicationCanonicalSampleData(data))
}

metabPublicationFixtureAssays <- function(profile_id) {
    unlist(
        metabPublicationAssayProfiles()[[profile_id]]$assays,
        use.names = FALSE
    )
}

metabPublicationBuildFixturePayload <- function(route, profile_id, root) {
    source_paths <- metabPublicationFixtureSourcePaths(route, profile_id)
    assays <- metabPublicationFixtureAssays(profile_id)
    members <- unlist(
        metabPublicationOutputMembers(route, profile_id),
        use.names = FALSE
    )
    payload_root <- file.path(root, "payload")
    dir.create(payload_root, recursive = TRUE, showWarnings = FALSE)
    offset <- 1L
    outputs <- list()
    for (index in seq_along(source_paths)) {
        source <- metabPublicationReadFixtureSource(source_paths[[index]])
        canonical <- if (identical(route, "custom")) {
            metabPublicationCanonicalCustom(source, assays[[index]], offset)
        } else {
            metabPublicationCanonicalMsdial(source, offset)
        }
        offset <- offset + nrow(canonical)
        if (identical(route, "custom") && length(source_paths) > 1L) {
            current <- if (length(outputs)) outputs[[1L]] else canonical[0L, ]
            outputs[[1L]] <- rbind(current, canonical)
        } else {
            outputs[[length(outputs) + 1L]] <- canonical
        }
    }
    if (length(outputs) != length(members)) {
        metabPublicationAbort("fixture output member count differs")
    }
    paths <- character(length(members))
    for (index in seq_along(members)) {
        paths[[index]] <- file.path(payload_root, members[[index]])
        metabPublicationWriteNegativeTable(outputs[[index]], paths[[index]])
    }
    list(
        payload = metabPublicationPayloadBinding(paths),
        paths = paths,
        source_bindings = metabPublicationFixtureBindings(source_paths),
        review_bindings = metabPublicationFixtureBindings(
            metabPublicationFixtureReviewPaths(profile_id)
        )
    )
}

metabPublicationFixtureMemberAssay <- function(filename, data, profile_id) {
    if ("assay" %in% names(data)) return(as.character(data$assay))
    assays <- metabPublicationFixtureAssays(profile_id)
    if (length(assays) == 1L) return(rep(assays[[1L]], nrow(data)))
    match <- assays[startsWith(filename, tolower(assays))]
    if (length(match) != 1L) {
        metabPublicationAbort("fixture member assay identity differs")
    }
    rep(match[[1L]], nrow(data))
}

metabPublicationFixtureDirectTruth <- function(route, profile_id, prepared) {
    sample_ids <- sprintf("METAB_S%03d", seq_len(if (identical(
        route,
        "custom"
    )) 6L else 4L))
    values <- numeric()
    ids <- integer()
    assays <- character()
    log2fc <- numeric()
    for (path in prepared$paths) {
        data <- metabPublicationReadFixtureSource(path)
        matrix <- as.matrix(data[sample_ids])
        storage.mode(matrix) <- "double"
        control <- seq_len(length(sample_ids) / 2L)
        treatment <- seq.int(length(control) + 1L, length(sample_ids))
        values <- c(values, as.numeric(matrix))
        ids <- c(ids, if (identical(route, "custom")) {
            as.integer(sub("^SYNMETAB", "", data$feature_id))
        } else {
            as.integer(data[["Alignment ID"]])
        })
        assays <- c(
            assays,
            metabPublicationFixtureMemberAssay(basename(path), data, profile_id)
        )
        log2fc <- c(log2fc, rowMeans(
            log2(matrix[, treatment, drop = FALSE] + 1),
            na.rm = TRUE
        ) - rowMeans(
            log2(matrix[, control, drop = FALSE] + 1),
            na.rm = TRUE
        ))
    }
    up <- ids[log2fc > 0.25]
    down <- ids[log2fc < -0.25]
    if (!length(up) || !length(down)) {
        metabPublicationAbort("fixture direction authority is incomplete")
    }
    counts <- table(factor(
        assays,
        levels = metabPublicationFixtureAssays(profile_id)
    ))
    list(
        sample_ids = sample_ids,
        values = values,
        ids = ids,
        assays = assays,
        log2fc = log2fc,
        up = up,
        down = down,
        assay_feature_counts = as.list(counts)
    )
}

metabPublicationBuildFixtureTruth <- function(route, profile_id, prepared) {
    direct <- metabPublicationFixtureDirectTruth(route, profile_id, prepared)
    capability <- metabPublicationExpectedCapability(route)
    member_count <- length(prepared$paths)
    dimensions <- list(
        aggregate_feature_count = as.integer(length(direct$ids)),
        assay_feature_counts = direct$assay_feature_counts,
        sample_count = as.integer(length(direct$sample_ids)),
        assay_count = as.integer(length(unique(direct$assays))),
        quantity_count = as.numeric(length(direct$values)),
        payload_member_count = as.integer(member_count)
    )
    profile <- list(
        profile_id = profile_id,
        assays = metabPublicationAssayProfiles()[[profile_id]]$assays,
        payload_mode = if (member_count == 3L) "bundle" else "single",
        member_count = as.integer(member_count)
    )
    workload_id <- metabPublicationWorkloadId(
        route,
        profile_id,
        "fixture_correctness"
    )
    list(
        schema = "multischolar.omics_publication_metabolomics_fixture_truth",
        schema_version = .METAB_PUBLICATION_VERSION,
        owner_ticket_id = .METAB_PUBLICATION_OWNER,
        truth_id = paste0(workload_id, ".truth"),
        workload_id = workload_id,
        workload_class = "fixture_correctness",
        capability = capability,
        assay_profile = profile,
        payload = prepared$payload,
        source_review_bindings = prepared$review_bindings,
        dimensions = dimensions,
        expected_import = list(
            member_count = as.integer(member_count),
            aggregate_feature_count = dimensions$aggregate_feature_count,
            assay_feature_counts = dimensions$assay_feature_counts,
            sample_count = dimensions$sample_count,
            quantity_count = dimensions$quantity_count,
            quantity_na_count = as.numeric(sum(is.na(direct$values))),
            quantity_sum = sum(direct$values, na.rm = TRUE),
            numerical_tolerance = max(1, abs(sum(direct$values, na.rm = TRUE))) *
                1e-12
        ),
        effects = list(
            effect_log2 = min(abs(direct$log2fc[abs(direct$log2fc) > 0.25])),
            minimum_sign_agreement = 1,
            median_margin_fraction = 0,
            up_count = as.integer(length(direct$up)),
            down_count = as.integer(length(direct$down)),
            unaffected_count = as.integer(
                length(direct$ids) - length(direct$up) - length(direct$down)
            ),
            up_feature_ids = as.list(direct$up),
            down_feature_ids = as.list(direct$down),
            assignment_sha256 = publicationObjectDigest(data.frame(
                feature_id = direct$ids,
                log2fc = direct$log2fc
            ))
        ),
        design = list(
            sample_ids = as.list(direct$sample_ids),
            control_ids = as.list(direct$sample_ids[
                seq_len(length(direct$sample_ids) / 2L)
            ]),
            treatment_ids = as.list(direct$sample_ids[
                seq.int(
                    length(direct$sample_ids) / 2L + 1L,
                    length(direct$sample_ids)
                )
            ])
        ),
        mapping = list(
            route = route,
            source = if (identical(route, "custom")) {
                "explicit_user_mapping_contract"
            } else {
                "msdial_schema_autodetection"
            },
            fallback_allowed = FALSE
        ),
        oracle_method =
            "direct_raw_table_arithmetic_and_reviewed_e2e_authority",
        verified_stages = capability$verified_stages,
        limitations = list(
            "Fixture values establish reviewed correctness, not scale.",
            "MS-DIAL evidence cannot widen its reader-characterized boundary."
        ),
        publication_authority = FALSE
    )
}

metabPublicationFreezeFixture <- function(route, profile_id, output_root) {
    if (file.exists(output_root) || dir.exists(output_root)) {
        metabPublicationAbort("fixture freeze output already exists")
    }
    prepared <- metabPublicationBuildFixturePayload(route, profile_id, output_root)
    truth <- metabPublicationBuildFixtureTruth(route, profile_id, prepared)
    truth_path <- file.path(output_root, "truth.json")
    publicationWriteJson(truth, truth_path)
    list(
        route = route,
        profile_id = profile_id,
        payload = prepared$payload,
        payload_paths = prepared$paths,
        source_bindings = prepared$source_bindings,
        review_bindings = prepared$review_bindings,
        truth = list(
            path = truth_path,
            sha256 = publicationFileDigest(truth_path),
            size_bytes = as.numeric(file.info(truth_path)$size),
            record = truth
        )
    )
}
