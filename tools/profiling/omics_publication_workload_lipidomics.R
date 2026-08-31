.LIPID_PUBLICATION_ADAPTER_ID <- "multischolar.omics_publication.lipidomics.v1"
.LIPID_PUBLICATION_ADAPTER_VERSION <- "1.0.0"

lipidPublicationReadContract <- function(path) {
    contract <- publicationReadJson(path)
    lipidPublicationValidateWorkload(contract)
    contract
}

lipidPublicationPreparationReceipt <- function(contract, payload, truth) {
    list(
        schema = "multischolar.omics_publication_lipidomics_preparation",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-067",
        status = "passed",
        workload_id = contract$workload_id,
        workload_class = contract$workload_class,
        contract_basis_sha256 = lipidPublicationContractBasis(contract),
        contract_sha256 = publicationObjectDigest(contract),
        payload = payload,
        truth = list(
            filename = basename(truth$path),
            sha256 = truth$sha256,
            size_bytes = truth$size_bytes
        ),
        route = contract$capability$input_format,
        profile_id = contract$assay_profile$profile_id,
        private_source_opened = FALSE,
        candidate_loaded = FALSE,
        promotion_authority = FALSE,
        publication_authority = FALSE
    )
}

lipidPublicationValidatePreparation <- function(record, contract) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id", "status",
        "workload_id", "workload_class", "contract_basis_sha256",
        "contract_sha256", "payload", "truth", "route", "profile_id",
        "private_source_opened", "candidate_loaded", "promotion_authority",
        "publication_authority"
    ), "Lipidomics preparation receipt")
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_lipidomics_preparation"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-067") &&
        identical(record$status, "passed") &&
        identical(record$workload_id, contract$workload_id) &&
        identical(record$workload_class, contract$workload_class) &&
        identical(
            record$contract_basis_sha256,
            lipidPublicationContractBasis(contract)
        ) && identical(record$contract_sha256, publicationObjectDigest(contract)) &&
        identical(record$route, contract$capability$input_format) &&
        identical(record$profile_id, contract$assay_profile$profile_id) &&
        !isTRUE(record$private_source_opened) && !isTRUE(record$candidate_loaded) &&
        !isTRUE(record$promotion_authority) &&
        !isTRUE(record$publication_authority)
    lipidPublicationRequireDigest(
        record$payload$payload_set_sha256,
        "Preparation payload"
    )
    lipidPublicationRequireDigest(record$truth$sha256, "Preparation truth")
    if (!valid) lipidPublicationAbort("lipidomics preparation receipt differs")
    invisible(record)
}

lipidPublicationPrepareGenerated <- function(
    contract,
    output_root,
    verify_expected = TRUE,
    allow_test_contract = FALSE
) {
    lipidPublicationValidateWorkload(
        contract,
        allow_test_contract = allow_test_contract,
        validate_authorities = !isTRUE(allow_test_contract)
    )
    if (identical(contract$workload_class, "fixture_correctness")) {
        lipidPublicationAbort("fixture workload requires fixture replay")
    }
    if (file.exists(output_root) || dir.exists(output_root)) {
        lipidPublicationAbort("lipidomics preparation root already exists")
    }
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    parameters <- publicationReadJson(contract$parameter_authority$path)
    lipidPublicationValidateParameters(parameters)
    plan <- lipidPublicationModelPlan(
        contract$dimensions$assay_feature_counts,
        contract$dimensions$sample_count,
        parameters,
        contract$rng$streams,
        contract$generator$chunk_rows
    )
    state <- lipidPublicationTruthState(plan)
    payload <- lipidPublicationSerialize(
        plan,
        contract$capability$input_format,
        contract$assay_profile$profile_id,
        file.path(output_root, "payload"),
        contract$generator$chunk_rows,
        observer = \(block, index) {
            lipidPublicationObserveTruth(state, block, index)
        }
    )
    truth_record <- lipidPublicationFinalizeTruth(
        state,
        plan,
        contract,
        payload
    )
    lipidPublicationValidateTruth(truth_record, contract)
    truth <- lipidPublicationWriteTruth(
        truth_record,
        file.path(output_root, contract$generator$truth_filename)
    )
    receipt <- lipidPublicationPreparationReceipt(contract, payload, truth)
    lipidPublicationValidatePreparation(receipt, contract)
    receipt_path <- file.path(output_root, "preparation-receipt.json")
    publicationWriteJson(receipt, receipt_path)
    if (isTRUE(verify_expected)) {
        valid <- identical(
            payload$payload_set_sha256,
            contract$expected_digests$payload_set_sha256
        ) && identical(truth$sha256, contract$expected_digests$truth_sha256)
        if (!valid) lipidPublicationAbort("prepared lipidomics digest differs")
    }
    list(
        payload = payload,
        truth = truth,
        receipt = list(
            path = receipt_path,
            sha256 = publicationFileDigest(receipt_path),
            record = receipt
        )
    )
}

lipidPublicationCopyFixtureBinding <- function(binding, destination) {
    source <- publicationPath(binding$path)
    if (!isTRUE(file.copy(source, destination, copy.mode = TRUE))) {
        lipidPublicationAbort("lipidomics fixture replay copy failed")
    }
    if (!identical(publicationFileDigest(destination), binding$sha256)) {
        lipidPublicationAbort("lipidomics fixture replay bytes differ")
    }
    destination
}

lipidPublicationPrepareFixture <- function(contract, output_root) {
    lipidPublicationValidateWorkload(contract)
    if (!identical(contract$workload_class, "fixture_correctness")) {
        lipidPublicationAbort("generated contract requires generation")
    }
    if (file.exists(output_root) || dir.exists(output_root)) {
        lipidPublicationAbort("lipidomics fixture output already exists")
    }
    payload_root <- file.path(output_root, "payload")
    dir.create(payload_root, recursive = TRUE, showWarnings = FALSE)
    members <- unlist(contract$generator$output_members, use.names = FALSE)
    copied <- vapply(seq_along(members), function(index) {
        lipidPublicationCopyFixtureBinding(
            contract$generator$fixture_payloads[[index]],
            file.path(payload_root, members[[index]])
        )
    }, character(1))
    payload <- lipidPublicationPayloadBinding(copied)
    truth_path <- file.path(output_root, contract$generator$truth_filename)
    lipidPublicationCopyFixtureBinding(
        contract$generator$fixture_truth,
        truth_path
    )
    truth <- list(
        path = truth_path,
        sha256 = publicationFileDigest(truth_path),
        size_bytes = as.numeric(file.info(truth_path)$size)
    )
    truth_record <- publicationReadJson(truth_path)
    lipidPublicationValidateFixtureTruth(truth_record, contract)
    receipt <- lipidPublicationPreparationReceipt(contract, payload, truth)
    lipidPublicationValidatePreparation(receipt, contract)
    receipt_path <- file.path(output_root, "preparation-receipt.json")
    publicationWriteJson(receipt, receipt_path)
    valid <- identical(
        payload$payload_set_sha256,
        contract$expected_digests$payload_set_sha256
    ) && identical(truth$sha256, contract$expected_digests$truth_sha256)
    if (!valid) lipidPublicationAbort("fixture payload or truth digest differs")
    list(
        payload = payload,
        truth = truth,
        receipt = list(
            path = receipt_path,
            sha256 = publicationFileDigest(receipt_path),
            record = receipt
        )
    )
}

lipidPublicationImporter <- function(route) {
    function_name <- if (identical(route, "lipidsearch")) {
        "importLipidSearchData"
    } else {
        "importLipidMSDIALData"
    }
    get(
        function_name,
        envir = asNamespace("MultiScholaR"),
        inherits = FALSE
    )
}

lipidPublicationReadMember <- function(path, route) {
    importer <- lipidPublicationImporter(route)
    arguments <- list(filepath = path)
    if (identical(route, "custom")) {
        arguments$lipid_id_column <- "lipid_id"
        arguments$annotation_column <- "annotation"
    }
    suppressMessages(suppressWarnings(do.call(importer, arguments)))
}

lipidPublicationMemberAssay <- function(filename, contract) {
    assays <- unlist(contract$assay_profile$assays, use.names = FALSE)
    if (length(assays) == 1L) return(assays[[1L]])
    matches <- assays[vapply(assays, function(assay) {
        grepl(paste0("_", tolower(assay), "[.]"), filename)
    }, logical(1))]
    if (length(matches) != 1L) {
        lipidPublicationAbort("payload member assay identity differs")
    }
    matches[[1L]]
}

lipidPublicationImportMembers <- function(contract, payload_root) {
    members <- unlist(contract$generator$output_members, use.names = FALSE)
    lapply(members, \(filename) {
        path <- file.path(payload_root, filename)
        list(
            filename = filename,
            assay = lipidPublicationMemberAssay(filename, contract),
            imported = lipidPublicationReadMember(
                path,
                contract$capability$input_format
            )
        )
    })
}

lipidPublicationImportSummary <- function(contract, members) {
    sample_ids <- sprintf(
        "LIPID_S%03d",
        seq_len(contract$dimensions$sample_count)
    )
    rows <- 0L
    quantities <- numeric()
    assays <- character()
    mapping_valid <- TRUE
    extra_detected_samples <- 0L
    effect_ids <- integer()
    effects <- numeric()
    lipid_classes <- character()
    family_ids <- character()
    pair_ids <- character()
    annotation_nonempty <- TRUE
    row_offset <- 0L
    route <- contract$capability$input_format
    for (member in members) {
        data <- as.data.frame(member$imported$data)
        missing_samples <- setdiff(sample_ids, names(data))
        if (length(missing_samples)) mapping_valid <- FALSE
        values <- as.matrix(data[, intersect(sample_ids, names(data)), drop = FALSE])
        quantities <- c(quantities, as.numeric(values))
        rows <- rows + nrow(data)
        member_assays <- rep(member$assay, nrow(data))
        assays <- c(assays, member_assays)
        ids <- if (identical(route, "custom")) {
            as.integer(sub("^SYNLIPID", "", data$lipid_id))
        } else if (identical(route, "msdial")) {
            as.integer(data[["Alignment ID"]])
        } else if ("Idx" %in% names(data)) {
            as.integer(data$Idx)
        } else {
            row_offset + seq_len(nrow(data))
        }
        route_classes <- switch(
            route,
            custom = as.character(data$lipid_class),
            msdial = as.character(data$Ontology),
            lipidsearch = as.character(data$LipidClass)
        )
        route_families <- switch(
            route,
            custom = as.character(data$composition_family_id),
            msdial = sub(";.*$", "", as.character(data$Comment)),
            lipidsearch = as.character(data$FattyAcid)
        )
        route_pairs <- switch(
            route,
            custom = as.character(data$isomer_pair_id),
            msdial = sub("^[^;]*;", "", as.character(data$Comment)),
            lipidsearch = as.character(data$LipidGroup)
        )
        route_pairs[route_pairs %in% c("", "none", route_families)] <- NA_character_
        annotation <- switch(
            route,
            custom = as.character(data$annotation),
            msdial = as.character(data$Name),
            lipidsearch = as.character(data$LipidName)
        )
        lipid_classes <- c(lipid_classes, route_classes)
        family_ids <- c(family_ids, route_families)
        pair_ids <- c(pair_ids, route_pairs)
        annotation_nonempty <- annotation_nonempty &&
            all(!is.na(annotation) & nzchar(annotation))
        control <- seq_len(length(sample_ids) / 2L)
        treatment <- seq.int(length(control) + 1L, length(sample_ids))
        row_mean <- function(columns) {
            rowMeans(log2(values[, columns, drop = FALSE] + 1), na.rm = TRUE)
        }
        effect_ids <- c(effect_ids, ids)
        effects <- c(effects, row_mean(treatment) - row_mean(control))
        if (!identical(route, "custom")) {
            mapping_valid <- mapping_valid && setequal(
                member$imported$sample_columns,
                sample_ids
            )
        }
        extra_detected_samples <- extra_detected_samples + length(setdiff(
            member$imported$sample_columns,
            sample_ids
        ))
        row_offset <- row_offset + nrow(data)
    }
    list(
        aggregate_feature_count = as.integer(rows),
        assay_feature_counts = as.list(table(factor(
            assays,
            levels = unlist(contract$assay_profile$assays, use.names = FALSE)
        ))),
        sample_count = as.integer(length(sample_ids)),
        quantity_count = as.numeric(length(quantities)),
        quantity_na_count = as.numeric(sum(is.na(quantities))),
        quantity_sum = sum(quantities, na.rm = TRUE),
        mapping_valid = mapping_valid,
        mapping_source = switch(
            route,
            lipidsearch = "lipidsearch_schema_autodetection",
            msdial = "msdial_schema_autodetection",
            custom = "explicit_user_mapping_contract"
        ),
        reader_extra_sample_column_count = as.integer(extra_detected_samples),
        lipid_class_count = as.integer(length(unique(lipid_classes))),
        composition_family_count = as.integer(length(unique(family_ids))),
        isomer_like_pair_count = as.integer(length(unique(stats::na.omit(
            pair_ids
        )))),
        annotation_nonempty = annotation_nonempty,
        effect_ids = effect_ids,
        effects = effects,
        member_count = as.integer(length(members))
    )
}

lipidPublicationDirectionEvidence <- function(summary, truth) {
    order <- order(summary$effect_ids, method = "radix")
    effects <- summary$effects[order]
    effect_ids <- summary$effect_ids[order]
    up_ids <- unlist(truth$effects$up_feature_ids, use.names = FALSE)
    down_ids <- unlist(truth$effects$down_feature_ids, use.names = FALSE)
    up_values <- effects[match(up_ids, effect_ids)]
    down_values <- effects[match(down_ids, effect_ids)]
    up_values <- up_values[is.finite(up_values)]
    down_values <- down_values[is.finite(down_values)]
    up_agreement <- if (length(up_values)) mean(up_values > 0) else NULL
    down_agreement <- if (length(down_values)) mean(down_values < 0) else NULL
    margin <- truth$effects$effect_log2 * truth$effects$median_margin_fraction
    up_valid <- !length(up_values) || (
        up_agreement >= truth$effects$minimum_sign_agreement &&
            stats::median(up_values) > margin
    )
    down_valid <- !length(down_values) || (
        down_agreement >= truth$effects$minimum_sign_agreement &&
            stats::median(down_values) < -margin
    )
    valid <- length(c(up_values, down_values)) > 0L && up_valid && down_valid
    list(
        valid = valid,
        observed_up_sign_agreement = up_agreement,
        observed_down_sign_agreement = down_agreement,
        observed_up_median_log2fc = if (length(up_values)) {
            stats::median(up_values)
        } else {
            NULL
        },
        observed_down_median_log2fc = if (length(down_values)) {
            stats::median(down_values)
        } else {
            NULL
        },
        minimum_sign_agreement = truth$effects$minimum_sign_agreement
    )
}

lipidPublicationValidateImportSummary <- function(summary, truth) {
    expected <- truth$expected_import
    valid <- identical(
        summary$aggregate_feature_count,
        expected$aggregate_feature_count
    ) && identical(
        publicationObjectDigest(summary$assay_feature_counts),
        publicationObjectDigest(expected$assay_feature_counts)
    ) && identical(summary$sample_count, expected$sample_count) &&
        summary$quantity_count == expected$quantity_count &&
        summary$quantity_na_count == expected$quantity_na_count &&
        abs(summary$quantity_sum - expected$quantity_sum) <=
            expected$numerical_tolerance && isTRUE(summary$mapping_valid) &&
        isTRUE(summary$annotation_nonempty) &&
        identical(summary$lipid_class_count, expected$lipid_class_count) &&
        identical(
            summary$composition_family_count,
            expected$composition_family_count
        ) && identical(
            summary$isomer_like_pair_count,
            expected$isomer_like_pair_count
        ) &&
        identical(summary$member_count, expected$member_count)
    if (!valid) lipidPublicationAbort("imported lipidomics truth differs")
    invisible(summary)
}

lipidPublicationRunImported <- function(contract, payload_root, truth_path) {
    truth <- publicationReadJson(truth_path)
    if (identical(contract$workload_class, "fixture_correctness")) {
        lipidPublicationValidateFixtureTruth(truth, contract)
    } else {
        lipidPublicationValidateTruth(truth, contract)
    }
    members <- lipidPublicationImportMembers(contract, payload_root)
    summary <- lipidPublicationImportSummary(contract, members)
    lipidPublicationValidateImportSummary(summary, truth)
    direction <- lipidPublicationDirectionEvidence(summary, truth)
    if (!isTRUE(direction$valid)) {
        lipidPublicationAbort("lipidomics effect direction differs")
    }
    retained_summary <- summary[setdiff(names(summary), c(
        "effect_ids", "effects"
    ))]
    list(
        status = "passed",
        workflow_evidence = list(
            truth_valid = TRUE,
            route = contract$capability$input_format,
            profile_id = contract$assay_profile$profile_id,
            verified_stages = contract$capability$verified_stages,
            import = retained_summary,
            differential_direction = direction,
            promotion_authority = FALSE
        ),
        query_evidence = list(),
        session_evidence = list(status = "memory_only", resource_count = 0L),
        report_evidence = list(status = "outside_workload", file_count = 0L),
        retained = lapply(members, \(member) member$imported$data)
    )
}

omicsWorkloadAdapter <- function() {
    list(
        adapter_id = .LIPID_PUBLICATION_ADAPTER_ID,
        adapter_version = .LIPID_PUBLICATION_ADAPTER_VERSION,
        supported_omics = "lipidomics",
        prepare = \(context) {
            output_root <- file.path(context$run_dir, "prepared")
            if (identical(
                context$contract$workload_class,
                "fixture_correctness"
            )) {
                lipidPublicationPrepareFixture(context$contract, output_root)
            } else {
                lipidPublicationPrepareGenerated(
                    context$contract,
                    output_root,
                    verify_expected = TRUE
                )
            }
        },
        run = \(context) lipidPublicationRunImported(
            context$contract,
            context$payload_root,
            context$truth_path
        )
    )
}
