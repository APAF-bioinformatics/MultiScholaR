.METAB_PUBLICATION_ADAPTER_ID <- "multischolar.omics_publication.metabolomics.v1"
.METAB_PUBLICATION_ADAPTER_VERSION <- "1.0.0"

metabPublicationReadContract <- function(path) {
    contract <- publicationReadJson(path)
    metabPublicationValidateWorkload(contract)
    contract
}

metabPublicationPreparationReceipt <- function(contract, payload, truth) {
    list(
        schema = "multischolar.omics_publication_metabolomics_preparation",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-066",
        status = "passed",
        workload_id = contract$workload_id,
        workload_class = contract$workload_class,
        contract_basis_sha256 = metabPublicationContractBasis(contract),
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

metabPublicationValidatePreparation <- function(record, contract) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id", "status",
        "workload_id", "workload_class", "contract_basis_sha256",
        "contract_sha256", "payload", "truth", "route", "profile_id",
        "private_source_opened", "candidate_loaded", "promotion_authority",
        "publication_authority"
    ), "Metabolomics preparation receipt")
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_metabolomics_preparation"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-066") &&
        identical(record$status, "passed") &&
        identical(record$workload_id, contract$workload_id) &&
        identical(record$workload_class, contract$workload_class) &&
        identical(
            record$contract_basis_sha256,
            metabPublicationContractBasis(contract)
        ) && identical(record$contract_sha256, publicationObjectDigest(contract)) &&
        identical(record$route, contract$capability$input_format) &&
        identical(record$profile_id, contract$assay_profile$profile_id) &&
        !isTRUE(record$private_source_opened) && !isTRUE(record$candidate_loaded) &&
        !isTRUE(record$promotion_authority) &&
        !isTRUE(record$publication_authority)
    metabPublicationRequireDigest(
        record$payload$payload_set_sha256,
        "Preparation payload"
    )
    metabPublicationRequireDigest(record$truth$sha256, "Preparation truth")
    if (!valid) metabPublicationAbort("metabolomics preparation receipt differs")
    invisible(record)
}

metabPublicationPrepareGenerated <- function(
    contract,
    output_root,
    verify_expected = TRUE,
    allow_test_contract = FALSE
) {
    metabPublicationValidateWorkload(
        contract,
        allow_test_contract = allow_test_contract,
        validate_authorities = !isTRUE(allow_test_contract)
    )
    if (identical(contract$workload_class, "fixture_correctness")) {
        metabPublicationAbort("fixture workload requires fixture replay")
    }
    if (file.exists(output_root) || dir.exists(output_root)) {
        metabPublicationAbort("metabolomics preparation root already exists")
    }
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    parameters <- publicationReadJson(contract$parameter_authority$path)
    metabPublicationValidateParameters(parameters)
    plan <- metabPublicationModelPlan(
        contract$dimensions$assay_feature_counts,
        contract$dimensions$sample_count,
        parameters,
        contract$rng$streams,
        contract$generator$chunk_rows
    )
    state <- metabPublicationTruthState(plan)
    payload <- metabPublicationSerialize(
        plan,
        contract$capability$input_format,
        contract$assay_profile$profile_id,
        file.path(output_root, "payload"),
        contract$generator$chunk_rows,
        observer = \(block, index) {
            metabPublicationObserveTruth(state, block, index)
        }
    )
    truth_record <- metabPublicationFinalizeTruth(
        state,
        plan,
        contract,
        payload
    )
    metabPublicationValidateTruth(truth_record, contract)
    truth <- metabPublicationWriteTruth(
        truth_record,
        file.path(output_root, contract$generator$truth_filename)
    )
    receipt <- metabPublicationPreparationReceipt(contract, payload, truth)
    metabPublicationValidatePreparation(receipt, contract)
    receipt_path <- file.path(output_root, "preparation-receipt.json")
    publicationWriteJson(receipt, receipt_path)
    if (isTRUE(verify_expected)) {
        valid <- identical(
            payload$payload_set_sha256,
            contract$expected_digests$payload_set_sha256
        ) && identical(truth$sha256, contract$expected_digests$truth_sha256)
        if (!valid) metabPublicationAbort("prepared metabolomics digest differs")
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

metabPublicationCopyFixtureBinding <- function(binding, destination) {
    source <- publicationPath(binding$path)
    if (!isTRUE(file.copy(source, destination, copy.mode = TRUE))) {
        metabPublicationAbort("metabolomics fixture replay copy failed")
    }
    if (!identical(publicationFileDigest(destination), binding$sha256)) {
        metabPublicationAbort("metabolomics fixture replay bytes differ")
    }
    destination
}

metabPublicationPrepareFixture <- function(contract, output_root) {
    metabPublicationValidateWorkload(contract)
    if (!identical(contract$workload_class, "fixture_correctness")) {
        metabPublicationAbort("generated contract requires generation")
    }
    if (file.exists(output_root) || dir.exists(output_root)) {
        metabPublicationAbort("metabolomics fixture output already exists")
    }
    payload_root <- file.path(output_root, "payload")
    dir.create(payload_root, recursive = TRUE, showWarnings = FALSE)
    members <- unlist(contract$generator$output_members, use.names = FALSE)
    copied <- vapply(seq_along(members), function(index) {
        metabPublicationCopyFixtureBinding(
            contract$generator$fixture_payloads[[index]],
            file.path(payload_root, members[[index]])
        )
    }, character(1))
    payload <- metabPublicationPayloadBinding(copied)
    truth_path <- file.path(output_root, contract$generator$truth_filename)
    metabPublicationCopyFixtureBinding(
        contract$generator$fixture_truth,
        truth_path
    )
    truth <- list(
        path = truth_path,
        sha256 = publicationFileDigest(truth_path),
        size_bytes = as.numeric(file.info(truth_path)$size)
    )
    truth_record <- publicationReadJson(truth_path)
    metabPublicationValidateFixtureTruth(truth_record, contract)
    receipt <- metabPublicationPreparationReceipt(contract, payload, truth)
    metabPublicationValidatePreparation(receipt, contract)
    receipt_path <- file.path(output_root, "preparation-receipt.json")
    publicationWriteJson(receipt, receipt_path)
    valid <- identical(
        payload$payload_set_sha256,
        contract$expected_digests$payload_set_sha256
    ) && identical(truth$sha256, contract$expected_digests$truth_sha256)
    if (!valid) metabPublicationAbort("fixture payload or truth digest differs")
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

metabPublicationImporter <- function() {
    get(
        "importMetabMSDIALData",
        envir = asNamespace("MultiScholaR"),
        inherits = FALSE
    )
}

metabPublicationReadMember <- function(path, route) {
    importer <- metabPublicationImporter()
    arguments <- list(filepath = path)
    if (identical(route, "custom")) {
        arguments$metabolite_id_column <- "feature_id"
        arguments$annotation_column <- "annotation"
    }
    suppressMessages(suppressWarnings(do.call(importer, arguments)))
}

metabPublicationMemberAssay <- function(filename, contract) {
    assays <- unlist(contract$assay_profile$assays, use.names = FALSE)
    if (length(assays) == 1L) return(assays[[1L]])
    if (identical(contract$capability$input_format, "custom")) {
        return("embedded_assay_column")
    }
    matches <- assays[startsWith(filename, tolower(assays))]
    if (length(matches) != 1L) {
        metabPublicationAbort("payload member assay identity differs")
    }
    matches[[1L]]
}

metabPublicationImportMembers <- function(contract, payload_root) {
    members <- unlist(contract$generator$output_members, use.names = FALSE)
    lapply(members, \(filename) {
        path <- file.path(payload_root, filename)
        list(
            filename = filename,
            assay = metabPublicationMemberAssay(filename, contract),
            imported = metabPublicationReadMember(
                path,
                contract$capability$input_format
            )
        )
    })
}

metabPublicationImportSummary <- function(contract, members) {
    sample_ids <- sprintf(
        "METAB_S%03d",
        seq_len(contract$dimensions$sample_count)
    )
    rows <- 0L
    quantities <- numeric()
    assays <- character()
    mapping_valid <- TRUE
    extra_detected_samples <- 0L
    effect_ids <- integer()
    effects <- numeric()
    for (member in members) {
        data <- as.data.frame(member$imported$data)
        missing_samples <- setdiff(sample_ids, names(data))
        if (length(missing_samples)) mapping_valid <- FALSE
        values <- as.matrix(data[, intersect(sample_ids, names(data)), drop = FALSE])
        quantities <- c(quantities, as.numeric(values))
        rows <- rows + nrow(data)
        member_assays <- if (identical(
            member$assay,
            "embedded_assay_column"
        )) {
            as.character(data$assay)
        } else {
            rep(member$assay, nrow(data))
        }
        assays <- c(assays, member_assays)
        ids <- if ("feature_id" %in% names(data)) {
            as.integer(sub("^SYNMETAB", "", data$feature_id))
        } else {
            as.integer(data[["Alignment ID"]])
        }
        control <- seq_len(length(sample_ids) / 2L)
        treatment <- seq.int(length(control) + 1L, length(sample_ids))
        row_mean <- function(columns) {
            rowMeans(log2(values[, columns, drop = FALSE] + 1), na.rm = TRUE)
        }
        effect_ids <- c(effect_ids, ids)
        effects <- c(effects, row_mean(treatment) - row_mean(control))
        if (identical(contract$capability$input_format, "msdial")) {
            mapping_valid <- mapping_valid && setequal(
                member$imported$sample_columns,
                sample_ids
            )
        }
        extra_detected_samples <- extra_detected_samples + length(setdiff(
            member$imported$sample_columns,
            sample_ids
        ))
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
        mapping_source = if (identical(
            contract$capability$input_format,
            "custom"
        )) {
            "explicit_user_mapping_contract"
        } else {
            "msdial_schema_autodetection"
        },
        reader_extra_sample_column_count = as.integer(extra_detected_samples),
        effect_ids = effect_ids,
        effects = effects,
        member_count = as.integer(length(members))
    )
}

metabPublicationDirectionEvidence <- function(summary, truth) {
    order <- order(summary$effect_ids, method = "radix")
    effects <- summary$effects[order]
    effect_ids <- summary$effect_ids[order]
    up_ids <- unlist(truth$effects$up_feature_ids, use.names = FALSE)
    down_ids <- unlist(truth$effects$down_feature_ids, use.names = FALSE)
    up_values <- effects[match(up_ids, effect_ids)]
    down_values <- effects[match(down_ids, effect_ids)]
    up_values <- up_values[is.finite(up_values)]
    down_values <- down_values[is.finite(down_values)]
    up_agreement <- mean(up_values > 0)
    down_agreement <- mean(down_values < 0)
    margin <- truth$effects$effect_log2 * truth$effects$median_margin_fraction
    valid <- length(up_values) > 0L && length(down_values) > 0L &&
        up_agreement >= truth$effects$minimum_sign_agreement &&
        down_agreement >= truth$effects$minimum_sign_agreement &&
        stats::median(up_values) > margin &&
        stats::median(down_values) < -margin
    list(
        valid = valid,
        observed_up_sign_agreement = up_agreement,
        observed_down_sign_agreement = down_agreement,
        observed_up_median_log2fc = stats::median(up_values),
        observed_down_median_log2fc = stats::median(down_values),
        minimum_sign_agreement = truth$effects$minimum_sign_agreement
    )
}

metabPublicationValidateImportSummary <- function(summary, truth) {
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
        identical(summary$member_count, expected$member_count)
    if (!valid) metabPublicationAbort("imported metabolomics truth differs")
    invisible(summary)
}

metabPublicationRunImported <- function(contract, payload_root, truth_path) {
    truth <- publicationReadJson(truth_path)
    if (identical(contract$workload_class, "fixture_correctness")) {
        metabPublicationValidateFixtureTruth(truth, contract)
    } else {
        metabPublicationValidateTruth(truth, contract)
    }
    members <- metabPublicationImportMembers(contract, payload_root)
    summary <- metabPublicationImportSummary(contract, members)
    metabPublicationValidateImportSummary(summary, truth)
    direction <- metabPublicationDirectionEvidence(summary, truth)
    if (!isTRUE(direction$valid)) {
        metabPublicationAbort("metabolomics effect direction differs")
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
        adapter_id = .METAB_PUBLICATION_ADAPTER_ID,
        adapter_version = .METAB_PUBLICATION_ADAPTER_VERSION,
        supported_omics = "metabolomics",
        prepare = \(context) {
            output_root <- file.path(context$run_dir, "prepared")
            if (identical(
                context$contract$workload_class,
                "fixture_correctness"
            )) {
                metabPublicationPrepareFixture(context$contract, output_root)
            } else {
                metabPublicationPrepareGenerated(
                    context$contract,
                    output_root,
                    verify_expected = TRUE
                )
            }
        },
        run = \(context) metabPublicationRunImported(
            context$contract,
            context$payload_root,
            context$truth_path
        )
    )
}
