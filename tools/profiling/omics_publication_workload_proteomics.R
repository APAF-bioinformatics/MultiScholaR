.PROTEOMICS_PUBLICATION_ADAPTER_ID <-
    "multischolar.omics_publication.proteomics.v1"
.PROTEOMICS_PUBLICATION_ADAPTER_VERSION <- "1.0.0"

proteomicsPublicationReadContract <- function(path) {
    contract <- publicationReadJson(path)
    proteomicsPublicationValidateWorkload(contract)
    contract
}

proteomicsPublicationReadParameters <- function(contract) {
    parameters <- publicationReadJson(contract$parameter_authority$path)
    proteomicsPublicationValidateParameters(parameters)
    parameters
}

proteomicsPublicationPayloadBinding <- function(binding) {
    list(
        filename = basename(binding$path),
        sha256 = binding$sha256,
        size_bytes = binding$size_bytes
    )
}

proteomicsPublicationPreparationReceipt <- function(
    contract,
    payload,
    truth,
    truth_record
) {
    list(
        schema = "multischolar.omics_publication_proteomics_preparation",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-065",
        status = "passed",
        workload_id = contract$workload_id,
        workload_class = contract$workload_class,
        capability = contract$capability,
        contract_basis_sha256 = proteomicsPublicationContractBasis(contract),
        contract_sha256 = publicationObjectDigest(contract),
        generator = contract$generator,
        rng = contract$rng,
        parameter_authority = contract$parameter_authority,
        source_authority = contract$source_authority,
        split_authority = contract$split_authority,
        payload = proteomicsPublicationPayloadBinding(payload),
        truth = list(
            filename = basename(truth$path),
            sha256 = truth$sha256,
            size_bytes = truth$size_bytes,
            record_sha256 = publicationObjectDigest(truth_record)
        ),
        private_source_opened = FALSE,
        candidate_loaded = FALSE,
        publication_authority = FALSE
    )
}

proteomicsPublicationValidatePreparation <- function(record, contract) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id", "status",
        "workload_id", "workload_class", "capability",
        "contract_basis_sha256", "contract_sha256", "generator", "rng",
        "parameter_authority", "source_authority", "split_authority",
        "payload", "truth", "private_source_opened", "candidate_loaded",
        "publication_authority"
    ), "Proteomics preparation receipt")
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_proteomics_preparation"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-065") &&
        identical(record$status, "passed") &&
        identical(record$workload_id, contract$workload_id) &&
        identical(record$workload_class, contract$workload_class) &&
        identical(record$capability, contract$capability) &&
        identical(
            record$contract_basis_sha256,
            proteomicsPublicationContractBasis(contract)
        ) && identical(record$contract_sha256, publicationObjectDigest(contract)) &&
        identical(record$generator, contract$generator) &&
        identical(record$rng, contract$rng) &&
        identical(record$parameter_authority, contract$parameter_authority) &&
        identical(record$source_authority, contract$source_authority) &&
        identical(record$split_authority, contract$split_authority) &&
        !isTRUE(record$private_source_opened) &&
        !isTRUE(record$candidate_loaded) &&
        !isTRUE(record$publication_authority)
    for (binding in c(list(record$payload), list(record$truth))) {
        proteomicsPublicationRequireDigest(binding$sha256, "Preparation digest")
    }
    if (!valid) proteomicsPublicationAbort("preparation receipt differs")
    invisible(record)
}

proteomicsPublicationPrepareGenerated <- function(
    contract,
    output_root,
    verify_expected = TRUE
) {
    proteomicsPublicationValidateWorkload(contract)
    if (identical(contract$workload_class, "fixture_correctness")) {
        proteomicsPublicationAbort("fixture contracts require fixture replay")
    }
    if (file.exists(output_root) || dir.exists(output_root)) {
        proteomicsPublicationAbort("preparation output root already exists")
    }
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    parameters <- proteomicsPublicationReadParameters(contract)
    plan <- proteomicsPublicationModelPlan(contract, parameters)
    state <- proteomicsPublicationTruthState(plan)
    payload_path <- file.path(output_root, contract$generator$output_filename)
    payload <- proteomicsPublicationSerialize(
        plan,
        payload_path,
        observer = \(block, index) {
            proteomicsPublicationObserveTruth(state, block, index)
        }
    )
    truth_record <- proteomicsPublicationFinalizeTruth(
        state,
        plan,
        proteomicsPublicationPayloadBinding(payload)
    )
    proteomicsPublicationValidateTruth(truth_record, contract)
    truth_path <- file.path(output_root, contract$generator$truth_filename)
    truth <- proteomicsPublicationWriteTruth(truth_record, truth_path)
    receipt <- proteomicsPublicationPreparationReceipt(
        contract,
        payload,
        truth,
        truth_record
    )
    proteomicsPublicationValidatePreparation(receipt, contract)
    receipt_path <- file.path(output_root, "preparation-receipt.json")
    publicationWriteJson(receipt, receipt_path)
    if (isTRUE(verify_expected)) {
        valid <- identical(
            payload$sha256,
            contract$expected_digests$payload_sha256
        ) && identical(truth$sha256, contract$expected_digests$truth_sha256)
        if (!valid) {
            proteomicsPublicationAbort("prepared payload or truth digest differs")
        }
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

proteomicsPublicationCopyBinding <- function(binding, destination) {
    source <- publicationPath(binding$path)
    if (!isTRUE(file.copy(source, destination, copy.mode = TRUE))) {
        proteomicsPublicationAbort("fixture replay copy failed")
    }
    if (!identical(publicationFileDigest(destination), binding$sha256)) {
        proteomicsPublicationAbort("fixture replay bytes differ")
    }
    list(
        path = destination,
        sha256 = binding$sha256,
        size_bytes = as.numeric(file.info(destination)$size)
    )
}

proteomicsPublicationPrepareFixture <- function(contract, output_root) {
    proteomicsPublicationValidateWorkload(contract)
    if (!identical(contract$workload_class, "fixture_correctness")) {
        proteomicsPublicationAbort("generated contracts require generation")
    }
    if (file.exists(output_root) || dir.exists(output_root)) {
        proteomicsPublicationAbort("fixture output root already exists")
    }
    dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
    payload <- proteomicsPublicationCopyBinding(
        contract$generator$fixture_payload,
        file.path(output_root, contract$generator$output_filename)
    )
    truth <- proteomicsPublicationCopyBinding(
        contract$generator$fixture_truth,
        file.path(output_root, contract$generator$truth_filename)
    )
    truth_record <- publicationReadJson(truth$path)
    proteomicsPublicationValidateFixtureTruth(truth_record, contract)
    receipt <- proteomicsPublicationPreparationReceipt(
        contract,
        payload,
        truth,
        truth_record
    )
    proteomicsPublicationValidatePreparation(receipt, contract)
    receipt_path <- file.path(output_root, "preparation-receipt.json")
    publicationWriteJson(receipt, receipt_path)
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

proteomicsPublicationImporter <- function(format) {
    namespace <- asNamespace("MultiScholaR")
    function_name <- switch(
        format,
        diann = "importDIANNData",
        maxquant = "importMaxQuantData",
        fragpipe = "importFragPipeData",
        pd_tmt = "importProteomeDiscovererTMTData",
        proteomicsPublicationAbort("importer format is unsupported")
    )
    get(function_name, envir = namespace, inherits = FALSE)
}

proteomicsPublicationImportSummary <- function(imported) {
    data <- as.data.frame(imported$data)
    mapping <- imported$column_mapping
    protein_ids <- as.character(data[[mapping$protein_col]])
    samples <- as.character(data[[mapping$run_col]])
    quantities <- as.numeric(data[[mapping$quantity_col]])
    peptide_ids <- if (!is.null(mapping$peptide_col)) {
        as.character(data[[mapping$peptide_col]])
    } else if ("Annotated.Sequence" %in% names(data)) {
        as.character(data[["Annotated.Sequence"]])
    } else {
        NULL
    }
    list(
        data = data,
        protein_ids = protein_ids,
        peptide_ids = peptide_ids,
        samples = samples,
        quantities = quantities,
        data_type = imported$data_type,
        row_count = as.numeric(nrow(data)),
        protein_count = as.integer(length(unique(protein_ids))),
        peptide_count = as.integer(if (is.null(peptide_ids)) {
            0L
        } else {
            length(unique(peptide_ids))
        }),
        sample_count = as.integer(length(unique(samples))),
        quantity_na_count = as.numeric(sum(is.na(quantities))),
        quantity_sum = sum(quantities, na.rm = TRUE),
        quantity_min = min(quantities, na.rm = TRUE),
        quantity_max = max(quantities, na.rm = TRUE)
    )
}

proteomicsPublicationValidateImportSummary <- function(summary, truth) {
    expected <- truth$expected_import
    tolerance <- expected$numerical_tolerance
    checks <- c(
        data_type = identical(summary$data_type, expected$data_type),
        row_count = summary$row_count == expected$row_count,
        protein_count = summary$protein_count == expected$protein_count,
        peptide_count = summary$peptide_count == expected$peptide_count,
        sample_count = summary$sample_count == expected$sample_count,
        quantity_na_count = summary$quantity_na_count ==
            expected$quantity_na_count,
        quantity_sum = abs(summary$quantity_sum - expected$quantity_sum) <=
            tolerance,
        quantity_min = abs(summary$quantity_min - expected$quantity_min) <= 1e-6,
        quantity_max = abs(summary$quantity_max - expected$quantity_max) <= 1e-6
    )
    if (!all(checks)) {
        proteomicsPublicationAbort(paste(
            "imported payload differs from truth:",
            paste(names(checks)[!checks], collapse = ", ")
        ))
    }
    invisible(summary)
}

proteomicsPublicationDirectionEvidence <- function(summary, truth) {
    values <- log2(summary$quantities + 1)
    control <- grepl("CTRL_", summary$samples, fixed = TRUE)
    treatment <- grepl("TREAT_", summary$samples, fixed = TRUE)
    means <- function(selected) {
        ids <- summary$protein_ids[selected]
        sums <- rowsum(values[selected], ids, reorder = FALSE, na.rm = TRUE)
        counts <- rowsum(
            as.numeric(is.finite(values[selected])),
            ids,
            reorder = FALSE,
            na.rm = TRUE
        )
        stats::setNames(as.numeric(sums / counts), rownames(sums))
    }
    control_mean <- means(control)
    treatment_mean <- means(treatment)
    proteins <- sprintf(
        "SYNPROT%08d",
        seq_len(truth$hierarchy$protein_count)
    )
    effects <- treatment_mean[proteins] - control_mean[proteins]
    up <- seq_len(truth$effects$up_count)
    down <- seq.int(
        truth$effects$up_count + 1L,
        truth$effects$up_count + truth$effects$down_count
    )
    up_values <- effects[up]
    down_values <- effects[down]
    up_finite <- is.finite(up_values)
    down_finite <- is.finite(down_values)
    up_agreement <- mean(up_values[up_finite] > 0)
    down_agreement <- mean(down_values[down_finite] < 0)
    margin <- truth$effects$effect_log2 * truth$effects$median_margin_fraction
    valid <- length(up_values[up_finite]) > 0L &&
        length(down_values[down_finite]) > 0L &&
        up_agreement >= truth$effects$minimum_sign_agreement &&
        down_agreement >= truth$effects$minimum_sign_agreement &&
        stats::median(up_values[up_finite]) > margin &&
        stats::median(down_values[down_finite]) < -margin
    list(
        valid = valid,
        expected_up_count = truth$effects$up_count,
        expected_down_count = truth$effects$down_count,
        minimum_sign_agreement = truth$effects$minimum_sign_agreement,
        observed_up_sign_agreement = up_agreement,
        observed_down_sign_agreement = down_agreement,
        observed_up_median_log2fc = stats::median(up_values[up_finite]),
        observed_down_median_log2fc = stats::median(down_values[down_finite]),
        observed_up_count = sum(effects > margin, na.rm = TRUE),
        observed_down_count = sum(effects < -margin, na.rm = TRUE),
        finite_effect_count = sum(is.finite(effects))
    )
}

proteomicsPublicationRunImported <- function(contract, payload_path, truth_path) {
    truth <- publicationReadJson(truth_path)
    fixture <- identical(
        truth$schema,
        "multischolar.omics_publication_proteomics_fixture_truth"
    )
    if (fixture) {
        proteomicsPublicationValidateFixtureTruth(truth, contract)
    } else {
        proteomicsPublicationValidateTruth(truth, contract)
    }
    importer <- proteomicsPublicationImporter(contract$capability$input_format)
    imported <- suppressMessages(suppressWarnings(importer(payload_path)))
    summary <- proteomicsPublicationImportSummary(imported)
    proteomicsPublicationValidateImportSummary(summary, truth)
    direction <- if (fixture) {
        list(
            valid = TRUE,
            oracle = "hand_reviewed_fixture_e2e",
            expected = truth$effects
        )
    } else {
        proteomicsPublicationDirectionEvidence(summary, truth)
    }
    if (!isTRUE(direction$valid)) {
        proteomicsPublicationAbort("imported effect directions differ from truth")
    }
    list(
        status = "passed",
        workflow_evidence = list(
            truth_valid = TRUE,
            evidence_class = contract$claim_scope$evidence_class,
            format = contract$capability$input_format,
            import = summary[setdiff(names(summary), c(
                "data", "protein_ids", "peptide_ids", "samples", "quantities"
            ))],
            differential_direction = direction,
            reachable_stage_claim = "import_only",
            promotion_authority = FALSE
        ),
        query_evidence = list(),
        session_evidence = list(status = "memory_only", resource_count = 0L),
        report_evidence = list(status = "not_in_workload", file_count = 0L),
        retained = imported$data
    )
}

proteomicsPublicationAdapterPrepare <- function(context) {
    prepared <- if (identical(
        context$contract$workload_class,
        "fixture_correctness"
    )) {
        proteomicsPublicationPrepareFixture(
            context$contract,
            file.path(context$run_dir, "prepared")
        )
    } else {
        proteomicsPublicationPrepareGenerated(
            context$contract,
            file.path(context$run_dir, "prepared"),
            verify_expected = TRUE
        )
    }
    list(
        payload_path = prepared$payload$path,
        truth_path = prepared$truth$path,
        metadata = list(
            workload_id = context$contract$workload_id,
            preparation_sha256 = prepared$receipt$sha256,
            evidence_class = context$contract$claim_scope$evidence_class
        )
    )
}

proteomicsPublicationAdapterRun <- function(context) {
    proteomicsPublicationRunImported(
        context$contract,
        context$payload_path,
        context$truth_path
    )
}

omicsWorkloadAdapter <- function() {
    list(
        adapter_id = .PROTEOMICS_PUBLICATION_ADAPTER_ID,
        adapter_version = .PROTEOMICS_PUBLICATION_ADAPTER_VERSION,
        supported_omics = "proteomics",
        prepare = proteomicsPublicationAdapterPrepare,
        run = proteomicsPublicationAdapterRun
    )
}
