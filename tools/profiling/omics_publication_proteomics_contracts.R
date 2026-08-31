.PROTEOMICS_PUBLICATION_OWNER <- "OMICS-ART-065"
.PROTEOMICS_PUBLICATION_VERSION <- "1.0.0"

proteomicsPublicationAbort <- function(message, class = "contract_error") {
    publicationAbort(
        message,
        paste0("multischolar_publication_proteomics_", class)
    )
}

proteomicsPublicationCapabilities <- function() {
    list(
        diann = list(
            capability_id = "proteomics.diann.peptide.dia.v1",
            input_format = "diann",
            data_level = "peptide",
            acquisition_mode = "dia"
        ),
        maxquant = list(
            capability_id = "proteomics.maxquant.protein.lfq.v1",
            input_format = "maxquant",
            data_level = "protein",
            acquisition_mode = "dda"
        ),
        fragpipe = list(
            capability_id = "proteomics.fragpipe.protein.lfq.v1",
            input_format = "fragpipe",
            data_level = "protein",
            acquisition_mode = "dda"
        ),
        pd_tmt = list(
            capability_id = "proteomics.pd_tmt.protein.tmt.v1",
            input_format = "pd_tmt",
            data_level = "protein",
            acquisition_mode = "not_recorded"
        )
    )
}

proteomicsPublicationClasses <- function() {
    c(
        "fixture_correctness",
        "representative",
        "operational_heavy",
        "stress"
    )
}

proteomicsPublicationRequireFlag <- function(value, label) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
        proteomicsPublicationAbort(
            paste(label, "must be one non-missing flag")
        )
    }
    invisible(value)
}

proteomicsPublicationRequireNumber <- function(
    value,
    label,
    minimum = -Inf,
    maximum = Inf
) {
    valid <- is.numeric(value) && length(value) == 1L && is.finite(value) &&
        value >= minimum && value <= maximum
    if (!valid) {
        proteomicsPublicationAbort(paste(label, "is outside its domain"))
    }
    invisible(value)
}

proteomicsPublicationRequireInteger <- function(
    value,
    label,
    minimum = 0L
) {
    proteomicsPublicationRequireNumber(value, label, minimum = minimum)
    if (value != as.integer(value)) {
        proteomicsPublicationAbort(paste(label, "must be an integer"))
    }
    invisible(value)
}

proteomicsPublicationRequireDigest <- function(value, label) {
    if (!publicationScalarString(value) || !grepl("^[0-9a-f]{64}$", value)) {
        proteomicsPublicationAbort(paste(label, "must be a SHA-256 digest"))
    }
    invisible(value)
}

proteomicsPublicationValidateBinding <- function(
    binding,
    label,
    require_current = TRUE
) {
    publicationRequireNames(binding, c("path", "sha256"), label)
    proteomicsPublicationRequireDigest(binding$sha256, paste(label, "sha256"))
    if (isTRUE(require_current)) {
        path <- publicationPath(binding$path)
        valid <- file.exists(path) && identical(
            publicationFileDigest(path),
            binding$sha256
        )
        if (!valid) proteomicsPublicationAbort(paste(label, "differs"))
    }
    invisible(binding)
}

proteomicsPublicationExpectedCapability <- function(format) {
    capability <- proteomicsPublicationCapabilities()[[format]]
    if (is.null(capability)) {
        proteomicsPublicationAbort("workload format is unsupported")
    }
    c(list(omic_type = "proteomics"), capability)
}

proteomicsPublicationValidateCapability <- function(capability) {
    publicationRequireNames(capability, c(
        "omic_type", "capability_id", "input_format", "data_level",
        "acquisition_mode"
    ), "Proteomics publication capability")
    expected <- proteomicsPublicationExpectedCapability(capability$input_format)
    if (!identical(capability, expected)) {
        proteomicsPublicationAbort("workload capability differs")
    }
    invisible(capability)
}

proteomicsPublicationValidateDimensions <- function(dimensions, format) {
    fields <- c(
        "protein_count", "peptide_count", "precursor_count", "sample_count",
        "quantity_count", "input_row_count", "input_column_count",
        "assay_count"
    )
    publicationRequireNames(dimensions, fields, "Proteomics dimensions")
    for (field in fields) {
        proteomicsPublicationRequireInteger(
            dimensions[[field]],
            paste("dimensions", field)
        )
    }
    valid <- dimensions$protein_count > 0L && dimensions$sample_count > 0L &&
        identical(dimensions$assay_count, 1L)
    if (identical(format, "diann")) {
        valid <- valid && dimensions$peptide_count > 0L &&
            dimensions$precursor_count >= dimensions$peptide_count &&
            dimensions$quantity_count ==
                dimensions$precursor_count * dimensions$sample_count &&
            dimensions$input_row_count == dimensions$quantity_count
    } else if (format %in% c("maxquant", "fragpipe")) {
        valid <- valid && dimensions$peptide_count == 0L &&
            dimensions$precursor_count == 0L &&
            dimensions$quantity_count ==
                dimensions$protein_count * dimensions$sample_count &&
            dimensions$input_row_count == dimensions$protein_count
    } else if (identical(format, "pd_tmt")) {
        valid <- valid && dimensions$peptide_count > 0L &&
            dimensions$precursor_count == 0L &&
            dimensions$quantity_count ==
                dimensions$peptide_count * dimensions$sample_count &&
            dimensions$input_row_count == dimensions$peptide_count &&
            dimensions$sample_count <= 18L
    }
    if (!valid || dimensions$input_column_count < 2L) {
        proteomicsPublicationAbort("workload dimensions are inconsistent")
    }
    invisible(dimensions)
}

proteomicsPublicationValidateFormatSchema <- function(schema, format) {
    publicationRequireNames(schema, c(
        "schema_id", "schema_version", "orientation", "delimiter",
        "required_columns", "quantity_column_rule", "line_ending",
        "encoding"
    ), "Proteomics format schema")
    orientations <- c(
        diann = "long",
        maxquant = "wide",
        fragpipe = "wide",
        pd_tmt = "wide_peptide_input"
    )
    valid <- publicationScalarString(schema$schema_id) &&
        identical(schema$schema_version, "1.0.0") &&
        identical(schema$orientation, unname(orientations[[format]])) &&
        identical(schema$delimiter, "tab") &&
        length(schema$required_columns) > 0L &&
        publicationScalarString(schema$quantity_column_rule) &&
        identical(schema$line_ending, "LF") &&
        identical(schema$encoding, "UTF-8")
    if (!valid) proteomicsPublicationAbort("format schema differs")
    invisible(schema)
}

proteomicsPublicationValidateRng <- function(rng, workload_class) {
    publicationRequireNames(rng, c(
        "kind", "normal_kind", "sample_kind", "seed_family_id", "seed",
        "streams"
    ), "Proteomics RNG")
    valid <- identical(rng$kind, "L'Ecuyer-CMRG") &&
        identical(rng$normal_kind, "Inversion") &&
        identical(rng$sample_kind, "Rejection") &&
        publicationScalarString(rng$seed_family_id) &&
        is.list(rng$streams) && !anyDuplicated(names(rng$streams))
    if (identical(workload_class, "fixture_correctness")) {
        valid <- valid && is.null(rng$seed) && !length(rng$streams)
    } else {
        valid <- valid && is.numeric(rng$seed) && length(rng$seed) == 1L &&
            rng$seed == as.integer(rng$seed) && length(rng$streams) >= 8L &&
            all(vapply(rng$streams, \(value) {
                is.numeric(value) && length(value) == 1L &&
                    value == as.integer(value)
            }, logical(1)))
    }
    if (!valid) proteomicsPublicationAbort("RNG contract differs")
    invisible(rng)
}

proteomicsPublicationValidateClaimScope <- function(scope) {
    publicationRequireNames(scope, c(
        "evidence_class", "scientific_authority", "performance_authority",
        "cross_project_authority", "promotion_authority", "limitations"
    ), "Proteomics claim scope")
    flags <- c(
        "scientific_authority", "performance_authority",
        "cross_project_authority", "promotion_authority"
    )
    lapply(flags, \(field) {
        proteomicsPublicationRequireFlag(scope[[field]], paste("claim", field))
    })
    valid <- scope$evidence_class %in% c(
        "fixture_correctness", "public_generated", "private_calibrated"
    ) && !isTRUE(scope$promotion_authority) &&
        !isTRUE(scope$cross_project_authority) && length(scope$limitations) > 0L
    if (!valid) proteomicsPublicationAbort("claim scope differs")
    invisible(scope)
}

proteomicsPublicationValidateTruthContract <- function(
    truth,
    workload_class
) {
    publicationRequireNames(truth, c(
        "schema_id", "oracle_type", "independent_from_importer",
        "source_bindings"
    ), "Proteomics truth contract")
    expected_oracle <- if (identical(workload_class, "fixture_correctness")) {
        "hand_reviewed_fixture"
    } else {
        "generator_derived_latent_truth"
    }
    valid <- publicationScalarString(truth$schema_id) &&
        identical(truth$oracle_type, expected_oracle) &&
        isTRUE(truth$independent_from_importer) &&
        length(truth$source_bindings) > 0L && all(vapply(
            truth$source_bindings,
            \(binding) tryCatch({
                proteomicsPublicationValidateBinding(binding, "Truth source")
                TRUE
            }, error = \(error) FALSE),
            logical(1)
        ))
    if (!valid) proteomicsPublicationAbort("truth contract differs")
    invisible(truth)
}

proteomicsPublicationValidateExecution <- function(execution) {
    publicationRequireNames(execution, c(
        "preparation_processes", "generation_inside_measured_worker",
        "generated_payload_committed", "timeout_seconds",
        "maximum_temporary_bytes"
    ), "Proteomics workload execution")
    valid <- identical(execution$preparation_processes, 2L) &&
        !isTRUE(execution$generation_inside_measured_worker) &&
        !isTRUE(execution$generated_payload_committed) &&
        is.numeric(execution$timeout_seconds) && execution$timeout_seconds > 0 &&
        is.numeric(execution$maximum_temporary_bytes) &&
        execution$maximum_temporary_bytes > 0
    if (!valid) proteomicsPublicationAbort("execution contract differs")
    invisible(execution)
}

proteomicsPublicationValidatePrivacy <- function(privacy, workload_class) {
    publicationRequireNames(privacy, c(
        "classification", "private_source", "private_values_retained",
        "source_project_identity_reused"
    ), "Proteomics workload privacy")
    expected <- if (identical(workload_class, "fixture_correctness")) {
        "public_fixture"
    } else {
        "public_generated"
    }
    valid <- identical(privacy$classification, expected) &&
        !isTRUE(privacy$private_source) &&
        !isTRUE(privacy$private_values_retained) &&
        !isTRUE(privacy$source_project_identity_reused)
    if (!valid) proteomicsPublicationAbort("privacy contract differs")
    invisible(privacy)
}

proteomicsPublicationContractFields <- function() {
    c(
        "schema", "schema_version", "contract_id", "owner_ticket_id",
        "status", "workload_id", "workload_class", "capability",
        "format_schema", "dimensions", "model_profile_id",
        "parameter_authority", "source_authority", "split_authority",
        "generator", "rng", "truth_contract", "execution", "privacy",
        "claim_scope", "expected_digests", "publication_authority"
    )
}

proteomicsPublicationValidateGenerator <- function(generator, workload_class) {
    publicationRequireNames(generator, c(
        "mode", "source_bindings", "chunk_rows", "output_filename",
        "truth_filename", "fixture_payload", "fixture_truth"
    ), "Proteomics generator")
    expected_mode <- if (identical(workload_class, "fixture_correctness")) {
        "fixture_replay"
    } else {
        "generated"
    }
    valid <- identical(generator$mode, expected_mode) &&
        length(generator$source_bindings) > 0L &&
        all(vapply(generator$source_bindings, \(binding) {
            tryCatch({
                proteomicsPublicationValidateBinding(binding, "Generator source")
                TRUE
            }, error = \(error) FALSE)
        }, logical(1))) && publicationScalarString(generator$output_filename) &&
        publicationScalarString(generator$truth_filename)
    if (identical(expected_mode, "generated")) {
        valid <- valid && is.numeric(generator$chunk_rows) &&
            generator$chunk_rows == as.integer(generator$chunk_rows) &&
            generator$chunk_rows >= 1L && is.null(generator$fixture_payload) &&
            is.null(generator$fixture_truth)
    } else {
        valid <- valid && is.null(generator$chunk_rows) &&
            is.list(generator$fixture_payload) &&
            is.list(generator$fixture_truth)
        if (valid) {
            proteomicsPublicationValidateBinding(
                generator$fixture_payload,
                "Fixture payload"
            )
            proteomicsPublicationValidateBinding(
                generator$fixture_truth,
                "Fixture truth"
            )
        }
    }
    if (!valid) proteomicsPublicationAbort("generator contract differs")
    invisible(generator)
}

proteomicsPublicationValidateWorkload <- function(contract) {
    publicationRequireNames(
        contract,
        proteomicsPublicationContractFields(),
        "Proteomics workload contract"
    )
    valid <- identical(
        contract$schema,
        "multischolar.omics_publication_proteomics_workload"
    ) && identical(contract$schema_version, .PROTEOMICS_PUBLICATION_VERSION) &&
        identical(contract$owner_ticket_id, .PROTEOMICS_PUBLICATION_OWNER) &&
        identical(contract$status, "frozen") &&
        publicationScalarString(contract$contract_id) &&
        publicationScalarString(contract$workload_id) &&
        contract$workload_class %in% proteomicsPublicationClasses() &&
        publicationScalarString(contract$model_profile_id) &&
        !isTRUE(contract$publication_authority)
    if (!valid) proteomicsPublicationAbort("workload header differs")
    proteomicsPublicationValidateCapability(contract$capability)
    format <- contract$capability$input_format
    proteomicsPublicationValidateDimensions(contract$dimensions, format)
    proteomicsPublicationValidateFormatSchema(contract$format_schema, format)
    for (field in c(
        "parameter_authority", "source_authority", "split_authority"
    )) {
        proteomicsPublicationValidateBinding(contract[[field]], field)
    }
    proteomicsPublicationValidateGenerator(contract$generator, contract$workload_class)
    proteomicsPublicationValidateRng(contract$rng, contract$workload_class)
    proteomicsPublicationValidateTruthContract(
        contract$truth_contract,
        contract$workload_class
    )
    proteomicsPublicationValidateExecution(contract$execution)
    proteomicsPublicationValidatePrivacy(contract$privacy, contract$workload_class)
    proteomicsPublicationValidateClaimScope(contract$claim_scope)
    publicationRequireNames(contract$expected_digests, c(
        "payload_sha256", "truth_sha256"
    ), "Proteomics expected digests")
    for (field in c("payload_sha256", "truth_sha256")) {
        proteomicsPublicationRequireDigest(
            contract$expected_digests[[field]],
            paste("expected", field)
        )
    }
    invisible(contract)
}

proteomicsPublicationContractBasis <- function(contract) {
    basis <- contract
    basis$expected_digests <- list(
        payload_sha256 = strrep("0", 64L),
        truth_sha256 = strrep("0", 64L)
    )
    publicationObjectDigest(basis)
}

proteomicsPublicationParameterOrigins <- function() {
    c(
        "methodological_design", "protocol_authority",
        "public_source_calibrated", "private_aggregate_calibrated"
    )
}

proteomicsPublicationValidateParameterDomain <- function(value, domain) {
    if (!is.list(domain) || !publicationScalarString(domain$type)) {
        proteomicsPublicationAbort("parameter domain is malformed")
    }
    if (domain$type %in% c("numeric_interval", "integer_interval")) {
        publicationRequireNames(
            domain,
            c("type", "minimum", "maximum"),
            "Parameter interval"
        )
        valid <- is.numeric(value) && length(value) == 1L && is.finite(value) &&
            value >= domain$minimum && value <= domain$maximum
        if (identical(domain$type, "integer_interval")) {
            valid <- valid && value == as.integer(value)
        }
    } else if (identical(domain$type, "probability_simplex")) {
        publicationRequireNames(
            domain,
            c("type", "length", "tolerance"),
            "Parameter simplex"
        )
        values <- unlist(value, use.names = FALSE)
        valid <- is.numeric(values) && length(values) == domain$length &&
            all(is.finite(values)) && all(values >= 0) &&
            abs(sum(values) - 1) <= domain$tolerance
    } else {
        valid <- FALSE
    }
    if (!valid) proteomicsPublicationAbort("parameter value is outside its domain")
    invisible(value)
}

proteomicsPublicationAllowedVocabulary <- function() {
    c(
        "declared_synthetic", "mechanistic_simulation", "protocol_minimum",
        "public_source_calibrated", "private_aggregate_calibrated"
    )
}

proteomicsPublicationValidateParameter <- function(parameter) {
    publicationRequireNames(parameter, c(
        "parameter_id", "value", "unit", "domain", "origin",
        "source_binding", "applicable_formats", "applicable_classes",
        "allowed_claim_vocabulary", "limitations"
    ), "Proteomics parameter")
    calibrated <- parameter$origin %in% c(
        "public_source_calibrated", "private_aggregate_calibrated",
        "protocol_authority"
    )
    valid <- publicationScalarString(parameter$parameter_id) &&
        publicationScalarString(parameter$unit) && is.list(parameter$domain) &&
        parameter$origin %in% proteomicsPublicationParameterOrigins() &&
        length(parameter$applicable_formats) > 0L &&
        all(parameter$applicable_formats %in% names(
            proteomicsPublicationCapabilities()
        )) && length(parameter$applicable_classes) > 0L &&
        all(parameter$applicable_classes %in% proteomicsPublicationClasses()) &&
        length(parameter$allowed_claim_vocabulary) > 0L &&
        all(parameter$allowed_claim_vocabulary %in%
            proteomicsPublicationAllowedVocabulary()) &&
        length(parameter$limitations) > 0L
    proteomicsPublicationValidateParameterDomain(parameter$value, parameter$domain)
    if (calibrated) {
        valid <- valid && is.list(parameter$source_binding)
        if (valid) {
            proteomicsPublicationValidateBinding(
                parameter$source_binding,
                "Parameter source"
            )
        }
    } else {
        valid <- valid && is.null(parameter$source_binding) &&
            !any(parameter$allowed_claim_vocabulary %in% c(
                "realistic", "empirical", "instrument_derived",
                "cohort_representative"
            ))
    }
    if (!valid) proteomicsPublicationAbort("parameter authority differs")
    invisible(parameter)
}

proteomicsPublicationValidateParameters <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "parameters_id", "owner_ticket_id",
        "status", "origin_classes", "forbidden_unbound_vocabulary",
        "parameters", "publication_authority"
    ), "Proteomics parameter authority")
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_proteomics_parameters"
    ) && identical(record$schema_version, .PROTEOMICS_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .PROTEOMICS_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        identical(record$origin_classes, as.list(
            proteomicsPublicationParameterOrigins()
        )) && identical(
            record$forbidden_unbound_vocabulary,
            as.list(c(
                "realistic", "empirical", "instrument_derived",
                "cohort_representative"
            ))
        ) && !isTRUE(record$publication_authority)
    ids <- vapply(record$parameters, `[[`, character(1), "parameter_id")
    lapply(record$parameters, proteomicsPublicationValidateParameter)
    if (!valid || anyDuplicated(ids)) {
        proteomicsPublicationAbort("parameter authority differs")
    }
    invisible(record)
}

proteomicsPublicationValidateExclusions <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "exclusions_id", "owner_ticket_id",
        "coverage_binding", "capability_binding", "excluded_capabilities",
        "unknown_policy", "publication_authority"
    ), "Proteomics exclusions")
    lapply(c("coverage_binding", "capability_binding"), \(field) {
        proteomicsPublicationValidateBinding(record[[field]], field)
    })
    expected_ids <- c(
        "proteomics.spectronaut.protein.lfq.v1",
        "proteomics.spectronaut.peptide.lfq.v1"
    )
    ids <- vapply(
        record$excluded_capabilities,
        `[[`,
        character(1),
        "capability_id"
    )
    valid_entries <- all(vapply(record$excluded_capabilities, \(entry) {
        setequal(names(entry), c(
            "capability_id", "format_id", "support_status", "reason",
            "positive_contract_count", "serializer_dispatch",
            "performance_authority", "promotion_authority"
        )) && identical(entry$format_id, "proteomics.spectronaut") &&
            identical(entry$support_status, "advertised_unverified") &&
            publicationScalarString(entry$reason) &&
            identical(entry$positive_contract_count, 0L) &&
            !isTRUE(entry$serializer_dispatch) &&
            !isTRUE(entry$performance_authority) &&
            !isTRUE(entry$promotion_authority)
    }, logical(1)))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_proteomics_exclusions"
    ) && identical(record$schema_version, .PROTEOMICS_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .PROTEOMICS_PUBLICATION_OWNER) &&
        setequal(ids, expected_ids) && !anyDuplicated(ids) && valid_entries &&
        identical(record$unknown_policy, "reject") &&
        !isTRUE(record$publication_authority)
    if (!valid) proteomicsPublicationAbort("exclusion authority differs")
    invisible(record)
}

proteomicsPublicationWorkloadId <- function(format, workload_class) {
    capability <- proteomicsPublicationCapabilities()[[format]]
    if (is.null(capability) || !workload_class %in% proteomicsPublicationClasses()) {
        proteomicsPublicationAbort("workload identity is unsupported")
    }
    paste(capability$capability_id, workload_class, sep = ".")
}

proteomicsPublicationSplitUnitFields <- function() {
    c(
        "source_id", "project_id", "independence_owner_id",
        "aggregate_profile_id", "generator_parameter_profile_id",
        "seed_family_id"
    )
}

proteomicsPublicationValidateSplitAssignment <- function(
    assignment,
    seed_families,
    calibration = FALSE
) {
    fields <- c(
        "assignment_id", "capability_id", "format", "workload_class",
        "split_role", proteomicsPublicationSplitUnitFields(), "seed"
    )
    publicationRequireNames(assignment, fields, "Proteomics split assignment")
    capability <- proteomicsPublicationCapabilities()[[assignment$format]]
    workload_id <- proteomicsPublicationWorkloadId(
        assignment$format,
        assignment$workload_class
    )
    expected_assignment_id <- if (isTRUE(calibration)) {
        paste0(workload_id, ".pilot_calibration")
    } else {
        workload_id
    }
    valid <- !is.null(capability) &&
        identical(assignment$capability_id, capability$capability_id) &&
        identical(assignment$assignment_id, expected_assignment_id) &&
        assignment$workload_class %in% proteomicsPublicationClasses() &&
        all(vapply(
            assignment[proteomicsPublicationSplitUnitFields()],
            publicationScalarString,
            logical(1)
        ))
    if (isTRUE(calibration)) {
        valid <- valid && identical(assignment$split_role, "pilot") &&
            identical(assignment$workload_class, "operational_heavy")
    } else if (identical(assignment$workload_class, "fixture_correctness")) {
        valid <- valid && identical(assignment$split_role, "fixture") &&
            is.null(assignment$seed)
    } else if (identical(assignment$workload_class, "stress")) {
        valid <- valid && identical(assignment$split_role, "stress")
    } else {
        valid <- valid && identical(assignment$split_role, "holdout")
    }
    if (!identical(assignment$split_role, "fixture")) {
        family <- seed_families[[assignment$seed_family_id]]
        valid <- valid && !is.null(family) && is.numeric(assignment$seed) &&
            length(assignment$seed) == 1L && assignment$seed == as.integer(
                assignment$seed
            ) && assignment$seed >= family$minimum_seed &&
            assignment$seed <= family$maximum_seed
    }
    if (!valid) proteomicsPublicationAbort("split assignment differs")
    invisible(assignment)
}

proteomicsPublicationSplitValues <- function(assignments, field) {
    vapply(assignments, `[[`, character(1), field)
}

proteomicsPublicationValidateSplitDisjointness <- function(
    pilot,
    final,
    rules
) {
    expected_rules <- stats::setNames(
        as.list(rep(FALSE, length(proteomicsPublicationSplitUnitFields()))),
        paste0(proteomicsPublicationSplitUnitFields(), "_overlap")
    )
    expected_rules$candidate_result_may_reassign_split <- FALSE
    if (!identical(rules, expected_rules)) {
        proteomicsPublicationAbort("split disjointness rules differ")
    }
    for (field in proteomicsPublicationSplitUnitFields()) {
        pilot_values <- proteomicsPublicationSplitValues(pilot, field)
        final_values <- proteomicsPublicationSplitValues(final, field)
        if (length(intersect(pilot_values, final_values))) {
            proteomicsPublicationAbort(paste("pilot/holdout overlap:", field))
        }
    }
    invisible(TRUE)
}

proteomicsPublicationValidateSplits <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "splits_id", "owner_ticket_id",
        "status", "splits_predecessor", "sources_binding", "split_units",
        "disjointness_rules", "seed_families", "pilot_calibration_assignments",
        "assignments", "readiness", "publication_authority"
    ), "Proteomics split authority")
    proteomicsPublicationValidateBinding(
        record$splits_predecessor,
        "Splits predecessor"
    )
    proteomicsPublicationValidateBinding(record$sources_binding, "Sources")
    predecessor <- publicationReadJson(record$splits_predecessor$path)
    expected_families <- predecessor$generated_seed_families
    names(expected_families) <- vapply(
        expected_families,
        `[[`,
        character(1),
        "seed_family_id"
    )
    families <- record$seed_families
    names(families) <- vapply(
        families,
        `[[`,
        character(1),
        "seed_family_id"
    )
    if (!identical(families, expected_families)) {
        proteomicsPublicationAbort("seed family authority differs")
    }
    lapply(
        record$pilot_calibration_assignments,
        proteomicsPublicationValidateSplitAssignment,
        seed_families = families,
        calibration = TRUE
    )
    lapply(
        record$assignments,
        proteomicsPublicationValidateSplitAssignment,
        seed_families = families,
        calibration = FALSE
    )
    assignment_ids <- vapply(
        c(record$pilot_calibration_assignments, record$assignments),
        `[[`,
        character(1),
        "assignment_id"
    )
    final_ids <- vapply(record$assignments, \(assignment) {
        proteomicsPublicationWorkloadId(
            assignment$format,
            assignment$workload_class
        )
    }, character(1))
    expected_ids <- unlist(lapply(
        names(proteomicsPublicationCapabilities()),
        \(format) vapply(proteomicsPublicationClasses(), \(workload_class) {
            proteomicsPublicationWorkloadId(format, workload_class)
        }, character(1))
    ), use.names = FALSE)
    proteomicsPublicationValidateSplitDisjointness(
        record$pilot_calibration_assignments,
        record$assignments,
        record$disjointness_rules
    )
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_proteomics_splits"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-065") &&
        identical(record$status, "frozen_pre_candidate") &&
        identical(
            record$split_units,
            as.list(proteomicsPublicationSplitUnitFields())
        ) && length(record$pilot_calibration_assignments) == 4L &&
        length(record$assignments) == 16L && !anyDuplicated(assignment_ids) &&
        setequal(final_ids, expected_ids) && !anyDuplicated(final_ids) &&
        !isTRUE(record$readiness$real_project_sources_ready) &&
        isTRUE(record$readiness$generated_assignments_complete) &&
        isTRUE(record$readiness$successor_required_before_candidate) &&
        !isTRUE(record$readiness$candidate_access_allowed) &&
        !isTRUE(record$publication_authority)
    if (!valid) proteomicsPublicationAbort("split authority differs")
    invisible(record)
}

proteomicsPublicationNegativeMutations <- function() {
    list(
        diann = c(
            "missing_protein_group", "missing_quantity", "nonrectangular",
            "zero_and_nonfinite_quantity"
        ),
        maxquant = c(
            "missing_protein_id", "missing_intensity",
            "native_space_lfq_columns", "duplicate_and_contaminant_rows"
        ),
        fragpipe = c(
            "missing_protein_id", "missing_intensity",
            "regular_intensity_without_maxlfq", "duplicate_protein_rows"
        ),
        pd_tmt = c(
            "missing_accession", "missing_abundance",
            "duplicate_normalized_channel", "native_grouped_abundance_columns"
        )
    )
}

proteomicsPublicationValidateNegativeCase <- function(case) {
    publicationRequireNames(case, c(
        "case_id", "format", "capability_id", "mutation", "expected_outcome",
        "expected_condition_class", "scientific_authority",
        "performance_authority", "promotion_authority"
    ), "Proteomics negative case")
    capability <- proteomicsPublicationCapabilities()[[case$format]]
    allowed <- proteomicsPublicationNegativeMutations()[[case$format]]
    rejection_mutations <- setdiff(allowed, c(
        "nonrectangular", "zero_and_nonfinite_quantity",
        "duplicate_and_contaminant_rows",
        "regular_intensity_without_maxlfq", "duplicate_protein_rows"
    ))
    expected <- if (case$mutation %in% rejection_mutations) {
        "classed_rejection"
    } else {
        "accepted_edge_characterization"
    }
    valid <- publicationScalarString(case$case_id) && !is.null(capability) &&
        identical(case$capability_id, capability$capability_id) &&
        case$mutation %in% allowed && identical(case$expected_outcome, expected) &&
        publicationScalarString(case$expected_condition_class) &&
        !isTRUE(case$scientific_authority) &&
        !isTRUE(case$performance_authority) &&
        !isTRUE(case$promotion_authority)
    if (!valid) proteomicsPublicationAbort("negative case differs")
    invisible(case)
}

proteomicsPublicationValidateNegativeAuthority <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "negative_id", "owner_ticket_id",
        "status", "cases", "unknown_mutation_policy", "publication_authority"
    ), "Proteomics negative authority")
    lapply(record$cases, proteomicsPublicationValidateNegativeCase)
    ids <- vapply(record$cases, `[[`, character(1), "case_id")
    formats <- vapply(record$cases, `[[`, character(1), "format")
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_proteomics_negative_contracts"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-065") &&
        identical(record$status, "frozen_pre_candidate") &&
        length(record$cases) == 16L && !anyDuplicated(ids) &&
        all(vapply(names(proteomicsPublicationCapabilities()), \(format) {
            sum(formats == format) == 4L
        }, logical(1))) &&
        identical(record$unknown_mutation_policy, "reject") &&
        !isTRUE(record$publication_authority)
    if (!valid) proteomicsPublicationAbort("negative authority differs")
    invisible(record)
}
