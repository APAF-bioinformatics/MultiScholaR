.METAB_PUBLICATION_OWNER <- "OMICS-ART-066"
.METAB_PUBLICATION_VERSION <- "1.0.0"

metabPublicationAbort <- function(message, class = "contract_error") {
    publicationAbort(
        message,
        paste0("multischolar_publication_metabolomics_", class)
    )
}

metabPublicationCapabilities <- function() {
    list(
        custom = list(
            capability_id = "metabolomics.custom.metabolite.standard.v1",
            input_format = "custom",
            support_tier = "full_workflow",
            verified_stages = list(
                "import", "design", "qc", "normalization_ruv", "session",
                "differential_abundance", "summary_report_export"
            )
        ),
        msdial = list(
            capability_id = "metabolomics.msdial.metabolite.standard.v1",
            input_format = "msdial",
            support_tier = "reader_boundary",
            verified_stages = list("import", "design")
        )
    )
}

metabPublicationAssayProfiles <- function() {
    list(
        lcms_pos = list(assays = list("LCMS_Pos")),
        lcms_neg = list(assays = list("LCMS_Neg")),
        gcms = list(assays = list("GCMS")),
        mixed = list(assays = list("LCMS_Pos", "LCMS_Neg", "GCMS"))
    )
}

metabPublicationClasses <- function() {
    c(
        "fixture_correctness", "representative", "operational_heavy", "stress"
    )
}

metabPublicationRequireDigest <- function(value, label) {
    if (!publicationScalarString(value) || !grepl("^[0-9a-f]{64}$", value)) {
        metabPublicationAbort(paste(label, "must be a SHA-256 digest"))
    }
    invisible(value)
}

metabPublicationScalarFlag <- function(value) {
    is.logical(value) && length(value) == 1L && !is.na(value)
}

metabPublicationValidateBinding <- function(
    binding,
    label,
    require_current = TRUE
) {
    publicationRequireNames(binding, c("path", "sha256"), label)
    metabPublicationRequireDigest(binding$sha256, paste(label, "sha256"))
    if (isTRUE(require_current)) {
        path <- publicationPath(binding$path)
        if (!file.exists(path) || !identical(
            publicationFileDigest(path),
            binding$sha256
        )) {
            metabPublicationAbort(paste(label, "differs"))
        }
    }
    invisible(binding)
}

metabPublicationExpectedCapability <- function(route) {
    expected <- metabPublicationCapabilities()[[route]]
    if (is.null(expected)) metabPublicationAbort("metabolomics route is unsupported")
    c(
        list(
            omic_type = "metabolomics",
            data_level = "metabolite",
            acquisition_mode = "not_recorded"
        ),
        expected
    )
}

metabPublicationValidateCapability <- function(capability) {
    publicationRequireNames(capability, c(
        "omic_type", "data_level", "acquisition_mode", "capability_id",
        "input_format", "support_tier", "verified_stages"
    ), "Metabolomics capability")
    if (!identical(
        capability,
        metabPublicationExpectedCapability(capability$input_format)
    )) {
        metabPublicationAbort("metabolomics capability differs")
    }
    invisible(capability)
}

metabPublicationValidateAssayProfile <- function(profile, route) {
    publicationRequireNames(profile, c(
        "profile_id", "assays", "payload_mode", "member_count"
    ), "Metabolomics assay profile")
    expected <- metabPublicationAssayProfiles()[[profile$profile_id]]
    expected_mode <- if (identical(profile$profile_id, "mixed") &&
        identical(route, "msdial")) {
        "bundle"
    } else {
        "single"
    }
    expected_count <- if (identical(expected_mode, "bundle")) 3L else 1L
    valid <- !is.null(expected) && identical(profile$assays, expected$assays) &&
        identical(profile$payload_mode, expected_mode) &&
        identical(profile$member_count, expected_count)
    if (!valid) metabPublicationAbort("metabolomics assay profile differs")
    invisible(profile)
}

metabPublicationValidateDimensions <- function(dimensions, profile) {
    publicationRequireNames(dimensions, c(
        "aggregate_feature_count", "assay_feature_counts", "sample_count",
        "assay_count", "quantity_count", "payload_member_count"
    ), "Metabolomics dimensions")
    counts <- unlist(dimensions$assay_feature_counts, use.names = TRUE)
    expected_assays <- unlist(profile$assays, use.names = FALSE)
    valid <- is.numeric(counts) && all(counts > 0) &&
        identical(sort(names(counts)), sort(expected_assays)) &&
        dimensions$aggregate_feature_count == sum(counts) &&
        dimensions$sample_count >= 4L && dimensions$sample_count %% 2L == 0L &&
        identical(dimensions$assay_count, as.integer(length(expected_assays))) &&
        dimensions$quantity_count ==
            dimensions$aggregate_feature_count * dimensions$sample_count &&
        identical(dimensions$payload_member_count, profile$member_count)
    if (!valid) metabPublicationAbort("metabolomics dimensions differ")
    invisible(dimensions)
}

metabPublicationExpectedScale <- function(profile_id, workload_class) {
    single <- !identical(profile_id, "mixed")
    switch(
        workload_class,
        fixture_correctness = list(
            minimum_features = if (single) 3L else 9L,
            sample_count = NULL
        ),
        representative = list(
            minimum_features = if (single) 20000L else 30000L,
            sample_count = 24L
        ),
        operational_heavy = list(
            minimum_features = 80000L,
            sample_count = 48L
        ),
        stress = list(
            minimum_features = if (single) 160000L else 180000L,
            sample_count = 96L
        ),
        metabPublicationAbort("metabolomics workload class is unsupported")
    )
}

metabPublicationValidateScale <- function(contract) {
    scale <- metabPublicationExpectedScale(
        contract$assay_profile$profile_id,
        contract$workload_class
    )
    dimensions <- contract$dimensions
    valid <- dimensions$aggregate_feature_count >= scale$minimum_features
    if (!is.null(scale$sample_count)) {
        valid <- valid && dimensions$sample_count >= scale$sample_count
    }
    if (identical(contract$assay_profile$profile_id, "mixed")) {
        counts <- unlist(dimensions$assay_feature_counts, use.names = FALSE)
        valid <- valid && length(counts) == 3L && all(counts > 0L)
    }
    if (!valid) metabPublicationAbort("metabolomics workload scale is below floor")
    invisible(contract)
}

metabPublicationWorkloadId <- function(route, profile, workload_class) {
    capability <- metabPublicationCapabilities()[[route]]
    if (is.null(capability) || is.null(metabPublicationAssayProfiles()[[profile]]) ||
        !workload_class %in% metabPublicationClasses()) {
        metabPublicationAbort("metabolomics workload identity is unsupported")
    }
    paste(capability$capability_id, profile, workload_class, sep = ".")
}

metabPublicationValidateMapping <- function(mapping) {
    publicationRequireNames(mapping, c(
        "mapping_id", "route", "profile_id", "requested_format", "observed_format",
        "mapping_source", "fallback_allowed", "bundle_member_mappings",
        "verified_stages", "promotion_authority"
    ), "Metabolomics mapping authority")
    capability <- metabPublicationCapabilities()[[mapping$route]]
    profile <- metabPublicationAssayProfiles()[[mapping$profile_id]]
    expected_source <- if (identical(mapping$route, "custom")) {
        "explicit_user_mapping_contract"
    } else {
        "msdial_schema_autodetection"
    }
    member_assays <- vapply(
        mapping$bundle_member_mappings,
        `[[`,
        character(1),
        "assay"
    )
    valid <- !is.null(capability) && !is.null(profile) &&
        identical(mapping$requested_format, mapping$route) &&
        identical(mapping$observed_format, mapping$route) &&
        identical(mapping$mapping_source, expected_source) &&
        !isTRUE(mapping$fallback_allowed) &&
        setequal(member_assays, unlist(profile$assays, use.names = FALSE)) &&
        !anyDuplicated(member_assays) &&
        identical(mapping$verified_stages, capability$verified_stages) &&
        !isTRUE(mapping$promotion_authority)
    if (!valid) metabPublicationAbort("metabolomics mapping differs")
    invisible(mapping)
}

metabPublicationValidateMappingAuthority <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "mapping_authority_id", "owner_ticket_id",
        "status", "mappings", "identity_alias_allowed", "unknown_policy",
        "publication_authority"
    ), "Metabolomics mapping authority")
    lapply(record$mappings, metabPublicationValidateMapping)
    ids <- vapply(record$mappings, `[[`, character(1), "mapping_id")
    pairs <- vapply(record$mappings, \(mapping) {
        paste(mapping$route, mapping$profile_id, sep = "::")
    }, character(1))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_metabolomics_mapping_authority"
    ) && identical(record$schema_version, .METAB_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .METAB_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        length(record$mappings) == 8L && !anyDuplicated(ids) &&
        !anyDuplicated(pairs) && !isTRUE(record$identity_alias_allowed) &&
        identical(record$unknown_policy, "reject") &&
        !isTRUE(record$publication_authority)
    if (!valid) metabPublicationAbort("mapping authority differs")
    invisible(record)
}

metabPublicationValidateSupportBoundaries <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "boundary_id", "owner_ticket_id",
        "coverage_binding", "capability_binding", "routes",
        "reachable_stage_claim_allowed", "unknown_policy",
        "publication_authority"
    ), "Metabolomics support boundaries")
    metabPublicationValidateBinding(record$coverage_binding, "Coverage")
    metabPublicationValidateBinding(record$capability_binding, "Capabilities")
    route_ids <- vapply(record$routes, `[[`, character(1), "route")
    valid_routes <- all(vapply(record$routes, \(route) {
        expected <- metabPublicationCapabilities()[[route$route]]
        setequal(names(route), c(
            "route", "capability_id", "support_tier", "verified_stages",
            "full_workflow_claim_allowed", "generated_scope"
        )) && !is.null(expected) &&
            identical(route$capability_id, expected$capability_id) &&
            identical(route$support_tier, expected$support_tier) &&
            identical(route$verified_stages, expected$verified_stages) &&
            identical(
                route$full_workflow_claim_allowed,
                identical(route$route, "custom")
            ) && publicationScalarString(route$generated_scope)
    }, logical(1)))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_metabolomics_support_boundaries"
    ) && identical(record$schema_version, .METAB_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .METAB_PUBLICATION_OWNER) &&
        setequal(route_ids, names(metabPublicationCapabilities())) &&
        !anyDuplicated(route_ids) && valid_routes &&
        !isTRUE(record$reachable_stage_claim_allowed) &&
        identical(record$unknown_policy, "reject") &&
        !isTRUE(record$publication_authority)
    if (!valid) metabPublicationAbort("support boundary differs")
    invisible(record)
}

metabPublicationValidateExclusions <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "exclusions_id", "owner_ticket_id",
        "capability_binding", "formats", "reader_invocation_allowed",
        "fallback_allowed", "promotion_authority", "publication_authority"
    ), "Metabolomics exclusions")
    metabPublicationValidateBinding(record$capability_binding, "Capabilities")
    expected <- c("progenesis", "xcms", "compound_discoverer", "unknown")
    ids <- vapply(record$formats, `[[`, character(1), "format")
    valid_formats <- all(vapply(record$formats, \(format) {
        expected_status <- if (identical(format$format, "unknown")) {
            "unknown_rejected"
        } else {
            "detection_only"
        }
        setequal(names(format), c(
            "format", "support_status", "capability_count", "reason"
        )) && identical(format$support_status, expected_status) &&
            identical(format$capability_count, 0L) &&
            publicationScalarString(format$reason)
    }, logical(1)))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_metabolomics_exclusions"
    ) && identical(record$schema_version, .METAB_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .METAB_PUBLICATION_OWNER) &&
        setequal(ids, expected) && !anyDuplicated(ids) && valid_formats &&
        !isTRUE(record$reader_invocation_allowed) &&
        !isTRUE(record$fallback_allowed) &&
        !isTRUE(record$promotion_authority) &&
        !isTRUE(record$publication_authority)
    if (!valid) metabPublicationAbort("metabolomics exclusions differ")
    invisible(record)
}

metabPublicationValidateMsdialBundle <- function(bundle) {
    publicationRequireNames(bundle, c(
        "bundle_id", "profile_id", "member_assays", "member_count",
        "member_schema_ids", "single_custom_table_substitution_allowed",
        "bundle_digest"
    ), "MS-DIAL bundle")
    expected <- metabPublicationAssayProfiles()[[bundle$profile_id]]
    valid <- identical(bundle$profile_id, "mixed") && !is.null(expected) &&
        identical(bundle$member_assays, expected$assays) &&
        identical(bundle$member_count, 3L) &&
        length(bundle$member_schema_ids) == 3L &&
        !anyDuplicated(unlist(bundle$member_schema_ids, use.names = FALSE)) &&
        !isTRUE(bundle$single_custom_table_substitution_allowed)
    metabPublicationRequireDigest(bundle$bundle_digest, "MS-DIAL bundle")
    basis <- bundle
    basis$bundle_digest <- strrep("0", 64L)
    valid <- valid && identical(
        bundle$bundle_digest,
        publicationObjectDigest(basis)
    )
    if (!valid) metabPublicationAbort("MS-DIAL bundle differs")
    invisible(bundle)
}

metabPublicationValidateMsdialBundles <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "bundles_id", "owner_ticket_id",
        "status", "mapping_authority", "bundles", "member_order_semantic",
        "missing_or_duplicate_member_policy", "custom_substitution_allowed",
        "publication_authority"
    ), "MS-DIAL bundle authority")
    metabPublicationValidateBinding(
        record$mapping_authority,
        "MS-DIAL mapping authority"
    )
    lapply(record$bundles, metabPublicationValidateMsdialBundle)
    ids <- vapply(record$bundles, `[[`, character(1), "bundle_id")
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_metabolomics_msdial_bundles"
    ) && identical(record$schema_version, .METAB_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .METAB_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        length(record$bundles) == length(metabPublicationClasses()) &&
        !anyDuplicated(ids) && !isTRUE(record$member_order_semantic) &&
        identical(record$missing_or_duplicate_member_policy, "reject") &&
        !isTRUE(record$custom_substitution_allowed) &&
        !isTRUE(record$publication_authority)
    if (!valid) metabPublicationAbort("MS-DIAL bundle authority differs")
    invisible(record)
}

metabPublicationParameterOrigins <- function() {
    c(
        "methodological_design", "protocol_authority",
        "public_source_calibrated", "private_scale_calibrated"
    )
}

metabPublicationParameterConsumers <- function() {
    c(
        treatment_fraction = "model.design",
        batch_count = "model.design",
        technical_replicate_size = "model.design",
        sample_correlation = "model.correlation",
        technical_replicate_correlation = "model.correlation",
        effect_fraction = "model.effects",
        effect_log2 = "model.effects_and_truth",
        lcms_pos_mean_log2 = "model.assay_means",
        lcms_neg_mean_log2 = "model.assay_means",
        gcms_mean_log2 = "model.assay_means",
        base_log2_sd = "model.correlated_group_baseline",
        correlated_group_size = "model.correlated_groups",
        group_offset_log2_sd = "model.correlated_group_offset",
        feature_offset_log2_sd = "model.feature_offset",
        batch_log2_sd = "model.batch_offset",
        residual_sigma_floor = "model.heteroscedastic_residual",
        residual_sigma_slope = "model.heteroscedastic_residual",
        heteroscedastic_reference_log2 = "model.heteroscedastic_residual",
        mcar_probability = "model.mcar",
        mar_logit_intercept = "model.mar",
        mar_batch_log_odds = "model.mar",
        detection_midpoint_log2 = "model.mnar",
        detection_slope_log2 = "model.mnar",
        lcms_pos_censor_log2 = "model.left_censoring",
        lcms_neg_censor_log2 = "model.left_censoring",
        gcms_censor_log2 = "model.left_censoring",
        internal_standard_fraction = "model.internal_standards",
        duplicate_fraction = "negative.duplicate_identifiers",
        annotation_width = "serializer.annotations",
        minimum_effect_sign_agreement = "truth.effect_direction",
        effect_median_margin_fraction = "truth.effect_direction"
    )
}

metabPublicationValidateParameterDomain <- function(value, domain) {
    publicationRequireNames(
        domain,
        c("type", "minimum", "maximum"),
        "Metabolomics parameter domain"
    )
    valid <- domain$type %in% c("numeric_interval", "integer_interval") &&
        is.numeric(value) && length(value) == 1L && is.finite(value) &&
        value >= domain$minimum && value <= domain$maximum
    if (identical(domain$type, "integer_interval")) {
        valid <- valid && value == as.integer(value)
    }
    if (!valid) metabPublicationAbort("parameter value is outside its domain")
    invisible(value)
}

metabPublicationValidateParameter <- function(parameter) {
    publicationRequireNames(parameter, c(
        "parameter_id", "value", "unit", "domain", "origin",
        "source_binding", "applicable_routes", "applicable_profiles",
        "applicable_classes", "allowed_claim_vocabulary", "limitations"
    ), "Metabolomics parameter")
    metabPublicationValidateParameterDomain(parameter$value, parameter$domain)
    design <- identical(parameter$origin, "methodological_design")
    calibrated <- parameter$origin %in% setdiff(
        metabPublicationParameterOrigins(),
        "methodological_design"
    )
    valid <- publicationScalarString(parameter$parameter_id) &&
        publicationScalarString(parameter$unit) &&
        parameter$origin %in% metabPublicationParameterOrigins() &&
        length(parameter$applicable_routes) > 0L &&
        all(parameter$applicable_routes %in% names(
            metabPublicationCapabilities()
        )) && length(parameter$applicable_profiles) > 0L &&
        all(parameter$applicable_profiles %in% names(
            metabPublicationAssayProfiles()
        )) && length(parameter$applicable_classes) > 0L &&
        all(parameter$applicable_classes %in% metabPublicationClasses()) &&
        length(parameter$allowed_claim_vocabulary) > 0L &&
        length(parameter$limitations) > 0L
    if (design) {
        valid <- valid && is.null(parameter$source_binding) &&
            !any(parameter$allowed_claim_vocabulary %in% c(
                "realistic", "empirical", "instrument_derived"
            ))
    }
    if (calibrated) {
        valid <- valid && is.list(parameter$source_binding)
        if (valid) {
            metabPublicationValidateBinding(
                parameter$source_binding,
                "Metabolomics parameter source"
            )
        }
    }
    if (!valid) metabPublicationAbort("metabolomics parameter differs")
    invisible(parameter)
}

metabPublicationValidateParameters <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "parameters_id", "owner_ticket_id",
        "status", "origin_classes", "forbidden_unbound_vocabulary",
        "parameters", "publication_authority"
    ), "Metabolomics parameters")
    lapply(record$parameters, metabPublicationValidateParameter)
    ids <- vapply(record$parameters, `[[`, character(1), "parameter_id")
    consumers <- metabPublicationParameterConsumers()
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_metabolomics_parameters"
    ) && identical(record$schema_version, .METAB_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .METAB_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        identical(
            record$origin_classes,
            as.list(metabPublicationParameterOrigins())
        ) && identical(
            record$forbidden_unbound_vocabulary,
            as.list(c("realistic", "empirical", "instrument_derived"))
        ) && !anyDuplicated(ids) && setequal(ids, names(consumers)) &&
        all(nzchar(unname(consumers))) && !isTRUE(record$publication_authority)
    if (!valid) metabPublicationAbort("metabolomics parameters differ")
    invisible(record)
}

metabPublicationParameterValues <- function(record) {
    metabPublicationValidateParameters(record)
    values <- lapply(record$parameters, `[[`, "value")
    names(values) <- vapply(
        record$parameters,
        `[[`,
        character(1),
        "parameter_id"
    )
    values
}

metabPublicationWorkloadFields <- function() {
    c(
        "schema", "schema_version", "contract_id", "owner_ticket_id",
        "status", "workload_id", "workload_class", "capability",
        "assay_profile", "dimensions", "model_profile_id",
        "parameter_authority", "source_authority", "split_authority",
        "mapping_authority", "contract_mapping_id", "support_boundary",
        "generator", "rng",
        "truth_contract", "execution", "privacy", "claim_scope",
        "expected_digests", "publication_authority"
    )
}

metabPublicationValidateGenerator <- function(generator, contract) {
    publicationRequireNames(generator, c(
        "mode", "source_bindings", "chunk_rows", "output_members",
        "truth_filename", "fixture_payloads", "fixture_truth"
    ), "Metabolomics generator")
    fixture <- identical(contract$workload_class, "fixture_correctness")
    expected_mode <- if (fixture) "fixture_replay" else "generated"
    member_count <- contract$dimensions$payload_member_count
    valid <- identical(generator$mode, expected_mode) &&
        length(generator$source_bindings) > 0L &&
        all(vapply(generator$source_bindings, \(binding) {
            tryCatch({
                metabPublicationValidateBinding(binding, "Generator source")
                TRUE
            }, error = \(error) FALSE)
        }, logical(1))) && length(generator$output_members) == member_count &&
        !anyDuplicated(unlist(generator$output_members, use.names = FALSE)) &&
        publicationScalarString(generator$truth_filename)
    if (fixture) {
        valid <- valid && is.null(generator$chunk_rows) &&
            length(generator$fixture_payloads) == member_count &&
            is.list(generator$fixture_truth)
        if (valid) {
            lapply(generator$fixture_payloads, \(binding) {
                metabPublicationValidateBinding(binding, "Fixture payload")
            })
            metabPublicationValidateBinding(generator$fixture_truth, "Fixture truth")
        }
    } else {
        valid <- valid && is.numeric(generator$chunk_rows) &&
            generator$chunk_rows == as.integer(generator$chunk_rows) &&
            generator$chunk_rows > 0L && !length(generator$fixture_payloads) &&
            is.null(generator$fixture_truth)
    }
    if (!valid) metabPublicationAbort("metabolomics generator differs")
    invisible(generator)
}

metabPublicationValidateRng <- function(rng, workload_class) {
    publicationRequireNames(rng, c(
        "kind", "normal_kind", "sample_kind", "seed_family_id", "seed",
        "streams"
    ), "Metabolomics RNG")
    fixture <- identical(workload_class, "fixture_correctness")
    valid <- identical(rng$kind, "L'Ecuyer-CMRG") &&
        identical(rng$normal_kind, "Inversion") &&
        identical(rng$sample_kind, "Rejection") &&
        publicationScalarString(rng$seed_family_id)
    if (fixture) {
        valid <- valid && is.null(rng$seed) && !length(rng$streams)
    } else {
        valid <- valid && is.numeric(rng$seed) && rng$seed == as.integer(rng$seed) &&
            length(rng$streams) >= 7L && !anyDuplicated(names(rng$streams))
    }
    if (!valid) metabPublicationAbort("metabolomics RNG differs")
    invisible(rng)
}

metabPublicationValidateClaimScope <- function(scope, capability, workload_class) {
    publicationRequireNames(scope, c(
        "evidence_class", "verified_stages", "scientific_authority",
        "performance_authority", "cross_project_authority",
        "promotion_authority", "limitations"
    ), "Metabolomics claim scope")
    fixture <- identical(workload_class, "fixture_correctness")
    performance <- workload_class %in% c(
        "representative",
        "operational_heavy"
    )
    expected_class <- switch(
        workload_class,
        fixture_correctness = "fixture_correctness",
        stress = "stress_characterization",
        "public_generated"
    )
    valid <- scope$evidence_class %in% c(
        "fixture_correctness", "public_generated", "stress_characterization"
    ) && identical(scope$evidence_class, expected_class) &&
        identical(scope$verified_stages, capability$verified_stages) &&
        identical(scope$scientific_authority, fixture) &&
        identical(scope$performance_authority, performance) &&
        !isTRUE(scope$cross_project_authority) &&
        !isTRUE(scope$promotion_authority) && length(scope$limitations) > 0L
    if (!valid) metabPublicationAbort("metabolomics claim scope differs")
    invisible(scope)
}

metabPublicationValidateTruthContract <- function(value, contract) {
    publicationRequireNames(value, c(
        "schema_id", "oracle_method", "validated_domains",
        "independent_of_package_reader", "support_widening_allowed"
    ), "Metabolomics truth contract")
    fixture <- identical(contract$workload_class, "fixture_correctness")
    expected_method <- if (fixture) {
        "direct_raw_table_arithmetic_and_reviewed_e2e_authority"
    } else {
        "independent_online_generated_truth"
    }
    required_domains <- c(
        "payload", "mapping", "assay_identity", "dimensions", "quantities",
        "missingness", "effects", "design"
    )
    valid <- identical(value$schema_id, "metabolomics.truth.v1") &&
        identical(value$oracle_method, expected_method) &&
        setequal(unlist(value$validated_domains), required_domains) &&
        !anyDuplicated(unlist(value$validated_domains)) &&
        isTRUE(value$independent_of_package_reader) &&
        !isTRUE(value$support_widening_allowed)
    if (!valid) metabPublicationAbort("metabolomics truth contract differs")
    invisible(value)
}

metabPublicationValidateExecution <- function(value, contract) {
    publicationRequireNames(value, c(
        "preparation_processes", "process_isolation", "candidate_loaded",
        "historical_pilot_required", "swap_allowed", "timeout_seconds",
        "expected_peak_bytes"
    ), "Metabolomics execution contract")
    expected_peak <- as.numeric(contract$dimensions$quantity_count * 16)
    valid <- identical(value$preparation_processes, 2L) &&
        identical(value$process_isolation, "fresh_vanilla_R") &&
        !isTRUE(value$candidate_loaded) && identical(
            value$historical_pilot_required,
            identical(contract$workload_class, "operational_heavy")
        ) && !isTRUE(value$swap_allowed) &&
        is.numeric(value$timeout_seconds) && value$timeout_seconds >= 60 &&
        is.numeric(value$expected_peak_bytes) &&
        value$expected_peak_bytes >= expected_peak
    if (!valid) metabPublicationAbort("metabolomics execution contract differs")
    invisible(value)
}

metabPublicationValidatePrivacy <- function(value, contract) {
    publicationRequireNames(value, c(
        "classification", "private_source_opened", "payload_tracked",
        "direct_identifiers_present", "cross_omic_values_reused"
    ), "Metabolomics privacy contract")
    expected <- if (identical(
        contract$workload_class,
        "fixture_correctness"
    )) {
        "public_fixture"
    } else {
        "public_generated"
    }
    valid <- identical(value$classification, expected) &&
        !isTRUE(value$private_source_opened) && !isTRUE(value$payload_tracked) &&
        !isTRUE(value$direct_identifiers_present) &&
        !isTRUE(value$cross_omic_values_reused)
    if (!valid) metabPublicationAbort("metabolomics privacy contract differs")
    invisible(value)
}

metabPublicationValidateWorkloadAuthorities <- function(contract) {
    mappings <- publicationReadJson(contract$mapping_authority$path)
    boundaries <- publicationReadJson(contract$support_boundary$path)
    sources <- publicationReadJson(contract$source_authority$path)
    splits <- publicationReadJson(contract$split_authority$path)
    metabPublicationValidateMappingAuthority(mappings)
    metabPublicationValidateSupportBoundaries(boundaries)
    metabPublicationValidateSources(sources)
    metabPublicationValidateSplits(splits)
    route <- contract$capability$input_format
    profile_id <- contract$assay_profile$profile_id
    mapping <- Filter(
        \(value) identical(value$route, route) &&
            identical(value$profile_id, profile_id),
        mappings$mappings
    )
    boundary <- Filter(\(value) identical(value$route, route), boundaries$routes)
    assignment <- Filter(
        \(value) identical(value$assignment_id, contract$workload_id),
        splits$assignments
    )
    source <- Filter(\(value) identical(value$route, route), sources$decisions)
    valid <- length(mapping) == 1L && length(boundary) == 1L &&
        length(assignment) == 1L && length(source) == 1L &&
        identical(mapping[[1L]]$mapping_id, contract$contract_mapping_id) &&
        identical(assignment[[1L]]$seed_family_id, contract$rng$seed_family_id) &&
        identical(assignment[[1L]]$seed, contract$rng$seed) &&
        identical(source[[1L]]$claim_scope, "project_specific_nonpromotional")
    if (!valid) metabPublicationAbort("metabolomics workload authority differs")
    invisible(contract)
}

metabPublicationValidateWorkload <- function(
    contract,
    allow_test_contract = FALSE,
    validate_authorities = TRUE
) {
    publicationRequireNames(
        contract,
        metabPublicationWorkloadFields(),
        "Metabolomics workload"
    )
    valid <- identical(
        contract$schema,
        "multischolar.omics_publication_metabolomics_workload"
    ) && identical(contract$schema_version, .METAB_PUBLICATION_VERSION) &&
        identical(contract$owner_ticket_id, .METAB_PUBLICATION_OWNER) &&
        contract$status %in% c("frozen", "test_only") &&
        (identical(contract$status, "frozen") || isTRUE(allow_test_contract)) &&
        publicationScalarString(contract$contract_id) &&
        publicationScalarString(contract$workload_id) &&
        contract$workload_class %in% metabPublicationClasses() &&
        publicationScalarString(contract$model_profile_id) &&
        !isTRUE(contract$publication_authority)
    if (!valid) metabPublicationAbort("metabolomics workload header differs")
    metabPublicationValidateCapability(contract$capability)
    metabPublicationValidateAssayProfile(
        contract$assay_profile,
        contract$capability$input_format
    )
    metabPublicationValidateDimensions(
        contract$dimensions,
        contract$assay_profile
    )
    if (identical(contract$status, "frozen")) metabPublicationValidateScale(contract)
    for (field in c(
        "parameter_authority", "source_authority", "split_authority",
        "mapping_authority", "support_boundary"
    )) {
        metabPublicationValidateBinding(contract[[field]], field)
    }
    metabPublicationValidateGenerator(contract$generator, contract)
    metabPublicationValidateRng(contract$rng, contract$workload_class)
    metabPublicationValidateTruthContract(contract$truth_contract, contract)
    metabPublicationValidateExecution(contract$execution, contract)
    metabPublicationValidatePrivacy(contract$privacy, contract)
    metabPublicationValidateClaimScope(
        contract$claim_scope,
        contract$capability,
        contract$workload_class
    )
    publicationRequireNames(contract$expected_digests, c(
        "payload_set_sha256", "truth_sha256"
    ), "Metabolomics expected digests")
    lapply(contract$expected_digests, metabPublicationRequireDigest, "expected")
    expected_id <- metabPublicationWorkloadId(
        contract$capability$input_format,
        contract$assay_profile$profile_id,
        contract$workload_class
    )
    if (!identical(contract$workload_id, expected_id)) {
        metabPublicationAbort("metabolomics workload id differs")
    }
    if (isTRUE(validate_authorities)) {
        metabPublicationValidateWorkloadAuthorities(contract)
    }
    invisible(contract)
}

metabPublicationContractBasis <- function(contract) {
    basis <- contract
    basis$expected_digests <- list(
        payload_set_sha256 = strrep("0", 64L),
        truth_sha256 = strrep("0", 64L)
    )
    publicationObjectDigest(basis)
}
