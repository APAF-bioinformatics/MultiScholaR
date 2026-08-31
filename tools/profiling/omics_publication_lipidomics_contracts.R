.LIPID_PUBLICATION_OWNER <- "OMICS-ART-067"
.LIPID_PUBLICATION_VERSION <- "1.0.0"

lipidPublicationAbort <- function(message, class = "contract_error") {
    publicationAbort(
        message,
        paste0("multischolar_publication_lipidomics_", class)
    )
}

lipidPublicationCapabilities <- function() {
    list(
        lipidsearch = list(
            capability_id = "lipidomics.lipidsearch.lipid.standard.v1",
            input_format = "lipidsearch",
            support_tier = "full_workflow",
            verified_stages = list(
                "import", "design", "qc", "normalization_ruv", "session",
                "differential_abundance", "summary_report_export"
            )
        ),
        msdial = list(
            capability_id = "lipidomics.msdial.lipid.standard.v1",
            input_format = "msdial",
            support_tier = "reader_boundary",
            verified_stages = list("import", "design")
        ),
        custom = list(
            capability_id = "lipidomics.custom.lipid.standard.v1",
            input_format = "custom",
            support_tier = "reader_boundary",
            verified_stages = list("import", "design")
        )
    )
}

lipidPublicationAssayProfiles <- function() {
    list(
        lcms_pos = list(assays = list("LCMS_Pos")),
        lcms_neg = list(assays = list("LCMS_Neg")),
        gcms = list(assays = list("GCMS")),
        mixed_lc = list(assays = list("LCMS_Pos", "LCMS_Neg")),
        mixed_lc_gcms = list(assays = list("LCMS_Pos", "GCMS"))
    )
}

lipidPublicationClasses <- function() {
    c(
        "fixture_correctness", "representative", "operational_heavy", "stress"
    )
}

lipidPublicationRequireDigest <- function(value, label) {
    if (!publicationScalarString(value) || !grepl("^[0-9a-f]{64}$", value)) {
        lipidPublicationAbort(paste(label, "must be a SHA-256 digest"))
    }
    invisible(value)
}

lipidPublicationScalarFlag <- function(value) {
    is.logical(value) && length(value) == 1L && !is.na(value)
}

lipidPublicationValidateBinding <- function(
    binding,
    label,
    require_current = TRUE
) {
    publicationRequireNames(binding, c("path", "sha256"), label)
    lipidPublicationRequireDigest(binding$sha256, paste(label, "sha256"))
    if (isTRUE(require_current)) {
        path <- publicationPath(binding$path)
        if (!file.exists(path) || !identical(
            publicationFileDigest(path),
            binding$sha256
        )) {
            lipidPublicationAbort(paste(label, "differs"))
        }
    }
    invisible(binding)
}

lipidPublicationExpectedCapability <- function(route) {
    expected <- lipidPublicationCapabilities()[[route]]
    if (is.null(expected)) lipidPublicationAbort("lipidomics route is unsupported")
    c(
        list(
            omic_type = "lipidomics",
            data_level = "lipid",
            acquisition_mode = "not_recorded"
        ),
        expected
    )
}

lipidPublicationValidateCapability <- function(capability) {
    publicationRequireNames(capability, c(
        "omic_type", "data_level", "acquisition_mode", "capability_id",
        "input_format", "support_tier", "verified_stages"
    ), "Lipidomics capability")
    if (!identical(
        capability,
        lipidPublicationExpectedCapability(capability$input_format)
    )) {
        lipidPublicationAbort("lipidomics capability differs")
    }
    invisible(capability)
}

lipidPublicationValidateAssayProfile <- function(profile, route) {
    publicationRequireNames(profile, c(
        "profile_id", "assays", "payload_mode", "member_count"
    ), "Lipidomics assay profile")
    expected <- lipidPublicationAssayProfiles()[[profile$profile_id]]
    mixed <- startsWith(profile$profile_id, "mixed_")
    expected_mode <- if (mixed) "bundle" else "single"
    expected_count <- if (mixed) 2L else 1L
    valid <- !is.null(expected) && identical(profile$assays, expected$assays) &&
        identical(profile$payload_mode, expected_mode) &&
        identical(profile$member_count, expected_count)
    if (!valid) lipidPublicationAbort("lipidomics assay profile differs")
    invisible(profile)
}

lipidPublicationValidateDimensions <- function(dimensions, profile) {
    publicationRequireNames(dimensions, c(
        "aggregate_feature_count", "assay_feature_counts", "sample_count",
        "assay_count", "quantity_count", "payload_member_count"
    ), "Lipidomics dimensions")
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
    if (!valid) lipidPublicationAbort("lipidomics dimensions differ")
    invisible(dimensions)
}

lipidPublicationExpectedScale <- function(profile_id, workload_class) {
    single <- !startsWith(profile_id, "mixed_")
    switch(
        workload_class,
        fixture_correctness = list(
            minimum_features = if (single) 3L else 6L,
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
        lipidPublicationAbort("lipidomics workload class is unsupported")
    )
}

lipidPublicationValidateScale <- function(contract) {
    scale <- lipidPublicationExpectedScale(
        contract$assay_profile$profile_id,
        contract$workload_class
    )
    dimensions <- contract$dimensions
    valid <- dimensions$aggregate_feature_count >= scale$minimum_features
    if (!is.null(scale$sample_count)) {
        valid <- valid && dimensions$sample_count >= scale$sample_count
    }
    if (startsWith(contract$assay_profile$profile_id, "mixed_")) {
        counts <- unlist(dimensions$assay_feature_counts, use.names = FALSE)
        valid <- valid && length(counts) == 2L && all(counts > 0L)
    }
    if (!valid) lipidPublicationAbort("lipidomics workload scale is below floor")
    invisible(contract)
}

lipidPublicationWorkloadId <- function(route, profile, workload_class) {
    capability <- lipidPublicationCapabilities()[[route]]
    if (is.null(capability) || is.null(lipidPublicationAssayProfiles()[[profile]]) ||
        !workload_class %in% lipidPublicationClasses()) {
        lipidPublicationAbort("lipidomics workload identity is unsupported")
    }
    paste(capability$capability_id, profile, workload_class, sep = ".")
}

lipidPublicationValidateMapping <- function(mapping) {
    publicationRequireNames(mapping, c(
        "mapping_id", "route", "profile_id", "requested_format", "observed_format",
        "mapping_source", "fallback_allowed", "bundle_member_mappings",
        "verified_stages", "promotion_authority"
    ), "Lipidomics mapping authority")
    capability <- lipidPublicationCapabilities()[[mapping$route]]
    profile <- lipidPublicationAssayProfiles()[[mapping$profile_id]]
    expected_source <- switch(
        mapping$route,
        lipidsearch = "lipidsearch_schema_autodetection",
        msdial = "msdial_schema_autodetection",
        custom = "explicit_user_mapping_contract"
    )
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
    if (!valid) lipidPublicationAbort("lipidomics mapping differs")
    invisible(mapping)
}

lipidPublicationValidateMappingAuthority <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "mapping_authority_id", "owner_ticket_id",
        "status", "mappings", "identity_alias_allowed", "unknown_policy",
        "publication_authority"
    ), "Lipidomics mapping authority")
    lapply(record$mappings, lipidPublicationValidateMapping)
    ids <- vapply(record$mappings, `[[`, character(1), "mapping_id")
    pairs <- vapply(record$mappings, \(mapping) {
        paste(mapping$route, mapping$profile_id, sep = "::")
    }, character(1))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_lipidomics_mapping_authority"
    ) && identical(record$schema_version, .LIPID_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .LIPID_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        length(record$mappings) == 15L && !anyDuplicated(ids) &&
        !anyDuplicated(pairs) && !isTRUE(record$identity_alias_allowed) &&
        identical(record$unknown_policy, "reject") &&
        !isTRUE(record$publication_authority)
    if (!valid) lipidPublicationAbort("mapping authority differs")
    invisible(record)
}

lipidPublicationValidateSupportBoundaries <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "boundary_id", "owner_ticket_id",
        "coverage_binding", "capability_binding", "routes",
        "reachable_stage_claim_allowed", "unknown_policy",
        "publication_authority"
    ), "Lipidomics support boundaries")
    lipidPublicationValidateBinding(record$coverage_binding, "Coverage")
    lipidPublicationValidateBinding(record$capability_binding, "Capabilities")
    route_ids <- vapply(record$routes, `[[`, character(1), "route")
    valid_routes <- all(vapply(record$routes, \(route) {
        expected <- lipidPublicationCapabilities()[[route$route]]
        setequal(names(route), c(
            "route", "capability_id", "support_tier", "verified_stages",
            "full_workflow_claim_allowed", "generated_scope"
        )) && !is.null(expected) &&
            identical(route$capability_id, expected$capability_id) &&
            identical(route$support_tier, expected$support_tier) &&
            identical(route$verified_stages, expected$verified_stages) &&
            identical(
                route$full_workflow_claim_allowed,
                identical(route$route, "lipidsearch")
            ) && publicationScalarString(route$generated_scope)
    }, logical(1)))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_lipidomics_support_boundaries"
    ) && identical(record$schema_version, .LIPID_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .LIPID_PUBLICATION_OWNER) &&
        setequal(route_ids, names(lipidPublicationCapabilities())) &&
        !anyDuplicated(route_ids) && valid_routes &&
        !isTRUE(record$reachable_stage_claim_allowed) &&
        identical(record$unknown_policy, "reject") &&
        !isTRUE(record$publication_authority)
    if (!valid) lipidPublicationAbort("support boundary differs")
    invisible(record)
}

lipidPublicationValidateExclusions <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "exclusions_id", "owner_ticket_id",
        "capability_binding", "formats", "reader_invocation_allowed",
        "fallback_allowed", "promotion_authority", "publication_authority"
    ), "Lipidomics exclusions")
    lipidPublicationValidateBinding(record$capability_binding, "Capabilities")
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
        "multischolar.omics_publication_lipidomics_exclusions"
    ) && identical(record$schema_version, .LIPID_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .LIPID_PUBLICATION_OWNER) &&
        setequal(ids, expected) && !anyDuplicated(ids) && valid_formats &&
        !isTRUE(record$reader_invocation_allowed) &&
        !isTRUE(record$fallback_allowed) &&
        !isTRUE(record$promotion_authority) &&
        !isTRUE(record$publication_authority)
    if (!valid) lipidPublicationAbort("lipidomics exclusions differ")
    invisible(record)
}

lipidPublicationValidateBundle <- function(bundle) {
    publicationRequireNames(bundle, c(
        "bundle_id", "route", "profile_id", "member_assays", "member_count",
        "member_schema_ids", "single_table_substitution_allowed",
        "bundle_digest"
    ), "Lipidomics bundle")
    expected <- lipidPublicationAssayProfiles()[[bundle$profile_id]]
    valid <- bundle$route %in% names(lipidPublicationCapabilities()) &&
        startsWith(bundle$profile_id, "mixed_") && !is.null(expected) &&
        identical(bundle$member_assays, expected$assays) &&
        identical(bundle$member_count, 2L) &&
        length(bundle$member_schema_ids) == 2L &&
        !anyDuplicated(unlist(bundle$member_schema_ids, use.names = FALSE)) &&
        !isTRUE(bundle$single_table_substitution_allowed)
    lipidPublicationRequireDigest(bundle$bundle_digest, "Lipidomics bundle")
    basis <- bundle
    basis$bundle_digest <- strrep("0", 64L)
    valid <- valid && identical(
        bundle$bundle_digest,
        publicationObjectDigest(basis)
    )
    if (!valid) lipidPublicationAbort("lipidomics bundle differs")
    invisible(bundle)
}

lipidPublicationValidateBundles <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "bundles_id", "owner_ticket_id",
        "status", "mapping_authority", "bundles", "member_order_semantic",
        "missing_or_duplicate_member_policy", "single_table_substitution_allowed",
        "publication_authority"
    ), "Lipidomics bundle authority")
    lipidPublicationValidateBinding(
        record$mapping_authority,
        "Lipidomics mapping authority"
    )
    lapply(record$bundles, lipidPublicationValidateBundle)
    ids <- vapply(record$bundles, `[[`, character(1), "bundle_id")
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_lipidomics_bundles"
    ) && identical(record$schema_version, .LIPID_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .LIPID_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        length(record$bundles) == 24L &&
        !anyDuplicated(ids) && !isTRUE(record$member_order_semantic) &&
        identical(record$missing_or_duplicate_member_policy, "reject") &&
        !isTRUE(record$single_table_substitution_allowed) &&
        !isTRUE(record$publication_authority)
    if (!valid) lipidPublicationAbort("lipidomics bundle authority differs")
    invisible(record)
}

lipidPublicationParameterOrigins <- function() {
    c(
        "methodological_design", "protocol_authority",
        "public_source_calibrated", "private_scale_calibrated"
    )
}

lipidPublicationParameterConsumers <- function() {
    c(
        treatment_fraction = "model.design",
        batch_count = "model.design",
        technical_replicate_size = "model.design",
        sample_correlation = "model.correlation",
        technical_replicate_correlation = "model.correlation",
        lipid_class_correlation = "model.correlation",
        effect_fraction = "model.effects",
        effect_log2 = "model.effects_and_truth",
        lcms_pos_mean_log2 = "model.assay_means",
        lcms_neg_mean_log2 = "model.assay_means",
        gcms_mean_log2 = "model.assay_means",
        base_log2_sd = "model.class_baseline",
        class_offset_log2_sd = "model.class_offset",
        family_size = "model.composition_families",
        family_offset_log2_sd = "model.composition_family_offset",
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
        isomer_like_fraction = "model.isomer_like_pairs",
        isomer_mass_delta_ppm_max = "model.isomer_like_pairs",
        isomer_rt_delta_minutes_max = "model.isomer_like_pairs",
        duplicate_fraction = "negative.duplicate_identifiers",
        annotation_width = "serializer.annotations",
        minimum_effect_sign_agreement = "truth.effect_direction",
        effect_median_margin_fraction = "truth.effect_direction"
    )
}

lipidPublicationValidateParameterDomain <- function(value, domain) {
    publicationRequireNames(
        domain,
        c("type", "minimum", "maximum"),
        "Lipidomics parameter domain"
    )
    valid <- domain$type %in% c("numeric_interval", "integer_interval") &&
        is.numeric(value) && length(value) == 1L && is.finite(value) &&
        value >= domain$minimum && value <= domain$maximum
    if (identical(domain$type, "integer_interval")) {
        valid <- valid && value == as.integer(value)
    }
    if (!valid) lipidPublicationAbort("parameter value is outside its domain")
    invisible(value)
}

lipidPublicationValidateParameter <- function(parameter) {
    publicationRequireNames(parameter, c(
        "parameter_id", "value", "unit", "domain", "origin",
        "source_binding", "applicable_routes", "applicable_profiles",
        "applicable_classes", "allowed_claim_vocabulary", "limitations"
    ), "Lipidomics parameter")
    lipidPublicationValidateParameterDomain(parameter$value, parameter$domain)
    design <- identical(parameter$origin, "methodological_design")
    calibrated <- parameter$origin %in% setdiff(
        lipidPublicationParameterOrigins(),
        "methodological_design"
    )
    valid <- publicationScalarString(parameter$parameter_id) &&
        publicationScalarString(parameter$unit) &&
        parameter$origin %in% lipidPublicationParameterOrigins() &&
        length(parameter$applicable_routes) > 0L &&
        all(parameter$applicable_routes %in% names(
            lipidPublicationCapabilities()
        )) && length(parameter$applicable_profiles) > 0L &&
        all(parameter$applicable_profiles %in% names(
            lipidPublicationAssayProfiles()
        )) && length(parameter$applicable_classes) > 0L &&
        all(parameter$applicable_classes %in% lipidPublicationClasses()) &&
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
            lipidPublicationValidateBinding(
                parameter$source_binding,
                "Lipidomics parameter source"
            )
        }
    }
    if (!valid) lipidPublicationAbort("lipidomics parameter differs")
    invisible(parameter)
}

lipidPublicationValidateParameters <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "parameters_id", "owner_ticket_id",
        "status", "origin_classes", "forbidden_unbound_vocabulary",
        "parameters", "publication_authority"
    ), "Lipidomics parameters")
    lapply(record$parameters, lipidPublicationValidateParameter)
    ids <- vapply(record$parameters, `[[`, character(1), "parameter_id")
    consumers <- lipidPublicationParameterConsumers()
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_lipidomics_parameters"
    ) && identical(record$schema_version, .LIPID_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .LIPID_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        identical(
            record$origin_classes,
            as.list(lipidPublicationParameterOrigins())
        ) && identical(
            record$forbidden_unbound_vocabulary,
            as.list(c("realistic", "empirical", "instrument_derived"))
        ) && !anyDuplicated(ids) && setequal(ids, names(consumers)) &&
        all(nzchar(unname(consumers))) && !isTRUE(record$publication_authority)
    if (!valid) lipidPublicationAbort("lipidomics parameters differ")
    invisible(record)
}

lipidPublicationParameterValues <- function(record) {
    lipidPublicationValidateParameters(record)
    values <- lapply(record$parameters, `[[`, "value")
    names(values) <- vapply(
        record$parameters,
        `[[`,
        character(1),
        "parameter_id"
    )
    values
}

lipidPublicationWorkloadFields <- function() {
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

lipidPublicationValidateGenerator <- function(generator, contract) {
    publicationRequireNames(generator, c(
        "mode", "source_bindings", "chunk_rows", "output_members",
        "truth_filename", "fixture_payloads", "fixture_truth"
    ), "Lipidomics generator")
    fixture <- identical(contract$workload_class, "fixture_correctness")
    expected_mode <- if (fixture) "fixture_replay" else "generated"
    member_count <- contract$dimensions$payload_member_count
    valid <- identical(generator$mode, expected_mode) &&
        length(generator$source_bindings) > 0L &&
        all(vapply(generator$source_bindings, \(binding) {
            tryCatch({
                lipidPublicationValidateBinding(binding, "Generator source")
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
                lipidPublicationValidateBinding(binding, "Fixture payload")
            })
            lipidPublicationValidateBinding(generator$fixture_truth, "Fixture truth")
        }
    } else {
        valid <- valid && is.numeric(generator$chunk_rows) &&
            generator$chunk_rows == as.integer(generator$chunk_rows) &&
            generator$chunk_rows > 0L && !length(generator$fixture_payloads) &&
            is.null(generator$fixture_truth)
    }
    if (!valid) lipidPublicationAbort("lipidomics generator differs")
    invisible(generator)
}

lipidPublicationValidateRng <- function(rng, workload_class) {
    publicationRequireNames(rng, c(
        "kind", "normal_kind", "sample_kind", "seed_family_id", "seed",
        "streams"
    ), "Lipidomics RNG")
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
    if (!valid) lipidPublicationAbort("lipidomics RNG differs")
    invisible(rng)
}

lipidPublicationValidateClaimScope <- function(
    scope,
    capability,
    profile_id,
    workload_class
) {
    publicationRequireNames(scope, c(
        "evidence_class", "verified_stages", "scientific_authority",
        "performance_authority", "cross_project_authority",
        "vendor_authority", "gc_vendor_authority",
        "three_file_workflow_authority", "promotion_authority", "limitations"
    ), "Lipidomics claim scope")
    fixture <- identical(workload_class, "fixture_correctness")
    route <- capability$input_format
    reviewed <- fixture && identical(route, "lipidsearch")
    gc_smoke <- reviewed && profile_id %in% c("gcms", "mixed_lc_gcms")
    performance <- workload_class %in% c(
        "representative",
        "operational_heavy"
    )
    expected_class <- if (fixture) {
        if (!reviewed) {
            "generated_reader_characterization"
        } else if (gc_smoke) {
            "application_assay_name_smoke"
        } else {
            "reviewed_fixture_correctness"
        }
    } else if (identical(workload_class, "stress")) {
        "stress_characterization"
    } else {
        "public_generated"
    }
    valid <- scope$evidence_class %in% c(
        "reviewed_fixture_correctness", "application_assay_name_smoke",
        "generated_reader_characterization", "public_generated",
        "stress_characterization"
    ) && identical(scope$evidence_class, expected_class) &&
        identical(scope$verified_stages, capability$verified_stages) &&
        identical(scope$scientific_authority, reviewed) &&
        identical(scope$performance_authority, performance) &&
        !isTRUE(scope$cross_project_authority) &&
        identical(scope$vendor_authority, reviewed && !gc_smoke) &&
        !isTRUE(scope$gc_vendor_authority) &&
        !isTRUE(scope$three_file_workflow_authority) &&
        !isTRUE(scope$promotion_authority) && length(scope$limitations) > 0L
    if (!valid) lipidPublicationAbort("lipidomics claim scope differs")
    invisible(scope)
}

lipidPublicationValidateTruthContract <- function(value, contract) {
    publicationRequireNames(value, c(
        "schema_id", "oracle_method", "validated_domains",
        "independent_of_package_reader", "support_widening_allowed"
    ), "Lipidomics truth contract")
    fixture <- identical(contract$workload_class, "fixture_correctness")
    route <- contract$capability$input_format
    profile_id <- contract$assay_profile$profile_id
    expected_method <- if (fixture && identical(route, "lipidsearch")) {
        if (profile_id %in% c("gcms", "mixed_lc_gcms")) {
            "direct_arithmetic_application_assay_name_smoke"
        } else {
            "direct_raw_table_arithmetic_and_reviewed_e2e_authority"
        }
    } else if (fixture) {
        "generated_reader_schema_truth"
    } else {
        "independent_online_generated_truth"
    }
    required_domains <- c(
        "payload", "mapping", "assay_identity", "dimensions", "quantities",
        "missingness", "effects", "design", "lipid_classes",
        "isomer_like_relations"
    )
    valid <- identical(value$schema_id, "lipidomics.truth.v1") &&
        identical(value$oracle_method, expected_method) &&
        setequal(unlist(value$validated_domains), required_domains) &&
        !anyDuplicated(unlist(value$validated_domains)) &&
        isTRUE(value$independent_of_package_reader) &&
        !isTRUE(value$support_widening_allowed)
    if (!valid) lipidPublicationAbort("lipidomics truth contract differs")
    invisible(value)
}

lipidPublicationValidateExecution <- function(value, contract) {
    publicationRequireNames(value, c(
        "preparation_processes", "process_isolation", "candidate_loaded",
        "historical_pilot_required", "swap_allowed", "timeout_seconds",
        "expected_peak_bytes"
    ), "Lipidomics execution contract")
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
    if (!valid) lipidPublicationAbort("lipidomics execution contract differs")
    invisible(value)
}

lipidPublicationValidatePrivacy <- function(value, contract) {
    publicationRequireNames(value, c(
        "classification", "private_source_opened", "payload_tracked",
        "direct_identifiers_present", "cross_omic_values_reused"
    ), "Lipidomics privacy contract")
    expected <- if (identical(
        contract$workload_class,
        "fixture_correctness"
    ) && identical(contract$capability$input_format, "lipidsearch")) {
        "public_fixture"
    } else {
        "public_generated"
    }
    valid <- identical(value$classification, expected) &&
        !isTRUE(value$private_source_opened) && !isTRUE(value$payload_tracked) &&
        !isTRUE(value$direct_identifiers_present) &&
        !isTRUE(value$cross_omic_values_reused)
    if (!valid) lipidPublicationAbort("lipidomics privacy contract differs")
    invisible(value)
}

lipidPublicationValidateWorkloadAuthorities <- function(contract) {
    mappings <- publicationReadJson(contract$mapping_authority$path)
    boundaries <- publicationReadJson(contract$support_boundary$path)
    sources <- publicationReadJson(contract$source_authority$path)
    splits <- publicationReadJson(contract$split_authority$path)
    lipidPublicationValidateMappingAuthority(mappings)
    lipidPublicationValidateSupportBoundaries(boundaries)
    lipidPublicationValidateSources(sources)
    lipidPublicationValidateSplits(splits)
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
    if (!valid) lipidPublicationAbort("lipidomics workload authority differs")
    invisible(contract)
}

lipidPublicationValidateWorkload <- function(
    contract,
    allow_test_contract = FALSE,
    validate_authorities = TRUE
) {
    publicationRequireNames(
        contract,
        lipidPublicationWorkloadFields(),
        "Lipidomics workload"
    )
    valid <- identical(
        contract$schema,
        "multischolar.omics_publication_lipidomics_workload"
    ) && identical(contract$schema_version, .LIPID_PUBLICATION_VERSION) &&
        identical(contract$owner_ticket_id, .LIPID_PUBLICATION_OWNER) &&
        contract$status %in% c("frozen", "test_only") &&
        (identical(contract$status, "frozen") || isTRUE(allow_test_contract)) &&
        publicationScalarString(contract$contract_id) &&
        publicationScalarString(contract$workload_id) &&
        contract$workload_class %in% lipidPublicationClasses() &&
        publicationScalarString(contract$model_profile_id) &&
        !isTRUE(contract$publication_authority)
    if (!valid) lipidPublicationAbort("lipidomics workload header differs")
    lipidPublicationValidateCapability(contract$capability)
    lipidPublicationValidateAssayProfile(
        contract$assay_profile,
        contract$capability$input_format
    )
    lipidPublicationValidateDimensions(
        contract$dimensions,
        contract$assay_profile
    )
    if (identical(contract$status, "frozen")) lipidPublicationValidateScale(contract)
    for (field in c(
        "parameter_authority", "source_authority", "split_authority",
        "mapping_authority", "support_boundary"
    )) {
        lipidPublicationValidateBinding(contract[[field]], field)
    }
    lipidPublicationValidateGenerator(contract$generator, contract)
    lipidPublicationValidateRng(contract$rng, contract$workload_class)
    lipidPublicationValidateTruthContract(contract$truth_contract, contract)
    lipidPublicationValidateExecution(contract$execution, contract)
    lipidPublicationValidatePrivacy(contract$privacy, contract)
    lipidPublicationValidateClaimScope(
        contract$claim_scope,
        contract$capability,
        contract$assay_profile$profile_id,
        contract$workload_class
    )
    publicationRequireNames(contract$expected_digests, c(
        "payload_set_sha256", "truth_sha256"
    ), "Lipidomics expected digests")
    lapply(contract$expected_digests, lipidPublicationRequireDigest, "expected")
    expected_id <- lipidPublicationWorkloadId(
        contract$capability$input_format,
        contract$assay_profile$profile_id,
        contract$workload_class
    )
    if (!identical(contract$workload_id, expected_id)) {
        lipidPublicationAbort("lipidomics workload id differs")
    }
    if (isTRUE(validate_authorities)) {
        lipidPublicationValidateWorkloadAuthorities(contract)
    }
    invisible(contract)
}

lipidPublicationContractBasis <- function(contract) {
    basis <- contract
    basis$expected_digests <- list(
        payload_set_sha256 = strrep("0", 64L),
        truth_sha256 = strrep("0", 64L)
    )
    publicationObjectDigest(basis)
}
