lipidPublicationParameterRecord <- function(
    parameter_id,
    value,
    unit,
    minimum,
    maximum,
    integer = FALSE
) {
    list(
        parameter_id = parameter_id,
        value = value,
        unit = unit,
        domain = list(
            type = if (integer) "integer_interval" else "numeric_interval",
            minimum = minimum,
            maximum = maximum
        ),
        origin = "methodological_design",
        source_binding = NULL,
        applicable_routes = as.list(names(lipidPublicationCapabilities())),
        applicable_profiles = as.list(names(lipidPublicationAssayProfiles())),
        applicable_classes = as.list(lipidPublicationClasses()),
        allowed_claim_vocabulary = list(
            "declared_synthetic",
            "mechanistic_simulation"
        ),
        limitations = list(
            paste(
                "Declared simulation parameter; no empirical instrument,",
                "cohort, or chemical identity claim."
            )
        )
    )
}

lipidPublicationBuildParameters <- function() {
    specs <- list(
        c("treatment_fraction", 0.5, "fraction_of_samples", 0.1, 0.9, FALSE),
        c("batch_count", 4, "batches", 1, 32, TRUE),
        c("technical_replicate_size", 2, "samples_per_group", 1, 8, TRUE),
        c("sample_correlation", 0.15, "pearson_correlation", 0, 0.9, FALSE),
        c("technical_replicate_correlation", 0.65, "pearson_correlation", 0, 0.98, FALSE),
        c("lipid_class_correlation", 0.3, "pearson_correlation", 0, 0.9, FALSE),
        c("effect_fraction", 0.05, "fraction_per_direction", 0.001, 0.24, FALSE),
        c("effect_log2", 1.2, "log2_fold_change", 0.1, 5, FALSE),
        c("lcms_pos_mean_log2", 18, "log2_intensity", 1, 40, FALSE),
        c("lcms_neg_mean_log2", 17.5, "log2_intensity", 1, 40, FALSE),
        c("gcms_mean_log2", 16, "log2_intensity", 1, 40, FALSE),
        c("base_log2_sd", 1.25, "log2_intensity_sd", 0.01, 8, FALSE),
        c("class_offset_log2_sd", 0.35, "log2_intensity_sd", 0, 4, FALSE),
        c("family_size", 4, "features_per_family", 2, 32, TRUE),
        c("family_offset_log2_sd", 0.25, "log2_intensity_sd", 0, 4, FALSE),
        c("feature_offset_log2_sd", 0.12, "log2_intensity_sd", 0, 4, FALSE),
        c("batch_log2_sd", 0.18, "log2_intensity_sd", 0, 4, FALSE),
        c("residual_sigma_floor", 0.1, "log2_intensity_sd", 0.001, 4, FALSE),
        c("residual_sigma_slope", 0.08, "sd_per_log2_unit", 0, 2, FALSE),
        c("heteroscedastic_reference_log2", 17, "log2_intensity", 1, 40, FALSE),
        c("mcar_probability", 0.015, "probability", 0, 0.5, FALSE),
        c("mar_logit_intercept", -4, "log_odds", -20, 20, FALSE),
        c("mar_batch_log_odds", 0.7, "log_odds", -10, 10, FALSE),
        c("detection_midpoint_log2", 14, "log2_intensity", 1, 40, FALSE),
        c("detection_slope_log2", 1, "log2_intensity", 0.01, 10, FALSE),
        c("lcms_pos_censor_log2", 15.5, "log2_intensity", 1, 40, FALSE),
        c("lcms_neg_censor_log2", 15, "log2_intensity", 1, 40, FALSE),
        c("gcms_censor_log2", 13.5, "log2_intensity", 1, 40, FALSE),
        c("internal_standard_fraction", 0.05, "fraction_of_features", 0.001, 0.5, FALSE),
        c("isomer_like_fraction", 0.02, "fraction_of_features", 0, 0.5, FALSE),
        c("isomer_mass_delta_ppm_max", 5, "parts_per_million", 0.01, 50, FALSE),
        c("isomer_rt_delta_minutes_max", 0.05, "minutes", 0.0001, 2, FALSE),
        c("duplicate_fraction", 0.01, "fraction_of_features", 0, 0.2, FALSE),
        c("annotation_width", 64, "characters", 8, 512, TRUE),
        c("minimum_effect_sign_agreement", 0.8, "fraction", 0.5, 1, FALSE),
        c("effect_median_margin_fraction", 0.25, "fraction", 0, 1, FALSE)
    )
    parameters <- lapply(specs, function(spec) {
        lipidPublicationParameterRecord(
            spec[[1L]],
            as.numeric(spec[[2L]]),
            spec[[3L]],
            as.numeric(spec[[4L]]),
            as.numeric(spec[[5L]]),
            identical(spec[[6L]], "TRUE")
        )
    })
    list(
        schema = "multischolar.omics_publication_lipidomics_parameters",
        schema_version = .LIPID_PUBLICATION_VERSION,
        parameters_id = paste0(
            "multischolar.omics_publication_lipidomics_parameters.2026-08-28.v1"
        ),
        owner_ticket_id = .LIPID_PUBLICATION_OWNER,
        status = "frozen_pre_candidate",
        origin_classes = as.list(lipidPublicationParameterOrigins()),
        forbidden_unbound_vocabulary = as.list(c(
            "realistic", "empirical", "instrument_derived"
        )),
        parameters = parameters,
        publication_authority = FALSE
    )
}

lipidPublicationMappingSource <- function(route) {
    switch(
        route,
        lipidsearch = "lipidsearch_schema_autodetection",
        msdial = "msdial_schema_autodetection",
        custom = "explicit_user_mapping_contract"
    )
}

lipidPublicationBuildMappings <- function() {
    mappings <- list()
    for (route in names(lipidPublicationCapabilities())) {
        for (profile_id in names(lipidPublicationAssayProfiles())) {
            assays <- unlist(
                lipidPublicationAssayProfiles()[[profile_id]]$assays,
                use.names = FALSE
            )
            mappings[[length(mappings) + 1L]] <- list(
                mapping_id = paste(
                    "lipidomics.mapping", route, profile_id, sep = "."
                ),
                route = route,
                profile_id = profile_id,
                requested_format = route,
                observed_format = route,
                mapping_source = lipidPublicationMappingSource(route),
                fallback_allowed = FALSE,
                bundle_member_mappings = lapply(as.list(assays), function(assay) {
                    list(
                        assay = assay,
                        mapping_id = paste(
                            route,
                            profile_id,
                            tolower(assay),
                            "v1",
                            sep = "."
                        )
                    )
                }),
                verified_stages =
                    lipidPublicationCapabilities()[[route]]$verified_stages,
                promotion_authority = FALSE
            )
        }
    }
    list(
        schema = "multischolar.omics_publication_lipidomics_mapping_authority",
        schema_version = .LIPID_PUBLICATION_VERSION,
        mapping_authority_id = paste0(
            "multischolar.omics_publication_lipidomics_mapping.2026-08-28.v1"
        ),
        owner_ticket_id = .LIPID_PUBLICATION_OWNER,
        status = "frozen_pre_candidate",
        mappings = mappings,
        identity_alias_allowed = FALSE,
        unknown_policy = "reject",
        publication_authority = FALSE
    )
}

lipidPublicationBuildSupportBoundaries <- function() {
    coverage_path <- "tests/testdata/omics-performance/coverage-v1.json"
    capability_path <- "tests/testdata/omics-capabilities.json"
    routes <- lapply(names(lipidPublicationCapabilities()), function(route) {
        capability <- lipidPublicationCapabilities()[[route]]
        list(
            route = route,
            capability_id = capability$capability_id,
            support_tier = capability$support_tier,
            verified_stages = capability$verified_stages,
            full_workflow_claim_allowed = identical(route, "lipidsearch"),
            generated_scope = if (identical(route, "lipidsearch")) {
                "generated_scale_truth_with_reviewed_fixture_science"
            } else {
                "exact_reader_schema_import_design_truth_only"
            }
        )
    })
    list(
        schema = "multischolar.omics_publication_lipidomics_support_boundaries",
        schema_version = .LIPID_PUBLICATION_VERSION,
        boundary_id = paste0(
            "multischolar.omics_publication_lipidomics_boundaries.2026-08-28.v1"
        ),
        owner_ticket_id = .LIPID_PUBLICATION_OWNER,
        coverage_binding = lipidPublicationAuthorityBinding(coverage_path),
        capability_binding = lipidPublicationAuthorityBinding(capability_path),
        routes = routes,
        reachable_stage_claim_allowed = FALSE,
        unknown_policy = "reject",
        publication_authority = FALSE
    )
}

lipidPublicationBuildExclusions <- function() {
    capability_path <- "tests/testdata/omics-capabilities.json"
    formats <- c("progenesis", "xcms", "compound_discoverer", "unknown")
    list(
        schema = "multischolar.omics_publication_lipidomics_exclusions",
        schema_version = .LIPID_PUBLICATION_VERSION,
        exclusions_id = paste0(
            "multischolar.omics_publication_lipidomics_exclusions.2026-08-28.v1"
        ),
        owner_ticket_id = .LIPID_PUBLICATION_OWNER,
        capability_binding = lipidPublicationAuthorityBinding(capability_path),
        formats = lapply(as.list(formats), function(format) {
            list(
                format = format,
                support_status = if (identical(format, "unknown")) {
                    "unknown_rejected"
                } else {
                    "detection_only"
                },
                capability_count = 0L,
                reason = "No dedicated reviewed reader and fixture authority."
            )
        }),
        reader_invocation_allowed = FALSE,
        fallback_allowed = FALSE,
        promotion_authority = FALSE,
        publication_authority = FALSE
    )
}

lipidPublicationBuildSources <- function() {
    projects_path <- "tests/testdata/omics-performance/projects-v1.json"
    private_path <- paste0(
        "tests/testdata/omics-performance/proteomics/private-envelope-v1.json"
    )
    private <- publicationReadJson(private_path)$report
    routes <- names(lipidPublicationCapabilities())
    list(
        schema = "multischolar.omics_publication_lipidomics_sources",
        schema_version = .LIPID_PUBLICATION_VERSION,
        sources_id = paste0(
            "multischolar.omics_publication_lipidomics_sources.2026-08-28.v1"
        ),
        owner_ticket_id = .LIPID_PUBLICATION_OWNER,
        status = "frozen_nonpromotional",
        projects_predecessor = lipidPublicationAuthorityBinding(projects_path),
        private_scale_predecessor = lipidPublicationAuthorityBinding(private_path),
        allowed_scale_fields = list(
            "row_count", "column_count", "byte_size",
            "salted_source_fingerprint"
        ),
        scale_receipt = list(
            row_count = private$row_count,
            column_count = private$column_count,
            byte_size = private$byte_size,
            salted_source_fingerprint = private$salted_source_fingerprint
        ),
        scale_mapping_policy =
            "proteomics_shape_to_bounded_lipidomics_dimensions_only_v1",
        minimum_real_projects = 3L,
        sources = list(),
        decisions = lapply(as.list(routes), function(route) {
            list(
                route = route,
                capability_id =
                    lipidPublicationCapabilities()[[route]]$capability_id,
                required_real_project_count = 3L,
                verified_real_project_count = 0L,
                claim_scope = "project_specific_nonpromotional",
                cross_project_source_ready = FALSE,
                promotion_eligible = FALSE
            )
        }),
        repository_searches = lapply(as.list(routes), function(route) {
            list(
                route = route,
                repository = "Metabolomics Workbench and MetaboLights",
                query = paste(route, "processed lipidomics exact result schema"),
                admissible_source_count = 0L,
                outcome = "no_admissible_format_exact_source_frozen"
            )
        }),
        generated_counts_as_real = FALSE,
        fixtures_count_as_real = FALSE,
        unknown_policy = "reject",
        publication_authority = FALSE
    )
}
