.publicationLipidRepoPath <- function(...) {
    file.path(
        normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
        ...
    )
}

for (.publication_lipid_source in c(
    "omics_workload_contract.R",
    "omics_publication_protocol.R",
    "omics_publication_lipidomics_contracts.R",
    "omics_publication_lipidomics_model.R",
    "omics_publication_lipidomics_serializers.R",
    "omics_publication_lipidomics_truth.R",
    "omics_publication_lipidomics_authorities.R",
    "omics_publication_lipidomics_sources.R",
    "omics_publication_lipidomics_governance.R",
    "omics_publication_lipidomics_negative.R",
    "omics_publication_lipidomics_fixtures.R",
    "omics_publication_lipidomics_freeze.R",
    "omics_publication_lipidomics_pilot.R",
    "omics_publication_lipidomics_manifest.R",
    "omics_publication_lipidomics_cleanup.R",
    "omics_publication_workload_lipidomics.R"
)) {
    sys.source(
        .publicationLipidRepoPath(
            "tools",
            "profiling",
            .publication_lipid_source
        ),
        envir = environment()
    )
}
rm(.publication_lipid_source)

.publicationLipidRecord <- function(name) {
    publicationReadJson(.publicationLipidRepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "lipidomics",
        name
    ))
}

.publicationLipidBinding <- function(path) {
    list(path = path, sha256 = publicationFileDigest(path))
}

.publicationLipidContract <- function(route, profile_id, seed = 310001L) {
    profile <- lipidPublicationAssayProfiles()[[profile_id]]
    assays <- unlist(profile$assays, use.names = FALSE)
    counts <- stats::setNames(as.list(rep(12L, length(assays))), assays)
    member_count <- if (startsWith(profile_id, "mixed_")) 2L else 1L
    parameter_path <- paste0(
        "tests/testdata/omics-performance/lipidomics/parameters-v1.json"
    )
    mappings <- "tests/testdata/omics-performance/lipidomics/mapping-authorities-v1.json"
    boundaries <- "tests/testdata/omics-performance/lipidomics/support-boundaries-v1.json"
    source_files <- c(
        "tools/profiling/omics_publication_lipidomics_contracts.R",
        "tools/profiling/omics_publication_lipidomics_model.R",
        "tools/profiling/omics_publication_lipidomics_serializers.R",
        "tools/profiling/omics_publication_lipidomics_truth.R",
        "tools/profiling/omics_publication_workload_lipidomics.R"
    )
    streams <- as.list(stats::setNames(
        seed + seq_len(8L) * 100L,
        c(
            "hierarchy", "feature_offsets", "class_residual", "batch",
            "residual", "mcar", "mar", "mnar"
        )
    ))
    list(
        schema = "multischolar.omics_publication_lipidomics_workload",
        schema_version = "1.0.0",
        contract_id = paste("test", route, profile_id, seed, sep = "."),
        owner_ticket_id = "OMICS-ART-067",
        status = "test_only",
        workload_id = lipidPublicationWorkloadId(
            route,
            profile_id,
            "representative"
        ),
        workload_class = "representative",
        capability = lipidPublicationExpectedCapability(route),
        assay_profile = list(
            profile_id = profile_id,
            assays = profile$assays,
            payload_mode = if (member_count == 2L) "bundle" else "single",
            member_count = member_count
        ),
        dimensions = list(
            aggregate_feature_count = as.integer(sum(unlist(counts))),
            assay_feature_counts = counts,
            sample_count = 8L,
            assay_count = as.integer(length(assays)),
            quantity_count = as.numeric(sum(unlist(counts)) * 8L),
            payload_member_count = member_count
        ),
        model_profile_id = "lipidomics.declared.test.v1",
        parameter_authority = .publicationLipidBinding(parameter_path),
        source_authority = .publicationLipidBinding(
            "tests/testdata/omics-performance/lipidomics/sources-v1.json"
        ),
        split_authority = .publicationLipidBinding(
            "tests/testdata/omics-performance/lipidomics/splits-v1.json"
        ),
        mapping_authority = .publicationLipidBinding(mappings),
        contract_mapping_id = paste(
            "lipidomics.mapping",
            route,
            profile_id,
            sep = "."
        ),
        support_boundary = .publicationLipidBinding(boundaries),
        generator = list(
            mode = "generated",
            source_bindings = lapply(as.list(source_files), .publicationLipidBinding),
            chunk_rows = 100L,
            output_members = lipidPublicationOutputMembers(route, profile_id),
            truth_filename = paste0(route, "-", profile_id, "-truth.json"),
            fixture_payloads = list(),
            fixture_truth = NULL
        ),
        rng = list(
            kind = "L'Ecuyer-CMRG",
            normal_kind = "Inversion",
            sample_kind = "Rejection",
            seed_family_id = "generated.holdout.300000-399999.v1",
            seed = seed,
            streams = streams
        ),
        truth_contract = list(
            schema_id = "lipidomics.truth.v1",
            oracle_method = "independent_online_generated_truth",
            validated_domains = as.list(c(
                "payload", "mapping", "assay_identity", "dimensions",
                "quantities", "missingness", "effects", "design",
                "lipid_classes", "isomer_like_relations"
            )),
            independent_of_package_reader = TRUE,
            support_widening_allowed = FALSE
        ),
        execution = list(
            preparation_processes = 2L,
            process_isolation = "fresh_vanilla_R",
            candidate_loaded = FALSE,
            historical_pilot_required = FALSE,
            swap_allowed = FALSE,
            timeout_seconds = 600,
            expected_peak_bytes = as.numeric(sum(unlist(counts)) * 8L * 16L)
        ),
        privacy = list(
            classification = "public_generated",
            private_source_opened = FALSE,
            payload_tracked = FALSE,
            direct_identifiers_present = FALSE,
            cross_omic_values_reused = FALSE
        ),
        claim_scope = list(
            evidence_class = "public_generated",
            verified_stages = lipidPublicationCapabilities()[[route]]$verified_stages,
            scientific_authority = FALSE,
            performance_authority = TRUE,
            cross_project_authority = FALSE,
            vendor_authority = FALSE,
            gc_vendor_authority = FALSE,
            three_file_workflow_authority = FALSE,
            promotion_authority = FALSE,
            limitations = list("Synthetic route-scoped truth only")
        ),
        expected_digests = list(
            payload_set_sha256 = strrep("a", 64L),
            truth_sha256 = strrep("b", 64L)
        ),
        publication_authority = FALSE
    )
}

test_that("lipidomics route support tiers are exact and non-aliased", {
    boundaries <- .publicationLipidRecord("support-boundaries-v1.json")
    mappings <- .publicationLipidRecord("mapping-authorities-v1.json")
    exclusions <- .publicationLipidRecord("exclusions-v1.json")
    expect_silent(lipidPublicationValidateSupportBoundaries(boundaries))
    expect_silent(lipidPublicationValidateMappingAuthority(mappings))
    expect_silent(lipidPublicationValidateExclusions(exclusions))
    expect_length(mappings$mappings, 15L)

    msdial <- Filter(\(mapping) identical(mapping$route, "msdial"), mappings$mappings)
    custom <- Filter(\(mapping) identical(mapping$route, "custom"), mappings$mappings)
    lipidsearch <- Filter(
        \(mapping) identical(mapping$route, "lipidsearch"),
        mappings$mappings
    )
    expect_true(all(vapply(msdial, \(mapping) {
        identical(mapping$verified_stages, list("import", "design")) &&
            identical(mapping$mapping_source, "msdial_schema_autodetection")
    }, logical(1))))
    expect_true(all(vapply(custom, \(mapping) {
        identical(mapping$mapping_source, "explicit_user_mapping_contract")
    }, logical(1))))
    expect_true(all(vapply(lipidsearch, \(mapping) {
        identical(
            mapping$mapping_source,
            "lipidsearch_schema_autodetection"
        ) && length(mapping$verified_stages) == 7L
    }, logical(1))))

    widened <- publicationGovernanceCopy(boundaries)
    widened$routes[[2L]]$verified_stages <- list("import", "design", "qc")
    expect_error(
        lipidPublicationValidateSupportBoundaries(widened),
        class = "multischolar_publication_lipidomics_contract_error"
    )

    aliased <- publicationGovernanceCopy(mappings)
    aliased$mappings[[1L]]$mapping_source <- "msdial_schema_autodetection"
    expect_error(
        lipidPublicationValidateMappingAuthority(aliased),
        class = "multischolar_publication_lipidomics_contract_error"
    )

    fallback <- publicationGovernanceCopy(mappings)
    fallback$mappings[[5L]]$fallback_allowed <- TRUE
    expect_error(
        lipidPublicationValidateMappingAuthority(fallback),
        class = "multischolar_publication_lipidomics_contract_error"
    )
})

test_that("mixed identities require exact two-file bundles", {
    authority <- .publicationLipidRecord("bundle-authorities-v1.json")
    expect_silent(lipidPublicationValidateBundles(authority))
    expect_length(authority$bundles, 24L)
    bundle <- authority$bundles[[1L]]
    expect_silent(lipidPublicationValidateBundle(bundle))

    missing <- publicationGovernanceCopy(bundle)
    missing$member_assays <- missing$member_assays[-1L]
    missing$member_count <- 1L
    expect_error(
        lipidPublicationValidateBundle(missing),
        class = "multischolar_publication_lipidomics_contract_error"
    )

    substituted <- publicationGovernanceCopy(bundle)
    substituted$single_table_substitution_allowed <- TRUE
    expect_error(
        lipidPublicationValidateBundle(substituted),
        class = "multischolar_publication_lipidomics_contract_error"
    )
})

test_that("detection-only lipidomics formats have no reader authority", {
    exclusions <- .publicationLipidRecord("exclusions-v1.json")
    expect_setequal(
        vapply(exclusions$formats, `[[`, character(1), "format"),
        c("progenesis", "xcms", "compound_discoverer", "unknown")
    )
    invoked <- publicationGovernanceCopy(exclusions)
    invoked$reader_invocation_allowed <- TRUE
    expect_error(
        lipidPublicationValidateExclusions(invoked),
        class = "multischolar_publication_lipidomics_contract_error"
    )
})

test_that("lipidomics workload identities preserve route assay and class", {
    ids <- unlist(lapply(names(lipidPublicationCapabilities()), \(route) {
        unlist(lapply(names(lipidPublicationAssayProfiles()), \(profile) {
            vapply(lipidPublicationClasses(), \(workload_class) {
                lipidPublicationWorkloadId(route, profile, workload_class)
            }, character(1))
        }), use.names = FALSE)
    }), use.names = FALSE)
    expect_length(ids, 60L)
    expect_identical(anyDuplicated(ids), 0L)
    expect_error(
        lipidPublicationWorkloadId(
            "progenesis",
            "mixed_lc",
            "representative"
        ),
        class = "multischolar_publication_lipidomics_contract_error"
    )
})

test_that("assay-aware latent model is deterministic and mask-attributable", {
    parameters <- .publicationLipidRecord("parameters-v1.json")
    expect_silent(lipidPublicationValidateParameters(parameters))
    expect_length(parameters$parameters, 36L)
    streams <- list(
        hierarchy = 310101L,
        feature_offsets = 310201L,
        class_residual = 310301L,
        batch = 310401L,
        residual = 310501L,
        mcar = 310601L,
        mar = 310701L,
        mnar = 310801L
    )
    for (profile_id in names(lipidPublicationAssayProfiles())) {
        assays <- unlist(
            lipidPublicationAssayProfiles()[[profile_id]]$assays,
            use.names = FALSE
        )
        counts <- stats::setNames(as.list(rep(12L, length(assays))), assays)
        first <- lipidPublicationModelPlan(
            counts,
            sample_count = 8L,
            parameter_authority = parameters,
            streams = streams,
            chunk_rows = 100L
        )
        second <- lipidPublicationModelPlan(
            counts,
            sample_count = 8L,
            parameter_authority = parameters,
            streams = streams,
            chunk_rows = 100L
        )
        block_first <- lipidPublicationGenerateBlock(
            first,
            seq_len(nrow(first$features))
        )
        block_second <- lipidPublicationGenerateBlock(
            second,
            seq_len(nrow(second$features))
        )
        expect_identical(first$features, second$features, info = profile_id)
        expect_identical(block_first, block_second, info = profile_id)
        expect_gt(first$correlation$minimum_eigenvalue, 0)
        expect_setequal(unique(first$features$assay), assays)
        family_assays <- tapply(
            first$features$assay,
            first$features$composition_family_id,
            \(value) length(unique(value))
        )
        family_classes <- tapply(
            first$features$lipid_class,
            first$features$composition_family_id,
            \(value) length(unique(value))
        )
        expect_true(all(family_assays == 1L), info = profile_id)
        expect_true(all(family_classes == 1L), info = profile_id)
        expect_true(all(
            first$features$effect_log2[first$features$internal_standard] == 0
        ), info = profile_id)
        masks <- list(
            block_first$mcar_missing,
            block_first$mar_missing,
            block_first$mnar_missing,
            block_first$censored_missing
        )
        overlap <- Reduce(`+`, lapply(masks, as.integer))
        expect_true(all(overlap <= 1L), info = profile_id)
        expect_true(all(
            is.na(block_first$values) == (overlap == 1L)
        ), info = profile_id)
        expect_true(all(block_first$residual_sigma > 0), info = profile_id)
    }
})

test_that("isomer-like relationships remain synthetic and domain-contained", {
    parameters <- .publicationLipidRecord("parameters-v1.json")
    counts <- list(LCMS_Pos = 600L, LCMS_Neg = 600L)
    plan <- lipidPublicationModelPlan(
        counts,
        sample_count = 8L,
        parameter_authority = parameters,
        streams = lipidPublicationFixtureStreams("custom", "mixed_lc"),
        chunk_rows = 2000L
    )
    paired <- !is.na(plan$features$isomer_pair_id)
    expect_gt(sum(paired), 0L)
    expect_true(all(startsWith(plan$features$lipid_class, "SYN_")))
    expect_identical(anyDuplicated(plan$features$feature_id), 0L)
    pair_assays <- tapply(
        plan$features$assay[paired],
        plan$features$isomer_pair_id[paired],
        function(value) length(unique(value))
    )
    pair_families <- tapply(
        plan$features$composition_family_id[paired],
        plan$features$isomer_pair_id[paired],
        function(value) length(unique(value))
    )
    expect_true(all(pair_assays == 1L))
    expect_true(all(pair_families == 1L))
    coordinates <- lipidPublicationCoordinates(plan$features, plan$parameters)
    pair_mz_span <- tapply(
        coordinates$mz[paired],
        plan$features$isomer_pair_id[paired],
        function(value) diff(range(value)) / mean(value) * 1e6
    )
    pair_rt_span <- tapply(
        coordinates$retention_time[paired],
        plan$features$isomer_pair_id[paired],
        function(value) diff(range(value))
    )
    expect_true(all(
        pair_mz_span <= plan$parameters$isomer_mass_delta_ppm_max + 1e-3
    ))
    expect_true(all(
        pair_rt_span <= plan$parameters$isomer_rt_delta_minutes_max + 1e-6
    ))
})

test_that("parameter provenance rejects unbound empirical vocabulary", {
    parameters <- .publicationLipidRecord("parameters-v1.json")
    expect_setequal(
        vapply(parameters$parameters, `[[`, character(1), "parameter_id"),
        names(lipidPublicationParameterConsumers())
    )
    changed <- publicationGovernanceCopy(parameters)
    changed$parameters[[1L]]$allowed_claim_vocabulary <- list("realistic")
    expect_error(
        lipidPublicationValidateParameters(changed),
        class = "multischolar_publication_lipidomics_contract_error"
    )

    missing <- publicationGovernanceCopy(parameters)
    missing$parameters <- missing$parameters[-1L]
    expect_error(
        lipidPublicationValidateParameters(missing),
        class = "multischolar_publication_lipidomics_contract_error"
    )
})

test_that("route serializers preserve mapping assay and lipid truth", {
    for (route in names(lipidPublicationCapabilities())) {
        for (profile_id in names(lipidPublicationAssayProfiles())) {
            contract <- .publicationLipidContract(route, profile_id)
            expect_silent(lipidPublicationValidateWorkload(
                contract,
                allow_test_contract = TRUE,
                validate_authorities = FALSE
            ))
            first <- lipidPublicationPrepareGenerated(
                contract,
                tempfile("lipid-publication-a-"),
                verify_expected = FALSE,
                allow_test_contract = TRUE
            )
            second <- lipidPublicationPrepareGenerated(
                contract,
                tempfile("lipid-publication-b-"),
                verify_expected = FALSE,
                allow_test_contract = TRUE
            )
            expect_identical(
                first$payload$payload_set_sha256,
                second$payload$payload_set_sha256
            )
            expect_identical(first$truth$sha256, second$truth$sha256)
            expect_identical(first$receipt$sha256, second$receipt$sha256)
            expected_members <- if (startsWith(profile_id, "mixed_")) 2L else 1L
            expect_identical(first$payload$member_count, expected_members)
            result <- lipidPublicationRunImported(
                contract,
                file.path(dirname(first$truth$path), "payload"),
                first$truth$path
            )
            expect_true(result$workflow_evidence$truth_valid)
            expect_identical(
                result$workflow_evidence$verified_stages,
                lipidPublicationCapabilities()[[route]]$verified_stages
            )
            expect_false(result$workflow_evidence$promotion_authority)
            expected_mapping <- switch(
                route,
                lipidsearch = "lipidsearch_schema_autodetection",
                msdial = "msdial_schema_autodetection",
                custom = "explicit_user_mapping_contract"
            )
            expect_identical(
                result$workflow_evidence$import$mapping_source,
                expected_mapping
            )
            if (!identical(route, "custom")) {
                expect_identical(
                    result$workflow_evidence$import$reader_extra_sample_column_count,
                    0L
                )
            } else {
                expect_gt(
                    result$workflow_evidence$import$reader_extra_sample_column_count,
                    0L
                )
            }
            expect_gt(
                result$workflow_evidence$import$lipid_class_count,
                0L
            )
            expect_gt(
                result$workflow_evidence$import$composition_family_count,
                0L
            )
        }
    }
})

test_that("lipidomics sources reject real-project and cross-omic substitution", {
    sources <- .publicationLipidRecord("sources-v1.json")
    expect_silent(lipidPublicationValidateSources(sources))
    expect_true(all(vapply(sources$decisions, \(decision) {
        !isTRUE(decision$promotion_eligible)
    }, logical(1))))
    expect_setequal(
        names(sources$scale_receipt),
        c("row_count", "column_count", "byte_size", "salted_source_fingerprint")
    )

    leaked <- publicationGovernanceCopy(sources)
    leaked$scale_receipt$unique_run_count <- 12L
    expect_error(
        lipidPublicationValidateSources(leaked),
        class = "multischolar_publication_schema_error"
    )

    promoted <- publicationGovernanceCopy(sources)
    promoted$decisions[[1L]]$promotion_eligible <- TRUE
    expect_error(
        lipidPublicationValidateSources(promoted),
        class = "multischolar_publication_lipidomics_contract_error"
    )
})

test_that("frozen splits are complete and pilot-holdout disjoint", {
    splits <- .publicationLipidRecord("splits-v1.json")
    bundles <- .publicationLipidRecord("bundle-authorities-v1.json")
    governance <- .publicationLipidRecord("governance-successor-v1.json")
    expect_silent(lipidPublicationValidateSplits(splits))
    expect_silent(lipidPublicationValidateBundles(bundles))
    expect_silent(lipidPublicationValidateGovernanceSuccessor(governance))
    expect_length(splits$pilot_calibration_assignments, 3L)
    expect_length(splits$assignments, 60L)
    expect_length(bundles$bundles, 24L)

    leaked <- publicationGovernanceCopy(splits)
    leaked$assignments[[1L]]$source_id <-
        leaked$pilot_calibration_assignments[[1L]]$source_id
    expect_error(
        lipidPublicationValidateSplits(leaked),
        class = "multischolar_publication_lipidomics_contract_error"
    )

    aliased <- publicationGovernanceCopy(bundles)
    aliased$single_table_substitution_allowed <- TRUE
    expect_error(
        lipidPublicationValidateBundles(aliased),
        class = "multischolar_publication_lipidomics_contract_error"
    )

    stale <- publicationGovernanceCopy(governance)
    stale$sources$sha256 <- strrep("0", 64L)
    expect_error(
        lipidPublicationValidateGovernanceSuccessor(stale),
        class = "multischolar_publication_lipidomics_contract_error"
    )
})

test_that("test-only contracts cannot masquerade as frozen workloads", {
    contract <- .publicationLipidContract("custom", "lcms_pos")
    expect_error(
        lipidPublicationValidateWorkload(contract),
        class = "multischolar_publication_lipidomics_contract_error"
    )
    frozen <- publicationGovernanceCopy(contract)
    frozen$status <- "frozen"
    expect_error(
        lipidPublicationValidateWorkload(
            frozen,
            validate_authorities = FALSE
        ),
        class = "multischolar_publication_lipidomics_contract_error"
    )

    widened <- publicationGovernanceCopy(contract)
    widened$truth_contract$support_widening_allowed <- TRUE
    expect_error(
        lipidPublicationValidateWorkload(
            widened,
            allow_test_contract = TRUE,
            validate_authorities = FALSE
        ),
        class = "multischolar_publication_lipidomics_contract_error"
    )

    private <- publicationGovernanceCopy(contract)
    private$privacy$private_source_opened <- TRUE
    expect_error(
        lipidPublicationValidateWorkload(
            private,
            allow_test_contract = TRUE,
            validate_authorities = FALSE
        ),
        class = "multischolar_publication_lipidomics_contract_error"
    )
})

test_that("negative contracts preserve rejection and invocation boundaries", {
    authority <- publicationReadJson(.publicationLipidRepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "lipidomics",
        "negative",
        "contracts-v1.json"
    ))
    expect_silent(lipidPublicationValidateNegativeAuthority(authority))
    expect_length(authority$cases, 40L)
    for (case in authority$cases) {
        suffix <- if (startsWith(case$mutation, "mixed_")) "" else ".tsv"
        path <- paste0(tempfile("lipidomics-negative-"), suffix)
        prepared <- lipidPublicationPrepareNegative(case, path)
        expect_true(file.exists(prepared$path) || dir.exists(prepared$path))
        evidence <- lipidPublicationEvaluateNegative(case, prepared$path)
        expect_identical(evidence$observed_outcome, case$expected_outcome)
        expect_identical(
            evidence$reader_invocation_count,
            as.integer(case$reader_invocation_expected)
        )
        expect_identical(evidence$fallback_invocation_count, 0L)
    }
})

test_that("reviewed fixture contracts replay with direct arithmetic truth", {
    paths <- sort(list.files(
        .publicationLipidRepoPath(
            "tests",
            "testdata",
            "omics-performance",
            "lipidomics",
            "workloads"
        ),
        pattern = "^fixture-correctness-.*[.]json$",
        full.names = TRUE
    ))
    expect_length(paths, 15L)
    for (path in paths) {
        contract <- lipidPublicationReadContract(path)
        first <- lipidPublicationPrepareFixture(
            contract,
            tempfile("lipidomics-fixture-a-")
        )
        second <- lipidPublicationPrepareFixture(
            contract,
            tempfile("lipidomics-fixture-b-")
        )
        expect_identical(
            first$payload$payload_set_sha256,
            second$payload$payload_set_sha256
        )
        expect_identical(first$truth$sha256, second$truth$sha256)
        expect_identical(first$receipt$sha256, second$receipt$sha256)
        result <- lipidPublicationRunImported(
            contract,
            file.path(dirname(first$truth$path), "payload"),
            first$truth$path
        )
        expect_true(result$workflow_evidence$truth_valid)
        expect_true(result$workflow_evidence$differential_direction$valid)
        expect_false(result$workflow_evidence$promotion_authority)
        if (contract$assay_profile$profile_id %in% c(
            "gcms",
            "mixed_lc_gcms"
        )) {
            expect_false(contract$claim_scope$gc_vendor_authority)
            expect_false(contract$claim_scope$vendor_authority)
        }
    }
})

test_that("heavy pilot contracts are historical calibration only", {
    contracts <- lapply(names(lipidPublicationCapabilities()), function(route) {
        lipidPublicationBuildPilotContract(
            route,
            aggregate_feature_count = 80000L,
            sample_count = 48L,
            grid_id = "floor"
        )
    })
    expect_identical(
        vapply(
            contracts,
            function(contract) contract$dimensions$payload_member_count,
            integer(1)
        ),
        c(2L, 2L, 2L)
    )
    expect_true(all(vapply(contracts, function(contract) {
        !isTRUE(contract$claim_scope$performance_authority) &&
            !isTRUE(contract$claim_scope$promotion_authority)
    }, logical(1))))

    promoted <- publicationGovernanceCopy(contracts[[1L]])
    promoted$claim_scope$performance_authority <- TRUE
    expect_error(
        lipidPublicationValidatePilotContract(promoted),
        class = "multischolar_publication_lipidomics_contract_error"
    )

    qualified <- lipidPublicationPilotQualification(list(
        status = "passed",
        publication_certifiable = TRUE,
        timed_out = FALSE,
        safety_aborted = FALSE,
        phase_evidence = list(valid = TRUE),
        safety_evidence = publicationTestSafetyEvidence(),
        metrics = list(
            peak_charged_memory_bytes = 4 * 1024^3,
            elapsed_seconds = 10
        )
    ))
    expect_true(qualified$qualified)
    aborted <- list(
        status = "failed",
        safety_aborted = TRUE,
        timed_out = FALSE
    )
    expect_identical(
        lipidPublicationPilotStatus(aborted, list(qualified = FALSE)),
        "safety_aborted_no_dimension_decision"
    )

    heavy <- lipidPublicationBuildHeavyContract(
        "msdial",
        "mixed_lc",
        aggregate_feature_count = 80000L,
        sample_count = 48L,
        qualification_id = "test_qualification"
    )
    expect_silent(lipidPublicationValidateWorkload(heavy))
    expect_identical(heavy$dimensions$payload_member_count, 2L)
    expect_true(heavy$claim_scope$performance_authority)

    incomplete_manifest <- list(
        schema = "multischolar.omics_publication_lipidomics_manifest"
    )
    expect_error(
        lipidPublicationValidateFinalManifest(incomplete_manifest),
        class = "multischolar_publication_schema_error"
    )
})

test_that("lipidomics payload cleanup is archived and symlink safe", {
    sandbox <- withr::local_tempdir(pattern = "lipidomics-cleanup-")
    withr::local_dir(sandbox)
    root <- file.path(".omics-publication-workloads", "lipidomics")
    generated <- file.path(root, "generated", "attempt-1")
    protected <- file.path(root, "authorities")
    dir.create(file.path(generated, "logs"), recursive = TRUE)
    dir.create(file.path(generated, "payload"), recursive = TRUE)
    dir.create(protected, recursive = TRUE)
    writeLines("log", file.path(generated, "logs", "run.log"))
    writeLines("payload", file.path(generated, "payload", "large.tsv"))
    writeLines("authority", file.path(protected, "receipt.json"))
    expect_true(file.symlink(
        normalizePath(file.path(protected, "receipt.json")),
        file.path(generated, "payload", "protected-link")
    ))
    plan <- list(
        schema = "multischolar.omics_publication_lipidomics_cleanup_plan",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-067",
        status = "approved_for_execution",
        root = root,
        archive_root = "cleanup-evidence/test",
        protected_paths = list("authorities"),
        removals = list(list(
            path = "generated/attempt-1",
            reason = "test generated payload"
        )),
        publication_authority = FALSE
    )
    dry_run <- lipidPublicationRunCleanup(plan, execute = FALSE)
    expect_identical(dry_run$status, "dry_run")
    expect_true(dir.exists(generated))
    expect_identical(dry_run$removals[[1L]]$before$symlink_count, 1L)

    result <- lipidPublicationRunCleanup(plan, execute = TRUE)
    expect_identical(result$status, "passed")
    expect_false(dir.exists(generated))
    expect_true(file.exists(file.path(protected, "receipt.json")))
    expect_true(file.exists(file.path(
        root,
        "cleanup-evidence/test/generated__attempt-1/retained/logs/run.log"
    )))
})

test_that("proteomics and metabolomics adapters reject lipidomics identity", {
    contract <- omicsWorkloadReadContract(.publicationLipidRepoPath(
        "tests",
        "testdata",
        "omics-parity",
        "lipidomics",
        "workloads",
        "lcms-pos-public-ci-v1.json"
    ))
    for (name in c("proteomics", "metabolomics")) {
        path <- .publicationLipidRepoPath(
            "tools",
            "profiling",
            paste0("omics_workload_", name, ".R")
        )
        environment <- new.env(parent = globalenv())
        sys.source(path, envir = environment)
        adapter <- environment$omicsWorkloadAdapter()
        changed <- contract
        changed$adapter <- list(
            adapter_id = adapter$adapter_id,
            adapter_version = adapter$adapter_version,
            source_sha256 = omicsWorkloadFileDigest(path)
        )
        expect_error(
            omicsWorkloadLoadAdapter(path, changed),
            class = "omics_workload_binding_error",
            info = name
        )
    }
})
