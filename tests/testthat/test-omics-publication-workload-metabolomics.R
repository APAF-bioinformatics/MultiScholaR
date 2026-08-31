.publicationMetabRepoPath <- function(...) {
    file.path(
        normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
        ...
    )
}

for (.publication_metab_source in c(
    "omics_workload_contract.R",
    "omics_publication_protocol.R",
    "omics_publication_metabolomics_contracts.R",
    "omics_publication_metabolomics_model.R",
    "omics_publication_metabolomics_serializers.R",
    "omics_publication_metabolomics_truth.R",
    "omics_publication_metabolomics_sources.R",
    "omics_publication_metabolomics_governance.R",
    "omics_publication_metabolomics_negative.R",
    "omics_publication_metabolomics_fixtures.R",
    "omics_publication_metabolomics_freeze.R",
    "omics_publication_metabolomics_pilot.R",
    "omics_publication_metabolomics_manifest.R",
    "omics_publication_metabolomics_cleanup.R",
    "omics_publication_workload_metabolomics.R"
)) {
    sys.source(
        .publicationMetabRepoPath(
            "tools",
            "profiling",
            .publication_metab_source
        ),
        envir = environment()
    )
}
rm(.publication_metab_source)

.publicationMetabRecord <- function(name) {
    publicationReadJson(.publicationMetabRepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "metabolomics",
        name
    ))
}

.publicationMetabBinding <- function(path) {
    list(path = path, sha256 = publicationFileDigest(path))
}

.publicationMetabContract <- function(route, profile_id, seed = 310001L) {
    profile <- metabPublicationAssayProfiles()[[profile_id]]
    assays <- unlist(profile$assays, use.names = FALSE)
    counts <- stats::setNames(as.list(rep(12L, length(assays))), assays)
    member_count <- if (identical(route, "msdial") &&
        identical(profile_id, "mixed")) 3L else 1L
    parameter_path <- paste0(
        "tests/testdata/omics-performance/metabolomics/parameters-v1.json"
    )
    mappings <- "tests/testdata/omics-performance/metabolomics/mapping-authorities-v1.json"
    boundaries <- "tests/testdata/omics-performance/metabolomics/support-boundaries-v1.json"
    source_files <- c(
        "tools/profiling/omics_publication_metabolomics_contracts.R",
        "tools/profiling/omics_publication_metabolomics_model.R",
        "tools/profiling/omics_publication_metabolomics_serializers.R",
        "tools/profiling/omics_publication_metabolomics_truth.R",
        "tools/profiling/omics_publication_workload_metabolomics.R"
    )
    streams <- as.list(stats::setNames(
        seed + seq_len(7L) * 100L,
        c("hierarchy", "feature_offsets", "batch", "residual", "mcar", "mar", "mnar")
    ))
    list(
        schema = "multischolar.omics_publication_metabolomics_workload",
        schema_version = "1.0.0",
        contract_id = paste("test", route, profile_id, seed, sep = "."),
        owner_ticket_id = "OMICS-ART-066",
        status = "test_only",
        workload_id = metabPublicationWorkloadId(
            route,
            profile_id,
            "representative"
        ),
        workload_class = "representative",
        capability = metabPublicationExpectedCapability(route),
        assay_profile = list(
            profile_id = profile_id,
            assays = profile$assays,
            payload_mode = if (member_count == 3L) "bundle" else "single",
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
        model_profile_id = "metabolomics.declared.test.v1",
        parameter_authority = .publicationMetabBinding(parameter_path),
        source_authority = .publicationMetabBinding(
            "tests/testdata/omics-performance/metabolomics/sources-v1.json"
        ),
        split_authority = .publicationMetabBinding(
            "tests/testdata/omics-performance/metabolomics/splits-v1.json"
        ),
        mapping_authority = .publicationMetabBinding(mappings),
        contract_mapping_id = paste(
            "metabolomics.mapping",
            route,
            profile_id,
            sep = "."
        ),
        support_boundary = .publicationMetabBinding(boundaries),
        generator = list(
            mode = "generated",
            source_bindings = lapply(as.list(source_files), .publicationMetabBinding),
            chunk_rows = 100L,
            output_members = metabPublicationOutputMembers(route, profile_id),
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
            schema_id = "metabolomics.truth.v1",
            oracle_method = "independent_online_generated_truth",
            validated_domains = as.list(c(
                "payload", "mapping", "assay_identity", "dimensions",
                "quantities", "missingness", "effects", "design"
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
            verified_stages = metabPublicationCapabilities()[[route]]$verified_stages,
            scientific_authority = FALSE,
            performance_authority = TRUE,
            cross_project_authority = FALSE,
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

test_that("custom and MS-DIAL support tiers are exact and non-aliased", {
    boundaries <- .publicationMetabRecord("support-boundaries-v1.json")
    mappings <- .publicationMetabRecord("mapping-authorities-v1.json")
    exclusions <- .publicationMetabRecord("exclusions-v1.json")
    expect_silent(metabPublicationValidateSupportBoundaries(boundaries))
    expect_silent(metabPublicationValidateMappingAuthority(mappings))
    expect_silent(metabPublicationValidateExclusions(exclusions))
    expect_length(mappings$mappings, 8L)

    msdial <- Filter(\(mapping) identical(mapping$route, "msdial"), mappings$mappings)
    custom <- Filter(\(mapping) identical(mapping$route, "custom"), mappings$mappings)
    expect_true(all(vapply(msdial, \(mapping) {
        identical(mapping$verified_stages, list("import", "design")) &&
            identical(mapping$mapping_source, "msdial_schema_autodetection")
    }, logical(1))))
    expect_true(all(vapply(custom, \(mapping) {
        identical(mapping$mapping_source, "explicit_user_mapping_contract")
    }, logical(1))))

    widened <- publicationGovernanceCopy(boundaries)
    widened$routes[[2L]]$verified_stages <- list("import", "design", "qc")
    expect_error(
        metabPublicationValidateSupportBoundaries(widened),
        class = "multischolar_publication_metabolomics_contract_error"
    )

    aliased <- publicationGovernanceCopy(mappings)
    aliased$mappings[[1L]]$mapping_source <- "msdial_schema_autodetection"
    expect_error(
        metabPublicationValidateMappingAuthority(aliased),
        class = "multischolar_publication_metabolomics_contract_error"
    )

    fallback <- publicationGovernanceCopy(mappings)
    fallback$mappings[[5L]]$fallback_allowed <- TRUE
    expect_error(
        metabPublicationValidateMappingAuthority(fallback),
        class = "multischolar_publication_metabolomics_contract_error"
    )
})

test_that("mixed MS-DIAL identity requires a three-file bundle", {
    bundle <- list(
        bundle_id = "metabolomics.msdial.mixed.fixture.v1",
        profile_id = "mixed",
        member_assays = list("LCMS_Pos", "LCMS_Neg", "GCMS"),
        member_count = 3L,
        member_schema_ids = list(
            "msdial.lcms_pos.v1",
            "msdial.lcms_neg.v1",
            "msdial.gcms.v1"
        ),
        single_custom_table_substitution_allowed = FALSE,
        bundle_digest = strrep("0", 64L)
    )
    bundle$bundle_digest <- publicationObjectDigest(bundle)
    expect_silent(metabPublicationValidateMsdialBundle(bundle))

    missing <- publicationGovernanceCopy(bundle)
    missing$member_assays <- list("LCMS_Pos", "GCMS")
    missing$member_count <- 2L
    expect_error(
        metabPublicationValidateMsdialBundle(missing),
        class = "multischolar_publication_metabolomics_contract_error"
    )

    substituted <- publicationGovernanceCopy(bundle)
    substituted$single_custom_table_substitution_allowed <- TRUE
    expect_error(
        metabPublicationValidateMsdialBundle(substituted),
        class = "multischolar_publication_metabolomics_contract_error"
    )
})

test_that("detection-only metabolomics formats have no reader authority", {
    exclusions <- .publicationMetabRecord("exclusions-v1.json")
    expect_setequal(
        vapply(exclusions$formats, `[[`, character(1), "format"),
        c("progenesis", "xcms", "compound_discoverer", "unknown")
    )
    invoked <- publicationGovernanceCopy(exclusions)
    invoked$reader_invocation_allowed <- TRUE
    expect_error(
        metabPublicationValidateExclusions(invoked),
        class = "multischolar_publication_metabolomics_contract_error"
    )
})

test_that("metabolomics workload identities preserve route assay and class", {
    ids <- unlist(lapply(names(metabPublicationCapabilities()), \(route) {
        unlist(lapply(names(metabPublicationAssayProfiles()), \(profile) {
            vapply(metabPublicationClasses(), \(workload_class) {
                metabPublicationWorkloadId(route, profile, workload_class)
            }, character(1))
        }), use.names = FALSE)
    }), use.names = FALSE)
    expect_length(ids, 32L)
    expect_identical(anyDuplicated(ids), 0L)
    expect_error(
        metabPublicationWorkloadId("progenesis", "mixed", "representative"),
        class = "multischolar_publication_metabolomics_contract_error"
    )
})

test_that("assay-aware latent model is deterministic and mask-attributable", {
    parameters <- .publicationMetabRecord("parameters-v1.json")
    expect_silent(metabPublicationValidateParameters(parameters))
    expect_length(parameters$parameters, 31L)
    streams <- list(
        hierarchy = 310101L,
        feature_offsets = 310201L,
        batch = 310301L,
        residual = 310401L,
        mcar = 310501L,
        mar = 310601L,
        mnar = 310701L
    )
    for (profile_id in names(metabPublicationAssayProfiles())) {
        assays <- unlist(
            metabPublicationAssayProfiles()[[profile_id]]$assays,
            use.names = FALSE
        )
        counts <- stats::setNames(as.list(rep(12L, length(assays))), assays)
        first <- metabPublicationModelPlan(
            counts,
            sample_count = 8L,
            parameter_authority = parameters,
            streams = streams,
            chunk_rows = 100L
        )
        second <- metabPublicationModelPlan(
            counts,
            sample_count = 8L,
            parameter_authority = parameters,
            streams = streams,
            chunk_rows = 100L
        )
        block_first <- metabPublicationGenerateBlock(
            first,
            seq_len(nrow(first$features))
        )
        block_second <- metabPublicationGenerateBlock(
            second,
            seq_len(nrow(second$features))
        )
        expect_identical(first$features, second$features, info = profile_id)
        expect_identical(block_first, block_second, info = profile_id)
        expect_gt(first$correlation$minimum_eigenvalue, 0)
        expect_setequal(unique(first$features$assay), assays)
        group_assays <- tapply(
            first$features$assay,
            first$features$correlated_group_id,
            \(value) length(unique(value))
        )
        expect_true(all(group_assays == 1L), info = profile_id)
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

test_that("parameter provenance rejects unbound empirical vocabulary", {
    parameters <- .publicationMetabRecord("parameters-v1.json")
    expect_setequal(
        vapply(parameters$parameters, `[[`, character(1), "parameter_id"),
        names(metabPublicationParameterConsumers())
    )
    changed <- publicationGovernanceCopy(parameters)
    changed$parameters[[1L]]$allowed_claim_vocabulary <- list("realistic")
    expect_error(
        metabPublicationValidateParameters(changed),
        class = "multischolar_publication_metabolomics_contract_error"
    )

    missing <- publicationGovernanceCopy(parameters)
    missing$parameters <- missing$parameters[-1L]
    expect_error(
        metabPublicationValidateParameters(missing),
        class = "multischolar_publication_metabolomics_contract_error"
    )
})

test_that("custom and MS-DIAL serializers preserve route and assay truth", {
    for (route in names(metabPublicationCapabilities())) {
        for (profile_id in names(metabPublicationAssayProfiles())) {
            contract <- .publicationMetabContract(route, profile_id)
            expect_silent(metabPublicationValidateWorkload(
                contract,
                allow_test_contract = TRUE,
                validate_authorities = FALSE
            ))
            first <- metabPublicationPrepareGenerated(
                contract,
                tempfile("metab-publication-a-"),
                verify_expected = FALSE,
                allow_test_contract = TRUE
            )
            second <- metabPublicationPrepareGenerated(
                contract,
                tempfile("metab-publication-b-"),
                verify_expected = FALSE,
                allow_test_contract = TRUE
            )
            expect_identical(
                first$payload$payload_set_sha256,
                second$payload$payload_set_sha256
            )
            expect_identical(first$truth$sha256, second$truth$sha256)
            expect_identical(first$receipt$sha256, second$receipt$sha256)
            expected_members <- if (identical(route, "msdial") &&
                identical(profile_id, "mixed")) 3L else 1L
            expect_identical(first$payload$member_count, expected_members)
            result <- metabPublicationRunImported(
                contract,
                file.path(dirname(first$truth$path), "payload"),
                first$truth$path
            )
            expect_true(result$workflow_evidence$truth_valid)
            expect_identical(
                result$workflow_evidence$verified_stages,
                metabPublicationCapabilities()[[route]]$verified_stages
            )
            expect_false(result$workflow_evidence$promotion_authority)
            expected_mapping <- if (identical(route, "custom")) {
                "explicit_user_mapping_contract"
            } else {
                "msdial_schema_autodetection"
            }
            expect_identical(
                result$workflow_evidence$import$mapping_source,
                expected_mapping
            )
            if (identical(route, "msdial")) {
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
        }
    }
})

test_that("metabolomics sources reject real-project and cross-omic substitution", {
    sources <- .publicationMetabRecord("sources-v1.json")
    expect_silent(metabPublicationValidateSources(sources))
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
        metabPublicationValidateSources(leaked),
        class = "multischolar_publication_schema_error"
    )

    promoted <- publicationGovernanceCopy(sources)
    promoted$decisions[[1L]]$promotion_eligible <- TRUE
    expect_error(
        metabPublicationValidateSources(promoted),
        class = "multischolar_publication_metabolomics_contract_error"
    )
})

test_that("frozen splits are complete and pilot-holdout disjoint", {
    splits <- .publicationMetabRecord("splits-v1.json")
    bundles <- .publicationMetabRecord("msdial-bundles-v1.json")
    governance <- .publicationMetabRecord("governance-successor-v1.json")
    expect_silent(metabPublicationValidateSplits(splits))
    expect_silent(metabPublicationValidateMsdialBundles(bundles))
    expect_silent(metabPublicationValidateGovernanceSuccessor(governance))
    expect_length(splits$pilot_calibration_assignments, 2L)
    expect_length(splits$assignments, 32L)
    expect_length(bundles$bundles, 4L)

    leaked <- publicationGovernanceCopy(splits)
    leaked$assignments[[1L]]$source_id <-
        leaked$pilot_calibration_assignments[[1L]]$source_id
    expect_error(
        metabPublicationValidateSplits(leaked),
        class = "multischolar_publication_metabolomics_contract_error"
    )

    aliased <- publicationGovernanceCopy(bundles)
    aliased$custom_substitution_allowed <- TRUE
    expect_error(
        metabPublicationValidateMsdialBundles(aliased),
        class = "multischolar_publication_metabolomics_contract_error"
    )

    stale <- publicationGovernanceCopy(governance)
    stale$sources$sha256 <- strrep("0", 64L)
    expect_error(
        metabPublicationValidateGovernanceSuccessor(stale),
        class = "multischolar_publication_metabolomics_contract_error"
    )
})

test_that("test-only contracts cannot masquerade as frozen workloads", {
    contract <- .publicationMetabContract("custom", "lcms_pos")
    expect_error(
        metabPublicationValidateWorkload(contract),
        class = "multischolar_publication_metabolomics_contract_error"
    )
    frozen <- publicationGovernanceCopy(contract)
    frozen$status <- "frozen"
    expect_error(
        metabPublicationValidateWorkload(
            frozen,
            validate_authorities = FALSE
        ),
        class = "multischolar_publication_metabolomics_contract_error"
    )

    widened <- publicationGovernanceCopy(contract)
    widened$truth_contract$support_widening_allowed <- TRUE
    expect_error(
        metabPublicationValidateWorkload(
            widened,
            allow_test_contract = TRUE,
            validate_authorities = FALSE
        ),
        class = "multischolar_publication_metabolomics_contract_error"
    )

    private <- publicationGovernanceCopy(contract)
    private$privacy$private_source_opened <- TRUE
    expect_error(
        metabPublicationValidateWorkload(
            private,
            allow_test_contract = TRUE,
            validate_authorities = FALSE
        ),
        class = "multischolar_publication_metabolomics_contract_error"
    )
})

test_that("negative contracts preserve rejection and invocation boundaries", {
    authority <- publicationReadJson(.publicationMetabRepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "metabolomics",
        "negative",
        "contracts-v1.json"
    ))
    expect_silent(metabPublicationValidateNegativeAuthority(authority))
    expect_length(authority$cases, 21L)
    for (case in authority$cases) {
        suffix <- if (startsWith(case$mutation, "mixed_")) "" else ".tsv"
        path <- paste0(tempfile("metabolomics-negative-"), suffix)
        prepared <- metabPublicationPrepareNegative(case, path)
        expect_true(file.exists(prepared$path) || dir.exists(prepared$path))
        evidence <- metabPublicationEvaluateNegative(case, prepared$path)
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
        .publicationMetabRepoPath(
            "tests",
            "testdata",
            "omics-performance",
            "metabolomics",
            "workloads"
        ),
        pattern = "^fixture-correctness-.*[.]json$",
        full.names = TRUE
    ))
    expect_length(paths, 8L)
    for (path in paths) {
        contract <- metabPublicationReadContract(path)
        first <- metabPublicationPrepareFixture(
            contract,
            tempfile("metabolomics-fixture-a-")
        )
        second <- metabPublicationPrepareFixture(
            contract,
            tempfile("metabolomics-fixture-b-")
        )
        expect_identical(
            first$payload$payload_set_sha256,
            second$payload$payload_set_sha256
        )
        expect_identical(first$truth$sha256, second$truth$sha256)
        expect_identical(first$receipt$sha256, second$receipt$sha256)
        result <- metabPublicationRunImported(
            contract,
            file.path(dirname(first$truth$path), "payload"),
            first$truth$path
        )
        expect_true(result$workflow_evidence$truth_valid)
        expect_true(result$workflow_evidence$differential_direction$valid)
        expect_false(result$workflow_evidence$promotion_authority)
    }
})

test_that("heavy pilot contracts are historical calibration only", {
    contracts <- lapply(names(metabPublicationCapabilities()), function(route) {
        metabPublicationBuildPilotContract(
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
        c(1L, 3L)
    )
    expect_true(all(vapply(contracts, function(contract) {
        !isTRUE(contract$claim_scope$performance_authority) &&
            !isTRUE(contract$claim_scope$promotion_authority)
    }, logical(1))))

    promoted <- publicationGovernanceCopy(contracts[[1L]])
    promoted$claim_scope$performance_authority <- TRUE
    expect_error(
        metabPublicationValidatePilotContract(promoted),
        class = "multischolar_publication_metabolomics_contract_error"
    )

    qualified <- metabPublicationPilotQualification(list(
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
        metabPublicationPilotStatus(aborted, list(qualified = FALSE)),
        "safety_aborted_no_dimension_decision"
    )

    heavy <- metabPublicationBuildHeavyContract(
        "msdial",
        "mixed",
        aggregate_feature_count = 80000L,
        sample_count = 48L,
        qualification_id = "test_qualification"
    )
    expect_silent(metabPublicationValidateWorkload(heavy))
    expect_identical(heavy$dimensions$payload_member_count, 3L)
    expect_true(heavy$claim_scope$performance_authority)

    incomplete_manifest <- list(
        schema = "multischolar.omics_publication_metabolomics_manifest"
    )
    expect_error(
        metabPublicationValidateFinalManifest(incomplete_manifest),
        class = "multischolar_publication_schema_error"
    )
})

test_that("metabolomics payload cleanup is archived and symlink safe", {
    sandbox <- withr::local_tempdir(pattern = "metabolomics-cleanup-")
    withr::local_dir(sandbox)
    root <- file.path(".omics-publication-workloads", "metabolomics")
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
        schema = "multischolar.omics_publication_metabolomics_cleanup_plan",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-066",
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
    dry_run <- metabPublicationRunCleanup(plan, execute = FALSE)
    expect_identical(dry_run$status, "dry_run")
    expect_true(dir.exists(generated))
    expect_identical(dry_run$removals[[1L]]$before$symlink_count, 1L)

    result <- metabPublicationRunCleanup(plan, execute = TRUE)
    expect_identical(result$status, "passed")
    expect_false(dir.exists(generated))
    expect_true(file.exists(file.path(protected, "receipt.json")))
    expect_true(file.exists(file.path(
        root,
        "cleanup-evidence/test/generated__attempt-1/retained/logs/run.log"
    )))
})

test_that("proteomics and lipidomics adapters reject metabolomics identity", {
    contract <- omicsWorkloadReadContract(.publicationMetabRepoPath(
        "tests",
        "testdata",
        "omics-parity",
        "metabolomics",
        "workloads",
        "lcms-pos-public-ci-v1.json"
    ))
    for (name in c("proteomics", "lipidomics")) {
        path <- .publicationMetabRepoPath(
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
