operationalArtifactDecisionMap <- function() {
    c(
        artifactProteomicsNonDiaCloseoutDecisions(),
        artifactMetabolomicsCloseoutDecisions(),
        artifactLipidomicsCloseoutDecisions()
    )
}

test_that("all-omics scorecard is exact complete and fail closed", {
    closeout <- operationalArtifactReadCloseout()
    expect_identical(
        closeout$schema,
        "multischolar.all_omics_operational_closeout"
    )
    expect_identical(closeout$schema_version, "1.0.0")
    expect_identical(closeout$decision, "withheld")
    expect_identical(closeout$effective_default_backend, "memory")
    expect_length(closeout$promoted_capability_ids, 0L)

    capability_ids <- names(workflowCapabilityCatalogue())
    recorded_ids <- unlist(closeout$capability_ids, use.names = FALSE)
    expect_setequal(recorded_ids, capability_ids)
    expect_identical(anyDuplicated(recorded_ids), 0L)

    descriptors <- artifactWorkflowDescriptorCatalogue()$descriptors
    canaries <- closeout$descriptor_canaries
    canary_ids <- vapply(canaries, `[[`, character(1), "capability_id")
    expect_setequal(canary_ids, names(descriptors))
    expect_identical(anyDuplicated(canary_ids), 0L)
    for (canary in canaries) {
        descriptor <- descriptors[[canary$capability_id]]
        expect_identical(canary$descriptor_version, descriptor$descriptor_version)
        expect_identical(canary$descriptor_digest, descriptor$descriptor_digest)
        expect_identical(canary$promotion_status, "withheld")
        expect_identical(canary$maximum_forced_rollout, "dual_write")
        expect_false(canary$auto_eligible)
        expect_false(canary$public_private_strata_pooled)
    }
})

test_that("track decisions independently reproduce the empty promoted set", {
    closeout <- operationalArtifactReadCloseout()
    canaries <- closeout$descriptor_canaries
    canaries <- stats::setNames(
        canaries,
        vapply(canaries, `[[`, character(1), "capability_id")
    )
    decisions <- operationalArtifactDecisionMap()
    nondia_ids <- names(artifactProteomicsNonDiaWorkflowDescriptors())
    for (capability_id in nondia_ids) {
        expect_identical(
            canaries[[capability_id]]$reason_code,
            decisions[[capability_id]]$reason_code
        )
    }
    for (capability_id in c(
        "metabolomics.custom.metabolite.standard.v1",
        "lipidomics.lipidsearch.lipid.standard.v1"
    )) {
        expect_identical(
            canaries[[capability_id]]$reason_code,
            decisions[[capability_id]]$reason_code
        )
    }
    dia <- artifactDiaWorkflowDescriptor()
    expect_false(dia$certification$auto_eligible)
    expect_identical(dia$certification$status, "dual_write")
    expect_identical(
        canaries[[dia$descriptor_id]]$reason_code,
        "one_or_more_performance_gates_failed"
    )
    expect_true(all(vapply(decisions, \(decision) {
        identical(decision$promotion_status, "withheld") &&
            identical(decision$effective_default_backend, "memory")
    }, logical(1))))
})

test_that("representative generated workload contracts stay deterministic", {
    closeout <- operationalArtifactReadCloseout()
    relative_paths <- unlist(closeout$workload_contracts, use.names = FALSE)
    contracts <- lapply(
        operationalArtifactRepoPath(relative_paths),
        jsonlite::read_json,
        simplifyVector = FALSE
    )
    for (contract in contracts) {
        expect_identical(
            contract$schema,
            "multischolar.omics_workload_contract"
        )
        expect_identical(contract$schema_version, "1.0.0")
        expect_identical(contract$privacy$classification, "public_synthetic")
        expect_null(contract$privacy$scale_metadata)
        expect_identical(contract$execution$repetitions, 5L)
        expect_identical(contract$execution$threads, 1L)
        expect_true(length(contract$rng$kind) == 3L)
        expect_true(is.numeric(contract$rng$seed))
        expect_true(contract$dimensions$feature_count > 0L)
        expect_true(contract$dimensions$sample_count > 0L)
        expect_identical(contract$dimensions$assay_count, 3L)
        expect_setequal(
            names(contract$assay_mix),
            c("LCMS_Pos", "LCMS_Neg", "GCMS")
        )
        expect_match(contract$adapter$source_sha256, "^[a-f0-9]{64}$")
        expect_match(contract$expected_digests$payload_sha256, "^[a-f0-9]{64}$")
        expect_match(contract$expected_digests$truth_sha256, "^[a-f0-9]{64}$")
        expect_identical(contract$adapter_parameters$generator_version, "1.0.0")
    }
    expect_identical(
        anyDuplicated(vapply(contracts, `[[`, character(1), "workload_id")),
        0L
    )
})

test_that("portability gaps and private evidence remain fail closed", {
    closeout <- operationalArtifactReadCloseout()
    contract <- closeout$operational_contract
    expect_identical(
        contract$closed_project_portability,
        "passed_per_track"
    )
    expect_identical(
        contract$live_or_incomplete_copy_rejection,
        "absent_required_for_promotion"
    )
    expect_identical(
        contract$supported_os_promotion_matrix,
        "absent_required_for_promotion"
    )
    expect_identical(contract$schema_migration_changed, FALSE)
    boundaries <- closeout$evidence_boundaries
    expect_true(boundaries$generated_evidence_authorizes_operational_gates_only)
    expect_true(boundaries$scientific_parity_requires_hand_reviewed_fixtures)
    expect_false(boundaries$public_private_strata_pooled)
    expect_false(boundaries$private_source_path_retained)
    expect_false(boundaries$unsalted_private_fingerprint_retained)
    expect_false(boundaries$private_headers_identifiers_values_sequences_retained)
    retained_values <- unname(unlist(
        closeout,
        recursive = TRUE,
        use.names = FALSE
    ))
    encoded <- paste(as.character(retained_values), collapse = "\n")
    expect_false(grepl(
        "(/home/|/Users/|[A-Za-z]:\\\\|cotton|fasta|sequence|header)",
        encoded,
        ignore.case = TRUE
    ))
    performance_paths <- unlist(
        closeout$performance_evidence,
        use.names = FALSE
    )
    expect_true(all(file.exists(operationalArtifactRepoPath(
        performance_paths
    ))))
    expect_setequal(basename(performance_paths), c(
        "test-omics-artifact-baseline.R",
        "test-omics-artifact-prot-nondia-baseline.R",
        "test-omics-artifact-metab-baseline.R",
        "test-omics-artifact-lipid-baseline.R"
    ))
})

test_that("closed project artifact paths remain move portable", {
    operationalArtifactSkipDependencies()
    original <- withr::local_tempdir(pattern = "omics-art-048-original-")
    fixture <- operationalArtifactFixture(
        original,
        "proteomics",
        "portable-study",
        "prot_dia",
        "portable-generation"
    )
    ref <- artifactStoreWriteParquet(
        fixture$store,
        fixture$encoded,
        fixture$logical_key
    )
    moved <- withr::local_tempdir(pattern = "omics-art-048-moved-")
    entries <- list.files(original, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    expect_true(all(file.copy(entries, moved, recursive = TRUE)))
    moved_fixture <- operationalArtifactFixture(
        moved,
        "proteomics",
        "portable-study",
        "prot_dia",
        "portable-generation"
    )
    sidecar <- artifactStoreManagedPaths(
        moved_fixture$store,
        moved_fixture$logical_key,
        ref$artifact_id
    )$sidecar
    restored <- artifactStoreReadSidecar(
        moved_fixture$store,
        sidecar,
        validate_payload = TRUE
    )
    expect_identical(restored$artifact_ref, ref)
    expect_false(grepl(original, workflowSessionJson(restored), fixed = TRUE))
})
