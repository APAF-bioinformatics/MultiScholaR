operationalArtifactLoadAssignments(
    "test-omics-artifact-dia-normalization.R",
    c(
        "omics014Capability", "omics014Context", "omics014Manager",
        "omics014Protein"
    )
)

test_that("all-omics store reconciliation is idempotent at every write boundary", {
    operationalArtifactSkipDependencies()
    specifications <- list(
        dia = c("proteomics", "proteomics-study", "prot_dia"),
        nondia = c("proteomics", "nondia-study", "prot_lfq"),
        metabolomics = c("metabolomics", "metabolomics-study", "metab_standard"),
        lipidomics = c("lipidomics", "lipidomics-study", "lipid_standard")
    )
    stages <- c(
        "before_write", "after_temp_write", "after_validation",
        "after_payload_rename", "after_sidecar_rename"
    )
    expected <- c(
        before_write = "intent_only",
        after_temp_write = "unvalidated_temporary",
        after_validation = "validated_temporary",
        after_payload_rename = "payload_published_sidecar_pending",
        after_sidecar_rename = "committed"
    )

    for (index in seq_along(stages)) {
        stage <- stages[[index]]
        specification <- specifications[[((index - 1L) %% 4L) + 1L]]
        root <- withr::local_tempdir(
            pattern = paste0("omics-art-048-recovery-", index, "-")
        )
        fixture <- operationalArtifactFixture(
            root,
            specification[[1L]],
            specification[[2L]],
            specification[[3L]],
            paste0("generation-", index)
        )
        expect_error(
            artifactStoreWriteParquet(
                fixture$store,
                fixture$encoded,
                fixture$logical_key,
                failure_injector = operationalArtifactFailAt(stage)
            ),
            class = "multischolar_test_operational_failure",
            info = stage
        )
        inventory <- reconcileArtifactStore(fixture$store)
        expect_identical(inventory$state, unname(expected[[stage]]), info = stage)
        first <- reconcileArtifactStore(fixture$store, repair = TRUE)
        second <- reconcileArtifactStore(fixture$store, repair = TRUE)
        expect_identical(first, second, info = stage)
        if (!identical(stage, "before_write")) {
            expect_identical(first$state, "committed", info = stage)
        }
    }
})

test_that("registry failure preserves the exact current generation", {
    operationalArtifactSkipDependencies()
    root <- withr::local_tempdir(pattern = "omics-art-048-current-")
    context <- omics014Context(root, "omics-art-048-current")
    manager <- omics014Manager(context)
    withr::defer(manager$close())
    object <- omics014Protein()
    manager$saveState(
        "stable_state",
        object,
        object@args,
        "stable operational state"
    )
    generation_id <- manager$getCurrentGenerationId()
    snapshot <- manager$exportState()

    expect_error(
        manager$saveState(
            "must_not_commit",
            object,
            object@args,
            "injected registry failure",
            failure_injector = operationalArtifactFailAt(
                "before_state_registry_commit"
            )
        ),
        class = "multischolar_test_operational_failure"
    )
    expect_identical(manager$getCurrentGenerationId(), generation_id)
    expect_identical(manager$exportState(), snapshot)
    expect_identical(manager$getState(), object)
    expect_false(manager$hasState("must_not_commit"))
})

test_that("recovery evidence covers every shared destructive boundary owner", {
    closeout <- operationalArtifactReadCloseout()
    paths <- unlist(closeout$boundary_evidence, use.names = FALSE)
    expect_true(all(file.exists(operationalArtifactRepoPath(paths))))
    expect_true(any(grepl("artifact-store", paths, fixed = TRUE)))
    expect_true(any(grepl("artifact-duckdb", paths, fixed = TRUE)))
    expect_true(any(grepl("workflow-state", paths, fixed = TRUE)))
    expect_true(any(grepl("session-lifecycle", paths, fixed = TRUE)))
    expect_true(any(grepl("artifact-hybrid", paths, fixed = TRUE)))
    expect_true(all(c(
        "passed_shared_kernel",
        "passed_per_track"
    ) %in% unlist(closeout$operational_contract, use.names = FALSE)))
})
