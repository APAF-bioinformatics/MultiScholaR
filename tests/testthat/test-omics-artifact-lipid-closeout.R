lipid046SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

lipid046Repo <- function(...) {
    file.path(normalizePath(testthat::test_path("..", "..")), ...)
}

lipid046Context <- function(root, project_id) {
    capability <- workflowCapabilityCatalogue()[[
        "lipidomics.lipidsearch.lipid.standard.v1"
    ]]
    capability$artifact_eligible <- TRUE
    capability$maximum_artifact_rollout <- "dual_write"
    capability$explicit_maximum_artifact_rollout <- "dual_write"
    context <- createWorkflowContext(
        list(
            base_dir = root,
            project_id = project_id,
            omic_type = "lipidomics",
            omic_label = "lipidomics-study"
        ),
        "lipidomics",
        "lipidomics-study",
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = project_id
        )
    )
    bindWorkflowContextFromImport(
        context,
        workflow_type = "lipidomics_standard",
        input_format = "lipidsearch",
        data_level = "lipid",
        capabilities = list(capability)
    )
    context
}

lipid046Manager <- function(context, object) {
    manager <- ArtifactWorkflowState$new(
        workflow_context = context,
        dehydrate_fn = dehydrateLipidomicsS4Artifact,
        validate_bundle_fn = validateLipidomicsS4Bundle,
        hydrate_fn = hydrateLipidomicsS4Artifact,
        descriptor_contract = artifactStageDescriptorContract(
            artifactLipidomicsWorkflowDescriptor()
        )
    )
    manager$setWorkflowType("lipidomics_standard")
    manager$saveState(
        "lipid_norm_complete",
        object,
        object@args,
        "closeout lifecycle"
    )
    manager
}

test_that("lipidomics closeout withholds every tuple fail closed", {
    decisions <- artifactLipidomicsCloseoutDecisions()
    expect_true(all(vapply(
        decisions,
        function(decision) identical(decision$promotion_status, "withheld"),
        logical(1)
    )))
    expect_true(all(vapply(
        decisions,
        function(decision) identical(decision$effective_default_backend, "memory"),
        logical(1)
    )))
    custom <- decisions[["lipidomics.lipidsearch.lipid.standard.v1"]]
    expect_identical(custom$support_status, "scientifically_supported")
    expect_identical(
        custom$reason_code,
        "eviction_and_performance_gates_failed"
    )
    expect_false(custom$certification$auto_eligible)
    expect_identical(custom$maximum_forced_rollout, "dual_write")
    expect_true(custom$performance_evidence$checks[["enough_pairs"]])
    expect_false(custom$performance_evidence$checks[["enough_runtime_pairs"]])
    expect_false(custom$performance_evidence$checks[["eviction_passed"]])
    expect_false(custom$performance_evidence$checks[["peak_rss"]])
    expect_false(custom$performance_evidence$checks[["runtime"]])
    expect_lt(
        custom$performance_evidence$measured$peak_rss_reduction_fraction,
        0
    )
    expect_gt(custom$performance_evidence$measured$runtime_ratio, 1)
    msdial <- decisions[["lipidomics.msdial.lipid.standard.v1"]]
    expect_identical(msdial$support_status, "reader_characterized")
    expect_identical(msdial$maximum_forced_rollout, "none")
})

test_that("closeout validation rejects forged promotion and complete checks", {
    decision <- artifactLipidomicsCloseoutDecisions()[[1L]]
    promoted <- decision
    promoted$promotion_status <- "promoted"
    expect_error(
        validateArtifactLipidomicsCloseout(promoted),
        class = "multischolar_invalid_lipid_closeout"
    )
    forged <- decision
    forged$performance_evidence$checks[] <- TRUE
    expect_error(
        validateArtifactLipidomicsCloseout(forged),
        class = "multischolar_incomplete_lipid_closeout_evidence"
    )
    non_finite <- decision
    non_finite$performance_evidence$artifact$settled_rss_ratio <- Inf
    non_finite$performance_evidence$checks[["finite"]] <- FALSE
    expect_identical(
        validateArtifactLipidomicsCloseout(non_finite)$promotion_status,
        "withheld"
    )
})

test_that("E2E manifest freezes exact no-promotion lipidomics scope", {
    manifest <- jsonlite::read_json(
        lipid046Repo("tests", "testdata", "e2e", "manifest.json"),
        simplifyVector = FALSE
    )
    lanes <- manifest$lanes
    names(lanes) <- vapply(lanes, `[[`, character(1), "lane_id")
    closeout <- lanes[["lipid_canonical"]]$artifact_closeout
    expect_identical(
        closeout$capability_id,
        "lipidomics.lipidsearch.lipid.standard.v1"
    )
    expect_false(closeout$promotion_candidate)
    expect_identical(
        closeout$reason_code,
        "eviction_and_performance_gates_failed"
    )
    expect_false(closeout$generated_evidence_authorizes_science)
})

test_that("frozen representative contract reproduces exact input digests", {
    repo <- lipid046Repo()
    adapter_env <- new.env(parent = globalenv())
    sys.source(
        file.path(repo, "tools", "profiling", "omics_workload_lipidomics.R"),
        envir = adapter_env
    )
    contract <- jsonlite::read_json(
        file.path(
            repo, "tests", "testdata", "omics-parity", "lipidomics",
            "workloads", "mixed-public-representative-v1.json"
        ),
        simplifyVector = FALSE
    )
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    prepared <- adapter_env$lipidWorkloadPrepare(list(
        contract = contract,
        run_dir = withr::local_tempdir()
    ))
    expect_identical(
        digest::digest(
            prepared$payload_path,
            algo = "sha256",
            serialize = FALSE,
            file = TRUE
        ),
        contract$expected_digests$payload_sha256
    )
    expect_identical(
        digest::digest(
            prepared$truth_path,
            algo = "sha256",
            serialize = FALSE,
            file = TRUE
        ),
        contract$expected_digests$truth_sha256
    )
    evidence <- artifactLipidomicsCloseoutEvidence()
    expect_identical(evidence$workload_id, contract$workload_id)
    expect_identical(evidence$generator_version, contract$adapter_parameters$generator_version)
    expect_identical(evidence$assay_order, names(contract$assay_mix))
    expect_identical(evidence$dimensions$features, contract$dimensions$feature_count)
    expect_true(evidence$comparable)
})

test_that("closeout keeps private calibration isolated and data-free", {
    evidence <- artifactLipidomicsCloseoutEvidence()
    private <- evidence$private_evidence
    expect_false(private$pooled_with_public)
    expect_false(private$source_path_retained)
    expect_false(private$unsalted_fingerprint_retained)
    expect_false(private$headers_or_values_retained)
    encoded <- workflowSessionJson(private)
    expect_false(grepl("/home|cotton|fasta|sequence", encoded, ignore.case = TRUE))
})

test_that("repeated artifact lifecycle closes writers locks and temporary paths", {
    lipid046SkipDependencies()
    global_names <- c(
        "project_dirs", "filtering_progress_lipidomics", "config_list"
    )
    before_globals <- vapply(
        global_names,
        exists,
        logical(1),
        envir = .GlobalEnv,
        inherits = FALSE
    )
    for (cycle in seq_len(5L)) {
        root <- withr::local_tempdir(pattern = paste0("lipid046-", cycle, "-"))
        context <- lipid046Context(root, paste0("omics-art-046-cycle-", cycle))
        object <- module_ci_lipid_norm_object(
            layout = if (cycle %% 2L) "combined" else "gc",
            positive_only = TRUE
        )
        manager <- lipid046Manager(context, object)
        expect_identical(methods::validObject(manager$getState(), test = TRUE), TRUE)
        expect_true(manager$close())
        lock_dir <- artifactPath(context$getPaths(), "locks")
        lock_path <- file.path(lock_dir, "project-registry.lock")
        lock <- filelock::lock(lock_path, timeout = 0)
        expect_false(is.null(lock), info = paste("cycle", cycle))
        filelock::unlock(lock)
    }
    after_globals <- vapply(
        global_names,
        exists,
        logical(1),
        envir = .GlobalEnv,
        inherits = FALSE
    )
    expect_identical(after_globals, before_globals)
})

test_that("closeout rollback remains memory default and explicit migration", {
    decisions <- artifactLipidomicsCloseoutDecisions()
    for (decision in decisions) {
        expect_identical(decision$rollback$target_backend, "memory")
        expect_true(decision$existing_projects_require_explicit_migration)
        expect_identical(decision$effective_default_backend, "memory")
    }
    capabilities <- workflowCapabilityCatalogue()
    custom <- capabilities[["lipidomics.lipidsearch.lipid.standard.v1"]]
    expect_false(custom$auto_eligible)
    expect_false(custom$artifact_eligible)
})
