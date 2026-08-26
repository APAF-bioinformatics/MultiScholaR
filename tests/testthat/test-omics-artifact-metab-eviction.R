.metab051LoadResumeHelpers <- function() {
    path <- testthat::test_path("test-omics-artifact-metab-resume.R")
    expressions <- parse(path)
    target <- parent.frame()
    for (expression in expressions) {
        assignment <- is.call(expression) &&
            identical(expression[[1L]], as.name("<-"))
        if (isTRUE(assignment)) eval(expression, envir = target)
    }
    invisible(TRUE)
}

.metab051LoadResumeHelpers()

.metab051Resumed <- function(kind = "mixed") {
    root <- tempfile("metab051-")
    dir.create(root, recursive = TRUE)
    built <- .metab032PersistProject(root, kind)
    workflow <- .metab032FreshWorkflow(built$paths)
    result <- resumeMetabArtifactWorkflowSafely(
        workflow,
        built$paths,
        "metabolomics-study",
        log_warn = \(...) invisible(NULL)
    )
    stopifnot(isTRUE(result$resumed))
    list(root = root, built = built, workflow = workflow, resume = result)
}

test_that("metabolomics eviction inventory and readiness are exact", {
    resumed <- .metab051Resumed("mixed")
    workflow <- resumed$workflow
    inventory <- metabEvictionConsumerInventory()
    expect_setequal(unique(inventory$category), c(
        "preview", "compatibility", "qc", "normalization", "session", "da",
        "report"
    ))
    expect_true(all(inventory$verified))
    expect_true(all(nzchar(inventory$source_after_eviction)))
    contract <- metabPayloadEvictionContract()
    expect_identical(contract$source_fields, c("data_tbl", "data_cln"))
    readiness <- metabEvictionReadiness(workflow, \(...) "evict")
    expect_true(all(readiness))
    data_before <- workflow$data_tbl
    cleaned_before <- workflow$data_cln
    for (name in names(readiness)) {
        candidate <- readiness
        candidate[[name]] <- FALSE
        result <- evictMetabArtifactWorkflowPayloads(
            workflow,
            rollout_fn = \(...) "evict",
            readiness_fn = \(...) candidate
        )
        expect_false(result$evicted, info = name)
        expect_contains(result$failed_prerequisites, name)
        expect_identical(workflow$data_tbl, data_before)
        expect_identical(workflow$data_cln, cleaned_before)
    }
    expect_true(workflow$state_manager$close())
})

test_that("successful eviction releases lists and reconstructs exact assays", {
    for (kind in c("lc", "gc", "mixed")) {
        resumed <- .metab051Resumed(kind)
        workflow <- resumed$workflow
        data_before <- workflow$data_tbl
        cleaned_before <- workflow$data_cln
        source_bytes <- sum(artifactPayloadSourceFieldBytes(
            workflow,
            METAB_EVICT_FIELDS
        ))
        result <- evictMetabArtifactWorkflowPayloads(
            workflow,
            rollout_fn = \(...) "evict"
        )
        expect_true(result$evicted, info = kind)
        expect_null(workflow$data_tbl)
        expect_null(workflow$data_cln)
        expect_identical(workflow$state_manager$getCacheInfo()$entries, 0L)
        gate <- evaluateMetabEvictionStageGate(
            source_bytes,
            result,
            workflow$state_manager
        )
        expect_true(gate$passed, info = kind)
        expect_identical(
            resolveMetabWorkflowAssays(workflow, "data_tbl"),
            data_before,
            info = kind
        )
        expect_identical(
            resolveMetabWorkflowAssays(workflow, "data_cln"),
            cleaned_before,
            info = kind
        )
        expect_null(workflow$data_tbl)
        expect_null(workflow$data_cln)
        expect_true(metabWorkflowPayloadAvailable(workflow, "data_tbl"))
        expect_true(metabWorkflowPayloadAvailable(workflow, "data_cln"))
        state <- workflow$state_manager$getState()
        expect_identical(state@metabolite_data, cleaned_before)
        expect_identical(methods::validObject(state), TRUE)
        expect_true(workflow$state_manager$close())
    }
})

test_that("settled resume avoids source and S4 hydration until explicit demand", {
    for (kind in c("lc", "gc", "mixed")) {
        root <- withr::local_tempdir(pattern = paste0("metab051-settled-", kind, "-"))
        built <- .metab032PersistProject(root, kind)
        expect_true(built$design$settled_resume_snapshot$written, info = kind)
        prepared <- createMetabResumeContext(
            built$paths,
            "metabolomics-study"
        )
        bundle <- hydrateMetabResumeBundle(
            prepared$context,
            retain_source_payloads = FALSE
        )
        workflow <- .metab032FreshWorkflow(built$paths)
        applyMetabResumeBundle(workflow, bundle)
        expect_false(bundle$source_payloads_retained, info = kind)
        expect_identical(bundle$readthrough_mode, "settled", info = kind)
        expect_null(bundle$state_object, info = kind)
        expect_null(workflow$data_tbl, info = kind)
        expect_null(workflow$data_cln, info = kind)
        expect_identical(workflow$state_manager$getCacheInfo()$entries, 0L)
        expect_false(
            workflow$state_manager$getResourceInfo()$registry_connection,
            info = kind
        )
        expect_identical(
            metabWorkflowAssayNames(workflow),
            names(built$payload$assayList),
            info = kind
        )
        expect_identical(workflow$state_manager$getCacheInfo()$entries, 0L)
        result <- settleMetabArtifactWorkflowSafely(
            workflow,
            rollout_fn = \(...) "evict",
            log_warn = \(...) invisible(NULL)
        )
        expect_true(result$evicted, info = kind)
        expect_true(result$source_hydration_avoided, info = kind)
        expect_identical(result$reason, "artifact_payload_hydration_avoided")
        gate <- evaluateMetabEvictionStageGate(1, result, workflow$state_manager)
        expect_true(gate$passed, info = kind)
        expect_identical(
            resolveMetabWorkflowAssays(workflow, "data_cln"),
            built$object@metabolite_data,
            info = kind
        )
        expect_identical(workflow$state_manager$getCacheInfo()$entries, 0L)
        expect_identical(workflow$state_manager$getState(), built$object, info = kind)
        expect_identical(workflow$state_manager$getCacheInfo()$entries, 1L)
        expect_false(
            workflow$state_manager$getResourceInfo()$registry_connection,
            info = kind
        )
        expect_true(workflowStateReleaseHydrationCache(workflow$state_manager))
        expect_identical(workflow$state_manager$getCacheInfo()$entries, 0L)
        expect_true(workflow$state_manager$close())
    }
})

test_that("settled resume validates payload bytes before changing workflow state", {
    root <- withr::local_tempdir(pattern = "metab051-settled-corrupt-")
    built <- .metab032PersistProject(root, "gc")
    ref <- built$import$refs[["assay_0001"]]
    path <- file.path(root, ref$relative_path)
    connection <- file(path, open = "r+b")
    byte <- readBin(connection, what = "raw", n = 1L)
    seek(connection, where = 0L, origin = "start")
    writeBin(as.raw(bitwXor(as.integer(byte), 1L)), connection)
    close(connection)
    prepared <- createMetabResumeContext(built$paths, "metabolomics-study")
    bundle <- hydrateMetabResumeBundle(
        prepared$context,
        retain_source_payloads = FALSE
    )
    workflow <- .metab032FreshWorkflow(built$paths)
    applyMetabResumeBundle(workflow, bundle)
    expect_error(
        resolveMetabWorkflowAssays(workflow, "data_tbl"),
        class = "multischolar_artifact_hash_mismatch"
    )
    expect_true(workflow$state_manager$close())
})

test_that("clear, cache, checkpoint, and ownership failures are atomic", {
    resumed <- .metab051Resumed("mixed")
    workflow <- resumed$workflow
    data_before <- workflow$data_tbl
    cleaned_before <- workflow$data_cln
    clear_result <- evictMetabArtifactWorkflowPayloadsSafely(
        workflow,
        rollout_fn = \(...) "evict",
        clear_fn = \(owner, name, value) {
            if (name == "data_cln") stop("injected clear failure")
            setArtifactWorkflowField(owner, name, value)
        },
        log_warn = \(...) invisible(NULL)
    )
    expect_false(clear_result$evicted)
    expect_identical(workflow$data_tbl, data_before)
    expect_identical(workflow$data_cln, cleaned_before)
    cache_result <- evictMetabArtifactWorkflowPayloadsSafely(
        workflow,
        rollout_fn = \(...) "evict",
        release_cache_fn = \(...) stop("injected cache failure"),
        log_warn = \(...) invisible(NULL)
    )
    expect_false(cache_result$evicted)
    expect_identical(workflow$data_tbl, data_before)
    expect_identical(workflow$data_cln, cleaned_before)
    checkpoint <- workflow$artifact_compatibility_checkpoint
    workflow$artifact_compatibility_checkpoint$descriptor_id <- "wrong"
    blocked <- evictMetabArtifactWorkflowPayloads(
        workflow,
        rollout_fn = \(...) "evict"
    )
    expect_false(blocked$evicted)
    expect_contains(
        blocked$failed_prerequisites,
        "compatibility_checkpoint_verified"
    )
    workflow$artifact_compatibility_checkpoint <- checkpoint
    expect_true(workflow$state_manager$close())

    memory_root <- withr::local_tempdir(pattern = "metab051-memory-")
    memory_paths <- .metab032Paths(memory_root, "metab051-memory")
    memory <- .metab032Workflow(memory_paths, backend = "memory")
    expect_false(evictMetabArtifactWorkflowPayloads(
        memory,
        rollout_fn = \(...) "evict"
    )$evicted)
    expect_false(metabArtifactCoordinatorOwned(memory))
})

test_that("settlement defaults retain payloads and explicit evict releases them", {
    resumed <- .metab051Resumed("gc")
    workflow <- resumed$workflow
    default <- settleMetabArtifactWorkflowSafely(workflow)
    expect_false(default$evicted)
    expect_identical(default$reason, "rollout_below_evict")
    expect_false(is.null(workflow$data_tbl))
    result <- settleMetabArtifactWorkflowSafely(
        workflow,
        rollout_fn = \(...) "evict",
        log_warn = \(...) invisible(NULL)
    )
    expect_true(result$evicted)
    expect_null(workflow$data_tbl)
    expect_null(workflow$data_cln)
    expect_true(workflow$state_manager$close())
})

test_that("evict design persistence hands the current session to artifacts", {
    root <- withr::local_tempdir(pattern = "metab051-current-handoff-")
    built <- .metab032PersistProject(root, "mixed", rollout = "evict")
    workflow <- built$workflow
    expect_true(built$design$settlement$evicted)
    expect_null(workflow$data_tbl)
    expect_null(workflow$data_cln)
    expect_s3_class(workflow$state_manager, "ArtifactWorkflowState")
    expect_identical(workflow$state_manager$getCacheInfo()$entries, 0L)
    expect_false(workflow$state_manager$getResourceInfo()$registry_connection)
    expect_identical(workflow$state_manager$getState(), built$object)
    expect_true(workflow$state_manager$close())
})

test_that("metabolomics eviction cannot own proteomics or lipidomics contexts", {
    descriptor <- artifactMetabolomicsWorkflowDescriptor()
    expect_identical(descriptor$identity$omic_type, "metabolomics")
    expect_identical(
        metabPayloadEvictionContract()$descriptor$descriptor_id,
        "metabolomics.custom.metabolite.standard.v1"
    )
    source <- paste(readLines(
        testthat::test_path("..", "..", "R", "mod_metab_eviction_helpers.R"),
        warn = FALSE
    ), collapse = "\n")
    expect_false(grepl("proteomics|lipidomics", source))
})

test_that("auto promotion is size-aware for supported assay workflows", {
    scenarios <- list(
        metabolomics = list(
            label = "metabolomics-study",
            data_format = "custom",
            data_type = "metabolite",
            prepare = prepareMetabArtifactContext
        ),
        lipidomics = list(
            label = "lipidomics-study",
            data_format = "lipidsearch",
            data_type = "lipid",
            prepare = prepareLipidArtifactContext
        )
    )
    for (omic in names(scenarios)) {
        scenario <- scenarios[[omic]]
        make_workflow <- function(rows) {
            root <- tempfile(paste0("adaptive-", omic, "-"))
            dir.create(root)
            workflow <- new.env(parent = emptyenv())
            workflow$workflow_context <- createWorkflowContext(
                list(base_dir = root, omic_label = scenario$label),
                omic,
                scenario$label,
                storage_policy = list()
            )
            workflow$data_tbl <- list(assay = as.data.frame(matrix(
                1,
                nrow = rows,
                ncol = 8L
            )))
            workflow$data_format <- scenario$data_format
            workflow$data_type <- scenario$data_type
            workflow
        }
        small <- make_workflow(10L)
        small_result <- scenario$prepare(small)
        expect_false(small_result$enabled, info = omic)
        expect_identical(
            small$workflow_context$getStorageDecision()$effective_backend,
            "memory",
            info = omic
        )
        large <- make_workflow(300000L)
        expect_gte(
            artifactPayloadProjectedSourceBytes(large$data_tbl),
            artifactPayloadAutoPromotionGate()$minimum_projected_source_bytes
        )
        large_result <- scenario$prepare(large)
        expect_true(large_result$enabled, info = omic)
        expect_identical(large_result$decision$effective_backend, "artifact")
        expect_identical(large_result$decision$effective_rollout, "evict")
        expect_identical(
            large_result$decision$reason_code,
            "auto_promoted_new_project"
        )
    }
})
