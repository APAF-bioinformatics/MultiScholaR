.lipid052LoadResumeHelpers <- function() {
    path <- testthat::test_path("test-omics-artifact-lipid-resume.R")
    expressions <- parse(path)
    target <- parent.frame()
    for (expression in expressions) {
        assignment <- is.call(expression) &&
            identical(expression[[1L]], as.name("<-"))
        if (isTRUE(assignment)) eval(expression, envir = target)
    }
    invisible(TRUE)
}

.lipid052LoadResumeHelpers()

.lipid052Resumed <- function(kind = "mixed") {
    root <- tempfile("lipid052-")
    dir.create(root, recursive = TRUE)
    built <- .lipid041PersistProject(root, kind)
    workflow <- .lipid041FreshWorkflow(built$paths)
    result <- resumeLipidArtifactWorkflowSafely(
        workflow,
        built$paths,
        "lipidomics-study",
        log_warn = \(...) invisible(NULL)
    )
    stopifnot(isTRUE(result$resumed))
    list(root = root, built = built, workflow = workflow, resume = result)
}

test_that("lipidomics eviction inventory and readiness are exact", {
    resumed <- .lipid052Resumed("mixed")
    workflow <- resumed$workflow
    inventory <- lipidEvictionConsumerInventory()
    expect_setequal(unique(inventory$category), c(
        "preview", "compatibility", "qc", "normalization", "session", "da",
        "report"
    ))
    expect_true(all(inventory$verified))
    expect_true(all(nzchar(inventory$source_after_eviction)))
    contract <- lipidPayloadEvictionContract()
    expect_identical(contract$source_fields, c("data_tbl", "data_cln"))
    readiness <- lipidEvictionReadiness(workflow, \(...) "evict")
    expect_true(all(readiness))
    data_before <- workflow$data_tbl
    cleaned_before <- workflow$data_cln
    for (name in names(readiness)) {
        candidate <- readiness
        candidate[[name]] <- FALSE
        result <- evictLipidArtifactWorkflowPayloads(
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
        resumed <- .lipid052Resumed(kind)
        workflow <- resumed$workflow
        data_before <- workflow$data_tbl
        cleaned_before <- workflow$data_cln
        source_bytes <- sum(artifactPayloadSourceFieldBytes(
            workflow,
            LIPID_EVICT_FIELDS
        ))
        result <- evictLipidArtifactWorkflowPayloads(
            workflow,
            rollout_fn = \(...) "evict"
        )
        expect_true(result$evicted, info = kind)
        expect_null(workflow$data_tbl)
        expect_null(workflow$data_cln)
        expect_identical(workflow$state_manager$getCacheInfo()$entries, 0L)
        gate <- evaluateLipidEvictionStageGate(
            source_bytes,
            result,
            workflow$state_manager
        )
        expect_true(gate$passed, info = kind)
        expect_identical(
            resolveLipidWorkflowAssays(workflow, "data_tbl"),
            data_before,
            info = kind
        )
        expect_identical(
            resolveLipidWorkflowAssays(workflow, "data_cln"),
            cleaned_before,
            info = kind
        )
        expect_null(workflow$data_tbl)
        expect_null(workflow$data_cln)
        expect_true(lipidWorkflowPayloadAvailable(workflow, "data_tbl"))
        expect_true(lipidWorkflowPayloadAvailable(workflow, "data_cln"))
        state <- workflow$state_manager$getState()
        expect_identical(state@lipid_data, cleaned_before)
        expect_identical(methods::validObject(state), TRUE)
        expect_true(workflow$state_manager$close())
    }
})

test_that("clear, cache, checkpoint, and ownership failures are atomic", {
    resumed <- .lipid052Resumed("mixed")
    workflow <- resumed$workflow
    data_before <- workflow$data_tbl
    cleaned_before <- workflow$data_cln
    clear_result <- evictLipidArtifactWorkflowPayloadsSafely(
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
    cache_result <- evictLipidArtifactWorkflowPayloadsSafely(
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
    blocked <- evictLipidArtifactWorkflowPayloads(
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

    memory_root <- withr::local_tempdir(pattern = "lipid052-memory-")
    memory_paths <- .lipid041Paths(memory_root, "lipid052-memory")
    memory <- .lipid041Workflow(memory_paths, backend = "memory")
    expect_false(evictLipidArtifactWorkflowPayloads(
        memory,
        rollout_fn = \(...) "evict"
    )$evicted)
    expect_false(lipidArtifactCoordinatorOwned(memory))
})

test_that("settlement defaults retain payloads and explicit evict releases them", {
    resumed <- .lipid052Resumed("gc")
    workflow <- resumed$workflow
    default <- settleLipidArtifactWorkflowSafely(workflow)
    expect_false(default$evicted)
    expect_identical(default$reason, "rollout_below_evict")
    expect_false(is.null(workflow$data_tbl))
    result <- settleLipidArtifactWorkflowSafely(
        workflow,
        rollout_fn = \(...) "evict",
        log_warn = \(...) invisible(NULL)
    )
    expect_true(result$evicted)
    expect_null(workflow$data_tbl)
    expect_null(workflow$data_cln)
    expect_true(workflow$state_manager$close())
})

test_that("lipidomics eviction cannot own proteomics or metabolomics contexts", {
    descriptor <- artifactLipidomicsWorkflowDescriptor()
    expect_identical(descriptor$identity$omic_type, "lipidomics")
    expect_identical(
        lipidPayloadEvictionContract()$descriptor$descriptor_id,
        "lipidomics.lipidsearch.lipid.standard.v1"
    )
    source <- paste(readLines(
        testthat::test_path("..", "..", "R", "mod_lipid_eviction_helpers.R"),
        warn = FALSE
    ), collapse = "\n")
    expect_false(grepl("proteomics|metabolomics", source))
})
