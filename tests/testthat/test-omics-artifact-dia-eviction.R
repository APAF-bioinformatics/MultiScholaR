diaArtifact011SkipDependencies <- function() {
    for (package in c("arrow", "DBI", "duckdb", "filelock", "lobstr")) {
        testthat::skip_if_not_installed(package)
    }
}

diaArtifact011Paths <- function(root) {
    paths <- list(
        base_dir = root,
        project_id = "dia-artifact-011",
        omic_type = "proteomics",
        omic_label = "dia-eviction-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE)
    dir.create(paths$results_dir, recursive = TRUE)
    paths
}

diaArtifact011Workflow <- function(paths) {
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        "proteomics",
        "dia-eviction-study",
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = paths$project_id
        )
    )
    workflow$state_manager <- WorkflowState$new()
    workflow$artifact_stage_results <- NULL
    workflow$tab_status <- list(
        setup_import = "complete",
        design_matrix = "complete",
        quality_control = "pending",
        normalization = "disabled",
        differential_expression = "disabled",
        enrichment_analysis = "disabled",
        session_summary = "disabled"
    )
    workflow
}

diaArtifact011Rss <- function() {
    if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
    as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]])
}

diaArtifact011StateBytes <- function(workflow, hydrate = TRUE) {
    sources <- lapply(PROT_DIA_EVICT_FIELDS, \(name) workflow[[name]])
    state <- if (isTRUE(hydrate)) workflow$state_manager$getState() else NULL
    as.numeric(lobstr::obj_size(list(sources = sources, state = state)))
}

diaArtifact011Build <- function(root, copies = 1L, annotations = TRUE) {
    paths <- diaArtifact011Paths(root)
    source <- file.path(paths$source_dir, "report.tsv")
    fixture <- testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    stopifnot(file.copy(fixture, source))
    imported <- suppressMessages(importDIANNData(source))
    if (copies > 1L) {
        imported$data <- imported$data[rep(seq_len(nrow(imported$data)), copies), ]
        rownames(imported$data) <- NULL
    }
    workflow <- diaArtifact011Workflow(paths)
    workflow$data_tbl <- imported$data
    workflow$data_cln <- imported$data
    workflow$data_format <- "diann"
    workflow$data_type <- imported$data_type
    workflow$column_mapping <- imported$column_mapping
    workflow$state_manager$setWorkflowType("DIA")
    before_write <- list(
        object_bytes = diaArtifact011StateBytes(workflow, hydrate = FALSE),
        rss_bytes = diaArtifact011Rss()
    )
    import_result <- persistProtDiaImportArtifacts(workflow, imported, source)
    stopifnot(isTRUE(import_result$ok))

    runs <- unique(workflow$data_cln$Run)
    workflow$design_matrix <- data.frame(
        Run = runs,
        group = sub("_.*", "", runs),
        replicates = seq_along(runs),
        tech_rep_group = runs,
        stringsAsFactors = FALSE
    )
    workflow$config_list <- list(globalParameters = list(workflow_type = "DIA"))
    workflow$contrasts_tbl <- data.frame(
        contrasts = "groupKO-groupWT",
        friendly_names = "KO_vs_WT",
        full_format = "KO_vs_WT=groupKO-groupWT",
        stringsAsFactors = FALSE
    )
    proteins <- unique(workflow$data_cln$Protein.Group)
    if (isTRUE(annotations)) {
        workflow$uniprot_dat_cln <- data.frame(
            Protein.Group = proteins,
            Gene = paste0("GENE", seq_along(proteins)),
            stringsAsFactors = FALSE
        )
        workflow$aa_seq_tbl_final <- data.frame(
            accession = proteins,
            sequence = rep("PEPTIDE", length(proteins)),
            stringsAsFactors = FALSE
        )
    } else {
        workflow$uniprot_dat_cln <- NULL
        workflow$aa_seq_tbl_final <- NULL
    }
    object <- PeptideQuantitativeDataDiann(
        workflow$data_cln,
        workflow$design_matrix,
        technical_replicate_id = "tech_rep_group",
        args = workflow$config_list
    )
    workflow$state_manager$saveState(
        "raw_data_s4",
        object,
        workflow$config_list,
        "DIA-NN design memory checkpoint."
    )
    design_result <- persistProtDiaDesignArtifacts(workflow)
    stopifnot(isTRUE(design_result$ok))
    list(
        paths = paths,
        workflow = workflow,
        object = object,
        import_data = imported$data,
        before_write = before_write
    )
}

diaArtifact011BlankWorkflow <- function() {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new()
    workflow$artifact_stage_results <- NULL
    workflow$tab_status <- list(
        setup_import = "pending",
        design_matrix = "disabled",
        quality_control = "disabled",
        normalization = "disabled",
        differential_expression = "disabled",
        enrichment_analysis = "disabled",
        session_summary = "disabled"
    )
    workflow
}

expectDiaArtifact011StateExact <- function(expected, actual) {
    expect_identical(class(actual), class(expected))
    for (slot_name in methods::slotNames(expected)) {
        expect_identical(
            methods::slot(actual, slot_name),
            methods::slot(expected, slot_name),
            info = slot_name
        )
    }
    expect_identical(methods::validObject(actual, test = TRUE), TRUE)
}

diaArtifact011Resume <- function(built) {
    workflow <- diaArtifact011BlankWorkflow()
    result <- resumeProtDiaArtifactWorkflowSafely(
        workflow,
        built$paths,
        "dia-eviction-study"
    )
    stopifnot(isTRUE(result$resumed))
    list(workflow = workflow, result = result)
}

test_that("DIA eviction inventory covers every current reader class", {
    inventory <- protDiaEvictionConsumerInventory()
    expect_setequal(
        unique(inventory$category),
        c(
            "preview", "annotation", "design", "qc", "session", "report",
            "compatibility"
        )
    )
    expect_false(anyDuplicated(inventory$reader) > 0L)
    expect_true(all(inventory$verified))
    expect_true(all(nzchar(inventory$reader)))
    expect_true(all(nzchar(inventory$source_after_eviction)))
    expect_identical(PROT_DIA_EVICT_FIELDS, c("data_tbl", "data_cln"))
})

test_that("DIA eviction releases duplicate tables and preserves every reader", {
    diaArtifact011SkipDependencies()
    root <- withr::local_tempdir(pattern = "dia-artifact-011-evict-")
    built <- diaArtifact011Build(root, copies = 20L)
    current_settlement <- settleProtDiaArtifactWorkflowSafely(
        built$workflow,
        built$paths,
        "dia-eviction-study",
        storage_policy = built$workflow$workflow_context$getStoragePolicy(),
        rollout_fn = \(...) "evict"
    )
    expect_true(current_settlement$evicted)
    expect_null(built$workflow$data_tbl)
    expect_null(built$workflow$data_cln)
    expectDiaArtifact011StateExact(
        built$object,
        built$workflow$state_manager$getState()
    )

    resumed <- diaArtifact011Resume(built)
    workflow <- resumed$workflow
    manager <- workflow$state_manager
    withr::defer(manager$close())
    source_before <- lapply(PROT_DIA_EVICT_FIELDS, \(name) workflow[[name]])
    names(source_before) <- PROT_DIA_EVICT_FIELDS
    stages <- list(
        before_write = built$before_write,
        after_readthrough = list(
            object_bytes = diaArtifact011StateBytes(workflow),
            rss_bytes = diaArtifact011Rss()
        )
    )

    below <- evictProtDiaWorkflowPayloads(workflow)
    expect_false(below$evicted)
    expect_contains(below$failed_prerequisites, "evict_rollout")
    for (name in PROT_DIA_EVICT_FIELDS) {
        expect_identical(workflow[[name]], source_before[[name]])
    }

    evicted <- evictProtDiaWorkflowPayloads(
        workflow,
        rollout_fn = \(...) "evict"
    )
    expect_true(evicted$evicted)
    expect_true(all(evicted$released_source_field_bytes > 0))
    expect_gt(evicted$released_source_bytes_upper_bound, 0)
    expect_true(all(vapply(
        PROT_DIA_EVICT_FIELDS,
        \(name) is.null(workflow[[name]]),
        logical(1)
    )))
    expect_identical(manager$getCacheInfo()$entries, 0L)
    stages$after_eviction <- list(
        object_bytes = diaArtifact011StateBytes(workflow, hydrate = FALSE),
        rss_bytes = diaArtifact011Rss()
    )

    expect_true(protDiaWorkflowPayloadAvailable(workflow, "data_tbl"))
    expect_true(protDiaWorkflowPayloadAvailable(workflow, "data_cln"))
    expect_true(protDiaWorkflowPayloadMarker(workflow, "data_tbl"))
    rebuilt_import <- resolveProtDiaWorkflowTable(workflow, "data_tbl")
    expect_identical(rebuilt_import, built$import_data)
    expect_null(workflow$data_tbl)
    expect_null(workflow$data_cln)
    expect_identical(workflow$design_matrix, built$object@design_matrix)
    expect_false(is.null(workflow$uniprot_dat_cln))
    expect_false(is.null(workflow$aa_seq_tbl_final))

    downstream <- manager$getState()
    expectDiaArtifact011StateExact(built$object, downstream)
    expect_identical(manager$getCacheInfo()$entries, 1L)
    expect_null(workflow$data_tbl)
    expect_null(workflow$data_cln)
    stages$after_downstream_hydration <- list(
        object_bytes = diaArtifact011StateBytes(workflow),
        rss_bytes = diaArtifact011Rss()
    )
    reduction <- 1 - (
        stages$after_downstream_hydration$object_bytes /
            stages$after_readthrough$object_bytes
    )
    expect_gt(reduction, 0)
    expect_identical(
        names(stages),
        c(
            "before_write", "after_readthrough", "after_eviction",
            "after_downstream_hydration"
        )
    )
    rss <- vapply(stages, `[[`, numeric(1), "rss_bytes")
    expect_true(all(is.na(rss) | is.finite(rss)))
    expect_identical(
        evictProtDiaWorkflowPayloads(workflow, rollout_fn = \(...) "evict"),
        evicted
    )

    manager$close()
    rollback <- diaArtifact011Resume(built)
    withr::defer(rollback$workflow$state_manager$close())
    expect_identical(rollback$workflow$data_tbl, built$import_data)
    expectDiaArtifact011StateExact(
        built$object,
        rollback$workflow$state_manager$getState()
    )
    rollback$workflow$state_manager$close()

    prepared <- createProtDiaResumeContext(
        built$paths,
        "dia-eviction-study",
        built$workflow$workflow_context$getStoragePolicy()
    )
    settled_bundle <- hydrateProtDiaResumeBundle(
        prepared$context,
        retain_source_payloads = FALSE
    )
    settled_workflow <- diaArtifact011BlankWorkflow()
    applyProtDiaResumeBundle(settled_workflow, settled_bundle)
    withr::defer(settled_workflow$state_manager$close())
    expect_null(settled_workflow$data_tbl)
    expect_null(settled_workflow$data_cln)
    expect_identical(settled_bundle$readthrough_mode, "settled")
    expect_identical(
        settled_workflow$artifact_readthrough_proof$readthrough_mode,
        "settled"
    )
    expect_identical(settled_workflow$state_manager$getCacheInfo()$entries, 1L)
    expectDiaArtifact011StateExact(
        built$object,
        settled_workflow$state_manager$getState()
    )
    settled_workflow$state_manager$close()
    import_ref <- settled_workflow$artifact_readthrough_refs$import$canonical_data
    store <- newArtifactStore(
        prepared$context$getPaths(),
        prepared$context$getIdentity()$project_id
    )
    payload_path <- artifactStoreResolveFile(
        store,
        import_ref$relative_path,
        must_exist = TRUE
    )
    connection <- file(payload_path, open = "r+b")
    original_byte <- readBin(connection, what = "raw", n = 1L)
    seek(connection, where = 0L, origin = "start")
    writeBin(as.raw(bitwXor(as.integer(original_byte), 1L)), connection)
    close(connection)
    expect_error(
        hydrateProtDiaResumeBundle(
            prepared$context,
            retain_source_payloads = FALSE
        ),
        class = "multischolar_artifact_hash_mismatch"
    )
})

test_that("completed DIA annotation may be represented by an explicit absence", {
    diaArtifact011SkipDependencies()
    root <- withr::local_tempdir(pattern = "dia-artifact-011-no-annotation-")
    built <- diaArtifact011Build(root, annotations = FALSE)
    resumed <- diaArtifact011Resume(built)
    workflow <- resumed$workflow
    withr::defer(workflow$state_manager$close())

    expect_null(workflow$uniprot_dat_cln)
    expect_null(workflow$aa_seq_tbl_final)
    expect_true(workflow$artifact_readthrough_proof$annotation_completed)
    evicted <- evictProtDiaWorkflowPayloads(
        workflow,
        rollout_fn = \(...) "evict"
    )
    expect_true(evicted$evicted)
})

test_that("every failed DIA eviction prerequisite preserves source values", {
    diaArtifact011SkipDependencies()
    root <- withr::local_tempdir(pattern = "dia-artifact-011-gates-")
    built <- diaArtifact011Build(root)
    resumed <- diaArtifact011Resume(built)
    workflow <- resumed$workflow
    withr::defer(workflow$state_manager$close())
    source_before <- lapply(PROT_DIA_EVICT_FIELDS, \(name) workflow[[name]])
    names(source_before) <- PROT_DIA_EVICT_FIELDS
    baseline <- protDiaEvictionReadiness(
        workflow,
        rollout_fn = \(...) "evict"
    )
    expect_true(all(baseline))

    for (failed_name in names(baseline)) {
        candidate <- baseline
        candidate[[failed_name]] <- FALSE
        result <- evictProtDiaWorkflowPayloads(
            workflow,
            rollout_fn = \(...) "evict",
            readiness_fn = \(...) candidate,
            gc_fn = \() invisible(NULL)
        )
        expect_false(result$evicted, info = failed_name)
        expect_identical(
            result$failed_prerequisites,
            failed_name,
            info = failed_name
        )
        for (name in PROT_DIA_EVICT_FIELDS) {
            expect_identical(
                workflow[[name]],
                source_before[[name]],
                info = paste(failed_name, name)
            )
        }
    }

    calls <- 0L
    failed <- evictProtDiaWorkflowPayloadsSafely(
        workflow,
        rollout_fn = \(...) "evict",
        readiness_fn = \(...) baseline,
        clear_fn = function(workflow_data, name, value) {
            calls <<- calls + 1L
            if (calls == 2L) stop("injected clear failure")
            setProtDiaWorkflowField(workflow_data, name, value)
        },
        gc_fn = \() invisible(NULL),
        log_warn = \(...) invisible(NULL)
    )
    expect_false(failed$ok)
    expect_false(failed$evicted)
    for (name in PROT_DIA_EVICT_FIELDS) {
        expect_identical(workflow[[name]], source_before[[name]])
    }
})

test_that("memory and non-DIA workflows cannot execute DIA eviction", {
    workflow <- new.env(parent = emptyenv())
    workflow$data_tbl <- data.frame(value = 1:3)
    workflow$data_cln <- data.frame(value = 1:3)
    workflow$state_manager <- WorkflowState$new()
    workflow$state_manager$setWorkflowType("TMT")
    workflow$artifact_stage_results <- NULL
    before <- lapply(PROT_DIA_EVICT_FIELDS, \(name) workflow[[name]])
    result <- evictProtDiaWorkflowPayloads(
        workflow,
        rollout_fn = \(...) "evict"
    )
    expect_false(result$enabled)
    expect_false(result$evicted)
    expect_contains(result$failed_prerequisites, "exact_dia_canary")
    expect_identical(workflow$data_tbl, before[[1L]])
    expect_identical(workflow$data_cln, before[[2L]])
})
