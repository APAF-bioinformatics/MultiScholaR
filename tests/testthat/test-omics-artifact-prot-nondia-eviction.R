nondiaArtifact050SkipDependencies <- function() {
    for (package in c("arrow", "DBI", "duckdb", "filelock", "ps")) {
        testthat::skip_if_not_installed(package)
    }
}

nondiaArtifact050RepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

source(
    nondiaArtifact050RepoPath("tools", "profiling", "omics_workload_contract.R"),
    local = TRUE
)

nondiaArtifact050WorkloadAdapterPath <- nondiaArtifact050RepoPath(
    "tools", "profiling", "omics_workload_proteomics.R"
)

nondiaArtifact050Scenarios <- function() {
    manifest <- jsonlite::read_json(
        nondiaArtifact050RepoPath(
            "tests", "testdata", "omics-parity", "proteomics", "manifest.json"
        ),
        simplifyVector = FALSE
    )
    manifest$fixture_scenarios
}

nondiaArtifact050Importer <- function(format) {
    switch(
        format,
        maxquant = importMaxQuantData,
        fragpipe = importFragPipeData,
        pd_tmt = importProteomeDiscovererTMTData
    )
}

nondiaArtifact050WorkflowType <- function(format) {
    if (identical(format, "pd_tmt")) "TMT" else "LFQ"
}

nondiaArtifact050CapabilityId <- function(format) {
    paste0(
        "proteomics.",
        switch(format, maxquant = "maxquant", fragpipe = "fragpipe", pd_tmt = "pd_tmt"),
        ".protein.",
        if (identical(format, "pd_tmt")) "tmt" else "lfq",
        ".v1"
    )
}

nondiaArtifact050Paths <- function(root, project_id, label = "eviction-study") {
    paths <- list(
        base_dir = root,
        project_id = project_id,
        omic_type = "proteomics",
        omic_label = label,
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    paths
}

nondiaArtifact050BlankWorkflow <- function() {
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- NULL
    workflow$state_manager <- WorkflowState$new()
    workflow$data_tbl <- NULL
    workflow$data_cln <- NULL
    workflow$design_matrix <- NULL
    workflow$contrasts_tbl <- NULL
    workflow$config_list <- NULL
    workflow$column_mapping <- NULL
    workflow$uniprot_dat_cln <- NULL
    workflow$aa_seq_tbl_final <- NULL
    workflow$data_format <- NULL
    workflow$data_type <- NULL
    workflow$artifact_stage_results <- NULL
    workflow$artifact_readthrough_refs <- NULL
    workflow$artifact_readthrough_proof <- NULL
    workflow$artifact_compatibility_checkpoint <- NULL
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

nondiaArtifact050Build <- function(
    root,
    format,
    source_path,
    project_id = paste0("nondia-050-", format),
    label = "eviction-study",
    worker_owned = FALSE
) {
    paths <- nondiaArtifact050Paths(root, project_id, label)
    source <- file.path(paths$source_dir, basename(source_path))
    stopifnot(file.copy(source_path, source))
    workflow <- nondiaArtifact050BlankWorkflow()
    workflow$workflow_context <- createWorkflowContext(
        paths,
        "proteomics",
        label,
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = project_id
        )
    )
    staged <- if (isTRUE(worker_owned)) {
        suppressMessages(stageProtNonDiaImportArtifacts(
            workflow,
            source,
            format
        ))
    } else {
        NULL
    }
    imported <- if (is.list(staged) && isTRUE(staged$ok)) {
        staged$result
    } else {
        suppressMessages(nondiaArtifact050Importer(format)(source))
    }
    workflow$data_tbl <- imported$data
    workflow$data_cln <- imported$data
    workflow$data_format <- format
    workflow$data_type <- imported$data_type
    workflow$column_mapping <- imported$column_mapping
    workflow_type <- nondiaArtifact050WorkflowType(format)
    workflow$state_manager$setWorkflowType(workflow_type)
    import_result <- persistProtNonDiaImportArtifacts(
        workflow,
        imported,
        source,
        pending_stage = staged$pending_stage %||% NULL,
        worker_attempted = isTRUE(staged$attempted)
    )
    stopifnot(isTRUE(import_result$ok))
    runs <- unique(as.character(
        workflow$data_cln[[workflow$column_mapping$run_col]]
    ))
    groups <- ifelse(grepl("KO|treatment", runs), "KO", "WT")
    replicates <- ave(seq_along(runs), groups, FUN = seq_along)
    workflow$design_matrix <- data.frame(
        Run = runs,
        group = groups,
        replicates = paste0("R", replicates),
        stringsAsFactors = FALSE
    )
    workflow$config_list <- list(
        globalParameters = list(workflow_type = workflow_type)
    )
    workflow$contrasts_tbl <- data.frame(
        contrasts = "groupKO-groupWT",
        friendly_names = "KO_vs_WT",
        full_format = "KO_vs_WT=groupKO-groupWT",
        stringsAsFactors = FALSE
    )
    proteins <- unique(as.character(
        workflow$data_cln[[workflow$column_mapping$protein_col]]
    ))
    workflow$uniprot_dat_cln <- data.frame(
        Protein.Ids = proteins,
        Gene = paste0("GENE", seq_along(proteins)),
        stringsAsFactors = FALSE
    )
    workflow$aa_seq_tbl_final <- data.frame(
        accession = proteins,
        sequence = rep("PEPTIDE", length(proteins)),
        stringsAsFactors = FALSE
    )
    suppressMessages(buildProtDesignStateCheckpoint(
        workflow,
        workflow_type,
        "non-DIA eviction fixture",
        validateColumnMapping = TRUE
    ))
    object <- workflow$state_manager$getState("protein_s4_initial")
    design_result <- persistProtNonDiaDesignArtifacts(workflow)
    stopifnot(isTRUE(design_result$ok))
    list(
        paths = paths,
        source = source,
        workflow = workflow,
        imported = imported,
        object = object,
        import_result = import_result,
        design_result = design_result,
        capability_id = nondiaArtifact050CapabilityId(format),
        descriptor = protNonDiaReadthroughDescriptor(
            nondiaArtifact050CapabilityId(format)
        ),
        workflow_type = workflow_type,
        format = format
    )
}

nondiaArtifact050BuildScenario <- function(root, scenario) {
    nondiaArtifact050Build(
        root,
        scenario$format,
        nondiaArtifact050RepoPath(scenario$fixture_path)
    )
}

nondiaArtifact050Resume <- function(built) {
    workflow <- nondiaArtifact050BlankWorkflow()
    result <- resumeProtNonDiaArtifactWorkflow(
        workflow,
        built$paths,
        "eviction-study"
    )
    stopifnot(isTRUE(result$resumed))
    list(workflow = workflow, result = result)
}

expectNondiaArtifact050StateExact <- function(expected, actual) {
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

nondiaArtifact050Rss <- function() {
    unname(as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]]))
}

nondiaArtifact050ArrowBytes <- function() {
    tryCatch(
        as.numeric(arrow::default_memory_pool()$bytes_allocated),
        error = \(...) NA_real_
    )
}

nondiaArtifact050Normalize <- function(object) {
    invisible(utils::capture.output(
        normalized <- suppressMessages(normaliseBetweenSamples(
            object,
            normalisation_method = "cyclicloess"
        ))
    ))
    normalized
}

nondiaArtifact050RunDa <- function(object) {
    matrix <- as.matrix(object@protein_quant_table[
        ,
        setdiff(names(object@protein_quant_table), object@protein_id_column),
        drop = FALSE
    ])
    rownames(matrix) <- object@protein_quant_table[[object@protein_id_column]]
    design <- object@design_matrix
    rownames(design) <- design$Run
    suppressMessages(runTestsContrasts(
        data = matrix,
        contrast_strings = "KO_vs_WT=groupKO-groupWT",
        design_matrix = design,
        formula_string = "~ 0 + group",
        treat_lfc_cutoff = 0,
        eBayes_trend = TRUE,
        eBayes_robust = TRUE
    ))$results[[1L]]
}

nondiaArtifact050WorkloadPaths <- function() {
    sort(list.files(
        nondiaArtifact050RepoPath(
            "tests", "testdata", "omics-parity", "proteomics", "workloads"
        ),
        pattern = "[.]json$",
        full.names = TRUE
    ))
}

nondiaArtifact050PrepareWorkload <- function(contract) {
    adapter <- omicsWorkloadLoadAdapter(
        nondiaArtifact050WorkloadAdapterPath,
        contract
    )
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    run_dir <- tempfile("nondia-050-workload-")
    dir.create(run_dir, recursive = TRUE)
    prepared <- adapter$prepare(list(
        contract = contract,
        run_dir = run_dir,
        repetition = 1L
    ))
    list(run_dir = run_dir, prepared = prepared)
}

test_that("non-DIA eviction inventories and gates are exact-tuple bound", {
    descriptors <- protNonDiaReadthroughDescriptors()
    gates <- protNonDiaEvictionStageGates()
    expect_setequal(names(gates), names(descriptors))
    for (descriptor in descriptors) {
        inventory <- protNonDiaEvictionConsumerInventory(descriptor)
        expect_true(all(inventory$descriptor_id == descriptor$descriptor_id))
        expect_true(all(inventory$input_format == descriptor$identity$input_format))
        expect_true(all(inventory$verified))
        expect_setequal(unique(inventory$category), c(
            "preview", "normalization", "session", "da", "report",
            "compatibility"
        ))
        gate <- protNonDiaEvictionStageGate(descriptor$descriptor_id)
        expect_identical(gate$minimum_source_graph_reduction, 0.25)
        expect_identical(gate$maximum_post_eviction_cache_entries, 0L)
        expect_identical(
            artifactDescriptorMaximumRollout(descriptor),
            "dual_write"
        )
    }
})

test_that("existing-session eviction releases duplicates and preserves consumers", {
    nondiaArtifact050SkipDependencies()
    for (scenario in nondiaArtifact050Scenarios()) {
        root <- withr::local_tempdir(
            pattern = paste0("nondia-050-evict-", scenario$format, "-")
        )
        built <- nondiaArtifact050BuildScenario(root, scenario)
        resumed <- nondiaArtifact050Resume(built)
        workflow <- resumed$workflow
        withr::defer(workflow$state_manager$close())
        data_oracle <- workflow$data_tbl
        cleaned_oracle <- workflow$data_cln
        before_source_graph <- sum(vapply(
            PROT_NONDIA_EVICT_FIELDS,
            \(field) as.numeric(utils::object.size(workflow[[field]])),
            numeric(1)
        ))
        resources <- list(
            before_rss = nondiaArtifact050Rss(),
            before_arrow = nondiaArtifact050ArrowBytes()
        )
        result <- evictProtNonDiaWorkflowPayloads(
            workflow,
            built$descriptor,
            rollout_fn = \(...) "evict"
        )
        resources$after_rss <- nondiaArtifact050Rss()
        resources$after_arrow <- nondiaArtifact050ArrowBytes()

        expect_true(result$evicted, info = scenario$format)
        expect_null(workflow$data_tbl)
        expect_null(workflow$data_cln)
        expect_identical(
            workflow$state_manager$getCacheInfo()$entries,
            0L
        )
        gate <- evaluateProtNonDiaEvictionStageGate(
            built$capability_id,
            before_source_graph,
            result,
            workflow$state_manager
        )
        expect_true(gate$passed)
        expect_true(all(gate$checks))
        expect_true(all(is.finite(unlist(resources[!grepl("arrow", names(resources))]))))
        expect_true(all(vapply(resources[grepl("arrow", names(resources))], \(x) {
            is.na(x) || is.finite(x)
        }, logical(1))))
        expect_true(protNonDiaWorkflowPayloadAvailable(workflow, "data_tbl"))
        expect_true(protNonDiaWorkflowPayloadAvailable(workflow, "data_cln"))
        expect_identical(
            resolveProtNonDiaWorkflowTable(workflow, "data_tbl"),
            data_oracle
        )
        expect_identical(
            resolveProtNonDiaWorkflowTable(workflow, "data_cln"),
            cleaned_oracle
        )
        expectNondiaArtifact050StateExact(
            built$object,
            workflow$state_manager$getState()
        )
        normalized <- nondiaArtifact050Normalize(
            workflow$state_manager$getState()
        )
        expected_normalized <- nondiaArtifact050Normalize(built$object)
        expectNondiaArtifact050StateExact(expected_normalized, normalized)
        expect_equal(
            nondiaArtifact050RunDa(normalized)$logFC,
            nondiaArtifact050RunDa(expected_normalized)$logFC,
            tolerance = 1e-7
        )
        session_data <- collectProtNormExportSessionData(
            workflow,
            list(
                correlation_filtered_obj = NULL,
                best_k = NA_integer_,
                correlation_threshold = NA_real_
            ),
            list(norm_method = "cyclicloess", ruv_mode = "skip"),
            messageFn = \(...) invisible(NULL)
        )
        expectNondiaArtifact050StateExact(
            built$object,
            session_data$current_s4_object
        )
        preview <- previewProtNonDiaImportArtifact(
            workflow$workflow_context,
            resumed$result$import_ref,
            limit = 2L
        )
        expect_identical(nrow(preview), 2L)
        copy_inputs <- prepareProtSummaryCopyInputs(
            workflow,
            list(proteomics = built$paths),
            prepareSummaryDependenciesFn = \(...) NULL,
            catFn = \(...) invisible(NULL)
        )
        expect_identical(copy_inputs$designMatrix, workflow$design_matrix)
        expect_identical(copy_inputs$contrastsTbl, workflow$contrasts_tbl)
        expect_identical(
            resolveProtSummaryExpectedTemplate(built$workflow_type),
            if (identical(built$workflow_type, "TMT")) {
                "TMT_report.rmd"
            } else {
                "LFQ_report.rmd"
            }
        )
    }
})

test_that("settled resume validates sources without retaining duplicate tables", {
    nondiaArtifact050SkipDependencies()
    for (scenario in nondiaArtifact050Scenarios()) {
        original <- withr::local_tempdir(
            pattern = paste0("nondia-050-settled-original-", scenario$format, "-")
        )
        moved <- withr::local_tempdir(
            pattern = paste0("nondia-050-settled-moved-", scenario$format, "-")
        )
        built <- nondiaArtifact050BuildScenario(original, scenario)
        entries <- list.files(
            original,
            all.files = TRUE,
            no.. = TRUE,
            full.names = TRUE
        )
        expect_true(all(file.copy(entries, moved, recursive = TRUE)))
        unlink(file.path(moved, "source"), recursive = TRUE, force = TRUE)
        moved_paths <- list(
            base_dir = moved,
            omic_type = "proteomics",
            omic_label = "eviction-study",
            source_dir = file.path(moved, "source"),
            results_dir = file.path(moved, "results")
        )
        prepared <- createProtNonDiaResumeContext(
            moved_paths,
            "eviction-study",
            built$capability_id
        )
        expect_true(prepared$enabled)
        bundle <- hydrateProtNonDiaResumeBundle(
            prepared$context,
            built$descriptor,
            retain_source_payloads = FALSE
        )
        withr::defer(bundle$state_manager$close())
        expect_false(bundle$source_payloads_retained)
        expect_identical(bundle$readthrough_mode, "settled")
        expect_null(bundle$data_tbl)
        expect_null(bundle$data_cln)
        expectNondiaArtifact050StateExact(built$object, bundle$state_object)
        expect_identical(
            bundle$context$getIdentity()$project_id,
            built$paths$project_id
        )
        expect_identical(
            bundle$context$getProjectRoot(),
            normalizePath(moved, winslash = "/")
        )
    }
})

test_that("every non-DIA eviction prerequisite and mutation failure fails closed", {
    nondiaArtifact050SkipDependencies()
    for (scenario in nondiaArtifact050Scenarios()) {
        root <- withr::local_tempdir(
            pattern = paste0("nondia-050-failure-", scenario$format, "-")
        )
        built <- nondiaArtifact050BuildScenario(root, scenario)
        resumed <- nondiaArtifact050Resume(built)
        workflow <- resumed$workflow
        withr::defer(workflow$state_manager$close())
        data_before <- workflow$data_tbl
        cleaned_before <- workflow$data_cln
        readiness <- protNonDiaEvictionReadiness(
            workflow,
            built$descriptor,
            rollout_fn = \(...) "evict"
        )
        expect_true(all(readiness))
        for (name in names(readiness)) {
            candidate <- readiness
            candidate[[name]] <- FALSE
            result <- evictProtNonDiaWorkflowPayloads(
                workflow,
                built$descriptor,
                rollout_fn = \(...) "evict",
                readiness_fn = \(...) candidate
            )
            expect_false(result$evicted, info = paste(scenario$format, name))
            expect_contains(result$failed_prerequisites, name)
            expect_identical(workflow$data_tbl, data_before)
            expect_identical(workflow$data_cln, cleaned_before)
        }

        clear_result <- evictProtNonDiaWorkflowPayloadsSafely(
            workflow,
            built$descriptor,
            rollout_fn = \(...) "evict",
            clear_fn = \(owner, name, value) {
                if (identical(name, "data_cln")) stop("injected clear failure")
                setArtifactWorkflowField(owner, name, value)
            },
            log_warn = \(...) invisible(NULL)
        )
        expect_false(clear_result$evicted)
        expect_identical(workflow$data_tbl, data_before)
        expect_identical(workflow$data_cln, cleaned_before)

        cache_result <- evictProtNonDiaWorkflowPayloadsSafely(
            workflow,
            built$descriptor,
            rollout_fn = \(...) "evict",
            release_cache_fn = \(...) stop("injected cache failure"),
            log_warn = \(...) invisible(NULL)
        )
        expect_false(cache_result$evicted)
        expect_identical(workflow$data_tbl, data_before)
        expect_identical(workflow$data_cln, cleaned_before)

        checkpoint <- workflow$artifact_compatibility_checkpoint
        workflow$artifact_compatibility_checkpoint$descriptor_id <- "wrong"
        blocked <- evictProtNonDiaWorkflowPayloads(
            workflow,
            built$descriptor,
            rollout_fn = \(...) "evict"
        )
        expect_false(blocked$evicted)
        expect_contains(
            blocked$failed_prerequisites,
            "compatibility_checkpoint_verified"
        )
        workflow$artifact_compatibility_checkpoint <- checkpoint

        ref <- artifactStoreNormalizeRef(
            built$import_result$refs$canonical_data
        )
        store <- newArtifactStore(
            built$workflow$workflow_context$getPaths(),
            built$paths$project_id
        )
        payload <- artifactStoreResolveFile(
            store,
            ref$relative_path,
            must_exist = TRUE
        )
        expect_true(file.rename(payload, paste0(payload, ".missing")))
        hydration_error <- rlang::catch_cnd(hydrateProtNonDiaResumeBundle(
            workflow$workflow_context,
            built$descriptor,
            retain_source_payloads = FALSE
        ))
        expect_true(inherits(hydration_error, "error"))
        expect_identical(workflow$data_tbl, data_before)
        expect_identical(workflow$data_cln, cleaned_before)
        expectNondiaArtifact050StateExact(
            built$object,
            workflow$state_manager$getState()
        )
    }
})

test_that("current-session settlement is tuple-bound and defaults retain memory", {
    nondiaArtifact050SkipDependencies()
    for (scenario in nondiaArtifact050Scenarios()) {
        root <- withr::local_tempdir(
            pattern = paste0("nondia-050-settle-", scenario$format, "-")
        )
        built <- nondiaArtifact050BuildScenario(root, scenario)
        default <- settleProtArtifactWorkflowSafely(
            built$workflow,
            built$paths,
            "eviction-study"
        )
        expect_false(default$evicted)
        expect_identical(default$reason, "rollout_below_evict")
        expect_false(is.null(built$workflow$data_tbl))
        result <- settleProtArtifactWorkflowSafely(
            built$workflow,
            built$paths,
            "eviction-study",
            rollout_fn = \(...) "evict"
        )
        expect_true(result$evicted)
        expect_identical(result$descriptor_id, built$capability_id)
        expect_true(result$state_manager_replaced)
        expect_false(result$complete_payload_returned)
        expect_null(built$workflow$data_tbl)
        expect_null(built$workflow$data_cln)
        expect_s3_class(
            built$workflow$state_manager,
            "ArtifactWorkflowState"
        )
        expect_identical(
            built$workflow$state_manager$getCacheInfo()$entries,
            0L
        )
        expectNondiaArtifact050StateExact(
            built$object,
            built$workflow$state_manager$getState()
        )
    }
})

test_that("public generated workloads meet independent source-graph stage gates", {
    nondiaArtifact050SkipDependencies()
    for (path in nondiaArtifact050WorkloadPaths()) {
        contract <- omicsWorkloadReadContract(path)
        prepared <- nondiaArtifact050PrepareWorkload(contract)
        withr::defer(unlink(prepared$run_dir, recursive = TRUE, force = TRUE))
        format <- contract$capability$input_format
        root <- withr::local_tempdir(pattern = paste0("nondia-050-gate-", format, "-"))
        built <- nondiaArtifact050Build(
            root,
            format,
            prepared$prepared$payload_path,
            project_id = paste0("nondia-050-gate-", format),
            worker_owned = format %in% c("maxquant", "fragpipe")
        )
        if (format %in% c("maxquant", "fragpipe")) {
            expect_true(
                isTRUE(
                    built$import_result$process_evidence$distinct_workers
                ) || isTRUE(
                    built$import_result$process_evidence$sequential_bounded_inline
                )
            )
            expect_false(
                built$import_result$process_evidence$complete_payload_returned
            )
        }
        resumed <- nondiaArtifact050Resume(built)
        workflow <- resumed$workflow
        withr::defer(workflow$state_manager$close())
        before_source_graph <- sum(vapply(
            PROT_NONDIA_EVICT_FIELDS,
            \(field) as.numeric(utils::object.size(workflow[[field]])),
            numeric(1)
        ))
        before_rss <- nondiaArtifact050Rss()
        result <- evictProtNonDiaWorkflowPayloads(
            workflow,
            built$descriptor,
            rollout_fn = \(...) "evict"
        )
        after_rss <- nondiaArtifact050Rss()
        decision <- evaluateProtNonDiaEvictionStageGate(
            built$capability_id,
            before_source_graph,
            result,
            workflow$state_manager
        )
        expect_true(decision$passed, info = built$capability_id)
        expect_gte(
            decision$measured$source_graph_reduction,
            decision$gate$minimum_source_graph_reduction
        )
        expect_gt(decision$measured$released_source_bytes, 0)
        expect_true(is.finite(before_rss) && is.finite(after_rss))
        expect_identical(
            artifactDescriptorMaximumRollout(built$descriptor),
            "dual_write"
        )
    }
})
