nondiaArtifact024SkipDependencies <- function() {
    for (package in c("arrow", "DBI", "duckdb", "filelock", "ps")) {
        testthat::skip_if_not_installed(package)
    }
}

nondiaArtifact024RepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

nondiaArtifact024Scenarios <- function() {
    manifest <- jsonlite::read_json(
        nondiaArtifact024RepoPath(
            "tests", "testdata", "omics-parity", "proteomics", "manifest.json"
        ),
        simplifyVector = FALSE
    )
    manifest$fixture_scenarios
}

nondiaArtifact024Importer <- function(format) {
    switch(
        format,
        maxquant = importMaxQuantData,
        fragpipe = importFragPipeData,
        pd_tmt = importProteomeDiscovererTMTData
    )
}

nondiaArtifact024WorkflowType <- function(format) {
    if (identical(format, "pd_tmt")) "TMT" else "LFQ"
}

nondiaArtifact024CapabilityId <- function(format) {
    paste0(
        "proteomics.",
        switch(format, maxquant = "maxquant", fragpipe = "fragpipe", pd_tmt = "pd_tmt"),
        ".protein.",
        if (identical(format, "pd_tmt")) "tmt" else "lfq",
        ".v1"
    )
}

nondiaArtifact024Paths <- function(root, project_id, label = "non-dia-study") {
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

nondiaArtifact024BlankWorkflow <- function() {
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

nondiaArtifact024Build <- function(
    root,
    scenario,
    project_id = paste0("nondia-024-", scenario$scenario_id),
    label = "non-dia-study"
) {
    paths <- nondiaArtifact024Paths(root, project_id, label)
    fixture <- nondiaArtifact024RepoPath(scenario$fixture_path)
    source <- file.path(paths$source_dir, basename(fixture))
    stopifnot(file.copy(fixture, source))
    workflow <- nondiaArtifact024BlankWorkflow()
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
    imported <- suppressMessages(nondiaArtifact024Importer(scenario$format)(source))
    workflow$data_tbl <- imported$data
    workflow$data_cln <- imported$data
    workflow$data_format <- scenario$format
    workflow$data_type <- imported$data_type
    workflow$column_mapping <- imported$column_mapping
    workflow_type <- nondiaArtifact024WorkflowType(scenario$format)
    workflow$state_manager$setWorkflowType(workflow_type)
    import_result <- persistProtNonDiaImportArtifacts(
        workflow,
        imported,
        source
    )
    stopifnot(isTRUE(import_result$ok))
    runs <- unique(as.character(
        workflow$data_cln[[workflow$column_mapping$run_col]]
    ))
    groups <- ifelse(grepl("KO", runs, fixed = TRUE), "KO", "WT")
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
        "non-DIA read-through fixture",
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
        capability_id = nondiaArtifact024CapabilityId(scenario$format),
        workflow_type = workflow_type,
        scenario = scenario
    )
}

expectNondiaArtifact024StateExact <- function(expected, actual) {
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

nondiaArtifact024Normalize <- function(object) {
    invisible(utils::capture.output(
        normalized <- suppressMessages(normaliseBetweenSamples(
            object,
            normalisation_method = "cyclicloess"
        ))
    ))
    normalized
}

nondiaArtifact024RunDa <- function(object) {
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

test_that("every supported non-DIA tuple resumes exact retained payloads", {
    nondiaArtifact024SkipDependencies()
    prior_globals <- captureProtContextLegacyGlobals(c(
        "design_matrix", "contrasts", "config", "annotations", "sequences"
    ))
    withr::defer(restoreProtContextLegacyGlobals(prior_globals))
    inventory <- validateProtNonDiaReadthroughConsumerInventory(
        protNonDiaReadthroughConsumerInventory()
    )
    expect_setequal(unique(inventory$category), c(
        "preview", "normalization", "session", "da", "report",
        "compatibility"
    ))

    for (scenario in nondiaArtifact024Scenarios()) {
        root <- withr::local_tempdir(
            pattern = paste0("nondia-024-resume-", scenario$format, "-")
        )
        built <- nondiaArtifact024Build(root, scenario)
        oracle <- lapply(c(
            "data_tbl", "data_cln", "design_matrix", "contrasts_tbl",
            "config_list", "column_mapping", "uniprot_dat_cln",
            "aa_seq_tbl_final"
        ), \(name) built$workflow[[name]])
        names(oracle) <- c(
            "data_tbl", "data_cln", "design_matrix", "contrasts_tbl",
            "config_list", "column_mapping", "uniprot_dat_cln",
            "aa_seq_tbl_final"
        )
        expect_identical(unlink(built$source, force = TRUE), 0L)
        resumed <- nondiaArtifact024BlankWorkflow()
        result <- resumeProtNonDiaArtifactWorkflow(
            resumed,
            built$paths,
            "non-dia-study"
        )
        withr::defer(resumed$state_manager$close())

        expect_true(result$resumed, info = scenario$format)
        expect_identical(result$capability_id, built$capability_id)
        expect_true(result$source_payloads_retained)
        expect_identical(result$readthrough_mode, "full")
        expect_null(resumed$artifact_stage_results$eviction)
        for (name in names(oracle)) {
            expect_identical(resumed[[name]], oracle[[name]], info = paste(
                scenario$format,
                name
            ))
        }
        expectNondiaArtifact024StateExact(
            built$object,
            resumed$state_manager$getState()
        )
        expect_identical(workflowStateType(resumed$state_manager), built$workflow_type)
        expect_identical(resumed$tab_status$setup_import, "complete")
        expect_identical(resumed$tab_status$design_matrix, "complete")
        expect_identical(resumed$tab_status$quality_control, "complete")
        expect_identical(resumed$tab_status$normalization, "pending")
        expect_true(all(vapply(
            result$resource_evidence,
            \(phase) is.numeric(phase$rss_bytes) && phase$rss_bytes > 0,
            logical(1)
        )))
        expect_gt(sum(result$resource_evidence$after$object_bytes), 0)
        expect_gt(result$resource_evidence$after$state_object_bytes, 0)

        preview <- previewProtNonDiaImportArtifact(
            resumed$workflow_context,
            result$import_ref,
            limit = 3L
        )
        expect_identical(nrow(preview), 3L)
        state_input <- hydrateArtifactStageInput(
            workflow_state = resumed$state_manager
        )
        expectNondiaArtifact024StateExact(built$object, state_input$value)
        expected_norm <- nondiaArtifact024Normalize(built$object)
        actual_norm <- nondiaArtifact024Normalize(state_input$value)
        expectNondiaArtifact024StateExact(expected_norm, actual_norm)
        expect_equal(
            nondiaArtifact024RunDa(actual_norm)$logFC,
            nondiaArtifact024RunDa(expected_norm)$logFC,
            tolerance = 1e-7
        )
        contrasts <- resolveProtDaAvailableContrasts(resumed)
        expect_setequal(contrasts$values, resumed$contrasts_tbl$contrasts)
        session_data <- collectProtNormExportSessionData(
            resumed,
            list(
                correlation_filtered_obj = NULL,
                best_k = NA_integer_,
                correlation_threshold = NA_real_
            ),
            list(norm_method = "cyclicloess", ruv_mode = "skip"),
            messageFn = \(...) invisible(NULL)
        )
        expectNondiaArtifact024StateExact(
            built$object,
            session_data$current_s4_object
        )
        expect_identical(session_data$workflow_type, built$workflow_type)
        expect_identical(
            resolveProtSummaryExpectedTemplate(built$workflow_type),
            if (identical(built$workflow_type, "TMT")) {
                "TMT_report.rmd"
            } else {
                "LFQ_report.rmd"
            }
        )
        copy_inputs <- prepareProtSummaryCopyInputs(
            resumed,
            list(proteomics = built$paths),
            prepareSummaryDependenciesFn = \(...) NULL,
            catFn = \(...) invisible(NULL)
        )
        expect_identical(copy_inputs$designMatrix, resumed$design_matrix)
        expect_identical(copy_inputs$contrastsTbl, resumed$contrasts_tbl)
    }
})

test_that("moved non-DIA projects resume in fresh processes without sources", {
    nondiaArtifact024SkipDependencies()
    repository <- normalizePath(testthat::test_path("..", ".."), winslash = "/")
    for (scenario in nondiaArtifact024Scenarios()) {
        original <- withr::local_tempdir(
            pattern = paste0("nondia-024-original-", scenario$format, "-")
        )
        moved <- withr::local_tempdir(
            pattern = paste0("nondia-024-moved-", scenario$format, "-")
        )
        built <- nondiaArtifact024Build(original, scenario)
        entries <- list.files(
            original,
            all.files = TRUE,
            no.. = TRUE,
            full.names = TRUE
        )
        expect_true(all(file.copy(entries, moved, recursive = TRUE)))
        unlink(file.path(moved, "source"), recursive = TRUE, force = TRUE)
        unlink(original, recursive = TRUE, force = TRUE)
        moved_paths <- list(
            base_dir = moved,
            omic_type = "proteomics",
            omic_label = "non-dia-study",
            source_dir = file.path(moved, "source"),
            results_dir = file.path(moved, "results")
        )
        script <- tempfile("nondia-024-fresh-", fileext = ".R")
        input <- tempfile("nondia-024-fresh-", fileext = ".rds")
        output <- tempfile("nondia-024-fresh-", fileext = ".rds")
        saveRDS(moved_paths, input)
        withr::defer(unlink(c(script, input, output), force = TRUE))
        writeLines(c(
            sprintf("devtools::load_all(%s, quiet = TRUE)", deparse(repository)),
            sprintf("paths <- readRDS(%s)", deparse(input)),
            "workflow <- new.env(parent = emptyenv())",
            "workflow$state_manager <- WorkflowState$new()",
            "workflow$artifact_stage_results <- NULL",
            paste0(
                "workflow$tab_status <- list(setup_import='pending', ",
                "design_matrix='disabled', quality_control='disabled', ",
                "normalization='disabled', differential_expression='disabled', ",
                "enrichment_analysis='disabled', session_summary='disabled')"
            ),
            paste0(
                "result <- resumeProtNonDiaArtifactWorkflow(",
                "workflow, paths, 'non-dia-study')"
            ),
            paste0(
                "saveRDS(list(result=result, data=workflow$data_tbl, ",
                "state=workflow$state_manager$getState(), ",
                "identity=workflow$workflow_context$getIdentity(), ",
                "root=workflow$workflow_context$getProjectRoot()), ",
                deparse(output), ")"
            ),
            "workflow$state_manager$close()"
        ), script)
        status <- system2(
            file.path(R.home("bin"), "Rscript"),
            c("--vanilla", shQuote(script)),
            stdout = TRUE,
            stderr = TRUE
        )
        expect_identical(
            attr(status, "status"),
            NULL,
            info = paste(status, collapse = "\n")
        )
        fresh <- readRDS(output)
        expect_true(fresh$result$resumed)
        expect_identical(fresh$result$capability_id, built$capability_id)
        expect_identical(fresh$identity$project_id, built$paths$project_id)
        expect_identical(fresh$root, normalizePath(moved, winslash = "/"))
        expect_identical(fresh$data, built$workflow$data_tbl)
        expectNondiaArtifact024StateExact(built$object, fresh$state)
    }
})

test_that("non-DIA resume failures leave the previous memory workflow current", {
    nondiaArtifact024SkipDependencies()
    for (scenario in nondiaArtifact024Scenarios()) {
        root <- withr::local_tempdir(
            pattern = paste0("nondia-024-failure-", scenario$format, "-")
        )
        built <- nondiaArtifact024Build(root, scenario)
        previous <- nondiaArtifact024BlankWorkflow()
        previous$data_tbl <- data.frame(sentinel = 1)
        previous$data_cln <- data.frame(sentinel = 2)
        previous$state_manager$setWorkflowType("LFQ")
        previous_values <- list(
            data_tbl = previous$data_tbl,
            data_cln = previous$data_cln,
            state_manager = previous$state_manager,
            tab_status = previous$tab_status
        )
        error <- rlang::catch_cnd(resumeProtNonDiaArtifactWorkflow(
            previous,
            built$paths,
            "non-dia-study",
            failure_injector = \(stage, ...) {
                if (identical(stage, "after_resume_apply")) {
                    stop("injected apply failure")
                }
            }
        ))
        expect_s3_class(error, "simpleError")
        expect_identical(previous$data_tbl, previous_values$data_tbl)
        expect_identical(previous$data_cln, previous_values$data_cln)
        expect_identical(previous$state_manager, previous_values$state_manager)
        expect_identical(previous$tab_status, previous_values$tab_status)

        inventory_error <- rlang::catch_cnd(resumeProtNonDiaArtifactWorkflow(
            previous,
            built$paths,
            "non-dia-study",
            inventory_fn = \() data.frame()
        ))
        expect_s3_class(
            inventory_error,
            "multischolar_incomplete_prot_nondia_consumer_inventory"
        )
        expect_identical(previous$data_tbl, previous_values$data_tbl)
        expect_identical(previous$state_manager, previous_values$state_manager)

        ref <- artifactStoreNormalizeRef(
            built$import_result$refs$canonical_data
        )
        payload <- artifactStoreResolveFile(
            newArtifactStore(
                built$workflow$workflow_context$getPaths(),
                built$paths$project_id
            ),
            ref$relative_path,
            must_exist = TRUE
        )
        expect_true(file.rename(payload, paste0(payload, ".missing")))
        hydration_error <- rlang::catch_cnd(resumeProtNonDiaArtifactWorkflow(
            previous,
            built$paths,
            "non-dia-study"
        ))
        expect_true(inherits(hydration_error, "error"))
        expect_identical(previous$data_tbl, previous_values$data_tbl)
        expect_identical(previous$state_manager, previous_values$state_manager)
    }
})

test_that("non-DIA resume and revert preserve immutable state generations", {
    nondiaArtifact024SkipDependencies()
    for (scenario in nondiaArtifact024Scenarios()) {
        root <- withr::local_tempdir(
            pattern = paste0("nondia-024-lineage-", scenario$format, "-")
        )
        built <- nondiaArtifact024Build(root, scenario)
        resumed <- nondiaArtifact024BlankWorkflow()
        resumeProtNonDiaArtifactWorkflow(
            resumed,
            built$paths,
            "non-dia-study"
        )
        manager <- resumed$state_manager
        withr::defer(manager$close())
        parent_name <- "protein_s4_initial"
        parent_generation <- manager$getCurrentGenerationId()
        parent_refs <- manager$getStateMetadata()$artifact_refs
        child <- manager$getState()
        child@args$readthrough_child <- scenario$format
        manager$saveState("readthrough_child", child, child@args, "resume child")
        child_generation <- manager$getCurrentGenerationId()
        child_refs <- manager$getStateMetadata()$artifact_refs
        child_paths <- vapply(child_refs, `[[`, character(1), "relative_path")

        expectNondiaArtifact024StateExact(
            built$object,
            manager$revertToState(parent_name)
        )
        expect_identical(manager$getCurrentGenerationId(), parent_generation)
        expect_identical(manager$getStateMetadata()$artifact_refs, parent_refs)
        expect_identical(manager$states$readthrough_child$status, "stale")
        expectNondiaArtifact024StateExact(
            child,
            manager$revertToState("readthrough_child")
        )
        expect_identical(manager$getCurrentGenerationId(), child_generation)
        expect_identical(manager$getStateMetadata()$artifact_refs, child_refs)
        expect_true(all(file.exists(file.path(root, child_paths))))
        expect_identical(
            tail(manager$getEvents(), 1L)[[1L]]$event_type,
            "state_resumed"
        )
    }
})

test_that("tuple policy, explicit memory, and read-through disable fail closed", {
    nondiaArtifact024SkipDependencies()
    scenario <- nondiaArtifact024Scenarios()[[1L]]
    root <- withr::local_tempdir(pattern = "nondia-024-policy-")
    built <- nondiaArtifact024Build(root, scenario)

    dispatched <- nondiaArtifact024BlankWorkflow()
    dispatched_result <- resumeProtArtifactWorkflowSafely(
        dispatched,
        built$paths,
        "non-dia-study"
    )
    withr::defer(dispatched$state_manager$close())
    expect_true(dispatched_result$resumed)
    expect_identical(dispatched_result$capability_id, built$capability_id)
    expect_true(dispatched_result$source_payloads_retained)
    expect_false(is.null(dispatched$data_tbl))
    expect_false(is.null(dispatched$data_cln))
    expect_null(dispatched$artifact_stage_results$eviction)

    disabled <- nondiaArtifact024BlankWorkflow()
    result <- resumeProtArtifactWorkflowSafely(
        disabled,
        built$paths,
        "non-dia-study",
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            readthrough_enabled = FALSE
        )
    )
    expect_false(result$resumed)
    expect_identical(result$reason, "readthrough_disabled")
    expect_null(disabled$data_tbl)
    expect_true(file.exists(built$source))
    expect_true(length(list.files(file.path(root, "state"), recursive = TRUE)) > 0L)

    memory_root <- withr::local_tempdir(pattern = "nondia-024-memory-")
    memory_paths <- nondiaArtifact024Paths(memory_root, "memory-project")
    memory <- nondiaArtifact024BlankWorkflow()
    memory_result <- resumeProtArtifactWorkflowSafely(
        memory,
        memory_paths,
        "non-dia-study",
        storage_policy = list(requested_backend = "memory")
    )
    expect_false(memory_result$resumed)
    expect_false(memory_result$artifact_project)
    expect_false(dir.exists(file.path(memory_root, "state")))

    wrong <- setdiff(
        names(protNonDiaReadthroughDescriptors()),
        built$capability_id
    )[[1L]]
    cross_tuple <- createProtNonDiaResumeContext(
        built$paths,
        "non-dia-study",
        capability_id = wrong
    )
    expect_false(cross_tuple$enabled)
    expect_identical(cross_tuple$reason, "artifact_manifest_absent")
    descriptor_ids <- names(protNonDiaReadthroughDescriptors())
    expect_false(any(grepl("spectronaut", descriptor_ids, fixed = TRUE)))
    expect_true(all(vapply(
        protNonDiaReadthroughDescriptors(),
        \(descriptor) identical(descriptor$identity$data_level, "protein"),
        logical(1)
    )))
    expect_true(grepl(
        "resumeProtArtifactWorkflowSafely",
        paste(deparse(body(mod_proteomics_server)), collapse = "\n"),
        fixed = TRUE
    ))

    ambiguous_root <- withr::local_tempdir(pattern = "nondia-024-ambiguous-")
    shared_project_id <- "nondia-024-ambiguous"
    nondiaArtifact024Build(
        ambiguous_root,
        nondiaArtifact024Scenarios()[[1L]],
        project_id = shared_project_id
    )
    nondiaArtifact024Build(
        ambiguous_root,
        nondiaArtifact024Scenarios()[[2L]],
        project_id = shared_project_id
    )
    ambiguous_paths <- nondiaArtifact024Paths(
        ambiguous_root,
        shared_project_id
    )
    ambiguous <- nondiaArtifact024BlankWorkflow()
    ambiguous_result <- resumeProtArtifactWorkflowSafely(
        ambiguous,
        ambiguous_paths,
        "non-dia-study",
        log_warn = \(...) invisible(NULL)
    )
    expect_false(ambiguous_result$resumed)
    expect_identical(ambiguous_result$reason, "artifact_resume_probe_failed")
    expect_true(
        "multischolar_ambiguous_prot_nondia_artifact_project" %in%
            ambiguous_result$error_class
    )
    expect_null(ambiguous$data_tbl)
})
