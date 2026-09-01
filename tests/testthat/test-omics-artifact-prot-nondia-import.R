nondiaArtifact023SkipDependencies <- function() {
    for (package in c("arrow", "DBI", "duckdb", "filelock")) {
        testthat::skip_if_not_installed(package)
    }
}

nondiaArtifact023RepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

nondiaArtifact023Scenarios <- function() {
    manifest <- jsonlite::read_json(
        nondiaArtifact023RepoPath(
            "tests", "testdata", "omics-parity", "proteomics", "manifest.json"
        ),
        simplifyVector = FALSE
    )
    manifest$fixture_scenarios
}

nondiaArtifact023Importer <- function(format) {
    switch(
        format,
        maxquant = importMaxQuantData,
        fragpipe = importFragPipeData,
        pd_tmt = importProteomeDiscovererTMTData
    )
}

nondiaArtifact023Paths <- function(root, project_id, omic_type = "proteomics") {
    paths <- list(
        base_dir = root,
        project_id = project_id,
        omic_type = omic_type,
        omic_label = "non-dia-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE)
    dir.create(paths$results_dir, recursive = TRUE)
    paths
}

nondiaArtifact023Workflow <- function(
    root,
    scenario,
    backend = "artifact",
    omic_type = "proteomics"
) {
    project_id <- paste0("nondia-023-", scenario$scenario_id)
    paths <- nondiaArtifact023Paths(root, project_id, omic_type)
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        omic_type,
        "non-dia-study",
        storage_policy = list(
            requested_backend = backend,
            requested_rollout = "dual_write",
            migration_requested = identical(backend, "artifact"),
            project_id = project_id
        )
    )
    workflow$state_manager <- WorkflowState$new()
    workflow$artifact_stage_results <- NULL
    workflow$tab_status <- list(
        setup_import = "complete",
        design_matrix = "pending",
        protein_qc = "disabled"
    )
    workflow
}

nondiaArtifact023Import <- function(workflow, scenario, source_path = NULL) {
    if (is.null(source_path)) {
        source_path <- nondiaArtifact023RepoPath(scenario$fixture_path)
    }
    imported <- suppressMessages(
        nondiaArtifact023Importer(scenario$format)(source_path)
    )
    workflow$data_tbl <- imported$data
    workflow$data_cln <- imported$data
    workflow$data_format <- scenario$format
    workflow$data_type <- imported$data_type
    workflow$column_mapping <- imported$column_mapping
    workflow_type <- if (identical(scenario$format, "pd_tmt")) "TMT" else "LFQ"
    workflow$state_manager$setWorkflowType(workflow_type)
    list(result = imported, source_path = source_path, workflow_type = workflow_type)
}

nondiaArtifact023Design <- function(workflow, workflow_type) {
    runs <- unique(as.character(workflow$data_cln[[workflow$column_mapping$run_col]]))
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
    proteins <- unique(
        as.character(workflow$data_cln[[workflow$column_mapping$protein_col]])
    )
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
        "non-DIA artifact fixture",
        validateColumnMapping = TRUE
    ))
    workflow$state_manager$getState("protein_s4_initial")
}

nondiaArtifact023ReadRef <- function(context, ref) {
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    ref <- artifactStoreNormalizeRef(ref)
    managed <- artifactStoreManagedPaths(store, ref$logical_key, ref$artifact_id)
    sidecar <- artifactStoreReadSidecar(store, managed$sidecar, TRUE)
    payload <- arrow::read_parquet(
        artifactStoreResolveFile(store, ref$relative_path, TRUE),
        as_data_frame = FALSE
    )
    decodeArtifactRectangular(payload, sidecar$codec_metadata)
}

nondiaArtifact023Registry <- function(context) {
    registry <- projectRegistryForContext(context)
    openProjectRegistryReadOnly(registry)
}

nondiaArtifact023CorruptingManager <- function(
    workflow_context,
    workflow_descriptor,
    descriptor_catalogue,
    codec_catalogue
) {
    descriptor <- findArtifactWorkflowDescriptor(
        workflow_context$getIdentity(),
        descriptor_catalogue
    )
    expect_identical(descriptor, workflow_descriptor)
    adapter <- artifactCodecAdapter(descriptor, codec_catalogue)
    ArtifactWorkflowState$new(
        workflow_context,
        dehydrate_fn = adapter$dehydrate,
        validate_bundle_fn = adapter$validate,
        hydrate_fn = function(bundle) {
            value <- adapter$hydrate(bundle)
            value@args$changed_after_hydration <- TRUE
            value
        },
        descriptor_contract = artifactStageDescriptorContract(descriptor)
    )
}

test_that("only supported non-DIA protein tuples receive dual-write descriptors", {
    descriptors <- artifactProteomicsNonDiaWorkflowDescriptors()
    expect_setequal(names(descriptors), c(
        "proteomics.maxquant.protein.lfq.v1",
        "proteomics.fragpipe.protein.lfq.v1",
        "proteomics.pd_tmt.protein.tmt.v1"
    ))
    expect_true(all(vapply(descriptors, function(descriptor) {
        identical(descriptor$identity$data_level, "protein") &&
            identical(descriptor$certification$status, "dual_write") &&
            identical(descriptor$certification$auto_eligible, FALSE) &&
            identical(
                artifactDescriptorMaximumRollout(descriptor),
                "dual_write"
            )
    }, logical(1))))
    catalogue <- artifactWorkflowDescriptorCatalogue()
    expect_setequal(
        names(catalogue$descriptors),
        c(
            "proteomics.diann.peptide.dia.v1",
            names(descriptors),
            "metabolomics.custom.metabolite.standard.v1",
            "lipidomics.lipidsearch.lipid.standard.v1"
        )
    )
    expect_false(any(grepl("spectronaut", names(catalogue$descriptors))))
    expect_false(any(vapply(descriptors, function(descriptor) {
        identical(descriptor$identity$data_level, "peptide")
    }, logical(1))))
    specs <- artifactProteomicsNonDiaDescriptorSpecs()
    for (capability_id in names(specs)) {
        resource <- jsonlite::read_json(
            nondiaArtifact023RepoPath(specs[[capability_id]]$resource_id),
            simplifyVector = TRUE
        )
        thresholds <- descriptors[[capability_id]]$evidence$performance_thresholds
        observed <- unlist(
            resource$release_gates[names(thresholds)],
            use.names = TRUE
        )
        expect_equal(observed, thresholds, info = capability_id)
        expect_true(all(c(
            specs[[capability_id]]$resource_id,
            specs[[capability_id]]$workload_id
        ) %in% descriptors[[capability_id]]$evidence$inventory_ids))
    }
})

test_that("supported non-DIA imports dual-write exact portable provenance", {
    nondiaArtifact023SkipDependencies()
    for (scenario in nondiaArtifact023Scenarios()) {
        root <- withr::local_tempdir()
        workflow <- nondiaArtifact023Workflow(root, scenario)
        source <- file.path(root, basename(scenario$fixture_path))
        expect_true(file.copy(
            nondiaArtifact023RepoPath(scenario$fixture_path),
            source
        ))
        imported <- nondiaArtifact023Import(workflow, scenario, source)
        memory_before <- imported$result$data
        output <- persistProtImportArtifacts(
            workflow,
            imported$result,
            source,
            log_warn = function(...) NULL
        )
        expect_true(output$enabled, info = scenario$scenario_id)
        expect_true(output$ok, info = scenario$scenario_id)
        expect_true(output$committed, info = scenario$scenario_id)
        expect_identical(workflow$data_tbl, memory_before)
        expect_identical(
            nondiaArtifact023ReadRef(
                workflow$workflow_context,
                output$refs$canonical_data
            ),
            memory_before
        )

        session <- nondiaArtifact023Registry(workflow$workflow_context)
        identity <- workflow$workflow_context$getIdentity()
        sources <- projectRegistryQuery(
            session,
            "sources",
            filters = list(
                workflow_id = identity$workflow_id,
                run_id = output$run_id
            )
        )
        parameters <- projectRegistryQuery(
            session,
            "parameters",
            filters = list(
                workflow_id = identity$workflow_id,
                run_id = output$run_id
            )
        )
        spec <- protNonDiaArtifactImportSpec(scenario$format)
        expect_identical(sources$source_digest, artifactByteDigest(source))
        expect_true(is.na(sources$source_uri))
        expect_identical(sources$parser_id, spec$parser_id)
        expect_identical(sources$format_id, spec$format_id)
        expect_identical(sources$data_level, "protein")
        expect_true(all(c(
            "capability_id", "column_mapping", "parser_parameters",
            "workflow_type", "input_format", "data_level"
        ) %in% parameters$parameter_name))
        closeProjectRegistry(session)
        expect_identical(unlink(source, force = TRUE), 0L)
        expect_identical(
            nondiaArtifact023ReadRef(
                workflow$workflow_context,
                output$refs$canonical_data
            ),
            memory_before
        )
    }
})

test_that("MaxQuant and FragPipe imports stage in payload-free workers", {
    nondiaArtifact023SkipDependencies()
    testthat::skip_if_not_installed("processx")
    scenarios <- Filter(
        \(scenario) scenario$format %in% c("maxquant", "fragpipe"),
        nondiaArtifact023Scenarios()
    )

    for (scenario in scenarios) {
        root <- tempfile(paste0("nondia-worker-", scenario$format, "-"))
        dir.create(root)
        withr::defer(unlink(root, recursive = TRUE, force = TRUE))
        workflow <- nondiaArtifact023Workflow(root, scenario)
        source <- nondiaArtifact023RepoPath(scenario$fixture_path)
        expected <- suppressMessages(
            nondiaArtifact023Importer(scenario$format)(source)
        )

        staged <- suppressMessages(stageProtNonDiaImportArtifacts(
            workflow,
            source,
            scenario$format
        ))
        expect_true(staged$ok, info = scenario$format)
        expect_true(staged$process_evidence$distinct_workers)
        expect_false(staged$process_evidence$complete_payload_returned)
        expect_identical(staged$result, expected)
        expect_null(workflow$artifact_stage_results)

        workflow$data_tbl <- staged$result$data
        workflow$data_cln <- staged$result$data
        workflow$data_format <- scenario$format
        workflow$data_type <- staged$result$data_type
        workflow$column_mapping <- staged$result$column_mapping
        workflow$state_manager$setWorkflowType(
            if (identical(scenario$format, "pd_tmt")) "TMT" else "LFQ"
        )
        committed <- persistProtNonDiaImportArtifacts(
            workflow,
            staged$result,
            source,
            pending_stage = staged$pending_stage,
            worker_attempted = TRUE
        )
        expect_true(committed$committed)
        expect_identical(
            workflow$artifact_stage_results$import$process_evidence,
            staged$process_evidence
        )
        expect_identical(
            nondiaArtifact023ReadRef(
                workflow$workflow_context,
                committed$refs$canonical_data
            ),
            expected$data
        )
    }
})

test_that("non-DIA import worker failures discard unpublished generations", {
    nondiaArtifact023SkipDependencies()
    testthat::skip_if_not_installed("processx")
    scenario <- Filter(
        \(candidate) identical(candidate$format, "maxquant"),
        nondiaArtifact023Scenarios()
    )[[1L]]

    for (failure in list(
        list(writer_failure_stage = "after_write"),
        list(verifier_failure_stage = "after_verify")
    )) {
        root <- tempfile("nondia-worker-failure-")
        dir.create(root)
        withr::defer(unlink(root, recursive = TRUE, force = TRUE))
        workflow <- nondiaArtifact023Workflow(root, scenario)
        source <- nondiaArtifact023RepoPath(scenario$fixture_path)
        result <- suppressMessages(do.call(
            stageProtNonDiaImportArtifactsSafely,
            c(list(
                workflow_data = workflow,
                source_path = source,
                format = scenario$format,
                log_warn = \(...) invisible(NULL)
            ), failure)
        ))

        expect_false(result$ok)
        expect_true(result$attempted)
        expect_null(result$result)
        audit <- workflow$artifact_stage_results$import_worker
        expect_null(audit$result)
        expect_false(isTRUE(audit$ok))
        expect_length(list.files(
            file.path(root, "artifacts", "intents"),
            all.files = TRUE,
            no.. = TRUE
        ), 0L)
    }
})

test_that("supported non-DIA design checkpoints dual-write and hydrate exactly", {
    nondiaArtifact023SkipDependencies()
    for (scenario in nondiaArtifact023Scenarios()) {
        root <- withr::local_tempdir()
        workflow <- nondiaArtifact023Workflow(root, scenario)
        imported <- nondiaArtifact023Import(workflow, scenario)
        expect_true(persistProtImportArtifacts(
            workflow,
            imported$result,
            imported$source_path,
            log_warn = function(...) NULL
        )$ok)
        object <- nondiaArtifact023Design(workflow, imported$workflow_type)
        memory_before <- workflow$state_manager$exportState()
        output <- persistProtDesignArtifacts(
            workflow,
            log_warn = function(...) NULL
        )
        expect_true(output$enabled, info = scenario$scenario_id)
        expect_true(output$ok, info = scenario$scenario_id)
        expect_true(output$committed, info = scenario$scenario_id)
        if (scenario$format %in% c("maxquant", "fragpipe")) {
            expect_true(output$process_evidence$distinct_workers)
            expect_true(output$process_evidence$valid_s4)
            expect_false(output$process_evidence$complete_payload_returned)
        } else {
            expect_null(output$process_evidence)
        }
        expect_setequal(names(output$refs), c(
            "cleaned_data", "design_matrix", "contrasts", "args",
            "annotations", "sequences"
        ))
        expect_identical(workflow$state_manager$exportState(), memory_before)
        expect_identical(
            workflow$state_manager$getState("protein_s4_initial"),
            object
        )
        expect_identical(
            nondiaArtifact023ReadRef(
                workflow$workflow_context,
                output$refs$design_matrix
            ),
            workflow$design_matrix
        )
        descriptor <- artifactProteomicsNonDiaWorkflowDescriptor(
            scenario$capability_id
        )
        manager <- newWorkflowState(
            workflow_context = workflow$workflow_context,
            workflow_descriptor = descriptor,
            descriptor_catalogue = artifactWorkflowDescriptorCatalogue(),
            codec_catalogue = artifactS4CodecCatalogue()
        )
        expect_identical(manager$getState("protein_s4_initial"), object)
        expect_identical(manager$getCurrentStateName(), "protein_s4_initial")
        expect_true(manager$close())
    }
})

test_that("non-DIA write and hydration failures preserve current memory state", {
    nondiaArtifact023SkipDependencies()
    for (scenario in nondiaArtifact023Scenarios()) {
        for (failure_stage in c("before_write", "before_registry_commit")) {
            root <- withr::local_tempdir()
            workflow <- nondiaArtifact023Workflow(root, scenario)
            imported <- nondiaArtifact023Import(workflow, scenario)
            memory_before <- workflow$data_tbl
            output <- persistProtImportArtifacts(
                workflow,
                imported$result,
                imported$source_path,
                failure_injector = function(stage, context) {
                    if (identical(stage, failure_stage)) {
                        rlang::abort("injected non-DIA import failure")
                    }
                    invisible(context)
                },
                log_warn = function(...) NULL
            )
            expect_true(output$enabled)
            expect_false(output$ok)
            expect_identical(workflow$data_tbl, memory_before)
            expect_identical(workflow$data_cln, memory_before)
            expect_identical(
                workflow$state_manager$getWorkflowType(),
                imported$workflow_type
            )
        }

        root <- withr::local_tempdir()
        workflow <- nondiaArtifact023Workflow(root, scenario)
        imported <- nondiaArtifact023Import(workflow, scenario)
        expect_true(persistProtImportArtifacts(
            workflow,
            imported$result,
            imported$source_path,
            log_warn = function(...) NULL
        )$ok)
        object <- nondiaArtifact023Design(workflow, imported$workflow_type)
        memory_before <- workflow$state_manager$exportState()
        output <- persistProtNonDiaDesignArtifacts(
            workflow,
            manager_factory = nondiaArtifact023CorruptingManager,
            log_warn = function(...) NULL
        )
        expect_false(output$ok, info = scenario$scenario_id)
        expect_identical(workflow$state_manager$exportState(), memory_before)
        expect_identical(
            workflow$state_manager$getState("protein_s4_initial"),
            object
        )

        codec_root <- withr::local_tempdir()
        codec_workflow <- nondiaArtifact023Workflow(codec_root, scenario)
        codec_import <- nondiaArtifact023Import(codec_workflow, scenario)
        expect_true(persistProtImportArtifacts(
            codec_workflow,
            codec_import$result,
            codec_import$source_path,
            log_warn = function(...) NULL
        )$ok)
        invalid <- nondiaArtifact023Design(
            codec_workflow,
            codec_import$workflow_type
        )
        invalid@args$globalParameters$workflow_type <- "INVALID"
        codec_workflow$state_manager$states$protein_s4_initial$data <- invalid
        memory_before <- codec_workflow$state_manager$exportState()
        output <- persistProtNonDiaDesignArtifacts(
            codec_workflow,
            log_warn = function(...) NULL
        )
        expect_false(output$ok, info = scenario$scenario_id)
        expect_true(
            "multischolar_proteomics_codec_shape_mismatch" %in%
                output$error_class
        )
        expect_identical(codec_workflow$state_manager$exportState(), memory_before)
    }
})

test_that("non-DIA artifacts use coordinator data instead of legacy globals", {
    nondiaArtifact023SkipDependencies()
    scenario <- nondiaArtifact023Scenarios()[[1L]]
    root <- withr::local_tempdir()
    workflow <- nondiaArtifact023Workflow(root, scenario)
    imported <- nondiaArtifact023Import(workflow, scenario)
    expect_true(persistProtImportArtifacts(
        workflow,
        imported$result,
        imported$source_path,
        log_warn = function(...) NULL
    )$ok)

    global_names <- c(
        "config_list", "contrasts_tbl", "uniprot_dat_cln", "aa_seq_tbl_final"
    )
    had_global <- vapply(
        global_names,
        exists,
        logical(1),
        envir = .GlobalEnv,
        inherits = FALSE
    )
    previous <- mget(
        global_names[had_global],
        envir = .GlobalEnv,
        inherits = FALSE
    )
    withr::defer({
        present <- global_names[vapply(
            global_names,
            exists,
            logical(1),
            envir = .GlobalEnv,
            inherits = FALSE
        )]
        if (length(present) > 0L) rm(list = present, envir = .GlobalEnv)
        list2env(previous, envir = .GlobalEnv)
    })
    poison <- stats::setNames(
        lapply(global_names, function(name) list(poison = name)),
        global_names
    )
    list2env(poison, envir = .GlobalEnv)

    nondiaArtifact023Design(workflow, imported$workflow_type)
    output <- persistProtDesignArtifacts(
        workflow,
        log_warn = function(...) NULL
    )
    expect_true(output$ok)
    expect_identical(
        nondiaArtifact023ReadRef(
            workflow$workflow_context,
            output$refs$contrasts
        ),
        workflow$contrasts_tbl
    )
    expect_identical(
        nondiaArtifact023ReadRef(
            workflow$workflow_context,
            output$refs$annotations
        ),
        workflow$uniprot_dat_cln
    )
    expect_identical(
        nondiaArtifact023ReadRef(
            workflow$workflow_context,
            output$refs$sequences
        ),
        workflow$aa_seq_tbl_final
    )
    for (name in global_names) {
        expect_identical(get(name, envir = .GlobalEnv), poison[[name]])
    }
})

test_that("memory, blocked, cross-omic, and cross-tuple paths create no artifacts", {
    nondiaArtifact023SkipDependencies()
    scenario <- nondiaArtifact023Scenarios()[[1L]]
    root <- withr::local_tempdir()
    memory <- nondiaArtifact023Workflow(root, scenario, backend = "memory")
    imported <- nondiaArtifact023Import(memory, scenario)
    output <- persistProtImportArtifacts(
        memory,
        imported$result,
        imported$source_path,
        log_warn = function(...) NULL
    )
    expect_false(output$enabled)
    expect_true(output$ok)
    expect_false(dir.exists(file.path(root, "state")))

    for (blocked in list(
        list(format = "spectronaut", data_type = "protein"),
        list(format = "spectronaut", data_type = "peptide"),
        list(format = "maxquant", data_type = "peptide")
    )) {
        blocked_root <- withr::local_tempdir()
        workflow <- nondiaArtifact023Workflow(blocked_root, scenario)
        workflow$data_format <- blocked$format
        workflow$data_type <- blocked$data_type
        prepared <- prepareProtNonDiaArtifactContext(workflow)
        expect_false(prepared$enabled)
        expect_false(workflow$workflow_context$isBound())
        expect_false(dir.exists(file.path(blocked_root, "state")))
    }

    other_root <- withr::local_tempdir()
    other <- nondiaArtifact023Workflow(
        other_root,
        scenario,
        omic_type = "metabolomics"
    )
    other$data_format <- "maxquant"
    other$data_type <- "protein"
    prepared <- prepareProtNonDiaArtifactContext(other)
    expect_false(prepared$enabled)
    expect_identical(prepared$reason, "not_proteomics_context")
    expect_false(other$workflow_context$isBound())
    expect_false(dir.exists(file.path(other_root, "state")))

    cross_root <- withr::local_tempdir()
    cross <- nondiaArtifact023Workflow(cross_root, scenario)
    nondiaArtifact023Import(cross, scenario)
    expect_true(prepareProtNonDiaArtifactContext(cross)$enabled)
    cross$data_format <- "fragpipe"
    prepared <- prepareProtNonDiaArtifactContext(cross)
    expect_false(prepared$enabled)
    expect_identical(
        cross$workflow_context$getIdentity()$input_format,
        "maxquant"
    )
})

test_that("production defaults dispatch exact DIA and non-DIA artifact adapters", {
    expect_identical(
        formals(applyProtImportResultToWorkflow)$prepareArtifactContext,
        quote(prepareProtArtifactContextSafely)
    )
    expect_identical(
        formals(registerProtImportObservers)$persistImportArtifactsFn,
        quote(persistProtImportArtifacts)
    )
    expect_identical(
        formals(completeProtDesignPostCheckpoint)$persistArtifactFn,
        quote(persistProtDesignArtifacts)
    )
})
