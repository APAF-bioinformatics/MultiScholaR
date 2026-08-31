diaArtifact009SkipDependencies <- function() {
    for (package in c("arrow", "DBI", "duckdb", "filelock")) {
        testthat::skip_if_not_installed(package)
    }
}

diaArtifact009Paths <- function(root, project_id = "dia-artifact-009") {
    paths <- list(
        base_dir = root,
        project_id = project_id,
        omic_type = "proteomics",
        omic_label = "dia-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    paths
}

diaArtifact009Workflow <- function(
    root,
    backend = "artifact",
    rollout = "dual_write"
) {
    paths <- diaArtifact009Paths(root)
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        "proteomics",
        "dia-study",
        storage_policy = list(
            requested_backend = backend,
            requested_rollout = rollout,
            migration_requested = identical(backend, "artifact"),
            project_id = paths$project_id
        )
    )
    workflow$state_manager <- WorkflowState$new()
    workflow$artifact_stage_results <- NULL
    workflow$tab_status <- list(setup_import = "complete", design_matrix = "pending")
    workflow$processing_log <- list()
    workflow
}

diaArtifact009Import <- function(workflow, source_path) {
    result <- suppressMessages(importDIANNData(source_path))
    workflow$data_tbl <- result$data
    workflow$data_cln <- result$data
    workflow$data_format <- "diann"
    workflow$data_type <- result$data_type
    workflow$column_mapping <- result$column_mapping
    workflow$state_manager$setWorkflowType("DIA")
    list(result = result, memory = result$data)
}

diaArtifact009Design <- function(workflow) {
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
    object <- PeptideQuantitativeDataDiann(
        peptide_data = workflow$data_cln,
        design_matrix = workflow$design_matrix,
        technical_replicate_id = "tech_rep_group",
        args = workflow$config_list
    )
    workflow$state_manager$saveState(
        "raw_data_s4",
        object,
        workflow$config_list,
        "DIA-NN design memory checkpoint."
    )
    object
}

diaArtifact009ReadRef <- function(context, ref) {
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

diaArtifact009Registry <- function(context) {
    registry <- projectRegistryForContext(context)
    session <- openProjectRegistryReadOnly(registry)
    withr::defer(closeProjectRegistry(session), envir = parent.frame())
    session
}

test_that("DIA-NN TSV and Parquet imports have one portable scientific shape", {
    diaArtifact009SkipDependencies()
    source <- testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    tsv <- suppressMessages(importDIANNData(source))
    parquet_path <- tempfile(fileext = ".parquet")
    withr::defer(unlink(parquet_path, force = TRUE))
    arrow::write_parquet(vroom::vroom(source, show_col_types = FALSE), parquet_path)
    parquet <- suppressMessages(importDIANNData(parquet_path))

    expect_identical(tsv, parquet)
    expect_s3_class(tsv$data, "tbl_df")
    expect_false(inherits(tsv$data, "spec_tbl_df"))
    expect_identical(
        setdiff(names(attributes(tsv$data)), c("names", "row.names", "class")),
        character()
    )
})

test_that("DIA evict ingress publishes refs without parent table hydration", {
    diaArtifact009SkipDependencies()
    root <- withr::local_tempdir(pattern = "dia-artifact-071-ingress-")
    source <- testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    workflow <- diaArtifact009Workflow(root, rollout = "evict")
    memory <- suppressMessages(importDIANNData(source))

    staged <- stageProtDiaImportArtifacts(
        workflow,
        source,
        "diann",
        hydrate_parent = NULL
    )
    expect_true(staged$ok)
    expect_true(staged$pending_stage$bounded_streaming)
    expect_identical(
        staged$pending_stage$stage$refs$canonical_data$codec$id,
        "multischolar.rectangular_streaming"
    )
    expect_false(staged$preflight$complete_payload_materialized)
    expect_identical(
        staged$preflight$source_size_bytes,
        as.numeric(file.info(source)$size)
    )
    expect_true(staged$result$parent_hydration_deferred)
    expect_null(staged$result$data)
    summary <- protDiaValidateImportSummary(staged$result$import_summary)
    expect_identical(summary$rows, as.numeric(nrow(memory$data)))
    expect_identical(summary$columns, as.integer(ncol(memory$data)))
    expect_identical(
        unlist(summary$column_names, use.names = FALSE),
        names(memory$data)
    )
    expect_setequal(
        unlist(summary$run_values, use.names = FALSE),
        unique(memory$data$Run)
    )

    applyProtImportResultToWorkflow(
        workflow,
        staged$result,
        "diann",
        fastaPath = "fixture.fasta",
        prepareArtifactContext = function(...) invisible(NULL)
    )
    expect_null(workflow$data_tbl)
    expect_null(workflow$data_cln)
    expect_identical(workflow$artifact_import_summary, summary)
    published <- persistProtDiaImportArtifacts(
        workflow,
        staged$result,
        source,
        pending_stage = staged$pending_stage,
        worker_attempted = TRUE
    )
    expect_true(published$committed)
    expect_true(protDiaWorkflowPayloadAvailable(workflow, "data_tbl"))
    expect_identical(
        resolveProtDiaWorkflowTable(workflow, "data_tbl"),
        memory$data
    )
})

test_that("DIA ingress rejects invalid schema before worker staging", {
    diaArtifact009SkipDependencies()
    source <- tempfile(fileext = ".tsv")
    withr::defer(unlink(source, force = TRUE))
    writeLines("Protein.Group\tRun\tPrecursor.Quantity", source)

    expect_error(
        protDiaIngressPreflight(source),
        class = "multischolar_invalid_prot_dia_ingress_schema"
    )
})

test_that("deferred DIA design runs in a payload-isolated worker", {
    diaArtifact009SkipDependencies()
    root <- withr::local_tempdir(pattern = "dia-artifact-072-design-worker-")
    paths <- diaArtifact009Paths(root)
    source <- testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    workflow <- diaArtifact009Workflow(root, rollout = "evict")
    staged <- stageProtDiaImportArtifacts(workflow, source, "diann")
    applyProtImportResultToWorkflow(
        workflow,
        staged$result,
        "diann",
        "fixture.fasta",
        prepareArtifactContext = function(...) invisible(NULL)
    )
    expect_true(persistProtDiaImportArtifacts(
        workflow,
        staged$result,
        source,
        pending_stage = staged$pending_stage,
        worker_attempted = TRUE
    )$committed)
    workflow$config_list <- list(
        globalParameters = list(workflow_type = "DIA"),
        srlQvalueProteotypicPeptideClean = list(
            input_matrix_column_ids = unlist(
                staged$result$import_summary$column_names,
                use.names = FALSE
            )
        )
    )
    input <- protDiaDeferredDesignInput(workflow)
    runs <- input$Run
    results <- list(
        data_cln = input,
        design_matrix = data.frame(
            Run = runs,
            group = sub("_.*", "", runs),
            replicates = seq_along(runs),
            tech_rep_group = runs,
            stringsAsFactors = FALSE
        ),
        contrasts_tbl = data.frame(
            contrasts = "groupKO-groupWT",
            friendly_names = "KO_vs_WT",
            full_format = "KO_vs_WT=groupKO-groupWT",
            stringsAsFactors = FALSE
        ),
        config_list = workflow$config_list
    )
    worker <- runProtDiaDeferredDesignProcess(protDiaDeferredDesignSpec(
        workflow,
        results,
        paths,
        "dia-study"
    ))
    applyProtDiaDeferredDesignResult(workflow, worker)

    expect_false(identical(worker$worker_pid, as.integer(Sys.getpid())))
    expect_false(worker$complete_payload_returned)
    expect_null(workflow$data_tbl)
    expect_null(workflow$data_cln)
    expect_s3_class(workflow$state_manager, "ArtifactWorkflowState")
    expect_identical(workflow$state_manager$getCacheInfo()$entries, 0L)
    qc <- runPeptideQvalueApplyStep(workflow, 0.01, 0.01, TRUE)
    expect_true(nzchar(qc$resultText))
    expect_true(is.raw(qc$plot_png))
    expect_gt(length(qc$plot_png), 0L)
    expect_identical(workflow$state_manager$getCurrentStateName(), "qvalue_filtered")
    expect_identical(workflow$state_manager$getCacheInfo()$entries, 0L)

    peptide <- workflow$state_manager$getState()
    protein_groups <- c("PG1", "PG1", "PG2")
    protein_table <- data.frame(
        Protein.Group = protein_groups,
        matrix(
            seq_len(length(protein_groups) * length(runs)),
            nrow = length(protein_groups),
            dimnames = list(NULL, runs)
        ),
        check.names = FALSE
    )
    protein <- ProteinQuantitativeData(
        protein_quant_table = protein_table,
        protein_id_column = "Protein.Group",
        protein_id_table = .proteinIdTableFromPeptideLineage(
            peptide,
            protein_groups
        ),
        design_matrix = peptide@design_matrix,
        sample_id = peptide@sample_id,
        group_id = peptide@group_id,
        technical_replicate_id = peptide@technical_replicate_id,
        args = peptide@args
    )
    protein <- saveProtProteinQcState(
        workflow,
        workflow$state_manager,
        peptide,
        protein,
        "protein_rollup",
        "protein_s4_created",
        list(rollup_method = "test_fixture"),
        "Created test protein state",
        transformation_type = "aggregation"
    )
    previous_manager <- workflow$state_manager
    state_export <- previous_manager$exportState()
    workflow$state_manager <- newProtDiaSettledStateManager(
        workflow,
        state_export = state_export
    )
    previous_manager$close()

    deduplicated <- runProteinDuplicateRemovalStep(workflow, "mean")
    expect_true(nzchar(deduplicated$resultText))
    expect_true(is.raw(deduplicated$plot_png))
    expect_false("deduplicatedS4" %in% names(deduplicated))
    expect_identical(
        workflow$state_manager$getCurrentStateName(),
        "duplicates_removed"
    )
    expect_identical(workflow$state_manager$getCacheInfo()$entries, 0L)
})

test_that("DIA-NN import dual-write records provenance without owning source paths", {
    diaArtifact009SkipDependencies()
    root <- tempfile("dia-artifact-009-import-")
    dir.create(root)
    withr::defer(unlink(root, recursive = TRUE, force = TRUE))
    source <- file.path(root, "input.tsv")
    fixture <- testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    expect_true(file.copy(fixture, source))
    workflow <- diaArtifact009Workflow(root)
    imported <- diaArtifact009Import(workflow, source)
    memory_before <- imported$memory

    output <- persistProtDiaImportArtifacts(
        workflow,
        imported$result,
        source,
        use_precursor_norm = TRUE,
        sanitize_names = FALSE
    )
    expect_true(output$enabled)
    expect_true(output$ok)
    expect_true(output$committed)
    expect_identical(workflow$data_tbl, memory_before)
    expect_identical(
        diaArtifact009ReadRef(workflow$workflow_context, output$refs$canonical_data),
        memory_before
    )

    session <- diaArtifact009Registry(workflow$workflow_context)
    identity <- workflow$workflow_context$getIdentity()
    sources <- projectRegistryQuery(
        session,
        "sources",
        filters = list(workflow_id = identity$workflow_id, run_id = output$run_id)
    )
    parameters <- projectRegistryQuery(
        session,
        "parameters",
        filters = list(workflow_id = identity$workflow_id, run_id = output$run_id)
    )
    artifacts <- projectRegistryQuery(
        session,
        "artifacts",
        filters = list(workflow_id = identity$workflow_id, stage_id = "import")
    )
    artifacts <- artifacts[artifacts$run_id == output$run_id, , drop = FALSE]
    expect_identical(sources$source_digest, artifactByteDigest(source))
    expect_true(is.na(sources$source_uri))
    expect_identical(sources$parser_id, "MultiScholaR::importDIANNData")
    expect_identical(sources$parser_version, "1.0.0")
    expect_identical(sources$format_id, "diann.tsv")
    expect_identical(sources$data_level, "peptide")
    expect_setequal(
        parameters$parameter_name,
        c(
            "column_mapping", "column_mapping_serialized",
            "readthrough_contract_version", "use_precursor_norm",
            "sanitize_names"
        )
    )
    expect_identical(artifacts$state_role, "canonical_data")
    expect_identical(artifacts$status, "committed")
    closeProjectRegistry(session)
    expect_true(unlink(source, force = TRUE) == 0L)
    expect_identical(
        diaArtifact009ReadRef(workflow$workflow_context, output$refs$canonical_data),
        memory_before
    )
})

test_that("DIA-NN design dual-write captures refs and exact S4 independently", {
    diaArtifact009SkipDependencies()
    root <- tempfile("dia-artifact-009-design-")
    dir.create(root)
    withr::defer(unlink(root, recursive = TRUE, force = TRUE))
    source <- testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    workflow <- diaArtifact009Workflow(root)
    imported <- diaArtifact009Import(workflow, source)
    expect_true(persistProtDiaImportArtifacts(
        workflow,
        imported$result,
        source
    )$ok)
    object <- diaArtifact009Design(workflow)
    state_before <- workflow$state_manager$exportState()

    output <- persistProtDiaDesignArtifacts(workflow)
    expect_true(output$enabled)
    expect_true(output$ok)
    expect_true(output$committed)
    expect_setequal(names(output$refs), c(
        "cleaned_data", "design_matrix", "contrasts", "args",
        "annotations", "sequences"
    ))
    expect_identical(workflow$state_manager$exportState(), state_before)
    expect_identical(workflow$state_manager$getState("raw_data_s4"), object)
    expect_identical(output$state_manifest$current_state, "raw_data_s4")
    expect_true(length(output$state_metadata$artifact_refs) > 0L)
    verification <- output$hydration_verification
    expect_true(verification$valid)
    expect_identical(verification$mode, "inline_exact")
    expect_true(identical(
        verification$verifier_pid,
        as.integer(Sys.getpid())
    ))
    expect_identical(
        verification$expected_digest,
        verification$hydrated_digest
    )
    expect_match(verification$manifest_semantic_digest, "^[0-9a-f]{64}$")
    expect_false(verification$complete_payload_returned)
    expect_identical(
        diaArtifact009ReadRef(workflow$workflow_context, output$refs$cleaned_data),
        workflow$data_cln
    )
    expect_identical(
        diaArtifact009ReadRef(workflow$workflow_context, output$refs$design_matrix),
        workflow$design_matrix
    )
    expect_identical(
        diaArtifact009ReadRef(workflow$workflow_context, output$refs$contrasts),
        workflow$contrasts_tbl
    )
    expect_identical(
        diaArtifact009ReadRef(workflow$workflow_context, output$refs$annotations),
        workflow$uniprot_dat_cln
    )
    expect_identical(
        diaArtifact009ReadRef(workflow$workflow_context, output$refs$sequences),
        workflow$aa_seq_tbl_final
    )
})

test_that("DIA exact S4 proof is independent of ALTREP representation", {
    diaArtifact009SkipDependencies()
    root <- withr::local_tempdir(pattern = "dia-artifact-009-digest-")
    workflow <- diaArtifact009Workflow(root)
    imported <- diaArtifact009Import(
        workflow,
        testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    )
    object <- diaArtifact009Design(workflow)
    hydrated <- hydrateDiaS4Artifact(dehydrateDiaS4Artifact(object))

    expect_identical(hydrated, object)
    expect_identical(
        artifactExactS4HydrationDigest(hydrated),
        artifactExactS4HydrationDigest(object)
    )
    expect_identical(imported$memory, workflow$data_cln)
})

test_that("artifact-mode import and design state cannot be overridden by globals", {
    diaArtifact009SkipDependencies()
    root <- tempfile("dia-artifact-009-globals-")
    dir.create(root)
    withr::defer(unlink(root, recursive = TRUE, force = TRUE))
    source <- testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    workflow <- diaArtifact009Workflow(root)
    imported <- diaArtifact009Import(workflow, source)
    prepareProtDiaArtifactContext(workflow)
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
        lapply(global_names, \(name) list(poison = name)),
        global_names
    )
    list2env(poison, envir = .GlobalEnv)
    assignments <- 0L
    storeProtImportConfiguration(
        workflow,
        configList = list(globalParameters = list(workflow_type = "DIA")),
        assignFn = \(...) assignments <<- assignments + 1L,
        assignEnv = .GlobalEnv,
        logInfo = \(...) NULL
    )
    expect_identical(assignments, 0L)
    expect_identical(workflow$config_list$globalParameters$workflow_type, "DIA")
    expect_identical(workflow$data_tbl, imported$memory)
    diaArtifact009Design(workflow)
    output <- persistProtDiaDesignArtifacts(workflow)
    expect_true(output$ok)
    expect_identical(
        diaArtifact009ReadRef(workflow$workflow_context, output$refs$contrasts),
        workflow$contrasts_tbl
    )
    expect_identical(
        diaArtifact009ReadRef(workflow$workflow_context, output$refs$annotations),
        workflow$uniprot_dat_cln
    )
    expect_identical(
        diaArtifact009ReadRef(workflow$workflow_context, output$refs$sequences),
        workflow$aa_seq_tbl_final
    )
    for (name in global_names) {
        expect_identical(get(name, envir = .GlobalEnv), poison[[name]])
    }
})

test_that("store and registry failures leave successful memory import current", {
    diaArtifact009SkipDependencies()
    source <- testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    for (failure_stage in c("before_write", "before_registry_commit")) {
        root <- tempfile(paste0("dia-artifact-009-failure-", failure_stage, "-"))
        dir.create(root)
        withr::defer(unlink(root, recursive = TRUE, force = TRUE))
        workflow <- diaArtifact009Workflow(root)
        imported <- diaArtifact009Import(workflow, source)
        memory_before <- imported$memory
        output <- persistProtDiaImportArtifacts(
            workflow,
            imported$result,
            source,
            failure_injector = \(stage, context) {
                if (identical(stage, failure_stage)) {
                    rlang::abort("injected import artifact failure")
                }
                invisible(context)
            },
            log_warn = \(...) NULL
        )
        expect_true(output$enabled)
        expect_false(output$ok)
        expect_identical(workflow$data_tbl, memory_before)
        expect_identical(workflow$data_cln, memory_before)
        expect_identical(workflow$state_manager$getWorkflowType(), "DIA")
        expect_identical(workflow$tab_status$setup_import, "complete")
    }
})

test_that("failed DIA-NN workers are not retried by in-process persistence", {
    diaArtifact009SkipDependencies()
    root <- withr::local_tempdir()
    source <- testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    workflow <- diaArtifact009Workflow(root)
    imported <- diaArtifact009Import(workflow, source)
    output <- persistProtDiaImportArtifacts(
        workflow,
        imported$result,
        source,
        pending_stage = NULL,
        worker_attempted = TRUE,
        log_warn = \(...) invisible(NULL)
    )

    expect_true(output$enabled)
    expect_false(output$ok)
    expect_false(output$committed)
    expect_identical(output$reason, "artifact_worker_failed_no_retry")
    expect_length(list.files(root, pattern = "[.]parquet$", recursive = TRUE), 0L)
})

test_that("failed design artifact state cannot replace its memory checkpoint", {
    diaArtifact009SkipDependencies()
    root <- tempfile("dia-artifact-009-design-failure-")
    dir.create(root)
    withr::defer(unlink(root, recursive = TRUE, force = TRUE))
    source <- testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    workflow <- diaArtifact009Workflow(root)
    imported <- diaArtifact009Import(workflow, source)
    expect_true(persistProtDiaImportArtifacts(
        workflow,
        imported$result,
        source
    )$ok)
    object <- diaArtifact009Design(workflow)
    state_before <- workflow$state_manager$exportState()
    output <- persistProtDiaDesignArtifacts(
        workflow,
        save_state_fn = \(...) rlang::abort("injected state failure"),
        log_warn = \(...) NULL
    )
    expect_false(output$ok)
    expect_identical(workflow$state_manager$exportState(), state_before)
    expect_identical(workflow$state_manager$getState("raw_data_s4"), object)

    session <- diaArtifact009Registry(workflow$workflow_context)
    identity <- workflow$workflow_context$getIdentity()
    runs <- projectRegistryQuery(
        session,
        "runs",
        filters = list(workflow_id = identity$workflow_id, status = "failed")
    )
    expect_identical(nrow(runs), 1L)
})

test_that("independent hydration failure precedes artifact state registration", {
    diaArtifact009SkipDependencies()
    root <- tempfile("dia-artifact-009-hydration-failure-")
    dir.create(root)
    withr::defer(unlink(root, recursive = TRUE, force = TRUE))
    workflow <- diaArtifact009Workflow(root)
    imported <- diaArtifact009Import(
        workflow,
        testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    )
    prepareProtDiaArtifactContext(workflow)
    object <- diaArtifact009Design(workflow)
    manager <- ArtifactWorkflowState$new(
        workflow$workflow_context,
        hydrate_fn = \(bundle) {
            hydrated <- hydrateDiaS4Artifact(bundle)
            hydrated@args$changed_after_hydration <- TRUE
            hydrated
        }
    )
    withr::defer(manager$close())
    initial <- manager$exportState()
    expect_error(
        manager$saveState(
            "raw_data_s4",
            object,
            workflow$config_list,
            "must not register"
        ),
        class = "multischolar_inexact_artifact_state_hydration"
    )
    expect_identical(manager$exportState(), initial)
    expect_false(manager$hasState("raw_data_s4"))
})

test_that("DIA worker proof failure precedes state registration", {
    diaArtifact009SkipDependencies()
    root <- withr::local_tempdir(pattern = "dia-artifact-009-worker-failure-")
    source <- testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    workflow <- diaArtifact009Workflow(root)
    imported <- diaArtifact009Import(workflow, source)
    expect_true(persistProtDiaImportArtifacts(
        workflow,
        imported$result,
        source
    )$ok)
    object <- diaArtifact009Design(workflow)
    state_before <- workflow$state_manager$exportState()
    failing_factory <- function(..., verify_hydration_fn) {
        newWorkflowState(
            ...,
            verify_hydration_fn = function(
                store,
                manifest,
                expected_object,
                hydrate_fn
            ) {
                if (is.null(expected_object)) {
                    return(verifyProtDiaArtifactStateInWorker(
                        store,
                        manifest,
                        expected_object,
                        hydrate_fn
                    ))
                }
                rlang::abort("injected process-bound parity failure")
            }
        )
    }
    output <- persistProtDiaDesignArtifacts(
        workflow,
        manager_factory = failing_factory,
        log_warn = \(...) invisible(NULL)
    )

    expect_false(output$ok)
    expect_identical(workflow$state_manager$exportState(), state_before)
    expect_identical(workflow$state_manager$getState("raw_data_s4"), object)
    session <- diaArtifact009Registry(workflow$workflow_context)
    identity <- workflow$workflow_context$getIdentity()
    states <- projectRegistryQuery(
        session,
        "states",
        filters = list(workflow_id = identity$workflow_id)
    )
    expect_false(any(states$logical_name == "raw_data_s4"))
})

test_that("memory and non-DIA lanes cannot initialize the DIA canary registry", {
    diaArtifact009SkipDependencies()
    source <- testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    root <- tempfile("dia-artifact-009-memory-")
    dir.create(root)
    withr::defer(unlink(root, recursive = TRUE, force = TRUE))
    workflow <- diaArtifact009Workflow(root, backend = "memory")
    imported <- diaArtifact009Import(workflow, source)
    output <- persistProtDiaImportArtifacts(workflow, imported$result, source)
    expect_false(output$enabled)
    expect_true(output$ok)
    expect_false(dir.exists(file.path(root, "state")))

    descriptor_ids <- names(artifactWorkflowDescriptorCatalogue()$descriptors)
    expect_setequal(descriptor_ids, c(
        "proteomics.diann.peptide.dia.v1",
        "proteomics.maxquant.protein.lfq.v1",
        "proteomics.fragpipe.protein.lfq.v1",
        "proteomics.pd_tmt.protein.tmt.v1",
        "metabolomics.custom.metabolite.standard.v1",
        "lipidomics.lipidsearch.lipid.standard.v1"
    ))
    expect_false(any(grepl("spectronaut", descriptor_ids, fixed = TRUE)))

    unsupported_root <- tempfile("dia-artifact-009-unsupported-")
    dir.create(unsupported_root)
    withr::defer(unlink(unsupported_root, recursive = TRUE, force = TRUE))
    unsupported <- diaArtifact009Workflow(unsupported_root)
    unsupported$data_format <- "fragpipe"
    unsupported$data_type <- "protein"
    assignments <- 0L
    expect_false(protDiaArtifactCoordinatorOwned(unsupported))
    expect_false(prepareProtDiaArtifactContext(unsupported)$enabled)
    storeProtImportConfiguration(
        unsupported,
        configList = list(globalParameters = list(workflow_type = "LFQ")),
        assignFn = \(...) assignments <<- assignments + 1L,
        assignEnv = new.env(parent = emptyenv()),
        logInfo = \(...) NULL
    )
    expect_identical(assignments, 1L)
    expect_false(dir.exists(file.path(unsupported_root, "state")))
})

test_that("DIA canary helpers are collated once before import and design callers", {
    description <- read.dcf(testthat::test_path("..", "..", "DESCRIPTION"))
    collate <- strsplit(description[[1L, "Collate"]], "[[:space:]]+")[[1L]]
    collate <- gsub("^'|'$", "", collate)
    import_helper <- match("mod_prot_import_artifact_helpers.R", collate)
    design_helper <- match("mod_prot_design_artifact_helpers.R", collate)
    expect_identical(sum(collate == "mod_prot_import_artifact_helpers.R"), 1L)
    expect_identical(sum(collate == "mod_prot_design_artifact_helpers.R"), 1L)
    expect_lt(import_helper, design_helper)
    expect_lt(design_helper, match("mod_prot_design_state_helpers.R", collate))
    expect_lt(import_helper, match("mod_prot_import_processing_helpers.R", collate))
})
