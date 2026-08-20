diaArtifact010SkipDependencies <- function() {
    for (package in c("arrow", "DBI", "duckdb", "filelock")) {
        testthat::skip_if_not_installed(package)
    }
}

diaArtifact010Paths <- function(root, project_id = "dia-artifact-010") {
    paths <- list(
        base_dir = root,
        project_id = project_id,
        omic_type = "proteomics",
        omic_label = "dia-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE)
    dir.create(paths$results_dir, recursive = TRUE)
    paths
}

diaArtifact010Workflow <- function(paths) {
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        "proteomics",
        "dia-study",
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

diaArtifact010Build <- function(root, link_import = TRUE) {
    paths <- diaArtifact010Paths(root)
    source <- file.path(paths$source_dir, "report.tsv")
    fixture <- testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    stopifnot(file.copy(fixture, source))
    workflow <- diaArtifact010Workflow(paths)
    imported <- suppressMessages(importDIANNData(source))
    workflow$data_tbl <- imported$data
    workflow$data_cln <- imported$data
    workflow$data_format <- "diann"
    workflow$data_type <- imported$data_type
    workflow$column_mapping <- imported$column_mapping
    workflow$state_manager$setWorkflowType("DIA")
    import_result <- persistProtDiaImportArtifacts(
        workflow,
        imported,
        source
    )
    stopifnot(isTRUE(import_result$ok))
    import_oracle <- workflow$data_tbl
    if (!isTRUE(link_import)) workflow$artifact_stage_results$import <- NULL

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
        source = source,
        workflow = workflow,
        object = object,
        import_oracle = import_oracle,
        import_result = import_result,
        design_result = design_result
    )
}

diaArtifact010BlankWorkflow <- function() {
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

expectDiaArtifact010StateExact <- function(expected, actual) {
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

test_that("DIA read-through restores one parent-linked scientific workflow", {
    diaArtifact010SkipDependencies()
    root <- withr::local_tempdir(pattern = "dia-artifact-010-resume-")
    built <- diaArtifact010Build(root)
    memory_oracle <- list(
        data_tbl = built$workflow$data_tbl,
        data_cln = built$workflow$data_cln,
        design_matrix = built$workflow$design_matrix,
        contrasts_tbl = built$workflow$contrasts_tbl,
        config_list = built$workflow$config_list,
        column_mapping = built$workflow$column_mapping,
        uniprot_dat_cln = built$workflow$uniprot_dat_cln,
        aa_seq_tbl_final = built$workflow$aa_seq_tbl_final
    )
    expect_true(unlink(built$source, force = TRUE) == 0L)

    newer <- built$workflow$data_tbl
    newer$Precursor.Normalised <- newer$Precursor.Normalised + 1000
    built$workflow$data_tbl <- newer
    newer_import <- persistProtDiaImportArtifacts(
        built$workflow,
        list(column_mapping = built$workflow$column_mapping),
        testthat::test_path("..", "testdata", "e2e", "prot_dia", "report.tsv")
    )
    expect_true(newer_import$ok)
    resumed <- diaArtifact010BlankWorkflow()
    result <- resumeProtDiaArtifactWorkflow(
        resumed,
        built$paths,
        "dia-study"
    )
    withr::defer(resumed$state_manager$close())

    expect_true(result$resumed)
    expect_identical(result$import_run_id, built$import_result$run_id)
    expect_false(identical(result$import_run_id, newer_import$run_id))
    for (name in names(memory_oracle)) {
        expect_identical(resumed[[name]], memory_oracle[[name]], info = name)
    }
    expect_s3_class(resumed$state_manager, "ArtifactWorkflowState")
    expectDiaArtifact010StateExact(built$object, resumed$state_manager$getState())
    expect_identical(resumed$tab_status$setup_import, "complete")
    expect_identical(resumed$tab_status$design_matrix, "complete")
    expect_identical(resumed$tab_status$quality_control, "pending")

    preview <- previewProtDiaImportArtifact(
        resumed$workflow_context,
        result$import_ref,
        limit = 5L
    )
    descriptor <- artifactDiaWorkflowDescriptor()
    operation <- descriptor$queries[[1L]]
    expected <- memory_oracle$data_tbl[
        do.call(order, memory_oracle$data_tbl[operation$order_by]),
        operation$projections,
        drop = FALSE
    ]
    expected <- as.data.frame(expected, stringsAsFactors = FALSE)
    rownames(expected) <- NULL
    expect_identical(preview, head(expected, 5L))
    qc_input <- hydrateArtifactStageInput(
        workflow_state = resumed$state_manager
    )
    expectDiaArtifact010StateExact(built$object, qc_input$value)
})

test_that("artifact revert resumes an immutable child and preflights failures", {
    diaArtifact010SkipDependencies()
    root <- withr::local_tempdir(pattern = "dia-artifact-010-lineage-")
    built <- diaArtifact010Build(root)
    resumed <- diaArtifact010BlankWorkflow()
    resumeProtDiaArtifactWorkflow(resumed, built$paths, "dia-study")
    manager <- resumed$state_manager
    child <- manager$getState()
    child@args$qc_child <- TRUE
    manager$saveState("qc_child", child, child@args, "QC child")
    child_generation <- manager$getCurrentGenerationId()
    child_refs <- manager$getStateMetadata()$artifact_refs
    child_paths <- vapply(child_refs, `[[`, character(1), "relative_path")
    parent_name <- "raw_data_s4"

    expectDiaArtifact010StateExact(
        built$object,
        manager$revertToState(parent_name)
    )
    parent_generation <- manager$getCurrentGenerationId()
    expect_identical(manager$states$qc_child$status, "stale")
    expectDiaArtifact010StateExact(child, manager$revertToState("qc_child"))
    expect_identical(manager$getCurrentGenerationId(), child_generation)
    expect_identical(manager$getStateMetadata()$artifact_refs, child_refs)
    expect_true(all(file.exists(file.path(root, child_paths))))
    expect_identical(tail(manager$getEvents(), 1L)[[1L]]$event_type, "state_resumed")

    manager$revertToState(parent_name)
    events_before <- manager$getEvents()
    manager$close()
    failing <- ArtifactWorkflowState$new(
        resumed$workflow_context,
        hydrate_fn = \(...) rlang::abort("injected hydration failure"),
        descriptor_contract = protDiaArtifactDescriptorContract()
    )
    withr::defer(failing$close())
    error <- rlang::catch_cnd(failing$revertToState("qc_child"))
    expect_s3_class(
        error,
        "multischolar_artifact_state_selection_hydration_failed"
    )
    expect_identical(error$project_id, built$paths$project_id)
    expect_identical(error$workflow_id, "proteomics.gui")
    expect_identical(error$generation_id, child_generation)
    expect_identical(error$state_name, "qc_child")
    expect_true(workflowCapabilityScalarString(error$remediation))
    expect_identical(failing$getCurrentGenerationId(), parent_generation)
    expect_identical(failing$getEvents(), events_before)
    expect_identical(failing$getStateMetadata("qc_child")$artifact_refs, child_refs)
})

test_that("moved DIA projects resume in a fresh process without source data", {
    diaArtifact010SkipDependencies()
    original <- withr::local_tempdir(pattern = "dia-artifact-010-original-")
    moved <- withr::local_tempdir(pattern = "dia-artifact-010-moved-")
    built <- diaArtifact010Build(original)
    oracle <- list(data_tbl = built$import_oracle, state = built$object)
    entries <- list.files(original, all.files = TRUE, no.. = TRUE, full.names = TRUE)
    expect_true(all(file.copy(entries, moved, recursive = TRUE)))
    unlink(file.path(moved, "source", "report.tsv"), force = TRUE)
    unlink(original, recursive = TRUE, force = TRUE)
    moved_paths <- list(
        base_dir = moved,
        omic_type = "proteomics",
        omic_label = "dia-study",
        source_dir = file.path(moved, "source"),
        results_dir = file.path(moved, "results")
    )
    expect_false(file.exists(file.path(moved_paths$source_dir, "report.tsv")))
    script <- tempfile("dia-artifact-010-fresh-", fileext = ".R")
    output <- tempfile("dia-artifact-010-fresh-", fileext = ".rds")
    input <- tempfile("dia-artifact-010-fresh-", fileext = ".rds")
    saveRDS(moved_paths, input)
    withr::defer(unlink(c(script, output, input), force = TRUE))
    repository <- normalizePath(testthat::test_path("..", ".."), winslash = "/")
    writeLines(c(
        sprintf("devtools::load_all(%s, quiet = TRUE)", deparse(repository)),
        sprintf("paths <- readRDS(%s)", deparse(input)),
        "workflow <- new.env(parent = emptyenv())",
        "workflow$state_manager <- WorkflowState$new()",
        "workflow$artifact_stage_results <- NULL",
        paste0(
            "workflow$tab_status <- list(setup_import = 'pending', ",
            "design_matrix = 'disabled', quality_control = 'disabled', ",
            "normalization = 'disabled', differential_expression = 'disabled', ",
            "enrichment_analysis = 'disabled', session_summary = 'disabled')"
        ),
        "result <- resumeProtDiaArtifactWorkflow(workflow, paths, 'dia-study')",
        paste0(
            "preview <- previewProtDiaImportArtifact(workflow$workflow_context, ",
            "result$import_ref, limit = 3L)"
        ),
        paste0(
            "saveRDS(list(result = result, data_tbl = workflow$data_tbl, ",
            "state = workflow$state_manager$getState(), preview = preview, ",
            "identity = workflow$workflow_context$getIdentity(), ",
            "root = workflow$workflow_context$getProjectRoot()), ",
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
    expect_identical(attr(status, "status"), NULL, info = paste(status, collapse = "\n"))
    fresh <- readRDS(output)
    expect_true(fresh$result$resumed)
    expect_identical(fresh$identity$project_id, built$paths$project_id)
    expect_identical(fresh$root, normalizePath(moved, winslash = "/"))
    expect_identical(fresh$data_tbl, oracle$data_tbl)
    expectDiaArtifact010StateExact(oracle$state, fresh$state)
    expect_identical(nrow(fresh$preview), 3L)
})

test_that("legacy fallback is explicit and memory roots remain registry-free", {
    diaArtifact010SkipDependencies()
    root <- withr::local_tempdir(pattern = "dia-artifact-010-legacy-")
    paths <- diaArtifact010Paths(root)
    workflow <- diaArtifact010BlankWorkflow()
    disabled <- resumeProtDiaArtifactWorkflow(
        workflow,
        paths,
        "dia-study",
        storage_policy = list(requested_backend = "memory")
    )
    expect_false(disabled$enabled)
    expect_false(disabled$resumed)
    expect_false(dir.exists(file.path(root, "state")))

    artifact_root <- withr::local_tempdir(pattern = "dia-artifact-010-unlinked-")
    built <- diaArtifact010Build(artifact_root, link_import = FALSE)
    candidate <- diaArtifact010BlankWorkflow()
    calls <- 0L
    result <- resumeProtDiaArtifactWorkflow(
        candidate,
        built$paths,
        "dia-study",
        legacy_adapter = function(...) {
            calls <<- calls + 1L
            list(source = "legacy-session")
        }
    )
    expect_identical(calls, 1L)
    expect_identical(result$reason, "explicit_legacy_adapter")
    expect_false(result$resumed)
    expect_identical(candidate$state_manager$getCurrentStateName(), "initial")
    expect_false(exists("data_tbl", envir = candidate, inherits = FALSE))
})
