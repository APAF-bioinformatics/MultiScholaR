nondiaArtifact055SkipDependencies <- function() {
    for (package in c("arrow", "DBI", "duckdb", "filelock")) {
        testthat::skip_if_not_installed(package)
    }
}

nondiaArtifact055RepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

nondiaArtifact055Scenarios <- function() {
    manifest <- jsonlite::read_json(
        nondiaArtifact055RepoPath(
            "tests", "testdata", "omics-parity", "proteomics", "manifest.json"
        ),
        simplifyVector = FALSE
    )
    manifest$fixture_scenarios
}

nondiaArtifact055Importer <- function(format) {
    switch(
        format,
        maxquant = importMaxQuantData,
        fragpipe = importFragPipeData,
        pd_tmt = importProteomeDiscovererTMTData
    )
}

nondiaArtifact055WorkflowType <- function(format) {
    if (identical(format, "pd_tmt")) "TMT" else "LFQ"
}

nondiaArtifact055CapabilityId <- function(format) {
    paste0(
        "proteomics.",
        switch(format, maxquant = "maxquant", fragpipe = "fragpipe", pd_tmt = "pd_tmt"),
        ".protein.",
        if (identical(format, "pd_tmt")) "tmt" else "lfq",
        ".v1"
    )
}

nondiaArtifact055Paths <- function(root, project_id) {
    paths <- list(
        base_dir = root,
        project_id = project_id,
        omic_type = "proteomics",
        omic_label = "session-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    paths
}

nondiaArtifact055BlankWorkflow <- function() {
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- NULL
    workflow$state_manager <- WorkflowState$new()
    workflow$artifact_stage_results <- list()
    workflow$artifact_readthrough_refs <- NULL
    workflow$artifact_readthrough_proof <- NULL
    workflow$artifact_compatibility_checkpoint <- NULL
    workflow$data_tbl <- NULL
    workflow$data_cln <- NULL
    workflow$design_matrix <- NULL
    workflow$contrasts_tbl <- NULL
    workflow$config_list <- NULL
    workflow$column_mapping <- NULL
    workflow$uniprot_dat_cln <- NULL
    workflow$aa_seq_tbl_final <- NULL
    workflow$normalization_state_refs <- list()
    workflow$final_for_da_ref <- NULL
    workflow$ruv_optimization_result <- NULL
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

nondiaArtifact055Build <- function(root, scenario) {
    format <- scenario$format
    capability_id <- nondiaArtifact055CapabilityId(format)
    descriptor <- protNonDiaReadthroughDescriptor(capability_id)
    project_id <- paste0("nondia-055-", format)
    paths <- nondiaArtifact055Paths(root, project_id)
    fixture <- nondiaArtifact055RepoPath(scenario$fixture_path)
    source <- file.path(paths$source_dir, basename(fixture))
    stopifnot(file.copy(fixture, source))
    workflow <- nondiaArtifact055BlankWorkflow()
    workflow$workflow_context <- createWorkflowContext(
        paths,
        "proteomics",
        "session-study",
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = project_id
        )
    )
    imported <- suppressMessages(nondiaArtifact055Importer(format)(source))
    workflow$data_tbl <- imported$data
    workflow$data_cln <- imported$data
    workflow$data_format <- format
    workflow$data_type <- imported$data_type
    workflow$column_mapping <- imported$column_mapping
    workflow_type <- nondiaArtifact055WorkflowType(format)
    report_template <- if (identical(workflow_type, "TMT")) {
        "TMT_report.rmd"
    } else {
        "LFQ_report.rmd"
    }
    workflow$state_manager$setWorkflowType(workflow_type)
    import_result <- persistProtNonDiaImportArtifacts(workflow, imported, source)
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
    workflow$config_list <- list(globalParameters = list(
        workflow_type = workflow_type,
        report_template = report_template
    ))
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
        "non-DIA session fixture",
        validateColumnMapping = TRUE
    ))
    design_result <- persistProtNonDiaDesignArtifacts(workflow)
    stopifnot(isTRUE(design_result$ok))

    resumed <- nondiaArtifact055BlankWorkflow()
    resume_result <- resumeProtNonDiaArtifactWorkflow(
        resumed,
        paths,
        "session-study"
    )
    stopifnot(isTRUE(resume_result$resumed))
    eviction <- evictProtNonDiaWorkflowPayloads(
        resumed,
        descriptor,
        rollout_fn = \(...) "evict"
    )
    stopifnot(isTRUE(eviction$evicted))
    before <- resumed$state_manager$getState()
    normalized <- module_ci_prot_norm_s4_methods()$normalise(
        before,
        normalisation_method = "none"
    )
    norm_data <- new.env(parent = emptyenv())
    norm_data$state_manager <- resumed$state_manager
    norm_data$state_refs <- list()
    norm_data$normalized_protein_obj <- NULL
    norm_data$ruv_normalized_obj <- NULL
    norm_data$correlation_filtered_obj <- NULL
    normalized <- saveProtNormState(
        resumed,
        resumed$state_manager,
        before,
        normalized,
        "normalization",
        "normalized",
        list(normalization_method = "none", ruv_mode = "skip"),
        "Applied session normalization",
        list(normalization_method = "none", ruv_mode = "skip"),
        transformation_type = "normalization"
    )
    settleProtNormArtifactState(
        resumed,
        norm_data,
        "normalization",
        "normalized",
        normalized
    )
    resumed$normalization_context <- list(
        ruv_mode = "skip",
        ruv_grouping_variable = "group"
    )
    resumed$ruv_optimization_result <- buildProtNormSkippedRuvResult()
    norm_data$correlation_filtered_obj <- normalized
    finalizeProtNormCorrelationWorkflowState(
        normalized,
        resumed,
        normData = norm_data,
        correlationThreshold = 0,
        skipped = TRUE,
        timeFn = \() as.POSIXct("2026-08-24 00:00:00", tz = "UTC"),
        messageFn = \(...) invisible(NULL),
        catFn = \(...) invisible(NULL)
    )
    final <- resumed$state_manager$getState("correlation_filtered")
    resumed$config_list <- final@args
    resumed$tab_status$normalization <- "complete"
    resumed$tab_status$differential_expression <- "pending"
    plot_path <- file.path(paths$results_dir, "post_normalization_qc.png")
    writeLines("portable plot fixture", plot_path)
    norm_data$qc_plot_paths <- c(post_normalization = plot_path)
    session_data <- collectProtNormExportSessionData(
        resumed,
        list2env(list(
            correlation_filtered_obj = final,
            best_k = NA_integer_,
            correlation_threshold = 0,
            qc_plot_paths = norm_data$qc_plot_paths
        ), parent = emptyenv()),
        list(norm_method = "none", ruv_mode = "skip"),
        timeFn = \() as.POSIXct("2026-08-24 00:00:00", tz = "UTC"),
        messageFn = \(...) invisible(NULL)
    )
    session_data$normalization_state_refs <- resumed$normalization_state_refs
    session_data$final_for_da_ref <- resumed$final_for_da_ref
    export_artifacts <- saveProtNormExportArtifacts(
        session_data,
        paths$source_dir,
        timeFn = \() as.POSIXct("2026-08-24 00:00:00", tz = "UTC"),
        messageFn = \(...) invisible(NULL)
    )
    artifact_session <- saveProtArtifactSessionManifest(
        resumed,
        norm_data,
        session_data,
        export_artifacts,
        paths$source_dir
    )
    list(
        paths = paths,
        source = source,
        workflow = resumed,
        norm_data = norm_data,
        session_data = session_data,
        export_artifacts = export_artifacts,
        artifact_session = artifact_session,
        final = final,
        descriptor = descriptor,
        capability_id = capability_id,
        workflow_type = workflow_type,
        report_template = report_template,
        scenario = scenario
    )
}

expectNondiaArtifact055StateExact <- function(expected, actual) {
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

nondiaArtifact055Move <- function(built) {
    moved <- tempfile(
        pattern = paste0("nondia-055-moved-", built$scenario$format, "-")
    )
    dir.create(moved, recursive = TRUE)
    entries <- list.files(
        built$paths$base_dir,
        all.files = TRUE,
        no.. = TRUE,
        full.names = TRUE
    )
    expect_true(all(file.copy(entries, moved, recursive = TRUE)))
    built$workflow$state_manager$close()
    unlink(built$paths$base_dir, recursive = TRUE, force = TRUE)
    list(
        base_dir = moved,
        omic_type = "proteomics",
        omic_label = "session-study",
        source_dir = file.path(moved, "source"),
        results_dir = file.path(moved, "results")
    )
}

test_that("non-DIA manifests are bounded, relative, and generation pinned", {
    nondiaArtifact055SkipDependencies()
    for (scenario in nondiaArtifact055Scenarios()) {
        root <- withr::local_tempdir(
            pattern = paste0("nondia-055-export-", scenario$format, "-")
        )
        built <- nondiaArtifact055Build(root, scenario)
        manifest <- built$artifact_session$manifest
        expect_true(built$artifact_session$ok)
        expect_identical(manifest$descriptor$descriptor_id, built$capability_id)
        expect_identical(
            manifest$workflow_state$current_state,
            "correlation_filtered"
        )
        expect_identical(
            manifest$workflow_state$current_generation_id,
            built$workflow$state_manager$getCurrentGenerationId()
        )
        expect_setequal(names(manifest$stage_refs), c("import", "design"))
        expect_true(all(vapply(
            manifest$workflow_state$active_lineage,
            \(entry) !grepl("^/", entry$manifest_relative_path),
            logical(1)
        )))
        expect_false(grepl(
            normalizePath(root, winslash = "/"),
            workflowSessionJson(manifest),
            fixed = TRUE
        ))
        metadata <- workflowSessionDecodeMetadata(manifest$metadata_json)
        expect_null(metadata$r6_complete_states)
        expect_null(metadata$state_manager)
        expect_null(metadata$workflow_context)
        expect_lt(nchar(manifest$metadata_json, type = "bytes"), 1024^2)
        expect_identical(
            manifest$compatibility$byte_digest,
            artifactByteDigest(built$export_artifacts$sessionFilepath)
        )
        built$workflow$state_manager$close()
    }
})

test_that("moved artifact and legacy sessions restore equivalent current inputs", {
    nondiaArtifact055SkipDependencies()
    for (scenario in nondiaArtifact055Scenarios()) {
        root <- withr::local_tempdir(
            pattern = paste0("nondia-055-restore-", scenario$format, "-")
        )
        built <- nondiaArtifact055Build(root, scenario)
        oracle <- built$final
        moved_paths <- nondiaArtifact055Move(built)
        manifest_path <- file.path(
            moved_paths$source_dir,
            "filtered_session_artifact_latest.json"
        )
        bundle <- restoreProtArtifactSessionManifest(
            manifest_path,
            moved_paths
        )
        withr::defer(bundle$state_manager$close())
        expectNondiaArtifact055StateExact(
            oracle,
            bundle$session_data$current_s4_object
        )
        expect_identical(bundle$session_data$design_matrix, oracle@design_matrix)
        expect_identical(
            bundle$session_data$contrasts_tbl,
            built$session_data$contrasts_tbl
        )
        expect_identical(bundle$session_data$config_list, oracle@args)
        expect_identical(bundle$uniprot_dat_cln, built$workflow$uniprot_dat_cln)
        expect_identical(
            resolveProtSummaryExpectedTemplate(bundle$session_data$workflow_type),
            built$report_template
        )
        da_data <- new.env(parent = emptyenv())
        da_data$current_s4_object <- NULL
        da_data$contrasts_available <- NULL
        da_data$formula_from_s4 <- NULL
        workflow <- nondiaArtifact055BlankWorkflow()
        applied <- applyProtDaSessionBundle(workflow, da_data, bundle)
        expectNondiaArtifact055StateExact(oracle, da_data$current_s4_object)
        expect_setequal(
            da_data$contrasts_available,
            built$session_data$contrasts_tbl$contrasts
        )
        expect_identical(applied$state_snapshot$r6_current_state_name,
            "correlation_filtered")

        latest_rds <- file.path(
            moved_paths$source_dir,
            "filtered_session_data_latest.rds"
        )
        unlink(latest_rds, force = TRUE)
        compatibility <- ensureProtDaCompatibilitySession(
            bundle,
            moved_paths$source_dir
        )
        expect_true(file.exists(compatibility$path))
        reconstructed <- readRDS(compatibility$path)
        expectNondiaArtifact055StateExact(
            oracle,
            reconstructed$current_s4_object
        )
        expect_identical(
            reconstructed$artifact_source_generation,
            bundle$manifest$workflow_state$current_generation_id
        )

        withr::with_options(list(
            multischolar.prot_nondia.artifact_session_restore = "disabled"
        ), {
            legacy <- resolveProtDaSessionBundle(moved_paths)
            expect_identical(legacy$format, "legacy")
            expectNondiaArtifact055StateExact(
                oracle,
                legacy$session_data$current_s4_object
            )
        })
    }
})

test_that("non-DIA session publication, restore, and apply failures are atomic", {
    nondiaArtifact055SkipDependencies()
    scenario <- nondiaArtifact055Scenarios()[[1L]]
    root <- withr::local_tempdir(pattern = "nondia-055-failure-")
    built <- nondiaArtifact055Build(root, scenario)
    latest <- built$artifact_session$latestFilepath
    digest_before <- artifactByteDigest(latest)
    expect_error(
        saveProtArtifactSessionManifest(
            built$workflow,
            built$norm_data,
            built$session_data,
            list(
                sessionFilename = "filtered_session_data_second.rds",
                sessionFilepath = built$export_artifacts$sessionFilepath
            ),
            built$paths$source_dir,
            failure_injector = \(stage, ...) {
                if (identical(stage, "before_latest_publication")) {
                    stop("injected latest failure")
                }
            }
        ),
        "injected latest failure",
        fixed = TRUE
    )
    expect_identical(artifactByteDigest(latest), digest_before)

    moved_paths <- nondiaArtifact055Move(built)
    manifest_path <- file.path(
        moved_paths$source_dir,
        "filtered_session_artifact_latest.json"
    )
    manifest_before <- readWorkflowSessionManifest(manifest_path)
    plot_path <- artifactResolveContainedPath(
        moved_paths$base_dir,
        manifest_before$plot_refs$post_normalization$relative_path,
        must_exist = TRUE
    )
    plot_contents <- readLines(plot_path, warn = FALSE)
    writeLines("corrupt plot", plot_path)
    expect_error(
        restoreProtArtifactSessionManifest(manifest_path, moved_paths),
        class = "multischolar_prot_nondia_session_plot_mismatch"
    )
    writeLines(plot_contents, plot_path)
    bundle <- restoreProtArtifactSessionManifest(manifest_path, moved_paths)
    workflow <- nondiaArtifact055BlankWorkflow()
    sentinel <- module_ci_prot_norm_object()
    workflow$state_manager$saveState("sentinel", sentinel, list(), "sentinel")
    workflow$contrasts_tbl <- data.frame(contrasts = "sentinel")
    da_data <- new.env(parent = emptyenv())
    da_data$current_s4_object <- sentinel
    da_data$contrasts_available <- "sentinel"
    da_data$formula_from_s4 <- "~ sentinel"
    expect_error(
        applyProtDaSessionBundle(
            workflow,
            da_data,
            bundle,
            failure_injector = \(stage, ...) {
                if (identical(stage, "before_session_apply_commit")) {
                    stop("injected apply failure")
                }
            }
        ),
        "injected apply failure",
        fixed = TRUE
    )
    expect_identical(workflow$contrasts_tbl$contrasts, "sentinel")
    expect_identical(da_data$current_s4_object, sentinel)

    legacy_path <- file.path(
        moved_paths$source_dir,
        "filtered_session_data_latest.rds"
    )
    expect_true(file.exists(legacy_path))
    manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
    manifest$descriptor$descriptor_id <-
        "proteomics.spectronaut.protein.lfq.v1"
    manifest$manifest_digest <- workflowSessionContentDigest(manifest)
    jsonlite::write_json(
        manifest,
        manifest_path,
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null"
    )
    fallback <- resolveProtDaSessionBundle(moved_paths)
    expect_identical(fallback$format, "legacy")
    expect_true(inherits(fallback$artifact_error, "error"))
})
