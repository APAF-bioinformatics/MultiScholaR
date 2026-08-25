metab053SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

metab053Context <- function(root, project_id) {
    paths <- list(
        base_dir = root,
        project_id = project_id,
        omic_type = "metabolomics",
        omic_label = "metabolomics-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    capability <- workflowCapabilityCatalogue()[[
        "metabolomics.custom.metabolite.standard.v1"
    ]]
    capability$artifact_eligible <- TRUE
    capability$auto_eligible <- FALSE
    capability$maximum_artifact_rollout <- "dual_write"
    context <- createWorkflowContext(
        paths,
        "metabolomics",
        "metabolomics-study",
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = project_id
        )
    )
    bindWorkflowContextFromImport(
        context,
        workflow_type = "metabolomics_standard",
        input_format = "custom",
        data_level = "metabolite",
        capabilities = list(capability)
    )
    context
}

metab053Manager <- function(context) {
    manager <- ArtifactWorkflowState$new(
        workflow_context = context,
        dehydrate_fn = dehydrateMetabolomicsS4Artifact,
        validate_bundle_fn = validateMetabolomicsS4Bundle,
        hydrate_fn = hydrateMetabolomicsS4Artifact,
        descriptor_contract = artifactStageDescriptorContract(
            artifactMetabolomicsWorkflowDescriptor()
        )
    )
    manager$setWorkflowType("metabolomics_standard")
    manager
}

metab053Workflow <- function(manager, context, object) {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$workflow_context <- context
    workflow$config_list <- object@args
    workflow$contrasts_tbl <- data.frame(
        friendly_names = "KO_vs_WT",
        contrasts = "groupKO-groupWT",
        stringsAsFactors = FALSE
    )
    workflow$design_matrix <- object@design_matrix
    workflow$metabolite_counts <- module_ci_metab_norm_counts(object)
    workflow$artifact_stage_results <- list()
    workflow$data_format <- "custom"
    workflow$data_type <- "metabolite"
    workflow
}

metab053SessionData <- function(workflow, object, timestamp) {
    list(
        r6_current_state_name = "metab_norm_complete",
        current_s4_object = object,
        contrasts_tbl = workflow$contrasts_tbl,
        design_matrix = object@design_matrix,
        config_list = workflow$config_list,
        itsd_selections = list(),
        ruv_optimization_results = list(),
        correlation_results = list(),
        assay_names = names(object@metabolite_data),
        export_timestamp = timestamp,
        omic_type = "metabolomics",
        experiment_label = "portable-session",
        normalization_method = "quantile",
        ruv_mode = "skip",
        itsd_applied = FALSE,
        itsd_aggregation = NA_character_,
        log_offset = 1,
        correlation_threshold = 0.75,
        ruv_grouping_variable = "tech_rep_group",
        feature_counts = module_ci_metab_norm_counts(object),
        metabolite_counts = module_ci_metab_norm_counts(object),
        normalization_complete = TRUE,
        ruv_complete = FALSE,
        correlation_filtering_complete = TRUE
    )
}

metab053SaveProject <- function(root, layout) {
    project_id <- paste0("omics-art-053-", layout)
    context <- metab053Context(root, project_id)
    manager <- metab053Manager(context)
    object <- module_ci_metab_norm_object(
        layout = layout,
        positive_only = TRUE
    )
    object@metabolite_data <- lapply(object@metabolite_data, function(assay) {
        rownames(assay) <- NULL
        assay
    })
    manager$saveState(
        "metab_norm_complete",
        object,
        object@args,
        "portable metabolomics session"
    )
    workflow <- metab053Workflow(manager, context, object)
    session_data <- metab053SessionData(
        workflow,
        object,
        as.POSIXct("2026-08-25 12:00:00", tz = "UTC")
    )
    source_dir <- file.path(root, "source")
    export_files <- saveMetabNormExportSessionRdsFiles(
        session_data,
        source_dir,
        timeFn = function() as.POSIXct("2026-08-25 12:00:00", tz = "UTC"),
        logInfoFn = function(...) invisible(NULL)
    )
    artifact <- saveMetabSessionManifest(
        workflow,
        session_data,
        export_files,
        source_dir
    )
    list(
        context = context,
        manager = manager,
        workflow = workflow,
        object = object,
        session_data = session_data,
        export_files = export_files,
        artifact = artifact,
        project_id = project_id
    )
}

test_that("LC GC and mixed artifact sessions move and restore exactly", {
    metab053SkipDependencies()
    for (layout in c("lc", "gc", "combined")) {
        root <- withr::local_tempdir(pattern = paste0("metab053-", layout, "-"))
        built <- metab053SaveProject(root, layout)
        manifest <- built$artifact$manifest
        encoded <- workflowSessionJson(manifest)
        expect_lt(nchar(encoded, type = "bytes"), 64 * 1024)
        expect_false(grepl(normalizePath(root), encoded, fixed = TRUE), info = layout)
        expect_false(grepl("connection|cache", encoded, ignore.case = TRUE), info = layout)
        expect_true(all(vapply(
            manifest$workflow_state$active_lineage,
            function(entry) !grepl("^/", entry$manifest_relative_path),
            logical(1)
        )), info = layout)
        built$manager$close()

        moved <- withr::local_tempdir(pattern = paste0("metab053-moved-", layout, "-"))
        fs::dir_copy(root, moved, overwrite = TRUE)
        moved_context <- metab053Context(moved, built$project_id)
        restored <- restoreMetabSessionManifest(
            file.path(
                moved,
                "source",
                "metab_filtered_session_artifact_latest.json"
            ),
            moved_context
        )
        withr::defer(restored$manager$close())

        expect_identical(restored$session_data$current_s4_object, built$object)
        expect_identical(
            restored$session_data$current_s4_object@design_matrix,
            built$session_data$design_matrix
        )
        expect_identical(
            restored$session_data$contrasts_tbl,
            built$session_data$contrasts_tbl
        )
        expect_identical(
            restored$session_data$assay_names,
            names(built$object@metabolite_data)
        )
        expected_plots <- suppressMessages(plotPca(
            built$object,
            grouping_variable = "group"
        ))
        restored_plots <- suppressMessages(plotPca(
            restored$session_data$current_s4_object,
            grouping_variable = "group"
        ))
        expect_identical(
            lapply(expected_plots, function(plot) ggplot2::ggplot_build(plot)$data),
            lapply(restored_plots, function(plot) ggplot2::ggplot_build(plot)$data)
        )
        expect_true(restored$compatibility$available)

        reconstructed <- file.path(moved, "source", "reconstructed.rds")
        writeMetabCompatibilitySession(restored, reconstructed)
        legacy <- readRDS(reconstructed)
        expect_identical(legacy$current_s4_object, built$object)
        expect_identical(legacy$design_matrix, built$object@design_matrix)
    }
})

test_that("artifact and legacy readers produce equivalent DA session inputs", {
    metab053SkipDependencies()
    root <- withr::local_tempdir()
    built <- metab053SaveProject(root, "combined")
    artifact_data <- readMetabArtifactOrLegacySession(
        built$export_files$latestFilepath,
        built$workflow
    )
    withr::local_options(
        multischolar.metabolomics.artifact_session_restore = "disabled"
    )
    legacy_data <- readMetabArtifactOrLegacySession(
        built$export_files$latestFilepath,
        built$workflow
    )

    expect_identical(artifact_data$current_s4_object, legacy_data$current_s4_object)
    expect_identical(artifact_data$design_matrix, legacy_data$design_matrix)
    expect_identical(artifact_data$contrasts_tbl, legacy_data$contrasts_tbl)
    expect_identical(artifact_data$assay_names, legacy_data$assay_names)
    expect_identical(
        artifact_data$current_s4_object@args,
        legacy_data$current_s4_object@args
    )
    built$manager$close()
})

test_that("missing or corrupt artifact sessions preserve valid legacy reload", {
    metab053SkipDependencies()
    root <- withr::local_tempdir()
    built <- metab053SaveProject(root, "gc")
    manifest_path <- built$artifact$latestFilepath
    file.remove(manifest_path)
    expect_identical(
        readMetabArtifactOrLegacySession(
            built$export_files$latestFilepath,
            built$workflow
        ),
        readRDS(built$export_files$latestFilepath)
    )
    writeLines("{not-json", manifest_path)
    before <- readRDS(built$export_files$latestFilepath)
    expect_identical(
        readMetabArtifactOrLegacySession(
            built$export_files$latestFilepath,
            built$workflow
        ),
        before
    )
    expect_identical(readRDS(built$export_files$latestFilepath), before)
    built$manager$close()
})

test_that("failed manifest publication leaves prior artifact and RDS usable", {
    metab053SkipDependencies()
    root <- withr::local_tempdir()
    built <- metab053SaveProject(root, "lc")
    previous <- readWorkflowSessionManifest(built$artifact$latestFilepath)
    result <- saveMetabSessionManifestSafely(
        workflow_data = built$workflow,
        session_data = built$session_data,
        export_files = built$export_files,
        source_dir = file.path(root, "source"),
        failure_injector = function(stage, context) {
            if (identical(stage, "before_latest_publication")) {
                stop("injected manifest failure")
            }
            invisible(context)
        }
    )

    expect_false(result$ok)
    expect_identical(
        readWorkflowSessionManifest(built$artifact$latestFilepath),
        previous
    )
    expect_identical(
        readRDS(built$export_files$latestFilepath)$current_s4_object,
        built$object
    )
    built$manager$close()
})

test_that("manifest size is independent of assay payload volume", {
    metab053SkipDependencies()
    small_root <- withr::local_tempdir(pattern = "metab053-small-")
    large_root <- withr::local_tempdir(pattern = "metab053-large-")
    small <- metab053SaveProject(small_root, "combined")
    large_object <- module_ci_metab_norm_object(
        layout = "combined",
        positive_only = TRUE
    )
    large_object@metabolite_data <- lapply(
        large_object@metabolite_data,
        function(assay) {
            value <- assay[rep(seq_len(nrow(assay)), 200L), , drop = FALSE]
            value$database_identifier <- paste0(
                value$database_identifier,
                "_",
                seq_len(nrow(value))
            )
            value$annotation_id <- value$database_identifier
            value$Name <- value$database_identifier
            value$S1 <- value$S1 + seq_len(nrow(value)) / 1000
            value
        }
    )
    rownames(large_object@metabolite_data[[1L]]) <- NULL
    rownames(large_object@metabolite_data[[2L]]) <- NULL
    context <- metab053Context(large_root, "omics-art-053-large")
    manager <- metab053Manager(context)
    manager$saveState(
        "metab_norm_complete",
        large_object,
        large_object@args,
        "large portable session"
    )
    workflow <- metab053Workflow(manager, context, large_object)
    session_data <- metab053SessionData(
        workflow,
        large_object,
        as.POSIXct("2026-08-25 12:00:00", tz = "UTC")
    )
    files <- saveMetabNormExportSessionRdsFiles(
        session_data,
        file.path(large_root, "source"),
        timeFn = function() as.POSIXct("2026-08-25 12:00:00", tz = "UTC"),
        logInfoFn = function(...) invisible(NULL)
    )
    large <- saveMetabSessionManifest(
        workflow,
        session_data,
        files,
        file.path(large_root, "source")
    )
    small_bytes <- file.info(small$artifact$latestFilepath)$size
    large_bytes <- file.info(large$latestFilepath)$size

    expect_lt(abs(large_bytes - small_bytes), 2048)
    expect_gt(file.info(files$latestFilepath)$size, large_bytes * 10)
    small$manager$close()
    manager$close()
})

test_that("frozen mixed CI workload digests survive artifact session restore", {
    metab053SkipDependencies()
    repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    adapter_env <- new.env(parent = globalenv())
    sys.source(
        file.path(repo, "tools", "profiling", "omics_workload_metabolomics.R"),
        envir = adapter_env
    )
    contract <- jsonlite::read_json(
        file.path(
            repo, "tests", "testdata", "omics-parity", "metabolomics",
            "workloads", "mixed-public-ci-v1.json"
        ),
        simplifyVector = FALSE
    )
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    root <- withr::local_tempdir()
    prepared <- adapter_env$metabWorkloadPrepare(list(
        contract = contract,
        run_dir = root
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
    generated <- adapter_env$metabWorkloadRead(prepared$payload_path)
    truth <- jsonlite::read_json(prepared$truth_path, simplifyVector = TRUE)
    samples <- unlist(truth$sample_ids, use.names = FALSE)
    assay_names <- unlist(truth$assays, use.names = FALSE)
    assays <- lapply(assay_names, function(assay_name) {
        value <- generated[generated$assay == assay_name, , drop = FALSE]
        value$assay <- NULL
        rownames(value) <- NULL
        value
    })
    names(assays) <- assay_names
    object <- createMetaboliteAssayData(
        assays,
        data.frame(
            Run = samples,
            group = unlist(truth$group_assignments, use.names = FALSE),
            batch = unlist(truth$batch_assignments, use.names = FALSE),
            tech_rep_group = samples,
            stringsAsFactors = FALSE
        ),
        metabolite_id_column = "feature_id",
        annotation_id_column = "annotation",
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        args = list(evidence_class = "generated_scaling")
    )
    context <- metab053Context(root, "omics-art-053-frozen-ci")
    manager <- metab053Manager(context)
    manager$saveState(
        "metab_norm_complete",
        object,
        object@args,
        "frozen CI session"
    )
    workflow <- metab053Workflow(manager, context, object)
    session_data <- metab053SessionData(
        workflow,
        object,
        as.POSIXct("2026-08-25 12:00:00", tz = "UTC")
    )
    files <- saveMetabNormExportSessionRdsFiles(
        session_data,
        file.path(root, "source"),
        timeFn = function() as.POSIXct("2026-08-25 12:00:00", tz = "UTC"),
        logInfoFn = function(...) invisible(NULL)
    )
    artifact <- saveMetabSessionManifest(
        workflow,
        session_data,
        files,
        file.path(root, "source")
    )
    manager$close()
    restored <- restoreMetabSessionManifest(artifact$latestFilepath, context)

    expect_identical(restored$session_data$current_s4_object, object)
    expect_identical(restored$session_data$assay_names, assay_names)
    expect_lt(file.info(artifact$latestFilepath)$size, 64 * 1024)
    restored$manager$close()
})
