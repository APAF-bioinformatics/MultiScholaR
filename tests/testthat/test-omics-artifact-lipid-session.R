lipid054SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

lipid054Context <- function(root, project_id) {
    paths <- list(
        base_dir = root,
        project_id = project_id,
        omic_type = "lipidomics",
        omic_label = "lipidomics-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    capability <- workflowCapabilityCatalogue()[[
        "lipidomics.lipidsearch.lipid.standard.v1"
    ]]
    capability$artifact_eligible <- TRUE
    capability$auto_eligible <- FALSE
    capability$maximum_artifact_rollout <- "dual_write"
    context <- createWorkflowContext(
        paths,
        "lipidomics",
        "lipidomics-study",
        storage_policy = list(
            requested_backend = "artifact",
            requested_rollout = "dual_write",
            migration_requested = TRUE,
            project_id = project_id
        )
    )
    bindWorkflowContextFromImport(
        context,
        workflow_type = "lipidomics_standard",
        input_format = "lipidsearch",
        data_level = "lipid",
        capabilities = list(capability)
    )
    context
}

lipid054Manager <- function(context) {
    manager <- ArtifactWorkflowState$new(
        workflow_context = context,
        dehydrate_fn = dehydrateLipidomicsS4Artifact,
        validate_bundle_fn = validateLipidomicsS4Bundle,
        hydrate_fn = hydrateLipidomicsS4Artifact,
        descriptor_contract = artifactStageDescriptorContract(
            artifactLipidomicsWorkflowDescriptor()
        )
    )
    manager$setWorkflowType("lipidomics_standard")
    manager
}

lipid054Workflow <- function(manager, context, object) {
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
    workflow$lipid_counts <- module_ci_lipid_norm_counts(object)
    workflow$artifact_stage_results <- list()
    workflow$data_format <- "lipidsearch"
    workflow$data_type <- "lipid"
    workflow
}

lipid054SessionData <- function(workflow, object, timestamp) {
    list(
        r6_current_state_name = "lipid_norm_complete",
        current_s4_object = object,
        contrasts_tbl = workflow$contrasts_tbl,
        design_matrix = object@design_matrix,
        config_list = workflow$config_list,
        itsd_selections = list(),
        ruv_optimization_results = list(),
        correlation_results = list(),
        assay_names = names(object@lipid_data),
        export_timestamp = timestamp,
        omic_type = "lipidomics",
        experiment_label = "portable-session",
        normalization_method = "quantile",
        ruv_mode = "skip",
        itsd_applied = FALSE,
        itsd_aggregation = NA_character_,
        log_offset = 1,
        correlation_threshold = 0.75,
        ruv_grouping_variable = "tech_rep_group",
        feature_counts = module_ci_lipid_norm_counts(object),
        lipid_counts = module_ci_lipid_norm_counts(object),
        normalization_complete = TRUE,
        ruv_complete = FALSE,
        correlation_filtering_complete = TRUE
    )
}

lipid054SaveRdsFiles <- function(
    session_data,
    source_dir,
    timeFn = Sys.time,
    logInfoFn = logger::log_info
) {
    dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
    timestamp <- format(timeFn(), "%Y%m%d_%H%M%S")
    session_filename <- sprintf(
        "lipid_filtered_session_data_%s.rds",
        timestamp
    )
    session_filepath <- file.path(source_dir, session_filename)
    latest_filepath <- file.path(
        source_dir,
        "lipid_filtered_session_data_latest.rds"
    )
    saveRDS(session_data, session_filepath)
    saveRDS(session_data, latest_filepath)
    logInfoFn("saved lipid session fixtures")
    list(
        session_filename = session_filename,
        session_filepath = session_filepath,
        latest_filepath = latest_filepath
    )
}

lipid054SaveProject <- function(root, layout) {
    project_id <- paste0("omics-art-054-", layout)
    context <- lipid054Context(root, project_id)
    manager <- lipid054Manager(context)
    object <- module_ci_lipid_norm_object(
        layout = layout,
        positive_only = TRUE
    )
    object@lipid_data <- lapply(object@lipid_data, function(assay) {
        rownames(assay) <- NULL
        assay
    })
    manager$saveState(
        "lipid_norm_complete",
        object,
        object@args,
        "portable lipidomics session"
    )
    workflow <- lipid054Workflow(manager, context, object)
    session_data <- lipid054SessionData(
        workflow,
        object,
        as.POSIXct("2026-08-25 12:00:00", tz = "UTC")
    )
    source_dir <- file.path(root, "source")
    export_files <- lipid054SaveRdsFiles(
        session_data,
        source_dir,
        timeFn = function() as.POSIXct("2026-08-25 12:00:00", tz = "UTC"),
        logInfoFn = function(...) invisible(NULL)
    )
    artifact <- saveLipidSessionManifest(
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
    lipid054SkipDependencies()
    for (layout in c("lc", "gc", "combined")) {
        root <- withr::local_tempdir(pattern = paste0("lipid054-", layout, "-"))
        built <- lipid054SaveProject(root, layout)
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

        moved <- withr::local_tempdir(pattern = paste0("lipid054-moved-", layout, "-"))
        fs::dir_copy(root, moved, overwrite = TRUE)
        moved_context <- lipid054Context(moved, built$project_id)
        restored <- restoreLipidSessionManifest(
            file.path(
                moved,
                "source",
                "lipid_filtered_session_artifact_latest.json"
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
            names(built$object@lipid_data)
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
        writeLipidCompatibilitySession(restored, reconstructed)
        legacy <- readRDS(reconstructed)
        expect_identical(legacy$current_s4_object, built$object)
        expect_identical(legacy$design_matrix, built$object@design_matrix)
    }
})

test_that("artifact and legacy readers produce equivalent DA session inputs", {
    lipid054SkipDependencies()
    root <- withr::local_tempdir()
    built <- lipid054SaveProject(root, "combined")
    artifact_data <- readLipidArtifactOrLegacySession(
        built$export_files$latest_filepath,
        built$workflow
    )
    withr::local_options(
        multischolar.lipidomics.artifact_session_restore = "disabled"
    )
    legacy_data <- readLipidArtifactOrLegacySession(
        built$export_files$latest_filepath,
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

test_that("production export and DA reader support artifact-only reload", {
    lipid054SkipDependencies()
    root <- withr::local_tempdir()
    context <- lipid054Context(root, "omics-art-054-production")
    manager <- lipid054Manager(context)
    withr::defer(manager$close())
    object <- module_ci_lipid_norm_object(
        layout = "combined",
        positive_only = TRUE
    )
    manager$saveState(
        "lipid_norm_complete",
        object,
        object@args,
        "production session parent"
    )
    workflow <- lipid054Workflow(manager, context, object)
    workflow$qc_params <- list(correlation_threshold = 0.75)
    norm_data <- module_ci_lipid_norm_data(object)
    norm_data$normalization_complete <- TRUE
    norm_data$ruv_complete <- TRUE
    norm_data$correlation_filtering_complete <- TRUE
    source_dir <- file.path(root, "source")
    result <- handleLipidNormExportSession(
        input = module_ci_lipid_norm_input(
            norm_method = "quantile",
            ruv_mode = "skip",
            apply_itsd = FALSE
        ),
        workflowData = workflow,
        experimentPaths = list(source_dir = source_dir),
        experimentLabel = "production",
        normData = norm_data,
        addLog = function(...) invisible(NULL),
        reqFn = function(value) {
            if (is.null(value)) stop("required value missing")
            invisible(value)
        },
        showNotificationFn = function(...) invisible(NULL),
        withProgressFn = function(message, value, expr) force(expr),
        incProgressFn = function(...) invisible(NULL),
        getTimeFn = function() {
            as.POSIXct("2026-08-25 12:00:00", tz = "UTC")
        },
        logInfoFn = function(...) invisible(NULL),
        logWarnFn = function(...) invisible(NULL),
        logErrorFn = function(message) stop(message)
    )

    expect_identical(
        result$session_filename,
        "lipid_filtered_session_data_20260825_120000.rds"
    )
    expect_true(file.exists(result$session_filepath))
    expect_true(file.exists(result$artifact_session_filepath))
    expect_true(file.exists(result$artifact_session_latest_filepath))
    latest_rds <- file.path(source_dir, "lipid_filtered_session_data_latest.rds")
    file.remove(result$session_filepath, latest_rds)
    loaded <- readLipidDaSessionData(
        latest_rds,
        workflowData = workflow,
        logInfo = function(...) invisible(NULL)
    )
    expect_identical(loaded$sessionData$current_s4_object, object)
    expect_identical(loaded$sessionData$contrasts_tbl, workflow$contrasts_tbl)
    expect_identical(loaded$sessionData$design_matrix, object@design_matrix)
})

test_that("missing or corrupt artifact sessions preserve valid legacy reload", {
    lipid054SkipDependencies()
    root <- withr::local_tempdir()
    built <- lipid054SaveProject(root, "gc")
    manifest_path <- built$artifact$latest_filepath
    file.remove(manifest_path)
    expect_identical(
        readLipidArtifactOrLegacySession(
            built$export_files$latest_filepath,
            built$workflow
        ),
        readRDS(built$export_files$latest_filepath)
    )
    writeLines("{not-json", manifest_path)
    before <- readRDS(built$export_files$latest_filepath)
    expect_identical(
        readLipidArtifactOrLegacySession(
            built$export_files$latest_filepath,
            built$workflow
        ),
        before
    )
    expect_identical(readRDS(built$export_files$latest_filepath), before)
    built$manager$close()
})

test_that("failed manifest publication leaves prior artifact and RDS usable", {
    lipid054SkipDependencies()
    root <- withr::local_tempdir()
    built <- lipid054SaveProject(root, "lc")
    previous <- readWorkflowSessionManifest(built$artifact$latest_filepath)
    result <- saveLipidSessionManifestSafely(
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
        readWorkflowSessionManifest(built$artifact$latest_filepath),
        previous
    )
    expect_identical(
        readRDS(built$export_files$latest_filepath)$current_s4_object,
        built$object
    )
    built$manager$close()
})

test_that("manifest size is independent of assay payload volume", {
    lipid054SkipDependencies()
    small_root <- withr::local_tempdir(pattern = "lipid054-small-")
    large_root <- withr::local_tempdir(pattern = "lipid054-large-")
    small <- lipid054SaveProject(small_root, "combined")
    large_object <- module_ci_lipid_norm_object(
        layout = "combined",
        positive_only = TRUE
    )
    large_object@lipid_data <- lapply(
        large_object@lipid_data,
        function(assay) {
            value <- assay[rep(seq_len(nrow(assay)), 200L), , drop = FALSE]
            value$lipid_id <- paste0(
                value$lipid_id,
                "_",
                seq_len(nrow(value))
            )
            value$lipid <- value$lipid_id
            value$LipidClass <- ifelse(
                grepl("^IS_", value$lipid_id),
                "ISTD",
                "PC"
            )
            value$S1 <- value$S1 + seq_len(nrow(value)) / 1000
            value
        }
    )
    rownames(large_object@lipid_data[[1L]]) <- NULL
    rownames(large_object@lipid_data[[2L]]) <- NULL
    context <- lipid054Context(large_root, "omics-art-054-large")
    manager <- lipid054Manager(context)
    manager$saveState(
        "lipid_norm_complete",
        large_object,
        large_object@args,
        "large portable session"
    )
    workflow <- lipid054Workflow(manager, context, large_object)
    session_data <- lipid054SessionData(
        workflow,
        large_object,
        as.POSIXct("2026-08-25 12:00:00", tz = "UTC")
    )
    files <- lipid054SaveRdsFiles(
        session_data,
        file.path(large_root, "source"),
        timeFn = function() as.POSIXct("2026-08-25 12:00:00", tz = "UTC"),
        logInfoFn = function(...) invisible(NULL)
    )
    large <- saveLipidSessionManifest(
        workflow,
        session_data,
        files,
        file.path(large_root, "source")
    )
    small_bytes <- file.info(small$artifact$latest_filepath)$size
    large_bytes <- file.info(large$latest_filepath)$size

    expect_lt(abs(large_bytes - small_bytes), 2048)
    expect_gt(file.info(files$latest_filepath)$size, large_bytes * 10)
    small$manager$close()
    manager$close()
})

test_that("frozen mixed CI workload digests survive artifact session restore", {
    lipid054SkipDependencies()
    repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    adapter_env <- new.env(parent = globalenv())
    sys.source(
        file.path(repo, "tools", "profiling", "omics_workload_lipidomics.R"),
        envir = adapter_env
    )
    contract <- jsonlite::read_json(
        file.path(
            repo, "tests", "testdata", "omics-parity", "lipidomics",
            "workloads", "mixed-public-ci-v1.json"
        ),
        simplifyVector = FALSE
    )
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    root <- withr::local_tempdir()
    prepared <- adapter_env$lipidWorkloadPrepare(list(
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
    generated <- utils::read.delim(
        prepared$payload_path,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
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
    object <- createLipidomicsAssayData(
        assays,
        data.frame(
            Run = samples,
            group = unlist(truth$group_assignments, use.names = FALSE),
            batch = unlist(truth$batch_assignments, use.names = FALSE),
            tech_rep_group = samples,
            stringsAsFactors = FALSE
        ),
        lipid_id_column = "lipid_id",
        annotation_id_column = "lipid_class",
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        args = list(evidence_class = "generated_scaling")
    )
    context <- lipid054Context(root, "omics-art-054-frozen-ci")
    manager <- lipid054Manager(context)
    manager$saveState(
        "lipid_norm_complete",
        object,
        object@args,
        "frozen CI session"
    )
    workflow <- lipid054Workflow(manager, context, object)
    session_data <- lipid054SessionData(
        workflow,
        object,
        as.POSIXct("2026-08-25 12:00:00", tz = "UTC")
    )
    files <- lipid054SaveRdsFiles(
        session_data,
        file.path(root, "source"),
        timeFn = function() as.POSIXct("2026-08-25 12:00:00", tz = "UTC"),
        logInfoFn = function(...) invisible(NULL)
    )
    artifact <- saveLipidSessionManifest(
        workflow,
        session_data,
        files,
        file.path(root, "source")
    )
    manager$close()
    restored <- restoreLipidSessionManifest(artifact$latest_filepath, context)

    expect_identical(restored$session_data$current_s4_object, object)
    expect_identical(restored$session_data$assay_names, assay_names)
    expect_lt(file.info(artifact$latest_filepath)$size, 64 * 1024)
    restored$manager$close()
})
