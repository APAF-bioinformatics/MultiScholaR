lipid045SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

lipid045Context <- function(root) {
    project_id <- "omics-art-045"
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

lipid045ProjectDirs <- function(root) {
    paths <- list(
        base_dir = root,
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results"),
        results_summary_dir = file.path(root, "results_summary"),
        publication_graphs_dir = file.path(root, "publication_graphs"),
        time_dir = file.path(root, "time"),
        qc_dir = file.path(root, "qc"),
        da_output_dir = file.path(root, "da_output"),
        pathway_dir = file.path(root, "pathway"),
        feature_qc_dir = file.path(root, "feature_qc"),
        integration_dir = file.path(root, "integration")
    )
    invisible(lapply(paths[names(paths) != "base_dir"], dir.create,
        recursive = TRUE,
        showWarnings = FALSE
    ))
    list(lipidomics = paths)
}

lipid045Workflow <- function(root) {
    context <- lipid045Context(root)
    results <- module_ci_lipid_da_results_list("mixed")
    object <- results$theObject
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
    manager$saveState(
        "lipid_norm_complete",
        object,
        object@args,
        "summary final state"
    )
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- context
    workflow$state_manager <- manager
    workflow$artifact_stage_results <- list()
    workflow$contrasts_tbl <- module_ci_lipid_da_contrasts("multi_group")
    workflow$design_matrix <- object@design_matrix
    workflow$config_list <- object@args
    workflow$da_ui_params <- list(q_cutoff = 0.05, formula_string = "~ 0 + group")
    workflow$normalization_ui_params <- list(normalisation_method = "quantile")
    workflow$itsd_ui_params <- list(regex = "^IS_")
    workflow$ruv_optimization_result <- list(best_k = 2L)
    persistLipidDaArtifacts(
        workflow,
        results,
        workflow$contrasts_tbl,
        parameters = list(
            formula_string = "~ 0 + group",
            da_q_val_thresh = 0.05,
            treat_lfc_cutoff = 1
        )
    )
    list(workflow = workflow, manager = manager, object = object)
}

test_that("artifact summary dependencies are exact and releasable", {
    lipid045SkipDependencies()
    root <- withr::local_tempdir()
    built <- lipid045Workflow(root)
    withr::defer(built$manager$close())
    project_dirs <- lipid045ProjectDirs(root)
    dependencies <- prepareLipidSummaryDependencies(
        built$workflow,
        project_dirs,
        "lipidomics"
    )

    expect_true(dependencies$enabled)
    expect_identical(dependencies$object, built$object)
    expect_identical(
        dependencies$evidence$state_generation_id,
        built$manager$getCurrentGenerationId()
    )
    expect_identical(
        dependencies$evidence$da_run_id,
        built$workflow$artifact_stage_results$differential_abundance$run_id
    )
    expect_identical(dependencies$assay_order, c("LCMS_Pos", "LCMS_Neg", "GCMS"))
    releaseLipidSummaryDependencies(dependencies)
    expect_identical(ls(dependencies, all.names = TRUE), character())
})

test_that("artifact summary fails clearly without required DA dependencies", {
    lipid045SkipDependencies()
    root <- withr::local_tempdir()
    built <- lipid045Workflow(root)
    withr::defer(built$manager$close())
    built$workflow$artifact_stage_results$differential_abundance <- NULL
    expect_error(
        prepareLipidSummaryDependencies(
            built$workflow,
            lipid045ProjectDirs(root),
            "lipidomics"
        ),
        "requires current DA artifacts",
        fixed = TRUE
    )
})

test_that("artifact summary saves final products without GlobalEnv payloads", {
    lipid045SkipDependencies()
    root <- withr::local_tempdir()
    built <- lipid045Workflow(root)
    withr::defer(built$manager$close())
    project_dirs <- lipid045ProjectDirs(root)
    values <- new.env(parent = emptyenv())
    values$workflow_args_saved <- FALSE
    output <- new.env(parent = emptyenv())
    isolated <- new.env(parent = emptyenv())
    result <- handleLipidSummarySaveWorkflowArgs(
        input = list(
            experiment_label = "artifact-summary",
            description = "artifact summary test"
        ),
        projectDirs = project_dirs,
        omicType = "lipidomics",
        workflowData = built$workflow,
        values = values,
        output = output,
        collectContextFn = function(workflowData, catFn) {
            collectLipidSummaryWorkflowArgsContext(
                workflowData,
                globalEnv = isolated,
                catFn = catFn
            )
        },
        createWorkflowArgsFn = function(...) {
            path <- file.path(root, "source", "study_parameters.txt")
            writeLines("artifact summary", path)
            path
        },
        renderTextFn = function(expr) eval(substitute(expr), parent.frame()),
        showNotificationFn = function(...) invisible(NULL),
        catFn = function(...) invisible(NULL)
    )

    expect_true(is.list(result))
    expect_false(exists("config_list", envir = isolated, inherits = FALSE))
    expect_identical(readRDS(result$s4Filepath), built$object)
    metadata_dir <- file.path(root, "results_summary", "artifact_metadata")
    expect_true(file.exists(file.path(metadata_dir, "final_s4.json")))
    expect_true(file.exists(file.path(metadata_dir, "study_parameters.json")))
})

test_that("publication copy receives explicit project context in artifact mode", {
    lipid045SkipDependencies()
    root <- withr::local_tempdir()
    built <- lipid045Workflow(root)
    withr::defer(built$manager$close())
    project_dirs <- lipid045ProjectDirs(root)
    values <- new.env(parent = emptyenv())
    values$workflow_args_saved <- TRUE
    values$files_copied <- FALSE
    output <- new.env(parent = emptyenv())
    isolated <- new.env(parent = emptyenv())
    captured <- new.env(parent = emptyenv())
    result <- handleLipidSummaryCopyToPublication(
        input = list(
            experiment_label = "artifact-copy",
            description = "explicit context"
        ),
        projectDirs = project_dirs,
        omicType = "lipidomics",
        workflowData = built$workflow,
        values = values,
        output = output,
        copyResultsSummaryFn = function(
            omic_type,
            experiment_label,
            contrasts_tbl,
            design_matrix,
            force,
            project_dirs
        ) {
            captured$project_dirs <- project_dirs
            list()
        },
        renderTextFn = function(expr) eval(substitute(expr), parent.frame()),
        showNotificationFn = function(...) invisible(NULL),
        withProgressFn = function(message = NULL, expr, ...) {
            eval(substitute(expr), parent.frame())
        },
        catFn = function(...) invisible(NULL),
        globalEnv = isolated
    )

    expect_true(is.list(result))
    expect_identical(captured$project_dirs, project_dirs)
    expect_false(exists("project_dirs", envir = isolated, inherits = FALSE))
})

test_that("artifact report and session products retain declared dependencies", {
    lipid045SkipDependencies()
    root <- withr::local_tempdir()
    built <- lipid045Workflow(root)
    withr::defer(built$manager$close())
    project_dirs <- lipid045ProjectDirs(root)
    template <- file.path(
        root,
        "scripts",
        "lipidomics",
        "lipidomics_report.rmd"
    )
    dir.create(dirname(template), recursive = TRUE)
    writeLines("---\ntitle: test\n---", template)
    values <- new.env(parent = emptyenv())
    values$files_copied <- TRUE
    values$report_generated <- FALSE
    output <- new.env(parent = emptyenv())
    report <- handleLipidSummaryGenerateReport(
        input = list(
            experiment_label = "artifact-report",
            description = "artifact report test"
        ),
        projectDirs = project_dirs,
        omicType = "lipidomics",
        workflowData = built$workflow,
        values = values,
        output = output,
        renderReportFn = function(...) {
            path <- file.path(root, "results_summary", "lipidomics_report.html")
            writeLines("<html>artifact report</html>", path)
            path
        },
        reactiveFn = function(expr) eval(substitute(expr), parent.frame()),
        outputOptionsFn = function(...) invisible(NULL),
        downloadHandlerFn = function(filename, content) list(
            filename = filename,
            content = content
        ),
        renderTextFn = function(expr) eval(substitute(expr), parent.frame()),
        withProgressFn = function(message = NULL, expr, ...) {
            eval(substitute(expr), parent.frame())
        },
        incProgressFn = function(...) invisible(NULL),
        showNotificationFn = function(...) invisible(NULL),
        logInfoFn = function(...) invisible(NULL),
        catFn = function(...) invisible(NULL)
    )
    session_values <- list(
        workflow_args_saved = TRUE,
        files_copied = TRUE,
        report_generated = TRUE,
        report_path = values$report_path
    )
    session <- handleLipidSummaryExportSessionState(
        input = list(
            experiment_label = "artifact-report",
            description = "artifact report test"
        ),
        projectDirs = project_dirs,
        omicType = "lipidomics",
        values = session_values,
        workflowData = built$workflow,
        showNotificationFn = function(...) invisible(NULL),
        logInfoFn = function(...) invisible(NULL)
    )

    expect_null(report)
    expect_true(file.exists(session))
    metadata_dir <- file.path(root, "results_summary", "artifact_metadata")
    expect_true(file.exists(file.path(metadata_dir, "report_html.json")))
    expect_true(file.exists(file.path(metadata_dir, "session_state.json")))
    report_record <- jsonlite::read_json(
        file.path(metadata_dir, "report_html.json"),
        simplifyVector = FALSE
    )
    expect_identical(
        report_record$dependencies$assay_order,
        as.list(c("LCMS_Pos", "LCMS_Neg", "GCMS"))
    )
})

test_that("artifact summary mode disables independently", {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new()
    withr::local_options(
        multischolar.lipidomics.summary_artifacts = "disabled"
    )
    expect_false(lipidSummaryArtifactEligible(workflow))
    expect_true(lipidSummaryGlobalOwnershipAllowed(workflow))
})
