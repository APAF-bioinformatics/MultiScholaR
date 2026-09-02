metab036SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

metab036Context <- function(root) {
    project_id <- "omics-art-036"
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
    capability$explicit_maximum_artifact_rollout <- "dual_write"
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

metab036ProjectDirs <- function(root) {
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
    list(metabolomics = paths)
}

metab036Workflow <- function(root) {
    context <- metab036Context(root)
    results <- module_ci_metab_da_results_list("mixed")
    object <- results$theObject
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
    manager$saveState(
        "metab_norm_complete",
        object,
        object@args,
        "summary final state"
    )
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- context
    workflow$state_manager <- manager
    workflow$artifact_stage_results <- list()
    workflow$contrasts_tbl <- module_ci_metab_da_contrasts("multi_group")
    workflow$design_matrix <- object@design_matrix
    workflow$config_list <- object@args
    workflow$da_ui_params <- list(q_cutoff = 0.05, formula_string = "~ 0 + group")
    workflow$normalization_ui_params <- list(normalisation_method = "quantile")
    workflow$itsd_ui_params <- list(regex = "^IS_")
    workflow$ruv_optimization_result <- list(best_k = 2L)
    persistMetabDaArtifacts(
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
    metab036SkipDependencies()
    root <- withr::local_tempdir()
    built <- metab036Workflow(root)
    withr::defer(built$manager$close())
    project_dirs <- metab036ProjectDirs(root)
    dependencies <- prepareMetabSummaryDependencies(
        built$workflow,
        project_dirs,
        "metabolomics"
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
    expect_identical(dependencies$assay_order, c("LCMS_Pos", "GCMS"))
    releaseMetabSummaryDependencies(dependencies)
    expect_identical(ls(dependencies, all.names = TRUE), character())
})

test_that("artifact summary fails clearly without required DA dependencies", {
    metab036SkipDependencies()
    root <- withr::local_tempdir()
    built <- metab036Workflow(root)
    withr::defer(built$manager$close())
    built$workflow$artifact_stage_results$differential_abundance <- NULL
    expect_error(
        prepareMetabSummaryDependencies(
            built$workflow,
            metab036ProjectDirs(root),
            "metabolomics"
        ),
        "requires current DA artifacts",
        fixed = TRUE
    )
})

test_that("artifact summary saves final products without GlobalEnv payloads", {
    metab036SkipDependencies()
    root <- withr::local_tempdir()
    built <- metab036Workflow(root)
    withr::defer(built$manager$close())
    project_dirs <- metab036ProjectDirs(root)
    values <- new.env(parent = emptyenv())
    values$workflow_args_saved <- FALSE
    output <- new.env(parent = emptyenv())
    isolated <- new.env(parent = emptyenv())
    result <- runMetabSummarySaveWorkflowArgsObserverShell(
        inputValues = list(
            experiment_label = "artifact-summary",
            description = "artifact summary test"
        ),
        projectDirs = project_dirs,
        omicType = "metabolomics",
        workflowData = built$workflow,
        values = values,
        output = output,
        createWorkflowArgsFromConfigFn = function(...) {
            path <- file.path(root, "source", "study_parameters.txt")
            writeLines("artifact summary", path)
            path
        },
        renderTextFn = function(expr) eval(substitute(expr), parent.frame()),
        showNotificationFn = function(...) invisible(NULL),
        catFn = function(...) invisible(NULL),
        globalEnv = isolated
    )

    expect_identical(result$status, "success")
    expect_false(result$configListAssigned)
    expect_false(exists("config_list", envir = isolated, inherits = FALSE))
    expect_identical(readRDS(result$s4Filepath), built$object)
    metadata_dir <- file.path(root, "results_summary", "artifact_metadata")
    expect_true(file.exists(file.path(metadata_dir, "final_s4.json")))
    expect_true(file.exists(file.path(metadata_dir, "study_parameters.json")))
})

test_that("publication copy receives explicit project context in artifact mode", {
    metab036SkipDependencies()
    root <- withr::local_tempdir()
    built <- metab036Workflow(root)
    withr::defer(built$manager$close())
    project_dirs <- metab036ProjectDirs(root)
    values <- new.env(parent = emptyenv())
    values$workflow_args_saved <- TRUE
    values$files_copied <- FALSE
    output <- new.env(parent = emptyenv())
    isolated <- new.env(parent = emptyenv())
    captured <- new.env(parent = emptyenv())
    result <- runMetabSummaryCopyToPublicationObserverShell(
        inputValues = list(
            experiment_label = "artifact-copy",
            description = "explicit context"
        ),
        projectDirs = project_dirs,
        omicType = "metabolomics",
        workflowData = built$workflow,
        values = values,
        output = output,
        copyToResultsSummaryFn = function(
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

    expect_identical(result$status, "success")
    expect_identical(captured$project_dirs, project_dirs)
    expect_false(result$projectDirsAssigned)
    expect_false(exists("project_dirs", envir = isolated, inherits = FALSE))
})

test_that("artifact report and session products retain declared dependencies", {
    metab036SkipDependencies()
    root <- withr::local_tempdir()
    built <- metab036Workflow(root)
    withr::defer(built$manager$close())
    project_dirs <- metab036ProjectDirs(root)
    template <- file.path(
        root,
        "scripts",
        "metabolomics",
        "metabolomics_report.rmd"
    )
    dir.create(dirname(template), recursive = TRUE)
    writeLines("---\ntitle: test\n---", template)
    values <- new.env(parent = emptyenv())
    values$files_copied <- TRUE
    values$report_generated <- FALSE
    output <- new.env(parent = emptyenv())
    report <- runMetabSummaryGenerateReportObserverShell(
        inputValues = list(
            experiment_label = "artifact-report",
            description = "artifact report test"
        ),
        projectDirs = project_dirs,
        omicType = "metabolomics",
        workflowData = built$workflow,
        values = values,
        output = output,
        renderReportFn = function(...) {
            path <- file.path(root, "results_summary", "metabolomics_report.html")
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
        report_path = report$renderedPath
    )
    session <- runMetabSummaryExportSessionObserverShell(
        inputValues = list(
            experiment_label = "artifact-report",
            description = "artifact report test"
        ),
        projectDirs = project_dirs,
        omicType = "metabolomics",
        values = session_values,
        workflowData = built$workflow,
        showNotificationFn = function(...) invisible(NULL),
        logInfoFn = function(...) invisible(NULL)
    )

    expect_identical(report$status, "success")
    expect_identical(session$status, "success")
    metadata_dir <- file.path(root, "results_summary", "artifact_metadata")
    expect_true(file.exists(file.path(metadata_dir, "report.json")))
    expect_true(file.exists(file.path(metadata_dir, "session_state.json")))
    report_record <- jsonlite::read_json(
        file.path(metadata_dir, "report.json"),
        simplifyVector = FALSE
    )
    expect_identical(
        report_record$dependencies$assay_order,
        as.list(c("LCMS_Pos", "GCMS"))
    )
})

test_that("artifact summary mode disables independently", {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new()
    withr::local_options(
        multischolar.metabolomics.summary_artifacts = "disabled"
    )
    expect_false(metabSummaryArtifactEligible(workflow))
    expect_true(metabSummaryGlobalOwnershipAllowed(workflow))
})
