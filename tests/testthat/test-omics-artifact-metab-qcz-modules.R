metab033ModuleSkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

metab033ModuleContext <- function(root) {
    project_id <- "omics-art-033-modules"
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

metab033ModuleManager <- function(context) {
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

metab033ModuleWorkflow <- function(manager, context) {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$workflow_context <- context
    workflow$config_list <- list(
        globalParameters = list(workflow_type = "metabolomics_standard"),
        qc = list(mode = "artifact-module")
    )
    workflow$artifact_stage_results <- list()
    workflow$processing_log <- list()
    workflow$tab_status <- list(
        setup_import = "complete",
        design_matrix = "complete",
        quality_control = "pending",
        normalization = "disabled"
    )
    workflow$data_format <- "custom"
    workflow$data_type <- "metabolite"
    workflow
}

test_that("artifact-mode metabolomics QC modules reach final handoff", {
    metab033ModuleSkipDependencies()
    context <- metab033ModuleContext(withr::local_tempdir())
    manager <- metab033ModuleManager(context)
    withr::defer(manager$close())
    workflow <- metab033ModuleWorkflow(manager, context)
    parent <- module_ci_metab_qc_object(
        layout = "combined",
        include_duplicates = TRUE,
        include_itsd = TRUE
    )
    manager$saveState(
        "metab_raw_data_s4",
        parent,
        list(stage = "design"),
        "metabolomics QC module parent"
    )
    testthat::local_mocked_bindings(
        updateMetaboliteFiltering = function(...) NULL,
        .package = "MultiScholaR"
    )
    testthat::local_mocked_bindings(
        showNotification = function(...) invisible(NULL),
        removeNotification = function(...) invisible(NULL),
        .package = "shiny"
    )
    testthat::local_mocked_bindings(
        log_info = function(...) invisible(NULL),
        log_warn = function(...) invisible(NULL),
        log_error = function(...) invisible(NULL),
        .package = "logger"
    )

    shiny::testServer(
        mod_metab_qc_intensity_server,
        args = list(
            workflow_data = workflow,
            omic_type = "metabolomics",
            experiment_label = "artifact-module"
        ),
        {
            session$setInputs(
                intensity_cutoff_percentile = 10,
                proportion_below_cutoff = 0.5,
                apply_filter = 1
            )
        }
    )
    expect_identical(manager$getCurrentStateName(), "metab_intensity_filtered")

    shiny::testServer(
        mod_metab_qc_duplicates_server,
        args = list(
            workflow_data = workflow,
            omic_type = "metabolomics",
            experiment_label = "artifact-module"
        ),
        {
            session$setInputs(resolve_duplicates = 1)
        }
    )
    expect_identical(manager$getCurrentStateName(), "metab_duplicates_resolved")

    shiny::testServer(
        mod_metab_qc_itsd_server,
        args = list(
            workflow_data = workflow,
            omic_type = "metabolomics",
            experiment_label = "artifact-module"
        ),
        {
            session$setInputs(is_pattern = "^IS_", analyze_is = 1)
        }
    )
    expect_identical(
        manager$getCurrentStateName(),
        "metab_internal_standards_analyzed"
    )

    shiny::testServer(
        mod_metab_qc_s4_server,
        args = list(
            workflow_data = workflow,
            omic_type = "metabolomics",
            experiment_label = "artifact-module"
        ),
        {
            session$setInputs(finalize_qc = 1)
        }
    )
    expect_identical(manager$getCurrentStateName(), "metab_qc_complete")
    expect_identical(workflow$tab_status$quality_control, "complete")
    expect_identical(methods::validObject(manager$getState(), test = TRUE), TRUE)
    expect_true(all(c(
        "intensity_filter", "duplicate_resolution",
        "internal_standard_analysis", "qc_finalization"
    ) %in% names(workflow$artifact_stage_results)))
})
