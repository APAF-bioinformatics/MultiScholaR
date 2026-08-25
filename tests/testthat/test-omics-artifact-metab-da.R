metab035SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

metab035Context <- function(root, project_id = "omics-art-035") {
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

metab035Manager <- function(context, object) {
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
        "metabolomics DA parent"
    )
    manager
}

metab035Workflow <- function(context, manager, object) {
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- context
    workflow$state_manager <- manager
    workflow$artifact_stage_results <- list()
    workflow$contrasts_tbl <- module_ci_metab_da_contrasts("multi_group")
    workflow$config_list <- object@args
    workflow
}

metab035Parameters <- function() {
    list(
        formula_string = "~ 0 + group",
        da_q_val_thresh = 0.05,
        treat_lfc_cutoff = 1,
        eBayes_trend = TRUE,
        eBayes_robust = TRUE
    )
}

test_that("metabolomics DA artifacts retain complete assay-aware results", {
    metab035SkipDependencies()
    context <- metab035Context(withr::local_tempdir())
    results <- module_ci_metab_da_results_list("mixed")
    manager <- metab035Manager(context, results$theObject)
    withr::defer(manager$close())
    workflow <- metab035Workflow(context, manager, results$theObject)
    persisted <- persistMetabDaArtifacts(
        workflow,
        results,
        workflow$contrasts_tbl,
        metab035Parameters()
    )

    expect_true(persisted$committed)
    expect_identical(
        persisted$manifest$assays,
        as.list(names(results$theObject@metabolite_data))
    )
    expect_identical(
        persisted$manifest$state_generation_id,
        manager$getCurrentGenerationId()
    )
    expect_identical(
        persisted$manifest$design_digest,
        artifactSemanticDigest(results$theObject@design_matrix)
    )
    expect_true(file.exists(artifactStoreResolveFile(
        metabDaStore(context),
        persisted$manifest$model_ref$relative_path,
        must_exist = TRUE
    )))
    expect_true(all(c(
        "logFC", "raw_pvalue", "fdr_qvalue", "significant",
        "metabolite_name", "assay", "comparison"
    ) %in% persisted$query_specification$projections))
    expect_true(file.exists(artifactStoreResolveFile(
        metabDaStore(context),
        metabDaPaths(context)$current,
        must_exist = TRUE
    )))
})

test_that("metabolomics DA queries are bounded deterministic and assay-safe", {
    metab035SkipDependencies()
    context <- metab035Context(
        withr::local_tempdir(),
        "omics-art-035-query"
    )
    results <- module_ci_metab_da_results_list("mixed")
    manager <- metab035Manager(context, results$theObject)
    withr::defer(manager$close())
    workflow <- metab035Workflow(context, manager, results$theObject)
    persistMetabDaArtifacts(
        workflow,
        results,
        workflow$contrasts_tbl,
        metab035Parameters()
    )
    projections <- c(
        "assay", "comparison", "metabolite_id", "logFC", "fdr_qvalue"
    )
    first <- queryMetabDaPage(
        workflow,
        projections = projections,
        limit = 5L
    )
    second <- queryMetabDaPage(
        workflow,
        projections = projections,
        limit = 5L,
        cursor = first$next_cursor
    )
    repeated <- queryMetabDaPage(
        workflow,
        projections = projections,
        limit = 5L
    )

    expect_identical(first$data, repeated$data)
    expect_equal(nrow(first$data), 5L)
    expect_false(any(paste(
        first$data$assay,
        first$data$comparison,
        first$data$metabolite_id
    ) %in% paste(
        second$data$assay,
        second$data$comparison,
        second$data$metabolite_id
    )))
    gc <- queryMetabDaPage(
        workflow,
        projections = projections,
        filters = list(
            assay = list(operator = "equal", value = "GCMS"),
            metabolite_id = list(operator = "contains", value = "M_")
        ),
        sort_id = "effect",
        direction = "desc",
        limit = 20L
    )
    expect_true(all(gc$data$assay == "GCMS"))
    expect_true(all(grepl("M_", gc$data$metabolite_id, fixed = TRUE)))
    expect_error(
        queryMetabDaPage(workflow, projections = "unknown", limit = 1L),
        class = "multischolar_invalid_artifact_query_projection"
    )
    expect_error(
        queryMetabDaPage(workflow, limit = 5001L),
        class = "multischolar_artifact_query_row_limit_exceeded"
    )
})

test_that("overlapping metabolite IDs remain distinct across assays", {
    metab035SkipDependencies()
    context <- metab035Context(
        withr::local_tempdir(),
        "omics-art-035-overlap"
    )
    results <- module_ci_metab_da_results_list("mixed")
    results$da_metabolites_long$metabolite_id[
        results$da_metabolites_long$assay == "GCMS"
    ] <- results$da_metabolites_long$metabolite_id[
        results$da_metabolites_long$assay == "LCMS_Pos"
    ]
    manager <- metab035Manager(context, results$theObject)
    withr::defer(manager$close())
    workflow <- metab035Workflow(context, manager, results$theObject)
    persistMetabDaArtifacts(
        workflow,
        results,
        workflow$contrasts_tbl,
        metab035Parameters()
    )
    page <- queryMetabDaPage(
        workflow,
        projections = c("assay", "comparison", "metabolite_id"),
        filters = list(metabolite_id = list(
            operator = "equal",
            value = results$da_metabolites_long$metabolite_id[[1L]]
        )),
        limit = 20L
    )

    expect_setequal(page$data$assay, c("LCMS_Pos", "GCMS"))
})

test_that("partial and interrupted DA runs cannot replace current", {
    metab035SkipDependencies()
    context <- metab035Context(
        withr::local_tempdir(),
        "omics-art-035-failure"
    )
    results <- module_ci_metab_da_results_list("mixed")
    manager <- metab035Manager(context, results$theObject)
    withr::defer(manager$close())
    workflow <- metab035Workflow(context, manager, results$theObject)
    partial <- results
    partial$da_metabolites_long <- partial$da_metabolites_long[
        partial$da_metabolites_long$assay == "LCMS_Pos",
        ,
        drop = FALSE
    ]
    expect_error(
        persistMetabDaArtifacts(
            workflow,
            partial,
            workflow$contrasts_tbl,
            metab035Parameters()
        ),
        class = "multischolar_incomplete_metab_da_run"
    )
    current_path <- artifactStoreResolveFile(
        metabDaStore(context),
        metabDaPaths(context)$current
    )
    expect_false(file.exists(current_path))
    failed <- persistMetabDaArtifactsSafely(
        workflow_data = workflow,
        results = results,
        contrasts_tbl = workflow$contrasts_tbl,
        parameters = metab035Parameters(),
        failure_injector = function(stage, context) {
            if (identical(stage, "before_metab_da_current_publication")) {
                stop("injected current failure")
            }
            invisible(context)
        }
    )
    expect_false(failed$ok)
    expect_false(file.exists(current_path))
    expect_null(workflow$artifact_stage_results$differential_abundance)
})

test_that("DA observer publishes artifacts before marking analysis complete", {
    metab035SkipDependencies()
    object <- module_ci_metab_da_object("combined")
    workflow <- module_ci_metab_da_workflow(object)
    da_data <- module_ci_metab_da_data_env()
    captured <- new.env(parent = emptyenv())
    result <- runMetabDaAnalysisObserverShell(
        object,
        module_ci_metab_da_contrasts("multi_group"),
        "~ 0 + group",
        0.05,
        1,
        workflow,
        da_data,
        session = list(),
        experimentPaths = list(),
        runAnalysis = function(...) module_ci_metab_da_results_list("mixed"),
        persistResults = function(...) {
            captured$persisted_before_complete <- !isTRUE(da_data$analysis_complete)
            list(enabled = TRUE, ok = TRUE, committed = TRUE)
        },
        updateSelectors = function(...) NULL,
        writeArtifacts = function(...) NULL,
        showNotification = function(...) invisible(NULL),
        removeNotification = function(...) invisible(NULL),
        logInfo = function(...) invisible(NULL),
        logError = function(...) invisible(NULL)
    )

    expect_identical(result$status, "success")
    expect_true(captured$persisted_before_complete)
    expect_true(da_data$analysis_complete)
})
