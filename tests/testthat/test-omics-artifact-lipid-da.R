lipid044SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

lipid044Context <- function(root, project_id = "omics-art-044") {
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
    capability$explicit_maximum_artifact_rollout <- "dual_write"
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

lipid044Manager <- function(context, object) {
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
        "lipidomics DA parent"
    )
    manager
}

lipid044Workflow <- function(context, manager, object) {
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- context
    workflow$state_manager <- manager
    workflow$artifact_stage_results <- list()
    workflow$contrasts_tbl <- module_ci_lipid_da_contrasts("multi_group")
    workflow$config_list <- object@args
    workflow
}

lipid044Parameters <- function() {
    list(
        formula_string = "~ 0 + group",
        da_q_val_thresh = 0.05,
        treat_lfc_cutoff = 1,
        eBayes_trend = TRUE,
        eBayes_robust = TRUE
    )
}

test_that("lipidomics DA artifacts retain complete assay-aware results", {
    lipid044SkipDependencies()
    context <- lipid044Context(withr::local_tempdir())
    results <- module_ci_lipid_da_results_list("mixed")
    manager <- lipid044Manager(context, results$theObject)
    withr::defer(manager$close())
    workflow <- lipid044Workflow(context, manager, results$theObject)
    persisted <- persistLipidDaArtifacts(
        workflow,
        results,
        workflow$contrasts_tbl,
        lipid044Parameters()
    )

    expect_true(persisted$committed)
    expect_identical(
        persisted$manifest$assays,
        as.list(names(results$theObject@lipid_data))
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
        lipidDaStore(context),
        persisted$manifest$model_ref$relative_path,
        must_exist = TRUE
    )))
    expect_true(all(c(
        "logFC", "raw_pvalue", "fdr_qvalue", "significant",
        "lipid_name", "assay", "comparison"
    ) %in% persisted$query_specification$projections))
    expect_true(file.exists(artifactStoreResolveFile(
        lipidDaStore(context),
        lipidDaPaths(context)$current,
        must_exist = TRUE
    )))
})

test_that("lipidomics DA queries are bounded deterministic and assay-safe", {
    lipid044SkipDependencies()
    context <- lipid044Context(
        withr::local_tempdir(),
        "omics-art-044-query"
    )
    results <- module_ci_lipid_da_results_list("mixed")
    manager <- lipid044Manager(context, results$theObject)
    withr::defer(manager$close())
    workflow <- lipid044Workflow(context, manager, results$theObject)
    persistLipidDaArtifacts(
        workflow,
        results,
        workflow$contrasts_tbl,
        lipid044Parameters()
    )
    projections <- c(
        "assay", "comparison", "lipid_id", "logFC", "fdr_qvalue"
    )
    first <- queryLipidDaPage(
        workflow,
        projections = projections,
        limit = 5L
    )
    second <- queryLipidDaPage(
        workflow,
        projections = projections,
        limit = 5L,
        cursor = first$next_cursor
    )
    repeated <- queryLipidDaPage(
        workflow,
        projections = projections,
        limit = 5L
    )

    expect_identical(first$data, repeated$data)
    expect_equal(nrow(first$data), 5L)
    expect_false(any(paste(
        first$data$assay,
        first$data$comparison,
        first$data$lipid_id
    ) %in% paste(
        second$data$assay,
        second$data$comparison,
        second$data$lipid_id
    )))
    gc <- queryLipidDaPage(
        workflow,
        projections = projections,
        filters = list(
            assay = list(operator = "equal", value = "GCMS"),
            lipid_id = list(operator = "contains", value = "L_")
        ),
        sort_id = "effect",
        direction = "desc",
        limit = 20L
    )
    expect_true(all(gc$data$assay == "GCMS"))
    expect_true(all(grepl("L_", gc$data$lipid_id, fixed = TRUE)))
    projection_sets <- list(
        selector = c("assay", "comparison"),
        volcano = c(
            "assay", "comparison", "lipid_id", "logFC", "fdr_qvalue",
            "significant"
        ),
        heatmap = c(
            "assay", "comparison", "lipid_id",
            grep(
                "^intensity\\.",
                names(results$da_lipids_long),
                value = TRUE
            )[[1L]]
        ),
        summary = c("assay", "comparison", "significant"),
        table = names(results$da_lipids_long)
    )
    for (projection in projection_sets) {
        page <- queryLipidDaPage(
            workflow,
            projections = projection,
            limit = 10L
        )
        expect_lte(nrow(page$data), 10L)
        expect_identical(names(page$data), projection)
    }
    expect_error(
        queryLipidDaPage(workflow, projections = "unknown", limit = 1L),
        class = "multischolar_invalid_artifact_query_projection"
    )
    expect_error(
        queryLipidDaPage(workflow, limit = 5001L),
        class = "multischolar_artifact_query_row_limit_exceeded"
    )
})

test_that("overlapping lipid IDs remain distinct across assays", {
    lipid044SkipDependencies()
    context <- lipid044Context(
        withr::local_tempdir(),
        "omics-art-044-overlap"
    )
    results <- module_ci_lipid_da_results_list("mixed")
    results$da_lipids_long$lipid_id[
        results$da_lipids_long$assay == "GCMS"
    ] <- results$da_lipids_long$lipid_id[
        results$da_lipids_long$assay == "LCMS_Pos"
    ]
    manager <- lipid044Manager(context, results$theObject)
    withr::defer(manager$close())
    workflow <- lipid044Workflow(context, manager, results$theObject)
    persistLipidDaArtifacts(
        workflow,
        results,
        workflow$contrasts_tbl,
        lipid044Parameters()
    )
    page <- queryLipidDaPage(
        workflow,
        projections = c("assay", "comparison", "lipid_id"),
        filters = list(lipid_id = list(
            operator = "equal",
            value = results$da_lipids_long$lipid_id[[1L]]
        )),
        limit = 20L
    )

    expect_setequal(page$data$assay, c("LCMS_Pos", "LCMS_Neg", "GCMS"))
})

test_that("partial and interrupted DA runs cannot replace current", {
    lipid044SkipDependencies()
    context <- lipid044Context(
        withr::local_tempdir(),
        "omics-art-044-failure"
    )
    results <- module_ci_lipid_da_results_list("mixed")
    manager <- lipid044Manager(context, results$theObject)
    withr::defer(manager$close())
    workflow <- lipid044Workflow(context, manager, results$theObject)
    partial <- results
    partial$da_lipids_long <- partial$da_lipids_long[
        partial$da_lipids_long$assay == "LCMS_Pos",
        ,
        drop = FALSE
    ]
    expect_error(
        persistLipidDaArtifacts(
            workflow,
            partial,
            workflow$contrasts_tbl,
            lipid044Parameters()
        ),
        class = "multischolar_incomplete_lipid_da_run"
    )
    current_path <- artifactStoreResolveFile(
        lipidDaStore(context),
        lipidDaPaths(context)$current
    )
    expect_false(file.exists(current_path))
    failed <- persistLipidDaArtifactsSafely(
        workflow_data = workflow,
        results = results,
        contrasts_tbl = workflow$contrasts_tbl,
        parameters = lipid044Parameters(),
        failure_injector = function(stage, context) {
            if (identical(stage, "before_lipid_da_current_publication")) {
                stop("injected current failure")
            }
            invisible(context)
        }
    )
    expect_false(failed$ok)
    expect_false(file.exists(current_path))
    expect_null(workflow$artifact_stage_results$differential_abundance)
})

test_that("frozen CI workload persists complete bounded DA projections", {
    lipid044SkipDependencies()
    repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    adapter <- new.env(parent = globalenv())
    sys.source(
        file.path(repo, "tools", "profiling", "omics_workload_lipidomics.R"),
        envir = adapter
    )
    contract <- jsonlite::read_json(
        file.path(
            repo,
            "tests",
            "testdata",
            "omics-parity",
            "lipidomics",
            "workloads",
            "mixed-public-ci-v1.json"
        ),
        simplifyVector = FALSE
    )
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    root <- withr::local_tempdir()
    prepared <- adapter$lipidWorkloadPrepare(list(
        contract = contract,
        run_dir = root
    ))
    expect_identical(
        digest::digest(file = prepared$payload_path, algo = "sha256"),
        contract$expected_digests$payload_sha256
    )
    expect_identical(
        digest::digest(file = prepared$truth_path, algo = "sha256"),
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
        assay <- generated[generated$assay == assay_name, , drop = FALSE]
        assay$assay <- NULL
        assay$lipid_id <- make.unique(as.character(assay$lipid_id))
        rownames(assay) <- NULL
        assay
    })
    names(assays) <- assay_names
    design <- data.frame(
        Run = samples,
        group = unlist(truth$group_assignments, use.names = FALSE),
        batch = unlist(truth$batch_assignments, use.names = FALSE),
        tech_rep_group = rep(
            paste0("pair_", seq_len(length(samples) / 2L)),
            each = 2L
        ),
        stringsAsFactors = FALSE
    )
    object <- createLipidomicsAssayData(
        assays,
        design,
        lipid_id_column = "lipid_id",
        annotation_id_column = "lipid_class",
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        args = list(evidence_class = "generated_scaling")
    )
    contrasts <- module_ci_lipid_da_contrasts("multi_group")
    rows <- lapply(names(assays), function(assay_name) {
        assay <- assays[[assay_name]]
        do.call(rbind, lapply(seq_len(nrow(contrasts)), function(index) {
            raw <- seq_len(nrow(assay)) / (nrow(assay) + 1)
            data.frame(
                lipid_id = assay$lipid_id,
                lipid_name = assay$lipid_class,
                assay = assay_name,
                comparison = contrasts$contrasts[[index]],
                friendly_name = contrasts$friendly_names[[index]],
                logFC = seq(-2, 2, length.out = nrow(assay)),
                raw_pvalue = raw,
                fdr_qvalue = stats::p.adjust(raw, method = "BH"),
                significant = "NS",
                numerator = sub("-.*$", "", contrasts$contrasts[[index]]),
                denominator = sub("^.*-", "", contrasts$contrasts[[index]]),
                stringsAsFactors = FALSE
            )
        }))
    })
    results <- list(
        theObject = object,
        contrasts_results = list(),
        da_lipids_long = do.call(rbind, rows)
    )
    context <- lipid044Context(root, "omics-art-044-frozen-ci")
    manager <- lipid044Manager(context, object)
    withr::defer(manager$close())
    workflow <- lipid044Workflow(context, manager, object)
    persisted <- persistLipidDaArtifacts(
        workflow,
        results,
        contrasts,
        lipid044Parameters()
    )
    projections <- c(
        "assay", "comparison", "lipid_id", "logFC", "fdr_qvalue"
    )
    for (limit in c(1L, 7L, 50L)) {
        page <- queryLipidDaPage(
            workflow,
            projections = projections,
            sort_id = "effect",
            direction = "desc",
            limit = limit
        )
        expect_lte(nrow(page$data), limit)
    }
    expect_true(persisted$committed)
    expect_identical(
        sort(unique(results$da_lipids_long$assay)),
        sort(assay_names)
    )
})

test_that("DA observer publishes artifacts before marking analysis complete", {
    lipid044SkipDependencies()
    object <- module_ci_lipid_da_object("triple")
    workflow <- module_ci_lipid_da_workflow(object)
    da_data <- module_ci_lipid_da_data_env()
    captured <- new.env(parent = emptyenv())
    result <- executeLipidDaRunAnalysis(
        currentS4 = object,
        contrastsTbl = module_ci_lipid_da_contrasts("multi_group"),
        formulaString = "~ 0 + group",
        daQValThresh = 0.05,
        treatLfcCutoff = 1,
        workflowData = workflow,
        daData = da_data,
        session = list(),
        experimentPaths = list(),
        runLipidsDaFn = function(...) module_ci_lipid_da_results_list("mixed"),
        persistResultsFn = function(...) {
            captured$persisted_before_complete <- !isTRUE(da_data$analysis_complete)
            list(enabled = TRUE, ok = TRUE, committed = TRUE)
        },
        successFinalizerFn = function(...) {
            expect_true(captured$persisted_before_complete)
            da_data$analysis_complete <- TRUE
            list(status = "finalized")
        },
        errorFinalizerFn = function(errorMessage) stop(errorMessage)
    )

    expect_identical(result$status, "success")
    expect_true(captured$persisted_before_complete)
    expect_true(da_data$analysis_complete)
})
