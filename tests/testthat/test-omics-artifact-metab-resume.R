.metab032RepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

.metab032Paths <- function(root, project_id) {
    list(
        base_dir = root,
        project_id = project_id,
        omic_type = "metabolomics",
        omic_label = "metabolomics-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
}

.metab032Workflow <- function(paths, backend = "artifact") {
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        "metabolomics",
        "metabolomics-study",
        storage_policy = list(
            requested_backend = backend,
            requested_rollout = "dual_write",
            migration_requested = identical(backend, "artifact"),
            project_id = paths$project_id
        )
    )
    workflow$state_manager <- WorkflowState$new()
    workflow$artifact_stage_results <- list()
    workflow$processing_log <- list()
    workflow$tab_status <- list(
        setup_import = "pending",
        design_matrix = "disabled",
        quality_control = "disabled",
        normalization = "disabled"
    )
    workflow
}

.metab032FixtureSpec <- function(kind) {
    switch(kind,
        lc = c(
            LCMS_Pos = "metab_lc/lcms_pos_features.tsv",
            LCMS_Neg = "metab_lc/lcms_neg_features.tsv"
        ),
        gc = c(GCMS = "metab_gc/gcms_features.tsv"),
        mixed = c(
            GCMS = "metab_combined/combined_gcms_features.tsv",
            LCMS_Pos = "metab_combined/combined_lcms_features.tsv"
        )
    )
}

.metab032FixturePayload <- function(root, kind) {
    specifications <- .metab032FixtureSpec(kind)
    input_dir <- file.path(root, "inputs")
    dir.create(input_dir, recursive = TRUE)
    source_paths <- vapply(names(specifications), \(assay_name) {
        source <- .metab032RepoPath(
            "tests", "testdata", "e2e", specifications[[assay_name]]
        )
        target <- file.path(input_dir, paste0(assay_name, ".tsv"))
        stopifnot(file.copy(source, target))
        target
    }, character(1))
    assays <- lapply(source_paths, \(path) utils::read.delim(
        path,
        check.names = FALSE,
        stringsAsFactors = FALSE
    ))
    sample_columns <- grep("^(WT|KO)_", names(assays[[1L]]), value = TRUE)
    list(
        assayList = assays,
        sampleCols = sample_columns,
        sampleNamesSanitized = FALSE,
        dataFormat = "custom",
        columnMapping = list(
            metabolite_id_col = "Feature.Name",
            annotation_col = "Feature.Name",
            sample_columns = sample_columns,
            is_pattern = NA_character_
        ),
        processingLog = list(evidence_class = "independent_reviewed_fixture"),
        assaySourceFiles = as.list(source_paths),
        sourceFiles = as.list(source_paths)
    )
}

.metab032CleanAssays <- function(payload) {
    mapping <- payload$columnMapping
    metadata <- unique(c(
        mapping$metabolite_id_col,
        mapping$annotation_col
    ))
    lapply(payload$assayList, \(assay) {
        assay[unique(c(metadata, mapping$sample_columns))]
    })
}

.metab032PersistProject <- function(root, kind = "mixed", payload = NULL) {
    project_id <- paste0("metab032-", basename(root))
    paths <- .metab032Paths(root, project_id)
    workflow <- .metab032Workflow(paths)
    if (is.null(payload)) payload <- .metab032FixturePayload(root, kind)
    applyMetabImportWorkflowPayload(
        workflow,
        payload,
        logInfo = \(...) invisible(NULL)
    )
    import <- persistMetabImportArtifacts(
        workflow,
        payload,
        log_warn = \(...) invisible(NULL)
    )
    stopifnot(isTRUE(import$ok))
    samples <- payload$columnMapping$sample_columns
    groups <- payload$designGroups
    if (is.null(groups)) groups <- ifelse(grepl("^KO", samples), "KO", "WT")
    batches <- payload$designBatches
    if (is.null(batches)) batches <- rep(NA_character_, length(samples))
    replicates <- ave(seq_along(samples), groups, FUN = seq_along)
    workflow$data_cln <- .metab032CleanAssays(payload)
    workflow$design_matrix <- data.frame(
        Run = samples,
        group = groups,
        batch = batches,
        replicates = as.integer(replicates),
        tech_rep_group = paste(groups, replicates, sep = "_"),
        stringsAsFactors = FALSE
    )
    levels <- unique(groups)
    contrast <- paste0("group", levels[[2L]], "-group", levels[[1L]])
    friendly <- paste0(levels[[2L]], "_vs_", levels[[1L]])
    workflow$contrasts_tbl <- data.frame(
        contrasts = contrast,
        friendly_names = friendly,
        full_format = paste0(friendly, "=", contrast),
        stringsAsFactors = FALSE
    )
    workflow$config_list <- list(
        globalParameters = list(workflow_type = "metabolomics_standard"),
        deAnalysisParameters = list(formula_string = "~ 0 + group"),
        artifact_evidence = payload$artifactEvidence
    )
    object <- createMetaboliteAssayData(
        workflow$data_cln,
        workflow$design_matrix,
        metabolite_id_column = payload$columnMapping$metabolite_id_col,
        annotation_id_column = payload$columnMapping$annotation_col,
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        args = workflow$config_list
    )
    workflow$state_manager$setWorkflowType("metabolomics_standard")
    workflow$state_manager$saveState(
        "metab_raw_data_s4",
        object,
        workflow$config_list,
        "metabolomics resume fixture"
    )
    design <- persistMetabDesignArtifacts(
        workflow,
        log_warn = \(...) invisible(NULL)
    )
    stopifnot(isTRUE(design$ok))
    list(
        paths = paths,
        workflow = workflow,
        payload = payload,
        object = object,
        import = import,
        design = design
    )
}

.metab032FreshWorkflow <- function(paths) {
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        "metabolomics",
        "metabolomics-study",
        storage_policy = list()
    )
    workflow$state_manager <- WorkflowState$new()
    workflow$artifact_stage_results <- list()
    workflow$tab_status <- list(
        setup_import = "pending",
        design_matrix = "disabled",
        quality_control = "disabled",
        normalization = "disabled"
    )
    workflow
}

.metab032MovedPaths <- function(built) {
    destination <- tempfile("metab032-moved-")
    stopifnot(file.rename(built$paths$base_dir, destination))
    paths <- built$paths
    paths$base_dir <- destination
    paths$source_dir <- file.path(destination, "source")
    paths$results_dir <- file.path(destination, "results")
    paths
}

test_that("LC, GC, and mixed projects resume after move without sources", {
    for (kind in c("lc", "gc", "mixed")) {
        root <- withr::local_tempdir(pattern = paste0("metab032-", kind, "-"))
        built <- .metab032PersistProject(root, kind)
        artifact_files <- unlist(lapply(
            c(built$import$refs, built$design$refs),
            \(ref) file.path(root, ref$relative_path)
        ), use.names = FALSE)
        unlink(file.path(root, "inputs"), recursive = TRUE)
        moved_paths <- .metab032MovedPaths(built)
        expect_false(dir.exists(file.path(moved_paths$base_dir, "inputs")))
        workflow <- .metab032FreshWorkflow(moved_paths)
        result <- resumeMetabArtifactWorkflowSafely(
            workflow,
            moved_paths,
            "metabolomics-study",
            log_warn = \(...) invisible(NULL)
        )
        expect_true(result$resumed, info = kind)
        expect_identical(result$assay_order, names(built$object@metabolite_data))
        expect_identical(workflow$state_manager$getState(), built$object)
        expect_identical(workflow$data_cln, built$object@metabolite_data)
        expect_identical(workflow$design_matrix, built$object@design_matrix)
        expect_identical(workflow$contrasts_tbl, built$workflow$contrasts_tbl)
        expect_identical(methods::validObject(workflow$state_manager$getState()), TRUE)
        moved_files <- file.path(
            moved_paths$base_dir,
            substring(artifact_files, nchar(root) + 2L)
        )
        expect_true(all(file.exists(moved_files)))
        expect_true(result$source_payloads_retained)
        expect_false(result$eviction_performed)
        expect_true(workflow$state_manager$close())
    }
})

test_that("fresh R process reopens a moved mixed project", {
    root <- withr::local_tempdir(pattern = "metab032-fresh-")
    built <- .metab032PersistProject(root, "mixed")
    unlink(file.path(root, "inputs"), recursive = TRUE)
    paths <- .metab032MovedPaths(built)
    input <- tempfile(fileext = ".rds")
    output <- tempfile(fileext = ".rds")
    saveRDS(paths, input)
    expression <- paste(
        sprintf("devtools::load_all(%s, quiet=TRUE)", dQuote(.metab032RepoPath())),
        sprintf("paths <- readRDS(%s)", dQuote(input)),
        "w <- new.env(parent=emptyenv())",
        paste0(
            "w$workflow_context <- createWorkflowContext(paths, 'metabolomics', ",
            "'metabolomics-study', storage_policy=list())"
        ),
        "w$state_manager <- WorkflowState$new()",
        "w$artifact_stage_results <- list()",
        paste0(
            "w$tab_status <- list(setup_import='pending', ",
            "design_matrix='disabled', quality_control='disabled', ",
            "normalization='disabled')"
        ),
        paste0(
            "r <- resumeMetabArtifactWorkflowSafely(w, paths, ",
            "'metabolomics-study', log_warn=function(...) NULL)"
        ),
        paste0(
            "saveRDS(list(result=r, assays=names(w$data_cln), ",
            "valid=methods::validObject(w$state_manager$getState()), ",
            "state=w$state_manager$getCurrentStateName()), ",
            dQuote(output), ")"
        ),
        sep = ";"
    )
    child <- processx::run(
        file.path(R.home("bin"), "Rscript"),
        c("--vanilla", "-e", expression),
        wd = .metab032RepoPath(),
        error_on_status = FALSE,
        timeout = 120000
    )
    expect_identical(child$status, 0L, info = child$stderr)
    restored <- readRDS(output)
    expect_true(restored$result$resumed)
    expect_identical(restored$assays, names(built$object@metabolite_data))
    expect_identical(restored$valid, TRUE)
    expect_identical(restored$state, "metab_raw_data_s4")
})

test_that("resume lineage supports revert and branch without import eviction", {
    root <- withr::local_tempdir(pattern = "metab032-lineage-")
    built <- .metab032PersistProject(root, "mixed")
    workflow <- .metab032FreshWorkflow(built$paths)
    result <- resumeMetabArtifactWorkflowSafely(
        workflow,
        built$paths,
        "metabolomics-study",
        log_warn = \(...) invisible(NULL)
    )
    expect_true(result$resumed)
    manager <- workflow$state_manager
    original <- manager$getState()
    import_paths <- vapply(
        workflow$artifact_readthrough_refs$import,
        \(ref) file.path(built$paths$base_dir, ref$relative_path),
        character(1)
    )
    normalized <- original
    normalized@args$lineage_marker <- "normalized"
    manager$saveState(
        "metab_normalized_test",
        normalized,
        normalized@args,
        "lineage normalization branch"
    )
    expect_identical(manager$getState(), normalized)
    expect_identical(manager$revertToState("metab_raw_data_s4"), original)
    branch <- original
    branch@args$lineage_marker <- "branch"
    manager$saveState(
        "metab_branch_test",
        branch,
        branch@args,
        "lineage branch"
    )
    expect_identical(manager$getState(), branch)
    expect_true(all(file.exists(import_paths)))
    expect_true(result$source_payloads_retained)
    expect_false(result$eviction_performed)
    expect_gt(result$resource_evidence$after$rss_bytes, 0)
    expect_gt(result$resource_evidence$after$state_object_bytes, 0)
    expect_true(manager$close())
})

.metab032WorkloadPayload <- function(contract_name) {
    environment <- new.env(parent = .GlobalEnv)
    sys.source(
        .metab032RepoPath("tools", "profiling", "omics_workload_contract.R"),
        envir = environment
    )
    adapter_path <- .metab032RepoPath(
        "tools", "profiling", "omics_workload_metabolomics.R"
    )
    contract_path <- .metab032RepoPath(
        "tests", "testdata", "omics-parity", "metabolomics", "workloads",
        contract_name
    )
    contract <- environment$omicsWorkloadReadContract(contract_path)
    adapter <- environment$omicsWorkloadLoadAdapter(adapter_path, contract)
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    run_dir <- tempfile("metab032-workload-")
    dir.create(run_dir, recursive = TRUE)
    prepared <- adapter$prepare(list(contract = contract, run_dir = run_dir))
    payload_digest <- digest::digest(
        file = prepared$payload_path,
        algo = "sha256",
        serialize = FALSE
    )
    truth_digest <- digest::digest(
        file = prepared$truth_path,
        algo = "sha256",
        serialize = FALSE
    )
    stopifnot(
        identical(payload_digest, contract$expected_digests$payload_sha256),
        identical(truth_digest, contract$expected_digests$truth_sha256)
    )
    data <- utils::read.delim(
        prepared$payload_path,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        na.strings = "NA"
    )
    truth <- jsonlite::read_json(prepared$truth_path, simplifyVector = TRUE)
    assay_names <- unlist(truth$assays, use.names = FALSE)
    assays <- lapply(assay_names, \(assay_name) {
        assay <- data[data$assay == assay_name, , drop = FALSE]
        assay$assay <- NULL
        row.names(assay) <- NULL
        assay
    })
    names(assays) <- assay_names
    evidence <- list(
        workload_id = contract$workload_id,
        workload_digest = environment$omicsWorkloadDigest(contract),
        payload_sha256 = payload_digest,
        truth_sha256 = truth_digest,
        evidence_class = "generated_scaling_not_biological_validation"
    )
    list(
        assayList = assays,
        sampleCols = unlist(truth$sample_ids, use.names = FALSE),
        sampleNamesSanitized = FALSE,
        dataFormat = "custom",
        columnMapping = list(
            metabolite_id_col = "feature_id",
            annotation_col = "annotation",
            sample_columns = unlist(truth$sample_ids, use.names = FALSE),
            is_pattern = NA_character_
        ),
        processingLog = evidence,
        assaySourceFiles = stats::setNames(
            rep(list(prepared$payload_path), length(assay_names)),
            assay_names
        ),
        sourceFiles = list(assay_set = prepared$payload_path),
        designGroups = unlist(truth$group_assignments, use.names = FALSE),
        designBatches = unlist(truth$batch_assignments, use.names = FALSE),
        artifactEvidence = evidence
    )
}

test_that("frozen CI and representative workloads resume with bound evidence", {
    contracts <- c(
        "mixed-public-ci-v1.json",
        "mixed-public-representative-v1.json"
    )
    for (contract_name in contracts) {
        payload <- .metab032WorkloadPayload(contract_name)
        root <- withr::local_tempdir(pattern = "metab032-scale-")
        built <- .metab032PersistProject(root, payload = payload)
        unlink(dirname(payload$sourceFiles[[1L]]), recursive = TRUE)
        workflow <- .metab032FreshWorkflow(built$paths)
        result <- resumeMetabArtifactWorkflowSafely(
            workflow,
            built$paths,
            "metabolomics-study",
            log_warn = \(...) invisible(NULL)
        )
        expect_true(result$resumed, info = contract_name)
        expect_identical(
            workflow$state_manager$getState()@args$artifact_evidence,
            payload$artifactEvidence,
            info = contract_name
        )
        expect_identical(names(workflow$data_cln), c(
            "LCMS_Pos", "LCMS_Neg", "GCMS"
        ))
        expect_true(workflow$state_manager$close())
    }
})

test_that("consumer and apply failures preserve the prior workflow state", {
    root <- withr::local_tempdir(pattern = "metab032-failures-")
    built <- .metab032PersistProject(root, "mixed")
    categories <- unique(metabReadthroughConsumerInventory()$category)
    for (category in categories) {
        workflow <- .metab032FreshWorkflow(built$paths)
        sentinel <- workflow$state_manager
        result <- resumeMetabArtifactWorkflowSafely(
            workflow,
            built$paths,
            "metabolomics-study",
            inventory_fn = \() {
                inventory <- metabReadthroughConsumerInventory()
                inventory$verified[inventory$category == category] <- FALSE
                inventory
            },
            log_warn = \(...) invisible(NULL)
        )
        expect_false(result$ok, info = category)
        expect_identical(workflow$state_manager, sentinel, info = category)
        expect_null(workflow$data_tbl, info = category)
    }
    workflow <- .metab032FreshWorkflow(built$paths)
    sentinel <- workflow$state_manager
    result <- resumeMetabArtifactWorkflowSafely(
        workflow,
        built$paths,
        "metabolomics-study",
        failure_injector = \(stage, context) {
            if (stage == "after_resume_apply") stop("injected apply failure")
            invisible(context)
        },
        log_warn = \(...) invisible(NULL)
    )
    expect_false(result$ok)
    expect_identical(workflow$state_manager, sentinel)
    expect_null(workflow$data_tbl)

    resumed <- .metab032FreshWorkflow(built$paths)
    expect_true(resumeMetabArtifactWorkflowSafely(
        resumed,
        built$paths,
        "metabolomics-study",
        log_warn = \(...) invisible(NULL)
    )$resumed)
    state_before <- resumed$state_manager$getState()
    expect_error(
        previewMetabImportArtifact(
            resumed,
            "GCMS",
            projections = "missing_column"
        ),
        class = "multischolar_invalid_metabolomics_preview"
    )
    expect_identical(resumed$state_manager$getState(), state_before)
    invalid_design <- resumed$design_matrix[-1L, , drop = FALSE]
    design_result <- validateMetabDesignDaPreflight(
        invalid_design,
        resumed$data_cln,
        resumed$contrasts_tbl,
        "~ 0 + group",
        resumed$column_mapping
    )
    expect_false(design_result$valid)
    expect_identical(resumed$state_manager$getState(), state_before)
    expect_error(
        validateMetabQcS4FinalizeState(list()),
        "Current state is not a MetaboliteAssayData object"
    )
    expect_identical(resumed$state_manager$getState(), state_before)
    expect_error(normaliseBetweenSamples(list()))
    expect_identical(resumed$state_manager$getState(), state_before)
    invalid_session <- new.env(parent = emptyenv())
    invalid_session$state_manager <- NULL
    expect_error(buildMetabNormExportSessionData(
        invalid_session,
        list(),
        list(),
        "invalid-session"
    ))
    expect_identical(resumed$state_manager$getState(), state_before)
    expect_true(resumed$state_manager$close())
})

test_that("memory and disabled read-through leave dual-written projects untouched", {
    root <- withr::local_tempdir(pattern = "metab032-disabled-")
    built <- .metab032PersistProject(root, "gc")
    for (policy in list(
        list(requested_backend = "memory"),
        list(readthrough_enabled = FALSE)
    )) {
        workflow <- .metab032FreshWorkflow(built$paths)
        sentinel <- workflow$state_manager
        result <- resumeMetabArtifactWorkflowSafely(
            workflow,
            built$paths,
            "metabolomics-study",
            storage_policy = policy,
            log_warn = \(...) invisible(NULL)
        )
        expect_false(result$resumed)
        expect_true(result$ok)
        expect_identical(workflow$state_manager, sentinel)
        expect_null(workflow$data_tbl)
    }
    legacy_root <- withr::local_tempdir(pattern = "metab032-legacy-")
    legacy_paths <- .metab032Paths(legacy_root, "legacy-project")
    dir.create(legacy_paths$source_dir, recursive = TRUE)
    legacy <- .metab032FreshWorkflow(legacy_paths)
    result <- resumeMetabArtifactWorkflowSafely(
        legacy,
        legacy_paths,
        "metabolomics-study",
        log_warn = \(...) invisible(NULL)
    )
    expect_false(result$resumed)
    expect_identical(result$reason, "artifact_manifest_absent")
    expect_false(result$artifact_project)
    expect_false(exists("evictMetabArtifactWorkflowPayloads", mode = "function"))
})

test_that("resume defaults and consumer inventory are complete", {
    inventory <- validateMetabReadthroughConsumerInventory(
        metabReadthroughConsumerInventory()
    )
    expect_setequal(unique(inventory$category), c(
        "preview", "compatibility", "qc", "normalization", "session", "da",
        "report"
    ))
    expect_true(all(inventory$verified))
    source <- paste(readLines(
        .metab032RepoPath("R", "mod_metabolomics.R"),
        warn = FALSE
    ), collapse = "\n")
    expect_match(source, "resumeMetabArtifactWorkflowSafely", fixed = TRUE)
    expect_false(grepl("evictMetabArtifact", source, fixed = TRUE))
})
