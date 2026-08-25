.lipid041RepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

.lipid041Paths <- function(root, project_id) {
    list(
        base_dir = root,
        project_id = project_id,
        omic_type = "lipidomics",
        omic_label = "lipidomics-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
}

.lipid041Workflow <- function(paths, backend = "artifact") {
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(paths$results_dir, recursive = TRUE, showWarnings = FALSE)
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        "lipidomics",
        "lipidomics-study",
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

.lipid041FixtureSpec <- function(kind) {
    switch(kind,
        lc = c(
            LCMS_Pos = "lipid_canonical/lipidsearch_lcms_pos.txt",
            LCMS_Neg = "lipid_canonical/lipidsearch_lcms_neg.txt"
        ),
        gc = c(GCMS = "lipid_canonical/lipidsearch_gcms.txt"),
        mixed = c(
            LCMS_Pos = "lipid_canonical/lipidsearch_lcms_pos.txt",
            LCMS_Neg = "lipid_canonical/lipidsearch_lcms_neg.txt",
            GCMS = "lipid_canonical/lipidsearch_gcms.txt"
        )
    )
}

.lipid041FixturePayload <- function(root, kind) {
    specifications <- .lipid041FixtureSpec(kind)
    input_dir <- file.path(root, "inputs")
    dir.create(input_dir, recursive = TRUE)
    source_paths <- vapply(names(specifications), \(assay_name) {
        source <- .lipid041RepoPath(
            "tests", "testdata", "e2e", specifications[[assay_name]]
        )
        target <- file.path(input_dir, paste0(assay_name, ".txt"))
        stopifnot(file.copy(source, target))
        target
    }, character(1))
    assays <- lapply(source_paths, \(path) suppressMessages(
        importLipidSearchData(path)$data
    ))
    sample_columns <- grep("^(WT|KO)_", names(assays[[1L]]), value = TRUE)
    list(
        assayList = assays,
        sampleCols = sample_columns,
        sampleNamesSanitized = FALSE,
        dataFormat = "lipidsearch",
        columnMapping = list(
            lipid_id_col = "LipidName",
            annotation_col = "LipidClass",
            sample_columns = sample_columns,
            is_pattern = NA_character_
        ),
        processingLog = list(evidence_class = "independent_reviewed_fixture"),
        assaySourceFiles = as.list(source_paths),
        sourceFiles = as.list(source_paths)
    )
}

.lipid041ApplyMemory <- function(workflow, payload, ...) {
    applyLipidImportResultToWorkflow(
        workflowData = workflow,
        assayList = payload$assayList,
        dataFormat = payload$dataFormat,
        lipidIdCol = payload$columnMapping$lipid_id_col,
        annotationCol = payload$columnMapping$annotation_col,
        sampleColumns = payload$columnMapping$sample_columns,
        isPattern = payload$columnMapping$is_pattern,
        logInfo = \(...) invisible(NULL)
    )
    workflow$processing_log$setup_import <- payload$processingLog
    invisible(workflow)
}

.lipid041CleanAssays <- function(payload) {
    mapping <- payload$columnMapping
    metadata <- unique(c(
        mapping$lipid_id_col,
        mapping$annotation_col
    ))
    lapply(payload$assayList, \(assay) {
        assay[unique(c(metadata, mapping$sample_columns))]
    })
}

.lipid041PersistProject <- function(root, kind = "mixed", payload = NULL) {
    project_id <- paste0("lipid041-", basename(root))
    paths <- .lipid041Paths(root, project_id)
    workflow <- .lipid041Workflow(paths)
    if (is.null(payload)) payload <- .lipid041FixturePayload(root, kind)
    .lipid041ApplyMemory(
        workflow,
        payload,
        logInfo = \(...) invisible(NULL)
    )
    import <- persistLipidImportArtifacts(
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
    workflow$data_cln <- .lipid041CleanAssays(payload)
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
        globalParameters = list(workflow_type = "lipidomics_standard"),
        deAnalysisParameters = list(formula_string = "~ 0 + group"),
        artifact_evidence = payload$artifactEvidence
    )
    object <- createLipidomicsAssayData(
        workflow$data_cln,
        workflow$design_matrix,
        lipid_id_column = payload$columnMapping$lipid_id_col,
        annotation_id_column = payload$columnMapping$annotation_col,
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        args = workflow$config_list
    )
    workflow$state_manager$setWorkflowType("lipidomics_standard")
    workflow$state_manager$saveState(
        "lipid_raw_data_s4",
        object,
        workflow$config_list,
        "lipidomics resume fixture"
    )
    design <- persistLipidDesignArtifacts(
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

.lipid041FreshWorkflow <- function(paths) {
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        "lipidomics",
        "lipidomics-study",
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

.lipid041MovedPaths <- function(built) {
    destination <- tempfile("lipid041-moved-")
    stopifnot(file.rename(built$paths$base_dir, destination))
    paths <- built$paths
    paths$base_dir <- destination
    paths$source_dir <- file.path(destination, "source")
    paths$results_dir <- file.path(destination, "results")
    paths
}

test_that("LC, GC, and mixed projects resume after move without sources", {
    for (kind in c("lc", "gc", "mixed")) {
        root <- withr::local_tempdir(pattern = paste0("lipid041-", kind, "-"))
        built <- .lipid041PersistProject(root, kind)
        artifact_files <- unlist(lapply(
            c(built$import$refs, built$design$refs),
            \(ref) file.path(root, ref$relative_path)
        ), use.names = FALSE)
        unlink(file.path(root, "inputs"), recursive = TRUE)
        moved_paths <- .lipid041MovedPaths(built)
        expect_false(dir.exists(file.path(moved_paths$base_dir, "inputs")))
        workflow <- .lipid041FreshWorkflow(moved_paths)
        result <- resumeLipidArtifactWorkflowSafely(
            workflow,
            moved_paths,
            "lipidomics-study",
            log_warn = \(...) invisible(NULL)
        )
        expect_true(result$resumed, info = kind)
        expect_identical(result$assay_order, names(built$object@lipid_data))
        expect_identical(workflow$state_manager$getState(), built$object)
        expect_identical(workflow$data_cln, built$object@lipid_data)
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
    root <- withr::local_tempdir(pattern = "lipid041-fresh-")
    built <- .lipid041PersistProject(root, "mixed")
    unlink(file.path(root, "inputs"), recursive = TRUE)
    paths <- .lipid041MovedPaths(built)
    input <- tempfile(fileext = ".rds")
    output <- tempfile(fileext = ".rds")
    saveRDS(paths, input)
    expression <- paste(
        sprintf("devtools::load_all(%s, quiet=TRUE)", dQuote(.lipid041RepoPath())),
        sprintf("paths <- readRDS(%s)", dQuote(input)),
        "w <- new.env(parent=emptyenv())",
        paste0(
            "w$workflow_context <- createWorkflowContext(paths, 'lipidomics', ",
            "'lipidomics-study', storage_policy=list())"
        ),
        "w$state_manager <- WorkflowState$new()",
        "w$artifact_stage_results <- list()",
        paste0(
            "w$tab_status <- list(setup_import='pending', ",
            "design_matrix='disabled', quality_control='disabled', ",
            "normalization='disabled')"
        ),
        paste0(
            "r <- resumeLipidArtifactWorkflowSafely(w, paths, ",
            "'lipidomics-study', log_warn=function(...) NULL)"
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
        wd = .lipid041RepoPath(),
        error_on_status = FALSE,
        timeout = 120000
    )
    expect_identical(child$status, 0L, info = child$stderr)
    restored <- readRDS(output)
    expect_true(restored$result$resumed)
    expect_identical(restored$assays, names(built$object@lipid_data))
    expect_identical(restored$valid, TRUE)
    expect_identical(restored$state, "lipid_raw_data_s4")
})

test_that("resume lineage supports revert and branch without import eviction", {
    root <- withr::local_tempdir(pattern = "lipid041-lineage-")
    built <- .lipid041PersistProject(root, "mixed")
    workflow <- .lipid041FreshWorkflow(built$paths)
    result <- resumeLipidArtifactWorkflowSafely(
        workflow,
        built$paths,
        "lipidomics-study",
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
        "lipid_normalized_test",
        normalized,
        normalized@args,
        "lineage normalization branch"
    )
    expect_identical(manager$getState(), normalized)
    expect_identical(manager$revertToState("lipid_raw_data_s4"), original)
    branch <- original
    branch@args$lineage_marker <- "branch"
    manager$saveState(
        "lipid_branch_test",
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

.lipid041WorkloadPayload <- function(contract_name) {
    environment <- new.env(parent = .GlobalEnv)
    sys.source(
        .lipid041RepoPath("tools", "profiling", "omics_workload_contract.R"),
        envir = environment
    )
    adapter_path <- .lipid041RepoPath(
        "tools", "profiling", "omics_workload_lipidomics.R"
    )
    contract_path <- .lipid041RepoPath(
        "tests", "testdata", "omics-parity", "lipidomics", "workloads",
        contract_name
    )
    contract <- environment$omicsWorkloadReadContract(contract_path)
    adapter <- environment$omicsWorkloadLoadAdapter(adapter_path, contract)
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    run_dir <- tempfile("lipid041-workload-")
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
            lipid_id_col = "lipid_id",
            annotation_col = "lipid_class",
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

test_that("generated workloads hydrate exactly but cannot claim stage resume", {
    contracts <- c(
        "mixed-public-ci-v1.json",
        "mixed-public-representative-v1.json"
    )
    for (contract_name in contracts) {
        payload <- .lipid041WorkloadPayload(contract_name)
        root <- withr::local_tempdir(pattern = "lipid041-scale-")
        built <- .lipid041PersistProject(root, payload = payload)
        expect_false(built$import$enabled, info = contract_name)
        expect_false(built$design$enabled, info = contract_name)
        expect_null(lipidArtifactImportSpec("custom"))
        key <- list(
            project_id = built$paths$project_id,
            omic_type = "lipidomics",
            workflow_slug = "lipid_standard",
            stage_id = "design",
            state_role = "generated_scaling_smoke",
            generation_id = contract_name
        )
        restored <- hydrateLipidomicsS4Artifact(
            dehydrateLipidomicsS4Artifact(built$object, key)
        )
        expect_identical(
            restored@args$artifact_evidence,
            payload$artifactEvidence,
            info = contract_name
        )
        expect_identical(names(restored@lipid_data), c(
            "LCMS_Pos", "LCMS_Neg", "GCMS"
        ))
        workflow <- .lipid041FreshWorkflow(built$paths)
        result <- resumeLipidArtifactWorkflowSafely(
            workflow,
            built$paths,
            "lipidomics-study",
            log_warn = \(...) invisible(NULL)
        )
        expect_false(result$resumed, info = contract_name)
        expect_identical(result$reason, "artifact_manifest_absent")
    }
})

test_that("consumer and apply failures preserve the prior workflow state", {
    root <- withr::local_tempdir(pattern = "lipid041-failures-")
    built <- .lipid041PersistProject(root, "mixed")
    categories <- unique(lipidReadthroughConsumerInventory()$category)
    for (category in categories) {
        workflow <- .lipid041FreshWorkflow(built$paths)
        sentinel <- workflow$state_manager
        result <- resumeLipidArtifactWorkflowSafely(
            workflow,
            built$paths,
            "lipidomics-study",
            inventory_fn = \() {
                inventory <- lipidReadthroughConsumerInventory()
                inventory$verified[inventory$category == category] <- FALSE
                inventory
            },
            log_warn = \(...) invisible(NULL)
        )
        expect_false(result$ok, info = category)
        expect_identical(workflow$state_manager, sentinel, info = category)
        expect_null(workflow$data_tbl, info = category)
    }
    workflow <- .lipid041FreshWorkflow(built$paths)
    sentinel <- workflow$state_manager
    result <- resumeLipidArtifactWorkflowSafely(
        workflow,
        built$paths,
        "lipidomics-study",
        failure_injector = \(stage, context) {
            if (stage == "after_resume_apply") stop("injected apply failure")
            invisible(context)
        },
        log_warn = \(...) invisible(NULL)
    )
    expect_false(result$ok)
    expect_identical(workflow$state_manager, sentinel)
    expect_null(workflow$data_tbl)

    resumed <- .lipid041FreshWorkflow(built$paths)
    expect_true(resumeLipidArtifactWorkflowSafely(
        resumed,
        built$paths,
        "lipidomics-study",
        log_warn = \(...) invisible(NULL)
    )$resumed)
    state_before <- resumed$state_manager$getState()
    expect_error(
        previewLipidImportArtifact(
            resumed,
            "GCMS",
            projections = "missing_column"
        ),
        class = "multischolar_invalid_lipidomics_preview"
    )
    expect_identical(resumed$state_manager$getState(), state_before)
    invalid_design <- resumed$design_matrix[-1L, , drop = FALSE]
    design_result <- validateLipidDesignDaPreflight(
        invalid_design,
        resumed$data_cln,
        resumed$contrasts_tbl,
        "~ 0 + group",
        resumed$column_mapping
    )
    expect_false(design_result$valid)
    expect_identical(resumed$state_manager$getState(), state_before)
    expect_error(
        validateLipidQcS4FinalizeState(list()),
        "Current state is not a LipidomicsAssayData object"
    )
    expect_identical(resumed$state_manager$getState(), state_before)
    expect_error(normaliseBetweenSamples(list()))
    expect_identical(resumed$state_manager$getState(), state_before)
    invalid_session <- new.env(parent = emptyenv())
    invalid_session$state_manager <- NULL
    expect_error(handleLipidNormExportSession(
        input = list(),
        workflowData = invalid_session,
        experimentPaths = list(),
        experimentLabel = "invalid-session",
        normData = list(normalization_complete = TRUE),
        addLog = \(...) invisible(NULL),
        reqFn = \(...) stop("missing state manager")
    ))
    expect_identical(resumed$state_manager$getState(), state_before)
    expect_true(resumed$state_manager$close())
})

test_that("memory and disabled read-through leave dual-written projects untouched", {
    root <- withr::local_tempdir(pattern = "lipid041-disabled-")
    built <- .lipid041PersistProject(root, "gc")
    for (policy in list(
        list(requested_backend = "memory"),
        list(readthrough_enabled = FALSE)
    )) {
        workflow <- .lipid041FreshWorkflow(built$paths)
        sentinel <- workflow$state_manager
        result <- resumeLipidArtifactWorkflowSafely(
            workflow,
            built$paths,
            "lipidomics-study",
            storage_policy = policy,
            log_warn = \(...) invisible(NULL)
        )
        expect_false(result$resumed)
        expect_true(result$ok)
        expect_identical(workflow$state_manager, sentinel)
        expect_null(workflow$data_tbl)
    }
    legacy_root <- withr::local_tempdir(pattern = "lipid041-legacy-")
    legacy_paths <- .lipid041Paths(legacy_root, "legacy-project")
    dir.create(legacy_paths$source_dir, recursive = TRUE)
    legacy <- .lipid041FreshWorkflow(legacy_paths)
    result <- resumeLipidArtifactWorkflowSafely(
        legacy,
        legacy_paths,
        "lipidomics-study",
        log_warn = \(...) invisible(NULL)
    )
    expect_false(result$resumed)
    expect_identical(result$reason, "artifact_manifest_absent")
    expect_false(result$artifact_project)
    expect_true(exists("evictLipidArtifactWorkflowPayloads", mode = "function"))
    descriptor <- artifactLipidomicsWorkflowDescriptor()
    expect_identical(descriptor$certification$status, "dual_write")
    expect_false(descriptor$certification$auto_eligible)
})

test_that("resume defaults and consumer inventory are complete", {
    inventory <- validateLipidReadthroughConsumerInventory(
        lipidReadthroughConsumerInventory()
    )
    expect_setequal(unique(inventory$category), c(
        "preview", "compatibility", "qc", "normalization", "session", "da",
        "report"
    ))
    expect_true(all(inventory$verified))
    source <- paste(readLines(
        .lipid041RepoPath("R", "mod_lipidomics.R"),
        warn = FALSE
    ), collapse = "\n")
    expect_match(source, "resumeLipidArtifactWorkflowSafely", fixed = TRUE)
    expect_false(grepl("evictLipidArtifact", source, fixed = TRUE))
    description <- read.dcf(.lipid041RepoPath("DESCRIPTION"))
    collate <- strsplit(description[1L, "Collate"], "[[:space:]]+")[[1L]]
    collate <- gsub("^'|'$", "", collate[nzchar(collate)])
    helpers <- c(
        "mod_lipid_readthrough_registry.R",
        "mod_lipid_readthrough_helpers.R"
    )
    expect_true(all(vapply(
        helpers,
        \(helper) sum(collate == helper) == 1L,
        logical(1)
    )))
    expect_true(all(
        match(helpers, collate) < match("mod_lipidomics.R", collate)
    ))
})
