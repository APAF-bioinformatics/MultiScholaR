.lipid040RepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

.lipid040Paths <- function(root, project_id, omic_type = "lipidomics") {
    paths <- list(
        base_dir = root,
        project_id = project_id,
        omic_type = omic_type,
        omic_label = "lipidomics-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE)
    dir.create(paths$results_dir, recursive = TRUE)
    paths
}

.lipid040Workflow <- function(
    root,
    backend = "artifact",
    omic_type = "lipidomics",
    project_id = paste0("lipid040-", basename(root))
) {
    paths <- .lipid040Paths(root, project_id, omic_type)
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        omic_type,
        "lipidomics-study",
        storage_policy = list(
            requested_backend = backend,
            requested_rollout = "dual_write",
            migration_requested = identical(backend, "artifact"),
            project_id = project_id
        )
    )
    workflow$state_manager <- WorkflowState$new()
    workflow$artifact_stage_results <- list()
    workflow$processing_log <- list()
    workflow$tab_status <- list(
        setup_import = "pending",
        design_matrix = "disabled",
        quality_control = "disabled"
    )
    list(workflow = workflow, paths = paths)
}

.lipid040ReadTable <- function(context, ref) {
    identity <- context$getIdentity()
    store <- newArtifactStore(context$getPaths(), identity$project_id)
    ref <- artifactStoreNormalizeRef(ref)
    managed <- artifactStoreManagedPaths(store, ref$logical_key, ref$artifact_id)
    sidecar <- artifactStoreReadSidecar(store, managed$sidecar, TRUE)
    payload <- arrow::read_parquet(
        artifactStoreResolveFile(store, ref$relative_path, TRUE),
        as_data_frame = FALSE
    )
    decodeArtifactRectangular(payload, sidecar$codec_metadata)
}

.lipid076ExactEnvelope <- function() {
    envelope <- jsonlite::read_json(
        .lipid040RepoPath(
            "inst", "extdata", "omics-auto-policy-receipts-v2.json"
        ),
        simplifyVector = FALSE
    )
    receipt <- Filter(\(candidate) identical(
        candidate$capability_id,
        "lipidomics.lipidsearch.lipid.standard.v1"
    ), envelope$receipts)[[1L]]
    receipt$receipt_id <- "test.lipidomics.lipidsearch.preingress.v1"
    receipt$receipt_kind <- "proposed_pilot"
    receipt$owner_ticket_id <- "OMICS-ART-076"
    receipt$decision <- "proposed_pilot"
    receipt$size_measure <- list(
        measure_id = "total_uncompressed_input_bytes_v1",
        unit = "byte",
        exact = TRUE,
        available_before_full_parse = TRUE
    )
    receipt$threshold_bytes <- 1
    receipt$receipt_digest <- NULL
    receipt$receipt_digest <- workflowPolicyObjectDigest(receipt)
    envelope$receipts <- list(receipt)
    envelope$envelope_digest <- NULL
    envelope$envelope_digest <- workflowPolicyObjectDigest(envelope)
    envelope
}

.lipid040FixturePayload <- function(kind = c("lc", "gc", "mixed")) {
    kind <- match.arg(kind)
    specifications <- switch(kind,
        lc = c(
            LCMS_Pos = paste0(
                "tests/testdata/e2e/lipid_canonical/",
                "lipidsearch_lcms_pos.txt"
            ),
            LCMS_Neg = paste0(
                "tests/testdata/e2e/lipid_canonical/",
                "lipidsearch_lcms_neg.txt"
            )
        ),
        gc = c(
            GCMS = paste0(
                "tests/testdata/e2e/lipid_canonical/",
                "lipidsearch_gcms.txt"
            )
        ),
        mixed = c(
            LCMS_Pos = paste0(
                "tests/testdata/e2e/lipid_canonical/",
                "lipidsearch_lcms_pos.txt"
            ),
            LCMS_Neg = paste0(
                "tests/testdata/e2e/lipid_canonical/",
                "lipidsearch_lcms_neg.txt"
            ),
            GCMS = paste0(
                "tests/testdata/e2e/lipid_canonical/",
                "lipidsearch_gcms.txt"
            )
        )
    )
    source_paths <- vapply(
        specifications,
        .lipid040RepoPath,
        character(1)
    )
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
        processingLog = list(
            evidence_class = "independent_reviewed_fixture",
            fixture_kind = kind
        ),
        assaySourceFiles = as.list(source_paths),
        sourceFiles = as.list(source_paths)
    )
}

.lipid040ApplyMemory <- function(workflow, payload) {
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

.lipid040ApplyImport <- function(workflow, payload) {
    .lipid040ApplyMemory(workflow, payload)
    persistLipidImportArtifacts(
        workflow,
        payload,
        log_warn = \(...) invisible(NULL)
    )
}

.lipid040CleanAssays <- function(payload) {
    mapping <- payload$columnMapping
    metadata <- unique(c(
        mapping$lipid_id_col,
        mapping$annotation_col
    ))
    metadata <- metadata[!is.na(metadata) & nzchar(metadata)]
    lapply(payload$assayList, \(assay) {
        assay[c(metadata, mapping$sample_columns)]
    })
}

.lipid040Design <- function(workflow, payload, evidence = NULL) {
    samples <- payload$columnMapping$sample_columns
    groups <- payload$designGroups
    if (is.null(groups)) groups <- ifelse(grepl("^KO", samples), "KO", "WT")
    batches <- payload$designBatches
    if (is.null(batches)) batches <- rep(NA_character_, length(samples))
    replicates <- ave(seq_along(samples), groups, FUN = seq_along)
    workflow$data_cln <- .lipid040CleanAssays(payload)
    workflow$design_matrix <- data.frame(
        Run = samples,
        group = groups,
        batch = batches,
        replicates = as.integer(replicates),
        tech_rep_group = paste(groups, replicates, sep = "_"),
        stringsAsFactors = FALSE
    )
    group_levels <- unique(groups)
    stopifnot(length(group_levels) >= 2L)
    contrast <- paste0("group", group_levels[[2L]], "-group", group_levels[[1L]])
    friendly <- paste0(group_levels[[2L]], "_vs_", group_levels[[1L]])
    workflow$contrasts_tbl <- data.frame(
        contrasts = contrast,
        friendly_names = friendly,
        full_format = paste0(friendly, "=", contrast),
        stringsAsFactors = FALSE
    )
    workflow$config_list <- list(
        globalParameters = list(workflow_type = "lipidomics_standard"),
        deAnalysisParameters = list(formula_string = "~ 0 + group"),
        artifact_evidence = evidence
    )
    object <- createLipidomicsAssayData(
        lipid_data = workflow$data_cln,
        design_matrix = workflow$design_matrix,
        lipid_id_column = payload$columnMapping$lipid_id_col,
        annotation_id_column = payload$columnMapping$annotation_col,
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        database_identifier_type = "FixtureFeatureName",
        internal_standard_regex = NA_character_,
        args = workflow$config_list
    )
    workflow$state_manager$setWorkflowType("lipidomics_standard")
    workflow$state_manager$saveState(
        "lipid_raw_data_s4",
        object,
        workflow$config_list,
        "lipidomics design memory checkpoint"
    )
    object
}

.lipid040Registry <- function(context) {
    registry <- projectRegistryForContext(context)
    openProjectRegistryReadOnly(registry)
}

test_that("only reviewed LipidSearch receives an evict descriptor", {
    descriptor <- artifactLipidomicsWorkflowDescriptor()
    expect_identical(
        descriptor$descriptor_id,
        "lipidomics.lipidsearch.lipid.standard.v1"
    )
    expect_identical(descriptor$certification$status, "evict")
    expect_false(descriptor$certification$auto_eligible)
    expect_identical(
        names(descriptor$codecs),
        "multischolar.s4.lipidomics_assay_data"
    )
    expect_false(any(vapply(
        artifactDescriptorCatalogueValues(artifactWorkflowDescriptorCatalogue()),
        \(candidate) candidate$identity$input_format %in% c("msdial", "custom") &&
            identical(candidate$identity$omic_type, "lipidomics"),
        logical(1)
    )))
    support <- workflowFormatSupportCatalogue()
    lipidsearch <- support[[which(vapply(
        support,
        \(entry) identical(entry$format_id, "lipidomics.lipidsearch"),
        logical(1)
    ))]]
    custom <- support[[which(vapply(
        support,
        \(entry) identical(entry$format_id, "lipidomics.custom"),
        logical(1)
    ))]]
    msdial <- support[[which(vapply(
        support,
        \(entry) identical(entry$format_id, "lipidomics.msdial"),
        logical(1)
    ))]]
    expect_identical(lipidsearch$support_status, "scientifically_supported")
    expect_identical(custom$support_status, "reader_characterized")
    expect_identical(msdial$support_status, "reader_characterized")
})

test_that("exact LipidSearch auto routing binds before the full reader", {
    root <- withr::local_tempdir(pattern = "lipid076-preingress-")
    built <- .lipid040Workflow(root, backend = "auto")
    source <- file.path(root, "source", "lipidsearch.txt")
    expect_true(file.copy(
        .lipid040RepoPath(paste0(
            "tests/testdata/e2e/lipid_canonical/",
            "lipidsearch_lcms_pos.txt"
        )),
        source
    ))
    observed_n_max <- NULL
    withr::local_options(list(
        MultiScholaR.lipidomics.preingress_envelope =
            .lipid076ExactEnvelope()
    ))
    preview <- loadLipidImportAssayPreview(
        source,
        vendorFormat = "lipidsearch",
        workflowData = built$workflow,
        importLipidSearch = \(path, n_max = Inf) {
            observed_n_max <<- n_max
            importLipidSearchData(path, n_max = n_max)
        }
    )

    expect_true(preview$deferred)
    expect_identical(observed_n_max, 1000L)
    expect_true(built$workflow$workflow_context$isBound())
    expect_identical(
        built$workflow$workflow_context$getStorageDecision()$effective_backend,
        "artifact"
    )
    expect_identical(
        preview$preIngress$outcome$token$measure$bytes,
        unname(as.numeric(file.info(source)$size))
    )
    expect_false(preview$preIngress$outcome$token$probe_evidence[[
        "complete_payload_materialized"
    ]])
})

test_that("installed LipidSearch auto policy remains post-parse compatible", {
    root <- withr::local_tempdir(pattern = "lipid076-legacy-auto-")
    built <- .lipid040Workflow(root, backend = "auto")
    source <- .lipid040RepoPath(paste0(
        "tests/testdata/e2e/lipid_canonical/",
        "lipidsearch_lcms_pos.txt"
    ))
    preview <- loadLipidImportAssayPreview(
        source,
        vendorFormat = "lipidsearch",
        workflowData = built$workflow
    )
    expect_false(preview$deferred)
    expect_identical(
        preview$preIngress$status,
        "installed_legacy_policy_deferred"
    )
    expect_false(built$workflow$workflow_context$isBound())
    expect_gt(nrow(preview$assayData), 0L)
})

test_that("deferred LipidSearch processing preserves the public payload", {
    root <- withr::local_tempdir(pattern = "lipid076-processing-")
    built <- .lipid040Workflow(root)
    source <- .lipid040RepoPath(paste0(
        "tests/testdata/e2e/lipid_canonical/",
        "lipidsearch_lcms_pos.txt"
    ))
    preview <- loadLipidImportAssayPreview(
        source,
        vendorFormat = "lipidsearch",
        workflowData = built$workflow
    )
    result <- runLipidImportProcessing(
        workflowData = built$workflow,
        assay1Name = "LCMS_Pos",
        assay1Data = preview$assayData,
        assay1File = source,
        vendorFormat = "lipidsearch",
        detectedFormat = "lipidsearch",
        lipidIdCol = "LipidName",
        annotationCol = "LipidClass",
        sampleColumns = preview$importResult$sample_columns,
        isPattern = "",
        sanitizeNames = FALSE,
        experimentPaths = built$paths,
        writeImportArtifacts = \(...) list(written = FALSE),
        notify = \(...) invisible(NULL),
        removeNotify = \(...) invisible(NULL),
        assay1Deferred = preview$deferred
    )
    expected <- importLipidSearchData(source)$data

    expect_identical(result$status, "success")
    expect_true(result$artifactStageResult$committed)
    expect_identical(result$assayList$LCMS_Pos, expected)
    expect_identical(built$workflow$data_tbl$LCMS_Pos, expected)
})

test_that("reviewed LipidSearch assays dual-write exact independent provenance", {
    for (kind in c("lc", "gc", "mixed")) {
        root <- withr::local_tempdir(pattern = paste0("lipid040-", kind, "-"))
        built <- .lipid040Workflow(root)
        payload <- .lipid040FixturePayload(kind)
        memory_before <- unserialize(serialize(payload$assayList, NULL))
        output <- .lipid040ApplyImport(built$workflow, payload)
        expect_true(output$enabled, info = kind)
        expect_true(output$ok, info = kind)
        expect_true(output$committed, info = kind)
        expect_identical(built$workflow$data_tbl, memory_before, info = kind)
        expect_identical(names(output$assay_refs), names(memory_before))
        expect_identical(output$assay_manifest$assay_name, names(memory_before))
        expect_identical(output$assay_manifest$assay_order, seq_along(memory_before))
        expect_true(all(c(
            "assay_identity", "artifact_role", "source_digest", "parser_id",
            "parser_version", "format_id", "mapping_digest", "table_digest",
            "feature_count", "sample_count", "lipid_id_column",
            "lipid_class_column", "sample_columns_json", "column_schema_json"
        ) %in% names(output$assay_manifest)))
        expect_identical(anyDuplicated(
            output$assay_manifest$assay_identity
        ), 0L)
        expect_identical(anyDuplicated(vapply(
            output$assay_refs,
            `[[`,
            character(1),
            "artifact_id"
        )), 0L)
        for (assay_name in names(memory_before)) {
            expect_identical(
                .lipid040ReadTable(
                    built$workflow$workflow_context,
                    output$assay_refs[[assay_name]]
                ),
                memory_before[[assay_name]],
                info = paste(kind, assay_name)
            )
        }
        session <- .lipid040Registry(built$workflow$workflow_context)
        identity <- built$workflow$workflow_context$getIdentity()
        sources <- projectRegistryQuery(
            session,
            "sources",
            filters = list(workflow_id = identity$workflow_id, run_id = output$run_id)
        )
        expect_true(is.na(sources$source_uri))
        expect_identical(sources$parser_id, unique(
            output$assay_manifest$parser_id
        ))
        expect_identical(sources$format_id, "lipidsearch.direct_file")
        closeProjectRegistry(session)
    }
})

test_that("import coordinator invokes additive dual-write after memory handoff", {
    root <- withr::local_tempdir(pattern = "lipid040-import-hook-")
    built <- .lipid040Workflow(root)
    payload <- .lipid040FixturePayload("lc")
    sources <- unlist(payload$assaySourceFiles, use.names = TRUE)
    result <- runLipidImportProcessing(
        workflowData = built$workflow,
        assay1Name = "LCMS_Pos",
        assay1Data = payload$assayList$LCMS_Pos,
        assay1File = sources[["LCMS_Pos"]],
        assay2File = sources[["LCMS_Neg"]],
        assay2Name = "LCMS_Neg",
        vendorFormat = "lipidsearch",
        detectedFormat = "lipidsearch",
        lipidIdCol = "LipidName",
        annotationCol = "LipidClass",
        sampleColumns = payload$sampleCols,
        isPattern = NA_character_,
        sanitizeNames = FALSE,
        experimentPaths = built$paths,
        writeImportArtifacts = \(...) list(written = TRUE, paths = list()),
        notify = \(...) invisible(NULL),
        removeNotify = \(...) invisible(NULL),
        logError = \(...) invisible(NULL),
        formatConfidence = 1
    )
    expect_identical(result$status, "success")
    expect_true(result$artifactStageResult$ok)
    expect_true(result$artifactStageResult$committed)
    expect_identical(
        built$workflow$data_tbl,
        payload$assayList
    )
    expect_identical(
        built$workflow$processing_log$setup_import$artifact_stage$run_id,
        result$artifactStageResult$run_id
    )
})

test_that("lipidomics design dual-write hydrates exact S4 and dependencies", {
    root <- withr::local_tempdir(pattern = "lipid040-design-")
    built <- .lipid040Workflow(root)
    payload <- .lipid040FixturePayload("mixed")
    import <- .lipid040ApplyImport(built$workflow, payload)
    expect_true(import$ok)
    object <- .lipid040Design(built$workflow, payload)
    memory_before <- built$workflow$state_manager$exportState()
    output <- persistLipidDesignArtifacts(
        built$workflow,
        log_warn = \(...) invisible(NULL)
    )
    expect_true(output$enabled)
    expect_true(output$ok)
    expect_true(output$committed)
    expect_identical(built$workflow$state_manager$exportState(), memory_before)
    expect_identical(
        built$workflow$state_manager$getState("lipid_raw_data_s4"),
        object
    )
    expected_roles <- c(
        sprintf("cleaned_assay_%04d", seq_along(payload$assayList)),
        "design_matrix", "contrasts", "args", "column_mapping",
        "assay_alignment", "raw_s4_dependencies"
    )
    expect_setequal(names(output$refs), expected_roles)
    expect_identical(
        .lipid040ReadTable(
            built$workflow$workflow_context,
            output$refs$design_matrix
        ),
        built$workflow$design_matrix
    )
    dependencies <- .lipid040ReadTable(
        built$workflow$workflow_context,
        output$refs$raw_s4_dependencies
    )
    expect_identical(dependencies$assay_name, names(payload$assayList))
    expect_identical(
        dependencies$parent_import_run_id,
        rep(import$run_id, nrow(dependencies))
    )
    alignment <- .lipid040ReadTable(
        built$workflow$workflow_context,
        output$refs$assay_alignment
    )
    expect_identical(alignment$assay_name, names(payload$assayList))
    expect_identical(alignment$assay_order, seq_along(payload$assayList))
    expect_identical(
        .lipid040ReadTable(
            built$workflow$workflow_context,
            output$refs$args
        ),
        artifactStageMetadataTable(built$workflow$config_list)
    )
    expect_identical(
        .lipid040ReadTable(
            built$workflow$workflow_context,
            output$refs$column_mapping
        ),
        artifactStageMetadataTable(built$workflow$column_mapping)
    )
    session <- .lipid040Registry(built$workflow$workflow_context)
    parameters <- projectRegistryQuery(
        session,
        "parameters",
        filters = list(
            workflow_id = built$workflow$workflow_context$getIdentity()$workflow_id,
            run_id = output$run_id
        )
    )
    expect_true(all(c(
        "formula_string", "contrasts_kind", "assay_order",
        "parent_import_run_id", "parent_import_generation_id"
    ) %in% parameters$parameter_name))
    closeProjectRegistry(session)
    prepared <- prepareLipidArtifactContext(built$workflow)
    manager <- newLipidArtifactStateManager(prepared)
    expect_identical(manager$getState("lipid_raw_data_s4"), object)
    expect_true(manager$close())
})

test_that("alignment and malformed inputs fail before artifact state advances", {
    payload <- .lipid040FixturePayload("mixed")
    samples <- payload$columnMapping$sample_columns
    groups <- ifelse(grepl("^KO", samples), "KO", "WT")
    design <- data.frame(
        Run = samples,
        group = groups,
        replicates = ave(seq_along(samples), groups, FUN = seq_along),
        stringsAsFactors = FALSE
    )
    reordered <- .lipid040CleanAssays(payload)
    metadata <- setdiff(names(reordered[[2L]]), samples)
    reordered[[2L]] <- reordered[[2L]][c(metadata, rev(samples))]
    expect_true(validateLipidDesignAlignment(
        design,
        reordered,
        payload$columnMapping
    )$valid)
    missing <- reordered
    missing[[2L]][[samples[[1L]]]] <- NULL
    expect_false(validateLipidDesignAlignment(
        design,
        missing,
        payload$columnMapping
    )$valid)
    extra <- reordered
    extra[[2L]]$EXTRA_1 <- seq_len(nrow(extra[[2L]]))
    expect_false(validateLipidDesignAlignment(
        design,
        extra,
        payload$columnMapping
    )$valid)

    root <- withr::local_tempdir(pattern = "lipid040-malformed-")
    built <- .lipid040Workflow(root)
    .lipid040ApplyMemory(built$workflow, payload)
    invalid_source <- payload
    bad_path <- file.path(root, "assay.rds")
    saveRDS(payload$assayList[[1L]], bad_path)
    invalid_source$assaySourceFiles[[1L]] <- bad_path
    output <- persistLipidImportArtifacts(
        built$workflow,
        invalid_source,
        log_warn = \(...) invisible(NULL)
    )
    expect_false(output$ok)
    expect_true(
        "multischolar_unverified_lipidomics_source_encoding" %in%
            output$error_class
    )
    expect_identical(built$workflow$data_tbl, payload$assayList)
})

.lipid040CorruptingManager <- function(
    workflow_context,
    dehydrate_fn,
    validate_bundle_fn,
    hydrate_fn,
    descriptor_contract
) {
    ArtifactWorkflowState$new(
        workflow_context = workflow_context,
        dehydrate_fn = dehydrate_fn,
        validate_bundle_fn = validate_bundle_fn,
        hydrate_fn = \(bundle) {
            value <- hydrate_fn(bundle)
            value@args$changed_after_hydration <- TRUE
            value
        },
        descriptor_contract = descriptor_contract
    )
}

test_that("store, registry, and hydration failures preserve memory state", {
    for (failure_stage in c("before_write", "before_registry_commit")) {
        root <- withr::local_tempdir(pattern = "lipid040-failure-")
        built <- .lipid040Workflow(root)
        payload <- .lipid040FixturePayload("lc")
            .lipid040ApplyMemory(built$workflow, payload)
        memory_before <- built$workflow$data_tbl
        output <- persistLipidImportArtifacts(
            built$workflow,
            payload,
            failure_injector = \(stage, context) {
                if (identical(stage, failure_stage)) {
                    rlang::abort("injected lipidomics import failure")
                }
                invisible(context)
            },
            log_warn = \(...) invisible(NULL)
        )
        expect_true(output$enabled)
        expect_false(output$ok)
        expect_identical(built$workflow$data_tbl, memory_before)
    }

    root <- withr::local_tempdir(pattern = "lipid040-hydration-failure-")
    built <- .lipid040Workflow(root)
    payload <- .lipid040FixturePayload("mixed")
    expect_true(.lipid040ApplyImport(built$workflow, payload)$ok)
    object <- .lipid040Design(built$workflow, payload)
    memory_before <- built$workflow$state_manager$exportState()
    output <- persistLipidDesignArtifacts(
        built$workflow,
        manager_factory = .lipid040CorruptingManager,
        log_warn = \(...) invisible(NULL)
    )
    expect_false(output$ok)
    expect_identical(built$workflow$state_manager$exportState(), memory_before)
    expect_identical(
        built$workflow$state_manager$getState("lipid_raw_data_s4"),
        object
    )
})

test_that("memory, unsupported, legacy, and cross-omic paths stay isolated", {
    memory_root <- withr::local_tempdir(pattern = "lipid040-memory-")
    memory <- .lipid040Workflow(memory_root, backend = "memory")
    payload <- .lipid040FixturePayload("gc")
    output <- .lipid040ApplyImport(memory$workflow, payload)
    expect_false(output$enabled)
    expect_true(output$ok)
    expect_identical(memory$workflow$data_tbl, payload$assayList)
    memory_object <- .lipid040Design(memory$workflow, payload)
    memory_design <- persistLipidDesignArtifacts(
        memory$workflow,
        log_warn = \(...) invisible(NULL)
    )
    expect_false(memory_design$enabled)
    expect_identical(
        memory$workflow$state_manager$getState("lipid_raw_data_s4"),
        memory_object
    )
    expect_false(dir.exists(file.path(memory_root, "state")))

    unsupported_root <- withr::local_tempdir(pattern = "lipid040-msdial-")
    unsupported <- .lipid040Workflow(unsupported_root)
    unsupported$workflow$data_format <- "msdial"
    unsupported$workflow$data_type <- "lipid"
    prepared <- prepareLipidArtifactContext(unsupported$workflow)
    expect_false(prepared$enabled)
    expect_false(unsupported$workflow$workflow_context$isBound())
    expect_false(dir.exists(file.path(unsupported_root, "state")))

    cross_root <- withr::local_tempdir(pattern = "lipid040-cross-")
    cross <- .lipid040Workflow(cross_root, omic_type = "metabolomics")
    cross$workflow$data_format <- "lipidsearch"
    cross$workflow$data_type <- "lipid"
    prepared <- prepareLipidArtifactContext(cross$workflow)
    expect_false(prepared$enabled)
    expect_identical(prepared$reason, "not_lipidomics_context")
    expect_false(cross$workflow$workflow_context$isBound())

    custom_root <- withr::local_tempdir(pattern = "lipid040-custom-")
    custom <- .lipid040Workflow(custom_root)
    custom$workflow$data_format <- "custom"
    custom$workflow$data_type <- "lipid"
    prepared <- prepareLipidArtifactContext(custom$workflow)
    expect_false(prepared$enabled)
    expect_identical(prepared$reason, "not_supported_lipidomics_tuple")
    expect_false(custom$workflow$workflow_context$isBound())
})

test_that("artifact mode uses coordinator state and ignores package globals", {
    root <- withr::local_tempdir(pattern = "lipid040-globals-")
    built <- .lipid040Workflow(root)
    payload <- .lipid040FixturePayload("lc")
    expect_true(.lipid040ApplyImport(built$workflow, payload)$ok)
    object <- .lipid040Design(built$workflow, payload)
    global_names <- c("config_list", "contrasts_tbl")
    had_global <- vapply(
        global_names,
        exists,
        logical(1),
        envir = .GlobalEnv,
        inherits = FALSE
    )
    previous <- mget(
        global_names[had_global],
        envir = .GlobalEnv,
        inherits = FALSE
    )
    withr::defer({
        present <- global_names[vapply(
            global_names,
            exists,
            logical(1),
            envir = .GlobalEnv,
            inherits = FALSE
        )]
        if (length(present)) rm(list = present, envir = .GlobalEnv)
        list2env(previous, envir = .GlobalEnv)
    })
    poison <- list(
        config_list = list(poison = "config"),
        contrasts_tbl = list(poison = "contrasts")
    )
    list2env(poison, envir = .GlobalEnv)
    builder_results <- list(
        design_matrix = built$workflow$design_matrix,
        data_cln = built$workflow$data_cln,
        contrasts_tbl = built$workflow$contrasts_tbl,
        config_list = built$workflow$config_list
    )
    expect_silent(hydrateLipidDesignBuilderWorkflowState(
        results = builder_results,
        workflowData = built$workflow,
        assignFn = \(...) stop("artifact mode must not assign package globals")
    ))
    output <- persistLipidDesignArtifacts(
        built$workflow,
        log_warn = \(...) invisible(NULL)
    )
    expect_true(output$ok)
    expect_identical(
        built$workflow$state_manager$getState("lipid_raw_data_s4"),
        object
    )
    expect_identical(get("config_list", envir = .GlobalEnv), poison$config_list)
    expect_identical(
        get("contrasts_tbl", envir = .GlobalEnv),
        poison$contrasts_tbl
    )
    expect_identical(
        .lipid040ReadTable(
            built$workflow$workflow_context,
            output$refs$contrasts
        ),
        built$workflow$contrasts_tbl
    )
    expect_silent(hydrateLipidDesignImportWorkflowState(
        workflowData = built$workflow,
        importedConfig = built$workflow$config_list,
        importedDesign = built$workflow$design_matrix,
        assayList = built$workflow$data_cln,
        importedContrasts = built$workflow$contrasts_tbl,
        assignFn = \(...) stop("artifact mode must not assign package globals")
    ))
    expect_true(lipidArtifactCoordinatorOwned(built$workflow))
})

.lipid040FrozenWorkload <- function() {
    environment <- new.env(parent = .GlobalEnv)
    sys.source(
        .lipid040RepoPath("tools", "profiling", "omics_workload_contract.R"),
        envir = environment
    )
    adapter_path <- .lipid040RepoPath(
        "tools", "profiling", "omics_workload_lipidomics.R"
    )
    contract_path <- .lipid040RepoPath(
        "tests", "testdata", "omics-parity", "lipidomics", "workloads",
        "mixed-public-ci-v1.json"
    )
    contract <- environment$omicsWorkloadReadContract(contract_path)
    adapter <- environment$omicsWorkloadLoadAdapter(adapter_path, contract)
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    run_dir <- tempfile("lipid040-workload-")
    dir.create(run_dir, recursive = TRUE)
    prepared <- adapter$prepare(list(
        contract = contract,
        run_dir = run_dir
    ))
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
    expect_identical(
        payload_digest,
        contract$expected_digests$payload_sha256
    )
    expect_identical(truth_digest, contract$expected_digests$truth_sha256)
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
        payload = list(
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
            designBatches = unlist(truth$batch_assignments, use.names = FALSE)
        ),
        evidence = evidence
    )
}

test_that("generated custom workload remains memory-only and path-free", {
    frozen <- .lipid040FrozenWorkload()
    root <- withr::local_tempdir(pattern = "lipid040-frozen-")
    built <- .lipid040Workflow(root)
    import <- .lipid040ApplyImport(built$workflow, frozen$payload)
    expect_true(import$ok)
    expect_false(import$enabled)
    expect_identical(import$reason, "not_supported_lipidomics_tuple")
    expect_identical(names(built$workflow$data_tbl), c(
        "LCMS_Pos", "LCMS_Neg", "GCMS"
    ))
    object <- .lipid040Design(
        built$workflow,
        frozen$payload,
        frozen$evidence
    )
    design <- persistLipidDesignArtifacts(
        built$workflow,
        log_warn = \(...) invisible(NULL)
    )
    expect_true(design$ok)
    expect_false(design$enabled)
    expect_identical(
        built$workflow$state_manager$getState("lipid_raw_data_s4"),
        object
    )
    expect_identical(
        object@args$artifact_evidence,
        frozen$evidence
    )
    retained <- paste(capture.output(str(list(import, design))), collapse = "\n")
    expect_false(grepl("private_path|source_path", retained))
    expect_false(grepl("[.](fasta|fa|crdownload)", retained, ignore.case = TRUE))
})

test_that("artifact helpers collate before callers and preserve package defaults", {
    description <- read.dcf(.lipid040RepoPath("DESCRIPTION"))
    collate <- strsplit(description[1L, "Collate"], "[[:space:]]+")[[1L]]
    collate <- gsub("^'|'$", "", collate[nzchar(collate)])
    helpers <- c(
        "utils_artifact_lipid_descriptors.R",
        "mod_lipid_import_artifact_helpers.R",
        "mod_lipid_design_artifact_helpers.R"
    )
    expect_true(all(vapply(
        helpers,
        \(helper) sum(collate == helper) == 1L,
        logical(1)
    )))
    expect_lt(
        match("mod_lipid_import_artifact_helpers.R", collate),
        match("mod_lipid_import_processing_helpers.R", collate)
    )
    expect_lt(
        match("mod_lipid_design_artifact_helpers.R", collate),
        match("mod_lipid_design_builder_helpers.R", collate)
    )
    expect_lt(
        match("mod_lipid_design_artifact_helpers.R", collate),
        match("mod_lipid_design_import_helpers.R", collate)
    )
    expect_identical(
        formals(runLipidImportProcessing)$persistArtifactFn,
        quote(persistLipidImportArtifacts)
    )
    expect_identical(
        formals(registerLipidDesignBuilderResultsObserver)$persistArtifactFn,
        quote(persistLipidDesignArtifacts)
    )
    expect_identical(
        formals(runLipidDesignImportConfirmationShell)$persistArtifactFn,
        quote(persistLipidDesignArtifacts)
    )
})
