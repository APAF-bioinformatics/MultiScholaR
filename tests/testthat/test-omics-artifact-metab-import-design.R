.metab031RepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

.metab031Paths <- function(root, project_id, omic_type = "metabolomics") {
    paths <- list(
        base_dir = root,
        project_id = project_id,
        omic_type = omic_type,
        omic_label = "metabolomics-study",
        source_dir = file.path(root, "source"),
        results_dir = file.path(root, "results")
    )
    dir.create(paths$source_dir, recursive = TRUE)
    dir.create(paths$results_dir, recursive = TRUE)
    paths
}

.metab031Workflow <- function(
    root,
    backend = "artifact",
    omic_type = "metabolomics",
    project_id = paste0("metab031-", basename(root))
) {
    paths <- .metab031Paths(root, project_id, omic_type)
    workflow <- new.env(parent = emptyenv())
    workflow$workflow_context <- createWorkflowContext(
        paths,
        omic_type,
        "metabolomics-study",
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

.metab031ReadTable <- function(context, ref) {
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

.metab031FixturePayload <- function(kind = c("lc", "gc", "mixed")) {
    kind <- match.arg(kind)
    specifications <- switch(kind,
        lc = c(
            LCMS_Pos = "tests/testdata/e2e/metab_lc/lcms_pos_features.tsv",
            LCMS_Neg = "tests/testdata/e2e/metab_lc/lcms_neg_features.tsv"
        ),
        gc = c(
            GCMS = "tests/testdata/e2e/metab_gc/gcms_features.tsv"
        ),
        mixed = c(
            GCMS = "tests/testdata/e2e/metab_combined/combined_gcms_features.tsv",
            LCMS_Pos = paste0(
                "tests/testdata/e2e/metab_combined/",
                "combined_lcms_features.tsv"
            )
        )
    )
    source_paths <- vapply(
        specifications,
        .metab031RepoPath,
        character(1)
    )
    assays <- lapply(source_paths, \(path) {
        utils::read.delim(
            path,
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
    })
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
        processingLog = list(
            evidence_class = "independent_reviewed_fixture",
            fixture_kind = kind
        ),
        assaySourceFiles = as.list(source_paths),
        sourceFiles = as.list(source_paths)
    )
}

.metab031ApplyImport <- function(workflow, payload) {
    applyMetabImportWorkflowPayload(
        workflow,
        payload,
        logInfo = \(...) invisible(NULL)
    )
    persistMetabImportArtifacts(
        workflow,
        payload,
        log_warn = \(...) invisible(NULL)
    )
}

.metab031CleanAssays <- function(payload) {
    mapping <- payload$columnMapping
    metadata <- unique(c(
        mapping$metabolite_id_col,
        mapping$annotation_col
    ))
    metadata <- metadata[!is.na(metadata) & nzchar(metadata)]
    lapply(payload$assayList, \(assay) {
        assay[c(metadata, mapping$sample_columns)]
    })
}

.metab031Design <- function(workflow, payload, evidence = NULL) {
    samples <- payload$columnMapping$sample_columns
    groups <- payload$designGroups
    if (is.null(groups)) groups <- ifelse(grepl("^KO", samples), "KO", "WT")
    batches <- payload$designBatches
    if (is.null(batches)) batches <- rep(NA_character_, length(samples))
    replicates <- ave(seq_along(samples), groups, FUN = seq_along)
    workflow$data_cln <- .metab031CleanAssays(payload)
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
        globalParameters = list(workflow_type = "metabolomics_standard"),
        deAnalysisParameters = list(formula_string = "~ 0 + group"),
        artifact_evidence = evidence
    )
    object <- createMetaboliteAssayData(
        metabolite_data = workflow$data_cln,
        design_matrix = workflow$design_matrix,
        metabolite_id_column = payload$columnMapping$metabolite_id_col,
        annotation_id_column = payload$columnMapping$annotation_col,
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        database_identifier_type = "FixtureFeatureName",
        internal_standard_regex = NA_character_,
        args = workflow$config_list
    )
    workflow$state_manager$setWorkflowType("metabolomics_standard")
    workflow$state_manager$saveState(
        "metab_raw_data_s4",
        object,
        workflow$config_list,
        "metabolomics design memory checkpoint"
    )
    object
}

.metab031Registry <- function(context) {
    registry <- projectRegistryForContext(context)
    openProjectRegistryReadOnly(registry)
}

test_that("only custom metabolomics receives a dual-write descriptor", {
    descriptor <- artifactMetabolomicsWorkflowDescriptor()
    expect_identical(
        descriptor$descriptor_id,
        "metabolomics.custom.metabolite.standard.v1"
    )
    expect_identical(descriptor$certification$status, "dual_write")
    expect_false(descriptor$certification$auto_eligible)
    expect_identical(
        names(descriptor$codecs),
        "multischolar.s4.metabolite_assay_data"
    )
    expect_false(any(vapply(
        artifactDescriptorCatalogueValues(artifactWorkflowDescriptorCatalogue()),
        \(candidate) identical(candidate$identity$input_format, "msdial") &&
            identical(candidate$identity$omic_type, "metabolomics"),
        logical(1)
    )))
    support <- workflowFormatSupportCatalogue()
    custom <- support[[which(vapply(
        support,
        \(entry) identical(entry$format_id, "metabolomics.custom"),
        logical(1)
    ))]]
    msdial <- support[[which(vapply(
        support,
        \(entry) identical(entry$format_id, "metabolomics.msdial"),
        logical(1)
    ))]]
    expect_identical(custom$support_status, "scientifically_supported")
    expect_identical(msdial$support_status, "reader_characterized")
})

test_that("reviewed custom assays dual-write exact independent provenance", {
    for (kind in c("lc", "gc", "mixed")) {
        root <- withr::local_tempdir(pattern = paste0("metab031-", kind, "-"))
        built <- .metab031Workflow(root)
        payload <- .metab031FixturePayload(kind)
        memory_before <- unserialize(serialize(payload$assayList, NULL))
        output <- .metab031ApplyImport(built$workflow, payload)
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
            "feature_count", "sample_count", "feature_id_column",
            "annotation_id_column", "sample_columns_json", "column_schema_json"
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
                .metab031ReadTable(
                    built$workflow$workflow_context,
                    output$assay_refs[[assay_name]]
                ),
                memory_before[[assay_name]],
                info = paste(kind, assay_name)
            )
        }
        session <- .metab031Registry(built$workflow$workflow_context)
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
        expect_identical(sources$format_id, "custom.column_mapped_tabular")
        closeProjectRegistry(session)
    }
})

test_that("metabolomics design dual-write hydrates exact S4 and dependencies", {
    root <- withr::local_tempdir(pattern = "metab031-design-")
    built <- .metab031Workflow(root)
    payload <- .metab031FixturePayload("mixed")
    import <- .metab031ApplyImport(built$workflow, payload)
    expect_true(import$ok)
    object <- .metab031Design(built$workflow, payload)
    memory_before <- built$workflow$state_manager$exportState()
    output <- persistMetabDesignArtifacts(
        built$workflow,
        log_warn = \(...) invisible(NULL)
    )
    expect_true(output$enabled)
    expect_true(output$ok)
    expect_true(output$committed)
    expect_identical(built$workflow$state_manager$exportState(), memory_before)
    expect_identical(
        built$workflow$state_manager$getState("metab_raw_data_s4"),
        object
    )
    expected_roles <- c(
        sprintf("cleaned_assay_%04d", seq_along(payload$assayList)),
        "design_matrix", "contrasts", "args", "column_mapping",
        "assay_alignment", "raw_s4_dependencies"
    )
    expect_setequal(names(output$refs), expected_roles)
    expect_identical(
        .metab031ReadTable(
            built$workflow$workflow_context,
            output$refs$design_matrix
        ),
        built$workflow$design_matrix
    )
    dependencies <- .metab031ReadTable(
        built$workflow$workflow_context,
        output$refs$raw_s4_dependencies
    )
    expect_identical(dependencies$assay_name, names(payload$assayList))
    expect_identical(
        dependencies$parent_import_run_id,
        rep(import$run_id, nrow(dependencies))
    )
    alignment <- .metab031ReadTable(
        built$workflow$workflow_context,
        output$refs$assay_alignment
    )
    expect_identical(alignment$assay_name, names(payload$assayList))
    expect_identical(alignment$assay_order, seq_along(payload$assayList))
    expect_identical(
        .metab031ReadTable(
            built$workflow$workflow_context,
            output$refs$args
        ),
        artifactStageMetadataTable(built$workflow$config_list)
    )
    expect_identical(
        .metab031ReadTable(
            built$workflow$workflow_context,
            output$refs$column_mapping
        ),
        artifactStageMetadataTable(built$workflow$column_mapping)
    )
    session <- .metab031Registry(built$workflow$workflow_context)
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
    prepared <- prepareMetabArtifactContext(built$workflow)
    manager <- newMetabArtifactStateManager(prepared)
    expect_identical(manager$getState("metab_raw_data_s4"), object)
    expect_true(manager$close())
})

test_that("alignment and malformed inputs fail before artifact state advances", {
    payload <- .metab031FixturePayload("mixed")
    samples <- payload$columnMapping$sample_columns
    groups <- ifelse(grepl("^KO", samples), "KO", "WT")
    design <- data.frame(
        Run = samples,
        group = groups,
        replicates = ave(seq_along(samples), groups, FUN = seq_along),
        stringsAsFactors = FALSE
    )
    reordered <- .metab031CleanAssays(payload)
    metadata <- setdiff(names(reordered[[2L]]), samples)
    reordered[[2L]] <- reordered[[2L]][c(metadata, rev(samples))]
    expect_true(validateMetabDesignAlignment(
        design,
        reordered,
        payload$columnMapping
    )$valid)
    missing <- reordered
    missing[[2L]][[samples[[1L]]]] <- NULL
    expect_false(validateMetabDesignAlignment(
        design,
        missing,
        payload$columnMapping
    )$valid)
    extra <- reordered
    extra[[2L]]$EXTRA_1 <- seq_len(nrow(extra[[2L]]))
    expect_false(validateMetabDesignAlignment(
        design,
        extra,
        payload$columnMapping
    )$valid)

    root <- withr::local_tempdir(pattern = "metab031-malformed-")
    built <- .metab031Workflow(root)
    applyMetabImportWorkflowPayload(
        built$workflow,
        payload,
        logInfo = \(...) invisible(NULL)
    )
    invalid_source <- payload
    bad_path <- file.path(root, "assay.rds")
    saveRDS(payload$assayList[[1L]], bad_path)
    invalid_source$assaySourceFiles[[1L]] <- bad_path
    output <- persistMetabImportArtifacts(
        built$workflow,
        invalid_source,
        log_warn = \(...) invisible(NULL)
    )
    expect_false(output$ok)
    expect_true(
        "multischolar_unverified_metabolomics_source_encoding" %in%
            output$error_class
    )
    expect_identical(built$workflow$data_tbl, payload$assayList)
})

.metab031CorruptingManager <- function(
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
        root <- withr::local_tempdir(pattern = "metab031-failure-")
        built <- .metab031Workflow(root)
        payload <- .metab031FixturePayload("lc")
        applyMetabImportWorkflowPayload(
            built$workflow,
            payload,
            logInfo = \(...) invisible(NULL)
        )
        memory_before <- built$workflow$data_tbl
        output <- persistMetabImportArtifacts(
            built$workflow,
            payload,
            failure_injector = \(stage, context) {
                if (identical(stage, failure_stage)) {
                    rlang::abort("injected metabolomics import failure")
                }
                invisible(context)
            },
            log_warn = \(...) invisible(NULL)
        )
        expect_true(output$enabled)
        expect_false(output$ok)
        expect_identical(built$workflow$data_tbl, memory_before)
    }

    root <- withr::local_tempdir(pattern = "metab031-hydration-failure-")
    built <- .metab031Workflow(root)
    payload <- .metab031FixturePayload("mixed")
    expect_true(.metab031ApplyImport(built$workflow, payload)$ok)
    object <- .metab031Design(built$workflow, payload)
    memory_before <- built$workflow$state_manager$exportState()
    output <- persistMetabDesignArtifacts(
        built$workflow,
        manager_factory = .metab031CorruptingManager,
        log_warn = \(...) invisible(NULL)
    )
    expect_false(output$ok)
    expect_identical(built$workflow$state_manager$exportState(), memory_before)
    expect_identical(
        built$workflow$state_manager$getState("metab_raw_data_s4"),
        object
    )
})

test_that("memory, unsupported, legacy, and cross-omic paths stay isolated", {
    memory_root <- withr::local_tempdir(pattern = "metab031-memory-")
    memory <- .metab031Workflow(memory_root, backend = "memory")
    payload <- .metab031FixturePayload("gc")
    output <- .metab031ApplyImport(memory$workflow, payload)
    expect_false(output$enabled)
    expect_true(output$ok)
    expect_identical(memory$workflow$data_tbl, payload$assayList)
    memory_object <- .metab031Design(memory$workflow, payload)
    memory_design <- persistMetabDesignArtifacts(
        memory$workflow,
        log_warn = \(...) invisible(NULL)
    )
    expect_false(memory_design$enabled)
    expect_identical(
        memory$workflow$state_manager$getState("metab_raw_data_s4"),
        memory_object
    )
    expect_false(dir.exists(file.path(memory_root, "state")))

    unsupported_root <- withr::local_tempdir(pattern = "metab031-msdial-")
    unsupported <- .metab031Workflow(unsupported_root)
    unsupported$workflow$data_format <- "msdial"
    unsupported$workflow$data_type <- "metabolite"
    prepared <- prepareMetabArtifactContext(unsupported$workflow)
    expect_false(prepared$enabled)
    expect_false(unsupported$workflow$workflow_context$isBound())
    expect_false(dir.exists(file.path(unsupported_root, "state")))

    cross_root <- withr::local_tempdir(pattern = "metab031-cross-")
    cross <- .metab031Workflow(cross_root, omic_type = "lipidomics")
    cross$workflow$data_format <- "custom"
    cross$workflow$data_type <- "metabolite"
    prepared <- prepareMetabArtifactContext(cross$workflow)
    expect_false(prepared$enabled)
    expect_identical(prepared$reason, "not_metabolomics_context")
    expect_false(cross$workflow$workflow_context$isBound())
})

test_that("artifact mode uses coordinator state and ignores package globals", {
    root <- withr::local_tempdir(pattern = "metab031-globals-")
    built <- .metab031Workflow(root)
    payload <- .metab031FixturePayload("lc")
    expect_true(.metab031ApplyImport(built$workflow, payload)$ok)
    object <- .metab031Design(built$workflow, payload)
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
    output <- persistMetabDesignArtifacts(
        built$workflow,
        log_warn = \(...) invisible(NULL)
    )
    expect_true(output$ok)
    expect_identical(
        built$workflow$state_manager$getState("metab_raw_data_s4"),
        object
    )
    expect_identical(get("config_list", envir = .GlobalEnv), poison$config_list)
    expect_identical(
        get("contrasts_tbl", envir = .GlobalEnv),
        poison$contrasts_tbl
    )
    expect_identical(
        .metab031ReadTable(
            built$workflow$workflow_context,
            output$refs$contrasts
        ),
        built$workflow$contrasts_tbl
    )
    expect_silent(hydrateMetabDesignImportWorkflowState(
        workflowData = built$workflow,
        importedDesign = built$workflow$design_matrix,
        assayList = built$workflow$data_cln,
        importedContrasts = built$workflow$contrasts_tbl,
        assignFn = \(...) stop("artifact mode must not assign package globals")
    ))
    expect_true(metabArtifactCoordinatorOwned(built$workflow))
})

.metab031FrozenWorkload <- function() {
    environment <- new.env(parent = .GlobalEnv)
    sys.source(
        .metab031RepoPath("tools", "profiling", "omics_workload_contract.R"),
        envir = environment
    )
    adapter_path <- .metab031RepoPath(
        "tools", "profiling", "omics_workload_metabolomics.R"
    )
    contract_path <- .metab031RepoPath(
        "tests", "testdata", "omics-parity", "metabolomics", "workloads",
        "mixed-public-ci-v1.json"
    )
    contract <- environment$omicsWorkloadReadContract(contract_path)
    adapter <- environment$omicsWorkloadLoadAdapter(adapter_path, contract)
    do.call(RNGkind, as.list(unlist(contract$rng$kind, use.names = FALSE)))
    set.seed(as.integer(contract$rng$seed))
    run_dir <- tempfile("metab031-workload-")
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
            designBatches = unlist(truth$batch_assignments, use.names = FALSE)
        ),
        evidence = evidence
    )
}

test_that("frozen three-assay workload dual-writes without private inputs", {
    frozen <- .metab031FrozenWorkload()
    root <- withr::local_tempdir(pattern = "metab031-frozen-")
    built <- .metab031Workflow(root)
    import <- .metab031ApplyImport(built$workflow, frozen$payload)
    expect_true(import$ok)
    expect_identical(names(import$assay_refs), c(
        "LCMS_Pos", "LCMS_Neg", "GCMS"
    ))
    object <- .metab031Design(
        built$workflow,
        frozen$payload,
        frozen$evidence
    )
    design <- persistMetabDesignArtifacts(
        built$workflow,
        log_warn = \(...) invisible(NULL)
    )
    expect_true(design$ok)
    prepared <- prepareMetabArtifactContext(built$workflow)
    manager <- newMetabArtifactStateManager(prepared)
    expect_identical(manager$getState("metab_raw_data_s4"), object)
    expect_identical(
        manager$getState("metab_raw_data_s4")@args$artifact_evidence,
        frozen$evidence
    )
    expect_true(manager$close())
    retained <- paste(capture.output(str(list(import, design))), collapse = "\n")
    expect_false(grepl("/home/doktersmol", retained, fixed = TRUE))
    expect_false(grepl("cotton_report", retained, fixed = TRUE))
})

test_that("artifact helpers collate before callers and preserve package defaults", {
    description <- read.dcf(.metab031RepoPath("DESCRIPTION"))
    collate <- strsplit(description[1L, "Collate"], "[[:space:]]+")[[1L]]
    collate <- gsub("^'|'$", "", collate[nzchar(collate)])
    helpers <- c(
        "utils_artifact_metab_descriptors.R",
        "mod_metab_import_artifact_helpers.R",
        "mod_metab_design_artifact_helpers.R"
    )
    expect_true(all(vapply(
        helpers,
        \(helper) sum(collate == helper) == 1L,
        logical(1)
    )))
    expect_lt(
        match("mod_metab_import_artifact_helpers.R", collate),
        match("mod_metab_import_processing_helpers.R", collate)
    )
    expect_lt(
        match("mod_metab_design_artifact_helpers.R", collate),
        match("mod_metab_design.R", collate)
    )
    expect_identical(
        formals(runMetabImportProcessing)$persistArtifactFn,
        quote(persistMetabImportArtifacts)
    )
    expect_identical(
        formals(registerMetabDesignBuilderResultsObserver)$persistArtifactFn,
        quote(persistMetabDesignArtifacts)
    )
    expect_identical(
        formals(registerMetabDesignImportObserverShell)$persistArtifactFn,
        quote(persistMetabDesignArtifacts)
    )
})
