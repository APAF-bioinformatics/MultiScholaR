if (!methods::isClass("SyntheticOmicData")) {
    methods::setClass(
        "SyntheticOmicData",
        slots = c(assay = "data.frame", design = "data.frame", metadata = "list"),
        validity = function(object) {
            assay_columns <- c(
                "feature", "sample", "abundance", "quality", "score"
            )
            design_columns <- c("sample", "group")
            if (!identical(names(object@assay), assay_columns)) {
                return("synthetic assay columns are invalid")
            }
            if (!identical(names(object@design), design_columns)) {
                return("synthetic design columns are invalid")
            }
            if (!all(unique(object@assay$sample) %in% object@design$sample)) {
                return("synthetic assay contains an unknown sample")
            }
            TRUE
        }
    )
}

syntheticArtifactCodecId <- function() {
    "multischolar.s4.synthetic_omic_data.v1"
}

syntheticArtifactQueryId <- function() {
    "synthetic.assay.lookup.v1"
}

syntheticArtifactDescriptorId <- function() {
    "syntheticomics.synthetic.feature.v1"
}

syntheticArtifactIdentity <- function() {
    workflowCapabilityIdentity(
        omic_type = "syntheticomics",
        workflow_id = "syntheticomics.gui",
        workflow_type = "SYNTHETIC",
        workflow_slug = "synthetic_standard",
        input_format = "synthetic",
        data_level = "feature",
        acquisition_mode = "synthetic"
    )
}

syntheticArtifactCodec <- function(version = 1L) {
    list(
        codec_id = syntheticArtifactCodecId(),
        codec_version = as.integer(version),
        class_name = "SyntheticOmicData",
        payload_schema_id = "multischolar.rectangular",
        payload_schema_version = 1L
    )
}

syntheticArtifactDescriptor <- function(
    certification = "evict",
    auto_eligible = TRUE,
    descriptor_version = "1.0.0"
) {
    owner <- syntheticArtifactDescriptorId()
    codec_id <- syntheticArtifactCodecId()
    query_id <- syntheticArtifactQueryId()
    newArtifactWorkflowDescriptor(
        descriptor_id = owner,
        descriptor_version = descriptor_version,
        identity = syntheticArtifactIdentity(),
        stages = stats::setNames(list(list(
            stage_id = "import",
            state_roles = c("payload_0001", "payload_0002"),
            codec_ids = codec_id,
            query_operation_ids = query_id,
            maximum_rollout = "evict"
        )), "import"),
        codecs = stats::setNames(list(syntheticArtifactCodec()), codec_id),
        queries = stats::setNames(list(list(
            operation_id = query_id,
            state_role = "payload_0001",
            projections = c(
                "feature", "sample", "abundance", "quality", "score"
            ),
            filters = list(
                feature = list(
                    column = "feature",
                    type = "character",
                    operators = c("equal", "in")
                ),
                abundance = list(
                    column = "abundance",
                    type = "double",
                    operators = c("gte", "lte", "between", "is_missing")
                ),
                score = list(
                    column = "score",
                    type = "double",
                    operators = "is_missing"
                )
            ),
            order_by = c("feature", "sample"),
            max_rows = 100L,
            max_bytes = 1024L * 1024L
        )), query_id),
        fixtures = list(owner_id = owner, fixture_ids = "synthetic.fixture.v1"),
        scientific_oracle = list(
            owner_id = owner,
            oracle_id = "synthetic.oracle.identity.v1",
            oracle_version = "1.0.0",
            tolerances = c(abundance = 0)
        ),
        compatibility_products = list(
            owner_id = owner,
            product_ids = c("synthetic.s4", "synthetic.design")
        ),
        evidence = list(
            owner_id = owner,
            inventory_ids = "synthetic.inventory.v1",
            codec_ids = codec_id,
            stage_ids = "import",
            lifecycle_ids = "synthetic.lifecycle.v1",
            performance_thresholds = c(
                max_rss_bytes = 256 * 1024^2,
                max_query_p95_seconds = 0.5
            )
        ),
        migration = list(
            owner_id = owner,
            strategy_id = "synthetic.memory_to_artifact.v1",
            from_backend = "memory",
            to_backend = "artifact"
        ),
        rollback = list(
            owner_id = owner,
            strategy_id = "synthetic.force_memory.v1",
            target_backend = "memory"
        ),
        certification = list(
            owner_id = owner,
            status = certification,
            auto_eligible = isTRUE(auto_eligible)
        )
    )
}

syntheticArtifactObject <- function(offset = 0) {
    methods::new(
        "SyntheticOmicData",
        assay = data.frame(
            feature = c("B", "A", "A", "C"),
            sample = c("S2", "S1", "S2", "S1"),
            abundance = c(4, 1, 2, 3) + offset,
            quality = factor(
                c("high", "low", "high", "low"),
                levels = c("low", "high")
            ),
            score = c(NA_real_, NaN, Inf, -Inf),
            stringsAsFactors = FALSE
        ),
        design = data.frame(
            sample = c("S1", "S2"),
            group = factor(c("control", "case"), levels = c("case", "control")),
            stringsAsFactors = FALSE
        ),
        metadata = list(source = "synthetic", offset = offset)
    )
}

syntheticArtifactCatalogues <- function(descriptor = syntheticArtifactDescriptor()) {
    list(
        descriptors = newArtifactDescriptorCatalogue(list(descriptor)),
        codecs = newArtifactCodecCatalogue(list(syntheticArtifactCodec()))
    )
}

syntheticArtifactContext <- function(
    project_root,
    descriptor_catalogue,
    backend = "artifact",
    rollout = "evict"
) {
    context <- createWorkflowContext(
        experiment_paths = list(
            base_dir = project_root,
            omic_label = "synthetic_study"
        ),
        omic_type = "syntheticomics",
        storage_policy = list(
            requested_backend = backend,
            requested_rollout = rollout,
            migration_requested = TRUE,
            project_id = "synthetic-project"
        )
    )
    bindWorkflowContextFromImport(
        context,
        workflow_type = "SYNTHETIC",
        input_format = "synthetic",
        data_level = "feature",
        descriptor_catalogue = descriptor_catalogue
    )
    context
}

syntheticArtifactSkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}

syntheticArtifactAssayRef <- function(manager, store) {
    references <- manager$getStateMetadata()$artifact_refs
    matches <- vapply(references, function(ref) {
        sidecar <- artifactStoreReadSidecar(
            store,
            artifactStoreManagedPaths(store, ref$logical_key, ref$artifact_id)$sidecar,
            validate_payload = FALSE
        )
        identical(sidecar$codec_metadata$owner, "SyntheticOmicData@assay")
    }, logical(1))
    stopifnot(sum(matches) == 1L)
    references[[which(matches)]]
}

expectSyntheticArtifactExact <- function(before, after) {
    expect_identical(class(after), class(before))
    for (slot_name in methods::slotNames(before)) {
        expect_identical(
            methods::slot(after, slot_name),
            methods::slot(before, slot_name),
            info = slot_name
        )
    }
    expect_identical(methods::validObject(after, test = TRUE), TRUE)
}

test_that("artifact extension catalogues are complete, owned, and immutable", {
    descriptor <- syntheticArtifactDescriptor()
    catalogues <- syntheticArtifactCatalogues(descriptor)

    expect_identical(
        names(descriptor),
        c(
            "schema", "schema_version", "descriptor_id", "descriptor_version",
            "identity", "stages", "codecs", "queries", "fixtures",
            "scientific_oracle", "compatibility_products", "evidence", "migration",
            "rollback", "certification", "descriptor_digest"
        )
    )
    expect_true(environmentIsLocked(catalogues$descriptors))
    expect_true(environmentIsLocked(catalogues$codecs))
    expect_error(
        catalogues$descriptors$descriptors <- list(),
        regexp = "locked binding"
    )
    expect_error(
        catalogues$codecs$codecs <- list(),
        regexp = "locked binding"
    )
    copy <- artifactDescriptorCatalogueValues(catalogues$descriptors)
    copy[[1L]]$certification$status <- "uncertified"
    expect_identical(
        findArtifactWorkflowDescriptor(
            syntheticArtifactIdentity(), catalogues$descriptors
        )$certification$status,
        "evict"
    )

    expect_error(
        newArtifactDescriptorCatalogue(list(descriptor, descriptor)),
        class = "multischolar_duplicate_artifact_descriptor"
    )
    expect_error(
        newArtifactCodecCatalogue(list(syntheticArtifactCodec(), syntheticArtifactCodec())),
        class = "multischolar_duplicate_artifact_codec"
    )
    future_codec <- syntheticArtifactCodec(version = 2L)
    expect_error(
        newArtifactCodecCatalogue(list(future_codec)),
        class = "multischolar_unsupported_artifact_codec_version"
    )
    injected <- descriptor
    injected$migration$runtime <- function() NULL
    expect_error(
        validateArtifactWorkflowDescriptor(injected),
        class = "multischolar_artifact_descriptor_code_injection"
    )
    attributed <- descriptor
    attr(attributed$evidence$performance_thresholds, "runtime") <- function() NULL
    expect_error(
        validateArtifactWorkflowDescriptor(attributed),
        class = "multischolar_artifact_descriptor_code_injection"
    )
    future <- descriptor
    future$schema_version <- ARTIFACT_DESCRIPTOR_SCHEMA_VERSION + 1L
    expect_error(
        validateArtifactWorkflowDescriptor(future),
        class = "multischolar_future_artifact_descriptor_version"
    )
    expect_s3_class(
        newWorkflowState(workflow_descriptor = future),
        "WorkflowState"
    )
    wrong_owner <- descriptor
    wrong_owner$scientific_oracle$owner_id <- "another.omic.v1"
    expect_error(
        validateArtifactWorkflowDescriptor(wrong_owner),
        class = "multischolar_artifact_descriptor_evidence_mismatch"
    )
    incompatible <- descriptor
    incompatible$codecs[[syntheticArtifactCodecId()]]$codec_version <- 2L
    incompatible$descriptor_digest <- artifactDescriptorDigest(incompatible)
    expect_error(
        artifactCodecAdapter(incompatible, catalogues$codecs),
        class = "multischolar_incompatible_artifact_codec_catalogue"
    )
})

test_that("unknown and uncertified extension tuples fail before project resources", {
    syntheticArtifactSkipDependencies()
    empty_catalogue <- newArtifactDescriptorCatalogue()
    auto_root <- withr::local_tempdir()
    auto_context <- syntheticArtifactContext(
        auto_root,
        empty_catalogue,
        backend = "auto",
        rollout = "dual_write"
    )
    expect_identical(
        auto_context$getStorageDecision()$effective_backend,
        "memory"
    )
    expect_false(dir.exists(file.path(auto_root, "state")))

    forced_root <- withr::local_tempdir()
    expect_error(
        syntheticArtifactContext(
            forced_root,
            empty_catalogue,
            backend = "artifact",
            rollout = "dual_write"
        ),
        class = "multischolar_unknown_workflow_capability"
    )
    expect_false(dir.exists(file.path(forced_root, "state")))

    uncertified <- syntheticArtifactDescriptor(
        certification = "uncertified",
        auto_eligible = FALSE
    )
    uncertified_catalogue <- newArtifactDescriptorCatalogue(list(uncertified))
    uncertified_auto_root <- withr::local_tempdir()
    uncertified_auto <- syntheticArtifactContext(
        uncertified_auto_root,
        uncertified_catalogue,
        backend = "auto",
        rollout = "dual_write"
    )
    expect_identical(
        uncertified_auto$getStorageDecision()$effective_backend,
        "memory"
    )
    expect_false(dir.exists(file.path(uncertified_auto_root, "state")))
    uncertified_root <- withr::local_tempdir()
    expect_error(
        syntheticArtifactContext(
            uncertified_root,
            uncertified_catalogue,
            backend = "artifact",
            rollout = "dual_write"
        ),
        class = "multischolar_artifact_not_certified"
    )
    expect_false(dir.exists(file.path(uncertified_root, "state")))

    for (rollout in c("dual_write", "read_through", "evict")) {
        context <- syntheticArtifactContext(
            withr::local_tempdir(),
            syntheticArtifactCatalogues()$descriptors,
            backend = "artifact",
            rollout = rollout
        )
        expect_identical(context$getStorageDecision()$effective_rollout, rollout)
    }
    exact_context <- syntheticArtifactContext(
        withr::local_tempdir(),
        syntheticArtifactCatalogues()$descriptors,
        backend = "artifact",
        rollout = "dual_write"
    )
    expect_error(
        newWorkflowState(
            workflow_context = exact_context,
            descriptor_catalogue = empty_catalogue
        ),
        class = "multischolar_artifact_state_descriptor_mismatch"
    )
})

test_that("synthetic omic round trips, reopens, and rolls back through shared state", {
    syntheticArtifactSkipDependencies()
    project_root <- withr::local_tempdir()
    descriptor <- syntheticArtifactDescriptor()
    catalogues <- syntheticArtifactCatalogues(descriptor)
    context <- syntheticArtifactContext(project_root, catalogues$descriptors)
    manager <- newWorkflowState(
        workflow_context = context,
        workflow_descriptor = descriptor,
        descriptor_catalogue = catalogues$descriptors,
        codec_catalogue = catalogues$codecs
    )
    object <- syntheticArtifactObject()
    expect_identical(
        manager$saveState(
            "imported",
            object,
            list(stage = "import"),
            "Synthetic import"
        ),
        "imported"
    )
    expect_lte(manager$getCacheInfo()$entries, 1L)
    expectSyntheticArtifactExact(object, manager$getState())
    expect_identical(manager$getCacheInfo()$entries, 1L)
    expect_true(manager$close())

    pin <- jsonlite::read_json(
        file.path(
            project_root,
            context$getPaths()$relative_paths$workflow_state_root,
            "artifact-descriptor.json"
        ),
        simplifyVector = TRUE
    )
    expect_identical(pin$contract$descriptor_id, descriptor$descriptor_id)
    expect_identical(pin$contract$descriptor_version, descriptor$descriptor_version)
    expect_identical(pin$contract$descriptor_digest, descriptor$descriptor_digest)

    skip_if(Sys.which("Rscript") == "")
    descriptor_path <- tempfile(fileext = ".rds")
    script_path <- tempfile(fileext = ".R")
    output_path <- tempfile(fileext = ".json")
    log_path <- tempfile(fileext = ".log")
    repo_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    saveRDS(descriptor, descriptor_path, version = 3L)
    writeLines(c(
        "args <- commandArgs(trailingOnly = TRUE)",
        "devtools::load_all(args[[1L]], quiet = TRUE)",
        "methods::setClass(",
        "    'SyntheticOmicData',",
        "    slots = c(assay = 'data.frame', design = 'data.frame', metadata = 'list')",
        ")",
        "descriptor <- readRDS(args[[3L]])",
        "descriptors <- newArtifactDescriptorCatalogue(list(descriptor))",
        "codecs <- newArtifactCodecCatalogue(descriptor$codecs)",
        "context <- createWorkflowContext(",
        "    list(base_dir = args[[2L]], omic_label = 'synthetic_study'),",
        "    'syntheticomics',",
        "    storage_policy = list(",
        "        requested_backend = 'artifact',",
        "        requested_rollout = 'evict',",
        "        migration_requested = TRUE,",
        "        project_id = 'synthetic-project'",
        "    )",
        ")",
        "bindWorkflowContextFromImport(",
        "    context, 'SYNTHETIC', 'synthetic', 'feature',",
        "    descriptor_catalogue = descriptors",
        ")",
        "manager <- newWorkflowState(",
        "    workflow_context = context,",
        "    workflow_descriptor = descriptor,",
        "    descriptor_catalogue = descriptors,",
        "    codec_catalogue = codecs",
        ")",
        "object <- manager$getState()",
        "result <- list(",
        "    state = manager$getCurrentStateName(),",
        "    class = class(object)[[1L]],",
        "    valid = identical(methods::validObject(object, test = TRUE), TRUE),",
        "    rows = nrow(object@assay),",
        "    cache_entries = manager$getCacheInfo()$entries",
        ")",
        "jsonlite::write_json(result, args[[4L]], auto_unbox = TRUE)",
        "manager$close()"
    ), script_path)
    status <- system2(
        Sys.which("Rscript"),
        c(
            "--vanilla", script_path, repo_root, project_root,
            descriptor_path, output_path
        ),
        stdout = log_path,
        stderr = log_path
    )
    if (!identical(as.integer(status), 0L)) {
        testthat::fail(paste(readLines(log_path, warn = FALSE), collapse = "\n"))
    }
    fresh <- jsonlite::read_json(output_path, simplifyVector = TRUE)
    expect_identical(fresh$state, "imported")
    expect_identical(fresh$class, "SyntheticOmicData")
    expect_true(fresh$valid)
    expect_identical(fresh$rows, 4L)
    expect_identical(fresh$cache_entries, 1L)

    reopened_context <- syntheticArtifactContext(project_root, catalogues$descriptors)
    reopened <- newWorkflowState(
        workflow_context = reopened_context,
        workflow_descriptor = descriptor,
        descriptor_catalogue = catalogues$descriptors,
        codec_catalogue = catalogues$codecs
    )
    withr::defer(reopened$close())
    expect_identical(reopened$getCurrentStateName(), "imported")
    expect_identical(reopened$getCacheInfo()$entries, 0L)
    expectSyntheticArtifactExact(object, reopened$getState())
    expect_identical(reopened$getCacheInfo()$entries, 1L)

    expect_error(
        newWorkflowState(workflow_context = reopened_context),
        class = "multischolar_artifact_state_descriptor_required"
    )
    changed <- syntheticArtifactDescriptor(descriptor_version = "1.0.1")
    changed_catalogues <- syntheticArtifactCatalogues(changed)
    changed_context <- syntheticArtifactContext(
        project_root, changed_catalogues$descriptors
    )
    expect_error(
        newWorkflowState(
            workflow_context = changed_context,
            workflow_descriptor = changed,
            descriptor_catalogue = changed_catalogues$descriptors,
            codec_catalogue = changed_catalogues$codecs
        ),
        class = "multischolar_artifact_state_descriptor_pin_mismatch"
    )

    memory_context <- syntheticArtifactContext(
        project_root,
        catalogues$descriptors,
        backend = "memory",
        rollout = "dual_write"
    )
    memory <- newWorkflowState(
        workflow_context = memory_context,
        workflow_descriptor = descriptor,
        descriptor_catalogue = catalogues$descriptors,
        codec_catalogue = catalogues$codecs
    )
    expect_s3_class(memory, "WorkflowState")
    expect_identical(memory$getCurrentStateName(), "initial")
    expect_true(length(list.files(project_root, pattern = "[.]parquet$", recursive = TRUE)) > 0L)
})

test_that("synthetic omics can opt into generic parent-backed row selections", {
    syntheticArtifactSkipDependencies()
    project_root <- withr::local_tempdir()
    descriptor <- syntheticArtifactDescriptor()
    catalogues <- syntheticArtifactCatalogues(descriptor)
    context <- syntheticArtifactContext(project_root, catalogues$descriptors)
    manager <- newWorkflowState(
        workflow_context = context,
        workflow_descriptor = descriptor,
        descriptor_catalogue = catalogues$descriptors,
        codec_catalogue = catalogues$codecs
    )
    withr::defer(manager$close())
    parent <- syntheticArtifactObject()
    manager$saveState(
        "imported",
        parent,
        list(stage = "import"),
        "Synthetic import"
    )
    child <- parent
    child@assay <- child@assay[c(4L, 2L), , drop = FALSE]
    key_columns <- c("feature", "sample")
    parent_keys <- artifactWorkflowStateEntityKeys(parent@assay, key_columns)
    child_keys <- artifactWorkflowStateEntityKeys(child@assay, key_columns)
    rejected_keys <- parent_keys[!parent_keys %in% child_keys]
    hint <- newArtifactRowSelectionHint(
        slot_name = "assay",
        key_columns = key_columns,
        method = "synthetic_quality_filter",
        normalized_parameters = list(minimum_quality = "low"),
        software = list(
            name = "SyntheticExtension",
            version = "1.0.0",
            source = "testthat"
        ),
        lineage = list(
            audit_record_id = "synthetic-audit-1",
            state_name = "quality_filtered",
            parent_state = "imported",
            parent_record_id = "synthetic-audit-0"
        ),
        rejected_reasons = stats::setNames(
            rep("below_quality_threshold", length(rejected_keys)),
            rejected_keys
        )
    )
    result <- manager$commitState(
        "quality_filtered",
        child,
        list(stage = "quality_filter"),
        "Synthetic quality filter",
        persistence_hint = hint,
        expected_parent_generation_id = manager$getCurrentGenerationId()
    )

    expect_identical(result$status, "accepted")
    expect_identical(result$representation, "row_selection")
    expectSyntheticArtifactExact(child, manager$getState())
    expect_lte(manager$getCacheInfo()$entries, 1L)

    changed <- child
    changed@assay$abundance <- changed@assay$abundance + 1
    materialized <- manager$commitState(
        "normalized",
        changed,
        list(stage = "normalization"),
        "Synthetic normalization",
        expected_parent_generation_id = manager$getCurrentGenerationId()
    )
    expect_identical(materialized$representation, "materialized")
    expectSyntheticArtifactExact(changed, manager$getState())
})

test_that("descriptor queries are typed, deterministic, and bounded before collection", {
    syntheticArtifactSkipDependencies()
    project_root <- withr::local_tempdir()
    descriptor <- syntheticArtifactDescriptor()
    catalogues <- syntheticArtifactCatalogues(descriptor)
    context <- syntheticArtifactContext(project_root, catalogues$descriptors)
    manager <- newWorkflowState(
        workflow_context = context,
        workflow_descriptor = descriptor,
        descriptor_catalogue = catalogues$descriptors,
        codec_catalogue = catalogues$codecs
    )
    withr::defer(manager$close())
    manager$saveState(
        "imported",
        syntheticArtifactObject(),
        list(stage = "import"),
        "Synthetic import"
    )
    store <- newArtifactStore(context$getPaths(), "synthetic-project")
    ref <- syntheticArtifactAssayRef(manager, store)
    query_session <- newArtifactQuerySession(store)
    withr::defer(query_session$close())
    filters <- list(
        feature = list(operator = "in", value = c("A", "B")),
        abundance = list(operator = "gte", value = 2)
    )
    first <- queryArtifactRef(
        store,
        ref,
        descriptor,
        syntheticArtifactQueryId(),
        projections = c("feature", "sample", "abundance"),
        filters = filters,
        limit = 10L,
        query_session = query_session
    )
    second <- queryArtifactRef(
        store,
        ref,
        descriptor,
        syntheticArtifactQueryId(),
        projections = c("feature", "sample", "abundance"),
        filters = filters,
        limit = 10L,
        query_session = query_session
    )
    expect_identical(first, second)
    expect_identical(query_session$getInfo()$borrow_count, 2L)
    expect_identical(first$feature, c("A", "B"))
    expect_identical(first$sample, c("S2", "S2"))
    expect_identical(first$abundance, c(2, 4))
    special <- queryArtifactRef(
        store,
        ref,
        descriptor,
        syntheticArtifactQueryId(),
        projections = c("feature", "quality", "score"),
        limit = 4L
    )
    expect_s3_class(special$quality, "factor")
    expect_identical(levels(special$quality), c("low", "high"))
    expect_identical(
        special$score,
        c(NaN, Inf, NA_real_, -Inf)
    )
    missing <- queryArtifactRef(
        store,
        ref,
        descriptor,
        syntheticArtifactQueryId(),
        projections = c("feature", "score"),
        filters = list(score = list(operator = "is_missing", value = TRUE)),
        limit = 4L
    )
    expect_identical(missing$feature, "B")
    expect_true(is.na(missing$score) && !is.nan(missing$score))
    expect_error(
        queryArtifactRef(
            store, ref, descriptor, syntheticArtifactQueryId(),
            projections = "not_owned"
        ),
        class = "multischolar_invalid_artifact_query_projection"
    )
    expect_error(
        queryArtifactRef(
            store, ref, descriptor, syntheticArtifactQueryId(),
            filters = list(sql = list(operator = "equal", value = "DROP TABLE"))
        ),
        class = "multischolar_invalid_artifact_query_filter"
    )
    expect_error(
        queryArtifactRef(
            store, ref, descriptor, syntheticArtifactQueryId(), limit = 101L
        ),
        class = "multischolar_artifact_query_row_limit_exceeded"
    )
    changed <- syntheticArtifactDescriptor(descriptor_version = "1.0.1")
    expect_error(
        queryArtifactRef(
            store, ref, changed, syntheticArtifactQueryId(), limit = 4L
        ),
        class = "multischolar_artifact_query_descriptor_pin_mismatch"
    )
    expect_error(
        queryArtifactRef(
            store,
            ref,
            descriptor,
            syntheticArtifactQueryId(),
            limit = 4L,
            resource_policy = list(max_result_bytes = 1L)
        ),
        class = "multischolar_artifact_query_byte_limit_exceeded"
    )
})

test_that("stage adapters expose only validated ordinary scientific inputs", {
    syntheticArtifactSkipDependencies()
    object <- syntheticArtifactObject()
    legacy_calls <- 0L
    legacy <- function() {
        legacy_calls <<- legacy_calls + 1L
        syntheticArtifactObject(100)
    }
    scientific <- function(value, multiplier) {
        expect_s4_class(value, "SyntheticOmicData")
        expect_false(any(vapply(
            list(value),
            function(item) inherits(item, c("DBIConnection", "MultiScholaRArtifactRef")),
            logical(1)
        )))
        sum(value@assay$abundance) * multiplier
    }
    result <- runArtifactStageAdapter(
        scientific,
        multiplier = 2,
        explicit_input = object,
        legacy_provider = legacy,
        allow_legacy = TRUE
    )
    expect_identical(result, 20)
    expect_identical(legacy_calls, 0L)
    expect_error(
        runArtifactStageAdapter(
            scientific,
            multiplier = 1,
            legacy_provider = legacy,
            allow_legacy = FALSE
        ),
        class = "multischolar_missing_artifact_stage_input"
    )
    expect_identical(legacy_calls, 0L)
    expect_identical(
        runArtifactStageAdapter(
            scientific,
            multiplier = 1,
            legacy_provider = legacy,
            allow_legacy = TRUE
        ),
        410
    )
    expect_identical(legacy_calls, 1L)
})

test_that("only immutable descriptor-backed tuples are artifact certified", {
    production <- mergeWorkflowDescriptorCapabilities()
    descriptors <- artifactDescriptorCatalogueValues(
        artifactWorkflowDescriptorCatalogue()
    )
    descriptor_ids <- vapply(
        descriptors,
        `[[`,
        character(1),
        "descriptor_id"
    )
    expect_setequal(descriptor_ids, c(
        "proteomics.diann.peptide.dia.v1",
        "proteomics.maxquant.protein.lfq.v1",
        "proteomics.fragpipe.protein.lfq.v1",
        "proteomics.pd_tmt.protein.tmt.v1",
        "metabolomics.custom.metabolite.standard.v1",
        "lipidomics.lipidsearch.lipid.standard.v1"
    ))
    expect_false(any(grepl("spectronaut", descriptor_ids, fixed = TRUE)))
    eligible <- Filter(\(capability) capability$artifact_eligible, production)
    eligible_ids <- vapply(eligible, `[[`, character(1), "capability_id")
    expect_setequal(unname(eligible_ids), unname(descriptor_ids))
    expect_true(all(vapply(eligible, \(capability) {
        match(
            capability$maximum_artifact_rollout,
            .WORKFLOW_ARTIFACT_ROLLOUTS
        ) >= match("dual_write", .WORKFLOW_ARTIFACT_ROLLOUTS) &&
            !isTRUE(capability$auto_eligible)
    }, logical(1))))
    unsupported <- Filter(\(capability) !capability$artifact_eligible, production)
    expect_true(all(vapply(unsupported, function(capability) {
        !isTRUE(capability$auto_eligible) &&
            is.null(capability$maximum_artifact_rollout)
    }, logical(1))))
    sources <- vapply(
        c(
            "R/utils_artifact_store.R",
            "R/utils_project_duckdb.R",
            "R/utils_artifact_queries.R",
            "R/utils_workflow_state_artifact_backend.R"
        ),
        function(path) paste(readLines(
            testthat::test_path("..", "..", path),
            warn = FALSE
        ), collapse = "\n"),
        character(1)
    )
    expect_false(any(grepl(
        "syntheticomics|proteomics|metabolomics|lipidomics",
        sources
    )))
})
