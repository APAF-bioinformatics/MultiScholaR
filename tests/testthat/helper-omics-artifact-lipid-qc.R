lipid042SkipDependencies <- function() {
    for (package in c("DBI", "duckdb", "filelock", "arrow")) {
        testthat::skip_if_not_installed(package)
    }
}
lipid042Capability <- function() {
    capabilities <- workflowCapabilityCatalogue()
    capability <- capabilities[[
        "lipidomics.lipidsearch.lipid.standard.v1"
    ]]
    capability$artifact_eligible <- TRUE
    capability$auto_eligible <- FALSE
    capability$maximum_artifact_rollout <- "dual_write"
    capability$explicit_maximum_artifact_rollout <- "dual_write"
    capability
}

lipid042Context <- function(root, project_id = "omics-art-042") {
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
        capabilities = list(lipid042Capability())
    )
    context
}

lipid042Manager <- function(context) {
    descriptor <- artifactLipidomicsWorkflowDescriptor()
    manager <- ArtifactWorkflowState$new(
        workflow_context = context,
        dehydrate_fn = dehydrateLipidomicsS4Artifact,
        validate_bundle_fn = validateLipidomicsS4Bundle,
        hydrate_fn = hydrateLipidomicsS4Artifact,
        descriptor_contract = artifactStageDescriptorContract(descriptor)
    )
    manager$setWorkflowType("lipidomics_standard")
    manager
}

lipid042Object <- function(
    include_duplicates = FALSE,
    include_itsd = TRUE,
    layout = "combined"
) {
    module_ci_lipid_qc_object(
        layout = layout,
        include_duplicates = include_duplicates,
        include_itsd = include_itsd
    )
}

lipid042Workflow <- function(manager, context) {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$workflow_context <- context
    workflow$config_list <- list(
        globalParameters = list(workflow_type = "lipidomics_standard"),
        qc = list(mode = "artifact-ticket")
    )
    workflow$artifact_stage_results <- list()
    workflow$processing_log <- list()
    workflow$tab_status <- list(
        setup_import = "complete",
        design_matrix = "complete",
        quality_control = "pending",
        normalization = "disabled"
    )
    workflow$data_format <- "lipidsearch"
    workflow$data_type <- "lipid"
    workflow
}

lipid042SaveInitial <- function(manager, object) {
    manager$saveState(
        "lipid_raw_data_s4",
        object,
        list(stage = "design"),
        "lipidomics QC parent"
    )
    invisible(object)
}

lipid042Manifest <- function(context, manager, state_name = NULL) {
    if (is.null(state_name)) state_name <- manager$getCurrentStateName()
    row <- manager$states[[state_name]]
    store <- newArtifactStore(
        context$getPaths(),
        context$getIdentity()$project_id
    )
    artifactWorkflowStateReadManifest(store, row$manifest_relative_path)
}

lipid042Metadata <- function(context, manager, state_name = NULL) {
    manifest <- lipid042Manifest(context, manager, state_name)
    artifactWorkflowStateUnserializeMetadata(
        manifest$data$metadata_json,
        "OMICS-ART-042 metadata"
    )
}

lipid042ReadRef <- function(context, ref) {
    store <- newArtifactStore(
        context$getPaths(),
        context$getIdentity()$project_id
    )
    sidecar <- artifactStoreReadSidecar(
        store,
        artifactStoreManagedPaths(
            store,
            ref$logical_key,
            ref$artifact_id
        )$sidecar,
        validate_payload = TRUE
    )
    payload <- arrow::read_parquet(
        artifactStoreResolveFile(store, ref$relative_path, must_exist = TRUE),
        as_data_frame = FALSE
    )
    decodeArtifactRectangular(payload, sidecar$codec_metadata)
}

lipid042Filter <- function(object) {
    filtered <- object
    filtered@lipid_data <- lapply(
        object@lipid_data,
        function(assay) {
            value <- assay[c(4L, 1L, 3L), , drop = FALSE]
            rownames(value) <- NULL
            value
        }
    )
    methods::validObject(filtered)
    filtered
}

lipid042ResolveDuplicates <- function(object) {
    resolved <- lapply(object@lipid_data, \(assay) {
        sample_columns <- names(assay)[vapply(assay, is.numeric, logical(1))]
        resolveLipidDuplicateFeaturesByIntensity(
            assay,
            object@lipid_id_column,
            sample_columns
        )
    })
    stats <- Map(\(before, after) {
        list(
            original = nrow(before),
            resolved = nrow(after),
            removed = nrow(before) - nrow(after)
        )
    }, object@lipid_data, resolved)
    list(resolvedAssayList = resolved, statsList = stats)
}

lipid042ItsdAnalysis <- function(object, pattern, require_matches = FALSE) {
    assay_names <- names(object@lipid_data)
    sample_columns <- as.character(
        object@design_matrix[[object@sample_id]]
    )
    metrics <- lapply(object@lipid_data, \(assay) {
        getLipidInternalStandardMetrics(
            assay,
            pattern,
            object@lipid_id_column,
            object@sample_id,
            sample_columns
        )
    })
    total <- sum(vapply(metrics, nrow, integer(1)))
    if (isTRUE(require_matches) && total == 0L) {
        stop("No internal standards found", call. = FALSE)
    }
    matched_columns <- stats::setNames(
        rep(list(object@lipid_id_column), length(assay_names)),
        assay_names
    )
    searched_columns <- stats::setNames(
        rep(list(unique(c(
            object@lipid_id_column,
            object@annotation_id_column
        ))), length(assay_names)),
        assay_names
    )
    combined <- do.call(rbind, Map(\(value, assay_name) {
        value$assay <- rep(assay_name, nrow(value))
        value
    }, metrics, assay_names))
    list(
        metrics = combined,
        metricsByAssay = metrics,
        pattern = pattern,
        nIsTotal = total,
        matchedColumns = matched_columns,
        searchedColumns = searched_columns
    )
}

lipid042PlotBuildData <- function(object, plot_fn, ...) {
    plots <- suppressWarnings(suppressMessages(plot_fn(object, ...)))
    lapply(plots, function(plot) ggplot2::ggplot_build(plot)$data)
}

lipid042ReviewedObject <- function() {
    repo <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
    fixture_root <- file.path(
        repo,
        "tests",
        "testdata",
        "e2e",
        "lipid_canonical"
    )
    fixture_names <- c(
        LCMS_Pos = "lipidsearch_lcms_pos.txt",
        LCMS_Neg = "lipidsearch_lcms_neg.txt",
        GCMS = "lipidsearch_gcms.txt"
    )
    assays <- lapply(fixture_names, \(filename) {
        utils::read.delim(
            file.path(fixture_root, filename),
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
    })
    samples <- grep("^(WT|KO)_", names(assays[[1L]]), value = TRUE)
    design <- data.frame(
        Run = samples,
        group = sub("_.*$", "", samples),
        batch = rep(c("B1", "B2"), length.out = length(samples)),
        tech_rep_group = samples,
        stringsAsFactors = FALSE
    )
    createLipidomicsAssayData(
        assays,
        design,
        lipid_id_column = "LipidName",
        annotation_id_column = "LipidClass",
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        database_identifier_type = "lipid_name",
        internal_standard_regex = "^IS_",
        args = list(evidence_class = "independently_reviewed_fixture")
    )
}
