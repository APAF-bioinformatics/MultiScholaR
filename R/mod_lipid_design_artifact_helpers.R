# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

lipidArtifactImportParent <- function(workflow_data) {
    result <- workflow_data$artifact_stage_results$import
    assay_refs <- if (is.list(result)) result$assay_refs else NULL
    valid <- is.list(result) && isTRUE(result$enabled) && isTRUE(result$ok) &&
        isTRUE(result$committed) && workflowCapabilityScalarString(result$run_id) &&
        workflowCapabilityScalarString(result$generation_id) &&
        is.list(assay_refs) && length(assay_refs) > 0L
    if (!isTRUE(valid)) {
        return(list(
            run_id = NULL,
            generation_id = NULL,
            assay_refs = list(),
            assay_manifest = NULL
        ))
    }
    list(
        run_id = result$run_id,
        generation_id = result$generation_id,
        assay_refs = assay_refs,
        assay_manifest = result$assay_manifest
    )
}

lipidArtifactMemoryCheckpoint <- function(
    workflow_data,
    state_name = "lipid_raw_data_s4"
) {
    manager <- workflow_data$state_manager
    state_object <- tryCatch(manager$getState(state_name), error = \(error) NULL)
    valid <- inherits(manager, "WorkflowState") &&
        methods::is(state_object, "LipidomicsAssayData") &&
        identical(methods::validObject(state_object, test = TRUE), TRUE)
    if (!isTRUE(valid)) {
        lipidArtifactAbort(
            "lipidomics design requires its completed memory S4 checkpoint",
            "multischolar_missing_lipidomics_memory_checkpoint"
        )
    }
    exact <- identical(
        methods::slot(state_object, "lipid_data"),
        workflow_data$data_cln
    ) && identical(
        methods::slot(state_object, "design_matrix"),
        workflow_data$design_matrix
    ) && identical(
        methods::slot(state_object, "args"),
        workflow_data$config_list
    )
    if (!isTRUE(exact)) {
        lipidArtifactAbort(
            "lipidomics S4 checkpoint differs from coordinator design state",
            "multischolar_inexact_lipidomics_design_checkpoint"
        )
    }
    formula_string <- workflow_data$config_list$deAnalysisParameters$formula_string
    if (!workflowCapabilityScalarString(formula_string)) {
        formula_string <- "~ 0 + group"
    }
    preflight <- validateLipidDesignDaPreflight(
        designMatrix = workflow_data$design_matrix,
        assayList = workflow_data$data_cln,
        contrastsTbl = workflow_data$contrasts_tbl,
        formulaString = formula_string,
        columnMapping = workflow_data$column_mapping,
        requireContrasts = FALSE
    )
    if (!isTRUE(preflight$valid)) {
        lipidArtifactAbort(
            sprintf(
                "lipidomics artifact design preflight failed: %s",
                paste(preflight$errors, collapse = "; ")
            ),
            "multischolar_invalid_lipidomics_artifact_design"
        )
    }
    list(manager = manager, state_object = state_object)
}

lipidArtifactContrastsTable <- function(contrasts) {
    if (is.null(contrasts)) {
        return(data.frame(
            contrasts = character(),
            stringsAsFactors = FALSE
        ))
    }
    if (is.character(contrasts)) {
        return(data.frame(contrasts = contrasts, stringsAsFactors = FALSE))
    }
    artifactStagePortableTable(contrasts, "contrasts", lipidArtifactAbort)
}

lipidArtifactAssayAlignment <- function(workflow_data, state_object) {
    assays <- methods::slot(state_object, "lipid_data")
    mapping <- workflow_data$column_mapping
    sample_id <- methods::slot(state_object, "sample_id")
    design_samples <- as.character(
        methods::slot(state_object, "design_matrix")[[sample_id]]
    )
    rows <- lapply(seq_along(assays), \(index) {
        assay <- assays[[index]]
        sample_columns <- intersect(design_samples, names(assay))
        feature_column <- methods::slot(state_object, "lipid_id_column")
        feature_ids <- assay[[feature_column]]
        data.frame(
            assay_name = names(assays)[[index]],
            assay_order = as.integer(index),
            cleaned_artifact_role = sprintf("cleaned_assay_%04d", index),
            feature_count = as.integer(nrow(assay)),
            sample_count = as.integer(length(sample_columns)),
            lipid_id_column = feature_column,
            lipid_class_column = methods::slot(
                state_object,
                "annotation_id_column"
            ),
            feature_order_digest = artifactSemanticDigest(feature_ids),
            sample_order_digest = artifactSemanticDigest(sample_columns),
            design_sample_order_digest = artifactSemanticDigest(design_samples),
            sample_columns_json = artifactStageParameterJson(sample_columns),
            mapping_digest = artifactSemanticDigest(mapping),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

lipidArtifactRawDependencies <- function(parent, state_object) {
    assay_names <- names(methods::slot(state_object, "lipid_data"))
    refs <- parent$assay_refs[assay_names]
    rows <- lapply(seq_along(assay_names), \(index) {
        ref <- refs[[index]]
        data.frame(
            assay_name = assay_names[[index]],
            assay_order = as.integer(index),
            parent_import_run_id = parent$run_id %||% NA_character_,
            parent_import_generation_id = parent$generation_id %||% NA_character_,
            parent_artifact_id = if (is.list(ref)) {
                ref$artifact_id
            } else {
                NA_character_
            },
            parent_semantic_digest = if (is.list(ref)) {
                ref$hash_policy$semantic$digest
            } else {
                NA_character_
            },
            s4_class = class(state_object)[[1L]],
            s4_codec_id = "multischolar.s4.lipidomics_assay_data",
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, rows)
}

lipidArtifactDesignTables <- function(workflow_data, state_object, parent) {
    alignment <- lipidArtifactAssayAlignment(workflow_data, state_object)
    assays <- methods::slot(state_object, "lipid_data")
    assay_tables <- assays[alignment$assay_name]
    names(assay_tables) <- alignment$cleaned_artifact_role
    c(assay_tables, list(
        design_matrix = methods::slot(state_object, "design_matrix"),
        contrasts = lipidArtifactContrastsTable(workflow_data$contrasts_tbl),
        args = artifactStageMetadataTable(workflow_data$config_list),
        column_mapping = artifactStageMetadataTable(
            workflow_data$column_mapping
        ),
        assay_alignment = alignment,
        raw_s4_dependencies = lipidArtifactRawDependencies(parent, state_object)
    ))
}

lipidArtifactWriteDesignStage <- function(
    prepared,
    workflow_data,
    checkpoint,
    failure_injector
) {
    parent <- lipidArtifactImportParent(workflow_data)
    if (!workflowCapabilityScalarString(parent$run_id) ||
        length(parent$assay_refs) == 0L) {
        lipidArtifactAbort(
            "lipidomics design artifact has no committed import parent",
            "multischolar_missing_lipidomics_import_parent"
        )
    }
    formula_string <- workflow_data$config_list$deAnalysisParameters$formula_string
    if (!workflowCapabilityScalarString(formula_string)) {
        formula_string <- "~ 0 + group"
    }
    stage <- writeArtifactStage(
        prepared$context,
        prepared$descriptor,
        stage_id = "design",
        tables = lipidArtifactDesignTables(
            workflow_data,
            checkpoint$state_object,
            parent
        ),
        parameters = list(
            capability_id = prepared$descriptor$descriptor_id,
            state_name = "lipid_raw_data_s4",
            workflow_type = workflowStateType(checkpoint$manager),
            formula_string = formula_string,
            contrasts_kind = artifactStageContrastsKind(
                workflow_data$contrasts_tbl
            ),
            assay_order = names(methods::slot(
                checkpoint$state_object,
                "lipid_data"
            )),
            parent_import_run_id = parent$run_id,
            parent_import_generation_id = parent$generation_id,
            readthrough_contract_version = 1L
        ),
        deferred_commit = TRUE,
        failure_injector = failure_injector,
        abort_fn = lipidArtifactAbort
    )
    stage$parent_import <- parent
    stage$readthrough_contract_version <- 1L
    stage
}

lipidArtifactStateAudit <- function(stage) {
    list(
        record_id = artifactOpaqueId("audit"),
        provenance_status = "artifact_dual_write_validated",
        stage_id = stage$stage_id,
        run_id = stage$run_id,
        action_id = stage$action_id,
        readthrough_contract_version = stage$readthrough_contract_version,
        parent_import_run_id = stage$parent_import$run_id,
        parent_import_generation_id = stage$parent_import$generation_id,
        stage_artifact_refs = stage$refs
    )
}

newLipidArtifactStateManager <- function(
    prepared,
    manager_factory = ArtifactWorkflowState$new
) {
    manager_factory(
        workflow_context = prepared$context,
        dehydrate_fn = dehydrateLipidomicsS4Artifact,
        validate_bundle_fn = validateLipidomicsS4Bundle,
        hydrate_fn = hydrateLipidomicsS4Artifact,
        descriptor_contract = artifactStageDescriptorContract(
            prepared$descriptor
        )
    )
}

lipidArtifactSaveDesignState <- function(
    manager,
    workflow_data,
    state_object,
    stage,
    failure_injector = NULL
) {
    manager$setWorkflowType("lipidomics_standard")
    manager$saveState(
        state_name = "lipid_raw_data_s4",
        s4_data_object = state_object,
        config_object = workflow_data$config_list,
        description = "Lipidomics design dual-write checkpoint.",
        audit_metadata = lipidArtifactStateAudit(stage),
        failure_injector = failure_injector
    )
    hydrated <- manager$getState("lipid_raw_data_s4")
    if (!identical(hydrated, state_object) ||
        !identical(methods::validObject(hydrated, test = TRUE), TRUE)) {
        lipidArtifactAbort(
            "lipidomics design artifact differs from memory checkpoint",
            "multischolar_inexact_lipidomics_artifact_hydration"
        )
    }
    list(
        state_manifest = manager$exportState(),
        state_metadata = manager$getStateMetadata("lipid_raw_data_s4")
    )
}

lipidArtifactFailDesignStage <- function(context, stage, state_error) {
    try(
        artifactStageUpdateStatus(
            context,
            stage,
            completed = FALSE,
            abort_fn = lipidArtifactAbort
        ),
        silent = TRUE
    )
    stop(state_error)
}

persistLipidDesignArtifactsCore <- function(
    workflow_data,
    failure_injector = NULL,
    manager_factory = ArtifactWorkflowState$new
) {
    prepared <- prepareLipidArtifactContext(workflow_data)
    if (!isTRUE(prepared$enabled)) {
        return(list(
            enabled = FALSE,
            ok = TRUE,
            stage_id = "design",
            reason = prepared$reason
        ))
    }
    checkpoint <- lipidArtifactMemoryCheckpoint(workflow_data)
    stage <- lipidArtifactWriteDesignStage(
        prepared,
        workflow_data,
        checkpoint,
        failure_injector
    )
    manager <- NULL
    state_result <- tryCatch({
        manager <- newLipidArtifactStateManager(prepared, manager_factory)
        lipidArtifactSaveDesignState(
            manager,
            workflow_data,
            checkpoint$state_object,
            stage,
            failure_injector
        )
    }, error = \(error) error)
    if (inherits(state_result, "error")) {
        if (!is.null(manager)) try(manager$close(), silent = TRUE)
        lipidArtifactFailDesignStage(prepared$context, stage, state_result)
    }
    manager$close()
    manager <- NULL
    artifactStageUpdateStatus(
        prepared$context,
        stage,
        completed = TRUE,
        failure_injector = failure_injector,
        abort_fn = lipidArtifactAbort
    )
    stage$committed <- TRUE
    c(list(enabled = TRUE, ok = TRUE), stage, state_result)
}

persistLipidDesignArtifacts <- function(
    workflow_data,
    failure_injector = NULL,
    manager_factory = ArtifactWorkflowState$new,
    log_warn = logger::log_warn
) {
    runArtifactStageSafely(
        workflow_data,
        "design",
        \() persistLipidDesignArtifactsCore(
            workflow_data,
            failure_injector,
            manager_factory
        ),
        "lipidomics",
        log_warn
    )
}
