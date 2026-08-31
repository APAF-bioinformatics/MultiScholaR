# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_DIA_NORM_AUDIT_SCHEMA <- "multischolar.prot_dia_normalization_audit"
.PROT_DIA_NORM_AUDIT_VERSION <- 1L
.PROT_DIA_NORM_REF_SCHEMA <- "multischolar.prot_dia_normalization_state_ref"
.PROT_DIA_NORM_REF_VERSION <- 1L
.PROT_DIA_NORM_STAGES <- c(
    "normalization",
    "ruv_correction",
    "correlation_filter"
)

protDiaNormArtifactMode <- function(stage_id = NULL) {
    mode <- getOption("multischolar.prot_dia.normalization_artifacts", "enabled")
    if (!is.null(stage_id)) {
        stage_id <- match.arg(stage_id, .PROT_DIA_NORM_STAGES)
        mode <- getOption(
            paste0("multischolar.prot_dia.normalization_artifacts.", stage_id),
            mode
        )
    }
    match.arg(mode, c("enabled", "disabled"))
}

protDiaNormArgs <- function(object) {
    if (!isS4(object) || !"args" %in% methods::slotNames(object)) return(list())
    args <- object@args
    if (is.list(args)) args else list()
}

protDiaNormCurrentRecord <- function(object) {
    audit <- protDiaNormArgs(object)$normalization_artifact_audit
    if (!is.list(audit) || !is.list(audit$records) ||
        length(audit$records) == 0L) {
        return(NULL)
    }
    record <- tail(audit$records, 1L)[[1L]]
    if (!identical(record$record_id, audit$current_record_id)) {
        protDiaArtifactAbort(
            "DIA-NN normalization audit head differs from its latest record",
            "multischolar_invalid_prot_dia_normalization_audit"
        )
    }
    record
}

protDiaNormInheritedRecordId <- function(object, audit_name, fallback) {
    audit <- protDiaNormArgs(object)[[audit_name]]
    record_id <- audit$current_record_id
    if (workflowStateScalarString(record_id)) record_id else fallback
}

protDiaNormSummary <- function(object) {
    if (!methods::is(object, "ProteinQuantitativeData")) {
        return(list(class = class(object)[1L], rows = NA_integer_))
    }
    data <- object@protein_quant_table
    key <- resolveProteinQuantIdentityColumn(object)
    sample_columns <- setdiff(names(data), key)
    list(
        class = class(object)[1L],
        rows = as.integer(nrow(data)),
        columns = as.integer(ncol(data)),
        active_protein_key = key,
        distinct_proteins = as.integer(length(unique(data[[key]]))),
        sample_columns = sample_columns,
        data_digest = .peptideQcDigest(data),
        design_digest = .peptideQcDigest(object@design_matrix),
        protein_id_table_digest = .peptideQcDigest(object@protein_id_table)
    )
}

protDiaNormPackageVersion <- function(package) {
    tryCatch(
        as.character(utils::packageVersion(package)),
        error = \(error) "unavailable"
    )
}

protDiaNormSoftware <- function(stage_id, parameters) {
    engine <- switch(stage_id,
        normalization = {
            method <- parameters$normalization_method
            if (identical(method, "none")) {
                list(name = "MultiScholaR", version = .peptideQcImplementationVersion())
            } else {
                list(name = "limma", version = protDiaNormPackageVersion("limma"))
            }
        },
        ruv_correction = if (identical(parameters$ruv_mode, "skip")) {
            list(name = "MultiScholaR", version = .peptideQcImplementationVersion())
        } else {
            list(name = "RUVIIIC", version = protDiaNormPackageVersion("RUVIIIC"))
        },
        correlation_filter = list(
            name = "MultiScholaR",
            version = .peptideQcImplementationVersion()
        )
    )
    list(
        name = "MultiScholaR",
        version = .peptideQcImplementationVersion(),
        source = "R package",
        stage_id = stage_id,
        scientific_engine = engine
    )
}

protDiaNormBoundedValue <- function(value, owner) {
    if (is.null(value)) return(NULL)
    if (is.data.frame(value) || is.matrix(value)) {
        return(list(
            owner = owner,
            storage = "s4_args_named_payload",
            rows = as.integer(nrow(value)),
            columns = as.integer(ncol(value)),
            digest = .peptideQcDigest(value)
        ))
    }
    if (is.atomic(value) && as.numeric(object.size(value)) > 4096) {
        return(list(
            owner = owner,
            storage = "s4_args_named_payload",
            length = as.integer(length(value)),
            digest = .peptideQcDigest(value)
        ))
    }
    if (is.list(value)) {
        value_names <- names(value)
        bounded <- lapply(seq_along(value), \(index) {
            label <- if (!is.null(value_names) && nzchar(value_names[[index]])) {
                value_names[[index]]
            } else {
                as.character(index)
            }
            protDiaNormBoundedValue(value[[index]], paste0(owner, "[[", label, "]]"))
        })
        names(bounded) <- value_names
        return(bounded)
    }
    value
}

protDiaNormAuditParameters <- function(parameters) {
    bounded <- protDiaNormBoundedValue(
        parameters,
        "ProteinQuantitativeData@args$normalization_artifact_parameters"
    )
    .peptideQcCanonicalise(bounded)
}

protDiaNormAttachParameters <- function(object, parameters, stage_id) {
    if (!methods::is(object, "ProteinQuantitativeData")) return(object)
    stage_parameters <- object@args$normalization_artifact_parameters
    if (!is.list(stage_parameters)) stage_parameters <- list()
    stage_parameters[[stage_id]] <- parameters
    object@args$normalization_artifact_parameters <- stage_parameters
    object
}

protDiaNormEmitAudit <- function(
    before,
    after,
    stage_id,
    parent_state,
    current_state,
    parent_generation_id,
    parameters,
    status,
    transformation_type,
    software,
    now
) {
    if (!methods::is(after, "ProteinQuantitativeData")) return(after)
    prior <- protDiaNormCurrentRecord(before)
    parent_record_id <- if (is.null(prior)) {
        protDiaNormInheritedRecordId(
            before,
            "protein_qc_audit",
            "legacy_untracked_protein_qc"
        )
    } else {
        prior$record_id
    }
    bounded_parameters <- protDiaNormAuditParameters(parameters)
    substantive <- list(
        schema = .PROT_DIA_NORM_AUDIT_SCHEMA,
        schema_version = .PROT_DIA_NORM_AUDIT_VERSION,
        stage_id = stage_id,
        status = status,
        transformation_type = transformation_type,
        parent_state = parent_state,
        current_state = current_state,
        parent_generation_id = parent_generation_id,
        parent_record_id = parent_record_id,
        peptide_qc_record_id = protDiaNormInheritedRecordId(
            after,
            "peptide_qc_audit",
            "legacy_untracked_peptide_qc"
        ),
        protein_qc_record_id = protDiaNormInheritedRecordId(
            after,
            "protein_qc_audit",
            "legacy_untracked_protein_qc"
        ),
        resolved_parameters = bounded_parameters,
        software = software,
        before_summary = protDiaNormSummary(before),
        after_summary = protDiaNormSummary(after)
    )
    record_id <- paste0(
        "prot-norm:",
        substr(.peptideQcDigest(substantive), 1L, 24L)
    )
    record <- c(
        substantive,
        list(
            record_id = record_id,
            timestamp_utc = format(now, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
            canonical_digest = .peptideQcDigest(substantive)
        )
    )
    previous_audit <- protDiaNormArgs(before)$normalization_artifact_audit
    records <- previous_audit$records
    if (!is.list(records)) records <- list()
    records[[length(records) + 1L]] <- record
    after@args$normalization_artifact_audit <- list(
        schema = .PROT_DIA_NORM_AUDIT_SCHEMA,
        schema_version = .PROT_DIA_NORM_AUDIT_VERSION,
        records = records,
        current_record_id = record_id
    )
    after
}

protDiaNormTransition <- function(before, after, stage_id) {
    methods::is(before, "ProteinQuantitativeData") &&
        methods::is(after, "ProteinQuantitativeData") &&
        stage_id %in% .PROT_DIA_NORM_STAGES
}

protDiaNormWorkflowIsDia <- function(workflow_data, state_manager) {
    if (inherits(state_manager, "ArtifactWorkflowState")) {
        return(identical(workflowStateType(state_manager), "DIA"))
    }
    get_workflow_type <- tryCatch(
        state_manager$getWorkflowType,
        error = function(...) NULL
    )
    if (is.function(get_workflow_type)) {
        return(identical(get_workflow_type(), "DIA"))
    }
    workflow_type <- if (protDiaPeptideQcWorkflowData(workflow_data)) {
        workflow_data$config_list$globalParameters$workflow_type
    } else {
        NULL
    }
    identical(workflow_type, "DIA")
}

protDiaNormArtifactEligible <- function(
    workflow_data,
    state_manager,
    before,
    after,
    stage_id
) {
    if (!protDiaNormTransition(before, after, stage_id) ||
        identical(protDiaNormArtifactMode(stage_id), "disabled")) {
        return(FALSE)
    }
    if (inherits(state_manager, "ArtifactWorkflowState")) {
        return(identical(workflowStateType(state_manager), "DIA"))
    }
    protDiaPeptideQcWorkflowData(workflow_data) &&
        protDiaArtifactCoordinatorOwned(workflow_data)
}

protDiaNormValidateParent <- function(manager, before) {
    parent <- manager$getState()
    if (!identical(parent, before)) {
        protDiaArtifactAbort(
            "DIA-NN normalization parent differs from the active scientific state",
            "multischolar_prot_dia_normalization_parent_mismatch"
        )
    }
    invisible(parent)
}

protDiaNormSafeConfig <- function(config_object, software) {
    config <- config_object
    path_fields <- names(config)[grepl("(_dir|_path|_file)$", names(config))]
    for (field in path_fields) {
        value <- config[[field]]
        if (workflowStateScalarString(value)) config[[field]] <- basename(value)
    }
    config$artifact_provenance <- list(
        software = software,
        output_ownership = list(scope = "workflow_session")
    )
    config
}

protDiaNormSaveMemory <- function(
    state_manager,
    state_name,
    object,
    config_object,
    description,
    audit_metadata
) {
    save_state <- state_manager$saveState
    supported <- names(formals(save_state))
    args <- list(
        state_name = state_name,
        s4_data_object = object,
        config_object = config_object,
        description = description
    )
    if ("audit_metadata" %in% supported || "..." %in% supported) {
        args$audit_metadata <- audit_metadata
    }
    do.call(save_state, args)
}

saveProtDiaNormArtifactState <- function(
    workflow_data,
    state_manager,
    before,
    after,
    stage_id,
    state_name,
    config_object,
    description,
    parameters,
    status,
    transformation_type,
    now,
    failure_injector = NULL
) {
    if (!protDiaNormArtifactEligible(
        workflow_data,
        state_manager,
        before,
        after,
        stage_id
    )) {
        return(list(handled = FALSE, enabled = FALSE, object = after))
    }
    resources <- protDiaPeptideQcArtifactManager(workflow_data, state_manager)
    manager <- resources$manager
    if (isTRUE(resources$owned)) on.exit(manager$close(), add = TRUE)
    protDiaNormValidateParent(manager, before)
    parent_state <- manager$getCurrentStateName()
    memory_parent_state <- workflowStateCurrentName(state_manager)
    if (isTRUE(resources$owned) && !identical(parent_state, memory_parent_state)) {
        protDiaArtifactAbort(
            "DIA-NN normalization managers disagree on the active parent state",
            "multischolar_prot_dia_normalization_parent_mismatch"
        )
    }
    parent_generation_id <- manager$getCurrentGenerationId()
    after <- protDiaNormAttachParameters(after, parameters, stage_id)
    software <- protDiaNormSoftware(stage_id, parameters)
    after <- protDiaNormEmitAudit(
        before,
        after,
        stage_id,
        parent_state,
        state_name,
        parent_generation_id,
        parameters,
        status,
        transformation_type,
        software,
        now
    )
    record <- protDiaNormCurrentRecord(after)
    config_object$audit_record_id <- record$record_id
    artifactStoreInvokeFailure(
        failure_injector,
        "after_audit_creation",
        list(
            stage_id = stage_id,
            state_name = state_name,
            audit_record_id = record$record_id,
            parent_generation_id = parent_generation_id
        )
    )
    result <- manager$commitState(
        state_name = state_name,
        s4_data_object = after,
        config_object = protDiaNormSafeConfig(config_object, software),
        description = description,
        audit_metadata = record,
        failure_injector = failure_injector,
        action_id = artifactOpaqueId("action"),
        expected_parent_generation_id = parent_generation_id
    )
    if (!identical(result$status, "accepted")) {
        protDiaArtifactAbort(
            "DIA-NN normalization state did not advance its exact parent",
            "multischolar_prot_dia_normalization_commit_rejected",
            result = result
        )
    }
    if (isTRUE(resources$owned)) {
        memory_error <- tryCatch(
            {
                protDiaNormSaveMemory(
                    state_manager,
                    state_name,
                    after,
                    config_object,
                    description,
                    record
                )
                NULL
            },
            error = identity
        )
        if (inherits(memory_error, "error")) {
            restore_error <- tryCatch(
                {
                    manager$revertToState(parent_state)
                    NULL
                },
                error = identity
            )
            if (inherits(restore_error, "error")) {
                protDiaArtifactAbort(
                    "DIA-NN normalization compatibility save and rollback failed",
                    "multischolar_prot_dia_normalization_divergent_state",
                    memory_error = memory_error,
                    restore_error = restore_error
                )
            }
            stop(memory_error)
        }
    }
    hydrated <- manager$getState(state_name)
    if (!identical(hydrated, after) ||
        !isTRUE(methods::validObject(hydrated, test = TRUE))) {
        protDiaArtifactAbort(
            "DIA-NN normalization artifact hydration changed the scientific state",
            "multischolar_inexact_prot_dia_normalization_hydration"
        )
    }
    output <- list(
        handled = TRUE,
        enabled = TRUE,
        ok = TRUE,
        committed = TRUE,
        object = after,
        stage_id = stage_id,
        state_name = state_name,
        generation_id = result$generation_id,
        parent_generation_id = parent_generation_id,
        representation = result$representation,
        metrics = result$metrics,
        audit_record_id = record$record_id
    )
    recordProtDiaArtifactResult(workflow_data, stage_id, output)
    output
}

saveProtNormState <- function(
    workflow_data,
    state_manager,
    before,
    after,
    stage_id,
    state_name,
    config_object,
    description,
    parameters = config_object,
    status = "applied",
    transformation_type = "transformation",
    now = Sys.time(),
    failure_injector = NULL
) {
    if (!protDiaNormTransition(before, after, stage_id)) return(after)
    nondia_result <- saveProtNonDiaNormArtifactState(
        workflow_data,
        state_manager,
        before,
        after,
        stage_id,
        state_name,
        config_object,
        description,
        parameters,
        status,
        transformation_type,
        now,
        failure_injector
    )
    if (isTRUE(nondia_result$handled)) return(nondia_result$object)
    if (!protDiaNormWorkflowIsDia(workflow_data, state_manager)) {
        protDiaNormSaveMemory(
            state_manager,
            state_name,
            after,
            config_object,
            description,
            NULL
        )
        return(after)
    }
    artifact_result <- saveProtDiaNormArtifactState(
        workflow_data,
        state_manager,
        before,
        after,
        stage_id,
        state_name,
        config_object,
        description,
        parameters,
        status,
        transformation_type,
        now,
        failure_injector
    )
    if (isTRUE(artifact_result$handled)) return(artifact_result$object)
    parent_state <- workflowStateCurrentName(state_manager)
    after <- protDiaNormAttachParameters(after, parameters, stage_id)
    software <- protDiaNormSoftware(stage_id, parameters)
    after <- protDiaNormEmitAudit(
        before,
        after,
        stage_id,
        parent_state,
        state_name,
        paste0("memory:", parent_state),
        parameters,
        status,
        transformation_type,
        software,
        now
    )
    record <- protDiaNormCurrentRecord(after)
    config_object$audit_record_id <- record$record_id
    protDiaNormSaveMemory(
        state_manager,
        state_name,
        after,
        config_object,
        description,
        record
    )
    after
}

protDiaNormStateRef <- function(workflow_data, stage_id, state_name) {
    result <- workflow_data$artifact_stage_results[[stage_id]]
    if (!is.list(result) || !isTRUE(result$enabled) ||
        !workflowStateScalarString(result$generation_id)) {
        return(NULL)
    }
    list(
        schema = .PROT_DIA_NORM_REF_SCHEMA,
        schema_version = .PROT_DIA_NORM_REF_VERSION,
        stage_id = stage_id,
        state_name = state_name,
        generation_id = result$generation_id
    )
}

#' Resolve a normalization state ref by exact proteomics descriptor
#' @param workflow_data Mutable proteomics workflow state.
#' @param stage_id Normalization stage identifier.
#' @param state_name State name.
#' @return A payload-free generation ref, or `NULL`.
#' @noRd
protNormStateRef <- function(workflow_data, stage_id, state_name) {
    nondia <- protNonDiaNormStateRef(workflow_data, stage_id, state_name)
    if (!is.null(nondia)) return(nondia)
    protDiaNormStateRef(workflow_data, stage_id, state_name)
}

settleProtNormArtifactState <- function(
    workflow_data,
    norm_data,
    stage_id,
    state_name,
    object
) {
    ref <- protNormStateRef(workflow_data, stage_id, state_name)
    if (is.null(ref)) return(object)
    refs <- workflow_data$normalization_state_refs
    if (!is.list(refs)) refs <- list()
    refs[[stage_id]] <- ref
    workflow_data$normalization_state_refs <- refs
    norm_refs <- norm_data$state_refs
    if (!is.list(norm_refs)) norm_refs <- list()
    norm_refs[[stage_id]] <- ref
    norm_data$state_refs <- norm_refs
    if (identical(stage_id, "normalization")) {
        norm_data$normalized_protein_obj <- object
    } else if (identical(stage_id, "ruv_correction")) {
        norm_data$normalized_protein_obj <- NULL
        norm_data$ruv_normalized_obj <- object
    } else {
        norm_data$normalized_protein_obj <- NULL
        norm_data$ruv_normalized_obj <- NULL
        norm_data$correlation_filtered_obj <- NULL
        workflow_data$ruv_normalised_for_da_analysis_obj <- NULL
        workflow_data$final_for_da_ref <- ref
    }
    object
}

releaseProtNormArtifactStageObjects <- function(norm_data) {
    manager <- norm_data$state_manager
    if (!inherits(manager, "ArtifactWorkflowState")) return(invisible(FALSE))
    norm_data$normalized_protein_obj <- NULL
    norm_data$ruv_normalized_obj <- NULL
    workflowStateReleaseHydrationCache(manager)
}

resolveProtNormStateObject <- function(
    workflow_data = NULL,
    norm_data = NULL,
    state_names,
    legacy_object = NULL,
    stage_id = NULL
) {
    if (!is.null(legacy_object)) return(legacy_object)
    if (!is.null(stage_id) &&
        identical(protNormArtifactMode(workflow_data, stage_id), "disabled")) {
        return(NULL)
    }
    manager <- NULL
    if (protNormWorkflowDataAvailable(workflow_data)) {
        manager <- workflow_data$state_manager
    }
    if (is.null(manager) && !is.null(norm_data)) manager <- norm_data$state_manager
    get_state <- tryCatch(manager$getState, error = function(...) NULL)
    if (!is.function(get_state)) return(NULL)
    has_state <- tryCatch(manager$hasState, error = function(...) NULL)
    for (state_name in state_names) {
        available <- if (is.function(has_state)) {
            isTRUE(has_state(state_name))
        } else {
            identical(workflowStateCurrentName(manager), state_name)
        }
        if (isTRUE(available)) return(get_state(state_name))
    }
    NULL
}

protDiaNormDualManagerMode <- function(workflow_data, state_manager, stage_id) {
    !inherits(state_manager, "ArtifactWorkflowState") &&
        !identical(protDiaNormArtifactMode(stage_id), "disabled") &&
        protDiaPeptideQcWorkflowData(workflow_data) &&
        protDiaArtifactCoordinatorOwned(workflow_data)
}

revertProtDiaNormState <- function(workflow_data, state_name) {
    state_manager <- workflow_data$state_manager
    if (inherits(state_manager, "ArtifactWorkflowState")) {
        return(state_manager$revertToState(state_name))
    }
    if (!protDiaNormDualManagerMode(
        workflow_data,
        state_manager,
        "normalization"
    )) {
        return(state_manager$revertToState(state_name))
    }
    artifact_manager <- newProtDiaArtifactStateManager(workflow_data$workflow_context)
    on.exit(artifact_manager$close(), add = TRUE)
    memory_current <- state_manager$getCurrentStateName()
    artifact_current <- artifact_manager$getCurrentStateName()
    targets_match <- isTRUE(state_manager$hasState(state_name)) &&
        isTRUE(artifact_manager$hasState(state_name)) &&
        identical(memory_current, artifact_current) &&
        identical(state_manager$getState(), artifact_manager$getState()) &&
        identical(
            state_manager$getState(state_name),
            artifact_manager$getState(state_name)
        )
    if (!targets_match) {
        protDiaArtifactAbort(
            "DIA-NN normalization managers differ before revert",
            "multischolar_prot_dia_normalization_parent_mismatch"
        )
    }
    reverted_artifact <- artifact_manager$revertToState(state_name)
    memory_error <- tryCatch(
        {
            state_manager$revertToState(state_name)
            NULL
        },
        error = identity
    )
    if (inherits(memory_error, "error")) {
        protDiaPeptideQcRestoreRevert(
            artifact_manager,
            state_manager,
            artifact_current,
            memory_current
        )
        stop(memory_error)
    }
    reverted_memory <- state_manager$getState()
    if (!identical(reverted_artifact, reverted_memory)) {
        protDiaPeptideQcRestoreRevert(
            artifact_manager,
            state_manager,
            artifact_current,
            memory_current
        )
        protDiaArtifactAbort(
            "DIA-NN normalization revert hydrated different active states",
            "multischolar_inexact_prot_dia_normalization_hydration"
        )
    }
    reverted_memory
}

initializeProtNormWorkflowContext <- function(
    workflow_data,
    input,
    experiment_paths,
    omic_type,
    experiment_label,
    grouping_variable
) {
    if (!protNormWorkflowDataAvailable(workflow_data)) return(invisible(NULL))
    workflow_data$normalization_context <- list(
        workflow_context = workflow_data$workflow_context,
        experiment_paths = experiment_paths,
        omic_type = omic_type,
        experiment_label = experiment_label,
        normalization_method = input$norm_method,
        ruv_mode = input$ruv_mode,
        ruv_grouping_variable = grouping_variable,
        ruv_parameters = list(),
        correlation = list()
    )
    config <- workflow_data$config_list
    if (!is.list(config)) config <- list()
    if (!is.list(config$normaliseBetweenSamples)) {
        config$normaliseBetweenSamples <- list()
    }
    config$normaliseBetweenSamples$method <- input$norm_method
    workflow_data$config_list <- config
    invisible(workflow_data$normalization_context)
}

updateProtNormWorkflowRuvContext <- function(
    workflow_data,
    mode,
    grouping_variable,
    percentage,
    k,
    controls,
    optimization_result = NULL,
    input = NULL
) {
    if (!protNormWorkflowDataAvailable(workflow_data)) return(invisible(NULL))
    context <- workflow_data$normalization_context
    if (!is.list(context)) context <- list()
    trace <- optimization_result$optimization_results
    context$ruv_parameters <- list(
        ruv_mode = mode,
        ruv_grouping_variable = grouping_variable,
        percentage_as_neg_ctrl = percentage,
        ruv_k = k,
        control_genes_index = controls,
        optimization_results = trace,
        separation_metric = optimization_result$separation_metric_used,
        k_penalty_weight = optimization_result$k_penalty_weight,
        adaptive_k_penalty = optimization_result$adaptive_k_penalty_used,
        percentage_min = input$auto_percentage_min,
        percentage_max = input$auto_percentage_max,
        skip_reason = optimization_result$skip_reason
    )
    workflow_data$normalization_context <- context
    config <- workflow_data$config_list
    if (!is.list(config)) config <- list()
    config <- tryCatch(
        updateRuvParameters(config, k, controls, percentage),
        error = \(error) config
    )
    workflow_data$config_list <- config
    invisible(context$ruv_parameters)
}
