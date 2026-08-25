# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.LIPID_QC_AUDIT_SCHEMA <- "multischolar.lipidomics_qc_audit"
.LIPID_QC_AUDIT_VERSION <- 1L
.LIPID_QC_STAGES <- c(
    "intensity_filter",
    "duplicate_resolution",
    "internal_standard_analysis",
    "qc_finalization"
)

lipidQcArtifactMode <- function() {
    mode <- getOption("multischolar.lipidomics.qc_artifacts", "enabled")
    match.arg(mode, c("enabled", "disabled"))
}

lipidQcWorkflowData <- function(workflow_data) {
    is.environment(workflow_data) || inherits(workflow_data, "reactivevalues")
}

lipidQcPackageVersion <- function() {
    tryCatch(
        as.character(utils::packageVersion("MultiScholaR")),
        error = function(error) "development"
    )
}

lipidQcSoftware <- function(stage_id) {
    list(
        name = "MultiScholaR",
        version = lipidQcPackageVersion(),
        source = "R package",
        stage_id = stage_id
    )
}

lipidQcAssaySummary <- function(object) {
    assays <- object@lipid_data
    feature_column <- object@lipid_id_column
    lapply(seq_along(assays), function(index) {
        assay <- assays[[index]]
        feature_ids <- if (feature_column %in% names(assay)) {
            assay[[feature_column]]
        } else {
            character()
        }
        list(
            assay_name = names(assays)[[index]],
            assay_order = as.integer(index),
            rows = as.integer(nrow(assay)),
            columns = as.integer(ncol(assay)),
            distinct_features = as.integer(length(unique(feature_ids))),
            table_digest = artifactSemanticDigest(assay),
            feature_order_digest = artifactSemanticDigest(feature_ids)
        )
    })
}

lipidQcDuplicateProvenance <- function(before, after, stats) {
    assay_names <- names(before@lipid_data)
    per_assay <- lapply(assay_names, \(assay_name) {
        parent <- before@lipid_data[[assay_name]]
        retained <- after@lipid_data[[assay_name]]
        value_columns <- names(parent)[vapply(parent, is.numeric, logical(1))]
        class_columns <- grep(
            "class",
            names(parent),
            value = TRUE,
            ignore.case = TRUE
        )
        key_columns <- lipidQcStableKeyColumns(after, retained)
        key_data <- if (is.null(key_columns)) {
            data.frame()
        } else {
            retained[key_columns]
        }
        list(
            assay_name = assay_name,
            original_rows = as.integer(stats[[assay_name]]$original),
            retained_rows = as.integer(stats[[assay_name]]$resolved),
            removed_rows = as.integer(stats[[assay_name]]$removed),
            identity_columns = if (is.null(key_columns)) {
                character()
            } else {
                key_columns
            },
            identity_digest = artifactSemanticDigest(key_data),
            annotation_columns = intersect(
                unique(c(
                    before@annotation_id_column,
                    class_columns
                )),
                names(retained)
            ),
            annotation_digest = artifactSemanticDigest(
                retained[intersect(
                    unique(c(before@annotation_id_column, class_columns)),
                    names(retained)
                )]
            ),
            value_columns = value_columns,
            value_digest = artifactSemanticDigest(retained[value_columns]),
            retained_order_digest = artifactSemanticDigest(key_data)
        )
    })
    list(
        method = "highest_mean_intensity",
        representative = "maximum_row_mean",
        tie_break = "first_stable_input_row",
        aggregation = "none",
        missing_value_policy = "mean_na_rm_all_missing_to_negative_infinity",
        lipid_id_column = before@lipid_id_column,
        annotation_id_column = before@annotation_id_column,
        assays = per_assay
    )
}

lipidQcBuildAuditRecord <- function(
    before,
    after,
    stage_id,
    parent_state,
    current_state,
    parent_generation_id,
    parent_record_id,
    parameters,
    status,
    transformation_type,
    now = Sys.time()
) {
    if (!methods::is(before, "LipidomicsAssayData") ||
        !methods::is(after, "LipidomicsAssayData") ||
        !stage_id %in% .LIPID_QC_STAGES) {
        return(NULL)
    }
    substantive <- list(
        schema = .LIPID_QC_AUDIT_SCHEMA,
        schema_version = .LIPID_QC_AUDIT_VERSION,
        stage_id = stage_id,
        status = status,
        transformation_type = transformation_type,
        parent_state = parent_state,
        current_state = current_state,
        parent_generation_id = parent_generation_id,
        parent_record_id = parent_record_id,
        parameters = artifactCanonicalizeValue(parameters),
        software = lipidQcSoftware(stage_id),
        before_assays = lipidQcAssaySummary(before),
        after_assays = lipidQcAssaySummary(after),
        design_digest = artifactSemanticDigest(after@design_matrix),
        assay_order = names(after@lipid_data)
    )
    canonical_digest <- artifactSemanticDigest(substantive)
    record_id <- paste0("lipid-qc:", substr(canonical_digest, 1L, 24L))
    record <- c(substantive, list(
        record_id = record_id,
        timestamp_utc = format(now, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
        canonical_digest = canonical_digest
    ))
    record
}

lipidQcStateAuditRecord <- function(state_manager) {
    getter <- workflowStateMember(state_manager, "getStateAudit")
    if (!is.function(getter)) return(NULL)
    audit <- tryCatch(getter(), error = function(error) NULL)
    if (!is.list(audit) || !workflowStateScalarString(audit$record_id)) {
        return(NULL)
    }
    audit
}

lipidQcStableKeyColumns <- function(object, assay) {
    design_samples <- as.character(object@design_matrix[[object@sample_id]])
    metadata_columns <- setdiff(names(assay), design_samples)
    candidates <- unique(c(
        object@lipid_id_column,
        object@annotation_id_column,
        metadata_columns
    ))
    candidates <- candidates[
        !is.na(candidates) & nzchar(candidates) & candidates %in% names(assay)
    ]
    selected <- character()
    for (candidate in candidates) {
        selected <- c(selected, candidate)
        valid <- tryCatch({
            artifactWorkflowStateEntityKeys(assay, selected)
            TRUE
        }, error = function(error) FALSE)
        if (isTRUE(valid)) return(selected)
    }
    NULL
}

lipidQcAssayKeyColumns <- function(object) {
    keys <- lapply(object@lipid_data, function(assay) {
        lipidQcStableKeyColumns(object, assay)
    })
    if (any(vapply(keys, is.null, logical(1)))) return(NULL)
    keys
}

lipidQcRejectedReasons <- function(before, after, key_columns, stage_id) {
    reasons <- Map(function(parent, child, columns) {
        parent_keys <- artifactWorkflowStateEntityKeys(parent, columns)
        child_keys <- artifactWorkflowStateEntityKeys(child, columns)
        rejected <- parent_keys[!parent_keys %in% child_keys]
        stats::setNames(
            rep(paste0("removed_by_", stage_id), length(rejected)),
            rejected
        )
    }, before@lipid_data, after@lipid_data, key_columns)
    names(reasons) <- names(before@lipid_data)
    reasons
}

lipidQcSelectionHint <- function(
    before,
    after,
    stage_id,
    state_name,
    parent_state,
    parameters,
    transformation_type,
    record
) {
    if (!identical(transformation_type, "filter") ||
        !identical(stage_id, "intensity_filter") ||
        !identical(names(before@lipid_data), names(after@lipid_data))) {
        return(NULL)
    }
    other_slots <- setdiff(methods::slotNames(before), c("lipid_data", "args"))
    if (!all(vapply(other_slots, function(slot_name) {
        identical(
            methods::slot(before, slot_name),
            methods::slot(after, slot_name)
        )
    }, logical(1)))) {
        return(NULL)
    }
    key_columns <- lipidQcAssayKeyColumns(before)
    if (is.null(key_columns)) return(NULL)
    if (is.null(record)) {
        lipidArtifactAbort(
            "lipidomics QC selection requires its audit record",
            "multischolar_missing_lipidomics_qc_audit"
        )
    }
    hint <- newArtifactAssaySelectionHint(
        slot_name = "lipid_data",
        assay_key_columns = key_columns,
        method = stage_id,
        normalized_parameters = artifactCanonicalizeValue(parameters),
        software = lipidQcSoftware(stage_id)[c("name", "version", "source")],
        lineage = list(
            audit_record_id = record$record_id,
            state_name = state_name,
            parent_state = parent_state,
            parent_record_id = record$parent_record_id
        ),
        rejected_reasons = lipidQcRejectedReasons(
            before,
            after,
            key_columns,
            stage_id
        )
    )
    tryCatch(
        {
            artifactWorkflowStateAssaySelectionPlan(before, after, hint)
            hint
        },
        multischolar_artifact_selection_requires_materialization =
            function(error) NULL,
        multischolar_invalid_artifact_assay_selection = function(error) NULL,
        multischolar_invalid_artifact_row_selection = function(error) NULL,
        multischolar_ambiguous_artifact_row_selection = function(error) NULL
    )
}

lipidQcStateReferenceHint <- function(
    before,
    after,
    stage_id,
    state_name,
    parent_state,
    parameters,
    record
) {
    eligible <- stage_id %in% c(
        "internal_standard_analysis",
        "qc_finalization"
    ) && identical(before, after)
    if (!isTRUE(eligible)) return(NULL)
    newArtifactStateReferenceHint(
        method = stage_id,
        normalized_parameters = artifactCanonicalizeValue(parameters),
        software = lipidQcSoftware(stage_id)[c("name", "version", "source")],
        lineage = list(
            audit_record_id = record$record_id,
            state_name = state_name,
            parent_state = parent_state,
            parent_record_id = record$parent_record_id
        )
    )
}

lipidQcArtifactEligible <- function(workflow_data, state_manager, before, after) {
    valid_transition <- methods::is(before, "LipidomicsAssayData") &&
        methods::is(after, "LipidomicsAssayData")
    if (!isTRUE(valid_transition)) return(FALSE)
    if (inherits(state_manager, "ArtifactWorkflowState")) {
        return(identical(workflowStateType(state_manager), "lipidomics_standard"))
    }
    !identical(lipidQcArtifactMode(), "disabled") &&
        lipidQcWorkflowData(workflow_data) &&
        lipidArtifactCoordinatorOwned(workflow_data)
}

lipidQcArtifactManager <- function(workflow_data, state_manager) {
    if (inherits(state_manager, "ArtifactWorkflowState")) {
        return(list(manager = state_manager, owned = FALSE))
    }
    prepared <- prepareLipidArtifactContext(workflow_data)
    if (!isTRUE(prepared$enabled)) {
        lipidArtifactAbort(
            "lipidomics QC has no exact artifact context",
            "multischolar_invalid_lipidomics_artifact_context"
        )
    }
    manager <- newLipidArtifactStateManager(prepared)
    manager$setWorkflowType("lipidomics_standard")
    list(manager = manager, owned = TRUE)
}

lipidQcValidateParent <- function(manager, before) {
    parent <- manager$getState()
    if (!identical(parent, before)) {
        lipidArtifactAbort(
            "lipidomics QC parent differs from the active scientific state",
            "multischolar_lipidomics_qc_parent_mismatch"
        )
    }
    invisible(parent)
}

lipidQcSafeConfig <- function(config_object, software) {
    config <- config_object
    file_fields <- names(config)[grepl("_file$", names(config))]
    owned_files <- list()
    for (field in file_fields) {
        value <- config[[field]]
        if (workflowStateScalarString(value)) {
            owned_files[[field]] <- basename(value)
            config[[field]] <- basename(value)
        }
    }
    config$artifact_provenance <- list(
        software = software,
        output_ownership = list(scope = "workflow_session", files = owned_files)
    )
    config
}

lipidQcSaveMemory <- function(
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

lipidQcCurrentStateName <- function(state_manager) {
    state_name <- workflowStateCurrentName(state_manager)
    if (workflowStateScalarString(state_name)) return(state_name)
    history_getter <- workflowStateMember(state_manager, "getHistory")
    history <- if (is.function(history_getter)) history_getter() else character()
    if (length(history) > 0L && workflowStateScalarString(tail(history, 1L))) {
        return(tail(history, 1L))
    }
    "legacy_current_state"
}

saveLipidQcState <- function(
    workflow_data,
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
    state_manager <- workflow_data$state_manager
    parent_state <- lipidQcCurrentStateName(state_manager)
    artifact_eligible <- lipidQcArtifactEligible(
        workflow_data,
        state_manager,
        before,
        after
    )
    resources <- if (isTRUE(artifact_eligible)) {
        lipidQcArtifactManager(workflow_data, state_manager)
    } else {
        NULL
    }
    manager <- if (is.null(resources)) NULL else resources$manager
    if (!is.null(resources) && isTRUE(resources$owned)) {
        on.exit(manager$close(), add = TRUE)
    }
    parent_generation_id <- if (is.null(manager)) {
        paste0("memory:", parent_state)
    } else {
        lipidQcValidateParent(manager, before)
        artifact_parent <- manager$getCurrentStateName()
        if (isTRUE(resources$owned) && !identical(artifact_parent, parent_state)) {
            lipidArtifactAbort(
                "lipidomics QC managers disagree on the active parent",
                "multischolar_lipidomics_qc_parent_mismatch"
            )
        }
        manager$getCurrentGenerationId()
    }
    parent_audit <- if (is.null(manager)) {
        lipidQcStateAuditRecord(state_manager)
    } else {
        lipidQcStateAuditRecord(manager)
    }
    parent_record_id <- if (is.null(parent_audit)) {
        "legacy_or_design_checkpoint"
    } else {
        parent_audit$record_id
    }
    record <- lipidQcBuildAuditRecord(
        before,
        after,
        stage_id,
        parent_state,
        state_name,
        parent_generation_id,
        parent_record_id,
        parameters,
        status,
        transformation_type,
        now
    )
    if (is.null(manager)) {
        lipidQcSaveMemory(
            state_manager,
            state_name,
            after,
            config_object,
            description,
            record
        )
        return(after)
    }
    hint <- lipidQcStateReferenceHint(
        before,
        after,
        stage_id,
        state_name,
        parent_state,
        parameters,
        record
    )
    if (is.null(hint)) {
        hint <- lipidQcSelectionHint(
            before,
            after,
            stage_id,
            state_name,
            parent_state,
            parameters,
            transformation_type,
            record
        )
    }
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
        config_object = lipidQcSafeConfig(
            config_object,
            lipidQcSoftware(stage_id)
        ),
        description = description,
        audit_metadata = record,
        persistence_hint = hint,
        failure_injector = failure_injector,
        action_id = artifactOpaqueId("action"),
        expected_parent_generation_id = parent_generation_id
    )
    if (!identical(result$status, "accepted")) {
        lipidArtifactAbort(
            "lipidomics QC commit did not advance its exact parent",
            "multischolar_lipidomics_qc_commit_rejected",
            result = result
        )
    }
    if (isTRUE(resources$owned)) {
        memory_error <- tryCatch({
            lipidQcSaveMemory(
                state_manager,
                state_name,
                after,
                config_object,
                description,
                record
            )
            NULL
        }, error = identity)
        if (inherits(memory_error, "error")) {
            restore_error <- tryCatch({
                manager$revertToState(parent_state)
                NULL
            }, error = identity)
            if (inherits(restore_error, "error")) {
                lipidArtifactAbort(
                    "lipidomics QC memory save and artifact rollback both failed",
                    "multischolar_lipidomics_qc_divergent_state",
                    memory_error = memory_error,
                    restore_error = restore_error
                )
            }
            stop(memory_error)
        }
    }
    hydrated <- manager$getState(state_name)
    if (!identical(hydrated, after) ||
        !identical(methods::validObject(hydrated, test = TRUE), TRUE)) {
        lipidArtifactAbort(
            "lipidomics QC artifact hydration changed scientific state",
            "multischolar_inexact_lipidomics_qc_hydration"
        )
    }
    output <- list(
        enabled = TRUE,
        ok = TRUE,
        committed = TRUE,
        stage_id = stage_id,
        state_name = state_name,
        generation_id = result$generation_id,
        parent_generation_id = parent_generation_id,
        representation = result$representation,
        metrics = result$metrics,
        audit_record_id = record$record_id,
        source_payloads_retained = TRUE,
        eviction_performed = FALSE
    )
    recordArtifactStageResult(workflow_data, stage_id, output)
    after
}

lipidQcDualManagerMode <- function(workflow_data, state_manager) {
    !inherits(state_manager, "ArtifactWorkflowState") &&
        !identical(lipidQcArtifactMode(), "disabled") &&
        lipidQcWorkflowData(workflow_data) &&
        lipidArtifactCoordinatorOwned(workflow_data)
}

lipidQcRestoreRevert <- function(
    artifact_manager,
    state_manager,
    artifact_state,
    memory_state
) {
    errors <- list(
        artifact = tryCatch({
            artifact_manager$revertToState(artifact_state)
            NULL
        }, error = identity),
        memory = tryCatch({
            state_manager$revertToState(memory_state)
            NULL
        }, error = identity)
    )
    errors <- Filter(function(error) inherits(error, "error"), errors)
    if (length(errors) > 0L) {
        lipidArtifactAbort(
            "lipidomics QC revert rollback could not restore both managers",
            "multischolar_lipidomics_qc_divergent_state",
            restore_errors = errors
        )
    }
    invisible(TRUE)
}

revertLipidQcState <- function(workflow_data, state_name) {
    state_manager <- workflow_data$state_manager
    if (inherits(state_manager, "ArtifactWorkflowState") ||
        !lipidQcDualManagerMode(workflow_data, state_manager)) {
        return(state_manager$revertToState(state_name))
    }
    resources <- lipidQcArtifactManager(workflow_data, state_manager)
    artifact_manager <- resources$manager
    on.exit(artifact_manager$close(), add = TRUE)
    memory_current <- lipidQcCurrentStateName(state_manager)
    artifact_current <- artifact_manager$getCurrentStateName()
    valid <- identical(memory_current, artifact_current) &&
        isTRUE(state_manager$hasState(state_name)) &&
        isTRUE(artifact_manager$hasState(state_name)) &&
        identical(state_manager$getState(), artifact_manager$getState()) &&
        identical(
            state_manager$getState(state_name),
            artifact_manager$getState(state_name)
        )
    if (!isTRUE(valid)) {
        lipidArtifactAbort(
            "lipidomics QC managers differ before revert",
            "multischolar_lipidomics_qc_parent_mismatch"
        )
    }
    reverted_artifact <- artifact_manager$revertToState(state_name)
    memory_error <- tryCatch({
        state_manager$revertToState(state_name)
        NULL
    }, error = identity)
    if (inherits(memory_error, "error")) {
        lipidQcRestoreRevert(
            artifact_manager,
            state_manager,
            artifact_current,
            memory_current
        )
        stop(memory_error)
    }
    reverted_memory <- state_manager$getState()
    if (!identical(reverted_artifact, reverted_memory)) {
        lipidQcRestoreRevert(
            artifact_manager,
            state_manager,
            artifact_current,
            memory_current
        )
        lipidArtifactAbort(
            "lipidomics QC revert hydrated different active states",
            "multischolar_inexact_lipidomics_qc_hydration"
        )
    }
    recordArtifactStageResult(workflow_data, "qc_revert", list(
        enabled = TRUE,
        ok = TRUE,
        stage_id = "qc_revert",
        state_name = state_name,
        artifact_generation_id = artifact_manager$getCurrentGenerationId()
    ))
    reverted_memory
}

lipidQcItsdProvenance <- function(object, analysis, input_pattern) {
    metrics <- analysis$metricsByAssay
    if (is.null(metrics) && is.data.frame(analysis$metrics)) {
        metrics <- split(analysis$metrics, analysis$metrics$assay)
    }
    assay_names <- names(object@lipid_data)
    per_assay <- lapply(assay_names, function(assay_name) {
        value <- metrics[[assay_name]]
        matched_ids <- if (is.data.frame(value) && "is_id" %in% names(value)) {
            as.character(value$is_id)
        } else {
            character()
        }
        requested_match <- if (length(matched_ids) == 0L) {
            FALSE
        } else {
            tryCatch(
                all(grepl(analysis$pattern, matched_ids)),
                error = \(error) FALSE
            )
        }
        list(
            assay_name = assay_name,
            status = if (is.data.frame(value) && nrow(value) > 0L) {
                "matched"
            } else {
                "no_match"
            },
            matched_ids = matched_ids,
            match_count = if (is.data.frame(value)) {
                as.integer(nrow(value))
            } else {
                0L
            },
            matched_column = analysis$matchedColumns[[assay_name]],
            searched_columns = analysis$searchedColumns[[assay_name]],
            match_source = if (length(matched_ids) == 0L) {
                "none"
            } else if (isTRUE(requested_match)) {
                "resolved_regex"
            } else {
                "common_internal_standard_fallback_v1"
            }
        )
    })
    list(
        requested_regex = if (is.null(input_pattern)) "" else input_pattern,
        resolved_regex = analysis$pattern,
        regex_source = if (!is.null(input_pattern) && nzchar(input_pattern)) {
            "user_input"
        } else {
            "LipidomicsAssayData@internal_standard_regex"
        },
        assay_scope = assay_names,
        assays = per_assay,
        total_matches = as.integer(analysis$nIsTotal),
        failure_reason = if (as.integer(analysis$nIsTotal) == 0L) {
            "no_internal_standard_matches"
        } else {
            NULL
        },
        parameters = list(
            cv_threshold_percent = 30,
            matcher = "MultiScholaR::getLipidInternalStandardMetrics",
            fallback_policy_version = "1.0.0"
        )
    )
}

lipidQcItsdFailureProvenance <- function(
    object,
    input_pattern,
    error,
    searched_columns = list()
) {
    assay_names <- if (methods::is(object, "LipidomicsAssayData")) {
        names(object@lipid_data)
    } else {
        character()
    }
    list(
        requested_regex = if (is.null(input_pattern)) "" else input_pattern,
        resolved_regex = if (is.null(input_pattern)) "" else input_pattern,
        regex_source = if (!is.null(input_pattern) && nzchar(input_pattern)) {
            "user_input"
        } else {
            "LipidomicsAssayData@internal_standard_regex"
        },
        assay_scope = assay_names,
        assays = lapply(assay_names, function(assay_name) {
            list(
                assay_name = assay_name,
                status = "no_match",
                matched_ids = character(),
                match_count = 0L,
                matched_column = NULL,
                searched_columns = searched_columns[[assay_name]],
                match_source = "none"
            )
        }),
        total_matches = 0L,
        failure_reason = conditionMessage(error),
        parameters = list(
            cv_threshold_percent = 30,
            matcher = "MultiScholaR::getLipidInternalStandardMetrics",
            fallback_policy_version = "1.0.0"
        )
    )
}

recordLipidQcItsdAnalysis <- function(
    workflow_data,
    current_s4,
    analysis,
    input_pattern
) {
    manager <- workflow_data$state_manager
    if (is.null(manager) || !is.function(manager$saveState)) {
        return(invisible(list(enabled = FALSE, reason = "state_manager_not_writable")))
    }
    provenance <- lipidQcItsdProvenance(current_s4, analysis, input_pattern)
    skipped <- identical(provenance$total_matches, 0L)
    saved <- saveLipidQcState(
        workflow_data,
        current_s4,
        current_s4,
        stage_id = "internal_standard_analysis",
        state_name = if (skipped) {
            "lipid_internal_standards_skipped"
        } else {
            "lipid_internal_standards_analyzed"
        },
        config_object = workflow_data$config_list,
        description = "Recorded lipidomics internal-standard QC analysis",
        parameters = provenance,
        status = if (skipped) "skipped" else "analyzed",
        transformation_type = "analysis"
    )
    workflow_data$lipid_qc_itsd_provenance <- provenance
    invisible(list(enabled = TRUE, object = saved, provenance = provenance))
}

recordLipidQcItsdFailure <- function(
    workflow_data,
    current_s4,
    input_pattern,
    error,
    searched_columns = list()
) {
    provenance <- lipidQcItsdFailureProvenance(
        current_s4,
        input_pattern,
        error,
        searched_columns
    )
    workflow_data$lipid_qc_itsd_provenance <- provenance
    manager <- workflow_data$state_manager
    writable <- methods::is(current_s4, "LipidomicsAssayData") &&
        !is.null(manager) && is.function(manager$saveState)
    if (!isTRUE(writable)) {
        return(invisible(list(enabled = FALSE, provenance = provenance)))
    }
    saved <- saveLipidQcState(
        workflow_data,
        current_s4,
        current_s4,
        stage_id = "internal_standard_analysis",
        state_name = "lipid_internal_standards_skipped",
        config_object = workflow_data$config_list,
        description = "Recorded skipped lipidomics internal-standard QC analysis",
        parameters = provenance,
        status = "skipped",
        transformation_type = "analysis"
    )
    invisible(list(enabled = TRUE, object = saved, provenance = provenance))
}

lipidQcProgressOwner <- function(workflow_data) {
    if (!lipidQcWorkflowData(workflow_data)) return(NULL)
    manager <- workflow_data$state_manager
    artifact_mode <- inherits(manager, "ArtifactWorkflowState") || tryCatch(
        lipidArtifactCoordinatorOwned(workflow_data),
        error = function(error) FALSE
    )
    if (isTRUE(artifact_mode)) workflow_data else NULL
}

initializeLipidQcProjectContext <- function(
    workflow_data,
    experiment_paths,
    omic_type,
    experiment_label
) {
    scalar_path <- function(value) {
        if (workflowStateScalarString(value)) value else NULL
    }
    workflow_data$lipid_qc_project_context <- list(
        omic_type = omic_type,
        experiment_label = experiment_label,
        time_dir = scalar_path(experiment_paths$time_dir),
        publication_graphs_dir = scalar_path(
            experiment_paths$publication_graphs_dir
        )
    )
    invisible(workflow_data$lipid_qc_project_context)
}

lipidQcTrackingContextArgs <- function(workflow_data, supported) {
    if (!lipidQcWorkflowData(workflow_data)) return(list())
    context <- workflow_data$lipid_qc_project_context
    if (!is.list(context)) return(list())
    candidates <- list(
        omics_type = context$omic_type,
        time_dir = context$time_dir,
        publication_graphs_dir = context$publication_graphs_dir,
        use_explicit_paths = TRUE
    )
    candidates <- candidates[!vapply(candidates, is.null, logical(1))]
    if ("..." %in% supported) return(candidates)
    candidates[intersect(names(candidates), supported)]
}

invokeLipidQcTracking <- function(
    update_fn,
    args,
    workflow_data = NULL
) {
    owner <- lipidQcProgressOwner(workflow_data)
    supported <- names(formals(update_fn))
    if (!is.null(owner) &&
        ("progress_owner" %in% supported || "..." %in% supported)) {
        args$progress_owner <- owner
    }
    context_args <- if (is.null(owner)) {
        list()
    } else {
        lipidQcTrackingContextArgs(workflow_data, supported)
    }
    args[names(context_args)] <- context_args
    do.call(update_fn, args)
}

releaseLipidQcSessionOwnership <- function(workflow_data) {
    releaseFilteringProgressLipidomics(workflow_data)
    workflow_data$lipid_qc_project_context <- NULL
    workflow_data$lipid_qc_itsd_provenance <- NULL
    invisible(TRUE)
}

initializeLipidQcProgressOwnership <- function(workflow_data, session) {
    owner <- lipidQcProgressOwner(workflow_data)
    if (!is.null(owner)) getFilteringProgressLipidomics(owner)
    if (is.function(session$onSessionEnded)) {
        session$onSessionEnded(function() {
            releaseLipidQcSessionOwnership(workflow_data)
        })
    }
    invisible(!is.null(owner))
}
