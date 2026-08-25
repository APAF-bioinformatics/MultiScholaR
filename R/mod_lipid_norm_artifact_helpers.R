# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.LIPID_NORM_AUDIT_SCHEMA <- "multischolar.lipidomics_normalization_audit"
.LIPID_NORM_AUDIT_VERSION <- 1L
.LIPID_NORM_STAGES <- c(
    "itsd_normalization",
    "itsd_skip",
    "log2_transformation",
    "between_sample_normalization",
    "ruv_correction",
    "ruv_skip",
    "correlation_filter",
    "correlation_skip",
    "normalization_reset"
)

lipidNormArtifactMode <- function() {
    mode <- getOption("multischolar.lipidomics.norm_artifacts", "enabled")
    match.arg(mode, c("enabled", "disabled"))
}

lipidNormSafeDigest <- function(value) {
    tryCatch(
        artifactSemanticDigest(value),
        error = function(error) {
            digest::digest(
                serialize(value, NULL, version = 3L),
                algo = "sha256",
                serialize = FALSE
            )
        }
    )
}

lipidNormArtifactActive <- function(workflow_data) {
    manager <- workflow_data$state_manager
    coordinator_owned <- tryCatch(
        lipidArtifactCoordinatorOwned(workflow_data),
        error = \(error) FALSE
    )
    inherits(manager, "ArtifactWorkflowState") || (
        !identical(lipidNormArtifactMode(), "disabled") &&
        lipidQcWorkflowData(workflow_data) &&
        isTRUE(coordinator_owned)
    )
}

storeLipidNormIntermediate <- function(
    workflow_data,
    norm_data,
    field,
    object,
    state_name
) {
    if (!isTRUE(lipidNormArtifactActive(workflow_data))) {
        norm_data[[field]] <- object
        return(invisible(object))
    }
    refs <- norm_data$artifact_state_refs
    if (!is.list(refs)) refs <- list()
    refs[[field]] <- state_name
    norm_data$artifact_state_refs <- refs
    norm_data[[field]] <- NULL
    invisible(state_name)
}

resolveLipidNormIntermediate <- function(workflow_data, norm_data, field) {
    value <- norm_data[[field]]
    if (!is.null(value)) return(value)
    refs <- norm_data$artifact_state_refs
    state_name <- if (is.list(refs)) refs[[field]] else NULL
    if (!workflowStateScalarString(state_name)) return(NULL)
    workflow_data$state_manager$getState(state_name)
}

resolveLipidNormLatest <- function(workflow_data, norm_data, fields) {
    for (field in fields) {
        value <- resolveLipidNormIntermediate(workflow_data, norm_data, field)
        if (!is.null(value)) return(value)
    }
    NULL
}

lipidNormArgsSummary <- function(object) {
    args <- object@args
    arg_names <- names(args)
    if (is.null(arg_names)) arg_names <- character()
    list(
        names = sort(arg_names, method = "radix"),
        semantic_digest = lipidNormSafeDigest(args),
        object_bytes = unname(as.numeric(utils::object.size(args))),
        storage = "exact_s4_codec_payload_refs",
        compatible_getter = "ArtifactWorkflowState$getState"
    )
}

lipidNormBuildAuditRecord <- function(
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
        !stage_id %in% .LIPID_NORM_STAGES) {
        return(NULL)
    }
    substantive <- list(
        schema = .LIPID_NORM_AUDIT_SCHEMA,
        schema_version = .LIPID_NORM_AUDIT_VERSION,
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
        before_design_digest = artifactSemanticDigest(before@design_matrix),
        after_design_digest = artifactSemanticDigest(after@design_matrix),
        after_args = lipidNormArgsSummary(after),
        assay_order = names(after@lipid_data)
    )
    canonical_digest <- artifactSemanticDigest(substantive)
    c(substantive, list(
        record_id = paste0("lipid-norm:", substr(canonical_digest, 1L, 24L)),
        timestamp_utc = format(now, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
        canonical_digest = canonical_digest
    ))
}

lipidNormControlSummary <- function(best_k_per_assay, ctrl_per_assay) {
    assay_names <- union(names(best_k_per_assay), names(ctrl_per_assay))
    lapply(assay_names, function(assay_name) {
        controls <- ctrl_per_assay[[assay_name]]
        control_count <- if (is.null(controls)) {
            0L
        } else if (is.logical(controls)) {
            as.integer(sum(controls, na.rm = TRUE))
        } else {
            as.integer(length(controls))
        }
        list(
            assay_name = assay_name,
            k = best_k_per_assay[[assay_name]],
            control_kind = if (is.null(controls)) {
                "none"
            } else if (is.logical(controls)) {
                "logical_mask"
            } else {
                "indices"
            },
            control_count = control_count,
            control_length = if (is.null(controls)) 0L else length(controls),
            control_digest = lipidNormSafeDigest(controls)
        )
    })
}

lipidNormOptimizerSummary <- function(optimization_state) {
    if (is.null(optimization_state)) {
        return(list(status = "not_recorded", assays = list()))
    }
    results <- optimization_state$ruvResults
    if (!is.list(results)) results <- list()
    list(
        status = "recorded",
        parameters = optimization_state$ruvParams,
        result_digest = lipidNormSafeDigest(results),
        assays = lapply(names(results), function(assay_name) {
            result <- results[[assay_name]]
            list(
                assay_name = assay_name,
                success = isTRUE(result$success),
                best_k = result$best_k,
                best_percentage = result$best_percentage,
                error = result$error
            )
        })
    )
}

lipidNormFeatureSelectionSummary <- function(feature_ids, object = NULL) {
    if (is.null(feature_ids)) {
        assays <- if (!methods::is(object, "LipidomicsAssayData")) {
            list()
        } else {
            lapply(names(object@lipid_data), \(assay_name) {
                assay <- object@lipid_data[[assay_name]]
                annotation_column <- object@annotation_id_column
                matched <- if (annotation_column %in% names(assay) &&
                    workflowStateScalarString(object@internal_standard_regex)) {
                    grepl(
                        object@internal_standard_regex,
                        as.character(assay[[annotation_column]])
                    )
                } else {
                    rep(FALSE, nrow(assay))
                }
                ids <- as.character(assay[[object@lipid_id_column]][matched])
                list(
                    assay_name = assay_name,
                    searched_columns = annotation_column,
                    count = as.integer(length(ids)),
                    identity_digest = lipidNormSafeDigest(ids)
                )
            })
        }
        return(list(
            source = "object_internal_standard_regex",
            regex = if (methods::is(object, "LipidomicsAssayData")) {
                object@internal_standard_regex
            } else {
                NULL
            },
            assays = assays
        ))
    }
    values <- if (is.list(feature_ids)) feature_ids else list(all = feature_ids)
    list(
        source = "manual_assay_selection",
        assays = lapply(names(values), function(assay_name) {
            ids <- as.character(values[[assay_name]])
            list(
                assay_name = assay_name,
                count = as.integer(length(ids)),
                identity_digest = lipidNormSafeDigest(ids)
            )
        })
    )
}

lipidNormCorrelationSummary <- function(results) {
    values <- if (is.list(results)) results else list(all = results)
    list(
        assays = lapply(names(values), function(assay_name) {
            value <- values[[assay_name]]
            list(
                assay_name = assay_name,
                result_count = if (is.data.frame(value)) nrow(value) else length(value),
                semantic_digest = lipidNormSafeDigest(value)
            )
        })
    )
}

lipidNormStateReferenceHint <- function(
    before,
    after,
    stage_id,
    state_name,
    parent_state,
    parameters,
    transformation_type,
    record
) {
    if (!identical(transformation_type, "no_op") ||
        !identical(before, after)) {
        return(NULL)
    }
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

lipidNormArtifactEligible <- function(
    workflow_data,
    state_manager,
    before,
    after
) {
    valid <- methods::is(before, "LipidomicsAssayData") &&
        methods::is(after, "LipidomicsAssayData")
    if (!isTRUE(valid)) return(FALSE)
    if (inherits(state_manager, "ArtifactWorkflowState")) {
        return(identical(workflowStateType(state_manager), "lipidomics_standard"))
    }
    !identical(lipidNormArtifactMode(), "disabled") &&
        lipidQcWorkflowData(workflow_data) &&
        lipidArtifactCoordinatorOwned(workflow_data)
}

saveLipidNormState <- function(
    workflow_data,
    before,
    after,
    stage_id,
    state_name,
    config_object,
    description,
    parameters,
    status = "applied",
    transformation_type = "transformation",
    now = Sys.time(),
    failure_injector = NULL
) {
    state_manager <- workflow_data$state_manager
    parent_state <- lipidQcCurrentStateName(state_manager)
    enabled <- lipidNormArtifactEligible(
        workflow_data,
        state_manager,
        before,
        after
    )
    resources <- if (isTRUE(enabled)) {
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
                "lipidomics normalization managers disagree on the active parent",
                "multischolar_lipidomics_norm_parent_mismatch"
            )
        }
        manager$getCurrentGenerationId()
    }
    parent_audit <- lipidQcStateAuditRecord(
        if (is.null(manager)) state_manager else manager
    )
    parent_record_id <- if (is.null(parent_audit)) {
        "legacy_or_qc_checkpoint"
    } else {
        parent_audit$record_id
    }
    record <- lipidNormBuildAuditRecord(
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
    hint <- lipidNormStateReferenceHint(
        before,
        after,
        stage_id,
        state_name,
        parent_state,
        parameters,
        transformation_type,
        record
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
            "lipidomics normalization commit did not advance its exact parent",
            "multischolar_lipidomics_norm_commit_rejected",
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
                    paste(
                        "lipidomics normalization memory save and artifact",
                        "rollback both failed"
                    ),
                    "multischolar_lipidomics_norm_divergent_state",
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
            "lipidomics normalization artifact changed scientific state",
            "multischolar_inexact_lipidomics_norm_hydration"
        )
    }
    recordArtifactStageResult(workflow_data, stage_id, list(
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
    ))
    after
}

recordLipidNormNoOp <- function(
    workflow_data,
    current_s4,
    stage_id,
    parameters,
    status = "skipped"
) {
    if (!isTRUE(lipidNormArtifactActive(workflow_data))) return(current_s4)
    writable <- lipidQcWorkflowData(workflow_data) &&
        !is.null(workflow_data$state_manager) &&
        is.function(workflow_data$state_manager$saveState)
    if (!isTRUE(writable)) return(current_s4)
    state_name <- lipidQcCurrentStateName(workflow_data$state_manager)
    saveLipidNormState(
        workflow_data,
        current_s4,
        current_s4,
        stage_id,
        state_name,
        workflow_data$config_list,
        paste(stage_id, status),
        parameters,
        status = status,
        transformation_type = "no_op"
    )
}

persistLipidNormItsdState <- function(
    workflow_data,
    norm_data,
    before,
    after,
    aggregation,
    feature_ids,
    description
) {
    after <- saveLipidNormState(
        workflow_data,
        before,
        after,
        "itsd_normalization",
        "lipid_itsd_norm",
        workflow_data$config_list,
        description,
        parameters = list(
            method = "ITSD",
            aggregation = aggregation,
            feature_selection = lipidNormFeatureSelectionSummary(
                feature_ids,
                before
            ),
            internal_standard_regex = before@internal_standard_regex
        ),
        transformation_type = "normalization"
    )
    storeLipidNormIntermediate(
        workflow_data,
        norm_data,
        "post_itsd_obj",
        after,
        "lipid_itsd_norm"
    )
    after
}

persistLipidNormLog2State <- function(
    workflow_data,
    norm_data,
    before,
    after,
    offset,
    description
) {
    after <- saveLipidNormState(
        workflow_data,
        before,
        after,
        "log2_transformation",
        "lipid_log2",
        workflow_data$config_list,
        description,
        parameters = list(
            transform = "log2",
            offset = offset,
            non_finite_policy = "preserve_as_NA"
        ),
        transformation_type = "transformation"
    )
    storeLipidNormIntermediate(
        workflow_data,
        norm_data,
        "post_log2_obj",
        after,
        "lipid_log2"
    )
    after
}

persistLipidNormBetweenSampleState <- function(
    workflow_data,
    before,
    after,
    method,
    description
) {
    skipped <- identical(method, "none")
    saveLipidNormState(
        workflow_data,
        before,
        after,
        "between_sample_normalization",
        "lipid_normalized",
        workflow_data$config_list,
        description,
        parameters = list(method = method),
        status = if (skipped) "skipped" else "applied",
        transformation_type = if (skipped) "no_op" else "normalization"
    )
}

persistLipidNormRuvState <- function(
    workflow_data,
    before,
    after,
    grouping_variable,
    best_k_per_assay,
    ctrl_per_assay,
    ruv_mode,
    optimization_state,
    description
) {
    saveLipidNormState(
        workflow_data,
        before,
        after,
        "ruv_correction",
        "lipid_ruv_corrected",
        workflow_data$config_list,
        description,
        parameters = list(
            mode = if (is.null(ruv_mode)) "not_recorded" else ruv_mode,
            grouping_variable = grouping_variable,
            assays = lipidNormControlSummary(
                best_k_per_assay,
                ctrl_per_assay
            ),
            optimizer = lipidNormOptimizerSummary(optimization_state)
        ),
        transformation_type = "batch_correction"
    )
}
