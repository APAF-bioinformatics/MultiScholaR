# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.METAB_NORM_AUDIT_SCHEMA <- "multischolar.metabolomics_normalization_audit"
.METAB_NORM_AUDIT_VERSION <- 1L
.METAB_NORM_STAGES <- c(
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

metabNormArtifactMode <- function() {
    mode <- getOption("multischolar.metabolomics.norm_artifacts", "enabled")
    match.arg(mode, c("enabled", "disabled"))
}

metabNormSafeDigest <- function(value) {
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

metabNormArtifactActive <- function(workflow_data) {
    manager <- workflow_data$state_manager
    inherits(manager, "ArtifactWorkflowState") || (
        !identical(metabNormArtifactMode(), "disabled") &&
        metabQcWorkflowData(workflow_data) &&
        metabArtifactCoordinatorOwned(workflow_data)
    )
}

storeMetabNormIntermediate <- function(
    workflow_data,
    norm_data,
    field,
    object,
    state_name
) {
    if (!isTRUE(metabNormArtifactActive(workflow_data))) {
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

resolveMetabNormIntermediate <- function(workflow_data, norm_data, field) {
    value <- norm_data[[field]]
    if (!is.null(value)) return(value)
    refs <- norm_data$artifact_state_refs
    state_name <- if (is.list(refs)) refs[[field]] else NULL
    if (!workflowStateScalarString(state_name)) return(NULL)
    workflow_data$state_manager$getState(state_name)
}

metabNormArgsSummary <- function(object) {
    args <- object@args
    list(
        names = sort(names(args), method = "radix"),
        semantic_digest = metabNormSafeDigest(args),
        object_bytes = unname(as.numeric(utils::object.size(args))),
        storage = "exact_s4_codec_payload_refs",
        compatible_getter = "ArtifactWorkflowState$getState"
    )
}

metabNormBuildAuditRecord <- function(
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
    if (!methods::is(before, "MetaboliteAssayData") ||
        !methods::is(after, "MetaboliteAssayData") ||
        !stage_id %in% .METAB_NORM_STAGES) {
        return(NULL)
    }
    substantive <- list(
        schema = .METAB_NORM_AUDIT_SCHEMA,
        schema_version = .METAB_NORM_AUDIT_VERSION,
        stage_id = stage_id,
        status = status,
        transformation_type = transformation_type,
        parent_state = parent_state,
        current_state = current_state,
        parent_generation_id = parent_generation_id,
        parent_record_id = parent_record_id,
        parameters = artifactCanonicalizeValue(parameters),
        software = metabQcSoftware(stage_id),
        before_assays = metabQcAssaySummary(before),
        after_assays = metabQcAssaySummary(after),
        before_design_digest = artifactSemanticDigest(before@design_matrix),
        after_design_digest = artifactSemanticDigest(after@design_matrix),
        after_args = metabNormArgsSummary(after),
        assay_order = names(after@metabolite_data)
    )
    canonical_digest <- artifactSemanticDigest(substantive)
    c(substantive, list(
        record_id = paste0("metab-norm:", substr(canonical_digest, 1L, 24L)),
        timestamp_utc = format(now, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC"),
        canonical_digest = canonical_digest
    ))
}

metabNormControlSummary <- function(best_k_per_assay, ctrl_per_assay) {
    assay_names <- union(names(best_k_per_assay), names(ctrl_per_assay))
    lapply(assay_names, function(assay_name) {
        controls <- ctrl_per_assay[[assay_name]]
        list(
            assay_name = assay_name,
            k = best_k_per_assay[[assay_name]],
            control_count = if (is.null(controls)) {
                0L
            } else {
                as.integer(sum(controls, na.rm = TRUE))
            },
            control_length = if (is.null(controls)) 0L else length(controls),
            control_digest = metabNormSafeDigest(controls)
        )
    })
}

metabNormOptimizerSummary <- function(optimization_state) {
    if (is.null(optimization_state)) {
        return(list(status = "not_recorded", assays = list()))
    }
    results <- optimization_state$ruvResults
    if (!is.list(results)) results <- list()
    list(
        status = "recorded",
        parameters = optimization_state$ruvParams,
        result_digest = metabNormSafeDigest(results),
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

metabNormFeatureSelectionSummary <- function(feature_ids) {
    if (is.null(feature_ids)) {
        return(list(source = "object_internal_standard_regex", assays = list()))
    }
    values <- if (is.list(feature_ids)) feature_ids else list(all = feature_ids)
    list(
        source = "manual_assay_selection",
        assays = lapply(names(values), function(assay_name) {
            ids <- as.character(values[[assay_name]])
            list(
                assay_name = assay_name,
                count = as.integer(length(ids)),
                identity_digest = metabNormSafeDigest(ids)
            )
        })
    )
}

metabNormCorrelationSummary <- function(results) {
    values <- if (is.list(results)) results else list(all = results)
    list(
        assays = lapply(names(values), function(assay_name) {
            value <- values[[assay_name]]
            list(
                assay_name = assay_name,
                result_count = if (is.data.frame(value)) nrow(value) else length(value),
                semantic_digest = metabNormSafeDigest(value)
            )
        })
    )
}

metabNormStateReferenceHint <- function(
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
        software = metabQcSoftware(stage_id)[c("name", "version", "source")],
        lineage = list(
            audit_record_id = record$record_id,
            state_name = state_name,
            parent_state = parent_state,
            parent_record_id = record$parent_record_id
        )
    )
}

metabNormArtifactEligible <- function(
    workflow_data,
    state_manager,
    before,
    after
) {
    valid <- methods::is(before, "MetaboliteAssayData") &&
        methods::is(after, "MetaboliteAssayData")
    if (!isTRUE(valid)) return(FALSE)
    if (inherits(state_manager, "ArtifactWorkflowState")) {
        return(identical(workflowStateType(state_manager), "metabolomics_standard"))
    }
    !identical(metabNormArtifactMode(), "disabled") &&
        metabQcWorkflowData(workflow_data) &&
        metabArtifactCoordinatorOwned(workflow_data)
}

saveMetabNormState <- function(
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
    parent_state <- metabQcCurrentStateName(state_manager)
    enabled <- metabNormArtifactEligible(
        workflow_data,
        state_manager,
        before,
        after
    )
    resources <- if (isTRUE(enabled)) {
        metabQcArtifactManager(workflow_data, state_manager)
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
        metabQcValidateParent(manager, before)
        artifact_parent <- manager$getCurrentStateName()
        if (isTRUE(resources$owned) && !identical(artifact_parent, parent_state)) {
            metabArtifactAbort(
                "metabolomics normalization managers disagree on the active parent",
                "multischolar_metabolomics_norm_parent_mismatch"
            )
        }
        manager$getCurrentGenerationId()
    }
    parent_audit <- metabQcStateAuditRecord(
        if (is.null(manager)) state_manager else manager
    )
    parent_record_id <- if (is.null(parent_audit)) {
        "legacy_or_qc_checkpoint"
    } else {
        parent_audit$record_id
    }
    record <- metabNormBuildAuditRecord(
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
        metabQcSaveMemory(
            state_manager,
            state_name,
            after,
            config_object,
            description,
            record
        )
        return(after)
    }
    hint <- metabNormStateReferenceHint(
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
        config_object = metabQcSafeConfig(
            config_object,
            metabQcSoftware(stage_id)
        ),
        description = description,
        audit_metadata = record,
        persistence_hint = hint,
        failure_injector = failure_injector,
        action_id = artifactOpaqueId("action"),
        expected_parent_generation_id = parent_generation_id
    )
    if (!identical(result$status, "accepted")) {
        metabArtifactAbort(
            "metabolomics normalization commit did not advance its exact parent",
            "multischolar_metabolomics_norm_commit_rejected",
            result = result
        )
    }
    if (isTRUE(resources$owned)) {
        memory_error <- tryCatch({
            metabQcSaveMemory(
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
                metabArtifactAbort(
                    paste(
                        "metabolomics normalization memory save and artifact",
                        "rollback both failed"
                    ),
                    "multischolar_metabolomics_norm_divergent_state",
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
        metabArtifactAbort(
            "metabolomics normalization artifact changed scientific state",
            "multischolar_inexact_metabolomics_norm_hydration"
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

recordMetabNormNoOp <- function(
    workflow_data,
    current_s4,
    stage_id,
    parameters,
    status = "skipped"
) {
    writable <- metabQcWorkflowData(workflow_data) &&
        !is.null(workflow_data$state_manager) &&
        is.function(workflow_data$state_manager$saveState)
    if (!isTRUE(writable)) return(current_s4)
    state_name <- metabQcCurrentStateName(workflow_data$state_manager)
    saveMetabNormState(
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

persistMetabNormItsdState <- function(
    workflow_data,
    norm_data,
    before,
    after,
    aggregation,
    feature_ids,
    description
) {
    after <- saveMetabNormState(
        workflow_data,
        before,
        after,
        "itsd_normalization",
        "metab_itsd_norm",
        workflow_data$config_list,
        description,
        parameters = list(
            method = "ITSD",
            aggregation = aggregation,
            feature_selection = metabNormFeatureSelectionSummary(feature_ids),
            internal_standard_regex = before@internal_standard_regex
        ),
        transformation_type = "normalization"
    )
    storeMetabNormIntermediate(
        workflow_data,
        norm_data,
        "post_itsd_obj",
        after,
        "metab_itsd_norm"
    )
    after
}

persistMetabNormLog2State <- function(
    workflow_data,
    norm_data,
    before,
    after,
    offset,
    description
) {
    after <- saveMetabNormState(
        workflow_data,
        before,
        after,
        "log2_transformation",
        "metab_log2",
        workflow_data$config_list,
        description,
        parameters = list(
            transform = "log2",
            offset = offset,
            non_finite_policy = "preserve_as_NA"
        ),
        transformation_type = "transformation"
    )
    storeMetabNormIntermediate(
        workflow_data,
        norm_data,
        "post_log2_obj",
        after,
        "metab_log2"
    )
    after
}

persistMetabNormBetweenSampleState <- function(
    workflow_data,
    before,
    after,
    method,
    description
) {
    skipped <- identical(method, "none")
    saveMetabNormState(
        workflow_data,
        before,
        after,
        "between_sample_normalization",
        "metab_normalized",
        workflow_data$config_list,
        description,
        parameters = list(method = method),
        status = if (skipped) "skipped" else "applied",
        transformation_type = if (skipped) "no_op" else "normalization"
    )
}

persistMetabNormRuvState <- function(
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
    saveMetabNormState(
        workflow_data,
        before,
        after,
        "ruv_correction",
        "metab_ruv_corrected",
        workflow_data$config_list,
        description,
        parameters = list(
            mode = if (is.null(ruv_mode)) "not_recorded" else ruv_mode,
            grouping_variable = grouping_variable,
            assays = metabNormControlSummary(
                best_k_per_assay,
                ctrl_per_assay
            ),
            optimizer = metabNormOptimizerSummary(optimization_state)
        ),
        transformation_type = "batch_correction"
    )
}
