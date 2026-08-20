# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.ARTIFACT_ROW_SELECTION_HINT_SCHEMA <-
    "multischolar.artifact_row_selection_hint"
.ARTIFACT_ROW_SELECTION_HINT_VERSION <- 1L
.ARTIFACT_ROW_SELECTION_SCHEMA <- "multischolar.artifact_row_selection"
.ARTIFACT_ROW_SELECTION_VERSION <- 1L
.ARTIFACT_ROW_SELECTION_RECIPE_REF <- ".selection_recipe"
.ARTIFACT_ROW_SELECTION_PATCH_REF <- ".selection_patch"

newArtifactRowSelectionHint <- function(
    slot_name,
    key_columns,
    method,
    normalized_parameters,
    software,
    lineage,
    rejected_reasons = character(),
    compaction = list(
        enabled = FALSE,
        reason = "representative_measurement_not_available"
    )
) {
    hint <- list(
        schema = .ARTIFACT_ROW_SELECTION_HINT_SCHEMA,
        schema_version = .ARTIFACT_ROW_SELECTION_HINT_VERSION,
        representation = "row_selection",
        slot_name = slot_name,
        key_columns = key_columns,
        method = method,
        normalized_parameters = normalized_parameters,
        software = software,
        lineage = lineage,
        rejected_reasons = rejected_reasons,
        compaction = compaction
    )
    artifactWorkflowStateValidateSelectionHint(hint)
}

artifactWorkflowStateValidateSelectionHint <- function(hint) {
    required <- c(
        "schema", "schema_version", "representation", "slot_name",
        "key_columns", "method", "normalized_parameters", "software",
        "lineage", "rejected_reasons", "compaction"
    )
    valid <- is.list(hint) && identical(names(hint), required) &&
        identical(hint$schema, .ARTIFACT_ROW_SELECTION_HINT_SCHEMA) &&
        identical(
            as.integer(hint$schema_version),
            .ARTIFACT_ROW_SELECTION_HINT_VERSION
        ) && identical(hint$representation, "row_selection") &&
        workflowStateScalarString(hint$slot_name) &&
        is.character(hint$key_columns) && length(hint$key_columns) > 0L &&
        !anyNA(hint$key_columns) && !anyDuplicated(hint$key_columns) &&
        all(nzchar(hint$key_columns)) &&
        workflowStateScalarString(hint$method) &&
        is.list(hint$normalized_parameters) &&
        artifactWorkflowStateValidSelectionSoftware(hint$software) &&
        artifactWorkflowStateValidSelectionLineage(hint$lineage) &&
        artifactWorkflowStateValidRejectionReasons(hint$rejected_reasons) &&
        is.list(hint$compaction) && identical(hint$compaction$enabled, FALSE) &&
        workflowStateScalarString(hint$compaction$reason)
    if (!isTRUE(valid)) {
        artifactWorkflowStateAbort(
            "artifact row-selection hint is malformed",
            "multischolar_invalid_artifact_row_selection"
        )
    }
    artifactWorkflowStateAssertSafeMetadata(
        hint[setdiff(names(hint), "rejected_reasons")],
        "row-selection hint"
    )
    hint$schema_version <- as.integer(hint$schema_version)
    hint
}

artifactWorkflowStateValidSelectionSoftware <- function(software) {
    required <- c("name", "version", "source")
    is.list(software) && all(required %in% names(software)) &&
        all(vapply(software[required], workflowStateScalarString, logical(1)))
}

artifactWorkflowStateValidSelectionLineage <- function(lineage) {
    required <- c(
        "audit_record_id", "state_name", "parent_state", "parent_record_id"
    )
    is.list(lineage) && all(required %in% names(lineage)) &&
        all(vapply(lineage[required], workflowStateScalarString, logical(1)))
}

artifactWorkflowStateValidRejectionReasons <- function(rejected_reasons) {
    if (!is.character(rejected_reasons) || anyNA(rejected_reasons) ||
        any(!nzchar(rejected_reasons))) {
        return(FALSE)
    }
    if (length(rejected_reasons) == 0L) return(TRUE)
    reason_names <- names(rejected_reasons)
    !is.null(reason_names) && !anyNA(reason_names) &&
        all(nzchar(reason_names)) && !anyDuplicated(reason_names)
}

artifactWorkflowStateKeyColumn <- function(value, owner) {
    if (anyNA(value)) {
        artifactWorkflowStateAbort(
            sprintf("artifact selection key column '%s' contains missing values", owner),
            "multischolar_invalid_artifact_row_selection",
            key_column = owner
        )
    }
    if (is.factor(value)) value <- as.character(value)
    if (inherits(value, "Date")) {
        return(paste0("Date:", sprintf("%a", unclass(value))))
    }
    if (inherits(value, "POSIXt")) {
        return(paste0("POSIXct:", sprintf("%a", unclass(value))))
    }
    if (inherits(value, "integer64")) {
        return(paste0("integer64:", as.character(value)))
    }
    if (is.character(value)) {
        return(paste0("character:", nchar(value, type = "bytes"), ":", value))
    }
    if (is.integer(value)) return(paste0("integer:", value))
    if (is.double(value)) return(paste0("double:", sprintf("%a", value)))
    if (is.logical(value)) return(paste0("logical:", value))
    artifactWorkflowStateAbort(
        sprintf("artifact selection key column '%s' has an unsupported type", owner),
        "multischolar_invalid_artifact_row_selection",
        key_column = owner
    )
}

artifactWorkflowStateEntityKeys <- function(data, key_columns) {
    if (!is.data.frame(data) || !all(key_columns %in% names(data))) {
        artifactWorkflowStateAbort(
            "artifact selection key columns are absent from the scientific table",
            "multischolar_invalid_artifact_row_selection",
            key_columns = key_columns
        )
    }
    if (nrow(data) == 0L) return(character())
    components <- Map(
        artifactWorkflowStateKeyColumn,
        data[key_columns],
        key_columns
    )
    components <- Map(
        function(name, value) paste0(nchar(name, type = "bytes"), ":", name, "=", value),
        key_columns,
        components
    )
    keys <- do.call(paste, c(components, sep = "\u001f"))
    if (anyDuplicated(keys) > 0L) {
        artifactWorkflowStateAbort(
            "artifact selection keys are not unique",
            "multischolar_ambiguous_artifact_row_selection",
            key_columns = key_columns
        )
    }
    keys
}

artifactWorkflowStateSelectionPayloadKey <- function(metadata, slot_name) {
    node <- metadata$slot_values[[slot_name]]
    if (!is.list(node) || !identical(node$node_type, "rectangular") ||
        !workflowStateScalarString(node$payload_key)) {
        artifactWorkflowStateAbort(
            sprintf("artifact selection slot '%s' is not one rectangular payload", slot_name),
            "multischolar_invalid_artifact_row_selection",
            slot_name = slot_name
        )
    }
    node$payload_key
}

artifactWorkflowStateSelectionReason <- function(keys, rejected_reasons) {
    reasons <- rep(NA_character_, length(keys))
    if (length(rejected_reasons) == 0L) return(reasons)
    if (is.null(names(rejected_reasons)) || anyDuplicated(names(rejected_reasons))) {
        artifactWorkflowStateAbort(
            "artifact selection rejection reasons require unique entity-key names",
            "multischolar_invalid_artifact_row_selection"
        )
    }
    matched <- match(keys, names(rejected_reasons))
    present <- !is.na(matched)
    reasons[present] <- unname(rejected_reasons[matched[present]])
    reasons
}

artifactWorkflowStateValidateRejectedReasons <- function(
    parent_keys,
    child_keys,
    rejected_reasons
) {
    rejected_keys <- parent_keys[!parent_keys %in% child_keys]
    if (length(rejected_keys) == 0L && length(rejected_reasons) == 0L) {
        return(invisible(TRUE))
    }
    if (!identical(sort(names(rejected_reasons)), sort(rejected_keys))) {
        artifactWorkflowStateAbort(
            "artifact selection rejection reasons do not match rejected entities",
            "multischolar_invalid_artifact_row_selection"
        )
    }
    invisible(TRUE)
}

artifactWorkflowStateSelectionPlan <- function(parent_object, state_object, hint) {
    hint <- artifactWorkflowStateValidateSelectionHint(hint)
    if (!isS4(parent_object) || !isS4(state_object) ||
        !identical(class(parent_object), class(state_object)) ||
        !hint$slot_name %in% methods::slotNames(parent_object)) {
        artifactWorkflowStateAbort(
            "artifact selection requires matching S4 parent and child objects",
            "multischolar_invalid_artifact_row_selection"
        )
    }
    parent <- methods::slot(parent_object, hint$slot_name)
    child <- methods::slot(state_object, hint$slot_name)
    if (!is.data.frame(parent) || !is.data.frame(child)) {
        artifactWorkflowStateAbort(
            "artifact selection requires data-frame parent and child slots",
            "multischolar_invalid_artifact_row_selection"
        )
    }
    parent_keys <- artifactWorkflowStateEntityKeys(parent, hint$key_columns)
    child_keys <- artifactWorkflowStateEntityKeys(child, hint$key_columns)
    selected_parent <- match(child_keys, parent_keys)
    if (anyNA(selected_parent)) {
        artifactWorkflowStateAbort(
            "artifact selection child contains an entity absent from its parent",
            "multischolar_invalid_artifact_row_selection"
        )
    }
    removed_columns <- setdiff(names(parent), names(child))
    if (length(removed_columns) > 0L) {
        artifactWorkflowStateAbort(
            "artifact row selection removed existing scientific columns",
            "multischolar_artifact_selection_requires_materialization",
            removed_columns = removed_columns
        )
    }
    common_columns <- intersect(names(parent), names(child))
    changed_columns <- common_columns[!vapply(common_columns, function(column) {
        identical(child[[column]], parent[[column]][selected_parent])
    }, logical(1))]
    if (length(changed_columns) > 0L) {
        artifactWorkflowStateAbort(
            "artifact row selection changed existing scientific columns",
            "multischolar_artifact_selection_requires_materialization",
            changed_columns = changed_columns
        )
    }
    added_columns <- setdiff(names(child), names(parent))
    if ("entity_key" %in% added_columns) {
        artifactWorkflowStateAbort(
            "artifact row selection added the reserved recipe key column",
            "multischolar_artifact_selection_requires_materialization"
        )
    }
    artifactWorkflowStateValidateRejectedReasons(
        parent_keys,
        child_keys,
        hint$rejected_reasons
    )
    selected_order <- match(parent_keys, child_keys)
    selected <- !is.na(selected_order)
    recipe <- data.frame(
        entity_key = parent_keys,
        selection_status = ifelse(selected, "selected", "rejected"),
        parent_order = seq_along(parent_keys),
        selected_order = as.integer(selected_order),
        stringsAsFactors = FALSE
    )
    recipe <- cbind(
        recipe,
        parent[hint$key_columns],
        stringsAsFactors = FALSE
    )
    recipe$failure_reason <- artifactWorkflowStateSelectionReason(
        parent_keys,
        hint$rejected_reasons
    )
    patch <- if (length(added_columns) == 0L) {
        NULL
    } else {
        data.frame(
            entity_key = child_keys,
            child[added_columns],
            check.names = FALSE,
            stringsAsFactors = FALSE
        )
    }
    list(
        parent = parent,
        child = child,
        parent_keys = parent_keys,
        child_keys = child_keys,
        recipe = recipe,
        patch = patch,
        added_columns = added_columns,
        removed_columns = character(),
        selected_count = as.integer(length(child_keys)),
        rejected_count = as.integer(length(parent_keys) - length(child_keys))
    )
}

artifactWorkflowStateReusableReference <- function(
    previous_manifest,
    payload_key,
    payload_metadata
) {
    if (is.null(previous_manifest) || is.null(previous_manifest$data$metadata_json)) {
        return(NULL)
    }
    previous_metadata <- artifactWorkflowStateUnserializeMetadata(
        previous_manifest$data$metadata_json,
        "parent S4 codec metadata"
    )
    previous_payload <- previous_metadata$payloads[[payload_key]]
    previous_ref <- previous_manifest$data$artifact_refs[[payload_key]]
    if (is.null(previous_payload) || is.null(previous_ref) ||
        !identical(previous_payload, payload_metadata)) {
        return(NULL)
    }
    previous_ref
}

artifactWorkflowStateWriteEncodedPayload <- function(
    store,
    identity,
    generation_id,
    encoded,
    state_role,
    failure_injector = NULL
) {
    artifactStoreWriteParquet(
        store,
        encoded,
        logical_key = list(
            project_id = identity$project_id,
            omic_type = identity$omic_type,
            workflow_slug = identity$workflow_slug,
            stage_id = paste0("state_", substr(generation_id, 5L, 20L)),
            state_role = state_role,
            generation_id = generation_id
        ),
        failure_injector = failure_injector
    )
}

artifactWorkflowStateWriteSelectionData <- function(
    store,
    identity,
    generation_id,
    state_object,
    parent_object,
    previous_manifest,
    bundle,
    hint,
    failure_injector = NULL
) {
    plan <- artifactWorkflowStateSelectionPlan(
        parent_object,
        state_object,
        hint
    )
    target_key <- artifactWorkflowStateSelectionPayloadKey(
        bundle$metadata,
        hint$slot_name
    )
    references <- list()
    new_reference_names <- character()
    for (payload_key in setdiff(names(bundle$payloads), target_key)) {
        reusable <- artifactWorkflowStateReusableReference(
            previous_manifest,
            payload_key,
            bundle$metadata$payloads[[payload_key]]
        )
        if (!is.null(reusable)) {
            references[[payload_key]] <- reusable
            next
        }
        encoded <- structure(
            list(
                payload = bundle$payloads[[payload_key]],
                metadata = bundle$metadata$payloads[[payload_key]]
            ),
            class = c("MultiScholaRArtifactRectangular", "list")
        )
        references[[payload_key]] <- artifactWorkflowStateWriteEncodedPayload(
            store,
            identity,
            generation_id,
            encoded,
            sprintf("payload_%04d", match(payload_key, names(bundle$payloads))),
            failure_injector
        )
        new_reference_names <- c(new_reference_names, payload_key)
    }
    recipe_encoded <- encodeArtifactTable(
        plan$recipe,
        stable_key = "entity_key",
        owner = paste(identity$omic_type, hint$method, "selection_recipe", sep = ".")
    )
    references[[.ARTIFACT_ROW_SELECTION_RECIPE_REF]] <-
        artifactWorkflowStateWriteEncodedPayload(
            store,
            identity,
            generation_id,
            recipe_encoded,
            "selection_recipe",
            failure_injector
        )
    new_reference_names <- c(
        new_reference_names,
        .ARTIFACT_ROW_SELECTION_RECIPE_REF
    )
    patch_ref_name <- NULL
    if (!is.null(plan$patch)) {
        patch_encoded <- encodeArtifactTable(
            plan$patch,
            stable_key = "entity_key",
            owner = paste(identity$omic_type, hint$method, "selection_patch", sep = ".")
        )
        references[[.ARTIFACT_ROW_SELECTION_PATCH_REF]] <-
            artifactWorkflowStateWriteEncodedPayload(
                store,
                identity,
                generation_id,
                patch_encoded,
                "selection_patch",
                failure_injector
            )
        new_reference_names <- c(
            new_reference_names,
            .ARTIFACT_ROW_SELECTION_PATCH_REF
        )
        patch_ref_name <- .ARTIFACT_ROW_SELECTION_PATCH_REF
    }
    parent_artifact_ids <- vapply(
        previous_manifest$data$artifact_refs,
        `[[`,
        character(1),
        "artifact_id"
    )
    recipe_ref <- references[[.ARTIFACT_ROW_SELECTION_RECIPE_REF]]
    dependencies <- lapply(seq_along(parent_artifact_ids), function(index) {
        list(
            artifact_id = recipe_ref$artifact_id,
            depends_on_artifact_id = unname(parent_artifact_ids[[index]]),
            relationship_type = "selection_parent_generation",
            ordinal = as.integer(index - 1L)
        )
    })
    lineage <- hint$lineage
    lineage$parent_generation_id <- previous_manifest$generation_id
    lineage$parent_content_identity <- previous_manifest$data$semantic_digest
    bundle$metadata$derivation <- list(
        schema = .ARTIFACT_ROW_SELECTION_SCHEMA,
        schema_version = .ARTIFACT_ROW_SELECTION_VERSION,
        representation = "row_selection",
        target_slot = hint$slot_name,
        target_payload_key = target_key,
        parent_generation_id = previous_manifest$generation_id,
        parent_semantic_digest = previous_manifest$data$semantic_digest,
        parent_artifact_ids = unname(parent_artifact_ids),
        key_columns = hint$key_columns,
        method = hint$method,
        normalized_parameters = hint$normalized_parameters,
        software = hint$software,
        lineage = lineage,
        selected_count = plan$selected_count,
        rejected_count = plan$rejected_count,
        selected_order_digest = artifactSemanticDigest(plan$child_keys),
        rejected_order_digest = artifactSemanticDigest(
            plan$parent_keys[!plan$parent_keys %in% plan$child_keys]
        ),
        recipe_ref_name = .ARTIFACT_ROW_SELECTION_RECIPE_REF,
        patch_ref_name = patch_ref_name,
        patched_columns = plan$added_columns,
        removed_columns = plan$removed_columns,
        compaction = hint$compaction
    )
    list(
        s4_class = bundle$metadata$class_name,
        semantic_digest = bundle$metadata$semantic_digest,
        codec = bundle$metadata$codec,
        metadata_json = artifactWorkflowStateSerializeMetadata(
            bundle$metadata,
            "S4 selection codec metadata"
        ),
        artifact_refs = references,
        reused = FALSE,
        new_reference_names = new_reference_names,
        dependencies = dependencies
    )
}

artifactWorkflowStateValidateDerivation <- function(derivation, manifest) {
    required <- c(
        "schema", "schema_version", "representation", "target_slot",
        "target_payload_key", "parent_generation_id",
        "parent_semantic_digest", "parent_artifact_ids", "key_columns",
        "method", "normalized_parameters", "software", "lineage",
        "selected_count", "rejected_count", "selected_order_digest",
        "rejected_order_digest", "recipe_ref_name", "patch_ref_name",
        "patched_columns", "removed_columns", "compaction"
    )
    valid <- is.list(derivation) && identical(names(derivation), required) &&
        identical(derivation$schema, .ARTIFACT_ROW_SELECTION_SCHEMA) &&
        identical(
            as.integer(derivation$schema_version),
            .ARTIFACT_ROW_SELECTION_VERSION
        ) && identical(derivation$representation, "row_selection") &&
        workflowStateScalarString(derivation$target_slot) &&
        workflowStateScalarString(derivation$target_payload_key) &&
        identical(derivation$parent_generation_id, manifest$parent_generation_id) &&
        is.character(derivation$parent_artifact_ids) &&
        length(derivation$parent_artifact_ids) > 0L &&
        !anyNA(derivation$parent_artifact_ids) &&
        !anyDuplicated(derivation$parent_artifact_ids) &&
        identical(
            derivation$recipe_ref_name,
            .ARTIFACT_ROW_SELECTION_RECIPE_REF
        ) &&
        (is.null(derivation$patch_ref_name) || identical(
            derivation$patch_ref_name,
            .ARTIFACT_ROW_SELECTION_PATCH_REF
        )) &&
        is.character(derivation$key_columns) && length(derivation$key_columns) > 0L &&
        !anyNA(derivation$key_columns) && all(nzchar(derivation$key_columns)) &&
        !anyDuplicated(derivation$key_columns) &&
        workflowStateScalarString(derivation$method) &&
        is.list(derivation$normalized_parameters) &&
        artifactWorkflowStateValidSelectionSoftware(derivation$software) &&
        artifactWorkflowStateValidSelectionLineage(derivation$lineage) &&
        artifactWorkflowStateCount(derivation$selected_count) &&
        artifactWorkflowStateCount(derivation$rejected_count) &&
        artifactWorkflowStateOptionalScalarString(derivation$patch_ref_name) &&
        is.character(derivation$patched_columns) &&
        !anyNA(derivation$patched_columns) && all(nzchar(derivation$patched_columns)) &&
        !anyDuplicated(derivation$patched_columns) &&
        identical(derivation$removed_columns, character()) &&
        is.list(derivation$compaction) &&
        identical(derivation$compaction$enabled, FALSE) &&
        workflowStateScalarString(derivation$compaction$reason)
    if (!isTRUE(valid)) {
        artifactWorkflowStateAbort(
            "artifact row-selection derivation is malformed",
            "multischolar_invalid_artifact_row_selection"
        )
    }
    artifactRefValidateDigest(
        derivation$parent_semantic_digest,
        "parent_semantic_digest"
    )
    artifactRefValidateDigest(
        derivation$selected_order_digest,
        "selected_order_digest"
    )
    artifactRefValidateDigest(
        derivation$rejected_order_digest,
        "rejected_order_digest"
    )
    vapply(
        derivation$parent_artifact_ids,
        artifactRefValidateId,
        character(1),
        field = "parent_artifact_ids",
        prefix = "art"
    )
    artifactWorkflowStateAssertSafeMetadata(derivation, "row-selection derivation")
    derivation
}

artifactWorkflowStateCount <- function(value) {
    is.integer(value) && length(value) == 1L && !is.na(value) && value >= 0L
}

artifactWorkflowStateOptionalScalarString <- function(value) {
    is.null(value) || workflowStateScalarString(value)
}

artifactWorkflowStateApplyExpectedTableMetadata <- function(value, metadata) {
    value <- value[metadata$logical_names]
    if (identical(metadata$row_names$kind, "automatic")) {
        attr(value, "row.names") <- c(NA_integer_, -nrow(value))
    } else {
        attr(value, "row.names") <- metadata$row_names$value
    }
    class(value) <- metadata$data_frame_class
    value
}

artifactWorkflowStateHydrateSelectionPayloads <- function(
    store,
    manifest,
    metadata,
    payload_entries,
    hydrate_parent_fn
) {
    derivation <- artifactWorkflowStateValidateDerivation(
        metadata$derivation,
        manifest
    )
    target_key <- artifactWorkflowStateSelectionPayloadKey(
        metadata,
        derivation$target_slot
    )
    if (!identical(target_key, derivation$target_payload_key)) {
        artifactWorkflowStateAbort(
            "artifact row-selection target differs from its S4 slot metadata",
            "multischolar_invalid_artifact_row_selection"
        )
    }
    recipe_entry <- payload_entries[[derivation$recipe_ref_name]]
    if (is.null(recipe_entry)) {
        artifactWorkflowStateAbort(
            "artifact row-selection recipe is missing",
            "multischolar_invalid_artifact_row_selection"
        )
    }
    parent_manifest <- artifactWorkflowStateReadManifest(
        store,
        artifactWorkflowStateManifestRelativePath(
            store,
            derivation$parent_generation_id
        )
    )
    if (!identical(
        parent_manifest$data$semantic_digest,
        derivation$parent_semantic_digest
    )) {
        artifactWorkflowStateAbort(
            "artifact selection parent content identity changed",
            "multischolar_artifact_selection_parent_mismatch"
        )
    }
    actual_parent_ids <- unname(vapply(
        parent_manifest$data$artifact_refs,
        `[[`,
        character(1),
        "artifact_id"
    ))
    if (!identical(actual_parent_ids, unname(derivation$parent_artifact_ids))) {
        artifactWorkflowStateAbort(
            "artifact selection parent refs differ from the immutable recipe",
            "multischolar_artifact_selection_parent_mismatch"
        )
    }
    parent_object <- hydrate_parent_fn(parent_manifest)
    parent <- methods::slot(parent_object, derivation$target_slot)
    parent_keys <- artifactWorkflowStateEntityKeys(
        parent,
        derivation$key_columns
    )
    recipe <- decodeArtifactRectangular(
        recipe_entry$payload,
        recipe_entry$metadata
    )
    required_recipe <- c(
        "entity_key", "selection_status", "parent_order", "selected_order",
        derivation$key_columns, "failure_reason"
    )
    valid_recipe <- identical(names(recipe), required_recipe) &&
        identical(recipe$entity_key, parent_keys) &&
        identical(as.integer(recipe$parent_order), seq_along(parent_keys)) &&
        all(recipe$selection_status %in% c("selected", "rejected")) &&
        all(is.na(recipe$failure_reason[recipe$selection_status == "selected"])) &&
        all(!is.na(recipe$failure_reason[recipe$selection_status == "rejected"]))
    if (!isTRUE(valid_recipe)) {
        artifactWorkflowStateAbort(
            "artifact row-selection recipe does not match its parent",
            "multischolar_artifact_selection_parent_mismatch"
        )
    }
    selected <- recipe$selection_status == "selected"
    selected_order <- recipe$selected_order[selected]
    valid_order <- identical(
        sort(as.integer(selected_order)),
        seq_len(sum(selected))
    ) && all(is.na(recipe$selected_order[!selected]))
    selected_keys <- recipe$entity_key[selected][order(selected_order)]
    rejected_keys <- recipe$entity_key[!selected]
    valid_summary <- identical(sum(selected), derivation$selected_count) &&
        identical(sum(!selected), derivation$rejected_count) &&
        identical(
            artifactSemanticDigest(selected_keys),
            derivation$selected_order_digest
        ) && identical(
            artifactSemanticDigest(rejected_keys),
            derivation$rejected_order_digest
        )
    if (!isTRUE(valid_order) || !isTRUE(valid_summary)) {
        artifactWorkflowStateAbort(
            "artifact row-selection counts or ordering are invalid",
            "multischolar_invalid_artifact_row_selection"
        )
    }
    selected_parent <- match(selected_keys, parent_keys)
    value <- parent[selected_parent, , drop = FALSE]
    patch <- NULL
    if (!is.null(derivation$patch_ref_name)) {
        patch_entry <- payload_entries[[derivation$patch_ref_name]]
        if (is.null(patch_entry)) {
            artifactWorkflowStateAbort(
                "artifact row-selection patch is missing",
                "multischolar_invalid_artifact_row_selection"
            )
        }
        patch <- decodeArtifactRectangular(
            patch_entry$payload,
            patch_entry$metadata
        )
        if (!identical(patch$entity_key, selected_keys) ||
            !identical(setdiff(names(patch), "entity_key"), derivation$patched_columns)) {
            artifactWorkflowStateAbort(
                "artifact row-selection patch does not match selected entities",
                "multischolar_invalid_artifact_row_selection"
            )
        }
        for (column in derivation$patched_columns) {
            value[[column]] <- patch[[column]]
        }
    } else if (length(derivation$patched_columns) > 0L) {
        artifactWorkflowStateAbort(
            "artifact row-selection metadata requires an absent patch",
            "multischolar_invalid_artifact_row_selection"
        )
    }
    expected <- metadata$payloads[[derivation$target_payload_key]]
    value <- artifactWorkflowStateApplyExpectedTableMetadata(value, expected)
    encoded <- encodeArtifactTable(value, owner = expected$owner)
    if (!identical(encoded$metadata, expected)) {
        artifactWorkflowStateAbort(
            "artifact row-selection hydration changed the expected table",
            "multischolar_inexact_artifact_state_hydration"
        )
    }
    payloads <- lapply(names(metadata$payloads), function(payload_key) {
        if (identical(payload_key, derivation$target_payload_key)) {
            return(encoded$payload)
        }
        entry <- payload_entries[[payload_key]]
        if (is.null(entry)) {
            artifactWorkflowStateAbort(
                "artifact row-selection state is missing a required payload",
                "multischolar_invalid_artifact_row_selection",
                payload_key = payload_key
            )
        }
        entry$payload
    })
    names(payloads) <- names(metadata$payloads)
    payloads
}
