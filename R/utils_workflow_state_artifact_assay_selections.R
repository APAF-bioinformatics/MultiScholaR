# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.ARTIFACT_ASSAY_SELECTION_HINT_SCHEMA <-
    "multischolar.artifact_assay_selection_hint"
.ARTIFACT_ASSAY_SELECTION_HINT_VERSION <- 1L
.ARTIFACT_ASSAY_SELECTION_SCHEMA <-
    "multischolar.artifact_assay_selection"
.ARTIFACT_ASSAY_SELECTION_VERSION <- 1L
.ARTIFACT_ASSAY_SELECTION_RECIPE_PREFIX <- ".assay_selection_recipe_"
.ARTIFACT_ASSAY_SELECTION_PATCH_PREFIX <- ".assay_selection_patch_"

newArtifactAssaySelectionHint <- function(
    slot_name,
    assay_key_columns,
    method,
    normalized_parameters,
    software,
    lineage,
    rejected_reasons = list(),
    compaction = list(
        enabled = FALSE,
        reason = "representative_measurement_not_available"
    )
) {
    hint <- list(
        schema = .ARTIFACT_ASSAY_SELECTION_HINT_SCHEMA,
        schema_version = .ARTIFACT_ASSAY_SELECTION_HINT_VERSION,
        representation = "assay_selection",
        slot_name = slot_name,
        assay_key_columns = assay_key_columns,
        method = method,
        normalized_parameters = normalized_parameters,
        software = software,
        lineage = lineage,
        rejected_reasons = rejected_reasons,
        compaction = compaction
    )
    artifactWorkflowStateValidateAssaySelectionHint(hint)
}

artifactWorkflowStateValidateAssayKeyColumns <- function(value) {
    valid_names <- is.list(value) && length(value) > 0L &&
        !is.null(names(value)) && !anyNA(names(value)) &&
        all(nzchar(names(value))) && !anyDuplicated(names(value))
    valid_values <- isTRUE(valid_names) && all(vapply(value, function(columns) {
        is.character(columns) && length(columns) > 0L && !anyNA(columns) &&
            all(nzchar(columns)) && !anyDuplicated(columns)
    }, logical(1)))
    isTRUE(valid_values)
}

artifactWorkflowStateNormalizeAssayRejections <- function(
    rejected_reasons,
    assay_names
) {
    if (!is.list(rejected_reasons)) {
        artifactWorkflowStateAbort(
            "artifact assay-selection rejection reasons must be a list",
            "multischolar_invalid_artifact_assay_selection"
        )
    }
    if (length(rejected_reasons) == 0L) {
        rejected_reasons <- stats::setNames(
            rep(list(character()), length(assay_names)),
            assay_names
        )
    }
    valid <- identical(names(rejected_reasons), assay_names) &&
        all(vapply(
            rejected_reasons,
            artifactWorkflowStateValidRejectionReasons,
            logical(1)
        ))
    if (!isTRUE(valid)) {
        artifactWorkflowStateAbort(
            paste(
                "artifact assay-selection rejection reasons must match",
                "the ordered assay names"
            ),
            "multischolar_invalid_artifact_assay_selection"
        )
    }
    rejected_reasons
}

artifactWorkflowStateValidateAssaySelectionHint <- function(hint) {
    required <- c(
        "schema", "schema_version", "representation", "slot_name",
        "assay_key_columns", "method", "normalized_parameters", "software",
        "lineage", "rejected_reasons", "compaction"
    )
    valid <- is.list(hint) && identical(names(hint), required) &&
        identical(hint$schema, .ARTIFACT_ASSAY_SELECTION_HINT_SCHEMA) &&
        identical(
            as.integer(hint$schema_version),
            .ARTIFACT_ASSAY_SELECTION_HINT_VERSION
        ) && identical(hint$representation, "assay_selection") &&
        workflowStateScalarString(hint$slot_name) &&
        artifactWorkflowStateValidateAssayKeyColumns(
            hint$assay_key_columns
        ) && workflowStateScalarString(hint$method) &&
        is.list(hint$normalized_parameters) &&
        artifactWorkflowStateValidSelectionSoftware(hint$software) &&
        artifactWorkflowStateValidSelectionLineage(hint$lineage) &&
        is.list(hint$compaction) && identical(hint$compaction$enabled, FALSE) &&
        workflowStateScalarString(hint$compaction$reason)
    if (!isTRUE(valid)) {
        artifactWorkflowStateAbort(
            "artifact assay-selection hint is malformed",
            "multischolar_invalid_artifact_assay_selection"
        )
    }
    hint$rejected_reasons <- artifactWorkflowStateNormalizeAssayRejections(
        hint$rejected_reasons,
        names(hint$assay_key_columns)
    )
    artifactWorkflowStateAssertSafeMetadata(
        hint[setdiff(names(hint), "rejected_reasons")],
        "assay-selection hint"
    )
    hint$schema_version <- as.integer(hint$schema_version)
    hint
}

artifactWorkflowStateAssayPayloadKeys <- function(
    metadata,
    slot_name,
    assay_names = NULL
) {
    node <- metadata$slot_values[[slot_name]]
    valid <- is.list(node) && identical(node$node_type, "list") &&
        is.list(node$names) && is.character(node$names$values) &&
        identical(node$names$missing, integer()) && is.list(node$values) &&
        length(node$values) == length(node$names$values) &&
        all(vapply(node$values, function(value) {
            is.list(value) && identical(value$node_type, "rectangular") &&
                workflowStateScalarString(value$payload_key)
        }, logical(1)))
    if (!isTRUE(valid) || anyNA(node$names$values) ||
        any(!nzchar(node$names$values)) || anyDuplicated(node$names$values)) {
        artifactWorkflowStateAbort(
            sprintf(
                "artifact assay-selection slot '%s' is not a named table list",
                slot_name
            ),
            "multischolar_invalid_artifact_assay_selection",
            slot_name = slot_name
        )
    }
    names_in_metadata <- node$names$values
    if (!is.null(assay_names) && !identical(names_in_metadata, assay_names)) {
        artifactWorkflowStateAbort(
            "artifact assay-selection names differ from S4 codec metadata",
            "multischolar_invalid_artifact_assay_selection"
        )
    }
    keys <- vapply(node$values, `[[`, character(1), "payload_key")
    stats::setNames(keys, names_in_metadata)
}

artifactWorkflowStateTableSelectionPlan <- function(
    parent,
    child,
    key_columns,
    rejected_reasons
) {
    if (!is.data.frame(parent) || !is.data.frame(child)) {
        artifactWorkflowStateAbort(
            "artifact assay selection requires data-frame assay tables",
            "multischolar_invalid_artifact_assay_selection"
        )
    }
    parent_keys <- artifactWorkflowStateEntityKeys(parent, key_columns)
    child_keys <- artifactWorkflowStateEntityKeys(child, key_columns)
    selected_parent <- match(child_keys, parent_keys)
    if (anyNA(selected_parent)) {
        artifactWorkflowStateAbort(
            "artifact assay-selection child contains an entity absent from its parent",
            "multischolar_invalid_artifact_assay_selection"
        )
    }
    removed_columns <- setdiff(names(parent), names(child))
    if (length(removed_columns) > 0L) {
        artifactWorkflowStateAbort(
            "artifact assay selection removed scientific columns",
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
            "artifact assay selection changed scientific values",
            "multischolar_artifact_selection_requires_materialization",
            changed_columns = changed_columns
        )
    }
    added_columns <- setdiff(names(child), names(parent))
    if ("entity_key" %in% added_columns) {
        artifactWorkflowStateAbort(
            "artifact assay selection added the reserved recipe key column",
            "multischolar_artifact_selection_requires_materialization"
        )
    }
    artifactWorkflowStateValidateRejectedReasons(
        parent_keys,
        child_keys,
        rejected_reasons
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
        parent[key_columns],
        stringsAsFactors = FALSE
    )
    recipe$failure_reason <- artifactWorkflowStateSelectionReason(
        parent_keys,
        rejected_reasons
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
        parent_keys = parent_keys,
        child_keys = child_keys,
        recipe = recipe,
        patch = patch,
        added_columns = added_columns,
        selected_count = as.integer(length(child_keys)),
        rejected_count = as.integer(length(parent_keys) - length(child_keys))
    )
}

artifactWorkflowStateAssaySelectionPlan <- function(
    parent_object,
    state_object,
    hint
) {
    hint <- artifactWorkflowStateValidateAssaySelectionHint(hint)
    valid_objects <- isS4(parent_object) && isS4(state_object) &&
        identical(class(parent_object), class(state_object)) &&
        hint$slot_name %in% methods::slotNames(parent_object)
    if (!isTRUE(valid_objects)) {
        artifactWorkflowStateAbort(
            "artifact assay selection requires matching S4 parent and child objects",
            "multischolar_invalid_artifact_assay_selection"
        )
    }
    parent <- methods::slot(parent_object, hint$slot_name)
    child <- methods::slot(state_object, hint$slot_name)
    assay_names <- names(hint$assay_key_columns)
    valid_assays <- is.list(parent) && is.list(child) &&
        identical(names(parent), assay_names) && identical(names(child), assay_names)
    if (!isTRUE(valid_assays)) {
        artifactWorkflowStateAbort(
            "artifact assay selection requires the same ordered named assays",
            "multischolar_invalid_artifact_assay_selection"
        )
    }
    plans <- Map(
        artifactWorkflowStateTableSelectionPlan,
        parent,
        child,
        hint$assay_key_columns,
        hint$rejected_reasons
    )
    names(plans) <- assay_names
    plans
}

artifactWorkflowStateAssaySelectionRefName <- function(prefix, index) {
    paste0(prefix, sprintf("%04d", as.integer(index)))
}

artifactWorkflowStateWriteAssaySelectionReference <- function(
    store,
    identity,
    generation_id,
    value,
    owner,
    role,
    failure_injector
) {
    encoded <- encodeArtifactTable(value, stable_key = "entity_key", owner = owner)
    artifactWorkflowStateWriteEncodedPayload(
        store,
        identity,
        generation_id,
        encoded,
        role,
        failure_injector
    )
}

artifactWorkflowStateAssayParentReference <- function(
    manifest,
    metadata,
    slot_name,
    assay_name,
    payload_key
) {
    direct <- manifest$data$artifact_refs[[payload_key]]
    if (is.list(direct)) {
        return(list(ref_name = payload_key, ref = direct))
    }
    derivation <- metadata$derivation
    if (!is.list(derivation) ||
        !identical(derivation$schema, .ARTIFACT_ASSAY_SELECTION_SCHEMA) ||
        !identical(derivation$target_slot, slot_name)) {
        artifactWorkflowStateAbort(
            "artifact assay selection has no parent assay representation",
            "multischolar_artifact_selection_parent_mismatch",
            assay_name = assay_name
        )
    }
    derivation <- artifactWorkflowStateValidateAssaySelectionDerivation(
        derivation,
        manifest
    )
    matches <- Filter(
        function(entry) identical(entry$assay_name, assay_name),
        derivation$assays
    )
    if (length(matches) != 1L) {
        artifactWorkflowStateAbort(
            "artifact assay selection cannot resolve its parent assay recipe",
            "multischolar_artifact_selection_parent_mismatch",
            assay_name = assay_name
        )
    }
    ref_name <- matches[[1L]]$recipe_ref_name
    ref <- manifest$data$artifact_refs[[ref_name]]
    if (!is.list(ref)) {
        artifactWorkflowStateAbort(
            "artifact assay selection parent recipe reference is missing",
            "multischolar_artifact_selection_parent_mismatch",
            assay_name = assay_name
        )
    }
    list(ref_name = ref_name, ref = ref)
}

artifactWorkflowStateWriteAssaySelectionData <- function(
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
    hint <- artifactWorkflowStateValidateAssaySelectionHint(hint)
    plans <- artifactWorkflowStateAssaySelectionPlan(
        parent_object,
        state_object,
        hint
    )
    assay_names <- names(plans)
    target_keys <- artifactWorkflowStateAssayPayloadKeys(
        bundle$metadata,
        hint$slot_name,
        assay_names
    )
    parent_metadata <- artifactWorkflowStateUnserializeMetadata(
        previous_manifest$data$metadata_json,
        "parent S4 codec metadata"
    )
    parent_target_keys <- artifactWorkflowStateAssayPayloadKeys(
        parent_metadata,
        hint$slot_name,
        assay_names
    )
    references <- list()
    new_reference_names <- character()
    for (payload_key in setdiff(names(bundle$payloads), unname(target_keys))) {
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

    derivation_assays <- vector("list", length(assay_names))
    dependencies <- vector("list", length(assay_names))
    for (index in seq_along(assay_names)) {
        assay_name <- assay_names[[index]]
        plan <- plans[[assay_name]]
        recipe_name <- artifactWorkflowStateAssaySelectionRefName(
            .ARTIFACT_ASSAY_SELECTION_RECIPE_PREFIX,
            index
        )
        references[[recipe_name]] <-
            artifactWorkflowStateWriteAssaySelectionReference(
                store,
                identity,
                generation_id,
                plan$recipe,
                paste(identity$omic_type, hint$method, assay_name, sep = "."),
                sprintf("assay_selection_recipe_%04d", index),
                failure_injector
            )
        new_reference_names <- c(new_reference_names, recipe_name)
        patch_name <- NULL
        if (!is.null(plan$patch)) {
            patch_name <- artifactWorkflowStateAssaySelectionRefName(
                .ARTIFACT_ASSAY_SELECTION_PATCH_PREFIX,
                index
            )
            references[[patch_name]] <-
                artifactWorkflowStateWriteAssaySelectionReference(
                    store,
                    identity,
                    generation_id,
                    plan$patch,
                    paste(identity$omic_type, hint$method, assay_name, "patch", sep = "."),
                    sprintf("assay_selection_patch_%04d", index),
                    failure_injector
                )
            new_reference_names <- c(new_reference_names, patch_name)
        }
        parent_key <- parent_target_keys[[assay_name]]
        parent <- artifactWorkflowStateAssayParentReference(
            previous_manifest,
            parent_metadata,
            hint$slot_name,
            assay_name,
            parent_key
        )
        parent_ref <- parent$ref
        dependencies[[index]] <- list(
            artifact_id = references[[recipe_name]]$artifact_id,
            depends_on_artifact_id = parent_ref$artifact_id,
            relationship_type = "assay_selection_parent",
            ordinal = as.integer(index - 1L)
        )
        derivation_assays[[index]] <- list(
            assay_name = assay_name,
            target_payload_key = target_keys[[assay_name]],
            parent_payload_key = parent_key,
            parent_ref_name = parent$ref_name,
            parent_artifact_id = parent_ref$artifact_id,
            key_columns = hint$assay_key_columns[[assay_name]],
            selected_count = plan$selected_count,
            rejected_count = plan$rejected_count,
            selected_order_digest = artifactSemanticDigest(plan$child_keys),
            rejected_order_digest = artifactSemanticDigest(
                plan$parent_keys[!plan$parent_keys %in% plan$child_keys]
            ),
            recipe_ref_name = recipe_name,
            patch_ref_name = patch_name,
            patched_columns = plan$added_columns,
            removed_columns = character()
        )
    }
    lineage <- hint$lineage
    lineage$parent_generation_id <- previous_manifest$generation_id
    lineage$parent_content_identity <- previous_manifest$data$semantic_digest
    bundle$metadata$derivation <- list(
        schema = .ARTIFACT_ASSAY_SELECTION_SCHEMA,
        schema_version = .ARTIFACT_ASSAY_SELECTION_VERSION,
        representation = "assay_selection",
        target_slot = hint$slot_name,
        parent_generation_id = previous_manifest$generation_id,
        parent_semantic_digest = previous_manifest$data$semantic_digest,
        method = hint$method,
        normalized_parameters = hint$normalized_parameters,
        software = hint$software,
        lineage = lineage,
        assays = derivation_assays,
        compaction = hint$compaction
    )
    list(
        s4_class = bundle$metadata$class_name,
        semantic_digest = bundle$metadata$semantic_digest,
        codec = bundle$metadata$codec,
        metadata_json = artifactWorkflowStateSerializeMetadata(
            bundle$metadata,
            "S4 assay-selection codec metadata"
        ),
        artifact_refs = references,
        reused = FALSE,
        new_reference_names = new_reference_names,
        dependencies = dependencies
    )
}

artifactWorkflowStateValidateAssayDerivationEntry <- function(entry) {
    required <- c(
        "assay_name", "target_payload_key", "parent_payload_key",
        "parent_ref_name", "parent_artifact_id", "key_columns", "selected_count",
        "rejected_count", "selected_order_digest", "rejected_order_digest",
        "recipe_ref_name", "patch_ref_name", "patched_columns",
        "removed_columns"
    )
    valid <- is.list(entry) && identical(names(entry), required) &&
        all(vapply(entry[c(
            "assay_name", "target_payload_key", "parent_payload_key",
            "parent_ref_name", "parent_artifact_id", "recipe_ref_name"
        )], workflowStateScalarString, logical(1))) &&
        is.character(entry$key_columns) && length(entry$key_columns) > 0L &&
        !anyNA(entry$key_columns) && all(nzchar(entry$key_columns)) &&
        !anyDuplicated(entry$key_columns) &&
        artifactWorkflowStateCount(entry$selected_count) &&
        artifactWorkflowStateCount(entry$rejected_count) &&
        artifactWorkflowStateOptionalScalarString(entry$patch_ref_name) &&
        is.character(entry$patched_columns) && !anyNA(entry$patched_columns) &&
        all(nzchar(entry$patched_columns)) &&
        !anyDuplicated(entry$patched_columns) &&
        identical(entry$removed_columns, character())
    if (!isTRUE(valid)) {
        artifactWorkflowStateAbort(
            "artifact assay-selection derivation entry is malformed",
            "multischolar_invalid_artifact_assay_selection"
        )
    }
    artifactRefValidateId(
        entry$parent_artifact_id,
        "parent_artifact_id",
        "art"
    )
    artifactRefValidateDigest(
        entry$selected_order_digest,
        "selected_order_digest"
    )
    artifactRefValidateDigest(
        entry$rejected_order_digest,
        "rejected_order_digest"
    )
    entry
}

artifactWorkflowStateValidateAssaySelectionDerivation <- function(
    derivation,
    manifest
) {
    required <- c(
        "schema", "schema_version", "representation", "target_slot",
        "parent_generation_id", "parent_semantic_digest", "method",
        "normalized_parameters", "software", "lineage", "assays",
        "compaction"
    )
    valid <- is.list(derivation) && identical(names(derivation), required) &&
        identical(derivation$schema, .ARTIFACT_ASSAY_SELECTION_SCHEMA) &&
        identical(
            as.integer(derivation$schema_version),
            .ARTIFACT_ASSAY_SELECTION_VERSION
        ) && identical(derivation$representation, "assay_selection") &&
        workflowStateScalarString(derivation$target_slot) &&
        identical(
            derivation$parent_generation_id,
            manifest$parent_generation_id
        ) && workflowStateScalarString(derivation$method) &&
        is.list(derivation$normalized_parameters) &&
        artifactWorkflowStateValidSelectionSoftware(derivation$software) &&
        artifactWorkflowStateValidSelectionLineage(derivation$lineage) &&
        is.list(derivation$assays) && length(derivation$assays) > 0L &&
        is.list(derivation$compaction) &&
        identical(derivation$compaction$enabled, FALSE) &&
        workflowStateScalarString(derivation$compaction$reason)
    if (!isTRUE(valid)) {
        artifactWorkflowStateAbort(
            "artifact assay-selection derivation is malformed",
            "multischolar_invalid_artifact_assay_selection"
        )
    }
    artifactRefValidateDigest(
        derivation$parent_semantic_digest,
        "parent_semantic_digest"
    )
    derivation$assays <- lapply(
        derivation$assays,
        artifactWorkflowStateValidateAssayDerivationEntry
    )
    assay_names <- vapply(
        derivation$assays,
        `[[`,
        character(1),
        "assay_name"
    )
    if (anyDuplicated(assay_names)) {
        artifactWorkflowStateAbort(
            "artifact assay-selection derivation repeats an assay",
            "multischolar_invalid_artifact_assay_selection"
        )
    }
    artifactWorkflowStateAssertSafeMetadata(
        derivation,
        "assay-selection derivation"
    )
    derivation
}

artifactWorkflowStateDecodeSelectionRecipe <- function(
    payload_entries,
    entry,
    parent_keys
) {
    recipe_entry <- payload_entries[[entry$recipe_ref_name]]
    if (is.null(recipe_entry)) {
        artifactWorkflowStateAbort(
            "artifact assay-selection recipe is missing",
            "multischolar_invalid_artifact_assay_selection"
        )
    }
    recipe <- decodeArtifactRectangular(
        recipe_entry$payload,
        recipe_entry$metadata
    )
    required <- c(
        "entity_key", "selection_status", "parent_order", "selected_order",
        entry$key_columns, "failure_reason"
    )
    valid <- identical(names(recipe), required) &&
        identical(recipe$entity_key, parent_keys) &&
        identical(as.integer(recipe$parent_order), seq_along(parent_keys)) &&
        all(recipe$selection_status %in% c("selected", "rejected")) &&
        all(is.na(recipe$failure_reason[recipe$selection_status == "selected"])) &&
        all(!is.na(recipe$failure_reason[recipe$selection_status == "rejected"]))
    if (!isTRUE(valid)) {
        artifactWorkflowStateAbort(
            "artifact assay-selection recipe does not match its parent assay",
            "multischolar_artifact_selection_parent_mismatch"
        )
    }
    recipe
}

artifactWorkflowStateHydrateAssayValue <- function(
    entry,
    parent,
    metadata,
    payload_entries
) {
    parent_keys <- artifactWorkflowStateEntityKeys(parent, entry$key_columns)
    recipe <- artifactWorkflowStateDecodeSelectionRecipe(
        payload_entries,
        entry,
        parent_keys
    )
    selected <- recipe$selection_status == "selected"
    selected_order <- recipe$selected_order[selected]
    valid_order <- identical(
        sort(as.integer(selected_order)),
        seq_len(sum(selected))
    ) && all(is.na(recipe$selected_order[!selected]))
    selected_keys <- recipe$entity_key[selected][order(selected_order)]
    rejected_keys <- recipe$entity_key[!selected]
    valid_summary <- identical(sum(selected), entry$selected_count) &&
        identical(sum(!selected), entry$rejected_count) &&
        identical(
            artifactSemanticDigest(selected_keys),
            entry$selected_order_digest
        ) && identical(
            artifactSemanticDigest(rejected_keys),
            entry$rejected_order_digest
        )
    if (!isTRUE(valid_order) || !isTRUE(valid_summary)) {
        artifactWorkflowStateAbort(
            "artifact assay-selection counts or ordering are invalid",
            "multischolar_invalid_artifact_assay_selection"
        )
    }
    value <- parent[match(selected_keys, parent_keys), , drop = FALSE]
    if (!is.null(entry$patch_ref_name)) {
        patch_entry <- payload_entries[[entry$patch_ref_name]]
        if (is.null(patch_entry)) {
            artifactWorkflowStateAbort(
                "artifact assay-selection patch is missing",
                "multischolar_invalid_artifact_assay_selection"
            )
        }
        patch <- decodeArtifactRectangular(
            patch_entry$payload,
            patch_entry$metadata
        )
        valid_patch <- identical(patch$entity_key, selected_keys) &&
            identical(
                setdiff(names(patch), "entity_key"),
                entry$patched_columns
            )
        if (!isTRUE(valid_patch)) {
            artifactWorkflowStateAbort(
                "artifact assay-selection patch differs from selected entities",
                "multischolar_invalid_artifact_assay_selection"
            )
        }
        for (column in entry$patched_columns) {
            value[[column]] <- patch[[column]]
        }
    } else if (length(entry$patched_columns) > 0L) {
        artifactWorkflowStateAbort(
            "artifact assay-selection metadata requires an absent patch",
            "multischolar_invalid_artifact_assay_selection"
        )
    }
    expected <- metadata$payloads[[entry$target_payload_key]]
    value <- artifactWorkflowStateApplyExpectedTableMetadata(value, expected)
    encoded <- encodeArtifactTable(value, owner = expected$owner)
    if (!identical(encoded$metadata, expected)) {
        artifactWorkflowStateAbort(
            "artifact assay-selection hydration changed the expected table",
            "multischolar_inexact_artifact_state_hydration"
        )
    }
    encoded$payload
}

artifactWorkflowStateHydrateAssaySelectionPayloads <- function(
    store,
    manifest,
    metadata,
    payload_entries,
    hydrate_parent_fn
) {
    derivation <- artifactWorkflowStateValidateAssaySelectionDerivation(
        metadata$derivation,
        manifest
    )
    assay_names <- vapply(
        derivation$assays,
        `[[`,
        character(1),
        "assay_name"
    )
    target_keys <- artifactWorkflowStateAssayPayloadKeys(
        metadata,
        derivation$target_slot,
        assay_names
    )
    recorded_target_keys <- stats::setNames(vapply(
        derivation$assays,
        `[[`,
        character(1),
        "target_payload_key"
    ), assay_names)
    if (!identical(target_keys, recorded_target_keys)) {
        artifactWorkflowStateAbort(
            "artifact assay-selection targets differ from S4 codec metadata",
            "multischolar_invalid_artifact_assay_selection"
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
            "artifact assay-selection parent content identity changed",
            "multischolar_artifact_selection_parent_mismatch"
        )
    }
    parent_metadata <- artifactWorkflowStateUnserializeMetadata(
        parent_manifest$data$metadata_json,
        "parent S4 codec metadata"
    )
    parent_keys <- artifactWorkflowStateAssayPayloadKeys(
        parent_metadata,
        derivation$target_slot,
        assay_names
    )
    parent_object <- hydrate_parent_fn(parent_manifest)
    parent_assays <- methods::slot(parent_object, derivation$target_slot)
    derived_payloads <- list()
    for (index in seq_along(derivation$assays)) {
        entry <- derivation$assays[[index]]
        parent_key <- parent_keys[[entry$assay_name]]
        parent <- artifactWorkflowStateAssayParentReference(
            parent_manifest,
            parent_metadata,
            derivation$target_slot,
            entry$assay_name,
            parent_key
        )
        parent_ref <- parent$ref
        valid_parent <- identical(parent_key, entry$parent_payload_key) &&
            identical(parent$ref_name, entry$parent_ref_name) &&
            is.list(parent_ref) &&
            identical(parent_ref$artifact_id, entry$parent_artifact_id)
        if (!isTRUE(valid_parent)) {
            artifactWorkflowStateAbort(
                "artifact assay-selection parent ref differs from its recipe",
                "multischolar_artifact_selection_parent_mismatch",
                assay_name = entry$assay_name
            )
        }
        derived_payloads[[entry$target_payload_key]] <-
            artifactWorkflowStateHydrateAssayValue(
                entry,
                parent_assays[[entry$assay_name]],
                metadata,
                payload_entries
            )
    }
    payloads <- lapply(names(metadata$payloads), function(payload_key) {
        derived <- derived_payloads[[payload_key]]
        if (!is.null(derived)) return(derived)
        entry <- payload_entries[[payload_key]]
        if (is.null(entry)) {
            artifactWorkflowStateAbort(
                "artifact assay-selection state is missing a required payload",
                "multischolar_invalid_artifact_assay_selection",
                payload_key = payload_key
            )
        }
        entry$payload
    })
    names(payloads) <- names(metadata$payloads)
    payloads
}

artifactWorkflowStateWritePersistenceData <- function(
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
    representation <- hint$representation
    if (identical(representation, "row_selection")) {
        return(artifactWorkflowStateWriteSelectionData(
            store,
            identity,
            generation_id,
            state_object,
            parent_object,
            previous_manifest,
            bundle,
            hint,
            failure_injector
        ))
    }
    if (identical(representation, "assay_selection")) {
        return(artifactWorkflowStateWriteAssaySelectionData(
            store,
            identity,
            generation_id,
            state_object,
            parent_object,
            previous_manifest,
            bundle,
            hint,
            failure_injector
        ))
    }
    if (identical(representation, "state_reference")) {
        return(artifactWorkflowStateWriteStateReferenceData(
            previous_manifest,
            parent_object,
            state_object,
            hint
        ))
    }
    artifactWorkflowStateAbort(
        "artifact persistence hint has an unsupported representation",
        "multischolar_invalid_artifact_state_persistence_hint"
    )
}
