# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

artifactStoreInventoryRow <- function(artifact_id, state, payload_path,
                                      sidecar_path, registry_state, detail = NA_character_) {
    data.frame(
        artifact_id = artifact_id,
        state = state,
        payload_path = payload_path,
        sidecar_path = sidecar_path,
        registry_state = registry_state,
        detail = detail,
        stringsAsFactors = FALSE
    )
}

reconcileArtifactStore <- function(store, registered_artifact_ids = NULL,
                                   repair = FALSE) {
    store <- validateArtifactStore(store)
    if (!is.null(registered_artifact_ids) &&
        (!is.character(registered_artifact_ids) ||
            anyNA(registered_artifact_ids) ||
            anyDuplicated(registered_artifact_ids) > 0L)) {
        artifactStoreAbort(
            "registered artifact IDs must be a unique character vector",
            "multischolar_invalid_registry_artifact_ids"
        )
    }
    recovery_errors <- list()
    if (isTRUE(repair)) {
        for (intent_path in artifactStoreIntentPaths(store)) {
            tryCatch(
                artifactStoreRecoverIntent(store, intent_path),
                error = \(error) {
                    recovery_errors[[intent_path]] <<- conditionMessage(error)
                }
            )
        }
    }
    rows <- list()
    seen <- character()
    for (sidecar_path in artifactStoreSidecarPaths(store)) {
        result <- tryCatch(
            artifactStoreReadSidecar(store, sidecar_path, validate_payload = TRUE),
            error = identity
        )
        if (inherits(result, "error")) {
            future <- inherits(
                result,
                "multischolar_unsupported_artifact_sidecar_version"
            )
            rows[[length(rows) + 1L]] <- artifactStoreInventoryRow(
                NA_character_,
                if (future) "future_schema" else "corrupt_sidecar_or_payload",
                NA_character_,
                sidecar_path,
                "unknown",
                conditionMessage(result)
            )
            next
        }
        ref <- result$artifact_ref
        seen <- c(seen, ref$artifact_id)
        registry_state <- if (is.null(registered_artifact_ids)) {
            "not_compared"
        } else if (ref$artifact_id %in% registered_artifact_ids) {
            "registered"
        } else {
            "unregistered"
        }
        rows[[length(rows) + 1L]] <- artifactStoreInventoryRow(
            ref$artifact_id,
            "committed",
            ref$relative_path,
            sidecar_path,
            registry_state
        )
    }
    for (intent_path in artifactStoreIntentPaths(store)) {
        result <- tryCatch(artifactStoreReadIntent(store, intent_path), error = identity)
        if (inherits(result, "error")) {
            future <- inherits(
                result,
                "multischolar_unsupported_artifact_intent_version"
            )
            rows[[length(rows) + 1L]] <- artifactStoreInventoryRow(
                NA_character_,
                if (future) "future_schema" else "malformed_intent",
                NA_character_,
                NA_character_,
                "unknown",
                conditionMessage(result)
            )
            next
        }
        ref <- result$artifact_ref
        if (ref$artifact_id %in% seen) next
        exists_at <- function(path) {
            resolved <- artifactStoreResolveFile(store, path)
            file.exists(resolved) || dir.exists(resolved)
        }
        flags <- c(
            final_payload = exists_at(result$managed_paths$payload),
            final_sidecar = exists_at(result$managed_paths$sidecar),
            temporary_payload = exists_at(result$temporary_paths$payload),
            temporary_sidecar = exists_at(result$temporary_paths$sidecar)
        )
        state <- if (flags[["final_sidecar"]] && !flags[["final_payload"]]) {
            "sidecar_without_payload"
        } else if (flags[["final_payload"]]) {
            "payload_published_sidecar_pending"
        } else if (flags[["temporary_payload"]] && flags[["temporary_sidecar"]]) {
            "validated_temporary"
        } else if (flags[["temporary_payload"]]) {
            "unvalidated_temporary"
        } else {
            "intent_only"
        }
        rows[[length(rows) + 1L]] <- artifactStoreInventoryRow(
            ref$artifact_id,
            state,
            if (flags[["final_payload"]]) {
                result$managed_paths$payload
            } else {
                result$temporary_paths$payload
            },
            if (flags[["final_sidecar"]]) {
                result$managed_paths$sidecar
            } else {
                result$temporary_paths$sidecar
            },
            if (is.null(registered_artifact_ids)) "not_compared" else "unregistered",
            recovery_errors[[intent_path]] %||% NA_character_
        )
    }
    if (!is.null(registered_artifact_ids)) {
        missing_ids <- setdiff(registered_artifact_ids, seen)
        for (artifact_id in missing_ids) {
            rows[[length(rows) + 1L]] <- artifactStoreInventoryRow(
                artifact_id,
                "registry_missing_files",
                NA_character_,
                NA_character_,
                "registered"
            )
        }
    }
    if (length(rows) == 0L) {
        return(artifactStoreInventoryRow(
            character(), character(), character(), character(), character(), character()
        ))
    }
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}

artifactStoreTrash <- function(store, ref, current_generation_ids = character(),
                               referenced_generation_ids = character()) {
    store <- validateArtifactStore(store)
    ref <- validateArtifactRef(ref)
    artifactStoreValidateLogicalKey(store, ref$logical_key)
    protected <- unique(c(current_generation_ids, referenced_generation_ids))
    if (!is.character(protected) || anyNA(protected)) {
        artifactStoreAbort(
            "protected generation IDs must be character values",
            "multischolar_invalid_artifact_retention"
        )
    }
    if (ref$logical_key$generation_id %in% protected) {
        artifactStoreAbort(
            "current or referenced artifact generations cannot be trashed",
            "multischolar_protected_artifact_generation"
        )
    }
    managed <- artifactStoreManagedPaths(store, ref$logical_key, ref$artifact_id)
    sidecar <- artifactStoreReadSidecar(store, managed$sidecar, validate_payload = TRUE)
    if (!identical(sidecar$artifact_ref, ref)) {
        artifactStoreAbort(
            "trash reference does not match the immutable artifact sidecar",
            "multischolar_artifact_ref_mismatch"
        )
    }
    trash_relative <- artifactNormalizeRelativePath(file.path(
        store$relative_paths$trash,
        ref$artifact_id
    ))
    trash_path <- artifactStoreResolveFile(store, trash_relative)
    if (file.exists(trash_path) || dir.exists(trash_path)) {
        artifactStoreAbort(
            "artifact already has a trash destination",
            "multischolar_artifact_already_trashed"
        )
    }
    artifactStoreEnsureDirectory(store, store$relative_paths$trash)
    if (!dir.create(trash_path, recursive = FALSE, mode = "0700")) {
        artifactStoreAbort(
            "artifact trash destination could not be created",
            "multischolar_artifact_trash_failed"
        )
    }
    payload_directory <- artifactStoreResolveFile(
        store,
        managed$payload_directory,
        must_exist = TRUE
    )
    sidecar_path <- artifactStoreResolveFile(store, managed$sidecar, must_exist = TRUE)
    trash_payload <- file.path(trash_path, "payload")
    trash_sidecar <- file.path(trash_path, "artifact.json")
    payload_moved <- suppressWarnings(file.rename(payload_directory, trash_payload))
    if (!isTRUE(payload_moved)) {
        artifactStoreAbort(
            "artifact payload could not be moved to project-local trash",
            "multischolar_artifact_trash_failed"
        )
    }
    sidecar_moved <- suppressWarnings(file.rename(sidecar_path, trash_sidecar))
    if (!isTRUE(sidecar_moved)) {
        rolled_back <- suppressWarnings(file.rename(trash_payload, payload_directory))
        artifactStoreAbort(
            "artifact sidecar trash move failed",
            if (rolled_back) {
                "multischolar_artifact_trash_failed"
            } else {
                "multischolar_partial_artifact_trash"
            }
        )
    }
    list(
        artifact_id = ref$artifact_id,
        generation_id = ref$logical_key$generation_id,
        original_payload = ref$relative_path,
        trash_relative_path = trash_relative,
        status = "trashed",
        trashed_at = artifactRefUtcNow()
    )
}
