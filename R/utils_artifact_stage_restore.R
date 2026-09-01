# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

#' Test a serialized read-through contract version
#' @param value Serialized version value.
#' @param expected Expected integer contract version.
#' @return A scalar logical.
#' @noRd
artifactStageReadthroughContractVersion <- function(value, expected = 1L) {
    identical(workflowStateVersionValue(value), as.integer(expected))
}

#' Normalize a named artifact reference list
#' @param refs Named artifact references.
#' @return Named normalized artifact references.
#' @noRd
artifactStageNormalizedRefs <- function(refs) {
    normalized <- lapply(refs, artifactStoreNormalizeRef)
    names(normalized) <- names(refs)
    normalized
}

#' Restore one optional artifact table
#' @param adapter Validated read-through adapter.
#' @param value Hydrated optional table.
#' @param available Provenance availability flag.
#' @param role Artifact role.
#' @return The hydrated table or `NULL` when absence is proven.
#' @noRd
artifactStageRestoreOptional <- function(adapter, value, available, role) {
    if (!is.logical(available) || length(available) != 1L || is.na(available)) {
        artifactStageReadthroughAbort(
            adapter,
            "invalid_optional_payload",
            sprintf(
                "%s '%s' availability flag is invalid",
                adapter$owner_label,
                role
            )
        )
    }
    if (isTRUE(available)) return(value)
    expected <- artifactStageOptionalTable(
        NULL,
        role,
        adapter$abort_fn
    )
    if (!identical(value, expected)) {
        artifactStageReadthroughAbort(
            adapter,
            "invalid_optional_payload",
            sprintf(
                "%s '%s' absence marker is invalid",
                adapter$owner_label,
                role
            )
        )
    }
    NULL
}

#' Restore the exact contrasts representation from provenance
#' @param adapter Validated read-through adapter.
#' @param value Hydrated contrasts table.
#' @param kind Recorded contrasts representation.
#' @return A data frame, character vector, or `NULL`.
#' @noRd
artifactStageRestoreContrasts <- function(adapter, value, kind) {
    if (identical(kind, "data.frame")) return(value)
    if (identical(kind, "character") && identical(names(value), "contrasts")) {
        return(as.character(value$contrasts))
    }
    expected <- artifactStageOptionalTable(
        NULL,
        "contrasts",
        adapter$abort_fn
    )
    if (identical(kind, "null") && identical(value, expected)) return(NULL)
    artifactStageReadthroughAbort(
        adapter,
        "invalid_contrasts",
        sprintf(
            "%s contrast representation is incompatible with provenance",
            adapter$owner_label
        )
    )
}

#' Restore an exact serialized import column mapping
#' @param adapter Validated read-through adapter.
#' @param parameters Decoded import run parameters.
#' @return The restored column mapping list.
#' @noRd
artifactStageColumnMapping <- function(adapter, parameters) {
    encoded <- parameters$column_mapping_serialized
    mapping <- artifactWorkflowStateUnserializeMetadata(
        encoded,
        sprintf("%s import column mapping", adapter$owner_label)
    )
    if (!is.list(mapping)) {
        artifactStageReadthroughAbort(
            adapter,
            "invalid_column_mapping",
            sprintf(
                "%s import column mapping is not an R list",
                adapter$owner_label
            )
        )
    }
    mapping
}
