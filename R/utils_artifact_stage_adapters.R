# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

artifactStageAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_artifact_stage_adapter_error"),
        ...
    )
}

artifactStageValidateFlag <- function(value, owner) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
        artifactStageAbort(
            sprintf("artifact stage adapter '%s' flag is invalid", owner),
            "multischolar_invalid_artifact_stage_adapter"
        )
    }
    value
}

hydrateArtifactStageInput <- function(
    explicit_input = NULL,
    workflow_state = NULL,
    state_name = NULL,
    legacy_provider = NULL,
    allow_legacy = FALSE,
    validate_fn = methods::validObject
) {
    allow_legacy <- artifactStageValidateFlag(allow_legacy, "allow_legacy")
    if (!is.null(explicit_input)) {
        value <- explicit_input
        source <- "explicit"
    } else if (!is.null(workflow_state)) {
        state_container <- is.environment(workflow_state) || is.list(workflow_state)
        if (!isTRUE(state_container) || !is.function(workflow_state$getState)) {
            artifactStageAbort(
                "artifact stage workflow state does not expose getState",
                "multischolar_invalid_artifact_stage_adapter"
            )
        }
        value <- workflow_state$getState(state_name)
        source <- "workflow_state"
    } else if (isTRUE(allow_legacy) && is.function(legacy_provider)) {
        value <- legacy_provider()
        source <- "legacy"
    } else {
        artifactStageAbort(
            "artifact stage adapter has no permitted scientific input source",
            "multischolar_missing_artifact_stage_input"
        )
    }
    if (is.null(value) || !is.function(validate_fn)) {
        artifactStageAbort(
            "artifact stage adapter received an invalid input or validator",
            "multischolar_invalid_artifact_stage_adapter"
        )
    }
    validity <- tryCatch(
        validate_fn(value),
        error = function(error) error
    )
    if (inherits(validity, "error") ||
        !(isTRUE(validity) || identical(validity, character()))) {
        artifactStageAbort(
            "artifact stage input failed its scientific validity contract",
            "multischolar_invalid_artifact_stage_input",
            validity = validity
        )
    }
    structure(
        list(value = value, source = source),
        class = c("MultiScholaRArtifactStageInput", "list")
    )
}

runArtifactStageAdapter <- function(
    scientific_fn,
    ...,
    explicit_input = NULL,
    workflow_state = NULL,
    state_name = NULL,
    legacy_provider = NULL,
    allow_legacy = FALSE,
    validate_fn = methods::validObject
) {
    if (!is.function(scientific_fn)) {
        artifactStageAbort(
            "artifact stage scientific owner must be a function",
            "multischolar_invalid_artifact_stage_adapter"
        )
    }
    input <- hydrateArtifactStageInput(
        explicit_input = explicit_input,
        workflow_state = workflow_state,
        state_name = state_name,
        legacy_provider = legacy_provider,
        allow_legacy = allow_legacy,
        validate_fn = validate_fn
    )
    scientific_fn(input$value, ...)
}
