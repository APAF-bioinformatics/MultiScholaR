# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

artifactWorkflowStateInstallBootstrap <- function(
    private,
    public,
    bootstrap
) {
    bootstrap <- artifactWorkflowStateValidateSettledBootstrap(
        bootstrap,
        private$store,
        private$identity
    )
    private$state_rows <- bootstrap$state_rows
    private$settled_metadata <- bootstrap$state_metadata
    private$settled_workflow_type <- bootstrap$workflow_type
    private$settled_audit_enabled <- bootstrap$audit_enabled
    artifactWorkflowStateRefreshCached(private, public)
    invisible(TRUE)
}

artifactWorkflowStateAssertUsableInternal <- function(private) {
    if (isTRUE(private$closed)) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState is closed",
            "multischolar_closed_artifact_workflow_state"
        )
    }
    artifactResourceAssertCreatorProcess(
        private$process_id,
        "artifact WorkflowState"
    )
    invisible(TRUE)
}

artifactWorkflowStateEnsureOpen <- function(private) {
    artifactWorkflowStateAssertUsableInternal(private)
    if (is.null(private$session)) {
        if (is.null(private$registry)) {
            private$registry <- projectRegistryForContext(
                private$context,
                private$resource_policy
            )
        }
        private$session <- initializeProjectRegistry(private$registry)
        artifactWorkflowStateEnsureWorkflow(
            private$session,
            private$registry_identity
        )
        private$refresh()
    }
    invisible(TRUE)
}

artifactWorkflowStateReleaseRegistryInternal <- function(private) {
    if (isTRUE(private$closed) || is.null(private$session)) {
        return(invisible(FALSE))
    }
    if (!is.null(private$cache_generation_id) ||
        length(private$observers) > 0L) {
        artifactWorkflowStateAbort(
            "artifact WorkflowState cannot suspend with live owned resources",
            "multischolar_artifact_state_resources_live"
        )
    }
    cleanup <- closeArtifactWorkflowStateSession(private$session)
    if (isTRUE(cleanup$closed)) private$session <- NULL
    if (!is.null(cleanup$error)) stop(cleanup$error)
    invisible(TRUE)
}

artifactWorkflowStateCloseInternal <- function(private) {
    if (isTRUE(private$closed)) return(invisible(FALSE))
    private$clearCache()
    private$observers <- list()
    if (is.null(private$session)) {
        private$closed <- TRUE
        return(invisible(TRUE))
    }
    cleanup <- closeArtifactWorkflowStateSession(private$session)
    if (isTRUE(cleanup$closed)) {
        private$session <- NULL
        private$closed <- TRUE
    }
    if (!is.null(cleanup$error)) stop(cleanup$error)
    invisible(TRUE)
}

artifactWorkflowStateRefreshCached <- function(private, public) {
    current <- private$currentRow()
    private$current_generation_id <- if (is.null(current)) {
        NULL
    } else {
        current$generation_id
    }
    if (is.null(current)) {
        public$states <- list()
        public$state_history <- list()
        public$audit_records <- list()
        public$current_state <- NULL
        return(invisible(TRUE))
    }
    lineage <- private$activeLineageRows()
    history <- vapply(lineage, `[[`, character(1), "logical_name")
    history <- history[!duplicated(history, fromLast = TRUE)]
    public$state_history <- as.list(history)
    public$current_state <- current$logical_name
    latest <- list()
    audits <- list()
    for (row in private$state_rows) {
        metadata <- private$settled_metadata[[row$generation_id]]
        manifest <- if (is.list(metadata)) {
            NULL
        } else {
            private$manifestForRow(row)
        }
        latest[[row$logical_name]] <- list(
            generation_id = row$generation_id,
            parent_generation_id = row$parent_generation_id,
            manifest_relative_path = row$manifest_relative_path,
            status = row$status,
            artifact_refs = if (is.list(metadata)) {
                metadata$artifact_refs
            } else {
                manifest$data$artifact_refs
            }
        )
        audit <- if (is.list(metadata)) {
            metadata$audit_metadata
        } else {
            artifactWorkflowStateUnserializeMetadata(
                manifest$audit_json,
                "audit metadata"
            )
        }
        if (is.list(audit) && workflowStateScalarString(audit$record_id)) {
            existing <- audits[[audit$record_id]]
            if (!is.null(existing) && !identical(existing, audit)) {
                artifactWorkflowStateAbort(
                    "artifact WorkflowState audit record IDs are not immutable",
                    "multischolar_artifact_state_audit_mismatch"
                )
            }
            audits[[audit$record_id]] <- audit
        }
    }
    public$states <- latest
    public$audit_records <- audits
    if (!is.null(private$settled_workflow_type)) {
        public$audit_enabled <- private$settled_audit_enabled
        public$workflow_type <- private$settled_workflow_type
    } else {
        current_manifest <- private$manifestForRow(current)
        public$audit_enabled <- current_manifest$audit_enabled
        public$workflow_type <- current_manifest$workflow_type
    }
    if (!identical(private$cache_generation_id, current$generation_id)) {
        private$clearCache()
    }
    invisible(TRUE)
}
