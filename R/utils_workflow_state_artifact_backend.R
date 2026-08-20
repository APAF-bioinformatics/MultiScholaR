# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

ArtifactWorkflowState <- R6::R6Class(
    "ArtifactWorkflowState",
    public = list(
        states = list(),
        current_state = NULL,
        state_history = list(),
        audit_records = list(),
        audit_enabled = TRUE,
        workflow_type = "LFQ",

        initialize = function(
            workflow_context,
            audit_enabled = TRUE,
            resource_policy = NULL,
            dehydrate_fn = dehydrateDiaS4Artifact,
            validate_bundle_fn = validateDiaS4Bundle,
            hydrate_fn = hydrateDiaS4Artifact, descriptor_contract = NULL
        ) {
            context <- artifactWorkflowStateValidateContext(workflow_context)
            if (!is.function(dehydrate_fn) || !is.function(validate_bundle_fn) ||
                !is.function(hydrate_fn)) {
                artifactWorkflowStateAbort(
                    "artifact WorkflowState codec adapter is invalid",
                    "multischolar_invalid_artifact_state_codec"
                )
            }
            private$context <- context
            private$identity <- context$getIdentity()
            private$store <- newArtifactStore(
                context$getPaths(),
                private$identity$project_id
            )
            artifactWorkflowStateEnsureMetadata(private$store, private$identity,
                descriptor_contract)
            private$registry <- projectRegistryForContext(
                context,
                resource_policy = resource_policy
            )
            private$dehydrate_fn <- dehydrate_fn
            private$validate_bundle_fn <- validate_bundle_fn
            private$hydrate_fn <- hydrate_fn
            self$audit_enabled <- isTRUE(audit_enabled)
            private$observers <- list()
            private$next_observer_id <- 0L
            private$cache_generation_id <- NULL
            private$cache_object <- NULL
            private$hydration_count <- 0L
            private$closed <- FALSE
            private$process_id <- as.integer(Sys.getpid())
            private$session <- initializeProjectRegistry(private$registry)
            initialized <- FALSE
            on.exit({
                if (!initialized && !is.null(private$session)) {
                    try(closeProjectRegistry(private$session), silent = TRUE)
                }
            }, add = TRUE)
            artifactWorkflowStateEnsureWorkflow(private$session, private$identity)
            private$refresh()
            if (is.null(private$current_generation_id)) {
                private$createInitialGeneration()
            }
            initialized <- TRUE
        },
        saveState = function(
            state_name,
            s4_data_object,
            config_object,
            description,
            audit_metadata = NULL
        ) {
            result <- self$commitState(
                state_name = state_name,
                s4_data_object = s4_data_object,
                config_object = config_object,
                description = description,
                audit_metadata = audit_metadata,
                action_id = artifactOpaqueId("action"),
                expected_parent_generation_id = private$current_generation_id
            )
            if (!identical(result$status, "accepted")) {
                artifactWorkflowStateAbort(
                    "artifact WorkflowState save did not advance the current generation",
                    "multischolar_artifact_state_save_rejected",
                    result = result
                )
            }
            state_name
        },
        commitState = function(
            state_name,
            s4_data_object,
            config_object,
            description,
            audit_metadata = NULL,
            action_id = artifactOpaqueId("action"),
            expected_parent_generation_id = NULL
        ) {
            private$assertOpen()
            private$validateStateName(state_name)
            private$validateActionId(action_id)
            if (is.null(expected_parent_generation_id)) {
                expected_parent_generation_id <- private$current_generation_id
            }
            private$validateGenerationId(expected_parent_generation_id)
            duplicate <- private$actionResult(action_id)
            if (!is.null(duplicate)) {
                duplicate$idempotent <- TRUE
                return(duplicate)
            }
            audit_metadata <- private$normalizeAuditMetadata(
                s4_data_object,
                audit_metadata
            )
            artifactWorkflowStateValidateSaveMetadata(
                config_object, audit_metadata, description
            )
            actual_parent <- private$currentGenerationId()
            if (!identical(actual_parent, expected_parent_generation_id)) {
                return(private$recordRejectedAction(
                    action_id,
                    state_name,
                    expected_parent_generation_id,
                    actual_parent
                ))
            }
            private$commitStateInternal(
                state_name = state_name,
                state_object = s4_data_object,
                config = config_object,
                description = description,
                audit_metadata = audit_metadata,
                action_id = action_id,
                expected_parent = expected_parent_generation_id,
                event_type = "state_saved"
            )
        },
        getState = function(state_name = NULL) {
            private$assertOpen()
            row <- private$resolveStateRow(state_name)
            if (is.null(row)) return(NULL)
            is_current <- identical(row$generation_id, private$current_generation_id)
            if (is_current && identical(
                private$cache_generation_id,
                row$generation_id
            )) {
                return(private$cache_object)
            }
            manifest <- private$manifestForRow(row)
            object <- artifactWorkflowStateHydrateData(
                private$store,
                manifest,
                hydrate_fn = private$hydrate_fn
            )
            private$hydration_count <- private$hydration_count + 1L
            if (is_current) private$setCache(row$generation_id, object)
            object
        },
        getStateMetadata = function(state_name = NULL) {
            private$assertOpen()
            row <- private$resolveStateRow(state_name)
            if (is.null(row)) return(NULL)
            manifest <- private$manifestForRow(row)
            list(
                timestamp = as.POSIXct(manifest$created_at, tz = "UTC"),
                config = artifactWorkflowStateUnserializeMetadata(
                    manifest$config_json,
                    "config"
                ),
                description = manifest$description,
                previous_state = private$logicalNameForGeneration(
                    manifest$parent_generation_id
                ),
                s4_class = manifest$s4_class,
                audit_metadata = artifactWorkflowStateUnserializeMetadata(
                    manifest$audit_json,
                    "audit metadata"
                ),
                generation_id = manifest$generation_id,
                artifact_refs = manifest$data$artifact_refs,
                status = row$status
            )
        },
        getStateConfig = function(state_name = NULL) {
            metadata <- self$getStateMetadata(state_name)
            if (is.null(metadata)) NULL else metadata$config
        },
        hasState = function(state_name) {
            workflowStateScalarString(state_name) &&
                any(vapply(private$state_rows, function(row) {
                    identical(row$logical_name, state_name)
                }, logical(1)))
        },
        getCurrentStateName = function() self$current_state,
        getCurrentGenerationId = function() {
            private$current_generation_id
        },
        revertToState = function(state_name) {
            private$assertOpen()
            private$validateStateName(state_name)
            lineage <- private$activeLineageRows()
            selection <- artifactWorkflowStateSelectGeneration(
                private$state_rows,
                lineage,
                state_name
            )
            if (identical(selection$mode, "current")) {
                return(self$getState(state_name))
            }
            object <- private$preflightSelection(selection$target)
            private$activateSelection(selection, artifactOpaqueId("action"))
            private$hydration_count <- private$hydration_count + 1L
            private$setCache(selection$target$generation_id, object)
            object
        },
        getHistory = function() {
            unlist(self$state_history, use.names = FALSE)
        },
        getStateAudit = function(state_name = NULL) {
            metadata <- self$getStateMetadata(state_name)
            if (is.null(metadata) || is.null(metadata$audit_metadata)) {
                return(list(provenance_status = "legacy_untracked_state"))
            }
            metadata$audit_metadata
        },
        getAuditRecords = function() {
            self$audit_records
        },
        isAuditEnabled = function() isTRUE(self$audit_enabled),
        setWorkflowType = function(type) {
            private$assertOpen()
            if (!workflowStateScalarString(type) || !type %in% .WORKFLOW_STATE_TYPES) {
                stop(
                    paste(
                        "Invalid workflow type. Must be one of:",
                        paste(.WORKFLOW_STATE_TYPES, collapse = ", ")
                    )
                )
            }
            if (identical(type, self$workflow_type)) return(type)
            previous_type <- self$workflow_type
            current_manifest <- private$currentManifest()
            current_object <- self$getState()
            self$workflow_type <- type
            result <- tryCatch(
                private$commitStateInternal(
                    state_name = self$current_state,
                    state_object = current_object,
                    config = artifactWorkflowStateUnserializeMetadata(
                        current_manifest$config_json,
                        "config"
                    ),
                    description = current_manifest$description,
                    audit_metadata = artifactWorkflowStateUnserializeMetadata(
                        current_manifest$audit_json,
                        "audit metadata"
                    ),
                    action_id = artifactOpaqueId("action"),
                    expected_parent = private$current_generation_id,
                    event_type = "workflow_type_set",
                    event_metadata = list(previous_type = previous_type, type = type)
                ),
                error = function(error) {
                    self$workflow_type <- previous_type
                    stop(error)
                }
            )
            if (!identical(result$status, "accepted")) {
                self$workflow_type <- previous_type
                artifactWorkflowStateAbort(
                    "artifact WorkflowState workflow type update was rejected",
                    "multischolar_artifact_state_save_rejected"
                )
            }
            type
        },
        getWorkflowType = function() self$workflow_type,
        exportState = function() {
            private$assertOpen()
            artifactWorkflowStateExportSnapshot(
                private$identity,
                private$current_generation_id,
                self$current_state,
                private$activeLineageRows(),
                self$workflow_type,
                self$audit_enabled
            )
        },
        exportLegacyState = function() {
            list(
                r6_current_state_name = self$current_state,
                r6_complete_states = self$states,
                r6_state_history = self$getHistory()
            )
        },
        restoreState = function(manifest, schema_version = NULL) {
            private$assertOpen()
            artifactWorkflowStateValidateRestoreSnapshot(
                manifest,
                private$identity,
                private$currentGenerationId(),
                schema_version = schema_version
            )
            private$refresh()
            invisible(self)
        },
        getEvents = function() {
            private$assertOpen()
            artifactWorkflowStateEvents(
                private$session,
                private$identity$workflow_id
            )
        },
        getRevision = function() length(self$getEvents()),
        observeTransitions = function(callback) {
            if (!is.function(callback)) {
                stop("WorkflowState transition callback must be a function.")
            }
            private$next_observer_id <- private$next_observer_id + 1L
            observer_id <- as.character(private$next_observer_id)
            private$observers[[observer_id]] <- callback
            removed <- FALSE
            function() {
                if (!removed) {
                    private$observers[[observer_id]] <- NULL
                    removed <<- TRUE
                }
                invisible(NULL)
            }
        },

        getCacheInfo = function() {
            list(
                entries = if (is.null(private$cache_generation_id)) 0L else 1L,
                generation_id = private$cache_generation_id,
                hydration_count = private$hydration_count
            )
        },

        getResourceInfo = function() {
            artifactWorkflowStateResourceInfo(
                private$session,
                private$process_id,
                private$closed,
                private$cache_generation_id,
                length(private$observers)
            )
        },

        close = function() {
            if (isTRUE(private$closed)) return(invisible(FALSE))
            private$clearCache()
            private$observers <- list()
            cleanup <- closeArtifactWorkflowStateSession(private$session)
            if (isTRUE(cleanup$closed)) {
                private$session <- NULL
                private$closed <- TRUE
            }
            if (!is.null(cleanup$error)) stop(cleanup$error)
            invisible(TRUE)
        }
    ),
    private = list(
        context = NULL,
        identity = NULL,
        store = NULL,
        registry = NULL,
        session = NULL,
        state_rows = list(),
        current_generation_id = NULL,
        cache_generation_id = NULL,
        cache_object = NULL,
        hydration_count = 0L,
        observers = list(),
        next_observer_id = 0L,
        dehydrate_fn = NULL,
        validate_bundle_fn = NULL,
        hydrate_fn = NULL,
        closed = FALSE,
        process_id = NULL,

        finalize = function() {
            if (!isTRUE(private$closed) && !is.null(private$session)) {
                try(self$close(), silent = TRUE)
            }
        },

        assertOpen = function() {
            if (isTRUE(private$closed) || is.null(private$session)) {
                artifactWorkflowStateAbort(
                    "artifact WorkflowState is closed",
                    "multischolar_closed_artifact_workflow_state"
                )
            }
            invisible(TRUE)
        },

        validateStateName = function(state_name) {
            if (!workflowStateScalarString(state_name)) {
                stop("WorkflowState state name must be a non-empty string.")
            }
            invisible(state_name)
        },

        validateActionId = function(action_id) {
            if (!workflowStateScalarString(action_id)) {
                artifactWorkflowStateAbort(
                    "artifact WorkflowState action ID must be a non-empty string",
                    "multischolar_invalid_artifact_state_action"
                )
            }
            invisible(action_id)
        },

        validateGenerationId = function(generation_id) {
            if (!is.null(generation_id)) {
                artifactRefValidateId(generation_id, "generation_id", "gen")
            }
            invisible(generation_id)
        },

        stateRows = function() {
            rows <- projectRegistryQuery(
                private$session,
                "states",
                filters = list(workflow_id = private$identity$workflow_id)
            )
            if (nrow(rows) == 0L) return(list())
            lapply(seq_len(nrow(rows)), function(index) {
                row <- as.list(rows[index, , drop = FALSE])
                row <- lapply(row, `[[`, 1L)
                if (is.na(row$parent_generation_id)) {
                    row$parent_generation_id <- NULL
                }
                row
            })
        },

        currentRow = function() {
            current <- Filter(function(row) identical(row$status, "current"),
                private$state_rows
            )
            if (length(current) > 1L) {
                artifactWorkflowStateAbort(
                    "artifact WorkflowState has multiple current generations",
                    "multischolar_invalid_artifact_state_lineage"
                )
            }
            if (length(current) == 0L) NULL else current[[1L]]
        },

        currentGenerationId = function() {
            private$state_rows <- private$stateRows()
            row <- private$currentRow()
            if (is.null(row)) NULL else row$generation_id
        },

        rowByGeneration = function(generation_id) {
            if (is.null(generation_id)) return(NULL)
            matches <- Filter(function(row) {
                identical(row$generation_id, generation_id)
            }, private$state_rows)
            if (length(matches) != 1L) {
                artifactWorkflowStateAbort(
                    "artifact WorkflowState lineage references a missing generation",
                    "multischolar_invalid_artifact_state_lineage"
                )
            }
            matches[[1L]]
        },

        activeLineageRows = function() {
            current <- private$currentRow()
            if (is.null(current)) return(list())
            lineage <- list()
            seen <- character()
            row <- current
            while (!is.null(row)) {
                if (row$generation_id %in% seen) {
                    artifactWorkflowStateAbort(
                        "artifact WorkflowState lineage contains a cycle",
                        "multischolar_invalid_artifact_state_lineage"
                    )
                }
                seen <- c(seen, row$generation_id)
                lineage[[length(lineage) + 1L]] <- row
                row <- private$rowByGeneration(row$parent_generation_id)
            }
            rev(lineage)
        },

        logicalNameForGeneration = function(generation_id) {
            if (is.null(generation_id)) return(NULL)
            private$rowByGeneration(generation_id)$logical_name
        },

        manifestForRow = function(row) {
            manifest <- artifactWorkflowStateReadManifest(
                private$store,
                row$manifest_relative_path
            )
            if (!identical(manifest$project_id, private$identity$project_id) ||
                !identical(manifest$workflow_id, private$identity$workflow_id) ||
                !identical(manifest$generation_id, row$generation_id) ||
                !identical(manifest$logical_name, row$logical_name)) {
                artifactWorkflowStateAbort(
                    "artifact state manifest does not match its registry generation",
                    "multischolar_artifact_state_registry_mismatch"
                )
            }
            manifest
        },

        currentManifest = function() {
            private$manifestForRow(private$currentRow())
        },

        resolveStateRow = function(state_name) {
            if (is.null(state_name)) return(private$currentRow())
            private$validateStateName(state_name)
            lineage <- private$activeLineageRows()
            active <- Filter(function(row) {
                identical(row$logical_name, state_name)
            }, lineage)
            if (length(active) > 0L) return(tail(active, 1L)[[1L]])
            all_matches <- Filter(function(row) {
                identical(row$logical_name, state_name)
            }, private$state_rows)
            if (length(all_matches) == 0L) NULL else tail(all_matches, 1L)[[1L]]
        },

        refresh = function() {
            private$state_rows <- private$stateRows()
            current <- private$currentRow()
            private$current_generation_id <- if (is.null(current)) {
                NULL
            } else {
                current$generation_id
            }
            if (is.null(current)) {
                self$states <- list()
                self$state_history <- list()
                self$audit_records <- list()
                self$current_state <- NULL
                return(invisible(TRUE))
            }
            lineage <- private$activeLineageRows()
            history <- vapply(lineage, `[[`, character(1), "logical_name")
            history <- history[!duplicated(history, fromLast = TRUE)]
            self$state_history <- as.list(history)
            self$current_state <- current$logical_name
            latest <- list()
            audits <- list()
            for (row in private$state_rows) {
                manifest <- private$manifestForRow(row)
                summary <- list(
                    generation_id = row$generation_id,
                    parent_generation_id = row$parent_generation_id,
                    manifest_relative_path = row$manifest_relative_path,
                    status = row$status,
                    artifact_refs = manifest$data$artifact_refs
                )
                latest[[row$logical_name]] <- summary
                audit <- artifactWorkflowStateUnserializeMetadata(
                    manifest$audit_json,
                    "audit metadata"
                )
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
            self$states <- latest
            self$audit_records <- audits
            current_manifest <- private$manifestForRow(current)
            self$audit_enabled <- current_manifest$audit_enabled
            self$workflow_type <- current_manifest$workflow_type
            if (!identical(private$cache_generation_id, current$generation_id)) {
                private$clearCache()
            }
            invisible(TRUE)
        },

        createInitialGeneration = function() {
            generation_id <- artifactOpaqueId("gen")
            data <- list(
                s4_class = NULL,
                semantic_digest = NULL,
                codec = NULL,
                metadata_json = NULL,
                artifact_refs = list(),
                reused = FALSE
            )
            audit <- list(provenance_status = "legacy_or_not_applicable")
            manifest <- artifactWorkflowStateNewManifest(
                private$identity,
                generation_id,
                NULL,
                "initial",
                self$workflow_type,
                self$audit_enabled,
                "Initial empty state",
                data,
                NULL,
                audit
            )
            manifest_path <- artifactWorkflowStateWriteManifest(private$store, manifest)
            action_id <- artifactOpaqueId("action")
            timestamp <- manifest$created_at
            artifactWorkflowStateTransaction(private$session, function() {
                projectRegistryWrite(private$session, "state", list(
                    workflow_id = private$identity$workflow_id,
                    generation_id = generation_id,
                    parent_generation_id = NULL,
                    logical_name = "initial",
                    manifest_relative_path = manifest_path,
                    status = "current",
                    created_at = timestamp,
                    updated_at = timestamp
                ))
                private$writeRevision(
                    action_id,
                    generation_id,
                    NULL,
                    "accepted",
                    list(state_name = "initial")
                )
                private$writeEvent(
                    "initialized",
                    "accepted",
                    generation_id,
                    "initial",
                    NULL,
                    list(),
                    timestamp
                )
            })
            private$refresh()
            private$setCache(generation_id, NULL)
            invisible(TRUE)
        },

        normalizeAuditMetadata = function(state_object, audit_metadata) {
            if (isTRUE(self$audit_enabled) && is.null(audit_metadata) &&
                isS4(state_object) && "args" %in% methods::slotNames(state_object) &&
                is.list(state_object@args$peptide_qc_audit)) {
                lineage <- state_object@args$peptide_qc_audit
                record_id <- lineage$current_record_id
                if (!is.null(record_id) && length(record_id) == 1L) {
                    matching <- Filter(
                        function(record) identical(record$record_id, record_id),
                        lineage$records
                    )
                    if (length(matching) > 0L) audit_metadata <- matching[[1L]]
                }
            }
            if (is.null(audit_metadata)) {
                audit_metadata <- list(
                    provenance_status = "legacy_or_not_applicable"
                )
            }
            if (!is.list(audit_metadata)) {
                stop("WorkflowState audit metadata must be an R list.")
            }
            record_id <- audit_metadata$record_id
            if (!is.null(record_id)) {
                if (!workflowStateScalarString(record_id)) {
                    stop("WorkflowState audit record ID must be a non-empty string.")
                }
                existing <- self$audit_records[[record_id]]
                if (!is.null(existing) && !identical(existing, audit_metadata)) {
                    stop("WorkflowState audit record IDs must be immutable.")
                }
            }
            audit_metadata
        },

        actionResult = function(action_id) {
            rows <- projectRegistryQuery(
                private$session,
                "revisions",
                filters = list(
                    workflow_id = private$identity$workflow_id,
                    action_id = action_id
                ),
                limit = 1L
            )
            if (nrow(rows) == 0L) return(NULL)
            details <- jsonlite::fromJSON(rows$details_json[[1L]], simplifyVector = TRUE)
            generation_id <- rows$generation_id[[1L]]
            if (is.na(generation_id)) generation_id <- NULL
            list(
                status = rows$revision_status[[1L]],
                action_id = action_id,
                generation_id = generation_id,
                logical_name = details$state_name,
                idempotent = FALSE
            )
        },

        recordRejectedAction = function(action_id, state_name, expected, actual) {
            timestamp <- artifactRefUtcNow()
            artifactWorkflowStateTransaction(private$session, function() {
                private$writeRevision(
                    action_id,
                    NULL,
                    expected,
                    "rejected_stale_parent",
                    list(
                        state_name = state_name,
                        actual_parent_generation_id = actual
                    )
                )
                private$writeEvent(
                    "state_save_rejected",
                    "rejected_stale_parent",
                    NULL,
                    state_name,
                    self$current_state,
                    list(
                        expected_parent_generation_id = expected,
                        actual_parent_generation_id = actual
                    ),
                    timestamp
                )
            })
            private$notifyLastEvent()
            list(
                status = "rejected_stale_parent",
                action_id = action_id,
                generation_id = NULL,
                logical_name = state_name,
                idempotent = FALSE
            )
        },

        commitStateInternal = function(
            state_name,
            state_object,
            config,
            description,
            audit_metadata,
            action_id,
            expected_parent,
            event_type,
            event_metadata = list()
        ) {
            previous_manifest <- private$currentManifest()
            generation_id <- artifactOpaqueId("gen")
            data <- artifactWorkflowStateWriteData(
                private$store,
                private$identity,
                generation_id,
                state_object,
                previous_manifest = previous_manifest,
                dehydrate_fn = private$dehydrate_fn,
                validate_bundle_fn = private$validate_bundle_fn
            )
            manifest <- artifactWorkflowStateNewManifest(
                private$identity,
                generation_id,
                expected_parent,
                state_name,
                self$workflow_type,
                self$audit_enabled,
                description,
                data,
                config,
                audit_metadata
            )
            manifest_path <- artifactWorkflowStateWriteManifest(private$store, manifest)
            artifactWorkflowStateVerifyHydration(private$store, manifest, state_object,
                private$hydrate_fn)
            timestamp <- manifest$created_at
            artifactWorkflowStateTransaction(private$session, function() {
                actual_parent <- private$currentGenerationId()
                if (!identical(actual_parent, expected_parent)) {
                    artifactWorkflowStateAbort(
                        "artifact state parent changed before registry commit",
                        "multischolar_artifact_state_compare_and_set_failed"
                    )
                }
                if (!isTRUE(data$reused)) {
                    for (index in seq_along(data$artifact_refs)) {
                        projectRegistryWrite(
                            private$session,
                            "artifact",
                            artifactWorkflowStateArtifactRecord(
                                private$identity,
                                data$artifact_refs[[index]],
                                index - 1L
                            )
                        )
                    }
                }
                artifactWorkflowStateUpdateStatus(
                    private$session,
                    private$identity,
                    expected_parent,
                    "current",
                    "historical",
                    timestamp
                )
                projectRegistryWrite(private$session, "state", list(
                    workflow_id = private$identity$workflow_id,
                    generation_id = generation_id,
                    parent_generation_id = expected_parent,
                    logical_name = state_name,
                    manifest_relative_path = manifest_path,
                    status = "current",
                    created_at = timestamp,
                    updated_at = timestamp
                ))
                private$writeRevision(
                    action_id,
                    generation_id,
                    expected_parent,
                    "accepted",
                    list(state_name = state_name, artifacts_reused = isTRUE(data$reused))
                )
                private$writeEvent(
                    event_type,
                    "accepted",
                    generation_id,
                    state_name,
                    self$current_state,
                    event_metadata,
                    timestamp
                )
            })
            private$refresh()
            private$setCache(generation_id, state_object)
            private$notifyLastEvent()
            list(
                status = "accepted",
                action_id = action_id,
                generation_id = generation_id,
                logical_name = state_name,
                artifacts_reused = isTRUE(data$reused),
                idempotent = FALSE
            )
        },

        preflightSelection = function(target) {
            artifactWorkflowStatePreflightSelection(list(
                store = private$store,
                identity = private$identity,
                target = target,
                manifest = private$manifestForRow(target),
                hydrate_fn = private$hydrate_fn
            ))
        },

        activateSelection = function(selection, action_id) {
            request <- list(
                session = private$session,
                identity = private$identity,
                selection = selection,
                action_id = action_id,
                write_revision = private$writeRevision,
                write_event = private$writeEvent
            )
            if (identical(selection$mode, "revert")) {
                artifactWorkflowStateCommitRevert(request)
            } else {
                artifactWorkflowStateCommitResume(request)
            }
            private$refresh()
            private$notifyLastEvent()
            invisible(TRUE)
        },

        writeRevision = function(
            action_id,
            generation_id,
            expected_parent,
            status,
            details
        ) {
            projectRegistryWrite(private$session, "revision", list(
                workflow_id = private$identity$workflow_id,
                revision_id = artifactOpaqueId("revision"),
                generation_id = generation_id,
                action_id = action_id,
                expected_parent_generation_id = expected_parent,
                revision_status = status,
                details_json = artifactWorkflowStateJson(details),
                recorded_at = artifactRefUtcNow()
            ))
        },

        writeEvent = function(
            event_type,
            event_status,
            generation_id,
            state_name,
            previous_state,
            metadata,
            timestamp
        ) {
            projectRegistryWrite(private$session, "event", list(
                workflow_id = private$identity$workflow_id,
                event_id = artifactOpaqueId("event"),
                generation_id = generation_id,
                run_id = NULL,
                event_type = event_type,
                event_status = event_status,
                details_json = artifactWorkflowStateJson(list(
                    state_name = state_name,
                    previous_state = previous_state,
                    metadata = metadata
                )),
                recorded_at = timestamp
            ))
        },

        setCache = function(generation_id, object) {
            private$cache_generation_id <- generation_id
            private$cache_object <- object
            invisible(TRUE)
        },

        clearCache = function() {
            private$cache_generation_id <- NULL
            private$cache_object <- NULL
            invisible(TRUE)
        },

        notifyLastEvent = function() {
            events <- self$getEvents()
            if (length(events) == 0L) return(invisible(NULL))
            event <- tail(events, 1L)[[1L]]
            for (callback in private$observers) {
                tryCatch(
                    callback(event),
                    error = function(error) {
                        try(
                            logger::log_warn(paste(
                                "ArtifactWorkflowState transition observer failed:",
                                conditionMessage(error)
                            )),
                            silent = TRUE
                        )
                    }
                )
            }
            invisible(event)
        }
    )
)
