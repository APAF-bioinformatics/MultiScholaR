# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' @title WorkflowState
#' @description An R6 class to manage the state of a multi-step workflow.
#'
#' @details This class tracks data objects as they move through an interactive
#' workflow. Public fields remain available for serialized compatibility, while
#' named methods are the supported production access path.
#'
#' @importFrom R6 R6Class
#' @export
WorkflowState <- R6::R6Class(
    "WorkflowState",
    public = list(
        #' @field states Saved state records keyed by logical state name.
        states = list(),

        #' @field current_state Name of the current active state.
        current_state = NULL,

        #' @field state_history Active state lineage. Revert truncates this lineage.
        state_history = list(),

        #' @field audit_records Immutable audit records keyed by record ID.
        audit_records = list(),

        #' @field audit_enabled Whether peptide-QC audit capture is enabled.
        audit_enabled = TRUE,

        #' @field workflow_type Workflow route for this state manager.
        workflow_type = "LFQ",

        #' @description Initialize an empty memory-backed workflow state.
        #' @param audit_enabled Whether peptide-QC audit capture is enabled.
        initialize = function(audit_enabled = TRUE) {
            self$states <- list()
            self$state_history <- list()
            self$audit_records <- list()
            self$audit_enabled <- isTRUE(audit_enabled)
            self$current_state <- "initial"
            self$workflow_type <- "LFQ"
            private$events <- list()
            private$revision <- 0L
            private$observers <- list()
            private$nextObserverId <- 0L
            self$states[["initial"]] <- list(
                timestamp = Sys.time(),
                data = NULL,
                config = NULL,
                description = "Initial empty state",
                audit_metadata = list(
                    provenance_status = "legacy_or_not_applicable"
                )
            )
            self$state_history <- append(self$state_history, "initial")
            private$recordEvent("initialized", "initial", NULL)
        },

        #' @description Save a new state snapshot.
        #' @param state_name Logical state name.
        #' @param s4_data_object Data object to save.
        #' @param config_object Configuration used to produce the state.
        #' @param description Human-readable state description.
        #' @param audit_metadata Optional additive audit record.
        saveState = function(
            state_name,
            s4_data_object,
            config_object,
            description,
            audit_metadata = NULL
        ) {
            if (!workflowStateScalarString(state_name)) {
                stop("WorkflowState state name must be a non-empty string.")
            }
            previousState <- self$current_state
            if (isTRUE(self$audit_enabled) && is.null(audit_metadata) &&
                base::isS4(s4_data_object) &&
                "args" %in% methods::slotNames(s4_data_object) &&
                is.list(s4_data_object@args$peptide_qc_audit)) {
                auditLineage <- s4_data_object@args$peptide_qc_audit
                recordId <- auditLineage$current_record_id
                if (!is.null(recordId) && length(recordId) == 1L) {
                    matching <- Filter(
                        \(record) identical(record$record_id, recordId),
                        auditLineage$records
                    )
                    if (length(matching) > 0L) {
                        audit_metadata <- matching[[1L]]
                    }
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
            recordId <- audit_metadata$record_id
            if (!is.null(recordId)) {
                if (!workflowStateScalarString(recordId)) {
                    stop("WorkflowState audit record ID must be a non-empty string.")
                }
                existingRecord <- self$audit_records[[recordId]]
                if (!is.null(existingRecord) &&
                    !identical(existingRecord, audit_metadata)) {
                    stop("WorkflowState audit record IDs must be immutable.")
                }
            }
            self$states[[state_name]] <- list(
                timestamp = Sys.time(),
                data = s4_data_object,
                config = config_object,
                description = description,
                previous_state = previousState,
                s4_class = class(s4_data_object)[1L],
                audit_metadata = audit_metadata
            )
            if (!is.null(recordId)) {
                self$audit_records[[recordId]] <- audit_metadata
            }
            history <- self$getHistory()
            self$state_history <- as.list(c(history[history != state_name], state_name))
            self$current_state <- state_name
            private$recordEvent("state_saved", state_name, previousState)
            state_name
        },

        #' @description Retrieve data from a named or current state.
        #' @param state_name State name; defaults to the current state.
        #' @return The stored data object.
        getState = function(state_name = NULL) {
            if (is.null(state_name)) {
                state_name <- self$current_state
            }
            state <- self$states[[state_name]]
            if (is.list(state) && "data" %in% names(state)) {
                return(state$data)
            }
            state
        },

        #' @description Return metadata for a named or current state.
        #' @param state_name State name; defaults to the current state.
        #' @return State metadata without the data payload, or `NULL` if absent.
        getStateMetadata = function(state_name = NULL) {
            if (is.null(state_name)) {
                state_name <- self$current_state
            }
            state <- self$states[[state_name]]
            if (is.null(state)) {
                return(NULL)
            }
            if (!is.list(state)) {
                return(list())
            }
            state[setdiff(names(state), "data")]
        },

        #' @description Return a state's saved configuration.
        #' @param state_name State name; defaults to the current state.
        #' @return The saved configuration or `NULL`.
        getStateConfig = function(state_name = NULL) {
            self$getStateMetadata(state_name)$config
        },

        #' @description Test whether a logical state exists.
        #' @param state_name State name to test.
        #' @return A scalar logical.
        hasState = function(state_name) {
            workflowStateScalarString(state_name) && state_name %in% names(self$states)
        },

        #' @description Return the current logical state name.
        #' @return A scalar character value.
        getCurrentStateName = function() {
            self$current_state
        },

        #' @description Revert the active lineage to a saved state.
        #' @param state_name State name to activate.
        #' @return The data object from the newly active state.
        revertToState = function(state_name) {
            if (!workflowStateScalarString(state_name)) {
                stop("WorkflowState state name must be a non-empty string.")
            }
            if (state_name %in% names(self$states)) {
                previousState <- self$current_state
                history <- self$getHistory()
                revertIndex <- match(state_name, history, nomatch = 0L)
                if (revertIndex > 0L) {
                    history <- history[seq_len(revertIndex)]
                } else {
                    history <- c(history, state_name)
                }
                self$state_history <- as.list(history)
                self$current_state <- state_name
                private$recordEvent("state_reverted", state_name, previousState)
                return(self$getState(state_name))
            }
            stop("State not found: ", state_name)
        },

        #' @description Return the active state lineage.
        #' @return A character vector of state names.
        getHistory = function() {
            unlist(self$state_history, use.names = FALSE)
        },

        #' @description Return audit metadata associated with a state.
        #' @param state_name State name; defaults to the current state.
        getStateAudit = function(state_name = NULL) {
            metadata <- self$getStateMetadata(state_name)$audit_metadata
            if (is.null(metadata)) {
                return(list(provenance_status = "legacy_untracked_state"))
            }
            metadata
        },

        #' @description Return every immutable audit record.
        getAuditRecords = function() {
            self$audit_records
        },

        #' @description Return whether audit capture is enabled.
        isAuditEnabled = function() {
            isTRUE(self$audit_enabled)
        },

        #' @description Set the workflow route for this session.
        #' @param type One of the currently implemented workflow types.
        setWorkflowType = function(type) {
            if (!workflowStateScalarString(type) ||
                !type %in% .WORKFLOW_STATE_TYPES) {
                stop(
                    paste(
                        "Invalid workflow type. Must be one of:",
                        paste(.WORKFLOW_STATE_TYPES, collapse = ", ")
                    )
                )
            }
            previousType <- self$workflow_type
            self$workflow_type <- type
            private$recordEvent(
                "workflow_type_set",
                self$current_state,
                self$current_state,
                metadata = list(previous_type = previousType, type = type)
            )
            type
        },

        #' @description Return the configured workflow type.
        getWorkflowType = function() {
            self$workflow_type
        },

        #' @description Export the versioned memory-state manifest.
        #' @return A self-describing R list suitable for `restoreState()`.
        exportState = function() {
            list(
                schema = .WORKFLOW_STATE_SCHEMA,
                schema_version = .WORKFLOW_STATE_SCHEMA_VERSION,
                backend = "memory",
                workflow_type = self$workflow_type,
                audit_enabled = isTRUE(self$audit_enabled),
                current_state = self$current_state,
                active_lineage = self$getHistory(),
                states = self$states,
                audit_records = self$audit_records,
                events = private$events
            )
        },

        #' @description Export the legacy filtered-session state fields.
        #' @return A list preserving the existing session field names.
        exportLegacyState = function() {
            list(
                r6_current_state_name = self$current_state,
                r6_complete_states = self$states,
                r6_state_history = self$getHistory()
            )
        },

        #' @description Restore a validated modern or explicit legacy manifest.
        #' @param manifest State manifest or legacy filtered-session list.
        #' @param schema_version Explicit version for a legacy list. Modern
        #'   manifests carry their own version.
        restoreState = function(manifest, schema_version = NULL) {
            prepared <- private$prepareRestore(manifest, schema_version)
            self$states <- prepared$states
            self$state_history <- as.list(prepared$active_lineage)
            self$audit_records <- prepared$audit_records
            self$audit_enabled <- prepared$audit_enabled
            self$workflow_type <- prepared$workflow_type
            self$current_state <- prepared$current_state
            private$events <- prepared$events
            private$revision <- length(prepared$events)
            private$recordEvent(
                "state_restored",
                prepared$current_state,
                prepared$current_state,
                metadata = list(source_schema_version = prepared$source_version)
            )
            invisible(self)
        },

        #' @description Return immutable workflow transition events.
        #' @return A list ordered by monotonically increasing revision.
        getEvents = function() {
            private$events
        },

        #' @description Return the current transition revision.
        #' @return A non-negative integer.
        getRevision = function() {
            private$revision
        },

        #' @description Subscribe to subsequent workflow transitions.
        #' @param callback Function accepting one immutable event list.
        #' @return An idempotent function that removes the subscription.
        observeTransitions = function(callback) {
            if (!is.function(callback)) {
                stop("WorkflowState transition callback must be a function.")
            }
            private$nextObserverId <- private$nextObserverId + 1L
            observerId <- as.character(private$nextObserverId)
            private$observers[[observerId]] <- callback
            removed <- FALSE
            function() {
                if (!removed) {
                    private$observers[[observerId]] <- NULL
                    removed <<- TRUE
                }
                invisible(NULL)
            }
        }
    ),
    private = list(
        events = list(),
        revision = 0L,
        observers = list(),
        nextObserverId = 0L,

        recordEvent = function(
            eventType,
            stateName,
            previousState,
            metadata = list()
        ) {
            private$revision <- private$revision + 1L
            event <- list(
                revision = private$revision,
                event_type = eventType,
                state_name = stateName,
                previous_state = previousState,
                timestamp = Sys.time(),
                metadata = metadata
            )
            private$events[[length(private$events) + 1L]] <- event
            callbacks <- private$observers
            for (callback in callbacks) {
                tryCatch(
                    callback(event),
                    error = \(error) {
                        try(
                            logger::log_warn(paste(
                                "WorkflowState transition observer failed:",
                                conditionMessage(error)
                            )),
                            silent = TRUE
                        )
                        NULL
                    }
                )
            }
            invisible(event)
        },

        prepareRestore = function(manifest, schemaVersion) {
            if (!is.list(manifest)) {
                stop("WorkflowState manifest must be an R list.", call. = FALSE)
            }
            rawVersion <- if (is.null(schemaVersion)) {
                manifest$schema_version
            } else {
                schemaVersion
            }
            version <- workflowStateVersionValue(rawVersion)
            if (is.na(version)) {
                stop("WorkflowState manifest schema version is missing.", call. = FALSE)
            }
            if (version > .WORKFLOW_STATE_SCHEMA_VERSION) {
                stop(
                    sprintf("Unsupported future WorkflowState schema version: %d.", version),
                    call. = FALSE
                )
            }
            if (identical(version, .WORKFLOW_STATE_LEGACY_VERSION)) {
                return(private$prepareLegacyRestore(manifest))
            }
            if (!identical(version, .WORKFLOW_STATE_SCHEMA_VERSION)) {
                stop(sprintf("Unsupported WorkflowState schema version: %d.", version))
            }
            embeddedVersion <- workflowStateVersionValue(manifest$schema_version)
            if (is.na(embeddedVersion) || !identical(embeddedVersion, version)) {
                stop("WorkflowState manifest schema version is inconsistent.")
            }
            if (!identical(manifest$schema, .WORKFLOW_STATE_SCHEMA)) {
                stop("WorkflowState manifest schema identifier is invalid.", call. = FALSE)
            }
            if (!identical(manifest$backend, "memory")) {
                stop("WorkflowState manifest backend is invalid.", call. = FALSE)
            }
            private$validateRestoreValues(
                states = manifest$states,
                currentState = manifest$current_state,
                activeLineage = manifest$active_lineage,
                auditRecords = manifest$audit_records,
                auditEnabled = manifest$audit_enabled,
                workflowType = manifest$workflow_type,
                events = manifest$events,
                sourceVersion = version,
                requireRecords = TRUE
            )
        },

        prepareLegacyRestore = function(manifest) {
            states <- manifest$r6_complete_states
            if (is.null(states)) {
                states <- list()
            }
            currentState <- manifest$r6_current_state_name
            if (!is.null(currentState) && !is.null(manifest$current_s4_object) &&
                is.null(states[[currentState]])) {
                states[[currentState]] <- manifest$current_s4_object
            }
            activeLineage <- manifest$r6_state_history
            if (is.null(activeLineage)) {
                activeLineage <- names(states)
            }
            auditRecords <- private$collectLegacyAuditRecords(states)
            workflowType <- manifest$workflow_type
            if (is.null(workflowType)) {
                workflowType <- "LFQ"
            }
            private$validateRestoreValues(
                states = states,
                currentState = currentState,
                activeLineage = activeLineage,
                auditRecords = auditRecords,
                auditEnabled = TRUE,
                workflowType = workflowType,
                events = list(),
                sourceVersion = .WORKFLOW_STATE_LEGACY_VERSION,
                requireRecords = FALSE
            )
        },

        validateRestoreValues = function(
            states,
            currentState,
            activeLineage,
            auditRecords,
            auditEnabled,
            workflowType,
            events,
            sourceVersion,
            requireRecords
        ) {
            if (!is.list(states) || is.null(names(states)) ||
                any(!nzchar(names(states))) || anyDuplicated(names(states))) {
                stop("WorkflowState manifest states must be a uniquely named list.")
            }
            if (length(states) == 0L || !workflowStateScalarString(currentState) ||
                !currentState %in% names(states)) {
                stop("WorkflowState manifest current state is inconsistent.")
            }
            activeLineage <- unlist(activeLineage, use.names = FALSE)
            if (!is.character(activeLineage) || length(activeLineage) == 0L ||
                anyDuplicated(activeLineage) ||
                !identical(tail(activeLineage, 1L), currentState) ||
                any(!activeLineage %in% names(states))) {
                stop("WorkflowState manifest active lineage is inconsistent.")
            }
            if (isTRUE(requireRecords)) {
                validRecords <- vapply(
                    states,
                    \(state) {
                        required <- c(
                            "timestamp",
                            "data",
                            "config",
                            "description",
                            "audit_metadata"
                        )
                        is.list(state) && all(required %in% names(state)) &&
                            inherits(state$timestamp, "POSIXt") &&
                            is.list(state$audit_metadata)
                    },
                    logical(1)
                )
                if (!all(validRecords)) {
                    stop("WorkflowState manifest contains a malformed state record.")
                }
            }
            if (!is.list(auditRecords) || !is.logical(auditEnabled) ||
                length(auditEnabled) != 1L || is.na(auditEnabled)) {
                stop("WorkflowState manifest audit metadata is malformed.")
            }
            if (length(auditRecords) > 0L) {
                recordNames <- names(auditRecords)
                validAuditRecords <- !is.null(recordNames) &&
                    all(nzchar(recordNames)) && !anyDuplicated(recordNames) &&
                    all(vapply(seq_along(auditRecords), \(index) {
                        record <- auditRecords[[index]]
                        is.list(record) &&
                            identical(record$record_id, recordNames[[index]])
                    }, logical(1)))
                if (!validAuditRecords) {
                    stop("WorkflowState manifest audit records are inconsistent.")
                }
            }
            if (isTRUE(requireRecords)) {
                stateAuditsMatch <- all(vapply(states, \(state) {
                    recordId <- state$audit_metadata$record_id
                    is.null(recordId) ||
                        (workflowStateScalarString(recordId) &&
                            identical(auditRecords[[recordId]], state$audit_metadata))
                }, logical(1)))
                if (!stateAuditsMatch) {
                    stop("WorkflowState manifest state audits are inconsistent.")
                }
            }
            if (!workflowStateScalarString(workflowType) ||
                !workflowType %in% .WORKFLOW_STATE_TYPES) {
                stop("WorkflowState manifest workflow type is invalid.")
            }
            private$validateEvents(
                events,
                states,
                currentState,
                requireEvents = requireRecords
            )
            list(
                states = states,
                current_state = currentState,
                active_lineage = activeLineage,
                audit_records = auditRecords,
                audit_enabled = auditEnabled,
                workflow_type = workflowType,
                events = events,
                source_version = sourceVersion
            )
        },

        validateEvents = function(events, states, currentState, requireEvents) {
            if (!is.list(events)) {
                stop("WorkflowState manifest events must be a list.")
            }
            if (length(events) == 0L) {
                if (isTRUE(requireEvents)) {
                    stop("WorkflowState manifest transition events are missing.")
                }
                return(invisible(TRUE))
            }
            allowedTypes <- c(
                "initialized",
                "state_saved",
                "state_reverted",
                "workflow_type_set",
                "state_restored"
            )
            valid <- vapply(seq_along(events), \(index) {
                event <- events[[index]]
                previousState <- event$previous_state
                is.list(event) &&
                    identical(event$revision, as.integer(index)) &&
                    workflowStateScalarString(event$event_type) &&
                    event$event_type %in% allowedTypes &&
                    workflowStateScalarString(event$state_name) &&
                    event$state_name %in% names(states) &&
                    (is.null(previousState) ||
                        (workflowStateScalarString(previousState) &&
                            previousState %in% names(states))) &&
                    inherits(event$timestamp, "POSIXt") &&
                    is.list(event$metadata)
            }, logical(1))
            if (!all(valid)) {
                stop("WorkflowState manifest contains malformed transition events.")
            }
            if (!identical(tail(events, 1L)[[1L]]$state_name, currentState)) {
                stop("WorkflowState manifest transition head is inconsistent.")
            }
            invisible(TRUE)
        },

        collectLegacyAuditRecords = function(states) {
            records <- list()
            for (state in states) {
                if (!is.list(state)) {
                    next
                }
                auditMetadata <- state$audit_metadata
                if (is.list(auditMetadata) &&
                    workflowStateScalarString(auditMetadata$record_id)) {
                    records[[auditMetadata$record_id]] <- auditMetadata
                }
            }
            records
        }
    )
)
