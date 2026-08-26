# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.artifactResourceEventLimit <- 100L

artifactResourceAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_artifact_resource_error"),
        ...
    )
}

artifactResourceScalarString <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

artifactResourceOwnershipTable <- function() {
    data.frame(
        resource_type = c(
            "registry_connection", "dbi_result", "writer_guard",
            "hydration_cache", "arrow_query_object", "duckdb_query_handle",
            "project_temporary_path", "background_task", "artifact_observer",
            "artifact_workflow_state"
        ),
        creator = c(
            "projectRegistryConnect", "projectRegistryExecuteBound/FetchBound",
            "projectRegistryAcquireWriter", "ArtifactWorkflowState$getState",
            "artifact codec internals", "artifactQueryConnect",
            "projectRegistryConnect/artifactQueryConnect",
            "registerArtifactBackgroundTask", "registerArtifactObserverResource",
            "newWorkflowState"
        ),
        owner = c(
            "project registry session", "calling function", "project registry session",
            "artifact workflow state", "calling codec/query function",
            "artifact query session", "owning registry/query handle",
            "artifact resource scope", "artifact resource scope",
            "artifact resource scope"
        ),
        close_operation = c(
            "projectRegistryDisconnect", "DBI::dbClearResult",
            "projectRegistryReleaseWriter", "ArtifactWorkflowState$close",
            "convert to ordinary R and release local reference",
            "artifactQueryDisconnect", "artifactCleanupTemporaryPath",
            "cancel then join", "observer$destroy", "ArtifactWorkflowState$close"
        ),
        normal_close_hook = c(
            "ArtifactWorkflowState$close", "on.exit", "closeProjectRegistry",
            "ArtifactWorkflowState$close", "lexical function return",
            "ArtifactQuerySession$close", "owning handle close", "resource scope close",
            "resource scope close", "resource scope close"
        ),
        error_cancel_cleanup = c(
            "startup unwind/session close", "on.exit", "startup unwind/session close",
            "state close before registry close", "lexical unwind",
            "session/module close", "close after native handle release", "cancel and join",
            "destroy", "session/module close"
        ),
        process_bound = c(
            TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE
        ),
        stringsAsFactors = FALSE
    )
}

validateArtifactResourceOwnership <- function(table = artifactResourceOwnershipTable()) {
    required <- c(
        "resource_type", "creator", "owner", "close_operation",
        "normal_close_hook", "error_cancel_cleanup", "process_bound"
    )
    valid <- is.data.frame(table) && identical(names(table), required) &&
        nrow(table) > 0L && anyDuplicated(table$resource_type) == 0L &&
        all(vapply(table[setdiff(required, "process_bound")], \(column) {
            is.character(column) && !anyNA(column) && all(nzchar(column))
        }, logical(1))) &&
        is.logical(table$process_bound) && !anyNA(table$process_bound)
    if (!isTRUE(valid)) {
        artifactResourceAbort(
            "artifact resource ownership table is malformed",
            "multischolar_invalid_artifact_resource_ownership"
        )
    }
    table
}

artifactResourceAssertCreatorProcess <- function(process_id, owner) {
    if (!identical(as.integer(Sys.getpid()), as.integer(process_id))) {
        artifactResourceAbort(
            sprintf("artifact resource '%s' belongs to another R process", owner),
            "multischolar_cross_process_artifact_resource",
            creator_pid = as.integer(process_id),
            current_pid = as.integer(Sys.getpid())
        )
    }
    invisible(TRUE)
}

artifactResourceDataOnly <- function(value, owner = "worker payload") {
    unsafe_class <- inherits(value, c(
        "ArtifactResourceScope", "ArtifactWorkflowState", "ArtifactQuerySession",
        "MultiScholaRProjectRegistrySession", "DBIConnection", "DBIResult",
        "ArrowObject"
    ))
    if (unsafe_class || is.environment(value) || typeof(value) == "externalptr") {
        artifactResourceAbort(
            sprintf("%s contains a process-bound resource", owner),
            "multischolar_process_bound_worker_payload",
            owner = owner
        )
    }
    if (isS4(value)) {
        for (slot_name in methods::slotNames(value)) {
            artifactResourceDataOnly(
                methods::slot(value, slot_name),
                paste0(owner, "@", slot_name)
            )
        }
    } else if (is.list(value)) {
        value_names <- names(value)
        for (index in seq_along(value)) {
            child <- if (is.null(value_names) || !nzchar(value_names[[index]])) {
                paste0(owner, "[[", index, "]]")
            } else {
                paste0(owner, "$", value_names[[index]])
            }
            artifactResourceDataOnly(value[[index]], child)
        }
    }
    invisible(TRUE)
}

ArtifactResourceScope <- R6::R6Class(
    "ArtifactResourceScope",
    public = list(
        initialize = function(owner_id) {
            if (!artifactResourceScalarString(owner_id)) {
                artifactResourceAbort(
                    "artifact resource scope owner must be a non-empty string",
                    "multischolar_invalid_artifact_resource_owner"
                )
            }
            validateArtifactResourceOwnership()
            private$owner_id <- owner_id
            private$process_id <- as.integer(Sys.getpid())
            private$resources <- list()
            private$events <- list()
            private$next_ordinal <- 0L
            private$closed <- FALSE
            private$close_reason <- NULL
            private$cleanup_errors <- list()
        },

        register = function(resource_id, resource_type, close_fn, status_fn = NULL) {
            private$assertOpen()
            private$assertProcess()
            known_types <- validateArtifactResourceOwnership()$resource_type
            valid <- artifactResourceScalarString(resource_id) &&
                artifactResourceScalarString(resource_type) &&
                resource_type %in% known_types && is.function(close_fn) &&
                (is.null(status_fn) || is.function(status_fn))
            if (!isTRUE(valid)) {
                artifactResourceAbort(
                    "artifact resource registration is malformed",
                    "multischolar_invalid_artifact_resource_registration"
                )
            }
            if (resource_id %in% names(private$resources)) {
                artifactResourceAbort(
                    sprintf("artifact resource '%s' is already registered", resource_id),
                    "multischolar_duplicate_artifact_resource"
                )
            }
            private$next_ordinal <- private$next_ordinal + 1L
            private$resources[[resource_id]] <- list(
                resource_id = resource_id,
                resource_type = resource_type,
                close_fn = close_fn,
                status_fn = status_fn,
                ordinal = private$next_ordinal
            )
            private$recordEvent("registered", resource_id, resource_type, NULL)
            invisible(resource_id)
        },

        release = function(resource_id, reason = "explicit") {
            private$assertOpen()
            private$assertProcess()
            private$releaseEntry(resource_id, reason, raise = TRUE)
        },

        close = function(reason = "session_end") {
            private$assertProcess()
            if (isTRUE(private$closed)) return(invisible(FALSE))
            if (!artifactResourceScalarString(reason)) {
                artifactResourceAbort(
                    "artifact resource close reason must be a non-empty string",
                    "multischolar_invalid_artifact_resource_close"
                )
            }
            entries <- private$resources
            if (length(entries) > 0L) {
                ordinals <- vapply(entries, `[[`, integer(1), "ordinal")
                resource_ids <- names(entries)[order(ordinals, decreasing = TRUE)]
                for (resource_id in resource_ids) {
                    private$releaseEntry(resource_id, reason, raise = FALSE)
                }
            }
            private$resources <- list()
            private$closed <- TRUE
            private$close_reason <- reason
            errors <- private$cleanup_errors
            if (length(errors) > 0L) {
                artifactResourceAbort(
                    "one or more artifact resources failed to close",
                    "multischolar_artifact_resource_cleanup_failed",
                    cleanup_errors = unname(errors),
                    cleanup_report = self$getInfo()
                )
            }
            invisible(TRUE)
        },

        getInfo = function() {
            private$assertProcess()
            types <- if (length(private$resources) == 0L) {
                character()
            } else {
                vapply(private$resources, `[[`, character(1), "resource_type")
            }
            statuses <- lapply(private$resources, \(entry) {
                if (is.null(entry$status_fn)) return(NULL)
                tryCatch(
                    entry$status_fn(),
                    error = \(error) list(status_error = conditionMessage(error))
                )
            })
            list(
                owner_id = private$owner_id,
                creator_pid = private$process_id,
                closed = private$closed,
                close_reason = private$close_reason,
                resource_count = length(private$resources),
                resources_by_type = as.list(table(types)),
                resource_status = statuses,
                cleanup_errors = unname(private$cleanup_errors),
                events = private$events
            )
        },

        isClosed = function() isTRUE(private$closed)
    ),
    private = list(
        owner_id = NULL,
        process_id = NULL,
        resources = list(),
        events = list(),
        next_ordinal = 0L,
        closed = FALSE,
        close_reason = NULL,
        cleanup_errors = list(),

        assertProcess = function() {
            artifactResourceAssertCreatorProcess(private$process_id, private$owner_id)
        },

        assertOpen = function() {
            if (isTRUE(private$closed)) {
                artifactResourceAbort(
                    "artifact resource scope is closed",
                    "multischolar_closed_artifact_resource_scope"
                )
            }
            invisible(TRUE)
        },

        recordEvent = function(event, resource_id, resource_type, detail) {
            private$events[[length(private$events) + 1L]] <- list(
                event = event,
                resource_id = resource_id,
                resource_type = resource_type,
                detail = detail,
                occurred_at = artifactRefUtcNow()
            )
            if (length(private$events) > .artifactResourceEventLimit) {
                private$events <- tail(private$events, .artifactResourceEventLimit)
            }
            invisible(TRUE)
        },

        releaseEntry = function(resource_id, reason, raise) {
            if (!artifactResourceScalarString(resource_id)) {
                artifactResourceAbort(
                    "artifact resource ID must be a non-empty string",
                    "multischolar_invalid_artifact_resource_registration"
                )
            }
            entry <- private$resources[[resource_id]]
            if (is.null(entry)) return(invisible(FALSE))
            private$resources[[resource_id]] <- NULL
            cleanup_error <- tryCatch(
                {
                    entry$close_fn()
                    NULL
                },
                error = \(error) error
            )
            if (is.null(cleanup_error)) {
                private$recordEvent(
                    "released", resource_id, entry$resource_type, reason
                )
                return(invisible(TRUE))
            }
            detail <- conditionMessage(cleanup_error)
            private$cleanup_errors[[length(private$cleanup_errors) + 1L]] <- list(
                resource_id = resource_id,
                resource_type = entry$resource_type,
                reason = reason,
                message = detail
            )
            private$recordEvent("release_failed", resource_id, entry$resource_type, detail)
            if (isTRUE(raise)) {
                artifactResourceAbort(
                    sprintf("artifact resource '%s' failed to close", resource_id),
                    "multischolar_artifact_resource_cleanup_failed",
                    resource_id = resource_id,
                    resource_type = entry$resource_type,
                    parent = cleanup_error
                )
            }
            invisible(FALSE)
        },

        finalize = function() {
            if (!isTRUE(private$closed)) {
                try(self$close("finalizer_fallback"), silent = TRUE)
            }
        }
    )
)

registerArtifactObserverResource <- function(scope, resource_id, observer) {
    if (!inherits(scope, "ArtifactResourceScope") ||
        is.null(observer) || !is.function(observer$destroy)) {
        artifactResourceAbort(
            "artifact observer registration requires a resource scope and observer",
            "multischolar_invalid_artifact_observer"
        )
    }
    scope$register(
        resource_id,
        "artifact_observer",
        close_fn = \() observer$destroy()
    )
}

registerArtifactBackgroundTask <- function(
    scope,
    resource_id,
    task,
    cancel_fn,
    join_fn,
    worker_payload = NULL
) {
    if (!inherits(scope, "ArtifactResourceScope") ||
        !is.function(cancel_fn) || !is.function(join_fn)) {
        artifactResourceAbort(
            "artifact task registration is malformed",
            "multischolar_invalid_artifact_background_task"
        )
    }
    artifactResourceDataOnly(worker_payload)
    scope$register(
        resource_id,
        "background_task",
        close_fn = \() {
            cancel_error <- tryCatch({
                cancel_fn(task)
                NULL
            }, error = \(error) error)
            join_error <- tryCatch({
                join_fn(task)
                NULL
            }, error = \(error) error)
            errors <- Filter(Negate(is.null), list(cancel_error, join_error))
            if (length(errors) > 0L) {
                artifactResourceAbort(
                    "artifact background task cancellation or join failed",
                    "multischolar_artifact_task_cleanup_failed",
                    cleanup_messages = vapply(errors, conditionMessage, character(1))
                )
            }
            invisible(TRUE)
        }
    )
}

registerArtifactWorkflowStateResource <- function(scope, manager) {
    if (!inherits(scope, "ArtifactResourceScope") ||
        !inherits(manager, "ArtifactWorkflowState")) {
        artifactResourceAbort(
            "artifact workflow state registration is malformed",
            "multischolar_invalid_artifact_workflow_state_resource"
        )
    }
    scope$register(
        "workflow_state",
        "artifact_workflow_state",
        close_fn = \() manager$close(),
        status_fn = \() manager$getResourceInfo()
    )
}

registerArtifactQuerySessionResource <- function(
    scope,
    query_session,
    resource_id = "bounded_query_session"
) {
    if (!inherits(scope, "ArtifactResourceScope") ||
        !inherits(query_session, "ArtifactQuerySession")) {
        artifactResourceAbort(
            "artifact query session registration is malformed",
            "multischolar_invalid_artifact_query_session_resource"
        )
    }
    scope$register(
        resource_id,
        "duckdb_query_handle",
        close_fn = \() query_session$close(),
        status_fn = \() query_session$getInfo()
    )
}

artifactQuerySessionForWorkflow <- function(
    workflow_data,
    store,
    resource_policy = NULL
) {
    scope <- shiny::isolate(workflow_data$artifact_resource_scope)
    if (!inherits(scope, "ArtifactResourceScope") || scope$isClosed()) return(NULL)
    query_session <- shiny::isolate(workflow_data$artifact_query_session)
    if (inherits(query_session, "ArtifactQuerySession")) {
        query_session$assertCompatible(store, resource_policy)
        return(query_session)
    }
    query_session <- newArtifactQuerySession(store, resource_policy)
    registerArtifactQuerySessionResource(scope, query_session)
    workflow_data$artifact_query_session <- query_session
    query_session
}

artifactWorkflowStateResourceInfo <- function(
    session,
    process_id,
    closed,
    cache_generation_id,
    observer_count
) {
    artifactResourceAssertCreatorProcess(process_id, "artifact WorkflowState")
    connection_open <- !is.null(session) && !is.null(session$handle) &&
        isTRUE(tryCatch(
            DBI::dbIsValid(session$handle$connection),
            error = \(...) FALSE
        ))
    list(
        creator_pid = process_id,
        closed = closed,
        registry_connection = connection_open,
        writer_guard = !is.null(session) && !is.null(session$writer),
        hydration_cache_entries = if (is.null(cache_generation_id)) 0L else 1L,
        observer_count = as.integer(observer_count)
    )
}

artifactReleaseProcessAllocator <- function() {
    if (!identical(Sys.info()[["sysname"]], "Linux")) {
        return(invisible(FALSE))
    }
    released <- tryCatch(
        .Call(C_multischolar_malloc_trim),
        error = \(...) FALSE
    )
    invisible(isTRUE(released))
}

artifactReleaseTransientMemory <- function(full = FALSE) {
    invisible(gc(full = full))
    artifactReleaseProcessAllocator()
}

artifactScheduleTransientMemoryRelease <- function(delay = 0) {
    if (!requireNamespace("later", quietly = TRUE)) {
        return(invisible(FALSE))
    }
    later::later(
        \() artifactReleaseTransientMemory(),
        delay = delay
    )
    invisible(TRUE)
}

closeArtifactWorkflowStateSession <- function(session) {
    if (is.null(session)) return(list(closed = TRUE, error = NULL))
    cleanup_error <- tryCatch(
        {
            closeProjectRegistry(session)
            NULL
        },
        error = \(error) error
    )
    list(closed = isTRUE(session$closed), error = cleanup_error)
}

closeWorkflowArtifactLifecycle <- function(
    workflow_data,
    reason = "module_removed",
    strict = TRUE
) {
    scope <- shiny::isolate(workflow_data$artifact_resource_scope)
    if (!inherits(scope, "ArtifactResourceScope")) return(invisible(FALSE))
    cleanup_error <- tryCatch(
        {
            scope$close(reason)
            NULL
        },
        error = \(error) error
    )
    workflow_data$artifact_cleanup_report <- scope$getInfo()
    manager <- shiny::isolate(workflow_data$state_manager)
    if (inherits(manager, "ArtifactWorkflowState")) {
        workflow_data$state_manager <- NULL
    }
    workflow_data$artifact_resource_scope <- NULL
    workflow_data$artifact_query_session <- NULL
    if (!is.null(cleanup_error) && isTRUE(strict)) stop(cleanup_error)
    if (!is.null(cleanup_error)) {
        logger::log_error(paste(
            "Artifact session resource cleanup failed:",
            conditionMessage(cleanup_error)
        ))
    }
    invisible(is.null(cleanup_error))
}

registerWorkflowArtifactLifecycle <- function(
    workflow_data,
    session,
    context_observer = NULL
) {
    if (!is.function(session$onSessionEnded)) {
        artifactResourceAbort(
            "artifact lifecycle requires a Shiny session close hook",
            "multischolar_invalid_artifact_session"
        )
    }
    context <- workflow_data$workflow_context
    owner_id <- if (inherits(context, "WorkflowContext")) {
        context$getStaticIdentity()$workflow_id
    } else {
        "workflow.gui"
    }
    scope <- ArtifactResourceScope$new(owner_id)
    if (!is.null(context_observer)) {
        registerArtifactObserverResource(scope, "context_binding", context_observer)
    }
    tracked_manager <- NULL
    state_observer <- shiny::observe({
        manager <- workflow_data$state_manager
        if (identical(manager, tracked_manager)) return(invisible(NULL))
        if (inherits(tracked_manager, "ArtifactWorkflowState")) {
            scope$release("workflow_state", "state_manager_replaced")
        }
        tracked_manager <<- manager
        if (inherits(manager, "ArtifactWorkflowState")) {
            registerArtifactWorkflowStateResource(scope, manager)
        }
        invisible(NULL)
    })
    registerArtifactObserverResource(scope, "state_manager_tracking", state_observer)
    workflow_data$artifact_resource_scope <- scope
    workflow_data$artifact_query_session <- NULL
    workflow_data$artifact_cleanup_report <- NULL
    session$onSessionEnded(\() {
        closeWorkflowArtifactLifecycle(
            workflow_data,
            reason = "session_end",
            strict = FALSE
        )
    })
    scope
}
