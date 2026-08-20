# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

normalizeWorkflowStoragePolicy <- function(storage_policy = NULL) {
    if (is.null(storage_policy)) storage_policy <- list()
    if (!is.list(storage_policy)) {
        workflowCapabilityAbort(
            "workflow storage policy must be an R list",
            "multischolar_invalid_storage_policy"
        )
    }
    requested_backend <- storage_policy$requested_backend
    if (is.null(requested_backend)) requested_backend <- "auto"
    requested_rollout <- storage_policy$requested_rollout
    if (is.null(requested_rollout)) requested_rollout <- "dual_write"
    migration_requested <- storage_policy$migration_requested
    if (is.null(migration_requested)) migration_requested <- FALSE
    list(
        requested_backend = requested_backend,
        requested_rollout = requested_rollout,
        migration_requested = migration_requested,
        project_id = storage_policy$project_id
    )
}

workflowContextBindingSnapshot <- function(identity, resolution, paths) {
    list(
        identity = identity,
        resolution = unclass(resolution),
        path_metadata = if (is.null(paths)) NULL else artifactPathMetadata(paths)
    )
}

workflowNormalizeProjectRoot <- function(project_root) {
    if (!workflowCapabilityScalarString(project_root)) {
        workflowCapabilityAbort(
            "workflow project root must be one path",
            "multischolar_invalid_project_root"
        )
    }
    expanded <- path.expand(project_root)
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", expanded)) {
        expanded <- file.path(getwd(), expanded)
    }
    normalizePath(expanded, winslash = "/", mustWork = FALSE)
}

workflowProjectRootFromPaths <- function(experiment_paths) {
    if (workflowCapabilityScalarString(experiment_paths$base_dir)) {
        return(experiment_paths$base_dir)
    }
    legacy_keys <- c("data_dir", "results_dir", "source_dir")
    available <- Filter(
        workflowCapabilityScalarString,
        experiment_paths[intersect(legacy_keys, names(experiment_paths))]
    )
    if (length(available) == 0L) {
        workflowCapabilityAbort(
            "workflow paths require base_dir or a standard legacy directory",
            "multischolar_invalid_workflow_context_input"
        )
    }
    roots <- vapply(available, \(path) dirname(dirname(path)), character(1))
    normalized <- vapply(roots, workflowNormalizeProjectRoot, character(1))
    if (length(unique(normalized)) != 1L) {
        workflowCapabilityAbort(
            "legacy workflow directories do not share one project root",
            "multischolar_ambiguous_project_root"
        )
    }
    normalized[[1L]]
}

WorkflowContext <- R6::R6Class(
    "WorkflowContext",
    public = list(
        initialize = function(project_root, static_identity, storage_policy = NULL) {
            required <- c("project_id", "omic_type", "omic_label", "workflow_id")
            if (!all(required %in% names(static_identity)) ||
                !all(vapply(
                    static_identity[required],
                    workflowCapabilityScalarString,
                    logical(1)
                ))) {
                workflowCapabilityAbort(
                    "workflow context static identity is incomplete",
                    "multischolar_invalid_workflow_identity"
                )
            }
            private$project_root <- workflowNormalizeProjectRoot(project_root)
            private$static_identity <- static_identity[required]
            private$storage_policy <- normalizeWorkflowStoragePolicy(storage_policy)
            private$binding <- NULL
            private$observers <- list()
            private$next_observer_id <- 0L
        },

        isBound = function() {
            !is.null(private$binding)
        },

        getProjectRoot = function() {
            private$project_root
        },

        getStaticIdentity = function() {
            private$static_identity
        },

        getStoragePolicy = function() {
            private$storage_policy
        },

        getIdentity = function() {
            if (is.null(private$binding)) return(NULL)
            private$binding$identity
        },

        getStorageDecision = function() {
            if (is.null(private$binding)) return(NULL)
            private$binding$resolution
        },

        getPaths = function() {
            if (is.null(private$binding)) return(NULL)
            private$binding$paths
        },

        getSnapshot = function() {
            if (is.null(private$binding)) return(NULL)
            workflowContextBindingSnapshot(
                private$binding$identity,
                private$binding$resolution,
                private$binding$paths
            )
        },

        bind = function(identity, resolution, paths = NULL) {
            if (!all(.WORKFLOW_IDENTITY_FIELDS %in% names(identity)) ||
                !all(vapply(
                    identity[.WORKFLOW_IDENTITY_FIELDS],
                    workflowCapabilityScalarString,
                    logical(1)
                ))) {
                workflowCapabilityAbort(
                    "workflow context binding identity is incomplete",
                    "multischolar_invalid_workflow_identity"
                )
            }
            static_fields <- names(private$static_identity)
            if (!identical(
                identity[static_fields],
                private$static_identity[static_fields]
            )) {
                workflowCapabilityAbort(
                    "workflow context binding changes its static identity",
                    "multischolar_workflow_context_identity_mismatch"
                )
            }
            required_resolution <- c(
                "requested_backend",
                "effective_backend",
                "effective_rollout",
                "capability_id",
                "capability_version",
                "reason_code",
                "project_state"
            )
            if (!is.list(resolution) ||
                !all(required_resolution %in% names(resolution)) ||
                !resolution$effective_backend %in% c("memory", "artifact")) {
                workflowCapabilityAbort(
                    "workflow context storage resolution is invalid",
                    "multischolar_invalid_backend_resolution"
                )
            }
            if (identical(resolution$effective_backend, "memory") &&
                !is.null(paths)) {
                workflowCapabilityAbort(
                    "memory workflow context cannot own artifact paths",
                    "multischolar_memory_context_has_artifact_paths"
                )
            }
            if (identical(resolution$effective_backend, "artifact") &&
                !inherits(paths, "MultiScholaRArtifactPaths")) {
                workflowCapabilityAbort(
                    "artifact workflow context requires validated paths",
                    "multischolar_artifact_context_missing_paths"
                )
            }
            candidate <- list(identity = identity, resolution = resolution, paths = paths)
            if (!is.null(private$binding)) {
                if (identical(
                    workflowContextBindingSnapshot(
                        candidate$identity,
                        candidate$resolution,
                        candidate$paths
                    ),
                    self$getSnapshot()
                )) {
                    return(invisible(FALSE))
                }
                workflowCapabilityAbort(
                    "workflow context is immutable and is already bound",
                    "multischolar_workflow_context_already_bound"
                )
            }
            private$binding <- candidate
            callbacks <- private$observers
            snapshot <- self$getSnapshot()
            for (callback in callbacks) private$notifyObserver(callback, snapshot)
            invisible(TRUE)
        },

        observeBinding = function(callback, replay = TRUE) {
            if (!is.function(callback)) {
                workflowCapabilityAbort(
                    "workflow context observer must be a function",
                    "multischolar_invalid_workflow_context_observer"
                )
            }
            private$next_observer_id <- private$next_observer_id + 1L
            observer_id <- as.character(private$next_observer_id)
            private$observers[[observer_id]] <- callback
            if (isTRUE(replay) && self$isBound()) {
                private$notifyObserver(callback, self$getSnapshot())
            }
            removed <- FALSE
            function() {
                if (!removed) {
                    private$observers[[observer_id]] <- NULL
                    removed <<- TRUE
                }
                invisible(NULL)
            }
        }
    ),
    private = list(
        project_root = NULL,
        static_identity = NULL,
        storage_policy = NULL,
        binding = NULL,
        observers = list(),
        next_observer_id = 0L,

        notifyObserver = function(callback, snapshot) {
            tryCatch(
                callback(snapshot),
                error = \(error) {
                    try(
                        logger::log_warn(paste(
                            "WorkflowContext binding observer failed:",
                            conditionMessage(error)
                        )),
                        silent = TRUE
                    )
                    NULL
                }
            )
        }
    )
)

createWorkflowContext <- function(
    experiment_paths,
    omic_type,
    experiment_label = NULL,
    storage_policy = NULL
) {
    if (!is.list(experiment_paths) ||
        !workflowCapabilityScalarString(omic_type)) {
        workflowCapabilityAbort(
            "workflow context requires existing experiment paths and omic type",
            "multischolar_invalid_workflow_context_input"
        )
    }
    policy <- normalizeWorkflowStoragePolicy(storage_policy)
    root <- workflowNormalizeProjectRoot(
        workflowProjectRootFromPaths(experiment_paths)
    )
    project_id <- policy$project_id
    if (is.null(project_id)) project_id <- experiment_paths$project_id
    if (is.null(project_id)) project_id <- basename(root)
    omic_label <- experiment_paths$omic_label
    if (is.null(omic_label)) omic_label <- experiment_label
    if (is.null(omic_label)) omic_label <- omic_type
    context_omic_type <- experiment_paths$omic_type
    if (is.null(context_omic_type)) context_omic_type <- omic_type
    static_identity <- list(
        project_id = project_id,
        omic_type = context_omic_type,
        omic_label = omic_label,
        workflow_id = paste0(omic_type, ".gui")
    )
    WorkflowContext$new(
        project_root = root,
        static_identity = static_identity,
        storage_policy = policy
    )
}

bindWorkflowContextFromImport <- function(
    context,
    workflow_type,
    input_format,
    data_level,
    capabilities = NULL,
    path_builder = buildArtifactPaths,
    project_state_detector = detectWorkflowProjectState,
    descriptor_catalogue = NULL
) {
    if (!inherits(context, "WorkflowContext")) {
        workflowCapabilityAbort(
            "workflow context binding requires a WorkflowContext",
            "multischolar_invalid_workflow_context"
        )
    }
    if (!is.null(descriptor_catalogue)) {
        if (!is.null(capabilities)) {
            workflowCapabilityAbort(
                "workflow binding cannot combine capability and descriptor catalogues",
                "multischolar_ambiguous_workflow_capability"
            )
        }
        capabilities <- mergeWorkflowDescriptorCapabilities(
            descriptor_catalogue = descriptor_catalogue
        )
    }
    identity <- resolveImportedWorkflowIdentity(
        context$getStaticIdentity(),
        workflow_type,
        input_format,
        data_level,
        capabilities = capabilities
    )
    policy <- context$getStoragePolicy()
    if (identical(policy$requested_backend, "artifact")) {
        resolveWorkflowBackend(
            identity,
            requested_backend = policy$requested_backend,
            requested_rollout = policy$requested_rollout,
            project_state = "new",
            migration_requested = TRUE,
            capabilities = capabilities
        )
    }
    project_state <- if (identical(policy$requested_backend, "memory")) {
        list(status = "new", reason = "explicit_memory_without_artifact_probe")
    } else {
        project_state_detector(context$getProjectRoot(), identity)
    }
    resolution <- resolveWorkflowBackend(
        identity,
        requested_backend = policy$requested_backend,
        requested_rollout = policy$requested_rollout,
        project_state = project_state,
        migration_requested = policy$migration_requested,
        capabilities = capabilities
    )
    paths <- if (identical(resolution$effective_backend, "artifact")) {
        path_builder(context$getProjectRoot(), identity)
    } else {
        NULL
    }
    context$bind(identity, resolution, paths)
    invisible(context$getSnapshot())
}

registerWorkflowContextBindingObserver <- function(
    workflow_data,
    session,
    register_session_cleanup = TRUE,
    descriptor_catalogue_fn = artifactWorkflowDescriptorCatalogue
) {
    if (!is.logical(register_session_cleanup) ||
        length(register_session_cleanup) != 1L ||
        is.na(register_session_cleanup)) {
        workflowCapabilityAbort(
            "workflow context cleanup registration must be true or false",
            "multischolar_invalid_workflow_context_cleanup"
        )
    }
    context <- workflow_data$workflow_context
    if (!inherits(context, "WorkflowContext")) return(NULL)
    binding_inputs <- shiny::reactive({
        shiny::req(workflow_data$data_format, workflow_data$data_type)
        list(
            input_format = workflow_data$data_format,
            data_level = workflow_data$data_type
        )
    })
    observer <- shiny::observeEvent(
        binding_inputs(),
        {
            bindWorkflowContextFromImport(
                context,
                workflow_type = workflowStateType(workflow_data$state_manager),
                input_format = workflow_data$data_format,
                data_level = workflow_data$data_type,
                descriptor_catalogue = descriptor_catalogue_fn()
            )
        },
        ignoreNULL = TRUE,
        once = TRUE
    )
    if (isTRUE(register_session_cleanup)) {
        session$onSessionEnded(\() observer$destroy())
    }
    observer
}
