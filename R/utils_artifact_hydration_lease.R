# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

artifactHydrationLeaseAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_artifact_hydration_lease_error"),
        ...
    )
}

artifactHydrationLeaseTokenFields <- function() {
    c(
        "schema", "schema_version", "lease_id", "creator_pid", "project_id",
        "omic_type", "descriptor_id", "generation_id", "stage_id", "roles",
        "ref_digests", "state", "validated", "published"
    )
}

artifactHydrationLeaseNewToken <- function(
    project_id,
    omic_type,
    descriptor_id,
    generation_id,
    stage_id,
    refs
) {
    roles <- names(refs)
    ref_digests <- vapply(refs, function(ref) {
        digest <- ref$semantic_digest
        if (is.null(digest)) digest <- ref$digest
        digest
    }, character(1))
    valid <- all(vapply(
        c(project_id, omic_type, descriptor_id, generation_id, stage_id),
        workflowPolicyScalarString,
        logical(1)
    )) && is.list(refs) && length(refs) > 0L &&
        !is.null(roles) && all(nzchar(roles)) && !anyDuplicated(roles) &&
        all(vapply(ref_digests, workflowPolicyDigestValid, logical(1)))
    if (!valid) {
        artifactHydrationLeaseAbort(
            "artifact hydration lease identity or refs are invalid",
            "multischolar_invalid_artifact_hydration_lease"
        )
    }
    basis <- list(
        schema = "multischolar.artifact_hydration_lease",
        schema_version = "1.0.0",
        creator_pid = as.integer(Sys.getpid()),
        project_id = project_id,
        omic_type = omic_type,
        descriptor_id = descriptor_id,
        generation_id = generation_id,
        stage_id = stage_id,
        roles = as.list(roles),
        ref_digests = as.list(unname(ref_digests)),
        state = "active",
        validated = FALSE,
        published = FALSE
    )
    lease_id <- workflowPolicyObjectDigest(basis)
    token <- append(basis, list(lease_id = lease_id), after = 2L)
    token <- token[artifactHydrationLeaseTokenFields()]
    structure(token, class = "ArtifactHydrationLeaseToken")
}

validateArtifactHydrationLeaseToken <- function(token) {
    valid <- inherits(token, "ArtifactHydrationLeaseToken") &&
        is.list(token) && identical(
            names(token),
            artifactHydrationLeaseTokenFields()
        ) && workflowPolicyDigestValid(token$lease_id) &&
        identical(token$creator_pid, as.integer(Sys.getpid())) &&
        token$state %in% c("active", "validated", "published") &&
        identical(token$validated, token$state %in% c("validated", "published")) &&
        identical(token$published, identical(token$state, "published"))
    if (!valid) {
        artifactHydrationLeaseAbort(
            "artifact hydration lease token is invalid or process-stale",
            "multischolar_invalid_artifact_hydration_lease"
        )
    }
    invisible(token)
}

ArtifactHydrationLeaseManager <- R6::R6Class(
    "ArtifactHydrationLeaseManager",
    public = list(
        initialize = function(project_id, omic_type, descriptor_id, scope = NULL) {
            identity <- c(project_id, omic_type, descriptor_id)
            if (!all(vapply(
                identity,
                workflowPolicyScalarString,
                logical(1)
            ))) {
                artifactHydrationLeaseAbort(
                    "artifact hydration lease manager identity is invalid",
                    "multischolar_invalid_artifact_hydration_lease"
                )
            }
            if (is.null(scope)) {
                scope <- ArtifactResourceScope$new(
                    paste0("hydration-lease-", project_id)
                )
            }
            if (!inherits(scope, "ArtifactResourceScope")) {
                artifactHydrationLeaseAbort(
                    "artifact hydration lease requires a resource scope",
                    "multischolar_invalid_artifact_hydration_lease"
                )
            }
            private$project_id <- project_id
            private$omic_type <- omic_type
            private$descriptor_id <- descriptor_id
            private$creator_pid <- as.integer(Sys.getpid())
            private$scope <- scope
            private$active <- NULL
            private$closed <- FALSE
        },

        acquire = function(generation_id, stage_id, refs, release_fn) {
            private$assertUsable()
            if (!is.null(private$active)) {
                artifactHydrationLeaseAbort(
                    "only one complete artifact stage may be hydrated",
                    "multischolar_overlapping_artifact_hydration"
                )
            }
            if (!is.function(release_fn)) {
                artifactHydrationLeaseAbort(
                    "artifact hydration lease release callback is invalid",
                    "multischolar_invalid_artifact_hydration_lease"
                )
            }
            token <- artifactHydrationLeaseNewToken(
                private$project_id,
                private$omic_type,
                private$descriptor_id,
                generation_id,
                stage_id,
                refs
            )
            resource_id <- paste0("hydration-lease-", token$lease_id)
            private$active <- list(
                token = token,
                resource_id = resource_id,
                release_fn = release_fn
            )
            private$scope$register(
                resource_id,
                "hydration_cache",
                function() {
                    active <- private$active
                    if (!is.null(active) &&
                        identical(active$resource_id, resource_id)) {
                        on.exit(private$active <- NULL, add = TRUE)
                        active$release_fn()
                    }
                    invisible(TRUE)
                }
            )
            token
        },

        markValidated = function(token) {
            private$assertActiveToken(token, "active")
            token$state <- "validated"
            token$validated <- TRUE
            private$active$token <- token
            token
        },

        markPublished = function(token) {
            private$assertActiveToken(token, "validated")
            token$state <- "published"
            token$published <- TRUE
            private$active$token <- token
            token
        },

        release = function(token, reason = "stage_complete") {
            private$assertUsable()
            private$assertActiveToken(token, token$state)
            resource_id <- private$active$resource_id
            released <- private$scope$release(resource_id, reason)
            if (!isTRUE(released) || !is.null(private$active)) {
                artifactHydrationLeaseAbort(
                    "artifact hydration lease release is incomplete",
                    "multischolar_incomplete_artifact_hydration_release"
                )
            }
            invisible(TRUE)
        },

        isActive = function() {
            !is.null(private$active)
        },

        getInfo = function() {
            list(
                project_id = private$project_id,
                omic_type = private$omic_type,
                descriptor_id = private$descriptor_id,
                creator_pid = private$creator_pid,
                active = !is.null(private$active),
                active_token = if (is.null(private$active)) {
                    NULL
                } else {
                    private$active$token
                },
                resource_count = private$scope$getInfo()$resource_count,
                closed = private$closed
            )
        },

        close = function(reason = "manager_close") {
            if (private$closed) return(invisible(FALSE))
            on.exit({
                private$active <- NULL
                private$closed <- TRUE
            }, add = TRUE)
            private$scope$close(reason)
            invisible(TRUE)
        }
    ),
    private = list(
        project_id = NULL,
        omic_type = NULL,
        descriptor_id = NULL,
        creator_pid = NULL,
        scope = NULL,
        active = NULL,
        closed = FALSE,

        assertUsable = function() {
            artifactResourceAssertCreatorProcess(
                private$creator_pid,
                "artifact hydration lease manager"
            )
            if (private$closed) {
                artifactHydrationLeaseAbort(
                    "artifact hydration lease manager is closed",
                    "multischolar_closed_artifact_hydration_lease"
                )
            }
            invisible(TRUE)
        },

        assertActiveToken = function(token, expected_state) {
            private$assertUsable()
            validateArtifactHydrationLeaseToken(token)
            active <- private$active
            valid <- !is.null(active) &&
                identical(token$lease_id, active$token$lease_id) &&
                identical(token$state, active$token$state) &&
                identical(token$state, expected_state) &&
                identical(token$project_id, private$project_id) &&
                identical(token$omic_type, private$omic_type) &&
                identical(token$descriptor_id, private$descriptor_id) &&
                identical(
                    token$generation_id,
                    active$token$generation_id
                ) && identical(token$stage_id, active$token$stage_id) &&
                identical(token$roles, active$token$roles) &&
                identical(token$ref_digests, active$token$ref_digests)
            if (!valid) {
                artifactHydrationLeaseAbort(
                    "artifact hydration lease token is stale or mismatched",
                    "multischolar_stale_artifact_hydration_lease"
                )
            }
            invisible(TRUE)
        }
    )
)
