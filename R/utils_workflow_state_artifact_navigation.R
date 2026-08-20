# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

artifactWorkflowStateRowMatches <- function(row, state_name) {
    identical(row$logical_name, state_name)
}

artifactWorkflowStateSelectGeneration <- function(
    state_rows,
    active_lineage,
    state_name
) {
    active_matches <- which(vapply(
        active_lineage,
        artifactWorkflowStateRowMatches,
        logical(1),
        state_name = state_name
    ))
    current <- tail(active_lineage, 1L)[[1L]]
    if (length(active_matches) > 0L) {
        target_index <- tail(active_matches, 1L)
        target <- active_lineage[[target_index]]
        mode <- if (identical(target$generation_id, current$generation_id)) {
            "current"
        } else {
            "revert"
        }
        descendants <- if (target_index < length(active_lineage)) {
            active_lineage[seq.int(target_index + 1L, length(active_lineage))]
        } else {
            list()
        }
        return(list(
            mode = mode,
            target = target,
            previous = current,
            descendants = descendants
        ))
    }

    matching_rows <- Filter(
        \(row) artifactWorkflowStateRowMatches(row, state_name),
        state_rows
    )
    resumable <- Filter(function(row) {
        identical(row$status, "stale") &&
            identical(row$parent_generation_id, current$generation_id)
    }, matching_rows)
    if (length(resumable) > 0L) {
        return(list(
            mode = "resume",
            target = tail(resumable, 1L)[[1L]],
            previous = current,
            descendants = list()
        ))
    }
    if (length(matching_rows) == 0L) {
        stop("State not found: ", state_name)
    }

    target <- tail(matching_rows, 1L)[[1L]]
    artifactWorkflowStateAbort(
        sprintf(
            "artifact state '%s' cannot resume from the current generation",
            state_name
        ),
        "multischolar_artifact_state_resume_parent_required",
        generation_id = target$generation_id,
        required_parent_generation_id = target$parent_generation_id,
        current_generation_id = current$generation_id,
        remediation = paste(
            "Revert to the stale generation's recorded parent before",
            "resuming this state."
        )
    )
}

artifactWorkflowStatePreflightSelection <- function(request) {
    tryCatch(
        artifactWorkflowStateHydrateData(
            request$store,
            request$manifest,
            hydrate_fn = request$hydrate_fn
        ),
        error = \(error) artifactWorkflowStateAbort(
            sprintf(
                "artifact state '%s' could not be hydrated before activation",
                request$target$logical_name
            ),
            "multischolar_artifact_state_selection_hydration_failed",
            project_id = request$identity$project_id,
            workflow_id = request$identity$workflow_id,
            generation_id = request$target$generation_id,
            state_name = request$target$logical_name,
            remediation = paste(
                "Verify the generation manifest and immutable payloads, then",
                "retry without changing the current state."
            ),
            parent = error
        )
    )
}

artifactWorkflowStateGenerationIds <- function(rows) {
    if (length(rows) == 0L) return(character())
    vapply(rows, `[[`, character(1), "generation_id")
}

artifactWorkflowStateCommitRevert <- function(request) {
    timestamp <- artifactRefUtcNow()
    stale_generation_ids <- artifactWorkflowStateGenerationIds(
        request$selection$descendants
    )
    artifactWorkflowStateTransaction(request$session, function() {
        for (row in rev(request$selection$descendants)) {
            artifactWorkflowStateUpdateStatus(
                request$session,
                request$identity,
                row$generation_id,
                row$status,
                "stale",
                timestamp
            )
        }
        artifactWorkflowStateUpdateStatus(
            request$session,
            request$identity,
            request$selection$target$generation_id,
            request$selection$target$status,
            "current",
            timestamp
        )
        request$write_revision(
            request$action_id,
            request$selection$target$generation_id,
            request$selection$previous$generation_id,
            "reverted",
            list(
                state_name = request$selection$target$logical_name,
                stale_generation_ids = stale_generation_ids
            )
        )
        request$write_event(
            "state_reverted",
            "accepted",
            request$selection$target$generation_id,
            request$selection$target$logical_name,
            request$selection$previous$logical_name,
            list(stale_generation_ids = stale_generation_ids),
            timestamp
        )
    })
    invisible(TRUE)
}

artifactWorkflowStateCommitResume <- function(request) {
    timestamp <- artifactRefUtcNow()
    artifactWorkflowStateTransaction(request$session, function() {
        artifactWorkflowStateUpdateStatus(
            request$session,
            request$identity,
            request$selection$previous$generation_id,
            request$selection$previous$status,
            "historical",
            timestamp
        )
        artifactWorkflowStateUpdateStatus(
            request$session,
            request$identity,
            request$selection$target$generation_id,
            request$selection$target$status,
            "current",
            timestamp
        )
        request$write_revision(
            request$action_id,
            request$selection$target$generation_id,
            request$selection$previous$generation_id,
            "resumed",
            list(
                state_name = request$selection$target$logical_name,
                parent_generation_id = request$selection$target$parent_generation_id
            )
        )
        request$write_event(
            "state_resumed",
            "accepted",
            request$selection$target$generation_id,
            request$selection$target$logical_name,
            request$selection$previous$logical_name,
            list(
                parent_generation_id = request$selection$target$parent_generation_id
            ),
            timestamp
        )
    })
    invisible(TRUE)
}
