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

.WORKFLOW_STATE_SCHEMA <- "multischolar.workflow_state"
.WORKFLOW_STATE_SCHEMA_VERSION <- 1L
.WORKFLOW_STATE_LEGACY_VERSION <- 0L
.WORKFLOW_STATE_TYPES <- c(
    "LFQ",
    "TMT",
    "DIA",
    "metabolomics_standard",
    "lipidomics_standard"
)

workflowStateScalarString <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

workflowStateVersionValue <- function(value) {
    version <- suppressWarnings(as.integer(value))
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || length(version) != 1L || is.na(version) ||
        value != version) {
        return(NA_integer_)
    }
    version
}

workflowStateMember <- function(stateManager, member) {
    tryCatch(stateManager[[member]], error = \(error) NULL)
}

# These adapters keep legacy list/environment test doubles usable while all real
# WorkflowState instances pass through the validated owner methods.
workflowStateCurrentName <- function(stateManager) {
    getter <- workflowStateMember(stateManager, "getCurrentStateName")
    if (is.function(getter)) {
        return(getter())
    }
    workflowStateMember(stateManager, "current_state")
}

workflowStateType <- function(stateManager) {
    getter <- workflowStateMember(stateManager, "getWorkflowType")
    if (is.function(getter)) {
        return(getter())
    }
    workflowStateMember(stateManager, "workflow_type")
}

workflowStateAuditEnabled <- function(stateManager) {
    getter <- workflowStateMember(stateManager, "isAuditEnabled")
    if (is.function(getter)) {
        return(getter())
    }
    value <- workflowStateMember(stateManager, "audit_enabled")
    if (is.null(value)) {
        return(NULL)
    }
    isTRUE(value)
}

workflowStateMetadata <- function(stateManager, stateName = NULL) {
    getter <- workflowStateMember(stateManager, "getStateMetadata")
    if (is.function(getter)) {
        return(getter(stateName))
    }
    if (is.null(stateName)) {
        stateName <- workflowStateCurrentName(stateManager)
    }
    states <- workflowStateMember(stateManager, "states")
    state <- states[[stateName]]
    if (!is.list(state)) {
        return(list())
    }
    state[setdiff(names(state), "data")]
}

workflowStateHasState <- function(stateManager, stateName) {
    checker <- workflowStateMember(stateManager, "hasState")
    if (is.function(checker)) {
        return(checker(stateName))
    }
    stateName %in% names(workflowStateMember(stateManager, "states"))
}

workflowStateHistory <- function(stateManager) {
    getter <- workflowStateMember(stateManager, "getHistory")
    if (is.function(getter)) {
        return(as.character(getter()))
    }
    history <- workflowStateMember(stateManager, "state_history")
    if (is.null(history)) {
        return(NULL)
    }
    as.character(unlist(history, use.names = FALSE))
}

workflowStateNames <- function(stateManager) {
    names(workflowStateMember(stateManager, "states"))
}

workflowStateAuditRecords <- function(stateManager) {
    getter <- workflowStateMember(stateManager, "getAuditRecords")
    if (is.function(getter)) {
        return(getter())
    }
    workflowStateMember(stateManager, "audit_records")
}

workflowStateLegacySnapshot <- function(stateManager) {
    exporter <- workflowStateMember(stateManager, "exportLegacyState")
    if (is.function(exporter)) {
        return(exporter())
    }
    list(
        r6_current_state_name = workflowStateCurrentName(stateManager),
        r6_complete_states = workflowStateMember(stateManager, "states"),
        r6_state_history = workflowStateHistory(stateManager)
    )
}

workflowStateManifest <- function(stateManager) {
    exporter <- workflowStateMember(stateManager, "exportState")
    if (!is.function(exporter)) {
        return(NULL)
    }
    exporter()
}

restoreWorkflowStateFromSession <- function(stateManager, sessionData) {
    restorer <- workflowStateMember(stateManager, "restoreState")
    manifest <- sessionData$workflow_state_manifest
    if (is.function(restorer)) {
        if (!is.null(manifest)) {
            restorer(manifest)
        } else {
            restorer(sessionData, schema_version = .WORKFLOW_STATE_LEGACY_VERSION)
        }
        return(workflowStateLegacySnapshot(stateManager))
    }
    snapshot <- list(
        r6_current_state_name = sessionData$r6_current_state_name,
        r6_complete_states = sessionData$r6_complete_states,
        r6_state_history = sessionData$r6_state_history
    )
    if (is.null(snapshot$r6_complete_states)) {
        snapshot$r6_complete_states <- list()
    }
    currentName <- snapshot$r6_current_state_name
    if (!is.null(currentName) && !is.null(sessionData$current_s4_object) &&
        is.null(snapshot$r6_complete_states[[currentName]])) {
        snapshot$r6_complete_states[[currentName]] <- sessionData$current_s4_object
    }
    if (is.null(snapshot$r6_state_history)) {
        snapshot$r6_state_history <- names(snapshot$r6_complete_states)
    }
    stateManager$states <- snapshot$r6_complete_states
    stateManager$state_history <- snapshot$r6_state_history
    stateManager$current_state <- currentName
    snapshot
}

observeWorkflowStateTransitions <- function(stateManager, callback) {
    observer <- workflowStateMember(stateManager, "observeTransitions")
    if (!is.function(observer)) {
        return(function() invisible(NULL))
    }
    observer(callback)
}
