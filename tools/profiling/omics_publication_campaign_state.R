publicationCampaignAbort <- function(message) {
    publicationAbort(message, "multischolar_publication_campaign_error")
}

publicationCampaignUsageNames <- function() {
    c(
        "process_launches", "worker_wall_seconds", "worker_cpu_seconds",
        "maximum_temporary_bytes_observed", "generated_payload_bytes_retained",
        "raw_evidence_bytes_retained", "sanitized_publication_bytes"
    )
}

publicationCampaignEmptyUsage <- function() {
    stats::setNames(
        as.list(rep(0, length(publicationCampaignUsageNames()))),
        publicationCampaignUsageNames()
    )
}

publicationCampaignAttemptUsageNames <- function() {
    c(
        "worker_wall_seconds", "worker_cpu_seconds", "temporary_bytes",
        "generated_payload_bytes_retained", "raw_evidence_bytes_retained",
        "sanitized_publication_bytes"
    )
}

publicationCampaignStateDigest <- function(state) {
    candidate <- state
    candidate$state_digest <- NULL
    publicationObjectDigest(candidate)
}

publicationCampaignEventDigest <- function(event) {
    candidate <- event
    candidate$event_digest <- NULL
    publicationObjectDigest(candidate)
}

publicationCampaignEvent <- function(
    event_type,
    event_index,
    previous_event_digest = NULL,
    slot_id = NULL,
    details = list()
) {
    event <- list(
        event_index = as.integer(event_index),
        event_type = event_type,
        occurred_at = publicationUtcNow(),
        previous_event_digest = previous_event_digest,
        slot_id = slot_id,
        details = details,
        event_digest = NULL
    )
    event$event_digest <- publicationCampaignEventDigest(event)
    event
}

publicationCampaignAppendEvent <- function(
    state,
    event_type,
    slot_id = NULL,
    details = list()
) {
    previous <- if (length(state$events)) {
        state$events[[length(state$events)]]$event_digest
    } else {
        NULL
    }
    state$events[[length(state$events) + 1L]] <- publicationCampaignEvent(
        event_type,
        length(state$events) + 1L,
        previous,
        slot_id,
        details
    )
    state
}

publicationValidateCampaignEvents <- function(events) {
    if (!is.list(events) || !length(events)) {
        publicationCampaignAbort("Campaign event chain is empty")
    }
    previous <- NULL
    for (index in seq_along(events)) {
        event <- events[[index]]
        publicationRequireNames(event, c(
            "event_index", "event_type", "occurred_at", "previous_event_digest",
            "slot_id", "details", "event_digest"
        ), "Campaign event")
        valid <- identical(event$event_index, as.integer(index)) &&
            publicationScalarString(event$event_type) &&
            publicationScalarString(event$occurred_at) &&
            identical(event$previous_event_digest, previous) &&
            publicationScalarString(event$event_digest) &&
            grepl("^[0-9a-f]{64}$", event$event_digest) &&
            identical(event$event_digest, publicationCampaignEventDigest(event))
        if (!valid) publicationCampaignAbort("Campaign event chain differs")
        previous <- event$event_digest
    }
    invisible(TRUE)
}

publicationCampaignRefreshDigest <- function(state) {
    state$state_digest <- NULL
    state$state_digest <- publicationCampaignStateDigest(state)
    state
}

publicationValidateCampaignPreflight <- function(preflight, schedule) {
    publicationRequireNames(preflight, c(
        "preflight_id", "host_id", "sha256", "certified"
    ), "Campaign preflight")
    if (!publicationScalarString(preflight$preflight_id) ||
        !identical(preflight$host_id, schedule$host_id) ||
        !publicationScalarString(preflight$sha256) ||
        !grepl("^[0-9a-f]{64}$", preflight$sha256) ||
        !isTRUE(preflight$certified)) {
        publicationCampaignAbort("Campaign preflight is not certified or bound")
    }
    invisible(preflight)
}

publicationValidateCampaignEvidence <- function(evidence) {
    publicationRequireNames(evidence, c(
        "result_sha256", "cleanup_sha256", "worker_started",
        "candidate_work_executed", "usage"
    ), "Campaign attempt evidence")
    digests <- c(evidence$result_sha256, evidence$cleanup_sha256)
    usage <- evidence$usage
    if (any(!vapply(digests, publicationScalarString, logical(1))) ||
        any(!grepl("^[0-9a-f]{64}$", digests)) ||
        !is.logical(evidence$worker_started) ||
        length(evidence$worker_started) != 1L ||
        !is.logical(evidence$candidate_work_executed) ||
        length(evidence$candidate_work_executed) != 1L ||
        !is.list(usage) ||
        !setequal(names(usage), publicationCampaignAttemptUsageNames()) ||
        any(!vapply(usage, \(value) {
            publicationScalarNumber(value) && value >= 0
        }, logical(1)))) {
        publicationCampaignAbort("Campaign attempt evidence is malformed")
    }
    invisible(evidence)
}

publicationCampaignBudgetBreaches <- function(usage, attempt_usage, budget) {
    limits <- budget$hard_limits
    breaches <- character()
    if (attempt_usage$worker_wall_seconds >
        limits$maximum_measured_worker_seconds_per_attempt) {
        breaches <- c(breaches, "maximum_measured_worker_seconds_per_attempt")
    }
    if (attempt_usage$temporary_bytes >
        limits$maximum_temporary_bytes_per_run) {
        breaches <- c(breaches, "maximum_temporary_bytes_per_run")
    }
    comparisons <- list(
        maximum_process_launches = c(
            usage$process_launches,
            limits$maximum_process_launches
        ),
        maximum_total_worker_wall_seconds = c(
            usage$worker_wall_seconds,
            limits$maximum_total_worker_wall_seconds
        ),
        maximum_total_worker_cpu_seconds = c(
            usage$worker_cpu_seconds,
            limits$maximum_total_worker_cpu_seconds
        ),
        maximum_generated_payload_bytes_retained = c(
            usage$generated_payload_bytes_retained,
            limits$maximum_generated_payload_bytes_retained
        ),
        maximum_raw_evidence_bytes_retained = c(
            usage$raw_evidence_bytes_retained,
            limits$maximum_raw_evidence_bytes_retained
        ),
        maximum_sanitized_publication_bytes = c(
            usage$sanitized_publication_bytes,
            limits$maximum_sanitized_publication_bytes
        )
    )
    for (name in names(comparisons)) {
        values <- comparisons[[name]]
        if (values[[1L]] > values[[2L]]) breaches <- c(breaches, name)
    }
    unique(breaches)
}

publicationCampaignApplyUsage <- function(state, attempt_usage, budget) {
    usage <- state$budget_usage
    usage$worker_wall_seconds <- usage$worker_wall_seconds +
        attempt_usage$worker_wall_seconds
    usage$worker_cpu_seconds <- usage$worker_cpu_seconds +
        attempt_usage$worker_cpu_seconds
    usage$maximum_temporary_bytes_observed <- max(
        usage$maximum_temporary_bytes_observed,
        attempt_usage$temporary_bytes
    )
    usage$generated_payload_bytes_retained <-
        usage$generated_payload_bytes_retained +
            attempt_usage$generated_payload_bytes_retained
    usage$raw_evidence_bytes_retained <- usage$raw_evidence_bytes_retained +
        attempt_usage$raw_evidence_bytes_retained
    usage$sanitized_publication_bytes <- usage$sanitized_publication_bytes +
        attempt_usage$sanitized_publication_bytes
    breaches <- publicationCampaignBudgetBreaches(usage, attempt_usage, budget)
    state$budget_usage <- usage
    list(state = state, breaches = breaches)
}

publicationInitializeCampaign <- function(
    campaign_id,
    schedule,
    protocol_digest,
    candidate_digest,
    budget,
    retry_policy
) {
    publicationValidateSchedule(schedule)
    digests <- c(protocol_digest, candidate_digest)
    if (!publicationScalarString(campaign_id) ||
        any(!vapply(digests, publicationScalarString, logical(1))) ||
        any(!grepl("^[0-9a-f]{64}$", digests))) {
        publicationCampaignAbort("Campaign bindings are incomplete")
    }
    slot_states <- lapply(schedule$slots, \(slot) {
        list(
            slot_id = slot$slot_id,
            status = "pending",
            attempts = list(),
            terminal_outcome = NULL
        )
    })
    names(slot_states) <- vapply(schedule$slots, `[[`, character(1), "slot_id")
    state <- list(
        schema = "multischolar.omics_publication_campaign_state",
        schema_version = "1.0.0",
        campaign_id = campaign_id,
        schedule_id = schedule$schedule_id,
        schedule_digest = schedule$schedule_digest,
        protocol_digest = protocol_digest,
        candidate_digest = candidate_digest,
        retry_policy_id = retry_policy$retry_policy_id,
        budget_id = budget$campaign_budget_id,
        status = "initialized",
        active_slot_id = NULL,
        budget_usage = publicationCampaignEmptyUsage(),
        budget_breaches = list(),
        slot_states = slot_states,
        events = list(),
        state_digest = NULL
    )
    state <- publicationCampaignAppendEvent(state, "campaign_initialized")
    state <- publicationCampaignRefreshDigest(state)
    publicationValidateCampaignState(state, schedule, budget, retry_policy)
    state
}

publicationValidateCampaignAttempts <- function(state, schedule, retry_policy) {
    active <- character()
    attempt_ids <- character()
    for (slot in state$slot_states) {
        publicationRequireNames(slot, c(
            "slot_id", "status", "attempts", "terminal_outcome"
        ), "Campaign slot state")
        if (length(slot$attempts) >
            retry_policy$attempt_contract$maximum_attempts_per_schedule_slot) {
            publicationCampaignAbort("Campaign slot exceeds retry bound")
        }
        for (index in seq_along(slot$attempts)) {
            attempt <- slot$attempts[[index]]
            publicationRequireNames(attempt, c(
                "attempt_id", "attempt_index", "status", "started_at",
                "finished_at", "reported_outcome", "outcome",
                "preflight_binding", "evidence"
            ), "Campaign attempt")
            expected_id <- paste(slot$slot_id, paste0("attempt-", index), sep = "::")
            if (!identical(attempt$attempt_id, expected_id) ||
                !identical(attempt$attempt_index, as.integer(index)) ||
                !attempt$status %in% c("running", "completed")) {
                publicationCampaignAbort("Campaign attempt identity differs")
            }
            publicationValidateCampaignPreflight(
                attempt$preflight_binding,
                schedule
            )
            if (identical(attempt$status, "running")) {
                if (!is.null(attempt$finished_at) || !is.null(attempt$evidence)) {
                    publicationCampaignAbort("Running attempt has terminal evidence")
                }
                active <- c(active, slot$slot_id)
            } else {
                publicationValidateCampaignEvidence(attempt$evidence)
                if (!publicationScalarString(attempt$finished_at) ||
                    !publicationScalarString(attempt$outcome)) {
                    publicationCampaignAbort("Completed attempt lacks outcome")
                }
                publicationCampaignValidateOutcomeEvidence(
                    attempt$reported_outcome,
                    attempt$evidence,
                    retry_policy
                )
            }
            attempt_ids <- c(attempt_ids, attempt$attempt_id)
        }
        last <- if (length(slot$attempts)) {
            slot$attempts[[length(slot$attempts)]]
        } else {
            NULL
        }
        slot_valid <- if (identical(slot$status, "pending")) {
            is.null(last) && is.null(slot$terminal_outcome)
        } else if (identical(slot$status, "running")) {
            identical(last$status, "running") && is.null(slot$terminal_outcome)
        } else if (identical(slot$status, "retry_pending")) {
            identical(last$status, "completed") &&
                is.null(slot$terminal_outcome)
        } else if (identical(slot$status, "passed")) {
            identical(last$outcome, "passed") &&
                slot$terminal_outcome %in% c(
                    "passed_first_attempt", "passed_single_harness_retry"
                )
        } else if (slot$status %in% c("failed", "aborted")) {
            identical(last$status, "completed") &&
                publicationScalarString(slot$terminal_outcome)
        } else {
            FALSE
        }
        if (!slot_valid) publicationCampaignAbort("Campaign slot transition differs")
    }
    if (anyDuplicated(attempt_ids) || length(active) > 1L ||
        (!length(active) && !is.null(state$active_slot_id)) ||
        (length(active) == 1L && !identical(active, state$active_slot_id))) {
        publicationCampaignAbort("Campaign active attempt ownership differs")
    }
    invisible(TRUE)
}

publicationCampaignReplayBudget <- function(state, budget) {
    usage <- publicationCampaignEmptyUsage()
    breaches <- list()
    for (slot in state$slot_states) {
        for (attempt in slot$attempts) {
            usage$process_launches <- usage$process_launches + 1L
            if (!identical(attempt$status, "completed")) next
            replay_state <- list(budget_usage = usage)
            applied <- publicationCampaignApplyUsage(
                replay_state,
                attempt$evidence$usage,
                budget
            )
            usage <- applied$state$budget_usage
            if (length(applied$breaches)) {
                breaches[[length(breaches) + 1L]] <- list(
                    slot_id = slot$slot_id,
                    attempt_id = attempt$attempt_id,
                    limits = as.list(applied$breaches),
                    evidence_sha256 = attempt$evidence$result_sha256
                )
            }
        }
    }
    list(usage = usage, breaches = breaches)
}

publicationValidateCampaignUsage <- function(state, budget) {
    if (!is.list(state$budget_usage) ||
        !setequal(names(state$budget_usage), publicationCampaignUsageNames()) ||
        any(!vapply(state$budget_usage, \(value) {
            publicationScalarNumber(value) && value >= 0
        }, logical(1)))) {
        publicationCampaignAbort("Campaign budget usage is malformed")
    }
    replay <- publicationCampaignReplayBudget(state, budget)
    usage_matches <- identical(
        names(state$budget_usage),
        names(replay$usage)
    ) && all(
        as.numeric(unlist(state$budget_usage, use.names = FALSE)) ==
            as.numeric(unlist(replay$usage, use.names = FALSE))
    )
    if (!usage_matches || !identical(
        publicationObjectDigest(state$budget_breaches),
        publicationObjectDigest(replay$breaches)
    ) ||
        state$budget_usage$process_launches >
            budget$hard_limits$maximum_process_launches) {
        publicationCampaignAbort("Campaign budget replay differs")
    }
    if (length(state$budget_breaches) &&
        !state$status %in% c("aborted", "safely_aborted")) {
        publicationCampaignAbort("Campaign continued after a budget breach")
    }
    invisible(TRUE)
}

publicationValidateCampaignState <- function(state, schedule, budget, retry_policy) {
    expected <- c(
        "schema", "schema_version", "campaign_id", "schedule_id",
        "schedule_digest", "protocol_digest", "candidate_digest",
        "retry_policy_id", "budget_id", "status", "active_slot_id",
        "budget_usage", "budget_breaches", "slot_states", "events",
        "state_digest"
    )
    publicationRequireNames(state, expected, "Publication campaign state")
    digests <- c(state$protocol_digest, state$candidate_digest)
    valid <- identical(
        state$schema,
        "multischolar.omics_publication_campaign_state"
    ) && identical(state$schema_version, "1.0.0") &&
        identical(state$schedule_id, schedule$schedule_id) &&
        identical(state$schedule_digest, schedule$schedule_digest) &&
        identical(state$retry_policy_id, retry_policy$retry_policy_id) &&
        identical(state$budget_id, budget$campaign_budget_id) &&
        state$status %in% c(
            "initialized", "running", "completed", "failed", "invalid",
            "aborted", "safely_aborted"
        ) && all(grepl("^[0-9a-f]{64}$", digests))
    if (!valid) publicationCampaignAbort("Campaign state binding differs")
    schedule_slots <- vapply(schedule$slots, `[[`, character(1), "slot_id")
    if (!identical(names(state$slot_states), schedule_slots)) {
        publicationCampaignAbort("Campaign slot state differs from schedule")
    }
    publicationValidateCampaignAttempts(state, schedule, retry_policy)
    publicationValidateCampaignUsage(state, budget)
    publicationValidateCampaignEvents(state$events)
    statuses <- vapply(state$slot_states, `[[`, character(1), "status")
    state_consistent <- if (identical(state$status, "completed")) {
        all(statuses == "passed")
    } else if (identical(state$status, "failed")) {
        any(statuses == "failed")
    } else if (identical(state$status, "invalid")) {
        any(statuses == "failed")
    } else if (identical(state$status, "aborted")) {
        any(statuses == "aborted")
    } else {
        TRUE
    }
    if (!state_consistent) {
        publicationCampaignAbort("Campaign aggregate status differs")
    }
    if (!identical(publicationCampaignStateDigest(state), state$state_digest)) {
        publicationCampaignAbort("Campaign state digest mismatch")
    }
    invisible(state)
}

publicationCampaignNextSlotId <- function(state) {
    statuses <- vapply(state$slot_states, `[[`, character(1), "status")
    available <- which(statuses %in% c("pending", "retry_pending"))
    if (length(available)) names(state$slot_states)[available[[1L]]] else NULL
}

publicationCampaignBeginAttempt <- function(
    state,
    slot_id,
    preflight,
    schedule,
    budget,
    retry_policy
) {
    publicationValidateCampaignState(state, schedule, budget, retry_policy)
    publicationValidateCampaignPreflight(preflight, schedule)
    expected_slot <- publicationCampaignNextSlotId(state)
    if (!state$status %in% c("initialized", "running") ||
        !is.null(state$active_slot_id) || !identical(slot_id, expected_slot)) {
        publicationCampaignAbort("Campaign cannot start requested slot")
    }
    slot <- state$slot_states[[slot_id]]
    if (!slot$status %in% c("pending", "retry_pending") ||
        length(slot$attempts) >=
            retry_policy$attempt_contract$maximum_attempts_per_schedule_slot ||
        state$budget_usage$process_launches >=
            budget$hard_limits$maximum_process_launches) {
        publicationCampaignAbort("Campaign slot or launch budget is exhausted")
    }
    attempt_index <- length(slot$attempts) + 1L
    attempt <- list(
        attempt_id = paste(slot_id, paste0("attempt-", attempt_index), sep = "::"),
        attempt_index = as.integer(attempt_index),
        status = "running",
        started_at = publicationUtcNow(),
        finished_at = NULL,
        reported_outcome = NULL,
        outcome = NULL,
        preflight_binding = preflight,
        evidence = NULL
    )
    slot$attempts[[attempt_index]] <- attempt
    slot$status <- "running"
    state$slot_states[[slot_id]] <- slot
    state$active_slot_id <- slot_id
    state$status <- "running"
    state$budget_usage$process_launches <-
        state$budget_usage$process_launches + 1L
    state <- publicationCampaignAppendEvent(
        state,
        "attempt_started",
        slot_id,
        list(
            attempt_id = attempt$attempt_id,
            preflight_sha256 = preflight$sha256
        )
    )
    publicationCampaignRefreshDigest(state)
}

publicationCampaignOutcomeTransition <- function(slot, outcome, retry_policy) {
    if (identical(outcome, "passed")) {
        retried <- length(slot$attempts) == 2L
        return(list(
            slot_status = "passed",
            terminal_outcome = if (retried) {
                "passed_single_harness_retry"
            } else {
                "passed_first_attempt"
            },
            campaign_status = "running"
        ))
    }
    retryable <- vapply(
        retry_policy$retryable_failures,
        `[[`,
        character(1),
        "failure_class"
    )
    if (outcome %in% retryable) {
        retry_available <- length(slot$attempts) <
            retry_policy$attempt_contract$maximum_attempts_per_schedule_slot
        return(list(
            slot_status = if (retry_available) "retry_pending" else "failed",
            terminal_outcome = if (retry_available) NULL else
                "failed_harness_retry_exhausted",
            campaign_status = if (retry_available) "running" else "invalid"
        ))
    }
    if (outcome %in% unlist(
        retry_policy$non_retryable_candidate_failures,
        use.names = FALSE
    )) {
        return(list(
            slot_status = "failed",
            terminal_outcome = outcome,
            campaign_status = "failed"
        ))
    }
    if (outcome %in% unlist(
        retry_policy$non_retryable_campaign_failures,
        use.names = FALSE
    )) {
        return(list(
            slot_status = "aborted",
            terminal_outcome = outcome,
            campaign_status = "aborted"
        ))
    }
    publicationCampaignAbort(paste("Unknown campaign outcome:", outcome))
}

publicationCampaignValidateOutcomeEvidence <- function(
    outcome,
    evidence,
    retry_policy
) {
    retryable <- vapply(
        retry_policy$retryable_failures,
        `[[`,
        character(1),
        "failure_class"
    )
    candidate_failures <- unlist(
        retry_policy$non_retryable_candidate_failures,
        use.names = FALSE
    )
    valid <- if (identical(outcome, "passed")) {
        isTRUE(evidence$worker_started) &&
            isTRUE(evidence$candidate_work_executed)
    } else if (outcome %in% retryable) {
        !isTRUE(evidence$worker_started) &&
            !isTRUE(evidence$candidate_work_executed)
    } else if (outcome %in% candidate_failures) {
        isTRUE(evidence$worker_started) &&
            isTRUE(evidence$candidate_work_executed)
    } else {
        TRUE
    }
    if (!valid) {
        publicationCampaignAbort("Campaign outcome evidence ownership differs")
    }
    invisible(TRUE)
}

publicationCampaignCompleteAttempt <- function(
    state,
    reported_outcome,
    evidence,
    schedule,
    budget,
    retry_policy
) {
    publicationValidateCampaignState(state, schedule, budget, retry_policy)
    publicationValidateCampaignEvidence(evidence)
    publicationCampaignValidateOutcomeEvidence(
        reported_outcome,
        evidence,
        retry_policy
    )
    slot_id <- state$active_slot_id
    if (is.null(slot_id) || !publicationScalarString(reported_outcome)) {
        publicationCampaignAbort("Campaign has no completable active attempt")
    }
    applied <- publicationCampaignApplyUsage(state, evidence$usage, budget)
    state <- applied$state
    effective_outcome <- if (length(applied$breaches)) {
        "budget_exhausted"
    } else {
        reported_outcome
    }
    if (length(applied$breaches)) {
        state$budget_breaches[[length(state$budget_breaches) + 1L]] <- list(
            slot_id = slot_id,
            attempt_id = state$slot_states[[slot_id]]$attempts[[
                length(state$slot_states[[slot_id]]$attempts)
            ]]$attempt_id,
            limits = as.list(applied$breaches),
            evidence_sha256 = evidence$result_sha256
        )
    }
    slot <- state$slot_states[[slot_id]]
    attempt_index <- length(slot$attempts)
    attempt <- slot$attempts[[attempt_index]]
    attempt$status <- "completed"
    attempt$finished_at <- publicationUtcNow()
    attempt$reported_outcome <- reported_outcome
    attempt$outcome <- effective_outcome
    attempt$evidence <- evidence
    slot$attempts[[attempt_index]] <- attempt
    state["active_slot_id"] <- list(NULL)
    transition <- publicationCampaignOutcomeTransition(
        slot,
        effective_outcome,
        retry_policy
    )
    slot$status <- transition$slot_status
    slot["terminal_outcome"] <- list(transition$terminal_outcome)
    state$slot_states[[slot_id]] <- slot
    state$status <- transition$campaign_status
    state <- publicationCampaignAppendEvent(
        state,
        "attempt_completed",
        slot_id,
        list(
            attempt_id = attempt$attempt_id,
            reported_outcome = reported_outcome,
            outcome = effective_outcome,
            result_sha256 = evidence$result_sha256
        )
    )
    if (identical(state$status, "running") && all(vapply(
        state$slot_states,
        \(item) identical(item$status, "passed"),
        logical(1)
    ))) {
        state$status <- "completed"
        state <- publicationCampaignAppendEvent(state, "campaign_completed")
    }
    publicationCampaignRefreshDigest(state)
}

publicationCampaignSafeAbort <- function(
    state,
    reason,
    schedule,
    budget,
    retry_policy
) {
    publicationValidateCampaignState(state, schedule, budget, retry_policy)
    allowed <- unlist(
        retry_policy$non_retryable_campaign_failures,
        use.names = FALSE
    )
    if (!state$status %in% c("initialized", "running") ||
        !is.null(state$active_slot_id) || !reason %in% allowed) {
        publicationCampaignAbort("Campaign cannot safely abort")
    }
    state$status <- "safely_aborted"
    state <- publicationCampaignAppendEvent(
        state,
        "campaign_safely_aborted",
        details = list(reason = reason)
    )
    publicationCampaignRefreshDigest(state)
}

publicationCampaignCheckpoint <- function(
    state,
    path,
    schedule,
    budget,
    retry_policy
) {
    publicationValidateCampaignState(state, schedule, budget, retry_policy)
    if (file.exists(path)) {
        publicationCampaignAbort("Campaign checkpoint already exists")
    }
    publicationWriteJson(state, path)
}

publicationCampaignResume <- function(
    path,
    schedule,
    budget,
    retry_policy,
    protocol_digest,
    candidate_digest
) {
    state <- jsonlite::read_json(path, simplifyVector = FALSE)
    publicationValidateCampaignState(state, schedule, budget, retry_policy)
    if (!identical(state$protocol_digest, protocol_digest) ||
        !identical(state$candidate_digest, candidate_digest) ||
        !is.null(state$active_slot_id) ||
        state$status %in% c(
            "completed", "failed", "invalid", "aborted", "safely_aborted"
        )) {
        publicationCampaignAbort("Campaign checkpoint cannot be resumed")
    }
    state
}
