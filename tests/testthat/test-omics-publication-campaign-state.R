publicationCampaignFixture <- function(maximum_process_launches = 22000L) {
    schedule <- publicationBuildSchedule(
        projects = publicationScheduleFixtureProjects(),
        estimand_ids = "estimand-1",
        comparison_ids = "backend-effect",
        pair_count = 30L,
        host_id = "primary-host",
        work_bindings = publicationScheduleFixtureWorkBindings()
    )
    retry <- publicationGovernanceRecord("retry_policy")
    budget <- publicationGovernanceRecord("campaign_budget")
    budget$hard_limits$maximum_process_launches <- maximum_process_launches
    state <- publicationInitializeCampaign(
        "campaign-1",
        schedule,
        paste(rep("a", 64L), collapse = ""),
        paste(rep("b", 64L), collapse = ""),
        budget,
        retry
    )
    list(schedule = schedule, retry = retry, budget = budget, state = state)
}

publicationCampaignPreflightFixture <- function(schedule) {
    list(
        preflight_id = "preflight-1",
        host_id = schedule$host_id,
        sha256 = paste(rep("c", 64L), collapse = ""),
        certified = TRUE
    )
}

publicationCampaignUsageFixture <- function(
    worker_wall_seconds = 1,
    worker_cpu_seconds = 0.5,
    temporary_bytes = 100,
    generated_payload_bytes_retained = 10,
    raw_evidence_bytes_retained = 20,
    sanitized_publication_bytes = 5
) {
    list(
        worker_wall_seconds = worker_wall_seconds,
        worker_cpu_seconds = worker_cpu_seconds,
        temporary_bytes = temporary_bytes,
        generated_payload_bytes_retained = generated_payload_bytes_retained,
        raw_evidence_bytes_retained = raw_evidence_bytes_retained,
        sanitized_publication_bytes = sanitized_publication_bytes
    )
}

publicationCampaignEvidenceFixture <- function(
    worker_started = TRUE,
    candidate_work_executed = TRUE,
    usage = publicationCampaignUsageFixture(),
    digest_character = "d"
) {
    digest <- paste(rep(digest_character, 64L), collapse = "")
    list(
        result_sha256 = digest,
        cleanup_sha256 = paste(rep("e", 64L), collapse = ""),
        worker_started = worker_started,
        candidate_work_executed = candidate_work_executed,
        usage = usage
    )
}

test_that("campaign records ordered append-only successful attempts and usage", {
    fixture <- publicationCampaignFixture()
    first_id <- fixture$schedule$slots[[1L]]$slot_id
    second_id <- fixture$schedule$slots[[2L]]$slot_id
    preflight <- publicationCampaignPreflightFixture(fixture$schedule)

    expect_error(
        publicationCampaignBeginAttempt(
            fixture$state,
            second_id,
            preflight,
            fixture$schedule,
            fixture$budget,
            fixture$retry
        ),
        class = "multischolar_publication_campaign_error"
    )

    state <- publicationCampaignBeginAttempt(
        fixture$state,
        first_id,
        preflight,
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    expect_identical(state$active_slot_id, first_id)
    expect_identical(state$budget_usage$process_launches, 1)

    state <- publicationCampaignCompleteAttempt(
        state,
        "passed",
        publicationCampaignEvidenceFixture(),
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    expect_identical(state$slot_states[[first_id]]$status, "passed")
    expect_identical(
        state$slot_states[[first_id]]$terminal_outcome,
        "passed_first_attempt"
    )
    expect_identical(state$budget_usage$worker_wall_seconds, 1)
    expect_length(state$slot_states[[first_id]]$attempts, 1L)
    expect_null(state$active_slot_id)
    expect_silent(publicationValidateCampaignState(
        state,
        fixture$schedule,
        fixture$budget,
        fixture$retry
    ))
})

test_that("campaign retries only pre-start harness failures once", {
    fixture <- publicationCampaignFixture()
    slot_id <- fixture$schedule$slots[[1L]]$slot_id
    preflight <- publicationCampaignPreflightFixture(fixture$schedule)
    harness_evidence <- publicationCampaignEvidenceFixture(
        worker_started = FALSE,
        candidate_work_executed = FALSE,
        usage = publicationCampaignUsageFixture(
            worker_wall_seconds = 0,
            worker_cpu_seconds = 0,
            temporary_bytes = 0,
            generated_payload_bytes_retained = 0,
            raw_evidence_bytes_retained = 1,
            sanitized_publication_bytes = 0
        )
    )

    state <- publicationCampaignBeginAttempt(
        fixture$state,
        slot_id,
        preflight,
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    state <- publicationCampaignCompleteAttempt(
        state,
        "worker_not_started",
        harness_evidence,
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    expect_identical(state$slot_states[[slot_id]]$status, "retry_pending")

    state <- publicationCampaignBeginAttempt(
        state,
        slot_id,
        preflight,
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    state <- publicationCampaignCompleteAttempt(
        state,
        "passed",
        publicationCampaignEvidenceFixture(digest_character = "f"),
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    expect_identical(
        state$slot_states[[slot_id]]$terminal_outcome,
        "passed_single_harness_retry"
    )
    expect_length(state$slot_states[[slot_id]]$attempts, 2L)

    exhausted <- publicationCampaignBeginAttempt(
        fixture$state,
        slot_id,
        preflight,
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    exhausted <- publicationCampaignCompleteAttempt(
        exhausted,
        "worker_not_started",
        harness_evidence,
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    exhausted <- publicationCampaignBeginAttempt(
        exhausted,
        slot_id,
        preflight,
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    exhausted <- publicationCampaignCompleteAttempt(
        exhausted,
        "worker_not_started",
        harness_evidence,
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    expect_identical(exhausted$status, "invalid")
    expect_identical(
        exhausted$slot_states[[slot_id]]$terminal_outcome,
        "failed_harness_retry_exhausted"
    )
})

test_that("candidate failures and observed budget breaches are terminal", {
    fixture <- publicationCampaignFixture(maximum_process_launches = 1L)
    slot_id <- fixture$schedule$slots[[1L]]$slot_id
    preflight <- publicationCampaignPreflightFixture(fixture$schedule)
    state <- publicationCampaignBeginAttempt(
        fixture$state,
        slot_id,
        preflight,
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    state <- publicationCampaignCompleteAttempt(
        state,
        "oom",
        publicationCampaignEvidenceFixture(),
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    expect_identical(state$status, "failed")
    expect_identical(state$slot_states[[slot_id]]$terminal_outcome, "oom")

    over <- publicationCampaignFixture()
    over$budget$hard_limits$maximum_measured_worker_seconds_per_attempt <- 0.5
    over$state <- publicationInitializeCampaign(
        "campaign-over-budget",
        over$schedule,
        paste(rep("a", 64L), collapse = ""),
        paste(rep("b", 64L), collapse = ""),
        over$budget,
        over$retry
    )
    over_state <- publicationCampaignBeginAttempt(
        over$state,
        over$schedule$slots[[1L]]$slot_id,
        publicationCampaignPreflightFixture(over$schedule),
        over$schedule,
        over$budget,
        over$retry
    )
    over_state <- publicationCampaignCompleteAttempt(
        over_state,
        "passed",
        publicationCampaignEvidenceFixture(worker_started = TRUE),
        over$schedule,
        over$budget,
        over$retry
    )
    expect_identical(over_state$status, "aborted")
    expect_identical(
        over_state$slot_states[[1L]]$terminal_outcome,
        "budget_exhausted"
    )
    expect_length(over_state$budget_breaches, 1L)
    expect_true(
        "maximum_measured_worker_seconds_per_attempt" %in%
            unlist(over_state$budget_breaches[[1L]]$limits)
    )
})

test_that("checkpoint resume safe abort and event chains preserve bindings", {
    fixture <- publicationCampaignFixture()
    checkpoint <- tempfile(fileext = ".json")
    publicationCampaignCheckpoint(
        fixture$state,
        checkpoint,
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    resumed <- publicationCampaignResume(
        checkpoint,
        fixture$schedule,
        fixture$budget,
        fixture$retry,
        fixture$state$protocol_digest,
        fixture$state$candidate_digest
    )
    expect_identical(resumed$state_digest, fixture$state$state_digest)
    expect_error(
        publicationCampaignCheckpoint(
            fixture$state,
            checkpoint,
            fixture$schedule,
            fixture$budget,
            fixture$retry
        ),
        class = "multischolar_publication_campaign_error"
    )

    aborted <- publicationCampaignSafeAbort(
        fixture$state,
        "unsafe_host_state",
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    expect_identical(aborted$status, "safely_aborted")
    expect_error(
        publicationCampaignSafeAbort(
            aborted,
            "unsafe_host_state",
            fixture$schedule,
            fixture$budget,
            fixture$retry
        ),
        class = "multischolar_publication_campaign_error"
    )

    tampered <- publicationGovernanceCopy(fixture$state)
    tampered$events[[1L]]$details$extra <- "changed"
    tampered <- publicationCampaignRefreshDigest(tampered)
    expect_error(
        publicationValidateCampaignState(
            tampered,
            fixture$schedule,
            fixture$budget,
            fixture$retry
        ),
        class = "multischolar_publication_campaign_error"
    )

    usage_drift <- publicationGovernanceCopy(fixture$state)
    usage_drift$budget_usage$worker_cpu_seconds <- 1
    usage_drift <- publicationCampaignRefreshDigest(usage_drift)
    expect_error(
        publicationValidateCampaignState(
            usage_drift,
            fixture$schedule,
            fixture$budget,
            fixture$retry
        ),
        class = "multischolar_publication_campaign_error"
    )
})

test_that("campaign cannot misclassify candidate work as harness retry", {
    fixture <- publicationCampaignFixture()
    slot_id <- fixture$schedule$slots[[1L]]$slot_id
    state <- publicationCampaignBeginAttempt(
        fixture$state,
        slot_id,
        publicationCampaignPreflightFixture(fixture$schedule),
        fixture$schedule,
        fixture$budget,
        fixture$retry
    )
    expect_error(
        publicationCampaignCompleteAttempt(
            state,
            "worker_not_started",
            publicationCampaignEvidenceFixture(
                worker_started = TRUE,
                candidate_work_executed = TRUE
            ),
            fixture$schedule,
            fixture$budget,
            fixture$retry
        ),
        class = "multischolar_publication_campaign_error"
    )
})

test_that("campaign reaches completion only after every frozen slot", {
    fixture <- publicationCampaignFixture()
    state <- fixture$state
    preflight <- publicationCampaignPreflightFixture(fixture$schedule)
    evidence <- publicationCampaignEvidenceFixture(
        usage = publicationCampaignUsageFixture(
            worker_wall_seconds = 0.1,
            worker_cpu_seconds = 0.05,
            temporary_bytes = 1,
            generated_payload_bytes_retained = 0,
            raw_evidence_bytes_retained = 1,
            sanitized_publication_bytes = 0
        )
    )
    for (slot in fixture$schedule$slots) {
        state <- publicationCampaignBeginAttempt(
            state,
            slot$slot_id,
            preflight,
            fixture$schedule,
            fixture$budget,
            fixture$retry
        )
        state <- publicationCampaignCompleteAttempt(
            state,
            "passed",
            evidence,
            fixture$schedule,
            fixture$budget,
            fixture$retry
        )
    }

    expect_identical(state$status, "completed")
    expect_true(all(vapply(state$slot_states, \(slot) {
        identical(slot$terminal_outcome, "passed_first_attempt")
    }, logical(1))))
    expect_equal(
        state$budget_usage$process_launches,
        length(fixture$schedule$slots)
    )
})
