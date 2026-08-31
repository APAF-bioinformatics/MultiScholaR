test_that("retry policy is append-only bounded and candidate failures are terminal", {
    retry <- publicationGovernanceRead("retry-policy-v1.json")

    expect_silent(publicationGovernanceValidateRetry(retry))
    expect_true(all(c("oom", "scientific_parity_mismatch", "resource_leak") %in%
        unlist(retry$non_retryable_candidate_failures)))

    unbounded <- publicationGovernanceCopy(retry)
    unbounded$attempt_contract$maximum_retries_per_schedule_slot <- 2L
    expect_error(
        publicationGovernanceValidateRetry(unbounded),
        class = "multischolar_publication_governance_error"
    )

    selective <- publicationGovernanceCopy(retry)
    selective$retryable_failures[[1L]]$candidate_work_executed <- TRUE
    expect_error(
        publicationGovernanceValidateRetry(selective),
        class = "multischolar_publication_governance_error"
    )
})

test_that("campaign matrix has a finite reconciled process ceiling", {
    matrix <- publicationGovernanceRead("campaign-matrix-v1.json")

    expect_silent(publicationGovernanceValidateCampaignMatrix(matrix))
    expect_identical(matrix$totals$primary_maximum_pair_slots, 4975L)
    expect_identical(matrix$totals$computed_maximum_process_launches, 21868L)
    expect_identical(matrix$totals$governed_maximum_process_launches, 22000L)
    expect_false(matrix$current_readiness$exact_pair_count_selected)
    expect_false(matrix$current_readiness$campaign_execution_authorized)

    drift <- publicationGovernanceCopy(matrix)
    drift$totals$computed_maximum_process_launches <- 1L
    expect_error(
        publicationGovernanceValidateCampaignMatrix(drift),
        class = "multischolar_publication_governance_error"
    )
})

test_that("campaign budget binds matrix retries host safety and safe abort", {
    matrix <- publicationGovernanceRead("campaign-matrix-v1.json")
    retry <- publicationGovernanceRead("retry-policy-v1.json")
    budget <- publicationGovernanceRead("campaign-budget-v1.json")

    expect_silent(publicationGovernanceValidateBudget(budget, matrix, retry))
    expect_identical(
        budget$resume_contract$safe_abort_outcome,
        "campaign_safely_aborted_incomplete_non_promotional"
    )
    expect_false(budget$current_status$execution_authorized)

    too_small <- publicationGovernanceCopy(budget)
    too_small$hard_limits$maximum_process_launches <- 100L
    expect_error(
        publicationGovernanceValidateBudget(too_small, matrix, retry),
        class = "multischolar_publication_governance_error"
    )

    executable <- publicationGovernanceCopy(budget)
    executable$current_status$execution_authorized <- TRUE
    expect_error(
        publicationGovernanceValidateBudget(executable, matrix, retry),
        class = "multischolar_publication_governance_error"
    )
})
