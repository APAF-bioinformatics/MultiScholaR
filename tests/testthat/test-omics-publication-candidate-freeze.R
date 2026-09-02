.candidate077RepoPath <- function(...) {
    file.path(
        normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
        ...
    )
}

for (.candidate077_source in c(
    "omics_publication_protocol.R",
    "omics_publication_comparators.R"
)) {
    sys.source(
        .candidate077RepoPath("tools", "profiling", .candidate077_source),
        envir = environment()
    )
}
rm(.candidate077_source)

test_that("current candidate freeze readiness fails closed", {
    readiness <- publicationBuildCandidateFreezeReadiness()

    expect_silent(publicationValidateCandidateFreezeReadiness(readiness))
    expect_identical(readiness$status, "blocked")
    expect_false(readiness$candidate_freeze_allowed)
    expect_setequal(unlist(readiness$blockers), c(
        "role_principals_incomplete",
        "runtime_implementation_review_missing",
        "campaign_role_handoff_missing",
        "project_source_authority_incomplete",
        "split_successor_candidate_access_denied",
        "campaign_budget_execution_unauthorized",
        "auxiliary_candidate_access_denied",
        "proteomics_candidate_access_denied",
        "metabolomics_candidate_access_denied",
        "lipidomics_candidate_access_denied",
        "policy_receipts_authority_missing",
        "handoff_authority_missing",
        "blind_labels_authority_missing"
    ))
    unavailable <- names(readiness$authority_bindings)[!vapply(
        readiness$authority_bindings,
        `[[`,
        logical(1),
        "available"
    )]
    expect_setequal(
        unavailable,
        c("policy_receipts", "handoff", "blind_labels")
    )
})

test_that("candidate readiness invalidates every authority mutation", {
    readiness <- publicationBuildCandidateFreezeReadiness()
    available <- names(readiness$authority_bindings)[vapply(
        readiness$authority_bindings,
        `[[`,
        logical(1),
        "available"
    )]

    for (name in available) {
        changed <- publicationGovernanceCopy(readiness)
        changed$authority_bindings[[name]]$sha256 <- strrep("0", 64L)
        changed$readiness_digest <- publicationObjectDigest(
            changed[names(changed) != "readiness_digest"]
        )
        expect_error(
            publicationValidateCandidateFreezeReadiness(changed),
            class = "multischolar_publication_comparator_error",
            info = name
        )
    }

    blocker_drift <- publicationGovernanceCopy(readiness)
    blocker_drift$blockers <- list()
    expect_error(
        publicationValidateCandidateFreezeReadiness(blocker_drift),
        class = "multischolar_publication_comparator_error"
    )
})

test_that("candidate successor cannot bypass blocked readiness", {
    authority <- publicationReadJson(.candidate077RepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "comparators-v1.json"
    ))
    revision <- publicationComparatorRevisions()[[
        "pre_repair_performance"
    ]]
    digest <- strrep("a", 64L)
    bindings <- stats::setNames(
        as.list(rep(digest, length(publicationCandidateFreezeBindingNames()))),
        publicationCandidateFreezeBindingNames()
    )

    expect_error(
        publicationComparatorFreezeCandidate(
            authority,
            revision,
            publicationComparatorSourceIdentity(revision),
            digest,
            bindings,
            publicationBuildCandidateFreezeReadiness()
        ),
        class = "multischolar_publication_comparator_error"
    )
})
