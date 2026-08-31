.publicationExclusionRepoPath <- function(...) {
    file.path(
        normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
        ...
    )
}

for (.publication_exclusion_source in c(
    "omics_publication_protocol.R",
    "omics_publication_auxiliary_contracts.R",
    "omics_publication_auxiliary_exclusions.R"
)) {
    sys.source(
        .publicationExclusionRepoPath(
            "tools",
            "profiling",
            .publication_exclusion_source
        ),
        envir = environment()
    )
}
rm(.publication_exclusion_source)

.publicationExclusionRecord <- function() {
    publicationReadJson(.publicationExclusionRepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "auxiliary",
        "exclusions-v1.json"
    ))
}

test_that("auxiliary exclusions exactly reconcile both source authorities", {
    exclusions <- .publicationExclusionRecord()

    expect_silent(auxPublicationValidateExclusions(exclusions))
    expect_identical(
        publicationObjectDigest(exclusions),
        publicationObjectDigest(auxPublicationBuildExclusions())
    )
    expect_length(exclusions$cases, 13L)
    types <- vapply(exclusions$cases, `[[`, character(1), "subject_type")
    expect_identical(
        as.integer(table(types)[c(
            "advertised_unverified_capability",
            "detection_only_format",
            "ui_placeholder",
            "unknown_capability",
            "unknown_format",
            "unknown_surface"
        )]),
        c(2L, 6L, 2L, 1L, 1L, 1L)
    )
})

test_that("every auxiliary exclusion rejects before downstream invocation", {
    exclusions <- .publicationExclusionRecord()
    receipts <- lapply(exclusions$cases, auxPublicationEvaluateExclusion)

    expect_length(receipts, 13L)
    expect_true(all(vapply(receipts, function(receipt) {
        identical(
            receipt$condition_class,
            "multischolar_publication_auxiliary_surface_excluded"
        ) && identical(receipt$prepare_invocations, 0L) &&
            identical(receipt$runner_invocations, 0L) &&
            identical(receipt$artifact_invocations, 0L) &&
            identical(receipt$backend_invocations, 0L) &&
            !isTRUE(receipt$promotion_authority)
    }, logical(1))))
})

test_that("auxiliary exclusions reject identity and permission mutations", {
    exclusions <- .publicationExclusionRecord()

    renamed <- publicationGovernanceCopy(exclusions)
    renamed$cases[[1L]]$case_id <- "renamed.excluded.v1"
    expect_error(
        auxPublicationValidateExclusions(renamed),
        class = "multischolar_publication_auxiliary_contract_error"
    )

    invoked <- publicationGovernanceCopy(exclusions)
    invoked$cases[[1L]]$runner_invocation_allowed <- TRUE
    expect_error(
        auxPublicationValidateExclusions(invoked),
        class = "multischolar_publication_auxiliary_contract_error"
    )

    omitted <- publicationGovernanceCopy(exclusions)
    omitted$cases <- omitted$cases[-1L]
    expect_error(
        auxPublicationValidateExclusions(omitted),
        class = "multischolar_publication_auxiliary_contract_error"
    )

    duplicated <- publicationGovernanceCopy(exclusions)
    duplicated$cases[[2L]] <- duplicated$cases[[1L]]
    expect_error(
        auxPublicationValidateExclusions(duplicated),
        class = "multischolar_publication_auxiliary_contract_error"
    )

    stale <- publicationGovernanceCopy(exclusions)
    stale$coverage_binding$sha256 <- strrep("0", 64L)
    expect_error(
        auxPublicationValidateExclusions(stale),
        class = "multischolar_publication_auxiliary_contract_error"
    )
})
