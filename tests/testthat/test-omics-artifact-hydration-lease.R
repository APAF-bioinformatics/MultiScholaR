.hydrationLeaseRefs <- function(seed = "a") {
    list(
        data_tbl = list(semantic_digest = strrep(seed, 64L)),
        data_cln = list(semantic_digest = strrep(seed, 64L))
    )
}

.hydrationLeaseManager <- function(scope = NULL) {
    ArtifactHydrationLeaseManager$new(
        project_id = "lease-project",
        omic_type = "proteomics",
        descriptor_id = "proteomics.diann.peptide.dia.v1",
        scope = scope
    )
}

test_that("one lease validates publishes and releases explicit ownership", {
    owner <- new.env(parent = emptyenv())
    owner$value <- list(locator = "complete-stage-owner")
    releases <- 0L
    manager <- .hydrationLeaseManager()

    token <- manager$acquire(
        "generation-1",
        "import",
        .hydrationLeaseRefs(),
        function() {
            releases <<- releases + 1L
            owner$value <- NULL
        }
    )
    expect_true(manager$isActive())
    expect_identical(manager$getInfo()$resource_count, 1L)
    expect_silent(artifactResourceDataOnly(token))
    expect_error(
        manager$acquire(
            "generation-1",
            "design",
            .hydrationLeaseRefs(),
            function() NULL
        ),
        class = "multischolar_overlapping_artifact_hydration"
    )

    validated <- manager$markValidated(token)
    expect_identical(validated$state, "validated")
    expect_true(validated$validated)
    published <- manager$markPublished(validated)
    expect_identical(published$state, "published")
    expect_true(published$published)
    expect_true(manager$release(published))
    expect_false(manager$isActive())
    expect_identical(manager$getInfo()$resource_count, 0L)
    expect_identical(releases, 1L)
    expect_null(owner$value)
    expect_true(manager$close())
    expect_false(manager$close())
})

test_that("lease publication rejects unvalidated and stale tokens", {
    manager <- .hydrationLeaseManager()
    withr::defer(try(manager$close(), silent = TRUE))
    token <- manager$acquire(
        "generation-1",
        "import",
        .hydrationLeaseRefs(),
        function() NULL
    )

    expect_error(
        manager$markPublished(token),
        class = "multischolar_stale_artifact_hydration_lease"
    )
    validated <- manager$markValidated(token)
    expect_error(
        manager$markValidated(token),
        class = "multischolar_stale_artifact_hydration_lease"
    )

    cross_process <- validated
    cross_process$creator_pid <- cross_process$creator_pid + 1L
    expect_error(
        manager$markPublished(cross_process),
        class = "multischolar_invalid_artifact_hydration_lease"
    )

    wrong_generation <- validated
    wrong_generation$generation_id <- "generation-2"
    expect_error(
        manager$markPublished(wrong_generation),
        class = "multischolar_stale_artifact_hydration_lease"
    )

    published <- manager$markPublished(validated)
    expect_true(manager$release(published))
})

test_that("lease identity rejects every mixed ownership dimension", {
    manager <- .hydrationLeaseManager()
    withr::defer(try(manager$close(), silent = TRUE))
    token <- manager$acquire(
        "generation-1",
        "import",
        .hydrationLeaseRefs(),
        function() NULL
    )
    validated <- manager$markValidated(token)
    mutations <- list(
        project = function(value) {
            value$project_id <- "other-project"
            value
        },
        omic = function(value) {
            value$omic_type <- "lipidomics"
            value
        },
        descriptor = function(value) {
            value$descriptor_id <- "lipidomics.lipidsearch.lipid.standard.v1"
            value
        },
        generation = function(value) {
            value$generation_id <- "generation-2"
            value
        },
        stage = function(value) {
            value$stage_id <- "design"
            value
        },
        roles = function(value) {
            value$roles <- rev(value$roles)
            value
        },
        refs = function(value) {
            value$ref_digests[[1L]] <- strrep("f", 64L)
            value
        }
    )
    for (name in names(mutations)) {
        expect_error(
            manager$markPublished(mutations[[name]](validated)),
            class = "multischolar_stale_artifact_hydration_lease",
            info = name
        )
    }
    published <- manager$markPublished(validated)
    expect_true(manager$release(published))
    expect_error(
        manager$release(published),
        class = "multischolar_stale_artifact_hydration_lease"
    )
})

test_that("lease release clears owner cache query and registry references", {
    resources <- new.env(parent = emptyenv())
    resources$owner <- list(value = "hydrated")
    resources$cache <- list(value = "cached")
    resources$query <- list(value = "query")
    resources$registry <- list(value = "registry")
    manager <- .hydrationLeaseManager()
    token <- manager$acquire(
        "generation-1",
        "import",
        .hydrationLeaseRefs(),
        function() {
            resources$owner <- NULL
            resources$cache <- NULL
            resources$query <- NULL
            resources$registry <- NULL
        }
    )
    token <- manager$markValidated(token)
    expect_true(manager$release(token))
    expect_true(all(vapply(
        c("owner", "cache", "query", "registry"),
        function(name) is.null(resources[[name]]),
        logical(1)
    )))
    expect_identical(manager$getInfo()$resource_count, 0L)
    expect_true(manager$close())
})

test_that("resource scope cancellation releases an active lease", {
    scope <- ArtifactResourceScope$new("lease-cancellation")
    owner <- new.env(parent = emptyenv())
    owner$value <- "active"
    manager <- .hydrationLeaseManager(scope)
    manager$acquire(
        "generation-1",
        "import",
        .hydrationLeaseRefs(),
        function() owner$value <- NULL
    )

    expect_true(scope$close("session_cancelled"))
    expect_false(manager$isActive())
    expect_null(owner$value)
    expect_identical(manager$getInfo()$resource_count, 0L)
    expect_true(manager$close())
})

test_that("release callback failure clears lease state and reports cleanup", {
    manager <- .hydrationLeaseManager()
    token <- manager$acquire(
        "generation-1",
        "import",
        .hydrationLeaseRefs(),
        function() stop("injected release failure")
    )
    token <- manager$markValidated(token)

    expect_error(
        manager$release(token),
        class = "multischolar_artifact_resource_cleanup_failed"
    )
    expect_false(manager$isActive())
    expect_identical(manager$getInfo()$resource_count, 0L)
    expect_error(
        manager$close(),
        class = "multischolar_artifact_resource_cleanup_failed"
    )
    expect_false(manager$close())
})

test_that("lease identity rejects invalid refs and closed managers", {
    manager <- .hydrationLeaseManager()
    invalid <- .hydrationLeaseRefs()
    invalid$data_tbl$semantic_digest <- "not-a-digest"
    expect_error(
        manager$acquire(
            "generation-1",
            "import",
            invalid,
            function() NULL
        ),
        class = "multischolar_invalid_artifact_hydration_lease"
    )
    expect_true(manager$close())
    expect_error(
        manager$acquire(
            "generation-1",
            "import",
            .hydrationLeaseRefs(),
            function() NULL
        ),
        class = "multischolar_closed_artifact_hydration_lease"
    )
})

test_that("repeated lease cycles keep structural resource counts bounded", {
    manager <- .hydrationLeaseManager()
    withr::defer(try(manager$close(), silent = TRUE))
    counts <- integer(100L)
    releases <- 0L
    for (index in seq_len(100L)) {
        token <- manager$acquire(
            paste0("generation-", index),
            "import",
            .hydrationLeaseRefs(),
            function() releases <<- releases + 1L
        )
        expect_identical(manager$getInfo()$resource_count, 1L)
        token <- manager$markValidated(token)
        expect_true(manager$release(token))
        counts[[index]] <- manager$getInfo()$resource_count
    }
    expect_true(all(counts == 0L))
    expect_identical(releases, 100L)
    expect_false(manager$isActive())
})

test_that("hydration lease release does not call garbage collection", {
    source <- readLines(
        testthat::test_path(
            "..",
            "..",
            "R",
            "utils_artifact_hydration_lease.R"
        ),
        warn = FALSE
    )
    expect_false(any(grepl("gc(", source, fixed = TRUE)))
    expect_false(any(grepl("artifactReleaseTransientMemory", source, fixed = TRUE)))
})
