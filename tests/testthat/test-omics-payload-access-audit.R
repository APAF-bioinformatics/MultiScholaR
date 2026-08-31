.payloadAuditRepoPath <- function(...) {
    file.path(
        normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
        ...
    )
}

.payloadAuditEnvironment <- new.env(parent = globalenv())
sys.source(
    .payloadAuditRepoPath(
        "tools",
        "refactor",
        "audit_omics_payload_access.R"
    ),
    envir = .payloadAuditEnvironment
)

.payloadAuditCopy <- function(value) {
    unserialize(serialize(value, NULL, version = 3L))
}

test_that("current payload access keys have one owner and pending obligations", {
    authority_path <- .payloadAuditRepoPath(
        "tools",
        "refactor",
        "omics-payload-access-owners-v1.json"
    )
    audit <- .payloadAuditEnvironment$payloadAccessAudit(
        .payloadAuditRepoPath(),
        "R",
        authority_path
    )

    expect_identical(audit$status, "ownership_assigned_pending_adoption")
    expect_identical(audit$production_access_count, 75L)
    expect_identical(audit$unresolved_computed_count, 0L)
    expect_length(audit$parse_errors, 0L)
    expect_length(audit$reconciliation$unowned_keys, 0L)
    expect_length(audit$reconciliation$stale_owner_keys, 0L)
    expect_length(audit$reconciliation$duplicate_owner_keys, 0L)
    expect_length(audit$reconciliation$mismatched_owner_keys, 0L)
    expect_length(audit$reconciliation$unresolved_computed_keys, 0L)
    expect_false(audit$candidate_freeze_allowed)

    obligations <- jsonlite::read_json(
        .payloadAuditRepoPath(
            "tests",
            "testdata",
            "omics-performance",
            "adapter-obligations-v1.json"
        ),
        simplifyVector = FALSE
    )
    expect_true(obligations$all_production_keys_owned)
    expect_false(obligations$all_adapter_obligations_closed)
    expect_false(obligations$candidate_freeze_allowed)
    expect_identical(
        obligations$access_authority$sha256,
        digest::digest(
            file = authority_path,
            algo = "sha256",
            serialize = FALSE
        )
    )
})

test_that("AST scanner classifies direct computed and S4 access exactly", {
    root <- withr::local_tempdir()
    dir.create(file.path(root, "R"))
    path <- file.path(root, "R", "fixture.R")
    writeLines(c(
        "read_value <- workflow_data$data_tbl",
        "workflow_data$data_cln <- read_value",
        "exact_value <- workflow_data[[\"data_tbl\"]]",
        "computed_value <- workflow_data[[field]]",
        "ordinary <- .data[[column]]",
        "slot_value <- object@data_cln",
        "# workflow_data$data_tbl is a comment"
    ), path)

    result <- .payloadAuditEnvironment$payloadAccessScanFile(
        path,
        root,
        c("data_tbl", "data_cln")
    )
    accesses <- result$accesses
    expect_length(result$errors, 0L)
    expect_length(accesses, 5L)
    expect_setequal(
        vapply(accesses, `[[`, character(1), "operator"),
        c("$", "$", "[[", "[[", "@")
    )
    expect_identical(
        sum(vapply(accesses, `[[`, character(1), "access_kind") == "write"),
        1L
    )
    expect_identical(
        sum(vapply(accesses, `[[`, logical(1), "computed")),
        1L
    )
    expressions <- vapply(accesses, `[[`, character(1), "expression")
    expect_false(any(grepl("ordinary", expressions, fixed = TRUE)))
    expect_false(any(grepl("comment", expressions, fixed = TRUE)))
})

test_that("access keys are stable across line movement", {
    root <- withr::local_tempdir()
    dir.create(file.path(root, "R"))
    path <- file.path(root, "R", "fixture.R")
    writeLines("value <- workflow_data$data_tbl", path)
    first <- .payloadAuditEnvironment$payloadAccessScanFile(
        path,
        root,
        c("data_tbl", "data_cln")
    )$accesses[[1L]]
    writeLines(c("", "", "value <- workflow_data$data_tbl"), path)
    second <- .payloadAuditEnvironment$payloadAccessScanFile(
        path,
        root,
        c("data_tbl", "data_cln")
    )$accesses[[1L]]

    expect_identical(first$access_key, second$access_key)
    expect_false(identical(first$location$line, second$location$line))
})

test_that("parse failures and owner mutations fail closed", {
    root <- withr::local_tempdir()
    dir.create(file.path(root, "R"))
    bad <- file.path(root, "R", "bad.R")
    writeLines("value <- function(", bad)
    result <- .payloadAuditEnvironment$payloadAccessScanFile(
        bad,
        root,
        c("data_tbl", "data_cln")
    )
    expect_length(result$errors, 1L)

    authority <- jsonlite::read_json(
        .payloadAuditRepoPath(
            "tools",
            "refactor",
            "omics-payload-access-owners-v1.json"
        ),
        simplifyVector = FALSE
    )
    duplicate <- .payloadAuditCopy(authority)
    duplicate$owners[[2L]] <- duplicate$owners[[1L]]
    expect_error(
        .payloadAuditEnvironment$payloadAccessValidateAuthority(duplicate),
        "owner authority is invalid"
    )

    weakened <- .payloadAuditCopy(authority)
    computed <- which(vapply(
        weakened$owners,
        function(owner) isTRUE(owner$computed),
        logical(1)
    ))[[1L]]
    weakened$owners[[computed]]$computed_resolution <- "exact_field"
    expect_error(
        .payloadAuditEnvironment$payloadAccessValidateAuthority(weakened),
        "owner authority is invalid"
    )
})
