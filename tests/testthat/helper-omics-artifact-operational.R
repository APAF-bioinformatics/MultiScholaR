operationalArtifactRepoPath <- function(...) {
    file.path(
        normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
        ...
    )
}

operationalArtifactReadCloseout <- function() {
    jsonlite::read_json(
        operationalArtifactRepoPath(
            "tests", "testdata", "omics-parity",
            "all-omics-operational-closeout-v1.json"
        ),
        simplifyVector = FALSE
    )
}

operationalArtifactReadPolicyReconciliation <- function() {
    jsonlite::read_json(
        operationalArtifactRepoPath(
            "tests", "testdata", "omics-performance",
            "current-policy-reconciliation-v1.json"
        ),
        simplifyVector = FALSE
    )
}

operationalArtifactSkipDependencies <- function() {
    for (package in c("arrow", "DBI", "duckdb", "filelock")) {
        testthat::skip_if_not_installed(package)
    }
}

operationalArtifactFixture <- function(
    project_root,
    omic_type,
    omic_label,
    workflow_slug,
    generation_id
) {
    identity <- list(
        omic_type = omic_type,
        omic_label = omic_label,
        workflow_slug = workflow_slug
    )
    paths <- buildArtifactPaths(project_root, identity)
    project_id <- "all-omics-operational"
    store <- newArtifactStore(paths, project_id)
    logical_key <- list(
        project_id = project_id,
        omic_type = omic_type,
        workflow_slug = workflow_slug,
        stage_id = "operational_closeout",
        state_role = "recovery_probe",
        generation_id = generation_id
    )
    value <- data.frame(
        feature_id = c("F1", "F2", "F3"),
        sample_id = c("S1", "S2", "S3"),
        value = c(1, 2, 3),
        stringsAsFactors = FALSE
    )
    list(
        paths = paths,
        store = store,
        logical_key = logical_key,
        value = value,
        encoded = encodeArtifactTable(
            value,
            stable_key = "feature_id",
            owner = "all-omics-operational"
        )
    )
}

operationalArtifactFailAt <- function(target_stage) {
    force(target_stage)
    \(stage, context) {
        if (identical(stage, target_stage)) {
            rlang::abort(
                paste("injected operational failure at", stage),
                class = "multischolar_test_operational_failure"
            )
        }
        invisible(context)
    }
}

operationalArtifactLoadAssignments <- function(path, names) {
    expressions <- parse(testthat::test_path(path))
    target <- parent.frame()
    for (expression in expressions) {
        assignment <- is.call(expression) &&
            identical(expression[[1L]], as.name("<-"))
        if (isTRUE(assignment) && as.character(expression[[2L]]) %in% names) {
            eval(expression, envir = target)
        }
    }
    invisible(TRUE)
}
