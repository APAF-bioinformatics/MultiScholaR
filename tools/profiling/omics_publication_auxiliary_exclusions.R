auxPublicationExpectedExclusions <- function() {
    capability <- publicationReadJson("tests/testdata/omics-capabilities.json")
    coverage <- publicationReadJson(
        "tests/testdata/omics-performance/coverage-v1.json"
    )
    placeholders <- lapply(capability$non_workflow_surfaces, function(surface) {
        if (!identical(surface$surface_type, "ui_placeholder")) return(NULL)
        list(
            subject_id = surface$surface_id,
            subject_type = "ui_placeholder",
            support_status = "excluded_until_implemented"
        )
    })
    placeholders <- Filter(Negate(is.null), placeholders)
    capabilities <- lapply(coverage$excluded_capabilities, function(value) {
        list(
            subject_id = value$capability_id,
            subject_type = "advertised_unverified_capability",
            support_status = value$support_status
        )
    })
    formats <- lapply(coverage$excluded_formats, function(value) {
        list(
            subject_id = value$format_id,
            subject_type = "detection_only_format",
            support_status = value$support_status
        )
    })
    unknown <- list(
        list(
            subject_id = "unknown.surface",
            subject_type = "unknown_surface",
            support_status = "unknown"
        ),
        list(
            subject_id = "unknown.capability",
            subject_type = "unknown_capability",
            support_status = "unknown"
        ),
        list(
            subject_id = "unknown.format",
            subject_type = "unknown_format",
            support_status = "unknown"
        )
    )
    c(placeholders, capabilities, formats, unknown)
}

auxPublicationValidateExclusionCase <- function(case) {
    publicationRequireNames(case, c(
        "case_id", "subject_id", "subject_type", "support_status",
        "expected_condition_class", "prepare_invocation_allowed",
        "runner_invocation_allowed", "artifact_invocation_allowed",
        "backend_invocation_allowed", "positive_matrix_allowed",
        "promotion_authority"
    ), "Auxiliary exclusion case")
    expected <- Filter(function(value) {
        identical(value$subject_id, case$subject_id)
    }, auxPublicationExpectedExclusions())
    valid <- length(expected) == 1L &&
        identical(
            case$case_id,
            paste0(expected[[1L]]$subject_id, ".excluded.v1")
        ) &&
        identical(case$subject_type, expected[[1L]]$subject_type) &&
        identical(case$support_status, expected[[1L]]$support_status) &&
        identical(
            case$expected_condition_class,
            "multischolar_publication_auxiliary_surface_excluded"
        ) && !isTRUE(case$prepare_invocation_allowed) &&
        !isTRUE(case$runner_invocation_allowed) &&
        !isTRUE(case$artifact_invocation_allowed) &&
        !isTRUE(case$backend_invocation_allowed) &&
        !isTRUE(case$positive_matrix_allowed) &&
        !isTRUE(case$promotion_authority)
    if (!valid) auxPublicationAbort("auxiliary exclusion case differs")
    invisible(case)
}

auxPublicationValidateExclusions <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "exclusions_id", "owner_ticket_id", "status",
        "capability_binding", "coverage_binding", "cases", "unknown_policy",
        "publication_authority"
    ), "Auxiliary exclusions")
    auxPublicationValidateBinding(record$capability_binding, "Capability")
    auxPublicationValidateBinding(record$coverage_binding, "Coverage")
    lapply(record$cases, auxPublicationValidateExclusionCase)
    ids <- vapply(record$cases, `[[`, character(1), "subject_id")
    expected <- vapply(
        auxPublicationExpectedExclusions(),
        `[[`,
        character(1),
        "subject_id"
    )
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_auxiliary_exclusions"
    ) && identical(record$schema_version, .AUX_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .AUX_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        setequal(ids, expected) && length(ids) == length(expected) &&
        !anyDuplicated(ids) && identical(record$unknown_policy, "reject") &&
        !isTRUE(record$publication_authority)
    if (!valid) auxPublicationAbort("auxiliary exclusions differ")
    invisible(record)
}

auxPublicationBuildExclusions <- function() {
    cases <- lapply(auxPublicationExpectedExclusions(), function(value) {
        c(
            list(case_id = paste0(value$subject_id, ".excluded.v1")),
            value,
            list(
                expected_condition_class =
                    "multischolar_publication_auxiliary_surface_excluded",
                prepare_invocation_allowed = FALSE,
                runner_invocation_allowed = FALSE,
                artifact_invocation_allowed = FALSE,
                backend_invocation_allowed = FALSE,
                positive_matrix_allowed = FALSE,
                promotion_authority = FALSE
            )
        )
    })
    list(
        schema = "multischolar.omics_publication_auxiliary_exclusions",
        schema_version = .AUX_PUBLICATION_VERSION,
        exclusions_id = paste0(
            "multischolar.omics_publication_auxiliary_exclusions.2026-08-28.v1"
        ),
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "frozen_pre_candidate",
        capability_binding = auxPublicationBinding(
            "tests/testdata/omics-capabilities.json"
        ),
        coverage_binding = auxPublicationBinding(
            "tests/testdata/omics-performance/coverage-v1.json"
        ),
        cases = cases,
        unknown_policy = "reject",
        publication_authority = FALSE
    )
}

auxPublicationEvaluateExclusion <- function(case) {
    auxPublicationValidateExclusionCase(case)
    counts <- new.env(parent = emptyenv())
    counts$prepare <- 0L
    counts$runner <- 0L
    counts$artifact <- 0L
    counts$backend <- 0L
    spy <- function(name) function() {
        counts[[name]] <- counts[[name]] + 1L
    }
    error <- tryCatch({
        auxPublicationResolveSurface(
            case$subject_id,
            prepare = spy("prepare"),
            runner = spy("runner"),
            artifact = spy("artifact"),
            backend = spy("backend")
        )
        NULL
    }, error = function(error) error)
    valid <- inherits(error, case$expected_condition_class) &&
        identical(counts$prepare, 0L) && identical(counts$runner, 0L) &&
        identical(counts$artifact, 0L) && identical(counts$backend, 0L)
    if (!valid) auxPublicationAbort("auxiliary exclusion invocation differs")
    list(
        case_id = case$case_id,
        condition_class = class(error)[[1L]],
        prepare_invocations = counts$prepare,
        runner_invocations = counts$runner,
        artifact_invocations = counts$artifact,
        backend_invocations = counts$backend,
        promotion_authority = FALSE
    )
}
