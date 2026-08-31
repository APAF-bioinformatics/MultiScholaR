auxPublicationSeedFamilies <- function() {
    list(
        pilot = list(
            seed_family_id = "generated.pilot.200000-299999.v1",
            minimum_seed = 200000L,
            maximum_seed = 299999L
        ),
        holdout = list(
            seed_family_id = "generated.holdout.300000-399999.v1",
            minimum_seed = 300000L,
            maximum_seed = 399999L
        ),
        stress = list(
            seed_family_id = "generated.stress.400000-499999.v1",
            minimum_seed = 400000L,
            maximum_seed = 499999L
        )
    )
}

auxPublicationAssignment <- function(
    route_id,
    workload_class,
    split_role,
    seed_family_id,
    seed
) {
    route <- auxPublicationRoutes()[[route_id]]
    if (is.null(route)) auxPublicationAbort("auxiliary assignment route differs")
    list(
        assignment_id = paste(
            "auxiliary",
            route_id,
            workload_class,
            split_role,
            sep = "."
        ),
        route_id = route_id,
        surface_id = route$surface_id,
        workload_class = workload_class,
        split_role = split_role,
        source_id = paste(
            "generated.source",
            route_id,
            workload_class,
            split_role,
            "v1",
            sep = "."
        ),
        project_id = paste(
            "generated.project",
            route_id,
            workload_class,
            split_role,
            "v1",
            sep = "."
        ),
        independence_owner_id = paste(
            "generated.owner",
            route_id,
            workload_class,
            split_role,
            "v1",
            sep = "."
        ),
        response_profile_id = paste(
            "generated.response",
            route_id,
            workload_class,
            split_role,
            "v1",
            sep = "."
        ),
        generator_parameter_profile_id = paste(
            "generated.parameters",
            route_id,
            workload_class,
            split_role,
            "v1",
            sep = "."
        ),
        seed_family_id = seed_family_id,
        seed = seed,
        generated_counts_as_real = FALSE,
        promotion_authority = FALSE
    )
}

auxPublicationBuildAssignments <- function() {
    routes <- names(auxPublicationRoutes())
    assignments <- list()
    pilot <- list()
    for (route_index in seq_along(routes)) {
        route_id <- routes[[route_index]]
        assignments[[length(assignments) + 1L]] <- auxPublicationAssignment(
            route_id,
            "fixture_correctness",
            "fixture",
            paste0("fixture.no_rng.", route_id, ".v1"),
            NULL
        )
        assignments[[length(assignments) + 1L]] <- auxPublicationAssignment(
            route_id,
            "representative",
            "holdout",
            auxPublicationSeedFamilies()$holdout$seed_family_id,
            as.integer(300000L + route_index * 100L + 2L)
        )
        assignments[[length(assignments) + 1L]] <- auxPublicationAssignment(
            route_id,
            "operational_heavy",
            "holdout",
            auxPublicationSeedFamilies()$holdout$seed_family_id,
            as.integer(300000L + route_index * 100L + 3L)
        )
        assignments[[length(assignments) + 1L]] <- auxPublicationAssignment(
            route_id,
            "stress",
            "stress",
            auxPublicationSeedFamilies()$stress$seed_family_id,
            as.integer(400000L + route_index * 100L + 4L)
        )
        pilot[[length(pilot) + 1L]] <- auxPublicationAssignment(
            route_id,
            "operational_heavy",
            "pilot",
            auxPublicationSeedFamilies()$pilot$seed_family_id,
            as.integer(200000L + route_index * 100L + 1L)
        )
    }
    list(assignments = assignments, pilot = pilot)
}

auxPublicationBuildSplits <- function() {
    assignments <- auxPublicationBuildAssignments()
    list(
        schema = "multischolar.omics_publication_auxiliary_splits",
        schema_version = .AUX_PUBLICATION_VERSION,
        splits_id = paste0(
            "multischolar.omics_publication_auxiliary_splits.2026-08-28.v1"
        ),
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "frozen_pre_candidate",
        sources_binding = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/sources-v1.json"
        )),
        responses_binding = auxPublicationBinding(paste0(
            "tests/testdata/omics-performance/auxiliary/responses-v1.json"
        )),
        split_units = list(
            "source_id", "project_id", "independence_owner_id",
            "response_profile_id", "generator_parameter_profile_id",
            "seed_family_id"
        ),
        disjointness_rules = list(
            source_id_overlap = FALSE,
            project_id_overlap = FALSE,
            independence_owner_id_overlap = FALSE,
            response_profile_id_overlap = FALSE,
            generator_parameter_profile_id_overlap = FALSE,
            seed_family_id_overlap = FALSE,
            candidate_result_may_reassign_split = FALSE
        ),
        seed_families = auxPublicationSeedFamilies(),
        pilot_calibration_assignments = assignments$pilot,
        assignments = assignments$assignments,
        generated_counts_as_real = FALSE,
        candidate_access_allowed = FALSE,
        publication_authority = FALSE
    )
}

auxPublicationValidateAssignment <- function(assignment) {
    publicationRequireNames(assignment, c(
        "assignment_id", "route_id", "surface_id", "workload_class",
        "split_role", "source_id", "project_id", "independence_owner_id",
        "response_profile_id", "generator_parameter_profile_id",
        "seed_family_id", "seed", "generated_counts_as_real",
        "promotion_authority"
    ), "Auxiliary split assignment")
    route <- auxPublicationRoutes()[[assignment$route_id]]
    valid <- !is.null(route) &&
        identical(assignment$surface_id, route$surface_id) &&
        assignment$workload_class %in% c(
            "fixture_correctness", "representative", "operational_heavy",
            "stress"
        ) &&
        assignment$split_role %in% c("fixture", "pilot", "holdout", "stress") &&
        !isTRUE(assignment$generated_counts_as_real) &&
        !isTRUE(assignment$promotion_authority)
    if (identical(assignment$split_role, "fixture")) {
        valid <- valid && is.null(assignment$seed) &&
            startsWith(assignment$seed_family_id, "fixture.no_rng.")
    } else {
        families <- auxPublicationSeedFamilies()
        family <- families[[assignment$split_role]]
        valid <- valid && !is.null(family) &&
            identical(assignment$seed_family_id, family$seed_family_id) &&
            is.numeric(assignment$seed) && length(assignment$seed) == 1L &&
            assignment$seed >= family$minimum_seed &&
            assignment$seed <= family$maximum_seed
    }
    if (!valid) auxPublicationAbort("auxiliary split assignment differs")
    invisible(assignment)
}

auxPublicationSplitValues <- function(assignments, field) {
    vapply(assignments, `[[`, character(1), field)
}

auxPublicationValidateSplits <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "splits_id", "owner_ticket_id", "status",
        "sources_binding", "responses_binding", "split_units",
        "disjointness_rules", "seed_families", "pilot_calibration_assignments",
        "assignments", "generated_counts_as_real", "candidate_access_allowed",
        "publication_authority"
    ), "Auxiliary split authority")
    auxPublicationValidateBinding(record$sources_binding, "Sources")
    auxPublicationValidateBinding(record$responses_binding, "Responses")
    lapply(record$assignments, auxPublicationValidateAssignment)
    lapply(
        record$pilot_calibration_assignments,
        auxPublicationValidateAssignment
    )
    expected <- auxPublicationBuildAssignments()
    assignment_ids <- auxPublicationSplitValues(
        record$assignments,
        "assignment_id"
    )
    pilot_ids <- auxPublicationSplitValues(
        record$pilot_calibration_assignments,
        "assignment_id"
    )
    disjoint_fields <- c(
        "source_id", "project_id", "independence_owner_id",
        "response_profile_id", "generator_parameter_profile_id"
    )
    disjoint <- all(vapply(disjoint_fields, function(field) {
        !length(intersect(
            auxPublicationSplitValues(record$assignments, field),
            auxPublicationSplitValues(record$pilot_calibration_assignments, field)
        ))
    }, logical(1)))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_auxiliary_splits"
    ) && identical(record$schema_version, .AUX_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .AUX_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        identical(record$seed_families, auxPublicationSeedFamilies()) &&
        identical(record$assignments, expected$assignments) &&
        identical(record$pilot_calibration_assignments, expected$pilot) &&
        length(assignment_ids) == 20L && !anyDuplicated(assignment_ids) &&
        length(pilot_ids) == 5L && !anyDuplicated(pilot_ids) && disjoint &&
        !isTRUE(record$generated_counts_as_real) &&
        !isTRUE(record$candidate_access_allowed) &&
        !isTRUE(record$publication_authority)
    if (!valid) auxPublicationAbort("auxiliary split authority differs")
    invisible(record)
}

auxPublicationFindAssignment <- function(splits, route_id, workload_class) {
    matches <- Filter(function(assignment) {
        identical(assignment$route_id, route_id) &&
            identical(assignment$workload_class, workload_class) &&
            !identical(assignment$split_role, "pilot")
    }, splits$assignments)
    if (length(matches) != 1L) {
        auxPublicationAbort("auxiliary assignment is missing or duplicated")
    }
    matches[[1L]]
}
