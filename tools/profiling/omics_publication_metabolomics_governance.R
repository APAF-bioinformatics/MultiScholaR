metabPublicationSplitUnitFields <- function() {
    c(
        "source_id", "project_id", "independence_owner_id",
        "aggregate_profile_id", "generator_parameter_profile_id",
        "seed_family_id"
    )
}

metabPublicationValidateSplitAssignment <- function(
    assignment,
    seed_families,
    calibration = FALSE
) {
    fields <- c(
        "assignment_id", "capability_id", "route", "profile_id",
        "workload_class", "split_role", metabPublicationSplitUnitFields(),
        "seed"
    )
    publicationRequireNames(assignment, fields, "Metabolomics split assignment")
    capability <- metabPublicationCapabilities()[[assignment$route]]
    expected_id <- metabPublicationWorkloadId(
        assignment$route,
        assignment$profile_id,
        assignment$workload_class
    )
    if (isTRUE(calibration)) expected_id <- paste0(expected_id, ".pilot_calibration")
    valid <- !is.null(capability) &&
        identical(assignment$capability_id, capability$capability_id) &&
        identical(assignment$assignment_id, expected_id) &&
        assignment$profile_id %in% names(metabPublicationAssayProfiles()) &&
        assignment$workload_class %in% metabPublicationClasses() &&
        all(vapply(
            assignment[metabPublicationSplitUnitFields()],
            publicationScalarString,
            logical(1)
        ))
    if (isTRUE(calibration)) {
        valid <- valid && identical(assignment$profile_id, "mixed") &&
            identical(assignment$workload_class, "operational_heavy") &&
            identical(assignment$split_role, "pilot")
    } else if (identical(assignment$workload_class, "fixture_correctness")) {
        valid <- valid && identical(assignment$split_role, "fixture") &&
            is.null(assignment$seed)
    } else if (identical(assignment$workload_class, "stress")) {
        valid <- valid && identical(assignment$split_role, "stress")
    } else {
        valid <- valid && identical(assignment$split_role, "holdout")
    }
    if (!identical(assignment$split_role, "fixture")) {
        family <- seed_families[[assignment$seed_family_id]]
        valid <- valid && !is.null(family) && is.numeric(assignment$seed) &&
            length(assignment$seed) == 1L && is.finite(assignment$seed) &&
            assignment$seed == as.integer(assignment$seed) &&
            assignment$seed >= family$minimum_seed &&
            assignment$seed <= family$maximum_seed
    }
    if (!valid) metabPublicationAbort("metabolomics split assignment differs")
    invisible(assignment)
}

metabPublicationSplitValues <- function(assignments, field) {
    vapply(assignments, `[[`, character(1), field)
}

metabPublicationValidateSplitDisjointness <- function(pilot, final, rules) {
    expected <- stats::setNames(
        as.list(rep(FALSE, length(metabPublicationSplitUnitFields()))),
        paste0(metabPublicationSplitUnitFields(), "_overlap")
    )
    expected$candidate_result_may_reassign_split <- FALSE
    if (!identical(rules, expected)) {
        metabPublicationAbort("metabolomics split disjointness rules differ")
    }
    for (field in metabPublicationSplitUnitFields()) {
        if (length(intersect(
            metabPublicationSplitValues(pilot, field),
            metabPublicationSplitValues(final, field)
        ))) {
            metabPublicationAbort(paste("metabolomics split overlap:", field))
        }
    }
    invisible(TRUE)
}

metabPublicationValidateSplits <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "splits_id", "owner_ticket_id",
        "status", "splits_predecessor", "sources_binding", "split_units",
        "disjointness_rules", "seed_families", "pilot_calibration_assignments",
        "assignments", "readiness", "publication_authority"
    ), "Metabolomics split authority")
    metabPublicationValidateBinding(record$splits_predecessor, "Splits predecessor")
    metabPublicationValidateBinding(record$sources_binding, "Sources")
    predecessor <- publicationReadJson(record$splits_predecessor$path)
    expected_families <- predecessor$generated_seed_families
    names(expected_families) <- vapply(
        expected_families,
        `[[`,
        character(1),
        "seed_family_id"
    )
    families <- record$seed_families
    names(families) <- vapply(
        families,
        `[[`,
        character(1),
        "seed_family_id"
    )
    if (!identical(families, expected_families)) {
        metabPublicationAbort("metabolomics seed family authority differs")
    }
    lapply(
        record$pilot_calibration_assignments,
        metabPublicationValidateSplitAssignment,
        seed_families = families,
        calibration = TRUE
    )
    lapply(
        record$assignments,
        metabPublicationValidateSplitAssignment,
        seed_families = families,
        calibration = FALSE
    )
    all_assignments <- c(
        record$pilot_calibration_assignments,
        record$assignments
    )
    assignment_ids <- vapply(
        all_assignments,
        `[[`,
        character(1),
        "assignment_id"
    )
    final_ids <- vapply(record$assignments, function(assignment) {
        metabPublicationWorkloadId(
            assignment$route,
            assignment$profile_id,
            assignment$workload_class
        )
    }, character(1))
    expected_ids <- unlist(lapply(
        names(metabPublicationCapabilities()),
        function(route) unlist(lapply(
            names(metabPublicationAssayProfiles()),
            function(profile_id) vapply(
                metabPublicationClasses(),
                function(workload_class) metabPublicationWorkloadId(
                    route,
                    profile_id,
                    workload_class
                ),
                character(1)
            )
        ), use.names = FALSE)
    ), use.names = FALSE)
    metabPublicationValidateSplitDisjointness(
        record$pilot_calibration_assignments,
        record$assignments,
        record$disjointness_rules
    )
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_metabolomics_splits"
    ) && identical(record$schema_version, .METAB_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .METAB_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        identical(record$split_units, as.list(metabPublicationSplitUnitFields())) &&
        length(record$pilot_calibration_assignments) == 2L &&
        length(record$assignments) == 32L && !anyDuplicated(assignment_ids) &&
        setequal(final_ids, expected_ids) && !anyDuplicated(final_ids) &&
        !isTRUE(record$readiness$real_project_sources_ready) &&
        isTRUE(record$readiness$generated_assignments_complete) &&
        isTRUE(record$readiness$successor_required_before_candidate) &&
        !isTRUE(record$readiness$candidate_access_allowed) &&
        !isTRUE(record$publication_authority)
    if (!valid) metabPublicationAbort("metabolomics split authority differs")
    invisible(record)
}

metabPublicationSplitAssignment <- function(
    route,
    profile_id,
    workload_class,
    split_role,
    seed,
    calibration = FALSE
) {
    workload_id <- metabPublicationWorkloadId(route, profile_id, workload_class)
    assignment_id <- if (isTRUE(calibration)) {
        paste0(workload_id, ".pilot_calibration")
    } else {
        workload_id
    }
    label <- paste(route, profile_id, workload_class, split_role, sep = ".")
    family <- switch(
        split_role,
        fixture = paste("fixture.no_rng", route, profile_id, "v1", sep = "."),
        pilot = "generated.pilot.200000-299999.v1",
        holdout = "generated.holdout.300000-399999.v1",
        stress = "generated.stress.400000-499999.v1"
    )
    list(
        assignment_id = assignment_id,
        capability_id = metabPublicationCapabilities()[[route]]$capability_id,
        route = route,
        profile_id = profile_id,
        workload_class = workload_class,
        split_role = split_role,
        source_id = paste0("generated.source.", label, ".v1"),
        project_id = paste0("generated.project.", label, ".v1"),
        independence_owner_id = paste0("generated.owner.", label, ".v1"),
        aggregate_profile_id = paste0("generated.aggregate.", label, ".v1"),
        generator_parameter_profile_id = paste0(
            "generated.parameters.",
            label,
            ".v1"
        ),
        seed_family_id = family,
        seed = seed
    )
}

metabPublicationBuildSplits <- function() {
    predecessor_path <- "tests/testdata/omics-performance/splits-v1.json"
    sources_path <- "tests/testdata/omics-performance/metabolomics/sources-v1.json"
    predecessor <- publicationReadJson(predecessor_path)
    routes <- names(metabPublicationCapabilities())
    profiles <- names(metabPublicationAssayProfiles())
    pilot <- lapply(seq_along(routes), function(route_index) {
        metabPublicationSplitAssignment(
            routes[[route_index]],
            "mixed",
            "operational_heavy",
            "pilot",
            as.integer(200000L + route_index * 1000L + 401L),
            calibration = TRUE
        )
    })
    assignments <- list()
    for (route_index in seq_along(routes)) {
        for (profile_index in seq_along(profiles)) {
            for (workload_class in metabPublicationClasses()) {
                split_role <- switch(
                    workload_class,
                    fixture_correctness = "fixture",
                    stress = "stress",
                    "holdout"
                )
                base <- route_index * 10000L + profile_index * 100L
                offset <- match(workload_class, metabPublicationClasses())
                seed <- switch(
                    split_role,
                    fixture = NULL,
                    stress = as.integer(400000L + base + offset),
                    as.integer(300000L + base + offset)
                )
                assignments[[length(assignments) + 1L]] <-
                    metabPublicationSplitAssignment(
                        routes[[route_index]],
                        profiles[[profile_index]],
                        workload_class,
                        split_role,
                        seed
                    )
            }
        }
    }
    rules <- stats::setNames(
        as.list(rep(FALSE, length(metabPublicationSplitUnitFields()))),
        paste0(metabPublicationSplitUnitFields(), "_overlap")
    )
    rules$candidate_result_may_reassign_split <- FALSE
    list(
        schema = "multischolar.omics_publication_metabolomics_splits",
        schema_version = .METAB_PUBLICATION_VERSION,
        splits_id = paste0(
            "multischolar.omics_publication_metabolomics_splits.2026-08-28.v1"
        ),
        owner_ticket_id = .METAB_PUBLICATION_OWNER,
        status = "frozen_pre_candidate",
        splits_predecessor = list(
            path = predecessor_path,
            sha256 = publicationFileDigest(predecessor_path)
        ),
        sources_binding = list(
            path = sources_path,
            sha256 = publicationFileDigest(sources_path)
        ),
        split_units = as.list(metabPublicationSplitUnitFields()),
        disjointness_rules = rules,
        seed_families = predecessor$generated_seed_families,
        pilot_calibration_assignments = pilot,
        assignments = assignments,
        readiness = list(
            real_project_sources_ready = FALSE,
            generated_assignments_complete = TRUE,
            successor_required_before_candidate = TRUE,
            candidate_access_allowed = FALSE
        ),
        publication_authority = FALSE
    )
}

metabPublicationBuildMsdialBundles <- function() {
    bundles <- lapply(metabPublicationClasses(), function(workload_class) {
        bundle <- list(
            bundle_id = paste(
                "metabolomics.msdial.mixed",
                workload_class,
                "v1",
                sep = "."
            ),
            profile_id = "mixed",
            member_assays = metabPublicationAssayProfiles()$mixed$assays,
            member_count = 3L,
            member_schema_ids = list(
                "msdial.lcms_pos.v1", "msdial.lcms_neg.v1", "msdial.gcms.v1"
            ),
            single_custom_table_substitution_allowed = FALSE,
            bundle_digest = strrep("0", 64L)
        )
        bundle$bundle_digest <- publicationObjectDigest(bundle)
        bundle
    })
    mapping_path <- paste0(
        "tests/testdata/omics-performance/metabolomics/",
        "mapping-authorities-v1.json"
    )
    list(
        schema = "multischolar.omics_publication_metabolomics_msdial_bundles",
        schema_version = .METAB_PUBLICATION_VERSION,
        bundles_id = paste0(
            "multischolar.omics_publication_metabolomics_bundles.2026-08-28.v1"
        ),
        owner_ticket_id = .METAB_PUBLICATION_OWNER,
        status = "frozen_pre_candidate",
        mapping_authority = list(
            path = mapping_path,
            sha256 = publicationFileDigest(mapping_path)
        ),
        bundles = bundles,
        member_order_semantic = FALSE,
        missing_or_duplicate_member_policy = "reject",
        custom_substitution_allowed = FALSE,
        publication_authority = FALSE
    )
}

metabPublicationGovernanceBindingFields <- function() {
    c(
        "projects_predecessor", "splits_predecessor", "sources", "splits",
        "parameters", "support_boundaries", "exclusions",
        "mapping_authorities", "msdial_bundles"
    )
}

metabPublicationValidateGovernanceSuccessor <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "successor_id", "owner_ticket_id",
        "status", metabPublicationGovernanceBindingFields(), "supersession",
        "claim_scope", "candidate_access_allowed", "promotion_authority",
        "publication_authority"
    ), "Metabolomics governance successor")
    for (field in metabPublicationGovernanceBindingFields()) {
        metabPublicationValidateBinding(record[[field]], field)
    }
    metabPublicationValidateSources(publicationReadJson(record$sources$path))
    metabPublicationValidateSplits(publicationReadJson(record$splits$path))
    metabPublicationValidateParameters(publicationReadJson(record$parameters$path))
    metabPublicationValidateSupportBoundaries(
        publicationReadJson(record$support_boundaries$path)
    )
    metabPublicationValidateExclusions(publicationReadJson(record$exclusions$path))
    metabPublicationValidateMappingAuthority(
        publicationReadJson(record$mapping_authorities$path)
    )
    metabPublicationValidateMsdialBundles(
        publicationReadJson(record$msdial_bundles$path)
    )
    capability_ids <- vapply(
        metabPublicationCapabilities(),
        `[[`,
        character(1),
        "capability_id"
    )
    publicationRequireNames(record$supersession, c(
        "field_scope", "owned_capability_ids", "predecessors_mutated",
        "unrelated_subjects_changed", "generated_substitution_allowed",
        "successor_required_before_candidate"
    ), "Metabolomics governance supersession")
    projects <- publicationReadJson(record$projects_predecessor$path)
    splits <- publicationReadJson(record$splits_predecessor$path)
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_metabolomics_governance_successor"
    ) && identical(record$schema_version, .METAB_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .METAB_PUBLICATION_OWNER) &&
        identical(record$status, "metabolomics_scoped_nonpromotional") &&
        identical(projects$owner_ticket_id, "OMICS-ART-062") &&
        identical(splits$owner_ticket_id, "OMICS-ART-062") &&
        identical(
            record$supersession$field_scope,
            as.list(c("project_sources", "split_assignments", "route_boundaries"))
        ) && setequal(
            unlist(record$supersession$owned_capability_ids, use.names = FALSE),
            capability_ids
        ) && !isTRUE(record$supersession$predecessors_mutated) &&
        !isTRUE(record$supersession$unrelated_subjects_changed) &&
        !isTRUE(record$supersession$generated_substitution_allowed) &&
        isTRUE(record$supersession$successor_required_before_candidate) &&
        identical(record$claim_scope, "project_specific_nonpromotional") &&
        !isTRUE(record$candidate_access_allowed) &&
        !isTRUE(record$promotion_authority) &&
        !isTRUE(record$publication_authority)
    if (!valid) metabPublicationAbort("metabolomics governance successor differs")
    invisible(record)
}

metabPublicationHandoffConsumers <- function() {
    list(
        "OMICS-ART-075" = "metabolomics_runtime_repair",
        "OMICS-ART-080" = "metabolomics_confirmatory_benchmark"
    )
}

metabPublicationValidateHandoff <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "handoff_id", "owner_ticket_id", "status",
        "producer_bindings", "consumers", "route_alias_allowed",
        "stale_or_missing_policy", "promotion_authority", "publication_authority"
    ), "Metabolomics handoff")
    lapply(record$producer_bindings, function(binding) {
        metabPublicationValidateBinding(binding, "Metabolomics handoff producer")
    })
    expected <- metabPublicationHandoffConsumers()
    ids <- vapply(record$consumers, `[[`, character(1), "ticket_id")
    valid_consumers <- all(vapply(record$consumers, function(consumer) {
        setequal(names(consumer), c(
            "ticket_id", "purpose", "routes", "required_workload_classes",
            "required_assay_profiles", "unavailable_outcome"
        )) && identical(consumer$purpose, expected[[consumer$ticket_id]]) &&
            setequal(
                unlist(consumer$routes),
                names(metabPublicationCapabilities())
            ) && identical(
                consumer$required_workload_classes,
                as.list(metabPublicationClasses())
            ) && identical(
                consumer$required_assay_profiles,
                as.list(names(metabPublicationAssayProfiles()))
            ) && publicationScalarString(consumer$unavailable_outcome)
    }, logical(1)))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_metabolomics_handoff"
    ) && identical(record$schema_version, .METAB_PUBLICATION_VERSION) &&
        identical(record$owner_ticket_id, .METAB_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_pre_candidate") &&
        setequal(ids, names(expected)) && !anyDuplicated(ids) && valid_consumers &&
        !isTRUE(record$route_alias_allowed) &&
        identical(record$stale_or_missing_policy, "reject") &&
        !isTRUE(record$promotion_authority) &&
        !isTRUE(record$publication_authority)
    if (!valid) metabPublicationAbort("metabolomics handoff differs")
    invisible(record)
}
