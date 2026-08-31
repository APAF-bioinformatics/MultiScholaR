lipidPublicationForbiddenScaleFields <- function(value) {
    if (!is.list(value)) return(invisible(TRUE))
    forbidden <- c(
        "path", "values", "distribution", "quantiles", "identifiers",
        "sequences", "unique_run_count", "unique_protein_group_count",
        "unique_peptide_count", "assay_mix", "sample_count", "lipid_classes",
        "adducts", "annotations"
    )
    if (length(intersect(names(value), forbidden))) {
        lipidPublicationAbort("lipidomics source leaks cross-omic fields")
    }
    lapply(value, lipidPublicationForbiddenScaleFields)
    invisible(TRUE)
}

lipidPublicationValidateSources <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "sources_id", "owner_ticket_id", "status",
        "projects_predecessor", "private_scale_predecessor", "allowed_scale_fields",
        "scale_receipt", "scale_mapping_policy", "minimum_real_projects",
        "sources", "decisions",
        "repository_searches", "generated_counts_as_real",
        "fixtures_count_as_real", "unknown_policy", "publication_authority"
    ), "Lipidomics source authority")
    lipidPublicationValidateBinding(
        record$projects_predecessor,
        "Projects predecessor"
    )
    lipidPublicationValidateBinding(
        record$private_scale_predecessor,
        "Private scale predecessor"
    )
    expected_fields <- list(
        "row_count", "column_count", "byte_size",
        "salted_source_fingerprint"
    )
    private <- publicationReadJson(record$private_scale_predecessor$path)
    publicationRequireNames(
        record$scale_receipt,
        unlist(expected_fields, use.names = FALSE),
        "Lipidomics scale receipt"
    )
    scale_valid <- identical(record$scale_receipt$row_count, private$report$row_count) &&
        identical(record$scale_receipt$column_count, private$report$column_count) &&
        identical(record$scale_receipt$byte_size, private$report$byte_size) &&
        identical(
            record$scale_receipt$salted_source_fingerprint,
            private$report$salted_source_fingerprint
        )
    routes <- vapply(record$decisions, `[[`, character(1), "route")
    valid_decisions <- all(vapply(record$decisions, \(decision) {
        setequal(names(decision), c(
            "route", "capability_id", "required_real_project_count",
            "verified_real_project_count", "claim_scope",
            "cross_project_source_ready", "promotion_eligible"
        )) && decision$route %in% names(lipidPublicationCapabilities()) &&
            identical(
                decision$capability_id,
                lipidPublicationCapabilities()[[decision$route]]$capability_id
            ) && identical(decision$required_real_project_count, 3L) &&
            identical(decision$verified_real_project_count, 0L) &&
            identical(decision$claim_scope, "project_specific_nonpromotional") &&
            !isTRUE(decision$cross_project_source_ready) &&
            !isTRUE(decision$promotion_eligible)
    }, logical(1)))
    valid_searches <- all(vapply(record$repository_searches, \(search) {
        setequal(names(search), c(
            "route", "repository", "query", "admissible_source_count",
            "outcome"
        )) && search$route %in% names(lipidPublicationCapabilities()) &&
            publicationScalarString(search$repository) &&
            publicationScalarString(search$query) &&
            identical(search$admissible_source_count, 0L) &&
            identical(search$outcome, "no_admissible_format_exact_source_frozen")
    }, logical(1)))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_lipidomics_sources"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-067") &&
        identical(record$status, "frozen_nonpromotional") &&
        identical(record$allowed_scale_fields, expected_fields) &&
        scale_valid &&
        identical(
            record$scale_mapping_policy,
            "proteomics_shape_to_bounded_lipidomics_dimensions_only_v1"
        ) && identical(record$minimum_real_projects, 3L) &&
        !length(record$sources) && setequal(
            routes,
            names(lipidPublicationCapabilities())
        ) && !anyDuplicated(routes) && valid_decisions &&
        length(record$repository_searches) == 3L && valid_searches &&
        !isTRUE(record$generated_counts_as_real) &&
        !isTRUE(record$fixtures_count_as_real) &&
        identical(record$unknown_policy, "reject") &&
        !isTRUE(record$publication_authority)
    lipidPublicationForbiddenScaleFields(record$scale_receipt)
    if (!valid) lipidPublicationAbort("lipidomics source authority differs")
    invisible(record)
}
