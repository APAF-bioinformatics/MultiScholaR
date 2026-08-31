auxPublicationSourceSearches <- function() {
    list(
        list(
            surface_id = "phosphosite.api",
            repository = "PRIDE Archive",
            metadata_url = "https://www.ebi.ac.uk/pride/archive/",
            query_scope = paste0(
                "phosphoproteomics projects with processed evidence, ",
                "FASTA, design, license, and independent owner metadata"
            ),
            admissible_source_count = 0L,
            outcome = "no_exact_admissible_source_frozen"
        ),
        list(
            surface_id = "phosphosite.api",
            repository = "ProteomeXchange",
            metadata_url = "https://www.proteomexchange.org/",
            query_scope = paste0(
                "phosphosite processed-result projects with complete ",
                "redistribution and design authority"
            ),
            admissible_source_count = 0L,
            outcome = "no_exact_admissible_source_frozen"
        ),
        list(
            surface_id = "phosphosite.api",
            repository = "OmicsDI",
            metadata_url = "https://www.omicsdi.org/",
            query_scope = paste0(
                "independently owned phosphoproteomics projects with ",
                "processable evidence and FASTA pairs"
            ),
            admissible_source_count = 0L,
            outcome = "no_exact_admissible_source_frozen"
        ),
        list(
            surface_id = "multiomics.api",
            repository = "MetaboLights",
            metadata_url = "https://www.ebi.ac.uk/metabolights/",
            query_scope = paste0(
                "multiomics projects with proteomics, metabolomics, ",
                "sample linkage, layer identities, and redistribution authority"
            ),
            admissible_source_count = 0L,
            outcome = "no_exact_admissible_source_frozen"
        ),
        list(
            surface_id = "multiomics.api",
            repository = "OmicsDI",
            metadata_url = "https://www.omicsdi.org/",
            query_scope = paste0(
                "multiomics projects with at least three linked layers ",
                "and exact sample-level joins"
            ),
            admissible_source_count = 0L,
            outcome = "no_exact_admissible_source_frozen"
        ),
        list(
            surface_id = "multiomics.api",
            repository = "PRIDE Archive",
            metadata_url = "https://www.ebi.ac.uk/pride/archive/",
            query_scope = paste0(
                "proteomics projects with repository-linked metabolomics ",
                "and independently verifiable sample correspondence"
            ),
            admissible_source_count = 0L,
            outcome = "no_exact_admissible_source_frozen"
        )
    )
}

auxPublicationSourceDecisions <- function() {
    lapply(c("phosphosite.api", "multiomics.api"), function(surface_id) {
        list(
            surface_id = surface_id,
            required_real_project_count = 3L,
            verified_real_project_count = 0L,
            cross_project_claim_status = "insufficient_independent_projects",
            current_claim_scope = "generated_local_nonpromotional",
            cross_project_source_ready = FALSE,
            promotion_eligible = FALSE
        )
    })
}

auxPublicationBuildSourceAuthority <- function() {
    list(
        schema = "multischolar.omics_publication_auxiliary_sources",
        schema_version = .AUX_PUBLICATION_VERSION,
        sources_id = paste0(
            "multischolar.omics_publication_auxiliary_sources.2026-08-28.v1"
        ),
        owner_ticket_id = .AUX_PUBLICATION_OWNER,
        status = "frozen_nonpromotional",
        surface_binding = auxPublicationBinding(
            paste0(
                "tests/testdata/omics-performance/auxiliary/",
                "surfaces-v1.json"
            )
        ),
        minimum_real_projects_per_surface = 3L,
        sources = list(),
        decisions = auxPublicationSourceDecisions(),
        repository_searches = auxPublicationSourceSearches(),
        generated_counts_as_real = FALSE,
        fixtures_count_as_real = FALSE,
        metadata_search_counts_as_real = FALSE,
        private_source_opened = FALSE,
        unknown_policy = "reject",
        publication_authority = FALSE
    )
}

auxPublicationValidateSourceSearch <- function(search) {
    publicationRequireNames(search, c(
        "surface_id", "repository", "metadata_url", "query_scope",
        "admissible_source_count", "outcome"
    ), "Auxiliary source search")
    expected <- Filter(function(value) {
        identical(value$surface_id, search$surface_id) &&
            identical(value$repository, search$repository)
    }, auxPublicationSourceSearches())
    valid <- length(expected) == 1L && identical(search, expected[[1L]]) &&
        startsWith(search$metadata_url, "https://") &&
        identical(search$admissible_source_count, 0L) &&
        identical(search$outcome, "no_exact_admissible_source_frozen")
    if (!valid) auxPublicationAbort("auxiliary source search differs")
    invisible(search)
}

auxPublicationValidateSources <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "sources_id", "owner_ticket_id", "status",
        "surface_binding", "minimum_real_projects_per_surface", "sources",
        "decisions", "repository_searches", "generated_counts_as_real",
        "fixtures_count_as_real", "metadata_search_counts_as_real",
        "private_source_opened", "unknown_policy", "publication_authority"
    ), "Auxiliary source authority")
    auxPublicationValidateBinding(record$surface_binding, "Surface")
    lapply(record$repository_searches, auxPublicationValidateSourceSearch)
    expected_decisions <- auxPublicationSourceDecisions()
    ids <- vapply(record$decisions, `[[`, character(1), "surface_id")
    search_keys <- vapply(record$repository_searches, function(value) {
        paste(value$surface_id, value$repository, sep = "::")
    }, character(1))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_auxiliary_sources"
    ) && identical(record$schema_version, .AUX_PUBLICATION_VERSION) &&
        identical(
            record$sources_id,
            "multischolar.omics_publication_auxiliary_sources.2026-08-28.v1"
        ) &&
        identical(record$owner_ticket_id, .AUX_PUBLICATION_OWNER) &&
        identical(record$status, "frozen_nonpromotional") &&
        identical(record$minimum_real_projects_per_surface, 3L) &&
        !length(record$sources) &&
        identical(record$decisions, expected_decisions) &&
        identical(ids, c("phosphosite.api", "multiomics.api")) &&
        length(search_keys) == 6L && !anyDuplicated(search_keys) &&
        !isTRUE(record$generated_counts_as_real) &&
        !isTRUE(record$fixtures_count_as_real) &&
        !isTRUE(record$metadata_search_counts_as_real) &&
        !isTRUE(record$private_source_opened) &&
        identical(record$unknown_policy, "reject") &&
        !isTRUE(record$publication_authority)
    if (!valid) auxPublicationAbort("auxiliary source authority differs")
    invisible(record)
}
