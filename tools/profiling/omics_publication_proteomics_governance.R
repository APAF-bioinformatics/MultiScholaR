proteomicsPublicationGovernanceBindingFields <- function() {
    c(
        "projects_predecessor", "splits_predecessor", "sources", "splits",
        "parameters", "exclusions", "private_envelope", "public_acquisition"
    )
}

proteomicsPublicationValidateGovernanceSuccessor <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "successor_id", "owner_ticket_id",
        "status", proteomicsPublicationGovernanceBindingFields(),
        "supersession", "claim_scope", "candidate_access_allowed",
        "promotion_authority", "publication_authority"
    ), "Proteomics governance successor")
    for (field in proteomicsPublicationGovernanceBindingFields()) {
        proteomicsPublicationValidateBinding(record[[field]], field)
    }
    projects <- publicationReadJson(record$projects_predecessor$path)
    predecessor_splits <- publicationReadJson(record$splits_predecessor$path)
    sources <- publicationReadJson(record$sources$path)
    splits <- publicationReadJson(record$splits$path)
    parameters <- publicationReadJson(record$parameters$path)
    exclusions <- publicationReadJson(record$exclusions$path)
    private <- publicationReadJson(record$private_envelope$path)
    acquisition <- publicationReadJson(record$public_acquisition$path)
    proteomicsPublicationValidateSources(sources)
    proteomicsPublicationValidateSplits(splits)
    proteomicsPublicationValidateParameters(parameters)
    proteomicsPublicationValidateExclusions(exclusions)
    proteomicsPublicationValidatePrivateEnvelope(private)
    proteomicsPublicationValidateAcquisition(acquisition)
    capability_ids <- vapply(
        proteomicsPublicationCapabilities(),
        `[[`,
        character(1),
        "capability_id"
    )
    publicationRequireNames(record$supersession, c(
        "field_scope", "owned_capability_ids", "predecessors_mutated",
        "unrelated_subjects_changed", "generated_substitution_allowed",
        "successor_required_before_candidate"
    ), "Proteomics governance supersession")
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_proteomics_governance_successor"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-065") &&
        identical(record$status, "proteomics_scoped_nonpromotional") &&
        identical(projects$owner_ticket_id, "OMICS-ART-062") &&
        identical(predecessor_splits$owner_ticket_id, "OMICS-ART-062") &&
        identical(
            record$supersession$field_scope,
            as.list(c("project_sources", "split_assignments"))
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
    if (!valid) proteomicsPublicationAbort("governance successor differs")
    invisible(record)
}

proteomicsPublicationHandoffConsumers <- function() {
    list(
        "OMICS-ART-070" = list(
            formats = list("diann"),
            purpose = "dia_commit_repair"
        ),
        "OMICS-ART-071" = list(
            formats = list("diann"),
            purpose = "dia_bounded_ingress"
        ),
        "OMICS-ART-072" = list(
            formats = list("diann"),
            purpose = "dia_consumer_completion"
        ),
        "OMICS-ART-073" = list(
            formats = list("maxquant", "fragpipe"),
            purpose = "lfq_runtime_repair"
        ),
        "OMICS-ART-074" = list(
            formats = list("pd_tmt"),
            purpose = "tmt_runtime_repair"
        ),
        "OMICS-ART-078" = list(
            formats = list("diann"),
            purpose = "dia_confirmatory_benchmark"
        ),
        "OMICS-ART-079" = list(
            formats = list("maxquant", "fragpipe", "pd_tmt"),
            purpose = "non_dia_confirmatory_benchmark"
        )
    )
}

proteomicsPublicationValidateHandoff <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "handoff_id", "owner_ticket_id", "status",
        "producer_bindings", "consumers", "engineering_control_alias_allowed",
        "stale_or_missing_policy", "promotion_authority",
        "publication_authority"
    ), "Proteomics handoff")
    lapply(record$producer_bindings, \(binding) {
        proteomicsPublicationValidateBinding(binding, "Handoff producer")
    })
    expected <- proteomicsPublicationHandoffConsumers()
    ids <- vapply(record$consumers, `[[`, character(1), "ticket_id")
    valid_consumers <- all(vapply(record$consumers, \(consumer) {
        contract <- expected[[consumer$ticket_id]]
        setequal(names(consumer), c(
            "ticket_id", "formats", "purpose", "required_workload_classes",
            "unavailable_outcome"
        )) && !is.null(contract) &&
            identical(consumer$formats, contract$formats) &&
            identical(consumer$purpose, contract$purpose) &&
            identical(
                consumer$required_workload_classes,
                as.list(proteomicsPublicationClasses())
            ) && publicationScalarString(consumer$unavailable_outcome)
    }, logical(1)))
    valid <- identical(
        record$schema,
        "multischolar.omics_publication_proteomics_handoff"
    ) && identical(record$schema_version, "1.0.0") &&
        identical(record$owner_ticket_id, "OMICS-ART-065") &&
        identical(record$status, "frozen_pre_candidate") &&
        setequal(ids, names(expected)) && !anyDuplicated(ids) && valid_consumers &&
        !isTRUE(record$engineering_control_alias_allowed) &&
        identical(record$stale_or_missing_policy, "reject") &&
        !isTRUE(record$promotion_authority) &&
        !isTRUE(record$publication_authority)
    if (!valid) proteomicsPublicationAbort("handoff authority differs")
    invisible(record)
}
