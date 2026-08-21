# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

ARTIFACT_DESCRIPTOR_SCHEMA <- "multischolar.artifact_workflow_descriptor"
ARTIFACT_DESCRIPTOR_SCHEMA_VERSION <- 1L
ARTIFACT_DESCRIPTOR_CATALOGUE_SCHEMA <- "multischolar.artifact_descriptor_catalogue"
ARTIFACT_DESCRIPTOR_CATALOGUE_VERSION <- 1L
.ARTIFACT_CERTIFICATION_LEVELS <- c(
    "uncertified", "dual_write", "read_through", "evict"
)

artifactDescriptorAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_artifact_descriptor_error"),
        ...
    )
}

artifactDescriptorSemver <- function(value, field) {
    if (!workflowCapabilityScalarString(value) ||
        !grepl("^[0-9]+[.][0-9]+[.][0-9]+$", value)) {
        artifactDescriptorAbort(
            sprintf("artifact descriptor field '%s' must be a semantic version", field),
            "multischolar_invalid_artifact_descriptor"
        )
    }
    value
}

artifactDescriptorNamedList <- function(value, owner, allow_empty = FALSE) {
    if (isTRUE(allow_empty) && is.list(value) && length(value) == 0L) {
        return(value)
    }
    valid <- is.list(value) && !is.null(names(value)) &&
        !any(!nzchar(names(value))) && anyDuplicated(names(value)) == 0L
    if (!isTRUE(allow_empty) && length(value) == 0L) valid <- FALSE
    if (!isTRUE(valid)) {
        artifactDescriptorAbort(
            sprintf("artifact descriptor '%s' must be an unambiguous named list", owner),
            "multischolar_invalid_artifact_descriptor"
        )
    }
    value
}

artifactDescriptorAssertDataOnly <- function(value, owner = "descriptor") {
    unsafe <- is.environment(value) || is.function(value) || isS4(value) ||
        is.language(value) || is.expression(value) ||
        typeof(value) %in% c("externalptr", "weakref") || inherits(value, "connection")
    if (isTRUE(unsafe)) {
        artifactDescriptorAbort(
            sprintf("artifact %s contains runtime code or process state", owner),
            "multischolar_artifact_descriptor_code_injection"
        )
    }
    value_attributes <- attributes(value)
    inspected_attributes <- setdiff(names(value_attributes), c("names", "class"))
    for (attribute_name in inspected_attributes) {
        artifactDescriptorAssertDataOnly(
            value_attributes[[attribute_name]],
            paste0(owner, "@", attribute_name)
        )
    }
    if (!is.list(value)) {
        return(invisible(TRUE))
    }
    for (index in seq_along(value)) {
        artifactDescriptorAssertDataOnly(
            value[[index]],
            paste0(owner, "[[", index, "]]")
        )
    }
    invisible(TRUE)
}

artifactDescriptorScalarStrings <- function(value, owner, allow_empty = FALSE) {
    valid <- is.character(value) && !anyNA(value) &&
        all(nzchar(value)) && anyDuplicated(value) == 0L
    if (!isTRUE(allow_empty) && length(value) == 0L) valid <- FALSE
    if (!isTRUE(valid)) {
        artifactDescriptorAbort(
            sprintf("artifact descriptor '%s' must contain unique strings", owner),
            "multischolar_invalid_artifact_descriptor"
        )
    }
    unname(value)
}

artifactDescriptorPositiveInteger <- function(value, owner, maximum = 1000000L) {
    valid <- length(value) == 1L && is.numeric(value) && !is.na(value) &&
        is.finite(value) && value >= 1L && value == floor(value) && value <= maximum
    if (!isTRUE(valid)) {
        artifactDescriptorAbort(
            sprintf("artifact descriptor '%s' must be a bounded positive integer", owner),
            "multischolar_invalid_artifact_descriptor"
        )
    }
    as.integer(value)
}

artifactDescriptorValidateCodec <- function(codec, key) {
    required <- c(
        "codec_id", "codec_version", "class_name",
        "payload_schema_id", "payload_schema_version"
    )
    if (!is.list(codec) || !identical(names(codec), required) ||
        !all(vapply(
            codec[c("codec_id", "class_name", "payload_schema_id")],
            workflowCapabilityScalarString,
            logical(1)
        )) || !identical(codec$codec_id, key)) {
        artifactDescriptorAbort(
            sprintf("artifact codec declaration '%s' is malformed", key),
            "multischolar_invalid_artifact_descriptor"
        )
    }
    codec$codec_version <- artifactDescriptorPositiveInteger(
        codec$codec_version, paste0(key, ".codec_version")
    )
    codec$payload_schema_version <- artifactDescriptorPositiveInteger(
        codec$payload_schema_version, paste0(key, ".payload_schema_version")
    )
    codec
}

artifactDescriptorValidateQuery <- function(query, key) {
    required <- c(
        "operation_id", "state_role", "projections", "filters", "order_by",
        "max_rows", "max_bytes"
    )
    if (!is.list(query) || !identical(names(query), required) ||
        !identical(query$operation_id, key) ||
        !workflowCapabilityScalarString(query$state_role)) {
        artifactDescriptorAbort(
            sprintf("artifact query declaration '%s' is malformed", key),
            "multischolar_invalid_artifact_descriptor"
        )
    }
    query$projections <- artifactDescriptorScalarStrings(
        query$projections, paste0(key, ".projections")
    )
    query$order_by <- artifactDescriptorScalarStrings(
        query$order_by, paste0(key, ".order_by")
    )
    if (!all(query$order_by %in% query$projections)) {
        artifactDescriptorAbort(
            sprintf("artifact query '%s' orders by an undeclared projection", key),
            "multischolar_invalid_artifact_descriptor"
        )
    }
    query$filters <- artifactDescriptorNamedList(
        query$filters, paste0(key, ".filters"), allow_empty = TRUE
    )
    allowed_operators <- c(
        "equal", "in", "gt", "gte", "lt", "lte", "between",
        "is_missing", "contains"
    )
    filter_ids <- names(query$filters)
    query$filters <- Map(function(filter, filter_id) {
        if (!is.list(filter) ||
            !identical(names(filter), c("column", "type", "operators")) ||
            !workflowCapabilityScalarString(filter$column) ||
            !filter$column %in% query$projections ||
            !workflowCapabilityScalarString(filter$type) ||
            !filter$type %in% c("character", "integer", "double", "logical") ||
            !all(filter$operators %in% allowed_operators)) {
            artifactDescriptorAbort(
                sprintf("artifact query filter '%s' is malformed", filter_id),
                "multischolar_invalid_artifact_descriptor"
            )
        }
        filter$operators <- artifactDescriptorScalarStrings(
            filter$operators, paste0(filter_id, ".operators")
        )
        filter
    }, query$filters, filter_ids)
    names(query$filters) <- filter_ids
    query$max_rows <- artifactDescriptorPositiveInteger(query$max_rows, "max_rows")
    query$max_bytes <- artifactDescriptorPositiveInteger(
        query$max_bytes, "max_bytes", maximum = 2^31 - 1L
    )
    query
}

artifactDescriptorValidateOwnedBlock <- function(value, owner_id, fields, label) {
    if (!is.list(value) || !identical(names(value), c("owner_id", fields)) ||
        !identical(value$owner_id, owner_id)) {
        artifactDescriptorAbort(
            sprintf("artifact descriptor %s is missing tuple-owned evidence", label),
            "multischolar_artifact_descriptor_evidence_mismatch"
        )
    }
    value
}

artifactDescriptorDigest <- function(descriptor) {
    candidate <- descriptor
    candidate$descriptor_digest <- NULL
    artifactSemanticDigest(candidate)
}

validateArtifactWorkflowDescriptor <- function(descriptor) {
    required <- c(
        "schema", "schema_version", "descriptor_id", "descriptor_version",
        "identity", "stages", "codecs", "queries", "fixtures",
        "scientific_oracle", "compatibility_products", "evidence", "migration",
        "rollback", "certification", "descriptor_digest"
    )
    artifactDescriptorAssertDataOnly(descriptor)
    version <- if (is.list(descriptor)) descriptor$schema_version else NA_integer_
    if (length(version) == 1L && is.numeric(version) && !is.na(version) &&
        version > ARTIFACT_DESCRIPTOR_SCHEMA_VERSION) {
        artifactDescriptorAbort(
            "artifact workflow descriptor uses a future schema version",
            "multischolar_future_artifact_descriptor_version"
        )
    }
    if (!is.list(descriptor) || !identical(names(descriptor), required) ||
        !identical(descriptor$schema, ARTIFACT_DESCRIPTOR_SCHEMA) ||
        !identical(as.integer(version), ARTIFACT_DESCRIPTOR_SCHEMA_VERSION) ||
        !workflowCapabilityScalarString(descriptor$descriptor_id)) {
        artifactDescriptorAbort(
            "artifact workflow descriptor is incomplete or malformed",
            "multischolar_invalid_artifact_descriptor"
        )
    }
    artifactDescriptorSemver(descriptor$descriptor_version, "descriptor_version")
    if (!identical(names(descriptor$identity), .WORKFLOW_CAPABILITY_KEY_FIELDS) ||
        !all(vapply(descriptor$identity, workflowCapabilityScalarString, logical(1)))) {
        artifactDescriptorAbort(
            "artifact workflow descriptor identity is not an exact workflow tuple",
            "multischolar_invalid_artifact_descriptor"
        )
    }
    stages <- artifactDescriptorNamedList(descriptor$stages, "stages")
    codecs <- artifactDescriptorNamedList(descriptor$codecs, "codecs")
    queries <- artifactDescriptorNamedList(descriptor$queries, "queries")
    codecs <- Map(artifactDescriptorValidateCodec, codecs, names(codecs))
    names(codecs) <- names(descriptor$codecs)
    queries <- Map(artifactDescriptorValidateQuery, queries, names(queries))
    names(queries) <- names(descriptor$queries)
    stage_required <- c(
        "stage_id", "state_roles", "codec_ids", "query_operation_ids",
        "maximum_rollout"
    )
    stages <- Map(function(stage, key) {
        if (!is.list(stage) || !identical(names(stage), stage_required) ||
            !identical(stage$stage_id, key) ||
            !workflowCapabilityScalarString(stage$maximum_rollout) ||
            !stage$maximum_rollout %in% .WORKFLOW_ARTIFACT_ROLLOUTS) {
            artifactDescriptorAbort(
                sprintf("artifact stage declaration '%s' is malformed", key),
                "multischolar_invalid_artifact_descriptor"
            )
        }
        stage$state_roles <- artifactDescriptorScalarStrings(
            stage$state_roles, paste0(key, ".state_roles")
        )
        stage$codec_ids <- artifactDescriptorScalarStrings(
            stage$codec_ids, paste0(key, ".codec_ids")
        )
        stage$query_operation_ids <- artifactDescriptorScalarStrings(
            stage$query_operation_ids,
            paste0(key, ".query_operation_ids"),
            allow_empty = TRUE
        )
        if (!all(stage$codec_ids %in% names(codecs)) ||
            !all(stage$query_operation_ids %in% names(queries))) {
            artifactDescriptorAbort(
                sprintf("artifact stage '%s' references an unknown codec or query", key),
                "multischolar_invalid_artifact_descriptor"
            )
        }
        query_roles <- vapply(
            queries[stage$query_operation_ids], `[[`, character(1), "state_role"
        )
        if (!all(query_roles %in% stage$state_roles)) {
            artifactDescriptorAbort(
                sprintf("artifact stage '%s' queries an undeclared state role", key),
                "multischolar_invalid_artifact_descriptor"
            )
        }
        stage
    }, stages, names(stages))
    names(stages) <- names(descriptor$stages)
    owner <- descriptor$descriptor_id
    fixtures <- artifactDescriptorValidateOwnedBlock(
        descriptor$fixtures, owner, "fixture_ids", "fixtures"
    )
    fixtures$fixture_ids <- artifactDescriptorScalarStrings(
        fixtures$fixture_ids, "fixture_ids"
    )
    oracle <- artifactDescriptorValidateOwnedBlock(
        descriptor$scientific_oracle,
        owner,
        c("oracle_id", "oracle_version", "tolerances"),
        "scientific oracle"
    )
    artifactDescriptorScalarStrings(oracle$oracle_id, "oracle_id")
    artifactDescriptorSemver(oracle$oracle_version, "oracle_version")
    if (!is.numeric(oracle$tolerances) || is.null(names(oracle$tolerances)) ||
        any(!is.finite(oracle$tolerances)) || any(oracle$tolerances < 0)) {
        artifactDescriptorAbort(
            "artifact scientific oracle requires owned non-negative tolerances",
            "multischolar_invalid_artifact_descriptor"
        )
    }
    products <- artifactDescriptorValidateOwnedBlock(
        descriptor$compatibility_products, owner, "product_ids", "compatibility products"
    )
    products$product_ids <- artifactDescriptorScalarStrings(
        products$product_ids, "product_ids"
    )
    evidence_fields <- c(
        "inventory_ids", "codec_ids", "stage_ids", "lifecycle_ids",
        "performance_thresholds"
    )
    evidence <- artifactDescriptorValidateOwnedBlock(
        descriptor$evidence, owner, evidence_fields, "evidence"
    )
    for (field in evidence_fields[1:4]) {
        evidence[[field]] <- artifactDescriptorScalarStrings(evidence[[field]], field)
    }
    if (!setequal(evidence$codec_ids, names(codecs)) ||
        !setequal(evidence$stage_ids, names(stages)) ||
        !is.numeric(evidence$performance_thresholds) ||
        is.null(names(evidence$performance_thresholds)) ||
        any(!is.finite(evidence$performance_thresholds)) ||
        any(evidence$performance_thresholds <= 0)) {
        artifactDescriptorAbort(
            "artifact descriptor evidence does not cover its own codecs and stages",
            "multischolar_artifact_descriptor_evidence_mismatch"
        )
    }
    migration <- artifactDescriptorValidateOwnedBlock(
        descriptor$migration,
        owner,
        c("strategy_id", "from_backend", "to_backend"),
        "migration"
    )
    rollback <- artifactDescriptorValidateOwnedBlock(
        descriptor$rollback,
        owner,
        c("strategy_id", "target_backend"),
        "rollback"
    )
    certification <- artifactDescriptorValidateOwnedBlock(
        descriptor$certification,
        owner,
        c("status", "auto_eligible"),
        "certification"
    )
    scalar_fields <- c(
        migration[c("strategy_id", "from_backend", "to_backend")],
        rollback[c("strategy_id", "target_backend")]
    )
    if (!all(vapply(scalar_fields, workflowCapabilityScalarString, logical(1))) ||
        !workflowCapabilityScalarString(certification$status) ||
        !certification$status %in% .ARTIFACT_CERTIFICATION_LEVELS ||
        !is.logical(certification$auto_eligible) ||
        length(certification$auto_eligible) != 1L || is.na(certification$auto_eligible)) {
        artifactDescriptorAbort(
            "artifact descriptor migration, rollback, or certification is malformed",
            "multischolar_invalid_artifact_descriptor"
        )
    }
    if (isTRUE(certification$auto_eligible) &&
        identical(certification$status, "uncertified")) {
        artifactDescriptorAbort(
            "an uncertified artifact descriptor cannot be auto eligible",
            "multischolar_artifact_not_certified"
        )
    }
    if (!identical(descriptor$descriptor_digest, artifactDescriptorDigest(descriptor))) {
        artifactDescriptorAbort(
            "artifact workflow descriptor digest does not match its contents",
            "multischolar_artifact_descriptor_digest_mismatch"
        )
    }
    descriptor$stages <- stages
    descriptor$codecs <- codecs
    descriptor$queries <- queries
    descriptor$fixtures <- fixtures
    descriptor$compatibility_products <- products
    descriptor$evidence <- evidence
    descriptor
}

newArtifactWorkflowDescriptor <- function(
    descriptor_id,
    descriptor_version,
    identity,
    stages,
    codecs,
    queries,
    fixtures,
    scientific_oracle,
    compatibility_products,
    evidence,
    migration,
    rollback,
    certification
) {
    descriptor <- list(
        schema = ARTIFACT_DESCRIPTOR_SCHEMA,
        schema_version = ARTIFACT_DESCRIPTOR_SCHEMA_VERSION,
        descriptor_id = descriptor_id,
        descriptor_version = descriptor_version,
        identity = identity,
        stages = stages,
        codecs = codecs,
        queries = queries,
        fixtures = fixtures,
        scientific_oracle = scientific_oracle,
        compatibility_products = compatibility_products,
        evidence = evidence,
        migration = migration,
        rollback = rollback,
        certification = certification,
        descriptor_digest = NULL
    )
    artifactDescriptorAssertDataOnly(descriptor)
    descriptor$descriptor_digest <- artifactDescriptorDigest(descriptor)
    structure(
        validateArtifactWorkflowDescriptor(descriptor),
        class = c("MultiScholaRArtifactWorkflowDescriptor", "list")
    )
}

newArtifactDescriptorCatalogue <- function(descriptors = list()) {
    if (!is.list(descriptors)) {
        artifactDescriptorAbort(
            "artifact descriptor catalogue input must be a list",
            "multischolar_invalid_artifact_descriptor_catalogue"
        )
    }
    descriptors <- lapply(descriptors, validateArtifactWorkflowDescriptor)
    ids <- vapply(descriptors, `[[`, character(1), "descriptor_id")
    keys <- vapply(
        descriptors,
        function(value) workflowCapabilityKey(value$identity),
        character(1)
    )
    if (anyDuplicated(ids) > 0L || anyDuplicated(keys) > 0L) {
        artifactDescriptorAbort(
            "artifact descriptor catalogue contains duplicate IDs or workflow tuples",
            "multischolar_duplicate_artifact_descriptor"
        )
    }
    catalogue <- new.env(parent = emptyenv())
    catalogue$schema <- ARTIFACT_DESCRIPTOR_CATALOGUE_SCHEMA
    catalogue$schema_version <- ARTIFACT_DESCRIPTOR_CATALOGUE_VERSION
    catalogue$descriptors <- stats::setNames(descriptors, ids)
    catalogue$key_index <- stats::setNames(ids, keys)
    class(catalogue) <- c("MultiScholaRArtifactDescriptorCatalogue", "environment")
    lockEnvironment(catalogue, bindings = TRUE)
    catalogue
}

validateArtifactDescriptorCatalogue <- function(catalogue) {
    valid <- inherits(catalogue, "MultiScholaRArtifactDescriptorCatalogue") &&
        is.environment(catalogue) && environmentIsLocked(catalogue) &&
        identical(catalogue$schema, ARTIFACT_DESCRIPTOR_CATALOGUE_SCHEMA) &&
        identical(catalogue$schema_version, ARTIFACT_DESCRIPTOR_CATALOGUE_VERSION)
    if (!isTRUE(valid)) {
        artifactDescriptorAbort(
            "artifact descriptor catalogue is invalid or mutable",
            "multischolar_invalid_artifact_descriptor_catalogue"
        )
    }
    catalogue
}

artifactDescriptorCatalogueValues <- function(catalogue) {
    catalogue <- validateArtifactDescriptorCatalogue(catalogue)
    lapply(catalogue$descriptors, identity)
}

findArtifactWorkflowDescriptor <- function(identity, catalogue) {
    catalogue <- validateArtifactDescriptorCatalogue(catalogue)
    key <- workflowCapabilityKey(identity)
    descriptor_id <- unname(catalogue$key_index[key])
    if (length(descriptor_id) == 0L || is.na(descriptor_id)) return(NULL)
    catalogue$descriptors[[descriptor_id]]
}

artifactDescriptorMaximumRollout <- function(descriptor) {
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    if (identical(descriptor$certification$status, "uncertified")) return(NULL)
    stage_ranks <- match(
        vapply(descriptor$stages, `[[`, character(1), "maximum_rollout"),
        .WORKFLOW_ARTIFACT_ROLLOUTS
    )
    certification_rank <- match(
        descriptor$certification$status,
        .WORKFLOW_ARTIFACT_ROLLOUTS
    )
    .WORKFLOW_ARTIFACT_ROLLOUTS[[min(c(stage_ranks, certification_rank))]]
}

artifactDescriptorCapabilities <- function(catalogue) {
    lapply(artifactDescriptorCatalogueValues(catalogue), function(descriptor) {
        maximum <- artifactDescriptorMaximumRollout(descriptor)
        newWorkflowCapability(
            capability_id = descriptor$descriptor_id,
            identity = descriptor$identity,
            artifact_eligible = !is.null(maximum),
            auto_eligible = isTRUE(descriptor$certification$auto_eligible),
            maximum_artifact_rollout = maximum,
            capability_version = descriptor$descriptor_version
        )
    })
}

artifactDiaWorkflowCodecDeclarations <- function() {
    list(
        "multischolar.s4.peptide_quantitative_data.diann" = list(
            codec_id = "multischolar.s4.peptide_quantitative_data.diann",
            codec_version = 1L,
            class_name = "PeptideQuantitativeData",
            payload_schema_id = "multischolar.rectangular",
            payload_schema_version = 1L
        ),
        "multischolar.s4.protein_quantitative_data.diann" = list(
            codec_id = "multischolar.s4.protein_quantitative_data.diann",
            codec_version = 1L,
            class_name = "ProteinQuantitativeData",
            payload_schema_id = "multischolar.rectangular",
            payload_schema_version = 1L
        )
    )
}

artifactDiaWorkflowDescriptor <- function() {
    owner <- "proteomics.diann.peptide.dia.v1"
    codecs <- artifactDiaWorkflowCodecDeclarations()
    query <- "proteomics.diann.import.preview.v1"
    newArtifactWorkflowDescriptor(
        descriptor_id = owner,
        descriptor_version = "1.0.0",
        identity = workflowCapabilityIdentity(
            "proteomics", "proteomics.gui", "DIA", "prot_dia",
            "diann", "peptide", "dia"
        ),
        stages = list(
            import = list(
                stage_id = "import",
                state_roles = c("canonical_data"),
                codec_ids = names(codecs),
                query_operation_ids = query,
                maximum_rollout = "dual_write"
            ),
            design = list(
                stage_id = "design",
                state_roles = c(
                    "cleaned_data", "design_matrix", "contrasts", "args",
                    "annotations", "sequences", "raw_data_s4"
                ),
                codec_ids = names(codecs),
                query_operation_ids = character(),
                maximum_rollout = "dual_write"
            )
        ),
        codecs = codecs,
        queries = stats::setNames(list(list(
            operation_id = query,
            state_role = "canonical_data",
            projections = c(
                "Run", "Protein.Group", "Stripped.Sequence", "Precursor.Id",
                "Precursor.Normalised", "Q.Value"
            ),
            filters = list(
                run = list(
                    column = "Run",
                    type = "character",
                    operators = c("equal", "in")
                ),
                protein_group = list(
                    column = "Protein.Group",
                    type = "character",
                    operators = c("equal", "in")
                )
            ),
            order_by = c(
                "Run", "Protein.Group", "Stripped.Sequence", "Precursor.Id"
            ),
            max_rows = 10000L,
            max_bytes = 64L * 1024L * 1024L
        )), query),
        fixtures = list(
            owner_id = owner,
            fixture_ids = c(
                "tests/testdata/e2e/prot_dia/report.tsv",
                "tests/testdata/e2e/prot_dia_limpa/seed_report.tsv"
            )
        ),
        scientific_oracle = list(
            owner_id = owner,
            oracle_id = "tests/testdata/omics-parity/dia-memory-oracle.json",
            oracle_version = "1.0.0",
            tolerances = c(quantity_absolute = 1e-10, quantity_relative = 1e-10)
        ),
        compatibility_products = list(
            owner_id = owner,
            product_ids = c(
                "data_cln.tab", "design_matrix.tab", "contrast_strings.tab",
                "config.ini", "aa_seq_tbl_final.RDS", "uniprot_dat_cln.RDS"
            )
        ),
        evidence = list(
            owner_id = owner,
            inventory_ids = "tests/testdata/omics-capabilities.json",
            codec_ids = names(codecs),
            stage_ids = c("import", "design"),
            lifecycle_ids = "OMICS-ART-058",
            performance_thresholds = c(
                runtime_ratio = 1.25,
                committed_disk_ratio = 1.35,
                retained_rss_reduction = 0.40
            )
        ),
        migration = list(
            owner_id = owner,
            strategy_id = "explicit_dual_write_canary",
            from_backend = "memory",
            to_backend = "artifact"
        ),
        rollback = list(
            owner_id = owner,
            strategy_id = "force_memory_ignore_canary_generations",
            target_backend = "memory"
        ),
        certification = list(
            owner_id = owner,
            status = "dual_write",
            auto_eligible = FALSE
        )
    )
}

.ARTIFACT_WORKFLOW_DESCRIPTOR_CATALOGUE <- newArtifactDescriptorCatalogue(
    list(artifactDiaWorkflowDescriptor())
)

artifactWorkflowDescriptorCatalogue <- function() {
    .ARTIFACT_WORKFLOW_DESCRIPTOR_CATALOGUE
}
