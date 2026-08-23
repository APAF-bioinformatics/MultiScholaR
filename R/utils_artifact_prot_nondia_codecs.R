# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.artifactProteomicsNonDiaCodecVersion <- 1L
.artifactProteomicsNonDiaCapabilityVersion <- "1.0.0"
.artifactProteomicsNonDiaStateRoles <- c(
    "protein_s4_initial",
    "normalized",
    "ruv_corrected",
    "correlation_filtered"
)

newArtifactProteomicsNonDiaCodecRole <- function(
    capability_id,
    codec_id,
    workflow_type
) {
    list(
        capability_id = capability_id,
        capability_version = .artifactProteomicsNonDiaCapabilityVersion,
        data_level = "protein",
        s4_role = "protein_quantitative_state",
        state_roles = .artifactProteomicsNonDiaStateRoles,
        class_name = "ProteinQuantitativeData",
        workflow_type = workflow_type,
        codec_id = codec_id,
        codec_version = .artifactProteomicsNonDiaCodecVersion
    )
}

artifactProteomicsNonDiaCodecRoles <- function() {
    roles <- list(
        newArtifactProteomicsNonDiaCodecRole(
            "proteomics.maxquant.protein.lfq.v1",
            "multischolar.s4.protein_quantitative_data.maxquant.protein.lfq",
            "LFQ"
        ),
        newArtifactProteomicsNonDiaCodecRole(
            "proteomics.fragpipe.protein.lfq.v1",
            "multischolar.s4.protein_quantitative_data.fragpipe.protein.lfq",
            "LFQ"
        ),
        newArtifactProteomicsNonDiaCodecRole(
            "proteomics.pd_tmt.protein.tmt.v1",
            "multischolar.s4.protein_quantitative_data.pd_tmt.protein.tmt",
            "TMT"
        )
    )
    stats::setNames(
        roles,
        vapply(roles, `[[`, character(1), "capability_id")
    )
}

artifactProteomicsNonDiaCodecDeclarations <- function() {
    roles <- artifactProteomicsNonDiaCodecRoles()
    declarations <- lapply(roles, \(role) {
        list(
            codec_id = role$codec_id,
            codec_version = role$codec_version,
            class_name = role$class_name,
            payload_schema_id = "multischolar.rectangular",
            payload_schema_version = 1L
        )
    })
    stats::setNames(
        declarations,
        vapply(declarations, `[[`, character(1), "codec_id")
    )
}

artifactProteomicsNonDiaRoleRemediation <- function(
    capability_id,
    capability_version,
    roles
) {
    role <- if (workflowCapabilityScalarString(capability_id)) {
        roles[[capability_id]]
    } else {
        NULL
    }
    if (is.null(role)) {
        return(paste(
            "Keep this tuple memory-backed until a reviewed fixture, scientific oracle,",
            "and exact codec role are certified."
        ))
    }
    if (!identical(capability_version, role$capability_version)) {
        return(sprintf(
            "Use certified capability version '%s' or add a reviewed codec version.",
            role$capability_version
        ))
    }
    sprintf(
        "Use one of the certified state roles for this tuple: %s.",
        paste(role$state_roles, collapse = ", ")
    )
}

artifactProteomicsNonDiaCodecRole <- function(
    capability_id,
    state_role,
    capability_version = .artifactProteomicsNonDiaCapabilityVersion
) {
    roles <- artifactProteomicsNonDiaCodecRoles()
    valid_capability <- workflowCapabilityScalarString(capability_id)
    role <- if (valid_capability) roles[[capability_id]] else NULL
    supported <- valid_capability &&
        workflowCapabilityScalarString(capability_version) &&
        workflowCapabilityScalarString(state_role) &&
        !is.null(role) &&
        identical(capability_version, role$capability_version) &&
        state_role %in% role$state_roles
    if (!isTRUE(supported)) {
        artifactCodecAbort(
            "non-DIA proteomics capability or state role has no certified codec",
            "multischolar_unsupported_proteomics_codec_role",
            capability_id = capability_id,
            capability_version = capability_version,
            state_role = state_role,
            remediation = artifactProteomicsNonDiaRoleRemediation(
                capability_id,
                capability_version,
                roles
            )
        )
    }
    role
}

artifactProteomicsNonDiaShapeAbort <- function(
    role,
    state_role,
    slot_name,
    reason,
    remediation
) {
    artifactCodecAbort(
        sprintf(
            "non-DIA proteomics lane '%s' state '%s' has unsupported slot '%s': %s",
            role$capability_id,
            state_role,
            slot_name,
            reason
        ),
        "multischolar_proteomics_codec_shape_mismatch",
        capability_id = role$capability_id,
        capability_version = role$capability_version,
        state_role = state_role,
        slot_name = slot_name,
        remediation = remediation
    )
}

artifactProteomicsNonDiaScalarSlot <- function(value, slot_name, role, state_role) {
    slot_value <- methods::slot(value, slot_name)
    if (!workflowCapabilityScalarString(slot_value)) {
        artifactProteomicsNonDiaShapeAbort(
            role,
            state_role,
            slot_name,
            "the slot must contain one non-missing name",
            sprintf("Rebuild the state with a scalar '%s' slot.", slot_name)
        )
    }
    slot_value
}

artifactProteomicsNonDiaOpaqueArgPath <- function(value, owner = "args") {
    if (is.raw(value)) return(owner)
    if (!is.list(value)) return(NULL)
    value_names <- names(value)
    for (index in seq_along(value)) {
        label <- if (!is.null(value_names) && nzchar(value_names[[index]])) {
            value_names[[index]]
        } else {
            as.character(index)
        }
        found <- artifactProteomicsNonDiaOpaqueArgPath(
            value[[index]],
            paste0(owner, "[[", label, "]]")
        )
        if (!is.null(found)) return(found)
    }
    NULL
}

artifactProteomicsNonDiaValidateClass <- function(value, role, state_role) {
    class_name <- class(value)
    if (!isS4(value) || length(class_name) != 1L ||
        !identical(class_name[[1L]], role$class_name)) {
        artifactProteomicsNonDiaShapeAbort(
            role,
            state_role,
            ".class",
            sprintf("expected exactly '%s'", role$class_name),
            "Rebuild the state through the supported protein-level design workflow."
        )
    }
    validity <- tryCatch(
        methods::validObject(value, test = TRUE),
        error = \(error) conditionMessage(error)
    )
    if (!identical(validity, TRUE)) {
        artifactProteomicsNonDiaShapeAbort(
            role,
            state_role,
            ".validObject",
            paste(validity, collapse = "; "),
            "Repair or rebuild the S4 object in memory before artifact migration."
        )
    }
    invisible(value)
}

artifactProteomicsNonDiaScalarSlots <- function(value, role, state_role) {
    slot_names <- c(
        "protein_id_column", "sample_id", "group_id",
        "technical_replicate_id"
    )
    stats::setNames(lapply(
        slot_names,
        \(slot_name) artifactProteomicsNonDiaScalarSlot(
            value,
            slot_name,
            role,
            state_role
        )
    ), slot_names)
}

artifactProteomicsNonDiaValidateDesign <- function(
    value,
    scalar_slots,
    role,
    state_role
) {
    design <- value@design_matrix
    required <- unlist(scalar_slots[c(
        "sample_id", "group_id", "technical_replicate_id"
    )], use.names = FALSE)
    if (!all(required %in% names(design))) {
        artifactProteomicsNonDiaShapeAbort(
            role,
            state_role,
            "design_matrix",
            "required sample, group, or technical-replicate columns are missing",
            "Rebuild the design checkpoint with all declared design columns."
        )
    }
    sample_id <- scalar_slots$sample_id
    samples <- design[[sample_id]]
    if (!is.character(samples) || anyNA(samples) ||
        any(!nzchar(samples)) || anyDuplicated(samples) > 0L) {
        artifactProteomicsNonDiaShapeAbort(
            role,
            state_role,
            paste0("design_matrix$", sample_id),
            "sample identities must be unique, non-missing character values",
            "Normalize sample identities before artifact migration."
        )
    }
    metadata_columns <- unlist(scalar_slots[c(
        "group_id", "technical_replicate_id"
    )], use.names = FALSE)
    invalid_metadata <- vapply(design[metadata_columns], \(column) {
        values <- as.character(column)
        anyNA(values) || any(!nzchar(values))
    }, logical(1))
    if (any(invalid_metadata)) {
        artifactProteomicsNonDiaShapeAbort(
            role,
            state_role,
            paste0("design_matrix$", metadata_columns[invalid_metadata][[1L]]),
            "group and technical-replicate identities cannot be missing or empty",
            "Complete the design metadata before artifact migration."
        )
    }
    samples
}

artifactProteomicsNonDiaValidateQuant <- function(
    value,
    scalar_slots,
    samples,
    role,
    state_role
) {
    quant_table <- value@protein_quant_table
    protein_id <- scalar_slots$protein_id_column
    proteins <- quant_table[[protein_id]]
    if (!is.character(proteins) || anyNA(proteins) ||
        any(!nzchar(proteins)) || anyDuplicated(proteins) > 0L) {
        artifactProteomicsNonDiaShapeAbort(
            role,
            state_role,
            paste0("protein_quant_table$", protein_id),
            "protein identities must be unique, non-missing character values",
            "Resolve protein identity duplicates before artifact migration."
        )
    }
    sample_columns <- setdiff(names(quant_table), protein_id)
    numeric_samples <- vapply(quant_table[sample_columns], is.numeric, logical(1))
    if (!identical(sample_columns, samples) || !all(numeric_samples)) {
        artifactProteomicsNonDiaShapeAbort(
            role,
            state_role,
            "protein_quant_table",
            "ordered numeric sample columns must match the design sample order",
            "Rebuild the protein table from the declared design and column mapping."
        )
    }
    invisible(value)
}

artifactProteomicsNonDiaValidateMetadata <- function(
    value,
    scalar_slots,
    role,
    state_role
) {
    protein_id <- scalar_slots$protein_id_column
    if (!protein_id %in% names(value@protein_id_table)) {
        artifactProteomicsNonDiaShapeAbort(
            role,
            state_role,
            "protein_id_table",
            "the active protein identity column is missing",
            "Rebuild protein identity/provenance metadata before migration."
        )
    }
    provenance_ids <- value@protein_id_table[[protein_id]]
    quant_ids <- value@protein_quant_table[[protein_id]]
    if (!is.character(provenance_ids) || anyNA(provenance_ids) ||
        any(!nzchar(provenance_ids)) || !all(quant_ids %in% provenance_ids)) {
        artifactProteomicsNonDiaShapeAbort(
            role,
            state_role,
            paste0("protein_id_table$", protein_id),
            "protein provenance must cover every non-missing quantitative identity",
            "Rebuild protein identity/provenance metadata before migration."
        )
    }
    observed_workflow <- value@args$globalParameters$workflow_type
    if (!identical(observed_workflow, role$workflow_type)) {
        artifactProteomicsNonDiaShapeAbort(
            role,
            state_role,
            "args$globalParameters$workflow_type",
            sprintf("expected workflow type '%s'", role$workflow_type),
            "Rebuild the design checkpoint under the selected exact workflow tuple."
        )
    }
    opaque_arg <- artifactProteomicsNonDiaOpaqueArgPath(value@args)
    if (!is.null(opaque_arg)) {
        artifactProteomicsNonDiaShapeAbort(
            role,
            state_role,
            opaque_arg,
            "opaque raw/RDS bytes are not a certified state component",
            "Persist the value under an explicit immutable artifact role and reference it."
        )
    }
    invisible(value)
}

artifactValidateProteomicsNonDiaProteinState <- function(value, role, state_role) {
    artifactProteomicsNonDiaValidateClass(value, role, state_role)
    scalar_slots <- artifactProteomicsNonDiaScalarSlots(value, role, state_role)
    samples <- artifactProteomicsNonDiaValidateDesign(
        value,
        scalar_slots,
        role,
        state_role
    )
    artifactProteomicsNonDiaValidateQuant(
        value,
        scalar_slots,
        samples,
        role,
        state_role
    )
    artifactProteomicsNonDiaValidateMetadata(
        value,
        scalar_slots,
        role,
        state_role
    )
    invisible(value)
}

dehydrateProteomicsNonDiaS4Artifact <- function(
    value,
    capability_id,
    state_role,
    capability_version = .artifactProteomicsNonDiaCapabilityVersion,
    inline_limit_bytes = .artifactInlineLimitBytes
) {
    role <- artifactProteomicsNonDiaCodecRole(
        capability_id,
        state_role,
        capability_version
    )
    artifactValidateProteomicsNonDiaProteinState(value, role, state_role)
    codec <- artifactProteomicsNonDiaCodecDeclarations()[[role$codec_id]]
    dehydrateExactS4Artifact(value, codec, inline_limit_bytes)
}

validateProteomicsNonDiaS4Bundle <- function(
    bundle,
    capability_id,
    state_role,
    capability_version = .artifactProteomicsNonDiaCapabilityVersion
) {
    role <- artifactProteomicsNonDiaCodecRole(
        capability_id,
        state_role,
        capability_version
    )
    codec <- artifactProteomicsNonDiaCodecDeclarations()[[role$codec_id]]
    validateExactS4Bundle(bundle, codec)
}

hydrateProteomicsNonDiaS4Artifact <- function(
    bundle,
    capability_id,
    state_role,
    capability_version = .artifactProteomicsNonDiaCapabilityVersion
) {
    bundle <- validateProteomicsNonDiaS4Bundle(
        bundle,
        capability_id,
        state_role,
        capability_version
    )
    role <- artifactProteomicsNonDiaCodecRole(
        capability_id,
        state_role,
        capability_version
    )
    codec <- artifactProteomicsNonDiaCodecDeclarations()[[role$codec_id]]
    value <- hydrateExactS4Artifact(bundle, codec)
    artifactValidateProteomicsNonDiaProteinState(value, role, state_role)
    value
}
