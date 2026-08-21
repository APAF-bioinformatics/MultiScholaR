# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

protDiaLegacySessionRestoreMode <- function() {
    match.arg(
        getOption("multischolar.prot_dia.legacy_session_restore", "enabled"),
        c("enabled", "disabled")
    )
}

loadProtDaLegacySessionBundle <- function(
    source_dir,
    read_rds_fn = readRDS
) {
    if (!identical(protDiaLegacySessionRestoreMode(), "enabled")) {
        workflowSessionAbort(
            "legacy DIA-NN session restore is disabled",
            "multischolar_prot_dia_legacy_session_disabled"
        )
    }
    path <- file.path(source_dir, "filtered_session_data_latest.rds")
    if (!file.exists(path)) {
        workflowSessionAbort(
            paste(
                "No exported session data found. Please export session data",
                "from the Normalization tab first."
            ),
            "multischolar_missing_prot_dia_legacy_session"
        )
    }
    session_data <- read_rds_fn(path)
    validateProtDaFilteredSession(session_data, source_dir = source_dir)
    annotation_path <- file.path(source_dir, "uniprot_dat_cln.RDS")
    annotations <- if (file.exists(annotation_path)) {
        read_rds_fn(annotation_path)
    } else {
        NULL
    }
    list(
        format = "legacy",
        session_data = session_data,
        state_manager = NULL,
        workflow_context = NULL,
        uniprot_dat_cln = annotations,
        aa_seq_tbl_final = NULL,
        session_path = path,
        artifact_error = NULL
    )
}

resolveProtDaSessionBundle <- function(
    experiment_paths,
    resource_policy = NULL,
    restore_artifact_fn = restoreProtDiaSessionManifest,
    load_legacy_fn = loadProtDaLegacySessionBundle
) {
    source_dir <- experiment_paths$source_dir
    if (is.null(source_dir) || !dir.exists(source_dir)) {
        stop(
            "Could not find the source directory to load session data.",
            call. = FALSE
        )
    }
    artifact_path <- file.path(
        source_dir,
        "filtered_session_artifact_latest.json"
    )
    artifact_error <- NULL
    if (identical(protDiaSessionMode("restore"), "enabled") &&
        file.exists(artifact_path)) {
        bundle <- tryCatch(
            restore_artifact_fn(
                artifact_path,
                experiment_paths,
                resource_policy
            ),
            error = function(error) {
                artifact_error <<- error
                NULL
            }
        )
        if (!is.null(bundle)) {
            bundle$session_path <- artifact_path
            bundle$artifact_error <- NULL
            return(bundle)
        }
    }
    legacy <- tryCatch(
        load_legacy_fn(source_dir),
        error = function(error) {
            if (!is.null(artifact_error)) stop(artifact_error)
            stop(error)
        }
    )
    legacy$artifact_error <- artifact_error
    legacy
}

protDaSessionStateSnapshot <- function(workflow_data, bundle) {
    if (identical(bundle$format, "artifact")) {
        return(workflowStateLegacySnapshot(bundle$state_manager))
    }
    restoreWorkflowStateFromSession(
        workflow_data$state_manager,
        bundle$session_data
    )
}

protDaSessionFormula <- function(object) {
    if (is.null(object) || !isS4(object) ||
        !"args" %in% methods::slotNames(object)) {
        return(NULL)
    }
    parameters <- object@args$deAnalysisParameters
    formula <- parameters$formula_string
    if (workflowCapabilityScalarString(formula)) formula else NULL
}

protDaSessionContrasts <- function(contrasts_tbl) {
    if (is.null(contrasts_tbl)) return(NULL)
    if ("contrasts" %in% names(contrasts_tbl)) {
        return(contrasts_tbl$contrasts)
    }
    contrasts_tbl[, 1L]
}

protDaSessionFieldValues <- function(session_data) {
    fields <- c(
        "fasta_metadata", "accession_cleanup_results",
        "ruv_optimization_result", "qc_params", "protein_counts",
        "mixed_species_analysis"
    )
    session_data[intersect(fields, names(session_data))]
}

protDaSessionCaptureFields <- function(owner, fields) {
    values <- lapply(fields, function(name) {
        tryCatch(owner[[name]], error = function(error) NULL)
    })
    names(values) <- fields
    values
}

protDaSessionRestoreFields <- function(owner, values) {
    for (name in names(values)) owner[[name]] <- values[[name]]
    invisible(NULL)
}

protDaSessionCaptureGlobals <- function() {
    names <- c("contrasts_tbl", "config_list", "uniprot_dat_cln")
    values <- lapply(names, function(name) {
        if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
            return(list(exists = TRUE, value = get(
                name,
                envir = .GlobalEnv,
                inherits = FALSE
            )))
        }
        list(exists = FALSE, value = NULL)
    })
    stats::setNames(values, names)
}

protDaSessionRestoreGlobals <- function(values) {
    for (name in names(values)) {
        if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
            rm(list = name, envir = .GlobalEnv)
        }
        if (isTRUE(values[[name]]$exists)) {
            assign(name, values[[name]]$value, envir = .GlobalEnv)
        }
    }
    invisible(NULL)
}

protDaSessionCaptureManager <- function(manager) {
    list(
        manifest = tryCatch(
            workflowStateManifest(manager),
            error = function(error) NULL
        ),
        legacy = tryCatch(
            workflowStateLegacySnapshot(manager),
            error = function(error) NULL
        )
    )
}

protDaSessionRestoreManager <- function(manager, captured) {
    restorer <- workflowStateMember(manager, "restoreState")
    if (is.function(restorer) && !is.null(captured$manifest)) {
        restorer(captured$manifest)
        return(invisible(NULL))
    }
    if (is.list(captured$legacy)) {
        manager$states <- captured$legacy$r6_complete_states
        manager$state_history <- captured$legacy$r6_state_history
        manager$current_state <- captured$legacy$r6_current_state_name
    }
    invisible(NULL)
}

applyProtDaSessionBundle <- function(
    workflow_data,
    da_data,
    bundle,
    failure_injector = NULL
) {
    session_data <- bundle$session_data
    old_manager <- workflow_data$state_manager
    manager_before <- protDaSessionCaptureManager(old_manager)
    workflow_fields <- c(
        "workflow_context", "state_manager", "contrasts_tbl",
        "design_matrix", "config_list", "uniprot_dat_cln",
        "aa_seq_tbl_final", "artifact_session_manifest",
        names(protDaSessionFieldValues(session_data)),
        "tab_status", "state_update_trigger"
    )
    workflow_before <- protDaSessionCaptureFields(
        workflow_data,
        unique(workflow_fields)
    )
    da_before <- protDaSessionCaptureFields(
        da_data,
        c("current_s4_object", "contrasts_available", "formula_from_s4")
    )
    globals_before <- protDaSessionCaptureGlobals()
    result <- tryCatch({
        state_snapshot <- protDaSessionStateSnapshot(workflow_data, bundle)
        artifactStoreInvokeFailure(
            failure_injector,
            "after_session_state_prepare",
            list(format = bundle$format)
        )
        if (identical(bundle$format, "artifact")) {
            workflow_data$workflow_context <- bundle$workflow_context
            workflow_data$state_manager <- bundle$state_manager
            workflow_data$artifact_session_manifest <- bundle$manifest
        }
        workflow_data$contrasts_tbl <- session_data$contrasts_tbl
        workflow_data$design_matrix <- session_data$design_matrix
        workflow_data$config_list <- session_data$config_list
        workflow_data$uniprot_dat_cln <- bundle$uniprot_dat_cln
        if (!is.null(bundle$aa_seq_tbl_final)) {
            workflow_data$aa_seq_tbl_final <- bundle$aa_seq_tbl_final
        }
        for (name in names(protDaSessionFieldValues(session_data))) {
            workflow_data[[name]] <- session_data[[name]]
        }
        da_data$current_s4_object <- session_data$current_s4_object
        da_data$contrasts_available <- protDaSessionContrasts(
            session_data$contrasts_tbl
        )
        da_data$formula_from_s4 <- protDaSessionFormula(
            session_data$current_s4_object
        )
        if (!is.null(session_data$contrasts_tbl)) {
            assign(
                "contrasts_tbl",
                session_data$contrasts_tbl,
                envir = .GlobalEnv
            )
        }
        if (!is.null(session_data$config_list)) {
            assign("config_list", session_data$config_list, envir = .GlobalEnv)
        }
        if (!is.null(bundle$uniprot_dat_cln)) {
            assign(
                "uniprot_dat_cln",
                bundle$uniprot_dat_cln,
                envir = .GlobalEnv
            )
        }
        status <- workflow_data$tab_status
        status$normalization <- "complete"
        status$differential_expression <- "pending"
        workflow_data$tab_status <- status
        workflow_data$state_update_trigger <- Sys.time()
        artifactStoreInvokeFailure(
            failure_injector,
            "before_session_apply_commit",
            list(format = bundle$format)
        )
        list(
            state_snapshot = state_snapshot,
            formula = da_data$formula_from_s4,
            annotations = bundle$uniprot_dat_cln
        )
    }, error = function(error) {
        if (identical(bundle$format, "legacy")) {
            try(
                protDaSessionRestoreManager(old_manager, manager_before),
                silent = TRUE
            )
        }
        protDaSessionRestoreFields(workflow_data, workflow_before)
        protDaSessionRestoreFields(da_data, da_before)
        protDaSessionRestoreGlobals(globals_before)
        if (identical(bundle$format, "artifact") &&
            inherits(bundle$state_manager, "ArtifactWorkflowState")) {
            try(bundle$state_manager$close(), silent = TRUE)
        }
        stop(error)
    })
    if (identical(bundle$format, "artifact") &&
        inherits(old_manager, "ArtifactWorkflowState") &&
        !identical(old_manager, bundle$state_manager)) {
        try(old_manager$close(), silent = TRUE)
    }
    result
}

ensureProtDaCompatibilitySession <- function(bundle, source_dir) {
    if (!identical(bundle$format, "artifact")) return(NULL)
    path <- file.path(source_dir, "filtered_session_data_latest.rds")
    if (file.exists(path)) return(NULL)
    writeProtDiaCompatibilitySession(bundle, path)
}
