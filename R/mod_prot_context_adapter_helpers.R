# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

.PROT_CONTEXT_GLOBAL_OPTION <-
    "multischolar.prot.legacy_global_adapters"

.PROT_CONTEXT_DEPENDENCIES <- list(
    project_paths = list(
        field = NULL,
        global_name = "project_dirs",
        classification = "path",
        explicit_owner = "module project_dirs/experiment_paths argument",
        legacy_consumers = "report, QC, and limpa path fallbacks",
        removal_prerequisite = "all supported omics pass explicit path objects"
    ),
    design_matrix = list(
        field = "design_matrix",
        global_name = "design_matrix",
        classification = "scientific_payload",
        explicit_owner = "workflow_data$design_matrix and S4 design_matrix",
        legacy_consumers = "legacy enrichment design fallback",
        removal_prerequisite = "non-DIA proteomics enrichment migration"
    ),
    contrasts = list(
        field = "contrasts_tbl",
        global_name = "contrasts_tbl",
        classification = "scientific_payload",
        explicit_owner = "workflow_data$contrasts_tbl",
        legacy_consumers = "legacy DA and enrichment observers",
        removal_prerequisite = "non-DIA proteomics DA/enrichment migration"
    ),
    config = list(
        field = "config_list",
        global_name = "config_list",
        classification = "configuration",
        explicit_owner = "workflow_data$config_list and state generation config",
        legacy_consumers = "legacy config mutation and report fallbacks",
        removal_prerequisite = "non-DIA config mutation migration"
    ),
    annotations = list(
        field = "uniprot_dat_cln",
        global_name = "uniprot_dat_cln",
        classification = "scientific_payload",
        explicit_owner = "workflow_data$uniprot_dat_cln or artifact annotation ref",
        legacy_consumers = "legacy DA and enrichment annotation cache",
        removal_prerequisite = "non-DIA annotation migration"
    ),
    sequences = list(
        field = "aa_seq_tbl_final",
        global_name = "aa_seq_tbl_final",
        classification = "scientific_payload",
        explicit_owner = "workflow_data$aa_seq_tbl_final or artifact sequence ref",
        legacy_consumers = "legacy accession and phosphosite cleanup",
        removal_prerequisite = "non-DIA and phosphosite sequence migration"
    ),
    experiment_paths = list(
        field = NULL,
        global_name = "experiment_paths",
        classification = "path",
        explicit_owner = "module experiment_paths argument or QC session context",
        legacy_consumers = "legacy protein-QC output fallback",
        removal_prerequisite = "non-DIA protein-QC path migration"
    ),
    filtering_progress = list(
        field = NULL,
        global_name = "filtering_progress",
        classification = "helper_state",
        explicit_owner = "protein_qc_context$progress_env",
        legacy_consumers = "legacy protein-QC progress display",
        removal_prerequisite = "non-DIA protein-QC session-context migration"
    )
)

.PROT_CONTEXT_GLOBAL_USAGE <- list(
    list(
        file = "mod_prot_context_adapter_helpers.R",
        functions = paste(
            "resolveProtContextDependency, protContextLegacyGlobalExists,",
            "publishProtContextLegacyGlobal, captureProtContextLegacyGlobals,",
            "restoreProtContextLegacyGlobals"
        ),
        dependencies = names(.PROT_CONTEXT_DEPENDENCIES),
        artifact_policy = "central adapter rejects all artifact-owned access"
    ),
    list(
        file = "mod_prot_import_server_helpers.R",
        functions = "registerProtImportServerObservers",
        dependencies = c("config", "sequences"),
        artifact_policy = "injected environment reaches central publishers only"
    ),
    list(
        file = "mod_prot_import_orchestration_helpers.R",
        functions = "runProtImportPipeline",
        dependencies = c("config", "sequences"),
        artifact_policy = "injected environment reaches central publishers only"
    ),
    list(
        file = "mod_prot_import_config_helpers.R",
        functions = "storeProtImportConfiguration",
        dependencies = "config",
        artifact_policy = "central publisher rejects artifact-owned writes"
    ),
    list(
        file = "mod_prot_import_processing_helpers.R",
        functions = "processProtImportFastaData",
        dependencies = "sequences",
        artifact_policy = "central publisher rejects artifact-owned writes"
    ),
    list(
        file = "mod_prot_enrich_registration_helpers.R",
        functions = paste(
            "registerProtEnrichSelectedTabObserver,",
            "registerProtEnrichDaResultsObserver,",
            "registerProtEnrichSelectedContrastObserver"
        ),
        dependencies = "contrasts",
        artifact_policy = "artifact observers initialize from artifact indices"
    ),
    list(
        file = "mod_prot_enrich_input_helpers.R",
        functions = paste(
            "createProtEnrichRawContrastNameReactive,",
            "resolveProtEnrichRunDependencies,",
            "resolveProtEnrichUniprotAnnotations"
        ),
        dependencies = c("contrasts", "design_matrix", "annotations"),
        artifact_policy = "artifact enrichment uses explicit setup adapters"
    ),
    list(
        file = "mod_prot_enrich_analysis_helpers.R",
        functions = paste(
            "prepareProtEnrichAnalysisBodySetup, runProtEnrichAnalysisBody"
        ),
        dependencies = c("contrasts", "design_matrix", "annotations"),
        artifact_policy = "artifact enrichment bypasses legacy body setup"
    ),
    list(
        file = "mod_prot_norm_ruv_helpers.R",
        functions = "updateProtNormRuvAuditTrail",
        dependencies = "config",
        artifact_policy = "artifact DIA updates workflow_data$config_list"
    ),
    list(
        file = "mod_prot_norm_workflow_helpers.R",
        functions = "runProtNormBetweenSamplesStep",
        dependencies = "config",
        artifact_policy = "artifact DIA updates workflow_data$config_list"
    ),
    list(
        file = "mod_prot_summary_report_helpers.R",
        functions = "resolveProtSummaryReportTemplate",
        dependencies = "config",
        artifact_policy = "artifact summary disables the global fallback"
    ),
    list(
        file = "mod_prot_summary_publication_helpers.R",
        functions = "runProtSummaryPublicationCopy",
        dependencies = "project_paths",
        artifact_policy = "artifact dependencies disable the global fallback"
    ),
    list(
        file = "mod_prot_summary_session_helpers.R",
        functions = "completeProtSummaryWorkflowArgsSave",
        dependencies = "config",
        artifact_policy = "artifact dependencies disable the global write"
    ),
    list(
        file = "func_prot_qc_filtering_progress_helpers.R",
        functions = "updatePeptideFiltering",
        dependencies = c("project_paths", "filtering_progress"),
        artifact_policy = "artifact QC injects project paths and progress owner"
    ),
    list(
        file = "func_prot_limpa_qc_helpers.R",
        functions = "resolveLimpaQcSaveDir",
        dependencies = "project_paths",
        artifact_policy = "artifact QC supplies explicit project paths"
    ),
    list(
        file = "mod_prot_qc_protein_cleanup.R",
        functions = "runProteinAccessionCleanupStep",
        dependencies = c("sequences", "experiment_paths"),
        artifact_policy = "artifact-owned cleanup rejects both global fallbacks"
    )
)

protContextLegacyMode <- function() {
    match.arg(
        getOption(.PROT_CONTEXT_GLOBAL_OPTION, "enabled"),
        c("enabled", "disabled")
    )
}

protContextDependencyInventory <- function() {
    entries <- lapply(names(.PROT_CONTEXT_DEPENDENCIES), function(name) {
        specification <- .PROT_CONTEXT_DEPENDENCIES[[name]]
        data.frame(
            dependency = name,
            global_name = specification$global_name,
            classification = specification$classification,
            explicit_owner = specification$explicit_owner,
            legacy_consumers = specification$legacy_consumers,
            warning_policy = paste(
                "silent compatibility by default; disable with",
                .PROT_CONTEXT_GLOBAL_OPTION
            ),
            removal_prerequisite = specification$removal_prerequisite,
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, entries)
}

protContextGlobalUsageInventory <- function() {
    entries <- lapply(.PROT_CONTEXT_GLOBAL_USAGE, function(usage) {
        specifications <- lapply(
            usage$dependencies,
            protContextSpecification
        )
        data.frame(
            file = usage$file,
            functions = usage$functions,
            dependencies = paste(usage$dependencies, collapse = ", "),
            classification = paste(
                unique(vapply(
                    specifications,
                    `[[`,
                    character(1L),
                    "classification"
                )),
                collapse = ", "
            ),
            legacy_consumers = paste(
                unique(vapply(
                    specifications,
                    `[[`,
                    character(1L),
                    "legacy_consumers"
                )),
                collapse = "; "
            ),
            artifact_policy = usage$artifact_policy,
            warning_policy = paste(
                "silent compatibility by default; disable with",
                .PROT_CONTEXT_GLOBAL_OPTION
            ),
            removal_prerequisite = paste(
                unique(vapply(
                    specifications,
                    `[[`,
                    character(1L),
                    "removal_prerequisite"
                )),
                collapse = "; "
            ),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, entries)
}

protContextSpecification <- function(dependency) {
    specification <- .PROT_CONTEXT_DEPENDENCIES[[dependency]]
    if (is.null(specification)) {
        rlang::abort(
            sprintf("Unknown proteomics context dependency '%s'", dependency),
            class = "multischolar_unknown_prot_context_dependency"
        )
    }
    specification
}

protContextArtifactOwned <- function(workflow_data) {
    if (is.null(workflow_data)) return(FALSE)
    isTRUE(tryCatch({
        context <- workflow_data$workflow_context
        if (!inherits(context, "WorkflowContext") || !context$isBound()) {
            return(FALSE)
        }
        descriptor <- findArtifactWorkflowDescriptor(
            context$getIdentity(),
            artifactWorkflowDescriptorCatalogue()
        )
        !is.null(descriptor) &&
            artifactStageCoordinatorOwned(workflow_data, descriptor)
    }, error = \(...) FALSE))
}

protContextWorkflowValue <- function(workflow_data, field) {
    if (is.null(workflow_data) || is.null(field)) return(NULL)
    tryCatch(workflow_data[[field]], error = \(error) NULL)
}

protContextAbortMissing <- function(dependency, specification) {
    rlang::abort(
        sprintf(
            "Artifact DIA-NN dependency '%s' is missing from %s",
            dependency,
            specification$explicit_owner
        ),
        class = c(
            "multischolar_missing_prot_context_dependency",
            "multischolar_prot_context_error"
        ),
        dependency = dependency,
        explicit_owner = specification$explicit_owner
    )
}

resolveProtContextDependency <- function(
    dependency,
    workflow_data = NULL,
    explicit_value = NULL,
    required = FALSE,
    global_env = .GlobalEnv,
    exists_fn = exists,
    get_fn = get
) {
    specification <- protContextSpecification(dependency)
    value <- explicit_value %||%
        protContextWorkflowValue(workflow_data, specification$field)
    artifact_owned <- protContextArtifactOwned(workflow_data)

    if (artifact_owned) {
        if (is.null(value) && isTRUE(required)) {
            protContextAbortMissing(dependency, specification)
        }
        return(list(
            value = value,
            source = specification$explicit_owner,
            artifact_owned = TRUE,
            legacy_adapter_used = FALSE
        ))
    }

    global_available <- identical(protContextLegacyMode(), "enabled") &&
        exists_fn(
            specification$global_name,
            envir = global_env,
            inherits = FALSE
        )
    if (global_available) {
        return(list(
            value = get_fn(
                specification$global_name,
                envir = global_env,
                inherits = FALSE
            ),
            source = "legacy_global_adapter",
            artifact_owned = FALSE,
            legacy_adapter_used = TRUE
        ))
    }
    if (is.null(value) && isTRUE(required)) {
        protContextAbortMissing(dependency, specification)
    }
    list(
        value = value,
        source = if (is.null(value)) "missing" else specification$explicit_owner,
        artifact_owned = FALSE,
        legacy_adapter_used = FALSE
    )
}

protContextLegacyGlobalExists <- function(
    dependency,
    workflow_data = NULL,
    global_env = .GlobalEnv,
    exists_fn = exists
) {
    if (protContextArtifactOwned(workflow_data) ||
        identical(protContextLegacyMode(), "disabled")) {
        return(FALSE)
    }
    specification <- protContextSpecification(dependency)
    exists_fn(
        specification$global_name,
        envir = global_env,
        inherits = FALSE
    )
}

publishProtContextLegacyGlobal <- function(
    dependency,
    value,
    workflow_data = NULL,
    global_env = .GlobalEnv,
    assign_fn = assign
) {
    if (protContextArtifactOwned(workflow_data) ||
        identical(protContextLegacyMode(), "disabled")) {
        return(invisible(FALSE))
    }
    specification <- protContextSpecification(dependency)
    assign_fn(specification$global_name, value, envir = global_env)
    invisible(TRUE)
}

captureProtContextLegacyGlobals <- function(
    dependencies,
    global_env = .GlobalEnv
) {
    if (identical(protContextLegacyMode(), "disabled")) return(list())
    values <- lapply(dependencies, function(dependency) {
        specification <- protContextSpecification(dependency)
        present <- exists(
            specification$global_name,
            envir = global_env,
            inherits = FALSE
        )
        list(
            exists = present,
            value = if (present) {
                get(
                    specification$global_name,
                    envir = global_env,
                    inherits = FALSE
                )
            } else {
                NULL
            }
        )
    })
    names(values) <- dependencies
    values
}

restoreProtContextLegacyGlobals <- function(
    values,
    global_env = .GlobalEnv
) {
    if (identical(protContextLegacyMode(), "disabled")) {
        return(invisible(NULL))
    }
    for (dependency in names(values)) {
        specification <- protContextSpecification(dependency)
        name <- specification$global_name
        if (exists(name, envir = global_env, inherits = FALSE)) {
            rm(list = name, envir = global_env)
        }
        if (isTRUE(values[[dependency]]$exists)) {
            assign(name, values[[dependency]]$value, envir = global_env)
        }
    }
    invisible(NULL)
}
