# KIT-073: the shared declaration core of the lipidomics and metabolomics routes (KIT-031/032 wrote them as two
# near-identical files). One `newWorkflowOperation` per step, every one in process (`compute = "in_process"`, run
# through `runSavedOperation` on the `whole` variant), built from a per-omic binding: the class, the stored slots,
# the names and strings each omic records, and the small kernel adaptors where the two omics call their kernels
# differently (utils_workflow_lipid_declarations.R, utils_workflow_metab_declarations.R). The kernels are unchanged and
# called under the names the modules call them by, so the recorded `@args` keep their keys. The wrappers are what the
# modules and the record's replay call: the operation on the head (or a named source), its products on the checkpoint,
# the commit under the module's state name with the module's configuration and description. Nothing here reads GUI
# inputs, the global environment or a Shiny session.

#' One omics declaration: the binding's class in, both stored components hydrated, the merge rules of a partitioned
#' unit supplied.
#' @noRd
newOmicsOperation <- function(binding, name, unit, parameters, product_type, kernel, output_class = binding$class,
                              global_statistics = NULL, fixed = NULL) {
    partitioned <- unit %in% c("row", "complete_sample", "feature_group")
    newWorkflowOperation(name, binding$class, fields = binding$fields, unit = unit, kernel = kernel,
        output_class = output_class, parameters = parameters, product_type = product_type,
        global_statistics = global_statistics, merge = if (partitioned) binding$merge else NULL, fixed = fixed)
}

#' Unique feature ids per assay (the QC tabs' count), by the id column when the table has it.
#' @noRd
omicsFeatureCounts <- function(binding, object) {
    id_column <- methods::slot(object, binding$id_slot)
    vapply(methods::slot(object, binding$data_slot), function(assay) {
        if (id_column %in% names(assay)) length(unique(assay[[id_column]])) else nrow(assay)
    }, numeric(1))
}

#' The intensity filter (row, after the per-assay quantile pre-pass; the user's 2026-10-28 decision removed the
#' `ceiling()` from that quantile in both kernels before any oracle was recorded).
#' @noRd
omicsIntensityFilterOperation <- function(binding) {
    newOmicsOperation(binding, paste0(binding$prefix, "_intensity_filter"), "row",
        parameters = binding$intensity_parameters,
        product_type = paste0(binding$prefix, ".intensity_filter_result"),
        global_statistics = stats::setNames(
            "the percentile of every intensity of one assay table, computed once per assay before any row",
            binding$intensity_threshold),
        kernel = function(object, parameters, store) {
            before <- omicsFeatureCounts(binding, object)
            filtered <- binding$filter_intensity(object, parameters)
            object <- NULL
            list(final = store$checkpoint(filtered, paste0(binding$prefix, "-intensity-filtered")),
                counts = list(before = before, after = omicsFeatureCounts(binding, filtered)))
        })
}

#' Duplicate resolution (feature group by the id column), through the omic's adaptor (`resolve_duplicates` returns
#' the resolved object and the step's statistics under the name the module reads).
#' @noRd
omicsDuplicateResolutionOperation <- function(binding) {
    newOmicsOperation(binding, paste0(binding$prefix, "_duplicate_resolution"), "feature_group", parameters = character(),
        product_type = paste0(binding$prefix, ".duplicate_resolution_result"),
        kernel = function(object, parameters, store) {
            resolved <- binding$resolve_duplicates(object)
            c(list(final = store$checkpoint(resolved$object, paste0(binding$prefix, "-duplicates-resolved"))), resolved$value)
        })
}

#' A store step: the current object under a new name (the QC finalise, the correlation skip, the reset). No kernel,
#' no product; the declaration names the step for the record and the replay.
#' @noRd
omicsStoreOperation <- function(binding, name) {
    newOmicsOperation(binding, name, "store", parameters = character(), product_type = paste0(binding$prefix, ".store_step"),
        kernel = function(object, parameters, store) list())
}

#' Internal-standard normalisation (complete sample after the mean-factor pre-pass: one aggregate standard signal per
#' sample, the mean over all samples, each sample scaled by mean / own aggregate; the standards removed afterwards).
#' @noRd
omicsItsdNormalisationOperation <- function(binding) {
    newOmicsOperation(binding, paste0(binding$prefix, "_itsd_normalisation"), "complete_sample",
        parameters = c("itsd_aggregation", "itsd_feature_ids"), product_type = paste0(binding$prefix, ".itsd_normalisation_result"),
        global_statistics = c(mean_standard_factor =
            "the mean over all samples of each assay's per-sample aggregate standard signal, computed once before any sample is scaled"),
        fixed = list(method = "ITSD"),
        kernel = function(object, parameters, store) {
            normaliseUntransformedDataFn <- normaliseUntransformedData
            normalised <- normaliseUntransformedDataFn(theObject = object, method = "ITSD",
                itsd_aggregation = parameters$itsd_aggregation, itsd_feature_ids = parameters$itsd_feature_ids)
            object <- NULL
            list(final = store$checkpoint(normalised, paste0(binding$prefix, "-itsd-normalised")))
        })
}

#' Log2 transformation (row: `log2(x + offset)` per value).
#' @noRd
omicsLog2Operation <- function(binding) {
    newOmicsOperation(binding, paste0(binding$prefix, "_log2"), "row", parameters = "offset",
        product_type = paste0(binding$prefix, ".log2_result"),
        kernel = function(object, parameters, store) {
            logTransformAssaysFn <- logTransformAssays
            transformed <- logTransformAssaysFn(theObject = object, offset = parameters$offset)
            object <- NULL
            list(final = store$checkpoint(transformed, paste0(binding$prefix, "-log2")))
        })
}

#' Between-sample normalisation (whole matrix; `none` is the module's identity, the state still committed).
#' @noRd
omicsBetweenSampleNormalisationOperation <- function(binding) {
    newOmicsOperation(binding, paste0(binding$prefix, "_between_sample_normalisation"), "whole_matrix",
        parameters = "normalisation_method", product_type = paste0(binding$prefix, ".between_sample_normalisation_result"),
        kernel = function(object, parameters, store) {
            normaliseBetweenSamplesFn <- normaliseBetweenSamples
            normalised <- if (identical(parameters$normalisation_method, "none")) object else
                normaliseBetweenSamplesFn(theObject = object, normalisation_method = parameters$normalisation_method)
            object <- NULL
            list(final = store$checkpoint(normalised, paste0(binding$prefix, "-normalised")),
                applied = !identical(parameters$normalisation_method, "none"))
        })
}

omicsRuvParameterNames <- c("ruv_mode", "ruv_grouping_variable", "percentage_min", "percentage_max", "separation_metric",
    "k_penalty_weight", "adaptive_k_penalty", "manual_k", "manual_percentage")

#' RUV-III (whole matrix): the per-assay optimisation (automatic sweep or manual), then `ruvIII_C_Varying` with the
#' best k and the control features per assay (the module's optimisation and correction in one declaration). The
#' optimisation results are a product of the step.
#' @noRd
omicsRuvOperation <- function(binding) {
    newOmicsOperation(binding, paste0(binding$prefix, "_ruv"), "whole_matrix", parameters = omicsRuvParameterNames,
        product_type = paste0(binding$prefix, ".ruv_result"),
        kernel = function(object, parameters, store) {
            ruvIII_C_VaryingFn <- ruvIII_C_Varying
            ruv_results <- binding$optimise_ruv(object, parameters$ruv_mode, parameters[setdiff(omicsRuvParameterNames, "ruv_mode")])
            best_k <- binding$extract_best_k(ruv_results)
            controls <- binding$extract_controls(ruv_results)
            corrected <- ruvIII_C_VaryingFn(theObject = object, ruv_grouping_variable = parameters$ruv_grouping_variable,
                ruv_number_k = best_k, ctrl = controls)
            object <- NULL
            final <- store$checkpoint(corrected, paste0(binding$prefix, "-ruv-corrected"))
            corrected <- NULL
            saved_results <- lapply(ruv_results, function(result) {
                for (field in intersect(c("cancor_plot", "best_cancor_plot"), names(result))) {
                    if (inherits(result[[field]], c("gg", "ggplot"))) result[[field]] <- stripProtNormQcPlotOwners(result[[field]])
                }
                result
            })
            list(final = final, optimisation = store$product(saved_results, "optimisation", paste0(binding$prefix, "-ruv-optimisation")),
                best_k = best_k,
                control_count = lapply(controls, function(x) if (is.null(x)) NA_integer_ else if (is.logical(x)) sum(x, na.rm = TRUE) else length(x)),
                percentage = lapply(ruv_results, function(x) x$best_percentage %||% NA_real_))
        })
}

#' Pearson correlations per sample pair within the grouping variable (whole matrix; a product, no state).
#' @noRd
omicsPearsonPairsOperation <- function(binding) {
    newOmicsOperation(binding, paste0(binding$prefix, "_pearson_pairs"), "whole_matrix", parameters = "correlation_group",
        product_type = paste0(binding$prefix, "-pearson-pairs"), output_class = NULL,
        kernel = function(object, parameters, store) {
            pearsonCorForSamplePairsFn <- pearsonCorForSamplePairs
            pearsonCorForSamplePairsFn(theObject = object, correlation_group = parameters$correlation_group)
        })
}

#' The sample selection by the correlation threshold on the pairs product (whole matrix).
#' @noRd
omicsPearsonSelectionOperation <- function(binding) {
    newOmicsOperation(binding, paste0(binding$prefix, "_pearson_selection"), "whole_matrix",
        parameters = c("pairs", "min_pearson_correlation_threshold"),
        product_type = paste0(binding$prefix, ".pearson_selection_result"),
        kernel = function(object, parameters, store) {
            filtered <- binding$filter_correlation(object, parameters$correlations, parameters$min_pearson_correlation_threshold)
            object <- NULL
            list(final = store$checkpoint(filtered, paste0(binding$prefix, "-correlation-filtered")))
        })
}

#' Differential abundance (whole matrix, product only): the omic's DA kernel once over every contrast and assay, its
#' long tables saved per assay and combined, the kernel's list (less the object) saved as the display product. The
#' GUI never builds the S4 result class; the manifest (`<prefix>.da_results.v1`) is the step's output.
#' @noRd
omicsDaOperation <- function(binding) {
    # The declared names are the DA names every route's record and report overlay read (`formula`, `qvalue`, `lfc`,
    # the proteomics `DA-contrast-output`'s); the kernel receives them under its own argument names.
    long_table <- binding$long_table
    newOmicsOperation(binding, paste0(binding$prefix, "_da_output"), "whole_matrix",
        parameters = c("formula", "qvalue", "lfc", "contrasts"),
        product_type = paste0(binding$prefix, "-da-references"), output_class = NULL,
        fixed = list(eBayes_trend = TRUE, eBayes_robust = TRUE),
        kernel = function(object, parameters, store) {
            result <- binding$run_da(object, parameters)
            object <- NULL
            result$theObject <- NULL
            accepted <- list(schema = paste0(binding$prefix, ".da_references.v1"), source = store$source, tables = list(),
                results = list(), products = list(), qvalue_warnings = result$qvalue_warnings)
            if (is.data.frame(result[[long_table]])) accepted$tables[[long_table]] <- store$table(as.data.frame(result[[long_table]]))
            for (assay in names(result$per_assay_results)) {
                table <- result$per_assay_results[[assay]]
                accepted$results[[assay]] <- list(assay = assay,
                    tables = if (is.data.frame(table)) stats::setNames(list(store$table(as.data.frame(table))), long_table) else list(),
                    rows = if (is.data.frame(table)) nrow(table) else 0L,
                    significant_counts = result$significant_counts[[assay]])
            }
            accepted$products$da_results <- store$product(result, "da_results", paste0(binding$prefix, "-da-results"))
            accepted$rows <- if (is.data.frame(result[[long_table]])) nrow(result[[long_table]]) else 0L
            accepted$columns <- names(result[[long_table]])
            accepted$model <- list(model = "limma", formula_used = parameters$formula, eBayes_trend = TRUE,
                eBayes_robust = TRUE, treat_applied = isTRUE(parameters$lfc > 0))
            accepted
        })
}

#' Every declaration of one omic, by name.
#' @noRd
omicsWorkflowOperations <- function(binding) {
    prefix <- binding$prefix
    declarations <- list(omicsIntensityFilterOperation(binding), omicsDuplicateResolutionOperation(binding),
        omicsStoreOperation(binding, paste0(prefix, "_qc_finalize")), omicsItsdNormalisationOperation(binding),
        omicsLog2Operation(binding), omicsBetweenSampleNormalisationOperation(binding), omicsRuvOperation(binding),
        omicsPearsonPairsOperation(binding), omicsPearsonSelectionOperation(binding),
        omicsStoreOperation(binding, paste0(prefix, "_norm_complete")), omicsStoreOperation(binding, paste0(prefix, "_reset")),
        omicsDaOperation(binding))
    stats::setNames(declarations, vapply(declarations, `[[`, "", "name"))
}

#' The omics bindings the registry and the record read, by omic type.
#' @noRd
omicsDeclarationBindings <- function() {
    list(lipidomics = lipidomicsDeclarationBinding, metabolomics = metabolomicsDeclarationBinding)
}

# ---------------------------------------------------------------------------------------------------------------
# The wrappers: the operation on a source, its product on the checkpoint, the commit under the module's state name.
# ---------------------------------------------------------------------------------------------------------------

#' Run one state-producing omics declaration and commit its checkpoint (KIT-031; KIT-032: both omics).
#' @param manager The session owner (`DIAStateAdapter`).
#' @param declaration The declaration.
#' @param parameters Every declared parameter and any extra the kernel reads.
#' @param state_name,config,description The module's commit.
#' @param source Saved reference of the input state; the head when NULL.
#' @return `list(run, reference, value, object)`: the operation run, the committed state's reference, the kernel's
#'   value and the committed object hydrated for the module's display (the modules hydrate whole objects per read).
#' @noRd
runWorkflowOmicsReferenceOperation <- function(manager, declaration, parameters, state_name, config, description, source = NULL,
    dependencies = list(), hydrate = TRUE) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the step needs the session's state owner", class = "saved_entity_scope")
    withWorkflowStateScope(manager, function() {
        source <- source %||% getWorkflowStateReference(manager)
        run <- runSavedOperation(manager, declaration, source, parameters, dependencies)
        final <- readWorkflowStateMetadata(manager, run$value$final)$data
        final$products[[paste0(declaration$name, "_science")]] <- run$reference
        commitWorkflowReferenceState(manager, state_name, final, config, description)
        list(run = run, reference = getWorkflowStateReference(manager), value = run$value,
            object = if (hydrate) manager$getState(state_name) else NULL)
    })
}

#' @noRd
runWorkflowOmicsReferenceIntensityFilter <- function(binding, manager, percentile, proportion, config, description = NULL,
    state_name = NULL, hydrate = TRUE) {
    description <- description %||% sprintf("Applied %s intensity filter (percentile: %d%%, proportion: %.2f)",
        binding$feature_noun, as.integer(percentile), proportion)
    runWorkflowOmicsReferenceOperation(manager, omicsIntensityFilterOperation(binding),
        stats::setNames(list(percentile, proportion), binding$intensity_parameters),
        state_name %||% paste0(binding$prefix, "_intensity_filtered"), config, description, hydrate = hydrate)
}

#' @noRd
runWorkflowOmicsReferenceDuplicateResolution <- function(binding, manager, config, description = NULL, state_name = NULL,
    hydrate = TRUE) {
    description <- description %||% sprintf("Resolved duplicate %s features by keeping highest intensity", binding$feature_noun)
    runWorkflowOmicsReferenceOperation(manager, omicsDuplicateResolutionOperation(binding), list(),
        state_name %||% paste0(binding$prefix, "_duplicates_resolved"), config, description, hydrate = hydrate)
}

#' @noRd
runWorkflowOmicsReferenceItsdNormalisation <- function(binding, manager, aggregation, feature_ids, config, description = NULL,
    state_name = NULL, hydrate = TRUE) {
    parameters <- list(itsd_aggregation = aggregation)
    parameters["itsd_feature_ids"] <- list(feature_ids)
    runWorkflowOmicsReferenceOperation(manager, omicsItsdNormalisationOperation(binding), parameters,
        state_name %||% paste0(binding$prefix, "_itsd_norm"), config,
        description %||% paste("ITSD normalization (aggregation:", aggregation, ")"), hydrate = hydrate)
}

#' @noRd
runWorkflowOmicsReferenceLog2 <- function(binding, manager, offset, config, description = NULL, state_name = NULL, hydrate = TRUE) {
    runWorkflowOmicsReferenceOperation(manager, omicsLog2Operation(binding), list(offset = offset),
        state_name %||% paste0(binding$prefix, "_log2"), config,
        description %||% paste("Log2 transformation (offset:", offset, ")"), hydrate = hydrate)
}

#' @noRd
runWorkflowOmicsReferenceBetweenSampleNormalisation <- function(binding, manager, method, config, description = NULL,
    state_name = NULL, hydrate = TRUE) {
    runWorkflowOmicsReferenceOperation(manager, omicsBetweenSampleNormalisationOperation(binding), list(normalisation_method = method),
        state_name %||% paste0(binding$prefix, "_normalized"), config,
        description %||% paste("Between-sample normalization (method:", method, ")"), hydrate = hydrate)
}

#' @param parameters The nine RUV parameters by their declared names (the tab's `ruv_params` plus the mode and grouping).
#' @return The operation result with `optimisation`, the hydrated per-assay optimisation results.
#' @noRd
runWorkflowOmicsReferenceRuv <- function(binding, manager, parameters, config, description = NULL, state_name = NULL,
    hydrate = TRUE) {
    result <- runWorkflowOmicsReferenceOperation(manager, omicsRuvOperation(binding), parameters[omicsRuvParameterNames],
        state_name %||% paste0(binding$prefix, "_ruv_corrected"), config, description %||% "RUV-III batch correction complete",
        hydrate = hydrate)
    result$optimisation <- hydrateSavedWorkflowEntity(manager$.context$store, result$value$optimisation)$value
    result
}

#' The correlation filter: the pairs product, then the selection on it (the DIA two-operation pattern).
#' @return The selection's result with `correlations`, the per-assay correlation tables.
#' @noRd
runWorkflowOmicsReferenceCorrelation <- function(binding, manager, group, threshold, config, description = NULL,
    state_name = NULL, hydrate = TRUE) {
    if (!inherits(manager, "DIAStateAdapter")) {
        cli::cli_abort(paste(binding$step_scope, "needs the session's state owner"), class = "saved_entity_scope")
    }
    withWorkflowStateScope(manager, function() {
        source <- getWorkflowStateReference(manager)
        pairs <- runSavedOperation(manager, omicsPearsonPairsOperation(binding), source, list(correlation_group = group))
        result <- runWorkflowOmicsReferenceOperation(manager, omicsPearsonSelectionOperation(binding),
            list(pairs = pairs$reference$entity_id, min_pearson_correlation_threshold = threshold, correlations = pairs$value),
            state_name %||% paste0(binding$prefix, "_correlation_filtered"), config,
            description %||% paste("Correlation filtering (threshold:", threshold, ")"), source = source,
            dependencies = list(pairs$reference), hydrate = hydrate)
        result$correlations <- pairs$value
        result
    })
}

#' Name a store step's provenance before the module's own `saveState` (the reset): the commit is attributed to the
#' declaration instead of `legacy`. A non-owner manager is refused (KIT-072).
#' @noRd
markWorkflowOmicsStoreStep <- function(manager, operation) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the step needs the session's state owner", class = "saved_entity_scope")
    entry <- workflowOperationEntry(operation)
    setWorkflowStepProvenance(manager, operation, entry$default_variant %||% "whole", entry$parameters, entry$declaration$fixed)
    invisible(NULL)
}

#' Commit a store step on the owner (the QC finalise, the correlation skip): the head's record under the new name
#' with the step's provenance, the prior products kept and no component re-serialised.
#' @noRd
commitWorkflowOmicsStoreStep <- function(manager, operation, state_name, config, description) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the step needs the session's state owner", class = "saved_entity_scope")
    withWorkflowStateScope(manager, function() {
        markWorkflowOmicsStoreStep(manager, operation)
        commitWorkflowReferenceState(manager, state_name, readWorkflowStateMetadata(manager)$data, config, description)
        getWorkflowStateReference(manager)
    })
}

#' The design step on the owner (KIT-031/032): the `design_builder` provenance with the three design roles, and the
#' workflow's inputs (the import's source copies under `import_sources/` as `<input prefix><n>`, the design file).
#' @noRd
recordWorkflowOmicsDesignProvenance <- function(binding, manager, object, source_dir = NULL, design_file = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the design step needs the session's state owner", class = "saved_entity_scope")
    setWorkflowStepProvenance(manager, "design_builder",
        parameter_names = c("sample_id", "group_id", "technical_replicate_id"),
        values = list(sample_id = object@sample_id, group_id = object@group_id,
            technical_replicate_id = object@technical_replicate_id))
    paths <- character()
    if (!is.null(source_dir) && dir.exists(file.path(source_dir, "import_sources"))) {
        sources <- sort(list.files(file.path(source_dir, "import_sources"), full.names = TRUE))
        if (length(sources)) paths[paste0(binding$design_input_prefix, seq_along(sources))] <- sources
    }
    if (is.null(design_file) && !is.null(source_dir)) design_file <- file.path(source_dir, "design_matrix.tab")
    if (!is.null(design_file) && file.exists(design_file)) paths[["design_file"]] <- design_file
    if (length(paths)) recordDIAWorkflowInputs(manager, paths)
    invisible(NULL)
}

#' Run the omic's DA declaration on a state and accept its manifest on the head (the `runDIAReferenceDaAnalysis`
#' twin): one kernel call for every contrast and assay, the manifest `<prefix>.da_results.v1` saved as
#' `da-result-references`, the DA step's provenance cleared (it commits no state).
#' @param manager The session owner.
#' @param source Saved reference of the state the DA runs on.
#' @param settings `formula_string`, `da_q_val_thresh`, `treat_lfc_cutoff`.
#' @param contrasts The contrasts table (`contrasts`, `friendly_names`, `full_format`).
#' @return `list(manifest, results, reference)`: `results` is the kernel's list for the module's display and files,
#'   with `theObject` the hydrated source.
#' @noRd
runWorkflowOmicsReferenceDaAnalysis <- function(binding, manager, source, settings, contrasts) {
    if (!inherits(manager, "DIAStateAdapter")) {
        cli::cli_abort(paste(binding$da_scope, "needs the session's state owner"), class = "saved_entity_scope")
    }
    if (!isSavedWorkflowReference(source)) {
        cli::cli_abort(paste(binding$da_scope, "needs a saved state reference"), class = "saved_entity_scope")
    }
    withWorkflowStateScope(manager, function() {
        context <- manager$.context
        spec <- readWorkflowStateMetadata(manager, source)$data
        contrasts <- as.data.frame(contrasts, stringsAsFactors = FALSE)
        parameters <- list(formula = settings$formula_string, qvalue = settings$da_q_val_thresh,
            lfc = settings$treat_lfc_cutoff, contrasts = contrasts)
        run <- runSavedOperation(manager, omicsDaOperation(binding), source, parameters)
        if (identical(run$status, "computed")) {
            snapshot <- context$snapshot
            snapshot$da_dependencies[[run$reference$entity_id]] <- collectDIANativeReferences(run$value)
            publishDIAStateHead(context, snapshot)
        }
        clearDIAStepProvenance(manager)
        manifest <- c(list(schema = paste0(binding$prefix, ".da_results.v1"), source = source),
            stats::setNames(list(spec$slots[[binding$id_slot]]), binding$id_slot),
            list(sample_id = spec$slots$sample_id, contrasts = contrasts, results = run$value$results, tables = run$value$tables,
                products = run$value$products, parameters = parameters[c("formula", "qvalue", "lfc")],
                qvalue_warnings = run$value$qvalue_warnings, model = run$value$model))
        reference <- saveDIAStateValue(context, manifest, "da-result-references", codec = "rds")
        snapshot <- context$snapshot
        snapshot$da_results <- reference
        snapshot$da_dependencies[[reference$entity_id]] <- collectDIANativeReferences(list(source, run$value))
        publishDIAStateHead(context, snapshot)
        results <- hydrateSavedWorkflowEntity(context$store, run$value$products$da_results)$value
        results$theObject <- readWorkflowExecutionComponents(manager, source, fields = binding$fields)
        list(manifest = manifest, results = results, reference = reference)
    })
}

#' Replay one recorded omics step (KIT-031/032): the record's operation with its recorded values through the wrapper
#' the module used, committed under the recorded state name with the head's configuration and the recorded
#' description.
#' @noRd
replayWorkflowOmicsStep <- function(binding, manager, step, parameters) {
    operation <- step$operation
    config <- readWorkflowStateMetadata(manager)$config
    description <- step$description %||% operation
    state <- step$state
    unsupported <- function() cli::cli_abort("step {.val {state}} cannot be replayed: operation {.val {operation}}",
        class = "workflow_replay_unsupported")
    prefix <- paste0(binding$prefix, "_")
    if (!startsWith(operation, prefix)) unsupported()
    switch(substring(operation, nchar(prefix) + 1L),
        intensity_filter = runWorkflowOmicsReferenceIntensityFilter(binding, manager, parameters[[binding$intensity_parameters[[1L]]]],
            parameters[[binding$intensity_parameters[[2L]]]], config, description, state_name = state, hydrate = FALSE),
        duplicate_resolution = runWorkflowOmicsReferenceDuplicateResolution(binding, manager, config, description,
            state_name = state, hydrate = FALSE),
        qc_finalize = , norm_complete = commitWorkflowOmicsStoreStep(manager, operation, state, config, description),
        itsd_normalisation = runWorkflowOmicsReferenceItsdNormalisation(binding, manager, parameters$itsd_aggregation,
            parameters$itsd_feature_ids, config, description, state_name = state, hydrate = FALSE),
        log2 = runWorkflowOmicsReferenceLog2(binding, manager, parameters$offset, config, description, state_name = state,
            hydrate = FALSE),
        between_sample_normalisation = runWorkflowOmicsReferenceBetweenSampleNormalisation(binding, manager,
            parameters$normalisation_method, config, description, state_name = state, hydrate = FALSE),
        ruv = runWorkflowOmicsReferenceRuv(binding, manager, parameters, config, description, state_name = state, hydrate = FALSE),
        pearson_selection = {
            pairs_operation <- paste0(prefix, "pearson_pairs")
            previous <- Filter(function(p) identical(p$operation, pairs_operation), step$preceded_by %||% list())
            if (!length(previous)) unsupported()
            pairs <- workflowReplayParameters(previous[[1L]]$parameters, workflowOperationEntry(pairs_operation)$parameters)
            runWorkflowOmicsReferenceCorrelation(binding, manager, pairs$correlation_group, parameters$min_pearson_correlation_threshold,
                config, description, state_name = state, hydrate = FALSE)
        },
        unsupported())
    invisible(NULL)
}
