# KIT-031: the lipidomics declarations, the kit's completeness test. One `newWorkflowOperation` per step of the
# lipidomics route (the step map `kit-001/step-map-lipid.json`), every one in process (`compute = "in_process"`,
# run through `runSavedOperation` on the `whole` variant), the kernels called unchanged through their legacy
# `…Fn` aliases so the recorded `@args` keep their names. The wrappers below are what the lipid modules and the
# record's replay call: the operation on the head (or a named source), its products on the checkpoint, the
# commit under the module's state name with the module's configuration and description. Nothing here reads GUI
# inputs, the global environment or a Shiny session.

lipidomicsOperationFields <- c("lipid_data", "design_matrix")

#' The merge rules of a partitioned lipid declaration: the assay tables concatenate, the design is recomputed once.
#' @noRd
lipidomicsPartitionMerge <- c(lipid_data = "concatenate", design_matrix = "design")

#' One lipidomics declaration: `LipidomicsAssayData` in, both stored components hydrated, the merge rules of a
#' partitioned unit supplied.
#' @noRd
newLipidomicsOperation <- function(name, unit, parameters, product_type, kernel, output_class = "LipidomicsAssayData",
                                   global_statistics = NULL, fixed = NULL) {
    partitioned <- unit %in% c("row", "complete_sample", "feature_group")
    newWorkflowOperation(name, "LipidomicsAssayData", fields = lipidomicsOperationFields, unit = unit, kernel = kernel,
        output_class = output_class, parameters = parameters, product_type = product_type,
        global_statistics = global_statistics, merge = if (partitioned) lipidomicsPartitionMerge else NULL, fixed = fixed)
}

#' Unique feature ids per assay (the QC tabs' count), by the id column when the table has it.
#' @noRd
lipidomicsFeatureCounts <- function(object) {
    vapply(object@lipid_data, function(assay) {
        if (object@lipid_id_column %in% names(assay)) length(unique(assay[[object@lipid_id_column]])) else nrow(assay)
    }, numeric(1))
}

#' The intensity filter (row, after the per-assay quantile pre-pass; the user's 2026-10-28 decision removed the
#' `ceiling()` from that quantile in `lipidIntensityFiltering` before any oracle was recorded).
#' @noRd
lipidIntensityFilterOperation <- function() {
    newLipidomicsOperation("lipid_intensity_filter", "row",
        parameters = c("lipids_intensity_cutoff_percentile", "lipids_proportion_of_samples_below_cutoff"),
        product_type = "lipid.intensity_filter_result",
        global_statistics = c(min_lipid_intensity_threshold =
            "the percentile of every intensity of one assay table, computed once per assay before any row"),
        kernel = function(object, parameters, store) {
            lipidIntensityFilteringFn <- lipidIntensityFiltering
            before <- lipidomicsFeatureCounts(object)
            filtered <- lipidIntensityFilteringFn(object,
                lipids_intensity_cutoff_percentile = parameters$lipids_intensity_cutoff_percentile,
                lipids_proportion_of_samples_below_cutoff = parameters$lipids_proportion_of_samples_below_cutoff)
            object <- NULL
            list(final = store$checkpoint(filtered, "lipid-intensity-filtered"),
                counts = list(before = before, after = lipidomicsFeatureCounts(filtered)))
        })
}

#' Duplicate resolution (feature group by the id column): per assay, the row with the highest mean intensity is kept.
#' @noRd
lipidDuplicateResolutionOperation <- function() {
    newLipidomicsOperation("lipid_duplicate_resolution", "feature_group", parameters = character(),
        product_type = "lipid.duplicate_resolution_result",
        kernel = function(object, parameters, store) {
            resolveLipidDuplicateFeaturesByIntensityFn <- resolveLipidDuplicateFeaturesByIntensity
            assays <- object@lipid_data
            names(assays) <- names(assays) %||% paste0("Assay_", seq_along(assays))
            before <- vapply(assays, nrow, integer(1))
            resolved <- lapply(assays, function(assay) {
                sample_cols <- names(assay)[vapply(assay, is.numeric, logical(1))]
                if (!length(sample_cols)) return(assay)
                resolveLipidDuplicateFeaturesByIntensityFn(assay_tibble = assay, id_col = object@lipid_id_column,
                    sample_cols = sample_cols)
            })
            object@lipid_data <- resolved
            list(final = store$checkpoint(object, "lipid-duplicates-resolved"),
                counts = list(before = before, after = vapply(resolved, nrow, integer(1))))
        })
}

#' A store step: the current object under a new name (the QC finalise, the correlation skip, the reset). No kernel,
#' no product; the declaration names the step for the record and the replay.
#' @noRd
lipidStoreOperation <- function(name) {
    newLipidomicsOperation(name, "store", parameters = character(), product_type = "lipid.store_step",
        kernel = function(object, parameters, store) list())
}
lipidQcFinalizeOperation <- function() lipidStoreOperation("lipid_qc_finalize")
lipidNormCompleteOperation <- function() lipidStoreOperation("lipid_norm_complete")
lipidResetOperation <- function() lipidStoreOperation("lipid_reset")

#' Internal-standard normalisation (complete sample after the mean-factor pre-pass: one aggregate standard signal per
#' sample, the mean over all samples, each sample scaled by mean / own aggregate; the standards removed afterwards).
#' @noRd
lipidItsdNormalisationOperation <- function() {
    newLipidomicsOperation("lipid_itsd_normalisation", "complete_sample",
        parameters = c("itsd_aggregation", "itsd_feature_ids"), product_type = "lipid.itsd_normalisation_result",
        global_statistics = c(mean_standard_factor =
            "the mean over all samples of each assay's per-sample aggregate standard signal, computed once before any sample is scaled"),
        fixed = list(method = "ITSD"),
        kernel = function(object, parameters, store) {
            normaliseUntransformedDataFn <- normaliseUntransformedData
            normalised <- normaliseUntransformedDataFn(theObject = object, method = "ITSD",
                itsd_aggregation = parameters$itsd_aggregation, itsd_feature_ids = parameters$itsd_feature_ids)
            object <- NULL
            list(final = store$checkpoint(normalised, "lipid-itsd-normalised"))
        })
}

#' Log2 transformation (row: `log2(x + offset)` per value).
#' @noRd
lipidLog2Operation <- function() {
    newLipidomicsOperation("lipid_log2", "row", parameters = "offset", product_type = "lipid.log2_result",
        kernel = function(object, parameters, store) {
            logTransformAssaysFn <- logTransformAssays
            transformed <- logTransformAssaysFn(theObject = object, offset = parameters$offset)
            object <- NULL
            list(final = store$checkpoint(transformed, "lipid-log2"))
        })
}

#' Between-sample normalisation (whole matrix; `none` is the module's identity, the state still committed).
#' @noRd
lipidBetweenSampleNormalisationOperation <- function() {
    newLipidomicsOperation("lipid_between_sample_normalisation", "whole_matrix", parameters = "normalisation_method",
        product_type = "lipid.between_sample_normalisation_result",
        kernel = function(object, parameters, store) {
            normaliseBetweenSamplesFn <- normaliseBetweenSamples
            normalised <- if (identical(parameters$normalisation_method, "none")) object else
                normaliseBetweenSamplesFn(theObject = object, normalisation_method = parameters$normalisation_method)
            object <- NULL
            list(final = store$checkpoint(normalised, "lipid-normalised"),
                applied = !identical(parameters$normalisation_method, "none"))
        })
}

lipidRuvParameterNames <- c("ruv_mode", "ruv_grouping_variable", "percentage_min", "percentage_max", "separation_metric",
    "k_penalty_weight", "adaptive_k_penalty", "manual_k", "manual_percentage")

#' RUV-III (whole matrix): the per-assay optimisation (automatic sweep or manual), then `ruvIII_C_Varying` with the
#' best k and the control features per assay. The optimisation results are a product of the step.
#' @noRd
lipidRuvOperation <- function() {
    newLipidomicsOperation("lipid_ruv", "whole_matrix", parameters = lipidRuvParameterNames,
        product_type = "lipid.ruv_result",
        kernel = function(object, parameters, store) {
            runLipidPerAssayRuvOptimizationFn <- runLipidPerAssayRuvOptimization
            ruvIII_C_VaryingFn <- ruvIII_C_Varying
            ruv_params <- parameters[setdiff(lipidRuvParameterNames, "ruv_mode")]
            ruv_results <- runLipidPerAssayRuvOptimizationFn(theObject = object, ruv_mode = parameters$ruv_mode,
                params = ruv_params, experiment_paths = NULL)
            best_k <- extractLipidBestKPerAssay(ruv_results)
            controls <- extractLipidCtrlPerAssay(ruv_results)
            corrected <- ruvIII_C_VaryingFn(theObject = object, ruv_grouping_variable = parameters$ruv_grouping_variable,
                ruv_number_k = best_k, ctrl = controls)
            object <- NULL
            final <- store$checkpoint(corrected, "lipid-ruv-corrected")
            corrected <- NULL
            saved_results <- lapply(ruv_results, function(result) {
                for (field in intersect(c("cancor_plot", "best_cancor_plot"), names(result))) {
                    if (inherits(result[[field]], c("gg", "ggplot"))) result[[field]] <- stripProtNormQcPlotOwners(result[[field]])
                }
                result
            })
            list(final = final, optimisation = store$product(saved_results, "optimisation", "lipid-ruv-optimisation"),
                best_k = best_k,
                control_count = lapply(controls, function(x) if (is.null(x)) NA_integer_ else if (is.logical(x)) sum(x, na.rm = TRUE) else length(x)),
                percentage = lapply(ruv_results, function(x) x$best_percentage %||% NA_real_))
        })
}

#' Pearson correlations per sample pair within the grouping variable (whole matrix; a product, no state).
#' @noRd
lipidPearsonPairsOperation <- function() {
    newLipidomicsOperation("lipid_pearson_pairs", "whole_matrix", parameters = "correlation_group",
        product_type = "lipid-pearson-pairs", output_class = NULL,
        kernel = function(object, parameters, store) {
            pearsonCorForSamplePairsFn <- pearsonCorForSamplePairs
            pearsonCorForSamplePairsFn(theObject = object, correlation_group = parameters$correlation_group)
        })
}

#' The sample selection by the correlation threshold on the pairs product (whole matrix).
#' @noRd
lipidPearsonSelectionOperation <- function() {
    newLipidomicsOperation("lipid_pearson_selection", "whole_matrix",
        parameters = c("pairs", "min_pearson_correlation_threshold"), product_type = "lipid.pearson_selection_result",
        kernel = function(object, parameters, store) {
            filterSamplesByLipidCorrelationThresholdFn <- filterSamplesByLipidCorrelationThreshold
            filtered <- filterSamplesByLipidCorrelationThresholdFn(theObject = object,
                pearson_correlation_per_pair = parameters$correlations,
                min_pearson_correlation_threshold = parameters$min_pearson_correlation_threshold)
            object <- NULL
            list(final = store$checkpoint(filtered, "lipid-correlation-filtered"))
        })
}

#' Differential abundance (whole matrix, product only): `runLipidsDA` once over every contrast and assay, its long
#' tables saved per assay and combined, the kernel's list (less the object) saved as the display product. The GUI
#' never builds the S4 result class; the manifest `lipid.da_results.v1` is the step's output.
#' @noRd
lipidDaOperation <- function() {
    # The declared names are the DA names every route's record and report overlay read (`formula`, `qvalue`, `lfc`,
    # the proteomics `DA-contrast-output`'s); the kernel receives them under its own argument names.
    newLipidomicsOperation("lipid_da_output", "whole_matrix",
        parameters = c("formula", "qvalue", "lfc", "contrasts"),
        product_type = "lipid-da-references", output_class = NULL,
        fixed = list(eBayes_trend = TRUE, eBayes_robust = TRUE),
        kernel = function(object, parameters, store) {
            runLipidsDAFn <- runLipidsDA
            result <- runLipidsDAFn(theObject = object, contrasts_tbl = parameters$contrasts,
                formula_string = parameters$formula, da_q_val_thresh = parameters$qvalue,
                treat_lfc_cutoff = parameters$lfc, eBayes_trend = TRUE, eBayes_robust = TRUE)
            object <- NULL
            result$theObject <- NULL
            accepted <- list(schema = "lipid.da_references.v1", source = store$source, tables = list(), results = list(),
                products = list(), qvalue_warnings = result$qvalue_warnings)
            if (is.data.frame(result$da_lipids_long)) accepted$tables$da_lipids_long <- store$table(as.data.frame(result$da_lipids_long))
            for (assay in names(result$per_assay_results)) {
                table <- result$per_assay_results[[assay]]
                accepted$results[[assay]] <- list(assay = assay,
                    tables = if (is.data.frame(table)) list(da_lipids_long = store$table(as.data.frame(table))) else list(),
                    rows = if (is.data.frame(table)) nrow(table) else 0L,
                    significant_counts = result$significant_counts[[assay]])
            }
            accepted$products$da_results <- store$product(result, "da_results", "lipid-da-results")
            accepted$rows <- if (is.data.frame(result$da_lipids_long)) nrow(result$da_lipids_long) else 0L
            accepted$columns <- names(result$da_lipids_long)
            accepted$model <- list(model = "limma", formula_used = parameters$formula, eBayes_trend = TRUE,
                eBayes_robust = TRUE, treat_applied = isTRUE(parameters$lfc > 0))
            accepted
        })
}

#' Identify the lipidomics DA manifest.
#' @noRd
isLipidResultManifest <- function(value) {
    is.list(value) && identical(value[["schema", exact = TRUE]], "lipid.da_results.v1")
}

#' Every lipidomics declaration, by name, with the QC tracking label of the state each commits.
#' @noRd
lipidomicsWorkflowOperations <- function() {
    declarations <- list(lipidIntensityFilterOperation(), lipidDuplicateResolutionOperation(), lipidQcFinalizeOperation(),
        lipidItsdNormalisationOperation(), lipidLog2Operation(), lipidBetweenSampleNormalisationOperation(),
        lipidRuvOperation(), lipidPearsonPairsOperation(), lipidPearsonSelectionOperation(), lipidNormCompleteOperation(),
        lipidResetOperation(), lipidDaOperation())
    stats::setNames(declarations, vapply(declarations, `[[`, "", "name"))
}

lipidomicsStepLabels <- c(lipid_intensity_filter = "2_Intensity_Filtered", lipid_duplicate_resolution = "3_Duplicates_Resolved",
    lipid_qc_finalize = "4_QC_Complete", lipid_itsd_normalisation = "5_ITSD_Normalised", lipid_log2 = "6_Log2",
    lipid_between_sample_normalisation = "7_Normalised", lipid_ruv = "8_RUV_Corrected", lipid_pearson_selection = "9_Correlation_Filtered",
    lipid_norm_complete = "9_Normalisation_Complete", lipid_reset = "0_Reset")

# ---------------------------------------------------------------------------------------------------------------
# The wrappers: the operation on a source, its product on the checkpoint, the commit under the module's state name.
# ---------------------------------------------------------------------------------------------------------------

#' Run one state-producing omics declaration and commit its checkpoint (KIT-031; KIT-032: the metabolomics wrappers
#' call it too, under its omic-neutral name).
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

#' The lipid name of `runWorkflowOmicsReferenceOperation` (KIT-032; one body).
#' @noRd
runLipidReferenceOperation <- runWorkflowOmicsReferenceOperation

#' @noRd
runLipidReferenceIntensityFilter <- function(manager, percentile, proportion, config, description = sprintf(
    "Applied lipid intensity filter (percentile: %d%%, proportion: %.2f)", as.integer(percentile), proportion),
    state_name = "lipid_intensity_filtered", hydrate = TRUE) {
    runLipidReferenceOperation(manager, lipidIntensityFilterOperation(),
        list(lipids_intensity_cutoff_percentile = percentile, lipids_proportion_of_samples_below_cutoff = proportion),
        state_name, config, description, hydrate = hydrate)
}

#' @noRd
runLipidReferenceDuplicateResolution <- function(manager, config,
    description = "Resolved duplicate lipid features by keeping highest intensity", state_name = "lipid_duplicates_resolved",
    hydrate = TRUE) {
    runLipidReferenceOperation(manager, lipidDuplicateResolutionOperation(), list(), state_name, config, description, hydrate = hydrate)
}

#' @noRd
runLipidReferenceItsdNormalisation <- function(manager, aggregation, feature_ids, config,
    description = paste("ITSD normalization (aggregation:", aggregation, ")"), state_name = "lipid_itsd_norm", hydrate = TRUE) {
    parameters <- list(itsd_aggregation = aggregation)
    parameters["itsd_feature_ids"] <- list(feature_ids)
    runLipidReferenceOperation(manager, lipidItsdNormalisationOperation(), parameters, state_name, config, description, hydrate = hydrate)
}

#' @noRd
runLipidReferenceLog2 <- function(manager, offset, config, description = paste("Log2 transformation (offset:", offset, ")"),
    state_name = "lipid_log2", hydrate = TRUE) {
    runLipidReferenceOperation(manager, lipidLog2Operation(), list(offset = offset), state_name, config, description, hydrate = hydrate)
}

#' @noRd
runLipidReferenceBetweenSampleNormalisation <- function(manager, method, config,
    description = paste("Between-sample normalization (method:", method, ")"), state_name = "lipid_normalized", hydrate = TRUE) {
    runLipidReferenceOperation(manager, lipidBetweenSampleNormalisationOperation(), list(normalisation_method = method),
        state_name, config, description, hydrate = hydrate)
}

#' @param parameters The nine RUV parameters by their declared names (the tab's `ruv_params` plus the mode and grouping).
#' @return The operation result with `optimisation`, the hydrated per-assay optimisation results.
#' @noRd
runLipidReferenceRuv <- function(manager, parameters, config, description = "RUV-III batch correction complete",
    state_name = "lipid_ruv_corrected", hydrate = TRUE) {
    result <- runLipidReferenceOperation(manager, lipidRuvOperation(), parameters[lipidRuvParameterNames], state_name, config,
        description, hydrate = hydrate)
    result$optimisation <- hydrateSavedWorkflowEntity(manager$.context$store, result$value$optimisation)$value
    result
}

#' The correlation filter: the pairs product, then the selection on it (the DIA two-operation pattern).
#' @return The selection's result with `correlations`, the per-assay correlation tables.
#' @noRd
runLipidReferenceCorrelation <- function(manager, group, threshold, config,
    description = paste("Correlation filtering (threshold:", threshold, ")"), state_name = "lipid_correlation_filtered",
    hydrate = TRUE) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the lipid step needs the session's state owner", class = "saved_entity_scope")
    withWorkflowStateScope(manager, function() {
        source <- getWorkflowStateReference(manager)
        pairs <- runSavedOperation(manager, lipidPearsonPairsOperation(), source, list(correlation_group = group))
        result <- runLipidReferenceOperation(manager, lipidPearsonSelectionOperation(),
            list(pairs = pairs$reference$entity_id, min_pearson_correlation_threshold = threshold, correlations = pairs$value),
            state_name, config, description, source = source, dependencies = list(pairs$reference), hydrate = hydrate)
        result$correlations <- pairs$value
        result
    })
}

#' Name a store step's provenance before the module's own `saveState` (the QC finalise, the correlation skip, the
#' reset): the commit is attributed to the declaration instead of `legacy`. A non-owner manager is left alone.
#' @noRd
markWorkflowOmicsStoreStep <- function(manager, operation) {
    if (!inherits(manager, "DIAStateAdapter")) return(invisible(NULL))
    entry <- workflowOperationEntry(operation)
    setWorkflowStepProvenance(manager, operation, entry$default_variant %||% "whole", entry$parameters, entry$declaration$fixed)
    invisible(NULL)
}

#' The lipid name of `markWorkflowOmicsStoreStep` (KIT-032; one body).
#' @noRd
markLipidReferenceStoreStep <- markWorkflowOmicsStoreStep

#' Commit a store step on the owner (the QC finalise, the correlation skip): the head's record under the new name
#' with the step's provenance, the prior products kept and no component re-serialised. The module's double path
#' keeps its `saveState`.
#' @noRd
commitWorkflowOmicsStoreStep <- function(manager, operation, state_name, config, description) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the step needs the session's state owner", class = "saved_entity_scope")
    withWorkflowStateScope(manager, function() {
        markWorkflowOmicsStoreStep(manager, operation)
        commitWorkflowReferenceState(manager, state_name, readWorkflowStateMetadata(manager)$data, config, description)
        getWorkflowStateReference(manager)
    })
}

#' The lipid name of `commitWorkflowOmicsStoreStep` (KIT-032; one body).
#' @noRd
commitLipidReferenceStoreStep <- commitWorkflowOmicsStoreStep

#' The design step on the owner (KIT-031): the `design_builder` provenance with the three design roles, and the
#' workflow's inputs (the import's source copies under `import_sources/` as `lipid_assay_<n>`, the design file).
#' @noRd
recordLipidDesignProvenance <- function(manager, object, source_dir = NULL, design_file = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) return(invisible(NULL))
    setWorkflowStepProvenance(manager, "design_builder",
        parameter_names = c("sample_id", "group_id", "technical_replicate_id"),
        values = list(sample_id = object@sample_id, group_id = object@group_id,
            technical_replicate_id = object@technical_replicate_id))
    paths <- character()
    if (!is.null(source_dir) && dir.exists(file.path(source_dir, "import_sources"))) {
        sources <- sort(list.files(file.path(source_dir, "import_sources"), full.names = TRUE))
        if (length(sources)) paths[paste0("lipid_assay_", seq_along(sources))] <- sources
    }
    if (is.null(design_file) && !is.null(source_dir)) design_file <- file.path(source_dir, "design_matrix.tab")
    if (!is.null(design_file) && file.exists(design_file)) paths[["design_file"]] <- design_file
    if (length(paths)) recordDIAWorkflowInputs(manager, paths)
    invisible(NULL)
}

#' Run the lipid DA declaration on a state and accept its manifest on the head (the `runDIAReferenceDaAnalysis` twin):
#' one kernel call for every contrast and assay, the manifest `lipid.da_results.v1` saved as `da-result-references`,
#' the DA step's provenance cleared (it commits no state).
#' @param manager The lipid session owner.
#' @param source Saved reference of the state the DA runs on.
#' @param settings `formula_string`, `da_q_val_thresh`, `treat_lfc_cutoff`.
#' @param contrasts The contrasts table (`contrasts`, `friendly_names`, `full_format`).
#' @return `list(manifest, results, reference)`: `results` is the kernel's list for the module's display and files,
#'   with `theObject` the hydrated source.
#' @noRd
runLipidReferenceDaAnalysis <- function(manager, source, settings, contrasts) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the lipid DA needs the session's state owner", class = "saved_entity_scope")
    if (!isSavedWorkflowReference(source)) cli::cli_abort("the lipid DA needs a saved state reference", class = "saved_entity_scope")
    withWorkflowStateScope(manager, function() {
        context <- manager$.context
        spec <- readWorkflowStateMetadata(manager, source)$data
        contrasts <- as.data.frame(contrasts, stringsAsFactors = FALSE)
        parameters <- list(formula = settings$formula_string, qvalue = settings$da_q_val_thresh,
            lfc = settings$treat_lfc_cutoff, contrasts = contrasts)
        run <- runSavedOperation(manager, lipidDaOperation(), source, parameters)
        if (identical(run$status, "computed")) {
            snapshot <- context$snapshot
            snapshot$da_dependencies[[run$reference$entity_id]] <- collectDIANativeReferences(run$value)
            publishDIAStateHead(context, snapshot)
        }
        clearDIAStepProvenance(manager)
        manifest <- list(schema = "lipid.da_results.v1", source = source, lipid_id_column = spec$slots$lipid_id_column,
            sample_id = spec$slots$sample_id, contrasts = contrasts, results = run$value$results, tables = run$value$tables,
            products = run$value$products, parameters = parameters[c("formula", "qvalue", "lfc")],
            qvalue_warnings = run$value$qvalue_warnings, model = run$value$model)
        reference <- saveDIAStateValue(context, manifest, "da-result-references", codec = "rds")
        snapshot <- context$snapshot
        snapshot$da_results <- reference
        snapshot$da_dependencies[[reference$entity_id]] <- collectDIANativeReferences(list(source, run$value))
        publishDIAStateHead(context, snapshot)
        results <- hydrateSavedWorkflowEntity(context$store, run$value$products$da_results)$value
        results$theObject <- readWorkflowExecutionComponents(manager, source, fields = lipidomicsOperationFields)
        list(manifest = manifest, results = results, reference = reference)
    })
}

#' Replay one recorded lipid step (KIT-031): the record's operation with its recorded values through the wrapper the
#' module used, committed under the recorded state name with the head's configuration and the recorded description.
#' @noRd
replayLipidomicsStep <- function(manager, step, parameters) {
    operation <- step$operation
    config <- readWorkflowStateMetadata(manager)$config
    description <- step$description %||% operation
    state <- step$state
    unsupported <- function() cli::cli_abort("step {.val {state}} cannot be replayed: operation {.val {operation}}",
        class = "workflow_replay_unsupported")
    switch(operation,
        lipid_intensity_filter = runLipidReferenceIntensityFilter(manager, parameters$lipids_intensity_cutoff_percentile,
            parameters$lipids_proportion_of_samples_below_cutoff, config, description, state_name = state, hydrate = FALSE),
        lipid_duplicate_resolution = runLipidReferenceDuplicateResolution(manager, config, description, state_name = state, hydrate = FALSE),
        lipid_qc_finalize = , lipid_norm_complete = commitLipidReferenceStoreStep(manager, operation, state, config, description),
        lipid_itsd_normalisation = runLipidReferenceItsdNormalisation(manager, parameters$itsd_aggregation,
            parameters$itsd_feature_ids, config, description, state_name = state, hydrate = FALSE),
        lipid_log2 = runLipidReferenceLog2(manager, parameters$offset, config, description, state_name = state, hydrate = FALSE),
        lipid_between_sample_normalisation = runLipidReferenceBetweenSampleNormalisation(manager, parameters$normalisation_method,
            config, description, state_name = state, hydrate = FALSE),
        lipid_ruv = runLipidReferenceRuv(manager, parameters, config, description, state_name = state, hydrate = FALSE),
        lipid_pearson_selection = {
            previous <- Filter(function(p) identical(p$operation, "lipid_pearson_pairs"), step$preceded_by %||% list())
            if (!length(previous)) unsupported()
            pairs <- workflowReplayParameters(previous[[1L]]$parameters, workflowOperationEntry("lipid_pearson_pairs")$parameters)
            runLipidReferenceCorrelation(manager, pairs$correlation_group, parameters$min_pearson_correlation_threshold, config,
                description, state_name = state, hydrate = FALSE)
        },
        unsupported())
    invisible(NULL)
}
