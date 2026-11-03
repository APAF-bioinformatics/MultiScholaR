# KIT-073: the shared declaration core of the lipidomics and metabolomics routes (KIT-031/032 wrote them as two
# near-identical files). One `newWorkflowOperation` per step, every one in process (`compute = "in_process"`, run
# through `runSavedOperation` on the `whole` variant), built from a per-omic binding: the class, the stored slots,
# the names and strings each omic records, and the small kernel adaptors where the two omics call their kernels
# differently (utils_workflow_lipid_declarations.R, utils_workflow_metab_declarations.R). The kernels are unchanged and
# called under the names the modules call them by, so the recorded `@args` keep their keys. The wrappers are what the
# modules and the record's replay call: the operation on the head (or a named source), its products on the checkpoint,
# the commit under the module's state name with the module's configuration and description. Nothing here reads GUI
# inputs, the global environment or a Shiny session.

#' The declarations that run on the worker's `omics` role (KIT-078), by the name after the omic's prefix.
#' @noRd
omicsWorkerSteps <- c("intensity_filter", "duplicate_resolution", "itsd_normalisation", "log2", "between_sample_normalisation",
    "ruv", "pearson_pairs", "pearson_selection", "da_contrast")

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
            # KIT-078: the declared `pairs` entity, read through the store in either process (it was the undeclared
            # live `correlations`).
            manager <- store$manager
            pairs <- manager$.context$snapshot$products[[parameters$pairs]]
            if (is.null(pairs)) cli::cli_abort("the Pearson pairs product is not in the session", class = "workflow_operation_parameters")
            correlations <- hydrateSavedWorkflowEntity(manager$.context$store, pairs)$value
            filtered <- binding$filter_correlation(object, correlations, parameters$min_pearson_correlation_threshold)
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

#' Identify one assay's uncontrasted limma fit (KIT-076): the values in order, the design and formula, the software.
#' New code: the proteomics `newWorkflowDiaBaseFitIdentity` hashes protein-only inputs.
#' @param data Feature-by-sample matrix.
#' @param design Design table of the run.
#' @param formula_string The model formula.
#' @param assay The assay's name.
#' @param software Installed software identity.
#' @return The store's scientific identity.
#' @noRd
newWorkflowOmicsDaFitIdentity <- function(data, design, formula_string, assay, software) {
    newSavedWorkflowScientificIdentity(
        digest::digest(list(values = data, feature_order = rownames(data), sample_order = colnames(data)), algo = "sha256"),
        digest::digest(list(formula = formula_string, design = design), algo = "sha256"),
        digest::digest(list(layer = "omics_da_base_fit"), algo = "sha256"),
        digest::digest(list(assay = assay), algo = "sha256"), software)
}

#' One assay's fit for the contrast operations (KIT-076): the saved `da_fit` when the identity has one, else the fit half
#' of the omic's kernel, saved once. A fit failure drops the assay (the orchestrator's silent drop, kept for the
#' fit half only and recorded in the KIT-076 outcome); ids the store's envelope refuses keep the fit unsaved.
#' @return `list(fit, reference)`, or NULL when the fit failed.
#' @noRd
omicsDaAssayFit <- function(binding, store, data, design, formula_string, sample_id, assay) {
    manager <- store$manager
    identity <- newWorkflowOmicsDaFitIdentity(data, design, formula_string, assay, manager$.context$software)
    factors <- list(route = binding$prefix, layer = "base_fit")
    hit <- resolveSavedWorkflowEntity(manager$.context$store, "da_fit", identity, factors = factors)
    if (identical(hit$status, "hit")) return(list(fit = hit$value$fit, reference = hit$reference))
    fit <- tryCatch(binding$da_fit(data, design, formula_string, sample_id), error = function(e) {
        cat(sprintf("[ERROR] Processing assay %s: %s\n", assay, e$message))
        logger::log_error(sprintf("   ERROR processing assay %s: %s", assay, e$message))
        NULL
    })
    if (is.null(fit)) return(NULL)
    features <- rownames(fit$coefficients)
    if (anyDuplicated(features) || anyNA(features) || !all(nzchar(features))) return(list(fit = fit, reference = NULL))
    reference <- saveWorkflowConsumerEntity(manager$.context$store,
        newSavedWorkflowDaFit(fit, features, rownames(fit$design), colnames(fit$design)),
        basename(tempfile(paste0(binding$prefix, "-da-fit-"))), identity, "da_fit", factors = factors)
    rootWorkflowProduct(manager, reference, list(store$source))
    list(fit = fit, reference = reference)
}

#' Differential abundance for one contrast (KIT-076; whole matrix, product only): per assay the saved fit, the contrast
#' half of the omic's kernel, the orchestrator's table. A failure of this contrast on an assay is recorded under
#' `errors` and costs nothing else; store errors are not caught. The product keeps what the module's list is rebuilt
#' from (`readWorkflowOmicsDaResults`).
#' @noRd
omicsDaContrastOperation <- function(binding) {
    long_table <- binding$long_table
    newOmicsOperation(binding, paste0(binding$prefix, "_da_contrast"), "whole_matrix",
        parameters = c("formula", "qvalue", "lfc", "contrast", "friendly_name"),
        product_type = paste0(binding$prefix, "-da-contrast-references"), output_class = NULL,
        fixed = list(eBayes_trend = TRUE, eBayes_robust = TRUE),
        kernel = function(object, parameters, store) {
            contrast <- parameters$contrast
            accepted <- list(schema = paste0(binding$prefix, ".da_contrast_reference.v1"), source = store$source,
                contrast = contrast, comparison = parameters$friendly_name, assays = list(), errors = list(),
                qvalue_warnings = list(), products = list())
            assays <- methods::slot(object, binding$data_slot)
            design <- object@design_matrix
            preflight <- tryCatch({
                binding$da_preflight(object, data.frame(contrasts = contrast, friendly_names = parameters$friendly_name,
                    stringsAsFactors = FALSE), parameters$formula)
                NULL
            }, error = function(e) conditionMessage(e))
            if (!is.null(preflight)) {
                accepted$errors <- stats::setNames(rep(list(preflight), length(assays)), names(assays))
                return(accepted)
            }
            products <- list()
            for (assay in names(assays)) {
                input <- binding$da_assay_input(assays[[assay]], assay, design, object@sample_id,
                    methods::slot(object, binding$id_slot), object@annotation_id_column)
                if (is.null(input)) next
                fit <- omicsDaAssayFit(binding, store, input$expr_matrix, design, parameters$formula, object@sample_id, assay)
                if (is.null(fit)) next
                outcome <- tryCatch({
                    applied <- binding$da_apply(fit$fit, contrast, if (parameters$lfc > 0) parameters$lfc else NA, TRUE, TRUE)
                    list(applied = applied, table = binding$da_format(applied$results[[contrast]], contrast, assay,
                        parameters$friendly_name, input$annotation_info, input$expr_matrix, design, object@sample_id,
                        object@group_id, parameters$qvalue, parameters$lfc))
                }, error = function(e) {
                    logger::log_error(sprintf("   ERROR in contrast %s, assay %s: %s", contrast, assay, conditionMessage(e)))
                    list(error = conditionMessage(e))
                })
                if (!is.null(outcome$error)) {
                    accepted$errors[[assay]] <- outcome$error
                    next
                }
                table <- as.data.frame(outcome$table)
                accepted$assays[[assay]] <- list(assay = assay, tables = stats::setNames(list(store$table(table)), long_table),
                    rows = nrow(table), fit = fit$reference, significant_counts = binding$da_count(outcome$table))
                if (length(outcome$applied$qvalue_warnings)) accepted$qvalue_warnings[[assay]] <- outcome$applied$qvalue_warnings
                products[[assay]] <- list(results = outcome$applied$results[[contrast]], table = outcome$table,
                    fit.eb = outcome$applied$fit.eb)
            }
            object <- NULL
            accepted$products$da_contrast <- store$product(list(assays = products), "da_contrast", paste0(binding$prefix, "-da-contrast"))
            accepted$rows <- sum(vapply(accepted$assays, function(entry) entry$rows, integer(1)))
            accepted$columns <- if (length(products)) names(products[[1L]]$table)
            accepted
        })
}

#' The contrast strings and friendly names of a contrasts table, as the orchestrators read them.
#' @noRd
omicsDaContrastNames <- function(contrasts) {
    strings <- if ("contrasts" %in% colnames(contrasts)) {
        as.character(contrasts$contrasts)
    } else if ("contrast_string" %in% colnames(contrasts)) {
        as.character(contrasts$contrast_string)
    } else {
        stop("contrasts_tbl must have a 'contrasts' or 'contrast_string' column")
    }
    list(strings = strings, friendly = if ("friendly_names" %in% colnames(contrasts)) as.character(contrasts$friendly_names) else strings)
}

#' The kernel's list rebuilt from a v2 manifest (KIT-076), in the kernel's order: assay-major, then contrast. The
#' exports, the Glimma pages, the volcanoes, the heatmaps, the selectors and the tables read this list. `fit.eb` here
#' is a coefficient view (one column per contrast): its consumers read the names and the row order; the full
#' per-contrast fit stays in the saved product. `failed_contrasts` names what the run could not compute.
#' @param manifest A `<prefix>.da_results.v2` manifest.
#' @return The list `runLipidsDA` / `runMetabolitesDA` return, without the source object.
#' @noRd
readWorkflowOmicsDaResults <- function(binding, manager, manifest) {
    # A contrast the preflight refused saved no product.
    products <- lapply(manifest$results, function(entry) {
        if (is.null(entry$products$da_contrast)) return(list())
        hydrateSavedWorkflowEntity(manager$.context$store, entry$products$da_contrast)$value$assays
    })
    per_assay <- list()
    contrasts_results <- list()
    warnings <- list()
    for (assay in manifest$assays) {
        present <- Filter(Negate(is.null), lapply(products, `[[`, assay))
        if (!length(present)) next
        per_assay[[assay]] <- dplyr::bind_rows(lapply(present, `[[`, "table"))
        contrasts_results[[assay]] <- list(results = lapply(present, `[[`, "results"),
            fit.eb = list(coefficients = do.call(cbind, lapply(unname(present), function(entry) entry$fit.eb$coefficients))),
            qvalue_warnings = do.call(c, c(list(list()), lapply(unname(manifest$results), function(entry) entry$qvalue_warnings[[assay]]))))
        if (length(contrasts_results[[assay]]$qvalue_warnings)) warnings[[assay]] <- contrasts_results[[assay]]$qvalue_warnings
    }
    failed <- Filter(length, lapply(manifest$results, `[[`, "errors"))
    c(list(contrasts_results = contrasts_results), stats::setNames(list(dplyr::bind_rows(per_assay)), binding$long_table),
        list(per_assay_results = per_assay, significant_counts = lapply(per_assay, binding$da_count), qvalue_warnings = warnings,
            da_q_val_thresh = manifest$parameters$qvalue, treat_lfc_cutoff = manifest$parameters$lfc, failed_contrasts = failed))
}

#' Every declaration of one omic, by name. `<prefix>_da_output` is no longer run (KIT-076): it stays registered so v1
#' records and manifests keep resolving their operation.
#' @noRd
omicsWorkflowOperations <- function(binding) {
    prefix <- binding$prefix
    declarations <- list(omicsIntensityFilterOperation(binding), omicsDuplicateResolutionOperation(binding),
        omicsStoreOperation(binding, paste0(prefix, "_qc_finalize")), omicsItsdNormalisationOperation(binding),
        omicsLog2Operation(binding), omicsBetweenSampleNormalisationOperation(binding), omicsRuvOperation(binding),
        omicsPearsonPairsOperation(binding), omicsPearsonSelectionOperation(binding),
        omicsStoreOperation(binding, paste0(prefix, "_norm_complete")), omicsStoreOperation(binding, paste0(prefix, "_reset")),
        omicsDaOperation(binding), omicsDaContrastOperation(binding))
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
#'   value, and the committed object only when `hydrate` asks for it (KIT-075: the modules hold the reference and read
#'   per action; the argument stays for the harness's recorders).
#' @noRd
runWorkflowOmicsReferenceOperation <- function(manager, declaration, parameters, state_name, config, description, source = NULL,
    dependencies = list(), hydrate = FALSE) {
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

#' Read an omics state for one action (KIT-075): the binding's components of the state, read in the owner's scope.
#' The caller uses the object within the action and holds only references between actions.
#' @param reference The state's saved reference or name; the head when NULL.
#' @return The transient object, or NULL when the owner has no such state.
#' @noRd
readWorkflowOmicsState <- function(binding, manager, reference = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort(paste(binding$step_scope, "needs the session's state owner"), class = "saved_entity_scope")
    if (is.character(reference)) reference <- getWorkflowStateReference(manager, reference)
    reference <- reference %||% getWorkflowStateReference(manager)
    if (is.null(reference)) return(NULL)
    withWorkflowStateScope(manager, function() readWorkflowExecutionComponents(manager, reference, fields = binding$fields))
}

#' An omics state's metadata (KIT-075): the stored slots, the assay names and row counts (from the stored tables'
#' shapes) and the design (a sample table), read without any assay table.
#' @param reference The state's saved reference or name; the head when NULL.
#' @return `list(reference, class, slots, assay_names, assay_rows, design_matrix)`, or NULL when the owner has no such
#'   state.
#' @noRd
readWorkflowOmicsMetadata <- function(binding, manager, reference = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) return(NULL)
    if (is.character(reference)) reference <- getWorkflowStateReference(manager, reference)
    reference <- reference %||% getWorkflowStateReference(manager)
    if (is.null(reference)) return(NULL)
    withWorkflowStateScope(manager, function() {
        spec <- readWorkflowStateMetadata(manager, reference)$data
        assays <- spec$components[[binding$data_slot]]
        list(reference = reference, class = spec$class, slots = spec$slots, assay_names = names(assays),
            assay_rows = vapply(assays, function(table) as.integer(table$metadata$shape[[1L]]), integer(1)),
            design_matrix = materializeWorkflowTable(manager$.context, spec$components$design_matrix))
    })
}

#' Bind the omics session owner to its route at the import (KIT-075, the proteomics pattern): binding replaces the
#' pending store, so the import's tables can only be held by reference once the owner is bound. An owner that is
#' already bound, or a bus without an owner, is left as it is.
#' @noRd
bindWorkflowOmicsOwnerAtImport <- function(workflowData, workflow_type, format) {
    if (isDIASessionOwnerPending(workflowData$state_manager)) bindWorkflowOwnerTypeAtDesign(workflowData, workflow_type, format)
    invisible(workflowData$state_manager)
}

#' Hold an import or design table list on the bus by reference, in the omic's owner (KIT-075).
#' @noRd
setWorkflowOmicsBusTables <- function(workflowData, field, value) {
    setWorkflowOwnedPayload(workflowData, field, value, manager = workflowData$state_manager)
}

#' Read an import or design table list from the bus for one action (KIT-075).
#' @noRd
getWorkflowOmicsBusTables <- function(workflowData, field) {
    getWorkflowOwnedPayload(workflowData, field, manager = workflowData$state_manager)
}

#' A design builder's table-list holder (KIT-075): the get/set contract of a `reactiveVal`, holding a reference in the
#' owner `owner()` returns at each use (the builder starts before the import binds the owner); plain without an owner.
#' @noRd
omicsBuilderTableVal <- function(label, owner) {
    held <- shiny::reactiveVal()
    function(value) {
        if (missing(value)) return(resolveWorkflowOwnedValue(held(), owner()))
        held(saveWorkflowOwnedValue(value, label, owner()))
        invisible(NULL)
    }
}

#' One tab step, as a wrapper runs it in process and a managed action dispatches it (KIT-078): the declaration's name, its
#' parameters, and the module's state name and description for the commit.
#' @return `list(operation, parameters, state_name, description)`.
#' @noRd
omicsIntensityFilterStep <- function(binding, percentile, proportion, description = NULL, state_name = NULL) {
    list(operation = paste0(binding$prefix, "_intensity_filter"),
        parameters = stats::setNames(list(percentile, proportion), binding$intensity_parameters),
        state_name = state_name %||% paste0(binding$prefix, "_intensity_filtered"),
        description = description %||% sprintf("Applied %s intensity filter (percentile: %d%%, proportion: %.2f)",
            binding$feature_noun, as.integer(percentile), proportion))
}

#' @noRd
omicsDuplicateResolutionStep <- function(binding, description = NULL, state_name = NULL) {
    list(operation = paste0(binding$prefix, "_duplicate_resolution"), parameters = list(),
        state_name = state_name %||% paste0(binding$prefix, "_duplicates_resolved"),
        description = description %||% sprintf("Resolved duplicate %s features by keeping highest intensity", binding$feature_noun))
}

#' @noRd
omicsItsdNormalisationStep <- function(binding, aggregation, feature_ids, description = NULL, state_name = NULL) {
    parameters <- list(itsd_aggregation = aggregation)
    parameters["itsd_feature_ids"] <- list(feature_ids)
    list(operation = paste0(binding$prefix, "_itsd_normalisation"), parameters = parameters,
        state_name = state_name %||% paste0(binding$prefix, "_itsd_norm"),
        description = description %||% paste("ITSD normalization (aggregation:", aggregation, ")"))
}

#' @noRd
omicsLog2Step <- function(binding, offset, description = NULL, state_name = NULL) {
    list(operation = paste0(binding$prefix, "_log2"), parameters = list(offset = offset),
        state_name = state_name %||% paste0(binding$prefix, "_log2"),
        description = description %||% paste("Log2 transformation (offset:", offset, ")"))
}

#' @noRd
omicsBetweenSampleNormalisationStep <- function(binding, method, description = NULL, state_name = NULL) {
    list(operation = paste0(binding$prefix, "_between_sample_normalisation"), parameters = list(normalisation_method = method),
        state_name = state_name %||% paste0(binding$prefix, "_normalized"),
        description = description %||% paste("Between-sample normalization (method:", method, ")"))
}

#' @param parameters The nine RUV parameters by their declared names (the tab's `ruv_params` plus the mode and grouping).
#' @noRd
omicsRuvStep <- function(binding, parameters, description = NULL, state_name = NULL) {
    list(operation = paste0(binding$prefix, "_ruv"), parameters = parameters[omicsRuvParameterNames],
        state_name = state_name %||% paste0(binding$prefix, "_ruv_corrected"),
        description = description %||% "RUV-III batch correction complete")
}

#' @noRd
omicsPearsonSelectionStep <- function(binding, pairs, threshold, description = NULL, state_name = NULL) {
    list(operation = paste0(binding$prefix, "_pearson_selection"),
        parameters = list(pairs = pairs, min_pearson_correlation_threshold = threshold),
        state_name = state_name %||% paste0(binding$prefix, "_correlation_filtered"),
        description = description %||% paste("Correlation filtering (threshold:", threshold, ")"))
}

#' Run one tab step through its wrapper's commit.
#' @noRd
runWorkflowOmicsReferenceStep <- function(manager, step, config, source = NULL, dependencies = list(), hydrate = FALSE) {
    runWorkflowOmicsReferenceOperation(manager, workflowOperationEntry(step$operation)$declaration, step$parameters, step$state_name,
        config, step$description, source = source, dependencies = dependencies, hydrate = hydrate)
}

#' @noRd
runWorkflowOmicsReferenceIntensityFilter <- function(binding, manager, percentile, proportion, config, description = NULL,
    state_name = NULL, hydrate = FALSE) {
    runWorkflowOmicsReferenceStep(manager, omicsIntensityFilterStep(binding, percentile, proportion, description, state_name), config, hydrate = hydrate)
}

#' @noRd
runWorkflowOmicsReferenceDuplicateResolution <- function(binding, manager, config, description = NULL, state_name = NULL,
    hydrate = FALSE) {
    runWorkflowOmicsReferenceStep(manager, omicsDuplicateResolutionStep(binding, description, state_name), config, hydrate = hydrate)
}

#' @noRd
runWorkflowOmicsReferenceItsdNormalisation <- function(binding, manager, aggregation, feature_ids, config, description = NULL,
    state_name = NULL, hydrate = FALSE) {
    runWorkflowOmicsReferenceStep(manager, omicsItsdNormalisationStep(binding, aggregation, feature_ids, description, state_name), config, hydrate = hydrate)
}

#' @noRd
runWorkflowOmicsReferenceLog2 <- function(binding, manager, offset, config, description = NULL, state_name = NULL, hydrate = FALSE) {
    runWorkflowOmicsReferenceStep(manager, omicsLog2Step(binding, offset, description, state_name), config, hydrate = hydrate)
}

#' @noRd
runWorkflowOmicsReferenceBetweenSampleNormalisation <- function(binding, manager, method, config, description = NULL,
    state_name = NULL, hydrate = FALSE) {
    runWorkflowOmicsReferenceStep(manager, omicsBetweenSampleNormalisationStep(binding, method, description, state_name), config, hydrate = hydrate)
}

#' @param parameters The nine RUV parameters by their declared names (the tab's `ruv_params` plus the mode and grouping).
#' @return The operation result with `optimisation`, the hydrated per-assay optimisation results.
#' @noRd
runWorkflowOmicsReferenceRuv <- function(binding, manager, parameters, config, description = NULL, state_name = NULL,
    hydrate = FALSE) {
    result <- runWorkflowOmicsReferenceStep(manager, omicsRuvStep(binding, parameters, description, state_name), config, hydrate = hydrate)
    result$optimisation <- hydrateSavedWorkflowEntity(manager$.context$store, result$value$optimisation)$value
    result
}

#' The correlation filter: the pairs product, then the selection on it (the DIA two-operation pattern).
#' @return The selection's result with `correlations`, the per-assay correlation tables.
#' @noRd
runWorkflowOmicsReferenceCorrelation <- function(binding, manager, group, threshold, config, description = NULL,
    state_name = NULL, hydrate = FALSE) {
    if (!inherits(manager, "DIAStateAdapter")) {
        cli::cli_abort(paste(binding$step_scope, "needs the session's state owner"), class = "saved_entity_scope")
    }
    withWorkflowStateScope(manager, function() {
        source <- getWorkflowStateReference(manager)
        pairs <- runSavedOperation(manager, omicsPearsonPairsOperation(binding), source, list(correlation_group = group))
        result <- runWorkflowOmicsReferenceStep(manager, omicsPearsonSelectionStep(binding, pairs$reference$entity_id, threshold, description, state_name),
            config, source = source, dependencies = list(pairs$reference), hydrate = hydrate)
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

#' Commit a store step on the owner (the QC finalise, the correlation skip, the reset): the head's record, or the
#' `source` state's, under the new name with the step's provenance, the prior products kept and no component re-serialised.
#' @noRd
commitWorkflowOmicsStoreStep <- function(manager, operation, state_name, config, description, source = NULL) {
    if (!inherits(manager, "DIAStateAdapter")) cli::cli_abort("the step needs the session's state owner", class = "saved_entity_scope")
    withWorkflowStateScope(manager, function() {
        markWorkflowOmicsStoreStep(manager, operation)
        # KIT-075: `source` (a saved state reference) recommits that state's record, as the reset does; the head's otherwise.
        commitWorkflowReferenceState(manager, state_name, readWorkflowStateMetadata(manager, source)$data, config, description)
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
    if (length(paths)) recordWorkflowInputs(manager, paths)
    invisible(NULL)
}

#' Run the omic's DA on a state, one declaration per contrast on the saved per-assay fits, and accept the manifest on
#' the head (the `runDIAReferenceDaAnalysis` twin; KIT-076, `<prefix>.da_results.v2`, keyed by contrast), the DA step's
#' provenance cleared (it commits no state). A failed contrast is recorded and the others complete; when none
#' completes the run stops with their messages and the head keeps its previous result.
#' @param manager The session owner.
#' @param source Saved reference of the state the DA runs on.
#' @param settings `formula_string`, `da_q_val_thresh`, `treat_lfc_cutoff`.
#' @param contrasts The contrasts table (`contrasts`, `friendly_names`, `full_format`).
#' @return `list(manifest, results, reference)`: `results` is the kernel's list for the module's display and files
#'   (`readWorkflowOmicsDaResults`; without the source object, whose readers take its components per action).
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
        metadata <- readWorkflowOmicsMetadata(binding, manager, source)
        contrasts <- as.data.frame(contrasts, stringsAsFactors = FALSE)
        names_of <- omicsDaContrastNames(contrasts)
        parameters <- list(formula = settings$formula_string, qvalue = settings$da_q_val_thresh, lfc = settings$treat_lfc_cutoff)
        operation <- omicsDaContrastOperation(binding)
        results <- list()
        for (i in seq_along(names_of$strings)) {
            run <- runSavedOperation(manager, operation, source, c(parameters,
                list(contrast = names_of$strings[[i]], friendly_name = names_of$friendly[[i]])))
            if (identical(run$status, "computed")) {
                snapshot <- context$snapshot
                snapshot$da_dependencies[[run$reference$entity_id]] <- collectWorkflowNativeReferences(run$value)
                publishWorkflowStateHead(context, snapshot)
            }
            results[[names_of$strings[[i]]]] <- run$value
        }
        clearWorkflowStepProvenance(manager)
        if (length(results) && all(vapply(results, function(entry) !length(entry$assays) && length(entry$errors) > 0L, logical(1)))) {
            stop(paste(unique(unlist(lapply(results, `[[`, "errors"))), collapse = "; "), call. = FALSE)
        }
        manifest <- c(list(schema = paste0(binding$prefix, ".da_results.v2"), source = source),
            stats::setNames(list(metadata$slots[[binding$id_slot]]), binding$id_slot),
            list(sample_id = metadata$slots$sample_id, assays = metadata$assay_names, contrasts = contrasts, results = results,
                parameters = parameters, model = list(model = "limma", formula_used = parameters$formula, eBayes_trend = TRUE,
                    eBayes_robust = TRUE, treat_applied = isTRUE(parameters$lfc > 0))))
        reference <- saveWorkflowStateValue(context, manifest, "da-result-references", codec = "rds")
        snapshot <- context$snapshot
        snapshot$da_results <- reference
        snapshot$da_dependencies[[reference$entity_id]] <- collectWorkflowNativeReferences(list(source, results))
        publishWorkflowStateHead(context, snapshot)
        list(manifest = manifest, results = readWorkflowOmicsDaResults(binding, manager, manifest), reference = reference)
    })
}

# ---------------------------------------------------------------------------------------------------------------
# The worker (KIT-078): each action runs under `computeWorkflowOmicsAction`, in the owner's scope, and leaves the head
# to the controller's completion (`completeWorkflowOmicsAction`). Nothing here reads a Shiny session.
# ---------------------------------------------------------------------------------------------------------------

#' A directory the worker stages one role's files in; the controller publishes them (`publishWorkflowOmicsResult`).
#' @noRd
omicsStagedPublicationDir <- function(request, role) {
    path <- file.path(request$directory, "omics-publication", role)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
}

#' The filtering step's metrics, saved for the controller's figures (the metabolite tables can pass the result's limit).
#' A failure costs the figures only, as it does in process.
#' @return The saved value's reference, or NULL.
#' @noRd
omicsFilteringMetricsReference <- function(binding, manager, state) {
    tryCatch({
        metrics <- binding$filtering_metrics(readWorkflowExecutionComponents(manager, state, fields = binding$fields))
        saveWorkflowStateValue(manager$.context, metrics, "omics-filtering-metrics", codec = "rds")
    }, error = function(error) {
        logger::log_warn(paste("Could not compute the filtering metrics:", conditionMessage(error)))
        NULL
    })
}

#' The worker's one compute handler: the registered declaration on the request's state, by the request's variant. A
#' state-producing declaration returns the checkpoint the completion commits; a product-only one returns the source's.
#' @noRd
computeWorkflowOmicsOperationAction <- function(manager, request, variant) {
    entry <- workflowActionSpec(request$operation)
    declaration <- entry$declaration
    binding <- omicsDeclarationBindings()[[entry$omic_type]]
    source <- request$input_references$state
    parameters <- request$parameters
    parameters$variant <- variant
    dependencies <- if ("pairs" %in% declaration$parameters) list(manager$.context$snapshot$products[[parameters$pairs]]) else list()
    run <- runSavedOperation(manager, declaration, source, parameters, dependencies)
    # An earlier worker's product whose action never completed has no root yet.
    if (identical(run$status, "hit")) rootWorkflowProduct(manager, run$reference, c(list(source), dependencies, list(run$value)))
    if (is.null(declaration$output_class)) {
        return(list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = run$reference,
            native = if (identical(declaration$name, paste0(binding$prefix, "_da_contrast"))) collectWorkflowNativeReferences(run$value)))
    }
    final <- readWorkflowStateMetadata(manager, run$value$final)$data
    final$products[[paste0(declaration$name, "_science")]] <- run$reference
    filtering <- declaration$name %in% paste0(binding$prefix, c("_intensity_filter", "_duplicate_resolution"))
    list(checkpoint = final, reference = run$reference,
        display = if (filtering) list(filtering_metrics = omicsFilteringMetricsReference(binding, manager, run$value$final)))
}

#' One stage's per-assay QC figures from the request's state: the files staged, the assets saved.
#' @noRd
computeWorkflowOmicsQcFiguresAction <- function(manager, request, variant) {
    binding <- omicsDeclarationBindings()[[workflowActionSpec(request$operation)$omic_type]]
    source <- request$input_references$state
    binding$qc_plots(theObject = readWorkflowExecutionComponents(manager, source, fields = binding$fields),
        experiment_paths = stats::setNames(list(omicsStagedPublicationDir(request, binding$qc_dir_role)), binding$qc_dir_role),
        stage = request$parameters$stage, grouping_variable = request$parameters$grouping_variable,
        shape_variable = request$parameters$shape_variable)
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = source)
}

#' The DA tables and figures from the accepted manifest and its source state: the files staged, the assets saved. A
#' failure is reported to the tab, which warns as it does in process; the science is already accepted.
#' @noRd
computeWorkflowOmicsDaPublicationAction <- function(manager, request, variant) {
    binding <- omicsDeclarationBindings()[[workflowActionSpec(request$operation)$omic_type]]
    source <- request$input_references$state
    manifest_reference <- manager$.context$snapshot$da_results
    if (is.null(manifest_reference)) cli::cli_abort("no accepted DA manifest", class = "workflow_execution_invalid")
    manifest <- resolveWorkflowStateValue(manager$.context, manifest_reference)
    results <- readWorkflowOmicsDaResults(binding, manager, manifest)
    results$theObject <- readWorkflowOmicsState(binding, manager, manifest$source)
    export <- tryCatch(list(success = isTRUE(binding$da_export(da_results_list = results,
        da_output_dir = omicsStagedPublicationDir(request, "da_output_dir"),
        publication_graphs_dir = omicsStagedPublicationDir(request, "publication_graphs_dir"),
        da_q_val_thresh = manifest$parameters$qvalue, lfc_threshold = manifest$parameters$lfc, heatmap_top_n = 50,
        heatmap_clustering = "both", heatmap_color_scheme = "RdBu"))), error = function(error) {
            list(success = FALSE, message = substr(conditionMessage(error), 1L, 1024L))
        })
    list(checkpoint = readWorkflowStateMetadata(manager, source)$data, reference = source, display = export)
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
