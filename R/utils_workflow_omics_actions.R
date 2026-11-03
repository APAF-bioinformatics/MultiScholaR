# KIT-078: the lipidomics and metabolomics tabs on the managed worker. One runner per omic per session drives a tab's
# actions in order over `installWorkflowManagedStep`, the proteomics modules' state machine. Under a managed application
# (`workflowSessionApplication`) a tab asks for the runner and dispatches its declarations; with no runner (non-Linux, no
# managed application, a test) the same tab runs the same declarations in process through their wrappers.

#' Install the omic's runner on a managed session. An unmanaged session gets none.
#' @param session The omic's module session.
#' @param workflow_data The omic's bus.
#' @param omic_type `lipidomics` or `metabolomics`.
#' @param experiment_label,experiment_paths The publication context of every action.
#' @return The runner (also kept on the session for `workflowOmicsManagedActions`), or NULL.
#' @noRd
installWorkflowOmicsActions <- function(session, workflow_data, omic_type, experiment_label, experiment_paths) {
    if (is.null(workflowSessionApplication(session))) return(NULL)
    runner <- new.env(parent = emptyenv())
    run <- new.env(parent = emptyenv())
    start <- function() {
        action <- run$actions[[run$index]]
        if (is.function(action$before)) action$before()
        step$apply(action$parameters, action$operation, action$context)
    }
    begin <- function(actions, done, failed) {
        run$actions <- actions
        run$index <- 1L
        run$done <- done
        run$failed <- failed
        start()
    }
    delivered <- function() {
        run$index <- run$index + 1L
        if (run$index <= length(run$actions)) return(start())
        done <- run$done
        run$actions <- NULL
        if (is.function(done)) done()
        # A tab's own sequence (`done` starting the next) comes before a click that waited.
        queued <- run$queued
        if (!is.null(queued) && is.null(shiny::isolate(step$active()))) {
            run$queued <- NULL
            begin(queued$actions, queued$done, queued$failed)
        }
    }
    step <- installWorkflowManagedStep(session, NULL, workflow_data, NULL,
        publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
        ui = list(notification = paste0(omic_type, "_action_working"), label = paste(omic_type, "step"), delivered = delivered,
            failed = function(error) {
                run$queued <- NULL
                if (is.function(run$failed)) run$failed(error)
            }),
        deliver = function(result) shiny::isolate(run$actions[[run$index]]$deliver(result)),
        read_result = readManagedOmicsResult)
    runner$step <- step
    #' Run a tab's actions in order: `done()` after the last, `failed(error)` if one fails (the step has reported it).
    runner$run <- function(actions, done = NULL, failed = NULL) {
        # A publication that failed is retried from the action it stopped at; its science is never recomputed.
        if (!is.null(shiny::isolate(step$pending())) && length(run$actions)) return(start())
        # The tab's button can be pressed while another tab's action runs (the pre-normalisation figures, opening the tab).
        if (!is.null(shiny::isolate(step$active()))) {
            run$queued <- list(actions = actions, done = done, failed = failed)
            return(invisible(NULL))
        }
        begin(actions, done, failed)
    }
    runner$progress <- function(message) {
        progress <- shiny::Progress$new(session, min = 0, max = 1)
        progress$set(message = message, value = 0)
        progress
    }
    session$userData$workflow_omics_actions[[omic_type]] <- runner
    runner
}

#' The omic's runner on the current session, or NULL when the tab runs in process.
#' @noRd
workflowOmicsManagedActions <- function(omic_type) {
    session <- shiny::getDefaultReactiveDomain()
    if (is.null(session)) return(NULL)
    session$userData$workflow_omics_actions[[omic_type]]
}

#' A state-producing tab step as a managed action: the declaration's parameters and the module's commit.
#' @param step `omics*Step`.
#' @param config The state's configuration.
#' @param deliver `function(result)` with `readManagedOmicsResult`'s value, in the owner's scope.
#' @param before Optional function run just before the action is dispatched.
#' @noRd
omicsManagedCommit <- function(step, config, deliver, before = NULL) {
    list(operation = step$operation, parameters = step$parameters, before = before, deliver = deliver,
        context = list(state_name = step$state_name, config = config, description = step$description))
}

#' A managed action that commits nothing: a product-only declaration, the figures, the DA export.
#' @noRd
omicsManagedProduct <- function(operation, parameters, deliver, before = NULL) {
    list(operation = operation, parameters = parameters, before = before, deliver = deliver, context = NULL)
}

#' The filtering metrics a managed step's worker saved, for the figures the tab draws (`updateLipidFiltering`,
#' `updateMetaboliteFiltering`).
#' @return The metrics, or NULL when the worker could not compute them.
#' @noRd
readWorkflowOmicsFilteringMetrics <- function(manager, result) {
    reference <- result$display$filtering_metrics
    if (is.null(reference)) return(NULL)
    resolveWorkflowStateValue(manager$.context, reference)
}

#' The ITSD selection's feature ids from the candidates the tab's tables showed: the table's row indices per assay, read
#' without the matrix (`itsd_candidate_ids`, kept when the tables are built).
#' @return Feature ids per assay with a selection, or NULL with none.
#' @noRd
omicsManualItsdFeatureIds <- function(selections, candidates, addLog = function(entry) invisible(entry)) {
    selected <- list()
    for (assay in names(selections)) {
        rows <- selections[[assay]]
        if (is.null(rows) || !length(rows) || is.null(candidates[[assay]])) next
        selected[[assay]] <- candidates[[assay]][rows]
        addLog(paste("Assay", assay, ":", length(selected[[assay]]), "ITSD features selected"))
    }
    if (length(selected)) selected
}

#' The normalisation tab's pipeline on the worker (`handleLipidNormRunNormalization` and `runMetabNormNormalizationPipeline`
#' run it in process): the QC figures, ITSD, log2, between-sample normalisation, the figures, RUV-III and its figures, the
#' composite, one managed action each under one progress display. `normData` holds each step's state reference as in process.
#' @param binding The omic's binding.
#' @param runner The omic's runner.
#' @param manager The session owner.
#' @param settings `apply_itsd`, `itsd_aggregation`, `itsd_feature_ids`, `log_offset`, `norm_method`, `ruv_mode`,
#'   `ruv_parameters` (the nine by their declared names), `grouping_variable`, `shape_variable`, `assay_names`, `qc_dir`.
#' @param normData The tab's reactive values.
#' @param config The states' configuration.
#' @param addLog The tab's log.
#' @param finish Function run after the last action: the composite figure and the tab's completion.
#' @param failed `function(error)`, the tab's handling of a failure.
#' @noRd
runWorkflowOmicsManagedNormalisation <- function(binding, runner, manager, settings, normData, config, addLog, finish, failed) {
    prefix <- binding$prefix
    progress <- runner$progress("Running normalization pipeline...")
    total <- 6
    advance <- function(detail) progress$inc(1 / total, detail = detail)
    figures <- function(stage, done, before = NULL) omicsManagedProduct(paste0(prefix, "_qc_figures"),
        list(stage = stage, grouping_variable = settings$grouping_variable, shape_variable = settings$shape_variable),
        function(result) addLog(done), before)
    ruv <- !identical(settings$ruv_mode, "skip")
    actions <- list(figures("post_filter", "Pre-normalization QC plots generated", function() {
        advance("Capturing pre-normalization state...")
        addLog("Post-filtering state captured")
    }))
    if (isTRUE(settings$apply_itsd)) {
        actions <- c(actions, list(omicsManagedCommit(omicsItsdNormalisationStep(binding, settings$itsd_aggregation, settings$itsd_feature_ids), config,
            function(result) {
                normData$post_itsd_obj <- result$stateReference
                addLog("ITSD normalization complete")
            }, function() {
                advance("Applying ITSD normalization...")
                addLog(paste("Applying ITSD normalization (aggregation:", settings$itsd_aggregation, ")"))
            })))
    }
    actions <- c(actions, list(
        omicsManagedCommit(omicsLog2Step(binding, settings$log_offset), config, function(result) {
            normData$post_log2_obj <- result$stateReference
            addLog("Log2 transformation complete")
        }, function() {
            if (!isTRUE(settings$apply_itsd)) {
                advance("Applying ITSD normalization...")
                addLog("ITSD normalization skipped")
            }
            advance("Applying log2 transformation...")
            addLog(paste("Applying log2 transformation (offset:", settings$log_offset, ")"))
        }),
        omicsManagedCommit(omicsBetweenSampleNormalisationStep(binding, settings$norm_method), config, function(result) {
            normData$post_norm_obj <- result$stateReference
            addLog("Between-sample normalization complete")
            normData$normalization_complete <- TRUE
        }, function() {
            advance("Applying between-sample normalization...")
            if (!identical(settings$norm_method, "none")) addLog(paste("Applying between-sample normalization (method:", settings$norm_method, ")"))
        }),
        figures("post_norm", "Post-normalization QC plots generated")))
    if (ruv) {
        actions <- c(actions, list(omicsManagedCommit(omicsRuvStep(binding, settings$ruv_parameters), config, function(result) {
            normData$ruv_optimization_results <- hydrateSavedWorkflowEntity(manager$.context$store, result$value$optimisation)$value
            addLog(paste("RUV optimization complete. Best k per assay:",
                paste(names(result$value$best_k), "=", unlist(result$value$best_k), collapse = ", ")))
            normData$ruv_corrected_obj <- result$stateReference
            normData$ruv_complete <- TRUE
            addLog("RUV-III correction applied")
        }, function() {
            advance("Running RUV-III batch correction...")
            addLog(paste("Running RUV-III (mode:", settings$ruv_mode, ")"))
        }), figures("ruv_corrected", "RUV QC plots generated", function() advance("Generating RUV QC plots..."))))
    }
    warn_composite <- function(message) {
        addLog(paste("Warning: Could not generate composite QC figure:", message))
        logger::log_warn(paste("Composite QC generation failed:", message))
    }
    runner$run(actions, done = function() {
        if (!ruv) {
            addLog("RUV-III skipped")
            normData$ruv_corrected_obj <- normData$post_norm_obj
            normData$ruv_complete <- TRUE
        }
        close_and_finish <- function() {
            progress$close()
            finish()
        }
        if (is.null(settings$assay_names)) return(close_and_finish())
        addLog("Generating composite QC figure...")
        # The states are committed; the figure's failure warns, as it does in process.
        runner$run(list(omicsManagedProduct(paste0(prefix, "_composite_figure"), list(assay_names = settings$assay_names, ruv_mode = settings$ruv_mode),
            function(result) {
                if (isTRUE(result$display$saved)) addLog(sprintf("Composite QC figure saved to: %s", file.path(settings$qc_dir, "composite_QC_figure")))
                else if (!is.null(result$display$message)) warn_composite(result$display$message)
            })), done = close_and_finish, failed = function(error) {
                warn_composite(conditionMessage(error))
                close_and_finish()
            })
    }, failed = function(error) {
        progress$close()
        failed(error)
    })
    invisible(TRUE)
}

#' The correlation filter on the worker: the pairs, then the selection on them (`runWorkflowOmicsReferenceCorrelation`).
#' @param deliver `function(correlations, reference)`, the tab's completion.
#' @param failed `function(error)`, the tab's handling of a failure.
#' @noRd
runWorkflowOmicsManagedCorrelation <- function(binding, runner, group, threshold, config, deliver, failed) {
    runner$run(list(omicsManagedProduct(paste0(binding$prefix, "_pearson_pairs"), list(correlation_group = group), function(result) {
        runner$correlation_pairs <- result
    })), done = function() {
        pairs <- runner$correlation_pairs
        runner$correlation_pairs <- NULL
        runner$run(list(omicsManagedCommit(omicsPearsonSelectionStep(binding, pairs$reference$entity_id, threshold), config,
            function(result) deliver(pairs$value, result$stateReference))), failed = failed)
    }, failed = failed)
    invisible(TRUE)
}

#' The DA on the worker: one action per contrast on the shared fits, the manifest accepted on the head, then the tables and
#' figures in one action (`runWorkflowOmicsReferenceDaAnalysis`, and the tab's export, in process).
#' @param manager The session owner.
#' @param source The state the DA runs on.
#' @param settings `formula_string`, `da_q_val_thresh`, `treat_lfc_cutoff`.
#' @param contrasts The contrasts table.
#' @param deliver `function(results, export)`: the kernel's list and the export's outcome (`list(success, message)`).
#' @param failed `function(error)` for a failure before the results exist.
#' @noRd
runWorkflowOmicsManagedDa <- function(binding, runner, manager, source, settings, contrasts, deliver, failed) {
    # A worker's input is the head; a tab left on an earlier state reloads it (in process the DA could run on either).
    if (!identical(source, getWorkflowStateReference(manager))) {
        stop("the session has moved on from the state this tab loaded; load the session again", call. = FALSE)
    }
    values <- list()
    actions <- lapply(omicsDaContrastSteps(binding, settings, contrasts), function(step) {
        omicsManagedProduct(step$operation, step$parameters, function(result) values[[step$contrast]] <<- result$value)
    })
    runner$run(actions, done = function() {
        accepted <- withWorkflowStateScope(manager, function() acceptWorkflowOmicsDaManifest(binding, manager, source, settings, contrasts, values))
        export <- NULL
        finish <- function(outcome) deliver(accepted$results, outcome)
        runner$run(list(omicsManagedProduct(paste0(binding$prefix, "_da_publication"), list(), function(result) export <<- result$display)),
            done = function() finish(export),
            # The science is accepted; the tab warns about the files, as it does in process.
            failed = function(error) finish(list(success = FALSE, message = conditionMessage(error))))
    }, failed = failed)
    invisible(TRUE)
}

#' The tab's report of a worker's DA export, in the words of the in-process write (`list(status, success, message)`).
#' @noRd
reportWorkflowOmicsDaExport <- function(export, notify = shiny::showNotification, logInfo = logger::log_info, logWarn = logger::log_warn) {
    if (isTRUE(export$success)) {
        logInfo("   All DA results written to disk successfully")
        notify("DA results saved to disk (tables, volcano plots, heatmaps)", type = "message", duration = 5)
        return(list(status = "completed", success = TRUE))
    }
    message <- export$message %||% "the export reported no success"
    logWarn(paste("   Could not write DA results to disk:", message))
    notify(paste("Warning: Could not save results to disk:", message), type = "warning", duration = 8)
    list(status = "error", success = FALSE, message = message)
}
