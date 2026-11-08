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
    begin <- function(actions, done, failed, dropped = NULL) {
        run$actions <- actions
        run$index <- 1L
        run$done <- done
        run$failed <- failed
        run$dropped <- dropped
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
        if (length(queued) && is.null(shiny::isolate(step$active()))) {
            run$queued <- queued[-1L]
            begin(queued[[1L]]$actions, queued[[1L]]$done, queued[[1L]]$failed, queued[[1L]]$dropped)
        }
    }
    step <- installWorkflowManagedStep(session, NULL, workflow_data, NULL,
        publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
        ui = list(notification = paste0(omic_type, "_action_working"), label = paste(omic_type, "step"), delivered = delivered,
            failed = function(error) {
                queued <- run$queued
                run$queued <- NULL
                if (is.function(run$failed)) run$failed(error)
                # KIT-082: what was waiting is dropped and told, so that a display does not wait for ever. The failure is not
                # theirs, so a display asks again (`dropped`).
                for (waiting in queued) {
                    if (is.function(waiting$dropped)) waiting$dropped() else if (is.function(waiting$failed)) waiting$failed(error)
                }
            }),
        deliver = function(result) shiny::isolate(run$actions[[run$index]]$deliver(result)),
        read_result = readManagedOmicsResult)
    runner$step <- step
    #' Run a tab's actions in order: `done()` after the last, `failed(error)` if one fails (the step has reported it).
    #' `dropped()`, when given, replaces `failed` if another action's failure drops these before they start.
    runner$run <- function(actions, done = NULL, failed = NULL, dropped = NULL) {
        waiting <- list(list(actions = actions, done = done, failed = failed, dropped = dropped))
        if (runner$retrying()) {
            # A publication that failed is retried from the action it stopped at; its science is never recomputed. The same
            # actions are that retry, and the first request's callbacks finish it. Anything else waits its turn.
            operations <- function(list) vapply(list, function(action) action$operation, character(1))
            if (length(actions) && !identical(operations(actions), operations(run$actions))) run$queued <- c(run$queued, waiting)
            return(start())
        }
        # The tab's button can be pressed while another tab's action runs (the pre-normalisation figures, opening the tab).
        if (!is.null(shiny::isolate(step$active()))) {
            run$queued <- c(run$queued, waiting)
            return(invisible(NULL))
        }
        begin(actions, done, failed, dropped)
    }
    #' A failed publication waiting for its retry.
    runner$retrying <- function() !is.null(shiny::isolate(step$pending())) && length(run$actions) > 0L
    runner$progress <- function(message) {
        # A retry runs under the first request's progress, which its callbacks close.
        if (runner$retrying()) return(list(set = function(...) invisible(NULL), inc = function(...) invisible(NULL), close = function(...) invisible(NULL)))
        progress <- shiny::Progress$new(session, min = 0, max = 1)
        progress$set(message = message, value = 0)
        progress
    }
    runner$displays <- shiny::reactiveValues()
    runner$requested <- character()
    runner$failures <- list()
    session$userData$workflow_omics_actions[[omic_type]] <- runner
    runner
}

#' One display from the worker, asked for once per key by the render that shows it (KIT-082). The render reads the entry, so
#' it runs again when the action delivers it.
#' @param runner The omic's runner.
#' @param key What identifies this display.
#' @param operation,parameters The registered operation.
#' @param read `function(result)` giving what the render keeps.
#' @return The entry, or NULL while the worker has not delivered it. A failure gives `list(message =)`.
#' @noRd
workflowOmicsDisplay <- function(runner, key, operation, parameters, read = function(result) result$display) {
    # A failed or dropped request brings the render back through this.
    runner$displays[[".asked"]]
    entry <- runner$displays[[key]]
    if (!is.null(entry)) return(entry)
    # A failure is shown once: the next time the render runs it asks again, as an in-process render would compute again.
    failure <- runner$failures[[key]]
    if (!is.null(failure)) {
        runner$failures[[key]] <- NULL
        return(list(message = failure))
    }
    # The render comes back when the retry clears it.
    if (!is.null(runner$step$pending())) {
        return(list(message = "An earlier step's files could not be written. Press that step's button again; this display follows."))
    }
    if (key %in% runner$requested) return(NULL)
    runner$requested <- c(runner$requested, key)
    again <- function() shiny::isolate({
        runner$requested <- setdiff(runner$requested, key)
        runner$displays[[".asked"]] <- (runner$displays[[".asked"]] %||% 0L) + 1L
    })
    shiny::isolate(runner$run(list(omicsManagedProduct(operation, parameters, function(result) runner$displays[[key]] <- read(result))),
        failed = function(error) {
            runner$failures[[key]] <- conditionMessage(error)
            again()
        }, dropped = again))
    NULL
}

#' The final-QC PCA from the worker.
#' @return `list(plot =)`, NULL while it is being drawn, or `list(message =)` on a failure.
#' @noRd
workflowOmicsFinalPca <- function(binding, runner, manager, state, color, shape) {
    key <- paste0("final-pca-", digest::digest(list(state, color, shape), algo = "sha256"))
    workflowOmicsDisplay(runner, key, paste0(binding$prefix, "_final_pca"), list(state = state, grouping_variable = color, shape_variable = shape),
        function(result) if (is.null(result$display$plot)) result$display else list(plot = resolveWorkflowStateValue(manager$.context, result$display$plot)))
}

#' A volcano widget generator for a managed session: the tab's Glimma output builder keeps its messages, and this gives it the
#' worker's page in place of the widget.
#' @return A function with the generator's arguments.
#' @noRd
workflowOmicsVolcanoPageFactory <- function(binding, runner, manager) {
    function(da_results_list, selected_contrast, selected_assay, da_q_val_thresh) {
        key <- paste0("volcano-", digest::digest(list(manager$.context$snapshot$da_results$entity_id, selected_contrast, selected_assay, as.numeric(da_q_val_thresh)), algo = "sha256"))
        entry <- workflowOmicsDisplay(runner, key, paste0(binding$prefix, "_volcano_page"),
            list(contrast = selected_contrast, assay = selected_assay, qvalue = as.numeric(da_q_val_thresh)),
            # `list(page = NULL)` would read as an absent entry.
            function(result) if (is.null(result$display$page) && is.null(result$display$message)) list(unavailable = TRUE) else result$display)
        if (is.null(entry)) return(shiny::div(class = "alert alert-info", style = "margin: 20px;", "Preparing the interactive plot..."))
        if (!is.null(entry$message)) stop(entry$message, call. = FALSE)
        if (isTRUE(entry$unavailable)) return(NULL)
        shiny::tags$iframe(srcdoc = rawToChar(hydrateSavedWorkflowEntity(manager$.context$store, entry$page)$bytes),
            style = "width:100%;height:820px;border:0;", title = paste("Interactive volcano", selected_contrast))
    }
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

#' A QC step that commits a state and draws the filtering-progress figure (KIT-104 D6; the intensity filter and the
#' duplicate resolution, mod_lipid_qc_intensity_server_helpers.R:241-260 and mod_lipid_qc_duplicates_server_helpers.R
#' :256-268 at b423b764, the same in metab; those files were deleted by KIT-104 W5). On a runner: the managed commit, and the figure from the worker's metrics. In process:
#' the reference step, and the figure from the committed state read for this action only. The tab holds no object.
#' @param binding The omic's declaration binding.
#' @param workflowData The omic's bus (`state_manager`, `config_list`).
#' @param step `omicsIntensityFilterStep` or `omicsDuplicateResolutionStep`.
#' @param deliver `function(result, figure)`: `result` is `list(value, reference)`, the declaration's value and the
#'   committed state's reference; `figure()` draws the figure, when the tab draws it.
#' @param draw `function(arguments)` drawing the figure from `list(theObject =)` or `list(theObject = NULL, metrics =)`.
#' @param failed `function(error)` for a managed action's failure (the step has reported it). In process an error is
#'   raised to the caller, as the tabs' observers handle it.
#' @noRd
runWorkflowOmicsQcAction <- function(binding, workflowData, step, deliver, draw, failed = NULL) {
    runner <- workflowOmicsManagedActions(binding$omic_type)
    if (!is.null(runner)) {
        runner$run(list(omicsManagedCommit(step, workflowData$config_list, function(result) {
            deliver(list(value = result$value, reference = result$stateReference), function() {
                metrics <- readWorkflowOmicsFilteringMetrics(workflowData$state_manager, result)
                if (is.null(metrics)) stop("the worker computed no filtering metrics")
                draw(list(theObject = NULL, metrics = metrics))
            })
        })), failed = failed)
        return(invisible(NULL))
    }
    run <- runWorkflowOmicsReferenceStep(workflowData$state_manager, step, workflowData$config_list)
    deliver(list(value = run$value, reference = run$reference), function() {
        draw(list(theObject = readWorkflowOmicsState(binding, workflowData$state_manager, run$reference)))
    })
    invisible(NULL)
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
#' @param invalidate `function()` clearing the tab's results: called when the new manifest is on the head but its rows could
#'   not be read, so that the tab shows no results rather than the previous run's beside this run's volcano.
#' @noRd
runWorkflowOmicsManagedDa <- function(binding, runner, manager, source, settings, contrasts, deliver, failed, invalidate = NULL) {
    # A worker's input is the head; a tab left on an earlier state reloads it (in process the DA could run on either).
    if (!identical(source, getWorkflowStateReference(manager))) {
        stop("the session has moved on from the state this tab loaded; load the session again", call. = FALSE)
    }
    values <- list()
    actions <- lapply(omicsDaContrastSteps(binding, settings, contrasts), function(step) {
        omicsManagedProduct(step$operation, step$parameters, function(result) values[[step$contrast]] <<- result$value)
    })
    runner$run(actions, done = function() {
        # KIT-082: the worker reads the contrasts' tables, and the tab holds their rows without the per-sample columns.
        accepted <- withWorkflowStateScope(manager, function() acceptWorkflowOmicsDaManifest(binding, manager, source, settings, contrasts, values, read_results = FALSE))
        tab <- NULL
        export <- NULL
        finish <- function(outcome) {
            results <- resolveWorkflowStateValue(manager$.context, tab)
            columns <- results$download_columns
            results$write_download <- function(file) writeWorkflowOmicsDaDownload(binding, manager, accepted$manifest, columns, file)
            deliver(results, outcome)
        }
        runner$run(list(omicsManagedProduct(paste0(binding$prefix, "_da_tab_results"), list(), function(result) tab <<- result$display$tab),
            omicsManagedProduct(paste0(binding$prefix, "_da_publication"), list(), function(result) export <<- result$display)),
            done = function() finish(export),
            # The science is accepted; with the tab's rows the tab warns about the files, as it does in process.
            failed = function(error) {
                if (!is.null(tab)) return(finish(list(success = FALSE, message = conditionMessage(error))))
                if (is.function(invalidate)) invalidate()
                failed(error)
            })
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
