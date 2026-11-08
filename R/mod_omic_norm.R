# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# KIT-105 W4: the shared normalisation frame (kit-101 design 4.7), from the lipid tab's server at f8d4c5fc
# (mod_lipid_norm_server.R and its helper files, deleted by W4). `mod_omic_norm_ui(id, config)` renders the options
# column and the QC tabset of an omic's `norm` group from the panel library (mod_omic_norm_panels.R);
# `mod_omic_norm_server(...)` holds the tab's reactive state (state references only, KIT-075), its log, the pre-QC on
# tab entry and the run, reset, correlation apply and skip, and export actions. Each action calls the core: the
# managed function with the omic's runner, the in-process twin without one (D6). Functions a test mocks are called
# by name when they run. The module frame keeps no locals (the testServer masking rule; test-omic-design-module.R).
# Recorded divergence from conventions/R-shiny.md: the in-process or managed actions in place of ExtendedTask are
# kept, under the existing-code scope rule of conventions/R.md.

#' The normalisation tab's UI: the options column (`config$panels`, or a panel function) and the QC tabset.
#' @param id The module id.
#' @param config The omic's `norm` group (`omicGuiConfig(<omic>)$norm`), or a view with its UI fields (KIT-097: trx).
#' @return A tag list.
#' @noRd
mod_omic_norm_ui <- function(id, config) {
    ns <- shiny::NS(id)
    shiny::tagList(shiny::wellPanel(shiny::fluidRow(omicNormOptionsColumn(ns, config), omicNormQcTabset(ns, config))))
}

#' The normalisation tab's server.
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label,selected_tab The tab's arguments, as
#'   `mod_lipid_norm_server`'s.
#' @param config The omic's `norm` group (`omicGuiConfig(<omic>)$norm`); the binding and readers are the configuration
#'   of `omic_type`.
#' @return The tab's reactive values.
#' @noRd
mod_omic_norm_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL, config) {
    shiny::moduleServer(id, function(input, output, session) {
        runOmicNormServer(input, output, session, id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab,
            config)
    })
}

#' The server's body: the state, the log and the aesthetics, then the observers and outputs in the lipid tab's order
#' (mod_lipid_norm_runtime_helpers.R:203-297 at f8d4c5fc).
#' @noRd
runOmicNormServer <- function(input, output, session, id, workflowData, experimentPaths, omicType, experimentLabel, selectedTab,
    norm) {
    logger::log_info("=== METABOLOMICS NORMALIZATION MODULE STARTED ===")
    logger::log_info(paste("Module ID:", id))
    config <- omicGuiConfig(omicType)
    config$norm <- norm
    normData <- omicNormReactiveState()
    addLog <- omicNormAddLog(normData)
    getPlotAesthetics <- omicNormPlotAestheticsGetter(input)
    shiny::observe(omicNormDetectAssays(workflowData, normData, config))
    if (!is.null(selectedTab) && isTRUE(norm$pre_qc_on_entry)) {
        shiny::observeEvent(selectedTab(), runOmicNormTabEntry(selectedTab(), workflowData, experimentPaths, normData, addLog,
            getPlotAesthetics, config), ignoreInit = FALSE)
    }
    shiny::observe(omicNormDesignChoices(session, workflowData))
    registerOmicNormLogOutput(output, omicNormLogRenderer(normData))
    output$itsd_selection_ui <- omicNormItsdSelectionUi(normData, session$ns)
    registerOmicNormRuvQcOutput(output, omicNormRuvQcUiRenderer(normData, session$ns))
    registerOmicNormQcImageOutputs(output, omicNormQcImageRenderer(normData, experimentPaths, qcDirKey = config$binding$binding$qc_dir_role),
        assaySlots = norm$assay_slots)
    registerOmicNormAssayLabelOutputs(output, omicNormAssayLabelRenderer(normData), norm$assay_slots)
    shiny::observe(omicNormItsdTables(output, workflowData, normData, config))
    registerOmicNormItsdSelectionTracking(input, normData)
    shiny::observeEvent(input$run_normalization, runOmicNormPipeline(input, workflowData, experimentPaths, normData, addLog,
        getPlotAesthetics, config))
    shiny::observeEvent(input$reset_normalization, resetOmicNorm(workflowData, normData, addLog, config))
    registerOmicNormRuvCancorOutputs(output, normData)
    output$correlation_filter_summary <- omicNormCorrelationSummary(normData, workflowData, config$binding$read_metadata)
    output$final_qc_plot <- omicNormFinalQcPlot(normData, getPlotAesthetics, workflowData, config)
    shiny::observeEvent(input$apply_correlation_filter, applyOmicNormCorrelation(input, workflowData, normData, addLog, config))
    shiny::observeEvent(input$skip_correlation_filter, skipOmicNormCorrelation(workflowData, normData, addLog, config))
    shiny::observeEvent(input[[norm$panels$export$button]], exportOmicNormSession(input, workflowData, experimentPaths,
        experimentLabel, normData, addLog, config))
    invisible(normData)
}

# ------------------------------------------------------------------------------------------------
# The start: the assays, the design's choices and the pre-normalisation QC on tab entry.
# ------------------------------------------------------------------------------------------------

#' The assay names, from the head state's metadata (KIT-075), when its class is the binding's.
#' @param config The omic's GUI configuration.
#' @noRd
omicNormDetectAssays <- function(workflowData, normData, config) {
    shiny::req(workflowData$state_manager)
    tryCatch({
        current <- config$binding$read_metadata(workflowData$state_manager)
        if (!identical(current$class, config$binding$binding$class)) return(invisible(NULL))
        normData$assay_names <- current$assay_names
        logger::log_info(paste("Detected assays:", paste(current$assay_names, collapse = ", ")))
        for (assay_name in normData$assay_names) {
            if (!assay_name %in% names(normData$itsd_selections)) normData$itsd_selections[[assay_name]] <- NULL
        }
    }, error = function(e) {
        logger::log_warn(paste("Could not detect assay names:", e$message))
    })
}

#' The aesthetics' and the RUV grouping's choices: the design's columns among the known variables.
#' @noRd
omicNormDesignChoices <- function(session, workflowData) {
    if (is.null(workflowData$design_matrix)) return(invisible(NULL))
    design_cols <- colnames(workflowData$design_matrix)
    plot_vars <- intersect(design_cols, c("group", "factor1", "factor2", "batch", "technical_replicate_id", "sample_id"))
    if (length(plot_vars) > 0) {
        default_var <- if ("group" %in% plot_vars) "group" else plot_vars[1]
        shiny::updateSelectInput(session, "color_variable", choices = plot_vars, selected = default_var)
        shiny::updateSelectInput(session, "shape_variable", choices = plot_vars, selected = default_var)
    }
    ruv_vars <- intersect(design_cols, c("group", "factor1", "factor2", "batch"))
    if (length(ruv_vars) > 0) {
        shiny::updateSelectInput(session, "ruv_grouping_variable", choices = ruv_vars,
            selected = if ("group" %in% ruv_vars) "group" else ruv_vars[1])
    }
    invisible(NULL)
}

#' On entry to the tab (`norm`), the pre-normalisation QC figures once, when the head state is the binding's class.
#' @noRd
runOmicNormTabEntry <- function(selectedTabValue, workflowData, experimentPaths, normData, addLog, getPlotAesthetics, config) {
    if (is.null(selectedTabValue) || !identical(selectedTabValue, "norm")) return(invisible(FALSE))
    logger::log_info("Normalization tab selected - checking if pre-QC needed")
    if (normData$pre_norm_qc_generated) return(invisible(FALSE))
    shiny::req(workflowData$state_manager)
    # KIT-075: the class only.
    current <- tryCatch(config$binding$read_metadata(workflowData$state_manager), error = function(e) NULL)
    if (is.null(current) || !identical(current$class, config$binding$binding$class)) return(invisible(FALSE))
    logger::log_info("Auto-triggering pre-normalization QC plots")
    shiny::withProgress(message = "Generating Pre-Normalization QC...", value = 0.5, {
        runOmicNormPreQc(workflowData, experimentPaths, normData, addLog, getPlotAesthetics, config)
    })
    normData$pre_norm_qc_generated <- TRUE
    invisible(TRUE)
}

#' The pre-normalisation QC figures of the head state: on the worker under a managed session (KIT-078), else of the
#' state read for this action only (KIT-075).
#' @noRd
runOmicNormPreQc <- function(workflowData, experimentPaths, normData, addLog, getPlotAesthetics, config) {
    logger::log_info("=== GENERATING PRE-NORMALIZATION QC PLOTS ===")
    shiny::req(workflowData$state_manager)
    binding <- config$binding$binding
    done <- function() {
        normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1
        normData$pre_norm_qc_generated <- TRUE
        logger::log_info("Pre-normalization QC plots generated successfully")
    }
    runner <- workflowOmicsManagedActions(binding$omic_type)
    if (!is.null(runner)) {
        current <- config$binding$read_metadata(workflowData$state_manager)
        if (is.null(current)) return(invisible(FALSE))
        if (length(current$assay_names) > 0 && is.null(normData$assay_names)) normData$assay_names <- current$assay_names
        aesthetics <- getPlotAesthetics()
        runner$run(list(omicsManagedProduct(paste0(binding$prefix, "_qc_figures"), list(stage = "post_filter",
            grouping_variable = aesthetics$color_var, shape_variable = aesthetics$shape_var), function(result) done())),
            failed = function(error) addLog(paste("Error generating Pre-QC:", conditionMessage(error))))
        return(invisible(TRUE))
    }
    current_s4 <- config$binding$read_state(workflowData$state_manager)
    if (is.null(current_s4)) {
        logger::log_warn("No S4 object available for QC plot generation")
        return(invisible(FALSE))
    }
    detected <- if (inherits(current_s4, binding$class)) names(methods::slot(current_s4, binding$data_slot))
    if (length(detected) > 0 && is.null(normData$assay_names)) {
        normData$assay_names <- detected
        logger::log_info(paste("Set assay names:", paste(detected, collapse = ", ")))
    }
    aesthetics <- getPlotAesthetics()
    tryCatch({
        binding$qc_plots(theObject = current_s4, experiment_paths = experimentPaths, stage = "post_filter",
            grouping_variable = aesthetics$color_var, shape_variable = aesthetics$shape_var)
        done()
        invisible(TRUE)
    }, error = function(e) {
        logger::log_error(paste("Error generating pre-normalization QC:", e$message))
        addLog(paste("Error generating Pre-QC:", e$message))
        invisible(FALSE)
    })
}

# ------------------------------------------------------------------------------------------------
# The actions.
# ------------------------------------------------------------------------------------------------

#' The run's settings: the managed twin's list (`runWorkflowOmicsManagedNormalisation`).
#' @param itsdFeatureIds The manual ITSD selection's feature ids (managed: from the candidates; in process: NULL, D8).
#' @noRd
omicNormRunSettings <- function(input, aesthetics, normData, experimentPaths, binding, itsdFeatureIds = NULL) {
    list(apply_itsd = input$apply_itsd, itsd_aggregation = input$itsd_aggregation, itsd_feature_ids = itsdFeatureIds,
        log_offset = input$log_offset, norm_method = input$norm_method, ruv_mode = input$ruv_mode,
        ruv_parameters = list(ruv_mode = input$ruv_mode, percentage_min = input$auto_percentage_min,
            percentage_max = input$auto_percentage_max, ruv_grouping_variable = input$ruv_grouping_variable,
            separation_metric = input$separation_metric, k_penalty_weight = input$k_penalty_weight,
            adaptive_k_penalty = input$adaptive_k_penalty, manual_k = input$ruv_k, manual_percentage = input$ruv_percentage),
        grouping_variable = aesthetics$color_var, shape_variable = aesthetics$shape_var, assay_names = normData$assay_names,
        qc_dir = experimentPaths[[binding$qc_dir_role]])
}

#' The normalisation pipeline: the managed sequence with a runner (KIT-078), else the in-process twin under the progress
#' display (KIT-105 D6, D9); then the tab's completion.
#' @noRd
runOmicNormPipeline <- function(input, workflowData, experimentPaths, normData, addLog, getPlotAesthetics, config) {
    shiny::req(workflowData$state_manager)
    addLog("Starting normalization pipeline...")
    finish <- function() {
        normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1
        addLog("Normalization pipeline complete!")
        shiny::showNotification("Normalization pipeline complete!", type = "message")
    }
    binding <- config$binding$binding
    manager <- workflowData$state_manager
    runner <- workflowOmicsManagedActions(binding$omic_type)
    if (!is.null(runner)) {
        aesthetics <- getPlotAesthetics()
        normData$post_filter_obj <- getWorkflowStateReference(manager)
        itsd_feature_ids <- if (isTRUE(input$apply_itsd)) omicsManualItsdFeatureIds(normData$itsd_selections, normData$itsd_candidate_ids, addLog)
        settings <- omicNormRunSettings(input, aesthetics, normData, experimentPaths, binding, itsd_feature_ids)
        return(config$norm$run$managed(binding, runner, manager, settings, normData, workflowData$config_list, addLog,
            finish = finish, failed = function(error) {
                addLog(paste("ERROR:", conditionMessage(error)))
                logger::log_error(paste("Normalization pipeline error:", conditionMessage(error)))
            }))
    }
    settings <- c(omicNormRunSettings(input, getPlotAesthetics(), normData, experimentPaths, binding), list(experiment_paths = experimentPaths))
    shiny::withProgress(message = "Running normalization pipeline...", value = 0, {
        config$norm$run$in_process(binding, manager, settings, normData, workflowData$config_list, addLog,
            advance = function(detail) shiny::incProgress(1 / 6, detail = detail),
            notify = function(...) shiny::showNotification(...), finish = finish, itsd_table = config$norm$itsd_table)
    })
}

#' Reset to the state the run started from: the store step `<prefix>_reset` recommits its record (KIT-075), and the
#' tab's flags and results are cleared.
#' @noRd
resetOmicNorm <- function(workflowData, normData, addLog, config) {
    shiny::req(workflowData$state_manager)
    tryCatch({
        if (!is.null(normData$post_filter_obj)) {
            step <- paste0(config$binding$binding$prefix, "_reset")
            commitWorkflowOmicsStoreStep(workflowData$state_manager, step, step, workflowData$config_list,
                "Reset to pre-normalization state", source = normData$post_filter_obj)
        }
        normData$normalization_complete <- FALSE
        normData$ruv_complete <- FALSE
        normData$correlation_filtering_complete <- FALSE
        normData$post_norm_obj <- NULL
        normData$ruv_corrected_obj <- NULL
        normData$correlation_filtered_obj <- NULL
        normData$ruv_optimization_results <- list()
        addLog("Reset to pre-normalization state")
        shiny::showNotification("Reset to pre-normalization state", type = "message")
        invisible(TRUE)
    }, error = function(e) {
        addLog(paste("ERROR during reset:", e$message))
        shiny::showNotification(paste("Error:", e$message), type = "error")
        invisible(FALSE)
    })
}

#' The tab's completion after a correlation step: `quality_control` and `normalization` complete.
#' @noRd
completeOmicNormTab <- function(workflowData) {
    updated_status <- workflowData$tab_status
    updated_status$quality_control <- "complete"
    updated_status$normalization <- "complete"
    workflowData$tab_status <- updated_status
    invisible(updated_status)
}

#' The correlation filter on the head (the normalisation's last commit): the pairs and the selection, committed under
#' `<prefix>_correlation_filtered` (KIT-031, KIT-072); on the worker under a managed session (KIT-078).
#' @noRd
applyOmicNormCorrelation <- function(input, workflowData, normData, addLog, config) {
    shiny::req(workflowData$state_manager)
    shiny::req(normData$ruv_complete || normData$normalization_complete)
    threshold <- input$min_pearson_correlation_threshold
    addLog(paste("Applying correlation filter (threshold:", threshold, ")"))
    shiny::showNotification("Applying correlation filter...", id = "corr_working", duration = NULL)
    binding <- config$binding$binding
    tryCatch({
        shiny::req(normData$ruv_corrected_obj %||% normData$post_norm_obj)
        logger::log_info("Calculating Pearson correlations per sample pair...")
        finish <- function(corr_results, reference) {
            normData$correlation_results <- corr_results
            normData$correlation_filtered_obj <- reference # KIT-075: the state's reference
            normData$correlation_filtering_complete <- TRUE
            completeOmicNormTab(workflowData)
            addLog("Correlation filtering complete")
            shiny::removeNotification("corr_working")
            shiny::showNotification("Correlation filtering complete! Ready for DE analysis.", type = "message")
        }
        runner <- workflowOmicsManagedActions(binding$omic_type)
        if (!is.null(runner)) {
            config$norm$correlation$managed(binding, runner, input$ruv_grouping_variable, threshold, workflowData$config_list, finish,
                function(error) {
                    addLog(paste("ERROR in correlation filtering:", conditionMessage(error)))
                    logger::log_error(paste("Correlation filtering error:", conditionMessage(error)))
                    shiny::removeNotification("corr_working")
                })
            return(invisible(TRUE))
        }
        config$norm$correlation$in_process(binding, workflowData$state_manager, input$ruv_grouping_variable, threshold,
            workflowData$config_list, finish)
        invisible(TRUE)
    }, error = function(e) {
        addLog(paste("ERROR in correlation filtering:", e$message))
        logger::log_error(paste("Correlation filtering error:", e$message))
        shiny::removeNotification("corr_working")
        shiny::showNotification(paste("Error:", e$message), type = "error")
        invisible(FALSE)
    })
}

#' Skip the correlation filter: the head's record committed under `<prefix>_norm_complete` (a store step), and the tab
#' complete. The flags are set before the commit or after it, per omic (`norm$skip_correlation_flags`, D7). F105-02: a
#' failing commit is caught and notified, and the tab's flags are put back, so the session goes on with its state
#' unchanged (it was unhandled, and shiny ended the session).
#' @noRd
skipOmicNormCorrelation <- function(workflowData, normData, addLog, config) {
    shiny::req(workflowData$state_manager)
    shiny::req(normData$ruv_complete || normData$normalization_complete)
    reference <- normData$ruv_corrected_obj %||% normData$post_norm_obj
    if (is.null(reference)) return(invisible(FALSE))
    before <- identical(config$norm$skip_correlation_flags, "before_commit")
    held <- list(obj = normData$correlation_filtered_obj, complete = normData$correlation_filtering_complete)
    set_flags <- function(obj, complete) {
        normData$correlation_filtered_obj <- obj
        normData$correlation_filtering_complete <- complete
    }
    if (before) set_flags(reference, TRUE)
    step <- paste0(config$binding$binding$prefix, "_norm_complete")
    committed <- tryCatch({
        commitWorkflowOmicsStoreStep(workflowData$state_manager, step, step, workflowData$config_list,
            "Normalization complete (correlation filtering skipped)")
        TRUE
    }, error = function(e) {
        set_flags(held$obj, held$complete)
        addLog(paste("ERROR:", conditionMessage(e)))
        logger::log_error(paste("Skip correlation filtering error:", conditionMessage(e)))
        shiny::showNotification(paste("Error:", conditionMessage(e)), type = "error")
        FALSE
    })
    if (!committed) return(invisible(FALSE))
    completeOmicNormTab(workflowData)
    addLog("Correlation filtering skipped - ready for DE analysis")
    shiny::showNotification("Normalization complete! Proceeding to DE analysis.", type = "message")
    if (!before) set_flags(reference, TRUE)
    invisible(TRUE)
}

# ------------------------------------------------------------------------------------------------
# The export.
# ------------------------------------------------------------------------------------------------

#' The export's directory: the source directory, else the export directory (created when missing).
#' @noRd
omicNormExportDir <- function(experimentPaths) {
    source_dir <- experimentPaths$source_dir
    if (!is.null(source_dir) && dir.exists(source_dir)) return(source_dir)
    source_dir <- experimentPaths$export_dir
    if (is.null(source_dir)) stop("Could not find a valid directory to save session data.")
    if (!dir.exists(source_dir)) dir.create(source_dir, recursive = TRUE)
    source_dir
}

#' Export the normalised session: the owner's publication (KIT-071, KIT-042/043) of the session list the core collects
#' (`norm$export$collect`) and the files its saver writes (`norm$export$save`); on a managed session the worker
#' serialises the state and the owner's retry path publishes it (KIT-082).
#' @noRd
exportOmicNormSession <- function(input, workflowData, experimentPaths, experimentLabel, normData, addLog, config) {
    logger::log_info("=== EXPORT NORMALIZED SESSION BUTTON CLICKED ===")
    shiny::req(workflowData$state_manager)
    if (!normData$normalization_complete) {
        shiny::showNotification("Please complete normalization before exporting session data.", type = "warning", duration = 5)
        return(invisible(FALSE))
    }
    report_failure <- function(e) {
        message <- conditionMessage(e)
        logger::log_error(paste("*** ERROR in session export:", message, "***"))
        addLog(paste("Export error:", message))
        shiny::showNotification(paste("Export error:", message), type = "error", duration = 10)
        invisible(FALSE)
    }
    tryCatch({
        source_dir <- omicNormExportDir(experimentPaths)
        manager <- workflowData$state_manager
        if (!inherits(manager, "DIAStateAdapter")) {
            cli::cli_abort("the session export needs the session's state owner", class = "saved_entity_scope")
        }
        binding <- config$binding$binding
        export <- config$norm$export
        collect <- function(workflowData, normData, input) {
            session_data <- export$collect(binding, workflowData, normData, input, experimentLabel, export, exportTimestamp = Sys.time())
            logger::log_info("*** EXPORT: Gathered session data successfully ***")
            logger::log_info(sprintf("*** EXPORT: Assays: %s ***", paste(normData$assay_names, collapse = ", ")))
            logger::log_info(sprintf("*** EXPORT: Contrasts available: %d ***",
                ifelse(is.null(session_data$contrasts_tbl), 0, nrow(session_data$contrasts_tbl))))
            session_data
        }
        publish <- function() publishOmicNormExport(manager, workflowData, normData, input, source_dir, collect, addLog, config)
        runner <- workflowOmicsManagedActions(binding$omic_type)
        if (!is.null(runner)) {
            progress <- runner$progress("Exporting normalized session data...")
            return(startWorkflowManagedSessionExport(manager, runner, workflowData, normData, input, source_dir, collect, export$prefix,
                publish = function() tryCatch({
                    progress$close()
                    publish()
                }, error = report_failure),
                failed = function(error) {
                    progress$close()
                    report_failure(error)
                }))
        }
        publish()
    }, error = report_failure)
}

#' The export's publication under the owner's scope and the progress display, then the tab's report.
#' @noRd
publishOmicNormExport <- function(manager, workflowData, normData, input, sourceDir, collect, addLog, config) {
    binding <- config$binding$binding
    export <- config$norm$export
    artifacts <- withWorkflowStateScope(manager, function() {
        shiny::withProgress(message = "Exporting normalized session data...", value = 0, {
            shiny::incProgress(0.2, detail = "Gathering data and publishing outputs...")
            result <- runWorkflowOwnedSessionExport(manager, workflowData, normData, input, sourceDir, collectFn = collect,
                saveFn = function(sessionData, sourceDir, publicationRequest) {
                    export$save(binding, sessionData, sourceDir, publicationRequest, export)
                },
                serializeFn = saveRDS, prefix = export$prefix)
            shiny::incProgress(0.8, detail = "Export complete")
            result$exportArtifacts
        })
    })
    addLog(paste("Exported comprehensive session data to:", artifacts$session_filepath))
    shiny::showNotification(sprintf("Session data exported successfully!\nSaved as: %s\nSee summary file for details.",
        artifacts$session_filename), type = "message", duration = 10)
    logger::log_info("=== EXPORT NORMALIZED SESSION COMPLETED SUCCESSFULLY ===")
    invisible(artifacts)
}
