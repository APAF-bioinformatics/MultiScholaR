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
# `mod_omic_norm_server(...)` holds the tab's reactive state (state references only, KIT-075; an omic with its own
# server may add small display objects, as proteomics does: a RUV result and three plot objects), its log, the pre-QC on
# tab entry and the run, reset, correlation apply and skip, and export actions. Each action calls the core: the
# managed function with the omic's runner, the in-process twin without one (D6). KIT-110: trx's tab too; the frame
# registers only what an omic's tabs and panels use, and calls an omic's own run and reset handlers when its `norm`
# names them (D3, D4). Recorded divergence (KIT-110 D7): trx's tab now writes the two console start lines below, and
# its observers and outputs are registered in the frame's order. Functions a test mocks are called
# by name when they run. The module frame keeps no locals (the testServer masking rule; test-omic-design-module.R).
# Recorded divergence from conventions/R-shiny.md: the in-process or managed actions in place of ExtendedTask are
# kept, under the existing-code scope rule of conventions/R.md.

#' The normalisation tab's UI: the options column (`config$panels`) and the QC tabset; above them, with its gate
#' on, the RUV design section (KIT-129 W2).
#' @param id The module id.
#' @param config The omic's `norm` group (`omicGuiConfig(<omic>)$norm`).
#' @return A tag list.
#' @noRd
mod_omic_norm_ui <- function(id, config) {
    ns <- shiny::NS(id)
    shiny::tagList(shiny::wellPanel(omicNormRuvDesignSection(ns, config),
        shiny::fluidRow(omicNormOptionsColumn(ns, config), omicNormQcTabset(ns, config))))
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
#' (mod_lipid_norm_runtime_helpers.R:203-297 at f8d4c5fc). KIT-110 D3: each registration is made when the omic's tabs
#' or panels use it (for lipid and metab every one is, in this order); D4: the run and reset are the omic's own handlers
#' when its `norm` names them (trx), else the frame's.
#' @noRd
runOmicNormServer <- function(input, output, session, id, workflowData, experimentPaths, omicType, experimentLabel, selectedTab,
    norm) {
    # KIT-122 (audit F7): the line names the tab's omic (every omic's session logged "METABOLOMICS").
    logger::log_info(sprintf("=== %s NORMALIZATION MODULE STARTED ===", toupper(omicType)))
    logger::log_info(paste("Module ID:", id))
    # KIT-120 (audit F4): an omic with its own server needs only its `norm` group, so the exported proteomics wrapper
    # keeps accepting any `omic_type` label (it used it for file names only).
    config <- if (is.null(norm$server)) omicGuiConfig(omicType) else list()
    config$norm <- norm
    tabs <- norm$plot_tabs
    panels <- names(norm$panels)
    normData <- omicNormReactiveState()
    addLog <- omicNormAddLog(normData)
    getPlotAesthetics <- omicNormPlotAestheticsGetter(input)
    # Fixed assay names are set once (trx: its one count table); otherwise they are read from the head state.
    if (is.null(norm$assay_names)) {
        shiny::observe(omicNormDetectAssays(workflowData, normData, config))
    } else {
        normData$assay_names <- norm$assay_names
    }
    if (!is.null(selectedTab) && isTRUE(norm$pre_qc_on_entry) && is.null(norm$server)) {
        shiny::observeEvent(selectedTab(), runOmicNormTabEntry(selectedTab(), workflowData, experimentPaths, normData, addLog,
            getPlotAesthetics, config), ignoreInit = FALSE)
    }
    if (any(c("aesthetics", "ruv_iii") %in% panels)) shiny::observe(omicNormDesignChoices(session, workflowData))
    # KIT-129 W2: the RUV design section's observers, for every omic; none while its gate is off.
    registerOmicNormRuvDesign(input, session, workflowData, norm$ruv_design, output, normData, norm$assay_slots)
    # KIT-120: an omic with its own server (`norm$server`, a module function's name; proteomics) registers its actions and its other outputs
    # itself, on the frame's state. The frame keeps the design's choices above and the stage grid.
    if (!is.null(norm$server)) {
        registerOmicNormStageGrid(output, normData, experimentPaths, norm, norm$qc_dir %||% config$binding$binding$qc_dir_role)
        omicGuiFunction(norm$server)(input, output, session, list(workflow_data = workflowData, experiment_paths = experimentPaths,
            omic_type = omicType, experiment_label = experimentLabel, selected_tab = selectedTab, norm_data = normData,
            add_log = addLog, plot_aesthetics = getPlotAesthetics, config = config))
        return(invisible(normData))
    }
    registerOmicNormLogOutput(output, omicNormLogRenderer(normData))
    # KIT-172: the sentence of a project whose record used a log2 offset.
    if ("log2" %in% panels) output$log2_recorded_offset <- omicNormRecordedOffsetUi(workflowData, normData, config$binding$binding$prefix)
    if ("itsd" %in% tabs) output$itsd_selection_ui <- omicNormItsdSelectionUi(normData, session$ns)
    registerOmicNormRuvQcOutput(output, omicNormRuvQcUiRenderer(normData, session$ns))
    registerOmicNormStageGrid(output, normData, experimentPaths, norm, norm$qc_dir %||% config$binding$binding$qc_dir_role)
    if ("itsd" %in% tabs) {
        shiny::observe(omicNormItsdTables(output, workflowData, normData, config))
        registerOmicNormItsdSelectionTracking(input, normData)
    }
    shiny::observeEvent(input$run_normalization, if (is.null(norm$run$handler)) {
        runOmicNormPipeline(input, workflowData, experimentPaths, normData, addLog, getPlotAesthetics, config)
    } else {
        runOmicNormHandler(input, workflowData, experimentPaths, normData, addLog, config)
    })
    shiny::observeEvent(input$reset_normalization, if (is.null(norm$reset)) {
        resetOmicNorm(workflowData, normData, addLog, config)
    } else {
        resetOmicNormHandler(workflowData, normData, addLog, config)
    })
    registerOmicNormRuvCancorOutputs(output, normData)
    if ("correlation_filter" %in% tabs) {
        output$correlation_filter_summary <- omicNormCorrelationSummary(normData, workflowData, config$binding$read_metadata)
        output$final_qc_plot <- omicNormFinalQcPlot(normData, getPlotAesthetics, workflowData, config)
        shiny::observeEvent(input$apply_correlation_filter, applyOmicNormCorrelation(input, workflowData, normData, addLog, config))
        shiny::observeEvent(input$skip_correlation_filter, skipOmicNormCorrelation(workflowData, normData, addLog, config))
    }
    if ("export" %in% panels) {
        shiny::observeEvent(input[[norm$panels$export$button]], exportOmicNormSession(input, workflowData, experimentPaths,
            experimentLabel, normData, addLog, config))
    }
    invisible(normData)
}

#' The stage grid's image slots and labels for `norm$output_assay_slots` assay rows (else `norm$assay_slots`; trx
#' registers two rows although its tab shows one, P-21): the files `<stem>_<stage prefix>_<kind>.png` in the QC
#' directory, with the omic's fixed file stem `norm$image_stem` and stage prefixes `norm$stage_prefixes` when it has
#' them (trx: `trx`, and `raw`, `tmm`, `ruvs`; KIT-110 D3), else the assay's sanitised name and the registrar's.
#' @param qcDirKey The QC directory's key in `experimentPaths` (the binding's `qc_dir_role`).
#' @noRd
registerOmicNormStageGrid <- function(output, normData, experimentPaths, norm, qcDirKey) {
    slots <- norm$output_assay_slots %||% norm$assay_slots
    stem <- if (isTRUE(norm$image_stemless)) "" else norm$image_stem
    renderer <- if (is.null(stem)) {
        omicNormQcImageRenderer(normData, experimentPaths, qcDirKey = qcDirKey)
    } else {
        omicNormQcImageRenderer(normData, experimentPaths, sanitizeAssayNameFn = function(assayName) stem, qcDirKey = qcDirKey)
    }
    if (is.null(norm$stage_prefixes)) {
        registerOmicNormQcImageOutputs(output, renderer, assaySlots = slots)
    } else {
        registerOmicNormQcImageOutputs(output, renderer, stagePrefixes = norm$stage_prefixes, assaySlots = slots)
    }
    if (!isFALSE(norm$assay_labels)) registerOmicNormAssayLabelOutputs(output, omicNormAssayLabelRenderer(normData), slots)
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
#' @param ruvDesign KIT-129 fourth pass item 2: the RUV design section's choice (`omicNormRuvDesignChoice`); NULL (the
#'   gate off) or a legacy design keeps today's settings.
#' @param logOffset KIT-172: the log2 offset when no input supplies one (the tab no longer renders the field): the
#'   configuration's, 0 (plain log2). A supplied `log_offset` input (a test, a replayed setting) is passed through.
#' @noRd
omicNormRunSettings <- function(input, aesthetics, normData, experimentPaths, binding, itsdFeatureIds = NULL, ruvDesign = NULL,
    logOffset = 0) {
    settings <- list(apply_itsd = input$apply_itsd, itsd_aggregation = input$itsd_aggregation, itsd_feature_ids = itsdFeatureIds,
        log_offset = input$log_offset %||% logOffset, norm_method = input$norm_method, ruv_mode = input$ruv_mode,
        ruv_parameters = list(ruv_mode = input$ruv_mode, percentage_min = input$auto_percentage_min,
            percentage_max = input$auto_percentage_max, ruv_grouping_variable = input$ruv_grouping_variable,
            separation_metric = input$separation_metric, k_penalty_weight = input$k_penalty_weight,
            adaptive_k_penalty = input$adaptive_k_penalty, manual_k = input$ruv_k, manual_percentage = input$ruv_percentage),
        grouping_variable = aesthetics$color_var, shape_variable = aesthetics$shape_var, assay_names = normData$assay_names,
        qc_dir = experimentPaths[[binding$qc_dir_role]])
    omicNormApplyRuvDesign(settings, ruvDesign, input$ruv_grouping_variable)
}

#' The log2 offset the project's record used (KIT-172): the declared parameter of the latest `<prefix>_log2` step in the
#' state history, or NULL when the project has none (or no session owner).
#' @noRd
omicNormRecordedLog2Offset <- function(manager, prefix) {
    snapshot <- tryCatch(manager$.context$snapshot, error = function(e) NULL)
    if (is.null(snapshot)) return(NULL)
    for (state in rev(as.character(unlist(snapshot$state_history)))) {
        provenance <- snapshot$step_provenance[[state]]
        if (identical(provenance$operation, paste0(prefix, "_log2"))) return(provenance$values$offset)
    }
    NULL
}

#' The recorded-offset sentence's renderer (KIT-172): the sentence when the project's record used a log2 offset above 0,
#' else nothing. It is read again when the tab runs its log step or the project's tab status changes (a project opened).
#' @noRd
omicNormRecordedOffsetUi <- function(workflowData, normData, prefix) {
    shiny::renderUI({
        normData$post_log2_obj
        workflowData$tab_status
        offset <- omicNormRecordedLog2Offset(workflowData$state_manager, prefix)
        if (!is.numeric(offset) || length(offset) != 1L || is.na(offset) || offset <= 0) return(NULL)
        shiny::helpText(sprintf(omicGuiLog2Texts()[["recorded_offset"]], format(offset)))
    })
}

#' The new design's RUV request of lipidomics and metabolomics (KIT-129 fourth pass item 2): the RUV step runs (its mode
#' is the design's id, never `skip`) with the new design's parameters (`omicsRuvDesignParameterNames`). The
#' technical-replicate design takes the declared technical replicates as its replicate sets, so it names no column.
#' @param settings `omicNormRunSettings`' list.
#' @param choice `omicNormRuvDesignChoice`'s value, or NULL.
#' @param grouping Today's RUV grouping column (the tab's input).
#' @return `settings`.
#' @noRd
omicNormApplyRuvDesign <- function(settings, choice, grouping) {
    if (is.null(choice) || isTRUE(choice$legacy)) return(settings)
    settings$ruv_mode <- choice$ruv_design
    settings$ruv_parameters <- list(ruv_design = choice$ruv_design,
        ruv_grouping_variable = if (!identical(choice$ruv_design, "technical_replicates")) grouping,
        ruv_k = choice$ruv_k, ruv_k_rule = choice$ruv_k_rule, ruv_k_seed = choice$ruv_k_seed, ruv_control_ids = choice$ruv_control_ids)
    settings
}

#' The RUV design section's choice when a run starts (KIT-129 fourth pass item 2), or NULL; a refusal is the tab's
#' error notice and its log line, and the run does not start.
#' @param omic The omic (`binding$omic_type`).
#' @return `list(ok, choice)`.
#' @noRd
omicNormReadRuvDesign <- function(input, omic, addLog, route = NULL) {
    tryCatch(list(ok = TRUE, choice = omicNormRuvDesignChoice(input, omicNormRuvRunTable(omic, route))),
        omic_norm_ruv_design_invalid = function(error) {
            addLog(paste("ERROR:", conditionMessage(error)))
            shiny::showNotification(conditionMessage(error), type = "error", duration = 10)
            list(ok = FALSE, choice = NULL)
        })
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
    ruv_design <- omicNormReadRuvDesign(input, binding$omic_type, addLog)
    if (!ruv_design$ok) return(invisible(FALSE))
    runner <- workflowOmicsManagedActions(binding$omic_type)
    if (!is.null(runner)) {
        aesthetics <- getPlotAesthetics()
        normData$post_filter_obj <- getWorkflowStateReference(manager)
        itsd_feature_ids <- if (isTRUE(input$apply_itsd)) omicsManualItsdFeatureIds(normData$itsd_selections, normData$itsd_candidate_ids, addLog)
        settings <- omicNormRunSettings(input, aesthetics, normData, experimentPaths, binding, itsd_feature_ids, ruv_design$choice,
            logOffset = config$norm$panels$log2$offset %||% 0)
        return(config$norm$run$managed(binding, runner, manager, settings, normData, workflowData$config_list, addLog,
            finish = finish, failed = function(error) {
                addLog(paste("ERROR:", conditionMessage(error)))
                logger::log_error(paste("Normalization pipeline error:", conditionMessage(error)))
                # KIT-129 fourth pass item 3: under the gate a refusal of the worker's RUV step is the tab's error notice,
                # as in process.
                if (omicGuiRuvDesignEnabled()) shiny::showNotification(paste("Error:", conditionMessage(error)), type = "error")
            }))
    }
    settings <- c(omicNormRunSettings(input, getPlotAesthetics(), normData, experimentPaths, binding, ruvDesign = ruv_design$choice,
        logOffset = config$norm$panels$log2$offset %||% 0), list(experiment_paths = experimentPaths))
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

#' A `setStatus(tab, value)` callback writing one entry of the bus's `tab_status`.
#' @noRd
omicNormStatusSetter <- function(workflowData) {
    function(tab, value) {
        status <- workflowData$tab_status
        status[[tab]] <- value
        workflowData$tab_status <- status
    }
}

#' An omic's own run (`norm$run$handler`; KIT-110 D4, trx's `runTrxNormalisation`), with the session's state owner
#' required as trx's server did, the settings `norm$run$settings` takes from the tab's inputs, and the tab's notification
#' and status callbacks. The handler resolves the omic's runner itself; there is no progress display.
#' @noRd
runOmicNormHandler <- function(input, workflowData, experimentPaths, normData, addLog, config) {
    shiny::req(workflowData$state_manager)
    run <- config$norm$run
    # KIT-129 fourth pass item 2: a refusal of the RUV design section's choice is the tab's error notice (trx's settings
    # read it, `trxNormGridParameters`).
    settings <- tryCatch(run$settings(input), omic_norm_ruv_design_invalid = function(error) {
        addLog(paste("ERROR:", conditionMessage(error)))
        shiny::showNotification(conditionMessage(error), type = "error", duration = 10)
        NULL
    })
    if (is.null(settings)) return(invisible(FALSE))
    run$handler(run, config$binding$binding, workflowData$state_manager, settings, normData, workflowData$config_list,
        experimentPaths, addLog, notify = function(...) shiny::showNotification(...), setStatus = omicNormStatusSetter(workflowData))
}

#' An omic's own reset (`norm$reset`; KIT-110 D4, trx's `resetTrxNormalisation`), with no `req`, as trx's server called it.
#' @noRd
resetOmicNormHandler <- function(workflowData, normData, addLog, config) {
    config$norm$reset(workflowData$state_manager, normData, workflowData$config_list, addLog,
        notify = function(...) shiny::showNotification(...), setStatus = omicNormStatusSetter(workflowData))
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
                    # D7 (KIT-105 W5): metab's tab also shows the worker's failure.
                    if (isTRUE(config$norm$managed_correlation_notify)) {
                        shiny::showNotification(paste("Error:", conditionMessage(error)), type = "error")
                    }
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

#' The export's publication under the owner's scope and the progress display, then the tab's report. The display's
#' steps and the summary's time are the omic's (`norm$export$progress`, `norm$export$summary_time`; D7, KIT-105 W5).
#' @noRd
publishOmicNormExport <- function(manager, workflowData, normData, input, sourceDir, collect, addLog, config) {
    binding <- config$binding$binding
    export <- config$norm$export
    progress <- export$progress
    artifacts <- withWorkflowStateScope(manager, function() {
        shiny::withProgress(message = "Exporting normalized session data...", value = 0, {
            shiny::incProgress(0.2, detail = progress$gather)
            result <- runWorkflowOwnedSessionExport(manager, workflowData, normData, input, sourceDir, collectFn = collect,
                saveFn = function(sessionData, sourceDir, publicationRequest) {
                    # Metab's tab steps the display through the saver's files, and its summary carries the time it is
                    # written (`Sys.time()` is read when the summary is); lipid's carries the session's export time.
                    export$save(binding, sessionData, sourceDir, publicationRequest, export,
                        advance = if (isTRUE(progress$steps)) function(amount, detail) shiny::incProgress(amount, detail = detail),
                        summaryTime = if (identical(export$summary_time, "written")) Sys.time() else sessionData$export_timestamp)
                },
                serializeFn = saveRDS, prefix = export$prefix)
            if (length(progress$complete)) shiny::incProgress(0.8, detail = progress$complete)
            result$exportArtifacts
        })
    })
    addLog(paste("Exported comprehensive session data to:", artifacts$session_filepath))
    shiny::showNotification(sprintf("Session data exported successfully!\nSaved as: %s\nSee summary file for details.",
        artifacts$session_filename), type = "message", duration = 10)
    logger::log_info("=== EXPORT NORMALIZED SESSION COMPLETED SUCCESSFULLY ===")
    invisible(artifacts)
}
