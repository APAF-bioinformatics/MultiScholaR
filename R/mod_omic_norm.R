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

#' The values of the RUV grouping column whose samples belong to more than one group
#' (KIT-174 (a), VAL-29). NA is counted as a value.
#' @param design_matrix The design table.
#' @param column The RUV grouping column, or NULL (an omic without the selector).
#' @param group_column The design's group column.
#' @return The spanning values (character; NA for the missing value), empty when none or when
#'   `column` is NULL.
#' @noRd
omicNormRuvGroupingSpansGroups <- function(design_matrix, column, group_column = "group") {
    if (is.null(column) || !is.data.frame(design_matrix) ||
        !all(c(column, group_column) %in% names(design_matrix))) {
        return(character())
    }
    values <- factor(as.character(design_matrix[[column]]), exclude = NULL)
    groups <- as.character(design_matrix[[group_column]])
    spans <- vapply(split(groups, values), \(group) length(unique(group)) > 1L, logical(1))
    levels(values)[spans]
}

#' The RUV grouping warning's text (KIT-174 (a), VAL-29; proposed wording, implemented as
#' proposed, the owner edits before release, VAL-41): the column, the number of its values that
#' span groups, the number of its values.
#' @noRd
omicNormRuvGroupingWarningText <- function() {
    paste(
        "The RUV grouping column '%s' puts samples of different groups together (%d of its %d",
        "values contain samples of more than one group). RUV-III treats samples with the same",
        "value as replicates, so the correction would be built from the differences between your",
        "groups. Choose 'group' unless each value of the column lies inside one group."
    )
}

#' Whether the RUV grouping column gives RUV-III's replicate sets in this run (KIT-174 (a),
#' 3 item 2): the mode is not `skip`; under the RUV design gate the design is not `none` or
#' `technical_replicates`; on proteomics, declared technical replicates do not take the column's
#' place (they do on DIA with limpa under every design, and on the classic routes except under
#' `groups_all_features` and `supplied_controls`, which withdraw the declaration for the RUV
#' step, `protRuvGroupReplicateObject`).
#' @param mode The run's RUV mode.
#' @param design The gate's design id, or NULL with the gate off.
#' @param technical Whether the head state's design declares technical replicates (proteomics).
#' @param classic TRUE on the classic proteomics routes.
#' @return One logical.
#' @noRd
omicNormRuvGroupingIsReplicateSets <- function(mode, design = NULL, technical = FALSE,
                                               classic = FALSE) {
    if (is.null(mode) || identical(mode, "skip")) return(FALSE)
    if (isTRUE(design %in% c("none", "technical_replicates"))) return(FALSE)
    withdrawn <- classic && isTRUE(design %in% c("groups_all_features", "supplied_controls"))
    !(isTRUE(technical) && !withdrawn)
}

#' The head state's design table and group column (KIT-174 (a)): the table the run reads, from
#' the state's metadata (as `omicNormRuvHeadTechnicalReplicates`); NULL without a session state
#' or when it cannot be read.
#' @noRd
omicNormRuvHeadDesign <- function(manager) {
    if (!inherits(manager, "DIAStateAdapter")) return(NULL)
    tryCatch(withWorkflowStateScope(manager, function() {
        reference <- getWorkflowStateReference(manager)
        if (is.null(reference)) return(NULL)
        spec <- readWorkflowStateMetadata(manager, reference)$data
        design <- materializeWorkflowTable(manager$.context, spec$components$design_matrix)
        list(design = as.data.frame(design), group_id = spec$slots$group_id %||% "group")
    }), error = function(error) NULL)
}

#' The RUV grouping warning at a run action (KIT-174 (a), VAL-29): when the column gives
#' RUV-III's replicate sets (`omicNormRuvGroupingIsReplicateSets`) and one of its values holds
#' samples of more than one group, one warning notification and one log line with the same text
#' (and a line in the tab's log where the tab has one). The run proceeds; RUV is unchanged. The
#' design is the head state's, else the session's design table.
#' @param column The RUV grouping column (NULL: no warning).
#' @param addLog The tab's log, or NULL.
#' @return The spanning values, invisibly.
#' @noRd
omicNormWarnRuvGrouping <- function(workflowData, column, mode, design = NULL, technical = FALSE,
                                    classic = FALSE, showNotificationFn = shiny::showNotification,
                                    addLog = NULL) {
    if (is.null(column) || !omicNormRuvGroupingIsReplicateSets(mode, design, technical, classic)) {
        return(invisible(character()))
    }
    head <- omicNormRuvHeadDesign(shiny::isolate(workflowData$state_manager))
    table <- head$design %||% shiny::isolate(workflowData$design_matrix)
    spans <- omicNormRuvGroupingSpansGroups(table, column, head$group_id %||% "group")
    if (!length(spans)) return(invisible(spans))
    text <- sprintf(omicNormRuvGroupingWarningText(), column, length(spans),
        length(unique(as.character(table[[column]]))))
    showNotificationFn(text, type = "warning", duration = NULL)
    logger::log_warn(text)
    if (is.function(addLog)) addLog(paste("WARNING:", text))
    invisible(spans)
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
    # KIT-194: the run rule (a refusal, or the previous run's start selected), after every other refusal (7.3).
    started <- omicNormApplyRunStart(binding, manager, normData, addLog)
    if (inherits(started, "error")) {
        shiny::showNotification(conditionMessage(started), type = "error", duration = 15)
        return(invisible(FALSE))
    }
    # KIT-174 (a), VAL-29: the RUV grouping warning (a new design's mode is its id; a legacy
    # design keeps the tab's mode).
    choice <- ruv_design$choice
    new_design <- !is.null(choice) && !isTRUE(choice$legacy)
    ruv_mode <- if (new_design) choice$ruv_design else input$ruv_mode
    omicNormWarnRuvGrouping(workflowData, input$ruv_grouping_variable, ruv_mode, choice$ruv_design,
        addLog = addLog)
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
    manager <- workflowData$state_manager
    # KIT-254: whether this press returned the history to the run's result; a failure after it clears the tab's fields.
    pressed <- new.env(parent = emptyenv())
    pressed$cut <- FALSE
    tryCatch({
        # KIT-194: only the last run's result, or a correlation-step state made from it (7.4).
        assertOmicNormHeadIsResult(binding, manager, normData)
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
            # KIT-254: the check and the return to the run's result when the actions start (a press queued behind a running
            # action reads what that action left; KIT-194 A02's idiom).
            start <- function() {
                assertOmicNormHeadIsResult(binding, manager, normData)
                pressed$cut <- omicNormCorrelationStart(manager, normData, addLog)
            }
            config$norm$correlation$managed(binding, runner, input$ruv_grouping_variable, threshold, workflowData$config_list, finish,
                function(error) {
                    if (isTRUE(pressed$cut)) clearOmicNormCorrelation(normData)
                    addLog(paste("ERROR in correlation filtering:", conditionMessage(error)))
                    logger::log_error(paste("Correlation filtering error:", conditionMessage(error)))
                    shiny::removeNotification("corr_working")
                    # D7 (KIT-105 W5): metab's tab also shows the worker's failure.
                    if (isTRUE(config$norm$managed_correlation_notify)) {
                        shiny::showNotification(paste("Error:", conditionMessage(error)), type = "error")
                    }
                }, start)
            return(invisible(TRUE))
        }
        pressed$cut <- omicNormCorrelationStart(manager, normData, addLog)
        config$norm$correlation$in_process(binding, manager, input$ruv_grouping_variable, threshold,
            workflowData$config_list, finish)
        invisible(TRUE)
    }, error = function(e) {
        if (isTRUE(pressed$cut)) clearOmicNormCorrelation(normData)
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
    # KIT-194: only the last run's result, or a correlation-step state made from it (7.4); refused before any flag is set.
    # KIT-254: then, at a correlation-step state, the history returns to the run's result, so the Skip commits it.
    cut_back <- FALSE
    allowed <- tryCatch({
        assertOmicNormHeadIsResult(config$binding$binding, workflowData$state_manager, normData)
        cut_back <- omicNormCorrelationStart(workflowData$state_manager, normData, addLog)
        TRUE
    }, error = function(e) {
        addLog(paste("ERROR:", conditionMessage(e)))
        shiny::showNotification(paste("Error:", conditionMessage(e)), type = "error")
        FALSE
    })
    if (isFALSE(allowed)) return(invisible(FALSE))
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
        # KIT-254: after a return to the run's result the head is that table, so the tab's fields are cleared, not put back.
        if (cut_back) clearOmicNormCorrelation(normData) else set_flags(held$obj, held$complete)
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
    # KIT-194: the export writes the head as the normalised session, so only the last run's result or a correlation-step
    # state made from it (7.4).
    stale <- tryCatch(assertOmicNormHeadIsResult(config$binding$binding, workflowData$state_manager, normData, "export"),
        error = function(e) e)
    if (inherits(stale, "error")) return(report_failure(stale))
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

# ------------------------------------------------------------------------------------------------
# KIT-194 (VAL-99 (1), VAL-105): the normalisation's state guards of lipidomics, metabolomics and transcriptomics.
# Quality-control steps refuse a normalised table; "Run Normalization Pipeline" refuses a normalised input and, pressed
# again at the tab's own result, runs again from the table the previous run started from (D1 (b)); the correlation
# step, its Skip and the export refuse unless the head is the last run's result (or a correlation-step state made from
# it). The lineage is read from what the records already carry (the products a normalisation commit adds, the save labels
# of the states and values a state was made from); no field, schema or state is added. Replay and the wrappers are not
# guarded (scripts decide what they pass). Recorded divergence: the refusals use plain words for non-specialists, not
# conventions/R.md's "state the offending value" (as KIT-182).
# ------------------------------------------------------------------------------------------------

#' The user-facing texts of the guards (proposed; the owner edits them before release, VAL-41, VAL-99 (2)).
#' @return A named list of texts; `%s` slots are the sub-tab and the step.
#' @noRd
omicNormGuardTexts <- function() {
    list(
        qc_lead = "quality-control steps cannot be applied to a table that has already been normalised.",
        qc_revert = paste("To change quality control, open the %s sub-tab and press Revert until the current step is %s, the last",
            "quality-control result before normalisation; then apply the steps again and press \"Run Normalization Pipeline\"."),
        run_lead = "this table has already been normalised, so it cannot be normalised again.",
        run_unknown_lead = "this tab does not know which table the earlier normalisation started from, so it cannot run it again.",
        run_revert = paste("In the Quality Control tab, open the %s sub-tab and press Revert until the current step is %s, the last",
            "quality-control result before normalisation, then press \"Run Normalization Pipeline\"."),
        no_step = paste("This table was loaded from a saved normalised session, or no quality-control result before normalisation",
            "is left in this session, so it cannot be filtered or normalised again here. To change quality control or",
            "normalisation, import the data again."),
        busy = "A step is still running. Press \"Run Normalization Pipeline\" again when it has finished.",
        correlation = paste("the current table is not the one the last normalisation produced (Quality Control was reverted since),",
            "so the correlation filter would not use the normalised table the analysis reports. Press \"Run Normalization",
            "Pipeline\" again, then apply or skip the correlation filter."),
        export = paste("the current table is not the one the last normalisation produced (Quality Control was reverted since),",
            "so the correlation filter would not use the normalised table the analysis reports. Press \"Run Normalization",
            "Pipeline\" again, then export."),
        da = paste("The current table is not the result of the last normalisation: Quality Control was reverted, or the",
            "normalisation was reset, after it ran. Run the normalisation (and the correlation filter or its Skip) again, then",
            "run the analysis."),
        # KIT-207: the Summary's refusals (drafts for the owner's wording pass, VAL-41), by the guard's case.
        summary = paste("The current table is not the result of the last normalisation: Quality Control was reverted, or the",
            "normalisation was reset, after it ran. Run the normalisation (and the correlation filter or its Skip) again, or,",
            "if you have only just pressed Reset, press Revert once in a Quality Control sub-tab; then press this button again."),
        summary_trx = paste("The current table is not the result of the last normalisation: Quality Control was reverted, or the",
            "normalisation was reset, after it ran. Run the normalisation again, then press this button again."),
        summary_none = paste("No normalisation has been run in this session, so there is no final table to save. Run the",
            "normalisation (and the correlation filter or its Skip), then press this button again."),
        summary_loaded = paste("This session file was saved after Quality Control was reverted past the normalisation, so its",
            "table (\"%s\") was not normalised. Import the data again, run quality control and the normalisation, then press",
            "this button again."),
        # KIT-212 (VAL-172 (1), (5)): the Summary's refusals of a differential analysis not run on the current table, and of
        # DA files no analysis of this session made (drafts for the owner's wording pass, VAL-41); `%s` is the folder and
        # the files.
        summary_da_stale = paste("The differential analysis was run on an earlier table, not the current one: the table changed",
            "after the analysis ran (a new normalisation run, a Revert, or the correlation filter applied again). Run the",
            "differential analysis again, then press 'Save Workflow Arguments' before copying or generating the report."),
        summary_da_files_unmade = paste("The results folder %s holds differential analysis files that no differential analysis",
            "of this session made (%s): they are from an earlier session, an opened session file, or a contrast the last analysis",
            "did not produce. Run the differential analysis, remove any of these files it does not make again, then press 'Save",
            "Workflow Arguments' before copying or generating the report."),
        # KIT-243 (VAL-212 (2); part D): the Summary's refusals of a Save that does not describe the current analysis, and of
        # a differential analysis whose files were not all written (drafts for the owner's wording pass, VAL-41).
        summary_save_stale = paste("The saved workflow parameters do not describe the current analysis: a step or the",
            "differential analysis was run, reverted or changed after the analysis record (analysis_record.json) was last",
            "saved, or it was saved by another version of MultiScholaR. Press 'Save Workflow Arguments' again, then copy the",
            "results or generate the report. If this message returns after saving, the analysis record could not be written:",
            "see the warning shown when saving."),
        da_publication_pending = paste("The files of the last differential analysis were not all written to the project",
            "folder, so it may still hold files of an earlier analysis. Run the differential analysis again, then copy the",
            "results or generate the report."))
}

#' Whether the binding is transcriptomics' (one normalisation state, `trx_ruvs`).
#' @noRd
omicNormIsTrx <- function(binding) identical(binding$omic_type, "transcriptomics")

#' The normalisation's own states.
#' @param binding The omic's declaration binding.
#' @return State names.
#' @noRd
omicNormStates <- function(binding) {
    if (omicNormIsTrx(binding)) return("trx_ruvs")
    paste0(binding$prefix, c("_itsd_norm", "_log2", "_normalized", "_ruv_corrected", "_correlation_filtered", "_norm_complete"))
}

#' The save labels of normalised values: the normalisation's state names and its kernels' owned checkpoints.
#' @noRd
omicNormLineageLabels <- function(binding) {
    if (omicNormIsTrx(binding)) return(c("trx_ruvs", "owner:trx-ruvs"))
    c(omicNormStates(binding), paste0("owner:", binding$prefix, c("-itsd-normalised", "-log2", "-normalised", "-ruv-corrected")))
}

#' The products a normalisation commit adds to its state (`<operation>_science`).
#' @noRd
omicNormProducts <- function(binding) {
    if (omicNormIsTrx(binding)) return("trx_ruvs_science")
    paste0(binding$prefix, c("_itsd_normalisation", "_log2", "_between_sample_normalisation", "_ruv"), "_science")
}

#' A state's bounded record, or NULL (logged) when it cannot be read as a component checkpoint.
#' @noRd
omicNormReadSpec <- function(manager, reference) {
    read <- function(error) {
        logger::log_warn(paste("The normalisation guard could not read a state's record:", conditionMessage(error)))
        NULL
    }
    tryCatch(readWorkflowStateMetadata(manager, reference)$data, saved_entity_unsupported = read, saved_entity_corrupt = read)
}

#' Whether a saved state holds normalised values: its record carries a normalisation's product (i), or it was made from
#' a normalisation's state or owned value (ii; dependencies accumulate along a lineage).
#' @param binding The omic's declaration binding.
#' @param manager The session's state owner.
#' @param reference The state's saved reference.
#' @return TRUE or FALSE.
#' @noRd
omicNormalisedLineage <- function(binding, manager, reference) {
    if (!isSavedWorkflowReference(reference)) return(FALSE)
    spec <- omicNormReadSpec(manager, reference)
    if (is.null(spec)) return(FALSE)
    labels <- omicNormLineageLabels(binding)
    any(names(spec$products) %in% omicNormProducts(binding)) ||
        any(vapply(spec$dependencies, \(d) is.list(d) && isTRUE(d$identity$parameters %in% labels), logical(1)))
}

#' Whether a state was adopted from a session file by Load Filtered Session (its provenance marked `adopted`).
#' @noRd
omicNormAdopted <- function(manager, state) isTRUE(manager$.context$snapshot$step_provenance[[state]]$adopted)

#' Whether the head holds normalised values: a normalisation state, a state adopted from a loaded session (A08), or a
#' state with normalised lineage.
#' @noRd
omicNormHeadNormalised <- function(binding, manager) {
    head <- manager$current_state
    if (!is.character(head) || length(head) != 1L) return(FALSE)
    head %in% omicNormStates(binding) || omicNormAdopted(manager, head) ||
        omicNormalisedLineage(binding, manager, getWorkflowStateReference(manager))
}

#' Whether a state was made from the last run (A01): it is the states map's `<prefix>_normalized` (the latest run's), it
#' depends on that entity (RUV, correlation: dependencies accumulate), or it recommits its record (Skip at the result).
#' @return TRUE or FALSE; FALSE when the session holds no normalisation or a record cannot be read.
#' @noRd
omicNormDerivesFromCurrentRun <- function(binding, manager, reference) {
    normalised <- getWorkflowStateReference(manager, paste0(binding$prefix, "_normalized"))
    if (is.null(normalised) || !isSavedWorkflowReference(reference)) return(FALSE)
    if (identical(reference$entity_id, normalised$entity_id)) return(TRUE)
    spec <- omicNormReadSpec(manager, reference)
    if (is.null(spec)) return(FALSE)
    if (any(vapply(spec$dependencies, \(d) is.list(d) && identical(d$entity_id, normalised$entity_id), logical(1)))) return(TRUE)
    identical(spec$components, omicNormReadSpec(manager, normalised)$components)
}

#' The step to return to: walking the active history back from before its first normalisation state (from the head when
#' there is none), the first state that is not normalised, not adopted, not `initial` and has a reference.
#' @return The state name, or NULL when there is none (a table loaded from a session file).
#' @noRd
omicNormPreNormalisationStep <- function(binding, manager) {
    history <- as.character(unlist(manager$getHistory()))
    first <- match(TRUE, history %in% omicNormStates(binding))
    end <- if (is.na(first)) length(history) else first - 1L
    for (state in rev(history[seq_len(end)])) {
        if (identical(state, "initial") || omicNormAdopted(manager, state)) next
        reference <- getWorkflowStateReference(manager, state)
        if (is.null(reference)) next
        if (!omicNormalisedLineage(binding, manager, reference)) return(state)
    }
    NULL
}

#' A state named as the QC sub-tab that made it, with the name the Revert output shows (F-s).
#' @noRd
omicNormStepName <- function(binding, state) {
    titles <- c(raw_data_s4 = "the design step (before quality control)", intensity_filtered = "Intensity Filter",
        duplicates_resolved = "Duplicate Resolution", qc_complete = "Finalize QC", expression_filtered = "Expression Filter",
        reset = "the table saved by Reset to Pre-Normalization")
    own <- paste0(binding$prefix, "_")
    suffix <- if (startsWith(state, own)) substring(state, nchar(own) + 1L) else ""
    if (!suffix %in% names(titles)) return(sprintf("\"%s\"", state))
    sprintf("\"%s\" (shown as %s)", titles[[suffix]], state)
}

#' A refusal's text: its first sentence, then where to Revert to (or, with no step left, the A07 text).
#' @param lead The first sentence.
#' @param revert The Revert sentence's template (`%s` the sub-tab, then the step).
#' @noRd
omicNormRefusalText <- function(binding, manager, lead, revert) {
    texts <- omicNormGuardTexts()
    step <- omicNormPreNormalisationStep(binding, manager)
    if (is.null(step)) return(paste(lead, texts$no_step))
    sub_tab <- if (omicNormIsTrx(binding)) "Expression Filter" else "Intensity Filter or Duplicate Resolution"
    paste(lead, sprintf(revert, sub_tab, omicNormStepName(binding, step)))
}

#' The QC guard (7.2): a QC step refuses a normalised head. A no-op when the manager is not a session owner.
#' @param binding The omic's declaration binding.
#' @param manager The session's state manager.
#' @return `manager`, invisibly; otherwise a condition of class `omic_qc_input_normalised`.
#' @noRd
assertOmicQcHeadNotNormalised <- function(binding, manager) {
    if (!inherits(manager, "DIAStateAdapter") || !omicNormHeadNormalised(binding, manager)) return(invisible(manager))
    texts <- omicNormGuardTexts()
    rlang::abort(omicNormRefusalText(binding, manager, texts$qc_lead, texts$qc_revert), class = "omic_qc_input_normalised")
}

#' The run rule (7.3): refused while the omic's managed runner is busy (A02); at the tab's own result, the state the
#' previous run started from when this tab remembers it and it is unchanged in the active history and not normalised
#' (D1 (b)), else refused; refused on a normalised or adopted head; otherwise NULL (the run reads the head, as before).
#' @param binding The omic's declaration binding.
#' @param manager The session's state manager.
#' @param normData The tab's reactive values (`post_filter_obj`, the previous run's start).
#' @return The state name to start from, or NULL; refusals have class `omic_norm_runner_busy` or `omic_norm_input_normalised`.
#' @noRd
omicNormRunStart <- function(binding, manager, normData) {
    if (!inherits(manager, "DIAStateAdapter")) return(NULL)
    texts <- omicNormGuardTexts()
    runner <- workflowOmicsManagedActions(binding$omic_type)
    if (!is.null(runner) && !is.null(shiny::isolate(runner$step$active()))) {
        rlang::abort(texts$busy, class = "omic_norm_runner_busy")
    }
    head <- manager$current_state
    if (isTRUE(head %in% omicNormStates(binding))) {
        start <- shiny::isolate(normData$post_filter_obj)
        if (isSavedWorkflowReference(start)) {
            history <- as.character(unlist(manager$getHistory()))
            # The state in the active history whose current reference is the remembered start (its entity unchanged).
            name <- Find(\(state) identical(getWorkflowStateReference(manager, state)$entity_id, start$entity_id), history)
            if (!is.null(name) && !omicNormAdopted(manager, name) && !omicNormalisedLineage(binding, manager, start)) return(name)
        }
        rlang::abort(omicNormRefusalText(binding, manager, texts$run_unknown_lead, texts$run_revert), class = "omic_norm_input_normalised")
    }
    if (omicNormHeadNormalised(binding, manager)) {
        rlang::abort(omicNormRefusalText(binding, manager, texts$run_lead, texts$run_revert), class = "omic_norm_input_normalised")
    }
    NULL
}

#' The run rule applied before a run: a refusal is logged and returned as an error (the caller shows it); a remembered
#' start is selected (the history truncated there, as Revert to it would) and logged.
#' @return TRUE when the run may go on, else the condition.
#' @noRd
omicNormApplyRunStart <- function(binding, manager, normData, addLog) {
    tryCatch({
        start <- omicNormRunStart(binding, manager, normData)
        if (!is.null(start)) {
            selectWorkflowStateReference(manager, start)
            addLog(sprintf("Normalisation runs again from '%s', the table the previous run started from", start))
        }
        TRUE
    }, error = function(e) {
        addLog(paste("ERROR:", conditionMessage(e)))
        e
    })
}

#' KIT-254: before a correlation Apply or Skip at a correlation-step head (one `assertOmicNormHeadIsResult` accepted),
#' the active history returns to the state holding the run's result (selected as Revert to it would, as
#' `omicNormApplyRunStart` selects a remembered start), so the Apply reads, and the Skip commits, the table the filter
#' first received. The correlation-step states after it leave the active history; no state record is deleted. A no-op on
#' a manager that is not a session owner, on a tab whose final is not a saved reference, and at the run's result.
#' @return TRUE when the history was cut, else FALSE; a condition of class `omic_norm_correlation_input_stale` when no
#'   state of the active history holds the run's result.
#' @noRd
omicNormCorrelationStart <- function(manager, normData, addLog) {
    if (!inherits(manager, "DIAStateAdapter")) return(FALSE)
    final <- shiny::isolate(normData$ruv_corrected_obj %||% normData$post_norm_obj)
    if (!isSavedWorkflowReference(final)) return(FALSE)
    if (identical(getWorkflowStateReference(manager)$entity_id, final$entity_id)) return(FALSE)
    history <- as.character(unlist(manager$getHistory()))
    # The state in the active history whose current reference is the run's result (as `omicNormRunStart` finds a start).
    name <- Find(\(state) identical(getWorkflowStateReference(manager, state)$entity_id, final$entity_id), history)
    if (is.null(name)) rlang::abort(omicNormGuardTexts()$correlation, class = "omic_norm_correlation_input_stale")
    selectWorkflowStateReference(manager, name)
    addLog(sprintf("Correlation filter starts again from '%s', the table the normalisation produced", name))
    TRUE
}

#' KIT-254: the tab's correlation fields at their initial values (`omicNormReactiveState`) after a failure that followed
#' the return to the run's result, so the summary and the final PCA show the table the head and the analysis now hold.
#' @noRd
clearOmicNormCorrelation <- function(normData) {
    normData$correlation_results <- list()
    normData$correlation_filtered_obj <- NULL
    normData$correlation_filtering_complete <- FALSE
    invisible(normData)
}

#' The correlation and export guard (7.4, lipidomics and metabolomics): the head must be the last run's result, or a
#' correlation-step state made from it (A01). A no-op on a manager that is not a session owner and on a tab whose final
#' is not a saved reference (the KIT-075 harness form).
#' @param kind `correlation` or `export`: which text a refusal shows.
#' @return TRUE, invisibly; otherwise a condition of class `omic_norm_correlation_input_stale`.
#' @noRd
assertOmicNormHeadIsResult <- function(binding, manager, normData, kind = "correlation") {
    if (!inherits(manager, "DIAStateAdapter")) return(invisible(TRUE))
    final <- shiny::isolate(normData$ruv_corrected_obj %||% normData$post_norm_obj)
    if (!isSavedWorkflowReference(final)) return(invisible(TRUE))
    head <- getWorkflowStateReference(manager)
    if (identical(head$entity_id, final$entity_id)) return(invisible(TRUE))
    if (isTRUE(manager$current_state %in% paste0(binding$prefix, c("_correlation_filtered", "_norm_complete"))) &&
        omicNormDerivesFromCurrentRun(binding, manager, head)) {
        return(invisible(TRUE))
    }
    rlang::abort(omicNormGuardTexts()[[kind]], class = "omic_norm_correlation_input_stale")
}
