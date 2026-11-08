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

# KIT-106 W3 (D2, D3, D4, D7): the shared DA tab (kit-101 design 4.8), from lipid's DA module at 2269a45f
# (mod_lipid_da.R, mod_lipid_da_server.R and mod_lipid_da_{analysis,display,results,export,session}_helpers.R, deleted by
# W3) and trx's (mod_trx_da.R, deleted). `mod_omic_da_ui(id, config)` renders an omic's `da` group from the panel library
# (mod_omic_da_panels.R); `mod_omic_da_server(...)` holds the tab's reactive values (results and state names only:
# KIT-075, 078, 082) and registers the outputs and observers its views use (KIT-110 D3's precedent), in lipid's order.
# The run, the session load, the export and the interactive volcano call the core (utils_workflow_omics_da.R,
# utils_workflow_omics_actions.R): the managed function with the omic's runner, the in-process one without it. The
# interactive volcano is built on demand only (D3), as the binding's `<prefix>_volcano_page`. Functions a test mocks are
# called by name, or read from the configuration built when the server starts. The module frame keeps no locals (the
# testServer masking rule; test-omic-design-module.R). Recorded divergences: the tab returns its reactive values for
# every omic (D2; lipid's returned its outputs, which no caller read); the load's "[D66]" console trace and the export's
# two directory log lines are not written (console only). Recorded divergence from conventions/R-shiny.md: the managed
# runner or the in-process call in place of ExtendedTask is kept, under the existing-code scope rule of conventions/R.md.

#' The DA tab's UI: the layout of an omic's `da` group.
#' @param id The module id.
#' @param config The omic's `da` group (`omicGuiConfig(<omic>)$da`).
#' @return A tag list.
#' @noRd
mod_omic_da_ui <- function(id, config) {
    ns <- shiny::NS(id)
    shiny::tagList(omicDaRenderNode(ns, omicDaLayout(config)))
}

#' The DA tab's server.
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label The tab's arguments, as `mod_lipid_da_server`'s.
#' @param config The omic's `da` group (`omicGuiConfig(<omic>)$da`); the binding is the configuration of `omic_type`.
#' @param selected_tab The shell's selected tab (proteomics' plugins, W5); not read for lipid, metab and trx.
#' @return The tab's reactive values.
#' @noRd
mod_omic_da_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, config, selected_tab = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        if (!is.null(config$plugins)) {
            return(runOmicDaPluginServer(input, output, session, workflow_data, experiment_paths, omic_type, experiment_label,
                config, selected_tab))
        }
        runOmicDaServer(input, output, session, workflow_data, experiment_paths, omic_type, config)
    })
}

#' The tab's reactive values (mod_lipid_da.R:34-48 at 2269a45f): the results list without its source object, the state
#' the session load named, the contrasts and assays it restored and the heatmap's last drawing.
#' @noRd
omicDaTabState <- function() {
    shiny::reactiveValues(da_results_list = NULL, contrasts_available = NULL, assays_available = NULL, analysis_complete = FALSE,
        current_state_name = NULL, formula_from_s4 = NULL, current_row_clusters = NULL, current_col_clusters = NULL,
        current_heatmap_plot = NULL)
}

#' The server's body: the observers and outputs in lipid's order (mod_lipid_da_server.R:9-126 at 2269a45f), each made
#' when the omic's group uses it: the model's output when the model is shown, the session load with its test id, the
#' save note, heatmap, cluster summary and Save Heatmap with the heatmap tab.
#' @noRd
runOmicDaServer <- function(input, output, session, workflowData, experimentPaths, omicType, da) {
    tab <- omicDaTabState()
    binding <- omicGuiConfig(omicType)$binding$binding
    registerOmicDaTextOutputs(output, workflowData, tab, da, binding)
    if ("load_session" %in% names(da$testids)) {
        shiny::observeEvent(input$load_filtered_session, loadOmicDaSession(session, workflowData, tab, experimentPaths, da, binding))
    }
    shiny::observeEvent(input$run_da_analysis, runOmicDaAnalysis(input, session, workflowData, tab, experimentPaths, da, binding))
    if (isTRUE(da$views$heatmap$save_note)) {
        output$heatmap_manual_save_warning <- shiny::renderUI(omicDaHeatmapSaveNote(tab$analysis_complete))
    }
    registerOmicDaVolcanoGlimmaOutput(output, input, tab, workflowData, experimentPaths, da, binding)
    registerOmicDaVolcanoStaticOutput(output, input, tab, da)
    if ("heatmap" %in% names(da$tabs)) {
        registerOmicDaHeatmapPlotOutput(output, input, tab, workflowData, da, binding)
        output$cluster_summary <- shiny::renderPrint({
            shiny::req(input$heatmap_tree_cut_method != "none")
            cat(omicDaClusterSummaryText(tab$current_row_clusters, da$nouns[["cluster"]]))
        })
        shiny::observeEvent(input$save_heatmap, saveOmicDaHeatmap(input, tab, experimentPaths, da))
    }
    registerOmicDaTableOutputs(output, input, tab, da, binding)
    invisible(tab)
}

# ------------------------------------------------------------------------------------------------
# The outputs.
# ------------------------------------------------------------------------------------------------

#' The contrasts, the status and, when the model is shown, the model's text. The contrasts and the model read the bus
#' (an environment for the MCI fixtures: computed once per session).
#' @noRd
registerOmicDaTextOutputs <- function(output, workflowData, tab, da, binding) {
    output$contrasts_display <- shiny::renderPrint(cat(omicDaContrastsText(workflowData$contrasts_tbl, da$texts)))
    output$da_status <- shiny::renderPrint(cat(omicDaStatusText(tab$analysis_complete, tab$da_results_list,
        omicDaLoadedMetadata(tab, workflowData, binding), da$texts)))
    if (!isTRUE(da$model$editable)) {
        output[[da$model$output]] <- shiny::renderText({
            workflowData$tab_status
            da$model$text(workflowData)
        })
    }
    invisible(output)
}

#' The interactive volcano, built when its box is ticked (D3): the managed factory with the omic's runner, the in-process
#' page otherwise, each made at render time and given the contrast's entry when the omic keys its results by entry (trx).
#' No page is asked for "Combined".
#' @noRd
registerOmicDaVolcanoGlimmaOutput <- function(output, input, tab, workflowData, experimentPaths, da, binding) {
    output$volcano_glimma <- shiny::renderUI({
        shiny::req(isTRUE(input$volcano_interactive))
        omicDaVolcanoGlimmaUi(tab$da_results_list, input$volcano_contrast, input$volcano_assay, input$da_q_val_thresh,
            omicDaVolcanoPageFactory(workflowData, experimentPaths, da, binding))
    })
    invisible(output)
}

#' The volcano page factory of the session: the worker's page when managed, the in-process page otherwise.
#' @noRd
omicDaVolcanoPageFactory <- function(workflowData, experimentPaths, da, binding) {
    runner <- workflowOmicsManagedActions(binding$omic_type)
    factory <- if (is.null(runner)) {
        workflowOmicsInProcessVolcanoPageFactory(binding, workflowData$state_manager, experimentPaths)
    } else {
        workflowOmicsVolcanoPageFactory(binding, runner, workflowData$state_manager)
    }
    if (is.null(da$contrast_entry)) return(factory)
    function(da_results_list, selected_contrast, selected_assay, da_q_val_thresh) {
        factory(da_results_list, da$contrast_entry(workflowData, selected_contrast), selected_assay, da_q_val_thresh)
    }
}

#' The static volcano.
#' @noRd
registerOmicDaVolcanoStaticOutput <- function(output, input, tab, da) {
    output[[da$views$static_volcano$output]] <- shiny::renderPlot(omicDaStaticVolcano(tab$da_results_list, input$volcano_contrast,
        input$volcano_assay, input$da_q_val_thresh, input$treat_lfc_cutoff, da$views$static_volcano))
    invisible(output)
}

#' The heatmap, drawn on the tab's rows with the rows it chooses read from the state (KIT-082); its clusters and drawing
#' kept for the cluster summary and Save Heatmap.
#' @noRd
registerOmicDaHeatmapPlotOutput <- function(output, input, tab, workflowData, da, binding) {
    output$heatmap_plot <- shiny::renderPlot({
        state <- omicDaHeatmapState(omicsDaResultsWithRowReader(binding, tab$da_results_list, tab, workflowData),
            input$heatmap_contrast, input$heatmap_assay, input$heatmap_top_n, input$heatmap_clustering, input$heatmap_scaling,
            input$heatmap_cluster_method, input$heatmap_distance_method, input$heatmap_color_scheme, input$heatmap_show_labels,
            input$da_q_val_thresh, input$heatmap_tree_cut_method, input$heatmap_n_clusters, input$heatmap_cut_height,
            input$heatmap_min_cluster_size, da$views$heatmap)
        omicDaStoreHeatmap(state, tab)
    })
    invisible(output)
}

#' Keep a heatmap's clusters and drawing on the tab (mod_lipid_da_display_helpers.R:202-212 at 2269a45f).
#' @return The plot to show, or NULL.
#' @noRd
omicDaStoreHeatmap <- function(state, tab) {
    if (is.null(state)) return(NULL)
    tab$current_row_clusters <- state$rowClusters
    tab$current_col_clusters <- state$colClusters
    tab$current_heatmap_plot <- state$currentHeatmapPlot
    state$plot
}

#' The summary statistics, the results table and the download, on the binding's long table.
#' @noRd
registerOmicDaTableOutputs <- function(output, input, tab, da, binding) {
    output$da_summary_stats <- shiny::renderPrint({
        results <- tab$da_results_list
        shiny::req(results)
        cat(omicDaSummaryText(results[[binding$long_table]], input$table_contrast, input$table_assay, input$da_q_val_thresh,
            da$summary, da$table_filter))
    })
    output$da_results_table <- DT::renderDT({
        results <- tab$da_results_list
        shiny::req(results)
        omicDaTableWidget(results[[binding$long_table]], input$table_contrast, input$table_assay, input$table_significance,
            input$da_q_val_thresh, input$treat_lfc_cutoff, input$table_max_rows, da$columns, da$table_filter)
    })
    registerOmicDaDownloadOutput(output, tab, da, binding)
}

#' The download: the tab's results read at each click.
#' @noRd
registerOmicDaDownloadOutput <- function(output, tab, da, binding, ...) {
    output$download_da_results <- omicDaDownloadHandler(function() tab$da_results_list, binding$long_table, da$download, ...)
    invisible(output)
}

# ------------------------------------------------------------------------------------------------
# The source the DA runs on.
# ------------------------------------------------------------------------------------------------

#' The metadata of the state the session load named (KIT-075: never the object), or NULL when none is named.
#' @noRd
omicDaLoadedMetadata <- function(tab, workflowData, binding) {
    state <- tab$current_state_name
    if (is.null(state) || is.null(workflowData) || is.null(workflowData$state_manager)) return(NULL)
    readWorkflowOmicsMetadata(binding, workflowData$state_manager, state)
}

#' The metadata of the state the DA would run on when the session load named none: the group's reader, given the bus
#' (trx's `trxDaReadMetadata`) or the owner (the binding readers, on the head). F106-05: a head that cannot be read (the
#' owner's `initial` state) is no data, and the tab refuses as for no owner.
#' @noRd
omicDaReadSource <- function(workflowData, source) {
    read <- source$read_metadata
    tryCatch(if ("workflowData" %in% names(formals(read))) read(workflowData) else read(workflowData$state_manager),
        error = function(error) {
            logger::log_warn(sprintf("   The DA tab could not read the head state: %s", conditionMessage(error)))
            NULL
        })
}

# ------------------------------------------------------------------------------------------------
# The run.
# ------------------------------------------------------------------------------------------------

#' The run button (mod_lipid_da_analysis_helpers.R:12-138 at 2269a45f): the refusals (no data, no contrasts), the
#' running notice and the in-process preflight (before the notice when the group says `preflight_first`, D7), then the
#' run.
#' @noRd
runOmicDaAnalysis <- function(input, session, workflowData, tab, experimentPaths, da, binding) {
    settings <- list(formula_string = if (isTRUE(da$model$editable)) input$formula_string,
        da_q_val_thresh = input$da_q_val_thresh, treat_lfc_cutoff = input$treat_lfc_cutoff)
    logger::log_info("=== RUN DA ANALYSIS BUTTON CLICKED ===")
    current <- omicDaLoadedMetadata(tab, workflowData, binding) %||% omicDaReadSource(workflowData, da$source)
    if (!is.list(current) || !identical(current$class, da$source$class)) {
        shiny::showNotification(da$source$no_data, type = "error", duration = 5)
        return(invisible(NULL))
    }
    contrasts <- workflowData$contrasts_tbl
    if (is.null(contrasts) || nrow(contrasts) == 0) {
        shiny::showNotification("No contrasts defined. Please define contrasts in the design tab.", type = "error", duration = 5)
        return(invisible(NULL))
    }
    first <- isTRUE(da$preflight_first)
    refused <- if (first) omicDaPreflight(workflowData, current, contrasts, settings, da, binding, running = FALSE)
    if (!is.null(refused)) return(invisible(refused))
    shiny::showNotification("Running differential expression analysis...", id = "da_running", duration = NULL)
    refused <- if (!first) omicDaPreflight(workflowData, current, contrasts, settings, da, binding, running = TRUE)
    if (!is.null(refused)) return(invisible(refused))
    executeOmicDaRun(tab, workflowData, session, experimentPaths, settings, contrasts, da, binding)
}

#' The group's preflight in process (mod_lipid_da_analysis_helpers.R:84-115 at 2269a45f; under a managed session the
#' declaration's own preflight runs in the worker): on the state read for this action only (KIT-075). An invalid
#' preflight is refused with its errors, the running notice removed when it was shown.
#' @return NULL when the run may go on, else `list(status = "error", stage = "preflight", errorMessage)`.
#' @noRd
omicDaPreflight <- function(workflowData, current, contrasts, settings, da, binding, running) {
    if (!is.function(da$preflight) || !is.null(workflowOmicsManagedActions(binding$omic_type))) return(NULL)
    preflight <- da$preflight(currentS4 = readWorkflowOmicsState(binding, workflowData$state_manager, current$reference),
        contrastsTbl = contrasts, formulaString = settings$formula_string)
    if (isTRUE(preflight$valid)) return(NULL)
    message <- paste(preflight$errors, collapse = "; ")
    if (running) shiny::removeNotification("da_running")
    shiny::showNotification(message, type = "error", duration = 5)
    list(status = "error", stage = "preflight", errorMessage = message)
}

#' The DA on the tab's state (mod_lipid_da_analysis_helpers.R:140-214 at 2269a45f): the managed DA with the omic's runner
#' (the realised formula and the contrast entries when the group has them), else the in-process reference DA; then the
#' completion, or the error.
#' @param settings `formula_string`, `da_q_val_thresh`, `treat_lfc_cutoff`.
#' @return `list(status)`: "dispatched", "success" (with `results`) or "error" (with `errorMessage`).
#' @noRd
executeOmicDaRun <- function(tab, workflowData, session, experimentPaths, settings, contrastsTbl, da, binding) {
    tryCatch({
        contrasts <- if (is.null(da$contrast_entries)) contrastsTbl else da$contrast_entries(contrastsTbl)
        manager <- workflowData$state_manager
        source <- getWorkflowStateReference(manager, tab$current_state_name)
        runner <- workflowOmicsManagedActions(binding$omic_type)
        if (!is.null(runner)) {
            if (is.function(da$model$formula)) settings$formula_string <- da$model$formula(manager, source)
            da$run$managed(binding, runner, manager, source, settings, contrasts,
                deliver = function(results, export) {
                    completeOmicDaRun(results, tab, workflowData, session, experimentPaths, da, binding,
                        export = if (is.function(da$run$export_report)) da$run$export_report(export) else export)
                },
                failed = function(error) failOmicDaRun(conditionMessage(error)),
                invalidate = function() {
                    tab$da_results_list <- NULL
                    tab$analysis_complete <- FALSE
                    status <- workflowData$tab_status
                    status[[da$status_key]] <- "pending"
                    workflowData$tab_status <- status
                })
            return(list(status = "dispatched"))
        }
        results <- da$run$in_process(manager, source, settings, contrasts)$results
        completeOmicDaRun(results, tab, workflowData, session, experimentPaths, da, binding)
        list(status = "success", results = results)
    }, error = function(e) {
        failOmicDaRun(conditionMessage(e))
        list(status = "error", errorMessage = conditionMessage(e))
    })
}

#' A failed run: the running notice removed and the error shown.
#' @noRd
failOmicDaRun <- function(message) {
    logger::log_error(sprintf("   DA analysis error: %s", message))
    shiny::removeNotification("da_running")
    shiny::showNotification(sprintf("Analysis error: %s", message), type = "error", duration = 10)
}

#' The completion (mod_lipid_da_analysis_helpers.R:216-360 at 2269a45f): the results on the tab, the status, the notices,
#' the selectors on the long table's contrasts and assays (one assay: the volcano opens on it, D7), then the export: the
#' worker's when managed (`export`), else the binding's DA publication run here (D4).
#' @param export The managed export's outcome (after the group's `run$export_report`), or NULL in process.
#' @return The selectors' choices and the export's report.
#' @noRd
completeOmicDaRun <- function(results, tab, workflowData, session, experimentPaths, da, binding, export = NULL) {
    long_table <- binding$long_table
    # D2: an omic whose results are keyed by the contrast's name (trx's entries) names each row by its comparison.
    if (is.function(da$contrast_entry) && !is.null(results[[long_table]])) {
        results[[long_table]]$friendly_name <- results[[long_table]]$comparison
    }
    tab$da_results_list <- results
    tab$analysis_complete <- TRUE
    status <- workflowData$tab_status
    status[[da$status_key]] <- "complete"
    workflowData$tab_status <- status
    shiny::removeNotification("da_running")
    shiny::showNotification("Differential expression analysis complete!", type = "message", duration = 5)
    logger::log_info("   DA analysis completed successfully")
    # KIT-076: a failed contrast fails alone; the run says so instead of a table that is silently missing.
    if (length(results$failed_contrasts) > 0) {
        shiny::showNotification(sprintf("DA failed for: %s", paste(names(results$failed_contrasts), collapse = ", ")),
            type = "warning", duration = NULL)
    }
    choices <- list()
    long <- results[[long_table]]
    if (!is.null(long)) {
        contrasts <- unique(long$friendly_name)
        if (length(contrasts) == 0) contrasts <- unique(long$comparison)
        for (id in c("volcano_contrast", "heatmap_contrast", "table_contrast")) {
            shiny::updateSelectInput(session, id, choices = contrasts, selected = contrasts[1])
        }
        single <- identical(da$assays, "single")
        assays <- unique(long$assay)
        if (is.null(assays) && single) assays <- names(results$significant_counts)
        shiny::updateSelectInput(session, "volcano_assay", choices = c("Combined", assays),
            selected = if (single) assays[[1L]] else "Combined")
        shiny::updateSelectInput(session, "heatmap_assay", choices = c("Combined", assays), selected = "Combined")
        shiny::updateSelectInput(session, "table_assay", choices = c("All", unique(long$assay)), selected = "All")
        logger::log_info(sprintf("   Updated dropdowns: %d contrasts, %d assays", length(contrasts), length(assays)))
        choices <- list(contrastChoices = contrasts, assayChoices = c("Combined", assays),
            tableAssayChoices = c("All", unique(long$assay)))
    }
    logger::log_info("   Writing DA results to disk...")
    export <- if (!is.null(export)) reportWorkflowOmicsDaExport(export) else exportOmicDaInProcess(workflowData, experimentPaths, da, binding)
    invisible(c(choices, list(export = export)))
}

#' The in-process export (D4): the binding's `<prefix>_da_publication` run here on the head, reported as the worker's.
#' Unset output directories skip it (mod_lipid_da_analysis_helpers.R:347-349 at 2269a45f).
#' @noRd
exportOmicDaInProcess <- function(workflowData, experimentPaths, da, binding) {
    if (is.null(experimentPaths$da_output_dir) || is.null(experimentPaths$publication_graphs_dir)) {
        logger::log_warn("   Output directories not configured, skipping file output")
        return(list(status = "skipped", success = FALSE))
    }
    display <- tryCatch(runWorkflowOmicsPublicationInProcess(binding, workflowData$state_manager, paste0(binding$prefix, "_da_publication"),
        list(), experimentPaths), error = function(e) list(success = FALSE, message = conditionMessage(e)))
    if (is.function(da$run$export_report)) display <- da$run$export_report(display)
    reportWorkflowOmicsDaExport(display, shiny::showNotification)
}

# ------------------------------------------------------------------------------------------------
# The session load and Save Heatmap.
# ------------------------------------------------------------------------------------------------

#' The load button: the core's session load (`loadWorkflowOmicsDaSession`, in the module's reactive domain for the
#' legacy adoption), then the selectors, the formula and the notices (mod_lipid_da_session_helpers.R at 2269a45f).
#' @noRd
loadOmicDaSession <- function(session, workflowData, tab, experimentPaths, da, binding) {
    logger::log_info("=== LOAD FILTERED SESSION BUTTON CLICKED ===")
    res <- loadWorkflowOmicsDaSession(binding, workflowData, tab, experimentPaths, da$session_load,
        reading = function() shiny::showNotification("Loading filtered session...", id = "loading_session", duration = NULL))
    if (identical(res$status, "missing")) {
        shiny::showNotification(res$message, type = "error", duration = 5)
        return(invisible(res))
    }
    if (!is.null(res$contrast_choices)) {
        for (id in c("volcano_contrast", "heatmap_contrast", "table_contrast")) {
            shiny::updateSelectInput(session, id, choices = res$contrast_choices, selected = res$contrast_choices[1])
        }
    }
    if (!is.null(res$assay_choices)) {
        shiny::updateSelectInput(session, "volcano_assay", choices = res$assay_choices)
        shiny::updateSelectInput(session, "heatmap_assay", choices = res$assay_choices)
        shiny::updateSelectInput(session, "table_assay", choices = res$table_assay_choices)
    }
    if (!is.null(res$formula)) shiny::updateTextAreaInput(session, "formula_string", value = res$formula)
    shiny::removeNotification("loading_session")
    if (identical(res$status, "success")) {
        shiny::showNotification("Session loaded successfully!", type = "message", duration = 3)
    } else {
        shiny::showNotification(paste("Error loading session:", res$message), type = "error", duration = 10)
    }
    invisible(res)
}

#' Save Heatmap (mod_lipid_da_export_helpers.R at 2269a45f): the heatmap last drawn, its clusters and its parameters,
#' under the group's file prefix.
#' @noRd
saveOmicDaHeatmap <- function(input, tab, experimentPaths, da) {
    plot <- tab$current_heatmap_plot
    clusters <- tab$current_row_clusters
    shiny::req(plot, experimentPaths$publication_graphs_dir)
    logger::log_info("Save Heatmap button clicked")
    params <- list(contrast = input$heatmap_contrast, top_n = input$heatmap_top_n, clustering_method = input$heatmap_cluster_method,
        distance_method = input$heatmap_distance_method, cluster_rows = input$heatmap_clustering %in% c("both", "row"),
        cluster_cols = input$heatmap_clustering %in% c("both", "column"), scaling = input$heatmap_scaling,
        color_scheme = input$heatmap_color_scheme, tree_cut_method = input$heatmap_tree_cut_method,
        n_clusters = input$heatmap_n_clusters, cut_height = input$heatmap_cut_height,
        min_cluster_size = input$heatmap_min_cluster_size)
    prefix <- gsub("[^A-Za-z0-9_]", "_", paste0(da$views$heatmap$save_prefix, input$heatmap_contrast))
    save_heatmap_products(heatmap_obj = plot, row_clusters = clusters, params_list = params,
        output_dir = experimentPaths$publication_graphs_dir, file_prefix = prefix)
    shiny::showNotification("Heatmap and cluster info saved to publication_graphs/Heatmap", type = "message", duration = 5)
    invisible(list(params = params, prefix = prefix))
}

# ------------------------------------------------------------------------------------------------
# The plugin registry: an omic whose group names `plugins` (proteomics, KIT-106 W5, D5).
# ------------------------------------------------------------------------------------------------

#' The tab's reactive values when the group has plugins: proteomics' (mod_prot_da_server.R:21-30 at fe97aac1), which
#' its plugins read and write.
#' @noRd
omicDaPluginTabState <- function() {
    shiny::reactiveValues(da_results_list = NULL, contrasts_available = NULL, analysis_complete = FALSE,
        volcano_plot_identity = NULL, current_s4_object = NULL, formula_from_s4 = NULL, current_row_clusters = NULL,
        current_col_clusters = NULL)
}

#' The server of an omic whose group has plugins (mod_prot_da_server.R:9-57 at fe97aac1): its `init`, `load` and `run`
#' plugins, each looked up by name when called (`omicGuiFunction`), with their signatures unchanged; the frame's
#' registrations of the group's outputs between them, in the order the omic's server made them. Kernel placement is the
#' run plugin's (KIT-111 to 113). Recorded divergence: the server's six console `cat()` lines on entry are not written.
#' @return The tab's reactive values.
#' @noRd
runOmicDaPluginServer <- function(input, output, session, workflowData, experimentPaths, omicType, experimentLabel, da,
    selectedTab) {
    tab <- omicDaPluginTabState()
    omicGuiFunction(da$plugins[["init"]])(input, output, session, tab, workflowData, selectedTab)
    registerOmicDaPluginTextOutputs(output, tab, da)
    omicGuiFunction(da$plugins[["load"]])(input, output, session, tab, workflowData, experimentPaths)
    omicGuiFunction(da$plugins[["run"]])(input, output, session, session$ns, tab, workflowData, experimentPaths,
        omic_type = omicType, experiment_label = experimentLabel)
    registerOmicDaPluginVolcanoOutputs(output, input, tab, experimentPaths, da)
    registerOmicDaPluginHeatmapOutputs(output, input, tab, experimentPaths, da)
    registerOmicDaPluginTableOutputs(output, input, tab, da)
    tab
}

#' The contrasts and the status, as text (mod_prot_da_handlers_init.R:318-336 at fe97aac1).
#' @noRd
registerOmicDaPluginTextOutputs <- function(output, tab, da) {
    output$contrasts_display <- shiny::renderText(diaDaContrastsView(tab$contrasts_available, da$texts))
    output$da_status <- shiny::renderText(diaDaStatusView(tab$analysis_complete, tab$contrasts_available, da$texts))
    invisible(output)
}

#' The volcano (mod_prot_da_render_volcano.R at fe97aac1): the opt-in CP08 capture, the interactive volcano built when
#' its box is ticked (D3; the DA-time page), and the static volcano, whose alt text is the managed publication's plot
#' identity (KIT-019).
#' @noRd
registerOmicDaPluginVolcanoOutputs <- function(output, input, tab, experimentPaths, da) {
    # --- TESTTHAT CHECKPOINT CP08 (see test-prot-08-volcano.R) ---
    shiny::observeEvent(list(input$volcano_contrast, tab$da_results_list), {
        shiny::req(input$volcano_contrast, tab$da_results_list)
        .capture_checkpoint(list(
            da_results_list = tab$da_results_list,
            selected_contrast = input$volcano_contrast,
            da_q_val_thresh = input$da_q_val_thresh,
            args_row_id = diaGuiProteinKey(tab$current_s4_object)
        ), "cp08", "volcano_input")
    }, ignoreInit = TRUE)
    output$volcano_glimma <- shiny::renderUI({
        shiny::req(isTRUE(input$volcano_interactive))
        diaDaVolcanoGlimmaView(tab$da_results_list, input$volcano_contrast, input$da_q_val_thresh, tab$current_s4_object,
            experimentPaths$da_output_dir)
    })
    output[[da$views$static_volcano$output]] <- shiny::renderPlot(alt = function() tab$volcano_plot_identity %||% NA_character_,
        diaDaStaticVolcanoView(tab$da_results_list, input$volcano_contrast, input$da_q_val_thresh, input$treat_lfc_cutoff,
            input$volcano_show_labels, input$volcano_label_top_n))
    invisible(output)
}

#' The heatmap (mod_prot_da_render_heatmap.R at fe97aac1): the opt-in CP09 capture, the heatmap, the cluster summary
#' and the frame's Save Heatmap under the group's prefix.
#' @noRd
registerOmicDaPluginHeatmapOutputs <- function(output, input, tab, experimentPaths, da) {
    # --- TESTTHAT CHECKPOINT CP09 (see test-prot-09-heatmap.R) ---
    shiny::observeEvent(list(input$heatmap_contrast, tab$da_results_list), {
        shiny::req(input$heatmap_contrast, tab$da_results_list)
        plot_data_structure <- if (isDIAResultManifest(tab$da_results_list)) tab$da_results_list else list(
            da_proteins_long = tab$da_results_list$da_proteins_long,
            theObject = tab$da_results_list$theObject
        )
        .capture_checkpoint(list(
            da_results_list = plot_data_structure,
            selected_contrast = input$heatmap_contrast,
            top_n_genes = input$heatmap_top_n,
            clustering_method = input$heatmap_cluster_method,
            distance_method = input$heatmap_distance_method,
            cluster_rows = input$heatmap_clustering %in% c("both", "row"),
            cluster_cols = input$heatmap_clustering %in% c("both", "column"),
            scale_data = input$heatmap_scaling,
            color_scheme = input$heatmap_color_scheme,
            show_gene_names = input$heatmap_show_labels,
            da_q_val_thresh = input$da_q_val_thresh,
            tree_cut_method = input$heatmap_tree_cut_method,
            n_clusters = input$heatmap_n_clusters,
            cut_height = input$heatmap_cut_height,
            min_cluster_size = input$heatmap_min_cluster_size
        ), "cp09", "heatmap_input")
    }, ignoreInit = TRUE)
    output$heatmap_plot <- shiny::renderPlot(diaDaHeatmapView(tab, input$heatmap_contrast, input$heatmap_top_n,
        input$heatmap_cluster_method, input$heatmap_distance_method, input$heatmap_clustering, input$heatmap_scaling,
        input$heatmap_color_scheme, input$heatmap_show_labels, input$da_q_val_thresh, input$heatmap_tree_cut_method,
        input$heatmap_n_clusters, input$heatmap_cut_height, input$heatmap_min_cluster_size, input$heatmap_show_dendro))
    output$cluster_summary <- shiny::renderPrint(diaDaClusterSummaryView(input$heatmap_tree_cut_method, tab$current_row_clusters))
    shiny::observeEvent(input$save_heatmap, saveOmicDaHeatmap(input, tab, experimentPaths, da))
    invisible(output)
}

#' The results table paged from the store, the summary statistics and the download (mod_prot_da_render_table.R at
#' fe97aac1).
#' @noRd
registerOmicDaPluginTableOutputs <- function(output, input, tab, da) {
    pager <- new.env(parent = emptyenv())
    output$da_results_table <- DT::renderDT(diaDaResultsTableView(tab$da_results_list, input$table_contrast,
        input$table_significance, input$da_q_val_thresh, input$treat_lfc_cutoff, input$table_max_rows, pager),
        server = TRUE, funcFilter = diaDisplayTableFilter(pager))
    output$da_summary_stats <- shiny::renderText(diaDaSummaryStatsView(tab$da_results_list, input$table_contrast,
        input$da_q_val_thresh, input$treat_lfc_cutoff))
    output$download_da_results <- shiny::downloadHandler(
        filename = function() paste0(da$download$stem, Sys.Date(), ".", da$download$extension),
        content = diaDaDownloadContent)
    invisible(output)
}
