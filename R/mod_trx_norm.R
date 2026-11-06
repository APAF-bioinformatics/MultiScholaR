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

# KIT-097: the transcriptomics normalisation tab, the lipid norm module under the transcriptomics profile. The panels are
# the lipid tab's (RUV QC, PCA, RLE, Log; one assay row) with its renderers; the options and the sequence are the
# transcriptomics ones: TMM and logCPM (`trx_tmm_logcpm`), the RUV grid (`trx_ruv_grid`, automatic as the notebook or
# manual), RUVs (`trx_ruvs`, the DA head), then the PCA and RLE at raw, TMM and RUVs, the CCA panels and the stage grid.
# `runWorkflowOmicsManagedNormalisation` (ITSD, log2, between-sample, RUV-III) is not this sequence and is not reused.

#' The normalisation options of the transcriptomics profile (the lipid panel's layout and styling).
#' @noRd
trxNormOptionsPanel <- function(ns) {
    automatic <- shiny::conditionalPanel(condition = "input.ruv_mode == 'automatic'", ns = ns,
        shiny::h5("Automatic Optimization Parameters"),
        shiny::fluidRow(shiny::column(4, shiny::numericInput(ns("auto_percentage_min"), "Min %:", value = 5, min = 1, max = 50, width = "100%")),
            shiny::column(4, shiny::numericInput(ns("auto_percentage_max"), "Max %:", value = 30, min = 1, max = 50, width = "100%")),
            shiny::column(4, shiny::numericInput(ns("auto_percentage_step"), "Step %:", value = 5, min = 1, max = 25, width = "100%"))),
        shiny::numericInput(ns("k_penalty_weight"), "K Value Penalty Weight:", value = 0.05, min = 0, max = 0.9, step = 0.01, width = "100%"),
        shiny::helpText("The notebook's grid: 5 to 30 % by 5, k capped at 3 below 15 samples and 4 otherwise, penalty 0.05; the",
            "percentage and k with the maximum composite score are chosen."))
    manual <- shiny::conditionalPanel(condition = "input.ruv_mode == 'manual'", ns = ns,
        shiny::sliderInput(ns("ruv_percentage"), "% Genes as Inferred Negative Controls:", min = 1, max = 50, value = 10, step = 1, width = "100%"),
        shiny::numericInput(ns("ruv_k"), "Number of Factors (k):", value = 2, min = 1, max = 4, width = "100%"),
        shiny::helpText("k must not exceed the cap (3 below 15 samples, 4 otherwise); a larger k is refused."))
    shiny::column(3, shiny::wellPanel(
        shiny::h4("Normalization Options"), shiny::hr(),
        shiny::h4("TMM and logCPM"),
        shiny::helpText("edgeR TMM normalisation factors of the filtered counts; log2 counts per million with a prior count of 2."),
        shiny::hr(), shiny::h4("RUVs (RUVSeq)"),
        testid(shiny::radioButtons(ns("ruv_mode"), "RUV Parameter Tuning:", choices = list("Automatic (recommended)" = "automatic",
            "Manual (advanced users)" = "manual"), selected = "automatic", width = "100%"), "trx-norm-ruv-mode"),
        shiny::helpText("RUVs uses inferred negative controls (the least significant genes in a preliminary group GLM). Its factors",
            "W_1..W_k enter the DA model as covariates; the RUVs-corrected logCPM is for display only."),
        automatic, manual, shiny::hr(),
        testid(shiny::actionButton(ns("run_normalization"), "Run Normalization Pipeline", class = "btn-primary", width = "100%",
            icon = shiny::icon("play")), "trx-norm-run"),
        shiny::br(), shiny::br(),
        testid(shiny::actionButton(ns("reset_normalization"), "Reset to Pre-Normalization", class = "btn-warning", width = "100%",
            icon = shiny::icon("undo")), "trx-norm-reset")))
}

#' The grid's parameters from the tab's inputs (KIT-095's declared names).
#' @noRd
trxNormGridParameters <- function(input) {
    number <- function(value, default) if (is.null(value) || !length(value) || is.na(value)) default else as.numeric(value)
    if (identical(input$ruv_mode, "manual")) {
        return(trxRuvGridParameters("manual", manual_k = as.integer(number(input$ruv_k, NA)), manual_percentage = number(input$ruv_percentage, NA)))
    }
    trxRuvGridParameters("automatic", number(input$auto_percentage_min, 5), number(input$auto_percentage_max, 30),
        number(input$auto_percentage_step, 5), number(input$k_penalty_weight, 0.05))
}

#' What the lipid RUV QC renderers show for the one count table, from the grid product: the CCA panels, k, percentage,
#' separation, the number of inferred controls and the grid table. No count is held.
#' @noRd
trxNormRuvQcResult <- function(grid) {
    panels <- trxCcaPanels(grid$best, grid$grid_summary, length(grid$prelim_order))
    list(success = TRUE, cancor_plot = panels$grid$plot, best_k = grid$chosen$k, best_percentage = grid$chosen$percentage,
        separation_score = grid$best$separation_score, control_genes_index = rep(TRUE, length(grid$chosen$control_genes)),
        optimization_results = grid$grid_summary)
}

# ---------------------------------------------------------------------------------------------------------------
# The publication operations in process.
# ---------------------------------------------------------------------------------------------------------------

#' Run one of the binding's publication operations in process on the head and copy what it staged to the experiment's
#' directories (what `publishWorkflowOmicsResult` does for the worker).
#' @param operation The registered name (`trx_qc_figures`, ...).
#' @return The operation's display.
#' @noRd
runTrxPublicationInProcess <- function(manager, operation, parameters, experiment_paths) {
    entry <- Filter(function(p) identical(paste0("trx_", p[[1L]]), operation), trxRoutePublications)
    if (!length(entry)) cli::cli_abort("unknown transcriptomics publication {.val {operation}}", class = "workflow_operation_parameters")
    directory <- tempfile("trx-publication-")
    dir.create(directory, recursive = TRUE)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    request <- list(operation = operation, directory = directory, parameters = parameters,
        input_references = list(state = getWorkflowStateReference(manager)))
    done <- withWorkflowStateScope(manager, function() get(entry[[1L]][[2L]], mode = "function")(manager, request, NULL))
    staged <- file.path(directory, "omics-publication")
    for (role in list.dirs(staged, full.names = FALSE, recursive = FALSE)) {
        target <- experiment_paths[[role]]
        if (is.null(target)) next
        dir.create(target, recursive = TRUE, showWarnings = FALSE)
        file.copy(list.files(file.path(staged, role), full.names = TRUE, all.files = TRUE, no.. = TRUE), target, recursive = TRUE, overwrite = TRUE)
    }
    done$display
}

#' The figures and tables after RUVs (the stages' PCA and RLE, the CCA panels, the stage grid), in process.
#' @noRd
runTrxNormPublicationsInProcess <- function(manager, experiment_paths, addLog) {
    for (stage in c("raw", "tmm", "ruvs")) {
        runTrxPublicationInProcess(manager, "trx_qc_figures", list(stage = stage), experiment_paths)
        addLog(sprintf("PCA and RLE (%s) generated", stage))
    }
    runTrxPublicationInProcess(manager, "trx_ruv_publication", list(), experiment_paths)
    runTrxPublicationInProcess(manager, "trx_composite_figure", list(figure = "normalization_comparison"), experiment_paths)
    addLog("CCA panels, grid workbook and the normalisation comparison grid saved")
}

# ---------------------------------------------------------------------------------------------------------------
# The sequence.
# ---------------------------------------------------------------------------------------------------------------

#' The RUVs step on the products the tab ran (the commit of the managed sequence).
#' @noRd
trxNormRuvsStep <- function(logcpm, grid, k) {
    list(operation = "trx_ruvs", parameters = list(logcpm = logcpm, grid = grid), state_name = "trx_ruvs",
        description = sprintf("RUVs: W_1..W_%d added to the design", k))
}

#' The managed sequence: the raw figures and TMM/logCPM, then the grid, then the RUVs commit and its figures. Each later
#' action names the products the earlier ones delivered, so each runs once those have.
#' @noRd
runTrxNormManaged <- function(runner, grid_parameters, config, normData, addLog, finish, failed) {
    held <- new.env(parent = emptyenv())
    figures <- function(stage) omicsManagedProduct("trx_qc_figures", list(stage = stage), function(result) addLog(sprintf("PCA and RLE (%s) generated", stage)))
    after_grid <- function() runner$run(list(
        omicsManagedCommit(trxNormRuvsStep(held$logcpm, held$grid, held$k), config, function(result) {
            normData$ruv_corrected_obj <- result$stateReference
            addLog(sprintf("RUVs applied: W_1..W_%d added to the design (state 'trx_ruvs')", held$k))
        }),
        figures("tmm"), figures("ruvs"), omicsManagedProduct("trx_ruv_publication", list(), function(result) invisible(NULL)),
        omicsManagedProduct("trx_composite_figure", list(figure = "normalization_comparison"),
            function(result) addLog("CCA panels, grid workbook and the normalisation comparison grid saved"))), done = finish, failed = failed)
    after_logcpm <- function() runner$run(list(omicsManagedProduct("trx_ruv_grid", c(list(logcpm = held$logcpm), grid_parameters), function(result) {
        held$grid <- result$reference$entity_id
        held$k <- result$value$chosen$k
        normData$ruv_optimization_results <- list(counts = trxNormRuvQcResult(result$value))
        addLog(sprintf("RUV grid (%s): k = %d, %s%% inferred negative controls", result$value$mode, result$value$chosen$k, result$value$chosen$percentage))
    })), done = after_grid, failed = failed)
    runner$run(list(figures("raw"), omicsManagedProduct("trx_tmm_logcpm", list(), function(result) {
        held$logcpm <- result$reference$entity_id
        addLog("TMM factors and logCPM computed")
    })), done = after_logcpm, failed = failed)
    invisible(TRUE)
}

#' Run the tab's sequence on the worker or in process.
#' @noRd
runTrxNormalisation <- function(input, workflowData, experimentPaths, normData, addLog) {
    manager <- workflowData$state_manager
    grid_parameters <- trxNormGridParameters(input)
    normData$post_filter_obj <- getWorkflowStateReference(manager)
    addLog(sprintf("Running TMM, the RUV grid (%s) and RUVs", grid_parameters$ruv_mode))
    failed <- function(error) {
        message <- paste("Normalization error:", trxConditionText(error))
        addLog(message)
        logger::log_error(message)
        shiny::showNotification(message, type = "error", duration = 15)
    }
    finish <- function() {
        normData$normalization_complete <- TRUE
        normData$ruv_complete <- TRUE
        normData$plot_refresh_trigger <- (normData$plot_refresh_trigger %||% 0) + 1
        status <- workflowData$tab_status
        status$normalization <- "complete"
        workflowData$tab_status <- status
        addLog("Normalization complete - ready for differential expression")
        shiny::showNotification("Normalization complete! Proceeding to DE analysis.", type = "message")
    }
    runner <- workflowOmicsManagedActions("transcriptomics")
    if (!is.null(runner)) return(runTrxNormManaged(runner, grid_parameters, workflowData$config_list, normData, addLog, finish, failed))
    tryCatch({
        run <- runTrxReferenceRuvs(manager, grid_parameters, workflowData$config_list)
        normData$ruv_optimization_results <- list(counts = trxNormRuvQcResult(run$grid$value))
        normData$ruv_corrected_obj <- run$reference
        addLog(sprintf("RUVs applied: W_1..W_%d added to the design (state 'trx_ruvs')", run$grid$value$chosen$k))
        runTrxNormPublicationsInProcess(manager, experimentPaths, addLog)
        finish()
    }, error = failed)
}

#' Reset to the state the tab started from (the store step `trx_reset`, recommitting that state's record).
#' @noRd
resetTrxNormalisation <- function(workflowData, normData, addLog) {
    tryCatch({
        if (!is.null(normData$post_filter_obj)) {
            commitWorkflowOmicsStoreStep(workflowData$state_manager, "trx_reset", "trx_reset", workflowData$config_list,
                "Reset to pre-normalization state", source = normData$post_filter_obj)
        }
        normData$normalization_complete <- FALSE
        normData$ruv_complete <- FALSE
        normData$ruv_corrected_obj <- NULL
        normData$ruv_optimization_results <- list()
        status <- workflowData$tab_status
        status$normalization <- "pending"
        workflowData$tab_status <- status
        addLog("Reset to pre-normalization state")
        shiny::showNotification("Reset to pre-normalization state", type = "message")
    }, error = function(error) shiny::showNotification(paste("Error:", trxConditionText(error)), type = "error"))
}

#' @title Transcriptomics Normalisation Tab
#' @description The lipid norm module's server under the transcriptomics profile (`mod_lipid_norm_server(profile =)`).
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label,selected_tab As `mod_lipid_norm_server`.
#' @return The tab's reactive values.
#' @name mod_trx_norm
#' @noRd
mod_trx_norm_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        normData <- createLipidNormReactiveState()
        # One count table; its figures are the worker's trx_<stage>_<kind>.png in the gene QC directory.
        normData$assay_names <- "counts"
        addLog <- buildLipidNormAddLog(normData)
        registerLipidNormStaticQcImageOutputs(output, buildLipidNormQcImageRenderer(normData, experiment_paths,
            sanitizeAssayNameFn = function(assay_name) "trx", qcDirKey = "gene_qc_dir"),
            stagePrefixes = c(pre_norm = "raw", post_norm = "tmm", ruv_corrected = "ruvs"))
        registerLipidNormAssayLabelOutputs(output, buildLipidNormAssayLabelRenderer(normData))
        registerLipidNormLogOutput(output, buildLipidNormLogRenderer(normData))
        registerLipidNormRuvQcOutput(output, buildLipidNormRuvQcUiRenderer(normData, session$ns))
        registerLipidNormRuvCancorOutputs(output, normData)
        shiny::observeEvent(input$run_normalization, {
            shiny::req(workflow_data$state_manager)
            runTrxNormalisation(input, workflow_data, experiment_paths, normData, addLog)
        })
        shiny::observeEvent(input$reset_normalization, resetTrxNormalisation(workflow_data, normData, addLog))
        normData
    })
}
