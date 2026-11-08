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

# KIT-110 W2: the transcriptomics normalisation core, moved from R/mod_trx_norm.R (KIT-097) at d5244ffe. The tab's
# sequence: TMM and logCPM (`trx_tmm_logcpm`), the RUV grid (`trx_ruv_grid`, automatic as the notebook or manual), RUVs
# (`trx_ruvs`, the DA head), then the PCA and RLE at raw, TMM and RUVs, the CCA panels and the stage grid. The run and
# the reset are the omic's own handlers: they take the tab's settings and its `notify` and `setStatus` callbacks, and
# call no Shiny UI function (KIT-110 D4, D9). The managed sequence runs the same operations in the same order on the
# worker (KIT-075, 078, 082; D6). `runWorkflowOmicsManagedNormalisation` (ITSD, log2, between-sample, RUV-III) is not
# this sequence and is not reused.

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

#' The figures and tables after RUVs (the stages' PCA and RLE, the CCA panels, the stage grid), in process.
#' @noRd
runTrxNormPublicationsInProcess <- function(manager, experiment_paths, addLog) {
    for (stage in c("raw", "tmm", "ruvs")) {
        runWorkflowOmicsPublicationInProcess(transcriptomicsDeclarationBinding, manager, "trx_qc_figures", list(stage = stage),
            experiment_paths)
        addLog(sprintf("PCA and RLE (%s) generated", stage))
    }
    runWorkflowOmicsPublicationInProcess(transcriptomicsDeclarationBinding, manager, "trx_ruv_publication", list(), experiment_paths)
    runWorkflowOmicsPublicationInProcess(transcriptomicsDeclarationBinding, manager, "trx_composite_figure",
        list(figure = "normalization_comparison"), experiment_paths)
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
#' action names the products the earlier ones delivered, so each runs once those have. The arguments are the frame's
#' managed call (`runOmicNormPipeline`); `binding` and `manager` are not read.
#' @param settings The grid's parameters (`trxNormGridParameters`).
#' @noRd
runTrxNormManaged <- function(binding, runner, manager, settings, normData, config, addLog, finish, failed) {
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
    after_logcpm <- function() runner$run(list(omicsManagedProduct("trx_ruv_grid", c(list(logcpm = held$logcpm), settings), function(result) {
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

#' The sequence in process: RUVs on the head through the wrappers (`runTrxReferenceRuvs`), then the figures.
#' @param settings The grid's parameters (`trxNormGridParameters`).
#' @noRd
runTrxNormInProcess <- function(binding, manager, settings, normData, config, addLog, experimentPaths, finish, failed) {
    tryCatch({
        run <- runTrxReferenceRuvs(manager, settings, config)
        normData$ruv_optimization_results <- list(counts = trxNormRuvQcResult(run$grid$value))
        normData$ruv_corrected_obj <- run$reference
        addLog(sprintf("RUVs applied: W_1..W_%d added to the design (state 'trx_ruvs')", run$grid$value$chosen$k))
        runTrxNormPublicationsInProcess(manager, experimentPaths, addLog)
        finish()
    }, error = failed)
}

#' Run the tab's sequence on the worker (`run$managed`, with the omic's runner) or in process (`run$in_process`).
#' @param run `managed` and `in_process`, the sequence's two forms.
#' @param settings The grid's parameters (`trxNormGridParameters`).
#' @param notify `function(...)` with `shiny::showNotification`'s arguments.
#' @param setStatus `function(tab, value)` writing the tab's status.
#' @noRd
runTrxNormalisation <- function(run, binding, manager, settings, normData, config, experimentPaths, addLog, notify, setStatus) {
    force(manager)
    force(settings)
    normData$post_filter_obj <- getWorkflowStateReference(manager)
    addLog(sprintf("Running TMM, the RUV grid (%s) and RUVs", settings$ruv_mode))
    failed <- function(error) {
        message <- paste("Normalization error:", trxConditionText(error))
        addLog(message)
        logger::log_error(message)
        notify(message, type = "error", duration = 15)
    }
    finish <- function() {
        normData$normalization_complete <- TRUE
        normData$ruv_complete <- TRUE
        normData$plot_refresh_trigger <- (normData$plot_refresh_trigger %||% 0) + 1
        setStatus("normalization", "complete")
        addLog("Normalization complete - ready for differential expression")
        notify("Normalization complete! Proceeding to DE analysis.", type = "message")
    }
    runner <- workflowOmicsManagedActions(binding$omic_type)
    if (!is.null(runner)) return(run$managed(binding, runner, manager, settings, normData, config, addLog, finish, failed))
    run$in_process(binding, manager, settings, normData, config, addLog, experimentPaths, finish, failed)
}

#' Reset to the state the tab started from (the store step `trx_reset`, recommitting that state's record).
#' @param notify,setStatus As `runTrxNormalisation`'s.
#' @noRd
resetTrxNormalisation <- function(manager, normData, config, addLog, notify, setStatus) {
    tryCatch({
        if (!is.null(normData$post_filter_obj)) {
            commitWorkflowOmicsStoreStep(manager, "trx_reset", "trx_reset", config, "Reset to pre-normalization state",
                source = normData$post_filter_obj)
        }
        normData$normalization_complete <- FALSE
        normData$ruv_complete <- FALSE
        normData$ruv_corrected_obj <- NULL
        normData$ruv_optimization_results <- list()
        setStatus("normalization", "pending")
        addLog("Reset to pre-normalization state")
        notify("Reset to pre-normalization state", type = "message")
    }, error = function(error) notify(paste("Error:", trxConditionText(error)), type = "error"))
}
