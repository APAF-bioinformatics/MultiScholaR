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

# ============================================================================
# mod_lipid_qc_s4.R
# ============================================================================
# Purpose: Final S4 state confirmation module for lipidomics QC
#
# KIT-104 W5: the exported wrappers over the shared finalise sub-tab (mod_omic_qc_tabs.R) under the
# lipid configuration, or the profile's omic (trx) under its own. Until KIT-104 W6 moves trx onto
# the shared finalise, a profile that brings its own finalise collaborators (`profile$finalize`)
# keeps today's path: `runLipidQcS4FinalizeWorkflow` with them, and `reportLipidQcS4FinalizeSuccess`,
# which they call.
# ============================================================================

#' @title Lipidomics QC Finalization Module
#' @description A Shiny module for finalizing the QC step and confirming the S4 state.
#'              Displays QC summary statistics and state history.
#' @name mod_lipid_qc_s4
#' @param id,workflow_data,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
#' @param profile An omic profile (KIT-096): its summary and the finalise's collaborators. `NULL`, the default, is lipidomics.
NULL

#' @rdname mod_lipid_qc_s4
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 h5 p hr div actionButton verbatimTextOutput uiOutput icon tags plotOutput
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
mod_lipid_qc_s4_ui <- function(id, profile = NULL) {
    mod_omic_qc_finalize_ui(id, omicGuiConfig(profile$omic_type %||% "lipidomics"))
}

#' The finalise's success under a profile's collaborators (trx's `reportFinalizeSuccessFn` calls it with its text),
#' until KIT-104 W6: the results text, the log and the notification.
#' @noRd
reportLipidQcS4FinalizeSuccess <- function(
    currentS4,
    history,
    output,
    finalizeResultsOutputName = "finalize_results",
    buildResultsTextFn,
    renderTextFn = shiny::renderText,
    logInfoFn = logger::log_info,
    showNotificationFn = shiny::showNotification
) {
    resultText <- buildResultsTextFn(
        currentS4 = currentS4
        , history = history
    )
    # KIT-075: the output's closure keeps this frame, so the object leaves it first.
    rm(currentS4)

    output[[finalizeResultsOutputName]] <- renderTextFn(resultText)

    logInfoFn("Lipidomics QC finalized successfully")
    showNotificationFn(
        "QC complete! Proceed to Normalization."
        , type = "message"
        , duration = 5
    )

    invisible(resultText)
}

#' The finalise under a profile's collaborators (trx's `trxQcFinalizeCollaborators`), until KIT-104 W6: read,
#' check, commit, the figure (the worker's metrics on a managed session), the tab completed and the report.
#' @noRd
runLipidQcS4FinalizeWorkflow <- function(
    workflowData,
    omicType,
    filterPlot,
    output,
    getFinalizeStateFn,
    validateFinalizeStateFn,
    saveCompletedStateFn,
    updateTrackingPlotFn,
    completeTabStatusFn = function(workflowData) completeOmicQcTabStatus(workflowData, TRUE),
    getFinalizeHistoryFn = function(stateManager) stateManager$getHistory(),
    reportFinalizeSuccessFn,
    reportFinalizeErrorFn = reportOmicQcFinalizeError,
    profile = NULL
) {
    tryCatch({
        currentS4 <- getFinalizeStateFn(
            stateManager = workflowData$state_manager
        )
        currentS4 <- validateFinalizeStateFn(currentS4 = currentS4)

        saveCompletedStateFn(
            stateManager = workflowData$state_manager
            , currentS4 = currentS4
            , configObject = workflowData$config_list
        )

        finish <- function() {
            completeTabStatusFn(workflowData = workflowData)

            history <- getFinalizeHistoryFn(
                stateManager = workflowData$state_manager
            )

            reportFinalizeSuccessFn(
                currentS4 = currentS4
                , history = history
                , output = output
            )
        }

        # KIT-082: the tracking figure's metrics read the whole matrix, so a managed session's worker computes them.
        # KIT-096: a profile's runner, metrics operation and state reader.
        runner <- workflowOmicsManagedActions(profile$omic_type)
        if (!is.null(runner)) {
            runner$run(list(omicsManagedProduct(paste0(profile$prefix, "_filtering_metrics"), list(), function(result) {
                updateTrackingPlotFn(
                    currentS4 = NULL
                    , omicType = omicType
                    , setFilterPlotFn = filterPlot
                    , metrics = readWorkflowOmicsFilteringMetrics(workflowData$state_manager, result)
                )
            })),
                # In process a tracking figure's failure leaves the tab to complete, with no figure; the state is committed either way.
                done = function() tryCatch(finish(), error = function(error) reportFinalizeErrorFn(error = error)),
                failed = function(error) {
                    filterPlot(NULL)
                    tryCatch(finish(), error = function(error) reportFinalizeErrorFn(error = error))
                })
            return(invisible(NULL))
        }

        updateTrackingPlotFn(
            currentS4 = profile$read_state(workflowData$state_manager)
            , omicType = omicType
            , setFilterPlotFn = filterPlot
        )

        finish()
    }, error = function(e) {
        reportFinalizeErrorFn(error = e)
    })
}

#' @rdname mod_lipid_qc_s4
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification renderText renderUI observe tags renderPlot
#' @importFrom DT renderDT datatable
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_lipid_qc_s4_server <- function(id, workflow_data, omic_type, experiment_label, profile = NULL) {
    config <- omicGuiConfig(profile$omic_type %||% "lipidomics")
    if (is.null(profile$finalize)) {
        return(mod_omic_qc_finalize_server(id, workflow_data, omic_type, experiment_label, config))
    }
    shiny::moduleServer(id, function(input, output, session) {
        registerLipidQcS4ProfileFinalize(input, output, workflow_data, omic_type, profile)
    })
}

#' The finalise sub-tab under a profile with its own summary and collaborators (trx), until KIT-104 W6: the shared
#' history and figure outputs, the profile's summary, and Finalize QC through `runLipidQcS4FinalizeWorkflow`.
#' @noRd
registerLipidQcS4ProfileFinalize <- function(input, output, workflow_data, omic_type, profile) {
    plot <- shiny::reactiveVal(NULL)
    registerOmicQcStateHistory(output, workflow_data)
    # KIT-096: a profile's summary, read from the state's metadata.
    profile$register_qc_summary(output, workflow_data)
    finalize_after_start <- qcClickAfterStart(input, "finalize_qc") # KIT-101 F-24
    shiny::observeEvent(input$finalize_qc, {
        if (!finalize_after_start()) return(invisible(NULL))
        shiny::req(workflow_data$state_manager)
        do.call(runLipidQcS4FinalizeWorkflow, c(list(workflowData = workflow_data, omicType = omic_type,
            filterPlot = plot, output = output, profile = profile), profile$finalize))
    })
    registerOmicQcPlotOutput(output, plot)
    invisible(NULL)
}
