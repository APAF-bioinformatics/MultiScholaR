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
# This module provides a summary of QC processing and allows the user to
# finalize the QC step before proceeding to normalization.
# ============================================================================

#' @title Lipidomics QC Finalization Module
#' @description A Shiny module for finalizing the QC step and confirming the S4 state.
#'              Displays QC summary statistics and state history.
#' @name mod_lipid_qc_s4
NULL

#' @rdname mod_lipid_qc_s4
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 h5 p hr div actionButton verbatimTextOutput uiOutput icon tags plotOutput
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
mod_lipid_qc_s4_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tabPanel(
        "Finalize QC"
        , shiny::br()
        , shiny::fluidRow(
            shiny::column(6
                , shiny::wellPanel(
                    shiny::h4("QC Processing Summary")
                    , shiny::p("Review the QC processing steps before finalizing. You can revert to any previous state if needed.")
                    , shiny::hr()
                    
                    # State history
                    , shiny::h5("Processing History")
                    , shiny::uiOutput(ns("state_history"))
                    
                    , shiny::hr()
                    , shiny::div(
                        shiny::actionButton(
                            ns("finalize_qc")
                            , "Finalize QC"
                            , class = "btn-success"
                            , width = "100%"
                            , icon = shiny::icon("check")
                        )
                    )
                )
            )
            , shiny::column(6
                , shiny::wellPanel(
                    shiny::h4("Current Data Summary")
                    , shiny::hr()
                    , shiny::uiOutput(ns("data_summary"))
                )
                , shiny::wellPanel(
                    shiny::h4("Per-Assay Statistics")
                    , shiny::hr()
                    , DT::DTOutput(ns("assay_stats_table"))
                )
            )
        )
        , shiny::fluidRow(
            shiny::column(12
                , shiny::verbatimTextOutput(ns("finalize_results"))
                , shiny::br()
                # QC Progress Grid
                , shinyjqui::jqui_resizable(
                    shiny::plotOutput(ns("filter_plot"), height = "600px", width = "100%")
                )
            )
        )
    )
}

buildLipidQcS4StateHistoryUi <- function(history) {
    if (length(history) == 0) {
        return(shiny::p(
            shiny::icon("info-circle")
            , " No processing history available."
            , style = "color: #666;"
        ))
    }

    history_items <- lapply(seq_along(history), function(i) {
        state_name <- history[i]
        is_current <- i == length(history)

        icon_name <- if (is_current) "arrow-right" else "check"
        icon_color <- if (is_current) "blue" else "green"
        text_style <- if (is_current) "font-weight: bold;" else ""

        shiny::tags$li(
            shiny::icon(icon_name, style = paste("color:", icon_color))
            , shiny::tags$span(
                sprintf(" %d. %s", i, state_name)
                , style = text_style
            )
            , if (is_current) shiny::tags$span(" (current)", style = "color: blue;") else NULL
        )
    })

    shiny::tags$ol(history_items, style = "list-style: none; padding-left: 0;")
}

registerLipidQcS4StateHistoryOutput <- function(
    output,
    workflowData,
    renderUiFn = shiny::renderUI,
    reqFn = shiny::req,
    buildStateHistoryUiFn = buildLipidQcS4StateHistoryUi
) {
    output$state_history <- renderUiFn({
        reqFn(workflowData$state_manager)

        history <- tryCatch({
            workflowData$state_manager$getHistory()
        }, error = function(e) {
            character(0)
        })

        buildStateHistoryUiFn(history)
    })

    output
}

buildLipidQcS4DataSummaryUi <- function(currentS4) {
    if (is.null(currentS4) || !inherits(currentS4, "LipidomicsAssayData")) {
        return(shiny::p(
            shiny::icon("exclamation-triangle", style = "color: orange;")
            , " No LipidomicsAssayData object available."
        ))
    }

    assay_list <- currentS4@lipid_data
    n_assays <- length(assay_list)

    total_lipids <- sum(vapply(assay_list, function(a) {
        if (currentS4@lipid_id_column %in% names(a)) {
            length(unique(a[[currentS4@lipid_id_column]]))
        } else {
            nrow(a)
        }
    }, numeric(1)))

    sample_cols <- if (n_assays > 0) {
        names(assay_list[[1]])[vapply(assay_list[[1]], is.numeric, logical(1))]
    } else {
        character(0)
    }
    n_samples <- length(sample_cols)

    dm <- currentS4@design_matrix
    n_groups <- if (nrow(dm) > 0 && currentS4@group_id %in% names(dm)) {
        length(unique(dm[[currentS4@group_id]]))
    } else {
        NA
    }

    shiny::tagList(
        shiny::tags$table(
            class = "table table-condensed"
            , style = "margin-bottom: 0;"
            , shiny::tags$tbody(
                shiny::tags$tr(
                    shiny::tags$td(shiny::strong("Number of Assays:"))
                    , shiny::tags$td(n_assays)
                )
                , shiny::tags$tr(
                    shiny::tags$td(shiny::strong("Total Lipids:"))
                    , shiny::tags$td(format(total_lipids, big.mark = ","))
                )
                , shiny::tags$tr(
                    shiny::tags$td(shiny::strong("Number of Samples:"))
                    , shiny::tags$td(n_samples)
                )
                , shiny::tags$tr(
                    shiny::tags$td(shiny::strong("Experimental Groups:"))
                    , shiny::tags$td(if (!is.na(n_groups)) n_groups else "N/A")
                )
                , shiny::tags$tr(
                    shiny::tags$td(shiny::strong("Lipid ID Column:"))
                    , shiny::tags$td(currentS4@lipid_id_column)
                )
                , shiny::tags$tr(
                    shiny::tags$td(shiny::strong("Sample ID Column:"))
                    , shiny::tags$td(currentS4@sample_id)
                )
            )
        )
    )
}

registerLipidQcS4DataSummaryOutput <- function(
    output,
    workflowData,
    renderUiFn = shiny::renderUI,
    reqFn = shiny::req,
    buildDataSummaryUiFn = buildLipidQcS4DataSummaryUi
) {
    output$data_summary <- renderUiFn({
        reqFn(workflowData$state_manager)

        current_s4 <- tryCatch({
            workflowData$state_manager$getState()
        }, error = function(e) {
            NULL
        })

        buildDataSummaryUiFn(current_s4)
    })

    output
}

buildLipidQcS4AssayStatsTable <- function(
    currentS4,
    datatableFn = DT::datatable
) {
    if (is.null(currentS4) || !inherits(currentS4, "LipidomicsAssayData")) {
        return(NULL)
    }

    assay_list <- currentS4@lipid_data
    lipid_id_col <- currentS4@lipid_id_column

    stats_df <- data.frame(
        Assay = names(assay_list)
        , Lipids = sapply(assay_list, function(a) {
            if (lipid_id_col %in% names(a)) {
                length(unique(a[[lipid_id_col]]))
            } else {
                nrow(a)
            }
        })
        , Samples = sapply(assay_list, function(a) {
            sum(sapply(a, is.numeric))
        })
        , stringsAsFactors = FALSE
    )

    stats_df$Missingness <- sapply(assay_list, function(a) {
        quant_cols <- names(a)[sapply(a, is.numeric)]
        if (length(quant_cols) == 0) return(NA)

        total_cells <- nrow(a) * length(quant_cols)
        missing_cells <- sum(sapply(quant_cols, function(col) {
            sum(is.na(a[[col]]) | a[[col]] == 0)
        }))

        round((missing_cells / total_cells) * 100, 1)
    })

    datatableFn(
        stats_df
        , options = list(
            dom = "t"
            , paging = FALSE
            , ordering = FALSE
        )
        , rownames = FALSE
        , class = "compact stripe"
    )
}

registerLipidQcS4AssayStatsOutput <- function(
    output,
    workflowData,
    renderDtFn = DT::renderDT,
    reqFn = shiny::req,
    buildAssayStatsTableFn = buildLipidQcS4AssayStatsTable
) {
    output$assay_stats_table <- renderDtFn({
        reqFn(workflowData$state_manager)

        current_s4 <- tryCatch({
            workflowData$state_manager$getState()
        }, error = function(e) {
            NULL
        })

        buildAssayStatsTableFn(current_s4)
    })

    output
}

buildLipidQcS4FinalizeResultsText <- function(currentS4, history) {
    assay_list <- currentS4@lipid_data
    total_lipids <- sum(vapply(assay_list, function(a) {
        if (currentS4@lipid_id_column %in% names(a)) {
            length(unique(a[[currentS4@lipid_id_column]]))
        } else {
            nrow(a)
        }
    }, numeric(1)))

    paste(
        "QC Finalization Complete"
        , "========================"
        , ""
        , sprintf("Total lipids retained: %d", total_lipids)
        , sprintf("Number of assays: %d", length(assay_list))
        , sprintf("Processing steps completed: %d", length(history))
        , ""
        , "Processing History:"
        , paste(sprintf("  %d. %s", seq_along(history), history), collapse = "\n")
        , ""
        , "State saved as: 'lipid_qc_complete'"
        , ""
        , "You can now proceed to the Normalization tab."
        , sep = "\n"
    )
}

getLipidQcS4FinalizeState <- function(
    stateManager,
    stateGetterFn = function(manager) manager$getState()
) {
    stateGetterFn(stateManager)
}

validateLipidQcS4FinalizeState <- function(
    currentS4,
    reqFn = shiny::req,
    inheritsFn = inherits,
    expectedClass = "LipidomicsAssayData",
    errorMessage = paste("Current state is not a", expectedClass, "object"),
    requireAssays = TRUE,
    requireDesign = TRUE
) {
    reqFn(currentS4)

    if (!inheritsFn(currentS4, expectedClass)) {
        stop(errorMessage)
    }
    currentSlots <- methods::slotNames(currentS4)
    if (isTRUE(requireAssays) &&
        (!"lipid_data" %in% currentSlots ||
            length(currentS4@lipid_data) == 0L ||
            any(vapply(currentS4@lipid_data, nrow, integer(1)) == 0L))) {
        stop("Current LipidomicsAssayData object has no assay data to finalize.")
    }
    if (isTRUE(requireDesign) &&
        (!"design_matrix" %in% currentSlots ||
            !is.data.frame(currentS4@design_matrix) ||
            nrow(currentS4@design_matrix) == 0L)) {
        stop("Current LipidomicsAssayData object has no design matrix rows to finalize.")
    }

    if (isTRUE(requireDesign) && isTRUE(requireAssays)) {
        sample_id_col <- currentS4@sample_id
        group_id_col <- currentS4@group_id
        sample_columns <- if (sample_id_col %in% names(currentS4@design_matrix)) {
            as.character(currentS4@design_matrix[[sample_id_col]])
        } else {
            character()
        }
        alignment <- validateLipidDesignAlignment(
            designMatrix = currentS4@design_matrix,
            assayList = currentS4@lipid_data,
            columnMapping = list(
                lipid_id_col = currentS4@lipid_id_column,
                annotation_col = currentS4@annotation_id_column,
                sample_columns = sample_columns
            ),
            sampleIdCol = sample_id_col,
            groupCol = group_id_col,
            requireReplicates = FALSE
        )
        if (!isTRUE(alignment$valid)) {
            stop(sprintf(
                "Current LipidomicsAssayData object has assay/design mismatch: %s",
                paste(alignment$errors, collapse = "; ")
            ))
        }
    }

    currentS4
}

saveLipidQcS4CompletedState <- function(
    stateManager,
    currentS4,
    configObject,
    stateName = "lipid_qc_complete",
    description = "QC processing complete - ready for normalization"
) {
    stateManager$saveState(
        state_name = stateName
        , s4_data_object = currentS4
        , config_object = configObject
        , description = description
    )

    invisible(stateName)
}

completeLipidQcS4TabStatus <- function(
    workflowData,
    tabName = "quality_control",
    status = "complete",
    normalizationTabName = "normalization",
    normalizationReadyStatus = "pending"
) {
    updatedStatus <- workflowData$tab_status
    updatedStatus[[tabName]] <- status
    if (normalizationTabName %in% names(updatedStatus) &&
        updatedStatus[[normalizationTabName]] %in% c("disabled", "locked")) {
        updatedStatus[[normalizationTabName]] <- normalizationReadyStatus
    }
    workflowData$tab_status <- updatedStatus

    invisible(updatedStatus)
}

getLipidQcS4FinalizeHistory <- function(
    stateManager,
    historyGetterFn = function(manager) manager$getHistory()
) {
    historyGetterFn(stateManager)
}

updateLipidQcS4TrackingPlot <- function(
    currentS4,
    omicType,
    setFilterPlotFn,
    stepName = "4_QC_Complete",
    updateLipidFilteringFn = updateLipidFiltering
) {
    qcPlot <- tryCatch({
        updateLipidFilteringFn(
            theObject = currentS4
            , step_name = stepName
            , omics_type = omicType
            , return_grid = TRUE
            , overwrite = TRUE
        )
    }, error = function(e) {
        NULL
    })

    setFilterPlotFn(qcPlot)

    invisible(qcPlot)
}

reportLipidQcS4FinalizeSuccess <- function(
    currentS4,
    history,
    output,
    finalizeResultsOutputName = "finalize_results",
    buildResultsTextFn = buildLipidQcS4FinalizeResultsText,
    renderTextFn = shiny::renderText,
    logInfoFn = logger::log_info,
    showNotificationFn = shiny::showNotification
) {
    resultText <- buildResultsTextFn(
        currentS4 = currentS4
        , history = history
    )

    output[[finalizeResultsOutputName]] <- renderTextFn(resultText)

    logInfoFn("Lipidomics QC finalized successfully")
    showNotificationFn(
        "QC complete! Proceed to Normalization."
        , type = "message"
        , duration = 5
    )

    invisible(resultText)
}

reportLipidQcS4FinalizeError <- function(
    error,
    messagePrefix = "Error finalizing QC:",
    pasteFn = paste,
    logErrorFn = logger::log_error,
    showNotificationFn = shiny::showNotification,
    notificationType = "error"
) {
    messageText <- pasteFn(messagePrefix, error$message)

    logErrorFn(messageText)
    showNotificationFn(messageText, type = notificationType)

    invisible(messageText)
}

runLipidQcS4FinalizeWorkflow <- function(
    workflowData,
    omicType,
    filterPlot,
    output,
    getFinalizeStateFn = getLipidQcS4FinalizeState,
    validateFinalizeStateFn = validateLipidQcS4FinalizeState,
    saveCompletedStateFn = saveLipidQcS4CompletedState,
    updateTrackingPlotFn = updateLipidQcS4TrackingPlot,
    completeTabStatusFn = completeLipidQcS4TabStatus,
    getFinalizeHistoryFn = getLipidQcS4FinalizeHistory,
    reportFinalizeSuccessFn = reportLipidQcS4FinalizeSuccess,
    reportFinalizeErrorFn = reportLipidQcS4FinalizeError
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

        updateTrackingPlotFn(
            currentS4 = currentS4
            , omicType = omicType
            , setFilterPlotFn = filterPlot
        )

        completeTabStatusFn(workflowData = workflowData)

        history <- getFinalizeHistoryFn(
            stateManager = workflowData$state_manager
        )

        reportFinalizeSuccessFn(
            currentS4 = currentS4
            , history = history
            , output = output
        )
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
mod_lipid_qc_s4_server <- function(id, workflow_data, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        filter_plot <- shiny::reactiveVal(NULL)

        registerLipidQcS4StateHistoryOutput(
            output = output
            , workflowData = workflow_data
        )

        registerLipidQcS4DataSummaryOutput(
            output = output
            , workflowData = workflow_data
        )

        registerLipidQcS4AssayStatsOutput(
            output = output
            , workflowData = workflow_data
        )
        
        # Finalize QC
        shiny::observeEvent(input$finalize_qc, {
            shiny::req(workflow_data$state_manager)

            runLipidQcS4FinalizeWorkflow(
                workflowData = workflow_data
                , omicType = omic_type
                , filterPlot = filter_plot
                , output = output
            )
        })

        # Render QC progress plot
        output$filter_plot <- shiny::renderPlot({

            shiny::req(filter_plot())


            if (inherits(filter_plot(), "grob") || inherits(filter_plot(), "gtable")) {
                grid::grid.draw(filter_plot())
            } else if (inherits(filter_plot(), "ggplot")) {
                print(filter_plot())
            } else {
            }
        })
    })
}
