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
# mod_metab_qc_intensity.R
# ============================================================================
# Purpose: Metabolite intensity filtering Shiny module
#
# This module wraps the metaboliteIntensityFiltering() S4 method to provide
# an interactive UI for filtering metabolites based on intensity thresholds.
# ============================================================================

#' @title Metabolite Intensity Filter Module
#' @description A Shiny module for applying metabolite intensity and missing value filters.
#'              Wraps the metaboliteIntensityFiltering() S4 method.
#' @name mod_metab_qc_intensity
NULL

#' @rdname mod_metab_qc_intensity
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr numericInput sliderInput div actionButton verbatimTextOutput uiOutput
#' @importFrom shinyjqui jqui_resizable
mod_metab_qc_intensity_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tabPanel(
        "Intensity Filter"
        , shiny::br()
        , shiny::fluidRow(
            shiny::column(4
                , shiny::wellPanel(
                    shiny::h4("Metabolite Intensity Filter Parameters")
                    , shiny::p("Filter metabolites based on intensity thresholds. Metabolites with low intensities across many samples will be removed.")
                    , shiny::hr()
                    
                    # Intensity cutoff percentile
                    , shiny::sliderInput(
                        ns("intensity_cutoff_percentile")
                        , "Intensity Cutoff Percentile (%)"
                        , min = 1
                        , max = 50
                        , value = 10
                        , step = 1
                    )
                    , shiny::helpText("Intensities below this percentile are considered 'low' (default: 10%)")
                    
                    # Proportion threshold
                    , shiny::sliderInput(
                        ns("proportion_below_cutoff")
                        , "Max Proportion of Samples Below Threshold"
                        , min = 0.1
                        , max = 1.0
                        , value = 0.5
                        , step = 0.05
                    )
                    , shiny::helpText("Remove metabolites where this proportion of samples are below the intensity threshold (default: 0.5)")
                    
                    , shiny::hr()
                    , shiny::div(
                        shiny::actionButton(
                            ns("apply_filter")
                            , "Apply Filter"
                            , class = "btn-primary"
                            , width = "48%"
                        )
                        , shiny::actionButton(
                            ns("revert_filter")
                            , "Revert"
                            , class = "btn-warning"
                            , width = "48%"
                            , style = "margin-left: 4%"
                        )
                    )
                )
            )
            , shiny::column(8
                , shiny::verbatimTextOutput(ns("filter_results"))
                , shiny::br()
                # Per-assay results tabs
                , shiny::uiOutput(ns("assay_results_tabs"))
                , shiny::br()
                , shinyjqui::jqui_resizable(
                    shiny::plotOutput(ns("filter_plot"), height = "600px", width = "100%")
                )
            )
        )
    )
}

# Keep the per-assay tab assembly isolated so the wrapper can shed render logic
# without changing the public server entry point.
buildMetabIntensityAssayTabsUi <- function(
    stats,
    ns,
    tabPanelFn = shiny::tabPanel,
    brFn = shiny::br,
    tableTagFn = shiny::tags$table,
    tbodyTagFn = shiny::tags$tbody,
    trTagFn = shiny::tags$tr,
    tdTagFn = shiny::tags$td,
    strongFn = shiny::strong,
    tabsetPanelFn = shiny::tabsetPanel
) {
    if (is.null(stats) || nrow(stats) == 0) {
        return(NULL)
    }

    tab_list <- lapply(seq_len(nrow(stats)), function(i) {
        tabPanelFn(
            stats$Assay[i],
            brFn(),
            tableTagFn(
                class = "table table-striped table-condensed",
                style = "width: auto;",
                tbodyTagFn(
                    trTagFn(
                        tdTagFn(strongFn("Original metabolites:")),
                        tdTagFn(stats$Original[i])
                    ),
                    trTagFn(
                        tdTagFn(strongFn("After filtering:")),
                        tdTagFn(stats$Filtered[i])
                    ),
                    trTagFn(
                        tdTagFn(strongFn("Removed:")),
                        tdTagFn(stats$Removed[i])
                    ),
                    trTagFn(
                        tdTagFn(strongFn("Percent retained:")),
                        tdTagFn(paste0(stats$Percent_Retained[i], "%"))
                    )
                )
            )
        )
    })

    do.call(tabsetPanelFn, c(list(id = ns("assay_stats_tabs")), tab_list))
}

renderMetabIntensityFilterPlot <- function(
    filterPlot,
    reqFn = shiny::req,
    inheritsFn = inherits,
    gridDrawFn = grid::grid.draw,
    printFn = print
) {
    plotObject <- filterPlot()
    reqFn(plotObject)

    if (inheritsFn(plotObject, "grob") || inheritsFn(plotObject, "gtable")) {
        gridDrawFn(plotObject)
    } else if (inheritsFn(plotObject, "ggplot")) {
        printFn(plotObject)
    }

    invisible(NULL)
}

buildMetabIntensityFilterStats <- function(currentS4, filteredS4) {
    countMetabolites <- function(s4Object) {
        vapply(s4Object@metabolite_data, function(assay) {
            if (s4Object@metabolite_id_column %in% names(assay)) {
                length(unique(assay[[s4Object@metabolite_id_column]]))
            } else {
                nrow(assay)
            }
        }, numeric(1))
    }

    originalCounts <- countMetabolites(currentS4)
    filteredCounts <- countMetabolites(filteredS4)

    statsDf <- data.frame(
        Assay = names(originalCounts)
        , Original = as.numeric(originalCounts)
        , Filtered = as.numeric(filteredCounts)
        , Removed = as.numeric(originalCounts) - as.numeric(filteredCounts)
        , stringsAsFactors = FALSE
    )
    statsDf$Percent_Retained <- round(
        (statsDf$Filtered / statsDf$Original) * 100
        , 1
    )

    statsDf
}

updateMetabIntensityFilterQcPlot <- function(
    filteredS4,
    omicType,
    setFilterPlotFn,
    stepName = "2_Intensity_Filtered",
    updateMetaboliteFilteringFn = updateMetaboliteFiltering,
    logWarnFn = logger::log_warn
) {
    qcPlot <- tryCatch({
        updateMetaboliteFilteringFn(
            theObject = filteredS4
            , step_name = stepName
            , omics_type = omicType
            , return_grid = TRUE
            , overwrite = TRUE
        )
    }, error = function(e) {
        logWarnFn(paste("Could not generate QC plot:", e$message))
        NULL
    })

    setFilterPlotFn(qcPlot)

    invisible(qcPlot)
}

saveMetabIntensityFilterState <- function(
    stateManager,
    filteredS4,
    configObject,
    intensityCutoffPercentile,
    proportionBelowCutoff,
    stateName = "metab_intensity_filtered",
    sprintfFn = sprintf
) {
    stateManager$saveState(
        state_name = stateName
        , s4_data_object = filteredS4
        , config_object = configObject
        , description = sprintfFn(
            "Applied metabolite intensity filter (percentile: %d%%, proportion: %.2f)"
            , intensityCutoffPercentile
            , proportionBelowCutoff
        )
    )

    invisible(stateName)
}

buildMetabIntensityFilterSummary <- function(
    stats,
    intensityCutoffPercentile,
    proportionBelowCutoff,
    stateName = "metab_intensity_filtered"
) {
    totalOriginal <- sum(stats$Original)
    totalFiltered <- sum(stats$Filtered)
    totalRemoved <- sum(stats$Removed)

    perAssayText <- paste(vapply(seq_len(nrow(stats)), function(i) {
        sprintf(
            "  %s: %d -> %d (removed %d, %.1f%% retained)"
            , stats$Assay[i]
            , stats$Original[i]
            , stats$Filtered[i]
            , stats$Removed[i]
            , stats$Percent_Retained[i]
        )
    }, character(1)), collapse = "\n")

    list(
        totalOriginal = totalOriginal
        , totalFiltered = totalFiltered
        , totalRemoved = totalRemoved
        , resultText = paste(
            "Metabolite Intensity Filter Applied Successfully"
            , "============================================"
            , sprintf("Intensity cutoff percentile: %d%%", intensityCutoffPercentile)
            , sprintf("Max proportion below cutoff: %.2f", proportionBelowCutoff)
            , ""
            , "Per-Assay Results:"
            , perAssayText
            , ""
            , sprintf(
                "Total: %d -> %d metabolites (removed %d)"
                , totalOriginal
                , totalFiltered
                , totalRemoved
            )
            , ""
            , sprintf("State saved as: '%s'", stateName)
            , sep = "\n"
        )
    )
}

reportMetabIntensityFilterSuccess <- function(
    stats,
    intensityCutoffPercentile,
    proportionBelowCutoff,
    savedStateName,
    output,
    filterResultsOutputName = "filter_results",
    workingNotificationId = "metab_intensity_filter_working",
    buildSummaryFn = buildMetabIntensityFilterSummary,
    renderTextFn = shiny::renderText,
    logInfoFn = logger::log_info,
    removeNotificationFn = shiny::removeNotification,
    showNotificationFn = shiny::showNotification,
    sprintfFn = sprintf
) {
    filterSummary <- buildSummaryFn(
        stats = stats
        , intensityCutoffPercentile = intensityCutoffPercentile
        , proportionBelowCutoff = proportionBelowCutoff
        , stateName = savedStateName
    )

    output[[filterResultsOutputName]] <- renderTextFn(filterSummary$resultText)

    logInfoFn(paste(
        "Metabolite intensity filter applied: removed"
        , filterSummary$totalRemoved
        , "metabolites"
    ))

    removeNotificationFn(workingNotificationId)
    showNotificationFn(
        sprintfFn(
            "Intensity filter applied: %d metabolites retained"
            , filterSummary$totalFiltered
        )
        , type = "message"
    )

    invisible(filterSummary)
}

reportMetabIntensityFilterRevertSuccess <- function(
    prevStateName,
    output,
    setFilterStatsFn,
    setFilterPlotFn,
    filterResultsOutputName = "filter_results",
    renderTextFn = shiny::renderText,
    logInfoFn = logger::log_info,
    showNotificationFn = shiny::showNotification,
    sprintfFn = sprintf
) {
    revertMessage <- sprintfFn(
        "Reverted to previous state: %s"
        , prevStateName
    )

    output[[filterResultsOutputName]] <- renderTextFn(revertMessage)
    setFilterStatsFn(NULL)
    setFilterPlotFn(NULL)

    logInfoFn(paste("Reverted metabolite intensity filter to", prevStateName))
    showNotificationFn("Reverted successfully", type = "message")

    invisible(revertMessage)
}

#' @rdname mod_metab_qc_intensity
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot renderUI tabsetPanel tabPanel
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_metab_qc_intensity_server <- function(id, workflow_data, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        filter_plot <- shiny::reactiveVal(NULL)
        filter_stats <- shiny::reactiveVal(NULL)
        
        # Apply intensity filter
        shiny::observeEvent(input$apply_filter, {
            shiny::req(workflow_data$state_manager)
            
            shiny::showNotification(
                "Applying metabolite intensity filter..."
                , id = "metab_intensity_filter_working"
                , duration = NULL
            )
            
            tryCatch({
                # Get current MetaboliteAssayData S4 object from state
                current_s4 <- workflow_data$state_manager$getState()
                shiny::req(current_s4)
                
                # Verify it's a MetaboliteAssayData object
                if (!inherits(current_s4, "MetaboliteAssayData")) {
                    stop("Current state is not a MetaboliteAssayData object")
                }
                
                logger::log_info(paste(
                    "Applying metabolite intensity filter: percentile ="
                    , input$intensity_cutoff_percentile
                    , ", proportion ="
                    , input$proportion_below_cutoff
                ))
                
                # Apply S4 method for filtering
                filtered_s4 <- metaboliteIntensityFiltering(
                    current_s4
                    , metabolites_intensity_cutoff_percentile = input$intensity_cutoff_percentile
                    , metabolites_proportion_of_samples_below_cutoff = input$proportion_below_cutoff
                )
                
                stats_df <- buildMetabIntensityFilterStats(
                    currentS4 = current_s4
                    , filteredS4 = filtered_s4
                )
                filter_stats(stats_df)
                
                updateMetabIntensityFilterQcPlot(
                    filteredS4 = filtered_s4
                    , omicType = omic_type
                    , setFilterPlotFn = filter_plot
                )

                saved_state_name <- saveMetabIntensityFilterState(
                    stateManager = workflow_data$state_manager
                    , filteredS4 = filtered_s4
                    , configObject = workflow_data$config_list
                    , intensityCutoffPercentile = input$intensity_cutoff_percentile
                    , proportionBelowCutoff = input$proportion_below_cutoff
                )

                reportMetabIntensityFilterSuccess(
                    stats = stats_df
                    , intensityCutoffPercentile = input$intensity_cutoff_percentile
                    , proportionBelowCutoff = input$proportion_below_cutoff
                    , savedStateName = saved_state_name
                    , output = output
                )
                
            }, error = function(e) {
                msg <- paste("Error applying metabolite intensity filter:", e$message)
                logger::log_error(msg)
                shiny::showNotification(msg, type = "error", duration = 15)
                shiny::removeNotification("metab_intensity_filter_working")
            })
        })
        
        # Revert filter
        shiny::observeEvent(input$revert_filter, {
            tryCatch({
                history <- workflow_data$state_manager$getHistory()
                if (length(history) > 1) {
                    prev_state_name <- history[length(history) - 1]
                    workflow_data$state_manager$revertToState(prev_state_name)
                    reportMetabIntensityFilterRevertSuccess(
                        prevStateName = prev_state_name
                        , output = output
                        , setFilterStatsFn = filter_stats
                        , setFilterPlotFn = filter_plot
                    )
                } else {
                    stop("No previous state to revert to.")
                }
            }, error = function(e) {
                msg <- paste("Error reverting:", e$message)
                logger::log_error(msg)
                shiny::showNotification(msg, type = "error")
            })
        })
        
        # Render per-assay results tabs
        output$assay_results_tabs <- shiny::renderUI({
            buildMetabIntensityAssayTabsUi(
                stats = filter_stats(),
                ns = ns
            )
        })
        
        # Render filter plot
        output$filter_plot <- shiny::renderPlot({
            renderMetabIntensityFilterPlot(filterPlot = filter_plot)
        })
    })
}
