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

#' @title Peptide Intensity Filter Module
#'
#' @description A Shiny module for applying peptide intensity filtering.
#'
#' @name mod_prot_qc_peptide_intensity
NULL

#' @rdname mod_prot_qc_peptide_intensity
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr numericInput helpText div actionButton verbatimTextOutput plotOutput
mod_prot_qc_peptide_intensity_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Intensity Filter",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Missing Value Parameters & Intensity Filter"),
          shiny::p("Configure missing value thresholds and filter peptides on intensity and missing values."),
          shiny::hr(),
          
          # Strict Mode Toggle
          shiny::checkboxInput(ns("use_strict_mode"), 
            "Strict Mode (No Missing Values Allowed)", 
            value = FALSE,
            width = "100%"
          ),
          shiny::helpText("When enabled, removes any peptide with missing values in ANY sample across all groups. Overrides flexible threshold settings below."),
          
          shiny::hr(),
          
          # Flexible mode parameters
          shiny::conditionalPanel(
            condition = paste0("!input['", ns("use_strict_mode"), "']"),
            
            shiny::h5("Missing Value Parameters (Flexible Mode)"),
            shiny::numericInput(ns("min_reps_per_group"), 
              "Min Replicates per Group", 
              value = 2, min = 1, max = 10, step = 1,
              width = "100%"
            ),
            shiny::helpText("Minimum valid measurements required within each group (default: 2)"),
            
            shiny::numericInput(ns("min_groups"), 
              "Min Groups Required", 
              value = 2, min = 1, max = 10, step = 1,
              width = "100%"
            ),
            shiny::helpText("Minimum groups that must meet the replicate threshold (default: 2)"),
            
            shiny::hr(),
            
            # Calculated percentages (read-only display)
            shiny::h5("Calculated Thresholds"),
            shiny::p(style = "background-color: #f0f0f0; padding: 10px; border-radius: 5px; color: #333;",
              shiny::strong("These values are automatically calculated from your replicate settings above:"),
              shiny::br(),
              shiny::textOutput(ns("calculated_groupwise_percent"), inline = TRUE),
              shiny::br(),
              shiny::textOutput(ns("calculated_max_groups_percent"), inline = TRUE)
            ),
            
            shiny::hr()
          ),
          
          # Intensity cutoff
          shiny::h5("Intensity Threshold"),
          shiny::numericInput(ns("intensity_cutoff_percentile"), 
            "Peptide Intensity Cutoff Percentile (%)", 
            value = 1, min = 0.1, max = 10, step = 0.1,
            width = "100%"
          ),
          shiny::helpText("Intensity threshold percentile (default: 1%)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_intensity_filter"), "Apply Filter", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_intensity"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("intensity_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("intensity_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_peptide_intensity
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
buildPeptideIntensityThresholdPreview <- function(workflowData,
                                                 minRepsPerGroup,
                                                 minGroups,
                                                 updateMissingValueParametersFn = updateMissingValueParameters) {
  shiny::req(workflowData$state_manager)
  currentS4 <- workflowData$state_manager$getState()
  shiny::req(currentS4)

  previewS4 <- updateMissingValueParametersFn(
    theObject = currentS4,
    min_reps_per_group = minRepsPerGroup,
    min_groups = minGroups,
    function_name = "peptideIntensityFiltering",
    grouping_variable = "group"
  )

  list(
    groupwiseText = sprintf(
      "Groupwise %% cutoff: %.3f%%",
      previewS4@args$peptideIntensityFiltering$groupwise_percentage_cutoff
    ),
    maxGroupsText = sprintf(
      "Max groups %% cutoff: %.3f%%",
      previewS4@args$peptideIntensityFiltering$max_groups_percentage_cutoff
    )
  )
}

runPeptideIntensityRevertStep <- function(workflowData) {
  history <- workflowData$state_manager$getHistory()
  if (length(history) <= 1) {
    stop("No previous state to revert to.")
  }

  previousState <- history[length(history) - 1]
  revertedS4 <- workflowData$state_manager$revertToState(previousState)
  logger::log_info(paste("Reverted intensity filtering to", previousState))

  list(
    previousState = previousState,
    revertedS4 = revertedS4,
    resultText = paste("Reverted to previous state:", previousState)
  )
}

runPeptideIntensityRevertObserver <- function(workflowData,
                                              output,
                                              runRevertStepFn = runPeptideIntensityRevertStep,
                                              renderTextFn = shiny::renderText,
                                              showNotificationFn = shiny::showNotification,
                                              logErrorFn = logger::log_error) {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$intensity_results <- renderTextFn(revertResult$resultText)
    showNotificationFn("Reverted successfully", type = "message")

    list(
      status = "success",
      revertResult = revertResult
    )
  }, error = function(e) {
    errorMessage <- paste("Error reverting:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
}

runPeptideIntensityApplyStep <- function(workflowData,
                                         useStrictMode,
                                         minRepsPerGroup,
                                         minGroups,
                                         intensityCutoffPercentile,
                                         updateConfigParameterFn = updateConfigParameter,
                                         updateMissingValueParametersFn = updateMissingValueParameters,
                                         peptideIntensityFilteringFn = peptideIntensityFiltering,
                                         logInfoFn = logger::log_info,
                                         nowFn = Sys.time) {
  shiny::req(workflowData$state_manager)
  currentS4 <- workflowData$state_manager$getState()
  shiny::req(currentS4)

  if (useStrictMode) {
    logInfoFn("Peptide Processing: Using STRICT MODE")

    currentS4 <- updateConfigParameterFn(
      theObject = currentS4,
      function_name = "peptideIntensityFiltering",
      parameter_name = "groupwise_percentage_cutoff",
      new_value = 0
    )
    currentS4 <- updateConfigParameterFn(
      theObject = currentS4,
      function_name = "peptideIntensityFiltering",
      parameter_name = "max_groups_percentage_cutoff",
      new_value = 0
    )
  } else {
    logInfoFn("Peptide Processing: Using FLEXIBLE MODE")

    currentS4 <- updateMissingValueParametersFn(
      theObject = currentS4,
      min_reps_per_group = minRepsPerGroup,
      min_groups = minGroups,
      function_name = "peptideIntensityFiltering",
      grouping_variable = "group"
    )
  }

  currentS4 <- updateConfigParameterFn(
    theObject = currentS4,
    function_name = "peptideIntensityFiltering",
    parameter_name = "peptides_intensity_cutoff_percentile",
    new_value = intensityCutoffPercentile
  )

  filteredS4 <- peptideIntensityFilteringFn(theObject = currentS4)

  if (is.null(workflowData$qc_params)) {
    workflowData$qc_params <- list()
  }
  if (is.null(workflowData$qc_params$peptide_qc)) {
    workflowData$qc_params$peptide_qc <- list()
  }

  workflowData$qc_params$peptide_qc$intensity_filter <- list(
    strict_mode = useStrictMode,
    min_reps_per_group = if (!useStrictMode) minRepsPerGroup else NA,
    min_groups = if (!useStrictMode) minGroups else NA,
    intensity_cutoff_percentile = intensityCutoffPercentile,
    timestamp = nowFn()
  )

  workflowData$state_manager$saveState(
    state_name = "intensity_filtered",
    s4_data_object = filteredS4,
    config_object = list(
      strict_mode = useStrictMode,
      intensity_cutoff_percentile = intensityCutoffPercentile
    ),
    description = if (useStrictMode) {
      "Applied STRICT peptide intensity filter"
    } else {
      "Applied FLEXIBLE peptide intensity filter"
    }
  )

  proteinCount <- filteredS4@peptide_data |>
    dplyr::distinct(Protein.Ids) |>
    nrow()

  resultText <- paste(
    "Intensity Filter Applied Successfully\n",
    "====================================\n",
    sprintf("Mode: %s\n", if (useStrictMode) "STRICT" else "FLEXIBLE"),
    sprintf("Proteins remaining: %d\n", proteinCount),
    sprintf("Intensity cutoff percentile: %.1f%%\n", intensityCutoffPercentile),
    sprintf(
      "Groupwise %% cutoff: %.3f%%\n",
      currentS4@args$peptideIntensityFiltering$groupwise_percentage_cutoff
    ),
    sprintf(
      "Max groups %% cutoff: %.3f%%\n",
      currentS4@args$peptideIntensityFiltering$max_groups_percentage_cutoff
    ),
    "State saved as: 'intensity_filtered'\n"
  )

  list(
    filteredS4 = filteredS4,
    resultText = resultText
  )
}

updatePeptideIntensityOutputs <- function(output,
                                          intensityPlot,
                                          intensityResult,
                                          omicType,
                                          experimentLabel,
                                          renderTextFn = shiny::renderText,
                                          updateProteinFilteringFn = updateProteinFiltering) {
  output$intensity_results <- renderTextFn(intensityResult$resultText)

  plotGrid <- updateProteinFilteringFn(
    data = intensityResult$filteredS4@peptide_data,
    step_name = "4_intensity_filtered",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
  intensityPlot(plotGrid)

  invisible(plotGrid)
}

runPeptideIntensityApplyObserver <- function(workflowData,
                                             useStrictMode,
                                             minRepsPerGroup,
                                             minGroups,
                                             intensityCutoffPercentile,
                                             output,
                                             intensityPlot,
                                             omicType,
                                             experimentLabel,
                                             runApplyStepFn = runPeptideIntensityApplyStep,
                                             updateOutputsFn = updatePeptideIntensityOutputs,
                                             showNotificationFn = shiny::showNotification,
                                             removeNotificationFn = shiny::removeNotification,
                                             logErrorFn = logger::log_error) {
  showNotificationFn(
    "Applying intensity filter...",
    id = "intensity_working",
    duration = NULL
  )

  tryCatch({
    intensityResult <- runApplyStepFn(
      workflowData = workflowData,
      useStrictMode = useStrictMode,
      minRepsPerGroup = minRepsPerGroup,
      minGroups = minGroups,
      intensityCutoffPercentile = intensityCutoffPercentile
    )

    plotGrid <- updateOutputsFn(
      output = output,
      intensityPlot = intensityPlot,
      intensityResult = intensityResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )

    removeNotificationFn("intensity_working")
    showNotificationFn("Intensity filter applied successfully", type = "message")

    list(
      status = "success",
      intensityResult = intensityResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error applying intensity filter:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("intensity_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
}

setupPeptideIntensityServerBootstrap <- function(input,
                                                 output,
                                                 workflowData,
                                                 intensityPlot,
                                                 omicType,
                                                 experimentLabel,
                                                 buildThresholdPreviewFn = buildPeptideIntensityThresholdPreview,
                                                 runApplyObserverFn = runPeptideIntensityApplyObserver,
                                                 runRevertObserverFn = runPeptideIntensityRevertObserver,
                                                 reactiveFn = shiny::reactive,
                                                 renderTextFn = shiny::renderText,
                                                 observeEventFn = shiny::observeEvent,
                                                 renderPlotFn = shiny::renderPlot,
                                                 reqFn = shiny::req,
                                                 gridDrawFn = grid::grid.draw) {
  thresholdPreview <- reactiveFn({
    buildThresholdPreviewFn(
      workflowData = workflowData,
      minRepsPerGroup = input$min_reps_per_group,
      minGroups = input$min_groups
    )
  })

  output$calculated_groupwise_percent <- renderTextFn({
    thresholdPreview()$groupwiseText
  })

  output$calculated_max_groups_percent <- renderTextFn({
    thresholdPreview()$maxGroupsText
  })

  observeEventFn(input$apply_intensity_filter, {
    runApplyObserverFn(
      workflowData = workflowData,
      useStrictMode = isTRUE(input$use_strict_mode),
      minRepsPerGroup = input$min_reps_per_group,
      minGroups = input$min_groups,
      intensityCutoffPercentile = input$intensity_cutoff_percentile,
      output = output,
      intensityPlot = intensityPlot,
      omicType = omicType,
      experimentLabel = experimentLabel
    )
  })

  observeEventFn(input$revert_intensity, {
    runRevertObserverFn(
      workflowData = workflowData,
      output = output
    )
  })

  output$intensity_plot <- renderPlotFn({
    reqFn(intensityPlot())
    gridDrawFn(intensityPlot())
  })

  list(
    thresholdPreview = thresholdPreview,
    reason = "wired"
  )
}

mod_prot_qc_peptide_intensity_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    intensity_plot <- shiny::reactiveVal(NULL)

    setupPeptideIntensityServerBootstrap(
      input = input,
      output = output,
      workflowData = workflow_data,
      intensityPlot = intensity_plot,
      omicType = omic_type,
      experimentLabel = experiment_label
    )
  })
}
