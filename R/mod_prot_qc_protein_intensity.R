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

#' @title Protein Intensity Filter Module
#'
#' @description A Shiny module for applying protein intensity and missing value filters.
#'
#' @name mod_prot_qc_protein_intensity
NULL

#' @rdname mod_prot_qc_protein_intensity
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr checkboxInput helpText conditionalPanel h5 numericInput div actionButton verbatimTextOutput plotOutput strong
mod_prot_qc_protein_intensity_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Protein Intensity Filter",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Missing Value Parameters & Intensity Filter"),
          shiny::p("Configure missing value thresholds and filter proteins on intensity and missing values."),
          shiny::hr(),
          
          # Strict Mode Toggle
          shiny::checkboxInput(ns("use_strict_mode"), 
            "Strict Mode (No Missing Values Allowed)", 
            value = FALSE,
            width = "100%"
          ),
          shiny::helpText("When enabled, removes any protein with missing values in ANY sample across all groups. Overrides flexible threshold settings below."),
          
          shiny::hr(),
          
          # Flexible mode parameters (shown when strict mode is OFF)
          shiny::conditionalPanel(
            condition = paste0("!input['", ns("use_strict_mode"), "']"),
            
            # Missing value parameters (chunk 20)
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
            shiny::p(style = "background-color: #f0f0f0; padding: 10px; border-radius: 5px;",
              shiny::strong("These values are automatically calculated from your replicate settings above:"),
              shiny::br(),
              shiny::textOutput(ns("calculated_groupwise_percent"), inline = TRUE),
              shiny::br(),
              shiny::textOutput(ns("calculated_max_groups_percent"), inline = TRUE)
            ),
            
            shiny::hr()
          ),
          
          # Intensity cutoff (always shown)
          shiny::h5("Intensity Threshold"),
          shiny::numericInput(ns("proteins_intensity_cutoff_percentile"), 
            "Protein Intensity Cutoff Percentile (%)", 
            value = 1, min = 0.1, max = 10, step = 0.1,
            width = "100%"
          ),
          shiny::helpText("Intensity threshold percentile (default: 1%)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_protein_intensity_filter"), "Apply Filter", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_protein_intensity_filter"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("protein_intensity_filter_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("protein_intensity_filter_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_protein_intensity
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
runProteinIntensityFilterApplyStep <- function(workflowData,
                                               useStrictMode,
                                               minRepsPerGroup,
                                               minGroups,
                                               intensityCutoffPercentile,
                                               updateConfigParameterFn = updateConfigParameter,
                                               updateMissingValueParametersFn = updateMissingValueParameters,
                                               removeRowsWithMissingValuesPercentFn = removeRowsWithMissingValuesPercent,
                                               logInfoFn = logger::log_info,
                                               nowFn = Sys.time) {
  shiny::req(workflowData$state_manager)
  currentS4 <- workflowData$state_manager$getState()
  shiny::req(currentS4)

  if (useStrictMode) {
    logInfoFn("Protein Processing: Using STRICT MODE (no missing values allowed)")

    currentS4 <- updateConfigParameterFn(
      theObject = currentS4,
      function_name = "removeRowsWithMissingValuesPercent",
      parameter_name = "groupwise_percentage_cutoff",
      new_value = 0
    )
    currentS4 <- updateConfigParameterFn(
      theObject = currentS4,
      function_name = "removeRowsWithMissingValuesPercent",
      parameter_name = "max_groups_percentage_cutoff",
      new_value = 0
    )
  } else {
    logInfoFn("Protein Processing: Using FLEXIBLE MODE (adaptive thresholds)")

    currentS4 <- updateMissingValueParametersFn(
      theObject = currentS4,
      min_reps_per_group = minRepsPerGroup,
      min_groups = minGroups
    )
  }

  currentS4 <- updateConfigParameterFn(
    theObject = currentS4,
    function_name = "removeRowsWithMissingValuesPercent",
    parameter_name = "proteins_intensity_cutoff_percentile",
    new_value = intensityCutoffPercentile
  )

  filteredS4 <- removeRowsWithMissingValuesPercentFn(currentS4)

  if (is.null(workflowData$qc_params)) {
    workflowData$qc_params <- list()
  }
  if (is.null(workflowData$qc_params$protein_qc)) {
    workflowData$qc_params$protein_qc <- list()
  }

  workflowData$qc_params$protein_qc$intensity_filter <- list(
    strict_mode = useStrictMode,
    min_reps_per_group = if (!useStrictMode) minRepsPerGroup else NA,
    min_groups = if (!useStrictMode) minGroups else NA,
    groupwise_percentage_cutoff = currentS4@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff,
    max_groups_percentage_cutoff = currentS4@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff,
    proteins_intensity_cutoff_percentile = intensityCutoffPercentile,
    timestamp = nowFn()
  )

  workflowData$state_manager$saveState(
    state_name = "protein_intensity_filtered",
    s4_data_object = filteredS4,
    config_object = list(
      strict_mode = useStrictMode,
      min_reps_per_group = if (!useStrictMode) minRepsPerGroup else NA,
      min_groups = if (!useStrictMode) minGroups else NA,
      groupwise_percentage_cutoff = currentS4@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff,
      max_groups_percentage_cutoff = currentS4@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff,
      proteins_intensity_cutoff_percentile = intensityCutoffPercentile
    ),
    description = if (useStrictMode) {
      "Applied STRICT protein intensity filter (no missing values)"
    } else {
      "Applied FLEXIBLE protein intensity filter (adaptive thresholds)"
    }
  )

  proteinCount <- filteredS4@protein_quant_table |>
    dplyr::distinct(Protein.Ids) |>
    nrow()

  resultText <- if (useStrictMode) {
    paste(
      "Protein Intensity Filter Applied Successfully\n",
      "============================================\n",
      "Mode: STRICT (No Missing Values)\n",
      sprintf("Proteins remaining: %d\n", proteinCount),
      "Groupwise % cutoff: 0.000% (strict - no missing allowed)\n",
      "Max groups % cutoff: 0.000% (strict - all groups must pass)\n",
      sprintf("Intensity cutoff percentile: %.1f%%\n", intensityCutoffPercentile),
      "State saved as: 'protein_intensity_filtered'\n"
    )
  } else {
    paste(
      "Protein Intensity Filter Applied Successfully\n",
      "============================================\n",
      "Mode: FLEXIBLE (Adaptive Thresholds)\n",
      sprintf("Proteins remaining: %d\n", proteinCount),
      sprintf("Min replicates per group: %d\n", minRepsPerGroup),
      sprintf("Min groups required: %d\n", minGroups),
      sprintf(
        "Groupwise %% cutoff: %.3f%% (calculated)\n",
        currentS4@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff
      ),
      sprintf(
        "Max groups %% cutoff: %.3f%% (calculated)\n",
        currentS4@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff
      ),
      sprintf("Intensity cutoff percentile: %.1f%%\n", intensityCutoffPercentile),
      "State saved as: 'protein_intensity_filtered'\n"
    )
  }

  list(
    filteredS4 = filteredS4,
    resultText = resultText
  )
}

updateProteinIntensityFilterOutputs <- function(output,
                                                proteinIntensityFilterPlot,
                                                intensityFilterResult,
                                                omicType,
                                                experimentLabel,
                                                renderTextFn = shiny::renderText,
                                                updateProteinFilteringFn = updateProteinFiltering) {
  output$protein_intensity_filter_results <- renderTextFn(intensityFilterResult$resultText)

  plotGrid <- updateProteinFilteringFn(
    data = intensityFilterResult$filteredS4@protein_quant_table,
    step_name = "11_protein_intensity_filtered",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
  proteinIntensityFilterPlot(plotGrid)

  invisible(plotGrid)
}

runProteinIntensityFilterApplyObserver <- function(workflowData,
                                                   useStrictMode,
                                                   minRepsPerGroup,
                                                   minGroups,
                                                   intensityCutoffPercentile,
                                                   output,
                                                   proteinIntensityFilterPlot,
                                                   omicType,
                                                   experimentLabel,
                                                   runApplyStepFn = runProteinIntensityFilterApplyStep,
                                                   updateOutputsFn = updateProteinIntensityFilterOutputs,
                                                   showNotificationFn = shiny::showNotification,
                                                   removeNotificationFn = shiny::removeNotification,
                                                   logInfoFn = logger::log_info,
                                                   logErrorFn = logger::log_error) {
  showNotificationFn(
    "Applying protein intensity filter...",
    id = "protein_intensity_filter_working",
    duration = NULL
  )

  tryCatch({
    intensityFilterResult <- runApplyStepFn(
      workflowData = workflowData,
      useStrictMode = useStrictMode,
      minRepsPerGroup = minRepsPerGroup,
      minGroups = minGroups,
      intensityCutoffPercentile = intensityCutoffPercentile
    )

    plotGrid <- updateOutputsFn(
      output = output,
      proteinIntensityFilterPlot = proteinIntensityFilterPlot,
      intensityFilterResult = intensityFilterResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )

    logInfoFn("Protein intensity filter applied successfully")
    removeNotificationFn("protein_intensity_filter_working")
    showNotificationFn("Protein intensity filter applied successfully", type = "message")

    list(
      status = "success",
      intensityFilterResult = intensityFilterResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error applying protein intensity filter:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("protein_intensity_filter_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
}

runProteinIntensityFilterRevertStep <- function(workflowData) {
  shiny::req(workflowData$state_manager)
  history <- workflowData$state_manager$getHistory()

  if (length(history) <= 1) {
    stop("No previous state to revert to.")
  }

  previousState <- history[length(history) - 1]
  revertedS4 <- workflowData$state_manager$revertToState(previousState)

  list(
    previousState = previousState,
    revertedS4 = revertedS4,
    resultText = paste("Reverted to previous state:", previousState)
  )
}

runProteinIntensityFilterRevertObserver <- function(workflowData,
                                                    output,
                                                    runRevertStepFn = runProteinIntensityFilterRevertStep,
                                                    renderTextFn = shiny::renderText,
                                                    showNotificationFn = shiny::showNotification,
                                                    logInfoFn = logger::log_info,
                                                    logErrorFn = logger::log_error) {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$protein_intensity_filter_results <- renderTextFn(revertResult$resultText)
    logInfoFn(paste("Reverted protein intensity filter to", revertResult$previousState))
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

bindProteinIntensityFilterPlot <- function(output, proteinIntensityFilterPlot) {
  output$protein_intensity_filter_plot <- renderPlot({
    req(proteinIntensityFilterPlot())
    grid.draw(proteinIntensityFilterPlot())
  })
}

mod_prot_qc_protein_intensity_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    protein_intensity_filter_plot <- shiny::reactiveVal(NULL)
    
    # Step 3: Protein Intensity Filter (chunks 20+21 combined)
    shiny::observeEvent(input$apply_protein_intensity_filter, {
      runProteinIntensityFilterApplyObserver(
        workflowData = workflow_data,
        useStrictMode = isTRUE(input$use_strict_mode),
        minRepsPerGroup = input$min_reps_per_group,
        minGroups = input$min_groups,
        intensityCutoffPercentile = input$proteins_intensity_cutoff_percentile,
        output = output,
        proteinIntensityFilterPlot = protein_intensity_filter_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
    })
    
    # Revert Protein Intensity Filter
    shiny::observeEvent(input$revert_protein_intensity_filter, {
      runProteinIntensityFilterRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })

    bindProteinIntensityFilterPlot(
      output = output,
      proteinIntensityFilterPlot = protein_intensity_filter_plot
    )
  })
}
