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
#' @param id,workflowData,useStrictMode,minRepsPerGroup,minGroups,intensityCutoffPercentile,updateConfigParameterFn,updateMissingValueParametersFn,removeRowsWithMissingValuesPercentFn,logInfoFn,nowFn Runtime inputs used by this function; see the usage section for accepted values.
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
  if (usesDIAReferenceStates(workflowData$state_manager)) {
    pending <- if (identical(workflowData$state_manager$current_state, "protein_intensity_filtered"))
      resumeDIAProteinAction(workflowData, "protein_intensity_filtered") else NULL
    if (!is.null(pending)) return(pending)
    # KIT-015: an unmanaged launch runs the registered declaration in process with its production variant.
    entry <- workflowActionSpec("protein_intensity_filter")
    parameters <- list(strict_mode = isTRUE(useStrictMode), min_reps_per_group = minRepsPerGroup, min_groups = minGroups,
      proteins_intensity_cutoff_percentile = intensityCutoffPercentile)
    validateWorkflowActionParameters(entry$name, parameters, workflowExecutionPolicy(FALSE), resolved = FALSE)
    value <- computeDIAReferenceProteinIntensityFilter(workflowData$state_manager, getDIAStateReference(workflowData$state_manager),
      parameters, variant = entry$default_variant, now_fn = nowFn)
    return(finishDIAReferenceProteinStep(workflowData, value))
  }

  shiny::req(workflowData$state_manager)
  currentS4 <- readDIAProteinScience(workflowData$state_manager)
  predecessor <- getDIACurrentReference(workflowData$state_manager)
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

  filteredS4 <- inheritDIAExecutionSource(filteredS4, workflowData$state_manager, predecessor)
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

  groupwise_cutoff <- currentS4@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff
  max_groups_cutoff <- currentS4@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff
  releaseDIAWorkingReference(workflowData$state_manager, predecessor)
  currentS4 <- NULL
  proteinCount <- filteredS4@protein_quant_table |>
    dplyr::distinct(.data[[filteredS4@protein_id_column]]) |>
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
        groupwise_cutoff
      ),
      sprintf(
        "Max groups %% cutoff: %.3f%% (calculated)\n",
        max_groups_cutoff
      ),
      sprintf("Intensity cutoff percentile: %.1f%%\n", intensityCutoffPercentile),
      "State saved as: 'protein_intensity_filtered'\n"
    )
  }

  finishDIAProteinAction(workflowData, list(
    filteredS4 = filteredS4,
    resultText = resultText
  ), "filteredS4")
}

updateProteinIntensityFilterOutputs <- function(output,
                                                proteinIntensityFilterPlot,
                                                intensityFilterResult,
                                                omicType,
                                                experimentLabel,
                                                renderTextFn = shiny::renderText,
                                                updateProteinFilteringFn = updateProteinFiltering) {
  output$protein_intensity_filter_results <- renderPeptideQcResultText(intensityFilterResult$resultText, renderTextFn)

  plotGrid <- publishDIAProteinProgress(intensityFilterResult, "11_protein_intensity_filtered", omicType,
    experimentLabel, updateProteinFilteringFn, intensityFilterResult$filteredS4)
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
  withDIAStateScope(getDIAActiveManager(), function() {
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
  })
}

runProteinIntensityFilterRevertStep <- function(workflowData) {
  shiny::req(workflowData$state_manager)
  history <- workflowData$state_manager$getHistory()

  if (length(history) <= 1) {
    stop("No previous state to revert to.")
  }

  previousState <- history[length(history) - 1]
  revertedS4 <- selectDIAStateReference(workflowData$state_manager, previousState)

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
  withDIAStateScope(getDIAActiveManager(), function() {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$protein_intensity_filter_results <- renderPeptideQcResultText(revertResult$resultText, renderTextFn)
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
  })
}

bindProteinIntensityFilterPlot <- function(output, proteinIntensityFilterPlot, plotIdentity = NULL) {
  output$protein_intensity_filter_plot <- if (is.null(plotIdentity)) renderPlot({
    req(proteinIntensityFilterPlot())
    grid.draw(proteinIntensityFilterPlot())
  }) else renderPlot({
    req(proteinIntensityFilterPlot())
    grid.draw(proteinIntensityFilterPlot())
  }, alt = function() plotIdentity())
}

mod_prot_qc_protein_intensity_server <- function(id, workflow_data, experiment_paths = NULL, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    protein_intensity_filter_plot <- shiny::reactiveVal(NULL)
    protein_intensity_filter_plot_identity <- shiny::reactiveVal("Protein intensity filtering summary")
    protein_intensity_filter_plot_value <- protein_intensity_filter_plot

    # KIT-015: the managed worker when the session carries a managed application (KIT-008 pattern), else the
    # same registered declaration in process.
    managed_step <- installWorkflowManagedStep(session, input, workflow_data, "protein_intensity_filter",
      publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
      ui = list(text = "protein_intensity_filter_results", plot = "protein_intensity_filter_plot", ack = "workflow_ack", notification = "protein_intensity_filter_working",
        label = "protein intensity filter", step_label = "11_protein_intensity_filtered", plot_identity = protein_intensity_filter_plot_identity),
      deliver = function(result) updateProteinIntensityFilterOutputs(output = output, proteinIntensityFilterPlot = protein_intensity_filter_plot_value, intensityFilterResult = result,
        omicType = omic_type, experimentLabel = experiment_label),
      read_result = readManagedProteinResult)
    
    # Step 3: Protein Intensity Filter (chunks 20+21 combined)
    shiny::observeEvent(input$apply_protein_intensity_filter, {
      if (!is.null(workflowSessionApplication(session)) && usesDIAReferenceStates(workflow_data$state_manager)) {
        managed_step$apply(list(strict_mode = isTRUE(input$use_strict_mode), min_reps_per_group = input$min_reps_per_group,
          min_groups = input$min_groups, proteins_intensity_cutoff_percentile = input$proteins_intensity_cutoff_percentile))
        return(invisible(NULL))
      }
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
      managed_step$cancel("reverted")
      runProteinIntensityFilterRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })

    bindProteinIntensityFilterPlot(
      output = output,
      proteinIntensityFilterPlot = protein_intensity_filter_plot,
      plotIdentity = protein_intensity_filter_plot_identity
    )
  })
}
