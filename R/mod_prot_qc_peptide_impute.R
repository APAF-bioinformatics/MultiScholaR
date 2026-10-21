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

#' @title Peptide Imputation Module
#'
#' @description A Shiny module for applying missing value imputation.
#'
#' @name mod_prot_qc_peptide_impute
#' @param id,imputed_s4,proportion_missing_values,exclusion_column Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_prot_qc_peptide_impute
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr numericInput textInput helpText div actionButton verbatimTextOutput plotOutput
mod_prot_qc_peptide_impute_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Imputation",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Missing Value Imputation"),
          shiny::p("Impute missing values using technical replicate averages."),
          shiny::hr(),
          
          shiny::numericInput(ns("proportion_missing_values"), 
            "Max Proportion Missing", 
            value = 0.5, min = 0.1, max = 1.0, step = 0.1,
            width = "100%"
          ),
          # KIT-049: the replicate-mean imputation is optional; skipping keeps missing values for IQ.
          shiny::checkboxInput(ns("skip_imputation"), "Skip imputation (keep missing values for IQ)", value = FALSE),
          shiny::helpText("A peptide is eligible when its missing fraction is less than or equal to this value (default: 0.5)."),
          shiny::textInput(
            ns("exclusion_column"),
            "Design exclusion flag column (optional)",
            value = "",
            width = "100%"
          ),
          shiny::helpText("Optional design column containing explicit TRUE/FALSE flags. Sample names are never excluded by pattern by default."),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_imputation"), "Apply Imputation", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_imputation"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("imputation_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("imputation_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_peptide_impute
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
.peptideImputationAuditSummary <- function(imputed_s4,
                                           proportion_missing_values,
                                           exclusion_column = NULL) {
  component <- isDIAComponentCheckpoint(imputed_s4)
  section <- tryCatch(
    if (component) imputed_s4$slots$args$peptideMissingValueImputation else
      imputed_s4@args$peptideMissingValueImputation,
    error = function(error) NULL
  )
  if (!is.list(section)) {
    section <- list()
  }
  summary <- section$imputation_summary
  if (!is.list(summary)) {
    summary <- list()
  }
  value_or <- function(value, fallback) {
    if (is.null(value) || length(value) == 0L) fallback else value
  }
  imputed_count <- if (component) summary$imputed_rows %||% NA_integer_ else if ("is_imputed" %in% names(imputed_s4@peptide_data)) {
    sum(imputed_s4@peptide_data$is_imputed %in% TRUE)
  } else {
    NA_integer_
  }

  list(
    status = value_or(summary$status, "applied"),
    skip_reason = value_or(summary$skip_reason, NULL),
    quantity_column = value_or(summary$quantity_column, value_or(section$quantity_column, NA_character_)),
    imputed_value_column = value_or(summary$imputed_value_column, value_or(section$imputed_value_column, NA_character_)),
    eligibility_denominator = value_or(
      summary$eligibility_denominator,
      "distinct_design_runs_in_technical_replicate_group"
    ),
    eligibility_operator = value_or(summary$eligibility_operator, "<="),
    maximum_missing_fraction = value_or(summary$maximum_missing_fraction, proportion_missing_values),
    exclusion_source = value_or(
      summary$exclusion_source,
      if (is.null(exclusion_column)) "none" else paste0("design_column:", exclusion_column)
    ),
    excluded_runs = value_or(summary$excluded_runs, character(0)),
    observed_rows = value_or(summary$observed_rows, NA_integer_),
    missing_not_imputed_rows = value_or(summary$missing_not_imputed_rows, NA_integer_),
    imputed_rows = value_or(summary$imputed_rows, as.integer(imputed_count)),
    input_rows = value_or(summary$input_rows, NA_integer_),
    output_rows = value_or(summary$output_rows, if (component) imputed_s4$dimensions$rows else nrow(imputed_s4@peptide_data))
  )
}

runPeptideImputationStep <- function(workflow_data,
                                     proportionMissingValues,
                                     exclusionColumn = NULL,
                                     updateConfigParameterFn = updateConfigParameter,
                                     skipImputation = FALSE) {
  if (usesDIAReferenceStates(workflow_data$state_manager)) {
    parameters <- list(proportion_missing_values = proportionMissingValues, exclusion_column = exclusionColumn)
    if (isTRUE(skipImputation)) parameters$skip <- TRUE
    return(runDIAReferencePeptideGroups(workflow_data, "peptideMissingValueImputation", parameters, "imputed",
      if (isTRUE(skipImputation)) "Missing value imputation skipped by user (missing values kept for IQ)" else
        "Applied missing value imputation using technical replicates"))
  }

  current_s4 <- getDIAExecutionState(workflow_data$state_manager, "peptideMissingValueImputation")
  predecessor <- getDIACurrentReference(workflow_data$state_manager)
  shiny::req(current_s4)

  if (length(exclusionColumn) == 0L || is.na(exclusionColumn[[1L]]) ||
      !nzchar(trimws(as.character(exclusionColumn[[1L]])))) {
    exclusionColumn <- NULL
  } else {
    exclusionColumn <- trimws(as.character(exclusionColumn[[1L]]))
  }

  logger::log_info(sprintf(
    "QC Step: Applying missing value imputation (proportion: %s)",
    proportionMissingValues
  ))

  current_s4 <- updateConfigParameterFn(
    theObject = current_s4,
    function_name = "peptideMissingValueImputation",
    parameter_name = "proportion_missing_values",
    new_value = proportionMissingValues
  )
  if (!is.null(exclusionColumn)) {
    current_s4 <- updateConfigParameterFn(
      theObject = current_s4,
      function_name = "peptideMissingValueImputation",
      parameter_name = "exclusion_column",
      new_value = exclusionColumn
    )
  }

  imputed_s4 <- if (isTRUE(skipImputation)) {
    peptideMissingValueImputationSkipped(current_s4, proportionMissingValues, exclusionColumn)
  } else {
    peptideMissingValueImputation(
      theObject = current_s4,
      proportion_missing_values = proportionMissingValues,
      exclusion_column = exclusionColumn
    )
  }
  imputation_audit <- .peptideImputationAuditSummary(
    imputed_s4,
    proportion_missing_values = proportionMissingValues,
    exclusion_column = exclusionColumn
  )

  if (is.null(workflow_data$qc_params)) {
    workflow_data$qc_params <- list()
  }
  if (is.null(workflow_data$qc_params$peptide_qc)) {
    workflow_data$qc_params$peptide_qc <- list()
  }

  workflow_data$qc_params$peptide_qc$imputation <- list(
    proportion_missing_values = proportionMissingValues,
    exclusion_column = exclusionColumn,
    imputation_summary = imputation_audit,
    timestamp = Sys.time()
  )

  tryCatch({
    if (exists("experiment_paths") && !is.null(experiment_paths$source_dir)) {
      qc_params_file <- file.path(experiment_paths$source_dir, "qc_params.RDS")
      saveRDS(workflow_data$qc_params, qc_params_file)
      logger::log_info(sprintf("Saved QC parameters to: %s", qc_params_file))
    }
  }, error = function(e) {
    logger::log_warn(sprintf("Could not save QC parameters file: %s", e$message))
  })

  imputationConfig <- list(
    proportion_missing_values = proportionMissingValues,
    exclusion_column = exclusionColumn,
    imputation_summary = imputation_audit
  )
  imputed_s4 <- .savePeptideQcState(
    state_manager = workflow_data$state_manager,
    before = current_s4,
    after = imputed_s4,
    stage_id = "imputation",
    state_name = "imputed",
    config_object = imputationConfig,
    description = if (identical(imputation_audit$status, "skipped")) {
      paste("Skipped technical-replicate imputation:", imputation_audit$skip_reason)
    } else {
      "Applied missing value imputation using technical replicates"
    },
    status = if (identical(imputation_audit$status, "skipped")) "skipped" else "applied",
    decision_reason = imputation_audit$skip_reason %||% NA_character_,
    transformation_type = "imputation"
  )

  # --- TESTTHAT CHECKPOINT CP02 (see test-prot-02-qc-filtering.R) ---
  releaseDIAWorkingReference(workflow_data$state_manager, predecessor)
  current_s4 <- NULL
  .capture_checkpoint(imputed_s4, "cp02", "qc_filtered_peptide")
  # --- END CP02 ---

  protein_count <- .countPeptideProteinGroups(imputed_s4)

  status_heading <- if (identical(imputation_audit$status, "skipped")) {
    "Missing Value Imputation Skipped\n"
  } else {
    "Missing Value Imputation Applied Successfully\n"
  }
  audit_lines <- c(
    sprintf("Eligibility: missing fraction %s %.1f\n", imputation_audit$eligibility_operator, proportionMissingValues),
    sprintf("Denominator: %s\n", imputation_audit$eligibility_denominator),
    sprintf("Quantity column: %s\n", imputation_audit$quantity_column),
    sprintf("Exclusion source: %s\n", imputation_audit$exclusion_source)
  )
  if (!is.na(imputation_audit$imputed_rows)) {
    audit_lines <- c(audit_lines, sprintf("Rows imputed: %d\n", imputation_audit$imputed_rows))
  }
  if (identical(imputation_audit$status, "skipped")) {
    audit_lines <- c(audit_lines, sprintf("Skip reason: %s\n", imputation_audit$skip_reason))
  }

  result_text <- paste(
    status_heading,
    "============================================\n",
    sprintf("Proteins remaining: %d\n", protein_count),
    sprintf("Max proportion missing: %.1f\n", proportionMissingValues),
    paste(audit_lines, collapse = ""),
    "State saved as: 'imputed'\n",
    "\nReady for peptide-to-protein rollup step."
  )

  list(
    imputedS4 = imputed_s4,
    resultText = result_text
  )
}

runPeptideImputationRevertStep <- function(workflow_data) {
  history <- workflow_data$state_manager$getHistory()
  if (length(history) <= 1) {
    stop("No previous state to revert to.")
  }

  previousState <- history[length(history) - 1]
  revertedS4 <- selectDIAStateReference(workflow_data$state_manager, previousState)
  logger::log_info(paste("Reverted imputation to", previousState))

  list(
    previousState = previousState,
    revertedS4 = revertedS4,
    resultText = paste("Reverted to previous state:", previousState)
  )
}

updatePeptideImputationOutputs <- function(output, imputationPlot, imputationResult, omicType, experimentLabel, projectDirs = NULL) {
  output$imputation_results <- renderPeptideQcResultText(imputationResult$resultText, renderText)

  plot_arguments <- list(
    data = imputationResult$filteringSummary %||% imputationResult$imputedS4@peptide_data,
    step_name = "8_imputed",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
  if (!is.null(projectDirs)) plot_arguments$project_dirs <- projectDirs
  plotGrid <- if (!is.null(imputationResult$stateReference)) {
    publishDIAQvalueOutput(getDIAActiveManager(), imputationResult, plot_arguments, updateProteinFiltering)
  } else do.call(updateProteinFiltering, plot_arguments)
  imputationPlot(plotGrid)

  invisible(plotGrid)
}

runPeptideImputationApplyObserver <- function(workflow_data,
                                              proportionMissingValues,
                                              output,
                                              imputationPlot,
                                              omicType,
                                              experimentLabel,
                                              exclusionColumn = NULL,
                                              skipImputation = FALSE,
                                              runImputationStepFn = runPeptideImputationStep,
                                              updateOutputsFn = updatePeptideImputationOutputs,
                                              showNotificationFn = shiny::showNotification,
                                              removeNotificationFn = shiny::removeNotification,
                                              logInfoFn = logger::log_info,
                                              logErrorFn = logger::log_error, experimentPaths = NULL) {
  withDIAStateScope(getDIAActiveManager(), function() {
  showNotificationFn(
    "Applying missing value imputation...",
    id = "imputation_working",
    duration = NULL
  )

  tryCatch({
    run_arguments <- list(
      workflow_data = workflow_data,
      proportionMissingValues = proportionMissingValues
    )
    if (!is.null(exclusionColumn) && length(exclusionColumn) > 0L &&
        !is.na(exclusionColumn[[1L]]) && nzchar(trimws(as.character(exclusionColumn[[1L]])))) {
      run_arguments$exclusionColumn <- exclusionColumn
    }
    if (isTRUE(skipImputation)) run_arguments$skipImputation <- TRUE
    imputationResult <- do.call(runImputationStepFn, run_arguments)

    output_arguments <- list(
      output = output,
      imputationPlot = imputationPlot,
      imputationResult = imputationResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )
    if (!is.null(experimentPaths)) {
      output_arguments$projectDirs <- stats::setNames(
        list(experimentPaths), paste(omicType, experimentLabel, sep = "_")
      )
    }
    plotGrid <- do.call(updateOutputsFn, output_arguments)

    logInfoFn("Missing value imputation applied successfully")
    removeNotificationFn("imputation_working")
    showNotificationFn(
      "Missing value imputation applied successfully",
      type = "message"
    )

    list(
      status = "success",
      imputationResult = imputationResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error applying missing value imputation:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("imputation_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
  })
}

runPeptideImputationRevertObserver <- function(workflow_data,
                                               output,
                                               runRevertStepFn = runPeptideImputationRevertStep,
                                               renderTextFn = shiny::renderText,
                                               showNotificationFn = shiny::showNotification,
                                               logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  tryCatch({
    revertResult <- runRevertStepFn(workflow_data = workflow_data)
    output$imputation_results <- renderPeptideQcResultText(revertResult$resultText, renderTextFn)
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

bindPeptideImputationPlot <- function(output, imputationPlot, plotIdentity = NULL) {
  output$imputation_plot <- if (is.null(plotIdentity)) renderPlot({
    req(imputationPlot())
    grid.draw(imputationPlot())
  }) else renderPlot({
    req(imputationPlot())
    grid.draw(imputationPlot())
  }, alt = function() plotIdentity())
}

mod_prot_qc_peptide_impute_server <- function(id, workflow_data, omic_type, experiment_label, experiment_paths = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    
    imputation_plot <- shiny::reactiveVal(NULL)
    imputation_plot_identity <- shiny::reactiveVal("Missing value imputation summary")

    # KIT-014: the managed worker when the session carries a managed application (KIT-008 pattern), else the
    # same registered declaration in process.
    managed_step <- installWorkflowManagedStep(session, input, workflow_data, "peptide_imputation",
      publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
      ui = list(text = "imputation_results", plot = "imputation_plot", ack = "workflow_ack", notification = "imputation_working",
        label = "missing value imputation", step_label = "8_imputed", plot_identity = imputation_plot_identity),
      deliver = function(result) {
        arguments <- list(output = output, imputationPlot = imputation_plot, imputationResult = result, omicType = omic_type, experimentLabel = experiment_label)
        if (!is.null(experiment_paths)) arguments$projectDirs <- stats::setNames(list(experiment_paths),
          paste(omic_type, experiment_label, sep = "_"))
        do.call(updatePeptideImputationOutputs, arguments)
      })
    
    # Step 7: Missing Value Imputation (chunk 16)
    shiny::observeEvent(input$apply_imputation, {
      shiny::req(workflow_data$state_manager)
      if (!is.null(workflowSessionApplication(session)) && usesDIAReferenceStates(workflow_data$state_manager)) {
        parameters <- list(proportion_missing_values = input$proportion_missing_values)
        exclusion <- input$exclusion_column
        if (!is.null(exclusion) && length(exclusion) > 0L && !is.na(exclusion[[1L]]) &&
            nzchar(trimws(as.character(exclusion[[1L]])))) parameters$exclusion_column <- trimws(as.character(exclusion[[1L]]))
        if (isTRUE(input$skip_imputation)) parameters$skip <- TRUE
        managed_step$apply(parameters)
        return(invisible(NULL))
      }

      apply_arguments <- list(
        workflow_data = workflow_data,
        proportionMissingValues = input$proportion_missing_values,
        output = output,
        imputationPlot = imputation_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
      if (!is.null(input$exclusion_column) && length(input$exclusion_column) > 0L &&
          !is.na(input$exclusion_column[[1L]]) &&
          nzchar(trimws(as.character(input$exclusion_column[[1L]])))) {
        apply_arguments$exclusionColumn <- input$exclusion_column
      }
      if (!is.null(experiment_paths)) apply_arguments$experimentPaths <- experiment_paths
      if (isTRUE(input$skip_imputation)) apply_arguments$skipImputation <- TRUE
      do.call(runPeptideImputationApplyObserver, apply_arguments)
    })
    
    # Revert Imputation
    shiny::observeEvent(input$revert_imputation, {
      managed_step$cancel("reverted")
      runPeptideImputationRevertObserver(
        workflow_data = workflow_data,
        output = output
      )
    })
    
    bindPeptideImputationPlot(output = output, imputationPlot = imputation_plot, plotIdentity = imputation_plot_identity)
  })
}
