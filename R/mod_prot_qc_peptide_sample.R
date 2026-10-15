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

#' @title Peptide Sample Filter Module
#'
#' @description A Shiny module for applying the sample quality filter.
#'
#' @name mod_prot_qc_peptide_sample
#' @param id,workflowData,minPeptidesPerSample,updateConfigParameterFn,filterMinNumPeptidesPerSampleFn,logInfoFn,nowFn Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' Register the PERF-107 browser acknowledgement before dynamic QC UI insertion.
#' @return Static JavaScript handler for the application head.
#' @noRd
workflowPeptideSampleAckScript <- function() {
  "(function install() {\n  if (window.MultiScholaRPerf107Ack) return;\n  if (!window.Shiny || typeof window.Shiny.addCustomMessageHandler !== 'function') {\n    setTimeout(install, 25);\n    return;\n  }\n  window.MultiScholaRPerf107Ack = true;\n  var state = window.MultiScholaRPerf107AckState = {registered: true, received: 0};\n  var timer = null, active = null;\n  window.Shiny.addCustomMessageHandler('multischolar-perf107-ready', function(m) {\n    state.received += 1;\n    if (active && active.request_id === m.request_id && active.content_id === m.content_id) return;\n    if (timer !== null) clearTimeout(timer);\n    active = m;\n    state.request_id = m.request_id;\n    state.content_id = m.content_id;\n    state.expected_text = m.expected_text;\n    state.plot_identity = m.plot_identity;\n    state.sent = false;\n    var started = Date.now(), stable = null, last = null;\n    function check() {\n      var text = document.getElementById(m.text_id);\n      var no_plot = m.plot_id === null || m.plot_id === undefined;\n      var plot = no_plot ? null : document.getElementById(m.plot_id);\n      var img = plot && plot.querySelector('img');\n      state.text_ready = !!text && text.textContent === m.expected_text;\n      state.plot_ready = no_plot || (!!plot && !plot.classList.contains('recalculating'));\n      state.image_ready = no_plot || (!!img && img.complete && img.naturalWidth > 0 &&\n        img.naturalHeight > 0 && img.alt === m.plot_identity);\n      var signature = no_plot ? 'no-plot' : (img && (img.currentSrc || img.src));\n      var ready = state.text_ready && state.plot_ready && state.image_ready;\n      if (!ready || signature !== last) stable = null;\n      last = signature;\n      if (ready && stable === null) stable = Date.now();\n      if (ready && Date.now() - stable >= 200) {\n        state.sent = true;\n        window.Shiny.setInputValue(m.input_id, {\n          request_id: m.request_id, content_id: m.content_id,\n          plot_ready: true, stable_ms: Date.now() - stable\n        }, {priority: 'event'});\n        timer = null;\n        return;\n      }\n      if (Date.now() - started >= 30000) { state.timed_out = true; return; }\n      timer = setTimeout(check, 50);\n    }\n    check();\n  });\n})();"
}

#' @rdname mod_prot_qc_peptide_sample
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr numericInput helpText div actionButton verbatimTextOutput plotOutput
mod_prot_qc_peptide_sample_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Sample Quality",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
          shiny::wellPanel(
          shiny::h4("Minimum Peptides per Sample"),
          shiny::p("Remove samples with fewer than the required distinct protein-group/stripped-peptide identities."),
          shiny::hr(),
          
          shiny::numericInput(ns("min_peptides_per_sample"), 
            "Min Peptides per Sample", 
            value = 500, min = 100, max = 5000, step = 100,
            width = "100%"
          ),
          shiny::helpText("Repeated rows do not add evidence. Higher = stricter (packaged default: 500)."),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_sample_filter"), "Apply Filter", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_sample"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("sample_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("sample_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_peptide_sample
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
runPeptideSampleApplyStep <- function(workflowData,
                                      minPeptidesPerSample,
                                      updateConfigParameterFn = updateConfigParameter,
                                      filterMinNumPeptidesPerSampleFn = filterMinNumPeptidesPerSample,
                                      logInfoFn = logger::log_info,
                                      nowFn = Sys.time) {
  if (usesDIAReferenceStates(workflowData$state_manager)) {
    # KIT-008: unmanaged launches run the same declaration the worker runs, in process.
    return(runDIAReferenceSampleFilter(workflowData, list(peptides_per_sample_cutoff = minPeptidesPerSample)))
  }

  shiny::req(workflowData$state_manager)
  currentS4 <- getDIAExecutionState(workflowData$state_manager, "filterMinNumPeptidesPerSample")
  predecessor <- getDIACurrentReference(workflowData$state_manager)
  shiny::req(currentS4)

  logInfoFn(sprintf(
    "QC Step: Applying sample quality filter (min: %s)",
    minPeptidesPerSample
  ))

  currentS4 <- updateConfigParameterFn(
    theObject = currentS4,
    function_name = "filterMinNumPeptidesPerSample",
    parameter_name = "peptides_per_sample_cutoff",
    new_value = minPeptidesPerSample
  )

  samplesBefore <- currentS4@peptide_data |>
    dplyr::distinct(!!rlang::sym(currentS4@sample_id)) |>
    dplyr::pull(!!rlang::sym(currentS4@sample_id))

  filteredS4 <- filterMinNumPeptidesPerSampleFn(
    theObject = currentS4,
    peptides_per_sample_cutoff = minPeptidesPerSample
  )
  filterSummary <- filteredS4@args$filterMinNumPeptidesPerSample$filter_summary

  samplesAfter <- filteredS4@peptide_data |>
    dplyr::distinct(!!rlang::sym(filteredS4@sample_id)) |>
    dplyr::pull(!!rlang::sym(filteredS4@sample_id))

  samplesRemoved <- setdiff(samplesBefore, samplesAfter)
  samplesRemovedCount <- length(samplesRemoved)

  if (is.null(filteredS4@args$filterMinNumPeptidesPerSample)) {
    filteredS4@args$filterMinNumPeptidesPerSample <- list()
  }
  filteredS4@args$filterMinNumPeptidesPerSample$samples_removed <- samplesRemoved
  filteredS4@args$filterMinNumPeptidesPerSample$samples_removed_count <- samplesRemovedCount
  filteredS4@args$filterMinNumPeptidesPerSample$samples_before_count <- length(samplesBefore)
  filteredS4@args$filterMinNumPeptidesPerSample$samples_after_count <- length(samplesAfter)

  if (is.null(workflowData$qc_params)) {
    workflowData$qc_params <- list()
  }
  if (is.null(workflowData$qc_params$peptide_qc)) {
    workflowData$qc_params$peptide_qc <- list()
  }

  workflowData$qc_params$peptide_qc$sample_filter <- list(
    min_peptides_per_sample = minPeptidesPerSample,
    samples_removed = samplesRemoved,
    samples_removed_count = samplesRemovedCount,
    samples_before_count = length(samplesBefore),
    samples_after_count = length(samplesAfter),
    filter_summary = filterSummary,
    timestamp = nowFn()
  )

  filteredS4 <- .savePeptideQcState(
    state_manager = workflowData$state_manager,
    before = currentS4,
    after = filteredS4,
    stage_id = "sample_filter",
    state_name = "sample_filtered",
    config_object = list(
      min_peptides_per_sample = minPeptidesPerSample,
      samples_removed = samplesRemoved,
      samples_removed_count = samplesRemovedCount,
      filter_summary = filterSummary
    ),
    description = "Applied minimum peptides per sample filter",
    now = nowFn()
  )

  releaseDIAWorkingReference(workflowData$state_manager, predecessor)
  currentS4 <- NULL
  proteinCount <- .countPeptideProteinGroups(filteredS4)
  runCount <- filteredS4@peptide_data |>
    dplyr::distinct(.data[[filteredS4@sample_id]]) |>
    nrow()

  resultText <- paste(
    "Sample Quality Filter Applied Successfully\n",
    "=========================================\n",
    sprintf("Proteins remaining: %d\n", proteinCount),
    sprintf("Samples remaining: %d\n", runCount),
    sprintf("Samples removed: %d\n", samplesRemovedCount),
    sprintf("Min distinct protein-group/stripped-peptide identities per sample: %d\n", minPeptidesPerSample),
    "State saved as: 'sample_filtered'\n"
  )

  if (samplesRemovedCount > 0) {
    resultText <- paste0(
      resultText,
      "\nRemoved samples:\n",
      paste(samplesRemoved, collapse = ", ")
    )
  }

  list(
    filteredS4 = filteredS4,
    resultText = resultText
  )
}

updatePeptideSampleOutputs <- function(output,
                                       samplePlot,
                                       sampleResult,
                                       omicType,
                                       experimentLabel,
                                       renderTextFn = shiny::renderText,
                                       updateProteinFilteringFn = updateProteinFiltering, projectDirs = NULL) {
  output$sample_results <- renderPeptideQcResultText(sampleResult$resultText, renderTextFn)

  plot_arguments <- list(
    data = sampleResult$filteringSummary %||% sampleResult$filteredS4@peptide_data,
    step_name = "6_sample_filtered",
    omic_type = omicType,
    experiment_label = experimentLabel,
    return_grid = TRUE,
    overwrite = TRUE
  )
  if (!is.null(projectDirs)) plot_arguments$project_dirs <- projectDirs
  plotGrid <- if (!is.null(sampleResult$stateReference)) {
    publishDIAQvalueOutput(getDIAActiveManager(), sampleResult, plot_arguments, updateProteinFilteringFn)
  } else do.call(updateProteinFilteringFn, plot_arguments)
  samplePlot(plotGrid)

  invisible(plotGrid)
}

runPeptideSampleApplyObserver <- function(workflowData,
                                          minPeptidesPerSample,
                                          output,
                                          samplePlot,
                                          omicType,
                                          experimentLabel,
                                          runApplyStepFn = runPeptideSampleApplyStep,
                                          updateOutputsFn = updatePeptideSampleOutputs,
                                          showNotificationFn = shiny::showNotification,
                                          removeNotificationFn = shiny::removeNotification,
                                          logInfoFn = logger::log_info,
                                          logErrorFn = logger::log_error, experimentPaths = NULL) {
  withDIAStateScope(getDIAActiveManager(), function() {
  showNotificationFn(
    "Applying sample quality filter...",
    id = "sample_working",
    duration = NULL
  )

  tryCatch({
    sampleResult <- runApplyStepFn(
      workflowData = workflowData,
      minPeptidesPerSample = minPeptidesPerSample
    )

    output_arguments <- list(
      output = output,
      samplePlot = samplePlot,
      sampleResult = sampleResult,
      omicType = omicType,
      experimentLabel = experimentLabel
    )
    if (!is.null(experimentPaths)) {
      output_arguments$projectDirs <- stats::setNames(
        list(experimentPaths), paste(omicType, experimentLabel, sep = "_")
      )
    }
    plotGrid <- do.call(updateOutputsFn, output_arguments)

    logInfoFn("Sample quality filter applied successfully")
    removeNotificationFn("sample_working")
    showNotificationFn("Sample quality filter applied successfully", type = "message")

    list(
      status = "success",
      sampleResult = sampleResult,
      plotGrid = plotGrid
    )
  }, error = function(e) {
    errorMessage <- paste("Error applying sample quality filter:", e$message)
    logErrorFn(errorMessage)
    showNotificationFn(errorMessage, type = "error", duration = 15)
    removeNotificationFn("sample_working")

    list(
      status = "error",
      errorMessage = errorMessage
    )
  })
  })
}

runPeptideSampleRevertStep <- function(workflowData,
                                       logInfoFn = logger::log_info) {
  history <- workflowData$state_manager$getHistory()
  if (length(history) <= 1) {
    stop("No previous state to revert to.")
  }

  previousState <- history[length(history) - 1]
  revertedS4 <- selectDIAStateReference(workflowData$state_manager, previousState)
  logInfoFn(paste("Reverted sample filter to", previousState))

  list(
    previousState = previousState,
    revertedS4 = revertedS4,
    resultText = paste("Reverted to previous state:", previousState)
  )
}

runPeptideSampleRevertObserver <- function(workflowData,
                                           output,
                                           runRevertStepFn = runPeptideSampleRevertStep,
                                           renderTextFn = shiny::renderText,
                                           showNotificationFn = shiny::showNotification,
                                           logErrorFn = logger::log_error) {
  withDIAStateScope(getDIAActiveManager(), function() {
  tryCatch({
    revertResult <- runRevertStepFn(workflowData = workflowData)
    output$sample_results <- renderPeptideQcResultText(revertResult$resultText, renderTextFn)
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

mod_prot_qc_peptide_sample_server <- function(id, workflow_data, omic_type, experiment_label, experiment_paths = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    
    sample_plot <- shiny::reactiveVal(NULL)
    sample_plot_identity <- shiny::reactiveVal("Sample quality filtering summary")
    restored_cutoff <- shiny::reactiveVal(NULL)
    shiny::observeEvent(workflow_data$native_session_restored, {
      value <- workflow_data$qc_params$peptide_qc$sample_filter$min_peptides_per_sample
      restored_cutoff(NULL)
      if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
        restored_cutoff(value)
      }
      sample_plot(NULL)
      output$sample_results <- shiny::renderText(paste("Checkpoint restored:",
        workflow_data$state_manager$current_state))
    })

    shiny::observe({
      value <- restored_cutoff()
      # Fresh Open activates this module before its dynamic input is bound.
      # Send the saved value once that input exists, then leave user edits alone.
      if (is.null(value) || is.null(input$min_peptides_per_sample)) return()
      shiny::updateNumericInput(session, "min_peptides_per_sample", value = value)
      restored_cutoff(NULL)
    })

    # Entry hook (KIT-013 item 4): the module runs whatever hook the managed application carries
    # (installed by the launch from its record) and writes the hook's receipt. The module itself
    # knows no command and no qualification; without a hook the input is ignored.
    shiny::observeEvent(input$entry_hook, {
      application <- session$userData$workflow_execution
      hook <- if (inherits(application, "WorkflowExecutionApplication")) application$entry_hook else NULL
      command <- input$entry_hook
      if (!is.function(hook) || !is.character(command) || length(command) != 1L || is.na(command) ||
          !grepl("^[a-z][a-z0-9_]*$", command)) return()
      path <- file.path(application$config$root, paste0("entry-", command, ".rds"))
      if (file.exists(path)) return()
      result <- tryCatch({
        if (!is.null(managed_action())) stop("entry hook cannot overlap GUI completion")
        hook(command, workflow_data, session, experiment_paths)
      }, error = function(error) list(status = "failed", message = conditionMessage(error)))
      writeWorkflowExecutionMessage(path, result)
    })

    # KIT-008: the managed worker when the session carries a managed application, else in process.
    managed_step <- installWorkflowManagedStep(session, input, workflow_data, "minimum_peptides_per_sample",
      publication = list(experiment_label = experiment_label, experiment_paths = experiment_paths, omic_type = omic_type),
      ui = list(text = "sample_results", plot = "sample_plot", ack = "workflow_ack", notification = "sample_working",
        label = "sample quality filter", step_label = "6_sample_filtered", plot_identity = sample_plot_identity),
      deliver = function(result) {
        arguments <- list(output = output, samplePlot = sample_plot, sampleResult = result,
          omicType = omic_type, experimentLabel = experiment_label)
        if (!is.null(experiment_paths)) arguments$projectDirs <- stats::setNames(list(experiment_paths),
          paste(omic_type, experiment_label, sep = "_"))
        do.call(updatePeptideSampleOutputs, arguments)
      })
    managed_action <- managed_step$active

    # Step 5: Sample Quality Filter (chunk 14)
    shiny::observeEvent(input$apply_sample_filter, {
      if (!is.null(workflowSessionApplication(session)) && usesDIAReferenceStates(workflow_data$state_manager)) {
        managed_step$apply(list(peptides_per_sample_cutoff = input$min_peptides_per_sample))
        return(invisible(NULL))
      }
      apply_arguments <- list(
        workflowData = workflow_data,
        minPeptidesPerSample = input$min_peptides_per_sample,
        output = output,
        samplePlot = sample_plot,
        omicType = omic_type,
        experimentLabel = experiment_label
      )
      if (!is.null(experiment_paths)) apply_arguments$experimentPaths <- experiment_paths
      do.call(runPeptideSampleApplyObserver, apply_arguments)
    })

    # Revert Sample Filter
    shiny::observeEvent(input$revert_sample, {
      managed_step$cancel("reverted")
      runPeptideSampleRevertObserver(
        workflowData = workflow_data,
        output = output
      )
    })
    
    # Render sample filter plot
    output$sample_plot <- shiny::renderPlot({
      shiny::req(sample_plot())
      grid::grid.draw(sample_plot())
    }, alt = function() sample_plot_identity())

  })
}
