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
  "(function install() {\n  if (window.MultiScholaRPerf107Ack) return;\n  if (!window.Shiny || typeof window.Shiny.addCustomMessageHandler !== 'function') {\n    setTimeout(install, 25);\n    return;\n  }\n  window.MultiScholaRPerf107Ack = true;\n  var state = window.MultiScholaRPerf107AckState = {registered: true, received: 0};\n  var timer = null, active = null;\n  window.Shiny.addCustomMessageHandler('multischolar-perf107-ready', function(m) {\n    state.received += 1;\n    if (active && active.request_id === m.request_id && active.content_id === m.content_id) return;\n    if (timer !== null) clearTimeout(timer);\n    active = m;\n    state.request_id = m.request_id;\n    state.content_id = m.content_id;\n    state.expected_text = m.expected_text;\n    state.plot_identity = m.plot_identity;\n    state.sent = false;\n    var started = Date.now(), stable = null, last = null;\n    function check() {\n      var text = document.getElementById(m.text_id);\n      var plot = document.getElementById(m.plot_id);\n      var img = plot && plot.querySelector('img');\n      state.text_ready = !!text && text.textContent === m.expected_text;\n      state.plot_ready = !!plot && !plot.classList.contains('recalculating');\n      state.image_ready = !!img && img.complete && img.naturalWidth > 0 &&\n        img.naturalHeight > 0 && img.alt === m.plot_identity;\n      var signature = img && (img.currentSrc || img.src);\n      var ready = state.text_ready && state.plot_ready && state.image_ready;\n      if (!ready || signature !== last) stable = null;\n      last = signature;\n      if (ready && stable === null) stable = Date.now();\n      if (ready && Date.now() - stable >= 200) {\n        state.sent = true;\n        window.Shiny.setInputValue(m.input_id, {\n          request_id: m.request_id, content_id: m.content_id,\n          plot_ready: true, stable_ms: Date.now() - stable\n        }, {priority: 'event'});\n        timer = null;\n        return;\n      }\n      if (Date.now() - started >= 30000) { state.timed_out = true; return; }\n      timer = setTimeout(check, 50);\n    }\n    check();\n  });\n})();"
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

#' Resolve the one server-owned PERF-107 GUI route.
#' @param session Current Shiny session.
#' @return Frozen qualification binding or NULL for ordinary GUI behavior.
#' @noRd
workflowPeptideSampleQualification <- function(session) {
  application <- session$userData$workflow_execution
  if (!inherits(application, "WorkflowExecutionApplication")) return(NULL)
  qualification <- application$config$qualification
  if (!is.list(qualification) ||
      !identical(qualification$programme, "perf107-phase-a")) return(NULL)
  workflowQualificationVariant(application$config, "minimum_peptides_per_sample")
  qualification
}

#' Prepare or inspect only the explicitly bound small GUI entry fixture.
#' @param command Fixed setup or post-completion inspection command.
#' @param workflow_data,session,experiment_paths Existing GUI ownership and paths.
#' @return Metadata pointing at a saved verification store; no assay handoff.
#' @noRd
runWorkflowGUIEntryFixture <- function(command, workflow_data, session, experiment_paths) {
  qualification <- workflowPeptideSampleQualification(session)
  application <- session$userData$workflow_execution
  if (is.null(qualification) || !isTRUE(application$config$entry_fixture_test) ||
      !identical(as.numeric(qualification$deadline), 900) ||
      !is.character(command) || length(command) != 1L || is.na(command) ||
      !command %in% c("prepare_reference", "capture_sample_filtered", "capture_imputed")) {
    stop("entry fixture is not authorized")
  }
  checkWorkflowApplication(application)
  if (!is.null(application$active)) stop("entry fixture cannot overlap an action")
  manager <- workflow_data$state_manager
  expected_state <- switch(command, prepare_reference = "intensity_filtered",
    capture_sample_filtered = "sample_filtered", capture_imputed = "imputed")
  if (!identical(manager$current_state, expected_state)) stop("entry fixture state differs")
  rows <- if (usesDIAReferenceStates(manager)) {
    readDIAStateMetadata(manager)$data$dimensions$rows
  } else nrow(manager$getState()@peptide_data)
  if (!is.numeric(rows) || length(rows) != 1L || !is.finite(rows) || rows > 2000 || rows < 0) {
    stop("entry fixture exceeds the small verification boundary")
  }
  if (command == "prepare_reference") {
    if (!identical(qualification$entry, "reference") || usesDIAReferenceStates(manager)) {
      stop("entry fixture already prepared or wrong route")
    }
    manager <- stageDIAGUIReferenceOwner(manager, experiment_paths)
    manager <- adoptDIAGUIReferenceOwner(manager, workflow_data, session)
  } else if (!usesDIAReferenceStates(manager)) {
    # Snapshot only AFTER GUI completion; the native action retains its original
    # owner and scientific call. The existing codec strips worker handles.
    manager <- stageDIAGUIReferenceOwner(manager, experiment_paths)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
  }
  list(status = "complete", command = command,
    root = dirname(manager$.context$store$directory), current_state = manager$current_state,
    rows = rows, qualification = qualification,
    measured_action_active = !is.null(application$active), verification_only = TRUE,
    ownership = diaStateOwnership(manager))
}

#' Recover the bounded display product for one accepted managed action.
#' @param manager Reference-backed owner.
#' @param result Accepted worker/controller result.
#' @param workflowData Existing GUI workflow bus.
#' @return Existing sample-filter display handoff without assay hydration.
#' @noRd
readManagedPeptideSampleResult <- function(manager, result, workflowData) {
  record <- readDIAStateMetadata(manager, result$state_reference)
  product <- record$data$products$qc_progress
  envelope <- hydrateSavedWorkflowEntity(manager$.context$store, product)
  checkSavedWorkflowConsumerPayload("scientific_product", envelope)
  value <- c(list(stateReference = result$state_reference), envelope$value)
  workflowData$qc_params <- value$qcParams
  value
}

#' Invoke the ordinary native GUI helper for an explicit comparison arm.
#' @noRd
runWorkflowNativeSampleControl <- function(workflow_data, cutoff) {
  manager <- workflow_data$state_manager
  if (!usesDIAReferenceStates(manager)) {
    return(runPeptideSampleApplyStep(workflow_data, cutoff))
  }
  # Only the control deliberately hydrates its selected input. The facade keeps
  # public native compute/plot ownership and delegates durable history to its owner.
  native_workflow <- new.env(parent = emptyenv())
  native_workflow$qc_params <- workflow_data$qc_params
  native_workflow$state_manager <- list(
    current_state = manager$current_state, audit_enabled = manager$audit_enabled,
    getState = manager$getState, saveState = manager$saveState
  )
  result <- runPeptideSampleApplyStep(native_workflow, cutoff)
  workflow_data$qc_params <- native_workflow$qc_params
  result
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
    return(runDIAReferencePeptideGroups(workflowData, "filterMinNumPeptidesPerSample",
      list(peptides_per_sample_cutoff = minPeptidesPerSample), "sample_filtered", "Applied minimum peptides per sample filter"))
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
    managed_action <- shiny::reactiveVal(NULL)
    pending_publication <- shiny::reactiveVal(NULL)
    restored_cutoff <- shiny::reactiveVal(NULL)
    shiny::observeEvent(workflow_data$native_session_restored, {
      value <- workflow_data$qc_params$peptide_qc$sample_filter$min_peptides_per_sample
      restored_cutoff(NULL)
      if (is.numeric(value) && length(value) == 1L && is.finite(value)) {
        restored_cutoff(value)
      }
      pending_publication(NULL)
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

    # Private fixture setup/inspection, admitted only by a hashed server member.
    # No scientific GUI action is replaced by this observer.
    shiny::observeEvent(input$perf107_fixture, {
      qualification <- workflowPeptideSampleQualification(session)
      application <- session$userData$workflow_execution
      if (is.null(qualification) || !isTRUE(application$config$entry_fixture_test)) return()
      command <- input$perf107_fixture
      if (!is.character(command) || length(command) != 1L || is.na(command) ||
          !command %in% c("prepare_reference", "capture_sample_filtered", "capture_imputed")) return()
      path <- file.path(application$config$root, paste0("entry-", command, ".rds"))
      if (file.exists(path)) return()
      result <- tryCatch({
        if (!is.null(managed_action())) stop("entry fixture cannot overlap GUI completion")
        runWorkflowGUIEntryFixture(command, workflow_data, session, experiment_paths)
      }, error = function(error) list(status = "failed", message = conditionMessage(error)))
      writeWorkflowExecutionMessage(path, result)
    })

    fail_managed_action <- function(action, error, status = "error") {
      if (identical(status, "publication_pending")) pending_publication(action)
      if (!is.null(action$owner) && !is.null(action$owner$job)) {
        cancelWorkflowAction(action$owner, status)
      }
      if (!action$reservation$closed) {
        finishWorkflowGUIAction(action$reservation, status)
      }
      managed_action(NULL)
      shiny::removeNotification("sample_working")
      message <- paste("Error applying sample quality filter:", conditionMessage(error))
      logger::log_error(message)
      shiny::showNotification(message, type = "error", duration = 15)
      invisible(NULL)
    }

    deliver_managed_action <- function(action, sample_result) {
      plot_grid <- withDIAStateScope(workflow_data$state_manager, function() {
        arguments <- list(output = output, samplePlot = sample_plot,
          sampleResult = sample_result, omicType = omic_type,
          experimentLabel = experiment_label)
        if (!is.null(experiment_paths)) {
          arguments$projectDirs <- stats::setNames(list(experiment_paths),
            paste(omic_type, experiment_label, sep = "_"))
        }
        do.call(updatePeptideSampleOutputs, arguments)
      })
      content_id <- digest::digest(list(sample_result$resultText,
        sample_result$stateReference$entity_id %||% workflow_data$state_manager$current_state,
        "6_sample_filtered"), algo = "sha256")
      action$sample_result <- sample_result
      action$plot_grid <- plot_grid
      action$content_id <- content_id
      action$accepted_manager <- workflow_data$state_manager
      action$accepted_reference <- getDIAStateReference(workflow_data$state_manager)
      # Match renderText's output conversion exactly, including its treatment
      # of a final newline when no removed-sample list follows the summary.
      action$expected_text <- paste(utils::capture.output(cat(sample_result$resultText)), collapse = "\n")
      sample_plot_identity(content_id)
      action$ack_started <- proc.time()[["elapsed"]]
      action$phase <- "awaiting_browser"
      managed_action(action)
      session$onFlushed(function() {
        session$sendCustomMessage("multischolar-perf107-ready", list(
          request_id = action$reservation$id, content_id = content_id,
          expected_text = action$expected_text, plot_identity = content_id,
          input_id = session$ns("workflow_ack"),
          text_id = session$ns("sample_results"),
          plot_id = session$ns("sample_plot")
        ))
      }, once = TRUE)
      invisible(NULL)
    }
    
    # Step 5: Sample Quality Filter (chunk 14)
    shiny::observeEvent(input$apply_sample_filter, {
      qualification <- workflowPeptideSampleQualification(session)
      if (!is.null(qualification)) {
        if (!is.null(managed_action())) {
          shiny::showNotification("A sample-filter action is already running", type = "warning")
          return(invisible(NULL))
        }
        application <- session$userData$workflow_execution
        publication <- list(experiment_label = experiment_label,
          experiment_paths = experiment_paths, omic_type = omic_type)
        action <- list(owner = NULL, phase = "starting")
        action$reservation <- tryCatch(beginWorkflowGUIAction(application,
          workflow_data$state_manager, "minimum_peptides_per_sample", publication),
          error = function(error) error)
        if (inherits(action$reservation, "error")) {
          shiny::showNotification(conditionMessage(action$reservation),
            type = "error", duration = 15)
          return(invisible(NULL))
        }
        shiny::showNotification("Applying sample quality filter...",
          id = "sample_working", duration = NULL)
        tryCatch({
          manager <- workflow_data$state_manager
          pending <- pending_publication()
          if (!is.null(pending)) {
            if (!identical(pending$owner$manager, manager) ||
                !identical(pending$cutoff, input$min_peptides_per_sample)) {
              stop("retry pending publication with its original parameters first")
            }
            action$owner <- pending$owner
            action$worker_request_id <- pending$worker_request_id
            action$cutoff <- pending$cutoff
            accepted <- resumeWorkflowAction(action$owner, action$worker_request_id)
            display <- readManagedPeptideSampleResult(manager, accepted, workflow_data)
            deliver_managed_action(action, display)
            pending_publication(NULL)
            return(invisible(NULL))
          }
          action$cutoff <- input$min_peptides_per_sample
          is_reference <- inherits(manager, "DIAStateAdapter")
          if (!identical(qualification$entry,
              if (is_reference) "reference" else "resident") &&
              !(identical(qualification$entry, "resident") && is_reference &&
                isTRUE(manager$.context$gui_pilot_adopted))) {
            stop("GUI entry ownership differs from the frozen member")
          }
          if (identical(qualification$variant, "native")) {
            result <- runWorkflowNativeSampleControl(workflow_data,
              input$min_peptides_per_sample)
            deliver_managed_action(action, result)
          } else {
            if (!is_reference) {
              manager <- stageDIAGUIReferenceOwner(manager, experiment_paths)
              manager <- adoptDIAGUIReferenceOwner(manager, workflow_data, session)
              action$reservation$manager <- manager
            }
            owner <- manager$.context$store$execution_owner
            if (is.null(owner)) owner <- newWorkflowExecutionController(manager, application)
            action$owner <- owner
            action$phase <- "running"
            action$worker_request_id <- runWorkflowAction(
              "minimum_peptides_per_sample",
              list(state = getDIAStateReference(manager)),
              list(peptides_per_sample_cutoff = input$min_peptides_per_sample),
              owner, reservation = action$reservation
            )
            managed_action(action)
          }
        }, error = function(error) fail_managed_action(action, error))
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

    shiny::observe({
      action <- managed_action()
      if (is.null(action) || !identical(action$phase, "running")) return()
      if (!identical(workflow_data$state_manager, action$owner$manager)) {
        fail_managed_action(action, simpleError("project was replaced"), "project_replaced")
        return()
      }
      shiny::invalidateLater(50, session)
      result <- tryCatch(pollWorkflowAction(action$owner), error = function(error) error)
      if (inherits(result, "error")) {
        fail_managed_action(action, result)
        return()
      }
      if (result$status %in% c("admitting", "running")) return()
      if (!identical(result$status, "complete")) {
        fail_managed_action(action, simpleError(paste("managed action ended:",
          result$status)), result$status)
        return()
      }
      tryCatch({
        accepted <- resumeWorkflowAction(action$owner, action$worker_request_id)
        display <- readManagedPeptideSampleResult(workflow_data$state_manager,
          accepted, workflow_data)
        deliver_managed_action(action, display)
      }, error = function(error) fail_managed_action(action, error))
    })

    shiny::observe({
      action <- managed_action()
      if (is.null(action) || !identical(action$phase, "awaiting_browser")) return()
      if (!identical(workflow_data$state_manager, action$accepted_manager) ||
          !identical(getDIAStateReference(workflow_data$state_manager), action$accepted_reference)) {
        fail_managed_action(action, simpleError("result was superseded"), "stale")
        return()
      }
      elapsed <- proc.time()[["elapsed"]] - action$ack_started
      if (elapsed >= 30) {
        fail_managed_action(action,
          simpleError("browser acknowledgement timed out"),
          "browser_ack_timeout")
        return()
      }
      shiny::invalidateLater(1000, session)
      session$sendCustomMessage("multischolar-perf107-ready", list(
        request_id = action$reservation$id, content_id = action$content_id,
        expected_text = action$expected_text, plot_identity = action$content_id,
        input_id = session$ns("workflow_ack"),
        text_id = session$ns("sample_results"),
        plot_id = session$ns("sample_plot")
      ))
    })

    shiny::observeEvent(input$workflow_ack, {
      action <- managed_action()
      acknowledgement <- input$workflow_ack
      if (is.null(action) || !identical(action$phase, "awaiting_browser") ||
          !is.list(acknowledgement) ||
          !identical(acknowledgement$request_id, action$reservation$id) ||
          !identical(acknowledgement$content_id, action$content_id) ||
          !isTRUE(acknowledgement$plot_ready) ||
          !is.numeric(acknowledgement$stable_ms) ||
          length(acknowledgement$stable_ms) != 1L ||
          !is.finite(acknowledgement$stable_ms) ||
          acknowledgement$stable_ms < 200) return()
      if (!identical(workflow_data$state_manager, action$accepted_manager) ||
          !identical(getDIAStateReference(workflow_data$state_manager),
            action$accepted_reference)) {
        fail_managed_action(action, simpleError("result was superseded"), "stale")
        return()
      }
      receipt <- tryCatch(finishWorkflowGUIAction(action$reservation, "complete",
        list(request_id = acknowledgement$request_id,
          content_id = acknowledgement$content_id, plot_ready = TRUE,
          stable_seconds = acknowledgement$stable_ms / 1000)),
        error = function(error) error)
      if (inherits(receipt, "error")) {
        fail_managed_action(action, receipt)
        return()
      }
      managed_action(NULL)
      shiny::removeNotification("sample_working")
      logger::log_info("Sample quality filter applied successfully")
      shiny::showNotification("Sample quality filter applied successfully", type = "message")
    })
    
    # Revert Sample Filter
    shiny::observeEvent(input$revert_sample, {
      pending_publication(NULL)
      action <- managed_action()
      if (!is.null(action)) {
        if (!is.null(action$owner) && !is.null(action$owner$job)) {
          cancelWorkflowAction(action$owner, "reverted")
        }
        if (!action$reservation$closed) {
          finishWorkflowGUIAction(action$reservation, "reverted")
        }
        managed_action(NULL)
        shiny::removeNotification("sample_working")
      }
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

    if (is.function(session$onSessionEnded)) {
      session$onSessionEnded(function() {
        action <- shiny::isolate(managed_action())
        if (is.null(action)) return(invisible(NULL))
        if (!is.null(action$owner) && !is.null(action$owner$job)) {
          cancelWorkflowAction(action$owner, "session_ended")
        }
        if (!action$reservation$closed) {
          finishWorkflowGUIAction(action$reservation, "session_ended")
        }
        invisible(NULL)
      })
    }
  })
}
