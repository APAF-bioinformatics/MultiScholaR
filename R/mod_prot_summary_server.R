#' @rdname sessionSummaryModule
#' @export
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req renderText showNotification downloadHandler withProgress incProgress
#' @importFrom logger log_info log_error
#' @importFrom utils write.table
#' @importFrom purrr detect
#' @importFrom readr read_tsv
mod_prot_summary_server <- function(id, project_dirs, omic_type = "proteomics", experiment_label = NULL, workflow_data = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Auto-populate experiment label if provided
    if (!is.null(experiment_label)) {
      shiny::updateTextInput(session, "experiment_label", value = experiment_label)
    }
    
    # Reactive values for tracking state
    values <- shiny::reactiveValues(
      workflow_args_saved = FALSE,
      files_copied = FALSE,
      report_generated = FALSE,
      report_path = NULL
    )
    
    registerProtSummaryTemplateStatusOutput(
      output = output,
      projectDirs = project_dirs,
      omicType = omic_type,
      renderTextFn = shiny::renderText,
      reqFn = shiny::req,
      buildTemplateStatusFn = buildProtSummaryTemplateStatus
    )
    
    observeProtSummaryWorkflowArgsSave(
      input = input,
      output = output,
      values = values,
      projectDirs = project_dirs,
      workflowData = workflow_data,
      omicType = omic_type,
      completeWorkflowArgsSaveFn = completeProtSummaryWorkflowArgsSave,
      resolveFinalS4StateFn = resolveProtSummaryFinalS4State,
      assignFn = assign,
      createWorkflowArgsFn = createWorkflowArgsFromConfig,
      saveRDSFn = saveRDS,
      renderTextFn = shiny::renderText,
      showNotificationFn = shiny::showNotification,
      dirExistsFn = dir.exists,
      dirCreateFn = dir.create,
      fileExistsFn = file.exists,
      writeLinesFn = writeLines,
      timestampFn = Sys.time,
      catFn = cat
    )
    
    observeProtSummaryPublicationCopy(
      input = input,
      output = output,
      values = values,
      projectDirs = project_dirs,
      workflowData = workflow_data,
      omicType = omic_type,
      fallbackBootstrapFn = bootstrapProtSummaryCopyFallbackStudyParams,
      prepareCopyInputsFn = prepareProtSummaryCopyInputs,
      runPublicationCopyFn = runProtSummaryPublicationCopy,
      handleCopyErrorFn = handleProtSummaryPublicationCopyError,
      withProgressFn = shiny::withProgress,
      catFn = cat
    )
    
    observeProtSummaryReportGeneration(
      input = input,
      output = output,
      values = values,
      projectDirs = project_dirs,
      workflowData = workflow_data,
      omicType = omic_type,
      validateProjectDirsFn = validateProtSummaryProjectDirs,
      runReportProgressFn = runProtSummaryReportProgress,
      withProgressFn = shiny::withProgress
    )
    
    observeProtSummaryGithubPush(
      input = input,
      output = output,
      values = values,
      projectDirs = project_dirs,
      omicType = omic_type,
      completeGithubPushFn = completeProtSummaryGithubPush,
      withProgressFn = shiny::withProgress
    )
    
    observeProtSummarySessionStateExport(
      input = input,
      values = values,
      projectDirs = project_dirs,
      omicType = omic_type,
      completeSessionStateExportFn = completeProtSummarySessionStateExport
    )
    
    initializeProtSummaryDefaultOutputs(output = output)
    
  })
}

