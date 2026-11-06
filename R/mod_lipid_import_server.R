#' @rdname mod_lipid_import
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderText showNotification removeNotification updateSelectInput outputOptions
#' @importFrom shinyFiles shinyFileChoose parseFilePaths getVolumes
#' @importFrom logger log_info log_error log_warn
mod_lipid_import_server <- function(id, workflow_data, experiment_paths, volumes = NULL, profile = NULL) {
  emitLipidImportModuleServerEntryDiagnostics(
    volumesIsNull = is.null(volumes)
  )

  shiny::moduleServer(id, function(input, output, session) {
    emitLipidImportModuleServerEntryDiagnostics(
      insideModuleServer = TRUE
    )

    # Check if shinyFiles is available
    use_shiny_files <- probeLipidImportShinyFilesAvailability()

    # Local reactive values
    local_data <- initializeLipidImportLocalData()

    # KIT-096: a profile's preview loader and processing at the selection and process seams (the lipid ones otherwise).
    assay1_selected <- if (is.null(profile)) buildLipidImportAssay1SelectedCallback else function(session, localData) {
      buildLipidImportAssay1SelectedCallback(session, localData, runPreviewLoad = function(session, localData) {
        handleLipidImportDataPreviewLoad(session, localData, loadPreview = profile$load_preview)
      })
    }
    process_request <- if (is.null(profile)) handleLipidImportProcessRequest else function(...) {
      handleLipidImportProcessRequest(..., processImport = profile$process_import)
    }

    # Set up shinyFiles if available
    volumes <- setupLipidImportShinyFileInputs(
      useShinyFiles = use_shiny_files,
      input = input,
      output = output,
      session = session,
      localData = local_data,
      volumes = volumes,
      buildAssay1Selected = assay1_selected
    )

    setupLipidImportStandardFileInputs(
      useShinyFiles = use_shiny_files,
      input = input,
      output = output,
      session = session,
      localData = local_data,
      buildAssay1Selected = assay1_selected
    )

    column_selection_reactives <- buildLipidImportColumnSelectionReactives(
      input = input,
      localData = local_data
    )
    get_lipid_id_col <- column_selection_reactives$lipidIdCol
    get_annotation_col <- column_selection_reactives$annotationCol
    get_sample_columns <- column_selection_reactives$sampleColumns

    registerLipidImportModuleOutputs(
      output = output,
      input = input,
      localData = local_data,
      workflowData = workflow_data,
      getLipidIdCol = get_lipid_id_col,
      getSampleColumns = get_sample_columns
    )

    # Process import
    registerLipidImportProcessObserver(
      input = input,
      workflowData = workflow_data,
      localData = local_data,
      getLipidIdCol = get_lipid_id_col,
      getAnnotationCol = get_annotation_col,
      getSampleColumns = get_sample_columns,
      experimentPaths = experiment_paths,
      handleProcessRequest = process_request
    )
  })
}
