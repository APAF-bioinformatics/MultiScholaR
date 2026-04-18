#' @rdname mod_metab_import
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderText showNotification removeNotification updateSelectInput outputOptions
#' @importFrom shinyFiles shinyFileChoose parseFilePaths getVolumes
#' @importFrom logger log_info log_error log_warn
mod_metab_import_server <- function(id, workflow_data, experiment_paths, volumes = NULL) {
  message("--- Entering mod_metab_import_server ---")
  message(sprintf("   mod_metab_import_server: volumes param is NULL = %s", is.null(volumes)))

  shiny::moduleServer(id, function(input, output, session) {
    message("   mod_metab_import_server: Inside moduleServer function")

    # Check if shinyFiles is available
    use_shiny_files <- requireNamespace("shinyFiles", quietly = TRUE)
    message(sprintf("   mod_metab_import_server: shinyFiles available = %s", use_shiny_files))

    # Local reactive values
    local_data <- shiny::reactiveValues(
      assay1_file = NULL,
      assay1_data = NULL,
      assay1_import_result = NULL,
      assay2_file = NULL,
      assay2_data = NULL,
      assay2_import_result = NULL,
      detected_format = NULL,
      format_confidence = NULL,
      all_headers = NULL
    )

    importData <- setupMetabImportAssaySelectionCallback(
      localData = local_data,
      session = session
    )

    # Set up shinyFiles if available
    if (use_shiny_files) {
      volumes <- setupMetabImportShinyFiles(
        input = input,
        output = output,
        session = session,
        volumes = volumes,
        localData = local_data,
        importDataFn = importData
      )
    }

    columnAccessors <- setupMetabImportColumnAccessors(
      input = input,
      localData = local_data
    )

    setupMetabImportFileLoadedOutput(
      output = output,
      localData = local_data
    )

    setupMetabImportFormatDetectionStatusOutput(
      output = output,
      localData = local_data
    )

    setupMetabImportMetaboliteIdStatusOutput(
      output = output,
      input = input,
      localData = local_data
    )

    setupMetabImportAnnotationStatusOutput(
      output = output,
      input = input,
      localData = local_data
    )

    setupMetabImportSampleColumnsDisplayOutput(
      output = output,
      localData = local_data
    )

    setupMetabImportAvailableColumnsDisplayOutput(
      output = output,
      localData = local_data
    )

    setupMetabImportCustomMetaboliteIdStatusOutput(
      output = output,
      input = input,
      localData = local_data
    )

    setupMetabImportCustomAnnotationStatusOutput(
      output = output,
      input = input,
      localData = local_data
    )

    setupMetabImportValidationSummaryOutput(
      output = output,
      localData = local_data,
      columnAccessors = columnAccessors
    )

    setupMetabImportProcessingObserver(
      input = input,
      localData = local_data,
      columnAccessors = columnAccessors,
      workflowData = workflow_data
    )

    setupMetabImportStatusOutput(
      output = output,
      workflowData = workflow_data
    )
  })
}

