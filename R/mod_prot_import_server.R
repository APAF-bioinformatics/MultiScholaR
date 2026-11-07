#' Setup & Import Applet Server
#' 
#' Server logic for data import and initial setup with multi-format support
#' 
#' @param id Module ID
#' @param workflow_data Reactive values object to store workflow data
#' @param experiment_paths List of paths for this experiment
#' 
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req
#' @importFrom shiny renderUI showNotification removeNotification showModal modalDialog removeModal updateNumericInput updateTextInput
#' @importFrom shinyjs html
#' @importFrom logger log_info log_error log_warn
#' @importFrom vroom vroom
#' @importFrom ini read.ini
#' @importFrom stats setNames
#' @importFrom DT renderDT datatable formatStyle
#' @export
#' @param volumes Runtime inputs used by this function; see the usage section for accepted values.
mod_prot_import_server <- function(id, workflow_data, experiment_paths, volumes = NULL) {
  message(sprintf("--- Entering mod_prot_import_server ---"))
  message(sprintf("   mod_prot_import_server Arg: id = %s", id))
  # KIT-102 W6: the shared import frame runs the proteomics panel (protImportPanelServer) in its module session.
  mod_omic_import_server(id, workflow_data, experiment_paths, volumes = volumes, config = omicGuiConfig("proteomics"))
}

#' The proteomics import source panel's server: today's import module body, called by the import frame inside its
#' moduleServer (the frame's panel contract, R/mod_omic_import.R). Success writes the statuses through the frame's
#' `complete` (completeProtImportSuccessState).
#' @param input,output,session The frame's module session.
#' @param workflow_data The workflow bus.
#' @param experiment_paths The experiment's project paths.
#' @param volumes shinyFiles volumes.
#' @param config `omicGuiConfig("proteomics")`.
#' @param complete The frame's completion helper.
#' @noRd
protImportPanelServer <- function(input, output, session, workflow_data, experiment_paths, volumes, config,
                                  complete) {
    message(sprintf("   mod_prot_import_server Step: Inside moduleServer function"))

    local_data <- createProtImportLocalData()
    use_shiny_files <- requireNamespace("shinyFiles", quietly = TRUE) && !is_test_mode()

    if (use_shiny_files) {
      volumes <- resolveProtImportShinyFileVolumes(volumes = volumes)
      registerProtImportShinyFileInputs(
        input = input,
        output = output,
        session = session,
        localData = local_data,
        volumes = volumes
      )
    }
    import_paths <- createProtImportPathResolvers(
      input = input,
      localData = local_data,
      useShinyFiles = use_shiny_files
    )

    registerProtImportObservers(
      input = input,
      output = output,
      session = session,
      workflowData = workflow_data,
      localData = local_data,
      experimentPaths = experiment_paths,
      useShinyFiles = use_shiny_files,
      searchResultsPathReactive = import_paths$searchResultsPath,
      fastaPathReactive = import_paths$fastaPath,
      activeFormatReactive = import_paths$activeFormat,
      startProcessingFn = startProtImportProcessing,
      runProcessingSafelyFn = runProtImportProcessingSafely,
      readDataWithStatusFn = readProtImportDataWithStatus,
      applyWorkflowWithStatusFn = applyProtImportWorkflowWithStatus,
      runMixedSpeciesAnalysisIfNeededFn = runProtImportMixedSpeciesAnalysisIfNeeded,
      loadConfigurationResourcesFn = loadProtImportConfigurationResources,
      recordSetupLogFn = recordProtImportSetupLog,
      finalizeProcessingFn = finalizeProtImportProcessing,
      confirmOrganismSelectionFn = confirmProtImportOrganismSelection,
      updateStatusFn = updateProtImportProcessingStatus,
      importDataByFormatFn = importProtImportDataByFormat,
      recordImportedDataFn = recordProtImportImportedData,
      importDiannFn = importDIANNData,
      importSpectronautFn = importSpectronautData,
      importFragpipeFn = importFragPipeData,
      importMaxquantFn = importMaxQuantData,
      importPdTmtFn = importProteomeDiscovererTMTData,
      applyImportResultFn = applyProtImportResultToWorkflow,
      processFastaDataFn = processProtImportFastaData,
      sanitizeRunNamesFn = sanitizeProtImportRunNames,
      processFastaFn = processFastaFile,
      analyzeMixedSpeciesDataFn = analyzeProtImportMixedSpeciesData,
      extractOrganismsFn = extractOrganismsFromFasta,
      analyzeDistributionFn = analyzeOrganismDistribution,
      buildSelectionModalFn = buildProtImportOrganismSelectionModal,
      loadConfigurationFn = loadProtImportConfiguration,
      storeConfigurationFn = storeProtImportConfiguration,
      loadOptionalMappingsFn = loadProtImportOptionalMappings,
      readConfigFn = readProtGuiConfig,
      getDefaultConfigFn = getDefaultProteomicsConfig,
      readOptionalMappingFn = readProtImportOptionalMapping,
      fileExistsFn = file.exists,
      fileCopyFn = file.copy,
      downloadFileFn = download.file,
      systemFileFn = system.file,
      resolveFilenameFn = resolveProtImportInputFilename,
      finalizeSetupStateFn = finalizeProtImportSetupState,
      completeSuccessStateFn = function(...) completeProtImportSuccessState(..., complete = complete),
      logInfo = log_info,
      logError = log_error,
      logWarn = log_warn,
      messageFn = message
    )

    registerProtImportOutputs(
      output = output,
      workflowData = workflow_data,
      localData = local_data,
      session = session,
      activeFormatReactive = import_paths$activeFormat
    )
}
