createProtImportLocalData <- function(
    reactiveValuesFn = shiny::reactiveValues
) {
  reactiveValuesFn(
    processing = FALSE,
    detected_format = NULL,
    format_confidence = NULL,
    search_results_file = NULL,
    fasta_file_path = NULL,
    uniprot_mapping_file = NULL,
    uniparc_mapping_file = NULL,
    config_file_path = NULL,
    organism_distribution = NULL,
    organism_mapping = NULL,
    pending_import_completion = FALSE,
    waiting_for_organism_selection = FALSE
  )
}

registerProtImportShinyFileInputs <- function(
    input,
    output,
    session,
    localData,
    volumes,
    registerChooserFn = registerProtImportShinyFileChooser,
    handleSelectionFn = handleProtImportShinyFileSelection,
    observeEventFn = shiny::observeEvent
) {
  specs <- list(
    list(
      inputId = "search_results",
      filetypes = c("tsv", "txt", "tab", "csv", "xlsx", "zip", "parquet"),
      localField = "search_results_file",
      outputId = "search_results_path",
      catchErrors = TRUE
    ),
    list(
      inputId = "fasta_file",
      filetypes = c("fasta", "fa", "faa"),
      localField = "fasta_file_path",
      outputId = "fasta_file_path",
      catchErrors = FALSE
    ),
    list(
      inputId = "uniprot_mapping",
      filetypes = c("tsv", "txt", "tab"),
      localField = "uniprot_mapping_file",
      outputId = "uniprot_mapping_path",
      catchErrors = FALSE
    ),
    list(
      inputId = "uniparc_mapping",
      filetypes = c("tsv", "txt", "tab"),
      localField = "uniparc_mapping_file",
      outputId = "uniparc_mapping_path",
      catchErrors = FALSE
    ),
    list(
      inputId = "config_file",
      filetypes = c("ini", "cfg"),
      localField = "config_file_path",
      outputId = "config_file_path",
      catchErrors = FALSE
    )
  )

  for (spec in specs) {
    registerChooserFn(
      input = input,
      inputId = spec$inputId,
      volumes = volumes,
      session = session,
      filetypes = spec$filetypes
    )

    local({
      current_spec <- spec
      observeEventFn(input[[current_spec$inputId]], {
        handleSelectionFn(
          selectedInput = input[[current_spec$inputId]],
          volumes = volumes,
          localData = localData,
          localField = current_spec$localField,
          output = output,
          outputId = current_spec$outputId,
          catchErrors = current_spec$catchErrors
        )
      })
    })
  }

  invisible(volumes)
}

createProtImportPathResolvers <- function(
    input,
    localData,
    useShinyFiles,
    reactiveFn = shiny::reactive,
    resolvePrimaryUploadPathFn = resolveProtImportPrimaryUploadPath,
    resolveActiveFormatFn = resolveActiveProtImportFormat
) {
  list(
    searchResultsPath = reactiveFn({
      resolvePrimaryUploadPathFn(
        useShinyFiles = useShinyFiles,
        shinyPath = localData$search_results_file,
        standardInput = input$search_results_standard
      )
    }),
    fastaPath = reactiveFn({
      resolvePrimaryUploadPathFn(
        useShinyFiles = useShinyFiles,
        shinyPath = localData$fasta_file_path,
        standardInput = input$fasta_file_standard
      )
    }),
    activeFormat = reactiveFn({
      resolveActiveFormatFn(
        formatOverride = input$format_override,
        detectedFormat = localData$detected_format
      )
    })
  )
}

registerProtImportOutputs <- function(
    output,
    workflowData,
    localData,
    session,
    activeFormatReactive,
    renderUiFn = shiny::renderUI,
    renderDtFn = DT::renderDT,
    buildFormatDetectionUiFn = buildProtImportFormatDetectionUi,
    buildFormatSpecificOptionsUiFn = buildProtImportFormatSpecificOptionsUi,
    buildOrganismDistributionTableFn = buildProtImportOrganismDistributionTable,
    buildStatusUiFn = buildProtImportStatusUi,
    buildDataSummaryUiFn = buildProtImportDataSummaryUi,
    reqFn = shiny::req
) {
  output$format_detection <- renderUiFn({
    buildFormatDetectionUiFn(
      detectedFormat = localData$detected_format,
      formatConfidence = localData$format_confidence
    )
  })

  output$format_specific_options <- renderUiFn({
    buildFormatSpecificOptionsUiFn(
      format = activeFormatReactive(),
      ns = session$ns
    )
  })

  output$organism_dist_table <- renderDtFn({
    reqFn(localData$organism_distribution)
    buildOrganismDistributionTableFn(localData$organism_distribution)
  })

  output$import_status <- renderUiFn({
    buildStatusUiFn(
      isProcessing = localData$processing,
      dataTbl = workflowData$data_tbl,
      dataFormat = workflowData$data_format,
      dataType = workflowData$data_type
    )
  })

  output$data_summary <- renderUiFn({
    buildDataSummaryUiFn(
      dataTbl = workflowData$data_tbl,
      setupImportLog = workflowData$processing_log$setup_import
    )
  })

  invisible(output)
}

registerProtImportObservers <- function(
    input,
    output,
    session,
    workflowData,
    localData,
    experimentPaths,
    useShinyFiles,
    searchResultsPathReactive,
    fastaPathReactive,
    activeFormatReactive,
    updateCheckpointCaptureFn = updateProtImportCheckpointCaptureOption,
    toggleMixedSpeciesInputsFn = toggleProtImportMixedSpeciesInputs,
    runFormatDetectionFn = runProtImportFormatDetection,
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
    resolveFilenameFn = resolveProtImportInputFilename,
    finalizeSetupStateFn = finalizeProtImportSetupState,
    completeSuccessStateFn = completeProtImportSuccessState,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    showNotificationFn = shiny::showNotification,
    showModalFn = shiny::showModal,
    readConfigFn = readConfigFile,
    getDefaultConfigFn = getDefaultProteomicsConfig,
    readOptionalMappingFn = readProtImportOptionalMapping,
    fileExistsFn = file.exists,
    fileCopyFn = file.copy,
    downloadFileFn = download.file,
    systemFileFn = system.file,
    assignFn = assign,
    assignEnv = .GlobalEnv,
    logInfo = logger::log_info,
    logError = logger::log_error,
    logWarn = logger::log_warn,
    messageFn = message
) {
  observeEventFn(input$capture_checkpoints, {
    updateCheckpointCaptureFn(
      captureCheckpoints = input$capture_checkpoints,
      optionsFn = options,
      logInfo = logInfo
    )
  }, ignoreInit = FALSE)

  observeEventFn(input$mixed_species_fasta, {
    toggleMixedSpeciesInputsFn(
      mixedSpeciesFasta = input$mixed_species_fasta,
      disableInput = shinyjs::disable,
      enableInput = shinyjs::enable,
      messageFn = messageFn
    )
  }, ignoreInit = FALSE)

  observeEventFn(searchResultsPathReactive(), {
    file_path <- searchResultsPathReactive()
    reqFn(file_path)

    runFormatDetectionFn(
      filePath = file_path,
      useShinyFiles = useShinyFiles,
      localData = localData,
      searchResultsStandard = input$search_results_standard,
      logInfo = logInfo,
      logError = logError
    )
  })

  observeEventFn(input$process_data, {
    search_results_path <- searchResultsPathReactive()
    fasta_path <- fastaPathReactive()
    format <- activeFormatReactive()

    startProcessingFn(
      searchResultsPath = search_results_path,
      fastaPath = fasta_path,
      format = format,
      localData = localData,
      ns = session$ns
    )

    runProcessingSafelyFn(
      runProcessing = function() {
        data_import_result <- readDataWithStatusFn(
          format = format,
          searchResultsPath = search_results_path,
          input = input,
          updateStatus = updateStatusFn,
          logInfo = logInfo,
          importDataByFormat = importDataByFormatFn,
          recordImportedData = recordImportedDataFn,
          importDiann = importDiannFn,
          importSpectronaut = importSpectronautFn,
          importFragpipe = importFragpipeFn,
          importMaxquant = importMaxquantFn,
          importPdTmt = importPdTmtFn,
          logError = logError,
          captureCheckpoint = .capture_checkpoint,
          messageFn = messageFn
        )

        applyWorkflowWithStatusFn(
          workflowData = workflowData,
          dataImportResult = data_import_result,
          format = format,
          fastaPath = fasta_path,
          organismName = input$organism_name,
          experimentPaths = experimentPaths,
          sanitizeNames = input$sanitize_names,
          updateStatus = updateStatusFn,
          applyImportResult = applyImportResultFn,
          processFastaData = processFastaDataFn,
          logInfo = logInfo,
          logWarn = logWarn,
          showNotification = showNotificationFn,
          sanitizeRunNames = sanitizeRunNamesFn,
          processFasta = processFastaFn
        )

        runMixedSpeciesAnalysisIfNeededFn(
          workflowData = workflowData,
          localData = localData,
          dataImportResult = data_import_result,
          fastaPath = fasta_path,
          session = session,
          mixedSpeciesFasta = input$mixed_species_fasta,
          analyzeMixedSpeciesData = analyzeMixedSpeciesDataFn,
          updateStatus = updateStatusFn,
          messageFn = messageFn,
          logInfo = logInfo,
          logWarn = logWarn,
          showNotification = showNotificationFn,
          extractOrganisms = extractOrganismsFn,
          analyzeDistribution = analyzeDistributionFn,
          buildSelectionModal = buildSelectionModalFn,
          showModal = showModalFn
        )

        loadConfigurationResourcesFn(
          workflowData = workflowData,
          experimentPaths = experimentPaths,
          useShinyFiles = useShinyFiles,
          configFilePath = localData$config_file_path,
          configFileStandard = input$config_file_standard,
          uniprotMappingFile = localData$uniprot_mapping_file,
          uniprotMappingStandard = input$uniprot_mapping_standard,
          uniparcMappingFile = localData$uniparc_mapping_file,
          uniparcMappingStandard = input$uniparc_mapping_standard,
          loadConfiguration = loadConfigurationFn,
          storeConfiguration = storeConfigurationFn,
          loadOptionalMappings = loadOptionalMappingsFn,
          readConfig = readConfigFn,
          showNotification = showNotificationFn,
          logInfo = logInfo,
          logError = logError,
          getDefaultConfig = getDefaultConfigFn,
          fileExists = fileExistsFn,
          fileCopy = fileCopyFn,
          downloadFile = downloadFileFn,
          systemFileFn = systemFileFn,
          assignFn = assignFn,
          assignEnv = assignEnv,
          readOptionalMapping = readOptionalMappingFn
        )

        recordSetupLogFn(
          workflowData = workflowData,
          dataImportResult = data_import_result,
          format = format,
          useShinyFiles = useShinyFiles,
          searchResultsPath = search_results_path,
          fastaPath = fasta_path,
          searchResultsStandard = input$search_results_standard,
          fastaFileStandard = input$fasta_file_standard,
          taxonId = input$taxon_id,
          organismName = input$organism_name,
          mixedSpeciesFasta = input$mixed_species_fasta,
          resolveFilename = resolveFilenameFn,
          finalizeSetupState = finalizeSetupStateFn,
          logInfo = logInfo
        )

        finalizeProcessingFn(
          workflowData = workflowData,
          localData = localData,
          dataImportResult = data_import_result,
          format = format,
          completeSuccessState = completeSuccessStateFn
        )

        invisible(data_import_result)
      },
      workflowData = workflowData,
      localData = localData,
      logError = logError
    )
  })

  observeEventFn(input$confirm_organism, {
    reqFn(input$selected_organism)
    reqFn(localData$organism_distribution)

    confirmOrganismSelectionFn(
      selectedTaxon = as.integer(input$selected_organism),
      filterToOrganism = input$filter_to_organism,
      workflowData = workflowData,
      localData = localData,
      session = session
    )
  })

  invisible(TRUE)
}
