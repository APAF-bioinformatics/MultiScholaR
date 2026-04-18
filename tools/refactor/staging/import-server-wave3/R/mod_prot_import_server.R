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
mod_prot_import_server <- function(id, workflow_data, experiment_paths, volumes = NULL) {
  message(sprintf("--- Entering mod_prot_import_server ---"))
  message(sprintf("   mod_prot_import_server Arg: id = %s", id))
  
  shiny::moduleServer(id, function(input, output, session) {
    message(sprintf("   mod_prot_import_server Step: Inside moduleServer function"))
    
    # Local reactive values
    local_data <- shiny::reactiveValues(
      processing = FALSE,
      detected_format = NULL,
      format_confidence = NULL,
      # Store file paths for shinyFiles
      search_results_file = NULL,
      fasta_file_path = NULL,
      uniprot_mapping_file = NULL,
      uniparc_mapping_file = NULL,
      config_file_path = NULL,
      # Mixed species FASTA analysis
      organism_distribution = NULL,
      organism_mapping = NULL,
      pending_import_completion = FALSE,
      waiting_for_organism_selection = FALSE
    )
    
    # Check if shinyFiles is available
    use_shiny_files <- requireNamespace("shinyFiles", quietly = TRUE)
    
    # Set up shinyFiles if available
    if (use_shiny_files) {
      volumes <- resolveProtImportShinyFileVolumes(
        volumes = volumes
      )

      registerProtImportShinyFileChooser(
        input = input
        , inputId = "search_results"
        , volumes = volumes
        , session = session
        , filetypes = c("tsv", "txt", "tab", "csv", "xlsx", "zip", "parquet")
      )
      
      observeEvent(input$search_results, {
        handleProtImportShinyFileSelection(
          selectedInput = input$search_results
          , volumes = volumes
          , localData = local_data
          , localField = "search_results_file"
          , output = output
          , outputId = "search_results_path"
          , catchErrors = TRUE
        )
      })
      
      registerProtImportShinyFileChooser(
        input = input
        , inputId = "fasta_file"
        , volumes = volumes
        , session = session
        , filetypes = c("fasta", "fa", "faa")
      )
      
      observeEvent(input$fasta_file, {
        handleProtImportShinyFileSelection(
          selectedInput = input$fasta_file
          , volumes = volumes
          , localData = local_data
          , localField = "fasta_file_path"
          , output = output
          , outputId = "fasta_file_path"
        )
      })
      
      registerProtImportShinyFileChooser(
        input = input
        , inputId = "uniprot_mapping"
        , volumes = volumes
        , session = session
        , filetypes = c("tsv", "txt", "tab")
      )
      
      observeEvent(input$uniprot_mapping, {
        handleProtImportShinyFileSelection(
          selectedInput = input$uniprot_mapping
          , volumes = volumes
          , localData = local_data
          , localField = "uniprot_mapping_file"
          , output = output
          , outputId = "uniprot_mapping_path"
        )
      })
      
      registerProtImportShinyFileChooser(
        input = input
        , inputId = "uniparc_mapping"
        , volumes = volumes
        , session = session
        , filetypes = c("tsv", "txt", "tab")
      )
      
      observeEvent(input$uniparc_mapping, {
        handleProtImportShinyFileSelection(
          selectedInput = input$uniparc_mapping
          , volumes = volumes
          , localData = local_data
          , localField = "uniparc_mapping_file"
          , output = output
          , outputId = "uniparc_mapping_path"
        )
      })
      
      registerProtImportShinyFileChooser(
        input = input
        , inputId = "config_file"
        , volumes = volumes
        , session = session
        , filetypes = c("ini", "cfg")
      )
      
      observeEvent(input$config_file, {
        handleProtImportShinyFileSelection(
          selectedInput = input$config_file
          , volumes = volumes
          , localData = local_data
          , localField = "config_file_path"
          , output = output
          , outputId = "config_file_path"
        )
      })
    }
    
    # --- TESTTHAT CHECKPOINT CAPTURE (TEMPORARY) ---
    shiny::observeEvent(input$capture_checkpoints, {
      updateProtImportCheckpointCaptureOption(
        captureCheckpoints = input$capture_checkpoints
        , optionsFn = options
        , logInfo = logger::log_info
      )
    }, ignoreInit = FALSE)
    # --- END CHECKPOINT CAPTURE ---

    # --- Observer: Toggle Step 3 inputs when mixed species is checked ---
    shiny::observeEvent(input$mixed_species_fasta, {
      toggleProtImportMixedSpeciesInputs(
        mixedSpeciesFasta = input$mixed_species_fasta
        , disableInput = shinyjs::disable
        , enableInput = shinyjs::enable
        , messageFn = message
      )
    }, ignoreInit = FALSE)
    
    # Get the search results file path (works for both shinyFiles and standard)
    get_search_results_path <- reactive({
      resolveProtImportPrimaryUploadPath(
        useShinyFiles = use_shiny_files
        , shinyPath = local_data$search_results_file
        , standardInput = input$search_results_standard
      )
    })
    
    # Get the FASTA file path
    get_fasta_path <- reactive({
      resolveProtImportPrimaryUploadPath(
        useShinyFiles = use_shiny_files
        , shinyPath = local_data$fasta_file_path
        , standardInput = input$fasta_file_standard
      )
    })
    
    # Detect file format when file is uploaded
    shiny::observeEvent(get_search_results_path(), {
      file_path <- get_search_results_path()
      shiny::req(file_path)

      runProtImportFormatDetection(
        filePath = file_path
        , useShinyFiles = use_shiny_files
        , localData = local_data
        , searchResultsStandard = input$search_results_standard
        , logInfo = log_info
        , logError = log_error
      )
    })
    
    # Get active format (either detected or overridden)
    active_format <- shiny::reactive({
      resolveActiveProtImportFormat(
        formatOverride = input$format_override
        , detectedFormat = local_data$detected_format
      )
    })
    
    # Render format detection output
    output$format_detection <- shiny::renderUI({
      buildProtImportFormatDetectionUi(
        detectedFormat = local_data$detected_format,
        formatConfidence = local_data$format_confidence
      )
    })
    
    # Render format-specific options
    output$format_specific_options <- shiny::renderUI({
      buildProtImportFormatSpecificOptionsUi(
        format = active_format(),
        ns = session$ns
      )
    })
    
    # Process imported data
    shiny::observeEvent(input$process_data, {
      # Get file paths
      search_results_path <- get_search_results_path()
      fasta_path <- get_fasta_path()
      format <- active_format()

      startProtImportProcessing(
        searchResultsPath = search_results_path
        , fastaPath = fasta_path
        , format = format
        , localData = local_data
        , ns = session$ns
      )

      runProtImportProcessingSafely(
        runProcessing = function() {
          data_import_result <- readProtImportDataWithStatus(
            format = format
            , searchResultsPath = search_results_path
            , input = input
            , updateStatus = updateProtImportProcessingStatus
            , logInfo = log_info
            , importDataByFormat = importProtImportDataByFormat
            , recordImportedData = recordProtImportImportedData
            , importDiann = importDIANNData
            , importSpectronaut = importSpectronautData
            , importFragpipe = importFragPipeData
            , importMaxquant = importMaxQuantData
            , importPdTmt = importProteomeDiscovererTMTData
            , logError = log_error
            , captureCheckpoint = .capture_checkpoint
            , messageFn = message
          )

          applyProtImportWorkflowWithStatus(
            workflowData = workflow_data
            , dataImportResult = data_import_result
            , format = format
            , fastaPath = fasta_path
            , organismName = input$organism_name
            , experimentPaths = experiment_paths
            , sanitizeNames = input$sanitize_names
            , updateStatus = updateProtImportProcessingStatus
            , applyImportResult = applyProtImportResultToWorkflow
            , processFastaData = processProtImportFastaData
            , logInfo = log_info
            , logWarn = log_warn
            , showNotification = shiny::showNotification
            , sanitizeRunNames = sanitizeProtImportRunNames
            , processFasta = processFastaFile
          )

          runProtImportMixedSpeciesAnalysisIfNeeded(
            workflowData = workflow_data
            , localData = local_data
            , dataImportResult = data_import_result
            , fastaPath = fasta_path
            , session = session
            , mixedSpeciesFasta = input$mixed_species_fasta
            , analyzeMixedSpeciesData = analyzeProtImportMixedSpeciesData
            , updateStatus = updateProtImportProcessingStatus
            , messageFn = message
            , logInfo = log_info
            , logWarn = log_warn
            , showNotification = shiny::showNotification
            , extractOrganisms = extractOrganismsFromFasta
            , analyzeDistribution = analyzeOrganismDistribution
            , buildSelectionModal = buildProtImportOrganismSelectionModal
            , showModal = shiny::showModal
          )

          loadProtImportConfigurationResources(
            workflowData = workflow_data
            , experimentPaths = experiment_paths
            , useShinyFiles = use_shiny_files
            , configFilePath = local_data$config_file_path
            , configFileStandard = input$config_file_standard
            , uniprotMappingFile = local_data$uniprot_mapping_file
            , uniprotMappingStandard = input$uniprot_mapping_standard
            , uniparcMappingFile = local_data$uniparc_mapping_file
            , uniparcMappingStandard = input$uniparc_mapping_standard
            , loadConfiguration = loadProtImportConfiguration
            , storeConfiguration = storeProtImportConfiguration
            , loadOptionalMappings = loadProtImportOptionalMappings
            , readConfig = readConfigFile
            , showNotification = shiny::showNotification
            , logInfo = log_info
            , logError = log_error
            , getDefaultConfig = getDefaultProteomicsConfig
            , fileExists = file.exists
            , fileCopy = file.copy
            , downloadFile = download.file
            , systemFileFn = system.file
            , assignFn = assign
            , assignEnv = .GlobalEnv
            , readOptionalMapping = readProtImportOptionalMapping
          )

          recordProtImportSetupLog(
            workflowData = workflow_data
            , dataImportResult = data_import_result
            , format = format
            , useShinyFiles = use_shiny_files
            , searchResultsPath = search_results_path
            , fastaPath = fasta_path
            , searchResultsStandard = input$search_results_standard
            , fastaFileStandard = input$fasta_file_standard
            , taxonId = input$taxon_id
            , organismName = input$organism_name
            , mixedSpeciesFasta = input$mixed_species_fasta
            , resolveFilename = resolveProtImportInputFilename
            , finalizeSetupState = finalizeProtImportSetupState
            , logInfo = log_info
          )

          finalizeProtImportProcessing(
            workflowData = workflow_data
            , localData = local_data
            , dataImportResult = data_import_result
            , format = format
            , completeSuccessState = completeProtImportSuccessState
          )

          invisible(data_import_result)
        }
        , workflowData = workflow_data
        , localData = local_data
        , logError = log_error
      )
    })
    
    # --- Organism Distribution Table Output ---
    output$organism_dist_table <- DT::renderDT({
      shiny::req(local_data$organism_distribution)

      buildProtImportOrganismDistributionTable(local_data$organism_distribution)
    })
    
    # --- Observer for Organism Selection Confirmation ---
    shiny::observeEvent(input$confirm_organism, {
      shiny::req(input$selected_organism)
      shiny::req(local_data$organism_distribution)

      confirmProtImportOrganismSelection(
        selectedTaxon = as.integer(input$selected_organism)
        , filterToOrganism = input$filter_to_organism
        , workflowData = workflow_data
        , localData = local_data
        , session = session
      )
    })
    
    # Render import status
    output$import_status <- shiny::renderUI({
      buildProtImportStatusUi(
        isProcessing = local_data$processing,
        dataTbl = workflow_data$data_tbl,
        dataFormat = workflow_data$data_format,
        dataType = workflow_data$data_type
      )
    })
    
    # Render data summary
    output$data_summary <- shiny::renderUI({
      buildProtImportDataSummaryUi(
        dataTbl = workflow_data$data_tbl,
        setupImportLog = workflow_data$processing_log$setup_import
      )
    })
    
  })
}

