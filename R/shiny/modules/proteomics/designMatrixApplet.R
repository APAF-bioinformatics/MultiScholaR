#' @title designMatrixAppletModule
#' 
#' @description A Shiny module that serves as the main host for the design
#' matrix creation workflow step. It embeds the design matrix builder UI
#' and handles the logic for saving the results back to the main workflow.
#' 
#' @param id Module ID
#' @param workflow_data A reactive values object to store workflow data.
#' @param experiment_paths A list of paths for the current experiment.
#'
#' @name designMatrixAppletModule
NULL

#' @rdname designMatrixAppletModule
#' @export
#' @import shiny
#' @importFrom shiny NS wellPanel h3 p conditionalPanel div icon
#' @importFrom DT DTOutput
designMatrixAppletUI <- function(id) {
  ns <- NS(id)
  
  shiny::fluidRow(
    shiny::column(12,
      shiny::wellPanel(
        shiny::fluidRow(
          shiny::column(8,
            shiny::h3("Design Matrix Builder")
          ),
          shiny::column(4,
            shiny::actionButton(ns("show_import_modal"), "Import Existing Design", 
                                icon = shiny::icon("folder-open"), class = "btn-info pull-right")
          )
        ),
        shiny::p("Use the tools below to define your experimental groups and the contrasts for differential analysis. Alternatively, import an existing design from a previous analysis."),
        
        # Conditional panel to show the builder only when data is available
        shiny::conditionalPanel(
          condition = paste0("output['", ns("data_available"), "']"),
          # Embed the design matrix builder UI directly
          designMatrixBuilderUI(ns("builder")),
          
          # Add a section for previewing the saved results
          shiny::hr(),
          shiny::h3("Saved Results Preview"),
          shiny::p("This section shows the design matrix and contrasts that have been saved to the workflow."),
              shiny::conditionalPanel(
                condition = paste0("output['", ns("design_matrix_exists"), "']"),
                shiny::wellPanel(
                  shiny::h4("Current Design Matrix"),
                  DT::DTOutput(ns("design_matrix_preview")),
                  shiny::br(),
                  shiny::h4("Defined Contrasts"),
                  DT::DTOutput(ns("contrasts_preview"))
            )
          )
        ),
        
        # Conditional panel to show a message if data is not available
        shiny::conditionalPanel(
          condition = paste0("!output['", ns("data_available"), "']"),
          shiny::div(
            class = "alert alert-info",
            shiny::icon("info-circle"),
            " Please complete the 'Setup & Import' step first. The builder will appear here once data is available."
          )
        )
      )
    )
  )
}

#' @rdname designMatrixAppletModule
#' @export
#' @import shiny
#' @importFrom shiny moduleServer reactive observeEvent req renderDT
#' @importFrom shiny showNotification outputOptions
#' @importFrom logger log_info log_error
#' @importFrom utils write.table
#' @importFrom vroom vroom
#' @importFrom shinyFiles shinyDirButton shinyDirChoose parseDirPath getVolumes
designMatrixAppletServer <- function(id, workflow_data, experiment_paths, volumes = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # == Setup shinyFiles =======================================================
    # Resolve the volumes object to a static list to avoid reactive context issues.
    # This checks if 'volumes' is the getVolumes() function and calls it,
    # or creates it if it's NULL.
    resolved_volumes <- shiny::isolate({
      base_volumes <- if (is.function(volumes)) {
        volumes()
      } else if (is.null(volumes)) {
        shinyFiles::getVolumes()()
      } else {
        volumes
      }
      
      # ✅ Add base_dir as default starting location for import
      if (!is.null(experiment_paths) && !is.null(experiment_paths$base_dir) && 
          dir.exists(experiment_paths$base_dir)) {
        # Put base_dir first so it's the default
        enhanced_volumes <- c("Project Base Dir" = experiment_paths$base_dir, base_volumes)
        logger::log_info(paste("Added base_dir to volumes for easier navigation:", experiment_paths$base_dir))
        enhanced_volumes
      } else {
        base_volumes
      }
    })

    shinyFiles::shinyDirChoose(input, "import_dir", roots = resolved_volumes, session = session)

    # == Modal Logic for Import =================================================
    
    observeEvent(input$show_import_modal, {
      ns <- session$ns
      showModal(modalDialog(
        title = "Import Existing Design Matrix",
        shiny::p("Select the folder containing 'design_matrix.tab' and 'data_cln.tab' files."),
        shinyFiles::shinyDirButton(ns("import_dir"), "Select Folder", "Choose a directory"),
        verbatimTextOutput(ns("import_dir_path"), placeholder = TRUE),
        footer = tagList(
          modalButton("Cancel"),
          actionButton(ns("confirm_import"), "Import", class = "btn-primary")
        )
      ))
    })

    # Display selected directory path in modal
    output$import_dir_path <- renderText({
      req(input$import_dir)
      shinyFiles::parseDirPath(resolved_volumes, input$import_dir)
    })

    # == Handle Import Confirmation =============================================
    
    observeEvent(input$confirm_import, {
      req(input$import_dir)
      
      import_path <- shinyFiles::parseDirPath(resolved_volumes, input$import_dir)
      req(import_path)
      
      removeModal()
      
      design_file <- file.path(import_path, "design_matrix.tab")
      cln_data_file <- file.path(import_path, "data_cln.tab")
      contrast_file <- file.path(import_path, "contrast_strings.tab")
      
      # Validate that files exist
      if (!file.exists(design_file) || !file.exists(cln_data_file)) {
        msg <- "Import failed: 'design_matrix.tab' and/or 'data_cln.tab' not found in the selected directory."
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 10)
        return()
      }
      
      shiny::showNotification("Importing design files...", id = "importing_design", duration = NULL)
      
      tryCatch({
        # ---[ ADDED LOGIC FOR CONFIG FILE ] ---
        # Ensure a config file is available, downloading a default if necessary.
        config_path <- file.path(experiment_paths$source_dir, "config.ini")
        
        if (!file.exists(config_path)) {
          logger::log_info("config.ini not found in project. Downloading default config.")
          tryCatch({
            default_config_url <- "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/Workbooks/config.ini"
            download.file(default_config_url, destfile = config_path, quiet = TRUE)
            logger::log_info(paste("Default config.ini downloaded to:", config_path))
            shiny::showNotification("Using default config.ini.", type = "message")
          }, error = function(e_download) {
            msg <- paste("Failed to download default config.ini:", e_download$message)
            logger::log_error(msg)
            shiny::showNotification(msg, type = "error", duration = 10)
            # Stop execution if config is essential and cannot be obtained
            stop("Could not obtain a configuration file.")
          })
        }
        
        # Now, load the config_list into the workflow
        logger::log_info("Loading config.ini from project.")
        workflow_data$config_list <- readConfigFile(file = config_path)
        
        # ✅ FIXED: Create global config_list for updateConfigParameter compatibility
        assign("config_list", workflow_data$config_list, envir = .GlobalEnv)
        logger::log_info("Created global config_list for updateConfigParameter compatibility")
        # ---[ END OF ADDED LOGIC ] ---

        # Read files
        imported_design <- vroom::vroom(design_file, show_col_types = FALSE)
        imported_data_cln <- vroom::vroom(cln_data_file, show_col_types = FALSE)
        
        imported_contrasts <- if(file.exists(contrast_file)) {
          data.frame(contrasts = readLines(contrast_file))
        } else {
          NULL
        }
        
        # --- Update workflow_data with imported data ---
        workflow_data$design_matrix <- imported_design
        workflow_data$data_cln <- imported_data_cln
        workflow_data$contrasts_tbl <- imported_contrasts
        
        # --- Create S4 Object and Save State (same logic as builder) ---
        logger::log_info("Creating PeptideQuantitativeData S4 object from imported design.")
        peptide_data_s4 <- new(
          "PeptideQuantitativeData",
          peptide_data = workflow_data$data_cln,
          protein_id_column = "Protein.Ids",
          peptide_sequence_column = "Stripped.Sequence",
          q_value_column = "Q.Value",
          global_q_value_column = "Global.Q.Value",
          proteotypic_peptide_sequence_column = "Proteotypic",
          raw_quantity_column = "Precursor.Quantity",
          norm_quantity_column = "Precursor.Normalised",
          is_logged_data = FALSE,
          design_matrix = workflow_data$design_matrix,
          sample_id = "Run",
          group_id = "group",
          technical_replicate_id = "replicates",
          args = workflow_data$config_list
        )
        
        logger::log_info("Saving imported S4 object to R6 state manager as 'raw_data_s4'.")
        workflow_data$state_manager$saveState(
            state_name = "raw_data_s4",
            s4_data_object = peptide_data_s4,
            config_object = workflow_data$config_list,
            description = "S4 object created from imported design matrix."
        )
        
        workflow_data$tab_status$design_matrix <- "complete"
        
        shiny::removeNotification("importing_design")
        shiny::showNotification("Design imported successfully!", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error during import:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("importing_design")
      })
    })

    # == Reactivity Checks ======================================================
    
    # Check if data is available from the previous step
    output$data_available <- shiny::reactive({
      !is.null(workflow_data$data_tbl) && !is.null(workflow_data$config_list)
    })
    shiny::outputOptions(output, "data_available", suspendWhenHidden = FALSE)
    
    # Check if a design matrix has been saved to the workflow
    output$design_matrix_exists <- shiny::reactive({
      !is.null(workflow_data$design_matrix)
    })
    shiny::outputOptions(output, "design_matrix_exists", suspendWhenHidden = FALSE)
    
    # == Module Integration =====================================================
    
    # Call the builder module server and store its returned reactive
    # Pass the required data as reactive expressions
    builder_results_rv <- designMatrixBuilderServer(
      "builder",
      data_tbl = shiny::reactive(workflow_data$data_tbl),
      config_list = shiny::reactive(workflow_data$config_list)
    )
    
    # == Handle Builder Results =================================================
    
    # This observer fires when the user clicks "Save Design" in the builder module
    shiny::observeEvent(builder_results_rv(), {
      results <- builder_results_rv()
      shiny::req(results)
      
      logger::log_info("Received results from design matrix builder module. Now saving to workflow and disk.")
      
      # 1. Update the main workflow_data object
      workflow_data$design_matrix <- results$design_matrix
      workflow_data$data_cln <- results$data_cln
      workflow_data$contrasts_tbl <- results$contrasts_tbl
      workflow_data$config_list <- results$config_list
      
      # ✅ FIXED: Create global config_list for updateConfigParameter compatibility
      assign("config_list", workflow_data$config_list, envir = .GlobalEnv)
      logger::log_info("Updated global config_list for updateConfigParameter compatibility")
      
      # 2. Get the source directory for writing config files
      source_dir <- experiment_paths$source_dir
      if (is.null(source_dir) || !dir.exists(source_dir)) {
          msg <- "Could not find the source directory to save files. Check experiment paths."
          logger::log_error(msg)
          shiny::showNotification(msg, type = "error", duration = 15)
          return()
      }
      
      # 3. Perform the file-writing operations
      tryCatch({
          # --- Save Design Matrix ---
          design_matrix_path <- file.path(source_dir, "design_matrix.tab")
          logger::log_info(paste("Writing design matrix to:", design_matrix_path))
          utils::write.table(results$design_matrix, file = design_matrix_path, sep = "\t", row.names = FALSE, quote = FALSE)
          
          # --- Save Cleaned Data ---
          data_cln_path <- file.path(source_dir, "data_cln.tab")
          logger::log_info(paste("Writing cleaned data to:", data_cln_path))
          utils::write.table(results$data_cln, file = data_cln_path, sep = "\t", row.names = FALSE, quote = FALSE)
          
          # --- Save Contrasts ---
          if (!is.null(results$contrasts_tbl) && nrow(results$contrasts_tbl) > 0) {
              contrast_path <- file.path(source_dir, "contrast_strings.tab")
              logger::log_info(paste("Writing contrasts to:", contrast_path))
              writeLines(results$contrasts_tbl$contrasts, contrast_path)
          }

          # --- Create S4 Object and Save Initial State ---
          logger::log_info("Creating PeptideQuantitativeData S4 object.")
          peptide_data_s4 <- new(
            "PeptideQuantitativeData",
            peptide_data = workflow_data$data_cln,
            protein_id_column = "Protein.Ids",
            peptide_sequence_column = "Stripped.Sequence",
            q_value_column = "Q.Value",
            global_q_value_column = "Global.Q.Value",
            proteotypic_peptide_sequence_column = "Proteotypic",
            raw_quantity_column = "Precursor.Quantity",
            norm_quantity_column = "Precursor.Normalised",
            is_logged_data = FALSE,
            design_matrix = workflow_data$design_matrix,
            sample_id = "Run",
            group_id = "group",
            technical_replicate_id = "replicates",
            args = workflow_data$config_list
          )
          
          logger::log_info("Saving S4 object to R6 state manager as 'raw_data_s4'.")
          workflow_data$state_manager$saveState(
              state_name = "raw_data_s4",
              s4_data_object = peptide_data_s4,
              config_object = workflow_data$config_list,
              description = "Initial S4 object created after design matrix definition."
          )
          
          # 4. Update tab status to 'complete'
          workflow_data$tab_status$design_matrix <- "complete"
          
          shiny::showNotification("Design matrix and contrasts saved successfully!", type = "message")
          
      }, error = function(e) {
          msg <- paste("Error saving design matrix results:", e$message)
          logger::log_error(msg)
          shiny::showNotification(msg, type = "error", duration = 15)
      })
      
    }, ignoreNULL = TRUE) # ignoreNULL is crucial
    
    # == Previews of Saved Data =================================================
    
    # Display design matrix preview from the main workflow_data
    output$design_matrix_preview <- DT::renderDT({
      shiny::req(workflow_data$design_matrix)
      workflow_data$design_matrix
    }, options = list(pageLength = 5, scrollX = TRUE))
    
    # Display contrasts preview from the main workflow_data
    output$contrasts_preview <- DT::renderDT({
      shiny::req(workflow_data$contrasts_tbl)
      workflow_data$contrasts_tbl
    }, options = list(pageLength = 5, scrollX = TRUE))
    
  })
} 