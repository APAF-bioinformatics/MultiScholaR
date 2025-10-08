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
  message(sprintf("--- Entering designMatrixAppletServer ---"))
  message(sprintf("   designMatrixAppletServer Arg: id = %s", id))
  
  shiny::moduleServer(id, function(input, output, session) {
    message(sprintf("   designMatrixAppletServer Step: Inside moduleServer function"))
    
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
          # Read the raw contrast strings
          contrast_strings <- readLines(contrast_file)
          
          # Reconstruct the full contrasts_tbl format with all required columns
          if (length(contrast_strings) > 0) {
            contrast_info <- lapply(contrast_strings, function(contrast_string) {
              # Parse the contrast string to extract friendly name
              # Expected format: "groupGA_Elevated-groupGA_Control"
              
              # Remove "group" prefixes if present for friendly name
              clean_string <- gsub("^group", "", contrast_string)
              clean_string <- gsub("-group", "-", clean_string)
              
              # Create friendly name by replacing - with _vs_
              friendly_name <- gsub("-", "_vs_", clean_string)
              
              # Create full format: friendly_name=original_contrast_string
              full_format <- paste0(friendly_name, "=", contrast_string)
              
              list(
                contrast_string = contrast_string,
                friendly_name = friendly_name,
                full_format = full_format
              )
            })
            
            # Create properly formatted contrasts_tbl
            data.frame(
              contrasts = sapply(contrast_info, function(x) x$contrast_string),
              friendly_names = sapply(contrast_info, function(x) x$friendly_name),
              full_format = sapply(contrast_info, function(x) x$full_format),
              stringsAsFactors = FALSE
            )
          } else {
            NULL
          }
        } else {
          NULL
        }
        
        # ✅ FIXED: Load aa_seq_tbl_final if it exists in the import directory or scripts
        aa_seq_file_import <- file.path(import_path, "aa_seq_tbl_final.RDS")
        aa_seq_file_scripts <- file.path(experiment_paths$source_dir, "aa_seq_tbl_final.RDS")
        
        if (file.exists(aa_seq_file_import)) {
          logger::log_info("Loading aa_seq_tbl_final from import directory.")
          aa_seq_tbl_final <- readRDS(aa_seq_file_import)
          workflow_data$aa_seq_tbl_final <- aa_seq_tbl_final
          assign("aa_seq_tbl_final", aa_seq_tbl_final, envir = .GlobalEnv)
          
          # Copy to scripts directory for future use
          saveRDS(aa_seq_tbl_final, aa_seq_file_scripts)
          logger::log_info("Copied aa_seq_tbl_final to scripts directory for persistence.")
          
        } else if (file.exists(aa_seq_file_scripts)) {
          logger::log_info("Loading existing aa_seq_tbl_final from scripts directory.")
          aa_seq_tbl_final <- readRDS(aa_seq_file_scripts)
          workflow_data$aa_seq_tbl_final <- aa_seq_tbl_final
          assign("aa_seq_tbl_final", aa_seq_tbl_final, envir = .GlobalEnv)
        } else {
          logger::log_warn("No aa_seq_tbl_final found. Protein accession cleanup will be skipped.")
          workflow_data$aa_seq_tbl_final <- NULL
        }
        
        # ✅ NEW: Load uniprot_dat_cln if it exists (CRITICAL for DE analysis gene names!)
        uniprot_file_import <- file.path(import_path, "uniprot_dat_cln.RDS")
        uniprot_file_scripts <- file.path(experiment_paths$source_dir, "uniprot_dat_cln.RDS")
        
        if (file.exists(uniprot_file_import)) {
          logger::log_info("Loading uniprot_dat_cln from import directory.")
          uniprot_dat_cln <- readRDS(uniprot_file_import)
          workflow_data$uniprot_dat_cln <- uniprot_dat_cln
          assign("uniprot_dat_cln", uniprot_dat_cln, envir = .GlobalEnv)
          
          # Copy to scripts directory for future use
          saveRDS(uniprot_dat_cln, uniprot_file_scripts)
          logger::log_info(sprintf("Copied uniprot_dat_cln to scripts directory (%d annotations).", nrow(uniprot_dat_cln)))
          
          shiny::showNotification(
            sprintf("UniProt annotations loaded: %d protein annotations available", nrow(uniprot_dat_cln)),
            type = "message",
            duration = 5
          )
          
        } else if (file.exists(uniprot_file_scripts)) {
          logger::log_info("Loading existing uniprot_dat_cln from scripts directory.")
          uniprot_dat_cln <- readRDS(uniprot_file_scripts)
          workflow_data$uniprot_dat_cln <- uniprot_dat_cln
          assign("uniprot_dat_cln", uniprot_dat_cln, envir = .GlobalEnv)
          
          shiny::showNotification(
            sprintf("UniProt annotations loaded: %d protein annotations available", nrow(uniprot_dat_cln)),
            type = "message",
            duration = 5
          )
        } else {
          logger::log_warn("No uniprot_dat_cln found. DE analysis will use accession IDs as gene names.")
          workflow_data$uniprot_dat_cln <- NULL
          
          shiny::showNotification(
            "Warning: No UniProt annotations found. Gene names in DE results will be limited.",
            type = "warning",
            duration = 8
          )
        }
        
        # --- Update workflow_data with imported data ---
        workflow_data$design_matrix <- imported_design
        workflow_data$data_cln <- imported_data_cln
        workflow_data$contrasts_tbl <- imported_contrasts
        
        # ✅ FIXED: Save contrasts_tbl to global environment for DE analysis
        if (!is.null(imported_contrasts)) {
          assign("contrasts_tbl", imported_contrasts, envir = .GlobalEnv)
          logger::log_info("Saved contrasts_tbl to global environment for DE analysis.")
        }
        
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
      
      # ✅ FIXED: Save contrasts_tbl to global environment immediately when workflow_data is updated
      if (!is.null(results$contrasts_tbl)) {
        assign("contrasts_tbl", results$contrasts_tbl, envir = .GlobalEnv)
        logger::log_info("Updated contrasts_tbl in global environment from design builder.")
      }
      
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
              
              # ✅ FIXED: Save contrasts_tbl to global environment for DE analysis
              assign("contrasts_tbl", results$contrasts_tbl, envir = .GlobalEnv)
              logger::log_info("Saved contrasts_tbl to global environment for DE analysis.")
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