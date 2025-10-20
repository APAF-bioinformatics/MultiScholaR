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
  
  shiny::tagList(
    # JavaScript handler for UniProt progress updates
    shiny::tags$script(shiny::HTML("
      Shiny.addCustomMessageHandler('updateUniprotProgress', function(message) {
        $('#uniprot_progress_bar').css('width', message.percent + '%');
        $('#uniprot_progress_bar').text(message.percent + '%');
        $('#uniprot_progress_text').text(message.text);
      });
    ")),
    
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
#' @importFrom tidyr pivot_wider
#' @importFrom rlang sym
designMatrixAppletServer <- function(id, workflow_data, experiment_paths, volumes = NULL, qc_trigger = NULL) {
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
        
        # --- Detect or Load workflow_type ---
        # Try to get workflow_type from config.ini first
        workflow_type <- NULL
        if (!is.null(workflow_data$config_list$globalParameters$workflow_type)) {
          workflow_type <- workflow_data$config_list$globalParameters$workflow_type
          logger::log_info(paste("Loaded workflow_type from config.ini:", workflow_type))
        } else {
          # Fallback: detect from data structure
          logger::log_warn("workflow_type not found in config.ini. Attempting to detect from data structure.")
          
          # Check for peptide-specific columns (DIA/LFQ have these, TMT doesn't)
          has_precursor_cols <- any(c("Precursor.Id", "Precursor.Charge", "Precursor.Quantity") %in% names(imported_data_cln))
          has_peptide_cols <- any(c("Stripped.Sequence", "Modified.Sequence") %in% names(imported_data_cln))
          
          if (has_precursor_cols && has_peptide_cols) {
            # Peptide-level data (DIA or LFQ) - default to DIA as it's more common
            workflow_type <- "DIA"
            logger::log_info("Detected peptide-level data. Setting workflow_type to DIA.")
          } else {
            # Protein-level data (TMT)
            workflow_type <- "TMT"
            logger::log_info("Detected protein-level data. Setting workflow_type to TMT.")
          }
        }
        
        # Set workflow_type in state manager
        workflow_data$state_manager$setWorkflowType(workflow_type)
        logger::log_info(paste("Workflow type set to:", workflow_type))
        
        # Ensure column_mapping is available for pivot operations
        if (is.null(workflow_data$column_mapping)) {
          logger::log_warn("column_mapping not found in workflow_data. Inferring from data structure.")
          
          # Infer basic column mappings from data
          if (workflow_type == "TMT") {
            workflow_data$column_mapping <- list(
              protein_col = "Protein.Ids",
              run_col = "Run",
              quantity_col = "Abundance"
            )
          } else {
            workflow_data$column_mapping <- list(
              protein_col = "Protein.Ids",
              run_col = "Run",
              quantity_col = "Precursor.Quantity"
            )
          }
          logger::log_info("Inferred column_mapping for workflow operations.")
        }
        
        # --- Create S4 Object and Save State (Conditional Orchestration) ---
        if (workflow_type %in% c("DIA", "LFQ")) {
          # For peptide workflows, the data_cln is already in the correct long format.
          # Use the existing constructor for peptide-level data.
          s4_object <- PeptideQuantitativeDataDiann(
            peptide_data = workflow_data$data_cln,
            design_matrix = workflow_data$design_matrix,
            args = workflow_data$config_list
          )
          state_name <- "raw_data_s4"
          description <- "Initial Peptide S4 object created after design matrix."

        } else if (workflow_type == "TMT") {
          # For protein workflows (TMT), the data_cln is long, but the S4 constructor needs wide.
          # We must reshape the data first.
          
          logger::log_info("TMT workflow: Reshaping protein data from long to wide format for S4 creation.")
          
          # Validate column_mapping before pivot
          if (is.null(workflow_data$column_mapping$protein_col) || 
              is.null(workflow_data$column_mapping$run_col) ||
              is.null(workflow_data$column_mapping$quantity_col)) {
            stop("Missing required column mappings for TMT data reshape. Cannot proceed.")
          }
          
          # Validate columns exist in data
          required_cols <- c(workflow_data$column_mapping$protein_col, 
                             workflow_data$column_mapping$run_col,
                             workflow_data$column_mapping$quantity_col)
          missing_cols <- required_cols[!required_cols %in% names(workflow_data$data_cln)]
          if (length(missing_cols) > 0) {
            stop(paste("Missing required columns in data_cln:", paste(missing_cols, collapse = ", ")))
          }
          
          protein_quant_table_wide <- workflow_data$data_cln |>
            tidyr::pivot_wider(
              id_cols = !!sym(workflow_data$column_mapping$protein_col),
              names_from = !!sym(workflow_data$column_mapping$run_col),
              values_from = !!sym(workflow_data$column_mapping$quantity_col)
            )
          
          logger::log_info(sprintf("TMT Import: Pivoted data from %d rows (long) to %d rows (wide)", 
                                   nrow(workflow_data$data_cln), 
                                   nrow(protein_quant_table_wide)))
          logger::log_info(sprintf("TMT Import: Wide format has %d protein rows, %d sample columns",
                                   nrow(protein_quant_table_wide),
                                   ncol(protein_quant_table_wide) - 1))  # -1 for protein ID column
          
          # *** CRITICAL: Apply log2 transformation to TMT abundance values ***
          # TMT data from Proteome Discoverer is in raw abundance scale, not log2
          # Limma and normalization require log2-transformed data
          logger::log_info("TMT Import: Applying log2 transformation to abundance values...")
          
          protein_id_col <- workflow_data$column_mapping$protein_col
          protein_quant_table_wide <- protein_quant_table_wide |>
            dplyr::mutate(
              dplyr::across(
                -!!sym(protein_id_col),  # Transform all columns except protein ID
                ~ {
                  # Add pseudo-count to avoid log(0), then log2 transform
                  # Using 1 as pseudo-count (standard for proteomics)
                  log2(.x + 1)
                }
              )
            )
          
          logger::log_info("TMT Import: Log2 transformation completed")

          # Use the existing constructor for protein-level data with the now wide table.
          s4_object <- ProteinQuantitativeData(
            protein_quant_table = protein_quant_table_wide,
            protein_id_column = workflow_data$column_mapping$protein_col,
            protein_id_table = protein_quant_table_wide |> dplyr::distinct(!!sym(workflow_data$column_mapping$protein_col)),
            design_matrix = workflow_data$design_matrix,
            sample_id = "Run",
            group_id = "group",
            technical_replicate_id = "replicates",
            args = workflow_data$config_list
          )
          state_name <- "protein_s4_initial"
          description <- "Initial Protein S4 object created from TMT data after design matrix."

        } else {
          stop("Unknown workflow_type in designMatrixAppletServer: ", workflow_type)
        }
        
        # Save the dynamically created S4 object to the state manager
        logger::log_info(paste("Saving S4 object to R6 state manager as:", state_name))
        workflow_data$state_manager$saveState(
            state_name = state_name,
            s4_data_object = s4_object,
            config_object = workflow_data$config_list,
            description = description
        )
        logger::log_info(sprintf("Import: S4 object saved to R6 state manager as '%s'", state_name))
        logger::log_info("Import: This state is now ACTIVE - QC modules will read from it")
        logger::log_info("Import: User can proceed to QC → Accession Cleanup")
        
        # --- TRIGGER UNIPROT ANNOTATION ---
        log_info("Design Matrix complete. Triggering UniProt annotation.")
        
        # Show modal with progress bar
        shiny::showModal(shiny::modalDialog(
          title = "Retrieving UniProt Annotations",
          shiny::p("Please wait while we retrieve protein annotations from UniProt..."),
          shiny::div(id = "uniprot_progress_text", "Initializing..."),
          shiny::tags$div(class = "progress",
            shiny::tags$div(id = "uniprot_progress_bar", 
                     class = "progress-bar progress-bar-striped active",
                     role = "progressbar",
                     style = "width: 0%",
                     "0%")
          ),
          footer = NULL,
          easyClose = FALSE
        ))
        
        tryCatch({
          protein_column <- workflow_data$column_mapping$protein_col
          if (is.null(protein_column)) {
            stop("Protein column not found in column_mapping. Cannot get annotations.")
          }
          
          cache_dir <- if (!is.null(experiment_paths) && !is.null(experiment_paths$results_dir)) {
            file.path(experiment_paths$results_dir, "cache")
          } else {
            file.path(tempdir(), "proteomics_cache")
          }
          
          uniprot_cache_dir <- file.path(cache_dir, "uniprot_annotations")
          if (!dir.exists(uniprot_cache_dir)) {
            dir.create(uniprot_cache_dir, recursive = TRUE)
          }
          
          # Create progress callback function with error handling
          progress_updater <- function(current, total) {
            tryCatch({
              percent <- round((current / total) * 100)
              session$sendCustomMessage("updateUniprotProgress", list(
                percent = percent,
                text = sprintf("Processing chunk %d of %d (%d%%)", current, total, percent)
              ))
            }, error = function(e) {
              # Silently catch any errors to prevent disrupting the main process
              message(paste("Progress update failed:", e$message))
            })
          }

          uniprot_dat_cln <- getUniprotAnnotationsFull(
            data_tbl = workflow_data$data_cln,
            protein_id_column = protein_column,
            cache_dir = uniprot_cache_dir,
            taxon_id = workflow_data$taxon_id,
            progress_callback = progress_updater
          )
          
          workflow_data$uniprot_dat_cln <- uniprot_dat_cln
          assign("uniprot_dat_cln", uniprot_dat_cln, envir = .GlobalEnv)
          
          if (!is.null(experiment_paths) && !is.null(experiment_paths$source_dir)) {
            scripts_uniprot_path <- file.path(experiment_paths$source_dir, "uniprot_dat_cln.RDS")
            saveRDS(uniprot_dat_cln, scripts_uniprot_path)
            log_info(sprintf("Saved uniprot_dat_cln to scripts directory: %s", scripts_uniprot_path))
          }
          log_info(sprintf("UniProt annotations retrieved successfully. Found %d annotations", nrow(uniprot_dat_cln)))
          shiny::removeModal()
          shiny::showNotification("UniProt annotations retrieved successfully.", type = "message")
          
        }, error = function(e) {
          log_warn(paste("Error getting UniProt annotations:", e$message))
          workflow_data$uniprot_dat_cln <- NULL
          shiny::removeModal()
          shiny::showNotification(paste("Warning: Could not retrieve UniProt annotations:", e$message), type = "warning", duration = 8)
        })
        
        # Set the QC trigger to TRUE to signal that QC modules can now run
        if (!is.null(qc_trigger)) {
          qc_trigger(TRUE)
        }
        
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
      config_list = shiny::reactive(workflow_data$config_list),
      column_mapping = shiny::reactive(workflow_data$column_mapping)
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
          
          # --- Save config.ini for Import Compatibility ---
          if (!is.null(workflow_data$config_list)) {
            # Add workflow_type to config for import detection
            if (!is.null(workflow_data$state_manager$workflow_type)) {
              workflow_data$config_list$globalParameters$workflow_type <- workflow_data$state_manager$workflow_type
              logger::log_info(paste("Added workflow_type to config.ini:", workflow_data$state_manager$workflow_type))
            }
            
            config_path <- file.path(source_dir, "config.ini")
            logger::log_info(paste("Writing config.ini to:", config_path))
            tryCatch({
              # Write config list back to INI format
              ini::write.ini(workflow_data$config_list, config_path)
              logger::log_info("Saved config.ini for future import/export.")
            }, error = function(e) {
              logger::log_warn(paste("Could not save config.ini:", e$message))
              logger::log_warn("Import will use default config if this design is imported.")
            })
          } else {
            logger::log_warn("config_list is NULL, cannot save config.ini.")
          }

          # --- Create S4 Object and Save State (Conditional Orchestration) ---
          workflow_type <- shiny::isolate(workflow_data$state_manager$workflow_type)
          
          if (workflow_type %in% c("DIA", "LFQ")) {
            # For peptide workflows, the data_cln is already in the correct long format.
            # Use the existing constructor for peptide-level data.
            s4_object <- PeptideQuantitativeDataDiann(
              peptide_data = workflow_data$data_cln,
              design_matrix = workflow_data$design_matrix,
              args = workflow_data$config_list
            )
            state_name <- "raw_data_s4"
            description <- "Initial Peptide S4 object created after design matrix."

          } else if (workflow_type == "TMT") {
            # For protein workflows (TMT), the data_cln is long, but the S4 constructor needs wide.
            # We must reshape the data first.
            
            logger::log_info("TMT workflow: Reshaping protein data from long to wide format for S4 creation.")
            
            protein_quant_table_wide <- workflow_data$data_cln |>
              tidyr::pivot_wider(
                id_cols = !!sym(workflow_data$column_mapping$protein_col),
                names_from = !!sym(workflow_data$column_mapping$run_col),
                values_from = !!sym(workflow_data$column_mapping$quantity_col)
              )
            
            # *** CRITICAL: Apply log2 transformation to TMT abundance values ***
            # TMT data from Proteome Discoverer is in raw abundance scale, not log2
            # Limma and normalization require log2-transformed data
            logger::log_info("TMT Save Design: Applying log2 transformation to abundance values...")
            
            protein_id_col <- workflow_data$column_mapping$protein_col
            protein_quant_table_wide <- protein_quant_table_wide |>
              dplyr::mutate(
                dplyr::across(
                  -!!sym(protein_id_col),  # Transform all columns except protein ID
                  ~ {
                    # Add pseudo-count to avoid log(0), then log2 transform
                    # Using 1 as pseudo-count (standard for proteomics)
                    log2(.x + 1)
                  }
                )
              )
            
            logger::log_info("TMT Save Design: Log2 transformation completed")

            # Use the existing constructor for protein-level data with the now wide table.
            s4_object <- ProteinQuantitativeData(
              protein_quant_table = protein_quant_table_wide,
              protein_id_column = workflow_data$column_mapping$protein_col,
              protein_id_table = protein_quant_table_wide |> dplyr::distinct(!!sym(workflow_data$column_mapping$protein_col)),
              design_matrix = workflow_data$design_matrix,
              sample_id = "Run",
              group_id = "group",
              technical_replicate_id = "replicates",
              args = workflow_data$config_list
            )
            state_name <- "protein_s4_initial"
            description <- "Initial Protein S4 object created from TMT data after design matrix."

          } else {
            stop("Unknown workflow_type in designMatrixAppletServer: ", workflow_type)
          }
          
          # Save the dynamically created S4 object to the state manager
          logger::log_info(paste("Saving S4 object to R6 state manager as:", state_name))
          workflow_data$state_manager$saveState(
              state_name = state_name,
              s4_data_object = s4_object,
              config_object = workflow_data$config_list,
              description = description
          )
          
          # --- TRIGGER UNIPROT ANNOTATION ---
          log_info("Design Matrix complete. Triggering UniProt annotation.")
          
          # Show modal with progress bar
          shiny::showModal(shiny::modalDialog(
            title = "Retrieving UniProt Annotations",
            shiny::p("Please wait while we retrieve protein annotations from UniProt..."),
            shiny::div(id = "uniprot_progress_text", "Initializing..."),
            shiny::tags$div(class = "progress",
              shiny::tags$div(id = "uniprot_progress_bar", 
                       class = "progress-bar progress-bar-striped active",
                       role = "progressbar",
                       style = "width: 0%",
                       "0%")
            ),
            footer = NULL,
            easyClose = FALSE
          ))
          
          tryCatch({
            protein_column <- workflow_data$column_mapping$protein_col
            if (is.null(protein_column)) {
              stop("Protein column not found in column_mapping. Cannot get annotations.")
            }
            
            cache_dir <- if (!is.null(experiment_paths) && !is.null(experiment_paths$results_dir)) {
              file.path(experiment_paths$results_dir, "cache")
            } else {
              file.path(tempdir(), "proteomics_cache")
            }
            
            uniprot_cache_dir <- file.path(cache_dir, "uniprot_annotations")
            if (!dir.exists(uniprot_cache_dir)) {
              dir.create(uniprot_cache_dir, recursive = TRUE)
            }
            
            # Create progress callback function with error handling
            progress_updater <- function(current, total) {
              tryCatch({
                percent <- round((current / total) * 100)
                session$sendCustomMessage("updateUniprotProgress", list(
                  percent = percent,
                  text = sprintf("Processing chunk %d of %d (%d%%)", current, total, percent)
                ))
              }, error = function(e) {
                # Silently catch any errors to prevent disrupting the main process
                message(paste("Progress update failed:", e$message))
              })
            }

            uniprot_dat_cln <- getUniprotAnnotationsFull(
              data_tbl = workflow_data$data_cln,
              protein_id_column = protein_column,
              cache_dir = uniprot_cache_dir,
              taxon_id = workflow_data$taxon_id,
              progress_callback = progress_updater
            )
            
            workflow_data$uniprot_dat_cln <- uniprot_dat_cln
            assign("uniprot_dat_cln", uniprot_dat_cln, envir = .GlobalEnv)
            
            if (!is.null(experiment_paths) && !is.null(experiment_paths$source_dir)) {
              scripts_uniprot_path <- file.path(experiment_paths$source_dir, "uniprot_dat_cln.RDS")
              saveRDS(uniprot_dat_cln, scripts_uniprot_path)
              log_info(sprintf("Saved uniprot_dat_cln to scripts directory: %s", scripts_uniprot_path))
            }
            log_info(sprintf("UniProt annotations retrieved successfully. Found %d annotations", nrow(uniprot_dat_cln)))
            shiny::removeModal()
            shiny::showNotification("UniProt annotations retrieved successfully.", type = "message")
            
          }, error = function(e) {
            log_warn(paste("Error getting UniProt annotations:", e$message))
            workflow_data$uniprot_dat_cln <- NULL
            shiny::removeModal()
            shiny::showNotification(paste("Warning: Could not retrieve UniProt annotations:", e$message), type = "warning", duration = 8)
          })
          
          # 4. Set the QC trigger to TRUE to signal that QC modules can now run
          message("=== DEBUG66: designMatrixApplet Save Design - About to set qc_trigger ===")
          message(sprintf("   DEBUG66: qc_trigger is NULL = %s", is.null(qc_trigger)))
          if (!is.null(qc_trigger)) {
            message("   DEBUG66: Setting qc_trigger(TRUE)")
            qc_trigger(TRUE)
            message(sprintf("   DEBUG66: qc_trigger set. Current value = %s", qc_trigger()))
          } else {
            message("   DEBUG66: qc_trigger is NULL, cannot set")
          }
          message("=== DEBUG66: designMatrixApplet - qc_trigger setting complete ===")
          
          # 5. Update tab status to 'complete'
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