# Session Summary and Report Generation Module
# Implements the final workflow steps: save parameters, copy files, generate report

sessionSummaryUI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    h3("Session Summary & Report Generation"),
    
    fluidRow(
      # Left column: Workflow Parameters
      column(6,
        wellPanel(
          h4("Workflow Parameters"),
          textInput(ns("experiment_label"), "Experiment Label:", 
                   value = "", placeholder = "e.g., my_proteomics_analysis"),
          textAreaInput(ns("description"), "Description:", 
                       value = "Full protein analysis workflow with config parameters",
                       rows = 3, resize = "vertical"),
          br(),
          actionButton(ns("save_workflow_args"), "Save Workflow Arguments", 
                      class = "btn-primary", icon = icon("save"))
        )
      ),
      
      # Right column: File Management
      column(6,
        wellPanel(
          h4("File Management"),
          checkboxInput(ns("force_copy_publication"), 
                       "Force Overwrite Publication Files", FALSE),
          br(),
          actionButton(ns("copy_to_publication"), "Copy to Publication Directory", 
                      class = "btn-info", icon = icon("copy")),
          br(), br(),
          verbatimTextOutput(ns("copy_status"))
        )
      )
    ),
    
    fluidRow(
      # Left column: Report Generation
      column(6,
        wellPanel(
          h4("Report Generation"),
          textOutput(ns("template_status")),
          br(),
          actionButton(ns("generate_report"), "Generate Report", 
                      class = "btn-success", icon = icon("file-pdf")),
          br(), br(),
          conditionalPanel(
            condition = paste0("output['", ns("report_ready"), "']"),
            downloadButton(ns("download_report"), "Download Report", 
                          class = "btn-success")
          )
        )
      ),
      
      # Right column: GitHub Integration (Optional)
      column(6,
        wellPanel(
          h4("Version Control (Optional)"),
          checkboxInput(ns("enable_github"), "Enable GitHub Push", FALSE),
          conditionalPanel(
            condition = paste0("input['", ns("enable_github"), "']"),
            textInput(ns("github_org"), "GitHub Organization:", ""),
            textInput(ns("github_email"), "GitHub Email:", ""),
            textInput(ns("github_username"), "GitHub Username:", ""),
            textInput(ns("project_id"), "Project ID:", ""),
            br(),
            actionButton(ns("push_to_github"), "Push to GitHub", 
                        class = "btn-warning", icon = icon("github"))
          )
        )
      )
    ),
    
    fluidRow(
      column(12,
        wellPanel(
          h4("Session Summary"),
          verbatimTextOutput(ns("session_summary")),
          br(),
          actionButton(ns("export_session_state"), "Export Session State (.RDS)", 
                      class = "btn-secondary", icon = icon("download"))
        )
      )
    )
  )
}

sessionSummaryServer <- function(id, project_dirs, omic_type = "proteomics", experiment_label = NULL, workflow_data = NULL) {
  moduleServer(id, function(input, output, session) {
    
    # Auto-populate experiment label if provided
    if (!is.null(experiment_label)) {
      updateTextInput(session, "experiment_label", value = experiment_label)
    }
    
    # Reactive values for tracking state
    values <- reactiveValues(
      workflow_args_saved = FALSE,
      files_copied = FALSE,
      report_generated = FALSE,
      report_path = NULL
    )
    
    # Template status check
    output$template_status <- renderText({
      req(project_dirs)
      
      # Use just omic_type as key (not omic_type_experiment_label)
      if (!omic_type %in% names(project_dirs)) {
        return("⚠️ Project directories not available")
      }
      
      report_path <- file.path(
        project_dirs[[omic_type]]$base_dir, 
        "scripts", 
        omic_type, 
        "DIANN_report.rmd"
      )
      
      if (file.exists(report_path)) {
        "✅ Report template available"
      } else {
        "⚠️ Report template needs to be downloaded"
      }
    })
    
    # Save workflow arguments - NEW SIMPLE APPROACH
    observeEvent(input$save_workflow_args, {
      req(input$experiment_label)
      
      cat("SESSION SUMMARY: Starting workflow args save process\n")
      
      # Simple, direct approach - no S4 objects!
      tryCatch({
        # Prepare contrasts_tbl if available
        contrasts_tbl <- NULL
        if (!is.null(workflow_data) && !is.null(workflow_data$contrasts_tbl)) {
          contrasts_tbl <- workflow_data$contrasts_tbl
          cat("SESSION SUMMARY: Using contrasts_tbl from workflow_data\n")
        }
        
        # Ensure config_list is in global environment
        if (!is.null(workflow_data) && !is.null(workflow_data$config_list)) {
          assign("config_list", workflow_data$config_list, envir = .GlobalEnv)
          cat("SESSION SUMMARY: Config list available with", length(workflow_data$config_list), "items\n")
        }
        
        # Call the new simple function directly
        cat("SESSION SUMMARY: Creating study_parameters.txt file\n")
        study_params_file <- createStudyParametersFile(
          workflow_name = input$experiment_label,
          description = input$description,
          source_dir_path = project_dirs[[omic_type]]$source_dir,
          contrasts_tbl = contrasts_tbl
        )
        
        cat("SESSION SUMMARY: Successfully created study_parameters.txt at:", study_params_file, "\n")
        
        values$workflow_args_saved <- TRUE
        showNotification("Study parameters saved successfully", type = "success")
        
        output$session_summary <- renderText({
          paste("Study parameters created for:", input$experiment_label, 
                "\nDescription:", input$description,
                "\nTimestamp:", Sys.time(),
                "\nFile:", study_params_file,
                "\nStatus: Parameters saved ✅")
        })
        
      }, error = function(e) {
        cat("SESSION SUMMARY ERROR:", e$message, "\n")
        
        # Create basic file as fallback
        basic_params_file <- file.path(project_dirs[[omic_type]]$source_dir, "study_parameters.txt")
        if (!file.exists(basic_params_file)) {
          writeLines(c(
            "Study Parameters",
            "================",
            "",
            paste("Workflow Name:", input$experiment_label),
            paste("Description:", input$description),
            paste("Timestamp:", Sys.time()),
            paste("Error:", e$message)
          ), basic_params_file)
          cat("SESSION SUMMARY: Created basic fallback file at:", basic_params_file, "\n")
        }
        
        values$workflow_args_saved <- TRUE
        showNotification("Study parameters saved with warnings", type = "warning")
      })
    })
    
    # Copy files to publication directory
    observeEvent(input$copy_to_publication, {
      req(input$experiment_label)
      # ✅ FIX: Make the requirement more flexible - allow copy even if workflow args had issues
      if (!values$workflow_args_saved) {
        # Try to create a basic study_parameters.txt file if it doesn't exist
        basic_params_file <- file.path(project_dirs[[omic_type]]$source_dir, "study_parameters.txt")
        if (!file.exists(basic_params_file)) {
          cat("   SESSION SUMMARY: Creating basic study_parameters.txt as fallback\n")
          basic_content <- paste(
            "Study Parameters",
            "================",
            "",
            paste("Workflow Name:", input$experiment_label),
            paste("Description:", input$description),
            paste("Timestamp:", Sys.time()),
            paste("Note: Some parameters could not be saved due to serialization issues"),
            sep = "\n"
          )
          tryCatch({
            writeLines(basic_content, basic_params_file)
            values$workflow_args_saved <- TRUE
          }, error = function(e) {
            logger::log_error(paste("Failed to create basic study_parameters.txt:", e$message))
          })
        }
      }
      
      cat("   SESSION SUMMARY: Copy to Publication button clicked\n")
      cat("   SESSION SUMMARY: workflow_args_saved =", values$workflow_args_saved, "\n")
      cat("   SESSION SUMMARY: experiment_label =", input$experiment_label, "\n")
      
      withProgress(message = "Copying files to publication directory...", {
        tryCatch({
          # ✅ FIX: Get contrasts_tbl and design_matrix from workflow_data or files
          contrasts_tbl <- NULL
          design_matrix <- NULL
          
          if (!is.null(workflow_data)) {
            if (!is.null(workflow_data$contrasts_tbl)) {
              contrasts_tbl <- workflow_data$contrasts_tbl
              cat("   SESSION SUMMARY Step: Got contrasts_tbl from workflow_data\n")
            }
            if (!is.null(workflow_data$design_matrix)) {
              design_matrix <- workflow_data$design_matrix
              cat("   SESSION SUMMARY Step: Got design_matrix from workflow_data\n")
            }
          }
          
          # Fallback: Try to load from files if not in workflow_data
          if (is.null(design_matrix)) {
            design_matrix_file <- file.path(project_dirs[[omic_type]]$source_dir, "design_matrix.tab")
            if (file.exists(design_matrix_file)) {
              design_matrix <- readr::read_tsv(design_matrix_file, show_col_types = FALSE)
              cat("   SESSION SUMMARY Step: Loaded design_matrix from file\n")
            }
          }
          
          if (is.null(contrasts_tbl)) {
            contrasts_file <- file.path(project_dirs[[omic_type]]$source_dir, "contrasts_tbl.tab")
            if (file.exists(contrasts_file)) {
              contrasts_tbl <- readr::read_tsv(contrasts_file, show_col_types = FALSE)
              cat("   SESSION SUMMARY Step: Loaded contrasts_tbl from file\n")
            }
          }
          
          # Debug output
          cat("   SESSION SUMMARY: About to call copyToResultsSummary\n")
          cat("   SESSION SUMMARY: omic_type =", omic_type, "\n")
          cat("   SESSION SUMMARY: experiment_label =", input$experiment_label, "\n")
          cat("   SESSION SUMMARY: contrasts_tbl is", ifelse(is.null(contrasts_tbl), "NULL", "available"), "\n")
          cat("   SESSION SUMMARY: design_matrix is", ifelse(is.null(design_matrix), "NULL", "available"), "\n")
          
          # ✅ FIX: Ensure project_dirs is in global environment for copyToResultsSummary
          if (!exists("project_dirs", envir = .GlobalEnv)) {
            cat("   SESSION SUMMARY: Setting project_dirs in global environment\n")
            assign("project_dirs", project_dirs, envir = .GlobalEnv)
          }
          
          # Call copyToResultsSummary with explicit parameters
          copyToResultsSummary(
            omic_type = omic_type,
            experiment_label = input$experiment_label,
            contrasts_tbl = contrasts_tbl,
            design_matrix = design_matrix,
            force = input$force_copy_publication
          )
          
          values$files_copied <- TRUE
          
          output$copy_status <- renderText("Files copied to publication directory successfully ✅")
          showNotification("Publication files copied", type = "success")
          
          # Update session summary
          output$session_summary <- renderText({
            paste("Workflow args created for:", input$experiment_label,
                  "\nDescription:", input$description,
                  "\nTimestamp:", Sys.time(),
                  "\nStatus: Arguments saved ✅, Files copied ✅")
          })
          
        }, error = function(e) {
          output$copy_status <- renderText(paste("Error:", e$message))
          showNotification(paste("Copy error:", e$message), type = "error", duration = 10)
          # ✅ FIX: Use paste() instead of logger interpolation in error handler
          logger::log_error(paste("Failed to copy files:", e$message))
          cat("   SESSION SUMMARY ERROR:", e$message, "\n")
          # Print traceback for debugging
          cat("   SESSION SUMMARY Traceback:\n")
          traceback()
        })
      })
    }, ignoreInit = TRUE)  # ✅ FIX: Add ignoreInit to prevent triggering on module load
    
    # Generate report with template download logic
    observeEvent(input$generate_report, {
      req(input$experiment_label)
      req(values$files_copied)
      
      # Validate project_dirs structure
      if (!omic_type %in% names(project_dirs) || is.null(project_dirs[[omic_type]]$base_dir)) {
        showNotification("Error: Project directories not properly initialized", 
                        type = "error", duration = 10)
        return()
      }
      
      withProgress(message = "Generating report...", {
        
        # STEP 1: Ensure template exists (download if needed)
        report_template_path <- file.path(
          project_dirs[[omic_type]]$base_dir, 
          "scripts", 
          omic_type, 
          "DIANN_report.rmd"
        )
        
        if (!file.exists(report_template_path)) {
          incProgress(0.2, detail = "Downloading report template...")
          
          # Create directories if needed
          dir.create(dirname(report_template_path), recursive = TRUE, showWarnings = FALSE)
          
          tryCatch({
            template_url <- "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/Workbooks/proteomics/report/DIANN_report.rmd"
            download.file(template_url, destfile = report_template_path, quiet = TRUE)
            logger::log_info(paste("Template downloaded to:", report_template_path))
            showNotification("Report template downloaded", type = "message")
          }, error = function(e) {
            showNotification(paste("Template download failed:", e$message), 
                            type = "error", duration = 10)
            logger::log_error(paste("Failed to download template:", e$message))
            return()
          })
        }
        
        # STEP 2: Generate report
        incProgress(0.5, detail = "Rendering report...")
        
        tryCatch({
          # Check if RenderReport function exists
          if (!exists("RenderReport")) {
            showNotification("Error: RenderReport function not found. Please ensure MultiScholaR is properly loaded.", 
                            type = "error", duration = 15)
            return()
          }
          
          logger::log_info(paste("Calling RenderReport with omic_type:", omic_type, "experiment_label:", input$experiment_label))
          
          rendered_path <- RenderReport(
            omic_type = omic_type,
            experiment_label = input$experiment_label,
            rmd_filename = "DIANN_report.rmd"
          )
          
          logger::log_info(paste("RenderReport returned path:", rendered_path))
          
          if (!is.null(rendered_path) && file.exists(rendered_path)) {
            values$report_generated <- TRUE
            values$report_path <- rendered_path
            
            # Enable download button
            output$report_ready <- reactive({ TRUE })
            outputOptions(output, "report_ready", suspendWhenHidden = FALSE)
            
            # Setup download handler
            output$download_report <- downloadHandler(
              filename = function() basename(rendered_path),
              content = function(file) file.copy(rendered_path, file)
            )
            
            showNotification("Report generated successfully!", type = "success")
            
            # Update session summary
            output$session_summary <- renderText({
              paste("Workflow args created for:", input$experiment_label,
                    "\nDescription:", input$description,
                    "\nTimestamp:", Sys.time(),
                    "\nStatus: Arguments saved ✅, Files copied ✅, Report generated ✅",
                    "\nReport location:", rendered_path)
            })
            
          } else {
            showNotification("Report generation failed - no output file created", 
                            type = "error", duration = 10)
          }
          
        }, error = function(e) {
          error_msg <- paste("Report generation failed:", e$message)
          logger::log_error(paste("Failed to generate report:", e$message))
          logger::log_error("Error traceback:")
          logger::log_error(paste(capture.output(traceback()), collapse = "\n"))
          
          showNotification(error_msg, type = "error", duration = 15)
          
          # Additional troubleshooting info
          showNotification(
            "Troubleshooting: Check that all workflow steps are complete and data files exist.", 
            type = "warning", duration = 10
          )
        })
      })
    })
    
    # GitHub integration
    observeEvent(input$push_to_github, {
      req(input$enable_github, input$github_org, input$github_email, 
          input$github_username, input$project_id)
      req(values$report_generated)
      
      withProgress(message = "Pushing to GitHub...", {
        tryCatch({
          options(
            github_org = input$github_org,
            github_user_email = input$github_email,
            github_user_name = input$github_username
          )
          
          pushProjectToGithubFromDirs(
            project_dirs = project_dirs,
            omic_type = omic_type, 
            experiment_label = input$experiment_label,
            project_id = input$project_id
          )
          
          showNotification("Successfully pushed to GitHub", type = "success")
          
          # Update session summary
          output$session_summary <- renderText({
            paste("Workflow args created for:", input$experiment_label,
                  "\nDescription:", input$description,
                  "\nTimestamp:", Sys.time(),
                  "\nStatus: Arguments saved ✅, Files copied ✅, Report generated ✅, GitHub pushed ✅")
          })
          
        }, error = function(e) {
          showNotification(paste("GitHub push failed:", e$message), 
                          type = "error", duration = 10)
          logger::log_error(paste("Failed to push to GitHub:", e$message))
        })
      })
    })
    
    # Export session state
    observeEvent(input$export_session_state, {
      req(input$experiment_label)
      
      tryCatch({
        session_export_path <- file.path(
          project_dirs[[omic_type]]$source_dir,
          paste0("session_state_", Sys.Date(), ".RDS")
        )
        
        # Create session state object
        session_state <- list(
          experiment_label = input$experiment_label,
          description = input$description,
          timestamp = Sys.time(),
          omic_type = omic_type,
          workflow_args_saved = values$workflow_args_saved,
          files_copied = values$files_copied,
          report_generated = values$report_generated,
          report_path = values$report_path,
          project_dirs = project_dirs
        )
        
        saveRDS(session_state, session_export_path)
        
        showNotification(paste("Session state exported to:", session_export_path), 
                        type = "success")
        logger::log_info(paste("Session state exported to:", session_export_path))
        
      }, error = function(e) {
        showNotification(paste("Export failed:", e$message), type = "error", duration = 10)
        logger::log_error(paste("Failed to export session state:", e$message))
      })
    })
    
    # Initialize session summary
    output$session_summary <- renderText({
      "Ready to save workflow parameters and generate report"
    })
    
    # Initialize report_ready as FALSE
    output$report_ready <- reactive({ FALSE })
    outputOptions(output, "report_ready", suspendWhenHidden = FALSE)
    
  })
} 