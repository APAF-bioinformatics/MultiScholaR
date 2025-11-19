# Session Summary and Report Generation Module
# Implements the final workflow steps: save parameters, copy files, generate report

#' @title sessionSummaryModule
#'
#' @description A Shiny module for the Session Summary and Report Generation step.
#'
#' @name sessionSummaryModule
#' @export
NULL

#' @rdname sessionSummaryModule
#' @export
#' @import shiny
#' @import shinydashboard
mod_prot_summary_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::fluidPage(
    shiny::h3("Session Summary & Report Generation"),
    
    shiny::fluidRow(
      # Left column: Workflow Parameters
      shiny::column(6,
        shiny::wellPanel(
          shiny::h4("Workflow Parameters"),
          shiny::textInput(ns("experiment_label"), "Experiment Label:", 
                   value = "", placeholder = "e.g., my_proteomics_analysis"),
          shiny::textAreaInput(ns("description"), "Description:", 
                       value = "Full protein analysis workflow with config parameters",
                       rows = 3, resize = "vertical"),
          shiny::br(),
          shiny::actionButton(ns("save_workflow_args"), "Save Workflow Arguments", 
                      class = "btn-primary", icon = shiny::icon("save"))
        )
      ),
      
      # Right column: File Management
      shiny::column(6,
        shiny::wellPanel(
          shiny::h4("File Management"),
          shiny::br(),
          shiny::actionButton(ns("copy_to_publication"), "Copy to Publication Directory", 
                      class = "btn-info", icon = shiny::icon("copy")),
          shiny::br(), shiny::br(),
          shiny::verbatimTextOutput(ns("copy_status"))
        )
      )
    ),
    
    shiny::fluidRow(
      # Left column: Report Generation
      shiny::column(6,
        shiny::wellPanel(
          shiny::h4("Report Generation"),
          shiny::textOutput(ns("template_status")),
          shiny::br(),
          shiny::actionButton(ns("generate_report"), "Generate Report", 
                      class = "btn-success", icon = shiny::icon("file-pdf")),
          shiny::br(), shiny::br(),
          shiny::conditionalPanel(
            condition = paste0("output['", ns("report_ready"), "']"),
            shiny::downloadButton(ns("download_report"), "Download Report", 
                          class = "btn-success")
          )
        )
      ),
      
      # Right column: GitHub Integration (Optional)
      shiny::column(6,
        shiny::wellPanel(
          shiny::h4("Version Control (Optional)"),
          shiny::checkboxInput(ns("enable_github"), "Enable GitHub Push", FALSE),
          shiny::conditionalPanel(
            condition = paste0("input['", ns("enable_github"), "']"),
            shiny::textInput(ns("github_org"), "GitHub Organization:", ""),
            shiny::textInput(ns("github_email"), "GitHub Email:", ""),
            shiny::textInput(ns("github_username"), "GitHub Username:", ""),
            shiny::textInput(ns("project_id"), "Project ID:", ""),
            shiny::br(),
            shiny::actionButton(ns("push_to_github"), "Push to GitHub", 
                        class = "btn-warning", icon = shiny::icon("github"))
          )
        )
      )
    ),
    
    shiny::fluidRow(
      shiny::column(12,
        shiny::wellPanel(
          shiny::h4("Session Summary"),
          shiny::verbatimTextOutput(ns("session_summary")),
          shiny::br(),
          shiny::actionButton(ns("export_session_state"), "Export Session State (.RDS)", 
                      class = "btn-secondary", icon = shiny::icon("download"))
        )
      )
    )
  )
}

#' @rdname sessionSummaryModule
#' @export
#' @import shiny
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req renderText showNotification downloadHandler withProgress incProgress
#' @importFrom logger log_info log_error
#' @importFrom utils write.table
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
    
    # Template status check
    output$template_status <- shiny::renderText({
      shiny::req(project_dirs)
      
      # Use just omic_type as key (not omic_type_experiment_label)
      if (!omic_type %in% names(project_dirs)) {
        return("⚠️ Project directories not available")
      }
      
      # NEW: Check for both template types and show which is available
      template_dir <- file.path(
        project_dirs[[omic_type]]$base_dir, 
        "scripts", 
        omic_type
      )
      
      diann_template <- file.path(template_dir, "DIANN_report.rmd")
      tmt_template <- file.path(template_dir, "TMT_report.rmd")
      
      templates_status <- c()
      if (file.exists(diann_template)) {
        templates_status <- c(templates_status, "DIA-NN ✅")
      }
      if (file.exists(tmt_template)) {
        templates_status <- c(templates_status, "TMT ✅")
      }
      
      if (length(templates_status) > 0) {
        paste("Templates:", paste(templates_status, collapse = ", "))
      } else {
        "⚠️ Report templates will be downloaded when generating report"
      }
    })
    
    # Save workflow arguments - EXTRACT FROM FINAL S4 OBJECT
    shiny::observeEvent(input$save_workflow_args, {
      shiny::req(input$experiment_label)
      
      cat("SESSION SUMMARY: Starting workflow args save process\n")
      
      tryCatch({
        # Get the DATA S4 object (not enrichment results) from R6 state manager
        final_s4_object <- NULL
        if (!is.null(workflow_data$state_manager)) {
          # Get the data object from correlation_filtered state, NOT the current enrichment state
          data_states <- c("correlation_filtered", "ruv_corrected", "protein_replicate_filtered", "imputed")
          
          # Use functional approach to find first valid state
          available_states <- workflow_data$state_manager$getHistory()
          data_state_used <- purrr::detect(data_states, ~ .x %in% available_states)
          
          if (!is.null(data_state_used)) {
            final_s4_object <- workflow_data$state_manager$getState(data_state_used)
          }
          
          if (!is.null(final_s4_object)) {
            cat(sprintf("SESSION SUMMARY: Retrieved DATA S4 object from state '%s'\n", data_state_used))
            cat(sprintf("SESSION SUMMARY: S4 object class: %s\n", class(final_s4_object)))
            
            # Check if this S4 object has @args slot
            args_available <- tryCatch({
              !is.null(final_s4_object@args)
            }, error = function(e) {
              FALSE
            })
            
            if (args_available) {
              cat(sprintf("SESSION SUMMARY: S4 @args contains %d function groups\n", length(final_s4_object@args)))
            } else {
              cat("SESSION SUMMARY: S4 @args is NULL or slot doesn't exist\n")
            }
          } else {
            cat("SESSION SUMMARY: No data S4 object found in any valid state\n")
          }
        } else {
          cat("SESSION SUMMARY: No state manager available\n")
        }
        
        # Prepare contrasts_tbl if available
        contrasts_tbl <- NULL
        if (!is.null(workflow_data) && !is.null(workflow_data$contrasts_tbl)) {
          contrasts_tbl <- workflow_data$contrasts_tbl
          cat("SESSION SUMMARY: Using contrasts_tbl from workflow_data\n")
        }
        
        # Ensure config_list is in global environment for fallback
        if (!is.null(workflow_data) && !is.null(workflow_data$config_list)) {
          assign("config_list", workflow_data$config_list, envir = .GlobalEnv)
          cat("SESSION SUMMARY: Config list available with", length(workflow_data$config_list), "items\n")
        }
        
        # Call the updated function with S4 object and workflow_data
        cat("SESSION SUMMARY: Creating study_parameters.txt file with S4 parameters and RUV results\n")
        study_params_file <- createWorkflowArgsFromConfig(
          workflow_name = input$experiment_label,
          description = input$description,
          source_dir_path = project_dirs[[omic_type]]$source_dir,
          final_s4_object = final_s4_object,
          contrasts_tbl = contrasts_tbl,
          workflow_data = workflow_data
        )
        
        cat("SESSION SUMMARY: Successfully created study_parameters.txt at:", study_params_file, "\n")
        
        values$workflow_args_saved <- TRUE
        shiny::showNotification("Study parameters saved successfully", type = "message")
        
        output$session_summary <- shiny::renderText({
          paste("Study parameters created for:", input$experiment_label, 
                "\nDescription:", input$description,
                "\nTimestamp:", Sys.time(),
                "\nFile:", study_params_file,
                "\nSource: Final S4 object @args + config_list",
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
        shiny::showNotification("Study parameters saved with warnings", type = "warning")
      })
    })
    
    # Copy files to publication directory
    shiny::observeEvent(input$copy_to_publication, {
      shiny::req(input$experiment_label)
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
            logger::log_error(logger::skip_formatter(paste("Failed to create basic study_parameters.txt:", e$message)))
          })
        }
      }
      
      cat("   SESSION SUMMARY: Copy to Publication button clicked\n")
      cat("   SESSION SUMMARY: workflow_args_saved =", values$workflow_args_saved, "\n")
      cat("   SESSION SUMMARY: experiment_label =", input$experiment_label, "\n")
      
      shiny::withProgress(message = "Copying files to publication directory...", {
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
          
          # Debug: Log project_dirs structure
          cat("   SESSION SUMMARY: project_dirs keys:", paste(names(project_dirs), collapse = ", "), "\n")
          cat("   SESSION SUMMARY: Using omic_type =", omic_type, "experiment_label =", input$experiment_label, "\n")
          
          # Call copyToResultsSummary with explicit parameters
          copyToResultsSummary(
            omic_type = omic_type,
            experiment_label = input$experiment_label,
            contrasts_tbl = contrasts_tbl,
            design_matrix = design_matrix,
            force = TRUE  # Always force in Shiny - user can't interact with terminal prompts
          )
          
          values$files_copied <- TRUE
          
          output$copy_status <- shiny::renderText("Files copied to publication directory successfully ✅")
          shiny::showNotification("Publication files copied", type = "message")
          
          # Update session summary
          output$session_summary <- shiny::renderText({
            paste("Workflow args created for:", input$experiment_label,
                  "\nDescription:", input$description,
                  "\nTimestamp:", Sys.time(),
                  "\nStatus: Arguments saved ✅, Files copied ✅")
          })
          
        }, error = function(e) {
          output$copy_status <- shiny::renderText(paste("Error:", e$message))
          shiny::showNotification(paste("Copy error:", e$message), type = "error", duration = 10)
          # ✅ FIX: Use skip_formatter to avoid glue interpolation issues
          logger::log_error(logger::skip_formatter(paste("Failed to copy files:", e$message)))
          cat("   SESSION SUMMARY ERROR:", e$message, "\n")
          # Print traceback for debugging
          cat("   SESSION SUMMARY Traceback:\n")
          traceback()
        })
      })
    }, ignoreInit = TRUE)  # ✅ FIX: Add ignoreInit to prevent triggering on module load
    
    # Generate report with template download logic
    shiny::observeEvent(input$generate_report, {
      shiny::req(input$experiment_label)
      shiny::req(values$files_copied)
      
      # Validate project_dirs structure
      if (!omic_type %in% names(project_dirs) || is.null(project_dirs[[omic_type]]$base_dir)) {
        shiny::showNotification("Error: Project directories not properly initialized", 
                        type = "error", duration = 10)
        return()
      }
      
      shiny::withProgress(message = "Generating report...", {
        
        # STEP 1: Detect workflow type and download appropriate template
        shiny::incProgress(0.1, detail = "Detecting workflow type...")
        
        # Try to detect workflow type from multiple sources (in order of preference)
        workflow_type_detected <- NULL
        
        # 1. Check workflow_data first (most reliable)
        if (!is.null(workflow_data) && !is.null(workflow_data$config_list)) {
          workflow_type_detected <- workflow_data$config_list$globalParameters$workflow_type
          cat(sprintf("   REPORT: Detected workflow_type from workflow_data: %s\n", workflow_type_detected))
        }
        
        # 2. Fallback: Check S4 object from state manager
        if (is.null(workflow_type_detected) && !is.null(workflow_data$state_manager)) {
          # Get the DATA S4 object (not enrichment results)
          data_states <- c("correlation_filtered", "ruv_corrected", "protein_replicate_filtered", "imputed")
          available_states <- workflow_data$state_manager$getHistory()
          data_state_used <- purrr::detect(data_states, ~ .x %in% available_states)
          
          if (!is.null(data_state_used)) {
            current_s4 <- workflow_data$state_manager$getState(data_state_used)
            if (!is.null(current_s4) && !is.null(current_s4@args$globalParameters$workflow_type)) {
              workflow_type_detected <- current_s4@args$globalParameters$workflow_type
              cat(sprintf("   REPORT: Detected workflow_type from S4 object (state: %s): %s\n", data_state_used, workflow_type_detected))
            }
          }
        }
        
        # 3. Fallback: Check global environment
        if (is.null(workflow_type_detected) && exists("config_list", envir = .GlobalEnv)) {
          config_list <- get("config_list", envir = .GlobalEnv)
          if (!is.null(config_list$globalParameters$workflow_type)) {
            workflow_type_detected <- config_list$globalParameters$workflow_type
            cat(sprintf("   REPORT: Detected workflow_type from global config_list: %s\n", workflow_type_detected))
          }
        }
        
        # 4. Final fallback: Default to DIA for backward compatibility
        if (is.null(workflow_type_detected) || !nzchar(workflow_type_detected)) {
          workflow_type_detected <- "DIA"
          cat("   REPORT: Using default workflow_type: DIA (no workflow type found)\n")
        }
        
        # Determine template filename based on workflow type
        template_filename <- if (tolower(workflow_type_detected) == "tmt" || tolower(workflow_type_detected) == "tmt_pd") {
          "TMT_report.rmd"
        } else if (tolower(workflow_type_detected) == "lfq") {
          "LFQ_report.rmd"
        } else {
          "DIANN_report.rmd"  # Default for DIA and unknown types
        }
        
        cat(sprintf("   REPORT: Selected template: %s for workflow_type: %s\n", template_filename, workflow_type_detected))
        
        # Construct template path
        report_template_path <- file.path(
          project_dirs[[omic_type]]$base_dir, 
          "scripts", 
          omic_type, 
          template_filename
        )
        
        # Download template if it doesn't exist
        if (!file.exists(report_template_path)) {
          shiny::incProgress(0.2, detail = paste("Checking for", template_filename, "template..."))
          
          # Create directories if needed
          dir.create(dirname(report_template_path), recursive = TRUE, showWarnings = FALSE)
          
          tryCatch({
            # 1. Try to get from package installation first (system.file)
            pkg_file <- system.file("reports", "proteomics", template_filename, package = "MultiScholaR")
            
            if (file.exists(pkg_file) && pkg_file != "") {
              cat(sprintf("   REPORT: Found template in package at: %s\n", pkg_file))
              file.copy(pkg_file, report_template_path)
              logger::log_info("Template copied from package to: {report_template_path}")
              shiny::showNotification(paste(template_filename, "template copied from package"), type = "message")
            } else {
              # 2. Fallback to GitHub download
              cat(sprintf("   REPORT: Template not found in package, downloading from GitHub...\n"))
              
              # Construct URL based on template filename (using main branch and inst structure)
              template_url <- paste0(
                "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/reports/proteomics/",
                template_filename
              )
              
              cat(sprintf("   REPORT: Downloading template from: %s\n", template_url))
              download.file(template_url, destfile = report_template_path, quiet = TRUE)
              logger::log_info("Template downloaded to: {report_template_path}")
              shiny::showNotification(paste(template_filename, "template downloaded"), type = "message")
            }
          }, error = function(e) {
            shiny::showNotification(paste("Template retrieval failed:", e$message), 
                            type = "error", duration = 10)
            logger::log_error(logger::skip_formatter(paste("Failed to retrieve template:", e$message)))
            return()
          })
        } else {
          cat(sprintf("   REPORT: Template already exists at: %s\n", report_template_path))
        }
        
        # STEP 2: Generate report
        shiny::incProgress(0.5, detail = "Rendering report...")
        
        tryCatch({
          # Check if RenderReport function exists
          if (!exists("RenderReport")) {
            shiny::showNotification("Error: RenderReport function not found. Please ensure MultiScholaR is properly loaded.", 
                            type = "error", duration = 15)
            return()
          }
          
          logger::log_info("Calling RenderReport with omic_type: {omic_type}, experiment_label: {input$experiment_label}")
          
          rendered_path <- RenderReport(
            omic_type = omic_type,
            experiment_label = input$experiment_label,
            rmd_filename = template_filename  # Use detected template filename instead of hardcoded "DIANN_report.rmd"
          )
          
          logger::log_info("RenderReport returned path: {rendered_path}")
          
          if (!is.null(rendered_path) && file.exists(rendered_path)) {
            values$report_generated <- TRUE
            values$report_path <- rendered_path
            
            # Enable download button
            output$report_ready <- shiny::reactive({ TRUE })
            shiny::outputOptions(output, "report_ready", suspendWhenHidden = FALSE)
            
            # Setup download handler
            output$download_report <- shiny::downloadHandler(
              filename = function() basename(rendered_path),
              content = function(file) file.copy(rendered_path, file)
            )
            
            shiny::showNotification("Report generated successfully!", type = "message")
            
            # Update session summary
            output$session_summary <- shiny::renderText({
              paste("Workflow args created for:", input$experiment_label,
                    "\nDescription:", input$description,
                    "\nTimestamp:", Sys.time(),
                    "\nStatus: Arguments saved ✅, Files copied ✅, Report generated ✅",
                    "\nReport location:", rendered_path)
            })
            
          } else {
            shiny::showNotification("Report generation failed - no output file created", 
                            type = "error", duration = 10)
          }
          
        }, error = function(e) {
          error_msg <- paste("Report generation failed:", e$message)
          
          # Enhanced debug output
          cat("   DEBUG66: REPORT GENERATION ERROR\n")
          cat(sprintf("      Error class: %s\n", class(e)[1]))
          cat(sprintf("      Error message: %s\n", e$message))
          cat("      Full error object:\n")
          print(e)
          
          logger::log_error("Failed to generate report: {e$message}")
          logger::log_error("Error class: {class(e)[1]}")
          
          shiny::showNotification(error_msg, type = "error", duration = 15)
          
          # Show detailed error in console for user
          shiny::showNotification(
            paste("Debug info: Check R console for detailed error trace"), 
            type = "warning", duration = 10
          )
        })
      })
    })
    
    # GitHub integration
    shiny::observeEvent(input$push_to_github, {
      shiny::req(input$enable_github, input$github_org, input$github_email, 
          input$github_username, input$project_id)
      shiny::req(values$report_generated)
      
      shiny::withProgress(message = "Pushing to GitHub...", {
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
          
          shiny::showNotification("Successfully pushed to GitHub", type = "message")
          
          # Update session summary
          output$session_summary <- shiny::renderText({
            paste("Workflow args created for:", input$experiment_label,
                  "\nDescription:", input$description,
                  "\nTimestamp:", Sys.time(),
                  "\nStatus: Arguments saved ✅, Files copied ✅, Report generated ✅, GitHub pushed ✅")
          })
          
        }, error = function(e) {
          shiny::showNotification(paste("GitHub push failed:", e$message), 
                          type = "error", duration = 10)
          logger::log_error(logger::skip_formatter(paste("Failed to push to GitHub:", e$message)))
        })
      })
    })
    
    # Export session state
    shiny::observeEvent(input$export_session_state, {
      shiny::req(input$experiment_label)
      
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
        
        shiny::showNotification(paste("Session state exported to:", session_export_path), 
                        type = "message")
        logger::log_info(paste("Session state exported to:", session_export_path))
        
      }, error = function(e) {
        shiny::showNotification(paste("Export failed:", e$message), type = "error", duration = 10)
        logger::log_error(logger::skip_formatter(paste("Failed to export session state:", e$message)))
      })
    })
    
    # Initialize session summary
    output$session_summary <- shiny::renderText({
      "Ready to save workflow parameters and generate report"
    })
    
    # Initialize report_ready as FALSE
    output$report_ready <- shiny::reactive({ FALSE })
    shiny::outputOptions(output, "report_ready", suspendWhenHidden = FALSE)
    
  })
}

