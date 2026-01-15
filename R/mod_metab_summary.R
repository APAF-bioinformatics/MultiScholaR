# ============================================================================
# mod_metab_summary.R
# ============================================================================
# Purpose: Metabolomics session summary and integration export Shiny module
#
# This module provides a summary of the processing session and enables
# export of results for multi-omics integration.
# ============================================================================

#' @title Metabolomics Session Summary Module
#' @description A Shiny module for displaying processing summary and enabling
#'              data export for integration workflows.
#' @name mod_metab_summary
NULL

#' @rdname mod_metab_summary
#' @export
#' @importFrom shiny NS tagList fluidPage h3 fluidRow column wellPanel h4 textInput textAreaInput actionButton icon br hr downloadButton verbatimTextOutput checkboxInput conditionalPanel textOutput
mod_metab_summary_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::fluidPage(
            shiny::h3("Session Summary & Report Generation"),

            shiny::fluidRow(
                # Left column: Workflow Parameters
                shiny::column(6,
                    shiny::wellPanel(
                        shiny::h4("Workflow Parameters"),
                        shiny::textInput(ns("experiment_label"), "Experiment Label:",
                                        value = "", placeholder = "e.g., my_metabolomics_analysis"),
                        shiny::textAreaInput(ns("description"), "Description:",
                                            value = "Full metabolomics analysis workflow with normalization and DE",
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
    )
}

#' @rdname mod_metab_summary
#' @export
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req renderText showNotification downloadHandler withProgress incProgress
#' @importFrom logger log_info log_error
#' @importFrom utils write.table
#' @importFrom purrr detect
#' @importFrom readr read_tsv
mod_metab_summary_server <- function(id, project_dirs, omic_type = "metabolomics", experiment_label = NULL, workflow_data = NULL) {
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

            # Use just omic_type as key
            if (!omic_type %in% names(project_dirs)) {
                return("⚠️ Project directories not available")
            }

            # Check for metabolomics report template
            template_dir <- file.path(
                project_dirs[[omic_type]]$base_dir,
                "scripts",
                omic_type
            )

            metab_template <- file.path(template_dir, "metabolomics_report.rmd")

            if (file.exists(metab_template)) {
                "Template: Metabolomics Report ✅"
            } else {
                "⚠️ Report template will be downloaded when generating report"
            }
        })

        # Save workflow arguments - EXTRACT FROM FINAL S4 OBJECT
        shiny::observeEvent(input$save_workflow_args, {
            shiny::req(input$experiment_label)

            cat("SESSION SUMMARY: Starting workflow args save process\n")

            tryCatch({
                # Get the DATA S4 object from R6 state manager
                final_s4_object <- NULL
                if (!is.null(workflow_data$state_manager)) {
                    # Get the data object from latest valid state
                    data_states <- c("metab_correlation_filtered", "metab_norm_complete",
                                   "metab_ruv_corrected", "metab_normalized",
                                   "loaded_for_de", "metab_qc_complete")

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
                cat("SESSION SUMMARY: Creating study_parameters.txt file with S4 parameters\n")
                study_params_file <- createWorkflowArgsFromConfig(
                    workflow_name = input$experiment_label,
                    description = input$description,
                    source_dir_path = project_dirs[[omic_type]]$source_dir,
                    final_s4_object = final_s4_object,
                    contrasts_tbl = contrasts_tbl,
                    workflow_data = workflow_data
                )

                cat("SESSION SUMMARY: Successfully created study_parameters.txt at:", study_params_file, "\n")

                # --- SAVE INTEGRATION OBJECT ---
                cat("SESSION SUMMARY: Saving Integration S4 Object...\n")
                integration_dir <- project_dirs[[omic_type]]$integration_dir
                # Fallback if integration_dir not defined in project_dirs
                if (is.null(integration_dir)) {
                    integration_dir <- file.path(project_dirs[[omic_type]]$base_dir, "integration")
                }

                if (!dir.exists(integration_dir)) {
                    dir.create(integration_dir, recursive = TRUE, showWarnings = FALSE)
                }

                # Use standardized naming: [OmicType]_[ExperimentLabel]_final_s4.RDS
                s4_filename <- sprintf("%s_%s_final_s4.RDS", omic_type, input$experiment_label)
                s4_filepath <- file.path(integration_dir, s4_filename)

                saveRDS(final_s4_object, s4_filepath)
                cat(sprintf("SESSION SUMMARY: Saved Integration S4 object to: %s\n", s4_filepath))
                shiny::showNotification("Saved Integration S4 Object", type = "message")

                values$workflow_args_saved <- TRUE
                shiny::showNotification("Study parameters saved successfully", type = "message")

                output$session_summary <- shiny::renderText({
                    paste("Study parameters created for:", input$experiment_label,
                          "\nDescription:", input$description,
                          "\nTimestamp:", Sys.time(),
                          "\nFile:", study_params_file,
                          "\nSource: Final S4 object @args + config_list",
                          "\nIntegration Object:", s4_filename,
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

            # Make requirement flexible - allow copy even if workflow args had issues
            if (!values$workflow_args_saved) {
                basic_params_file <- file.path(project_dirs[[omic_type]]$source_dir, "study_parameters.txt")
                if (!file.exists(basic_params_file)) {
                    cat("SESSION SUMMARY: Creating basic study_parameters.txt as fallback\n")
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

            cat("SESSION SUMMARY: Copy to Publication button clicked\n")
            cat("SESSION SUMMARY: workflow_args_saved =", values$workflow_args_saved, "\n")
            cat("SESSION SUMMARY: experiment_label =", input$experiment_label, "\n")

            shiny::withProgress(message = "Copying files to publication directory...", {
                tryCatch({
                    # Get contrasts_tbl and design_matrix from workflow_data or files
                    contrasts_tbl <- NULL
                    design_matrix <- NULL

                    if (!is.null(workflow_data)) {
                        if (!is.null(workflow_data$contrasts_tbl)) {
                            contrasts_tbl <- workflow_data$contrasts_tbl
                            cat("SESSION SUMMARY: Got contrasts_tbl from workflow_data\n")
                        }
                        if (!is.null(workflow_data$design_matrix)) {
                            design_matrix <- workflow_data$design_matrix
                            cat("SESSION SUMMARY: Got design_matrix from workflow_data\n")
                        }
                    }

                    # Fallback: Try to load from files if not in workflow_data
                    if (is.null(design_matrix)) {
                        design_matrix_file <- file.path(project_dirs[[omic_type]]$source_dir, "design_matrix.tab")
                        if (file.exists(design_matrix_file)) {
                            design_matrix <- readr::read_tsv(design_matrix_file, show_col_types = FALSE)
                            cat("SESSION SUMMARY: Loaded design_matrix from file\n")
                        }
                    }

                    if (is.null(contrasts_tbl)) {
                        contrasts_file <- file.path(project_dirs[[omic_type]]$source_dir, "contrasts_tbl.tab")
                        if (file.exists(contrasts_file)) {
                            contrasts_tbl <- readr::read_tsv(contrasts_file, show_col_types = FALSE)
                            cat("SESSION SUMMARY: Loaded contrasts_tbl from file\n")
                        }
                    }

                    cat("SESSION SUMMARY: About to call copyToResultsSummary\n")
                    cat("SESSION SUMMARY: omic_type =", omic_type, "\n")
                    cat("SESSION SUMMARY: experiment_label =", input$experiment_label, "\n")
                    cat("SESSION SUMMARY: contrasts_tbl is", ifelse(is.null(contrasts_tbl), "NULL", "available"), "\n")
                    cat("SESSION SUMMARY: design_matrix is", ifelse(is.null(design_matrix), "NULL", "available"), "\n")

                    # Ensure project_dirs is in global environment for copyToResultsSummary
                    if (!exists("project_dirs", envir = .GlobalEnv)) {
                        cat("SESSION SUMMARY: Setting project_dirs in global environment\n")
                        assign("project_dirs", project_dirs, envir = .GlobalEnv)
                    }

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
                    logger::log_error(logger::skip_formatter(paste("Failed to copy files:", e$message)))
                    cat("SESSION SUMMARY ERROR:", e$message, "\n")
                    cat("SESSION SUMMARY Traceback:\n")
                    traceback()
                })
            })
        }, ignoreInit = TRUE)

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

                # Determine template filename for metabolomics
                template_filename <- "metabolomics_report.rmd"

                cat(sprintf("REPORT: Selected template: %s\n", template_filename))

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
                        # Try to get from package installation first
                        pkg_file <- system.file("reports", "metabolomics", template_filename, package = "MultiScholaR")

                        if (file.exists(pkg_file) && pkg_file != "") {
                            cat(sprintf("REPORT: Found template in package at: %s\n", pkg_file))
                            file.copy(pkg_file, report_template_path)
                            logger::log_info("Template copied from package to: {report_template_path}")
                            shiny::showNotification(paste(template_filename, "template copied from package"), type = "message")
                        } else {
                            # Fallback to GitHub download
                            cat(sprintf("REPORT: Template not found in package, downloading from GitHub...\n"))

                            template_url <- paste0(
                                "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/reports/metabolomics/",
                                template_filename
                            )

                            cat(sprintf("REPORT: Downloading template from: %s\n", template_url))
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
                    cat(sprintf("REPORT: Template already exists at: %s\n", report_template_path))
                }

                # Generate report
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
                        rmd_filename = template_filename
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

                    cat("DEBUG66: REPORT GENERATION ERROR\n")
                    cat(sprintf("Error class: %s\n", class(e)[1]))
                    cat(sprintf("Error message: %s\n", e$message))
                    cat("Full error object:\n")
                    print(e)

                    logger::log_error("Failed to generate report: {e$message}")
                    logger::log_error("Error class: {class(e)[1]}")

                    shiny::showNotification(error_msg, type = "error", duration = 15)

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

