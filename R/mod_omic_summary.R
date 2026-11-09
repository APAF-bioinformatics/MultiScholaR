# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# KIT-108 W3 (D2): the shared Summary & Export tab (kit-101 design 4.10), from the lipid and metab summary modules at
# 5be34ee5 (mod_lipid_summary.R's markup; the observers of mod_lipid_summary_*_helpers.R and
# mod_metab_summary_*_helpers.R, and trx's export panel, mod_trx_summary.R; deleted by W3). `mod_omic_summary_ui(id,
# config)` renders an omic's `summary` group: today's layout (byte-identical, W1 (a)) with the group's placeholder,
# description and test ids, and each export step's panel above the parameters. `mod_omic_summary_server(...)` holds the
# tab's reactive values, today's initial outputs, the label's prefill, the template status, the save, copy, report,
# GitHub and session-export observers with today's `req` gates and F-11's refusal notices, and one panel server per
# export step. Each observer calls one core entry (utils_workflow_summary.R), which notifies through the frame's
# callbacks; the frame renders. The texts a step renders are those of the values the observer passed (the label and
# description at the click). Functions a test mocks are called by name, or read from the configuration built when the
# server starts. The module frame keeps no locals (the testServer masking rule; test-omic-design-module.R). Recorded
# divergences from conventions/R-shiny.md, kept under conventions/R.md's existing-code scope rule: the copy and the
# render are not ExtendedTasks (D3 keeps them in process); the wellPanel/fluidRow markup is kept byte-identical; the
# moved functions keep their names.

#' The Summary & Export tab's UI: the layout of an omic's `summary` group.
#' @param id The module id.
#' @param config The omic's `summary` group (`omicGuiConfig(<omic>)$summary`).
#' @return A tag list.
#' @noRd
mod_omic_summary_ui <- function(id, config) {
    ns <- shiny::NS(id)
    tp <- config$testid

    shiny::tagList(
        shiny::fluidPage(
            shiny::h3("Session Summary & Report Generation"),
            lapply(config$export_steps, omicSummaryExportStepUi, ns = ns),

            shiny::fluidRow(
                # Left column: Workflow Parameters
                shiny::column(6,
                    shiny::wellPanel(
                        shiny::h4("Workflow Parameters"),
                        shiny::textInput(ns("experiment_label"), "Experiment Label:",
                                        value = "", placeholder = config$placeholder),
                        shiny::textAreaInput(ns("description"), "Description:",
                                            value = config$description,
                                            rows = 3, resize = "vertical"),
                        shiny::br(),
                        testid(
                            shiny::actionButton(ns("save_workflow_args"), "Save Workflow Arguments",
                                               class = "btn-primary", icon = shiny::icon("save")),
                            paste0(tp, "-summary-save-args")
                        )
                    )
                ),

                # Right column: File Management
                shiny::column(6,
                    shiny::wellPanel(
                        shiny::h4("File Management"),
                        shiny::br(),
                        testid(
                            shiny::actionButton(ns("copy_to_publication"), "Copy to Publication Directory",
                                               class = "btn-info", icon = shiny::icon("copy")),
                            paste0(tp, "-summary-copy-publication")
                        ),
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
                        testid(
                            shiny::actionButton(ns("generate_report"), "Generate Report",
                                               class = "btn-success", icon = shiny::icon("file-pdf")),
                            paste0(tp, "-summary-generate-report")
                        ),
                        shiny::br(), shiny::br(),
                        shiny::conditionalPanel(
                            condition = paste0("output['", ns("report_ready"), "']"),
                            testid(
                                shiny::downloadButton(ns("download_report"), "Download Report",
                                                     class = "btn-success"),
                                paste0(tp, "-summary-download-report")
                            )
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
                        testid(
                            shiny::actionButton(ns("export_session_state"), "Export Session State (.RDS)",
                                               class = "btn-secondary", icon = shiny::icon("download")),
                            paste0(tp, "-summary-export-state")
                        )
                    )
                )
            )
        )
    )
}

#' An export step's panel (D4; trx's `trxSummaryExportUi` at 5be34ee5, mod_trx_summary.R:43-54): its title, help text,
#' the choice's output, the Run button, the status and the figure.
#' @noRd
omicSummaryExportStepUi <- function(step, ns) {
    shiny::fluidRow(shiny::column(12, shiny::wellPanel(
        shiny::h4(step$title),
        do.call(shiny::helpText, as.list(step$help)),
        shiny::uiOutput(ns(step$choice$output)),
        testid(shiny::actionButton(ns(step$button$input), step$button$label, class = step$button$class,
            icon = shiny::icon(step$button$icon)), step$button$testid),
        shiny::br(), shiny::br(), shiny::verbatimTextOutput(ns(step$status$output)),
        shinyjqui::jqui_resizable(shiny::imageOutput(ns(step$figure$output), height = step$figure$height)))))
}

#' The Summary & Export tab's server.
#' @param id,project_dirs,omic_type,experiment_label,workflow_data The tab's arguments, as the shell's `summary` call
#'   passes them (`server(module_id, project_dirs, omic, label, data)`).
#' @param config The omic's `summary` group (`omicGuiConfig(<omic>)$summary`).
#' @return The tab's reactive values.
#' @noRd
mod_omic_summary_server <- function(id, project_dirs, omic_type, experiment_label, workflow_data, config) {
    shiny::moduleServer(id, function(input, output, session) {
        runOmicSummaryServer(input, output, session, project_dirs, omic_type, experiment_label, workflow_data, config)
    })
}

#' The tab's registrations, in today's order: the export steps' panels, the label's prefill, the reactive values, the
#' template status, the save, copy, report, GitHub and session-export observers, then the initial outputs.
#' @noRd
runOmicSummaryServer <- function(input, output, session, projectDirs, omicType, experimentLabel, workflowData, config) {
    lapply(config$export_steps, function(step) {
        registerOmicSummaryExportStep(input, output, session, step, workflowData, projectDirs[[omicType]], omicType)
    })
    if (!is.null(experimentLabel)) shiny::updateTextInput(session, "experiment_label", value = experimentLabel)
    values <- shiny::reactiveValues(workflow_args_saved = FALSE, files_copied = FALSE, report_generated = FALSE,
        report_path = NULL)

    output$template_status <- shiny::renderText({
        shiny::req(projectDirs)
        workflowSummaryTemplateStatus(config, projectDirs, omicType)
    })

    shiny::observeEvent(input$save_workflow_args, {
        shiny::req(input$experiment_label)
        omicSummarySave(output, values, input$experiment_label, input$description, projectDirs, omicType, workflowData,
            config)
    })

    shiny::observeEvent(input$copy_to_publication, {
        shiny::req(input$experiment_label)
        omicSummaryCopy(output, values, input$experiment_label, input$description, projectDirs, omicType, workflowData,
            config)
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$generate_report, {
        # KIT-101 F-11: a refused report says why before the gate stops it.
        if (!shiny::isTruthy(input$experiment_label)) {
            shiny::showNotification("Please enter an experiment label before generating the report.", type = "warning")
        } else if (!isTRUE(values$files_copied)) {
            shiny::showNotification(paste("Please copy the results to the publication directory before generating the report",
                "(Copy to Publication Directory)."), type = "warning")
        }
        shiny::req(input$experiment_label)
        shiny::req(values$files_copied)
        omicSummaryReport(output, values, input$experiment_label, input$description, projectDirs, omicType, workflowData,
            config)
    })

    shiny::observeEvent(input$push_to_github, {
        shiny::req(input$enable_github, input$github_org, input$github_email, input$github_username, input$project_id)
        shiny::req(values$report_generated)
        omicSummaryGithub(output, input$experiment_label, input$description, input$github_org, input$github_email,
            input$github_username, input$project_id, projectDirs, omicType)
    })

    shiny::observeEvent(input$export_session_state, {
        shiny::req(input$experiment_label)
        writeWorkflowSummarySessionState(config, workflowData, projectDirs, omicType, input$experiment_label,
            input$description, values = list(workflow_args_saved = values$workflow_args_saved,
                files_copied = values$files_copied, report_generated = values$report_generated,
                report_path = values$report_path), notify = omicSummaryNotify)
    })

    output$session_summary <- shiny::renderText({
        "Ready to save workflow parameters and generate report"
    })
    output$report_ready <- shiny::reactive({ FALSE })
    shiny::outputOptions(output, "report_ready", suspendWhenHidden = FALSE)

    invisible(values)
}

#' The core's notifications, shown.
#' @noRd
omicSummaryNotify <- function(message, type, duration) {
    shiny::showNotification(message, type = type, duration = duration)
}

#' The core's progress steps, shown.
#' @noRd
omicSummaryProgress <- function(value, detail) shiny::incProgress(value, detail = detail)

#' The save: the core's `saveWorkflowSummaryParameters`, then the session summary.
#' @noRd
omicSummarySave <- function(output, values, label, description, projectDirs, omicType, workflowData, config) {
    saveWorkflowSummaryParameters(config, workflowData, projectDirs, omicType, label, description,
        notify = omicSummaryNotify,
        done = function(saved) {
            values$workflow_args_saved <- TRUE
            output$session_summary <- shiny::renderText({
                paste(
                    "Study parameters created for:", label,
                    "\nDescription:", description,
                    "\nTimestamp:", Sys.time(),
                    "\nFile:", saved$studyParamsFile,
                    "\nSource: Final S4 object @args + config_list",
                    "\nIntegration Object:", saved$integrationFile,
                    "\nStatus: Parameters saved [OK]"
                )
            })
        },
        failed = function(error) values$workflow_args_saved <- TRUE)
}

#' The copy: the core's `copyWorkflowSummaryResults` under the progress display, then the copy status and the session
#' summary.
#' @noRd
omicSummaryCopy <- function(output, values, label, description, projectDirs, omicType, workflowData, config) {
    copied <- shiny::withProgress(message = "Copying files to publication directory...",
        copyWorkflowSummaryResults(config, workflowData, projectDirs, omicType, label, description,
            argsSaved = values$workflow_args_saved, notify = omicSummaryNotify))
    if (isTRUE(copied$fallbackWritten)) values$workflow_args_saved <- TRUE
    if (!is.null(copied$error)) {
        output$copy_status <- shiny::renderText(paste("Error:", copied$error))
        return(invisible(copied))
    }
    values$files_copied <- TRUE
    output$copy_status <- shiny::renderText(
        if (length(copied$optional) > 0L) {
            "Files copied to publication directory with optional missing artifacts [WARNING]"
        } else {
            "Files copied to publication directory successfully [OK]"
        }
    )
    output$session_summary <- shiny::renderText({
        paste(
            "Workflow args created for:", label,
            "\nDescription:", description,
            "\nTimestamp:", Sys.time(),
            "\nStatus: Arguments saved [OK], Files copied [OK]"
        )
    })
    invisible(copied)
}

#' The report: the core's `renderWorkflowSummaryReport` under the progress display, then the rendered report's values,
#' its download and the session summary.
#' @noRd
omicSummaryReport <- function(output, values, label, description, projectDirs, omicType, workflowData, config) {
    report <- renderWorkflowSummaryReport(config, workflowData, projectDirs, omicType, label, notify = omicSummaryNotify,
        progress = omicSummaryProgress, within = function(message, code) shiny::withProgress(message = message, code))
    if (is.null(report$path)) return(invisible(report))
    values$report_generated <- TRUE
    values$report_path <- report$path
    output$report_ready <- shiny::reactive({ TRUE })
    shiny::outputOptions(output, "report_ready", suspendWhenHidden = FALSE)
    output$download_report <- shiny::downloadHandler(
        filename = function() basename(report$path),
        content = function(file) file.copy(report$path, file)
    )
    output$session_summary <- shiny::renderText({
        paste(
            "Workflow args created for:", label,
            "\nDescription:", description,
            "\nTimestamp:", Sys.time(),
            "\nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK]",
            "\nReport location:", report$path
        )
    })
    invisible(report)
}

#' The GitHub push (`pushProjectToGithubFromDirs`), then the session summary.
#' @noRd
omicSummaryGithub <- function(output, label, description, org, email, username, projectId, projectDirs, omicType) {
    shiny::withProgress(message = "Pushing to GitHub...", tryCatch({
        options(github_org = org, github_user_email = email, github_user_name = username)
        pushProjectToGithubFromDirs(project_dirs = projectDirs, omic_type = omicType, experiment_label = label,
            project_id = projectId)
        shiny::showNotification("Successfully pushed to GitHub", type = "message")
        output$session_summary <- shiny::renderText({
            paste(
                "Workflow args created for:", label,
                "\nDescription:", description,
                "\nTimestamp:", Sys.time(),
                "\nStatus: Arguments saved [OK], Files copied [OK], Report generated [OK], GitHub pushed [OK]"
            )
        })
        invisible(TRUE)
    }, error = function(e) {
        shiny::showNotification(paste("GitHub push failed:", conditionMessage(e)), type = "error", duration = 10)
        logger::log_error(logger::skip_formatter(paste("Failed to push to GitHub:", conditionMessage(e))))
        invisible(NULL)
    }))
}

#' An export step's panel server (D4; trx's `trxSummaryServer` at 5be34ee5, mod_trx_summary.R:96-117): the choice's
#' radio buttons (with the unavailable text when the head's design does not offer the alternative), the Run (the core's
#' `runWorkflowSummaryExportStep`; a failure notified with the step's error prefix), the status and the figure.
#' @noRd
registerOmicSummaryExportStep <- function(input, output, session, step, workflowData, experimentPaths, omicType) {
    shown <- shiny::reactiveVal(NULL)
    choice <- step$choice
    output[[choice$output]] <- shiny::renderUI({
        workflowData$tab_status
        available <- choice$available(workflowData$state_manager)
        shiny::tagList(shiny::radioButtons(session$ns(choice$input), choice$label, choices = choice$choices(available),
            selected = choice$selected), if (!available) do.call(shiny::helpText, as.list(choice$unavailable)))
    })
    shiny::observeEvent(input[[step$button$input]], {
        shiny::req(workflowData$state_manager)
        runWorkflowSummaryExportStep(step, workflowData, experimentPaths, input[[choice$input]] %||% choice$selected,
            done = shown, binding = omicGuiConfig(omicType)$binding$binding, failed = function(error) {
                message <- paste(step$error$prefix, trxConditionText(error))
                logger::log_error(message)
                shiny::showNotification(message, type = step$error$type, duration = step$error$duration)
            })
    })
    output[[step$status$output]] <- shiny::renderText(if (is.null(shown())) step$status$idle else shown()$text)
    output[[step$figure$output]] <- shiny::renderImage({
        shown()
        path <- file.path(experimentPaths[[step$figure$dir]] %||% "", step$figure$file)
        if (!is.null(shown()) && file.exists(path)) list(src = path, contentType = step$figure$image[["content_type"]],
            width = step$figure$image[["width"]], height = step$figure$image[["height"]], alt = step$figure$alt) else
            list(src = "", alt = step$figure$alt_idle)
    }, deleteFile = FALSE)
    invisible(shown)
}
