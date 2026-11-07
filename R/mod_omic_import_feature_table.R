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

# KIT-102 W5: the feature-table source panel of the import frame (mod_omic_import.R), shared by
# lipidomics, metabolomics and transcriptomics. It renders and serves each omic's import from
# `config$import`: every DOM id, test id, label, choice, default, help text, notification and message
# is the omic's today (the lipid panel's layout, which metab takes). The computation is the core's
# (`previewOmicsImport()`, `importOmicsTables()`, utils_workflow_omics_import.R); the panel keeps the
# notifications, the logger lines and the commit to the bus, then calls the frame's completion
# helper. Recorded divergence from conventions/R-shiny.md: the wellPanel layout, the import's
# fileInput fallback in test mode and the in-process import are kept (the byte-identical contract;
# decision D4), under the existing-code scope rule of conventions/R.md.

# ------------------------------------------------------------------------------------------------
# The configuration's plugins.
# ------------------------------------------------------------------------------------------------

#' The configuration with its import plugins: `import$preview` (default `previewOmicsImport`) and,
#' for transcriptomics, `import$preview = previewTrxCountsImport` and the count policy
#' `import$check = checkTrxImportCounts` while the configuration does not carry them yet.
#' @noRd
omicImportFeatureTablePlugins <- function(config) {
    import <- config$import
    trx <- identical(config$shell$omic_type, "transcriptomics")
    if (is.null(import$preview)) import$preview <- if (trx) previewTrxCountsImport else previewOmicsImport
    if (is.null(import$check) && trx) import$check <- checkTrxImportCounts
    config$import <- import
    config
}

#' The preview of a first assay file, by the configuration's preview plugin.
#' @param path The assay file.
#' @param config An omic GUI configuration (`omicGuiConfig()`).
#' @return The preview (`previewOmicsImport()`'s shape).
#' @noRd
omicImportFeatureTablePreview <- function(path, config) {
    import <- omicImportFeatureTablePlugins(config)$import
    import$preview(path, import = import)
}

#' Whether the import offers native file pickers: shinyFiles is installed and the app is not in
#' test mode (shinytest2 uploads through `fileInput`).
#' @noRd
omicImportUsesShinyFiles <- function() {
    requireNamespace("shinyFiles", quietly = TRUE) && !is_test_mode()
}

#' A `fileInput` value's first upload path; NULL when there is none.
#' @noRd
omicImportStandardFilePath <- function(file_input) {
    if (is.null(file_input) || is.null(file_input$datapath) || length(file_input$datapath) == 0L) {
        return(NULL)
    }
    path <- as.character(file_input$datapath[1])
    if (!nzchar(path)) NULL else path
}

# ------------------------------------------------------------------------------------------------
# The UI.
# ------------------------------------------------------------------------------------------------

#' A tag with the configured test id `key`, when the omic has one.
#' @noRd
omicImportTestid <- function(tag, import, key) {
    id <- unname(import$testids[key])
    if (is.na(id)) tag else testid(tag, id)
}

#' The feature-table panel's UI (the frame's `ui(ns, config)`): the file import and column mapping
#' columns, a rule and the process footer.
#' @param ns The frame's namespace function.
#' @param config An omic GUI configuration (`omicGuiConfig()`).
#' @noRd
omicImportFeatureTableUi <- function(ns, config) {
    import <- config$import
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(6, omicImportFileSection(ns, import, omicImportUsesShinyFiles())),
            shiny::column(6, omicImportMappingSection(ns, import))
        ),
        shiny::hr(),
        shiny::fluidRow(shiny::column(
            12,
            omicImportTestid(
                shiny::actionButton(ns("process_import"), "Process Imported Data", class = "btn-success",
                    width = "100%", icon = shiny::icon("check")),
                import, "process"
            ),
            shiny::br(),
            shiny::br(),
            shiny::uiOutput(ns("import_status"))
        ))
    )
}

#' One assay's name and file inputs in a light well panel.
#' @param slot `assay1` or `assay2`.
#' @noRd
omicImportAssayPanel <- function(ns, import, slot, label, value, placeholder, use_shiny_files) {
    name_input <- omicImportTestid(
        shiny::textInput(ns(paste0(slot, "_name")), label, value = value, placeholder = placeholder),
        import, paste0(slot, "_name")
    )
    file_input <- if (use_shiny_files) {
        shiny::tagList(
            shinyFiles::shinyFilesButton(id = ns(paste0(slot, "_file")), label = "Select File",
                title = "Choose Data File", multiple = FALSE, icon = shiny::icon("file")),
            shiny::br(),
            shiny::verbatimTextOutput(outputId = ns(paste0(slot, "_path")), placeholder = TRUE)
        )
    } else {
        shiny::fileInput(inputId = ns(paste0(slot, "_file_std")), label = NULL,
            accept = c(".tsv", ".tab", ".txt", ".csv", ".xlsx", ".parquet"))
    }
    file_testid <- unname(import$testids[paste0(slot, "_file")])
    if (!is.na(file_testid)) file_input <- shiny::div("data-testid" = file_testid, file_input)
    shiny::wellPanel(
        style = "background-color: #f8f9fa;",
        shiny::fluidRow(shiny::column(4, name_input), shiny::column(8, file_input))
    )
}

#' Steps 1 and 2: the format choice and the assay files.
#' @noRd
omicImportFileSection <- function(ns, import, use_shiny_files) {
    shiny::tagList(
        shiny::h4(import$format_heading),
        omicImportTestid(
            shiny::radioButtons(ns("vendor_format"), NULL, choices = import$formats,
                selected = import$default_format, inline = TRUE),
            import, "vendor_format"
        ),
        shiny::hr(),
        shiny::h4("Step 2: Import Data Files"),
        shiny::p(import$help),
        omicImportAssayPanel(ns, import, "assay1", "Assay Name", import$assay_name, NULL,
            use_shiny_files),
        if (isTRUE(import$second_assay)) {
            omicImportAssayPanel(ns, import, "assay2", "Assay Name (Optional)", "", "e.g., LCMS_Neg",
                use_shiny_files)
        }
    )
}

#' A `conditionalPanel` shown for the custom format (`custom = TRUE`) or the others.
#' @noRd
omicImportFormatPanel <- function(ns, custom, ...) {
    shiny::conditionalPanel(
        condition = paste0("input['", ns("vendor_format"), "'] ", if (custom) "==" else "!=", " 'custom'"),
        ...
    )
}

#' One column choice (a select, or for the custom format a text input) beside its status output.
#' @noRd
omicImportColumnRow <- function(input, status) {
    shiny::fluidRow(shiny::column(6, input), shiny::column(6, status))
}

#' The omic's extra mapping inputs, in its order (`config$import$aux_inputs`).
#' @noRd
omicImportAuxInputs <- function(ns, import) {
    inputs <- list(
        exclude_norm = shiny::tagList(
            shiny::checkboxInput(ns("exclude_norm"), "Exclude Vendor-Normalized Columns", value = TRUE),
            shiny::helpText("Excludes columns ending in _Norm, _Normalized, etc.")
        ),
        sanitize_names = shiny::tagList(
            omicImportTestid(
                shiny::checkboxInput(ns("sanitize_names"), "Sanitize Sample Names", value = TRUE),
                import, "sanitize_names"
            ),
            shiny::helpText(paste("Clean sample IDs (e.g., '123-Sample!' -> 'x123_sample') for better",
                "compatibility with downstream analysis."))
        ),
        is_pattern = shiny::tagList(
            omicImportTestid(
                shiny::textInput(ns("is_pattern"), "Internal Standard Pattern (Regex)", value = "",
                    placeholder = "e.g., ^IS_|_d[0-9]+$|ISTD"),
                import, "is_pattern"
            ),
            shiny::helpText("Regular expression to identify internal standards")
        )
    )
    unname(inputs[import$aux_inputs])
}

#' Steps 3 and 4: the column mapping, the detected sample columns and the validation.
#' @noRd
omicImportMappingSection <- function(ns, import) {
    ids <- import$inputs
    outputs <- import$outputs
    shiny::tagList(
        shiny::h4("Step 3: Column Mapping"),
        shiny::p("Configure column mappings for your data. Dropdowns are auto-populated based on detected format."),
        shiny::conditionalPanel(
            condition = paste0("output['", ns("file_loaded"), "']"),
            shiny::wellPanel(
                omicImportFormatPanel(ns, FALSE, shiny::uiOutput(ns("format_detection_status")), shiny::hr()),
                omicImportFormatPanel(ns, TRUE,
                    shiny::div(class = "alert alert-info", shiny::icon("edit"), shiny::strong(" Custom Format"),
                        shiny::br(), "Type the exact column names from your file below."),
                    shiny::hr()),
                omicImportFormatPanel(ns, FALSE, omicImportColumnRow(
                    omicImportTestid(shiny::selectInput(ns(ids[["id_col"]]), import$id_label, choices = NULL),
                        import, "id_col"),
                    shiny::uiOutput(ns(outputs[["id_status"]])))),
                omicImportFormatPanel(ns, TRUE, omicImportColumnRow(
                    shiny::textInput(ns(ids[["id_col_custom"]]), import$id_custom_label, value = "",
                        placeholder = import$id_custom_placeholder),
                    shiny::uiOutput(ns(outputs[["id_status_custom"]])))),
                omicImportFormatPanel(ns, FALSE, omicImportColumnRow(
                    omicImportTestid(shiny::selectInput(ns("annotation_col"), import$annotation_label,
                        choices = NULL), import, "annotation_col"),
                    shiny::uiOutput(ns("annotation_status")))),
                omicImportFormatPanel(ns, TRUE, omicImportColumnRow(
                    shiny::textInput(ns("annotation_col_custom"), "Annotation Column (optional)", value = "",
                        placeholder = import$annotation_custom_placeholder),
                    shiny::uiOutput(ns("annotation_status_custom")))),
                omicImportFormatPanel(ns, TRUE,
                    shiny::textInput(ns("sample_cols_pattern"), "Sample Column Pattern (Regex)", value = "",
                        placeholder = "e.g., ^Sample_|_Rep[0-9]+$"),
                    shiny::helpText("Leave blank to auto-detect numeric columns as samples")),
                omicImportAuxInputs(ns, import),
                shiny::hr(),
                shiny::h5("Detected Sample Columns"),
                shiny::verbatimTextOutput(ns("sample_columns_display")),
                omicImportFormatPanel(ns, TRUE, shiny::hr(), shiny::h5("Available Columns in File"),
                    shiny::verbatimTextOutput(ns("available_columns_display")))
            ),
            shiny::hr(),
            shiny::h4("Step 4: Validation"),
            shiny::uiOutput(ns("validation_summary"))
        ),
        shiny::conditionalPanel(
            condition = paste0("!output['", ns("file_loaded"), "']"),
            shiny::div(class = "alert alert-info", shiny::icon("info-circle"),
                " Import a data file to configure column mappings.")
        )
    )
}

# ------------------------------------------------------------------------------------------------
# The displays.
# ------------------------------------------------------------------------------------------------

#' The detected format and its confidence, as an alert coloured by the confidence.
#' @param labels `config$import$format_labels`; a format without one is "Unknown".
#' @noRd
omicImportFormatStatus <- function(detected_format, format_confidence, labels) {
    confidence_pct <- round(format_confidence * 100)
    color_class <- if (confidence_pct >= 70) "success" else if (confidence_pct >= 40) "warning" else "danger"
    label <- if (detected_format %in% names(labels)) labels[[detected_format]] else "Unknown"
    shiny::tags$div(
        class = paste("alert", paste0("alert-", color_class)),
        shiny::tags$strong("Detected format: "),
        label,
        shiny::tags$br(),
        shiny::tags$small(sprintf("Confidence: %d%%", confidence_pct))
    )
}

#' A column choice's status: empty (optional, or a prompt), not found, or found with its message.
#' @param modes One of `config$import$status_modes`: `success` (`unique_ids`, `found_unique_ids`,
#'   `ok`, `found`), optional `empty` (`optional`, `prompt`) and `case` (`insensitive`).
#' @noRd
omicImportColumnStatus <- function(assay_data, column, modes) {
    status <- function(icon, color, text) {
        shiny::tags$span(shiny::icon(icon, style = sprintf("color: %s;", color)), text)
    }
    if (is.null(column) || !nzchar(column)) {
        if (modes["empty"] %in% "prompt") return(status("question-circle", "gray", " Enter column name"))
        return(status("minus-circle", "gray", " Optional"))
    }
    headers <- names(assay_data)
    actual <- if (modes["case"] %in% "insensitive") {
        headers[match(tolower(column), tolower(headers))]
    } else {
        headers[match(column, headers)]
    }
    if (length(actual) == 0L || is.na(actual)) return(status("times-circle", "red", " Column not found"))
    switch(modes[["success"]],
        unique_ids = status("check-circle", "green", sprintf(" %d unique IDs", length(unique(assay_data[[actual]])))),
        found_unique_ids = status("check-circle", "green",
            sprintf(" Found: %d unique IDs", length(unique(assay_data[[actual]])))),
        ok = status("check-circle", "green", " OK"),
        found = status("check-circle", "green", " Found"))
}

#' Column names as one line, the first `max_columns` and a count of the rest.
#' @noRd
omicImportColumnPreviewText <- function(columns, max_columns = Inf) {
    if (length(columns) > max_columns) {
        return(paste(paste(utils::head(columns, max_columns), collapse = ", "),
            sprintf("... and %d more", length(columns) - max_columns)))
    }
    paste(columns, collapse = ", ")
}

#' The first assay's mapping validation, passed (with the omic's feature count and any warnings) or
#' failed (with its errors).
#' @noRd
omicImportValidationSummary <- function(validation, import) {
    if (validation$valid) {
        return(shiny::tagList(
            shiny::tags$div(class = "alert alert-success", shiny::icon("check-circle"),
                shiny::tags$strong(" Validation Passed")),
            shiny::tags$ul(
                shiny::tags$li(sprintf("%s: %d", import$validation_noun, validation$summary[[import$count_key]])),
                shiny::tags$li(sprintf("Samples: %d", validation$summary$n_samples)),
                shiny::tags$li(sprintf("Missing values: %.1f%%", validation$summary$pct_missing))
            ),
            if (length(validation$warnings) > 0) {
                shiny::tags$div(class = "alert alert-warning", shiny::icon("exclamation-triangle"), " Warnings:",
                    shiny::tags$ul(lapply(validation$warnings, shiny::tags$li)))
            }
        ))
    }
    shiny::tags$div(class = "alert alert-danger", shiny::icon("times-circle"),
        shiny::tags$strong(" Validation Failed"), shiny::tags$ul(lapply(validation$errors, shiny::tags$li)))
}

#' The completed import's line: format, assays and samples; NULL before the import completes.
#' @noRd
omicImportStatusDisplay <- function(workflow_data) {
    if (!identical(workflow_data$tab_status$setup_import, "complete")) return(NULL)
    log_info <- workflow_data$processing_log$setup_import
    shiny::tags$div(
        class = "alert alert-success",
        shiny::icon("check-circle"),
        shiny::tags$strong(" Import Complete"),
        shiny::tags$br(),
        sprintf("Format: %s | Assays: %d | Samples: %d", toupper(log_info$detected_format), log_info$n_assays,
            log_info$n_samples)
    )
}

# ------------------------------------------------------------------------------------------------
# The processing and the commit.
# ------------------------------------------------------------------------------------------------

#' Process an import: the working notification, the core import (`importOmicsTables()`), the
#' sanitising notice, the commit to the bus, then the frame's completion helper. A refusal or a
#' reader error is logged and shown as `Error: <message>`; nothing is committed.
#' @param workflow_data The workflow bus.
#' @param request The core's request (`importOmicsTables()`).
#' @param config An omic GUI configuration (`omicGuiConfig()`).
#' @param complete The frame's completion helper (`complete(notify = TRUE)`).
#' @return `list(status = "success", payload)` or `list(status = "error", error)`, invisibly.
#' @noRd
processOmicImportFeatureTable <- function(workflow_data, request, config, complete,
                                          showNotification = shiny::showNotification,
                                          removeNotification = shiny::removeNotification,
                                          logError = logger::log_error, logInfo = logger::log_info) {
    config <- omicImportFeatureTablePlugins(config)
    working <- config$import$notification_id
    showNotification("Processing imported data...", id = working, duration = NULL)
    tryCatch({
        payload <- importOmicsTables(request, config)
        if (isTRUE(payload$sanitized)) {
            logInfo(sprintf("Sanitized %d sample column names.", length(payload$sample_columns)))
            showNotification("Sample names sanitized for R compatibility.", type = "message")
        }
        commitOmicImportFeatureTable(workflow_data, payload, config, logInfo)
        removeNotification(working)
        complete()
        invisible(list(status = "success", payload = payload))
    }, error = function(error) {
        logError(paste("Error processing import:", conditionMessage(error)))
        removeNotification(working)
        showNotification(paste("Error:", conditionMessage(error)), type = "error", duration = 10)
        invisible(list(status = "error", error = error))
    })
}

#' The import's commit (KIT-075): the owner bound to the omic's route, the tables held by reference,
#' the format, data type and column mapping, the owner's workflow type, and the processing log.
#' @noRd
commitOmicImportFeatureTable <- function(workflow_data, payload, config, logInfo = logger::log_info) {
    shell <- config$shell
    workflow_type <- shell$routes[[1L]]
    bindWorkflowOmicsOwnerAtImport(workflow_data, workflow_type, payload$data_format)
    setWorkflowOmicsBusTables(workflow_data, "data_tbl", payload$assays)
    workflow_data$data_format <- payload$data_format
    workflow_data$data_type <- shell$bus$data_type
    workflow_data$column_mapping <- payload$column_mapping
    if (!is.null(workflow_data$state_manager)) {
        workflow_data$state_manager$setWorkflowType(workflow_type)
        logInfo(sprintf("Workflow type set to: %s", workflow_type))
    }
    workflow_data$processing_log$setup_import <- payload$processing_log
    logInfo(sprintf("%s import complete: %d assays, %d total %s", shell$label, length(payload$assays),
        sum(vapply(payload$assays, nrow, integer(1))), tolower(config$import$validation_noun)))
    invisible(workflow_data)
}

# ------------------------------------------------------------------------------------------------
# The server.
# ------------------------------------------------------------------------------------------------

#' The feature-table panel's server (the frame's `server(...)`, inside its moduleServer): the file
#' inputs and the preview, the column choices, the outputs and the process observer.
#' @param input,output,session The frame's.
#' @param workflow_data The workflow bus.
#' @param experiment_paths The experiment's project paths (the artifacts' `source_dir`).
#' @param volumes shinyFiles volumes; NULL for the system's.
#' @param config An omic GUI configuration (`omicGuiConfig()`).
#' @param complete The frame's completion helper.
#' @noRd
omicImportFeatureTableServer <- function(input, output, session, workflow_data, experiment_paths, volumes,
                                         config, complete) {
    import <- config$import
    ids <- import$inputs
    outputs <- import$outputs
    local_data <- shiny::reactiveValues(assay1_file = NULL, assay1_data = NULL, assay1_import_result = NULL,
        assay2_file = NULL, detected_format = NULL, format_confidence = NULL, all_headers = NULL)
    slots <- if (isTRUE(import$second_assay)) c("assay1", "assay2") else "assay1"

    load_preview <- function() {
        shiny::req(local_data$assay1_file)
        tryCatch({
            preview <- omicImportFeatureTablePreview(local_data$assay1_file, config)
            local_data$all_headers <- preview$headers
            local_data$detected_format <- preview$detected_format
            local_data$format_confidence <- preview$format_confidence
            local_data$assay1_import_result <- preview$import_result
            local_data$assay1_data <- preview$assay_data
            shiny::updateSelectInput(session, ids[["id_col"]], choices = preview$updates$id$choices,
                selected = preview$updates$id$selected)
            shiny::updateSelectInput(session, "annotation_col", choices = preview$updates$annotation$choices,
                selected = preview$updates$annotation$selected)
            if (!is.null(preview$updates$is_pattern)) {
                shiny::updateTextInput(session, "is_pattern", value = preview$updates$is_pattern)
            }
            logger::log_info(sprintf("Imported assay 1: %d rows, %d columns, format: %s",
                nrow(preview$assay_data), ncol(preview$assay_data), preview$detected_format))
        }, error = function(error) {
            message <- paste("Error importing data:", conditionMessage(error))
            logger::log_error(message)
            shiny::showNotification(message, type = "error")
        })
    }
    select_file <- function(slot, path) {
        local_data[[paste0(slot, "_file")]] <- path
        output[[paste0(slot, "_path")]] <- shiny::renderText(path)
        if (identical(slot, "assay1")) load_preview()
    }

    if (omicImportUsesShinyFiles()) {
        if (is.null(volumes)) volumes <- shinyFiles::getVolumes()()
        for (slot in slots) local({
            slot <- slot
            file_id <- paste0(slot, "_file")
            shinyFiles::shinyFileChoose(input, file_id, roots = volumes, session = session,
                filetypes = c("tsv", "tab", "txt", "csv", "xlsx", "parquet"))
            shiny::observeEvent(input[[file_id]], {
                chosen <- input[[file_id]]
                if (is.null(chosen) || is.integer(chosen)) return()
                tryCatch({
                    files <- shinyFiles::parseFilePaths(volumes, chosen)
                    if (nrow(files) > 0) select_file(slot, as.character(files$datapath[1]))
                }, error = function(error) {
                    logger::log_error(paste("Error parsing file path:", conditionMessage(error)))
                })
            })
        })
    } else {
        for (slot in slots) local({
            slot <- slot
            shiny::observeEvent(input[[paste0(slot, "_file_std")]], {
                path <- omicImportStandardFilePath(input[[paste0(slot, "_file_std")]])
                if (!is.null(path)) select_file(slot, path)
            })
        })
    }

    id_col <- shiny::reactive(resolveOmicsImportColumn(names(local_data$assay1_data), input$vendor_format,
        input[[ids[["id_col"]]]], input[[ids[["id_col_custom"]]]]))
    annotation_col <- shiny::reactive(resolveOmicsImportColumn(names(local_data$assay1_data),
        input$vendor_format, input$annotation_col, input$annotation_col_custom))
    sample_columns <- shiny::reactive({
        shiny::req(local_data$assay1_data)
        resolveOmicsImportSampleColumns(local_data$assay1_data, input$vendor_format, input$sample_cols_pattern,
            local_data$assay1_import_result, if (isTRUE(input$exclude_norm)) OMICS_IMPORT_NORMALISED_COLUMNS)
    })

    output$file_loaded <- shiny::reactive(!is.null(local_data$assay1_data))
    shiny::outputOptions(output, "file_loaded", suspendWhenHidden = FALSE)
    output$format_detection_status <- shiny::renderUI({
        shiny::req(local_data$detected_format)
        omicImportFormatStatus(local_data$detected_format, local_data$format_confidence, import$format_labels)
    })
    output[[outputs[["id_status"]]]] <- shiny::renderUI({
        shiny::req(local_data$assay1_data)
        shiny::req(input[[ids[["id_col"]]]])
        omicImportColumnStatus(local_data$assay1_data, input[[ids[["id_col"]]]], import$status_modes$id)
    })
    output$annotation_status <- shiny::renderUI({
        shiny::req(local_data$assay1_data)
        omicImportColumnStatus(local_data$assay1_data, input$annotation_col, import$status_modes$annotation)
    })
    output$sample_columns_display <- shiny::renderText({
        shiny::req(local_data$assay1_import_result)
        omicImportColumnPreviewText(local_data$assay1_import_result$sample_columns, max_columns = 10)
    })
    output$available_columns_display <- shiny::renderText({
        shiny::req(local_data$all_headers)
        omicImportColumnPreviewText(local_data$all_headers)
    })
    output[[outputs[["id_status_custom"]]]] <- shiny::renderUI({
        shiny::req(local_data$assay1_data)
        omicImportColumnStatus(local_data$assay1_data, input[[ids[["id_col_custom"]]]],
            import$status_modes$id_custom)
    })
    output$annotation_status_custom <- shiny::renderUI({
        shiny::req(local_data$assay1_data)
        omicImportColumnStatus(local_data$assay1_data, input$annotation_col_custom,
            import$status_modes$annotation_custom)
    })
    output$validation_summary <- shiny::renderUI({
        shiny::req(local_data$assay1_data)
        shiny::req(id_col())
        omicImportValidationSummary(import$validate_mapping(local_data$assay1_data, id_col(), sample_columns()),
            import)
    })
    output$import_status <- shiny::renderUI(omicImportStatusDisplay(workflow_data))

    shiny::observeEvent(input$process_import, {
        shiny::req(local_data$assay1_data)
        shiny::req(id_col())
        request <- list(assay1_name = input$assay1_name, assay1_data = local_data$assay1_data,
            assay1_file = local_data$assay1_file, assay2_file = local_data$assay2_file,
            assay2_name = input$assay2_name, vendor_format = input$vendor_format,
            detected_format = local_data$detected_format, id_col = id_col(), annotation_col = annotation_col(),
            sample_columns = sample_columns(), sanitize_names = input$sanitize_names, is_pattern = input$is_pattern,
            experiment_paths = experiment_paths)
        processOmicImportFeatureTable(workflow_data, request, config, complete)
    })
    invisible(NULL)
}
