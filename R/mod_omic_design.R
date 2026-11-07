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

# KIT-103 W4: the shared design tab (kit-103 design 4.5), the frame of every omic's design step. It
# renders the lipid design tab's markup (mod_lipid_design_ui.R at 616e1147) with the omic's texts
# from `config$design`, and serves, in the lipid server's order (mod_lipid_design_server.R and its
# helper files at 616e1147): the import-existing-design plugin, the `data_available` and
# `design_matrix_exists` outputs, the builder slot (mod_omic_design_builder.R), the commit's GUI
# side and the saved-results previews. The commit and the import are the core's
# (`commitOmicsDesign()`, `readOmicsDesignImport()`, `commitOmicsDesignImport()`,
# utils_workflow_omics_design.R); the frame keeps the modal, the notifications, the QC trigger and
# the tab status (KIT-102's split). Every id, test id, text and notification is the omic's today.
# Proteomics' import plugin and commit hook (`proteomics`, `prot_hook`) are KIT-103 W6's. Recorded
# divergences from conventions/R-shiny.md: the wellPanel layout, the shinyFiles fallback to a text
# input in test mode, the in-process commit, and the two status outputs assigned a `reactive`
# (their client-side conditions read them; under testServer they read as an error, which the entry
# oracles record) are kept, under the existing-code scope rule of conventions/R.md.

# ------------------------------------------------------------------------------------------------
# The UI.
# ------------------------------------------------------------------------------------------------

#' A tag with the design's configured test id `key`, when the omic has one.
#' @noRd
omicsDesignTestid <- function(tag, design, key) {
    id <- unname(design$testids[key])
    if (is.na(id)) tag else testid(tag, id)
}

#' The builder slot: the shared builder, or proteomics' own builder UI (back-reference #4).
#' @noRd
omicsDesignBuilderSlotUi <- function(ns, config) {
    switch(config$design$builder_ui,
        shared = mod_omic_design_builder_ui(ns("builder"), config),
        proteomics = mod_prot_design_builder_ui(ns("builder")))
}

#' The design tab's UI: the heading with the import button (with an import plugin), the help, the
#' multi-assay banner (when configured), the builder and the saved-results previews while data is
#' available, else the no-data note.
#' @param id Module id (the design tab's).
#' @param config An omic GUI configuration (`omicGuiConfig()`).
#' @noRd
mod_omic_design_ui <- function(id, config) {
    validateOmicGuiConfig(config)
    ns <- shiny::NS(id)
    design <- config$design
    shiny::tagList(shiny::fluidRow(shiny::column(12, shiny::wellPanel(
        shiny::fluidRow(
            shiny::column(8, shiny::h3("Design Matrix Builder")),
            if (!is.null(design$design_import)) shiny::column(4, omicsDesignTestid(
                shiny::actionButton(ns("show_import_modal"), "Import Existing Design", icon = shiny::icon("folder-open"),
                    class = "btn-info pull-right"), design, "import_existing"))),
        shiny::p(design$help),
        if (isTRUE(design$multi_assay_banner)) shiny::div(class = "alert alert-info", style = "margin-top: 10px;",
            shiny::icon("info-circle"), shiny::tags$strong(" Multi-Assay Workflow: "),
            "Design applies uniformly across all assay modes (LCMS_Pos, LCMS_Neg, GCMS, etc.). ",
            "Sample assignments and contrasts are shared."),
        shiny::conditionalPanel(condition = paste0("output['", ns("data_available"), "']"),
            omicsDesignBuilderSlotUi(ns, config),
            shiny::hr(),
            shiny::h3("Saved Results Preview"),
            shiny::p(design$preview_text),
            shiny::conditionalPanel(condition = paste0("output['", ns("design_matrix_exists"), "']"),
                shiny::wellPanel(
                    shiny::h4("Current Design Matrix"),
                    DT::DTOutput(ns("design_matrix_preview")),
                    shiny::br(),
                    shiny::h4("Defined Contrasts"),
                    DT::DTOutput(ns("contrasts_preview")),
                    if (isTRUE(design$preview_assays)) shiny::tagList(shiny::br(), shiny::h4("Assays Included"),
                        shiny::verbatimTextOutput(ns("assays_preview")))))),
        shiny::conditionalPanel(condition = paste0("!output['", ns("data_available"), "']"),
            shiny::div(class = "alert alert-info", shiny::icon("info-circle"), design$no_data_text))))))
}

# ------------------------------------------------------------------------------------------------
# The server.
# ------------------------------------------------------------------------------------------------

#' The design tab's server, in today's order: the import plugin, the status outputs, the builder,
#' the commit of its results and the previews.
#' @param id Module id.
#' @param workflow_data The workflow bus.
#' @param experiment_paths The experiment's project paths (`source_dir`, `base_dir`).
#' @param volumes shinyFiles volumes (a function, a value, or NULL for the system's).
#' @param qc_trigger Optional QC trigger, set TRUE after a commit.
#' @param config An omic GUI configuration (`omicGuiConfig()`) with its design function fields.
#' @return The builder's results (a reactiveVal), invisibly.
#' @noRd
mod_omic_design_server <- function(id, workflow_data, experiment_paths, volumes = NULL, qc_trigger = NULL, config) {
    validateOmicGuiConfig(config)
    import_server <- omicsDesignImportPluginServer(config)
    shiny::moduleServer(id, function(input, output, session) {
        import_server(input, output, session, workflow_data, experiment_paths, volumes, qc_trigger, config)

        output$data_available <- shiny::reactive({
            !is.null(workflow_data$data_tbl) && !is.null(workflow_data$config_list)
        })
        shiny::outputOptions(output, "data_available", suspendWhenHidden = FALSE)
        output$design_matrix_exists <- shiny::reactive({
            !is.null(workflow_data$design_matrix)
        })
        shiny::outputOptions(output, "design_matrix_exists", suspendWhenHidden = FALSE)

        # KIT-075: the tables are read per use and never cached (a function, not a reactive).
        builder_results <- mod_omic_design_builder_server("builder",
            data_tbl = function() getWorkflowOmicsBusTables(workflow_data, "data_tbl"),
            config_list = shiny::reactive(workflow_data$config_list),
            column_mapping = shiny::reactive(workflow_data$column_mapping),
            existing_design_matrix = shiny::reactive(workflow_data$design_matrix),
            existing_contrasts = shiny::reactive(workflow_data$contrasts_tbl),
            state_owner = function() shiny::isolate(workflow_data$state_manager), config = config)

        shiny::observeEvent(builder_results(), {
            results <- builder_results()
            shiny::req(results)
            runOmicsDesignCommit(results, workflow_data, experiment_paths, qc_trigger, config)
        }, ignoreNULL = TRUE)

        registerOmicsDesignPreviewOutputs(output, workflow_data, config)
        invisible(builder_results)
    })
}

#' The saved-results previews: the design, the contrasts and (when configured) the assays.
#' @noRd
registerOmicsDesignPreviewOutputs <- function(output, workflowData, config) {
    output$design_matrix_preview <- DT::renderDT({
        shiny::req(workflowData$design_matrix)
        workflowData$design_matrix
    }, options = list(pageLength = 5, scrollX = TRUE))
    output$contrasts_preview <- DT::renderDT({
        shiny::req(workflowData$contrasts_tbl)
        workflowData$contrasts_tbl
    }, options = list(pageLength = 5, scrollX = TRUE))
    if (isTRUE(config$design$preview_assays)) {
        output$assays_preview <- shiny::renderText({
            shiny::req(workflowData$data_tbl)
            # KIT-075: read for this render.
            paste("Included assays:", paste(names(getWorkflowOmicsBusTables(workflowData, "data_tbl")), collapse = ", "))
        })
    }
    invisible(output)
}

# ------------------------------------------------------------------------------------------------
# The commit's GUI side.
# ------------------------------------------------------------------------------------------------

#' The modal shown while a design is committed.
#' @noRd
omicsDesignProcessingModal <- function() {
    shiny::modalDialog(title = "Processing Design Matrix",
        shiny::div(style = "text-align: center; padding: 20px;", shiny::icon("spinner", class = "fa-spin fa-3x"),
            shiny::br(), shiny::br(), shiny::p("Saving design matrix and preparing data...")),
        footer = NULL, easyClose = FALSE)
}

#' A failed save, as today: logged, the modal removed, then the error notification.
#' @noRd
omicsDesignSaveFailed <- function(message) {
    msg <- paste("Error saving design matrix results:", message)
    logger::log_error(msg)
    shiny::removeModal()
    shiny::showNotification(msg, type = "error", duration = 15)
    invisible(NULL)
}

#' Commit a builder result with the GUI around it: the processing modal, the core commit
#' (`commitOmicsDesign()`), then the QC trigger, the design tab complete and the success
#' notification. A refusal (`omics_design_refused`) shows its text; a failure after the bus write
#' (`omics_design_save_failed`) shows "Error saving design matrix results: <message>". Both remove
#' the modal and leave the QC trigger and the tab status alone.
#' @param results The builder's result.
#' @param workflowData The workflow bus.
#' @param experimentPaths The project paths (`source_dir`).
#' @param qcTrigger Optional QC trigger.
#' @param config An omic GUI configuration with its design function fields.
#' @return NULL, invisibly.
#' @noRd
runOmicsDesignCommit <- function(results, workflowData, experimentPaths, qcTrigger, config) {
    if (!identical(config$design$commit, "shared")) {
        cli::cli_abort("the {.val {config$design$commit}} design commit is not wired yet (KIT-103 W6).",
            class = "omic_gui_config_invalid")
    }
    shiny::showModal(omicsDesignProcessingModal())
    committed <- tryCatch(commitOmicsDesign(workflowData, results, experimentPaths$source_dir, config),
        omics_design_refused = function(e) {
            logger::log_error(e$message)
            shiny::showNotification(e$message, type = "error", duration = 15)
            shiny::removeModal()
            NULL
        },
        omics_design_save_failed = function(e) omicsDesignSaveFailed(e$message))
    if (is.null(committed)) return(invisible(NULL))
    tryCatch({
        if (!is.null(qcTrigger)) qcTrigger(TRUE)
        status <- workflowData$tab_status
        status$design_matrix <- "complete"
        workflowData$tab_status <- status
        shiny::removeModal()
        shiny::showNotification("Design matrix and contrasts saved successfully!", type = "message")
    }, error = function(e) omicsDesignSaveFailed(e$message))
    invisible(NULL)
}

# ------------------------------------------------------------------------------------------------
# The import-existing-design plugin.
# ------------------------------------------------------------------------------------------------

#' The import plugin's server by `design$design_import$plugin`: `feature_table` (lipid, metab).
#' A wide-table omic without an import (trx) registers it too, as today's lipid tab did under a
#' profile (its `import_dir_path` output is recorded by the entry oracles); with no button it
#' cannot be reached. Proteomics' own import (`proteomics`) is KIT-103 W6's.
#' @return A function `(input, output, session, workflow_data, experiment_paths, volumes,
#'   qc_trigger, config)`.
#' @noRd
omicsDesignImportPluginServer <- function(config) {
    plugin <- config$design$design_import$plugin
    if (is.null(plugin) && identical(config$design$data_adapter, "wide_assay_list")) plugin <- "feature_table"
    if (!identical(plugin, "feature_table")) {
        cli::cli_abort("the {.val {plugin}} design import plugin is not wired yet (KIT-103 W6).",
            class = "omic_gui_config_invalid")
    }
    registerOmicsDesignFeatureImport
}

#' The import's volumes: the given ones (or the system's), with the project's base directory first
#' when it exists.
#' @noRd
omicsDesignImportVolumes <- function(volumes, experimentPaths) {
    baseVolumes <- if (is.function(volumes)) volumes() else if (is.null(volumes)) shinyFiles::getVolumes()() else volumes
    baseDir <- experimentPaths$base_dir
    if (is.null(baseDir) || !dir.exists(baseDir)) return(baseVolumes)
    logger::log_info(paste("Added base_dir to volumes:", baseDir))
    c("Project Base Dir" = baseDir, baseVolumes)
}

#' The import modal: the directory chooser (a text input in test mode) and the confirm button, with
#' the omic's test ids.
#' @noRd
buildOmicsDesignImportModal <- function(ns, config, testMode = is_test_mode()) {
    design <- config$design
    dirChooser <- if (testMode) {
        omicsDesignTestid(shiny::textInput(ns("import_dir"), "Import Directory Path:", value = ""), design, "import_dir")
    } else {
        shinyFiles::shinyDirButton(ns("import_dir"), "Select Folder", "Choose a directory")
    }
    shiny::modalDialog(title = "Import Existing Design Matrix",
        shiny::p("Select the folder containing 'design_matrix.tab', assay data files, and optionally 'contrast_strings.tab'."),
        shiny::helpText("Required files: design_matrix.tab, assay_manifest.txt, data_cln_*.tab files"),
        dirChooser,
        shiny::verbatimTextOutput(ns("import_dir_path"), placeholder = TRUE),
        footer = shiny::tagList(shiny::modalButton("Cancel"), omicsDesignTestid(
            shiny::actionButton(ns("confirm_import"), "Import", class = "btn-primary"), design, "import_confirm")))
}

#' The feature-table import plugin's server: the volumes and directory chooser (shinyFiles; the
#' path as typed in test mode), the modal, the chosen path and the confirmation.
#' @noRd
registerOmicsDesignFeatureImport <- function(input, output, session, workflow_data, experiment_paths, volumes,
                                             qc_trigger, config) {
    testMode <- is_test_mode()
    resolvedVolumes <- shiny::isolate(omicsDesignImportVolumes(volumes, experiment_paths))
    if (!testMode) shinyFiles::shinyDirChoose(input, "import_dir", roots = resolvedVolumes, session = session)
    parseDir <- if (testMode) function(roots, selection) selection else shinyFiles::parseDirPath
    shiny::observeEvent(input$show_import_modal, {
        shiny::showModal(buildOmicsDesignImportModal(session$ns, config, testMode))
    })
    output$import_dir_path <- shiny::renderText({
        shiny::req(input$import_dir)
        parseDir(resolvedVolumes, input$import_dir)
    })
    shiny::observeEvent(input$confirm_import, {
        shiny::req(input$import_dir)
        importPath <- parseDir(resolvedVolumes, input$import_dir)
        shiny::req(importPath)
        runOmicsDesignImport(workflow_data, experiment_paths, importPath, qc_trigger, config)
    })
    invisible(NULL)
}

#' Import an exported design directory with the GUI around it (lipid's order): the modal closed,
#' the files checked (a missing one is today's error), the working notification, the core's read
#' and commit (`design_import$reader`, `commitOmicsDesignImport()`), the QC trigger, the design tab
#' complete and the success text; a failure is "Error during import: <message>". With
#' `design_import$check_files_first` (metab, P-11) the files are checked before the modal closes,
#' and a missing one leaves it open.
#' @param importPath The exported design directory.
#' @return NULL, invisibly.
#' @noRd
runOmicsDesignImport <- function(workflowData, experimentPaths, importPath, qcTrigger, config) {
    importConfig <- config$design$design_import
    shiny::req(importConfig)
    checkFirst <- isTRUE(importConfig$check_files_first)
    if (!checkFirst) shiny::removeModal()
    check <- checkOmicsDesignImportFiles(importPath)
    if (!isTRUE(check$ok)) {
        logger::log_error(check$message)
        shiny::showNotification(check$message, type = "error", duration = 10)
        return(invisible(NULL))
    }
    if (checkFirst) shiny::removeModal()
    shiny::showNotification("Importing design files...", id = "importing_design", duration = NULL)
    tryCatch({
        imported <- importConfig$reader(importPath, config, experimentPaths$source_dir)
        committed <- commitOmicsDesignImport(workflowData, imported, experimentPaths$source_dir, config)
        if (!is.null(qcTrigger)) qcTrigger(TRUE)
        status <- workflowData$tab_status
        status$design_matrix <- "complete"
        workflowData$tab_status <- status
        shiny::removeNotification("importing_design")
        shiny::showNotification(sprintf("Design imported successfully! Loaded %d assays with %d samples.",
            committed$n_assays, committed$n_samples), type = "message")
    }, error = function(e) {
        msg <- paste("Error during import:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("importing_design")
    })
    invisible(NULL)
}
