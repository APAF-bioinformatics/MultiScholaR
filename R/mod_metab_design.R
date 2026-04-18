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

# ============================================================================
# mod_metab_design.R
# ============================================================================
# Purpose: Metabolomics design matrix host module - 1:1 with proteomics
#
# This module embeds mod_metab_design_builder.R and handles orchestration,
# including Import Existing Design, S4 object creation, and state management.
# ============================================================================

#' @title Metabolomics Design Matrix Applet Module
#'
#' @description A Shiny module that serves as the main host for the design
#' matrix creation workflow step. Embeds mod_metab_design_builder_ui/server
#' and handles Import Existing Design, S4 object creation, and state saving.
#' Architecture is 1:1 with mod_prot_design.R.
#'
#' @param id Module ID
#' @param workflow_data A reactive values object to store workflow data.
#' @param experiment_paths A list of paths for the current experiment.
#' @param volumes A list of volumes for shinyFiles (optional).
#' @param qc_trigger A reactive trigger for QC execution (optional).
#'
#' @name metabolomicsDesignMatrixAppletModule
NULL

#' @rdname metabolomicsDesignMatrixAppletModule
#' @export
#' @importFrom shiny NS tagList wellPanel h3 h4 p conditionalPanel div icon tags HTML fluidRow column actionButton
#' @importFrom DT DTOutput renderDT
mod_metab_design_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::fluidRow(
            shiny::column(12
                , shiny::wellPanel(
                    shiny::fluidRow(
                        shiny::column(8
                            , shiny::h3("Design Matrix Builder")
                        )
                        , shiny::column(4
                            , shiny::actionButton(ns("show_import_modal"), "Import Existing Design"
                                , icon = shiny::icon("folder-open"), class = "btn-info pull-right")
                        )
                    )
                    , shiny::p("Use the tools below to define your experimental groups and contrasts for differential analysis. All assays (Pos/Neg modes) share the same design matrix.")

                    # Multi-assay info banner
                    , shiny::div(
                        class = "alert alert-info"
                        , style = "margin-top: 10px;"
                        , shiny::icon("info-circle")
                        , shiny::tags$strong(" Multi-Assay Workflow: ")
                        , "Design applies uniformly across all assay modes (LCMS_Pos, LCMS_Neg, GCMS, etc.). "
                        , "Sample assignments and contrasts are shared."
                    )

                    # Conditional panel: show builder when data is available
                    , shiny::conditionalPanel(
                        condition = paste0("output['", ns("data_available"), "']")
                        , if (exists("mod_metab_design_builder_ui")) {
                            mod_metab_design_builder_ui(ns("builder"))
                        } else {
                            shiny::div("Design builder module not loaded")
                        }

                        # Saved results preview
                        , shiny::hr()
                        , shiny::h3("Saved Results Preview")
                        , shiny::p("This section shows the design matrix and contrasts saved to the workflow.")
                        , shiny::conditionalPanel(
                            condition = paste0("output['", ns("design_matrix_exists"), "']")
                            , shiny::wellPanel(
                                shiny::h4("Current Design Matrix")
                                , DT::DTOutput(ns("design_matrix_preview"))
                                , shiny::br()
                                , shiny::h4("Defined Contrasts")
                                , DT::DTOutput(ns("contrasts_preview"))
                                , shiny::br()
                                , shiny::h4("Assays Included")
                                , shiny::verbatimTextOutput(ns("assays_preview"))
                            )
                        )
                    )

                    # Conditional panel: message if no data
                    , shiny::conditionalPanel(
                        condition = paste0("!output['", ns("data_available"), "']")
                        , shiny::div(
                            class = "alert alert-info"
                            , shiny::icon("info-circle")
                            , " Please complete the 'Import' step first. The builder will appear once data is available."
                        )
                    )
                )
            )
        )
    )
}

#' @rdname metabolomicsDesignMatrixAppletModule
#' @export
#' @importFrom shiny moduleServer reactive observeEvent req renderUI showNotification removeNotification outputOptions renderText
#' @importFrom logger log_info log_error log_warn
#' @importFrom utils write.table
#' @importFrom vroom vroom
#' @importFrom shinyFiles shinyDirButton shinyDirChoose parseDirPath getVolumes
#' @importFrom jsonlite write_json read_json
registerMetabDesignPreviewOutputs <- function(
    output,
    workflowData,
    renderDt = DT::renderDT,
    renderText = shiny::renderText,
    req = shiny::req
) {
    output$design_matrix_preview <- renderDt({
        req(workflowData$design_matrix)
        workflowData$design_matrix
    }, options = list(pageLength = 5, scrollX = TRUE))

    output$contrasts_preview <- renderDt({
        req(workflowData$contrasts_tbl)
        workflowData$contrasts_tbl
    }, options = list(pageLength = 5, scrollX = TRUE))

    output$assays_preview <- renderText({
        req(workflowData$data_tbl)
        assay_names <- names(workflowData$data_tbl)
        paste("Included assays:", paste(assay_names, collapse = ", "))
    })

    invisible(output)
}

registerMetabDesignStateOutputs <- function(
    output,
    workflowData,
    reactive = shiny::reactive,
    outputOptions = shiny::outputOptions
) {
    output$data_available <- reactive({
        !is.null(workflowData$data_tbl) && !is.null(workflowData$config_list)
    })
    outputOptions(output, "data_available", suspendWhenHidden = FALSE)

    output$design_matrix_exists <- reactive({
        !is.null(workflowData$design_matrix)
    })
    outputOptions(output, "design_matrix_exists", suspendWhenHidden = FALSE)

    invisible(output)
}

registerMetabDesignBuilderModule <- function(
    workflowData,
    moduleId = "builder",
    builderServerExists = exists("mod_metab_design_builder_server"),
    builderServerFn = NULL,
    reactiveFn = shiny::reactive,
    reactiveValFn = shiny::reactiveVal
) {
    if (isTRUE(builderServerExists)) {
        if (is.null(builderServerFn)) {
            builderServerFn <- mod_metab_design_builder_server
        }

        return(builderServerFn(
            moduleId,
            data_tbl = reactiveFn(workflowData$data_tbl),
            config_list = reactiveFn(workflowData$config_list),
            column_mapping = reactiveFn(workflowData$column_mapping),
            existing_design_matrix = reactiveFn(workflowData$design_matrix),
            existing_contrasts = reactiveFn(workflowData$contrasts_tbl)
        ))
    }

    reactiveValFn(NULL)
}

registerMetabDesignBuilderResultsObserver <- function(
    builderResultsRv,
    workflowData,
    experimentPaths,
    qcTrigger = NULL,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    runBuilderResultsFlow = NULL
) {
    if (is.null(runBuilderResultsFlow)) {
        runBuilderResultsFlow <- function(results, workflowData, experimentPaths, qcTrigger) {
            shiny::showModal(shiny::modalDialog(
                title = "Processing Design Matrix"
                , shiny::div(
                    style = "text-align: center; padding: 20px;"
                    , shiny::icon("spinner", class = "fa-spin fa-3x")
                    , shiny::br()
                    , shiny::br()
                    , shiny::p("Saving design matrix and preparing data...")
                )
                , footer = NULL
                , easyClose = FALSE
            ))

            logger::log_info("Received results from metabolomics design builder. Saving to workflow and disk.")

            workflowData$design_matrix <- results$design_matrix
            workflowData$data_cln <- results$data_cln
            workflowData$contrasts_tbl <- results$contrasts_tbl
            workflowData$config_list <- results$config_list

            if (!is.null(results$contrasts_tbl)) {
                assign("contrasts_tbl", results$contrasts_tbl, envir = .GlobalEnv)
                logger::log_info("Updated contrasts_tbl in global environment.")
            }

            assign("config_list", workflowData$config_list, envir = .GlobalEnv)
            logger::log_info("Updated global config_list.")

            source_dir <- experimentPaths$source_dir
            if (is.null(source_dir) || !dir.exists(source_dir)) {
                msg <- "Could not find source directory to save files."
                logger::log_error(msg)
                shiny::showNotification(msg, type = "error", duration = 15)
                shiny::removeModal()
                return()
            }

            tryCatch({
                design_matrix_path <- file.path(source_dir, "design_matrix.tab")
                logger::log_info(paste("Writing design matrix to:", design_matrix_path))
                utils::write.table(results$design_matrix, file = design_matrix_path
                    , sep = "\t", row.names = FALSE, quote = FALSE)

                if (!is.null(results$contrasts_tbl) && nrow(results$contrasts_tbl) > 0) {
                    contrast_path <- file.path(source_dir, "contrast_strings.tab")
                    logger::log_info(paste("Writing contrasts to:", contrast_path))
                    writeLines(results$contrasts_tbl$contrasts, contrast_path)
                }

                assay_names <- names(results$data_cln)
                for (assay_name in assay_names) {
                    assay_path <- file.path(source_dir, paste0("data_cln_", assay_name, ".tab"))
                    logger::log_info(paste("Writing assay data to:", assay_path))
                    utils::write.table(results$data_cln[[assay_name]], file = assay_path
                        , sep = "\t", row.names = FALSE, quote = FALSE)
                }

                manifest_path <- file.path(source_dir, "assay_manifest.txt")
                writeLines(assay_names, manifest_path)
                logger::log_info(sprintf("Saved assay manifest with %d assays: %s",
                    length(assay_names), paste(assay_names, collapse = ", ")))

                col_map <- workflowData$column_mapping
                if (!is.null(col_map)) {
                    col_map_path <- file.path(source_dir, "column_mapping.json")
                    jsonlite::write_json(col_map, col_map_path, auto_unbox = TRUE)
                    logger::log_info("Saved column_mapping.json")
                }

                manifest_json_path <- file.path(source_dir, "manifest.json")
                manifest_data <- list(
                    data_path = "assay_manifest.txt",
                    design_matrix_path = "design_matrix.tab",
                    contrast_strings_path = "contrast_strings.tab",
                    column_mapping_path = "column_mapping.json"
                )
                jsonlite::write_json(manifest_data, manifest_json_path, auto_unbox = TRUE, pretty = TRUE)
                logger::log_info(paste("Saved manifest.json to:", manifest_json_path))

                if (!is.null(workflowData$config_list)) {
                    config_path <- file.path(source_dir, "config.ini")
                    logger::log_info(paste("Writing config.ini to:", config_path))
                    tryCatch({
                        ini::write.ini(workflowData$config_list, config_path)
                        logger::log_info("Saved config.ini.")
                    }, error = function(e) {
                        logger::log_warn(paste("Could not save config.ini:", e$message))
                    })
                }

                s4_obj <- createMetaboliteAssayData(
                    metabolite_data = results$data_cln
                    , design_matrix = results$design_matrix
                    , metabolite_id_column = col_map$metabolite_id_col
                    , annotation_id_column = if (!is.null(col_map$annotation_col) && !is.na(col_map$annotation_col) && nzchar(col_map$annotation_col)) {
                        col_map$annotation_col
                    } else {
                        NA_character_
                    }
                    , sample_id = "Run"
                    , group_id = "group"
                    , technical_replicate_id = "tech_rep_group"
                    , database_identifier_type = "Unknown"
                    , internal_standard_regex = if (!is.null(col_map$is_pattern) && !is.na(col_map$is_pattern)) {
                        col_map$is_pattern
                    } else {
                        NA_character_
                    }
                    , args = results$config_list
                )

                if (is.null(workflowData$state_manager)) {
                    workflowData$state_manager <- WorkflowState$new("metabolomics")
                }

                logger::log_info("Saving MetaboliteAssayData S4 object to state manager as 'metab_raw_data_s4'")
                workflowData$state_manager$saveState(
                    state_name = "metab_raw_data_s4"
                    , s4_data_object = s4_obj
                    , config_object = results$config_list
                    , description = "Initial MetaboliteAssayData S4 object created after design matrix"
                )

                if (!is.null(qcTrigger)) {
                    qcTrigger(TRUE)
                    logger::log_info("QC trigger set to TRUE")
                }

                updated_status <- workflowData$tab_status
                updated_status$design_matrix <- "complete"
                workflowData$tab_status <- updated_status

                shiny::removeModal()
                shiny::showNotification("Design matrix and contrasts saved successfully!", type = "message")

                logger::log_info(sprintf("Design save complete. Files saved to: %s", source_dir))
                logger::log_info(sprintf("Saved: design_matrix.tab, %d assay files, assay_manifest.txt, column_mapping.json, config.ini",
                    length(assay_names)))
            }, error = function(e) {
                msg <- paste("Error saving design matrix results:", e$message)
                logger::log_error(msg)
                shiny::removeModal()
                shiny::showNotification(msg, type = "error", duration = 15)
            })
        }
    }

    observeEventFn(builderResultsRv(), {
        results <- builderResultsRv()
        reqFn(results)

        runBuilderResultsFlow(
            results = results,
            workflowData = workflowData,
            experimentPaths = experimentPaths,
            qcTrigger = qcTrigger
        )
    }, ignoreNULL = TRUE)
}












mod_metab_design_server <- function(id, workflow_data, experiment_paths, volumes = NULL, qc_trigger = NULL) {
    message(sprintf("--- Entering mod_metab_design_server ---"))
    message(sprintf("   mod_metab_design_server Arg: id = %s", id))

    shiny::moduleServer(id, function(input, output, session) {
        message(sprintf("   mod_metab_design_server Step: Inside moduleServer function"))

        # == Setup shinyFiles =======================================================
        resolved_volumes <- initializeMetabDesignImportBootstrap(
            input = input,
            session = session,
            experimentPaths = experiment_paths,
            volumes = volumes
        )

        # == Modal Logic for Import =================================================
        registerMetabDesignImportModalShell(
            input = input,
            output = output,
            session = session,
            resolvedVolumes = resolved_volumes
        )

        # == Handle Import Confirmation =============================================
        registerMetabDesignImportObserverShell(
            input = input,
            resolvedVolumes = resolved_volumes,
            workflowData = workflow_data,
            experimentPaths = experiment_paths,
            qcTrigger = qc_trigger
        )

        # == Reactivity Checks ======================================================

        registerMetabDesignStateOutputs(
            output = output,
            workflowData = workflow_data
        )

        # == Module Integration =====================================================

        builder_results_rv <- registerMetabDesignBuilderModule(
            workflowData = workflow_data
        )

        # == Handle Builder Results =================================================

        registerMetabDesignBuilderResultsObserver(
            builderResultsRv = builder_results_rv,
            workflowData = workflow_data,
            experimentPaths = experiment_paths,
            qcTrigger = qc_trigger
        )

        # == Previews of Saved Data =================================================

        registerMetabDesignPreviewOutputs(
            output = output,
            workflowData = workflow_data
        )
    })
}
