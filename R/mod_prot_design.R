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

#' @title designMatrixAppletModule
#' 
#' @description A Shiny module that serves as the main host for the design
#' matrix creation workflow step. It embeds the design matrix builder UI
#' and handles the logic for saving the results back to the main workflow.
#' 
#' @param id Module ID
#' @param workflow_data A reactive values object to store workflow data.
#' @param experiment_paths A list of paths for the current experiment.
#' @param volumes A list of volumes for shinyFiles (optional).
#' @param qc_trigger A reactive trigger for QC execution (optional).
#'
#' @name designMatrixAppletModule
NULL





















#' @rdname designMatrixAppletModule
#' @export
#' @importFrom shiny NS tagList wellPanel h3 p conditionalPanel div icon tags HTML
#' @importFrom DT DTOutput renderDT
mod_prot_design_ui <- function(id) {
  ns <- shiny::NS(id)
  
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
          # Using mod_prot_design_builder_ui
          if (exists("mod_prot_design_builder_ui")) {
            mod_prot_design_builder_ui(ns("builder"))
          } else {
            shiny::div("Design builder module not loaded")
          },
          
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
#' @importFrom shiny moduleServer reactive observeEvent req renderUI showNotification removeNotification outputOptions
#' @importFrom logger log_info log_error
#' @importFrom utils write.table
#' @importFrom vroom vroom
#' @importFrom shinyFiles shinyDirButton shinyDirChoose parseDirPath getVolumes
#' @importFrom tidyr pivot_wider
#' @importFrom rlang sym
registerProtDesignServerShells <- function(
    input,
    output,
    session,
    workflowData,
    experimentPaths,
    volumes = NULL,
    qcTrigger = NULL,
    initializeImportBootstrap = initializeProtDesignImportBootstrap,
    registerImportModalShell = registerProtDesignImportModalShell,
    registerImportConfirmationObserver = registerProtDesignImportConfirmationObserver,
    registerPreviewOutputs = registerProtDesignPreviewOutputs,
    registerBuilderModule = registerProtDesignBuilderModule,
    registerBuilderResultsObserver = registerProtDesignBuilderResultsObserver
) {
  importBootstrap <- initializeImportBootstrap(
    input = input,
    session = session,
    experimentPaths = experimentPaths,
    volumes = volumes
  )
  resolvedVolumes <- importBootstrap$resolvedVolumes
  importFastaPath <- importBootstrap$importFastaPath

  registerImportModalShell(
    input = input,
    output = output,
    session = session,
    resolvedVolumes = resolvedVolumes,
    importFastaPath = importFastaPath
  )

  registerImportConfirmationObserver(
    input = input,
    resolvedVolumes = resolvedVolumes,
    importFastaPath = importFastaPath,
    workflowData = workflowData,
    experimentPaths = experimentPaths,
    session = session,
    qcTrigger = qcTrigger
  )

  registerPreviewOutputs(
    output = output,
    workflowData = workflowData
  )

  builderResultsRv <- registerBuilderModule(
    workflowData = workflowData
  )

  registerBuilderResultsObserver(
    builderResultsRv = builderResultsRv,
    workflowData = workflowData,
    experimentPaths = experimentPaths,
    session = session,
    qcTrigger = qcTrigger
  )

  invisible(builderResultsRv)
}

mod_prot_design_server <- function(id, workflow_data, experiment_paths, volumes = NULL, qc_trigger = NULL) {
  message(sprintf("--- Entering mod_prot_design_server ---"))
  message(sprintf("   mod_prot_design_server Arg: id = %s", id))
  
  shiny::moduleServer(id, function(input, output, session) {
    message(sprintf("   mod_prot_design_server Step: Inside moduleServer function"))

    registerProtDesignServerShells(
      input = input,
      output = output,
      session = session,
      workflowData = workflow_data,
      experimentPaths = experiment_paths,
      volumes = volumes,
      qcTrigger = qc_trigger
    )
  })
}
