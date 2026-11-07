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
#' @param volumes A list of volumes for shinyFiles (optional).
#'
#' @name designMatrixAppletModule
#' @param input,output,session,workflowData,experimentPaths,qcTrigger,initializeImportBootstrap,registerImportModalShell,registerImportConfirmationObserver,registerPreviewOutputs,registerBuilderModule,registerBuilderResultsObserver Runtime inputs used by this function; see the usage section for accepted values.
NULL





















#' @rdname designMatrixAppletModule
#' @export
#' @importFrom shiny NS tagList wellPanel h3 p conditionalPanel div icon tags HTML
#' @importFrom DT DTOutput renderDT
mod_prot_design_ui <- function(id) {
  # KIT-103 W6: a wrapper over the shared design tab with the proteomics configuration (its own
  # builder UI as the `proteomics` plugin, back-reference #4).
  mod_omic_design_ui(id, omicGuiConfig("proteomics"))
}

#' The design tab's handler for the UniProt annotation progress (`completeProtDesignPostCheckpoint`
#' sends `updateUniprotProgress`); the shared tab renders it first for the `prot_hook` commit.
#' @noRd
protDesignUniprotProgressScript <- function() {
    shiny::tags$script(shiny::HTML("
      Shiny.addCustomMessageHandler('updateUniprotProgress', function(message) {
        $('#uniprot_progress_bar').css('width', message.percent + '%');
        $('#uniprot_progress_bar').text(message.percent + '%');
        $('#uniprot_progress_text').text(message.text);
      });
    "))
}

#' The proteomics plugins of the shared design tab's server (`registerOmicsDesignServer()`): its
#' own import panel (FASTA and organism sidecars; `runProtDesignImportConfirmationFlow`), the
#' builder over the DIA table's provider (`registerProtDesignBuilderModule`, the `dia_long`
#' adapter), its own commit (`registerProtDesignBuilderResultsObserver` ->
#' `runProtDesignBuilderSaveFlow`, the UniProt post-step and the DIA reclamation) and the previews.
#' The arguments are `registerProtDesignServerShells`' injectable steps.
#' @return A list of the frame's plugin functions (`import`, `builder`, `commit`, `previews`).
#' @noRd
protDesignServerPlugins <- function(
    initializeImportBootstrap = initializeProtDesignImportBootstrap,
    registerImportModalShell = registerProtDesignImportModalShell,
    registerImportConfirmationObserver = registerProtDesignImportConfirmationObserver,
    registerPreviewOutputs = registerProtDesignPreviewOutputs,
    registerBuilderModule = registerProtDesignBuilderModule,
    registerBuilderResultsObserver = registerProtDesignBuilderResultsObserver
) {
  list(
    import = function(input, output, session, workflow_data, experiment_paths, volumes, qc_trigger, config) {
      importBootstrap <- initializeImportBootstrap(input = input, session = session,
        experimentPaths = experiment_paths, volumes = volumes)
      registerImportModalShell(input = input, output = output, session = session,
        resolvedVolumes = importBootstrap$resolvedVolumes, importFastaPath = importBootstrap$importFastaPath)
      registerImportConfirmationObserver(input = input, resolvedVolumes = importBootstrap$resolvedVolumes,
        importFastaPath = importBootstrap$importFastaPath, workflowData = workflow_data,
        experimentPaths = experiment_paths, session = session, qcTrigger = qc_trigger)
    },
    builder = function(workflow_data, config) registerBuilderModule(workflowData = workflow_data),
    commit = function(builder_results, workflow_data, experiment_paths, session, qc_trigger, config) {
      registerBuilderResultsObserver(builderResultsRv = builder_results, workflowData = workflow_data,
        experimentPaths = experiment_paths, session = session, qcTrigger = qc_trigger)
    },
    previews = function(output, workflow_data, config) {
      registerPreviewOutputs(output = output, workflowData = workflow_data)
    })
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
  # KIT-103 W6: the shared design tab's server under the proteomics configuration, its plugins
  # built from these steps.
  registerOmicsDesignServer(input, output, session, workflowData, experimentPaths, volumes, qcTrigger,
    omicGuiConfig("proteomics"), plugins = protDesignServerPlugins(initializeImportBootstrap,
      registerImportModalShell, registerImportConfirmationObserver, registerPreviewOutputs,
      registerBuilderModule, registerBuilderResultsObserver))
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
