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

#' Proteomics Workflow UI Module
#'
#' @description Main UI for the proteomics workflow interface
#'
#' @param id Module ID
#'
#' @importFrom shiny NS tagList fluidRow column wellPanel h4 uiOutput tabsetPanel
#' @importFrom shiny tabPanel icon br conditionalPanel div hr
#' @importFrom shiny actionButton sliderInput numericInput plotOutput radioButtons
#' @importFrom shiny textAreaInput checkboxGroupInput textInput downloadButton
#' @importFrom shiny verbatimTextOutput h3 p
#' @importFrom DT DTOutput
#' @export
mod_proteomics_ui <- function(id) {
  # KIT-102: the shared shell (mod_omic_workflow.R) under the proteomics configuration.
  mod_omic_workflow_ui(id, omicGuiConfig("proteomics"))
}

#' Proteomics Workflow Server
#'
#' @description Server logic for the proteomics workflow interface
#'
#' @param id Module ID
#' @param project_dirs Project directory structure from setupDirectories
#' @param omic_type The omics type (should be "proteomics")
#' @param experiment_label The experiment label for this analysis
#' @param volumes Volumes for shinyFiles
#'
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req reactiveVal
#' @importFrom logger log_info log_error
#' @export
mod_proteomics_server <- function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
  # KIT-102: the shared shell (mod_omic_workflow.R) under the proteomics configuration.
  mod_omic_workflow_server(id, project_dirs, omic_type, experiment_label, volumes, omicGuiConfig("proteomics"))
}

# ------------------------------------------------------------------------------------------------
# The shell's proteomics extras (component `peptide_checkpoint`, omicGuiComponents): the DIA peptide-checkpoint panel
# with its save and open, the native-restore QC trigger and the aa_seq_tbl_final resume, moved unchanged from the
# pre-KIT-102 mod_proteomics_ui/server. They live in the shell's namespace.
# ------------------------------------------------------------------------------------------------

#' The peptide-checkpoint panel.
#' @param ns The shell's namespace function.
#' @noRd
protPeptideCheckpointUi <- function(ns) {
  shiny::wellPanel(
    shiny::h4("DIA peptide checkpoints"),
    shiny::actionButton(ns("save_peptide_checkpoint"), "Save peptide checkpoint"),
    if (requireNamespace("shinyFiles", quietly = TRUE) && !is_test_mode()) {
      shinyFiles::shinyFilesButton(ns("peptide_checkpoint_file"), "Choose checkpoint",
        "Choose a peptide checkpoint", multiple = FALSE)
    } else {
      shiny::fileInput(ns("peptide_checkpoint_upload"), "Checkpoint archive", accept = ".msdia.tar")
    },
    shiny::actionButton(ns("open_peptide_checkpoint"), "Open peptide checkpoint"),
    shiny::helpText("Save the selected peptide stage and its history. Open a checkpoint to resume peptide QC."),
    shiny::verbatimTextOutput(ns("peptide_checkpoint_status"))
  )
}

#' The extras' server: the native-restore QC trigger, the checkpoint save and open, and the aa_seq_tbl_final resume.
#' @param input,output,session The shell's.
#' @param context The shell's context: `workflow_data`, `experiment_paths` and `qc_trigger`.
#' @noRd
protPeptideCheckpointServer <- function(input, output, session, context) {
  workflow_data <- context$workflow_data
  experiment_paths <- context$experiment_paths
  qc_trigger <- context$qc_trigger
  shiny::observeEvent(workflow_data$native_session_restored, {
    qc_trigger(TRUE)
  }, ignoreNULL = TRUE)

  checkpoint_status <- shiny::reactiveVal("")
  output$peptide_checkpoint_status <- shiny::renderText(checkpoint_status())
  checkpoint_volumes <- NULL
  if (requireNamespace("shinyFiles", quietly = TRUE) && !is_test_mode()) {
    checkpoint_volumes <- shinyFiles::getVolumes()
    shinyFiles::shinyFileChoose(input, "peptide_checkpoint_file", roots = checkpoint_volumes,
      session = session, filetypes = "tar")
  }
  shiny::observeEvent(input$save_peptide_checkpoint, {
    tryCatch({
      published <- saveDIAPeptideCheckpoint(workflow_data, session, experiment_paths)
      checkpoint_status(paste("Checkpoint saved:", published$archive_path))
    }, error = function(error) {
      checkpoint_status(paste("Checkpoint save failed:", conditionMessage(error)))
      logger::log_error(paste("Peptide checkpoint save failed:", conditionMessage(error)))
    })
  })
  shiny::observeEvent(input$open_peptide_checkpoint, {
    tryCatch({
      path <- NULL
      if (!is.null(checkpoint_volumes) && !is.null(input$peptide_checkpoint_file)) {
        selected <- shinyFiles::parseFilePaths(checkpoint_volumes, input$peptide_checkpoint_file)
        if (nrow(selected) == 1L) path <- as.character(selected$datapath[[1L]])
      } else path <- input$peptide_checkpoint_upload$datapath
      if (is.null(path) || !file.exists(path)) stop("choose a peptide checkpoint first")
      state <- openDIAPeptideCheckpoint(path, workflow_data, session, experiment_paths)
      status <- workflow_data$tab_status
      status$setup_import <- "complete"
      status$design_matrix <- "complete"
      status$quality_control <- "in_progress"
      workflow_data$tab_status <- status
      shiny::updateTabsetPanel(session, "workflow_tabs", selected = "qc")
      checkpoint_status(paste("Checkpoint opened:", state))
    }, error = function(error) {
      checkpoint_status(paste("Checkpoint open failed:", conditionMessage(error)))
      logger::log_error(paste("Peptide checkpoint open failed:", conditionMessage(error)))
    })
  })
  # The on-disk DIA state manager is attached at import (initializeDIANativeImport, KIT-011); the
  # capacity-trial session of the PERF era is gone.

  # Load aa_seq_tbl_final from scripts directory if resuming session
  if (!is.null(experiment_paths) && !is.null(experiment_paths$source_dir)) {
    aa_seq_file_path <- file.path(experiment_paths$source_dir, "aa_seq_tbl_final.RDS")
    if (file.exists(aa_seq_file_path)) {
      log_info("Loading existing aa_seq_tbl_final from scripts directory for session resumption")
      tryCatch(
        {
          aa_seq_tbl_final <- readRDS(aa_seq_file_path)
          workflow_data$aa_seq_tbl_final <- aa_seq_tbl_final
          log_info(sprintf("Successfully loaded aa_seq_tbl_final with %d sequences", nrow(aa_seq_tbl_final)))
        },
        error = function(e) {
          log_warn(paste("Error loading aa_seq_tbl_final:", e$message))
        }
      )
    }
  }
  invisible(NULL)
}
