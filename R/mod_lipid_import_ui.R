#' Lipidomics Import Module
#'
#' A Shiny module for importing lipidomics data with vendor-format detection
#' and dynamic column mapping.
#'
#' @name mod_lipid_import
#' @param id,workflow_data,experiment_paths,volumes Runtime inputs used by this function; see the usage section for accepted values.
#' @param profile An omic profile (KIT-096, e.g. the transcriptomics profile); `NULL`, the default, is lipidomics.
NULL

#' @rdname mod_lipid_import
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br radioButtons selectInput textInput actionButton uiOutput verbatimTextOutput icon tags conditionalPanel helpText div
#' @importFrom shinyjs useShinyjs disable enable
mod_lipid_import_ui <- function(id, profile = NULL) {
  ns <- shiny::NS(id)

  # Check if shinyFiles is available and we are not in test mode
  use_shiny_files <- requireNamespace("shinyFiles", quietly = TRUE) && !is_test_mode()

  if (is.null(profile)) {
    return(buildLipidImportUiShell(
      ns = ns,
      useShinyFiles = use_shiny_files
    ))
  }

  # KIT-096: another omic's profile supplies its sections and header at the shell's seams.
  buildLipidImportUiShell(ns = ns, useShinyFiles = use_shiny_files, buildModulePanel = function(ns, useShinyFiles) {
    buildLipidImportModulePanel(ns, useShinyFiles,
      buildFileImportSection = function(ns, useShinyFiles) buildLipidImportFileImportSection(ns, useShinyFiles, profile = profile),
      buildColumnMappingSection = function(ns) buildLipidImportColumnMappingSection(ns, profile = profile),
      buildProcessFooterSection = function(ns) buildLipidImportProcessFooterSection(ns, profile = profile),
      buildHeader = function(...) shiny::h3(profile$import_title))
  })
}
