#' Lipidomics Design Matrix Applet Module
#'
#' A Shiny module that hosts the lipidomics design-matrix workflow, including
#' the design builder, existing-design import, S4 object creation, and state
#' persistence.
#'
#' @param id Module ID.
#' @param workflow_data Reactive values used to store workflow data.
#' @param experiment_paths Paths for the current experiment.
#' @param volumes Optional shinyFiles volumes.
#' @param qc_trigger Optional reactive QC trigger.
#' @param profile An omic profile (KIT-096): its help text, the builder's run order and level order, and the design commit's
#'   collaborators. `NULL`, the default, is lipidomics.
#'
#' @name lipidomicsDesignMatrixAppletModule
NULL

#' @rdname lipidomicsDesignMatrixAppletModule
#' @export
#' @importFrom shiny NS tagList wellPanel h3 h4 p conditionalPanel div icon tags HTML fluidRow column actionButton
#' @importFrom DT DTOutput renderDT
mod_lipid_design_ui <- function(id, profile = NULL) {
    # KIT-103 W4: the shared design tab (mod_omic_design.R) under the lipid configuration, or the
    # profile's omic (trx) under its own.
    mod_omic_design_ui(id, lipidDesignGuiConfig(profile))
}
