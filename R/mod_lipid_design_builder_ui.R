#' Lipidomics Design Matrix Builder Module
#'
#' A Shiny module to build a design matrix for lipidomics experiments. It
#' follows the proteomics design-builder UI and server modules, with the key
#' difference that lipidomics data is stored as a list of assays where samples
#' are column names rather than row values.
#'
#' @param id The module ID.
#' @param data_tbl A reactive expression returning the data-table list.
#' @param config_list A reactive expression returning the configuration list.
#' @param column_mapping A reactive expression returning the column mapping.
#'
#' @return A reactive expression containing the saved design matrix, filtered
#'   data, contrasts, and updated configuration.
#'
#' @name lipidomicsDesignMatrixBuilderModule
NULL

#' Lipidomics Design Matrix Builder UI Module
#'
#' @param id Module ID.
#' @return A fluidPage containing the UI for the design matrix builder.
#'
#' @note **ARCHITECTURAL EXCEPTION**: This module is designed as a standalone component
#' callable from R Markdown workflows, not integrated into the main Shiny app workflow.
#' Therefore, it returns `fluidPage()` instead of `tagList()`, and uses reactive expressions
#' instead of `workflow_data` for state management. This is an intentional deviation from
#' the standard Golem module pattern.
#'
#' @rdname lipidomicsDesignMatrixBuilderModule
#' @export
mod_lipid_design_builder_ui <- function(id) {
    # KIT-103 W4: the shared builder UI (mod_omic_design_builder.R) with the lipid texts.
    mod_omic_design_builder_ui(id, omicGuiConfig("lipidomics"))
}
