#' Design Matrix Builder Server Module
#' Design Matrix Builder Server Module
#'
#' @param id Module ID.
#' @param data_tbl A reactive expression that returns the data table.
#' @param config_list A reactive expression that returns the main configuration list.
#' @param column_mapping A reactive expression that returns the column mapping list.
#' 
#' @note **ARCHITECTURAL EXCEPTION**: This module uses reactive expressions instead of
#' `workflow_data` because it is designed as a standalone component callable from
#' R Markdown workflows. This is an intentional deviation from the standard Golem
#' module pattern.
#' 
#' @rdname designMatrixBuilderModule
#' @export
mod_prot_design_builder_server <- function(id, data_tbl, config_list, column_mapping) {
    # KIT-103 W6: back-reference #5, a wrapper over the shared builder server with the proteomics
    # configuration (the `dia_long` adapter); `.runDesignMatrixApplet` calls it with plain reactives.
    mod_omic_design_builder_server(id, data_tbl, config_list, column_mapping,
        config = omicGuiConfig("proteomics"))
}

