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
    runProtDesignBuilderServerEntryShell(
        id = id,
        dataTbl = data_tbl,
        configList = config_list,
        columnMapping = column_mapping
    )
}

